/*                                                                
**  Copyright (C) 2005,2007,2011  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */


#include <dslib.h>
#include <dsnan.h>
#include <histlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <stdlib.h>

#include "resample_img.h"

/* 
   These typedef's define a) an enumerated type to be used to know which trans
   forms to apply to do the reprojection; and b) to create a structure that
   will hold the Physical and World WCS dmDescriptors for the X and Y axes;
   for both the inpt and output images.
*/
typedef enum { coordLOGICAL, coordPHYSICAL, coordWORLD } CoordType;
typedef struct { dmDescriptor *xaxis; dmDescriptor *yaxis; } Axes;
typedef struct { Axes physical; Axes world; } WCS;
typedef struct { WCS ref; WCS image; } WCS_Descriptors;


/* Hold info for an input image */
typedef struct {
  void *data;        // pixel values
  dmDataType dt;     // pixel datatype
  long *lAxes;       // axis lenghts
  short *mask;        // mask of valid pixels
  dmDescriptor *xdesc;  // X (or sky) coordinate descriptor
  dmDescriptor *ydesc;  // Y coordinate descriptor
  dmBlock *block; // The block image came from
} Image;

/* Input parameters */
typedef struct {
  char instack[DS_SZ_PATHNAME];  /* Input stack of images */
  char reffile[DS_SZ_PATHNAME];  /* match file */
  char outfile[DS_SZ_PATHNAME];  /* output file */
  short clobber;
  short verbose;
  int subpix;
  char lookup[DS_SZ_PATHNAME];
  char csys[30];
  CoordType ctype;
  long quantum;
  long randseed;

} Parameters;


/* Buffer to hold areas */
typedef struct {
  long maxlen;
  long len;
  long *xx;
  long *yy;
  double *area;      
} Buffer;

typedef struct {
  Polygon *ref;
  Polygon *tmp;
  Polygon *clip;
} Polygons;



double get_contour_area( Polygons *poly, long xx, long yy );
int fill_polygon( long subpix, long xx, long yy, VertexList *refpixlist, WCS_Descriptors *descs, CoordType ctype);
int convert_coords( double *refpix, WCS_Descriptors *descs, double *imgpix, CoordType ctype);
int check_coords( CoordType ctype, WCS_Descriptors descs );
Image *load_image_file( dmBlock *inBlock );
Polygon *make_polygon(int subpix);
int find_bounding_box( Polygon *ref_poly, long *lAxes, long *xx_min, long *xx_max, long *yy_min, long*yy_max) ;
Parameters* get_parameters(void);
Image *load_ref_image( Parameters *pars, WCS_Descriptors *descs);
Image *load_infile_image(Parameters *pars, char *infile, Header_Type **hdr, WCS_Descriptors *descs  );
void free_image( Image *img );
Buffer *make_buffer(void);
void add_to_buffer( Buffer *buffer, long xx, long yy, double area );
void cummulative( Buffer *buffer );
long sample_distro( Buffer *buffer, Image *refImage);
int is_invalid_val( double val );
void distribute_counts( Parameters *pars, Buffer *buffer, Image *refImage, double val, double *out_data);
Polygons *make_polys( int subpix );
int process_infile( Image *inImage, Image *refImage, WCS_Descriptors *descs, Parameters *pars, Buffer *buffer, 
    Polygons *polys, double *out_data );
Header_Type *merge_headers( Parameters *pars, long num_infiles, Header_Type **hdr );
int write_output( Parameters *pars, Image *refImage, double *out_data, long num_infiles, Header_Type **hdr );


// --------------------------



/* Load image file.  This uses the routines defined in dmimgio.h */
#include <dmimgio.h>
#include <cxcregion.h>

Image *load_image_file( dmBlock *inBlock )
{
  dmDataType dt;
  void *data=NULL;
  long *lAxes=NULL;

  regRegion *dss=NULL;
  long null;
  short has_null;

  short *mask=NULL;
  dmDescriptor *xdesc=NULL;
  dmDescriptor *ydesc=NULL;

  /* Read input */

  dt = get_image_data( inBlock, &data, &lAxes, &dss, &null, &has_null );
  get_image_wcs( inBlock, &xdesc, &ydesc );
  mask = get_image_mask( inBlock, data, dt, lAxes, dss, null, has_null, 
                         xdesc, ydesc );

  Image *img = (Image*)calloc( 1,sizeof(Image));
  img->dt = dt;
  img->data = data;
  img->lAxes = lAxes;
  img->mask = mask;
  img->xdesc = xdesc;
  img->ydesc = ydesc;
  img->block = inBlock;

  return(img);
}


void free_image( Image *img )
{
    if ( NULL == img ) return;
    
    if ( NULL != img->data ) { 
        free( img->data );
        img->data = NULL;
    }
    if ( NULL != img->mask ) {
        free( img->mask );
        img->mask = NULL;        
    }
    if ( NULL != img->lAxes ) {
        free(img->lAxes);
        img->lAxes = NULL;
    }
    return;
}



/* Routine to convert coordinates */

/*
  Note:  There are cases where the dmCoordCalc/Invert will fail.
    We want to plow through these since the error returned may not be 
    fatal to what we want to do; so wait till end to return
    the error code. (Only caught in 1 place below).
*/
/*
 *  ** This code is different than reproject_image.
 * This code does coord calcs from image to ref.  The reproject_image
 * code is *intentionally* the other was around.
 */
int convert_coords( double *refpix,
                    WCS_Descriptors *descs,
                    double *imgpix,
                    CoordType ctype)
{

  int retval = 0;

  double refphys[2];
  double refwrld[2];
  double imgphys[2];
  switch (ctype) {

  case coordWORLD: {
    if ( dmSUCCESS != dmCoordCalc_d( descs->image.physical.xaxis, refpix, refphys ) ) retval=-1;
    if ( dmSUCCESS != dmCoordCalc_d( descs->image.world.xaxis, refphys, refwrld )) retval=-1;
    if ( dmSUCCESS != dmCoordInvert_d( descs->ref.world.xaxis, refwrld, imgphys )) retval=-1;
    if ( dmSUCCESS != dmCoordInvert_d( descs->ref.physical.xaxis, imgphys, imgpix )) retval=-1;
    
    if ( descs->image.physical.yaxis ) {
      if ( dmSUCCESS != dmCoordCalc_d( descs->image.physical.yaxis, refpix+1, refphys+1 )) retval=-1;
      if ( dmSUCCESS != dmCoordCalc_d( descs->image.world.yaxis, refphys+1, refwrld+1 )) retval=-1;
      if ( dmSUCCESS != dmCoordInvert_d(  descs->ref.world.yaxis, refwrld+1, imgphys+1 )) retval=-1;
      if ( dmSUCCESS != dmCoordInvert_d( descs->ref.physical.yaxis, imgphys+1, imgpix+1 )) retval=-1;
    }
    break;
  }

  case coordLOGICAL: {
    imgpix[0] = refpix[0];
    imgpix[1] = refpix[1];
    break;
  }

  case coordPHYSICAL: {
    if ( dmSUCCESS != dmCoordCalc_d( descs->image.physical.xaxis, refpix, refphys )) retval=-1;
    if ( dmSUCCESS != dmCoordInvert_d( descs->ref.physical.xaxis, refphys, imgpix )) retval=-1;
    if (  descs->image.physical.yaxis ) {
      if ( dmSUCCESS != dmCoordCalc_d( descs->image.physical.yaxis, refpix+1, refphys+1 )) retval=-1;
      if ( dmSUCCESS != dmCoordInvert_d( descs->ref.physical.yaxis, refphys+1, imgpix+1 )) retval=-1;
    }
    break;
  }
  } /* end switch */


  return(retval);
}
                    

/*
   Okay, now we are down to the real reproject_image* code.

  This routine is used to make a polygon around the square output pixel.  

  The points of this  polygon will then get mapped to the input image
  and the output pixel value will be constructed from those overlapping
  input pixels weighted by the area of the polygon that covers them.
*/

int fill_polygon( long subpix,
                  long xx, long yy,
                  VertexList *refpixlist,
                  WCS_Descriptors *descs,
                  CoordType ctype)
{

  long refidx;
  double refpix[2];
  double imgpix[2];
  long ii;
  
  double delta = 1.0/subpix;

  
  refidx = 0;
  refpix[0] = xx + 0.5; /* = + 1 - 0.5 */
  refpix[1] = yy + 0.5; /* = + 1 - 0.5 */
  for (ii=0;ii<subpix;ii++) { /* order matter */
    convert_coords( refpix, descs, imgpix,ctype);
    imgpix[0] -= 1;
    imgpix[1] -= 1;
    refpixlist->vertex[refidx].x=imgpix[0];
    refpixlist->vertex[refidx].y=imgpix[1];
    refidx++;
    refpix[1] += delta;
  }

  refpix[0] = xx + 0.5; /* = + 1 - 0.5 */
  refpix[1] = yy + 1.5; /* = + 1 - 0.5 */
  for (ii=0;ii<subpix;ii++) { /* order matter */
    convert_coords( refpix, descs, imgpix,ctype);
    imgpix[0] -= 1;
    imgpix[1] -= 1;
    refpixlist->vertex[refidx].x=imgpix[0];
    refpixlist->vertex[refidx].y=imgpix[1];
    refidx++;
    refpix[0] += delta;
  }
  
  refpix[0] = xx + 1.5; /* = + 1 - 0.5 */
  refpix[1] = yy + 1.5; /* = + 1 - 0.5 */
  for (ii=subpix;ii--;) { /* order matter */
    convert_coords( refpix, descs, imgpix,ctype);
    imgpix[0] -= 1;
    imgpix[1] -= 1;
    refpixlist->vertex[refidx].x=imgpix[0];
    refpixlist->vertex[refidx].y=imgpix[1];
    refidx++;
    refpix[1] -= delta;
  }
  
  refpix[0] = xx + 1.5; /* = + 1 + 0.5 */
  refpix[1] = yy + 0.5; /* = + 1 - 0.5 */
  for (ii=subpix;ii--;) { /* order matter */
    convert_coords( refpix, descs, imgpix,ctype);
    imgpix[0] -= 1;
    imgpix[1] -= 1;
    refpixlist->vertex[refidx].x=imgpix[0];
    refpixlist->vertex[refidx].y=imgpix[1];
    refidx++;
    refpix[0] -= delta;
  }
  

  return(0);
}


/*
  After the polygon has been clipped; we need to compute the area ==
  0.5 * sum ( x_i*y_(i+1) - x_(i+1)*y(i) )

*/
double get_contour_area( Polygons *polys, long xx, long yy)
{
  double area;
  double larea = 0;
  long zz;

  /* Create a polygon that is clipped by the current pixel */
  super_poly_clip( polys->ref, polys->tmp, polys->clip, xx, yy );

  /* Compute the area of that polygon */
  area = 0;    
  for (zz=0;zz<polys->clip->contour->num_vertices-1;zz++) {
      larea += (( polys->clip->contour->vertex[zz].x *
                  polys->clip->contour->vertex[zz+1].y ) -
                ( polys->clip->contour->vertex[zz+1].x *
                  polys->clip->contour->vertex[zz].y  ));
  }
    /* last point */
  larea += (( polys->clip->contour->vertex[zz].x *
              polys->clip->contour->vertex[0].y ) -
            ( polys->clip->contour->vertex[0].x *
              polys->clip->contour->vertex[zz].y ) );
    
  larea /= 2.0;
    
  area += fabs(larea);

  return(area);
}



int check_coords( CoordType ctype, WCS_Descriptors descs )
{
    /* Some error checking on the WORLD coord calcs */
    if ( coordWORLD == ctype ) {
      if ((NULL==descs.image.world.xaxis)||(NULL==descs.ref.world.xaxis)) {
        err_msg("ERROR: There must be a WCS on both input "
                "images if coord_sys=world\n");
        return(-1);
      } else {
        char iname[50];
        char rname[50];
        long npar;
        
        dmCoordGetTransformType( descs.image.world.xaxis, iname, &npar, 49);
        dmCoordGetTransformType( descs.ref.world.xaxis, rname, &npar, 49 ); 

        /* Not sure if this is really the check to be doing or
           if you should check some other names */
        if ( 0 != ds_strcmp_cis( iname, rname)) {
          err_msg("WARNING: The WCS on the image images should be the same "
                  "type when using coord_sys=wcs (ie TAN-P)");
        }                
      } /* end else coords */      
    } /* end if coord */

    return(0);
}


/* Alloc memory for a polygon struct */
Polygon *make_polygon(int subpix)
{
    VertexList *refpixlist = (VertexList*)calloc(1,sizeof(VertexList));   
    Polygon *ref_poly = (Polygon*)calloc(1,sizeof(Polygon));
        

    refpixlist->vertex = (Vertex*)calloc(4*subpix*4, sizeof(Vertex)); // why 4*4? 4 sides
    refpixlist->num_vertices = 4*subpix;    
    ref_poly->contour = refpixlist;
    return(ref_poly);
}


/* find the min/max x,y values from the polygon, must be w/i lAxes ranges */
int find_bounding_box( Polygon *ref_poly, long *lAxes, long *xx_min, long *xx_max, long *yy_min, long*yy_max) 
{
    long ii;
 
    *xx_min = ref_poly->contour->vertex[0].x-0.5;
    *yy_min = ref_poly->contour->vertex[0].y-0.5;
    *xx_max = ref_poly->contour->vertex[0].x+0.5;
    *yy_max = ref_poly->contour->vertex[0].y+0.5;
    for(ii=ref_poly->contour->num_vertices;--ii;) {
      *xx_min = MIN( *xx_min, (ref_poly->contour->vertex[ii].x-0.5));
      *yy_min = MIN( *yy_min, (ref_poly->contour->vertex[ii].y-0.5));
      *xx_max = MAX( *xx_max, (ref_poly->contour->vertex[ii].x+0.5));
      *yy_max = MAX( *yy_max, (ref_poly->contour->vertex[ii].y+0.5));
    }

    *xx_min = MIN( MAX( *xx_min, 0), lAxes[0]-1);
    *xx_max = MIN( MAX( *xx_max, 0), lAxes[0]-1);
    *yy_min = MIN( MAX( *yy_min, 0), lAxes[1]-1);
    *yy_max = MIN( MAX( *yy_max, 0), lAxes[1]-1);

    return(0);
}

/* Get parameters from .par file 
 * 
 * Also does some elementary checks and clobbering
 * */

Parameters *get_parameters(void)
{
  Parameters *pars = (Parameters*)calloc(1,sizeof(Parameters));
    
  /* Get the parameters */
  clgetstr( "infile", pars->instack, DS_SZ_FNAME );
  clgetstr( "matchfile", pars->reffile, DS_SZ_FNAME );
  clgetstr( "outfile", pars->outfile, DS_SZ_FNAME );
  pars->subpix = clgeti("resolution");
  pars->quantum = clgeti("quantum");
  clgetstr("coord_sys", pars->csys, 29);
  pars->randseed = clgeti("randseed");
  clgetstr( "lookupTab", pars->lookup, DS_SZ_FNAME);
  pars->clobber = clgetb( "clobber" );
  pars->verbose = clgeti( "verbose" );
  if ( pars->verbose ) {  
    printf("resample_image - parameters\n");
    printf("%15s = %-s\n", "infile", pars->instack );
    printf("%15s = %-s\n", "matchfile", pars->reffile );
    printf("%15s = %-s\n", "outfile", pars->outfile );
    printf("%15s = %d\n", "resolution", pars->subpix );
    printf("%15s = %ld\n", "quantum", pars->quantum );
    printf("%15s = %-s\n", "coord_sys", pars->csys );
    printf("%15s = %ld\n", "randseed", pars->randseed );
    printf("%15s = %-s\n", "lookupTab", pars->lookup );
    printf("%15s = %-s\n", "clobber", (pars->clobber ? "yes" : "no") );
    printf("%15s = %d\n", "verbose", pars->verbose );
  }

  if ( ( strlen( pars->reffile) == 0 ) ||
       ( ds_strcmp_cis(pars->reffile, "none" ) == 0 ) ) {
    err_msg("ERROR: Must supply a valid match file\n");
    return(NULL);
   }


  switch ( pars->csys[0] ) {
    case 'l': pars->ctype = coordLOGICAL; break;
    case 'p': pars->ctype = coordPHYSICAL; break;
    case 'w': pars->ctype = coordWORLD; break;
    default:
      err_msg("ERROR: Unknow coordinate type '%s'\n", pars->csys );
    return(NULL);
  }

  if ( ds_clobber( pars->outfile, pars->clobber, NULL ) != 0 ) {
    return(NULL);
  }

  Stack inStack;
  inStack = stk_build( pars->instack );
  if ( ( NULL == inStack ) || ( stk_count(inStack)==0 ) ||
       (( stk_count(inStack)==1 ) && ( strlen(stk_read_num(inStack,1))==0))) {
    err_msg("ERROR: problems opening stack '%s'\n", pars->instack );
    return(NULL);
  }    

  /* Set random seed */
  if (( -1 == pars->randseed ) || (INDEFL == pars->randseed)) {
        srand48( time(NULL));
  } else {
        srand48( pars->randseed );
  }

  return(pars);
    
}


/* Load the reference image.  Really only care about the WCS here */
Image *load_ref_image( Parameters *pars, WCS_Descriptors *descs)
{
  Image *refImage;
  dmBlock *blk;

  if ( NULL == ( blk = dmImageOpen( pars->reffile) ) ) {
      err_msg("ERROR: Cannot open image '%s'\n", pars->reffile );
      return(NULL);
  }
  
  if ( NULL == ( refImage = load_image_file( blk ))) {
      err_msg("ERROR: Cannot open image '%s'\n", pars->reffile );
      return(NULL);      
  }

  descs->ref.physical.xaxis = refImage->xdesc;
  descs->ref.physical.yaxis = refImage->ydesc;
  descs->ref.world.xaxis = dmDescriptorGetCoord( refImage->xdesc );
  descs->ref.world.yaxis = dmDescriptorGetCoord( refImage->ydesc );

  return(refImage);
}


/* Load the infile image. */
Image *load_infile_image(Parameters *pars, char *infile, Header_Type **hdr, WCS_Descriptors *descs  )
{
    Image *inImage;
    dmBlock *blk;

    if ( pars->verbose > 1 ) {
      printf("\nProcessing input file '%s'\n", infile );
    }
     
    /* Now load the image */
    if ( NULL == ( blk = dmImageOpen( infile) ) ) {
      err_msg("ERROR: Cannot open image '%s'\n", infile );
      return(NULL);
    }

    if ( NULL == ( inImage = load_image_file( blk ))) {
      err_msg("ERROR: Cannot open image '%s'\n", infile );
      return(NULL);        
    }

    /* The infile must be integer ! */
    if (( dmFLOAT == inImage->dt ) ||
        ( dmDOUBLE == inImage->dt ) ) {
        err_msg("ERROR: The input image must be an integer datatype.");
        return(NULL);
    }
    
    
    *hdr = getHdr( inImage->block, hdrDM_FILE );

    descs->image.physical.xaxis = inImage->xdesc;
    descs->image.physical.yaxis = inImage->ydesc;
    descs->image.world.xaxis = dmDescriptorGetCoord( inImage->xdesc );
    descs->image.world.yaxis = dmDescriptorGetCoord( inImage->ydesc );
    
    if ( 0 != check_coords( pars->ctype, *descs ) ) {
        return(NULL);
    }

    return(inImage);
}

 /* The buffer stuff then is keeping track of which pixels 
 * in the output map to the "current" input pixel and with what
 * area-fraction.
 */
 
Buffer *make_buffer(void)
{
    // Allocate the buffer.  It will grow if needed 

    Buffer *buffer = (Buffer*)calloc(1,sizeof(Buffer));
    buffer->maxlen = 100;
    buffer->len = 0;
    buffer->xx = (long*)calloc(buffer->maxlen, sizeof(long));
    buffer->yy = (long*)calloc(buffer->maxlen, sizeof(long));
    buffer->area = (double*)calloc(buffer->maxlen, sizeof(double));
    
    return(buffer);
    
}

void add_to_buffer( Buffer *buffer, long xx, long yy, double area )
{
    // Add an x, y, area tuple to the buffer.  Increase
    // buffer size if needed.

    buffer->xx[buffer->len] = xx;
    buffer->yy[buffer->len] = yy;
    buffer->area[buffer->len] = area;

    buffer->len += 1;
    if ( buffer->len >= buffer->maxlen ) {
        buffer->maxlen *= 2;
        buffer->xx = (long*)realloc( buffer->xx, sizeof(long)*buffer->maxlen);
        buffer->yy = (long*)realloc( buffer->yy, sizeof(long)*buffer->maxlen);
        buffer->area = (double*)realloc( buffer->area, sizeof(double)*buffer->maxlen);
    }
    
}


void cummulative( Buffer *buffer )
{
    // Turn the area array into a fractional (0 to 1), cummulative
    // probability distribution

    long bb;
    for (bb=1;bb<buffer->len; bb++) 
        buffer->area[bb] += buffer->area[bb-1];
    for (bb=0;bb<buffer->len;bb++) 
        buffer->area[bb] /= buffer->area[buffer->len-1];

}


long sample_distro( Buffer *buffer, Image *refImage)
{
    // Randomly sample the values in the buffer based on the
    // relative area fractions.
    // 
    double randval;
    long bb;
    randval = drand48();  // 0 to 1
    for (bb=0;bb<buffer->len;bb++) {
        // I'm using '<=' here so I don't have to special case "0" or "1"
        if (randval <= buffer->area[bb]) break;
    }
        
    long outpix;
    outpix = buffer->xx[bb] + buffer->yy[bb]*refImage->lAxes[0];
    return(outpix);

}

/* Check for 0 or NULL pixels */
int is_invalid_val( double val )
{
    static int show_less_than_zero_message = 1; // Only show message once

    if ( 0 == val )  return(1); // skip it.
    if ( ds_dNAN(val) ) return(1); // skip it.
    if ( 0 > val ) {
        if ( 1 == show_less_than_zero_message) {
          err_msg("WARNING: Pixel values less than 0 are ignored");
          show_less_than_zero_message=0;
        }
        return(1);
    }

    return(0);
}


void distribute_counts( Parameters *pars, Buffer *buffer, Image *refImage, double val, double *out_data)
{
    long vv;
    long outpix;
    long step = pars->quantum;
    /*
     * What's 'quantum' all about?  Let's say all the pixel
     * values are in the 3000 range.  Well, we may not care to
     * redistrbute all 3000 values individually.  Instead we may
     * go ahead and package them up into bundles and reproject them
     * "quantum" amount at a time to make things run faster.  
     * 
     * This progresses untill the last bundle where we force ourselfs
     * to go 1-by-1.  [Could think of a 'gracefully' algorithm like
     * divide step by 2, or go in sqrt() steps, but we'll leave that
     * for later]
     * 
     */

    for (vv=0;vv<val;vv=vv+step) {
        if (vv+step>val) step=1;
        outpix = sample_distro( buffer, refImage );
        out_data[outpix] += step;            
    } /* end for vv */


}

/* Allocate polygon structs */
Polygons *make_polys( int subpix )
{
    Polygons *polys = (Polygons*)calloc(1,sizeof(Polygons));
    polys->ref = make_polygon(subpix);
    polys->tmp = make_polygon(subpix);
    polys->clip = make_polygon(subpix);
    return(polys);
}


/*
 *  The way resample image works is that it determines the amount of
 * area each input covers in each of the output pixels. So for example
 * if the input image is binned by 4 and the reference/output image is 
 * binned by 1, then 1 pixel in the input will map to 16 pixels in the 
 * output.
 * 
 * The area is then used to randomly sample the **counts** (ie 'val')
 * from the input image, sending some fraction of the counts into
 * each of the output pixels based on the relative area coverage.
 * 
 * This preserves the integer nature of the data (input must be int
 * datatype).  Which means that the output can be used with
 * for wavdetect or any other tool expecting Possion stats.
 * 
 * The assumption is that the counts are evenly distributed within
 * the input pixel.  We don't know any better so this is the best,
 * unbiased assumption we can make.  We could define some kind of
 * weighting from the center kind of thing but that's messy.
 * 
 */
int process_infile( Image *inImage, Image *refImage, WCS_Descriptors *descs, Parameters *pars, Buffer *buffer, 
    Polygons *polys, double *out_data )
{

    /* For all pixels in the input image, we need to map them to the output, if any */
    long mm,nn;
    for (nn=0;nn<inImage->lAxes[1];nn++) { // infile Y axis
      for (mm=0;mm<inImage->lAxes[0];mm++ ) { // infile X axis

        /* get image value; */
        double val;
        val = get_image_value( inImage->data, inImage->dt, mm, nn, inImage->lAxes, inImage->mask );
        if (is_invalid_val(val)) continue;

        /* This creates a polygon around output pixel mm,nn                */        
        /* find the min and max pixels that need to intersect polygon with */
        long xx_min, xx_max, yy_min, yy_max;
        fill_polygon( pars->subpix, mm, nn, polys->ref->contour, descs, pars->ctype );
        find_bounding_box( polys->ref, refImage->lAxes, &xx_min, &xx_max, &yy_min, &yy_max );




        buffer->len = 0;
        long xx, yy;
        for(yy=yy_min;yy<=yy_max;yy++) { // refImage Y axis
          for (xx=xx_min;xx<=xx_max;xx++ ) {  // refImage X axis
            double area;
            /* clip polgon against the pixel and compute area */
            area = get_contour_area( polys, xx, yy );
            if ( 0 == area ) continue;
            add_to_buffer( buffer, xx, yy, area );

          } /* end for xx (x-axis in ref image ) */
        } /* end for yy (y-axis in ref image) */

        if (0 == buffer->len) continue; // no coverage, skip it

        cummulative( buffer );

        distribute_counts( pars, buffer, refImage, val, out_data );

      } /* end for mm (x-axis in input image) */
    } /* end for nn (y-axis in input image */
    
    return(0);

}


Header_Type *merge_headers( Parameters *pars, long num_infiles, Header_Type **hdr )
{

  if ( 1 == num_infiles ) 
    return(hdr[0]);
  if ( 0 == strlen(pars->lookup)) 
    return(hdr[0]);
  if ( 0 == ds_strcmp_cis(pars->lookup, "none")) 
    return(hdr[0]);

  return( mergeHdr( pars->lookup, hdr, num_infiles ));
}


int write_output( Parameters *pars, Image *refImage, double *out_data, long num_infiles, Header_Type **hdr )
{
  /* merge headers */
  Header_Type *outhdr;
  outhdr = merge_headers( pars, num_infiles, hdr );

  /* create output image  */
  dmBlock *outBlock = NULL;
  dmDescriptor *outDesc;
  
  char unit[100];
  memset(unit, 0, 100);
  dmGetUnit( dmImageGetDataDescriptor( refImage->block), unit, 99 );
    
  
  if ( NULL == ( outBlock = dmImageCreate( pars->outfile, dmLONG,refImage->lAxes,2 ))){
      err_msg("ERROR: Cannot create output image '%s'\n", pars->outfile );
      return(-1);
  }
  outDesc = dmImageGetDataDescriptor( outBlock );
  putHdr( outBlock, hdrDM_FILE, outhdr, ALL_STS, "resample_image");
  put_param_hist_info( outBlock, "resample_image", NULL, 0 );
  dmSetUnit( outDesc, unit );
  dmBlockCopyWCS( refImage->block, outBlock);
  dmSetArray_d( outDesc, out_data, refImage->lAxes[0]*refImage->lAxes[1]);
  dmImageClose( outBlock );
  dmImageClose( refImage->block );

  return(0);
}


//---------------------------------------------------------
/* Now onto the main routine */

int resample_img(void)
{

  Parameters *pars;
  Polygons *polys;
  Buffer *buffer;
  WCS_Descriptors *descs;
  
  if ( NULL == ( pars = get_parameters())) {
      return(-1); // error message internal
  }

  /* These polys are used to store regions around each point */
  polys = make_polys( pars->subpix);

  /* These buffer values contain the fraction of areas the input pixel covers
   * in the reference & output image */
  buffer = make_buffer();

  /* Store the coordinate descriptors */
  descs = (WCS_Descriptors*)calloc(1,sizeof( WCS_Descriptors));

  /* Load the ref image info*/
  Image *refImage;
  if ( NULL == (refImage = load_ref_image( pars, descs))) {
      return(-1);
  }

  /* Allocate output array */
  double *out_data;
  out_data = (double*)calloc(refImage->lAxes[0]*refImage->lAxes[1], sizeof(double));

  /* Now let's start on the input stack */
  Stack inStack;
  inStack = stk_build( pars->instack );

  long num_infiles;
  num_infiles = stk_count(inStack);

  Header_Type **hdr;
  hdr = (Header_Type**) calloc( num_infiles, sizeof(Header_Type*));

  /* Begin loop over input files */  
  long ii;
  for (ii=0;ii<num_infiles;ii++) {
    char *infile;                  /* individual image in stack */
    Image *inImage;

    infile = stk_read_num(inStack, ii+1);
    if ( NULL == ( inImage = load_infile_image( pars, infile, hdr+ii, descs ))) {
        return(-1);
    }    

    process_infile( inImage, refImage, descs, pars, buffer, polys, out_data );

    dmImageClose( inImage->block ); // Need to leave open to do coord x-forms
    free_image( inImage );

  }  /* end for ii */


  write_output( pars, refImage, out_data, num_infiles, hdr );

  return(0);


}
