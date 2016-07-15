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

#include "repro_img.h"

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
  short do_norm;
  char which_norm[30];
  char csys[30];
  CoordType ctype;
  long quantum;

} Parameters;


/* Buffer to hold areas */
typedef struct {
  long maxlen;
  long len;
  long *xx;
  long *yy;
  double *area;      
} Buffer;



double get_contour_area( Polygon *clip_poly );
int fill_polygon( long subpix, long xx, long yy, VertexList *refpixlist, WCS_Descriptors *descs, CoordType ctype);
int convert_coords( double *refpix, WCS_Descriptors *descs, double *imgpix, CoordType ctype);
int check_coords( CoordType ctype, WCS_Descriptors descs );
Image *load_image_file( dmBlock *inBlock );
Polygon *make_polygon(int subpix);
int find_bounding_box( Polygon *ref_poly, long *xx_min, long *xx_max, long *yy_min, long*yy_max) ;
Parameters* get_parameters(void);
Image *load_ref_image( Parameters *pars, WCS_Descriptors *descs, dmBlock **refBlock, double **out_data);
Image *load_infile_image(Parameters *pars, char *infile, dmBlock **inBlock, Header_Type **hdr, WCS_Descriptors *descs  );
void free_image( Image *img );
Buffer *make_buffer(void);
void add_to_buffer( Buffer *buffer, long xx, long yy, double area );
void cummulative( Buffer *buffer );
long sample_distro( Buffer *buffer, Image *refImage);
int is_invalid_val( double val );
void distribute_counts( Parameters *pars, Buffer *buffer, Image *refImage, double val, double *out_data);

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
    if ( dmSUCCESS != dmCoordCalc_d( descs->ref.physical.xaxis, refpix, refphys ) ) retval=-1;
    if ( dmSUCCESS != dmCoordCalc_d( descs->ref.world.xaxis, refphys, refwrld )) retval=-1;
    if ( dmSUCCESS != dmCoordInvert_d( descs->image.world.xaxis, refwrld, imgphys )) retval=-1;
    if ( dmSUCCESS != dmCoordInvert_d( descs->image.physical.xaxis, imgphys, imgpix )) retval=-1;
    
    if ( descs->ref.physical.yaxis ) {
      if ( dmSUCCESS != dmCoordCalc_d( descs->ref.physical.yaxis, refpix+1, refphys+1 )) retval=-1;
      if ( dmSUCCESS != dmCoordCalc_d( descs->ref.world.yaxis, refphys+1, refwrld+1 )) retval=-1;
      if ( dmSUCCESS != dmCoordInvert_d(  descs->image.world.yaxis, refwrld+1, imgphys+1 )) retval=-1;
      if ( dmSUCCESS != dmCoordInvert_d( descs->image.physical.yaxis, imgphys+1, imgpix+1 )) retval=-1;
    }
    break;
  }

  case coordLOGICAL: {
    imgpix[0] = refpix[0];
    imgpix[1] = refpix[1];
    break;
  }

  case coordPHYSICAL: {
    if ( dmSUCCESS != dmCoordCalc_d( descs->ref.physical.xaxis, refpix, refphys )) retval=-1;
    if ( dmSUCCESS != dmCoordInvert_d( descs->image.physical.xaxis, refphys, imgpix )) retval=-1;
    if (  descs->ref.physical.yaxis ) {
      if ( dmSUCCESS != dmCoordCalc_d( descs->ref.physical.yaxis, refpix+1, refphys+1 )) retval=-1;
      if ( dmSUCCESS != dmCoordInvert_d( descs->image.physical.yaxis, refphys+1, imgpix+1 )) retval=-1;
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


  /* This is new here ...
   * 
   * In reproject image, we take the output pixel and map it back
   * to the input pixels, accumulating the fraction of the area it covers
   * 
   * In resampling we need to take the input pixel and determine which
   * pixels in the output it will overlap.  I could re-write these routines,
   * but to keep kind of in sync with reproject_image I'll keep as are.
   * 
   * Basically need to swap 'ref' and 'image' in the above convert_coords
   * routine.
   * 
   */
  WCS tmp_wcs;
  tmp_wcs = descs->ref;
  descs->ref = descs->image;
  descs->image = tmp_wcs;

  
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
  
  /* And now restore these back to their original order */  
  tmp_wcs = descs->ref;
  descs->ref = descs->image;
  descs->image = tmp_wcs;

  return(0);
}


/*
  After the polygon has been clipped; we need to compute the area ==
  0.5 * sum ( x_i*y_(i+1) - x_(i+1)*y(i) )

*/
double get_contour_area( Polygon *clip_poly )
{
  double area;
  double larea = 0;
  long zz;

  area = 0;
    
  for (zz=0;zz<clip_poly->contour->num_vertices-1;zz++) {
      larea += (( clip_poly->contour->vertex[zz].x *
                  clip_poly->contour->vertex[zz+1].y ) -
                ( clip_poly->contour->vertex[zz+1].x *
                  clip_poly->contour->vertex[zz].y  ));
  }
    /* last point */
  larea += (( clip_poly->contour->vertex[zz].x *
              clip_poly->contour->vertex[0].y ) -
            ( clip_poly->contour->vertex[0].x *
              clip_poly->contour->vertex[zz].y ) );
    
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


/* find the min/max x,y values from the polygon */
int find_bounding_box( Polygon *ref_poly, long *xx_min, long *xx_max, long *yy_min, long*yy_max) 
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
  clgetstr("method", pars->which_norm, 29);
  clgetstr("coord_sys", pars->csys, 29);
  clgetstr( "lookupTab", pars->lookup, DS_SZ_FNAME);
  pars->clobber = clgetb( "clobber" );
  pars->verbose = clgeti( "verbose" );
  pars->quantum = 4;
  if ( pars->verbose ) {  
    printf("resample_image - parameters\n");
    printf("%15s = %-s\n", "infile", pars->instack );
    printf("%15s = %-s\n", "matchfile", pars->reffile );
    printf("%15s = %-s\n", "outfile", pars->outfile );
    printf("%15s = %d\n", "resolution", pars->subpix );
    printf("%15s = %-s\n", "method", pars->which_norm );
    printf("%15s = %-s\n", "coord_sys", pars->csys );
    printf("%15s = %-s\n", "lookupTab", pars->lookup );
    printf("%15s = %-s\n", "clobber", (pars->clobber ? "yes" : "no") );
    printf("%15s = %d\n", "verbose", pars->verbose );
  }

  if ( ( strlen( pars->reffile) == 0 ) ||
       ( ds_strcmp_cis(pars->reffile, "none" ) == 0 ) ) {

    err_msg("ERROR: Must supply a valid match file\n");
    return(NULL);
   }

  pars->do_norm =  ( strcmp( pars->which_norm, "sum" ) == 0 ) ? 0 : 1;

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
  srand48( time(NULL));


  return(pars);
    
}


/* Load the reference image.  Really only care about the WCS here */
Image *load_ref_image( Parameters *pars, WCS_Descriptors *descs, dmBlock **refBlock, double **out_data)
{
  Image *refImage;

  if ( NULL == ( *refBlock = dmImageOpen( pars->reffile) ) ) {
      err_msg("ERROR: Cannot open image '%s'\n", pars->reffile );
      return(NULL);
  }
  
  if ( NULL == ( refImage = load_image_file( *refBlock ))) {
      err_msg("ERROR: Cannot open image '%s'\n", pars->reffile );
      return(NULL);      
  }

  descs->ref.physical.xaxis = refImage->xdesc;
  descs->ref.physical.yaxis = refImage->ydesc;
  descs->ref.world.xaxis = dmDescriptorGetCoord( refImage->xdesc );
  descs->ref.world.yaxis = dmDescriptorGetCoord( refImage->ydesc );

  /* Alloc outptu data array */
  *out_data = ( double*)calloc( refImage->lAxes[0]*refImage->lAxes[1], sizeof(double));

  return(refImage);
}


/* Load the infile image. */
Image *load_infile_image(Parameters *pars, char *infile, dmBlock **inBlock, Header_Type **hdr, WCS_Descriptors *descs  )
{
    Image *inImage;

    if ( pars->verbose > 1 ) {
      printf("\nProcessing input file '%s'\n", infile );
    }
     
    /* Now load the image */
    if ( NULL == ( *inBlock = dmImageOpen( infile) ) ) {
      err_msg("ERROR: Cannot open image '%s'\n", infile );
      return(NULL);
    }

    if ( NULL == ( inImage = load_image_file( *inBlock ))) {
      err_msg("ERROR: Cannot open image '%s'\n", infile );
      return(NULL);        
    }

    /* The infile must be integer ! */
    if (( dmFLOAT == inImage->dt ) ||
        ( dmDOUBLE == inImage->dt ) ) {
        err_msg("ERROR: The input image must be an integer datatype.");
        return(NULL);
    }
    
    
    *hdr = getHdr( *inBlock, hdrDM_FILE );

    descs->image.physical.xaxis = inImage->xdesc;
    descs->image.physical.yaxis = inImage->ydesc;
    descs->image.world.xaxis = dmDescriptorGetCoord( inImage->xdesc );
    descs->image.world.yaxis = dmDescriptorGetCoord( inImage->ydesc );
    
    if ( 0 != check_coords( pars->ctype, *descs ) ) {
        return(NULL);
    }

    return(inImage);
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
 * The buffer stuff then is keeping track of which pixels 
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


int is_invalid_val( double val )
{
    static int show_less_than_zero_message = 1;

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


//---------------------------------------------------------
/* Now onto the main routine */

int resample_img(void)
{

  Parameters *pars;
  
  if ( NULL == ( pars = get_parameters())) {
      return(-1); // error message internal
  }

  /* These polys are used to store regions around each point */
  Polygon *ref_poly, *clip_poly, *tmp_clip_poly;    
  ref_poly = make_polygon(pars->subpix);
  clip_poly = make_polygon(pars->subpix);
  tmp_clip_poly = make_polygon(pars->subpix);

  /* These buffer values contain the fraction of areas the input pixel covers
   * in the reference & output image */
  Buffer *buffer;
  buffer = make_buffer();

  Image *refImage;
  WCS_Descriptors descs;
  double *out_data;
  dmBlock *refBlock;  
  if ( NULL == (refImage = load_ref_image( pars, &descs, &refBlock, &out_data))) {
      return(-1);
  }

  /* Now let's start on the input stack */
  Stack inStack;
  long num_infiles;
  Header_Type **hdr;

  inStack = stk_build( pars->instack );
  num_infiles = stk_count(inStack);
  hdr = (Header_Type**) calloc( num_infiles, sizeof(Header_Type*));
  stk_rewind(inStack);
  
  char *infile;                  /* individual image in stack */
  while ( NULL != (infile = stk_read_next(inStack))) {

    Image *inImage;
    dmBlock *inBlock;    
    num_infiles--;
    if ( NULL == ( inImage = load_infile_image( pars, infile, &inBlock, hdr+num_infiles, &descs ))) {
        return(-1);
    }    

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
        fill_polygon( pars->subpix, mm, nn, ref_poly->contour, &descs, pars->ctype );
        find_bounding_box( ref_poly, &xx_min, &xx_max, &yy_min, &yy_max );

        buffer->len = 0;
        long xx, yy;
        for(yy=yy_min;yy<=yy_max;yy++) { // refImage Y axis
          for (xx=xx_min;xx<=xx_max;xx++ ) {  // refImage X axis
            double area;

            /* clip polgon against the pixel and compute area */
            super_poly_clip( ref_poly, tmp_clip_poly, clip_poly, xx, yy );
            area = get_contour_area( clip_poly );
            if ( 0 == area ) continue;
            add_to_buffer( buffer, xx, yy, area );

          } /* end for xx (x-axis in ref image ) */
        } /* end for yy (y-axis in ref image) */

        if (0 == buffer->len) continue; // no coverage, skip it

        cummulative( buffer );

        distribute_counts( pars, buffer, refImage, val, out_data );

      } /* end for mm (x-axis in input image) */
    } /* end for nn (y-axis in input image */
    

    if ( pars->verbose > 2 )  printf("Processing %g%% complete        \n", 100.0 );
    free_image( inImage );
    dmImageClose( inBlock ); // Need to leave open to do coord x-forms

  }  /* end while infile loop over stack */


  /* merge headers */
  num_infiles = stk_count( inStack );
  Header_Type *outhdr;
  if (( (strlen( pars->lookup ) == 0 ) ||
        (ds_strcmp_cis(pars->lookup, "none") ==0 )) ||
      ( num_infiles == 1 ) ) {
    outhdr = hdr[0];
  } else {
    outhdr = mergeHdr( pars->lookup, hdr, num_infiles );
  }


  /* create output image  */
  dmBlock *outBlock = NULL;
  dmDescriptor *outDesc;
  
  if ( NULL == ( outBlock = dmImageCreate( pars->outfile, dmLONG,refImage->lAxes,2 ))){
      err_msg("ERROR: Cannot create output image '%s'\n", pars->outfile );
      return(-1);
  }
  outDesc = dmImageGetDataDescriptor( outBlock );
  putHdr( outBlock, hdrDM_FILE, outhdr, ALL_STS, "resample_image");
  put_param_hist_info( outBlock, "resample_image", NULL, 0 );
  dmBlockCopyWCS( refBlock, outBlock);
  dmSetArray_d( outDesc, out_data, refImage->lAxes[0]*refImage->lAxes[1]);
  dmImageClose(outBlock );
  dmImageClose( refBlock );

  return(0);


}
