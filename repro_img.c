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
#include <cxcregion.h>


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
typedef struct {
 WCS ref; 
 WCS image;
 } WCS_Descriptors;



/* Hold info for an input image */
typedef struct {
  void *data;        // pixel values
  dmDataType dt;     // pixel datatype
  long *lAxes;       // axis lenghts
  short *mask;        // mask of valid pixels
  dmDescriptor *xdesc;  // X (or sky) coordinate descriptor
  dmDescriptor *ydesc;  // Y coordinate descriptor
} Image;






/*
  All these get_ routines are pilfered from other (too be released) DM tools 
  that do some basic 2D image I/O.  Allows an arbitrary image data-type to be
  read in and a mask created to indicate which pixels are good & bad.
*/
int find_chip_corners( long *min_ref_x, long *min_ref_y,
                       long *max_ref_x, long *max_ref_y,
                       long *reflen, long *inlen,
                       WCS_Descriptors *desc );
double get_contour_area( Polygon *clip_poly );

int fill_polygon( long subpix,
                  long xx, long yy,
                  VertexList *refpixlist,
                  WCS_Descriptors *descs,
                  CoordType ctype);
                  
int convert_coords( double *refpix,
                    WCS_Descriptors *descs,
                    double *imgpix,
                    CoordType ctype);

int check_coords( CoordType ctype, WCS_Descriptors descs );
                       
Image *load_image_file( dmBlock *inBlock );
int make_polygon(int subpix, Polygon *poly);
int find_bounding_box( Polygon ref_poly, long *xx_min, long *xx_max, long *yy_min, long*yy_max) ;
void super_poly_clip( Polygon *ref_poly, Polygon *tmp_clip_poly, Polygon *clip_poly, long mm, long nn ) ;


// --------------------------



/* Load image file.  This uses the routines defined in dmimgio.h */
#include <dmimgio.h>
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

       
  case coordWORLD:
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

  case coordLOGICAL:
    imgpix[0] = refpix[0];
    imgpix[1] = refpix[1];
    break;
    
  case coordPHYSICAL:
    if ( dmSUCCESS != dmCoordCalc_d( descs->ref.physical.xaxis, refpix, refphys )) retval=-1;
    if ( dmSUCCESS != dmCoordInvert_d( descs->image.physical.xaxis, refphys, imgpix )) retval=-1;
    if (  descs->ref.physical.yaxis ) {
      if ( dmSUCCESS != dmCoordCalc_d( descs->ref.physical.yaxis, refpix+1, refphys+1 )) retval=-1;
      if ( dmSUCCESS != dmCoordInvert_d( descs->image.physical.yaxis, refphys+1, imgpix+1 )) retval=-1;
    }

    break;
  } /* end switch */


  return(retval);
}
                    

/*
   Okay, now we are down to the real reproject_image* code.

*/


/*

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
    int stat;
    stat = convert_coords( refpix, descs, imgpix,ctype);
    if ( 0 != stat ) {
    }
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
    int stat;
    stat = convert_coords( refpix, descs, imgpix,ctype);
    if ( 0 != stat ) {
    }
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
    int stat;
    stat = convert_coords( refpix, descs, imgpix,ctype);
    if ( 0 != stat ) {
    }
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
    int stat;
    stat = convert_coords( refpix, descs, imgpix,ctype);
    if ( 0 != stat ) {
    }
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

double get_contour_area( Polygon *clip_poly )
{
  double area;
  long kk=0;

  area = 0;
    double larea = 0;
    long zz;
    
    for (zz=0;zz<clip_poly->contour[kk].num_vertices-1;zz++) {
      larea += (( clip_poly->contour[kk].vertex[zz].x *
                  clip_poly->contour[kk].vertex[zz+1].y ) -
                ( clip_poly->contour[kk].vertex[zz+1].x *
                  clip_poly->contour[kk].vertex[zz].y  ));
    }
    /* last point */
    larea += (( clip_poly->contour[kk].vertex[zz].x *
                clip_poly->contour[kk].vertex[0].y ) -
              ( clip_poly->contour[kk].vertex[0].x *
                clip_poly->contour[kk].vertex[zz].y ) );
    
    larea /= 2.0;
    
      area += fabs(larea);

  return(area);
}




/*
  Okay has nothing to do wtih 'chip''s 

  This is a routine that does the inverse of the bulk of the rest of the
  processing.  It maps the input image corners to the output image.

  This gives use a potentially much smaller area to iterate over when working in 
  the other direction (from output pixels to input pixels).


*/

int find_chip_corners( long *min_ref_x, long *min_ref_y,
                       long *max_ref_x, long *max_ref_y,
                       long *reflen, long *inlen,
                       WCS_Descriptors *desc )
{

  double in_log_loc[2];
  double ref_log_loc[2];
  long vals_x[2], vals_y[2];


    /* In this case we are going from the input image to the reference image instead
       of how the other code goes from the ref. to the input.  So we just
       swap the pointers.  Have to be sure to return to other order before
       we get back to the main routine */

    WCS tmp_wcs;
    tmp_wcs = desc->ref;
    desc->ref = desc->image;
    desc->image = tmp_wcs;

    /* Why + and - 2?  1 to get from image to c-array index'es and 1 to get
       to the end of the pixel.  Could probably use -1 and +2.
       Ranges get clipped at the end so extra padding doesn't
       hurt
    */    

    long xlen, ylen, ilen, jlen;
    xlen = reflen[0];
    ylen = reflen[1];
    ilen = inlen[0];
    jlen = inlen[1];


    vals_x[0] = -2; 
    vals_y[0] = -2;
    vals_x[1] = ilen+2;
    vals_y[1] = jlen+2;
    
    *min_ref_x = xlen+1;
    *min_ref_y = ylen+1;
    *max_ref_x = 0;
    *max_ref_y = 0;

    long ii,jj;
    for (ii=0;ii<2;ii++) {
      for (jj=0;jj<2;jj++ ) {

        in_log_loc[0] = vals_x[ii];
        in_log_loc[1] = vals_y[jj];
        
        if ( 0 != convert_coords( in_log_loc, desc, ref_log_loc, coordWORLD ) ) {
          tmp_wcs = desc->ref;
          desc->ref = desc->image;
          desc->image = tmp_wcs;
          *min_ref_x = 0;
          *min_ref_y = 0;
          *max_ref_x = xlen;
          *max_ref_y = ylen;
          return(0);
        }

        *min_ref_x = MIN( *min_ref_x, (ref_log_loc[0]-1));
        *min_ref_y = MIN( *min_ref_y, (ref_log_loc[1]-1));
        *max_ref_x = MAX( *max_ref_x, (ref_log_loc[0]-1));
        *max_ref_y = MAX( *max_ref_y, (ref_log_loc[1]-1));
        
      } /* end for jj */
    } /* end for ii */
    
    if ( *min_ref_x < 0 ) *min_ref_x = 0;
    if ( *min_ref_y < 0 ) *min_ref_y = 0;
    if ( *max_ref_x > xlen ) *max_ref_x = xlen;
    if ( *max_ref_y > ylen ) *max_ref_y = ylen;
    
    if ( *max_ref_x < 0 ) *max_ref_x = 0; /* Well, not on image! */
    if ( *max_ref_y < 0 ) *max_ref_y = 0; /* Well, not on image! */
    if ( *min_ref_x > xlen ) *min_ref_x = xlen;
    if ( *min_ref_y > ylen ) *min_ref_y = ylen;
  
    tmp_wcs = desc->ref;
    desc->ref = desc->image;
    desc->image = tmp_wcs;

    return(0);

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


int make_polygon(int subpix, Polygon *poly)
{
    VertexList *refpixlist = (VertexList*)calloc(1,sizeof(VertexList));   
    Polygon *ref_poly = (Polygon*)calloc(1,sizeof(Polygon));
        

    refpixlist->vertex    = (Vertex*)calloc(4*subpix*4, sizeof(Vertex)); // why 4*4, 4 sides
    refpixlist->num_vertices = 4*subpix;    
    ref_poly->contour      = refpixlist;
    *poly = *ref_poly;

    return(0);
}



int find_bounding_box( Polygon ref_poly, long *xx_min, long *xx_max, long *yy_min, long*yy_max) 
{
    long ii;
 
    *xx_min = ref_poly.contour[0].vertex[0].x-0.5;
    *yy_min = ref_poly.contour[0].vertex[0].y-0.5;
    *xx_max = ref_poly.contour[0].vertex[0].x+0.5;
    *yy_max = ref_poly.contour[0].vertex[0].y+0.5;
    for(ii=ref_poly.contour[0].num_vertices;--ii;) {
      *xx_min = MIN( *xx_min, (ref_poly.contour[0].vertex[ii].x-0.5));
      *yy_min = MIN( *yy_min, (ref_poly.contour[0].vertex[ii].y-0.5));
      *xx_max = MAX( *xx_max, (ref_poly.contour[0].vertex[ii].x+0.5));
      *yy_max = MAX( *yy_max, (ref_poly.contour[0].vertex[ii].y+0.5));
    }

    return(0);
}


void super_poly_clip( Polygon *ref_poly, Polygon *tmp_clip_poly, Polygon *clip_poly, long mm, long nn ) 
{
    /* intersect polygons ; once for each side of the rectangle*/
    poly_clip( ref_poly,      mm, nn, tmp_clip_poly ,LEFT );
    poly_clip( tmp_clip_poly, mm, nn, clip_poly     ,RIGHT );
    poly_clip( clip_poly,     mm, nn, tmp_clip_poly ,BOTTOM );
    poly_clip( tmp_clip_poly, mm, nn, clip_poly     ,TOP );
}







/* Now onto the main routine */

int resample_img(void)
{

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
  CoordType ctype = coordWORLD ;


  /* Get the parameters */
  clgetstr( "infile", instack, DS_SZ_FNAME );
  clgetstr( "matchfile", reffile, DS_SZ_FNAME );
  clgetstr( "outfile", outfile, DS_SZ_FNAME );
  subpix = clgeti("resolution");
  clgetstr("method", which_norm, 29);
  clgetstr("coord_sys", csys, 29);
  clgetstr( "lookupTab", lookup, DS_SZ_FNAME);
  clobber = clgetb( "clobber" );
  verbose = clgeti( "verbose" );
  
  if ( verbose ) {  
    printf("resample_image - parameters\n");
    printf("%15s = %-s\n", "infile", instack );
    printf("%15s = %-s\n", "matchfile", reffile );
    printf("%15s = %-s\n", "outfile", outfile );
    printf("%15s = %d\n", "resolution", subpix );
    printf("%15s = %-s\n", "method", which_norm );
    printf("%15s = %-s\n", "coord_sys", csys );
    printf("%15s = %-s\n", "lookupTab", lookup );
    printf("%15s = %-s\n", "clobber", (clobber ? "yes" : "no") );
    printf("%15s = %d\n", "verbose", verbose );
  }

  if ( ( strlen( reffile) == 0 ) ||
       ( ds_strcmp_cis(reffile, "none" ) == 0 ) ) {

    err_msg("ERROR: Must supply a valid match file\n");
    return(-1);
   }

  do_norm =  ( strcmp( which_norm, "sum" ) == 0 ) ? 0 : 1;

  switch ( csys[0] ) {
    case 'l': ctype = coordLOGICAL; break;
    case 'p': ctype = coordPHYSICAL; break;
    case 'w': ctype = coordWORLD; break;
    default:
      err_msg("ERROR: Unknow coordinate type '%s'\n", csys );
    return(-1);
  }

  /* ----------------------------- */


  if ( ds_clobber( outfile, clobber, NULL ) != 0 ) {
    return(-1);
  }

    // ------------------------




  dmBlock *refBlock;
  Image *refImage;
  if ( NULL == ( refBlock = dmImageOpen( reffile) ) ) {
      err_msg("ERROR: Cannot open image '%s'\n", reffile );
      return(-1);
  }
  
  if ( NULL == ( refImage = load_image_file( refBlock ))) {
      err_msg("ERROR: Cannot open image '%s'\n", reffile );
      return(-1);      
  }


  WCS_Descriptors descs;
  descs.ref.physical.xaxis = refImage->xdesc;
  descs.ref.physical.yaxis = refImage->ydesc;
  descs.ref.world.xaxis = dmDescriptorGetCoord( refImage->xdesc );
  descs.ref.world.yaxis = dmDescriptorGetCoord( refImage->ydesc );


  /* Alloc outptu data array */
  double *out_data;
  out_data = ( double*)calloc( refImage->lAxes[0]*refImage->lAxes[1], sizeof(double));



  /* Now let's start on the input stack */
  Stack inStack;
  long num_infiles;
  Header_Type **hdr;

  inStack = stk_build( instack );
  if ( ( NULL == inStack ) || ( stk_count(inStack)==0 ) ||
       (( stk_count(inStack)==1 ) && ( strlen(stk_read_num(inStack,1))==0))) {
    err_msg("ERROR: problems opening stack '%s'\n", instack );
    return(-3);
  }    
  num_infiles = stk_count(inStack);
  hdr = (Header_Type**) calloc( num_infiles, sizeof(Header_Type*));
  stk_rewind(inStack);


    
  char *infile;                  /* individual image in stack */
  while ( NULL != (infile = stk_read_next(inStack))) {

    if ( verbose > 1 ) {
      printf("\nProcessing input file '%s'\n", infile );
    }
    
    
    /* Now load the image */
    dmBlock *inBlock;
    Image *inImage;
    if ( NULL == ( inBlock = dmImageOpen( infile) ) ) {
      err_msg("ERROR: Cannot open image '%s'\n", infile );
      return(-1);
    }

    if ( NULL == ( inImage = load_image_file( inBlock ))) {
      err_msg("ERROR: Cannot open image '%s'\n", infile );
      return(-1);        
    }
    
    hdr[--num_infiles] = getHdr( inBlock, hdrDM_FILE );

    descs.image.physical.xaxis = inImage->xdesc;
    descs.image.physical.yaxis = inImage->ydesc;
    descs.image.world.xaxis = dmDescriptorGetCoord( inImage->xdesc );
    descs.image.world.yaxis = dmDescriptorGetCoord( inImage->ydesc );
    
    if ( 0 != check_coords( ctype, descs ) ) {
        return(-1);
    }

    
    Polygon ref_poly, clip_poly, tmp_clip_poly;    
    make_polygon(subpix, &ref_poly);
    make_polygon(subpix, &clip_poly);
    make_polygon(subpix, &tmp_clip_poly);

    /* Okay, we'll look to the input image and map the corners
       to the image in the output.  This way we only have to process those
       pixels in the output which could conver the input.  Should
       save time especially when input image is small piece of the
       output image ... eg mosaics */
    long min_ref_x,min_ref_y,max_ref_x,max_ref_y;
    find_chip_corners( &min_ref_x, &min_ref_y, &max_ref_x, &max_ref_y, 
                 refImage->lAxes, inImage->lAxes, &descs );

    /* Begin loop through data */

    long xx, yy;
    for(yy=min_ref_y;yy<max_ref_y;yy++) {
      long ypix_off;
      ypix_off = yy*refImage->lAxes[0];

      if ( verbose > 2 ) {
        double perc = 100.0 * yy / max_ref_y;
        printf("Processing %g%% complete        \r", perc );
      }

      for (xx=min_ref_x;xx<max_ref_x;xx++ ) {
        
        /* This creates a polygon around output pixel xx,yy */
        fill_polygon( subpix, xx, yy,ref_poly.contour, &descs, ctype );
        
        /* find the min and max pixels that need to intersect polygon with */
        long xx_min, xx_max, yy_min, yy_max;
        find_bounding_box( ref_poly, &xx_min, &xx_max, &yy_min, &yy_max );
        

        double weight;  
        double sum;
        weight  = 0;
        sum = 0;

        /* for each pixel that could intersect poly check it and get area */
        long mm,nn;
        for (nn=yy_min;nn<=yy_max;nn++) {
          for (mm=xx_min;mm<=xx_max;mm++ ) {
                    
            double area;
            double val;

            super_poly_clip( &ref_poly, &tmp_clip_poly, &clip_poly, mm, nn );

            /* area of polygon */
            area = get_contour_area( &clip_poly );

            /* get image value; weight output pixel by it */
            val = get_image_value( inImage->data, inImage->dt, mm, nn, inImage->lAxes, inImage->mask );
            
            if ( !ds_dNAN(val) )
              sum += ( val * area );
            
            weight += area;
            
          } /* end for mm */
        } /* end for nn */
        
        
        /* Can either sum up values or compute the average value; this is done here */
        long outpix;
        outpix = ypix_off + xx;
        if ( weight > 0 ) {
          if ( do_norm )
            out_data[outpix] += (sum/weight);
          else
            out_data[outpix] += sum;
        }

        
      } /* end for xx */
    } /* end for yy */

    if ( verbose > 2 )  printf("Processing %g%% complete        \n", 100.0 );


    // TODO free(data);
    // TODO free(lAxes);
    // TODO if ( mask ) free(mask);
    dmImageClose( inBlock );

  }  /* end while infile loop over stack */


  /* merge headers */
  num_infiles = stk_count( inStack );
  Header_Type *outhdr;
  if (( (strlen( lookup ) == 0 ) ||
        (ds_strcmp_cis(lookup, "none") ==0 )) ||
      ( num_infiles == 1 ) ) {
    outhdr = hdr[0];
  } else {
    outhdr = mergeHdr( lookup, hdr, num_infiles );
  }



  /* create output image  */
  dmBlock *outBlock = NULL;
  dmDescriptor *outDesc;
  
  if ( NULL == ( outBlock = dmImageCreate( outfile, dmDOUBLE,refImage->lAxes,2 ))){
      err_msg("ERROR: Cannot create output image '%s'\n", outfile );
      return(-1);
  }
  outDesc = dmImageGetDataDescriptor( outBlock );
  putHdr( outBlock, hdrDM_FILE, outhdr, ALL_STS, "resample_image");
  put_param_hist_info( outBlock, "resample_image", NULL, 0 );
  dmBlockCopyWCS( refBlock, outBlock);
  dmSetArray_d( outDesc, out_data, refImage->lAxes[0]*refImage->lAxes[1]);
  dmImageClose(outBlock );
  if ( outBlock != refBlock ) dmImageClose( refBlock );

  return(0);


}
