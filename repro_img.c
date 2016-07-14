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


/*
  All these get_ routines are pilfered from other (too be released) DM tools 
  that do some basic 2D image I/O.  Allows an arbitrary image data-type to be
  read in and a mask created to indicate which pixels are good & bad.
*/
int find_chip_corners( long *min_ref_x, long *min_ref_y,
                       long *max_ref_x, long *max_ref_y,
                       long xlen, long ylen,
                       long ilen, long jlen,
                       WCS_Descriptors *desc );
double get_contour_area( Polygon *clip_poly );

int fill_polygon( long subpix,
                  double delta,
                  long xx, long yy,
                  VertexList *refpixlist,
                  WCS_Descriptors *descs,
                  CoordType ctype);
                  
int convert_coords( double *refpix,
                    WCS_Descriptors *descs,
                    double *imgpix,
                    CoordType ctype);

                       

// --------------------------

short *get_image_mask( dmBlock *inBlock, void *data, dmDataType dt, 
                       long *lAxes, regRegion *dss, long null, short has_null, 
                       dmDescriptor *xAxis, dmDescriptor *yAxis );


double get_image_value( void *data, dmDataType dt, 
                        long xx, long yy, long *lAxes, 
                        short *mask );

dmDataType get_image_data( dmBlock *inBlock, void **data, long **lAxes,
                           regRegion **dss, long *nullval, short *nullset );


short  get_image_wcs( dmBlock *imgBlock, dmDescriptor **xAxis, 
                      dmDescriptor **yAxis );

/*
  This routine is used to look at the pixel values to see if they are
  an integer NULL value, a floating point NaN, or if they happen to
  fall outside the data sub space.  If any of these are true
  the output mask will be 0 at that pixel location; otherwise it will be 1.

*/
short *get_image_mask( dmBlock *inBlock, void *data, dmDataType dt, 
                       long *lAxes, regRegion *dss, long null, short has_null, 
                       dmDescriptor *xAxis, dmDescriptor *yAxis )
{
  long npix = lAxes[0] * lAxes[1];
  short *mask;
  long xx, yy;
  mask = (short*)calloc( npix, sizeof(short));

  
  for ( xx=lAxes[0]; xx--; ) {
    for ( yy=lAxes[1]; yy--; ) {
      double dat;
      long idx;
      idx = xx + ( yy * lAxes[0] );
      
      dat = get_image_value( data, dt, xx, yy, lAxes, NULL );
      
      /* Now ... if it is an integer data type, it could possibly have a
         null value. Check for that */
      if ( ( has_null && ( dat == null ) ) ||
           ds_dNAN( dat ) ) {
        continue;
      }
            
      /* If the image has a data sub space (aka a region filter applied)
         then need to convert coords to physical and check */
      if ( dss && xAxis ) {
        double pos[2];
        double loc[2];
        pos[0]=xx+1;
        pos[1]=yy+1;
        
        if (yAxis) {  /* If no y axis, then xAxis has 2 components */
          dmCoordCalc_d( xAxis, pos, loc );
          dmCoordCalc_d( yAxis, pos+1, loc+1 );
        } else {
          dmCoordCalc_d( xAxis, pos, loc );
        }
        if ( !regInsideRegion( dss, loc[0], loc[1] ) )
          continue;
      }
      
      mask[idx] = 1;
    }
  }

  return(mask );
}



/*
  This routine is used to lookup the pixel values and check if the
  pixel is in the mask or not.

  If the pixel is off the image or outside DSS then it returns 0.

*/

double get_image_value( void *data, dmDataType dt, 
                        long xx, long yy, long *lAxes, 
                        short *mask )
{

  long npix = xx + (yy * lAxes[0] );
  double retval;

  /* Okay, first get all the data from the different data types.  
     Cast everything to doubles */


  if (( xx < 0 ) || ( xx >= lAxes[0] ) ||
      ( yy < 0 ) || ( yy >= lAxes[1] )) {
    return(0);
  }



  switch ( dt ) {
    
  case dmBYTE: {
    unsigned char *img = (unsigned char*)data;
    retval = img[npix];
    break;
  }
    
  case dmSHORT: {
    short *img = (short*)data;
    retval = img[npix];
    break;
  }
    
  case dmUSHORT: {
    unsigned short *img = (unsigned short*)data;
    retval = img[npix];
    break;
  }
    
  case dmLONG: {
    long *img = (long*)data;
    retval = img[npix];
    break;
  }
    
  case dmULONG: {
    unsigned long *img = (unsigned long*)data;
    retval = img[npix];
    break;
  }
    
  case dmFLOAT: {
    float *img = (float*)data;
    retval = img[npix];
    break;
  }
  case dmDOUBLE: {
    double *img = (double*)data;
    retval = img[npix];
    break;
  }
  default:
    ds_MAKE_DNAN( retval );

  }


  if ( mask ) {
    if ( !mask[npix] ) {
      ds_MAKE_DNAN( retval );
    }
  }


  return(retval);

}


/* Load the data into memory,  check for DSS, null values */
dmDataType get_image_data( dmBlock *inBlock, void **data, long **lAxes,
                           regRegion **dss, long *nullval, short *nullset )
{

  dmDescriptor *imgDesc;
  dmDataType dt;
  dmDescriptor *grp;
  dmDescriptor *imgdss;

  long naxes;
  long npix;
  char ems[1000];

  *nullval = INDEFL;
  *dss = NULL;
  *nullset = 0;
  
  imgDesc = dmImageGetDataDescriptor( inBlock );

  /* Sanity check, only 2D images */
  naxes = dmGetArrayDimensions( imgDesc, lAxes );
  if ( naxes != 2 ) {
    return( dmUNKNOWNTYPE );
  }
  npix = (*lAxes)[0] * (*lAxes)[1];
  dt = dmGetDataType( imgDesc );


  /* Okay, first lets get the image descriptor */
  grp = dmArrayGetAxisGroup( imgDesc, 1 );
  dmGetName( grp, ems, 1000);
  imgdss = dmSubspaceColOpen( inBlock, ems );
  if ( imgdss )
    *dss = dmSubspaceColGetRegion( imgdss);
  
  
  switch ( dt ) 
    {
    case dmBYTE:
      *data = ( void *)calloc( npix, sizeof(char ));
      dmGetArray_ub( imgDesc, (unsigned char*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmSHORT:
      *data = ( void *)calloc( npix, sizeof(short ));
      dmGetArray_s( imgDesc, (short*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmUSHORT:
      *data = ( void *)calloc( npix, sizeof(short ));
      dmGetArray_us( imgDesc, (unsigned short*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmLONG:
      *data = ( void *)calloc( npix, sizeof(long ));
      dmGetArray_l( imgDesc, (long*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmULONG:
      *data = ( void *)calloc( npix, sizeof(long ));
      dmGetArray_ul( imgDesc, (unsigned long*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmFLOAT:
      *data = ( void *)calloc( npix, sizeof(float ));
      dmGetArray_f( imgDesc, (float*) *data, npix );
      *nullset = 0;
      break;
      
    case dmDOUBLE:
      *data = ( void *)calloc( npix, sizeof(double ));
      dmGetArray_d( imgDesc, (double*) *data, npix );
      *nullset = 0;
      break;
      
    default:
      return( dmUNKNOWNTYPE );
    }

  return(dt);

}




/* Get the WCS descriptor */
short  get_image_wcs( dmBlock *imgBlock, dmDescriptor **xAxis, 
                      dmDescriptor **yAxis )
{
  

  dmDescriptor *imgData;
  long n_axis_groups;

  imgData = dmImageGetDataDescriptor( imgBlock );
  n_axis_groups = dmArrayGetNoAxisGroups( imgData );
  

  /* This is the usual trick ... can have 1 axis group w/ 
     dimensionality 2 (eg a vector column) or can have
     2 axis groups w/ dimensionaity 1 (eg 2 disjoint columns)*/
  if ( n_axis_groups == 1 ) {
    dmDescriptor *pos = dmArrayGetAxisGroup( imgData, 1 );
    dmDescriptor *xcol;
    long n_components;
    
    n_components = dmGetElementDim( pos );
    if ( n_components != 2 ) {
      err_msg("ERROR: could not find 2D image\n");
      return(-1);
    }
    
    xcol = dmGetCpt( pos, 1 );
    
    *xAxis = pos;
    *yAxis = NULL;
    
  } else if ( n_axis_groups == 2 ) {
    dmDescriptor *xcol;
    dmDescriptor *ycol;
  
    xcol = dmArrayGetAxisGroup( imgData, 1 );
    ycol = dmArrayGetAxisGroup( imgData, 2 );

    *xAxis = xcol;
    *yAxis = ycol;
    
  } else {
    err_msg("Invalid number of axis groups\n");
    *xAxis = NULL;
    *yAxis = NULL;
    return(-1);
  }

  return(0);

}





/* ----------------------------------------------------------- */

/*
   Okay, now we are down to the real reproject_image* code.

*/






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

  This routine is used to make a polygon around the square output pixel.  

  The points of this  polygon will then get mapped to the input image
  and the output pixel value will be constructed from those overlapping
  input pixels weighted by the area of the polygon that covers them.
*/

int fill_polygon( long subpix,
                  double delta,
                  long xx, long yy,
                  VertexList *refpixlist,
                  WCS_Descriptors *descs,
                  CoordType ctype)
{

  long refidx;
  double refpix[2];
  double imgpix[2];
  long ii;
  
  
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
                       long xlen, long ylen,
                       long ilen, long jlen,
                       WCS_Descriptors *desc )
{

  double in_log_loc[2];
  double ref_log_loc[2];

    long vals_x[2], vals_y[2];
    long ii,jj;


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
    vals_x[0] = -2; 
    vals_y[0] = -2;
    vals_x[1] = ilen+2;
    vals_y[1] = jlen+2;
    
    *min_ref_x = xlen+1;
    *min_ref_y = ylen+1;
    *max_ref_x = 0;
    *max_ref_y = 0;

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




  void *data;
  long *lAxes;
  regRegion *dss;
  long null;
  short has_null;
  short *mask = NULL;

  dmDataType dt;
  dmBlock *inBlock;
  dmDescriptor *xdesc, *ydesc;
  dmDescriptor *xdesc_w, *ydesc_w;


  
  double *out_data;
  short do_norm;
  char which_norm[30];
  char csys[30];

  CoordType ctype = coordWORLD ;

  WCS_Descriptors descs;



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
  dmDescriptor *refxdesc, *refydesc;
  dmDescriptor *refxdesc_w, *refydesc_w;
  long *refAxes;

  if ( NULL == ( refBlock = dmImageOpen( reffile) ) ) {
      err_msg("ERROR: Cannot open image '%s'\n", reffile );
      return(-1);
  }

  if ( dmUNKNOWNTYPE == ( dt = get_image_data( refBlock, &data,  &refAxes, 
                                                 &dss, &null, &has_null ) ) ) {
      err_msg("ERROR: Cannot get image data or unknown image data-type for "
              "file '%s'\n", reffile);
      return(-1);
  }

  if ( 0 != get_image_wcs( refBlock, &refxdesc, &refydesc ) ) {
      err_msg("ERROR: Cannot load WCS for file '%s'\n", reffile );
      return(-1);
  }

  refxdesc_w = dmDescriptorGetCoord( refxdesc );
  refydesc_w = dmDescriptorGetCoord( refydesc );
  descs.ref.physical.xaxis = refxdesc;
  descs.ref.physical.yaxis = refydesc;
  descs.ref.world.xaxis = refxdesc_w;
  descs.ref.world.yaxis = refydesc_w;

  free(data);  /* We don't need the data for the ref image, just the WCS */

  /* Alloc outptu data array */
  out_data = ( double*)calloc( refAxes[1]*refAxes[0], sizeof(double));



  /* Now let's start on the input stack */

  Stack inStack;

  inStack = stk_build( instack );
  if ( ( NULL == inStack ) ||
       ( stk_count(inStack)==0 ) ||
       (( stk_count(inStack)==1 ) && ( strlen(stk_read_num(inStack,1))==0))) {
    err_msg("ERROR: problems opening stack '%s'\n", instack );
    return(-3);
  }



  /* Okay, start to process the data */
    
  long num_infiles;
  Header_Type **hdr;
  num_infiles = stk_count(inStack);
  hdr = (Header_Type**) calloc( num_infiles, sizeof(Header_Type*));
  
  
  char *infile;                  /* individual image in stack */

  stk_rewind(inStack);
  while ( NULL != (infile = stk_read_next(inStack))) {

    if ( verbose > 1 ) {
      printf("\nProcessing input file '%s'\n", infile );
    }
    
    
    /* Now load the image */
    if ( NULL == ( inBlock = dmImageOpen( infile) ) ) {
      err_msg("ERROR: Cannot open image '%s'\n", infile );
      return(-1);
    }
    
    if ( dmUNKNOWNTYPE == ( dt = get_image_data( inBlock, &data,  &lAxes, 
                                                 &dss, &null, &has_null ) ) ) {
      err_msg("ERROR: Cannot get image data or unknown image data-type for "
              "file '%s'\n", infile);
      return(-1);
    }
    
    if ( 0 != get_image_wcs( inBlock, &xdesc, &ydesc ) ) {
      err_msg("ERROR: Cannot load WCS for file '%s'\n", infile );
      return(-1);
    }
    
    if ( NULL == ( mask = get_image_mask( inBlock, data, dt, 
                                          lAxes, dss, null, has_null, 
                                          xdesc, ydesc ))){
      /* ERROR ?*/
    }

    xdesc_w = dmDescriptorGetCoord( xdesc );
    ydesc_w = dmDescriptorGetCoord( ydesc );

    descs.image.physical.xaxis = xdesc;
    descs.image.physical.yaxis = ydesc;
    descs.image.world.xaxis = xdesc_w;
    descs.image.world.yaxis = ydesc_w;
    
    

    /* Some error checking on the WORLD coord calcs */
    if ( coordWORLD == ctype ) {
      if ((NULL==xdesc_w)||(NULL==refxdesc_w)) {
        err_msg("ERROR: There must be a WCS on both input "
                "images if coord_sys=world\n");
        return(-1);
      } else {
        char iname[50];
        char rname[50];
        long npar;
        
        dmCoordGetTransformType( xdesc_w, iname, &npar, 49);
        dmCoordGetTransformType( refxdesc_w, rname, &npar, 49 ); 

        /* Not sure if this is really the check to be doing or
           if you should check some other names */
        if ( 0 != ds_strcmp_cis( iname, rname)) {
          err_msg("WARNING: The WCS on the image images should be the same "
                  "type when using coord_sys=wcs (ie TAN-P)");
        }
                
      } /* end else coords */
      
    } /* end if coord */

    hdr[--num_infiles] = getHdr( inBlock, hdrDM_FILE );
  
  
    VertexList refpixlist;   
    VertexList clip_list;
    VertexList tmp_clip_list;

    Polygon ref_poly;
    Polygon clip_poly;
    Polygon tmp_clip_poly;

    double delta;
    
    delta = 1.0 / subpix;
    
    
    /* Set up some data structures and alloc memory */
    refpixlist.vertex = ( Vertex*)calloc(4*subpix, sizeof(Vertex)); /*leak*/
    refpixlist.num_vertices = 4*subpix;
    
    ref_poly.contour = &refpixlist;
    clip_poly.contour = &clip_list;
    tmp_clip_poly.contour = &tmp_clip_list;

    clip_list.vertex = (Vertex*)calloc(4*subpix*4,sizeof(Vertex));
    tmp_clip_list.vertex = (Vertex*)calloc(4*subpix*4,sizeof(Vertex));

    
    /* Okay, we'll look to the input image and map the corners
       to the image in the output.  This way we only have to process those
       pixels in the output which could conver the input.  Should
       save time especially when input image is small piece of the
       output image ... eg mosaics */
    long min_ref_x,min_ref_y;
    long max_ref_x,max_ref_y;
    find_chip_corners( &min_ref_x, &min_ref_y, &max_ref_x, &max_ref_y, 
                 refAxes[0], refAxes[1], lAxes[0], lAxes[1], &descs );

    /* Begin loop through data */

    long xx, yy;
    long ii;

    for(yy=min_ref_y;yy<max_ref_y;yy++) {
      long ypix_off;
      ypix_off = yy*refAxes[0];

      if ( verbose > 2 ) {
        double perc = 100.0 * yy / max_ref_y;
        printf("Processing %g%% complete        \r", perc );
      }

      for (xx=min_ref_x;xx<max_ref_x;xx++ ) {
        long outpix;
        outpix = ypix_off + xx;
        
        if ( subpix == 0 ) {  // TODO DELETE THIS, REMOVE FROM .t FILE
          double val;
          double refpix[2];
          double imgpix[2];        
          long mm,nn;

          refpix[0] = xx+1;
          refpix[1] = yy+1;
          convert_coords( refpix, &descs, imgpix,ctype);
          imgpix[0] -= 0.5; /* = -1 + 0.5 to round */
          imgpix[1] -= 0.5;
          
          mm = imgpix[0]; /* round it */
          nn = imgpix[1];
          
          /* get image value */
          val = get_image_value( data, dt, mm, nn, lAxes, mask );
          
          if ( !ds_dNAN(val) )
            out_data[outpix] += val;
          outpix--;
          continue;
        }


      
        /* This creates a polygon around output pixel xx,yy */
        long xx_min, xx_max, yy_min, yy_max;
        fill_polygon( subpix, delta, xx, yy,&refpixlist, &descs, ctype );
        
        /* find the min and max pixels that need to intersect polygon with */
        xx_min = refpixlist.vertex[0].x-0.5;
        yy_min = refpixlist.vertex[0].y-0.5;
        xx_max = refpixlist.vertex[0].x+0.5;
        yy_max = refpixlist.vertex[0].y+0.5;
        for (ii=4*subpix;--ii;) {
          xx_min = MIN( xx_min, (refpixlist.vertex[ii].x-0.5));
          yy_min = MIN( yy_min, (refpixlist.vertex[ii].y-0.5));
          xx_max = MAX( xx_max, (refpixlist.vertex[ii].x+0.5));
          yy_max = MAX( yy_max, (refpixlist.vertex[ii].y+0.5));
        }

        
        long mm,nn;
        double weight;  
        double sum;
        weight  = 0;
        sum = 0;
        /* for each pixel that could intersect poly check it and get area */
        for (nn=yy_min;nn<=yy_max;nn++) {

          for (mm=xx_min;mm<=xx_max;mm++ ) {
                    
            double area;
            double val;

            /* intersect polygons ; once for each side of the rectangle*/
            poly_clip( &ref_poly,      mm, nn, &tmp_clip_poly ,LEFT );
            poly_clip( &tmp_clip_poly, mm, nn, &clip_poly     ,RIGHT );
            poly_clip( &clip_poly,     mm, nn, &tmp_clip_poly ,BOTTOM );
            poly_clip( &tmp_clip_poly, mm, nn, &clip_poly     ,TOP );

            /* area of polygon */
            area = get_contour_area( &clip_poly );

            /* get image value; weight output pixel by it */
            val = get_image_value( data, dt, mm, nn, lAxes, mask );
            
            if ( !ds_dNAN(val) )
              sum += ( val * area );
            
            weight += area;
            
          } /* end for mm */
        } /* end for nn */
        
        
        /* Can either sum up values or compute the average value; this is done here */
        if ( weight > 0 ) {
          if ( do_norm )
            out_data[outpix] += (sum/weight);
          else
            out_data[outpix] += sum;
        }

        
      } /* end for xx */
    } /* end for yy */

    if ( verbose > 2 )  printf("Processing %g%% complete        \n", 100.0 );


    free(data);
    free(lAxes);
    if ( mask ) free(mask);
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
  
  if ( NULL == ( outBlock = dmImageCreate( outfile, dmDOUBLE,refAxes,2 ))){
      err_msg("ERROR: Cannot create output image '%s'\n", outfile );
      return(-1);
  }
  outDesc = dmImageGetDataDescriptor( outBlock );
  putHdr( outBlock, hdrDM_FILE, outhdr, ALL_STS, "resample_image");
  put_param_hist_info( outBlock, "resample_image", NULL, 0 );
  dmBlockCopyWCS( refBlock, outBlock);
  dmSetArray_d( outDesc, out_data, (refAxes[0]*refAxes[1]));
  dmImageClose(outBlock );
  if ( outBlock != refBlock ) dmImageClose( refBlock );

  return(0);


}
