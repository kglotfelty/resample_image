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

#include "gpc.h"

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
 WCS ref; WCS image; short swap_phys; short swap_wrld;
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
double convert_hms_dms( char *radec, double factor );
double convert_pixsize( char *psize );
double get_contour_area( gpc_polygon *clip_poly );
int fill_polygon( long subpix,
                  double delta,
                  long xx, long yy,
                  gpc_vertex_list *refpixlist,
                  WCS_Descriptors *descs,
                  CoordType ctype);
                  
int convert_coords( double *refpix,
                    WCS_Descriptors *descs,
                    double *imgpix,
                    CoordType ctype);

                       

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




extern int poly_clip( gpc_polygon *pp, long qx, long qy,
               gpc_polygon *ret, short edge);




/* Routine to convert coordinates */

/* Note: 
   Thereis a known bug in that if the input image has RA,DEC and the matchfile has DEC,RA
   the output image will be bogus.  It'll be in this routine where we'll need to
   swap the axes at the appropriate point in the transform.

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
                  gpc_vertex_list *refpixlist,
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

double get_contour_area( gpc_polygon *clip_poly )
{
  double area;
  long kk;

  area = 0;
  for (kk=0;kk<clip_poly->num_contours;kk++) {
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
    
    if ( clip_poly->hole[kk] == 0 )
      area += fabs(larea);
    else {
      area -= fabs(larea);
    }
  }

  return(area);
}


/*
  Convert arcmin/sec to degress based on the pixel size.
*/

double convert_pixsize( char *psize )
{

  char *ptr;
  char *ptr2;
  double retval;

  ptr = psize + strlen(psize);

  retval = strtod( psize, &ptr2 );
  if (( ptr2!=NULL) && ( strlen(ptr2)>0)) {
    if (ptr2==(ptr-1)) {
      if (*ptr2 == '"' )
        retval /= 3600.0;
      else if ( *ptr2 == '\'' )
        retval /= 60.0;
      else
        ds_MAKE_DNAN(retval);
    } else {
      ds_MAKE_DNAN(retval);
    }
  }

  return(retval);
}

/*
  Convert a string in HMS or DMS format to degrees. (adapted from dmcoords)
*/
double convert_hms_dms( char *radec, double factor )
{
  double retval;
  char *ptr;

  retval = strtod( radec, &ptr );
  if (( ptr != NULL )&&( strlen(ptr)>0)) {

    int hh, mm;
    float sec;
    char *ptr2;

    hh = strtol( radec, &ptr, 10 );

    if (( ptr == NULL ) || ( *ptr != ':' ) ) {
      ds_MAKE_DNAN( retval);
      return(retval);
    }
    ptr2 = ptr+1;
    mm = strtol( ptr2, &ptr, 10 );
    if (( ptr == NULL ) || ( *ptr != ':' ) ) {
      ds_MAKE_DNAN( retval);
      return(retval);
    }
    ptr2 = ptr+1;
    sec = strtod( ptr2, &ptr );
    if (( ptr != NULL ) && (strlen(ptr) > 0 )) {
      ds_MAKE_DNAN( retval);
      return(retval);
    }
    
    if ( ( sec < 0 ) || ( sec >= 60 ) ) {
      ds_MAKE_DNAN( retval);
      return(retval);
    }
    if ( ( mm < 0 ) || ( mm >= 60 ) ) {
      ds_MAKE_DNAN( retval);
      return(retval);
    }

    if ( factor == 1 ) {
      if ( (hh < -90) || ( hh > 90 ) ) {
        ds_MAKE_DNAN( retval);
        return(retval);
      }
    } else {
      if ( (hh < 0) || ( hh > 24)) {
        ds_MAKE_DNAN( retval);
        return(retval);
        
      }
      
    }
    
    /*retval = (((hh * 60.0 + mm) * 60.0) + sec ) / 3600.0;*/


    int sgn = ( hh < 0 ) ? -1 : 1; /* get sign of hour*/

    hh *= sgn;  /* make positive */
    
    retval = (hh + (mm/60.0) + ( sec/3600.0)) * sgn;
    
    retval *= factor;
    
  }

  return(retval);
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
    vals_x[0] = -2; vals_y[0] = -2;
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
          return(-1);
          
        }


        *min_ref_x = MIN( *min_ref_x, (ref_log_loc[0]-1));
        *min_ref_y = MIN( *min_ref_y, (ref_log_loc[1]-1));
        *max_ref_x = MAX( *max_ref_x, (ref_log_loc[0]-1));
        *max_ref_y = MAX( *max_ref_y, (ref_log_loc[1]-1));
        
      }
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

int repro_img(void);

int repro_img(void)
{

  char instack[DS_SZ_PATHNAME];  /* Input stack of images */
  char *infile;                  /* individual image in stack */
  char reffile[DS_SZ_PATHNAME];  /* match file */
  char outfile[DS_SZ_PATHNAME];  /* output file */
  short clobber;
  short verbose;

  Stack inStack;

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

  dmBlock *refBlock;
  dmDescriptor *refxdesc, *refydesc;
  dmDescriptor *refxdesc_w, *refydesc_w;
  long *refAxes;
  
  dmBlock *outBlock = NULL;
  dmDescriptor *outDesc;

  long xx, yy;
  long ii;

  int subpix;
  double delta;

  /* These gpc_* data structures are defined in gpc.h*/
  /* 
     They are a left over from the original prototype
     that used a general polygon clip algorithm but
     now we use a more specialized rectangle-clip-polygon
     algorithm.  However, they are somewhat convienient so
     we kept them as is.
  */

  gpc_vertex_list imgpixlist;
  gpc_vertex_list refpixlist;   

  gpc_vertex_list clip_list;
  gpc_vertex_list tmp_clip_list;

  gpc_polygon img_poly;
  gpc_polygon ref_poly;

  gpc_polygon clip_poly;
  gpc_polygon tmp_clip_poly;


  int hole = 0;

  
  long xx_min, xx_max, yy_min, yy_max;
  
  double *out_data;
  long outpix;
  short do_norm;
  char which_norm[30];
  char csys[30];

  CoordType ctype = coordWORLD ;

  short using_refimg=1;

  double crval[2];
  double crpix[2];
  double cdelt[2];  
  long   naxes[2];


  /* This .c file is used to compile two different tools:
     reproject_image and reproject_image_grid.

     The later will have -D_RI_GRID on the compile line to compile the 
     minimally different code between the two (all to
     do with the user interface
  */
  
  char   wrld_vals[2][100];
  char   psize[50];
  
  long num_infiles;
  Header_Type **hdr;
  Header_Type *outhdr;
  char lookup[DS_SZ_PATHNAME];

  WCS_Descriptors descs;


  /* Get the parameters */
  clgetstr( "infile", instack, DS_SZ_FNAME );
  clgetstr( "matchfile", reffile, DS_SZ_FNAME );
  clgetstr( "outfile", outfile, DS_SZ_FNAME );


  if ( ( strlen( reffile) == 0 ) ||
       ( ds_strcmp_cis(reffile, "none" ) == 0 ) ) {

    err_msg("ERROR: Must supply a valid match file\n");
    return(-1);
  } else {
    using_refimg = 1;
  }


  subpix = clgeti("resolution");

  clgetstr("method", which_norm, 29);
  if ( strcmp( which_norm, "sum" ) == 0 ) 
    do_norm=0;
  else
    do_norm=1;
  clgetstr("coord_sys", csys, 29);
  switch ( csys[0] ) {
  case 'l': ctype = coordLOGICAL; break;
  case 'p': ctype = coordPHYSICAL; break;
  case 'w': ctype = coordWORLD; break;
  default:
    err_msg("ERROR: Unknow coordinate type '%s'\n", csys );
    return(-1);
  }

  clgetstr( "lookupTab", lookup, DS_SZ_FNAME);
  clobber = clgetb( "clobber" );
  verbose = clgeti( "verbose" );
  

  if ( verbose ) {  
    printf("reproject_image - parameters\n");
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






  if ( ds_clobber( outfile, clobber, NULL ) != 0 ) {
    return(-1);
  }


  if (  0 == using_refimg ) {

    char *skyname = "SKY";
    char *xynames[2] = { "X", "Y" };
    char *wrldname = "EQPOS";
    char *radecnames[2] = { "RA", "DEC"};

    refAxes = (long*)calloc(2,sizeof(long));
    refAxes[0] = naxes[0];
    refAxes[1] = naxes[1];
    if ( NULL == ( outBlock = dmImageCreate( outfile, dmDOUBLE,refAxes,2 ))){
      err_msg("ERROR: Cannot create output image '%s'\n", outfile );
      return(-1);
    }
    
    cdelt[1] = convert_pixsize( psize );
    if ( ds_dNAN( (cdelt[1])) ) {
      err_msg("ERROR: Cannot convert pixelsize='%s'\n", psize );
      return(-2);
    }
    cdelt[0] = -cdelt[1]; /* East is to the left */
    
    
    crval[0] = convert_hms_dms( wrld_vals[0], 15.0 );
    if ( ds_dNAN( (crval[0]))) {
      err_msg("ERROR: cannot convert xcenter value='%s'\n", wrld_vals[0]);
      return(-2);
    }
    
    crval[1] = convert_hms_dms( wrld_vals[1], 1.0 );
    if ( ds_dNAN( (crval[1]))) {
      err_msg("ERROR: cannot convert xcenter value='%s'\n", wrld_vals[1]);
      return(-2);
    }
    
   
    refBlock = outBlock;
    outDesc = dmImageGetDataDescriptor( outBlock );
    
    refxdesc = dmArrayCreateAxisGroup( outDesc, skyname, dmDOUBLE,
                                       "pixel", xynames, 2 );
    refydesc = NULL;
    /* do not use lowercase 'tan' ... bad results */
    refxdesc_w = dmCoordCreate_d( refxdesc, wrldname, "degree",
                                  radecnames, 2, "TAN",
                                  crpix, crval, cdelt, NULL);
    refydesc_w = NULL;
    
  } else {
    /* do ref image first */
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

    free(data);


  } /* end else has refimg */


  /* Now let's start on the input stack */

  inStack = stk_build( instack );
  if ( ( NULL == inStack ) ||
       ( stk_count(inStack)==0 ) ||
       (( stk_count(inStack)==1 ) && ( strlen(stk_read_num(inStack,1))==0))) {
    err_msg("ERROR: problems opening stack '%s'\n", instack );
    return(-3);
  }


  num_infiles = stk_count(inStack);
  hdr = (Header_Type**) calloc( num_infiles, sizeof(Header_Type*));
  
  /* Okay, start to process the data */

  out_data = ( double*)calloc( refAxes[1]*refAxes[0], sizeof(double));
  

  stk_rewind(inStack);
  while ( NULL != (infile = stk_read_next(inStack))) {

    long min_ref_x,min_ref_y;
    long max_ref_x,max_ref_y;

    outpix = refAxes[1]*refAxes[0];

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
    
    

    /* 
       We want the output image to have a physical coordinate system
       that is _consistent_ with the 1st input file.  They will not be identical.
       Typically the physical pixel is defined at the bottom left corner.
       This maps the center of the field/tan-point/etc. to physical and uses
       that as the transform values.
    */
    if (( 0 == using_refimg) && ( 1 == stk_current( inStack )) ) {
      double ppix[2];
      double pval[2];
      double pdlt[2];
      short dim = ( ydesc ? 1 : 2 );
      double phyinp[2];
      
      
      dmCoordGetTransform_d( xdesc, ppix, pval, pdlt, dim );
      if ( ydesc )
        dmCoordGetTransform_d( ydesc, ppix+1, pval+1, pdlt+1, dim );

      cdelt[0] /= pdlt[0]; /* New cdelt values */
      cdelt[1] /= pdlt[1];
      
      /* Okay, we map the crval to the physical pixel in the input image */
      dmCoordInvert_d( xdesc_w, crval, phyinp); 
      if ( ydesc )
        dmCoordInvert_d( ydesc_w, crval+1, phyinp+1); 


      ppix[0] = crpix[0];
      ppix[1] = crpix[1];
      
      dmCoordSetTransform_d( refxdesc, ppix, phyinp, pdlt, 2 );
      dmCoordSetTransform_d( refxdesc_w, phyinp, crval, cdelt,2 );
      
    }

    descs.swap_phys = 0;
    descs.swap_wrld = 0;

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

        dmGetCptName( xdesc_w, 1, iname, 49 );
        dmGetCptName( refxdesc_w, 1, rname, 49 );
        
        if ( 0 != ds_strcmp_cis( iname, rname ) ) 
          descs.swap_wrld = 1;

        dmGetCptName( xdesc, 1, iname, 49 );
        dmGetCptName( refxdesc, 1, rname, 49 );
        
        if ( 0 != ds_strcmp_cis( iname, rname ) ) 
          descs.swap_phys = 1;


        
      } /* end else coords */
      
    } /* end if coord */

    hdr[--num_infiles] = getHdr( inBlock, hdrDM_FILE );
  
    
    
    delta = 1.0 / subpix;
    
    
    /* Set up some data structures and alloc memory */
    refpixlist.vertex = ( gpc_vertex*)calloc(4*subpix, sizeof(gpc_vertex)); /*leak*/
    refpixlist.num_vertices = 4*subpix;
    
    imgpixlist.vertex = ( gpc_vertex*)calloc(4,sizeof(gpc_vertex)); /* leak*/
    imgpixlist.num_vertices = 4;
    
    img_poly.contour = &imgpixlist;
    img_poly.hole = &hole;
    img_poly.num_contours = 1;
    
    ref_poly.contour = &refpixlist;
    ref_poly.hole = &hole;
    ref_poly.num_contours = 1;
    
    clip_poly.contour = &clip_list;
    clip_poly.hole = &hole;
    clip_poly.num_contours = 1;

    tmp_clip_poly.contour = &tmp_clip_list;
    tmp_clip_poly.hole = &hole;
    tmp_clip_poly.num_contours = 1;

    clip_list.vertex = (gpc_vertex*)calloc(4*subpix*4,sizeof(gpc_vertex));
    tmp_clip_list.vertex = (gpc_vertex*)calloc(4*subpix*4,sizeof(gpc_vertex));

    /*
      Put all this stuff into a struct to make it faster to pass around 
    */
    descs.ref.physical.xaxis = refxdesc;
    descs.ref.physical.yaxis = refydesc;
    descs.ref.world.xaxis = refxdesc_w;
    descs.ref.world.yaxis = refydesc_w;
    descs.image.physical.xaxis = xdesc;
    descs.image.physical.yaxis = ydesc;
    descs.image.world.xaxis = xdesc_w;
    descs.image.world.yaxis = ydesc_w;
    
    /* Okay, we'll look to the input image and map the corners
       to the image in the output.  This way we only have to process those
       pixels in the output which could conver the input.  Should
       save time especially when input image is small piece of the
       output image ... eg mosaics */
    if ( 0 != find_chip_corners( &min_ref_x, &min_ref_y, 
                                 &max_ref_x, &max_ref_y, 
                                 refAxes[0], refAxes[1],
                                 lAxes[0], lAxes[1],
                                 &descs )) {

      min_ref_x = 0;
      min_ref_y = 0;
      max_ref_x = refAxes[0];
      max_ref_y = refAxes[1];
      /* WARNING? */
    }

    /* Begin loop through data */
    outpix--; /* output pixel counter */
    for(yy=min_ref_y;yy<max_ref_y;yy++) {
      double ypix_off;
      ypix_off = yy*refAxes[0];

      if ( verbose > 2 ) {
        double perc = 100.0 * yy / max_ref_y;
        printf("Processing %g%% complete        \r", perc );
      }

      for (xx=min_ref_x;xx<max_ref_x;xx++ ) {
        double val;
        double refpix[2];
        double imgpix[2];
        
        long mm,nn;
        double weight;  
        double sum;

        outpix = ypix_off + xx;
        
        if ( subpix == 0 ) {
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
        
        
        weight  = 0;
        sum = 0;
        /* for each pixel that could intersect poly check it and get area */
        for (nn=yy_min;nn<=yy_max;nn++) {
          for (mm=xx_min;mm<=xx_max;mm++ ) {
                    
            double area;
            /* imgpixlist is data under img_poly */
            imgpixlist.vertex[0].x=mm-0.5;
            imgpixlist.vertex[1].x=mm-0.5;
            imgpixlist.vertex[2].x=mm+0.5;
            imgpixlist.vertex[3].x=mm+0.5;
            imgpixlist.vertex[0].y=nn-0.5;
            imgpixlist.vertex[1].y=nn+0.5;
            imgpixlist.vertex[2].y=nn+0.5;
            imgpixlist.vertex[3].y=nn-0.5;


            /* intersect polygons ; once for each side of the rectangle*/
            poly_clip( &ref_poly,      mm, nn, &tmp_clip_poly ,1 );
            poly_clip( &tmp_clip_poly, mm, nn, &clip_poly     ,2 );
            poly_clip( &clip_poly,     mm, nn, &tmp_clip_poly ,3 );
            poly_clip( &tmp_clip_poly, mm, nn, &clip_poly     ,4 );

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

        outpix--;
        
      } /* end for xx */
    } /* end for yy */

    if ( verbose > 2 )  printf("Processing %g%% complete        \n", 100.0 );


    free(data);
    free(lAxes);
    if ( mask ) free(mask);
    dmImageClose( inBlock );

  }  /* end while infile loop over stack */

  num_infiles = stk_count( inStack );

  /* merge headers */
  if (( (strlen( lookup ) == 0 ) ||
        (ds_strcmp_cis(lookup, "none") ==0 )) ||
      ( num_infiles == 1 ) ) {
    outhdr = hdr[num_infiles-1];
  } else {
    outhdr = mergeHdr( lookup, hdr, num_infiles );
  }
  

  /* create output image if necessary */
  if ( 1 == using_refimg ) {
    if ( NULL == ( outBlock = dmImageCreate( outfile, dmDOUBLE,refAxes,2 ))){
      err_msg("ERROR: Cannot create output image '%s'\n", outfile );
      return(-1);
    }
  }
  outDesc = dmImageGetDataDescriptor( outBlock );
  
  /*   dmBlockCopy( inBlock, outBlock, "HEADER"); */
  putHdr( outBlock, hdrDM_FILE, outhdr, ALL_STS, "reproject_image");
  put_param_hist_info( outBlock, "reproject_image", NULL, 0 );

  if ( 1 == using_refimg ) dmBlockCopyWCS( refBlock, outBlock);
  
  dmSetArray_d( outDesc, out_data, (refAxes[0]*refAxes[1]));
  dmImageClose(outBlock );
  if ( outBlock != refBlock ) dmImageClose( refBlock );

  return(0);


}
