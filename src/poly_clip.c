/*                                                                
**  Copyright (C) 2005,2007, 2016  Smithsonian Astrophysical Observatory 
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


#include <resample_img.h>
#include <stdlib.h>
     

/*

  Okay, the idea is this:

  For each edge of the rectange we see which points are inside and 
  outside.

  Those points inside are propagated to the next stage.

  Where we have a point inside and the next point in the
  polyon is outside then we have to find the point
  of intersection and that point get added.  (The point
  outside is never added.)

  We do this for all 4 edges and then we have a polygon
  that is the original clipped by a rectangle.

*/

typedef enum { LEFT=1, RIGHT=2, BOTTOM=3, TOP=4 } EdgeType;


static int poly_clip( Polygon *pp, long qx, long qy,
               Polygon *ret, EdgeType edge);

int poly_clip( Polygon *pp, long qx, long qy,
               Polygon *ret, EdgeType edge)
     
{
  
  double llx,lly, urx,ury;
  long ii;
  int *nout;

    /* p1 and p2 describe whether the point is relative to the edge
     * being EdgeType edge being checked (for example , to the left of the left-most
     * edge "MISSES" where as being to the right of the left-most edge "CROSSES".
     * */
     

  enum { MISSES, CROSSES } p1, p2;

  /* qx and qy are the center of the square pixel being checked 
   * 
   * This could be generalized to be a rectangle by passing in 
   * separate lengths for X and Y (or by passing in llx,lly 
   * and urx,ury directly)
   * 
   * */
  llx = qx -0.5;
  lly = qy -0.5;
  urx = qx +0.5;
  ury = qy +0.5;

  
  nout = &ret->contour->num_vertices;
  (*nout) = 0;
  
  for (ii=0;ii<(pp->contour->num_vertices);ii++) {
    short jj;    
    double ax, bx, cx, dx;
    double ay, by, cy, dy;

    p1=MISSES; p2=MISSES;
    cx =0 ; cy = 0; dx = 0; dy = 0;
    jj=(ii+1)%pp->contour->num_vertices; /* got back to 1st point at ii=num_vertices */
    ax = pp->contour->vertex[ii].x;
    ay = pp->contour->vertex[ii].y;
    bx = pp->contour->vertex[jj].x;
    by = pp->contour->vertex[jj].y;
    

    /* These are testing to see if  
       p1 (current point) and p2 (next point) are CROSSES or MISSES.
       
       We go ahead and save the coordinates of the pixel edge as well in
       the cx,cy & dx,dy parameter.
    */
    switch (edge) {

    case LEFT: /* left */
      p1 = ( ax >= llx ) ? CROSSES : MISSES;
      p2 = ( bx >= llx ) ? CROSSES : MISSES;
      cx = llx;     
      cy = lly;
      dx = llx;     
      dy = ury;
      break;
    case RIGHT:
      p1 = ( ax <= urx ) ? CROSSES : MISSES;
      p2 = ( bx <= urx ) ? CROSSES : MISSES;
      cx = urx;     
      cy = lly;
      dx = urx;     
      dy = ury;
      break;
    case BOTTOM: /* bottom */
      p1 = ( ay >= lly ) ? CROSSES : MISSES;
      p2 = ( by >= lly ) ? CROSSES : MISSES;
      cx = llx;
      cy = lly;
      dx = urx;     
      dy = lly;
      break;
    case TOP: /* top */
      p1 = ( ay <= ury ) ? CROSSES : MISSES;
      p2 = ( by <= ury ) ? CROSSES : MISSES;
      cx = llx;     
      cy = ury;
      dx = urx;     
      dy = ury;
      break;
    }

    
    /* check each side, see which it intersects */
    
    if ( (CROSSES==p1) && (CROSSES==p2) ) { /* both points CROSSES, add the jj-th one */
      ret->contour->vertex[(*nout)].x = bx;
      ret->contour->vertex[(*nout)].y = by;
      (*nout)++;
      continue;
    } else if ( (MISSES==p1) && (MISSES==p2) ) {/* neither point CROSSES */
      /* do nothing */
      continue;

    } else if (( (CROSSES==p1) && (MISSES==p2) ) || ( (MISSES==p1) && (CROSSES==p2))){  /* one or the other point CROSSES */

    /* If one point is CROSSES and the other is not, then we add 
     * the point where the line segment between the two points intersects
     * the edge.
     * */
     

      double dom, rr, ss;

      /* basically solving eqn for a line; speed up? */

      dom = (bx-ax)*(dy-cy)-(by-ay)*(dx-cx);
      rr  = (ay-cy)*(dx-cx)-(ax-cx)*(dy-cy);
      ss  = (ay-cy)*(bx-ax)-(ax-cx)*(by-ay);
        

      /* dom is demonimator, if 0, then lines are parallel, no intersection */
      if ( 0 == dom ) continue;
      rr /= dom;
      ss /= dom;
        
      ret->contour->vertex[(*nout)].x = ax+rr*(bx-ax); 
      ret->contour->vertex[(*nout)].y = ay+rr*(by-ay);
      (*nout)++;
      if (CROSSES==p2) { // if p2 is CROSSES (means p1 is MISSES) so add p2
        ret->contour->vertex[(*nout)].x = bx;
        ret->contour->vertex[(*nout)].y = by;
        (*nout)++;
      }


    } // end else if p1 || p2 but not both


  } /* end loop over ii */


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






#ifdef HAS_MAIN

// Test routine only

int main() 
{

  Polygon pp;
  Polygon tmp;
  Polygon ret;
  VertexList vv_list;
  Vertex *vv;
  int hole = 0;

  pp.num_contours = 1;
  pp.hole = &hole;
  pp.contour = &vv_list;

  vv = (Vertex*)calloc(4,sizeof(Vertex));
  vv_list.num_vertices = 4;
  vv_list.vertex = vv;
  
  vv[0].x = 9.0;  vv[0].y = 9.0;
  vv[1].x = 10.0; vv[1].y = 9.0;
  vv[2].x = 10.0; vv[2].y = 10.0;
  vv[3].x = 9.0;  vv[3].y = 10.0;
  
  poly_clip( &pp, 10, 10, &ret, 1 );
  poly_clip( &ret, 10, 10, &tmp, 2 );
  /* free */
  poly_clip( &tmp, 10, 10, &ret, 3 );
  poly_clip( &ret, 10, 10, &tmp, 4 );
  

  return(0);

}

#endif


