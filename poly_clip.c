/*                                                                
**  Copyright (C) 2005,2007  Smithsonian Astrophysical Observatory 
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


#include "gpc.h"
#include <stdlib.h>


int poly_clip( gpc_polygon *pp, long qx, long qy,
	       gpc_polygon *ret, short edge);


#ifdef HAS_MAIN

int main() 
{

  gpc_polygon pp;
  gpc_polygon tmp;
  gpc_polygon ret;
  gpc_vertex_list vv_list;
  gpc_vertex *vv;
  int hole = 0;

  pp.num_contours = 1;
  pp.hole = &hole;
  pp.contour = &vv_list;

  vv = (gpc_vertex*)calloc(4,sizeof(gpc_vertex));
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

int poly_clip( gpc_polygon *pp, long qx, long qy,
	       gpc_polygon *ret, short edge)
     
{
  
  double llx,lly, urx,ury;
  long ii;
  int *nout;


  llx = qx -0.5;
  lly = qy -0.5;
  urx = qx +0.5;
  ury = qy +0.5;

  
  nout = &ret->contour[0].num_vertices;
  (*nout) = 0;
  
  for (ii=0;ii<(pp->contour[0].num_vertices);ii++) {
    short p1, p2;
    short jj;
    
    double px_ii, px_jj, py_ii, py_jj;

    p1=0; p2=0;
    
    jj=(ii+1)%pp->contour[0].num_vertices; /* got back to 1st point at ii=num_vertices */
    px_ii = pp->contour[0].vertex[ii].x;
    py_ii = pp->contour[0].vertex[ii].y;
    px_jj = pp->contour[0].vertex[jj].x;
    py_jj = pp->contour[0].vertex[jj].y;
    

    /* These are testing to see if  
       p1 (current point) and p2 (next point) are inside or outside.
    */
    switch (edge) {

    case 1: /* left */
      if ( px_ii >= llx ) 
	p1=1;
      else
	p1=0;
      if ( px_jj >= llx )
	p2=1;
      else
	p2=0;

      break;
    case 2:
      if ( px_ii <= urx ) 
	p1=1;
      else
	p1=0;
      if ( px_jj <= urx )
	p2=1;
      else
	p2=0;
      

      break;
    case 3: /* bottom */
      if ( py_ii >= lly ) 
	p1=1;
      else
	p1=0;
      if ( py_jj >= lly )
	p2=1;
      else
	p2=0;

      break;
    case 4: /* top */
      if ( py_ii <= ury ) 
	p1=1;
      else
	p1=0;
      if ( py_jj <= ury )
	p2=1;
      else
	p2=0;


      break;



    }

    
    
    if ( (1==p1) && (1==p2) ) { /* both points inside, add the jj-th one */
      ret->contour[0].vertex[(*nout)].x = px_jj;
      ret->contour[0].vertex[(*nout)].y = py_jj;
      (*nout)++;
    } else if (( (1==p1) && (0==p2) ) || ( (0==p1) && (1==p2))){  /* one or the other point inside */
      /* check each side, see which it intersects */

      double dom, rr, ss;
      double ax, bx, cx, dx;
      double ay, by, cy, dy;
      
      ax = px_ii;  ay = py_ii;
      bx = px_jj;  by = py_jj;


      /* basically solving eqn for a line; speed up? */

      switch ( edge ) {

      case 1: /* left */
	
	cx = llx;     cy = lly;
	dx = llx;     dy = ury;
	
	dom = (bx-ax)*(dy-cy)-(by-ay)*(dx-cx);
	rr = (ay-cy)*(dx-cx)-(ax-cx)*(dy-cy);
	ss = (ay-cy)*(bx-ax)-(ax-cx)*(by-ay);
	
	if ( 0 == dom ) continue;
	rr /= dom;
	ss /= dom;
	
	ret->contour[0].vertex[(*nout)].x = ax+rr*(bx-ax); /* better be llx! */
	ret->contour[0].vertex[(*nout)].y = ay+rr*(by-ay);
	(*nout)++;
	if (0==p1) {
	  ret->contour[0].vertex[(*nout)].x = px_jj;
	  ret->contour[0].vertex[(*nout)].y = py_jj;
	  (*nout)++;
	}
	continue; /* next ii */

      
      case 2: /* right */
	cx = urx;     cy = lly;
	dx = urx;     dy = ury;
	
	dom = (bx-ax)*(dy-cy)-(by-ay)*(dx-cx);
	rr = (ay-cy)*(dx-cx)-(ax-cx)*(dy-cy);
	ss = (ay-cy)*(bx-ax)-(ax-cx)*(by-ay);
	
	if ( 0 == dom ) continue;
	rr /= dom;
	ss /= dom;
	
	ret->contour[0].vertex[(*nout)].x = ax+rr*(bx-ax); /* better be urx! */
	ret->contour[0].vertex[(*nout)].y = ay+rr*(by-ay);
	(*nout)++;
	if (0==p1) {
	  ret->contour[0].vertex[(*nout)].x = px_jj;
	  ret->contour[0].vertex[(*nout)].y = py_jj;
	  (*nout)++;
	}
	continue; /* next ii */
	
      case 3: /* bottom */
	cx = llx;     cy = lly;
	dx = urx;     dy = lly;
	
	dom = (bx-ax)*(dy-cy)-(by-ay)*(dx-cx);
	rr = (ay-cy)*(dx-cx)-(ax-cx)*(dy-cy);
	ss = (ay-cy)*(bx-ax)-(ax-cx)*(by-ay);
	
	if ( 0 == dom ) continue;
	rr /= dom;
	ss /= dom;
	
	ret->contour[0].vertex[(*nout)].x = ax+rr*(bx-ax);
	ret->contour[0].vertex[(*nout)].y = ay+rr*(by-ay); /* better be lly */
	(*nout)++;
	if (0==p1) {
	  ret->contour[0].vertex[(*nout)].x = px_jj;
	  ret->contour[0].vertex[(*nout)].y = py_jj;
	  (*nout)++;
	}
	continue; /* next ii */

      case 4: /* top */
	cx = llx;     cy = ury;
	dx = urx;     dy = ury;
	
	dom = (bx-ax)*(dy-cy)-(by-ay)*(dx-cx);
	rr = (ay-cy)*(dx-cx)-(ax-cx)*(dy-cy);
	ss = (ay-cy)*(bx-ax)-(ax-cx)*(by-ay);
	
	if ( 0 == dom ) continue;
	rr /= dom;
	ss /= dom;
	
	ret->contour[0].vertex[(*nout)].x = ax+rr*(bx-ax);
	ret->contour[0].vertex[(*nout)].y = ay+rr*(by-ay); /* bette be ury */
	(*nout)++;
	if (0==p1) {
	  ret->contour[0].vertex[(*nout)].x = px_jj;
	  ret->contour[0].vertex[(*nout)].y = py_jj;
	  (*nout)++;
	}
	continue; /* next ii */
      } /* end switch */

    } else if ( (0==p1) && (0==p2) ) {/* neither point inside */
      /* do nothing */
    } 



  } /* end loop over ii */




  return(0);


}
