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

 
/* H*****************************************************************
 *
 * FILE NAME: t_repro_img.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * IRAF wrapper for repro_img tool
 *  
 *
 *
H***************************************************************** */

#include "dslib.h"
extern int repro_img(void);


int main(int argc, char** argv)
{
    int fail_status = 0; 


#ifdef _RI_GRID_
    dsErrInitLib(dsPTGRPERR, "reproject_image_grid");
#else
    dsErrInitLib(dsPTGRPERR, "reproject_image");
#endif

    /* INIT THE IRAF ENVIRONMENT */
    /* IRAFHOST_INIT;		*/

    /* OPEN THE PARAMETER FILE */
    if(clinit(argv, argc, "rw") == NULL)
    {
       err_msg( "Problem opening parameter file %s.par\n", argv[0]);
       err_msg( "Parameter library error: %s.\n", paramerrstr());
       fail_status = -1;
    }
    else
    {    
       /* EXECUTE OUR PROGRAM */ 
       fail_status = repro_img();
    
       /* CLOSE PARAMETER FILE AND RETURN TO THE OS */
       clclose();
    } 
    
    exit(fail_status); 
}
