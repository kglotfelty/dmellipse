/*                                                                
**  Copyright (C) 2004-2008,2011,2014-2015  Smithsonian Astrophysical Observatory 
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

/* ------------------------------------ */

#include "ellipse.h"
#include <dslib.h>
#include "dmimgio.h"
#include <stack.h> 
#include <histlib.h>

#include <dsnan.h>
#include <math.h>
#include <float.h>



/*Convert input C array index i,j into physical coords */
int convert_coords( Image *img, long xx, long yy, double *xat, double *yat )
{
  
  if ( img->xdesc ) {
    double lgc[2];
    double phy[2];
    lgc[0] = xx+1;  // C index starts at 0, logical coords start at 1
    lgc[1] = yy+1;
    dmCoordCalc_d( img->xdesc, lgc, phy );
    if ( img->ydesc ) {
      dmCoordCalc_d( img->ydesc, lgc+1, phy+1 );
    }
    *xat = phy[0];
    *yat = phy[1];
  } else {
    *xat = xx;
    *yat = yy;
  }

  return(0); // not checked
}



/* Compute the moments from a 2D array (image) */
int get_moments_img( Image *img, Moments *mom )
{
  double sum_vals;
  int numvals;

  long xx, yy;

  numvals = 0;
  sum_vals = 0;

  /* init to 0 */
  mom->x_mu = 0;
  mom->y_mu = 0;
  mom->out[0][0] = mom->out[0][1] = mom->out[0][2] = 0;
  mom->out[1][0] = mom->out[1][1] = mom->out[1][2] = 0;
  mom->out[2][0] = mom->out[2][1] = mom->out[2][2] = 0;

  for (yy=img->lAxes[1]; yy--; ) {
    for ( xx=img->lAxes[0]; xx--; ) {
      double val;
      double xat, yat;
      val = get_image_value( img->data, img->dt, xx, yy, img->lAxes, img->mask );
      if ( ds_dNAN(val ) ) continue;
      
      val = fabs(val);
      numvals++;
      convert_coords( img, xx, yy, &xat, &yat );

      sum_vals += val;
      mom->x_mu += xat * val;
      mom->y_mu += yat * val;

    } /* end for loop over xx */
  } /* end for loop over yy */

  if ( 0 == sum_vals ) return(0);

  mom->x_mu /= sum_vals;  // X centroid
  mom->y_mu /= sum_vals;  // Y centroid


  for (yy=img->lAxes[1]; yy--; ) {
    double yy_qq;
    
    for ( xx=img->lAxes[0]; xx--; ) {
      double val;
      double xx_pp;
      double xat, yat;
      
      long ii, jj;

      val = get_image_value( img->data, img->dt, xx, yy, img->lAxes, img->mask );
      if ( ds_dNAN(val ) ) continue;
      val = fabs(val);

      convert_coords( img, xx, yy, &xat, &yat );

      // 2nd order moment matrix
      for ( ii=0;ii<3;ii++ ) {
        for ( jj=0;jj<3;jj++ ) {
          xx_pp = pow((xat-(mom->x_mu)),ii);
          yy_qq = pow((yat-(mom->y_mu)),jj);
          
          mom->out[ii][jj] += ( xx_pp * yy_qq * val ); 
        } // end for jj
      } // end for ii

    } /* end for loop over xx */
  } /* end for loop over yy */

  compute_ellipse_vals_from_moment( mom, sum_vals );

  return(numvals);
}


/* Similar to above; compute moments from and array of X, Y, and Z values */
int get_moments_table( Table *tab, Moments *mom )
{

  double sum_vals;
  int numvals;

  long nn;

  numvals = 0;
  sum_vals = 0;

  // Init to 0
  mom->x_mu = 0;
  mom->y_mu = 0;
  mom->out[0][0] = mom->out[0][1] = mom->out[0][2] = 0;
  mom->out[1][0] = mom->out[1][1] = mom->out[1][2] = 0;
  mom->out[2][0] = mom->out[2][1] = mom->out[2][2] = 0;

  for (nn=tab->nrows; nn--; ) {
        double val;

        val = tab->vals[nn]; 
        if ((0 == tab->mask[nn]) || (ds_dNAN(val) )) continue;
        val = fabs(val);
        numvals++;

        sum_vals += val;
        mom->x_mu += tab->xx[nn] * val;
        mom->y_mu += tab->yy[nn] * val;
  }

  if ( 0 == sum_vals ) return(0);

  mom->x_mu /= sum_vals;  // X centroid
  mom->y_mu /= sum_vals;  // Y centroid

  for ( nn = tab->nrows; nn--; ) {
      double val;
      double xx_pp, yy_qq;
      
      long ii, jj;

      val = tab->vals[nn];
      if ((0 == tab->mask[nn]) || (ds_dNAN(val) )) continue;
      val = fabs(val);
 
      // 2nd order 3x3 moments matrix
      for ( ii=0;ii<3;ii++ ) {
        for ( jj=0;jj<3;jj++ ) {
          xx_pp = pow((tab->xx[nn]-(mom->x_mu)),ii);
          yy_qq = pow((tab->yy[nn]-(mom->y_mu)),jj);
          mom->out[ii][jj] += ( xx_pp * yy_qq * val );
        }
      }

  } /* end for loop over nn*/

  compute_ellipse_vals_from_moment( mom, sum_vals );

  return(numvals);
}



/* Compute the ellipse parameters (major axis, minor axis, rotation angle) from
 * the 2nd order moments matrix */
void compute_ellipse_vals_from_moment( Moments *mom, double sum_vals )
{

  double mm_0_2, mm_2_0, mm_1_1;
  double mm_diag_diff;
  double mm_0_2_sq, mm_2_0_sq, mm_1_1_sq;
  
  mm_2_0 = mom->out[2][0]/sum_vals;
  mm_0_2 = mom->out[0][2]/sum_vals;
  mm_1_1 = mom->out[1][1]/sum_vals;

  mm_2_0_sq = mm_2_0 * mm_2_0;
  mm_0_2_sq = mm_0_2 * mm_0_2;
  mm_1_1_sq = mm_1_1 * mm_1_1;

  // "xsig" is major axis, "ysig" is minor axis
  mm_diag_diff = mm_2_0 - mm_0_2;
  mom->xsig = ( mm_2_0 + mm_0_2 ) + sqrt ( mm_diag_diff * mm_diag_diff + 
                                           4 * mm_1_1 * mm_1_1 );
  mom->xsig /= 2.0;
  mom->xsig = 2.0 * sqrt( mom->xsig );

  mom->ysig = ( mm_2_0 + mm_0_2 ) - sqrt ( mm_diag_diff * mm_diag_diff + 
                                           4 * mm_1_1 * mm_1_1 );
  mom->ysig /= 2.0;
  mom->ysig = 2.0 * sqrt( mom->ysig );

  // angle in degrees, CCW from +X axis to major axis
  if ( fabs(mm_diag_diff) > FLT_EPSILON ) {
    mom->phi = 0.5 * (180.0/3.141592)*atan2( (2.0 * mm_1_1) , mm_diag_diff );
  } else {
    mom->phi = 90;
  }

  if ( mom->ysig > FLT_EPSILON ) {
    mom->ell = mom->xsig / mom->ysig;
  } else {
    mom->ell = 1.0;
  }


}



/*
  Set the mask pixels that are outside the current ellipse/shape defined
  by the moments (mom) at the given radius (offset) to 0.
*/
void make_ellipse_mask_img( Image *img,
                        Moments *mom, 
                        float offset,
                        char *shape)
{
  long xx, yy;

  regRegion *ellipse;
  char ellstr[1000];
  double r0, r1;

  r0 = offset;          // major axis
  r1 = r0 / mom->ell;   // minor axis
  
  sprintf(ellstr, "%s(%g,%g,%g,%g,%g)", shape,
          mom->x_mu, mom->y_mu, r0, r1, mom->phi );
  ellipse = regParse(ellstr);

  for (yy=img->lAxes[1]; yy--; ) {
    for ( xx=img->lAxes[0]; xx--; ) {
      double xat, yat;

      if ( img->mask[xx+yy*img->lAxes[0]] ) {
        convert_coords( img, xx, yy, &xat, &yat );
        if ( !regInsideRegion( ellipse, xat, yat ) )
          img->mask[xx+yy*img->lAxes[0]] = 0;
      } // end if mask

    } // end for xx
  }  // end for yy

  regFree(ellipse);

}



/* Similar to above but for rows in a table */
void make_ellipse_mask_table( Table *tab,
                        Moments *mom, 
                        float offset,
                        char *shape)
{
  long nn;

  regRegion *ellipse;
  char ellstr[1000];
  double r0, r1;

  r0 = offset;
  r1 = r0 / mom->ell;
  
  sprintf(ellstr, "%s(%g,%g,%g,%g,%g)", shape,
          mom->x_mu, mom->y_mu, r0, r1, mom->phi );
  ellipse = regParse(ellstr);


  for (nn = tab->nrows; nn--; ) {
      if (tab->mask[nn]) {
         if ( !regInsideRegion( ellipse, tab->xx[nn], tab->yy[nn] ) )
            tab->mask[nn] = 0;
      }
  } //end for nn

  regFree(ellipse);

}

/* Set the mask image back to original values (image)*/
void reset_mask_img(Image *img, short *mask)
{
  long ii;
  for (ii=img->lAxes[1]*img->lAxes[0];ii--; ) {
    if (mask) {
      img->mask[ii] = mask[ii];
    } else {
      img->mask[ii] = 1;
    }
  }
}


/* Set the mask image back to original values (table)*/
void reset_mask_table(Table *tab, short *mask)
{
  long ii;
  for (ii=tab->nrows;ii--; ) {
    if (mask) {
      tab->mask[ii] = mask[ii];
    } else {
      tab->mask[ii] = 1;
    }
  }
}


/* Parse/stack-expand the grid of fractions to compute */
long parse_grid( char *grid, double **vals )
{
  long nval;
  Stack _gs;
  long ii;
  nval = 0;

  if (( NULL == grid ) || ( 0 == strlen(grid))) {
    err_msg("ERROR: Empty grid\n");
    return(nval);
  }

 
  _gs = stk_build(grid);
  if (( NULL == _gs ) ||
      ( 0 == stk_count(_gs) ) ||
      ( ( 1 == stk_count(_gs) && ( 0 == strlen(stk_read_num(_gs,1))))) ){
    err_msg("ERROR: Could not stack expand grid='%s'\n", grid );
    return(nval);
  }

  stk_rewind(_gs);
  nval = stk_count(_gs);
  *vals = (double*)calloc(nval, sizeof(double));
  for (ii=1;ii<=nval;ii++ ) {
    char *val = stk_read_num(_gs, ii);
    char *ptr;
    double dval;
    
    dval = strtod(val, &ptr );
    if ( (ptr ) && (ptr[0] != 0 )) {
      err_msg("ERROR: Could not parse '%s'\n", val );
      return(0);
    }
    free(val);
    (*vals)[ii-1] = dval;
    
  } /* end for ii */

  return(nval);
}



/* Load image file.  This uses the routines defined in dmimgio.h */
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


/* Load a table.  3 columns are loaded.  If the zcolumn is blank, then
 * use the value 1.0 for Z.  (Could also opt to use 1.0/N). */
Table *load_table_file( dmBlock *inBlock, char *xcol, char *ycol, char *zcol )
{

  long ii;

  Table *tab = (Table*)calloc(1,sizeof(Table));

  tab->nrows = dmTableGetNoRows(inBlock );

  dmDescriptor *xd = dmTableOpenColumn( inBlock, xcol );
  if ( NULL == xd ) {
    err_msg("ERROR: Cannot find column '%s' in input file", xcol );
    return(NULL);
  }
  tab->xx = (double*)calloc(tab->nrows, sizeof(double));
  dmGetScalars_d( xd, tab->xx, 1, tab->nrows);

  xd = dmTableOpenColumn( inBlock, ycol );
  if ( NULL == xd ) {
    err_msg("ERROR: Cannot find column '%s' in input file", ycol );
    return(NULL);
  }
  tab->yy = (double*)calloc(tab->nrows, sizeof(double));
  dmGetScalars_d( xd, tab->yy, 1, tab->nrows);


  tab->vals = (double*)calloc(tab->nrows, sizeof(double));
  tab->mask = (short*)calloc(tab->nrows, sizeof(short));

  if (strlen(zcol) > 0 ) {
       xd = dmTableOpenColumn( inBlock, zcol );
       if ( NULL == xd ) {
         err_msg("WARNING: Cannot find column '%s' in input file.  Using 1's", zcol );
         for (ii=tab->nrows; ii--; ) {
            tab->vals[ii] = 1.0;
         }
         for (ii=tab->nrows; ii--; ) {
             tab->mask[ii] = 1;
         }
       } else {
           dmGetScalars_d( xd, tab->vals, 1, tab->nrows);
           get_table_mask( xd, tab );
       }
   } else { /* No z column */
        xd = NULL;
        for (ii=tab->nrows; ii--; ) {
            tab->vals[ii] = 1.0;
        }
        for (ii=tab->nrows; ii--; ) {
             tab->mask[ii] = 1;
        }
   }


  purge_duplicates_in_tab( tab );

  return(tab);
}



 
void purge_duplicates_in_tab( Table *tab )
{    
    /*
     * Multiple rows with the same x,y in the table can slow things down.
     * 
     * Collapse the same x,y into a value with the sum() of the vals's
     * 
     * This has an interesting side-effect w.r.t. accumulated
     * error in float-point values when computing the moments -- it seems
     * to make things more stable.  
     * 
     * */
     
    
    long ii;
    long jj;
    
    long n_minus_1 = tab->nrows-1;
    long nout = 0;
    
    /* Identify duplicate (x,y) */
    
    for (ii=0;ii< n_minus_1; ii++) {        
        for (jj=ii+1; jj< tab->nrows; jj++) {
            /* Could do an epsilon check here */
            if (( tab->xx[ii] == tab->xx[jj] ) && ( tab->yy[ii] == tab->yy[jj] )) {

                tab->vals[ii] += tab->vals[jj];  // Increment old value
                tab->mask[jj] = 0;   // Set j-th to be bad.

            }
        } // end for jj
    } // end for ii
    


    /* Repack the table values purging duplicates */
    
    double *xx, *yy, *vals;
    short *mask;
    xx = (double*)calloc(tab->nrows, sizeof(double));
    yy = (double*)calloc(tab->nrows, sizeof(double));
    vals = (double*)calloc(tab->nrows, sizeof(double));
    mask = (short*)calloc(tab->nrows, sizeof(short));
    nout = 0;
    for (ii=0;ii<tab->nrows;ii++) {
        if ( 0 == tab->mask[ii] )
            continue;
        
        xx[nout] = tab->xx[ii];
        yy[nout] = tab->yy[ii];
        vals[nout] = tab->vals[ii];
        mask[nout] = tab->mask[ii];
        nout++;
                
    } // end for ii


    // printf("IN: %ld\tOUT:%ld\n", tab->nrows, nout );

    tab->nrows = nout;
    free( tab->xx );
    tab->xx = xx;    
    free(tab->yy);
    tab->yy = yy;
    free(tab->mask);
    tab->mask = mask;
    free(tab->vals);
    tab->vals = vals;

    return;
}





/* If zcolumn is an interger, and has a NUL value set then setup mask */
void get_table_mask( dmDescriptor *zcolDesc, Table *tab)
{
 
  dmDataType dt = dmGetDataType( zcolDesc );
  long null_value;
  long ii;
  
  switch (dt)
  {
    case dmBYTE: // fall thru intended
    case dmSHORT:
    case dmUSHORT:
    case dmLONG:
    case dmULONG:
        if ( 0 == dmDescriptorGetNull_l( zcolDesc, &null_value ) ) {
            for (ii=tab->nrows; ii--; ) {
                tab->mask[ii] = (tab->vals[ii] == null_value) ? 0 : 1;
            }
        } else {
            for (ii=tab->nrows; ii--; ) {
                 tab->mask[ii] = 1;
              }
        }
        break;
    case dmFLOAT:
    case dmDOUBLE:
        for (ii=tab->nrows; ii--; ) {
             tab->mask[ii] = 1;
          }
        break;
    default:
        err_msg("Unsupported z-column datatype");
        break;
  } // end switch
    
  return;    
}


/* Create the output file */
OutputFile *make_output_file( dmBlock *inBlock, char *outfile, double tol, short clobber ) 
{

  OutputFile *outtab = (OutputFile*)calloc(1,sizeof(OutputFile));

  if ( 0 != ds_clobber( outfile, clobber, NULL )) {
    return(NULL);
  }

  strcat( outfile, "[REGION]");  // Always want block named 'region'

  outtab->outBlock = dmTableCreate( outfile );
  if ( NULL == outtab->outBlock ) {
    err_msg("ERROR: Cannot create output file '%s'\n", outfile );
    return(NULL);
  }

  Header_Type *hdr;
  hdr = getHdr( (void*) inBlock, hdrDM_FILE );

  if ( 1 != dmBlockGetNo( outtab->outBlock ) ) { // ie FITS file 
    dmDataset *ds = dmBlockGetDataset( outtab->outBlock );
    dmBlock *pBlock = dmDatasetMoveToBlock( ds, 1 );
    putHdr( pBlock, hdrDM_FILE, hdr, PRIMARY_STS, "dmellipse") ;
  }

  char *colname[2] = { "x", "y"};

  outtab->skycol = dmColumnCreateVector( outtab->outBlock, "pos", dmDOUBLE, 0, "pixel", "sky coordinates",
                                 colname, 2 );
  outtab->scol = dmColumnCreate( outtab->outBlock, "shape", dmTEXT, 30, "", "Region Shape");
  outtab->rcol = dmColumnCreateArray( outtab->outBlock, "r", dmDOUBLE, 0, "pixels", "Radii",
                              2 );
  outtab->acol = dmColumnCreate( outtab->outBlock, "rotang", dmDOUBLE, 0, "deg", "Rotation angle");

  /* TODO: Uncomment out below to add RFE new column.  Only commented out
   * to make regresson tests PASS */

  outtab->fcol = dmColumnCreate( outtab->outBlock, "fraction", dmDOUBLE, 0, "", "Fraction");
  outtab->dfcol = dmColumnCreate( outtab->outBlock, "input_fraction", dmDOUBLE, 0, "", "Input Fraction");

  outtab->ccol = dmColumnCreate( outtab->outBlock, "component", dmSHORT, 0, "", "Component");
  
  outtab->Xcol = dmColumnCreate( outtab->outBlock, "state", dmSHORT, 0, "", "Bit encoded state");
  outtab->Scol = dmColumnCreate( outtab->outBlock, "converge", dmTEXT, 100, "", "Convergence criteria");

  /* HDU STUFF */

  putHdr( outtab->outBlock, hdrDM_FILE, hdr, BASIC_STS, "dmellipse" );
    
  dmKeyWrite_c( outtab->outBlock, "ORIGIN",  "ASC", "", "" );
  dmKeyWrite_c( outtab->outBlock, "CONTENT", "REGION", "", "");
  dmKeyWrite_c( outtab->outBlock, "HDUDOC", "ASC-FITS-REGION-1.2: Rots, McDowell: FITS REGION "
                "Binary Table Design", "", "" );
  dmKeyWrite_c( outtab->outBlock, "HDUVERS", "1.2.0", "", "");
  dmKeyWrite_c( outtab->outBlock, "HDUCLASS", "ASC", "", "" );
  dmKeyWrite_c( outtab->outBlock, "HDUCLAS1", "REGION", "", "" );
  dmKeyWrite_c( outtab->outBlock, "HDUCLAS2", "STANDARD", "", "" );
  dmKeyWrite_d( outtab->outBlock, "FRACTOL", tol, "", "Requested tolerance on the input_fraction" );

  put_param_hist_info( outtab->outBlock, "dmellipse", NULL, 0 );

  return(outtab);
}


/* Convert the state code into a human readable string */
void get_state_string( short state, char *comment )
{

    switch (state) {
      
    case Success:
      sprintf( comment, "OK" );
      break;
    case Converge:
      sprintf( comment, "CANNOT_CONVERGE" );
      break;
    case TooSmall:
      sprintf( comment, "ELLIPSE_TOO_SMALL" );
      break;
    case TooBig:
      sprintf( comment, "ELLIPSE_TOO_BIG" );
      break;
    case TooManySteps:
      sprintf( comment, "WALK_TOO_MANY_STEPS" );
      break;
    case Working:
    default:
      sprintf( comment, "ERROR: You shouldn't be here"); /* Should never get here */
      break;
    }

}



OutputValues *find_ellipse( ProcessingParams *proc, short is_image, Image *img, Table *tab, short *mask, float frac_in  )
{
  /*

  The algorithm is something like this:

   - Get the centroid, ellipticity, and size of moments for entire image
   - Use these as a starting point.  Draw an ellipse using
     those parameters and get the enclosed fraction.
   - Now start iteration.  If the fraction is smaller than
     desired, then increase the radii (keeping the same ellipticity)
     if fraction is greater then decrease radii.
   - If the previous iteration was increasing and the current is
     decreasing (or visa versa) then change the amount we increment
     the radii by half.  So that each time we 'change direction' of
     the search we make smaller and smaller steps until we
     converge.

   - During each iteration the previously determined centroid,
     angle, and ellipticity are used with the radii.

   - Once we find a possible solution, we repeat test to 
     see if the centroid, and moments of the solution 
     still give us the desired result.  If not then we
     don't change the step size we just change the centroid,
     etc. and let the algorithm continue.  We call this
     letting the parameters 'walk' (kind of a random walk
     thing)

   - There are two other 'exit' conditions:
     . if the step size gets too small then we stop.
     . if we 'walk' too many steps then we could be in some
       kind of infinite loop so we break out -- basically
       solution is not uniuqe.

  */
    OutputValues *outvals = (OutputValues*)calloc(1,sizeof(OutputValues));

    typedef enum { DECREASE=-1, INCREASE=1} Direction;
    Direction dir = INCREASE;

    Moments mom;
    float step = proc->step_start;
    State_Code state = Working;

    float radius;  /* current major axis radius */
    double total;  /* sum of pixels for normalization */
    double min_tol, max_tol; /* convergence limits */
    short walk = 0;  /* random walk counter */
    int all_nn;     /* total number of pixels/rows */
    double frac;    /* current fraction */
    
    min_tol = frac_in - proc->tol;
    max_tol = frac_in + proc->tol;
    
    memset(&mom, 0, sizeof(Moments));
 
    if ( is_image ) {
         img->mask = (short*)calloc(img->lAxes[0]*img->lAxes[1],sizeof(short) );
         reset_mask_img( img, mask );
         all_nn = get_moments_img( img, &mom);
    } else {
          tab->mask = (short*)calloc( tab->nrows, sizeof(short));
         reset_mask_table( tab, mask );
         all_nn = get_moments_table( tab, &mom);
    }

    if ( 0 == all_nn) {
      err_msg("ERROR: No data in image\n");
    }

    /* If user supplied default/starting values, set them now */
    if ( INDEFD != proc->def_x_mu ) mom.x_mu = proc->def_x_mu;
    if ( INDEFD != proc->def_y_mu ) mom.y_mu = proc->def_y_mu;
    if ( INDEFD != proc->def_ell ) mom.ell = proc->def_ell;
    if ( INDEFD != proc->def_phi ) mom.phi = proc->def_phi;

    /* Setup normalization */
    if ( proc->normalize )
      total = mom.out[0][0];  // [0][0] is sum of pixel values 
    else
      total = 1;

    /* First guess uses a radius = major axis from entire dataset */
    radius = mom.xsig;


    while ( Working == state ) {    

      Moments dad;
      int nn;
      memset(&dad,0,sizeof(Moments));

      if ( is_image ) {
          reset_mask_img( img, mask /*orignal mask*/);
          make_ellipse_mask_img(  img,  &mom,  radius, proc->shape );
          nn = get_moments_img( img, &dad );
      } else {
          reset_mask_table( tab, mask /*orignal mask*/);
          make_ellipse_mask_table(  tab,  &mom,  radius, proc->shape );
          nn = get_moments_table( tab, &dad );
      }

      /* The (0,0) moment is the sum() of the pixel values */
      frac = dad.out[0][0] / total;


      /* Check the number of pixels/values first */

      if ( 0 == nn ) {   /* There are no pixels in the centroid/etc */
        dir = INCREASE;
        step /= 2.0;
        if ( step < proc->step_tol ) {
          err_msg("WARNING: Step size limit reached\n");
          state=Converge;
        }
        else 
          radius += step;
        goto verb;
      }

      if ( ( 1 == nn ) &&  /* There is only 1 pixel left */
           ( frac > min_tol ) ) { /* But fraction is still too high */
        err_msg("WARNING: Could not make ellipse small enough\n");
        state=TooSmall;
        goto verb;
      }

      if (( nn == all_nn) && ( frac < min_tol ) ) {
        err_msg("WARNING: Could not make ellipse big enough\n");
        state=TooBig;
        goto verb;
      }

      /* --Now check the fraction -------------------------------*/

      if ( frac < min_tol ) {
        walk = 0;
        if (DECREASE == dir ) {
          dir = INCREASE;
          step /= 2.0; /* Change dir, so go half steps */
        }
        if ( step < proc->step_tol ) {
          err_msg("WARNING: Could not converge to within tolerance=%g "
                  "of value=%g%%; final value is %g%%\n", proc->tol, 
                  (frac_in*100.0), (frac*100.0) );
          state = Converge;
        } else {
          radius += (dir*step);
        }

      } else if ( frac > max_tol ) {
        walk = 0;
        if (INCREASE == dir ) {
          dir = DECREASE;
          step /= 2.0; /* Change dir, so go half steps */
        }
        if ( step < proc->step_tol ) {
          err_msg("WARNING: Could not converge to within tolerance=%g "
                  "of value=%g%%; final value is %g%%\n", proc->tol, 
                  (frac_in*100.0), (frac*100.0) );
          state = Converge;
        } else {
          radius += (dir*step);
          if ( radius < 0 ) {  
            /*This shouldn't happen, but floats are funny so we check just to make sure */  
            radius = fabs(radius);
            step /=2.0;
          }
        } 
      } //end frac > max_tol



      /* If user asked to freeze parameters, change the new values back
       * to the frozen values */
      if ( 0 != proc->freeze_centroid ) {
        dad.x_mu = mom.x_mu;
        dad.y_mu = mom.y_mu;
      }
      if ( 0 != proc->freeze_ell ) {
        dad.ell = mom.ell;
      }
      if ( 0 != proc->freeze_angle ) {
        dad.phi = mom.phi;
      }

      /* Save the new value */
      memcpy( &mom, &dad, sizeof(Moments));

      if (( frac >= min_tol ) && ( frac <= max_tol )) {
        /* have to double check that centroid & angle are still good */

        memset(&dad,0,sizeof(Moments));
        if ( is_image ) {
            reset_mask_img( img, mask /*original mask */); 
            make_ellipse_mask_img( img, &mom, radius, proc->shape );
            get_moments_img( img, &dad );
        } else {
            reset_mask_table( tab, mask /*original mask */); 
            make_ellipse_mask_table( tab, &mom, radius, proc->shape );
            get_moments_table( tab, &dad );
        }

        frac = dad.out[0][0] / total;
        
        if ( (frac >= min_tol ) && (frac <= max_tol )) {
          state = Success; /* We're good; keep radius, and mom values */
          /* 
             At this point we could use the new dad values; but
             we verified that 'mom' converged.  If we were to use
             'dad' then we'd have to check 'dad' ... which
             would cause us to loop over and over again
          */
        } else {
          walk++;
          if ( walk >= proc->walk_max ) {
            err_msg("WARNING: Did not find a unique solution after %d"
                    " iterations for %g%%\n", walk,
                    (frac_in*100.0) );
            state = TooManySteps;
          }
          /* Keep going; we haven't converged yet */
        }
      }
      

    /* Ok, I don't like "goto"'s either, but sometimes they can be useful.  
     * Could have used continue's above but then wouldn't get the verbose
     * output which can be useful when debugging (user debugging params, not
     * developer debugging code).*/
    verb:
      if ( proc->verbose > 1 ) {
        fprintf(stderr, "state=%d\t", state);
        fprintf(stderr, "centroid=(%g,%g)\t", mom.x_mu, mom.y_mu);
        fprintf(stderr, "ellipticity=%g\t", mom.ell);
        fprintf(stderr, "angle=%g\t", mom.phi);
        fprintf(stderr, "radius=%g\t", radius );
        fprintf(stderr, "frac=%g\t", frac);
        fprintf(stderr, "dir=%d\t", dir);
        fprintf(stderr, "step=%g\t", step);
        fprintf(stderr, "\n");
      }

      
    }  /* End while() */

     
    if (is_image) 
        free(img->mask);
    else
        free(tab->mask);


    // Save output values 
    outvals->xx = mom.x_mu;
    outvals->yy = mom.y_mu;
    strcpy(outvals->shape, proc->shape);
    outvals->angle = mom.phi;
    outvals->major_axis = radius;
    outvals->minor_axis = (radius/mom.ell);
    outvals->frac_in = frac_in;
    outvals->frac_out = frac;
    outvals->state = state;

    return(outvals);
}


/* Write the values to the output file */
void write_output( OutputFile *outtab, OutputValues *outvals, int ii )
{
    /* Done, write output */
    char comment[1000];
    get_state_string( outvals->state, comment );

    double tmpval[2];
    tmpval[0] = outvals->xx;
    tmpval[1] = outvals->yy;
    dmSetVector_d( outtab->skycol, tmpval, 2 );
    dmSetScalar_c( outtab->scol, outvals->shape );
    dmSetScalar_d( outtab->acol, outvals->angle );
    tmpval[0] =  outvals->major_axis;
    tmpval[1] =  outvals->minor_axis;
    dmSetArray_d( outtab->rcol, tmpval, 2 );
    dmSetScalar_f( outtab->fcol, outvals->frac_out );
    dmSetScalar_f( outtab->dfcol, outvals->frac_in );
    dmSetScalar_s( outtab->ccol, ii );
    dmSetScalar_s( outtab->Xcol, outvals->state );
    dmSetScalar_c( outtab->Scol, comment);

    dmTableNextRow( outtab->outBlock );

}



/* Main routine */
int ellipse(void)
{
  char infile[DS_SZ_PATHNAME];
  char outfile[DS_SZ_PATHNAME];
  char grid[DS_SZ_PATHNAME];
  dmBlock *inBlock;
  
  ProcessingParams proc;
  OutputFile *outtab;
  OutputValues *outvals;

  short clobber;
  short verbose;
  
  double *frac_lo;
  long nfrac;
  long ii;

  /* Load parameters */

  clgetstr( "infile", infile, DS_SZ_PATHNAME );
  clgetstr( "outfile", outfile, DS_SZ_PATHNAME );
  clgetstr( "fraction", grid, DS_SZ_PATHNAME );

  clgetstr( "shape", proc.shape, 100 );

  proc.def_x_mu = clgetd("x_centroid");
  proc.def_y_mu = clgetd("y_centroid");
  proc.def_phi = clgetd("angle");
  proc.def_ell = clgetd("ellipticity");

  proc.freeze_centroid = clgetb("fix_centroid");
  proc.freeze_angle =clgetb("fix_angle");
  proc.freeze_ell = clgetb("fix_ellipticity");

  proc.tol = clgetd("tolerance");
  proc.step_tol = clgetd("minstep");
  proc.walk_max = clgetd("maxwalk");
  proc.step_start = clgetd("step");
  proc.normalize = clgetb("normalize");

  clobber = clgetb("clobber");
  verbose = clgeti("verbose");

  proc.verbose = verbose;



  /* Load data */
  if ( NULL == (inBlock=dmBlockOpen(NULL, infile) )) {
    err_msg("ERROR: Cannot open input file '%s'\n", infile );
    return(-1);
  }

  Image *img = NULL;
  Table *tab = NULL;
  short *mask = NULL;
  short is_image;

  if ( dmBlockGetType(inBlock) == dmIMAGE ) {
    img = load_image_file( inBlock );
    mask = img->mask;
    is_image = 1;
  } else if ( dmBlockGetType(inBlock) == dmTABLE ) {
    tab = load_table_file( inBlock, "x", "y", "z"); // TODO: make col names parameters
    mask = tab->mask;
    is_image = 0;
  } else {
    err_msg("You cannot be here.  Something bad has happened");
    return(-334);
  }
  
  /* Setup output */ 
  if ( NULL == ( outtab = make_output_file( inBlock, outfile, proc.tol, clobber)) ) {
     return(-34);
  } 


  /* If INDEF, update initial step based on image size*/
  if ( INDEFD == proc.step_start ) {
    if ( is_image ) {
        proc.step_start = sqrt( 1.0*MAX(img->lAxes[0],img->lAxes[1])) ;
    } else {
        proc.step_start = 1.0; /* TODO: maybe take sqrt(range of values); OK for now*/
    }
  }


  /* Parse the stack of fractions */
  if ( 0 == (nfrac = parse_grid(grid, &frac_lo ))) {
    return(-3);
  }


  /* Start algorithm */
  for ( ii=0;ii<nfrac;ii++ ) {

    outvals = find_ellipse( &proc, is_image, img, tab, mask, frac_lo[ii]);
    write_output( outtab, outvals, ii );
 
   free(outvals);
  } /* end loop over fraction */


  dmTableClose( outtab->outBlock);

  return(0);


}

