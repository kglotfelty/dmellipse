/*                                                                
**  Copyright (C) 2015  Smithsonian Astrophysical Observatory 
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


#include <ascdm.h>

/* algorithm (not tool) exit codes.  */
typedef enum {
  Working = -1,  /* Algorithm is iterating */
  Success,        /* Finshed successfully, fraction within tolerance */
  Converge,       /* Could not converge to value within the tolerance specified */
  TooSmall,       /* Cannot make ellipse small enough (1 pixel remains) */
  TooBig,         /* Cannot make ellipse big enough, fills entire image */
  TooManySteps    /* Cannot find unique solution within specified 'walk' size*/
} State_Code;


/* 2D moments */
typedef struct {
  double x_mu;    // x-centroid
  double y_mu;     // y-centroid
  double out[3][3];  // matrix of 2nd order moments
  double phi;        // rotation angle (CCW +x to major axis)
  double ecc;        // eccenetricity
  double xsig;       // major axis
  double ysig;       // minor axis
  double ell;        // ellipticity
} Moments;


/* Hold info for an input image */
typedef struct {
  void *data;        // pixel values
  dmDataType dt;     // pixel datatype
  long *lAxes;       // axis lenghts
  short *mask;        // mask of valid pixels
  dmDescriptor *xdesc;  // X (or sky) coordinate descriptor
  dmDescriptor *ydesc;  // Y coordinate descriptor
} Image;


/* Hold info for a table */
typedef struct {
  double *xx;       // X values
  double *yy;       // Y values
  double *vals;     // Row values (optional, may set = 1)
  short *mask;      // Valid rows  /* TODO: deal w/ interger null values */
  long nrows;       // Number of rows (len of earlier arrays)
} Table;


/* Hold the parameters that control the algorithm processing */
typedef struct {

  double tol;           // Tolerance on the input fraction (absolute value)
  double step_tol;      // Minimum step size
  double step_start;    // Initial step isze
  long walk_max;        // Maximum number of equal-fraction "walks" that are allowed
  short freeze_centroid;    // Keep centroid at initial value?
  short freeze_ell;         // Keep ellipticity at initial value?
  short freeze_angle;       // Keep angle at initial value?
  char shape[100];          // Shape of region (ellipse or rotbox)

  double def_x_mu;          // initial x-centroid, INDEF == compute
  double def_y_mu;          // initial y-centroid, INDEF == compute
  double def_ell;           // initial ellipticity, INDEF == compute
  double def_phi;           // initial rotation angle, INDEF == compute
  short normalize;          // Should values be renormalized to sum to 1.0?
  short verbose;            // verbose/chatter setting
} ProcessingParams;


typedef struct {
  dmBlock *outBlock;
  dmDescriptor *skycol; // sky x,y column
  dmDescriptor *scol;   // shape column
  dmDescriptor *ccol;   // component number column
  dmDescriptor *rcol;   // radii column
  dmDescriptor *acol;   // angle column
  dmDescriptor *Xcol;   // state column
  dmDescriptor *fcol;   // fraction achieved column
  dmDescriptor *Scol;   // state description
  dmDescriptor *dfcol;  // desired fraction column
} OutputFile;


typedef struct {
  double xx;          // x centroid
  double yy;            // y centroid
  char shape[10];       // shape (ellipse or rotbox)
  double angle;         // rotation angle CCW from +x
  double major_axis;    // ellipse major axis
  double minor_axis;    // ellipse minor axis
  double frac_in;       // requested input fraction
  double frac_out;      // fraction achieved
  short state;          // convergence/exit criteria
} OutputValues;


/* Convert index/image coords to physical coordinates */
int convert_coords( Image *img,  // Image data
                    long xx,      // i: x value
                    long yy,      // i: y value
                    double *xat,  // o: physical x
                    double *yat   // o: physical y
                    );

/* Get the moments for an image or a table */
int get_moments_img( Image *img, Moments *mom );
int get_moments_table( Table *tab, Moments *mom );

/* Create a mask which includes the $shape defined by the moments, offset*/
void make_ellipse_mask_img( Image *img, Moments *mom, float offset, char *shape);
void make_ellipse_mask_table( Table *tab, Moments *mom, float offset, char *shape);

/* Reset mask to original values */
void reset_mask_img(Image *img, short *mask );
void reset_mask_table(Table *tab, short *mask );

/* Parse the input grid parameter */
long parse_grid( char *grid, double **vals );

/* Main routine */
int ellipse(void);

/* Load image data */
Image *load_image_file( dmBlock *inBlock );

/* Load table data */
Table *load_table_file( dmBlock *inBlock, char *xcol, char *ycol, char *zcol);

/* Compute the major, minor, ellipticity from the moment matrix */
void compute_ellipse_vals_from_moment( Moments *mom, double sum_vals );

/* Create output file */
OutputFile *make_output_file( dmBlock *inBlock, char *outfile, short clobber );

/* Turn state code into read words */
void get_state_string( short state, char *comment );

/* Write output, one row at a time */
void write_output( OutputFile *outtab, OutputValues *outvals, int ii );

/* Main ellipse finding routine */
OutputValues *find_ellipse( ProcessingParams *proc, short is_image, Image *img, Table *tab, short *mask, float frac_in  );


void get_table_mask( dmDescriptor *zcolDesc, Table *tab);
