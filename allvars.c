#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"

 GRID ***pg_ddm;                       /* this is for mother particles */
 GRID ***pg_dm;                        /* this is for daughter particles */
 POW_BIN *p_bin_ddm;
 POW_BIN *p_bin_dm;                    /* power bins for daughters */

 int ng;                               /* grid number in each dimension */
 int gadget_file_num;
 double boxsize;
 double grid_dis;
 double total_mass;
 double total_mass_ddm;
 double total_mass_dm;
 long int total_part_num;
 double total_grid_num;
 
 GADGET_HEAD g_head;
