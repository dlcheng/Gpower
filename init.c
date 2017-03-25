#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

void init_all()
{
  int i, j, k;
  ng = NG;
  
  pg_ddm = alloc_3d_array(); 
  pg_dm  = alloc_3d_array();

  p_bin_ddm = alloc_power_bin_array();
  p_bin_dm  = alloc_power_bin_array();
 
  for(i = 0; i < ng; i++)
    for(j = 0; j < ng; j++)
      for(k = 0; k < ng; k++)
        {
        pg_ddm[i][j][k].Re = 0;
        pg_ddm[i][j][k].Im = 0;
        pg_dm[i][j][k].Re  = 0;
        pg_dm[i][j][k].Im  = 0;
	      }
	    
  detect_and_link_gadget_file();
  
  grid_dis = boxsize / (double) ng;  
  total_grid_num = ng * ng * ng;
  total_mass =  (3e10 * MPCTOM /(8.0 * Pi * G0) / SUNTOKG) * pow(boxsize, 3) * g_head.Omega0;
  total_mass_ddm = 0.0;
  total_mass_dm  = 0.0;
  
  printf("Omega_m = %.2f, Boxsize = %.1f[Mpc/h]\n", g_head.Omega0, g_head.BoxSize/1e3);
  printf("total_mass = %.3e, Total_num = %ld\n", total_mass, total_part_num);
  printf("Original_num of particles = %ld\n", g_head.Original_num);

}    /* end init_all */
