#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


int main()
{
  omp_set_num_threads(N_thread); 
    
  init_all();
  state("Initilzation done ...");
  
  load_gedget_part_to_array();
  state("Loading Gadget done ...");
  
  grid_transfer(pg_ddm, total_mass_ddm);
  grid_transfer(pg_dm,  total_mass_dm);
  state("FFT done ...");
  
  collect_power_bin_data(p_bin_ddm, pg_ddm);
  collect_power_bin_data(p_bin_dm, pg_dm);
  state("Power done ...");
  
  write_file();
  state("Output done ...");
  
  free_all();

  return 1;
}   /* end main */
