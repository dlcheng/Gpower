#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"

void write_file()
{
  int i;
  char file_name[500];
  FILE *fp;

#ifndef NOFILE_OUTPUT  
  sprintf(file_name, "%s%s%.3f%s", OUTPUT_PATH, "g_power_z_", g_head.Redshift, ".txt");
  if((fp = fopen(file_name, "w+")) == NULL)    
    {
	  printf("Fail to open output file %s\n", file_name);
	  exit(0);
	  }
#endif
    
  for(i = 0; i < BIN_NUMBER; i++)
    {	
      double temp_k = pow(10, p_bin_ddm[i].log_k);
      double temp_pk_ddm = p_bin_ddm[i].p * pow(boxsize, 3) / pow(total_grid_num, 2);
      double temp_pk_dm = p_bin_dm[i].p   * pow(boxsize, 3) / pow(total_grid_num, 2);
      double fraction = total_mass_dm / (total_mass_dm + total_mass_ddm);
      int num = p_bin_ddm[i].n;

      printf("%.6e\t%.6e\t%.6e\t%.6e\t%d\n", temp_k, temp_pk_ddm, temp_pk_dm, fraction, num);
#ifndef NOFILE_OUTPUT      
      fprintf(fp, "%.6e\t%.6e\t%.6e\t%.6e\t%d\n", temp_k, temp_pk_ddm, temp_pk_dm, fraction, num);	
#endif      
    }
#ifndef NOFILE_OUTPUT
  fclose(fp);	
#endif
}    /* end write file */	
