#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


void collect_power_bin_data(POW_BIN *p_bin, GRID ***pg)
{
  int i, j, k;
  int i_n, j_n, k_n;
  double kx, ky, kz, wx, wy, wz, wk;
  int m;
  double k_min = 2 * Pi / boxsize;
  double k_max = Pi / boxsize * ng * sqrt(3) / 2.0;		/* half the length of the three-dimension Hyquist frequency */
  double log_k_dis = (log10(k_max) - log10(k_min)) / (double) BIN_NUMBER;
  double k_now;
  
  init_power_bin(p_bin, k_min, log_k_dis);
  
  for(i = 0 ; i < ng; i++)
    for(j = 0; j < ng; j++)
      for(k = 0; k < ng; k++)
        {
		 i_n = i - ng/2;
		 j_n = j - ng/2;
		 k_n = k - ng/2;
		 
#ifdef NGP		 
		 kx = 2 * Pi / boxsize * i_n;
		 ky = 2 * Pi / boxsize * j_n;
		 kz = 2 * Pi / boxsize * k_n;	
		 
		 if(i_n == 0)
		   wx = 1;
		 else 
		   wx = sin(0.5 * grid_dis * kx) / (0.5 * grid_dis * kx);

		 if(j_n == 0)
		   wy = 1;
		 else 
		   wy = sin(0.5 * grid_dis * ky) / (0.5 * grid_dis * ky);

		 if(k_n == 0)
		   wz = 1;
		 else 
		   wz = sin(0.5 * grid_dis * kz) / (0.5 * grid_dis * kz);
		   
		  wk = wx * wy * wz; 		   		     
#endif
		 		
	     k_now = 2 * Pi / boxsize * sqrt(i_n * i_n + j_n * j_n + k_n * k_n);	     
           	     
	     if(k_now >= k_min && k_now <= k_max)
	       {
			m = (int) (log10(k_now/k_min) / log_k_dis);
			if( m > BIN_NUMBER - 1)
			  m = BIN_NUMBER - 1;
			   
		    p_bin[m].p = p_bin[m].p + (pow(pg[i][j][k].Re, 2) + pow(pg[i][j][k].Im , 2)) / wk / wk;
		    p_bin[m].n++;	   
           }			   
	    }
	    
   for(i = 0 ; i < BIN_NUMBER; i++)	    
	 p_bin[i].p = p_bin[i].p / (double) p_bin[i].n;
	
}    /* end collect_power_bin_data */	


void init_power_bin(POW_BIN * p_bin, double k_min, double log_k_dis)
{
  int i;
	
  for(i = 0 ; i < BIN_NUMBER; i++)
	{
	p_bin[i].log_k = log10(k_min) + (i + 0.5) * log_k_dis;
	p_bin[i].p = 0;
	p_bin[i].n = 0;
    }
}   /* end init_power_bin */	
