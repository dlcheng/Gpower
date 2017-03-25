#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_complex.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


void grid_transfer(GRID ***pg, double local_total_mass)
{  
  
  flip_data(pg, local_total_mass);
  fft_3d(1, pg);
	
}     /* end grid_transfer */	

void flip_data(GRID ***pg, double local_total_mass)
{
  int i, j, k;
  int m;
  
/* OpenMP */
#pragma omp parallel shared(ng, pg, total_grid_num, local_total_mass) private(i, j, k, m)
{  
#pragma omp for  
  for(i = 0; i < ng; i++)
    for(j = 0; j < ng; j++)
      for(k = 0; k < ng; k++)
        {
		    m = i + j + k;	

		    pg[i][j][k].Re = pg[i][j][k].Re * total_grid_num / local_total_mass - 1;  /* overdensity */
		    pg[i][j][k].Im = pg[i][j][k].Im * total_grid_num / local_total_mass - 1;		 
         
		    if((m % 2) == 1)
		      {
		      pg[i][j][k].Re = -1 * pg[i][j][k].Re;
		      pg[i][j][k].Im = -1 * pg[i][j][k].Im;		
	        }	
	      }	
}      
}	        /* end flip_data */

void fft_3d(int flag, GRID ***f)
{                             /* flag = 0 means k->x, 
                                 flag = 1 means x->k 
                                 f is the fft cube
                                 only flag = 1 is allowed here
                              */


 int i,j,k;
 double *fft_data_local;
 
 #pragma omp parallel shared(ng, f, flag) private(i, j, k, fft_data_local)
{
 fft_data_local = (double *)malloc(2 * ng * sizeof(double));

 #pragma omp barrier
 #pragma omp for
 for(j=0; j<ng; j++)      
  {
   for(k=0; k<ng; k++)
     {
     copy_to_fft_array(-1,j,k,f,fft_data_local);

     if(flag == 0)
        gsl_fft_complex_radix2_inverse(fft_data_local,1,ng);
     if(flag == 1)	   
        gsl_fft_complex_radix2_forward(fft_data_local,1,ng);      

     copy_from_fft_array(-1,j,k,f,fft_data_local);     
     }
  }
  
  #pragma omp barrier
  #pragma omp for
  for(k=0; k<ng; k++)
  {
   for(i=0; i<ng; i++)
     {
     copy_to_fft_array(i,-1,k,f,fft_data_local);

     if(flag == 0)
        gsl_fft_complex_radix2_inverse(fft_data_local,1,ng);
     if(flag == 1)
        gsl_fft_complex_radix2_forward(fft_data_local,1,ng);

     copy_from_fft_array(i,-1,k,f,fft_data_local);     
     }
  }

 #pragma omp barrier
 #pragma omp for
 for(i=0; i<ng; i++)
   {
   for(j=0; j<ng; j++)
     {
     copy_to_fft_array(i,j,-1,f,fft_data_local);

     if(flag == 0)
        gsl_fft_complex_radix2_inverse(fft_data_local,1,ng);
     if(flag == 1)
        gsl_fft_complex_radix2_forward(fft_data_local,1,ng);

     copy_from_fft_array(i,j,-1,f,fft_data_local);     
     }
    }
 free(fft_data_local);    
}     /* end OpenMP */
}                             /* end three_d_fft*/

void copy_to_fft_array(int i, int j, int k, GRID ***f, double *fft_data_local)
{
 int m;
 if(i == -1)
   {
    for(m=0; m<ng; m++)
      {	  
      fft_data_local[2*m] = f[m][j][k].Re;
      fft_data_local[2*m+1] = f[m][j][k].Im;
      }
   }

 if(j == -1)
   {
    for(m=0; m<ng; m++)
      {
      fft_data_local[2*m] = f[i][m][k].Re;
      fft_data_local[2*m+1] = f[i][m][k].Im;
      }
   }

 if(k == -1)
   {
    for(m=0; m<ng; m++)
      {
      fft_data_local[2*m] = f[i][j][m].Re;
      fft_data_local[2*m+1] = f[i][j][m].Im;
      }
   }
}                             /* end copy_to_fft_array */

void copy_from_fft_array(int i, int j, int k, GRID *** f, double *fft_data_local)
{
 int m;
 if(i == -1)
   {
    for(m=0; m<ng; m++)
      {
      f[m][j][k].Re = fft_data_local[2*m];
      f[m][j][k].Im = fft_data_local[2*m+1];
      }
   }

 if(j == -1)
   {
    for(m=0; m<ng; m++)
      {
      f[i][m][k].Re = fft_data_local[2*m];
      f[i][m][k].Im = fft_data_local[2*m+1];
      }
   }

 if(k == -1)
   {
    for(m=0; m<ng; m++)
      {
      f[i][j][m].Re = fft_data_local[2*m];
      f[i][j][m].Im = fft_data_local[2*m+1];
      }
   }
}                             /* end copy_from_fft_array */
