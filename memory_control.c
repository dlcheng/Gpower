#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


GRID *** alloc_3d_array()
{
 int i,j;

  GRID *** pointer = (GRID ***) malloc(ng * sizeof(GRID **));

  for(i=0; i<ng ; i++)
    pointer[i] = (GRID **) malloc(ng * sizeof(GRID *));
   
  for(i=0; i<ng ; i++)
    for(j=0; j<ng ; j++)
       pointer[i][j] = (GRID *) malloc(ng * sizeof(GRID));

  return pointer;
}                              /* end create_3d_array */

void free_3d_array(GRID *** array)
{
 int i,j;
  
 for(i=0; i<ng; i++)
   for(j=0; j<ng ; j++)
      free(array[i][j]);

 for(i=0; i<ng; i++)
   free(array[i]);

 free(array);
}                               /* end free_3d_array */
  
POW_BIN * alloc_power_bin_array()
{
  POW_BIN *p_bin;

  if(( p_bin = (POW_BIN *) malloc(BIN_NUMBER * sizeof(POW_BIN)) ) == NULL)
    warn_and_end("Fail to allocate memory for the POWER_BIN data"); 		

  return p_bin;
}	


void free_all()
{  
  free_3d_array(pg_ddm);
  free_3d_array(pg_dm);
  
  free(p_bin_ddm);
  free(p_bin_dm);
		
}   /* end free_all */	
