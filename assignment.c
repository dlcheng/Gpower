#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


void assign_part(float *pos, float mass, unsigned int id)
{
#ifdef NGP
  int i, j, k;

  i = (int) (pos[0] / grid_dis);
  if(i > ng-1)
    i = 0;

  j = (int) (pos[1] / grid_dis);
  if(j > ng-1)
    j = 0;

  k = (int) (pos[2] / grid_dis);
  if(k > ng-1)
    k = 0;
  
  if(id < g_head.Original_num)
    {
     pg_ddm[i][j][k].Re += mass;    /* this is mother particle */
     total_mass_ddm += mass;
    }
  else
    {
     pg_dm[i][j][k].Re += mass;  
     total_mass_dm += mass;
    }
  
#endif
}    /* end assign_part */
