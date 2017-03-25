#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"
#include "define.h"


void detect_and_link_gadget_file()
{
  char file_name[500];
  unsigned int block_size;
  
  FILE *fp;
  
  sprintf(file_name, "%s%s", GADGET_INPUT_PATH, GADGET_FILE_BASE);	
	
  if((fp = fopen(file_name, "rb")) != NULL)	
    {
	 state("Found SINGLE gadget file.");	 		   		  
	 rewind(fp);
	 if(fread(&block_size, sizeof(int), 1, fp));
	 if(fread(&g_head, sizeof(GADGET_HEAD), 1, fp));
	 gadget_file_num = 1;	 
	 fclose(fp);
    }
  else
    {
	 sprintf(file_name, "%s%s%s%d", GADGET_INPUT_PATH, GADGET_FILE_BASE, ".", 0);
	 if((fp = fopen(file_name, "rb")) != NULL)	
           {
	   state("Found MULTI gadget files.");		   		  
	   rewind(fp);
	   if(fread(&block_size, sizeof(int), 1, fp));
	   if(fread(&g_head, sizeof(GADGET_HEAD), 1, fp));
	   gadget_file_num = g_head.NumFiles;	 
	   fclose(fp);
           }
         else
           warn_and_end("OPEN GADGET FILE");  	 							
    }		
	
/* Total number of dark matter particles used in the simulation */ 
  total_part_num = g_head.NallHW[1];
  total_part_num = total_part_num << 32;
  total_part_num += g_head.Nall[1];
  
  boxsize = g_head.BoxSize / 1e3;                  /* boxsize in unit of Mpc/h */
  
}                            /* end load_gadget_head */


/* This function fills the p_part array */
void load_gedget_part_to_array()
{
  int i;
  unsigned block_size, j;
  GADGET_HEAD local_head;
  char file_name[500];
  FILE *fp_mass;
  FILE *fp_id;
  FILE *fp_pos;
  long int disp;
  float pos[3];
  float mass;
  unsigned int id;
  
  for(i=0; i<gadget_file_num; i++)
     {		 
		  if(gadget_file_num == 1)	
	      {
        sprintf(file_name, "%s%s", GADGET_INPUT_PATH, GADGET_FILE_BASE);		   
	      fp_pos = fopen(file_name, "rb");
        fp_id  = fopen(file_name, "rb");
	      fp_mass = fopen(file_name, "rb");
	      rewind(fp_pos);
        rewind(fp_id);
	      rewind(fp_mass);	   	   
	   
	      check_gadget_file(file_name);		      		
  	    }
  	  else
  	    {
	      sprintf(file_name, "%s%s%s%d", GADGET_INPUT_PATH, GADGET_FILE_BASE, ".", i);		   		   
	      fp_pos = fopen(file_name, "rb");
        fp_id  = fopen(file_name, "rb");
	      fp_mass = fopen(file_name, "rb");
	      rewind(fp_pos);
        rewind(fp_id);
	      rewind(fp_mass);	   	   
	   
	      check_gadget_file(file_name);	   
	      }  		 
		 
      if(fread(&block_size, sizeof(block_size),  1, fp_pos));
      if(fread(&local_head, sizeof(GADGET_HEAD), 1, fp_pos));
      if(fread(&block_size, sizeof(block_size),  1, fp_pos));
      if(fread(&block_size, sizeof(block_size),  1, fp_pos));
      
      disp = sizeof(int) * 9 + 256 + sizeof(float) * 3 * local_head.Npart[1] * 2 + sizeof(int) * local_head.Npart[1];
      if(local_head.Massarr[1] == 0)
        fseek(fp_mass, disp, 1);                   /* now pointing to the first mass */

      disp = sizeof(int) * 7 + 256 + sizeof(float) * 3 * local_head.Npart[1] * 2;
      fseek(fp_id, disp, 1);
        
/* Read the POS and mass */     
      for(j=0; j<local_head.Npart[1]; j++)
        {
         if(fread(pos, 3 * sizeof(float), 1, fp_pos))
           {        
           pos[0] = pos[0] / 1e3;
           pos[1] = pos[1] / 1e3;
           pos[2] = pos[2] / 1e3;                /* all in unit of Mpc/h */
	         }
	       
         if(local_head.Massarr[1] == 0)
           {
			     if(fread(&mass, sizeof(float), 1, fp_mass));  
			        mass = mass * 1e10; 			   
		       }
		     else
		       mass = local_head.Massarr[1] * 1e10;          /* in unit of 1 solar Mass /h */

         if(fread(&id, sizeof(unsigned int), 1, fp_id));  /* load the ID of particle */
                
         wrap_pos(pos);
         assign_part(pos, mass, id);		   
         
//         printf("Pos:[%.3e, %.3e, %.3e], mass: %.3e\n", pos[0], pos[1], pos[2], mass);        
	    }
	    
      fclose(fp_pos);
      fclose(fp_mass);	             
      fclose(fp_id);
    }
}   /* end load_gadget_part */

/* translate particles to [0, boxsize) */
void wrap_pos(float *pos)
{
  int i;	
	
  for(i=0; i<3; i++)
   {
	while(pos[i] < 0)
	  {
	  pos[i] += boxsize;  
	  }	  
	while(pos[i] >= boxsize)
	  {
	  pos[i] -= boxsize;
	  }	    
   }
}    /* end wrap_pos */	

/* This function is used to check the consistence of the gadget file */

void check_gadget_file(char *file_name)
{ 
  FILE *fp = fopen(file_name, "rb");
  GADGET_HEAD H;
  unsigned temp1, temp2;  /* used to read the block 4-byte-head */
  unsigned cali;          /* block value include overflow */
  long int disp;          /* displacement to move the file point */

/* HEAD block check */
#ifdef VERBOSE 
  printf("Begin to check the file %s\n", file_name); 
  printf("Checking the HEAD ...\n");
#endif  
  if(fread(&temp1, sizeof(temp1), 1, fp));
  if(fread(&H, 256, 1, fp));
  if(fread(&temp2, sizeof(temp2), 1, fp));     
  if(temp2!=temp1)
    {
    printf("Eorror detected in the checking HEAD block.\n");
    exit(0);
    }
/* POS block check */   
#ifdef VERBOSE
  printf("Checking the POS ...\n");  
#endif  
  cali = 3 * sizeof(float) * H.Npart[1];
  disp = 3 * sizeof(float) * H.Npart[1];  /* disp != clai */
  if(fread(&temp1, sizeof(temp1), 1, fp));
  if(temp1 != cali)
    {
	 printf("Error detected in the checking POS block.\n");
	 exit(0);
	}
  fseek(fp, disp, 1);
  if(fread(&temp2, sizeof(temp2), 1, fp));
  if(temp1 != temp2)
    {
	 printf("Error detected in the checking POS block.\n");
	 exit(0);
	}
/* VEL block check */ 	
#ifdef VERBOSE
  printf("Checking the VEL structure...\n");
#endif  
  cali = 3 * sizeof(float) * H.Npart[1];
  disp = 3 * sizeof(float) * H.Npart[1];  /* disp != clai */  
  if(fread(&temp1, sizeof(temp1), 1, fp));
  if(temp1 != cali)
    {
	 printf("Error detected in the checking VEL block.\n");
	 exit(0);
	}
  fseek(fp, disp, 1);
  if(fread(&temp2, sizeof(temp2), 1, fp));
  if(temp1 != temp2)
    {
	 printf("Error detected in the checking VEL block.\n");
	 exit(0);
	}	
/* ID block check */
#ifdef VERBOSE
  printf("Checking the ID structure...\n"); 
#endif  
  cali = sizeof(int) * H.Npart[1];
  disp = sizeof(float) * H.Npart[1];  /* disp != clai */  
  if(fread(&temp1, sizeof(temp1), 1, fp));
  if(temp1 != cali)
    {
	 printf("Error detected in the checking ID block.\n");
	 exit(0);
	}
  fseek(fp, disp, 1);
  if(fread(&temp2, sizeof(temp2), 1, fp));
  if(temp1 != temp2)
    {
	 printf("Error detected in the checking ID block.\n");
	 exit(0);
	}	
/* Mass block check, if any */
  if(H.Massarr[1] == 0.0)
    { 
#ifdef VERBOSE			
    printf("Checking the MASS structure...\n");		
#endif    
    cali = sizeof(float) * H.Npart[1];
    disp = sizeof(float) * H.Npart[1];  /* disp != clai */    	
    if(fread(&temp1, sizeof(temp1), 1, fp));
    if(temp1 != cali)
      {
       printf("Error detected in the checking MASS block.\n");
	   exit(0);
	  }
    fseek(fp, disp, 1);
    if(fread(&temp2, sizeof(temp2), 1, fp));
    if(temp1 != temp2)
      {
	   printf("Error detected in the checking MASS block.\n");
	   exit(0);
      }		
    }

  fclose(fp);
/* Done with the check */
#ifdef VERBOSE
  printf("Gadget file check done, the file is consistent.\n");        
#endif  
		
}                            /* end check_gadget_file */	
