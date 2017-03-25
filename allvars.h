#ifndef ALLVAR_H
#define ALLVAR_H

typedef struct grid GRID;
typedef struct gadget_head GADGET_HEAD;
typedef struct pow_bin POW_BIN;

struct gadget_head
{  
  unsigned int Npart[6];
  double Massarr[6];
  double Time;
  double Redshift;
  int FlagSfr;
  int FlagFeedback;
  int Nall[6];
  int  FlagCooling;
  int NumFiles;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int FlagAge;
  int FlagMetals;
  int NallHW[6];
  int Flag_entr_ics;
#ifndef DECAY_DARK_MATTER
  char unused[60];         
#else
  unsigned int Original_num;                 /* the original particle number before decay */
  double Original_soften_length;             /* original particle's softening length  */
  char fill[48];                             /* new unused */
#endif
}; 

struct grid
{
  float Re;                                  /* real part of the complex number */
  float Im;                                  /* imaginary part of the complex number */
};

struct pow_bin
{ 
  double log_k;                              /* log mean k of the bin */
  double p;	
  unsigned int n;                            /* number of grids in the bin */
};	


extern GRID ***pg_ddm;                       /* this is for mother particles */
extern GRID ***pg_dm;                        /* this is for daughter particles */
extern POW_BIN *p_bin_ddm;
extern POW_BIN *p_bin_dm;                    /* power bins for daughters */

extern int ng;                               /* grid number in each dimension */
extern int gadget_file_num;
extern double boxsize;
extern double grid_dis;
extern double total_mass;
extern double total_mass_ddm;
extern double total_mass_dm;
extern long int total_part_num;
extern double total_grid_num;

extern GADGET_HEAD g_head;


#endif
