/* unit used here:
   1 mass = 1 soloar mass /h
   1 length = 1 Mpc/h
   k_max is the half of the nyquist frequency
*/

#define  N_thread                     10   /* number of cores for OpenMP */
#define  NG                           512  /* the dimension of the grid, should be power of 2, like 128, 256, 512...*/
#define  BIN_NUMBER                   20

#define  GADGET_INPUT_PATH            "./ICs/512-512/"
#define  OUTPUT_PATH                  "./Result/512-512/"
#define  GADGET_FILE_BASE             "snapshot_005" 
/*
  decay-dark-matter:
  snapshot_000  ->  z=4
  snapshot_001  ->  z=3  
  snapshot_002  ->  z=2
  snapshot_003  ->  z=1
  snapshot_004  ->  z=0.5
  snapshot_005  ->  z=0.4
  snapshot_006  ->  z=0.3
  snapshot_007  ->  z=0.2
  snapshot_008  ->  z=0.1
  snapshot_009  ->  z=0
*/

/* ************************************ */
#define  Pi          3.1415926535897932384626433832795028842              /* pi */

#define  SUNTOKG     1.9891e30
#define  MPCTOM      3.08567758e22
#define  G0          6.67384e-11
