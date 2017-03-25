#ifndef ALLVAR_H
 #include "allvars.h"
#endif

void assign_part(float *pos, float mass, unsigned int id);

void grid_transfer(GRID ***pg, double local_total_mass);
void fft_3d(int flag, GRID ***f);
void flip_data(GRID ***pg, double local_total_mass);
void copy_to_fft_array(int i, int j, int k, GRID ***f, double *fft_data_local);
void copy_from_fft_array(int i, int j, int k, GRID *** f, double *fft_data_local);

void init_all();

void detect_and_link_gadget_file();
void load_gedget_part_to_array();
void wrap_pos(float *pos);
void check_gadget_file(char *file_name);

GRID *** alloc_3d_array();
void free_3d_array(GRID *** array);
POW_BIN * alloc_power_bin_array();
void free_all();

void write_file();

void collect_power_bin_data(POW_BIN *p_bin, GRID ***pg);
void init_power_bin(POW_BIN *p_bin, double k_min, double log_k_dis);

void warn_and_end(char *s);
void state(char *s);
