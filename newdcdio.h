/*newdcdio.h*/
#ifndef NEWDCDIO_H
#define NEWDCDIO_H
#include "setup.h"




void dcd_read_header(FILE *dcd_file, char *dcd_filename, int* N, int* NFILE, int* NPRIV, int* NSAVC, float* DELTA);
int dcd_read_header_N(FILE *header_file, setup_data *setup);
int dcd_read_frame(FILE *dcd_file, int N, float *X, float *Y, float *Z);
int dcd_read_file(FILE *dcd_file, int N, float ***X, float ***Y, float ***Z, int start_frame, int end_frame);
int dcd_single_analysis(setup_data *setup);

#endif
