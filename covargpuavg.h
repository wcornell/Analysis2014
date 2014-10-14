/*covargpuavg.h*/
#ifndef COVARGPUAVG_H
#define COVARGPUAVG_H
#include "setup.h"

#ifdef __cplusplus
extern "C"  
#endif 
	int covargpuavg_setup(setup_data *setup);

#ifdef __cplusplus
extern "C"  
#endif 
	int covargpuavg(setup_data *setup, float ***X, float ***Y, float ***Z);

#ifdef __cplusplus
extern "C"  
#endif 
	int covargpuavg_post(setup_data *setup);


#endif
