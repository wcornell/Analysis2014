/*covargpu.h*/
#ifndef COVARGPU_H
#define COVARGPU_H
#include "setup.h"

#ifdef __cplusplus
extern "C"  
#endif 
	int covargpu_setup(setup_data *setup);

#ifdef __cplusplus
extern "C"  
#endif 
	int covargpu(setup_data *setup, float ***X, float ***Y, float ***Z);

#ifdef __cplusplus
extern "C"  
#endif 
	int covargpu_post(setup_data *setup);


#endif
