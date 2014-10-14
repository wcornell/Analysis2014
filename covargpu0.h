/*covargpu0.h*/
#ifndef COVARGPU0_H
#define COVARGPU0_H
#include "setup.h"

#ifdef __cplusplus
extern "C"  
#endif 
	int covargpu0_setup(setup_data *setup);

#ifdef __cplusplus
extern "C"  
#endif 
	int covargpu0(setup_data *setup, float ***X, float ***Y, float ***Z);

#ifdef __cplusplus
extern "C"  
#endif 
	int covargpu0_post(setup_data *setup);


#endif
