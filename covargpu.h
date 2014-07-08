/*covargpu.h*/
#ifndef COVARGPU_H
#define COVARGPU_H
#include "setup.h"

extern "C" int covargpu_setup(setup_data *setup);
extern "C" int covargpu(setup_data *setup, float ***X, float ***Y, float ***Z);
extern "C" int covargpu_post(setup_data *setup);


#endif
