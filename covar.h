/*covar.h*/
#ifndef COVAR_H
#define COVAR_H
#include "setup.h"

int covar_setup(setup_data *setup);
int covar(setup_data *setup, float ***X, float ***Y, float ***Z);
int covar_post(setup_data *setup);



#endif
