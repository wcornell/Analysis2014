/*pairdist.h*/
#ifndef PAIRDIST_H
#define PAIRDIST_H
#include "setup.h"

int pairdist(setup_data *setup, float ***X, float ***Y, float ***Z);
int pairdist_setup(setup_data *setup);
int pairdist_post(setup_data *setup);

#endif
