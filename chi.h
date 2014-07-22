/*chi.h*/
#ifndef CHI_H
#define CHI_H
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "setup.h"

int chi_setup(setup_data *setup);
int chi(setup_data *setup, float***X, float ***Y, float ***Z);
int chi_post(setup_data *setup);

#endif
