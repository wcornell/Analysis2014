/*pdbio.h*/
#ifndef PDBIO_H
#define PDBIO_H
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "setup.h"

int pdb_setup(setup_data *setup);
int pdb_write(setup_data *setup, float ***X, float ***Y, float ***Z);
int pdb_post(setup_data *setup);

#endif
