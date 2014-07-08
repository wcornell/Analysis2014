/*setup.h*/
#ifndef _SETUP_H
#define _SETUP_H

typedef struct{
	char *protein_name;
	char *sim_type;
	char function_name[20];
	char dcd_filename[50];
	int runstart, runcount, runnum;
	int start_frame, end_frame;
	int (*analysis_setup_funcptr)();
	int (*analysis_funcptr)();
	int (*analysis_post_funcptr)();
	int N;

	//pairdist variables
	int numpairs;
	int *RESID1;
	int *RESID2;
	int *RESID_map;
	
	//covar variable
	float **covar_vals;

	//covargpu variables
	float *dev_X, *dev_Y, *dev_Z;
	float *dev_covar_vals, *lin_covar_vals;

	//pdbio variable
	char **line_beginning;
	FILE *pdb_file;
	
} setup_data;


int read_setup_file(char *setup_filename, setup_data *setup);

#endif
