/*covar.c*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "setup.h"
#include "covar.h"
#include "newdcdio.h"



int covar_setup(setup_data *setup){
	setup->covar_vals = (float **) malloc(setup->N*sizeof(float *));
	for(int i = 0; i < setup->N; i++){
		setup->covar_vals[i] = (float*) malloc(setup->N*sizeof(float));
		for(int j = 0; j < setup->N; j++){
			setup->covar_vals[i][j] = 0;
		}
	}
	//printf("Memory allocated\n");
	return 0;
}



int covar(setup_data *setup, float ***X, float ***Y, float ***Z){

	float delX_i, delY_i, delZ_i;
	float delX_j, delY_j, delZ_j;

	float numerator, denominator;
	float percent_done;

	for(int frame = setup->start_frame; frame < setup->end_frame; frame++){
		for(int i = 0; i < setup->N; i++){
			for(int j = 0; j < setup->N; j++){
				if(i <= j){				
					delX_i = (*X)[frame+1][i]-(*X)[0][i];
					delY_i = (*Y)[frame+1][i]-(*Y)[0][i];
					delZ_i = (*Z)[frame+1][i]-(*Z)[0][i];
					delX_j = (*X)[frame+1][j]-(*X)[0][j];
					delY_j = (*Y)[frame+1][j]-(*Y)[0][j];
					delZ_j = (*Z)[frame+1][j]-(*Z)[0][j];

					if((delX_i == 0 && delY_i == 0 && delZ_i == 0)||(delX_j == 0 && delY_j == 0 && delZ_j == 0)){
						setup->covar_vals[i][j] = 0;
					}
					else{
						numerator = delX_i*delX_j + delY_i*delY_j + delZ_i*delZ_j;
						denominator = sqrt(delX_i*delX_i + delY_i*delY_i + delZ_i*delZ_i)*sqrt(delX_j*delX_j + delY_j*delY_j + delZ_j*delZ_j);
						//printf("%f\n", numerator/denominator);
						setup->covar_vals[i][j] += numerator/denominator;
						if(j!=i) setup->covar_vals[j][i] += numerator/denominator;
					}
				}
				//printf("%7.3f/%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f\n", numerator, denominator, delX_i, delY_i, delZ_i, delX_j, delY_j, delZ_j);
			}
		}
		percent_done = ((float) frame / (setup->end_frame-setup->start_frame)) * 100;
		printf("\r%02.0f%%  ", percent_done);
		fflush(stdout);
	}
	printf("\b\b\b\b");
	fflush(stdout);



	return 0;
}


/*int covar_post(setup_data *setup){
	char dat_filename[50];
	FILE *dat_file;

	sprintf(dat_filename, "%s/%s/dat/covar_%s_%d-%d_%d-%d.dat", setup->protein_name, setup->sim_type, setup->protein_name, setup->start_frame, setup->end_frame, setup->runstart, setup->runstart+setup->runnum);
	dat_file = fopen(dat_filename, "w");
	if(dat_file == NULL){
		printf("Failed to open file: %s", dat_filename);
		exit(0);
	}
	float iterations = (setup->runcount*(setup->end_frame-setup->start_frame));
	for(int i = 0; i < setup->N; i++){
		for(int j = 0; j < setup->N; j++){
			setup->covar_vals[i][j] /= iterations;
			fprintf(dat_file, "%d %d %f\n", i, j, setup->covar_vals[i][j]);
		}
	}
	
	fclose(dat_file);



	for(int i = 0; i < setup->N; i++){
		free(setup->covar_vals[i]);
	}
	free(setup->covar_vals);
	//printf("Memory freed\n");
	return 0;
}*/

int covar_post(setup_data *setup){
	char dat_filename[50];
	FILE *dat_file;

	sprintf(dat_filename, "%s/%s/dat/covar_%s_%d-%d_%d-%d.dat", setup->protein_name, setup->sim_type, setup->protein_name, setup->start_frame, setup->end_frame, setup->runstart, setup->runstart+setup->runnum);
	dat_file = fopen(dat_filename, "w");
	if(dat_file == NULL){
		printf("Failed to open file: %s", dat_filename);
		exit(0);
	}
	float iterations = (setup->runcount*(setup->end_frame-setup->start_frame));
	for(int i = 0; i < setup->N; i++){
		for(int j = 0; j < setup->N; j++){
			setup->covar_vals[i][j] /= iterations;
			fprintf(dat_file, "%10.6f", setup->covar_vals[i][j]);
		}
	fprintf(dat_file, "\n");
	}
	
	fclose(dat_file);



	for(int i = 0; i < setup->N; i++){
		free(setup->covar_vals[i]);
	}
	free(setup->covar_vals);
	//printf("Memory freed\n");
	return 0;
}
