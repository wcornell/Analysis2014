/*covargpu.c*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "setup.h"
#include "covargpu.h"
#include "newdcdio.h"


__global__ void covargpu_add(setup_data *setup){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	float delX_i, delY_i, delZ_i;
	float delX_j, delY_j, delZ_j;

	float numerator, denominator;
	for(int i = 0; i < setup->N; i++){
		for(int j = 0; j < setup->N; j++){
			while(tid < (setup->end_frame-setup->start_frame+1)){
				delX_i = (setup->dev_X)[tid][i]-(setup->dev_X)[0][i];
				delY_i = (setup->dev_Y)[tid][i]-(setup->dev_Y)[0][i];
				delZ_i = (setup->dev_Z)[tid][i]-(setup->dev_Z)[0][i];
				delX_j = (setup->dev_X)[tid][j]-(setup->dev_X)[0][j];
				delY_j = (setup->dev_Y)[tid][j]-(setup->dev_Y)[0][j];
				delZ_j = (setup->dev_Z)[tid][j]-(setup->dev_Z)[0][j];
				if((delX_i == 0 && delY_i == 0 && delZ_i == 0)||(delX_j == 0 && delY_j == 0 && delZ_j == 0)){
					setup->dev_covar_vals[i][j] += 0;
				}
				else{
					numerator = delX_i*delX_j + delY_i*delY_j + delZ_i*delZ_j;
					denominator = sqrt(delX_i*delX_i + delY_i*delY_i + delZ_i*delZ_i)*sqrt(delX_j*delX_j + delY_j*delY_j + delZ_j*delZ_j);
					//printf("%f\n", numerator/denominator);
					setup->dev_covar_vals[i][j] += numerator/denominator;
					if(j!=i) setup->dev_covar_vals[j][i] += numerator/denominator;
				}
				tid += blockDim.x * gridDim.x;
			}
		}
	}
	return;
}


extern "C" int covargpu_setup(setup_data *setup){
	int num_frames = setup->end_frame-setup->start_frame+1;	
	cudaMalloc((void**)setup->dev_X, num_frames * sizeof(float**));
	cudaMalloc((void**)setup->dev_Y, num_frames * sizeof(float**));
	cudaMalloc((void**)setup->dev_Z, num_frames * sizeof(float**));
	for(int i = 0; i < num_frames; i++) {
		//Within each frame, allocate space for each amino acid
		cudaMalloc((void**)setup->dev_X[i], setup->N*sizeof(float*));
		cudaMalloc((void**)setup->dev_Y[i], setup->N*sizeof(float*));
		cudaMalloc((void**)setup->dev_Z[i], setup->N*sizeof(float*)); 
	}

	setup->covar_vals = (float **) malloc(setup->N*sizeof(float *));
	for(int i = 0; i < setup->N; i++){
		setup->covar_vals[i] = (float*) malloc(setup->N*sizeof(float));
		for(int j = 0; j < setup->N; j++){
			setup->covar_vals[i][j] = 0;
		}
	}

	cudaMalloc((void**)&setup->dev_covar_vals, setup->N*sizeof(float*));
	for(int i = 0; i < setup->N; i++){
		cudaMalloc((void**)&setup->dev_covar_vals[i], setup->N*sizeof(float));
		for(int j = 0; j < setup->N; j++){
			cudaMemcpy(&setup->dev_covar_vals[i][j], &setup->covar_vals[i][j], setup->N*sizeof(float), cudaMemcpyHostToDevice);
		}
	}
	return 0;
}

extern "C" int covargpu(setup_data *setup, float ***X, float ***Y, float ***Z){
	//Copy XYZ coordinate arrays to device
	for(int frame = setup->start_frame; frame < setup->end_frame; frame++){
			cudaMemcpy(setup->dev_X[frame], (*X)[frame], setup->N*sizeof(float), cudaMemcpyHostToDevice);
			cudaMemcpy(setup->dev_Y[frame], (*Y)[frame], setup->N*sizeof(float), cudaMemcpyHostToDevice);
			cudaMemcpy(setup->dev_Z[frame], (*Z)[frame], setup->N*sizeof(float), cudaMemcpyHostToDevice);
	}
	//Copy running sum of covar values to device
	for(int i = 0; i < setup->N; i++){
		for(int j = 0; j < setup->N; j++){
			cudaMemcpy(&setup->dev_covar_vals[i][j], &setup->covar_vals[i][j], setup->N*sizeof(float), cudaMemcpyHostToDevice);
		}
	}

	//Add to running sum
	covargpu_add<<<128,128>>>(setup);

	//Copy running sum back to host
	for(int i = 0; i < setup->N; i++){
		for(int j = 0; j < setup->N; j++){
			cudaMemcpy(&setup->covar_vals[i][j], &setup->dev_covar_vals[i][j], setup->N*sizeof(float), cudaMemcpyDeviceToHost);
		}
	}
	return 0;
}


extern "C" int covargpu_post(setup_data *setup){
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
