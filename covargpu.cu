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
	int num_frames = setup->end_frame-setup->start_frame+1;	

	float delX_i, delY_i, delZ_i;
	float delX_j, delY_j, delZ_j;

	float numerator, denominator;
    for (int i = 0; i < setup->N; i++) {
		while(tid < setup->N){
        	for (int frame = 0; frame <= num_frames; frame++) {
				delX_i = (setup->dev_X)[frame*setup->N+i]-(setup->dev_X)[i];
				delY_i = (setup->dev_Y)[frame*setup->N+i]-(setup->dev_Y)[i];
				delZ_i = (setup->dev_Z)[frame*setup->N+i]-(setup->dev_Z)[i];
				delX_j = (setup->dev_X)[frame*setup->N+tid]-(setup->dev_X)[tid];
				delY_j = (setup->dev_Y)[frame*setup->N+tid]-(setup->dev_Y)[tid];
				delZ_j = (setup->dev_Z)[frame*setup->N+tid]-(setup->dev_Z)[tid];
				if((delX_i == 0 && delY_i == 0 && delZ_i == 0)||(delX_j == 0 && delY_j == 0 && delZ_j == 0)){
					setup->dev_covar_vals[i + tid*setup->N] += 0;
				}
				else{
					numerator = delX_i*delX_j + delY_i*delY_j + delZ_i*delZ_j;
					denominator = sqrt(delX_i*delX_i + delY_i*delY_i + delZ_i*delZ_i)*sqrt(delX_j*delX_j + delY_j*delY_j + delZ_j*delZ_j);
					setup->dev_covar_vals[i*setup->N + tid] += numerator/denominator;
					if(tid!=i) setup->dev_covar_vals[i*setup->N + tid] += numerator/denominator;
				}
			}
		tid += blockDim.x * gridDim.x;
		}
	}
	return;
}


extern "C" int covargpu_setup(setup_data *setup){
	//Get GPU from user
	int dev;
	printf("Enter GPU ID: "); fflush(stdout);
	scanf("%d", &dev);
	cudaSetDevice(dev);
	//printf("Setting up device memory\n");

	//Allocate space for linearized XYZ arrays
	int num_frames = setup->end_frame-setup->start_frame+1;	
	cudaMalloc((void**)&setup->dev_X, num_frames*sizeof(float)*setup->N);
	cudaMalloc((void**)&setup->dev_Y, num_frames*sizeof(float)*setup->N);
	cudaMalloc((void**)&setup->dev_Z, num_frames*sizeof(float)*setup->N); 
	//printf("XYZ allocated\n");

	//Allocate host covar vals, initialize to 0
	setup->lin_covar_vals = (float *) malloc(setup->N*setup->N*sizeof(float));
    for (int i = 0; i < setup->N*setup->N; ++i) {
        setup->lin_covar_vals[i] = 0;
    }

	//Allocate device covar vals
	cudaMalloc((void**)&setup->dev_covar_vals, setup->N*setup->N*sizeof(float*));
	//printf("Allocation complete\n");
	return 0;
}

extern "C" int covargpu(setup_data *setup, float ***X, float ***Y, float ***Z){
	//Allocate linearized, contiguous arrays for XYZ (num_frames x N)
	int num_frames = setup->end_frame-setup->start_frame+1;	
	float *linX = (float*) malloc(setup->N*num_frames*sizeof(float));	
	float *linY = (float*) malloc(setup->N*num_frames*sizeof(float));
	float *linZ = (float*) malloc(setup->N*num_frames*sizeof(float));
	//printf("lin declared\n");
	for(int i = setup->start_frame; i < setup->end_frame; i++){
		for(int j = 0; j < setup->N; j++){
			linX[setup->N*(i-setup->start_frame) + j] = (*X)[i][j];
			linY[setup->N*(i-setup->start_frame) + j] = (*Y)[i][j];
			linZ[setup->N*(i-setup->start_frame) + j] = (*Z)[i][j];
		}
	}	//printf("Linear matrices created\n");
	
	//Copy XYZ coordinate arrays to device
	cudaMemcpy(linX, &setup->dev_X, setup->N*num_frames*sizeof(float), cudaMemcpyHostToDevice);	
	cudaMemcpy(linY, &setup->dev_Y, setup->N*num_frames*sizeof(float), cudaMemcpyHostToDevice);	
	cudaMemcpy(linZ, &setup->dev_Z, setup->N*num_frames*sizeof(float), cudaMemcpyHostToDevice);
	//printf("XYZ copied\n");

	//Copy running sum of covar values to device
	cudaMemcpy(setup->lin_covar_vals, setup->dev_covar_vals, setup->N*setup->N*sizeof(float), cudaMemcpyDeviceToHost);
	//printf("Covar matrix copied\n");

	//Add to running sum
	covargpu_add<<<128,128>>>(setup);

	//Copy running sum back to host
	cudaMemcpy(setup->dev_covar_vals, setup->lin_covar_vals, setup->N*setup->N*sizeof(float), cudaMemcpyHostToDevice);
	//printf("Covar matrix updated\n");

	//Free host linear arrays
	free(linX);
	free(linY);
	free(linZ);

	return 0;
}


extern "C" int covargpu_post(setup_data *setup){

	//Create/Open .dat file
	char dat_filename[50];
	FILE *dat_file;
	sprintf(dat_filename, "%s/%s/dat/covar_%s_%d-%d_%d-%d.dat", setup->protein_name, setup->sim_type, setup->protein_name, setup->start_frame, setup->end_frame, setup->runstart, setup->runstart+setup->runnum);
	dat_file = fopen(dat_filename, "w");
	if(dat_file == NULL){
		printf("Failed to open file: %s", dat_filename);
		exit(0);
	}

	//Print matrix to .dat file
	float iterations = (setup->runcount*(setup->end_frame-setup->start_frame));
	for(int i = 0; i < setup->N; i++){
		for(int j = 0; j < setup->N; j++){
			setup->lin_covar_vals[i + j*setup->N] /= iterations;
			fprintf(dat_file, "%10.6f", setup->lin_covar_vals[i + j*setup->N]);
		}
	fprintf(dat_file, "\n");
	}
	fclose(dat_file);

	//Free memory
	cudaFree(setup->dev_X);
	cudaFree(setup->dev_Y);
	cudaFree(setup->dev_Z);
	cudaFree(setup->dev_covar_vals);
	free(setup->lin_covar_vals);
	return 0;
}
