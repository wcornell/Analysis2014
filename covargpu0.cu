/*covargpu0.cu*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "setup.h"
#include "covargpu0.h"
#include "newdcdio.h"

__global__ void covargpu0_add(float *covar_vals, float *X, float *Y, float *Z, int N, int num_frames){

	float delX_i, delY_i, delZ_i;
	float delX_j, delY_j, delZ_j;

	float numerator, denominator;
	int tidx = threadIdx.x + blockIdx.x * blockDim.x;
	while(tidx < N){
		int tidy = threadIdx.y + blockIdx.y * blockDim.y;
		while(tidy < N){
        	for (int frame = 0; frame < num_frames; frame++) {
				if(tidx <= tidy){					
					int index1 = frame*N+tidx, index2 = frame*N+tidy;
					delX_i = X[index1]-X[tidx];
					delY_i = Y[index1]-Y[tidx];
					delZ_i = Z[index1]-Z[tidx];
					delX_j = X[index2]-X[tidy];
					delY_j = Y[index2]-Y[tidy];
					delZ_j = Z[index2]-Z[tidy];
					if((delX_i == 0 && delY_i == 0 && delZ_i == 0)||(delX_j == 0 && delY_j == 0 && delZ_j == 0)){
						covar_vals[tidx*N + tidy] += 0;
					}
					else{
						numerator = delX_i*delX_j + delY_i*delY_j + delZ_i*delZ_j;
						denominator = sqrtf((delX_i*delX_i + delY_i*delY_i + delZ_i*delZ_i)*(delX_j*delX_j + delY_j*delY_j + delZ_j*delZ_j));
						covar_vals[tidx*N + tidy] += numerator/denominator;
						if(tidx!=tidy) covar_vals[tidy*N + tidx] += numerator/denominator;
					}
				}
			}
			tidy += blockDim.y * gridDim.y;
		}
		tidx += blockDim.x * gridDim.x;
	}
	return;
}



extern "C" int covargpu0_setup(setup_data *setup){
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

	sprintf(setup->dcd_filename, "%s/%s/dcd/%s_%d_%s.dcd", setup->protein_name, setup->sim_type, setup->protein_name, setup->runstart, setup->sim_type);
	FILE *dcd_file = fopen(setup->dcd_filename, "r");
	if(dcd_file == NULL){
		printf("Failed to open file: %s\n", setup->dcd_filename);
		exit(0);
	}
	float tempXYZ[setup->N*3];
	int iin;
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(tempXYZ, 4*setup->N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(&tempXYZ[setup->N], 4*setup->N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(&tempXYZ[setup->N*2], 4*setup->N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fclose(dcd_file);

	//Allocate host covar vals, initialize to 0
	setup->lin_covar_vals = (float *) malloc(setup->N*setup->N*sizeof(float));
	for(int i = 0; i < setup->N*setup->N; i++) setup->lin_covar_vals[i] = 0;

	//Allocate device covar vals, copy 0 initialized matrix
	cudaMalloc((void**)&setup->dev_covar_vals, setup->N*setup->N*sizeof(float*));
	cudaMemcpy(setup->dev_covar_vals, setup->lin_covar_vals, setup->N*setup->N*sizeof(float), cudaMemcpyHostToDevice);


	//printf("Allocation complete\n");
	return 0;
}


extern "C" int covargpu0(setup_data *setup, float ***X, float ***Y, float ***Z){
	//Allocate linearized, contiguous arrays for XYZ (num_frames x N)
	int num_frames = setup->end_frame-setup->start_frame+1;	
	float *linX = (float*) malloc(setup->N*num_frames*sizeof(float));	
	float *linY = (float*) malloc(setup->N*num_frames*sizeof(float));
	float *linZ = (float*) malloc(setup->N*num_frames*sizeof(float));
	//printf("lin declared\n");
	for(int i = setup->start_frame; i <= setup->end_frame; i++){
		for(int j = 0; j < setup->N; j++){
			linX[setup->N*(i-setup->start_frame) + j] = (*X)[i][j];
			linY[setup->N*(i-setup->start_frame) + j] = (*Y)[i][j];
			linZ[setup->N*(i-setup->start_frame) + j] = (*Z)[i][j];
		}
	}	//printf("Linear matrices created\n");


	//Copy XYZ coordinate arrays to device
	cudaMemcpy(setup->dev_X, linX, setup->N*num_frames*sizeof(float), cudaMemcpyHostToDevice);	
	cudaMemcpy(setup->dev_Y, linY, setup->N*num_frames*sizeof(float), cudaMemcpyHostToDevice);	
	cudaMemcpy(setup->dev_Z, linZ, setup->N*num_frames*sizeof(float), cudaMemcpyHostToDevice);
	//printf("XYZ copied\n");

	//for(int i = 0; i < setup->N*num_frames; i++) printf("%8.3f%8.3f%8.3f\n", linX[i], linY[i], linZ[i]);

	//setup timing variables
	cudaEvent_t start, stop;
	float elapsed_time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	//block and thread dimension setup
	dim3 threadsPerBlock(32, 32);
	dim3 numBlocks(32, 32);

	//Call kernel to add to running sum
	covargpu0_add<<<numBlocks,threadsPerBlock>>>(setup->dev_covar_vals, setup->dev_X, setup->dev_Y, setup->dev_Z, setup->N, num_frames);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsed_time, start, stop);
	printf("Covariance calculated in %3.2f s\n", elapsed_time/1000);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);


	//Free host linear arrays
	free(linX);
	free(linY);
	free(linZ);

	return 0;
}



extern "C" int covargpu0_post(setup_data *setup){

	//Create/Open .dat file
	char dat_filename[50];
	FILE *dat_file;
	sprintf(dat_filename, "%s/%s/dat/covargpu0_%s_%d-%d_%d-%d.dat", setup->protein_name, setup->sim_type, setup->protein_name, setup->start_frame, setup->end_frame, setup->runstart, setup->runstart+setup->runnum);
	dat_file = fopen(dat_filename, "w");
	if(dat_file == NULL){
		printf("Failed to open file: %s", dat_filename);
		exit(0);
	}

	//Copy running sum back to host
	cudaMemcpy(setup->lin_covar_vals, setup->dev_covar_vals, setup->N*setup->N*sizeof(float), cudaMemcpyDeviceToHost);

	//Print matrix to .dat file
	int iterations = (setup->runcount*(setup->end_frame-setup->start_frame));	
	for(int i = 0; i < setup->N; i++){
		for(int j = 0; j < setup->N; j++){
			fprintf(dat_file, "%10.6f", setup->lin_covar_vals[i*setup->N + j]/iterations);
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
