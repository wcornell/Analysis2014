/*newdcdio.c*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "setup.h"
#include "pairdist.h"
#include "newdcdio.h"



void dcd_read_header(FILE* dcd_file, char *dcd_filename, int* N, int* NFILE, int* NPRIV, int* NSAVC, float* DELTA){
	int iin; 
	char cin[4];
	fread(&iin, 4, 1, dcd_file); //check first four bytes
	if(iin != 84){
		printf("Error! Wrong DCD file: DCD supposed to have '84' in first 4 bytes, but it hasn't.\n");
		exit(0);
	}
	fread(&cin, 4, 1, dcd_file); //check next four bytes
	/*if(!((cin[0]=="C")&&(cin[1]=="O")&&(cin[2]=="R")&&(cin[3]=="D")){
	if(strcmp(cin, cord_char)){ //check that those four bytes are CORD*
		printf("Error! Wrong DCD file: no 'CORD' sequence at the beginning of the file. Found: %s \n", cin);
		exit(0);
	}*/

	fread(NFILE, 4, 1, dcd_file); 
	fread(NPRIV, 4, 1, dcd_file);
	fread(NSAVC, 4, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	for(int i = 0; i < 5; i++){
		fread(&iin, 4, 1, dcd_file);
	}
	fread(DELTA, 4, 1, dcd_file);
	for(int i = 0; i < 9; i++){
		fread(&iin, 4, 1, dcd_file);
	}
	//24
	fread(&iin, 4, 1, dcd_file);
	//84;
	fread(&iin, 4, 1, dcd_file);
	//164;
	fread(&iin, 4, 1, dcd_file);
	//2;
	fread(&iin, 4, 1, dcd_file);
	char title[80];
	fread(&title, 80, 1, dcd_file);
	//printf("Title1: %s\n", title);
	fread(&title, 80, 1, dcd_file);
	//printf("Title2: %s\n", title);
	//164
	fread(&iin, 4, 1, dcd_file);
	//4
	fread(&iin, 4, 1, dcd_file);
	//N
	fread(N, 4, 1, dcd_file);
	//4
	fread(&iin, 4, 1, dcd_file);
}

int dcd_read_header_N(FILE *header_file, setup_data *setup){
	int NFILE;
	int NPRIV;
	int NSAVC;
	float DELTA;
	char *dcd_filename;

	dcd_read_header(header_file, dcd_filename, &setup->N, &NFILE, &NPRIV, &NSAVC, &DELTA);

	return 0;
}

int dcd_read_frame(FILE* dcd_file, int N, float *X, float *Y, float *Z){
	int iin;
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(X, 4*N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(Y, 4*N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(Z, 4*N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	if(feof(dcd_file) == 0){
		return 0;
	} else {
		return -1;
	}
}

int dcd_read_file(FILE *dcd_file, int N, float ***X, float ***Y, float ***Z, int start_frame, int end_frame){
	int file_offset = (24 + 12*N)*start_frame; //How many bytes into file to start reading
	fseek(dcd_file, file_offset, SEEK_CUR); //set location in file to beginning of start frame

	for(int i = 0; i < end_frame-start_frame+1; i++){
		if(dcd_read_frame(dcd_file, N, (*X)[i], (*Y)[i], (*Z)[i])){
			printf("Not enough frames in file. Stopped at last frame, #%d\n", i+start_frame);
			exit(0);
		}
		//printf("Read and stored frame #%d\n", i+start_frame);
		
	}
	return 0;
}


int dcd_single_analysis(setup_data *setup){

	int num_frames = setup->end_frame-setup->start_frame+1;	
	int NFILE; //(not used in this implementation)
	int NPRIV; //Starting timestep of DCD file - NOT ZERO (not used in this implementation)
	int NSAVC; //Timesteps between DCD saves (not used in this implementation)
	float DELTA; //length of a timestep (not used in this implementation)

	float **X; 
	float **Y;
	float **Z;
	
	/*-------------------------GET HEADER----------------------------*/
	FILE* dcd_file = fopen(setup->dcd_filename, "r"); //Open DCD file
	if(dcd_file == NULL){
		printf("Unable to open DCD file: %s\n", setup->dcd_filename);
		exit(0);
	}
	dcd_read_header(dcd_file, setup->dcd_filename, &setup->N, &NFILE, &NPRIV, &NSAVC, &DELTA);	
	//printf("File header read\n");
	

	/*----------------------MEMORY ALLOCATION------------------------*/
	//allocate space for XYZ position for each frame
	X = (float**) malloc(num_frames*sizeof(float**));
	Y = (float**) malloc(num_frames*sizeof(float**));
	Z = (float**) malloc(num_frames*sizeof(float**)); 
	for(int i = 0; i < num_frames; i++) {
		//Within each frame, allocate space for each amino acid
		X[i] = (float*) malloc(setup->N*sizeof(float*));
		Y[i] = (float*) malloc(setup->N*sizeof(float*));
		Z[i] = (float*) malloc(setup->N*sizeof(float*)); 
	}

	//printf("Memory allocated for %d frames\n", num_frames);

	/*-------------------------READ DATA-----------------------------*/
	dcd_read_file(dcd_file, setup->N, &X, &Y, &Z, setup->start_frame, setup->end_frame);
	//printf("File read\n");
	//if(!(fflush(dcd_file))) printf("buffer flushed\n");
	if(fclose(dcd_file)){ //Close DCD file
		printf("Failed to close file: %s\n", setup->dcd_filename);
		exit(0);
	} 
	//printf("%d\n", fclose(dcd_file));
    printf("File closed\n");
	
	/*printf("Frame 0 X190: %f\n", X[0][190]);
	printf("Frame 0 Y190: %f\n", Y[0][190]);
	printf("Frame 0 Z190: %f\n", Z[0][190]);*/

	/*--------------------------ANALYSIS-----------------------------*/
	setup->analysis_funcptr(setup, &X, &Y, &Z);
	//printf("exited analysis_func\n");

	/*------------------------DEALLOCATION---------------------------*/
	for(int i = 0; i < num_frames; i++){
		free(X[i]);
		free(Y[i]);
		free(Z[i]);
	}
	free(X);
	free(Y);
	free(Z);
	
	//printf("Memory deallocated\n");	
	return 0;
}


























