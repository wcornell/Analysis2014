/*pdbio.c*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "setup.h"



int pdb_setup(setup_data *setup){
	setup->line_beginning = (char **) malloc(setup->N*sizeof(char *));
	for(int i = 0; i < setup->N; i++){
		setup->line_beginning[i] = (char *) malloc(27*sizeof(char));
	}
	//char pdb_filename[50];
	char buffer[100];
	//sprintf(pdb_filename, "%s/structures/%s.ref.pdb", setup->protein_name, setup->protein_name);
	setup->pdb_file = fopen(setup->ref_pdb_filename, "r");
	if(!setup->pdb_file){
		printf("Failed to open file: %s\n", setup->ref_pdb_filename);
		exit(0);
	}

	int counter = 0;	
	while(counter < setup->N){
		fgets(buffer, 100, setup->pdb_file);
		if(!strncmp(buffer, "ATOM", 4)){
			strncpy(setup->line_beginning[counter], buffer, 26);
			setup->line_beginning[counter][26] = 0;
			counter++;
		}
	}
	return 0;
}

int pdb_write(setup_data *setup, float ***X, float ***Y, float ***Z){
	FILE *new_pdb_file;
	char new_pdb_filename[100];
	for(int frame = setup->start_frame; frame <= setup->end_frame; frame++){	
		sprintf(new_pdb_filename, "%s/%s/dat/%s_%d_%d.pdb", setup->protein_name, setup->sim_type, setup->protein_name, frame, setup->runnum);
		new_pdb_file = fopen(new_pdb_filename, "w");
		if(new_pdb_file == NULL){
			printf("Failed to open pdb file: %s\n", new_pdb_filename);
			exit(0);
		}
		for(int i = 0; i < setup->N; i++){
			fprintf(new_pdb_file, "%s%11.3f%8.3f%8.3f\n", setup->line_beginning[i], (*X)[0][i], (*Y)[0][i], (*Z)[0][i]);
		}
	}
	return 0;
}

int pdb_post(setup_data *setup){
	for(int i = 0; i < setup->N; i++){
		free(setup->line_beginning[i]);
	}
	free(setup->line_beginning);
	return 0;
}
