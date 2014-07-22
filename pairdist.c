/*pairdist.c*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "setup.h"
#include "pairdist.h"
#include "newdcdio.h"

int pairdist_setup(setup_data *setup){

	setup->RESID_map = (int *) malloc(sizeof(int)*setup->N);
	char pdb_filename[50];	
	FILE *pdb_file;

	printf("Enter desired number of pairs:");
	scanf("%d", &(setup->numpairs));
	setup->RESID1 = (int *) malloc(sizeof(int)*(setup->numpairs));
	setup->RESID2 = (int *) malloc(sizeof(int)*(setup->numpairs));	
	
	for(int i = 0; i < setup->numpairs; i++){
		printf("Pair #%d: Enter RESID 1: ", i);
		scanf("%d", &(setup->RESID1[i]));
		
		printf("Pair #%d: Enter RESID 2: ", i);
		scanf("%d", &(setup->RESID2[i]));
	}

	sprintf(pdb_filename, "%s/structures/%s.ref.pdb", setup->protein_name, setup->protein_name);
	//printf("Opening file: %s\n", pdb_filename);
	pdb_file = fopen(pdb_filename, "r");
	if(pdb_file == NULL){
		printf("Failed to open pdb file: %s", pdb_filename);
		exit(0);
	}
	char buffer[100];
	int counter = 0;	
	while(counter < setup->N){
		fgets(buffer, 100, pdb_file);
		if(!strncmp(buffer, "ATOM", 4)){
			sscanf(&(buffer[22]), "%d", &setup->RESID_map[counter]);
			//printf("%4d%4d\n", counter, setup->RESID_map[counter]);
			counter++;
		}
	}
	return 0;
}

int pairdist(setup_data *setup, float ***X, float ***Y, float ***Z){
	
	int num_frames = setup->end_frame-setup->start_frame+1;
	char dat_filename[50];
	FILE *dat_file;
	float X1, X2;
	float Y1, Y2;
	float Z1, Z2;

	for(int i = 0; i < setup->numpairs; i++){
		sprintf(dat_filename, "%s/%s/%s/pairdist_%d-%d_%d.dat", setup->protein_name, setup->sim_type, "dat", setup->RESID1[i], setup->RESID2[i], setup->runnum);
		//printf("Creating .dat file: %s\n", dat_filename);
		dat_file = fopen(dat_filename, "w");
		if(!dat_file){ 
			printf("Failed to create file: %s\n", dat_filename);
			exit(0);
		}

		int RESINDEX1 = -1, RESINDEX2 = -1;
		int temp = setup->RESID1[i];
		while(RESINDEX1 == -1 && temp >= 0){
			if((setup->RESID_map)[temp] == setup->RESID1[i]){
				RESINDEX1 = temp;
				//printf("RESINDEX1: %d\n", RESINDEX1);
			}
			//printf("%d %d %d\n", temp, setup->RESID1[i], setup->RESID_map[temp]);			
			temp--;
		}

		if(temp < 0){
			printf("RESID %d does not exist\n", setup->RESID1[i]);
			exit(0);
		}

		temp = setup->RESID2[i];
		while(RESINDEX2 == -1 && temp >= 0){
			if(setup->RESID_map[temp]==setup->RESID2[i]){
				RESINDEX2 = temp;
				//printf("RESINDEX2: %d\n", RESINDEX2);
			}
			temp--;
		}

		if(temp < 0){
			printf("RESID %d does not exist\n", setup->RESID2[i]);
			exit(0);
		}

		for(int j = 0; j < num_frames; j++){
			X1 = (*X)[j][RESINDEX1];
			X2 = (*X)[j][RESINDEX2];
			Y1 = (*Y)[j][RESINDEX1];
			Y2 = (*Y)[j][RESINDEX2];
			Z1 = (*Z)[j][RESINDEX1];
			Z2 = (*Z)[j][RESINDEX2];


			/*printf("X1: %f\n", X1);
			printf("Y1: %f\n", Y1);
			printf("Z1: %f\n", Z1);
			printf("X2: %f\n", X2);
			printf("Y2: %f\n", Y2);
			printf("Z2: %f\n", Z2);*/

			float dist = sqrt((X2-X1)*(X2-X1) + (Y2-Y1)*(Y2-Y1) + (Z2-Z1)*(Z2-Z1));
			//printf("Pair #%d, Frame %d, Distance %f\n", i, j, dist);
			fprintf(dat_file, "%6d%12f\n", j+setup->start_frame, dist);
		}
		fclose(dat_file);
	}
	return 0;
}

int pairdist_post(setup_data *setup){
		if(setup->runcount > 1){	
		char buffer[20];
		//printf("%d %d", setup.RESID1[0], setup.RESID2[0]);
		int answer_given = 0;	
		do{
			printf("\nWould you like to average all data files produced? (y/n) ");
			scanf("%s", buffer);
			if(!strcmp(buffer, "y")) answer_given = 1;
			if(!strcmp(buffer, "n")) answer_given = 1;
		} while(!answer_given);	
		if(!strcmp(buffer, "y")){
			for(int pair = 0; pair < setup->numpairs; pair++){
				FILE *files[setup->runcount];
				char dat_filename[100];
				for(int run = setup->runstart; run < setup->runstart+setup->runnum; run++){
					sprintf(dat_filename, "%s/%s/dat/pairdist_%d-%d_%d.dat", setup->protein_name, setup->sim_type, setup->RESID1[pair], setup->RESID2[pair], run);
					files[run-setup->runstart] = fopen(dat_filename, "r");
					if(files[run-setup->runstart] == NULL){
						printf("Failed to open file: %s\n", dat_filename);
						exit(0);
					}
				}
				sprintf(dat_filename, "%s/%s/dat/pairdist_%d-%d_average_%d-%d.dat", setup->protein_name, setup->sim_type, 
						setup->RESID1[pair], setup->RESID2[pair], setup->runstart, setup->runstart+setup->runnum);
				FILE *avg_dat_file = fopen(dat_filename, "w");


				for(int frame = setup->start_frame; frame <= setup->end_frame; frame++){	
					int framenumtotal = 0;
					float valuetotal = 0;
	
					for(int run = setup->runstart; run < setup->runstart+setup->runnum; run++){
						char buffer[50];
						int framenum;
						float value;
						fgets(buffer, 49, files[run - setup->runstart]);
						sscanf(buffer, "%i %f", &framenum, &value);
						//printf("frame: %d, framenum: %d\n", frame, framenum);
						framenumtotal += framenum;
						valuetotal += value;
					}
					int framenum = framenumtotal/setup->runnum;
					float value = valuetotal/setup->runnum;
					fprintf(avg_dat_file, "%6d%12f\n", framenum, value);
				}
				for(int run = setup->runstart; run < setup->runstart+setup->runnum; run++){
					sprintf(dat_filename, "%s/%s/dat/pairdist_%d-%d_%d.dat", setup->protein_name, setup->sim_type, setup->RESID1[pair], setup->RESID2[pair], run);
					if(fclose(files[run-setup->runstart]) == EOF){
						printf("Failed to close file: %s\n", dat_filename);
						exit(0);
					}
				}
			}
			char buffer[20];
			answer_given = 0;
			do{
				printf("Would you like to erase previous files? (y/n) ");	
				scanf("%s", buffer);			
				if(!strcmp(buffer, "y")) answer_given = 1;
				if(!strcmp(buffer, "n")) answer_given = 1;		
			} while(!answer_given);
			if(!strcmp(buffer, "y")){
				char rm_dat_file_cmd[100];
				for(int pair = 0; pair < setup->numpairs; pair++){
					for(int run = setup->runstart; run <= setup->runstart+setup->runnum; run++){
						sprintf(rm_dat_file_cmd, "rm %s/%s/dat/pairdist_%d-%d_%d.dat", setup->protein_name, setup->sim_type, setup->RESID1[pair], setup->RESID2[pair], run);
						system(rm_dat_file_cmd);
					}
				}
			}
		} //else printf("Exiting\n");
	}
	free(setup->RESID1);
	free(setup->RESID2);
	free(setup->RESID_map);

	return 0;
}










