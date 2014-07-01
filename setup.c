/*setup.c*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "setup.h"



int read_setup_file(char *setup_filename, setup_data *setup){
	FILE *setup_file = fopen(setup_filename, "r");
	if(setup_file == NULL){
		printf("Failed to open setup file: %s", setup_filename);
		exit(0);
	}
	else{
		char buffer1[20];
		char buffer2[20];
		fgets(buffer1, 19, setup_file);
		sscanf(buffer1, "%s\t%d", buffer2, &(setup->runstart));
		//printf("runstart: %d\n", setup->runstart);
		fgets(buffer1, 19, setup_file);
		sscanf(buffer1, "%s\t%d", buffer2, &(setup->runcount));
		//printf("runcount: %d\n", setup->runcount);
		fgets(buffer1, 19, setup_file);
		sscanf(buffer1, "%s\t%d", buffer2, &(setup->start_frame));
		//printf("start_frame: %d\n", setup->start_frame);
		fgets(buffer1, 19, setup_file);
		sscanf(buffer1, "%s\t%d", buffer2, &(setup->end_frame));
		//printf("end_frame: %d\n", setup->end_frame);
	}
	return 0;
}
