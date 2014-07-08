/*=====================================*\
 * DCD and PDB Analysis Functions      *                            
 * Author: William Cornell             *
 * Date: 7/1/14                        *
 * Tehver Research Group               *
\*=====================================*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "setup.h"
#include "newdcdio.h"
#include "pairdist.h"
#include "covar.h"
#include "covargpu.h"
#include "pdbio.h"


int main(int argc, char *argv[]){

	/*==============================================COMMAND LINE============================================================================*/
	if(argc < 3){
		printf("Insufficient number of command line arguments given\n");
		exit(0);
	}
	//Command line arguments (4 and 5 optional): $ ./analysis <protein_name> <sim_type> <setup_file> <analysis_function>
	setup_data setup; //struct from setup.h
	setup.protein_name = argv[1];
	//printf("Protein Name: %s\n", setup.protein_name);
	setup.sim_type = argv[2];
	//printf("Sim Type: %s\n", setup.sim_type);
	if(argc < 5){ //if no fifth argument of function to run
		char buffer[20];
		int func_given = 0;
		while(!func_given){
			printf("Enter program to run (or type 'help' for more options): ");
			scanf("%s", buffer);
			if(!strcmp(buffer, "help")){
				printf("\nPrograms available:\n\n");
				printf("\tpairdist -- asks for pairs of RESIDs and creates .dat files of their relative positions as a function of time\n\n");
				printf("\tcovar    -- creates a .dat file of the covariance matrix to 6 decimal places\n\n");
				printf("\tcovargpu -- GPU accelerated version of covar\n\n");
				/*Add help documentation for new functions here*/
				printf("For more documentation, see the README included with the source code\n\n");
			} 
			if(!strcmp(buffer, "pairdist")){
				sprintf(setup.function_name, "pairdist");
				func_given = 1;
			}
			if(!strcmp(buffer, "covar")){
				sprintf(setup.function_name, "covar");
				func_given = 1;
			}
			if(!strcmp(buffer, "covargpu")){
				sprintf(setup.function_name, "covargpu");
				func_given = 1;
			}
			if(!strcmp(buffer, "pdbwrite")){
				sprintf(setup.function_name, "pdbwrite");
				func_given = 1;
			}
			/*Add conditions for additional functions here*/
		}
	} else sprintf(setup.function_name, "%s", argv[4]);

	if(!strcmp(setup.function_name, "pairdist")){
		setup.analysis_setup_funcptr = &pairdist_setup;
		setup.analysis_funcptr = &pairdist;
		setup.analysis_post_funcptr = &pairdist_post;
	} 
	if(!strcmp(setup.function_name, "covar")){
		setup.analysis_setup_funcptr = &covar_setup;
		setup.analysis_funcptr = &covar;
		setup.analysis_post_funcptr = &covar_post;
	}
	if(!strcmp(setup.function_name, "covargpu")){
		setup.analysis_setup_funcptr = &covargpu_setup;
		setup.analysis_funcptr = &covargpu;
		setup.analysis_post_funcptr = &covargpu_post;
	}
	if(!strcmp(setup.function_name, "pdbwrite")){
		setup.analysis_setup_funcptr = &pdb_setup;
		setup.analysis_funcptr = &pdb_write;
		setup.analysis_post_funcptr = &pdb_post;
	}
	/*Add additional conditionals and function pointer assignment for new functions here*/

	if(argc > 3) read_setup_file(argv[3], &setup);
	else{
		printf("Enter starting run number: ");
		scanf("%d", &(setup.runstart));
		printf("Enter number of runs: ");
		scanf("%d", &(setup.runcount));


		printf("Enter desired starting frame: ");
		scanf("%d", &(setup.start_frame));
		printf("Enter desired ending frame: ");
		scanf("%d", &(setup.end_frame));
	}

	/*===================================================READ HEADER================================================================*/

    //READS the header to obtain number of atoms, stored in setup.N
	sprintf(setup.dcd_filename, "%s/%s/dcd/%s_%d_%s.dcd", setup.protein_name, setup.sim_type, setup.protein_name, setup.runstart, setup.sim_type);
	FILE *header_file = fopen(setup.dcd_filename, "r");
	if(header_file == NULL){
		printf("Failed to open file: %s\n", setup.dcd_filename);
		exit(0);
	}
	dcd_read_header_N(header_file, &setup);
	fclose(header_file);
	//printf("N: %d\n", setup.N);

	/*====================================================ANALYSIS===================================================================*/

	setup.analysis_setup_funcptr(&setup); //call function specific one time setup
	//printf("%d %d", setup.RESID1[0], setup.RESID2[0]);
	clock_t start = clock(), diff; //inititalize a few variables to time analysis
	
	for(int i = setup.runstart; i < setup.runstart + setup.runcount; i++){ //loop over each run
		sprintf(setup.dcd_filename, "%s/%s/dcd/%s_%d_%s.dcd", setup.protein_name, setup.sim_type, setup.protein_name, i, setup.sim_type);
		setup.runnum = i;
		printf("\nReading file: %s\n", setup.dcd_filename);
		dcd_single_analysis(&setup);
	}

	
	diff = clock() - start; //difference in time
	int msec = diff*1000/CLOCKS_PER_SEC;

	printf("\nAnalysis completed in %d.%d seconds (CPU time)\n", msec/1000, msec%1000);

	setup.analysis_post_funcptr(&setup);

	/*===============================================================================================================================*/

	return 0;
}























