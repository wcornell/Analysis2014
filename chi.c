/*chi.c*/

/*
int chi_setup(setup_data *setup){
	FILE *ref_pdb;

	ref_pdb = fopen(setup->ref_pdb_filename, "r");
	if(ref_pdb == NULL){
		printf("Failed to open file: %s", ref_pdb_filename);
		exit(0);
	}

	setup->refX = (float*) malloc(setup->N*sizeof(float));
	setup->refY = (float*) malloc(setup->N*sizeof(float));
	setup->refZ = (float*) malloc(setup->N*sizeof(float));
	
	int counter = 0;
	char tempstr[50];
	int temp;
	while(counter < setup->N){
		fgets(buffer, 100, setup->pdb_file);
		if(!strncmp(buffer, "ATOM", 4)){
			sscanf(buffer, "%s %d %s %s %s %d %f %f %f %d %d", tempstr, temp, tempstr, tempstr, tempstr, temp,
				   setup->refX[counter], setup->refY[counter], setup->refZ[counter]);
			counter++;
		}
	}
	return 0;
	
}

int chi(double **dist0){   // dist0 is the reference (PDB-based) pair-distance matrix

	int diff=0, total=0;
	double dist, dx, dy, dz, sum, chi;

	for(int i = 1 ; i <= tot_amino - 1; i ++){
		for(int j = i + 1 ; j <= tot_amino ; j ++){
			dx = Amino[i].x - Amino[j].x;
			dy = Amino[i].y - Amino[j].y;
			dz = Amino[i].z - Amino[j].z;
			dist = sqrt(dx*dx+dy*dy+dz*dz);
			sum = fabs(dist - dist0[i][j]);
		if(sum <= R_limit_chi) diff ++;   // R_limit_chi is a paramater in def_param.h, its value is always 2.0 
		total ++;
		}
	}
	chi = 1.0 - ((double) diff)/((double) total);
	return 0;
}

int chi_post(setup_data *setup){
	free(setup->refX);
	free(setup->refY);
	free(setup->refZ);
	return 0;
}*/
