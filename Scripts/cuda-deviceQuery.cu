/*cuda-deviceQuery.cu*/
#include <stdio.h>
#include <stdlib.h>

int main(void){
	cudaDeviceProp prop;
	int count;
	
	cudaGetDeviceCount(&count);
	printf("There are %d CUDA compatible devices available\n\n", count);
	for(int i = 0; i < count; i++){
		cudaGetDeviceProperties(&prop, i);
		printf("--------------------  Device %d --------------------\n", i);
		printf("  General Information\n");
		printf("    Name: %s\n", prop.name);
		printf("    Compute capability: %d.%d\n", prop.major, prop.minor);
		printf("    Clock rate: %d\n", prop.clockRate);
		printf("  Memory Information\n");
		printf("    Total global memory: %ld\n", prop.totalGlobalMem);
		printf("    Total constant memory: %ld\n", prop.totalConstMem);
		printf("  Multi-processor information\n");
		printf("    Multiprocessor count: %d\n", prop.multiProcessorCount);
		printf("    Shared memory per mp: %ld\n", prop.sharedMemPerBlock);
		printf("    Threads in warp: %d\n", prop.warpSize);
		printf("    Max threads per block: %d\n", prop.maxThreadsPerBlock);
		printf("    Max thread dimensions: (%d, %d, %d)\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
		printf("    Max Grid dimensions: (%d, %d, %d)\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
		printf("----------------------------------------------------\n\n");
	}
	return 0;
}
