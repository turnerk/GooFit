#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char** argv) {

	MPI_Init(NULL, NULL); // Initialize the MPI environment
	
  int world_size;
  int world_rank;

	// Get the number of processes  
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  // Get the rank of the processes 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int num_elements = 10;
  int num_elements_per_proc = num_elements / world_size;

  double *nums;

  if (world_rank == 0) {
    nums = new double[num_elements];
    for (int i = 0; i < num_elements; i++) {
    	nums[i] = (double) i+1;
    }
  }
 
	int *sendcounts = new int[world_size];  // array describing how many elements to send to each process
  int *displs = new int[world_size];  // array describing the displacements where each segment begins

  //indexing for copying events over!
	for (int i = 0; i < world_size - 1; i++) {
  	sendcounts[i] = num_elements_per_proc;
	}
  sendcounts[world_size - 1] = num_elements - num_elements_per_proc*(world_size - 1); 

  //displacements into the array for indexing!
  displs[0] = 0;
  for (int i = 1; i < world_size; i++) {
    displs[i] = displs[i - 1] + sendcounts[i - 1];
	}

	double *sub_nums = new double[sendcounts[world_rank]];

  MPI_Scatterv(nums, sendcounts, displs, MPI_DOUBLE, sub_nums,
              sendcounts[world_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (world_rank == 0) {
    printf("\nNum array: ");
    for (int i = 0; i < world_size; i++) {
      for (int j = 0; j < 5; j++) {
        printf("%f ", nums[displs[i]+j]);
      }
    }
  }

  // Print the numbers scattered to each processor
  printf("\nProcessor rank %i out of %i processors: ", world_rank, world_size);

  for (int i = 0; i < sendcounts[world_rank]; i++) {
        printf("%f ", sub_nums[i]); 
  }

	if (world_rank == 0) 
    delete [] nums;
	delete [] sub_nums;

  // Finalize the MPI environment
  MPI_Finalize();

	delete [] sendcounts;
	delete [] displs;

	if (world_rank == world_size -1) 
    printf("\n");

  return 0;
}
