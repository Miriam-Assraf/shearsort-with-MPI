/*
 ============================================================================
 Name        : ex2.c
 Author      : Miriam Assraf
 Description : Shear sort using odd-even sort for cuboids matrix
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "sorting_funcs.h"
#include "input_output_funcs.h"
#include "constants.h"

MPI_Comm cartCreate(int n, int my_rank, int* upRank, int* bottomRank, int* leftRank, int* rightRank, int my_coords[]);

int main(int argc, char* argv[]) {
	int my_rank;
	int num_processes;
	int upRank, bottomRank, leftRank, rightRank, my_coords[2];
	int n = COLS;
	double cuboids[ROWS][COLS];
	double* my_cube = (double*)malloc(sizeof(double) * n);
	MPI_Comm comm2D;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

	if (num_processes != ROWS) {	// make sure number of processes matches number of lines = number of cubes
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (my_rank == ROOT) {	// root process read all cuboids from file
		readCuboids(cuboids);
	}

	comm2D = cartCreate(n, my_rank, &upRank, &bottomRank, &leftRank, &rightRank, my_coords);

	// scatter rows of matrix to different processes- each process get one line (one cube)
	MPI_Scatter(&cuboids[0][0], n, MPI_DOUBLE, &my_cube[0], n, MPI_DOUBLE, ROOT, comm2D);

	shearSort(n, my_rank, upRank, bottomRank, leftRank, rightRank, my_coords, my_cube, comm2D);

	// gather results to root process
	MPI_Gather(&my_cube[0], n, MPI_DOUBLE, &cuboids[0][0], n, MPI_DOUBLE, ROOT, comm2D);

	if (my_rank == ROOT) {	// root process writes results of sorted matrix
		printResults(n, cuboids);
	}

	MPI_Finalize();
	return 0;
}

MPI_Comm cartCreate(int n, int my_rank, int* upRank, int* bottomRank, int* leftRank, int* rightRank, int my_coords[]) {
	int dim[2], period[2];
	int reorder;
	MPI_Comm comm2D;

	dim[0] = n;
	dim[1] = n;
	period[0] = 0;
	period[1] = 0;
	reorder = 0;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm2D);

	// Each process displays its rank and cartesian coordinates
	MPI_Cart_coords(comm2D, my_rank, 2, my_coords);
	fflush(stdout);
	MPI_Cart_shift(comm2D, 1, -1, rightRank, leftRank);
	MPI_Cart_shift(comm2D, 0, -1, bottomRank, upRank);

	return comm2D;
}