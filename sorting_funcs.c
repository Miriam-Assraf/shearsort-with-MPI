/*
 ============================================================================
 Author      : Miriam Assraf
 Description : Shear sort functions
 ============================================================================
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "sorting_funcs.h"
#include "constants.h"

void shearSort(int n, int my_rank, int upRank, int bottomRank, int leftRank, int rightRank, int my_coords[], double *cuboid, MPI_Comm comm) {
	for (int i = 0; i <= ceil(log2(n) + 1); i++) {
		// even steps sort cols
		if (i % 2 == 0) {
			oddEvenSort(my_rank, cuboid, bottomRank, upRank, my_coords[0], LEFT_TO_RIGHT, n, comm);	// only increasing order
		}
		// odd steps sort rows
		else {
			if (my_coords[0] % 2 == 0) {	// even rows increasing order
				oddEvenSort(my_rank, cuboid, rightRank, leftRank, my_coords[1], LEFT_TO_RIGHT, n, comm);
			} else {	// odd rows decreasing order
				oddEvenSort(my_rank, cuboid, rightRank, leftRank, my_coords[1], RIGHT_TO_LEFT, n, comm);
			}
		}
	}
}

void oddEvenSort(int my_rank, double *my_cube, int next_proc, int previous_proc, int position, int order, int n, MPI_Comm comm) {
	for (int i = 0; i <= n; i++) {
		// even steps
		if (i % 2 == 0) {
			evenStep(my_cube, next_proc, previous_proc, position, order, n, comm);
		}
		// odd steps
		else {
			oddStep(my_cube, next_proc, previous_proc, position, order, n, comm);
		}
	}
}

void evenStep(double *my_cube, int next_proc, int previous_proc, int position, int order, int n, MPI_Comm comm) {
	double my_value = my_cube[1] * my_cube[2] * my_cube[3];
	double neighbor_value;
	double neighbor_cube[n];
	MPI_Status status;

	// even row/column process exchange values with right/bottom neighbor
	if (position % 2 == 0 && position != (n - 1)) {	// without last one- doesn't have right/bottom neighbor
		MPI_Recv(neighbor_cube, n, MPI_DOUBLE, next_proc, 0, comm, &status);	// get cube from right/bottom neighbor
		neighbor_value = neighbor_cube[1] * neighbor_cube[2] * neighbor_cube[3];
		if (my_value == neighbor_value) {	// if same volume, compare by height
			my_value = my_cube[3];
			neighbor_value = neighbor_cube[3];
		}
		my_value = exchangeWithNext(order, n, my_value, neighbor_value);	// compare by given order
		if (my_value == neighbor_value) {	// was exchanged
			MPI_Send(my_cube, n, MPI_DOUBLE, next_proc, 0, comm);	// send to neighbor my_cube to exchange
			memcpy(my_cube, neighbor_cube, sizeof(double) * n);
		} else {	// wasn't exchanged
			MPI_Send(neighbor_cube, n, MPI_DOUBLE, next_proc, 0, comm);	// send to neighbor it's own cube
		}
	}
	// odd row/column process exchange values with left/upper neighbor
	else {
		MPI_Send(my_cube, n, MPI_DOUBLE, previous_proc, 0, comm);	// send my_cube to left/upper neighbor
		MPI_Recv(neighbor_cube, n, MPI_DOUBLE, previous_proc, 0, comm, &status);	// receive cube for exchange from left/upper neighbor
		memcpy(my_cube, neighbor_cube, sizeof(double) * n);
	}
}

void oddStep(double *my_cube, int next_proc, int previous_proc, int position, int order, int n, MPI_Comm comm) {
	double my_value = my_cube[1] * my_cube[2] * my_cube[3];
	double neighbor_value;
	double neighbor_cube[n];
	MPI_Status status;

	// even row/column process exchange values with left/upper neighbor
	if (position % 2 == 0) {
		if (position != 0) {	// first diesn't have left/upper neighbor
			MPI_Recv(neighbor_cube, n, MPI_DOUBLE, previous_proc, 0, comm, &status);	// receive cube from left/upper neighbor
			neighbor_value = neighbor_cube[1] * neighbor_cube[2] * neighbor_cube[3];
			if (my_value == neighbor_value) {	// if same volume, compare by height
				my_value = my_cube[3];
				neighbor_value = neighbor_cube[3];
			}
			my_value = exchangeWithPrevious(order, n, my_value, neighbor_value);	// compare by given order
			if (my_value == neighbor_value) {	// was exchanged
				MPI_Send(my_cube, n, MPI_DOUBLE, previous_proc, 0, comm);	// send to neighbor my_cube to exchange
				memcpy(my_cube, neighbor_cube, sizeof(double) * n);
			} else {
				MPI_Send(neighbor_cube, n, MPI_DOUBLE, previous_proc, 0, comm);	// send to neighbor it's own cube
			}
		}
	}
	// odd row/column process exchange values with right/bottom neighbor
	else {
		if (position != (n - 1)) {
			MPI_Send(neighbor_cube, n, MPI_DOUBLE, next_proc, 0, comm);	// send my_cube to right/bottom neighbor
			MPI_Recv(neighbor_cube, n, MPI_DOUBLE, next_proc, 0, comm, &status);	// receive cube for exchange from right/bottom neighbor
			memcpy(my_cube, neighbor_cube, sizeof(double) * n);
		}
	}
}

double exchangeWithNext(int order, int n, double my_volume, double neighbor_volume) {
	if (order == LEFT_TO_RIGHT) { // if increasing order-gets the smaller value
		if (my_volume > neighbor_volume) {
			my_volume = neighbor_volume;
		}
	} else {	// if decreasing order- gets the greater value
		if (my_volume < neighbor_volume) {
			my_volume = neighbor_volume;
		}
	}
	return my_volume;
}

double exchangeWithPrevious(int order, int n, double my_volume, double neighbor_volume) {
	if (order == LEFT_TO_RIGHT) {	// if increasing order- gets greater value
		if (my_volume < neighbor_volume) {
			my_volume = neighbor_volume;
		}
	} else {	// if decreasing order- gets smaller value
		if (my_volume > neighbor_volume) {
			my_volume = neighbor_volume;
		}
	}
	return my_volume;
}
