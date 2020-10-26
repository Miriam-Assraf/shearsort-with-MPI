/*
 ============================================================================
 Name        : ex2.c
 Author      : Miriam Assraf
 Description : Shear sort using odd-even sort for cuboids matrix
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"

#define ROOT 0
#define LEFT_TO_RIGHT 0
#define RIGHT_TO_LEFT 1

const int ROWS = 16;
const int COLS = 4;

void readCuboids(double cuboids[][COLS]);
void printResults(int n, double cuboids[][COLS]);
MPI_Comm cartCreate(int n, int my_rank, int *upRank, int *bottomRank, int *leftRank, int *rightRank, int my_coords[]);
void shearSort(int n, int my_rank, int upRank, int bottomRank, int leftRank, int rightRank, int my_coords[], double *cuboid, MPI_Comm comm);
void oddEvenSort(int my_rank, double *my_cube, int next_proc, int previous_proc, int position, int order, int n, MPI_Comm comm);
void evenStep(double *my_cube, int next_proc, int previous_proc, int position, int order, int n, MPI_Comm comm);
void oddStep(double *my_cube, int next_proc, int previous_proc, int position, int order, int n, MPI_Comm comm);
double exchangeWithNext(int order, int n, double my_volume, double neighbor_volume);
double exchangeWithPrevious(int order, int n, double my_volume, double neighbor_volume);

int main(int argc, char *argv[]) {
	int my_rank;
	int num_processes;
	int upRank, bottomRank, leftRank, rightRank, my_coords[2];
	int n = 4;
	double cuboids[ROWS][COLS];
	double *my_cube = (double*) malloc(sizeof(double) * n);
	MPI_Comm comm2D;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

	if(num_processes!=ROWS){	// make sure number of processes matches number of lines = number of cubes
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

void readCuboids(double cuboids[][COLS]) {
	FILE *inputCuboiudsFile = fopen("cuboids.dat.txt", "r");
	int i = 0;

	if (inputCuboiudsFile == NULL) {
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	while (fscanf(inputCuboiudsFile, "%lf %lf %lf %lf", &cuboids[i][0], &cuboids[i][1], &cuboids[i][2], &cuboids[i][3]) != EOF) {
		i++;
	}
	fclose(inputCuboiudsFile);
}

void printResults(int n, double cuboids[][COLS]) {
	for (int i = 0; i < n * n; i++) {
		printf("cube %0.0lf and volume %0.2lf\n", cuboids[i][0],
				cuboids[i][1] * cuboids[i][2] * cuboids[i][3]);
	}
	FILE *outputResultsFile = fopen("result.dat.txt", "w");

	for (int i = 0; i < n; i++) {
		fprintf(outputResultsFile, "%0.0lf ", cuboids[i][0]);
	}
	for (int j = 2 * n - 1; j >= n; j--) {
		fprintf(outputResultsFile, "%0.0lf ", cuboids[j][0]);
	}

	for (int i = 2 * n; i < 3 * n; i++) {
		fprintf(outputResultsFile, "%0.0lf ", cuboids[i][0]);
	}

	for (int i = 4 * n - 1; i >= 3 * n; i--) {
		fprintf(outputResultsFile, "%0.0lf ", cuboids[i][0]);
	}

	fclose(outputResultsFile);
}

MPI_Comm cartCreate(int n, int my_rank, int *upRank, int *bottomRank, int *leftRank, int *rightRank, int my_coords[]) {
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
