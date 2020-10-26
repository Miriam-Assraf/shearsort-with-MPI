/*
 ============================================================================
 Author      : Miriam Assraf
 Description : Read input and write output
 ============================================================================
 */
#include <stdio.h>
#include "mpi.h"
#include "input_output_funcs.h"

void readCuboids(double cuboids[][COLS]) {
	FILE* inputCuboiudsFile = fopen("cuboids.dat.txt", "r");
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
	// print to console
	for (int i = 0; i < n * n; i++) {
		printf("cube %0.0lf and volume %0.2lf\n", cuboids[i][0],
			cuboids[i][1] * cuboids[i][2] * cuboids[i][3]);
	}

	// print to file in the right order - sorted in shearsort
	// even rows from left to right, odd rows from right to left
	FILE* outputResultsFile = fopen("result.dat.txt", "w");

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