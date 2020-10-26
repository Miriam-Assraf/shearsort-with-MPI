#pragma once
void shearSort(int n, int my_rank, int upRank, int bottomRank, int leftRank, int rightRank, int my_coords[], double* cuboid, MPI_Comm comm);
void oddEvenSort(int my_rank, double* my_cube, int next_proc, int previous_proc, int position, int order, int n, MPI_Comm comm);
void evenStep(double* my_cube, int next_proc, int previous_proc, int position, int order, int n, MPI_Comm comm);
void oddStep(double* my_cube, int next_proc, int previous_proc, int position, int order, int n, MPI_Comm comm);
double exchangeWithNext(int order, int n, double my_volume, double neighbor_volume);
double exchangeWithPrevious(int order, int n, double my_volume, double neighbor_volume);