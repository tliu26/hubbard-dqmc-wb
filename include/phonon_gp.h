#ifndef PHONON_GP_H
#define PHONON_GP_H
#include "sim_params.h"
#include <mkl.h>

lapack_int ComputeXUpdateMatrixLU(const int i, const int N, const int Ncell, const int n, const int *restrict phonon_neighbors,
    const double delta_i, const double delta_j, const double *restrict G, double *restrict partial_gg_mat, lapack_int *restrict ipiv);

double ComputeXUpdateDeterminant(const int n, const double *restrict partial_gg_mat, const lapack_int *restrict ipiv);

void GreenXWoodburyUpdate(const int i, const int N, const int Ncell, const int n, const int *restrict phonon_neighbors,
    const double delta_i, const double delta_j, double *restrict G, const double *restrict partial_g_mat, lapack_int *restrict ipiv);

void ConstructPhononNeigborsMap(const int Nx, const int Ny, int *restrict phonon_neighbors);

void PrintMatrix(const int m, const int n, const char *name, const double *restrict mat);

#endif