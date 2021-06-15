#include "phonon_gp.h"
#include "dupio.h"
#include <string.h>

lapack_int ComputeXUpdateMatrixLU(const int i, const int N, const int Ncell, const int n, const int *restrict phonon_neighbors,
    const double delta_i, const double delta_j, const double *restrict G, double *restrict partial_gg_mat, lapack_int *restrict ipiv)
{
    const int o = i / Ncell;
    const int i_cell = i % Ncell;
    int k1, k2;
    const int *phonon_i_neighbors = &phonon_neighbors[i_cell*n];
    for (k1 = 0; k1 < n; k1++)
    {
        for (k2 = 0; k2 < n; k2++)
        {
            int k2k1 = k2 + k1 * n;
            const double g_k2k1 = G[phonon_i_neighbors[k2] + o*Ncell + (phonon_i_neighbors[k1] + o*Ncell)*N];
            if (k2 == k1) { partial_gg_mat[k2k1] = 1 + (1 - g_k2k1) * (k2 == 0 ? delta_i : delta_j); }
            else          { partial_gg_mat[k2k1] =         -g_k2k1  * (k2 == 0 ? delta_i : delta_j); }
        }
    }
    lapack_int info;
    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, partial_gg_mat, n, ipiv);
    return info;
}

double ComputeXUpdateDeterminant(const int n, const double *restrict partial_gg_mat, const lapack_int *restrict ipiv)
{
    int j;
    double det = 1;
    for (j = 0; j <= NUM_NEIGHBORS; j++) { det *= partial_gg_mat[j + j * n] * (ipiv[j] == (j + 1) ? 1 : -1); }
    return det;
}

void GreenXWoodburyUpdate(const int i, const int N, const int Ncell, const int n, const int *restrict phonon_neighbors,
    const double delta_i, const double delta_j, double *restrict G, const double *restrict partial_g_mat, lapack_int *restrict ipiv)
{
    double *gU = (double *)MKL_malloc(N * n * sizeof(double), MEM_DATA_ALIGN);
    double *Vg = (double *)MKL_malloc(N * n * sizeof(double), MEM_DATA_ALIGN);
    double *pgm = (double *)MKL_malloc(n * n * sizeof(double), MEM_DATA_ALIGN);
    memcpy(pgm, partial_g_mat, n * n *sizeof(double));
    const int i_cell = i % Ncell;
    const int o = i / Ncell;
    const int *phonon_i_neighbors = &phonon_neighbors[i_cell*n];
    int k1;
    for (k1 = 0; k1 < n; k1++)
    {
        int j = phonon_i_neighbors[k1] + o*Ncell;
        int i_tmp;
        for (i_tmp = 0; i_tmp < N; i_tmp++)
        {
            gU[i_tmp + k1*N] = G[i_tmp + j*N];
            if (i_tmp == j) { Vg[k1 + i_tmp*n] = (k1 == 0 ? (1 - G[j + i_tmp*N]) * delta_i : (1 - G[j + i_tmp*N]) * delta_j); }
            else            { Vg[k1 + i_tmp*n] = (k1 == 0 ?    - G[j + i_tmp*N]  * delta_i :    - G[j + i_tmp*N]  * delta_j); }
        }
    }
    lapack_int info = LAPACKE_dgetri(LAPACK_COL_MAJOR, n, pgm, n, ipiv);
    if (info < 0) { duprintf("Intel MKL 'LAPACKE_dgetri' failed, error code: %i.\n", info); }
    double *T = (double *)MKL_malloc(N * n * sizeof(double), MEM_DATA_ALIGN);
    cblas_dgemm(LAPACK_COL_MAJOR, CblasNoTrans, CblasNoTrans, n, N, n,  1.0, pgm, n, Vg, n, 0.0, T, n);
    cblas_dgemm(LAPACK_COL_MAJOR, CblasNoTrans, CblasNoTrans, N, N, n, -1.0,  gU, N,  T, n, 1.0, G, N);
    MKL_free(T);
    MKL_free(gU);
    MKL_free(Vg);
    MKL_free(pgm);
}

void ConstructPhononNeigborsMap(const int Nx, const int Ny, int *restrict phonon_neighbors)
{
    const int n = NUM_NEIGHBORS + 1;
    int i;
    for (i = 0; i < Nx*Ny; i++)
    {
        const int ix = i % Nx;
        const int iy = i / Nx;
        phonon_neighbors[0 + i*n] = i;
        phonon_neighbors[1 + i*n] = (ix < Nx - 1 ? ix + 1 : 0     ) + iy * Nx;
        phonon_neighbors[2 + i*n] = (ix > 0      ? ix - 1 : Nx - 1) + iy * Nx;
        phonon_neighbors[3 + i*n] = ix + (iy < Ny - 1 ? iy + 1 : 0     ) * Nx;
        phonon_neighbors[4 + i*n] = ix + (iy > 0      ? iy - 1 : Ny - 1) * Nx;
    }
}

void PrintMatrix(const int m, const int n, const char *name, const double *restrict mat)
{
    duprintf("%s = \n", name);
    int i;
    for (i = 0; i < m; i++)
    {
        int j;
        for (j = 0; j < n; j++)
        {
            duprintf("%.12f ", mat[i + j*m]);
        }
        duprintf("\n");
    }
}