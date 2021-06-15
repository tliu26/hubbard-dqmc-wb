#include "phonon_g.h"
#include "dupio.h"
#include "util.h"
#include <mkl.h>

void CalculateElectronPhononCouplingMatrix(const sim_params_t *restrict params, double *restrict g_mat)
{
	const int Norb = params->Norb;
	const int Nx   = params->Nx;
	const int Ny   = params->Ny;

	const int Ncell = Nx*Ny;
	const int N     = Norb*Ncell;

    int i, j;
    int o;

    for (o = 0; o < Norb; o++)
    {
        for (i = 0; i < Ncell; i++)
        {
            const int ix = i % Nx;
            const int iy = i / Nx;
            for (j = 0; j < Ncell; j++)
            {
                const int jx = j % Nx;
                const int jy = j / Nx;
                const int ij_ind = i + o*Ncell + (j + o*Ncell) * N;
                if (i == j)                               {g_mat[ij_ind] = params->phonon_params.g[o];}
                else if (ix == jx && iy == (jy + 1) % Ny) {g_mat[ij_ind] = params->phonon_params.gp[o];}
                else if (ix == jx && jy == (iy + 1) % Ny) {g_mat[ij_ind] = params->phonon_params.gp[o];}
                else if (iy == jy && ix == (jx + 1) % Nx) {g_mat[ij_ind] = params->phonon_params.gp[o];}
                else if (iy == jy && jx == (ix + 1) % Nx) {g_mat[ij_ind] = params->phonon_params.gp[o];}
            }
        }
    }
}

int CalculateInverseElectronPhononCouplingMatrix(const double *restrict g_mat, const int N, double *restrict g_mat_inv)
{
    memcpy(g_mat_inv, g_mat, N * N * sizeof(double));
    lapack_int *ipiv = MKL_malloc(N * sizeof(lapack_int), MEM_DATA_ALIGN);
    lapack_int info;
    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, N, g_mat_inv, N, ipiv);
    if (info < 0)
    {
        duprintf("Intel MKL 'LAPACKE_dgetrf' failed, error code: %i.\n", info);
		return info;
    }
    info = LAPACKE_dgetri(LAPACK_COL_MAJOR, N, g_mat_inv, N, ipiv);
    if (info < 0)
    {
        duprintf("Intel MKL 'LAPACKE_dgetri' failed, error code: %i.\n", info);
		return info;
    }
    MKL_free(ipiv);
    return info;
}

void TransformPhononCoordinates(const double *restrict g_mat_inv, const double *restrict X, const int N, double *restrict X_old)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            X_old[i] += g_mat_inv[i + j*N] * X[j];
        }
    }
}

double CalculatePhononPotentialEnergyWithGP(const int i, const int Ncell, const int N, const phonon_params_t *restrict phonon_params,
    const double dx, const double *restrict X_l,
    const double *restrict g_mat_inv)
{
    double *X_l_old  = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);
    double *dx_l_old = (double *)MKL_malloc(N*sizeof(double), MEM_DATA_ALIGN);
    TransformPhononCoordinates(g_mat_inv, X_l, N, X_l_old);
    int j;
    for (j = 0; j < N; j++)
    {
        dx_l_old[j] = dx * g_mat_inv[j + i*N];
    }
    double dPEph = 0;
    for (j = 0; j < N; j++)
    {
        const int o = j / Ncell;
        dPEph += 0.5 * square(phonon_params->omega[o]) * dx_l_old[j] * (dx_l_old[j] + 2*X_l_old[j]);
    }
    MKL_free(X_l_old);
    MKL_free(dx_l_old);
    return dPEph;
}

double CalculatePhononEnergygp(const int i, const int Ncell, const int N, const phonon_params_t *restrict phonon_params,
    const double inv_dt_sq, const double dx, const double *restrict X_l, const double *restrict X_l_next_old, const double *restrict X_l_prev_old,
    const double *restrict g_mat_inv)
{
    double *X_l_old  = (double *)MKL_calloc(N, sizeof(double), MEM_DATA_ALIGN);
    double *dx_l_old = (double *)MKL_malloc(N*sizeof(double), MEM_DATA_ALIGN);
    TransformPhononCoordinates(g_mat_inv, X_l, N, X_l_old);
    int j;
    for (j = 0; j < N; j++)
    {
        dx_l_old[j] = dx * g_mat_inv[j + i*N];
    }
    double dEph = 0;
    for (j = 0; j < N; j++)
    {
        const int o = j / Ncell;
        dEph += dx_l_old[j] * (0.5*square(phonon_params->omega[o])*(dx_l_old[j] + 2*X_l_old[j]) + inv_dt_sq*(dx_l_old[j] - (X_l_next_old[j] - 2*X_l_old[j] + X_l_prev_old[j])));
    }
    MKL_free(X_l_old);
    MKL_free(dx_l_old);
    return dEph;
}