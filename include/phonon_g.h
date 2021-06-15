#ifndef PHONON_G_H
#define PHONON_G_H

#include "sim_params.h"

void CalculateElectronPhononCouplingMatrix(const sim_params_t *restrict params, double *restrict g_mat);

int CalculateInverseElectronPhononCouplingMatrix(const double *restrict g_mat, const int N, double *restrict g_mat_inv);

void TransformPhononCoordinates(const double *restrict g_mat_inv, const double *restrict X, const int N, double *restrict X_old);

double CalculatePhononEnergygp(const int i, const int Ncell, const int N, const phonon_params_t *restrict phonon_params,
    const double inv_dt_sq, const double dx, const double *restrict X_l, const double *restrict X_l_next_old, const double *restrict X_l_prev_old,
    const double *restrict g_mat_inv);

double CalculatePhononPotentialEnergyWithGP(const int i, const int Ncell, const int N, const phonon_params_t *restrict phonon_params,
    const double dx, const double *restrict X_l,
    const double *restrict g_mat_inv);


#endif