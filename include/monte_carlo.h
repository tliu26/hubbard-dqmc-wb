#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "kinetic.h"
#include "greens_func.h"
#include "random.h"
#include "sim_params.h"
#include "measurement.h"


void DQMCIteration(const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params,
	const int nwraps, randseed_t *restrict seed, spin_field_t *restrict s,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd,
	const int neqlt, measurement_data_t *restrict meas_data);

void DQMCPhononIteration(const double dt, const double mu, const kinetic_t *restrict kinetic, const int noHS, const stratonovich_params_t *restrict stratonovich_params, const phonon_params_t *restrict phonon_params,
	const int nwraps, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX, const int *restrict phonon_neighbors,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd,
	const int neqlt, measurement_data_t *restrict meas_data, measurement_data_phonon_t *restrict meas_data_phonon,
	const sim_params_t *restrict params, const char *fnbase);


void PhononBlockUpdates(const double dt, const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params,
	const phonon_params_t *restrict phonon_params, randseed_t *restrict seed, const spin_field_t *restrict s, double *restrict X, double *restrict expX, const int *restrict phonon_neighbors,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd,
	int *n_block_accept, int *n_block_total);


void PhononFlipUpdates(const double dt, const double mu, const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params,
	const phonon_params_t *restrict phonon_params, randseed_t *restrict seed, const spin_field_t *restrict s, double *restrict X, double *restrict expX, const int *restrict phonon_neighbors,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd,
	int *n_flip_accept, int *n_flip_total);


void DQMCSimulation(const sim_params_t *restrict params,
	measurement_data_t *restrict meas_data, measurement_data_unequal_time_t *restrict meas_data_uneqlt, measurement_data_phonon_t *restrict meas_data_phonon,
	int *restrict i_start, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX, const int *restrict phonon_neighbors, const char *fnbase);



#endif
