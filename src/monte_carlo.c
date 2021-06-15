#include "monte_carlo.h"
#include "greens_func.h"
#include "util.h"
#include "dupio.h"
#include "profiler.h"
#include "checkpoint.h"
#include "progress.h"
#include "phonon_gp.h"
#include <mkl.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <stdint.h>


//________________________________________________________________________________________________________________________
///
/// \brief Perform a Determinant Quantum Monte Carlo (DQMC) iteration
///
/// \param kinetic              matrix exponential of kinetic energy operator
/// \param stratonovich_params  precomputed Hubbard-Stratonovich parameters
/// \param nwraps               number of "time slice wraps" before recomputing the Green's function; must be a multiple of 'prodBlen'
/// \param seed                 random number "seed" structure, will be updated during function call
/// \param s                    Hubbard-Stratonovich field; will be updated
/// \param tsm_u                time step matrices for the spin-up Green's function, will be updated
/// \param tsm_d                time step matrices for the spin-down Green's function, will be updated
/// \param Gu                   spin-up   Green's function, must have been computed on input and will be updated
/// \param Gd                   spin-down Green's function, must have been computed on input and will be updated
/// \param neqlt                perform an equal time measurement every 'neqlt' time slices
/// \param meas_data            measurement data
///
void DQMCIteration(const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params, const int nwraps,
	randseed_t *restrict seed, spin_field_t *restrict s, time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d,
	greens_func_t *restrict Gu, greens_func_t *restrict Gd, const int neqlt, measurement_data_t *restrict meas_data)
{
	Profile_Begin("DQMCIter");
	__assume_aligned(s, MEM_DATA_ALIGN);

	// dimension consistency checks
	assert(tsm_u->N == kinetic->Ncell * kinetic->Norb);
	assert(tsm_u->N == tsm_d->N);
	assert(tsm_u->L == tsm_d->L);
	const int N = tsm_u->N;
	const int Ncell = kinetic->Ncell;
	const int L = tsm_u->L;

	assert(tsm_u->prodBlen == tsm_d->prodBlen);
	assert(nwraps % tsm_u->prodBlen == 0);  // must be a multiple of 'prodBlen'

	// store Green's functions before recomputing them to estimate error
	#if defined(DEBUG) | defined(_DEBUG)
	greens_func_t Gu_old, Gd_old;
	AllocateGreensFunction(N, &Gu_old);
	AllocateGreensFunction(N, &Gd_old);
	__assume_aligned(Gu_old.mat, MEM_DATA_ALIGN);
	__assume_aligned(Gd_old.mat, MEM_DATA_ALIGN);
	#endif

	// random shuffle of lattice cells and orbitals
	int *orb_cell_order = MKL_malloc(N * sizeof(int), MEM_DATA_ALIGN);
	__assume_aligned(orb_cell_order, MEM_DATA_ALIGN);

	// iterate over time slices
	int l;
	for (l = 0; l < L; l++)
	{
		Profile_Begin("DQMCIter_Wraps");
		#pragma omp parallel sections
		{
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_u->B[l], tsm_u->invB[l], Gu->mat);
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_d->B[l], tsm_d->invB[l], Gd->mat);
		}
		Profile_End("DQMCIter_Wraps");

		// iterate over lattice sites in random order, updating the Hubbard-Stratonovich field
		Profile_Begin("DQMCIter_SiteUpdate");
		Random_Shuffle(seed, N, orb_cell_order);
		int j;
		for (j = 0; j < N; j++)
		{
			const int i = orb_cell_order[j];
			const int o = i / Ncell;    // orbital index
			assert(0 <= i && i < N);
			assert(0 <= o && o < kinetic->Norb);

			// Eq. (13)
			// suggest flipping s_{i,l}
			const double du = 1 + (1 - Gu->mat[i + i*N]) * stratonovich_params->delta[  s[i + l*N]][o];
			const double dd = 1 + (1 - Gd->mat[i + i*N]) * stratonovich_params->delta[1-s[i + l*N]][o];
			if (Random_GetUniform(seed) < fabs(du*dd))
			{
				// Eq. (15)
				#pragma omp parallel sections
				{
					#pragma omp section
					GreenShermanMorrisonUpdate(stratonovich_params->delta[  s[i + l*N]][o], N, i, Gu->mat);
					#pragma omp section
					GreenShermanMorrisonUpdate(stratonovich_params->delta[1-s[i + l*N]][o], N, i, Gd->mat);
				}
				// correspondingly update determinants
				Gu->logdet -= log(fabs(du));
				Gd->logdet -= log(fabs(dd));
				if (du < 0) { Gu->sgndet = -Gu->sgndet; }
				if (dd < 0) { Gd->sgndet = -Gd->sgndet; }

				// actually flip spin of Hubbard-Stratonovich field entry
				s[i + l*N] = 1 - s[i + l*N];
			}
		}
		Profile_End("DQMCIter_SiteUpdate");

		// re-compute corresponding B matrices
		Profile_Begin("DQMCIter_Brecomp");
		#pragma omp parallel sections
		{
			#pragma omp section
			UpdateTimeStepMatrices(kinetic, stratonovich_params->expVu, s, l, tsm_u);
			#pragma omp section
			UpdateTimeStepMatrices(kinetic, stratonovich_params->expVd, s, l, tsm_d);
		}
		Profile_End("DQMCIter_Brecomp");

		// recompute Green's function after several time slice "wraps"
		if ((l + 1) % nwraps == 0)
		{
			Profile_Begin("DQMCIter_Grecomp");
			// store current Green's function matrices to compare with newly constructed ones
			#if defined(DEBUG) | defined(_DEBUG)
			CopyGreensFunction(Gu, &Gu_old);
			CopyGreensFunction(Gd, &Gd_old);
			#endif

			#pragma omp parallel sections
			{
				#pragma omp section
				GreenConstruct(tsm_u, (l + 1) % L, Gu);
				#pragma omp section
				GreenConstruct(tsm_d, (l + 1) % L, Gd);
			}

			#if defined(DEBUG) | defined(_DEBUG)
			// deviation of matrix entries
			double err_u = UniformDistance(N*N, Gu_old.mat, Gu->mat);
			double err_d = UniformDistance(N*N, Gd_old.mat, Gd->mat);
			if (err_u > 1e-8*N || err_d > 1e-8*N) {
				duprintf("Warning: after calling 'GreenConstruct()', largest entrywise distance between previous and new Green's functions is %g (up) and %g (down).\n", err_u, err_d);
			}
			// deviation of matrix determinants
			err_u = fabs(Gu_old.logdet - Gu->logdet);
			err_d = fabs(Gd_old.logdet - Gd->logdet);
			if (err_u > 1e-8 || err_d > 1e-8) {
				duprintf("Warning: after calling 'GreenConstruct()', largest distance between logarithm of previous and new Green's function determinants is %g (up) and %g (down).\n", err_u, err_d);
			}
			if (Gu_old.sgndet != Gu->sgndet || Gd_old.sgndet != Gd->sgndet) {
				duprintf("Warning: after calling 'GreenConstruct()', determinant sign has changed.\n");
			}
			#endif
			Profile_End("DQMCIter_Grecomp");
		}

		if (neqlt > 0 && (l + 1) % neqlt == 0)
		{
			// accumulate equal time "measurement" data
			Profile_Begin("DQMCIter_AccumulateEqMeas");
			AccumulateMeasurement(Gu, Gd, meas_data);
			Profile_End("DQMCIter_AccumulateEqMeas");
		}
	}

	// clean up
	MKL_free(orb_cell_order);
	#if defined(DEBUG) | defined(_DEBUG)
	DeleteGreensFunction(&Gd_old);
	DeleteGreensFunction(&Gu_old);
	#endif
	Profile_End("DQMCIter");
}



//________________________________________________________________________________________________________________________
///
/// \brief Perform block updates of the phonon field
///
/// \param dt                   imaginary-time step
/// \param kinetic              matrix exponential of kinetic energy operator
/// \param stratonovich_params  precomputed Hubbard-Stratonovich parameters
/// \param phonon_params        phonon field parameters
/// \param seed                 random number "seed" structure, will be updated during function call
/// \param s                    Hubbard-Stratonovich field, remains constant during phonon block updates
/// \param X                    phonon field, will be updated
/// \param expX                 entrywise exponential of the phonon field: exp(-dt g X)
/// \param tsm_u                time step matrices for the spin-up Green's function, will be updated
/// \param tsm_d                time step matrices for the spin-down Green's function, will be updated
/// \param Gu                   spin-up Green's function, will be updated
/// \param Gd                   spin-down Green's function, will be updated
/// \param n_block_accept       number of accepted block updates (updated by function)
/// \param n_block_total        number of total block updates (updated by function)
///
void PhononBlockUpdates(const double dt, const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params,
	const phonon_params_t *restrict phonon_params, randseed_t *restrict seed, const spin_field_t *restrict s, double *restrict X, double *restrict expX, const int *restrict phonon_neighbors,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd,
	int *n_block_accept, int *n_block_total)
{
	__assume_aligned(   X, MEM_DATA_ALIGN);
	__assume_aligned(expX, MEM_DATA_ALIGN);

	// dimension consistency checks
	assert(tsm_u->N == kinetic->Ncell * kinetic->Norb);
	assert(tsm_u->N == tsm_d->N);
	assert(tsm_u->L == tsm_d->L);
	const int N = tsm_u->N;
	const int Ncell = kinetic->Ncell;
	const int L = tsm_u->L;

	// fast return for zero block updates
	if (phonon_params->n_block_updates <= 0) {
		return;
	}

	// store X_{i,l} and corresponding exponential for all 'l'
	double *X_i    = (double *)MKL_malloc(L * sizeof(double), MEM_DATA_ALIGN);
	double *expX_i = (double *)MKL_malloc(L * sizeof(double), MEM_DATA_ALIGN);
	double *expX_js = (double *)MKL_malloc(L * NUM_NEIGHBORS * sizeof(double), MEM_DATA_ALIGN);
	__assume_aligned(X_i,    MEM_DATA_ALIGN);
	__assume_aligned(expX_i, MEM_DATA_ALIGN);
	__assume_aligned(expX_js, MEM_DATA_ALIGN);

	// storage for new Green's functions
	greens_func_t Gu_new, Gd_new;
	AllocateGreensFunction(N, &Gu_new);
	AllocateGreensFunction(N, &Gd_new);

	// new time step B matrices after the phonon field update
	time_step_matrices_t tsm_u_new;
	time_step_matrices_t tsm_d_new;
	AllocateTimeStepMatrices(N, L, tsm_u->prodBlen, &tsm_u_new);
	AllocateTimeStepMatrices(N, L, tsm_d->prodBlen, &tsm_d_new);

	int n;
	for (n = 0; n < phonon_params->n_block_updates; n++)
	{
		// randomly select a lattice site
		int i = (int)(Random_GetBoundedUint(seed, N));
		int o = i / Ncell;
		int j = i % Ncell;
		// ignore sites without phonon coupling
		if (phonon_params->g[o] == 0 && phonon_params->gp[o] == 0)
		{
			continue;
		}

		// backup X_{i,l} and corresponding exponential for all time slices
		int l;
		for (l = 0; l < L; l++)
		{
			   X_i[l] =    X[i + l*N];
			expX_i[l] = expX[i + l*N];
			int jn;
			for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
			{
				expX_js[l + L*(jn-1)] = expX[phonon_neighbors[jn + j*(NUM_NEIGHBORS+1)] + o*Ncell + l*N];
			}
		}

		// suggest a simultaneous shift of X_{i,l} for all 'l'
		const double dx = (Random_GetUniform(seed) - 0.5) * phonon_params->block_box_width;

		// calculate change of the phonon (lattice) energy
		double dEph = 0;
		double dEph1 = 0;
		for (l = 0; l < L; l++)
		{
			dEph += (dx + 2*X[i + l*N]);
			if (phonon_params->J != 0)
			{
				const double xpx = X[phonon_neighbors[1 + j * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				const double xmx = X[phonon_neighbors[2 + j * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				const double xpy = X[phonon_neighbors[3 + j * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				const double xmy = X[phonon_neighbors[4 + j * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				dEph1 += 2 * (4 * X[i + l*N] - xpx - xmx - xpy - xmy) + 4 * dx;
			}
		}
		dEph *= 0.5*square(phonon_params->omega[o]) * dx;
		if (phonon_params->J != 0)
		{
			dEph += (phonon_params->J / 2) * dx * dEph1;
		}

		// actually shift phonon field entries
		for (l = 0; l < L; l++)
		{
			   X[i + l*N] += dx;
			expX[i + l*N] *= exp(-dt*phonon_params->g[o] * dx);
			// expX[i + l*N] = exp(-dt * phonon_params->g[o] * X[i + l*N]);
			int jn;
			for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
			{
				expX[phonon_neighbors[jn + j*(NUM_NEIGHBORS+1)] + o*Ncell + l*N] *= exp(-dt * phonon_params->gp[o] * dx);
			}
		}

		// calculate new time step matrices
		#pragma omp parallel sections
		{
			#pragma omp section
			InitPhononTimeStepMatrices(kinetic, stratonovich_params->expVu, s, expX, &tsm_u_new);
			#pragma omp section
			InitPhononTimeStepMatrices(kinetic, stratonovich_params->expVd, s, expX, &tsm_d_new);
		}

		// calculate new Green's functions
		#pragma omp parallel sections
		{
			#pragma omp section
			GreenConstruct(&tsm_u_new, 0, &Gu_new);
			#pragma omp section
			GreenConstruct(&tsm_d_new, 0, &Gd_new);
		}

		// decide whether block update is accepted; note that det(G) = 1/det(M)
		(*n_block_total)++;
		if (Random_GetUniform(seed) < exp((Gu->logdet + Gd->logdet) - (Gu_new.logdet + Gd_new.logdet) - dt * dEph))
		{
			(*n_block_accept)++;

			// copy new stuff
			CopyGreensFunction(&Gu_new, Gu);
			CopyGreensFunction(&Gd_new, Gd);
			CopyTimeStepMatrices(&tsm_u_new, tsm_u);
			CopyTimeStepMatrices(&tsm_d_new, tsm_d);
		}
		else
		{
			// undo changes
			int l;
			for (l = 0; l < L; l++)
			{
				   X[i + l*N] =    X_i[l];
				expX[i + l*N] = expX_i[l];
				int jn;
				for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
				{
					expX[phonon_neighbors[jn + j*(NUM_NEIGHBORS+1)] + o*Ncell + l*N] = expX_js[l + L*(jn-1)];
				}
			}
		}
	}

	// clean up
	DeleteTimeStepMatrices(&tsm_d_new);
	DeleteTimeStepMatrices(&tsm_u_new);
	DeleteGreensFunction(&Gd_new);
	DeleteGreensFunction(&Gu_new);
	MKL_free(expX_i);
	MKL_free(X_i);
	MKL_free(expX_js);
}



//________________________________________________________________________________________________________________________
///
/// \brief Perform flip updates of the phonon field
///
/// \param dt                   imaginary-time step
/// \param mu                   chemical potential
/// \param kinetic              matrix exponential of kinetic energy operator
/// \param stratonovich_params  precomputed Hubbard-Stratonovich parameters
/// \param phonon_params        phonon field parameters
/// \param seed                 random number "seed" structure, will be updated during function call
/// \param s                    Hubbard-Stratonovich field, remains constant during phonon block updates
/// \param X                    phonon field, will be updated
/// \param expX                 entrywise exponential of the phonon field: exp(-dt g X)
/// \param tsm_u                time step matrices for the spin-up Green's function, will be updated
/// \param tsm_d                time step matrices for the spin-down Green's function, will be updated
/// \param Gu                   spin-up Green's function, will be updated
/// \param Gd                   spin-down Green's function, will be updated
/// \param n_flip_accept        number of accepted flip updates (updated by function)
/// \param n_flip_total         number of total flip updates (updated by function)
///
void PhononFlipUpdates(const double dt, const double mu, const kinetic_t *restrict kinetic, const stratonovich_params_t *restrict stratonovich_params,
	const phonon_params_t *restrict phonon_params, randseed_t *restrict seed, const spin_field_t *restrict s, double *restrict X, double *restrict expX, const int *restrict phonon_neighbors,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd,
	int *n_flip_accept, int *n_flip_total)
{
	__assume_aligned(   X, MEM_DATA_ALIGN);
	__assume_aligned(expX, MEM_DATA_ALIGN);

	// dimension consistency checks
	assert(tsm_u->N == kinetic->Ncell * kinetic->Norb);
	assert(tsm_u->N == tsm_d->N);
	assert(tsm_u->L == tsm_d->L);
	const int N = tsm_u->N;
	const int Ncell = kinetic->Ncell;
	const int L = tsm_u->L;

	// fast return for zero block updates
	if (phonon_params->n_block_updates <= 0) {
		return;
	}

	// store X_{i,l} and corresponding exponential for all 'l'
	double *X_ref    = (double *)MKL_malloc(L * sizeof(double), MEM_DATA_ALIGN);
	double *expX_ref = (double *)MKL_malloc(L * sizeof(double), MEM_DATA_ALIGN);
	double *expX_js = (double *)MKL_malloc(L * NUM_NEIGHBORS * sizeof(double), MEM_DATA_ALIGN);
	__assume_aligned(X_ref,    MEM_DATA_ALIGN);
	__assume_aligned(expX_ref, MEM_DATA_ALIGN);
	__assume_aligned(expX_js, MEM_DATA_ALIGN);

	// storage for new Green's functions
	greens_func_t Gu_new, Gd_new;
	AllocateGreensFunction(N, &Gu_new);
	AllocateGreensFunction(N, &Gd_new);

	// new time step B matrices after the phonon field update
	time_step_matrices_t tsm_u_new;
	time_step_matrices_t tsm_d_new;
	AllocateTimeStepMatrices(N, L, tsm_u->prodBlen, &tsm_u_new);
	AllocateTimeStepMatrices(N, L, tsm_d->prodBlen, &tsm_d_new);

	int n;
	for (n = 0; n < phonon_params->n_block_updates; n++)
	{
		// randomly select a lattice site
		const int i = (int)(Random_GetBoundedUint(seed, N));
		const int o = i / Ncell;
		int j = i % Ncell;
		// ignore sites without phonon coupling
		if (phonon_params->g[o] == 0 && phonon_params->gp[o] == 0)
		{
			continue;
		}

		// backup X_{i,l} and corresponding exponential for all time slices
		int l;
		for (l = 0; l < L; l++)
		{
			   X_ref[l] =    X[i + l*N];
			expX_ref[l] = expX[i + l*N];
			int jn;
			for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
			{
				expX_js[l + L*(jn-1)] = expX[phonon_neighbors[jn + j*(NUM_NEIGHBORS+1)] + o*Ncell + l*N];
			}
		}

		// calculate change of the phonon (lattice) energy
		// actually shift phonon field entries
		double dEph = 0;

		const double flip = 2*mu/(phonon_params->g[o] + NUM_NEIGHBORS*phonon_params->gp[o]);
		for (l = 0; l < L; l++)
		{
			dEph += 0.5*square(phonon_params->omega[o]) * (-flip) * (-flip + 2*X_ref[l]);
			X[i + l*N] = flip - X_ref[l];
			const double dx = flip - 2*X_ref[l];
			// expX[i + l*N] = exp(-dt*phonon_params->g[o] * X[i + l*N]);
			expX[i + l*N] *= exp(-dt * phonon_params->g[o] * dx);
			int jn;
			for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
			{
				expX[phonon_neighbors[jn + j*(NUM_NEIGHBORS+1)] + o*Ncell + l*N] *= exp(-dt * phonon_params->gp[o] * dx);
			}
			if (phonon_params->J != 0)
			{
				const double xpx = X[phonon_neighbors[1 + j * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				const double xmx = X[phonon_neighbors[2 + j * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				const double xpy = X[phonon_neighbors[3 + j * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				const double xmy = X[phonon_neighbors[4 + j * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				dEph += (phonon_params->J / 2) * dx * (2 * (4 * X_ref[l] - xpx - xmx - xpy - xmy) + 4 * dx);
			}
		}

		// calculate new time step matrices
		#pragma omp parallel sections
		{
			#pragma omp section
			InitPhononTimeStepMatrices(kinetic, stratonovich_params->expVu, s, expX, &tsm_u_new);
			#pragma omp section
			InitPhononTimeStepMatrices(kinetic, stratonovich_params->expVd, s, expX, &tsm_d_new);
		}

		// calculate new Green's functions
		#pragma omp parallel sections
		{
			#pragma omp section
			GreenConstruct(&tsm_u_new, 0, &Gu_new);
			#pragma omp section
			GreenConstruct(&tsm_d_new, 0, &Gd_new);
		}

		// decide whether block update is accepted; note that det(G) = 1/det(M)
		const double p = exp((Gu->logdet + Gd->logdet) - (Gu_new.logdet + Gd_new.logdet) - dt * dEph);
		(*n_flip_total)++;
		if (Random_GetUniform(seed) < p)
		{
			(*n_flip_accept)++;

			// copy new stuff
			CopyGreensFunction(&Gu_new, Gu);
			CopyGreensFunction(&Gd_new, Gd);
			CopyTimeStepMatrices(&tsm_u_new, tsm_u);
			CopyTimeStepMatrices(&tsm_d_new, tsm_d);
		}
		else
		{
			// undo changes
			for (l = 0; l < L; l++)
			{
				   X[i + l*N] =    X_ref[l];
				expX[i + l*N] = expX_ref[l];
				int jn;
				for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
				{
					expX[phonon_neighbors[jn + j*(NUM_NEIGHBORS+1)] + o*Ncell + l*N] = expX_js[l + L*(jn-1)];
				}
			}
		}
	}

	// clean up
	DeleteTimeStepMatrices(&tsm_d_new);
	DeleteTimeStepMatrices(&tsm_u_new);
	DeleteGreensFunction(&Gd_new);
	DeleteGreensFunction(&Gu_new);
	MKL_free(expX_ref);
	MKL_free(X_ref);
	MKL_free(expX_js);
}



//________________________________________________________________________________________________________________________
///
/// \brief Perform a Determinant Quantum Monte Carlo (DQMC) iteration, taking phonons into account
///
/// \param dt                   imaginary-time step
/// \param mu                   chemical potential
/// \param kinetic              matrix exponential of kinetic energy operator
/// \param noHS                 set to true to skip updating the Hubbard-Stratonovich field
/// \param stratonovich_params  precomputed Hubbard-Stratonovich parameters
/// \param phonon_params        phonon field parameters
/// \param nwraps               number of "time slice wraps" before recomputing the Green's function
/// \param seed                 random number "seed" structure, will be updated during function call
/// \param s                    Hubbard-Stratonovich field; will be updated
/// \param X                    phonon field; will be updated
/// \param expX                 entrywise exponential of the phonon field: exp(-dt g X)
/// \param tsm_u                time step matrices for the spin-up Green's function, will be updated
/// \param tsm_d                time step matrices for the spin-down Green's function, will be updated
/// \param Gu                   spin-up   Green's function, must have been computed on input and will be updated
/// \param Gd                   spin-down Green's function, must have been computed on input and will be updated
/// \param neqlt                perform an equal time measurement every 'neqlt' time slices
/// \param meas_data            basic measurement data
/// \param meas_data_phonon     phonon measurement data
///
void DQMCPhononIteration(const double dt, const double mu, const kinetic_t *restrict kinetic, const int noHS, const stratonovich_params_t *restrict stratonovich_params, const phonon_params_t *restrict phonon_params,
	const int nwraps, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX, const int *restrict phonon_neighbors,
	time_step_matrices_t *restrict tsm_u, time_step_matrices_t *restrict tsm_d, greens_func_t *restrict Gu, greens_func_t *restrict Gd,
	const int neqlt, measurement_data_t *restrict meas_data, measurement_data_phonon_t *restrict meas_data_phonon,
	const sim_params_t *restrict params, const char *fnbase)
{
	Profile_Begin("DQMCIter");
	__assume_aligned(s, MEM_DATA_ALIGN);

	// dimension consistency checks
	assert(tsm_u->N == kinetic->Ncell * kinetic->Norb);
	assert(tsm_u->N == tsm_d->N);
	assert(tsm_u->N == Gu->N);
	assert(tsm_u->N == Gd->N);
	assert(tsm_u->L == tsm_d->L);
	const int N = tsm_u->N;
	const int Ncell = kinetic->Ncell;
	const int L = tsm_u->L;

	assert(tsm_u->prodBlen == tsm_d->prodBlen);
	assert(nwraps % tsm_u->prodBlen == 0);  // must be a multiple of 'prodBlen'

	// store Green's functions before recomputing them to estimate error
	#if defined(DEBUG) | defined(_DEBUG)
	greens_func_t Gu_old, Gd_old;
	AllocateGreensFunction(N, &Gu_old);
	AllocateGreensFunction(N, &Gd_old);
	__assume_aligned(Gu_old.mat, MEM_DATA_ALIGN);
	__assume_aligned(Gd_old.mat, MEM_DATA_ALIGN);
	#endif

	// pre-compute 1/dt^2
	const double inv_dt_sq = 1.0 / square(dt);

	// random shuffle of lattice cells and orbitals
	int *orb_cell_order = MKL_malloc(N * sizeof(int), MEM_DATA_ALIGN);
	 __assume_aligned(orb_cell_order, MEM_DATA_ALIGN);

	// iterate over time slices
	int l;
	for (l = 0; l < L; l++)
	{
		Profile_Begin("DQMCIter_Wraps");
		#pragma omp parallel sections
		{
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_u->B[l], tsm_u->invB[l], Gu->mat);
			#pragma omp section
			GreenTimeSliceWrap(N, tsm_d->B[l], tsm_d->invB[l], Gd->mat);
		}
		Profile_End("DQMCIter_Wraps");

		if (!noHS)
		{
			// iterate over lattice sites randomly, updating the Hubbard-Stratonovich field
			Profile_Begin("DQMCIter_HSUpdate");
			Random_Shuffle(seed, N, orb_cell_order);
			int j;
			for (j = 0; j < N; j++)
			{
				const int i = orb_cell_order[j];
				const int o = i / Ncell;    // orbital index
				assert(0 <= i && i < N);
				assert(0 <= o && o < kinetic->Norb);

				// Eq. (13)
				// suggest flipping s_{i,l}
				const double du = 1 + (1 - Gu->mat[i + i*N]) * stratonovich_params->delta[  s[i + l*N]][o];
				const double dd = 1 + (1 - Gd->mat[i + i*N]) * stratonovich_params->delta[1-s[i + l*N]][o];
				if (Random_GetUniform(seed) < fabs(du*dd))
				{
					// Eq. (15)
					#pragma omp parallel sections
					{
						#pragma omp section
						GreenShermanMorrisonUpdate(stratonovich_params->delta[  s[i + l*N]][o], N, i, Gu->mat);
						#pragma omp section
						GreenShermanMorrisonUpdate(stratonovich_params->delta[1-s[i + l*N]][o], N, i, Gd->mat);
					}
					// correspondingly update determinants
					Gu->logdet -= log(fabs(du));
					Gd->logdet -= log(fabs(dd));
					if (du < 0) { Gu->sgndet = -Gu->sgndet; }
					if (dd < 0) { Gd->sgndet = -Gd->sgndet; }

					// actually flip spin
					s[i + l*N] = 1 - s[i + l*N];
				}
			}
			Profile_End("DQMCIter_HSUpdate");
		}

		// next and previous time slices; required for the phonon field with periodic boundary conditions
		const int l_next = (l + 1    ) % L;
		const int l_prev = (l + L - 1) % L;

		// iterate over lattice sites, updating the phonon field
		Profile_Begin("DQMCIter_XUpdate");
		Random_Shuffle(seed, N, orb_cell_order);
		int j;
		for (j = 0; j < phonon_params->n_local_updates; j++)
		{
			const int i = orb_cell_order[j % N];
			const int o = i / Ncell;    // orbital index
			assert(0 <= i && i < N);
			assert(0 <= o && o < kinetic->Norb);

			// skip orbitals without phonon coupling
			if (phonon_params->g[o] == 0 && phonon_params->gp[o] == 0)
			{
				continue;
			}

			// suggest a shift of X_{i,l}
			const double dx = (Random_GetUniform(seed) - 0.5) * phonon_params->local_box_width;

			// Eq. (19) in PRB 87, 235133 (2013)
			// const double delta = expm1(-dt*phonon_params->g[o] * dx);
			// const double du = 1 + (1 - Gu->mat[i + i*N]) * delta;
			// const double dd = 1 + (1 - Gd->mat[i + i*N]) * delta;

			// Woodbury update
			const double delta_i = expm1(-dt*phonon_params->g[o] * dx);
			const double delta_j = expm1(-dt*phonon_params->gp[o] * dx);
			const int n = NUM_NEIGHBORS + 1;
			double *partial_gu_mat = (double *)MKL_malloc(n * n * sizeof(double), MEM_DATA_ALIGN);
			lapack_int *ipiv_u = (lapack_int *)MKL_malloc(n * sizeof(lapack_int), MEM_DATA_ALIGN);
			double *partial_gd_mat = (double *)MKL_malloc(n * n * sizeof(double), MEM_DATA_ALIGN);
			lapack_int *ipiv_d = (lapack_int *)MKL_malloc(n * sizeof(lapack_int), MEM_DATA_ALIGN);
			lapack_int info1, info2;
			double du, dd;
			Profile_Begin("DQMCIter_XU_Determinant");
			// #pragma omp parallel sections
			{
				// #pragma omp section
				{
					info1 = ComputeXUpdateMatrixLU(i, N, Ncell, n, phonon_neighbors, delta_i, delta_j, Gu->mat, partial_gu_mat, ipiv_u);
					du = ComputeXUpdateDeterminant(n, partial_gu_mat, ipiv_u);
				}
				// #pragma omp section
				{
					info2 = ComputeXUpdateMatrixLU(i, N, Ncell, n, phonon_neighbors, delta_i, delta_j, Gd->mat, partial_gd_mat, ipiv_d);
					dd = ComputeXUpdateDeterminant(n, partial_gd_mat, ipiv_d);
				}
			}
			Profile_End("DQMCIter_XU_Determinant");
			if (info1 < 0 || info2 < 0) { duprintf("Intel MKL 'LAPACKE_dgetrf' failed, error code: %i, $i.\n", info1, info2); }

			#if defined(DEBUG) | defined(_DEBUG)
			if (l == nwraps - 1)
			{
				// check the determinant
				double *expX_tmp = (double *)MKL_malloc(L * N * sizeof(double), MEM_DATA_ALIGN);
				memcpy(expX_tmp, expX, L * N * sizeof(double));
				expX_tmp[i + l*N] *= exp(-dt*phonon_params->g[o] * dx);
				int jn;
				for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
				{
					expX_tmp[phonon_neighbors[jn + (i % Ncell)*(NUM_NEIGHBORS+1)] + o*Ncell + l*N] *= exp(-dt*phonon_params->gp[o] * dx);
				}
				time_step_matrices_t tsm_u_tmp;
				time_step_matrices_t tsm_d_tmp;
				AllocateTimeStepMatrices(N, L, params->prodBlen, &tsm_u_tmp);
				AllocateTimeStepMatrices(N, L, params->prodBlen, &tsm_d_tmp);
				#pragma omp parallel sections
				{
					#pragma omp section
					InitPhononTimeStepMatrices(kinetic, stratonovich_params->expVu, s, expX_tmp, &tsm_u_tmp);
					#pragma omp section
					InitPhononTimeStepMatrices(kinetic, stratonovich_params->expVd, s, expX_tmp, &tsm_d_tmp);
				}
				greens_func_t Gu_tmp, Gd_tmp;
				AllocateGreensFunction(N, &Gu_tmp);
				AllocateGreensFunction(N, &Gd_tmp);
				#pragma omp parallel sections
				{
					#pragma omp section
					GreenConstruct(&tsm_u_tmp, (l + 1) % L, &Gu_tmp);
					#pragma omp section
					GreenConstruct(&tsm_d_tmp, (l + 1) % L, &Gd_tmp);
				}
				double *Gu_w = (double *)MKL_malloc(N * N * sizeof(double), MEM_DATA_ALIGN);
				double *Gd_w = (double *)MKL_malloc(N * N * sizeof(double), MEM_DATA_ALIGN);
				double *Gu_s = (double *)MKL_malloc(N * N * sizeof(double), MEM_DATA_ALIGN);
				double *Gd_s = (double *)MKL_malloc(N * N * sizeof(double), MEM_DATA_ALIGN);
				memcpy(Gu_w, Gu->mat, N * N * sizeof(double));
				memcpy(Gd_w, Gd->mat, N * N * sizeof(double));
				memcpy(Gu_s, Gu->mat, N * N * sizeof(double));
				memcpy(Gd_s, Gd->mat, N * N * sizeof(double));
				GreenShermanMorrisonUpdate(delta_i, N, i, Gu_s);
				for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
				{
					GreenShermanMorrisonUpdate(delta_j, N, phonon_neighbors[jn + (i % Ncell)*(NUM_NEIGHBORS+1)], Gu_s);
				}
				GreenShermanMorrisonUpdate(delta_i, N, i, Gd_s);
				for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
				{
					GreenShermanMorrisonUpdate(delta_j, N, phonon_neighbors[jn + (i % Ncell)*(NUM_NEIGHBORS+1)], Gd_s);
				}
				GreenXWoodburyUpdate(i, N, Ncell, n, phonon_neighbors, delta_i, delta_j, Gu_w, partial_gu_mat, ipiv_u);
				GreenXWoodburyUpdate(i, N, Ncell, n, phonon_neighbors, delta_i, delta_j, Gd_w, partial_gd_mat, ipiv_d);

				char path[1024];
				sprintf(path, "%s_Gu_w_debug.dat", fnbase); WriteData(path, Gu_w, sizeof(double), N*N, true);
				sprintf(path, "%s_Gd_w_debug.dat", fnbase); WriteData(path, Gd_w, sizeof(double), N*N, true);
				sprintf(path, "%s_Gu_s_debug.dat", fnbase); WriteData(path, Gu_s, sizeof(double), N*N, true);
				sprintf(path, "%s_Gd_s_debug.dat", fnbase); WriteData(path, Gd_s, sizeof(double), N*N, true);

				// PrintMatrix(N, N, "Gu_tmp", Gu_tmp.mat);
				// PrintMatrix(N, N, "Gu_w  ", Gu_w);
				// PrintMatrix(N, N, "Gu_s  ", Gu_s);
				// PrintMatrix(N, N, "Gd_tmp", Gd_tmp.mat);
				// PrintMatrix(N, N, "Gd_w  ", Gd_w);
				// PrintMatrix(N, N, "Gd_s  ", Gd_s);
				double gu_err = UniformDistance(N*N, Gu_tmp.mat, Gu_w);
				double gd_err = UniformDistance(N*N, Gd_tmp.mat, Gd_w);
				if (gu_err > 1e-12 || gd_err > 1e-12) {
					duprintf("Woodbury warning: Largest entrywise distance between updated Green's function and from scratch is %g (up) and %g (down).\n", gu_err, gd_err);
				}
				// gu_err = UniformDistance(N*N, Gu_s, Gu_w);
				// gd_err = UniformDistance(N*N, Gd_s, Gd_w);
				// if (gu_err > 1e-16 || gd_err > 1e-16) {
				// 	duprintf("Woodbury warning: Largest entrywise distance between woodbury and sherman morrison is %g (up) and %g (down).\n", gu_err, gd_err);
				// }


 				double du_tmp = exp(Gu->logdet - Gu_tmp.logdet);
    			double dd_tmp = exp(Gd->logdet - Gd_tmp.logdet);
				double det_tmp = exp((Gu->logdet + Gd->logdet) - (Gu_tmp.logdet + Gd_tmp.logdet));
 				duprintf("du_gp =%.16f\n", du);
                duprintf("du_tmp=%.16f\n", du_tmp);
                duprintf("dd_gp =%.16f\n", dd);
                duprintf("dd_tmp=%.16f\n", dd_tmp);
				duprintf("det     = %.16f\n", du * dd);
				duprintf("det_tmp = %.16f\n", det_tmp);
				duprintf("Gu.logdet %g, gu.logdet %g, gu_tmp.logdet %g, gd_tmp.logdet %g\n\n", Gu->logdet, Gd->logdet, Gu_tmp.logdet, Gd_tmp.logdet);
				if (du * dd - det_tmp > 1e-12) {
					duprintf("Woodbury warning: Difference between determinant from woodbury and from scratch is %g.\n", du * dd - det_tmp);
				}


				// char path[1024];
				// sprintf(path, "%s_gu_logdet.dat", fnbase); WriteData(path, &Gu->logdet, sizeof(double), 1, true);
				// sprintf(path, "%s_gd_logdet.dat", fnbase); WriteData(path, &Gd->logdet, sizeof(double), 1, true);
				// sprintf(path, "%s_gu_tmp_logdet.dat", fnbase); WriteData(path, &Gu_tmp.logdet, sizeof(double), 1, true);
				// sprintf(path, "%s_gd_tmp_logdet.dat", fnbase); WriteData(path, &Gd_tmp.logdet, sizeof(double), 1, true);
				// sprintf(path, "%s_det_tmp_logdet.dat", fnbase); WriteData(path, &det_tmp, sizeof(double), 1, true);
				DeleteGreensFunction(&Gu_tmp);
				DeleteGreensFunction(&Gd_tmp);
				DeleteTimeStepMatrices(&tsm_d_tmp);
				DeleteTimeStepMatrices(&tsm_u_tmp);
				MKL_free(Gu_w);
				MKL_free(Gd_w);
				MKL_free(Gu_s);
				MKL_free(Gd_s);
				MKL_free(expX_tmp);
			}
			#endif


			// change of the phonon (lattice) energy
			double dEph = dx * (0.5*square(phonon_params->omega[o])*(dx + 2*X[i + l*N]) + inv_dt_sq*(dx - (X[i + l_next*N] - 2*X[i + l*N] + X[i + l_prev*N])));
			if (phonon_params->J != 0)
			{
				const double xpx = X[phonon_neighbors[1 + (i % Ncell) * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				const double xmx = X[phonon_neighbors[2 + (i % Ncell) * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				const double xpy = X[phonon_neighbors[3 + (i % Ncell) * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				const double xmy = X[phonon_neighbors[4 + (i % Ncell) * (NUM_NEIGHBORS + 1)] + o * Ncell + l*N];
				dEph += (phonon_params->J / 2) * dx * (2 * (4 * X[i + l*N] - xpx - xmx - xpy - xmy) + 4 * dx);
			}

			meas_data_phonon->n_local_total++;
			if (Random_GetUniform(seed) < fabs(du*dd) * exp(-dt * dEph))
			{
				meas_data_phonon->n_local_accept++;
				// Eq. (15)
				int jn;
				// Profile_Begin("DQMCIter_XUpdate_ShermanMorrisonUpdate");
				// GreenShermanMorrisonUpdate(delta_i, N, i, Gu->mat);
				// for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
				// {
				// 	GreenShermanMorrisonUpdate(delta_j, N, phonon_neighbors[jn + (i % Ncell)*(NUM_NEIGHBORS+1)], Gu->mat);
				// }
				// GreenShermanMorrisonUpdate(delta_i, N, i, Gd->mat);
				// for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
				// {
				// 	GreenShermanMorrisonUpdate(delta_j, N, phonon_neighbors[jn + (i % Ncell)*(NUM_NEIGHBORS+1)], Gd->mat);
				// }
				// Profile_End("DQMCIter_XUpdate_ShermanMorrisonUpdate");
				Profile_Begin("DQMCIter_XU_Woodbury");
				#pragma omp parallel sections
				{
					#pragma omp section
					GreenXWoodburyUpdate(i, N, Ncell, n, phonon_neighbors, delta_i, delta_j, Gu->mat, partial_gu_mat, ipiv_u);
					#pragma omp section
					GreenXWoodburyUpdate(i, N, Ncell, n, phonon_neighbors, delta_i, delta_j, Gd->mat, partial_gd_mat, ipiv_d);
				}
				Profile_End("DQMCIter_XU_Woodbury");
				// correspondingly update determinants
				Gu->logdet -= log(fabs(du));
				Gd->logdet -= log(fabs(dd));
				if (du < 0) { Gu->sgndet = -Gu->sgndet; }
				if (dd < 0) { Gd->sgndet = -Gd->sgndet; }

				// actually update the phonon field
				   X[i + l*N] += dx;
				expX[i + l*N] *= exp(-dt*phonon_params->g[o] * dx);

                for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
 				{
 					expX[phonon_neighbors[jn + (i % Ncell)*(NUM_NEIGHBORS+1)] + o*Ncell + l*N] *= exp(-dt*phonon_params->gp[o] * dx);
 				}
			}
			MKL_free(partial_gu_mat);
			MKL_free(partial_gd_mat);
			MKL_free(ipiv_u);
			MKL_free(ipiv_d);
		}
		Profile_End("DQMCIter_XUpdate");

		// re-compute corresponding B matrices
		Profile_Begin("DQMCIter_Brecomp");
		#pragma omp parallel sections
		{
			#pragma omp section
			UpdatePhononTimeStepMatrices(kinetic, stratonovich_params->expVu, s, expX, l, tsm_u);
			#pragma omp section
			UpdatePhononTimeStepMatrices(kinetic, stratonovich_params->expVd, s, expX, l, tsm_d);
		}
		Profile_End("DQMCIter_Brecomp");

		// recompute Green's function after several time slice "wraps"
		if ((l + 1) % nwraps == 0)
		{
			Profile_Begin("DQMCIter_Grecomp");
			// store current Green's function matrices to compare with newly constructed ones
			#if defined(DEBUG) | defined(_DEBUG)
			CopyGreensFunction(Gu, &Gu_old);
			CopyGreensFunction(Gd, &Gd_old);
			#endif

			#pragma omp parallel sections
			{
				#pragma omp section
				GreenConstruct(tsm_u, (l + 1) % L, Gu);
				#pragma omp section
				GreenConstruct(tsm_d, (l + 1) % L, Gd);
			}

			#if defined(DEBUG) | defined(_DEBUG)
			// deviation of matrix entries
			double err_u = UniformDistance(N*N, Gu_old.mat, Gu->mat);
			double err_d = UniformDistance(N*N, Gd_old.mat, Gd->mat);
			if (err_u > 1e-8*N || err_d > 1e-8*N) {
				duprintf("Warning: after calling 'GreenConstruct()', largest entrywise distance between previous and new Green's functions is %g (up) and %g (down).\n", err_u, err_d);
			}
			// deviation of matrix determinants
			err_u = fabs(Gu_old.logdet - Gu->logdet);
			err_d = fabs(Gd_old.logdet - Gd->logdet);
			if (err_u > 1e-8 || err_d > 1e-8) {
				duprintf("Warning: after calling 'GreenConstruct()', largest distance between logarithm of previous and new Green's function determinants is %g (up) and %g (down).\n", err_u, err_d);
			}
			if (Gu_old.sgndet != Gu->sgndet || Gd_old.sgndet != Gd->sgndet) {
				duprintf("Warning: after calling 'GreenConstruct()', determinant sign has changed.\n");
			}
			#endif
			Profile_End("DQMCIter_Grecomp");
		}

		if (neqlt > 0 && (l + 1) % neqlt == 0)
		{
			// accumulate equal time "measurement" data
			Profile_Begin("DQMCIter_AccumulateEqMeas");
			AccumulateMeasurement(Gu, Gd, meas_data);
			Profile_End("DQMCIter_AccumulateEqMeas");
			} else if (neqlt > 0 && l == 0) {
			// accumulate phonon data
			Profile_Begin("DQMCIter_AccumulatePhonon");
			AccumulatePhononData(Gu, Gd, l, X, dt, phonon_params->omega, meas_data_phonon);
			Profile_End("DQMCIter_AccumulatePhonon");
		}
	}

    /*
    Profile_Begin("DI_4");
	if (phonon_params->track_phonon_ite)
	{
		int Norb = kinetic->Norb;
		double *X0 = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
		double sign = (double)(Gu->sgndet * Gd->sgndet);
		int o;
		for (l = 0; l < L; l++)
		{
			for (o = 0; o < Norb; o++)
			{
				X0[o] += (sign/L) * X[o*Ncell + l*N];
			}
		}
		char path[1024];
		sprintf(path, "%s_phonon_iteration2_X0.dat",   fnbase); WriteData(path, X0,    sizeof(double), Norb, true);
		sprintf(path, "%s_phonon_iteration2_sign.dat", fnbase); WriteData(path, &sign, sizeof(double), 1,    true);
		MKL_free(X0);
	}
	Profile_End("DI_4");
	*/

    Profile_Begin("DI_4");
	if (phonon_params->track_phonon_ite)
	{
		int Norb = kinetic->Norb;
		double *PE = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
		double *KE = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
		double *X0 = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
		double *X_avg = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
		double *Xl0_avg = (double *)MKL_calloc(Norb, sizeof(double), MEM_DATA_ALIGN);
		double sign = (double)(Gu->sgndet * Gd->sgndet);
		double signfac = sign / L / Ncell;
		int o;
		for (l = 0; l < L; l++)
		{
			const int lplus = (l + 1) % L;
			for (o = 0; o < Norb; o++)
			{
				int i;
				for (i = 0; i < Ncell; i++)
				{
					PE[o] += 0.5*square(phonon_params->omega[o])*signfac*square(X[i + o*Ncell + l*N]);
					KE[o] += 0.5/(dt*dt)*signfac*square(X[i + o*Ncell + lplus*N] - X[i + o*Ncell + l*N]);
					X_avg[o] += signfac * X[i + o*Ncell + l*N];
					if (l == 0)
					{
						Xl0_avg[o] += (sign/Ncell) * X[i + o*Ncell];
					}
				}
				X0[o] += (sign/L) * X[o*Ncell + l*N];
			}
		}
		char path[1024];
		sprintf(path, "%s_phonon_iteration2_KE.dat",    fnbase); WriteData(path, KE,    sizeof(double), Norb, true);
		sprintf(path, "%s_phonon_iteration2_PE.dat",    fnbase); WriteData(path, PE,    sizeof(double), Norb, true);
		sprintf(path, "%s_phonon_iteration2_X0.dat",    fnbase); WriteData(path, X0,    sizeof(double), Norb, true);
		sprintf(path, "%s_phonon_iteration2_X_avg.dat", fnbase); WriteData(path, X_avg, sizeof(double), Norb, true);
		sprintf(path, "%s_phonon_iteration2_Xl0_avg.dat", fnbase); WriteData(path, Xl0_avg, sizeof(double), Norb, true);
		sprintf(path, "%s_phonon_iteration2_sign.dat",  fnbase); WriteData(path, &sign, sizeof(double), 1,    true);
		MKL_free(X0);
		MKL_free(X_avg);
		MKL_free(Xl0_avg);
		MKL_free(PE);
		MKL_free(KE);
	}
	Profile_End("DI_4");

	// perform block updates
	Profile_Begin("DQMCIter_PhononBlock");
	PhononBlockUpdates(dt, kinetic, stratonovich_params, phonon_params, seed, s, X, expX, phonon_neighbors, tsm_u, tsm_d, Gu, Gd,
	                   &meas_data_phonon->n_block_accept, &meas_data_phonon->n_block_total);
	Profile_End("DQMCIter_PhononBlock");

	Profile_Begin("DQMCIter_PhononFlip");
	PhononFlipUpdates(dt, mu, kinetic, stratonovich_params, phonon_params, seed, s, X, expX, phonon_neighbors, tsm_u, tsm_d, Gu, Gd,
	                  &meas_data_phonon->n_flip_accept, &meas_data_phonon->n_flip_total);
	Profile_End("DQMCIter_PhononFlip");

	// clean up
	MKL_free(orb_cell_order);
	#if defined(DEBUG) | defined(_DEBUG)
	DeleteGreensFunction(&Gd_old);
	DeleteGreensFunction(&Gu_old);
	#endif
	Profile_End("DQMCIter");
}


//________________________________________________________________________________________________________________________
///
/// \brief Perform a Determinant Quantum Monte Carlo (DQMC) simulation
///
/// \param params               simulation parameters
/// \param meas_data            measurement data structure for accumulating measurements
/// \param meas_data_uneqlt     unequal time measurement data structure
/// \param meas_data_phonon     phonon measurement data structure
/// \param iteration            pointer to iteration counter (i.e. number of iterations completed from previous runs)
/// \param seed                 random number generator seed
/// \param s                    Hubbard-Stratonovich field; will be updated
/// \param X                    phonon field; only accessed if params->use_phonons is true
/// \param expX                 entrywise exponential of the phonon field: exp(-dt g X)
///
void DQMCSimulation(const sim_params_t *restrict params,
	measurement_data_t *restrict meas_data, measurement_data_unequal_time_t *restrict meas_data_uneqlt, measurement_data_phonon_t *restrict meas_data_phonon,
	int *restrict iteration, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX, const int *restrict phonon_neighbors, const char *fnbase)
{
	const int Norb  = params->Norb;
	const int Ncell = params->Nx * params->Ny;
	const int N     = Norb * Ncell;

	// get the time (in ticks) for automatic stopping and checkpointing
	const uint64_t t_end = GetTicks() + params->max_time * GetTickRes();

	// calculate matrix exponential of the kinetic nearest neighbor hopping matrix
	kinetic_t kinetic;
	RectangularKineticExponential(params, &kinetic);

	// pre-calculate some stuff related to the Hubbard-Stratonovich field, for every orbital
	int noHS = 1; // flag to disable H-S updates if all U == 0
	int o;
	for (o = 0; o < Norb; o++)
	{
		if (params->U[o] != 0)
		{
			noHS = 0;
			break;
		}
	}
	stratonovich_params_t stratonovich_params;
	FillStratonovichParameters(Norb, params->U, params->dt, &stratonovich_params);

	// time step matrices
	time_step_matrices_t tsm_u;
	time_step_matrices_t tsm_d;
	AllocateTimeStepMatrices(N, params->L, params->prodBlen, &tsm_u);
	AllocateTimeStepMatrices(N, params->L, params->prodBlen, &tsm_d);
	if (params->use_phonons)
	{
		#pragma omp parallel sections
		{
			#pragma omp section
			InitPhononTimeStepMatrices(&kinetic, stratonovich_params.expVu, s, expX, &tsm_u);
			#pragma omp section
			InitPhononTimeStepMatrices(&kinetic, stratonovich_params.expVd, s, expX, &tsm_d);
		}
	}
	else
	{
		#pragma omp parallel sections
		{
			#pragma omp section
			InitTimeStepMatrices(&kinetic, stratonovich_params.expVu, s, &tsm_u);
			#pragma omp section
			InitTimeStepMatrices(&kinetic, stratonovich_params.expVd, s, &tsm_d);
		}
	}

	// allocate Green's functions
	greens_func_t Gu, Gd;
	AllocateGreensFunction(N, &Gu);
	AllocateGreensFunction(N, &Gd);

	// construct initial Green's functions
	#pragma omp parallel sections
	{
		#pragma omp section
		GreenConstruct(&tsm_u, 0, &Gu);
		#pragma omp section
		GreenConstruct(&tsm_d, 0, &Gd);
	}

	// register a signal handler to stop simulation on SIGINT
	StopOnSIGINT();

	// first call to UpdateProgress to indicate completion of initialization phase
	UpdateProgress();

	// perform dqmc iterations
	duprintf("Starting DQMC iterations...\n");
	for (; *iteration < params->nequil + params->nsampl; (*iteration)++)
	{
		if (params->max_time > 0 && GetTicks() >= t_end)
		{
			duprintf("Reached time limit of %d seconds.\n", params->max_time);
			stopped = 2;
		}

		if (stopped == 1 || stopped == 2) // either the above happened or SIGINT was received
		{
			duprintf("Stopping DQMC iterations early.\n");
			break;
		}

		// set neqlt to 0 in equilibration stage so that no measurements are made
		const int neqlt = (*iteration >= params->nequil) ? params->neqlt : 0;

		if (params->use_phonons)
		{
			DQMCPhononIteration(params->dt, params->mu, &kinetic, noHS, &stratonovich_params, &params->phonon_params, params->nwraps, seed, s, X, expX, phonon_neighbors, &tsm_u, &tsm_d, &Gu, &Gd, neqlt, meas_data, meas_data_phonon, params, fnbase);
		}
		else
		{
			DQMCIteration(&kinetic, &stratonovich_params, params->nwraps, seed, s, &tsm_u, &tsm_d, &Gu, &Gd, neqlt, meas_data);
		}

		// accumulate unequal time "measurement" data
		if (*iteration >= params->nequil && params->nuneqlt > 0 && (*iteration % params->nuneqlt) == 0)
		{
			Profile_Begin("DQMCSim_AccumulateUneqMeas");
			AccumulateUnequalTimeMeasurement((double)(Gu.sgndet * Gd.sgndet), &tsm_u, &tsm_d, meas_data_uneqlt);
			Profile_End("DQMCSim_AccumulateUneqMeas");
		}

		UpdateProgress();
	}

	// clean up
	DeleteGreensFunction(&Gd);
	DeleteGreensFunction(&Gu);
	DeleteTimeStepMatrices(&tsm_d);
	DeleteTimeStepMatrices(&tsm_u);
	DeleteStratonovichParameters(&stratonovich_params);
	DeleteKineticExponential(&kinetic);
}
