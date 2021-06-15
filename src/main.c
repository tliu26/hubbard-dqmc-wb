#include "monte_carlo.h"
#include "kinetic.h"
#include "sim_params.h"
#include "random.h"
#include "stratonovich.h"
#include "profiler.h"
#include "checkpoint.h"
#include "progress.h"
#include "util.h"
#include "dupio.h"
#include "phonon_gp.h"
#include <mkl.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h>

// for sleep function and creating directories
#ifdef _WIN32
#include <windows.h>
#include <direct.h>

static int makedir(const char *path)
{
	return _mkdir(path);
}

#else

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static int makedir(const char *path)
{
	return mkdir(path, 0755);
}

#endif

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif


int main(int argc, char *argv[])
{
	int status;

	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	// if compiled with OpenMP, use only 2 threads max for efficiency
	omp_set_num_threads(2);

	// initialize profiling
	Profile_Start();

	// make sure there's an input file
	if (argc < 2)
	{
		duprintf("No input file specified!\n");
		return -1;
	}

	// zero simulation parameters and use initial time as random seed (will be overwritten if explicitly specified in parameter file)
	sim_params_t params = { 0 };
	params.itime = GetTicks();

	// read parameters from input file
	duprintf("Reading simulation parameters from file '%s'...\n", argv[1]);
	status = ParseParameterFile(argv[1], &params);  // this also allocates memory for the arrays in params
	if (status < 0)
	{
		duprintf("Error parsing parameter file, exiting...\n");
		return -1;
	}

	// perform admissibility checks of the simulation parameters
	status = ValidateSimulationParameters(&params);
	if (status < 0) {
		duprintf("Parameters could not be validated, exiting...\n");
		return -2;
	}

	// create output directory or read from 2nd command line argument
	char path[1024];
	if (argc < 3)
	{
		duprintf("No output directory specified.\n");
		duprintf("Creating default directory based on parameters and initial seed...\n");
		makedir("output");
		sprintf(path, "output/N%ix%i_beta%g_mu%g_sim_%llu", params.Nx, params.Ny, params.L * params.dt, params.mu, params.itime);
	}
	else
	{
		strcpy(path, argv[2]);
	}
	makedir(path);
	duprintf("Using output directory '%s'.\n\n", path);

	// base output file name
	char fnbase[1024];
	sprintf(fnbase, "%s/sim_%llu", path, params.itime);

	// open simulation log file for writing
	sprintf(path, "%s_simulation.log", fnbase);
	fd_log = fopen(path, "a");
	if (fd_log == NULL)
	{
		duprintf("Cannot open log file '%s', exiting...\n", path);
		return -3;
	}

	// allocate and initialize equal time measurement data structure
	measurement_data_t meas_data;
	AllocateMeasurementData(params.Norb, params.Nx, params.Ny, params.pbc_shift, &meas_data);

	// allocate and initialize unequal time measurement data structure
	measurement_data_unequal_time_t meas_data_uneqlt;
	if (params.nuneqlt > 0)
	{
		status = AllocateUnequalTimeMeasurementData(params.Norb, params.Nx, params.Ny, params.pbc_shift, params.L, params.L / params.prodBlen, &meas_data_uneqlt);
		if (status < 0) {
			duprintf("Could not allocate unequal time measurement data structure (probably out of memory), exiting...\n");
			return -4;
		}
	}

	// allocate and initialize phonon measurement data structure
	measurement_data_phonon_t meas_data_phonon;
	if (params.use_phonons)
	{
		AllocatePhononData(params.Norb, params.Nx, params.Ny, params.pbc_shift, params.L,
			params.nsampl, &meas_data_phonon);
	}

	// dimensions
	const int Ncell = params.Nx * params.Ny;
	const int N     = params.Norb * Ncell;
	const int LxN   = params.L * N;

	// state variables for DQMC simulation; these are also the variables stored in a checkpoint
	// TODO: include phonons in checkpoint
	int iteration = 0;
	randseed_t seed;
	spin_field_t *s = (spin_field_t *)MKL_malloc(LxN * sizeof(spin_field_t), MEM_DATA_ALIGN);

	// initialize phonon field
	double *X = NULL, *expX = NULL;
	int *phonon_neighbors;
	phonon_neighbors = (int *)MKL_malloc(Ncell * (NUM_NEIGHBORS + 1) * sizeof(int), MEM_DATA_ALIGN);
	ConstructPhononNeigborsMap(params.Nx, params.Ny, phonon_neighbors);
	if (params.use_phonons)
	{
		X    = (double *)MKL_malloc(LxN * sizeof(double), MEM_DATA_ALIGN);
		expX = (double *)MKL_malloc(LxN * sizeof(double), MEM_DATA_ALIGN);
	}

	// check for previous data in output directory
	if (SearchCheckpoint(fnbase, params.use_phonons) != 0) // no previous data found, start from scratch
	{
		// random generator seed; multiplicative constant from Pierre L'Ecuyer's paper
		Random_SeedInit(1865811235122147685LL * params.itime, &seed);

		// random initial Hubbard-Stratonovich field
		int i;
		for (i = 0; i < LxN; i++)
		{
			s[i] = (Random_GetUniform(&seed) < 0.5 ? 0 : 1);
		}

		if (params.use_phonons)
		{
			int l;
			for (l = 0; l < params.L; l++)
			{
				for (i = 0; i < N; i++)
				{
					const int o = i / Ncell;    // orbital index
					if (params.phonon_params.g[o] == 0 && params.phonon_params.gp[o] == 0) // set X to 0 if coupling is zero
					{
						X[i + l*N] = 0;
					}
					else
					{
						X[i + l*N] = (Random_GetUniform(&seed) - 0.5) * params.phonon_params.local_box_width;
					}
					expX[i + l*N] = exp(-params.dt*params.phonon_params.g[o] * X[i + l*N]);
				}
				// phonon g'
				for (i = 0; i < N; i++)
				{
					const int o = i / Ncell;
					const int j = i % Ncell;
					int jn;
					for (jn = 1; jn <= NUM_NEIGHBORS; jn++)
					{
						expX[phonon_neighbors[jn + j*(NUM_NEIGHBORS+1)] + o*Ncell + l*N] *= exp(-params.dt * params.phonon_params.gp[o] * X[i + l*N]);
					}
				}
			}
		}

		// print a header and the parameters
		duprintf("Hubbard model DQMC\n  git commit id %s\n  compiled on %s\n", VERSION, __DATE__);
		duprintf("_______________________________________________________________________________\n");
		PrintSimulationParameters(&params);

		duprintf("Starting DQMC simulation...\n");
	}
	else // found a previous checkpoint, so load the previous state
	{
		status = LoadCheckpoint(fnbase, &iteration, &seed, s, X, expX, LxN, params.use_phonons);
		if (status < 0)
		{
			duprintf("Failed to load previous checkpoint, exiting...\n");
			return -5;
		}
		duprintf("Loaded checkpoint.\n");

		if (iteration == params.nequil + params.nsampl)
		{
			duprintf("All iterations already completed in previous checkpoint, exiting...\n");
			return -6;
		}

		// load previous measurement data
		LoadMeasurementData(fnbase, &meas_data);
		if (params.nuneqlt > 0)
		{
			LoadUnequalTimeMeasurementData(fnbase, &meas_data_uneqlt);
		}
		if (params.use_phonons)
		{
			LoadPhononData(fnbase, &meas_data_phonon);
		}

		duprintf("Resuming DQMC simulation at iteration %d.\n", iteration);
	}

	// enable progress tracking (progress of simulation is shown whenever a SIGUSR1 signal is received)
	InitProgressTracking(&iteration, params.nequil, params.nsampl);

	// start timer
	const clock_t t_start = clock();

	// perform simulation
	DQMCSimulation(&params, &meas_data, &meas_data_uneqlt, &meas_data_phonon, &iteration, &seed, s, X, expX, phonon_neighbors, fnbase);

	// stop timer
	const clock_t t_end = clock();
	double cpu_time = (double)(t_end - t_start) / CLOCKS_PER_SEC;
	duprintf("%d iterations completed, CPU time: %g\n", iteration, cpu_time);

	// at end of simulation, normalize measurement data and show summary of results
	if (iteration == params.nequil + params.nsampl)
	{
		duprintf("All iterations completed.\n");
		NormalizeMeasurementData(&meas_data);
		if (params.nuneqlt > 0)
		{
			NormalizeUnequalTimeMeasurementData(&meas_data_uneqlt);
		}
		// show some simulation results
		PrintMeasurementDataSummary(&meas_data);

		if (params.use_phonons)
		{
			NormalizePhononData(&meas_data_phonon);
			PrintPhononData(&meas_data_phonon);
		}
	}

	// save checkpoint for next run. even if simulation is finished, saving the HS field
	// is helpful if someone wants to extend the simulation further.
	duprintf("Saving checkpoint to disk...");
	status = SaveCheckpoint(fnbase, &iteration, &seed, s, X, expX, LxN, params.use_phonons);
	if (status < 0)
	{
		duprintf("\nFailed to save checkpoint, exiting...\n");
		return -7;
	}
	duprintf(" done.\n");

	// save simulation results as binary data to disk
	duprintf("Saving simulation results to disk...");
	SaveMeasurementData(fnbase, &meas_data);
	if (params.nuneqlt > 0)
	{
		SaveUnequalTimeMeasurementData(fnbase, &meas_data_uneqlt);
	}
	if (params.use_phonons)
	{
		SavePhononData(fnbase, &meas_data_phonon);
	}
	duprintf(" done.\n");

	// clean up
	Profile_Stop();
	fclose(fd_log);
	if (params.use_phonons)
	{
		MKL_free(expX);
		MKL_free(X);
		DeletePhononData(&meas_data_phonon);
	}
	MKL_free(phonon_neighbors);
	MKL_free(s);
	if (params.nuneqlt > 0)
	{
		DeleteUnequalTimeMeasurementData(&meas_data_uneqlt);
	}
	DeleteMeasurementData(&meas_data);
	DeleteSimulationParameters(&params);
	MKL_Free_Buffers();

	// check for MKL memory leaks
	#ifdef _DEBUG
	int nbuffers;
	MKL_INT64 nbytes_alloc;
	nbytes_alloc = MKL_Mem_Stat(&nbuffers);
	if (nbytes_alloc > 0)
	{
		printf("\nMKL memory leak detected! MKL still uses %lld bytes in %d buffer(s).\n", nbytes_alloc, nbuffers);
	}
	else
	{
		printf("\nMKL memory leak check appears to be fine.\n");
	}
	#endif

	if (stopped == 2)
	{
		return 0;
	}
	return stopped; // 0 if simulation ran to completion, 1 if stopped by SIGINT or reaching the max run time limit
}
