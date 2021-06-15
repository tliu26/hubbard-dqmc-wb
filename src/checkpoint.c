#include "checkpoint.h"
#include <dupio.h>
#include <util.h>
#include <signal.h>
#include <sys/stat.h>
#include <stdint.h>


int stopped;

static void StopSimulation(int signum)
{
	duprintf("Signal %d received.\n", signum);
	stopped = 1;
}

int StopOnSIGINT(void)
{
	return (signal(SIGINT, StopSimulation) == SIG_ERR) ? -1 : 0;
	//return sigaction(SIGINT, &(struct sigaction){.sa_handler = StopSimulation}, NULL);
}


int SearchCheckpoint(const char *fnbase, bool use_phonons)
{
	char path[1024];
	struct stat st;

	sprintf(path, "%s_iteration.dat", fnbase);
	if (stat(path, &st) != 0) { return -1; }

	sprintf(path, "%s_randseed.dat", fnbase);
	if (stat(path, &st) != 0) { return -1; }

	sprintf(path, "%s_hsfield.dat", fnbase);
	if (stat(path, &st) != 0) { return -1; }

	if (use_phonons)
	{
		sprintf(path, "%s_X.dat", fnbase);
		if (stat(path, &st) != 0) { return -1; }

		sprintf(path, "%s_expX.dat", fnbase);
		if (stat(path, &st) != 0) { return -1; }
	}

	return 0;
}


int LoadCheckpoint(const char *restrict fnbase, int *restrict iteration, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX, const int LxN, bool use_phonons)
{
	int status = 0;
	char path[1024];

	sprintf(path, "%s_iteration.dat", fnbase);
	status = ReadData(path, iteration, sizeof(int), 1);
	if (status < 0) { return status; }

	sprintf(path, "%s_randseed.dat", fnbase);
	status = ReadData(path, seed->s, sizeof(uint64_t), 16);
	if (status < 0) { return status; }
	seed->p = 0;

	sprintf(path, "%s_hsfield.dat", fnbase);
	status = ReadData(path, s, sizeof(spin_field_t), LxN);
	if (status < 0) { return status; }

	if (use_phonons)
	{
		sprintf(path, "%s_X.dat", fnbase);
		status = ReadData(path, X, sizeof(double), LxN);
		if (status < 0) { return status; }

		sprintf(path, "%s_expX.dat", fnbase);
		status = ReadData(path, expX, sizeof(double), LxN);
		if (status < 0) { return status; }
	}

	return 0;
}


int SaveCheckpoint(const char *restrict fnbase, const int *restrict iteration, const randseed_t *restrict seed, const spin_field_t *restrict s, const double *X, const double *expX, const int LxN, bool use_phonons)
{
	int status = 0;
	char path[1024];

	sprintf(path, "%s_iteration.dat", fnbase);
	status = WriteData(path, iteration, sizeof(int), 1, false);
	if (status < 0) { return status; }

	// instead of storing seed->p, just save a permutation of seed->s where the 0th element is the p'th element of seed->s.
	sprintf(path, "%s_randseed.dat", fnbase);
	uint64_t a[16];
	int i;
	for (i = 0; i < 16; i++)
	{
		a[i] = seed->s[(i + seed->p + 16) % 16];
	}
	status = WriteData(path, a, sizeof(uint64_t), 16, false);
	if (status < 0) { return status; }

	sprintf(path, "%s_hsfield.dat", fnbase);
	status = WriteData(path, s, sizeof(spin_field_t), LxN, false);
	if (status < 0) { return status; }

	if (use_phonons)
	{
		sprintf(path, "%s_X.dat", fnbase);
		status = WriteData(path, X, sizeof(double), LxN, false);
		if (status < 0) { return status; }

		sprintf(path, "%s_expX.dat", fnbase);
		status = WriteData(path, expX, sizeof(double), LxN, false);
		if (status < 0) { return status; }
	}

	return 0;
}
