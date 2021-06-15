#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "random.h"
#include "stratonovich.h"
#include <stdbool.h>


extern int stopped;


int StopOnSIGINT(void);


int InitCheckpointing(int nequil, int nsampl);


int SearchCheckpoint(const char *fnbase, bool use_phonons);


int LoadCheckpoint(const char *restrict fnbase, int *restrict iteration, randseed_t *restrict seed, spin_field_t *restrict s, double *restrict X, double *restrict expX, const int LxN, bool use_phonons);

int SaveCheckpoint(const char *restrict fnbase, const int *restrict iteration, const randseed_t *restrict seed, const spin_field_t *restrict s, const double *X, const double *expX, const int LxN, bool use_phonons);



#endif
