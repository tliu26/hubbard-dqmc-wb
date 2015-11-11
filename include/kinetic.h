#ifndef KINETIC_H
#define KINETIC_H

#include "param_parser.h"


//________________________________________________________________________________________________________________________
///
/// \brief Matrix exponential of the kinetic energy operator
///
typedef struct
{
	double *expK;			//!< exp(-dt K)
	double *inv_expK;		//!< exp( dt K)
	int Norb;				//!< number of orbitals per unit cell
	int Ncell;				//!< total number of unit cells
	int N;					//!< total number of orbitals
}
kinetic_t;


void RectangularKineticExponential(const sim_params_t *restrict params, kinetic_t *restrict kinetic);

void DeleteKineticExponential(kinetic_t *restrict kinetic);



#endif
