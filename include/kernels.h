///@endcond
//##############################################################################
/**
 *  @file    kernels.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief GPU kernel definitions
 *
 *  @section DESCRIPTION
 *  Kernel definitions for all CUDA-enabled routines for solving GPE.
 */
//##############################################################################

#ifndef KERNELS_H
#define KERNELS_H
#include<stdio.h>

/**
* @brief	subtraction operation for 2 double2 values
* @ingroup	gpu
*/
__device__ double2 subtract(double2 a, double2 b);

/**
* @brief	addition operation for 2 double2 values
* @ingroup	gpu
*/
__device__ double2 add(double2 a, double2 b);

/**
* @brief	power operation for a double2
* @param        base number
* @param        power
* @ingroup	gpu
*/
__device__ double2 pow(double2 a, int b);

/**
* @brief	multiplication operation for double2 and double values
* @ingroup	gpu
*/
__device__ double2 mult(double2 a, double b);

/**
* @brief	multiplication operation for 2 double2 values
* @ingroup	gpu
*/
__device__ double2 mult(double2 a, double2 b);


/**
* @brief	Indexing of threads on grid
* @ingroup	gpu
*/
__device__ unsigned int getGid3d3d();

//##############################################################################
/**
* Helper functions for complex numbers
*/
//##############################################################################

/**
* @brief        Sums double2* and double2* energies
* @ingroup      gpu
* @param        Array 1
* @param        Array 2
* @param        Output
*/
__global__ void energy_sum(double2 *in1, double2 *in2, double *out);

/**
* @brief        Sums double2* and double2* to an output double2
* @ingroup      gpu
* @param        Array 1
* @param        Array 2
* @param        Output
*/
__global__ void sum(double2 *in1, double2 *in2, double2 *out);


/**
* @brief	Return the squared magnitude of a complex number. $|(a+\textrm{i}b)*(a-\textrm{i}b)|$
* @ingroup	gpu
* @param	in Complex number
* @return	Absolute-squared complex number
*/
__host__ __device__ double complexMagnitudeSquared(double2 in);

/**
* @brief        Complex magnitude of a double2 array
* @ingroup      gpu
*/
__global__ void complexMagnitudeSquared(double2 *in, double *out);

/**
* @brief        Complex magnitude of a double2 array
* @ingroup      gpu
*/
__global__ void complexMagnitudeSquared(double2 *in, double2 *out);

//##############################################################################
/**
 * Multiplication for linear, non-linear and phase-imprinting of the condensate.
 */
//##############################################################################

/**
* @brief	Kernel for complex multiplication
* @ingroup	gpu
* @param	in1 Wavefunction input
* @param	in2 Evolution operator input
* @param	out Pass by reference output for multiplcation result
*/
__global__ void cMult(double2* in1, double2* in2, double2* out);


/**
* @brief	Kernel for complex multiplication with nonlinear density term
* @ingroup	gpu
* @param	in1 Wavefunction input
* @param	in2 Evolution operator input
* @param	out Pass by reference output for multiplication result
* @param	dt Timestep for evolution
* @param	gState If performing real (1) or imaginary (0) time evolution
* @param	gDenConst a constant for evolution
*/
__global__ void cMultDensity(double2* in1, double2* in2, double2* out, double dt, int gstate, double g);

/**
* @brief        Complex field scaling and renormalisation. Used mainly post-FFT.
* @ingroup      gpu
* @param        in Complex field to be scaled (divided, not multiplied)
* @param        factor Scaling vector to be used
* @param        out Pass by reference output for result
*/
__global__ void vecMult(double2 *in, double *factor, double2 *out);

/**
* @brief        Complex field Summation
* @ingroup      gpu
* @param        in Complex field to be scaled (divided, not multiplied)
* @param        factor Scaling vector to be used
* @param        out Pass by reference output for result
*/
__global__ void vecSum(double2 *in, double *factor, double2 *out);

/**
* @brief        field scaling
* @ingroup      gpu
* @param        in field to be scaled (divided, not multiplied)
* @param        factor Scaling vector to be used
* @param        out Pass by reference output for result
*/
__global__ void vecMult(double *in, double *factor, double *out);

/**
* @brief        field Summation
* @ingroup      gpu
* @param        in field to be scaled (divided, not multiplied)
* @param        factor Scaling vector to be used
* @param        out Pass by reference output for result
*/
__global__ void vecSum(double *in, double *factor, double *out);

/**
* @brief        Complex field scaling and renormalisation. Used mainly post-FFT.
* @ingroup      gpu
* @param        in Complex field to be scaled (multiplied, not divided)
* @param        scaling factor to be used
* @param        out Pass by reference output for result
*/
__global__ void scalarMult(double2* in, double factor, double2* out);

/**
* @brief        field scaling and renormalisation. Used mainly post-FFT.
* @ingroup      gpu
* @param        in field to be scaled (multiplied, not divided)
* @param        scalaing factor to be used
* @param        out Pass by reference output for result
*/
__global__ void scalarMult(double* in, double factor, double* out);

/**
* @brief        Complex field scaling and renormalisation. Used mainly post-FFT.
* @ingroup      gpu
* @param        in Complex field to be scaled (multiplied, not divided)
* @param        complex scaling factor to be used
* @param        out Pass by reference output for result
*/
__global__ void scalarMult(double2* in, double2 factor, double2* out);


/**
* @brief        Conjugate of double2*.
* @ingroup      gpu
* @param        in Complex field to be conjugated
* @param        out Pass by reference output for result
*/
__global__ void vecConjugate(double2 *in, double2 *out);

/**
* @brief	Used as part of multipass to renormalise the wavefucntion
* @ingroup	gpu
* @param	in Complex field to be renormalised
* @param	dr Smallest area element of grid (dx*dy)
* @param	pSum GPU array used to store intermediate results during parallel summation
*/
__global__ void scalarDiv_wfcNorm(double2* in, double dr, double* pSum, double2* out);

//##############################################################################

/**
* @brief	Performs wavefunction renormalisation using parallel summation and applying scalarDiv_wfcNorm
* @ingroup	gpu
* @param	input Wavefunction to be renormalised
* @param	output Pass by reference return of renormalised wavefunction
* @param	pass Number of passes performed by routine
*/
__global__ void multipass(double2* input, double2* output, int pass);

/**
* @brief        Performs parallel summation of double arrays
* @ingroup      gpu
* @param        input double array
* @param        input double array after summation
*/
__global__ void multipass(double* input, double* output);

#endif
