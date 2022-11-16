///@endcond
//##############################################################################
/**
 *  @file    split_op.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief Host and device declarations for simulations.
 *
 *  @section DESCRIPTION
 *  These functions and variables are necessary for carrying out the GPUE
 *	simulations. This file will be re-written in an improved form in some
 *	future release.
 */
//##############################################################################

#ifndef SPLIT_OP_H
#define SPLIT_OP_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <ctype.h>
#include <getopt.h>
#include "ds.h"

#ifdef __linux
	#include<omp.h>
#elif __APPLE__
	//printf("OpenMP support disabled due to Clang/LLVM being behind the trend.",);
#endif

/**
 * @brief	Checks if a CUDA operation has succeeded. Prints to stdout.
 * @ingroup	data
 * @param	result Result code of CUDA operation
 */
void cudaHandleError(cudaError_t result);

/**
 * @brief Checks if the last CUDA operation has succeeded. Prints to stdout.
 * @ingroup data
 */
void cudaCheckError();

/**
 * @brief Checks if a cuFFT operation has succeeded. Prints to stdout.
 * @ingroup data
 * @param result Result code of cuFFT operation
 */
void cufftHandleError(cufftResult result);

/**
 * @brief Performs parallel summation on a flat array of data
 * @ingroup gpu
 * @param data Gpu-allocated array to reduce
 * @param length Length of data
 * @param threadCount number of threads to use
 */
void gpuReduce(double* data, int length, int threadCount);

/**
* @brief	Performs parallel summation and renormalises the wavefunction
* @ingroup	data
* @param	gpuWfc GPU memory location for wavefunction
* @param	par Parameter class
*/
void parSum(double2* gpuWfc, Grid &par);

/**
* @brief	Calculates the energy of the condensate.
* @ingroup	data
* @param	gpuWfc Device wavefunction array
* @param	Parameter class
* @return	$\langle \Psi | H | \Psi \rangle$
*/
double energy_calc(Grid &par, double2* wfc, double interaction);

#endif
