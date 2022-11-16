///@endcond
//##############################################################################
/**
 *  @file    evolution.h
 *  @author  James Ryan Schloss (leios)
 *  @date    5/31/2016
 *  @version 0.1
 *
 *  @brief function for evolution.
 *
 *  @section DESCRIPTION
 *  These functions and variables are necessary for carrying out the GPUE
 *	simulations. This file will be re-written in an improved form in some
 *	future release.
 */
//##############################################################################

#ifndef EVOLUTION_H
#define EVOLUTION_H

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
#include "split_op.h"
#include "kernels.h"
#include "fileIO.h"


double gRamp(Grid &par, double time);
void write_ramp(Grid &par, std::string filename);

 /**
 * @brief       performs real or imaginary time evolution
 * @ingroup     data
 * @param       Parameter set
 * @param       Total number of steps
 * @param       Real (1) or imaginary (1) time evolution
 * @return      0 for success. See CUDA failure codes in cuda.h for other values
 */
void evolve(Grid &par, int numSteps, unsigned int gstate);
#endif
