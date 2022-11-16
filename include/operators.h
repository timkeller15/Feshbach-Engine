//##############################################################################
/**
 *  @file    operators.h
 *  @author  James R. Schloss (Leios)
 *  @date    1/1/2017
 *  @version 0.1
 *
 *  @brief File to hold all operators for finding fields on the GPU
 *
 *  @section DESCRIPTION
 *      This file holds all operators for finding fields on the GPU
 */
 //#############################################################################

#ifndef OPERATORS_H
#define OPERATORS_H

#include "../include/ds.h"
#include <sys/stat.h>
#include <unordered_map>
//#include <boost/math/special_functions.hpp>

 /**
 * @brief       Determines if file exists, requests new file if it does not
 * @ingroup     data
 */

std::string filecheck(std::string filename);

/*----------------------------------------------------------------------------//
* GPU KERNELS
*-----------------------------------------------------------------------------*/

 /**
 * @brief      Function to generate momentum grids
 */
void generate_p_space(Grid &par);

 /**
 * @brief       This function calls the appropriate K kernel
 */
void generate_K(Grid &par);

 /**
 * @brief       Simple kernel for generating K
 */
__global__ void simple_K(double *xp, double *yp, double *zp, double *K);

 /**
 * @brief       Function to generate V
 */
void generate_fields(Grid &par);

 /**
 * @brief       Kernel to generate harmonic V
 */
__global__ void kharmonic_V(double *x, double *y, double *z, double *items, double *V);

/**
* @brief       Kernel to generate all auxiliary fields
*/

__global__ void aux_fields(double *V, double *K, double dt,
                          double2* GV, double2* EV, double2* GK, double2* EK);
                          
// Function to generate grid and treads
void generate_grid(Grid &par);

#endif
