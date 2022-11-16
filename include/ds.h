///@endcond
//##############################################################################
/**
 *  @file    ds.h
 *  @author  James R. Schloss (leios) and Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief Dastructure for simulation runtime parameters
 *
 *  @section DESCRIPTION
 *      This file holds necessary classes and structs for all GPUE simulations.
 *      EqnNode and EqnNode_gpu are for dynamic parsing and Grid is for
 *      general use.
 */
 //#############################################################################

#ifndef DS_H
#define DS_H
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <typeinfo>
#include <cassert>
#include <iostream>
#include <iomanip>

/*----------------------------------------------------------------------------//
* CLASSES
*-----------------------------------------------------------------------------*/

/**
 * @brief       Struct for an x, y, z position
 * @ingroup     data
 */
struct pos{
    double x, y, z;
};

/**
 * @brief       function pointer type
 * @ingroup     data
 */
typedef double (*fnPtr) (double, double);

/**
 * @brief       Class to hold the variable map and grid information
 * @ingroup     data
 */
class Grid{
    // Here we keep our variable map (unordered for performance)
    // and also grid information. Note that dx dy, and dz are in param_double
    private:
        typedef void (*functionPtrV)(double*, double*, double*, double*, double*);
        typedef void (*functionPtrwfc)(double*, double*, double*, double*, double2*);
        std::unordered_map<std::string, int> param_int;
        std::unordered_map<std::string, double> param_double;
        std::unordered_map<std::string, double*> param_dstar;
        std::unordered_map<std::string, bool> param_bool;
        std::unordered_map<std::string, cufftDoubleComplex*> sobel;
        std::unordered_map<std::string, std::string> param_string;

        // List of all strings for parsing into the appropriate param map
        // 1 -> int, 2 -> double, 3 -> double*
        std::unordered_map<std::string, int> id_list;

    // Here we keep the functions to store variables and access grid data
    public:
        dim3 grid, threads;

        // Map for function pointers and keys K and V
        functionPtrV V_fn;
        functionPtrwfc wfc_fn;

        // placing grid parameters in public for now
        double *x, *y, *z, *xp, *yp, *zp;

        // Function to store sobel_fft operators into the sobel map
        void store(std::string id, cufftDoubleComplex* d2param);

        // Function to store integer into param_int
        void store(std::string id, int iparam);

        // Function to store double into param_double
        void store(std::string id, double dparam);

        // Function to store double* into param_dstar
        void store(std::string id, double *dsparam);

        // Function to store bool into param_bool
        void store(std::string id, bool bparam);

        // Function to store string into data_dir
        void store(std::string id, std::string sparam);

        // Function to retrieve integer value from param_int
        int ival(std::string id);

        // Function to retrieve double value from param_double
        double dval(std::string id);

        // Function to retrieve double star values from param_dstar
        double *dsval(std::string id);

        // Function to retrieve bool from param_bool
        bool bval(std::string id);

        // Fucntion to retrieve string from data_dir
        std::string sval(std::string id);

        // Function to call back the sobel operators
        cufftDoubleComplex *cufftDoubleComplexval(std::string id);

        // Function for file writing
        void write(std::string filename);
};
typedef class Grid Grid;

#endif
