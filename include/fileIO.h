///@endcond
//##############################################################################
/**
 *  @file    fileIO.h
 *  @author  Lee J. O'Riordan (mlxd)
 *  @date    12/11/2015
 *  @version 0.1
 *
 *  @brief Routines for input and output of simulation data.
 *
 *  @section DESCRIPTION
 *  The functions herein are used to write the simulation data to text-based
 *  files (HDF was planned, but for simplicity I removed it). Data from previous
 *  simulations can also be read into memory.
 */
 //##############################################################################

#ifndef FILEIO_H
#define FILEIO_H
#include "../include/ds.h"
#include <vector>
#include <string>

/** Check source file for further information on functions **/
namespace FileIO {

    /**
    * @brief	Reads in the real and imaginary components from text files
    * @ingroup	helper
    *
    * @param	*fileR Name of data file of real components
    * @param	*fileI Name of data file of imaginary components
    * @param	xDim Size of x-grid
    * @param	yDim Size of y-grid
    * @return	*double2 Memory address of read-in data. Complex only
    */
    double2 *readIn(std::string fileR, std::string fileI, int gSize);

    /**
    * @brief	Writes the specified double2 array to a text file
    * @ingroup	helper
    *
    * @param	*buffer Char buffer for use by function internals. char[100] usually
    * @param	*file Name of data file name for saving to
    * @param	*data double2 array to be written out
    * @param	length Overall length of the file to write out
    */
    void writeOut(std::string buffer, std::string file, double2 *data, int length);
}


#endif
