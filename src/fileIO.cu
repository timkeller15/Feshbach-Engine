
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include "../include/fileIO.h"

namespace FileIO{

    /*
     * Reads datafile into memory.
     */
    double2* readIn(std::string fileR, std::string fileI,
                        int gSize){
        FILE *f;
        f = fopen(fileR.c_str(),"r");
        int i = 0;
        double2 *arr = (double2*) malloc(sizeof(double2)*gSize);
        double line;
        while(fscanf(f,"%lE",&line) > 0){
            arr[i].x = line;
            ++i;
        }
        fclose(f);
        f = fopen(fileI.c_str(),"r");
        i = 0;
        while(fscanf(f,"%lE",&line) > 0){
            arr[i].y = line;
            ++i;
        }
        fclose(f);
        return arr;
    }

    /*
     * Writes out double2 complex data files.
     */
    void writeOut(std::string buffer, std::string file, double2 *data, int length){
        FILE *f;
        sprintf ((char *)buffer.c_str(), "%s_real.dat", file.c_str());
        f = fopen (buffer.c_str(),"w");
        int i;
        for (i = 0; i < length; i++)
            fprintf (f, "%.16e\n",data[i].x);
        fclose (f);

        sprintf ((char *)buffer.c_str(), "%s_imag.dat", file.c_str());
        f = fopen (buffer.c_str(),"w");
        for (i = 0; i < length; i++)
            fprintf (f, "%.16e\n",data[i].y);
        fclose (f);

    }
}
