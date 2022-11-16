#include "../include/split_op.h"
#include "../include/kernels.h"
#include "../include/fileIO.h"
#include "../include/parser.h"

#include "../include/evolution.h"
#include <string>
#include <iostream>


char buffer[100]; //Buffer for printing out. Need to replace by a better write-out procedure. Consider binary or HDF.
int device; //GPU ID choice.
double timeTotal;

void cudaHandleError(cudaError_t result) {
    if (result != cudaSuccess) {
        std::cout << "Cuda operation failed!\n Error code: " << result << '\n';
        exit(1);
    }
}

void cudaCheckError() {
    cudaError_t result = cudaGetLastError();
    if (result != cudaSuccess) {
        std::cout << "Cuda kernel failed!\n Error code: " << result << '\n';
        exit(1);
    }
}

void cufftHandleError(cufftResult result) {
    if (result != CUFFT_SUCCESS) {
        std::cout << "cufft operation failed!\n Error code: " << result << '\n';
    }
}

/*
 * General-purpose summation of an array on the gpu, storing the result in the first element
*/
void gpuReduce(double* data, int length, int threadCount) {
    dim3 block(length / threadCount, 1, 1);
    dim3 threads(threadCount, 1, 1);

    while((double)length/threadCount > 1.0){
        multipass<<<block,threads,threadCount*sizeof(double)>>>(&data[0],
                                                                &data[0]);
        cudaCheckError();
        length /= threadCount;
        block = (int) ceil((double)length/threadCount);
    }
    multipass<<<1,length,threadCount*sizeof(double)>>>(&data[0],
                                                       &data[0]);
    cudaCheckError();
}

/*
 * Used to perform parallel summation on WFC for normalisation.
 */
void parSum(double2* gpuWfc, Grid &par){
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double dz = par.dval("dz");
    dim3 threads = par.threads;
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    dim3 grid_tmp(xDim*yDim*zDim, 1, 1);
    int gsize = xDim*yDim*zDim;
    double dg = dx*dy*dz;

    dim3 block(grid_tmp.x/threads.x, 1, 1);

    double *density;
    cudaHandleError( cudaMalloc((void**) &density, sizeof(double)*gsize) );

    complexMagnitudeSquared<<<par.grid, par.threads>>>(gpuWfc, density);
    cudaCheckError();

    gpuReduce(density, grid_tmp.x, threads.x);

    scalarDiv_wfcNorm<<<par.grid,par.threads>>>(gpuWfc, dg, density, gpuWfc);
    cudaCheckError();

    cudaHandleError( cudaFree(density) );
}

double energy_calc(Grid &par, double2* wfc, double interaction){
    double* K = par.dsval("K_gpu");
    double* V = par.dsval("V_gpu");

    dim3 grid = par.grid;
    dim3 threads = par.threads;

    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int gsize = xDim*yDim*zDim;

    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double dz = par.dval("dz");
    double dg = dx*dy*dz;

    cufftHandle plan;
    plan = par.ival("plan_3d");

    double renorm_factor = 1.0/pow(gsize,0.5);

    double2 *wfc_c, *wfc_k;
    double2 *energy_r, *energy_k;
    double *energy;

    cudaHandleError( cudaMalloc((void **) &wfc_c, sizeof(double2)*gsize) );
    cudaHandleError( cudaMalloc((void **) &wfc_k, sizeof(double2)*gsize) );
    cudaHandleError( cudaMalloc((void **) &energy_r, sizeof(double2)*gsize) );
    cudaHandleError( cudaMalloc((void **) &energy_k, sizeof(double2)*gsize) );

    cudaHandleError( cudaMalloc((void **) &energy, sizeof(double)*gsize)  );

    // Finding conjugate
    vecConjugate<<<grid, threads>>>(wfc, wfc_c);
    cudaCheckError();

    // Momentum-space energy
    cufftHandleError( cufftExecZ2Z(plan, wfc, wfc_k, CUFFT_FORWARD) );
    scalarMult<<<grid, threads>>>(wfc_k, renorm_factor, wfc_k);
    cudaCheckError();

    vecMult<<<grid, threads>>>(wfc_k, K, energy_k);
    cudaCheckError();
    cudaHandleError( cudaFree(wfc_k) );

    cufftHandleError( cufftExecZ2Z(plan, energy_k, energy_k, CUFFT_INVERSE) );
    scalarMult<<<grid, threads>>>(energy_k, renorm_factor, energy_k);
    cudaCheckError();

    cMult<<<grid, threads>>>(wfc_c, energy_k, energy_k);
    cudaCheckError();

    // Position-space energy
    int N = par.ival("N");

    double *real_comp;
    cudaHandleError( cudaMalloc((void**) &real_comp, sizeof(double)*gsize) );
    complexMagnitudeSquared<<<grid, threads>>>(wfc, real_comp);
    cudaCheckError();
    scalarMult<<<grid, threads>>>(real_comp,
                                  0.5*N*interaction,
                                  real_comp);
    cudaCheckError();
    vecSum<<<grid, threads>>>(real_comp, V, real_comp);
    cudaCheckError();
    vecMult<<<grid, threads>>>(wfc, real_comp, energy_r);
    cudaCheckError();

    cudaHandleError( cudaFree(real_comp) );

    cMult<<<grid, threads>>>(wfc_c, energy_r, energy_r);
    cudaCheckError();

    energy_sum<<<grid, threads>>>(energy_r, energy_k, energy);
    cudaCheckError();

    cudaHandleError( cudaFree(energy_r) );
    cudaHandleError( cudaFree(energy_k) );

    gpuReduce(energy, gsize, threads.x);

    double sum = 0;

    cudaHandleError( cudaMemcpy(&sum, energy, sizeof(double),
                                cudaMemcpyDeviceToHost) );

    sum *= dg;

    cudaHandleError( cudaFree(energy) );
    cudaHandleError( cudaFree(wfc_c) );

    return sum;
}
