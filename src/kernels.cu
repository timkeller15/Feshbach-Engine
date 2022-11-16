#include "../include/ds.h"
#include "../include/kernels.h"
#include <stdio.h>

__device__ double2 subtract(double2 a, double2 b){
    return {a.x-b.x, a.y-b.y};
}
__device__ double2 add(double2 a, double2 b){
    return {a.x+b.x, a.y+b.y};
}
__device__ double2 pow(double2 a, double b){
    double r = sqrt(a.x*a.x + a.y*a.y);
    double theta = atan(a.y / a.x);
    return{pow(r,b)*cos(b*theta),pow(r,b)*sin(b*theta)};
}

__device__ double2 mult(double2 a, double b){
    return {a.x*b, a.y*b};
}
__device__ double2 mult(double2 a, double2 b){
    return {a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x};
}

__device__ unsigned int getGid3d3d(){
    int blockId = blockIdx.x + blockIdx.y * gridDim.x
                  + gridDim.x * gridDim.y * blockIdx.z;
    int threadId = blockId * (blockDim.x * blockDim.y * blockDim.z)
                   + (threadIdx.y * blockDim.x)
                   + (threadIdx.z * (blockDim.x * blockDim.y)) + threadIdx.x;
    return threadId;
}

__global__ void sum(double2 *in1, double2 *in2, double2 *out){
    int gid = getGid3d3d();
    out[gid].x = in1[gid].x + in2[gid].x;
    out[gid].y = in1[gid].y + in2[gid].y;
}

__global__ void energy_sum(double2 *in1, double2 *in2, double *out){
    int gid = getGid3d3d();
    out[gid] = in1[gid].x + in2[gid].x;
}

__host__ __device__ double complexMagnitudeSquared(double2 in){
    return in.x*in.x + in.y*in.y;
}

__global__ void complexMagnitudeSquared(double2 *in, double *out){
    int gid = getGid3d3d();
    out[gid] = in[gid].x*in[gid].x + in[gid].y*in[gid].y;
}

__global__ void complexMagnitudeSquared(double2 *in, double2 *out){
    int gid = getGid3d3d();
    out[gid].x = in[gid].x*in[gid].x + in[gid].y*in[gid].y;
    out[gid].y = 0;
}

/**
 * Performs complex multiplication of in1 and in2, giving result as out.
 */
__global__ void cMult(double2* in1, double2* in2, double2* out){
    unsigned int gid = getGid3d3d();
    double2 result;
    double2 tin1 = in1[gid];
    double2 tin2 = in2[gid];
    result.x = (tin1.x*tin2.x - tin1.y*tin2.y);
    result.y = (tin1.x*tin2.y + tin1.y*tin2.x);
    out[gid] = result;
}

/**
 * Performs multiplication of double* with double2*
 */
__global__ void vecMult(double2 *in, double *factor, double2 *out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = in[gid].x * factor[gid];
    result.y = in[gid].y * factor[gid];
    out[gid] = result;
}

__global__ void vecMult(double *in, double *factor, double *out){
    double result;
    unsigned int gid = getGid3d3d();
    result = in[gid] * factor[gid];
    out[gid] = result;
}


__global__ void vecSum(double2 *in, double *factor, double2 *out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = in[gid].x + factor[gid];
    result.y = in[gid].y + factor[gid];
    out[gid] = result;
}

__global__ void vecSum(double *in, double *factor, double *out){
    double result;
    unsigned int gid = getGid3d3d();
    result = in[gid] + factor[gid];
    out[gid] = result;
}

/**
 * Performs the non-linear evolution term of Gross--Pitaevskii equation.
 */
__global__ void cMultDensity(double2* in1, double2* in2, double2* out, double dt, int gstate, double g){
    double2 result;
    double gDensity;

    int gid = getGid3d3d();
    double2 tin1 = in1[gid];
    double2 tin2 = in2[gid];
    gDensity = g*complexMagnitudeSquared(in2[gid])*dt;

    if(gstate == 0){
        double tmp = in1[gid].x*exp(-gDensity);
        result.x = (tmp)*tin2.x - (tin1.y)*tin2.y;
        result.y = (tmp)*tin2.y + (tin1.y)*tin2.x;
    }
    else{
        double2 tmp;
        tmp.x = tin1.x*cos(-gDensity) - tin1.y*sin(-gDensity);
        tmp.y = tin1.y*cos(-gDensity) + tin1.x*sin(-gDensity);

        result.x = (tmp.x)*tin2.x - (tmp.y)*tin2.y;
        result.y = (tmp.x)*tin2.y + (tmp.y)*tin2.x;
    }
    out[gid] = result;
}

/**
 * Multiplies both components of vector type "in", by the value "factor".
 * Results given with "out".
 */
__global__ void scalarMult(double2* in, double factor, double2* out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = (in[gid].x * factor);
    result.y = (in[gid].y * factor);
    out[gid] = result;
}

__global__ void scalarMult(double* in, double factor, double* out){
    double result;
    unsigned int gid = getGid3d3d();
    result = (in[gid] * factor);
    out[gid] = result;
}

__global__ void scalarMult(double2* in, double2 factor, double2* out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = (in[gid].x * factor.x - in[gid].y*factor.y);
    result.y = (in[gid].x * factor.y + in[gid].y*factor.x);
    out[gid] = result;
}

/**
 * As above, but normalises for wfc
 */
__global__ void scalarDiv_wfcNorm(double2* in, double dr, double* pSum, double2* out){
    unsigned int gid = getGid3d3d();
    double2 result;
    double norm = sqrt((pSum[0])*dr);
    result.x = (in[gid].x/norm);
    result.y = (in[gid].y/norm);
    out[gid] = result;
}


/**
 * Finds conjugate for double2*
 */
__global__ void vecConjugate(double2 *in, double2 *out){
    double2 result;
    unsigned int gid = getGid3d3d();
    result.x = in[gid].x;
    result.y = -in[gid].y;
    out[gid] = result;
}

/**
 * Routine for parallel summation. Can be looped over from host.
 */
__global__ void multipass(double2* input, double2* output, int pass){
    unsigned int tid = threadIdx.x + threadIdx.y*blockDim.x
                       + threadIdx.z * blockDim.x * blockDim.y;
    unsigned int bid = blockIdx.x + blockIdx.y * gridDim.x
                       + gridDim.x * gridDim.y * blockIdx.z;

    unsigned int gid = getGid3d3d();
    extern __shared__ double2 sdata[];
    sdata[tid] = input[gid];
    __syncthreads();
    for(int i = blockDim.x>>1; i > 0; i>>=1){
        if(tid < i){
            sdata[tid].x += sdata[tid + i].x;
            sdata[tid].y += sdata[tid + i].y;
        }
        __syncthreads();
    }
    if(tid==0){
        output[bid] = sdata[0];
    }
}

/**
 * Routine for parallel summation. Can be looped over from host.
 */
__global__ void multipass(double* input, double* output){
    unsigned int tid = threadIdx.x + threadIdx.y*blockDim.x
                       + threadIdx.z * blockDim.x * blockDim.y;
    unsigned int bid = blockIdx.x + blockIdx.y * gridDim.x
                       + gridDim.x * gridDim.y * blockIdx.z;

    unsigned int gid = getGid3d3d();
    extern __shared__ double sdatad[];
    sdatad[tid] = input[gid];
    __syncthreads();

    for(int i = blockDim.x>>1; i > 0; i>>=1){
        if(tid < i){
            sdatad[tid] += sdatad[tid + i];
        }
        __syncthreads();
    }
    if(tid==0){
        output[bid] = sdatad[0];
    }
}
