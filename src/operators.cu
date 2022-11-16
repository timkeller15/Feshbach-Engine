#include "../include/operators.h"
#include "../include/split_op.h"
#include "../include/kernels.h"


// Function to check whether a file exists
std::string filecheck(std::string filename){

    struct stat buffer = {0};
    if (stat(filename.c_str(), &buffer) == -1){
        std::cout << "File " << filename << " does not exist!" << '\n';
        std::cout << "Please select a new file:" << '\n';
        std::cin >> filename;
        filename = filecheck(filename);
    }

    return filename;
}

/*----------------------------------------------------------------------------//
* GPU KERNELS
*-----------------------------------------------------------------------------*/

// Function to generate momentum grids
void generate_p_space(Grid &par){

    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    double xMax = par.dval("xMax");
    double yMax = par.dval("yMax");
    double zMax = par.dval("zMax");

    double pxMax = par.dval("pxMax");
    double pyMax = par.dval("pyMax");
    double pzMax = par.dval("pzMax");

    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double dz = par.dval("dz");

    double dpx = par.dval("dpx");
    double dpy = par.dval("dpy");
    double dpz = par.dval("dpz");


    double *x, *y, *z, *px, *py, *pz,
           *x_gpu, *y_gpu, *z_gpu,
           *px_gpu, *py_gpu, *pz_gpu;

    x = (double *) malloc(sizeof(double) * xDim);
    y = (double *) malloc(sizeof(double) * yDim);
    z = (double *) malloc(sizeof(double) * zDim);
    px = (double *) malloc(sizeof(double) * xDim);
    py = (double *) malloc(sizeof(double) * yDim);
    pz = (double *) malloc(sizeof(double) * zDim);


    for(int i=0; i<xDim/2; ++i){
        x[i] = -xMax + i*dx;
        x[i + (xDim/2)] = i*dx;

        px[i] = i*dpx;
        px[i + (xDim/2)] = -pxMax + i*dpx;

    }
    for(int i=0; i<yDim/2; ++i){
        y[i] = -yMax + i*dy;
        y[i + (yDim/2)] = i*dy;

        py[i] = i*dpy;
        py[i + (yDim/2)] = -pyMax + i*dpy;

    }
    for(int i=0; i<zDim/2; ++i){
        z[i] = -zMax + i*dz;
        z[i + (zDim/2)] = i*dz;

        pz[i] = i*dpz;
        pz[i + (zDim/2)] = -pzMax + i*dpz;

    }

    par.store("x",x);
    par.store("y",y);
    par.store("z",z);
    par.store("px",px);
    par.store("py",py);
    par.store("pz",pz);

    // Now move these items to the gpu
    cudaHandleError( cudaMalloc((void**) &x_gpu, sizeof(double) * xDim) );
    cudaHandleError( cudaMalloc((void**) &y_gpu, sizeof(double) * yDim) );
    cudaHandleError( cudaMalloc((void**) &z_gpu, sizeof(double) * zDim) );
    cudaHandleError( cudaMalloc((void**) &px_gpu, sizeof(double) * xDim) );
    cudaHandleError( cudaMalloc((void**) &py_gpu, sizeof(double) * yDim) );
    cudaHandleError( cudaMalloc((void**) &pz_gpu, sizeof(double) * zDim) );

    cudaHandleError( cudaMemcpy(x_gpu, x, sizeof(double)*xDim, cudaMemcpyHostToDevice) );
    cudaHandleError( cudaMemcpy(y_gpu, y, sizeof(double)*yDim, cudaMemcpyHostToDevice) );
    cudaHandleError( cudaMemcpy(z_gpu, z, sizeof(double)*zDim, cudaMemcpyHostToDevice) );
    cudaHandleError( cudaMemcpy(px_gpu, px, sizeof(double)*xDim, cudaMemcpyHostToDevice) );
    cudaHandleError( cudaMemcpy(py_gpu, py, sizeof(double)*yDim, cudaMemcpyHostToDevice) );
    cudaHandleError( cudaMemcpy(pz_gpu, pz, sizeof(double)*zDim, cudaMemcpyHostToDevice) );

    par.store("x_gpu",x_gpu);
    par.store("y_gpu",y_gpu);
    par.store("z_gpu",z_gpu);
    par.store("px_gpu",px_gpu);
    par.store("py_gpu",py_gpu);
    par.store("pz_gpu",pz_gpu);
}

// This function is basically a wrapper to call the appropriate K kernel
void generate_K(Grid &par){

    // For k, we need xp, yp, and zp. These will also be used in generating
    // pAxyz parameters, so it should already be stored in par.
    double *px_gpu = par.dsval("px_gpu");
    double *py_gpu = par.dsval("py_gpu");
    double *pz_gpu = par.dsval("pz_gpu");
    double gSize = par.ival("gSize");

    // Creating K to work with
    double *K, *K_gpu;
    K = (double*)malloc(sizeof(double)*gSize);
    cudaHandleError( cudaMalloc((void**) &K_gpu, sizeof(double)*gSize) );

    simple_K<<<par.grid, par.threads>>>(px_gpu, py_gpu, pz_gpu, K_gpu);
    cudaCheckError();

    cudaHandleError( cudaMemcpy(K, K_gpu, sizeof(double)*gSize, cudaMemcpyDeviceToHost) );
    par.store("K",K);
    par.store("K_gpu",K_gpu);

}

// Simple kernel for generating K
__global__ void simple_K(double *xp, double *yp, double *zp, double *K){

    unsigned int gid = getGid3d3d();
    //std::cout << "GID: " << gid << '\n';
    unsigned int xid = blockDim.x*blockIdx.x + threadIdx.x;
    unsigned int yid = blockDim.y*blockIdx.y + threadIdx.y;
    unsigned int zid = blockDim.z*blockIdx.z + threadIdx.z;
    //std::cout << "xid: " << xid << " yid: " << yid " zid: "<< zid << '\n';
    K[gid] = 0.5*(xp[xid]*xp[xid] + yp[yid]*yp[yid] + zp[zid]*zp[zid]);
}

// function to generate V
void generate_fields(Grid &par){

    generate_p_space(par);
    generate_K(par);

    int gSize = par.ival("gSize");

    double dt = par.dval("dt");
    double *x_gpu = par.dsval("x_gpu");
    double *y_gpu = par.dsval("y_gpu");
    double *z_gpu = par.dsval("z_gpu");
    double *px_gpu = par.dsval("px_gpu");
    double *py_gpu = par.dsval("py_gpu");
    double *pz_gpu = par.dsval("pz_gpu");
    double *K_gpu = par.dsval("K_gpu");

    // Creating items list for kernels

    double *items, *items_gpu;
    int item_size = 3;
    items = (double*)malloc(sizeof(double)*item_size);
    cudaHandleError( cudaMalloc((void**) &items_gpu, sizeof(double)*item_size) );

    for (int i = 0; i < item_size; ++i){
        items[i] = 0;
    }
    items[0] = par.dval("xMax");
    items[1] = par.dval("yMax");
    items[2] = par.dval("zMax");

    cudaHandleError( cudaMemcpy(items_gpu, items, sizeof(double)*item_size, cudaMemcpyHostToDevice) );

    // Generating V

    double *V, *V_gpu;

    V = (double *)malloc(sizeof(double)*gSize);

    cudaMalloc((void **) &V_gpu, sizeof(double)*gSize);

    kharmonic_V<<<par.grid, par.threads>>>(x_gpu, y_gpu, z_gpu, items_gpu, V_gpu);
    cudaCheckError();

    cudaHandleError( cudaMemcpy(V, V_gpu, sizeof(double)*gSize, cudaMemcpyDeviceToHost) );

    // Generating wfc

    double2 *wfc, *wfc_gpu;

    wfc = (double2 *)malloc(sizeof(double2)*gSize);

    cudaHandleError( cudaMalloc((void**) &wfc_gpu, sizeof(double2)*gSize) );
    wfc = par.cufftDoubleComplexval("wfc");
    cudaHandleError( cudaMemcpy(wfc_gpu, wfc, sizeof(double2)*gSize, cudaMemcpyHostToDevice) );


    // generating aux fields.
    double2 *GV, *EV, *GK, *EK;
    double2 *GV_gpu, *EV_gpu, *GK_gpu, *EK_gpu;

    GV = (double2 *)malloc(sizeof(double2)*gSize);
    EV = (double2 *)malloc(sizeof(double2)*gSize);
    GK = (double2 *)malloc(sizeof(double2)*gSize);
    EK = (double2 *)malloc(sizeof(double2)*gSize);

    cudaHandleError( cudaMalloc((void**) &GV_gpu, sizeof(double2)*gSize) );
    cudaHandleError( cudaMalloc((void**) &EV_gpu, sizeof(double2)*gSize) );
    cudaHandleError( cudaMalloc((void**) &GK_gpu, sizeof(double2)*gSize) );
    cudaHandleError( cudaMalloc((void**) &EK_gpu, sizeof(double2)*gSize) );

    aux_fields<<<par.grid, par.threads>>>(V_gpu, K_gpu, dt, GV_gpu, EV_gpu, GK_gpu, EK_gpu);
    cudaCheckError();

    cudaHandleError( cudaMemcpy(GV, GV_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost) );
    cudaHandleError( cudaMemcpy(EV, EV_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost) );
    cudaHandleError( cudaMemcpy(GK, GK_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost) );
    cudaHandleError( cudaMemcpy(EK, EK_gpu, sizeof(double2)*gSize, cudaMemcpyDeviceToHost) );

    // Storing variables
    cudaHandleError( cudaFree(items_gpu) );
    cudaHandleError( cudaFree(GV_gpu) );
    cudaHandleError( cudaFree(EV_gpu) );
    cudaHandleError( cudaFree(GK_gpu) );
    cudaHandleError( cudaFree(EK_gpu) );

    cudaHandleError( cudaFree(x_gpu) );
    cudaHandleError( cudaFree(y_gpu) );
    cudaHandleError( cudaFree(z_gpu) );

    cudaHandleError( cudaFree(px_gpu) );
    cudaHandleError( cudaFree(py_gpu) );
    cudaHandleError( cudaFree(pz_gpu) );


    par.store("V_gpu",V_gpu);

    par.store("V",V);
    par.store("items", items);
    par.store("wfc", wfc);
    par.store("wfc_gpu", wfc_gpu);

    par.store("GV",GV);
    par.store("EV",EV);
    par.store("GK",GK);
    par.store("EK",EK);
}

    __global__ void kharmonic_V(double *x, double *y, double *z, double* items, double *V){

    int gid = getGid3d3d();
    int xid = blockDim.x*blockIdx.x + threadIdx.x;
    int yid = blockDim.y*blockIdx.y + threadIdx.y;
    int zid = blockDim.z*blockIdx.z + threadIdx.z;

    double V_x = x[xid];
    double V_y = y[yid];
    double V_z = z[zid];

    V[gid] = 0.5*(V_x*V_x + V_y*V_y + V_z*V_z);
}

__global__ void aux_fields(double *V, double *K, double dt,
                           double2* GV, double2* EV, double2* GK, double2* EK){
    int gid = getGid3d3d();

    GV[gid].x = exp(-V[gid]*0.5*dt);
    GK[gid].x = exp(-K[gid]*dt);
    GV[gid].y = 0.0;
    GK[gid].y = 0.0;

    EV[gid].x=cos(-V[gid]*0.5*dt);
    EV[gid].y=sin(-V[gid]*0.5*dt);
    EK[gid].x=cos(-K[gid]*dt);
    EK[gid].y=sin(-K[gid]*dt);
}

// Function to generate grids and treads for 2d and 3d cases
void generate_grid(Grid& par){

    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int xD = 1, yD = 1, zD = 1;
    int max_threads = 256;
    if (xDim < max_threads){
        max_threads = xDim;
    }

    if (xDim <= max_threads){
        par.threads.x = xDim;
        par.threads.y = 1;
        par.threads.z = 1;

        xD = 1;
        yD = yDim;
        zD = zDim;
    }
    else{
        int count = 0;
        int dim_tmp = xDim;
        while (dim_tmp > max_threads){
            count++;
            dim_tmp /= 2;
        }

        std::cout << "count is: " << count << '\n';

        par.threads.x = dim_tmp;
        par.threads.y = 1;
        par.threads.z = 1;
        xD = pow(2,count);
        yD = yDim;
        zD = zDim;
    }

    par.grid.x=xD;
    par.grid.y=yD;
    par.grid.z=zD;
}
