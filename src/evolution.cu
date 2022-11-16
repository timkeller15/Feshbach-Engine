#include "../include/evolution.h"
#include "../include/fileIO.h"

// Returns the current interaction strength during time evolution
// Available ramp modes are TRA, STA or constant
double gRamp(Grid &par, double time){

    int rampmode = par.ival("rampmode");
    int N = par.ival("N");
    double gi = par.dval("gi");
    double gf = par.dval("gf");
    double Tf = par.dval("Tf");
    double r1, r2;
    r1 = (double)gf/gi;
    r2 = (double)time/Tf;

    double a = 1.0 + (pow(r1,0.2) - 1.0)*(10.0*pow(r2,3.0) - 15.0*pow(r2,4.0) + 6.0*pow(r2,5.0)); // smoother step polynomial for scaling function
    double att = (pow(r1,0.2) - 1.0)*(60.0*r2 - 180.0*pow(r2,2.0) + 120.0*pow(r2,3.0))/(Tf*Tf); // and its second derivative

    double g;

    switch (rampmode)
    {
        case 0: //TRA, time-rescaled adiabatic reference from setting att = 0
        {
            g = N*gi*pow(a,5.0);
            break;
        }
        case 1: //STA, shortcut to adiabaticity for Thomas-Fermi regime
        {
            g = N*gi*pow(a,4.0)*(att + a);
            break;
        }
        case 2: //constant
        {
            g = N*gi;
            break;
        }
        default:
        {
            g = N*gi;
        }
    }

    return g;
}

// Output interaction ramp if necessary during testing
void write_ramp(Grid &par, std::string filename){
    std::ofstream outputt, outputg;
    outputt.open(filename + "_t");
    outputg.open(filename + "_g");
    double t;
    double g;
    double Tf = par.dval("Tf");

    for(int i=0; i < 100; ++i){
        t = (i/99.0)*Tf;
        g = gRamp(par,t);
        outputt << t << '\n';
        outputg << g << '\n';
    }

    outputt.close();
    outputg.close();
}

void evolve(Grid &par, int numSteps, unsigned int gstate){

    // Re-establishing variables from parsed Grid class
    std::string data_dir = par.sval("data_dir");
    double dt = par.dval("dt");
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double dz = par.dval("dz");

    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");

    int N = par.ival("N");
    double gi = par.dval("gi");
    double gf = par.dval("gf");
    double Tf = par.dval("Tf");

    cufftDoubleComplex *wfc = par.cufftDoubleComplexval("wfc");
    cufftDoubleComplex *gpuWfc = par.cufftDoubleComplexval("wfc_gpu");
    cufftDoubleComplex *K_gpu = par.cufftDoubleComplexval("K_gpu");
    cufftDoubleComplex *V_gpu = par.cufftDoubleComplexval("V_gpu");

    int gridSize = xDim * yDim * zDim;

    // getting data from Cuda class
    cufftHandle plan_3d = par.ival("plan_3d");

    dim3 threads = par.threads;
    dim3 grid = par.grid;

    // Because no two operations are created equally.
    // Multiplication is faster than divisions.
    double renorm_factor_nd=1.0/pow(gridSize,0.5);
    double renorm_factor_x=1.0/pow(xDim,0.5);
    double renorm_factor_y=1.0/pow(yDim,0.5);
    double renorm_factor_z=1.0/pow(zDim,0.5);

    //write_ramp(par,data_dir + "ramp");

    double Ei = energy_calc(par, gpuWfc, gi);
    par.store("Ei",Ei);
    printf("Initial energy: %3.10f\n",Ei);

    // Fourier Split-Step Method
    // Iterating through all of the steps in either g or esteps.
    for(int i=0; i < numSteps; ++i){
        double time = dt*i;

        // U_r(dt/2)*wfc
        cMultDensity<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc,0.5*dt,gstate,gRamp(par,time));
        cudaCheckError();

        // U_p(dt)*fft2(wfc)
        cufftHandleError( cufftExecZ2Z(plan_3d,gpuWfc,gpuWfc,CUFFT_FORWARD) );

        // Normalise
        scalarMult<<<grid,threads>>>(gpuWfc,renorm_factor_nd,gpuWfc);
        cudaCheckError();

        cMult<<<grid,threads>>>(K_gpu,gpuWfc,gpuWfc);
        cudaCheckError();

        cufftHandleError( cufftExecZ2Z(plan_3d,gpuWfc,gpuWfc,CUFFT_INVERSE) );

        // Normalise
        scalarMult<<<grid,threads>>>(gpuWfc,renorm_factor_nd,gpuWfc);
        cudaCheckError();

        // U_r(dt/2)*wfc
        cMultDensity<<<grid,threads>>>(V_gpu,gpuWfc,gpuWfc,0.5*dt,gstate,gRamp(par,time));
        cudaCheckError();

        if(gstate==0){
            parSum(gpuWfc, par);
        }
    }

    par.store("wfc", wfc);
    par.store("wfc_gpu", gpuWfc);

    double Ef = energy_calc(par, gpuWfc, gRamp(par,Tf)/N);
    par.store("Ef",Ef);
    printf("Final energy: %3.10f\n",Ef);
}
