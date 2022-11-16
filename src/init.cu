#include "../include/init.h"
#define PI 3.141592653589793

void check_memory(Grid &par){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");

    bool energy_calc = true;

    int gSize = xDim*yDim*zDim;
    size_t free = 0;
    size_t total = 0;

    cudaHandleError( cudaMemGetInfo(&free, &total) );

    // Note that this check is specifically for the case where we need to keep
    // 8 double2* values on the GPU. This is not the case for dynamic fields
    // and the test should be updated accordingly as these are used more.
    size_t req_memory = 16*8*(size_t)gSize;
    if (energy_calc){
        req_memory += 4*16*(size_t)gSize;
    }
    if (free < req_memory){
        std::cout << "Not enough GPU memory for gridsize!\n";
        std::cout << "Free memory is: " << free << '\n';
        std::cout << "Required memory is: " << req_memory << '\n';
        if (energy_calc){
            std::cout << "Required memory for energy calc is: "
                      << 4*16*(size_t)gSize << '\n';
        }
        std::cout << "xDim is: " << xDim << '\n';
        std::cout << "yDim is: " << yDim << '\n';
        std::cout << "zDim is: " << zDim << '\n';
        std::cout << "gSize is: " << gSize << '\n';
        exit(1);
    }
}

int init(Grid &par){

    check_memory(par);

    // Re-establishing variables from parsed Grid class
    // Initializes uninitialized variables to 0 values
    std::string data_dir = par.sval("data_dir");
    int N = par.ival("N");
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    unsigned int gSize = xDim*yDim*zDim;

    double posmax = par.dval("posmax");

    cufftDoubleComplex *wfc;
    wfc = par.cufftDoubleComplexval("wfc");

    cufftHandle plan_3d;

    generate_grid(par);
    cudaCheckError();

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    double sum = 0.0;

    double xMax, yMax, zMax;
    if (posmax > 0){
        xMax = posmax;
        yMax = posmax;
        zMax = posmax;
    }
    else{
        xMax = 10.0;
        yMax = 10.0;
        zMax = 10.0;
    }
    par.store("xMax",xMax);
    par.store("yMax",yMax);
    par.store("zMax",zMax);

    double pxMax, pyMax, pzMax;
    pxMax = (PI/xMax)*(xDim>>1);
    pyMax = (PI/yMax)*(yDim>>1);
    pzMax = (PI/zMax)*(zDim>>1);
    par.store("pxMax",pxMax);
    par.store("pyMax",pyMax);
    par.store("pzMax",pzMax);

    double dx = xMax/(xDim>>1);
    double dy = yMax/(yDim>>1);
    double dz = zMax/(zDim>>1);

    par.store("dx",dx);
    par.store("dy",dy);
    par.store("dz",dz);

    double dpx, dpy, dpz;
    dpx = PI/(xMax);
    dpy = PI/(yMax);
    dpz = PI/(zMax);

    par.store("dpx",dpx);
    par.store("dpy",dpy);
    par.store("dpz",dpz);

    double dt = par.dval("dt");
    double Tf = par.dval("Tf");
    int steps = (int) (Tf/dt);
    if(!par.bval("batchmode")){
        std::cout << "Number of steps: " << steps << '\n';
    }
    par.store("steps", steps);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
    /* Initialise wavefunction, momentum, position,
       imaginary and real-time evolution operators . */
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    par.store("gSize", xDim*yDim*zDim);

    generate_fields(par);

    double *K = par.dsval("K");
    double *V = par.dsval("V");

    double *x = par.dsval("x");
    double *y = par.dsval("y");
    double *z = par.dsval("z");

    double2 *GV = par.cufftDoubleComplexval("GV");
    double2 *EV = par.cufftDoubleComplexval("EV");
    double2 *GK = par.cufftDoubleComplexval("GK");
    double2 *EK = par.cufftDoubleComplexval("EK");

    wfc = par.cufftDoubleComplexval("wfc");

    for(int i=0; i < gSize; i++ ){
        sum+=sqrt(wfc[i].x*wfc[i].x + wfc[i].y*wfc[i].y);
    }

    cufftHandleError( cufftPlan3d(&plan_3d, xDim, yDim, zDim, CUFFT_Z2Z) );

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    // Storing variables that have been initialized
    // Re-establishing variables from parsed Grid class
    // Initializes uninitialized variables to 0 values
    par.store("wfc", wfc);
    par.store("plan_3d", plan_3d);

    double gi = par.dval("gi");
    double gf = par.dval("gf");
    int rampmode = par.ival("rampmode");
    std::string mode, fname;
    switch (rampmode)
    {
        case 0:
        {
            mode = "TRA";
            break;
        }
        case 1:
        {
            mode = "STA";
            break;
        }
        case 2:
        {
            mode = "constant";
            break;
        }
        default:
            printf("Wrong rampmode!\n");
            abort();
    }

    char numstr[21]; // enough to hold all numbers up to 64-bits

    if(par.bval("imaginarytime")){
        fname = par.sval("input_wfc");
        sprintf(numstr, "_N%d_g%2.1f", N,gi);
        fname += numstr;
    } else {
        fname = par.sval("filename");
        if (par.bval("batchmode")){
            sprintf(numstr, "_N%d_gi%2.1f_gf%2.1f_batch_", N,gi,gf);
        } else {
            sprintf(numstr, "_N%d_gi%2.1f_gf%2.1f_Tf%2.1f_", N,gi,gf,Tf);
        }
        fname += numstr + mode;
    }

    par.store("filename",fname);
    std::cout << "Filename: " << fname << "\n";

    return 0;
}

void set_variables(Grid &par, bool ev_type){
    // Re-establishing variables from parsed Grid class
    // Note that 3d variables are set to nullptr's unless needed
    //      This might need to be fixed later
    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double2 *V_gpu;
    double2 *K_gpu;
    cufftDoubleComplex *wfc = par.cufftDoubleComplexval("wfc");
    cufftDoubleComplex *wfc_gpu = par.cufftDoubleComplexval("wfc_gpu");

    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int gsize = xDim*yDim*zDim;

    cudaHandleError( cudaMalloc((void**) &V_gpu, sizeof(double2)*gsize) );
    cudaHandleError( cudaMalloc((void**) &K_gpu, sizeof(double2)*gsize) );

    if (ev_type == 0){
        cufftDoubleComplex *GK = par.cufftDoubleComplexval("GK");
        cufftDoubleComplex *GV = par.cufftDoubleComplexval("GV");

        cudaHandleError( cudaMemcpy(K_gpu, GK, sizeof(cufftDoubleComplex)*gsize,
                                        cudaMemcpyHostToDevice) );
        cudaHandleError( cudaMemcpy(V_gpu, GV, sizeof(cufftDoubleComplex)*gsize,
                                        cudaMemcpyHostToDevice) );

        cudaHandleError( cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*gsize,
                                    cudaMemcpyHostToDevice) );
        par.store("K_gpu", K_gpu);
        par.store("V_gpu", V_gpu);
        par.store("wfc_gpu", wfc_gpu);

        free(GV); free(GK);
    }
    else if (ev_type == 1){

        cufftDoubleComplex *EV = par.cufftDoubleComplexval("EV");
        cufftDoubleComplex *EK = par.cufftDoubleComplexval("EK");

        cudaHandleError( cudaMemcpy(K_gpu, EK, sizeof(cufftDoubleComplex)*gsize,
                                    cudaMemcpyHostToDevice) );
        par.store("K_gpu", K_gpu);

        cudaHandleError( cudaMemcpy(V_gpu, EV, sizeof(cufftDoubleComplex)*gsize,
                                    cudaMemcpyHostToDevice) );
        par.store("V_gpu", V_gpu);

        cudaHandleError( cudaMemcpy(wfc_gpu, wfc, sizeof(cufftDoubleComplex)*gsize,
                                    cudaMemcpyHostToDevice) );

        par.store("wfc_gpu", wfc_gpu);

        free(EV);
        free(EK);
    }

}

int main(int argc, char **argv){

    Grid par = parseArgs(argc,argv);

    int device = 0;
    cudaHandleError( cudaSetDevice(device) );

    std::string buffer;
    time_t start,fin;
    time(&start);
    printf("Start: %s\n", ctime(&start));

    //************************************************************//
    /*
    * Initialise the Params data structure to track params and variables
    */
    //************************************************************//

    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int gSize = xDim * yDim * zDim;

    cufftDoubleComplex *wfc;

    int N = par.ival("N");
    double gi = par.dval("gi");
    double gf = par.dval("gf");
    std::string input_real, input_imag, target_real, target_imag;
    std::string data_dir = par.sval("data_dir");
    char numstr[21];
    sprintf(numstr, "_N%d_g%2.1f", N,gi);
    std::string input_wfc = par.sval("input_wfc") + numstr;
    sprintf(numstr, "_N%d_g%2.1f", N,gf);
    std::string target_wfc = par.sval("target_wfc") + numstr;

    input_real = data_dir + input_wfc + "_real.dat";
    input_imag = data_dir + input_wfc + "_imag.dat";
    target_real = data_dir + target_wfc + "_real.dat";
    target_imag = data_dir + target_wfc + "_imag.dat";

    printf("Loading input wavefunction...");
    wfc = FileIO::readIn(input_real,input_imag,gSize);
    par.store("wfc",wfc);
    printf("completed.\n");

    init(par);

    int runs = 1;
    double cycles[10000];
    cycles[0] = par.dval("Tf");

    // Perform the same ramp for several durations Tf loaded from file "ramp_durations.txt" in batchmode
    if(par.bval("batchmode")){

        runs = 0;
        std::string line;
        std::ifstream inputfile(data_dir + "ramp_durations.txt");

        while (std::getline(inputfile, line)){
                ++runs;
        }
        inputfile.clear();
        inputfile.seekg(0, std::ios::beg);

        int count = 0;
        while (std::getline(inputfile, line)){
                    cycles[count++] = std::stod(line);
        }
        inputfile.close();
    }

    std::string fname = par.sval("filename");

    if(par.bval("imaginarytime")){
        par.store("gf",gi);
        std::cout << "Imaginary-time evolution started..." << '\n';
        set_variables(par, 0);
        evolve(par, par.ival("steps"), 0);
    }
    else{
        std::cout << "Real-time evolution started..." << '\n';
        set_variables(par, 1);

        // Load target wave function from file for calculating fidelity
        cufftDoubleComplex *wf_ini, *wf_target;
        printf("Loading target wavefunctions...");
        wf_ini = FileIO::readIn(input_real,input_imag,gSize);
        wf_target = FileIO::readIn(target_real,target_imag,gSize);
        par.store("wf_ini",wf_ini);
        par.store("wf_target",wf_target);
        printf("completed!\n");

        for(int i=0; i<runs; i++){

            double Tf = cycles[i];
            double dt = par.dval("dt");
            int steps = (int) (Tf/dt);
            par.store("steps", steps);
            par.store("Tf",Tf);
            printf("Current Tf: %3.2f...\n",Tf);

            evolve(par, par.ival("steps"), 1);

            // Retrieve final wave function from GPU
            cufftDoubleComplex *wfc_gpu = par.cufftDoubleComplexval("wfc_gpu");
            cufftDoubleComplex *wf_final = par.cufftDoubleComplexval("wfc");
            cudaHandleError(cudaMemcpy(wf_final, wfc_gpu, sizeof(cufftDoubleComplex)*gSize,cudaMemcpyDeviceToHost));

            // Calculate fidelity after the ramp
            double fid = fidelity(par,wf_final,wf_target);
            std::ofstream file;
            file.open(data_dir + fname + ".dat", std::ios::out | std::ios::app);
            file << std::setprecision(12) << Tf << '\t' << par.dval("Ef") << '\t' << fid << '\n';
            file.close();
            printf("Fidelity: %1.10f\n",fid);

            // Reset wave function on GPU with initial one for next run
            cudaHandleError(cudaMemcpy(wfc_gpu, wf_ini, sizeof(cufftDoubleComplex)*xDim*yDim*zDim, cudaMemcpyHostToDevice));
        }
    }

    par.write(data_dir + fname + "_params.dat");

    if (par.bval("write_wfc")) {
        std::cout << "Writing final wave function" << '\n';
        cufftDoubleComplex *wfc_gpu = par.cufftDoubleComplexval("wfc_gpu");
        cufftDoubleComplex *wf_final = par.cufftDoubleComplexval("wfc");
        cudaHandleError(cudaMemcpy(wf_final, wfc_gpu, sizeof(cufftDoubleComplex)*gSize,cudaMemcpyDeviceToHost));
        FileIO::writeOut(buffer, data_dir + fname, wf_final, xDim*yDim*zDim);
    }

    time(&fin);
    printf("Finish: %s\n", ctime(&fin));
    printf("Total time: %ld seconds\n ",(long)fin-start);
    std::cout << '\n';
    return 0;
}


double fidelity(Grid &par, cufftDoubleComplex *in1, cufftDoubleComplex *in2){
    int xDim = par.ival("xDim");
    int yDim = par.ival("yDim");
    int zDim = par.ival("zDim");
    int gSize = xDim*yDim*zDim;

    double2 result;
    result.x = 0.0;
    result.y = 0.0;

    for(int i=0; i < gSize; i++ ){
        result.x += in1[i].x*in2[i].x + in1[i].y*in2[i].y;
        result.y += in1[i].x*in2[i].y - in1[i].y*in2[i].x;
    }

    double dx = par.dval("dx");
    double dy = par.dval("dy");
    double dz = par.dval("dz");

    result.x *= dx*dy*dz;
    result.y *= dx*dy*dz;

    return result.x*result.x + result.y*result.y;
}
