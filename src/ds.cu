
#include "../include/ds.h"
#include "../include/operators.h"
#include "../include/split_op.h"

/*----------------------------------------------------------------------------//
* GRID
*-----------------------------------------------------------------------------*/

// Function to store sobel_fft operators and stuff
void Grid::store(std::string id, cufftDoubleComplex *d2param){
    sobel[id] = d2param;
}

// Function to store integer into Grid->param_int
void Grid::store(std::string id, int iparam){
    param_int[id] = iparam;
}

// Function to store double into Grid->param_double
void Grid::store(std::string id, double dparam){
    param_double[id] = dparam;
}

// Function to store double* into param_dstar
void Grid::store(std::string id, double *dsparam){
    param_dstar[id] = dsparam;
}

// Function to store bool into param_bool
void Grid::store(std::string id, bool bparam){
    param_bool[id] = bparam;
}

// Function to store string into data_dir
void Grid::store(std::string id, std::string sparam){
    param_string[id] = sparam;
}

// Function to retrieve integer from Grid->param_int
int Grid::ival(std::string id){
    return param_int[id];
}

// Function to retrieve double from Grid->param_double
double Grid::dval(std::string id){
    auto it = param_double.find(id);
    if (it == param_double.end()){
        std::cout << "ERROR: could not find string " << id
                  << " in Grid::param_double." << '\n';
        assert(it != param_double.end());
    }
    return it->second;
}

// Function to retrieve double star values from param_dstar
double *Grid::dsval(std::string id){
    auto it = param_dstar.find(id);
    if (it == param_dstar.end()){
        std::cout << "ERROR: could not find string " << id
                  << " in Grid::param_dstar." << '\n';
        assert(it != param_dstar.end());
    }
    return it->second;
}

// Function to retrieve bool values from param_bool
bool Grid::bval(std::string id){
    auto it = param_bool.find(id);
    if (it == param_bool.end()){
        std::cout << "ERROR: could not find string " << id
                  << " in Grid::param_bool." << '\n';
        assert(it != param_bool.end());
    }
    return it->second;
}

// Function to retrieve string from data_dir
std::string Grid::sval(std::string id){
    auto it = param_string.find(id);
    if (it == param_string.end()){
        std::cout << "ERROR: could not find string " << id
                  << " in Grid::param_string." << '\n';
        assert(it != param_string.end());
    }
    return it->second;
}

// Function to call back the sobel operators
cufftDoubleComplex *Grid::cufftDoubleComplexval(std::string id){
    auto it = sobel.find(id);
    if (it == sobel.end()){
        std::cout << "ERROR: could not find string " << id
                  << " in Grid::sobel." << '\n';
        assert(it != sobel.end());
    }
    return it->second;
}

// Function for file writing
void Grid::write(std::string filename){
    std::ofstream output;
    output.open(filename);

    //Needed to recognise Params.dat as .ini format for python post processing
    output << "[Params]" <<"\n";

    // We simply iterate through the int and double param maps
    for (auto item : param_double){
        output << std::setprecision(12) << item.first << "=" << item.second << '\n';
    }

    for (auto item : param_int){
        output << item.first << "=" << item.second << '\n';
    }

    for (auto item : param_bool){
        output << item.first << "=" << item.second << '\n';
    }

    for (auto item : param_string){
        output << item.first << "=" << item.second << '\n';
    }

    output.close();
}
