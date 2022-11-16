#include "../include/parser.h"
#include "../include/operators.h"

struct stat st = {0};

Grid parseArgs(int argc, char** argv){

    // Creates initial grid class for the parameters
    Grid par;
    int opt;

    // Setting default values
    par.store("Ngrid",256);
    par.store("xDim",256);
    par.store("yDim",256);
    par.store("zDim",256);
    par.store("posmax",10.0);
    par.store("dt", 1e-4);
    par.store("Tf",0.5);
    par.store("imaginarytime",false);
    par.store("N", (int)1e4);
    par.store("input_wfc", (std::string)"wfc_in");
    par.store("target_wfc", (std::string)"wfc_tar");
    par.store("samplestep",1000);
    par.store("gi",1.0);
    par.store("gf",0.8);
    par.store("data_dir", (std::string)"data/");
    par.store("write_wfc", false);
    par.store("rampmode",0);
    par.store("batchmode",false);
    par.store("filename",(std::string)"results");


    static struct option long_options[] =
    {
        {"Ngrid", required_argument, NULL, 'a'},
        {"posmax", required_argument, NULL, 'b'},
        {"dt", required_argument, NULL, 'c'},
        {"Tf", required_argument, NULL, 'd'},
        {"imaginarytime", no_argument, NULL, 'e'},
        {"atoms", required_argument, NULL, 'f'},
        {"gi", required_argument, NULL, 'g'},
        {"gf", required_argument, NULL, 'h'},
        {"rampmode", required_argument, NULL, 'i'},
        {"input", required_argument, NULL, 'j'},
        {"samplestep", required_argument, NULL, 'k'},
        {"directory", required_argument, NULL, 'l'},
        {"target", required_argument, NULL, 'm'},
        {"writewfc", no_argument, NULL, 'n'},
        {"batchmode", no_argument, NULL, 'o'},
        {"filename", required_argument, NULL, 'p'},
        {NULL, 0, NULL, 0}
    };

    optind = 1;

    while ((opt = getopt_long(argc, argv,
           "a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:", long_options, NULL)) !=-1)
    {
        switch (opt)
        {
            case 'a':
            {
                int Ngrid = atoi(optarg);
                printf("Ngrid: %d\n",Ngrid);
                par.store("Ngrid",(int)Ngrid);
                par.store("xDim",(int)Ngrid);
                par.store("yDim",(int)Ngrid);
                par.store("zDim",(int)Ngrid);
                break;
            }
            case 'b':
            {
                double posmax = atof(optarg);
                printf("Posmax: %3.1f\n",posmax);
                par.store("posmax",(double)posmax);
                break;
            }
            case 'c':
            {
                double dt = atof(optarg);
                printf("Timestep: %g\n",dt);
                par.store("dt",(double)dt);
                break;
            }
            case 'd':
            {
                double Tf = atof(optarg);
                printf("Timespan: %3.3f\n",Tf);
                par.store("Tf",(double)Tf);
                break;
            }
            case 'e':
            {
                printf("Imaginary time-evolution engaged.\n");
                par.store("imaginarytime",true);
                break;
            }
            case 'f':
            {
                int N = atof(optarg);
                printf("Atoms: %d\n",N);
                par.store("N",(int)N);
                break;
            }
            case 'g':
            {
                double gi = atof(optarg);
                printf("Initial interaction: %3.3f\n",gi);
                par.store("gi",gi);
                break;
            }
            case 'h':
            {
                double gf = atof(optarg);
                printf("Final interaction: %3.3f\n",gf);
                par.store("gf",gf);
                break;
            }
            case 'i':
            {
                int rampmode = atoi(optarg);
                switch (rampmode)
                {
                    case 0:
                    {
                        printf("Ramp Mode: TRA\n");
                        break;
                    }
                    case 1:
                    {
                        printf("Ramp Mode: STA\n");
                        break;
                    }
                    case 2:
                    {
                        printf("Ramp Mode: constant\n");
                        break;
                    }
                    default:
                        abort();
                }
                par.store("rampmode",(int)rampmode);
                break;
            }
            case 'j':
            {
                std::string input_wfc = optarg;
                std::cout << "Reading input wave function " << input_wfc << " from file \n";
                par.store("input_wfc",input_wfc);
                break;
            }
            case 'k':
            {
                int samplestep = atoi(optarg);
                printf("Samplestep: %d\n",samplestep);
                par.store("samplestep",(int)samplestep);
                break;
            }
            case 'l':
            {
                std::string data_dir = optarg;
                std::cout << "Data directory: " << data_dir << '\n';
                par.store("data_dir", data_dir + "/");
                break;
            }
            case 'm':
            {
                std::string target_wfc = optarg;
                std::cout << "Reading target wave function " << target_wfc << " from file \n";
                par.store("target_wfc",target_wfc);
                break;
            }
            case 'n':
            {
                printf("Writing out the final wave function.\n");
                par.store("write_wfc",true);
                break;
            }
            case 'o':
            {
                printf("Batch mode engaged. Reading values of Tf from file ramp_durations.txt\n");
                par.store("batchmode",true);
                break;
            }
            case 'p':
            {
                std::string filename = optarg;
                std::cout << "Filename: " << filename << '\n';
                par.store("filename",filename);
                break;
            }
            case '?':
            {
                if (optopt == 'c') {
                    fprintf (stderr,
                             "Option -%c requires an argument.\n", optopt);
                }
                else if (isprint (optopt)) {
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                }
                else {
                    fprintf (stderr,
                             "Unknown option character `\\x%x'.\n",optopt);
                }
                return par;
            default:
                abort ();
            }
        }
    }

    std::string data_dir = par.sval("data_dir");

    // Setting variables
    if (stat(data_dir.c_str(), &st) == -1) {
        mkdir(data_dir.c_str(), 0700);
    }

    return par;
}
