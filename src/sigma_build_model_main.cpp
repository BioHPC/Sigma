//=============================================================================
// sigma_build_model_main.cpp
//   : main cpp code for building a Q-matrix
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified: 11/05/2013 - V1.0.1 is released.
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include "utils.h"
#include "parse_config.h"
#include "sigma_config.h"
#include "sigma_core.h"


//=============================================================================
// Help usage 
//=============================================================================
void usage(std::string & program_name);


//=============================================================================
// Help usage 
//=============================================================================
void usage(std::string & program_name)
{
    std::cout << std::endl 
              << "  [Usage]" << std::endl 
              << "    " << program_name << " [options] -c <config file path> -w <working directory>" << std::endl
              << std::endl
              << "  [Inputs]" << std::endl
              << "    1. config file path (default: sigma_config.cfg)" << std::endl
              << "      - if config file is not specified, the program will search it in the working directory" << std::endl
              << "    2. working directory (default: current running directory)" << std::endl
              << "      - if working_directory is not specified, the program will work in the current directory" << std::endl
              << "      - results will be generated in working directory" << std::endl
              << std::endl
              << "  [Options]" << std::endl
              << "    -h/--help" << std::endl
              << "    -v/--version" << std::endl
              << std::endl
              << "  [Outputs]" << std::endl
              << "    Q matrix file: sigma_out.qmatrix.txt will be generated in the working directory." << std::endl
              << std::endl;
}


//=============================================================================
// Initialize arguments
//=============================================================================
void initializeArguments(
        int argc, char **argv,              // (in) argv
        std::string & version,              // (in) program version
        std::string & working_directory,    // (out) working directory
        std::string & config_path)          // (out) config path
{
    // initialize variables
    int i;
    std::string config_filename, program_name, program_path;
    std::string working_directory_default, config_filename_default, config_path_default;

    // grab command line arguments
    std::vector<std::string> arguments_vector;

    // push to arguments_vector
    while(argc--) {
        arguments_vector.push_back(*argv++);
    }   

    // default option values
    working_directory_default = std::string(".") + Utils::getPathSeparator();
    config_filename_default   = "sigma_config.cfg";
    config_path_default       = working_directory_default + config_filename_default;

    // option values
    working_directory = working_directory_default;
    config_filename   = config_filename_default;
    config_path       = config_path_default;

    // get program name
    program_path = arguments_vector[0];
    program_name = Utils::getProgramName(program_path);

    // get argements
    for(i = 1; i <= (int)arguments_vector.size()-1; i++) 
    {   
        // working directory
        if(arguments_vector[i] == "-w") {
            working_directory = arguments_vector[++i];
            if (*working_directory.rbegin() != Utils::getPathSeparator())  
                working_directory += Utils::getPathSeparator();
        }
        // config path
        else if (arguments_vector[i] == "-c") {
            config_path = arguments_vector[++i];
        }
        // version
        else if (arguments_vector[i] == "-v" || 
                 arguments_vector[i] == "-V" || 
                 arguments_vector[i] == "--version") {
            std::cout << std::endl << program_name << " V" << version 
                << std::endl << std::endl;
            exit(0);
        }
        // help usage
        else if (arguments_vector[i] == "-h" || 
                 arguments_vector[i] == "--help") {
            usage(program_name);
            exit(0);
        }
        // unknown option
        else {
            std::cerr << "*** Error: Unknown option " << arguments_vector[i] 
                << std::endl << std::endl;
            usage(program_name);
            exit(1);
        }
    }

    // only -w is provided
    if ((working_directory != working_directory_default) && (config_path == config_path_default)) {
        config_path = working_directory + config_filename;
    }

    // check config path exist
    if (!Utils::isFileExist(config_path)) {
        usage(program_name);
        Utils::exitWithError("*** Error: cannot find a config file " + config_path + "\n");
    }
}


//=============================================================================
// Main
//=============================================================================
int main(int argc, char **argv)
{

    // version check
    std::string version = "1.0.1 (Beta)";

    // get program name
    std::string program_path = argv[0];
    std::string program_name = Utils::getProgramName(program_path);

    // initialize variables 
    std::string working_directory, config_path, output_parent_genome_directory;

    // initialize arguments
    initializeArguments(argc, argv, version, // (in)
        working_directory, config_path); // (out)
        
    // for elapsed time
    clock_t start_time, finish_time;
    time(&start_time);

    // display work start and time record
    std::cout << std::endl 
        << "********************************************************************************" << std::endl
        << Utils::currentDateTime() << " Beginning " << program_name << " V" << version <<  std::endl;
              
    //=============================================================================
    // gernerate Q matrix
    //=============================================================================
    generateQMatrix(config_path, working_directory);
   
    // for elapsed time
    time(&finish_time);
    double elapsed_time = difftime(finish_time, start_time);

    // display elapsed time
    std::cout << Utils::currentDateTime() << " Ending " << program_name << std::endl
        << "Total Elapsed Time =  " << elapsed_time << " [seconds]" << std::endl
        << "********************************************************************************" 
        << std::endl << std::endl;
              
    return 0;
}
