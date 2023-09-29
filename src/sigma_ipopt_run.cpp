//=============================================================================
// sigma_ipopt_run.cpp
//   : Initialize IPOPT with options.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified: 10/01/2016 - V1.1.0 is released.
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include <omp.h>

#include "IpIpoptApplication.hpp"
#include "sigma_ipopt_nlp.h"
#include "utils.h"
#include "parse_config.h"
#include "sigma_config.h"
#include "sigma_core.h"


//=============================================================================
// run IPOPT
//=============================================================================
int runIpopt(const std::string & config_path, 
             const std::string & working_directory, 
             const std::string & input_qmatrix, 
             const int number_threads)
{
    // [Step 1] prepare
    std::cout << "  ** Check the qmatrix model: Running -> " << std::flush;

    // initialize
    std::string qmatrix_filename, qmatrix_filebase, qmatrix_filepath;

    // get qmatrix filename and path
    if (input_qmatrix != "") 
        qmatrix_filename = Utils::getFilename(input_qmatrix);
    else
        qmatrix_filename = "sigma_out.qmatrix.txt";
    qmatrix_filebase = working_directory + Utils::getFilebase2(qmatrix_filename);
    qmatrix_filepath = working_directory + qmatrix_filename;

    // get minimum_relative_abundance
    SigmaConfig config(config_path);
    double minimum_relative_abundance = config.getMinimumRelativeAbundance();

    // Create a new instance of your nlp
    //  (use a SmartPtr, not raw)
    SmartPtr<TNLP> mynlp = new SIGMA_IPOPT_NLP(qmatrix_filepath, 
                                               qmatrix_filebase, 
                                               number_threads,
                                               minimum_relative_abundance);

    // for ipopt output
    std::string ipopt_output = qmatrix_filebase + ".ipopt.txt";

    // [Step 1] done
    std::cout << "Done!" << std::endl;
   
    // [Step 2] align reads to each genome
    std::cout << "  ** Solve the model using IPOPT: Running -> " << std::flush;

    // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    SmartPtr<IpoptApplication> app = new IpoptApplication();

    // Change some options
    // Note: The following choices are only examples, they might not be
    //       suitable for your optimization problem.
    //app->Options()->SetNumericValue("tol", 1e-8);
    app->Options()->SetStringValue("linear_solver", "mumps");
    //app->Options()->SetNumericValue("mumps pivtol", 1e-5);
    //app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetIntegerValue("max_iter", 100);
    app->Options()->SetStringValue("output_file", ipopt_output);

    // Intialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int) status;
    }

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);

    if (status == Solve_Succeeded) {
        std::cout << "  ** The problem solved!" << std::endl << std::endl;
    }
    else {
        std::cout << "  ** The problem FAILED!" << std::endl << std::endl;
    }

    return (int) status;
}
