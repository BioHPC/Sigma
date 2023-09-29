//=============================================================================
// sigma_bootstrap_mpi_main.cpp
//   : MPI main cpp code for statistical bootstrap method.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified: 11/05/2013 - V1.0.1 (Beta) is relaeased
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include <mpi.h>
#include "math.h"
#include "utils.h"
#include "parse_config.h"
#include "sigma_config.h"
#include "sigma_core.h"
#include "sigma_stat.h"

#define WORKTAG    1
#define DIETAG     2


//=============================================================================
// Help usage 
//=============================================================================
void usage(std::string & program_name);


//=============================================================================
// Help usage 
//=============================================================================
void usage(std::string & program_name)
{
    std::cout << "  [Usage]" << std::endl 
              << "    " << program_name << " [options] -c <config file path> -w <working directory>" << std::endl
              << std::endl
              << "  [Inputs]" << std::endl
              << "    1. config file path (default: sigma_config.cfg)" << std::endl
              << "      - if config file is not specified, the program will search it in the working directory" << std::endl
              << "    2. working directory (default: current running directory)" << std::endl
              << "      - if working_directory is not specified, the program will work in the current directory" << std::endl
              << "      - results will be generated in the working directory" << std::endl
              << std::endl
              << "  [Options]" << std::endl
              << "    -h/--help" << std::endl
              << "    -v/--version" << std::endl
              << "    -t/--multi-threads <int>    # number of threads (default: 1)" << std::endl
              << std::endl
              << "  [Outputs]" << std::endl
              << "    sigma_out.stat_bootstrap.txt" << std::endl << std::endl
              << std::endl;
}


//=============================================================================
// initialize arguments 
//=============================================================================
void initializeArguments(
        int argc, char **argv,              // (in) argv
        std::string & version,              // (in) program version
        std::string & working_directory,    // (out) working directory
        std::string & config_path,          // (out) config path
        int & number_threads)               // (out) number of multi-processes
{
    // initialize variables
    int i, number_threads_default;
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
    number_threads_default  = 1;

    // option values
    working_directory = working_directory_default;
    config_filename   = config_filename_default;
    config_path       = config_path_default;
    number_threads    = number_threads_default;

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
        // multi-processes
        else if (arguments_vector[i] == "-t" ||
                 arguments_vector[i] == "--multi-threads") {
            std::stringstream(arguments_vector[++i]) >> number_threads;
            if (number_threads < 1)
                Utils::exitWithError("*** Error: check -p option value.");
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
        Utils::exitWithError("*** Error: config file " + config_path + " couldn't be found!\n");
    }
}


//=============================================================================
// Master process
//=============================================================================
void masterProcess(int & proc_size,
                   const unsigned int & bootstrap_iteration_number)
{
    // MPI variables
    int rank;               // process rank
    int number_processes;   // total number of processes
    int number_workers;     // number of worker processes (number_processes - 1)
    int total_number_tasks; // total number of jobs
    int init_number_tasks;  // initial assigning number of jobs (depends of number of workers)
    int task_id;            // taks id based on index of vector
    int result;             // result of MPI call

    // MPI call
    MPI_Status status;

    // get number of initial assigning jobs
    number_processes   = proc_size;
    total_number_tasks = bootstrap_iteration_number;
    number_workers     = number_processes -1;
    init_number_tasks  = ((total_number_tasks <= number_workers) ? total_number_tasks : number_workers) ;

    // MPI send call for workers 
    for (rank=1; rank<=init_number_tasks; rank++)
    {
        task_id = rank-1;
        MPI_Send(&task_id,          // message buffer
                 1,                 // one data item
                 MPI_INT,           // data item is an integer
                 rank,              // destination process rank
                 WORKTAG,           // user chosen message tag
                 MPI_COMM_WORLD);   // default communicator
    }

    // remaining jobs
    if ((int)total_number_tasks > number_workers)
    {
        task_id = number_workers;
        while (task_id < (int) total_number_tasks)
        {
            // MPI receive call
            MPI_Recv(&result,           // message buffer
                     1,                 // one data item
                     MPI_INT,           // data item is an integer
                     MPI_ANY_SOURCE,    // receive from any sender
                     MPI_ANY_TAG,       // any type of message
                     MPI_COMM_WORLD,    // default communicator
                     &status);          //info about the received message

            // MPI send for reaming job
            MPI_Send(&task_id,          // message buffer
                     1,                 // one data item
                     MPI_INT,           // data item is an integer
                     status.MPI_SOURCE, // destination process rank
                     WORKTAG,           // user chosen message tag
                     MPI_COMM_WORLD);   // default communicator

            // task index ++
            task_id ++;
        }
    }

    // MPI receive call for workers 
    for (rank=1; rank<=init_number_tasks; rank++)
    {
        // MPI receive call (same as before)
        MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    // tell all the slaves to exit by sending an empty message with the DIETAG.
    for (rank=1; rank<= number_workers; rank++) 
    {
        // MPI send call with DIETAG
        MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    }
}


//=============================================================================
// Slave process
//=============================================================================
void slaveProcess(int & proc_id,    // (in) process (rank) ID
                  int & proc_size,  // (in) number of processes (cores)
                  const std::string & working_directory,    // (in) working directory
                  const std::string & program_directory,    // (in) program directory
                  const int & number_threads,               // (in) number of threads
                  const std::string & qmatrix_filename,     // (in) qmatrix filename
                  const std::string & qmatrix_filepath,     // (in) qmatrix filepath
                  const std::vector<std::string> & qmatrix_comments_list,   // (in) qmatrix comments
                  const std::vector<std::string> & qmatrix_data_list)       // (in) qmatrix data
{
    // MPI status
    MPI_Status status;

    // variables
    int task_id;
    std::string qmatrix_filebase; // qmatrix filebase (sigma_config.qmatrix.txt -> sigma_config)
    std::string bootstrap_qmatrix_filename; // new generated Q matrix filename
    std::string bootstrap_qmatrix_filepath; // new generated Q matrix filepath
    std::string bootstrap_gvector_filename; // new generated gvector filename
    std::string bootstrap_gvector_filepath; // new generated gvector filepath
    std::string system_command;

    // loop until getting DIETAG
    while(true)
    {
        // MPI_Recv receives a task
        MPI_Recv(&task_id, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG)
            break;

        // change random seed
        time_t timer;
        double seconds = time(&timer);
        //std::cout << "  ** seconds = " << seconds << std::endl;
        unsigned int seed = abs(int(fmod(((seconds*181)*((task_id)*359)), 104729)));
        //std::cout << "  ** seed = " << seed << std::endl;
        srand (seed);

        // bootstrap qmatrix filename and filepath
        qmatrix_filebase = Utils::getFilebase2(qmatrix_filename);
        bootstrap_qmatrix_filename = qmatrix_filebase + ".bootstrap." + Utils::intToString(task_id) + ".qmatrix.txt";
        bootstrap_qmatrix_filepath = working_directory + bootstrap_qmatrix_filename;

        // bootstrap gvector filename and filepath
        bootstrap_gvector_filename = qmatrix_filebase + ".bootstrap." + Utils::intToString(task_id) + ".gvector.txt";
        bootstrap_gvector_filepath = working_directory + bootstrap_gvector_filename;

        // do-while until get the bootstrap_qmatrix
        do
        {
            // write bootstrap qmatrix 
            writeBootstrapQmatrix(proc_id,
                                  bootstrap_qmatrix_filepath,
                                  qmatrix_comments_list,
                                  qmatrix_data_list);
        }
        while(!Utils::isFileExist(bootstrap_qmatrix_filepath));

        // debug
        //std::cout << "The task id is " << task_id << std::endl;

        // get system command
        system_command = program_directory + "sigma-solve-model"
                       + " -t " + Utils::intToString(number_threads) 
                       + " -i " + bootstrap_qmatrix_filepath  + " > /dev/null 2>&1";

        // do-while until finish the run of system call
        do
        {
            // debug
            std::cout << "System command: " << system_command << std::endl;

            // system call
            system(system_command.c_str());
        }
        while(!Utils::isFileExist(bootstrap_gvector_filepath));

        // remove bootstrap qmatrix
        Utils::ifFileExistRemove(bootstrap_qmatrix_filepath);

        // MPI_Send sends finish tag
        MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
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
    std::string program_directory = Utils::getProgramDir(program_path);

    // initialize variables 
    int proc_id, proc_size;    /* MPI rank -> proc_id, MPI size -> proc_size */
    int number_threads;
    std::string working_directory, config_path;
    double start_time, finish_time, elapsed_time;

    // config file variables
    unsigned int bootstrap_iteration_number;
    std::string qmatrix_filename, qmatrix_filepath;

    // initialize MPI
    MPI_Init(&argc,&argv);

    // get current process ID and size
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    // initialize arguments
    initializeArguments(argc, argv, version,                // (in)
        working_directory, config_path, number_threads);    // (out)

    // record start time
    if (proc_id == 0 ) {
        start_time = MPI_Wtime();
        // display work start and time record
        std::cout << std::endl 
            << "********************************************************************************" << std::endl
            << Utils::currentDateTime() << " Beginning " << program_name << " V" << version <<  std::endl;

        // [Step 1] prepare
        std::cout << "  [Step 1] Prepare bootstrapping with Q-matirx: Running -> " << std::flush;
    }

    // get config values
    SigmaConfig config(config_path);
    bootstrap_iteration_number = config.getBootstrapIterationNumber();
    qmatrix_filename = "sigma_out.qmatrix.txt";

    // get qmatrix_path
    qmatrix_filepath = working_directory + qmatrix_filename;

    // get qmatrix data
    std::vector<std::string> qmatrix_comments_list, qmatrix_data_list;
    // check path exist
    if (Utils::isFileExist(qmatrix_filepath))
    {   
        getBootstrapQmatrixData(qmatrix_filepath,       // (in)
                                qmatrix_comments_list,  // (out)
                                qmatrix_data_list);     // (out)
    }   
    else 
    {   
        Utils::exitWithError("*** Error: Failed to open " + qmatrix_filepath);
    }   

    // for mater process
    if (proc_id == 0 ) {

        // [Step 1] done
        std::cout << "Done!" << std::endl;

        // [Step 2] align reads to each genome
        std::cout << "  [Step 2] Run sigma-solve-model with resampled Q-matrix: Running -> " << std::flush;

        masterProcess(proc_size,                   // (in) number of processes (cores)
                      bootstrap_iteration_number); // (in) bootstrapping iteration number
    }
    // for slave processes
    else {
        slaveProcess(proc_id,            // (in) process (rank) ID
                     proc_size,          // (in) number of processes (cores)
                     working_directory,  // (in) working directory
                     program_directory,  // (in) program directory
                     number_threads,     // (in) number of threads
                     qmatrix_filename,   // (in) qmatrix filename
                     qmatrix_filepath,   // (in) qmatrix filepath
                     qmatrix_comments_list, // (in) qmatrix comments 
                     qmatrix_data_list);    // (in) qmatrix data
    }
        
    // get statistics
    if (proc_id == 0 ) {

        // [Step 2] done
        std::cout << "Done!" << std::endl;

        // [Step 3] calculate statistics 
        std::cout << "  [Step 3] Calculate statistics: Running -> " << std::flush;
        calculateBootstrapGvectorStat(working_directory,
                                      qmatrix_filename,
                                      bootstrap_iteration_number);

        // [Step 3] done
        std::cout << "Done!" << std::endl;

        // finish time and elpased time
        finish_time = MPI_Wtime();
        elapsed_time = double(finish_time - start_time);

        // display elapsed time
        std::cout << Utils::currentDateTime() << " Ending " << program_name << std::endl
            << "Total Elapsed Time =  " << elapsed_time << " [seconds]" << std::endl
            << "********************************************************************************"
            << std::endl << std::endl;
    }

    // finalize MPI
    MPI_Finalize();

    return 0;
}
