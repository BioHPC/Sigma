//=============================================================================
// sigma_jackknife_mpi_main.cpp
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
              << "    sigma_out.stat_jackknife.txt" << std::endl << std::endl
              << std::endl;
}


//=============================================================================
// initialize arguments 
//=============================================================================
void initializeArguments(int argc, char **argv, // (in) argv
        std::string & version,                  // (in) program version
        std::string & working_directory,        // (out) working directory
        std::string & config_path,              // (out) config path
        int & number_threads)                   // (out) number of multi-processes
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
    number_threads  = number_threads_default;

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
        // multi-threads
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
                   const std::vector<unsigned int> & genome_index_list)
{
    // MPI variables
    int rank;               // process rank
    int number_processes;   // total number of processes
    int number_workers;     // number of worker processes (number_threads - 1)
    int total_number_tasks; // total number of jobs
    int init_number_tasks;  // initial assigning number of jobs (depends of number of workers)
    int task_id;            // taks id based on index of vector
    int result;             // result of MPI call

    // MPI call
    MPI_Status status;

    // get number of initial assigning jobs
    number_processes   = proc_size;
    total_number_tasks = genome_index_list.size();
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
void slaveProcess(
        int & proc_id,                          // (in) process (rank) ID
        int & proc_size,                        // (in) number of processes (cores)
        const std::string & working_directory,  // (in) working directory
        const std::string & program_directory,  // (in) program directory
        const int & number_threads,             // (in) number of threads
        const std::string & qmatrix_filename,   // (in) qmatrix filename
        const std::string & qmatrix_filepath,   // (in) qmatrix filepath
        const std::vector<std::string> & qmatrix_comments_list,       // (in) qmatrix comments
        const std::vector<unsigned int> & genome_index_list,          // (in) genome_index_list
        const std::map<unsigned int, std::string> & genome_name_map,  // (in) genome name map
        const std::map<unsigned int, double> & align_rate_map,        // (in) align rate map
        const std::vector<std::string> & read_id_data,                // (in) read_id data
        const std::vector< std::vector<unsigned int> > & genome_index_data,   // (in) genome_index_data
        const std::vector< std::vector<double> > & qvalue_data)       // (in) qvalue_data
{
    // MPI status
    MPI_Status status;

    // variables
    int task_id;
    std::string qmatrix_filebase;           // qmatrix filebase (sigma_config.qmatrix.txt -> sigma_config)
    std::string jackknife_qmatrix_filename; // new generated Q matrix filename
    std::string jackknife_qmatrix_filepath; // new generated Q matrix filepath
    std::string jackknife_gvector_filename; // new generated gvector filename
    std::string jackknife_gvector_filepath; // new generated gvector filepath
    std::string system_command;

    // loop until getting DIETAG
    while(true)
    {
        // MPI_Recv receives a task
        MPI_Recv(&task_id, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG)
            break;

        // jackknife qmatrix filename and filepath
        qmatrix_filebase = Utils::getFilebase2(qmatrix_filename);
        jackknife_qmatrix_filename = qmatrix_filebase + ".jackknife." + Utils::intToString(task_id) + ".qmatrix.txt";
        jackknife_qmatrix_filepath = working_directory + jackknife_qmatrix_filename;

        // jackknife gvector filename and filepath
        jackknife_gvector_filename = qmatrix_filebase + ".jackknife." + Utils::intToString(task_id) + ".gvector.txt";
        jackknife_gvector_filepath = working_directory + jackknife_gvector_filename;

        // do-while until get the bootstrap_qmatrix
        do  
        {
            // write jackknife qmatrix 
            writeJackknifeQmatrix(proc_id,
                                  task_id,
                                  jackknife_qmatrix_filepath,
                                  qmatrix_comments_list,
                                  genome_index_list,
                                  genome_name_map,
                                  align_rate_map,
                                  read_id_data, 
                                  genome_index_data,
                                  qvalue_data);
        }   
        while(!Utils::isFileExist(jackknife_qmatrix_filepath));

        // std::cout << "The task id is " << task_id << std::endl;
        // get system command
        system_command = program_directory + "sigma-solve-model"
                       + " -t " + Utils::intToString(number_threads) 
                       + " -i " + jackknife_qmatrix_filepath + " > /dev/null 2>&1";

        // do-while until finish the run of system call
        do  
        {   
            // system call
            system(system_command.c_str());
        }   
        while(!Utils::isFileExist(jackknife_gvector_filepath));

        // remove jackknife qmatrix
        Utils::ifFileExistRemove(jackknife_qmatrix_filepath);

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
    std::string qmatrix_filename, qmatrix_filepath;
    unsigned int longest_read_length;
    double mismatch_probability;

    // initialize MPI
    MPI_Init(&argc,&argv);

    // get current process ID and size
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    // initialize arguments
    initializeArguments(argc, argv, version,  // (in)
        working_directory, config_path, number_threads);  // (out)

    // record start time
    if (proc_id == 0 ) 
    {
        // display work start and time record
        start_time = MPI_Wtime();
        std::cout << std::endl 
            << "********************************************************************************" << std::endl
            << Utils::currentDateTime() << " Beginning " << program_name << " V" << version <<  std::endl;

        // [Step 1] prepare
        std::cout << "  [Step 1] Prepare jackknifeping with Q-matirx: Running -> " << std::flush;
    }

    // get config values
    SigmaConfig config(config_path);
    qmatrix_filename = "sigma_out.qmatrix.txt";
    mismatch_probability = config.getMismatchProbability();
    longest_read_length = 100;

    // get qmatrix_path
    qmatrix_filepath = working_directory + qmatrix_filename;

    // get qmatrix data
    std::vector<std::string> qmatrix_comments_list;
    std::vector<unsigned int> genome_index_list;
    std::map<unsigned int, std::string> genome_name_map;
    std::map<unsigned int, double> align_rate_map;
    std::vector<std::string> read_id_data;
    std::vector< std::vector<unsigned int> > genome_index_data;
    std::vector< std::vector<double> > qvalue_data;
    // check path exist
    if (Utils::isFileExist(qmatrix_filepath))
    {   
        getJackknifeQmatrixData(qmatrix_filepath,       // (in)
                                qmatrix_comments_list,  // (out)
                                genome_index_list,      // (out)
                                genome_name_map,        // (out)
                                align_rate_map,         // (out)
                                read_id_data,           // (out)
                                genome_index_data,      // (out)
                                qvalue_data);           // (out)
    } 
    else 
    {
        Utils::exitWithError("*** Error: Failed to open " + qmatrix_filepath);
    }

    // for mater process
    if (proc_id == 0 ) 
    {
        // [Step 1] done
        std::cout << "Done!" << std::endl;

        // [Step 2] align reads to each genome
        std::cout << "  [Step 2] Run sigma-solve-model with resampled Q-matrix: Running -> " << std::flush;

        masterProcess(proc_size,          // (in) number of processes (cores)
                      genome_index_list); // (in) genome index list
    }
    // for slave processes
    else 
    {
        // return values for workers
        slaveProcess(proc_id,               // (in) process (rank) ID
                     proc_size,             // (in) number of processes (cores)
                     working_directory,     // (in) working directory
                     program_directory,     // (in) program directory
                     number_threads,      // (in) number of threads
                     qmatrix_filename,      // (in) qmatrix filename
                     qmatrix_filepath,      // (in) qmatrix filepath
                     qmatrix_comments_list, // (in) qmatrix comments 
                     genome_index_list,     // (in) qmatrix data
                     genome_name_map,       // (in) genome name map
                     align_rate_map,        // (in) align rate map
                     read_id_data,          // (in) read_id data
                     genome_index_data,     // (in) qmatrix data
                     qvalue_data);          // (in) qmatrix data
    }
        
    // get statistics
    if (proc_id == 0 ) 
    {
        // [Step 2] done
        std::cout << "Done!" << std::endl;

        // [Step 3] calculate statistics 
        std::cout << "  [Step 3] Calculate statistics: Running -> " << std::flush;
        calculateJackknifeGvectorStat(working_directory,     // (in) working directory
                                      qmatrix_filename,      // (in) qmatrix filename
                                      qmatrix_filepath,      // (in) qmatrix filepath
                                      longest_read_length,   // (in) longest read length
                                      mismatch_probability,  // (in) mismatch probability
                                      qmatrix_comments_list, // (in) qmatrix comments 
                                      genome_index_list,     // (in) qmatrix data
                                      genome_name_map,       // (in) genome name map
                                      align_rate_map,        // (in) align rate map
                                      read_id_data,          // (in) read_id data
                                      genome_index_data,     // (in) qmatrix data
                                      qvalue_data);          // (in) qmatrix data

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
