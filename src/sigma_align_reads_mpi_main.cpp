//=============================================================================
// sigma_align_reads_mpi_main.cpp
//   : MPI main cpp code for aligning reads to reference genome set.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified: 11/05/2013 - V1.0.1 (Beta) is relaeased
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include <mpi.h>
#include "utils.h"
#include "parse_config.h"
#include "sigma_config.h"
#include "sigma_core.h"

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
    std::cout << std::endl 
              << "  [Usage]" << std::endl 
              << "    " << program_name << " [options] -c <config file path> -w <working directory>" << std::endl
              << std::endl
              << "  [Inputs]" << std::endl
              << "    1. config file path (default: sigma_config.cfg)" << std::endl
              << "      - if config file is not specified, the program will search it in the working directory" << std::endl
              << "      - include bowtie search options and more" << std::endl
              << "    2. working directory (default: current running directory)" << std::endl
              << "      - if working_directory is not specified, the program will work in the current directory" << std::endl
              << "      - results will be generated in working directory" << std::endl
              << std::endl
              << "  [Options]" << std::endl
              << "    -h/--help" << std::endl
              << "    -v/--version" << std::endl
              << std::endl
              << "  [Outputs]" << std::endl
              << "    bam format aligned results are reported in each genome directory" << std::endl
              << std::endl;
}


//=============================================================================
// Initialize arguments 
//=============================================================================
void initializeArguments(
        int argc, char **argv,              // (in) argv
        std::string & version,              // (in) program version
        std::string & working_directory,    // (out) working directory
        std::string & config_path,          // (out) config path
        std::string & config_base)          // (out) config base
{
    // initialize variables
    int i;
    std::string config_filename, program_name, program_path,
        working_directory_default, config_filename_default, config_path_default;

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

    // get config_base 
    config_filename = Utils::getFilename(config_path);
    config_base = Utils::getFilebase(config_filename);
}


//=============================================================================
// Master process 
//=============================================================================
void masterProcess(std::string & output_parent_genome_directory,
        const std::vector<std::string> & output_genome_directory_list)
{
    // MPI variables
    int rank;               // process rank
    int number_processes;   // total number of processes
    int number_workers;     // number of worker processes (number_processes - 1)
    int total_number_tasks; // total number of jobs
    int init_number_tasks;  // initial assigning number of jobs (depends of number of workers)
    int task_id;            // taks id based on index of vector
    int result;             // result of MPI call
    int iter_num;
    int max_iter_num = 1000;

    // MPI call
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);

    // if output_parent_genome_directory doesn't exist, mkdir the directory
    if (!Utils::isDirectory(output_parent_genome_directory)) {
        std::string system_command = "mkdir " + output_parent_genome_directory;
        // system call
        const int res = system(system_command.c_str());
        if (res!=0) 
            Utils::exitWithError("*** Error: Failed command: " + system_command);
    }

    // create sub directories
    for(std::vector<int>::size_type i = 0; i != output_genome_directory_list.size(); i++) {
        std::string output_genome_directory = output_genome_directory_list[i];
        // if directory doesn't exist, mkdir the directory
        if (!Utils::isDirectory(output_genome_directory)) {
            std::string system_command = "mkdir " + output_genome_directory;
            // system call
            const int res = system(system_command.c_str());
            if (res!=0) 
                Utils::exitWithError("*** Error: Failed CMD: " + system_command);
        }
    }
    
    // get number of initial assigning jobs
    total_number_tasks = output_genome_directory_list.size();
    number_workers = number_processes -1;
    init_number_tasks = ((total_number_tasks <= number_workers) ? total_number_tasks : number_workers) ;

    // MPI send call for workers 
    for (rank=1; rank<=init_number_tasks; rank++)
    {
        task_id = rank-1;
        MPI_Send(&task_id,         // message buffer
                 1,                // one data item
                 MPI_INT,          // data item is an integer
                 rank,             // destination process rank
                 WORKTAG,          // user chosen message tag
                 MPI_COMM_WORLD);  // default communicator
    }

    // remaining jobs
    if ((int)total_number_tasks > number_workers)
    {
        task_id = number_workers;
        while (task_id < (int) total_number_tasks)
        {
            // MPI receive call
            MPI_Recv(&result,        // message buffer
                     1,              // one data item
                     MPI_INT,        // data item is an integer
                     MPI_ANY_SOURCE, // receive from any sender
                     MPI_ANY_TAG,    // any type of message
                     MPI_COMM_WORLD, // default communicator
                     &status);       //info about the received message

            // MPI send for reaming job
            MPI_Send(&task_id,           // message buffer
                     1,                  // one data item
                     MPI_INT,            // data item is an integer
                     status.MPI_SOURCE,  // destination process rank
                     WORKTAG,            // user chosen message tag
                     MPI_COMM_WORLD);    // default communicator

            // task index ++
            task_id ++;
        }
    }

    // MPI receive call for workers 
    for (rank=1; rank<=init_number_tasks; rank++)
    {
        // MPI receive call
        MPI_Recv(&result,        // message buffer
                 1,              // one data item
                 MPI_INT,        // data item is an integer
                 MPI_ANY_SOURCE, // receive from any sender
                 MPI_ANY_TAG,    // any type of message
                 MPI_COMM_WORLD, // default communicator
                 &status);       //info about the received message
    }

    // system() call sometimes makes error with OpenMPI, so check and run again until sucess
    iter_num = 1;
    while (iter_num < max_iter_num)
    {
        std::vector<int> remain_genome_index_list;

        // search remain genome directory list and save to vector
        remainGenomeDirectoryList(output_parent_genome_directory, // (in), output_parent_genome_directory genome directory
            remain_genome_index_list);  // (out) remain genome index list

        if (remain_genome_index_list.size() == 0) {
            break; 
        }
        else
        {
            std::cout << "  ** iter_num = " << Utils::intToString(iter_num) << std::endl;
            std::cout << "  ** remain_genome_index_list.size() = " << Utils::intToString(remain_genome_index_list.size()) << std::endl;

            // get number of initial assigning jobs
            total_number_tasks = remain_genome_index_list.size();
            number_workers = number_processes -1;
            init_number_tasks = ((total_number_tasks <= number_workers) ? total_number_tasks : number_workers) ;

            // MPI send call for workers 
            for (rank=1; rank<=init_number_tasks; rank++)
            {
                task_id = remain_genome_index_list.at(rank-1);
                MPI_Send(&task_id, 1, MPI_INT, rank, WORKTAG, MPI_COMM_WORLD);
            }

            // MPI receive call for workers 
            for (rank=1; rank<=init_number_tasks; rank++)
            {
                // MPI receive call
                MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
            }
        }
        iter_num++;
    }
        
    // tell all the slaves to exit by sending an empty message with the DIETAG.
    for (rank=1; rank<= number_workers; rank++) {
        MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    }
}


//=============================================================================
// Slave process 
//=============================================================================
void slaveProcess(const std::string & config_base,
                  const std::vector<std::string> & genome_index_base_list,
                  const std::vector<std::string> & output_genome_directory_list,
                  const std::string & bowtie_option_cmd,
                  const std::string & read_option_cmd,
                  const std::string & samtools_cmd)
{
    MPI_Status status;

    // variables
    int task_id;
    std::string genome_index_base, output_genome_directory, genome_name;
    std::string bamout_path, bamlog_path, bowtie_cmd;
    std::string system_run_log_path;   // system run log output path

    while(true)
    {
        MPI_Recv(&task_id, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIETAG)
            break;

        // get each genome info
        genome_index_base = genome_index_base_list.at(task_id);
        output_genome_directory = output_genome_directory_list.at(task_id);
        genome_name = Utils::getFilename(genome_index_base);
        bamout_path = output_genome_directory + Utils::getPathSeparator() 
                    + genome_name + ".align.bam";
        bamlog_path = output_genome_directory + Utils::getPathSeparator() 
                    + genome_name + ".align.log";

        // final bowtie command
        bowtie_cmd = "("  + bowtie_option_cmd + " -x " + genome_index_base + " " + read_option_cmd 
                   + " | " + samtools_cmd + " -o " + bamout_path + " - ) 2> " + bamlog_path;

        // system call
        system(bowtie_cmd.c_str());

        // MPI_Send
        MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}


//=============================================================================
// Main 
//=============================================================================
int main(int argc, char **argv)
{

    // version check
    std::string version = "0.0.1 (Alpha)";

    // get program name
    std::string program_path = argv[0];
    std::string program_name = Utils::getProgramName(program_path);

    // initialize variables 
    int proc_id, proc_size;    /* MPI rank -> proc_id, MPI size -> proc_size */
    std::string working_directory, config_path, config_base, output_parent_genome_directory;
    double start_time, finish_time, elapsed_time;

    // initialize MPI
    MPI_Init(&argc,&argv);

    // get current process ID and size
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    // initialize arguments
    initializeArguments(argc, argv, version, working_directory, config_path, config_base); 

    // record start time
    if (proc_id == 0 ) {
        start_time = MPI_Wtime();
        // display work start and time record
        std::cout << std::endl 
            << "********************************************************************************" << std::endl
            << Utils::currentDateTime() << " Beginning " << program_name << " V" << version <<  std::endl;

        // [Step 1] prepare
        std::cout << "  [Step 1] Prepare alignment: Running -> " << std::flush;
    }

    // parse config and get bowtie option
    std::string reference_genome_directory, bowtie_option_cmd, read_option_cmd, samtools_cmd;

    // parse config and get bowtie option
    getBowtieOption(
        config_path,                // (in) config file path
        reference_genome_directory, // (out) reference genome dir
        bowtie_option_cmd,          // (out) bowtie option command
        read_option_cmd,            // (out) bowtie read files and option 
        samtools_cmd);              // (out) samtools command

    // search genome directory list and save to vector
    std::vector<std::string> genome_directory_list, 
        genome_index_base_list, output_genome_directory_list;
    std::vector< std::vector<std::string> > genome_fasta_path_list;

    // search genome directory list and save to vector
    searchGenomeDirectoryList(
        working_directory,              // (in), working directory
        reference_genome_directory,     // (in), reference genome directory
        genome_directory_list,          // (out) genome directory list (vector)
        genome_fasta_path_list,         // (out) genome fasta path list (vec of vec))
        genome_index_base_list,         // (out) genome index base list (vector)
        output_parent_genome_directory, // (out) output_parent_genome_directory (string)
        output_genome_directory_list);  // (out) output genome directory list (vector)
   
    // for mater process
    if (proc_id == 0 ) {

        // [Step 1] done
        std::cout << "Done!" << std::endl;

        // [Step 2] align reads to each genome
        std::cout << "  [Step 2] Align reads to each genome: Running -> " << std::flush;

        masterProcess(output_parent_genome_directory, // (in) output_parent_genome_directory
                      output_genome_directory_list);  // (in) output_genome_directory_list
    }
    // for slave processes
    else {
        slaveProcess(
            config_base, 
            genome_index_base_list,       // (in) genome index base list
            output_genome_directory_list, // (in) output genome directory list
            bowtie_option_cmd,            // (in) bowtie option command
            read_option_cmd,              // (in) bowtie read files and option 
            samtools_cmd);                // (in) samtools command
    }
        
    // record start time
    if (proc_id == 0 ) {

        // [Step 2] done
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
