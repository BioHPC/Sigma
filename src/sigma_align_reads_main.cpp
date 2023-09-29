//=============================================================================
// sigma_align_reads_main.cpp
//   : main cpp code for aligning reads to reference genome set.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified: 11/05/2013 - V1.0.1 is released.
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include <unistd.h>
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
              << "      - include bowtie search options and more" << std::endl
              << "    2. working directory (default: current running directory)" << std::endl
              << "      - if working_directory is not specified, the program will work in the current directory" << std::endl
              << "      - results will be generated in working directory" << std::endl
              << std::endl
              << "  [Options]" << std::endl
              << "    -h/--help" << std::endl
              << "    -v/--version" << std::endl
              << "    -p/--multi-processes <int>    # number of multi-processes (default: 1)" << std::endl
              << std::endl
              << "  [Outputs]" << std::endl
              << "    bam format aligned results are reported in each output genome directory" << std::endl
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
        std::string & config_base,          // (out) config base
        int & number_processes)             // (out) number of multi-processes
{
    // initialize variables
    int i, number_processes_default;
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
    number_processes_default  = 1;

    // option values
    working_directory = working_directory_default;
    config_filename   = config_filename_default;
    config_path       = config_path_default;
    number_processes  = number_processes_default;

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
        else if (arguments_vector[i] == "-p" || 
                 arguments_vector[i] == "--multi-processes") {
            std::stringstream(arguments_vector[++i]) >> number_processes;
            if (number_processes < 1) 
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
        usage(program_name);
        Utils::exitWithError("*** Error: cannot find a config file " + config_path + "\n");
    }

    // get config_base 
    config_filename = Utils::getFilename(config_path);
    config_base = Utils::getFilebase(config_filename);
}


//=============================================================================
// align reads to genome index by bowtie2
//=============================================================================
void runMultiProcesses(
        const std::string & config_base,        // (in) config file basename (string)
        int & number_processes,                 // (in) number of multi-processes (int)
        const std::vector<std::string> & genome_index_base_list,        // (in) genome index base list (vector)
        const std::string & output_parent_genome_directory,             // (in) output_parent_genome_directory (string)
        const std::vector<std::string> & output_genome_directory_list,  // (in) output genome directory list (vector)
        const std::string & bowtie_option_cmd,  // (in) bowtie option command (string)
        const std::string & read_option_cmd,    // (in) bowtie read files and option  (string)
        const std::string & samtools_cmd)       // (in) samtools command (string)
{
    // variables
    int number_workers;         // number of worker processes (number_processes - 1)
    int total_number_tasks;     // total number of jobs
    int init_number_tasks;      // initial assigning number of jobs (depends of number of workers)
    std::string genome_index_base;          // genome index base
    std::string output_genome_directory;    // output genome directory name
    std::string genome_name;    // genome name
    std::string bamout_path;    // sam mapping results path
    std::string bamlog_path;    // sam mapping log path
    std::string bowtie_cmd;     // botiew system command

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
    number_workers = number_processes;
    init_number_tasks = ((total_number_tasks <= number_workers) ? total_number_tasks : number_workers) ;

    // for process list
    std::vector<FILE*> process_list;

    // MPI send call for workers 
    for (int task_id=0; task_id<init_number_tasks; task_id++)
    {
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

        // setup popen
        std::string popen_result;
        FILE *process;

        // call popen (non-blocking popen for multi-processing)
        process = popen(bowtie_cmd.c_str(), "r");
        if (!process)
            Utils::exitWithError("*** Error: bowtie system command error.");
        
        // push back the process into list
        process_list.push_back(process);
    }

    int task_id = init_number_tasks;

    // run until no more work
    while (process_list.size() > 0 ) 
    {
        // get process from the list
        std::string popen_result;
        char buffer[512];
        FILE *process = process_list.at(0);

        // wait until the process finish
        while ( fgets(buffer, 512, process) != NULL )
            popen_result.append(buffer);

        // close the process
        pclose(process);

        // delete the process from the list
        process_list.erase(process_list.begin(), process_list.begin()+1);
        sleep(1);

        // if job remains, run new job
        if (task_id < total_number_tasks) {

            // get each genome info
            genome_index_base = genome_index_base_list.at(task_id);
            output_genome_directory = output_genome_directory_list.at(task_id);
            genome_name = Utils::getFilename(genome_index_base);
            bamout_path = output_genome_directory + Utils::getPathSeparator() 
                        + genome_name +  ".align.bam";
            bamlog_path = output_genome_directory + Utils::getPathSeparator() 
                        + genome_name +  ".align.log";

            // final bowtie command
            bowtie_cmd = "("  + bowtie_option_cmd + " -x " + genome_index_base + " " + read_option_cmd 
                       + " | " + samtools_cmd + " -o " + bamout_path + " - ) 2> " + bamlog_path;

            // setup popen
            FILE *new_process;

            // call popen 
            new_process = popen(bowtie_cmd.c_str(), "r");
            if (!new_process)
                Utils::exitWithError("*** Error: bowtie system command error.");
        
            // push back to list
            process_list.push_back(new_process);

            // increase task_id
            task_id++;
        }
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
    int number_processes;
    std::string working_directory, config_path, config_base, output_parent_genome_directory;

    // initialize arguments
    initializeArguments(argc, argv, version, // (in)
        working_directory, config_path, config_base, number_processes); // (out)
        
    // for elapsed time
    clock_t start_time, finish_time;
    time(&start_time);

    // display work start and time record
    std::cout << std::endl 
        << "********************************************************************************" << std::endl
        << Utils::currentDateTime() << " Beginning " << program_name << " V" << version <<  std::endl;
              
    // [Step 1] prepare
    std::cout << "  [Step 1] Prepare alignment: Running -> " << std::flush;

    // parse config and get bowtie option
    std::string reference_genome_directory, bowtie_option_cmd, read_option_cmd, samtools_cmd;

    // parse config and get bowtie option
    getBowtieOption(
        config_path,                // (in) config file path
        reference_genome_directory, // (out) reference genome dir
        bowtie_option_cmd,          // (out) bowtie option command
        read_option_cmd,            // (out) bowtie read files and option 
        samtools_cmd);              // (out) samtools command

    // search genome directory list and save to vector set-up
    std::vector<std::string> genome_directory_list, 
        genome_index_base_list, output_genome_directory_list;
    std::vector< std::vector<std::string> > genome_fasta_path_list;

    // search genome directory list and save to vector
    searchGenomeDirectoryList(
        working_directory,              // (in) working directory
        reference_genome_directory,     // (in) reference genome directory
        genome_directory_list,          // (out) genome directory list (vector)
        genome_fasta_path_list,         // (out) genome fasta path list (vec of vec))
        genome_index_base_list,         // (out) genome index base list (vector)
        output_parent_genome_directory, // (out) output_parent_genome_directory (string)
        output_genome_directory_list);  // (out) output genome directory list (vector)

    // [Step 1] done
    std::cout << "Done!" << std::endl;
   
    // [Step 2] align reads to each genome
    std::cout << "  [Step 2] Align reads to each genome: Running -> " << std::flush;

    // [Step 2] align reads to genome index by bowtie2
    runMultiProcesses(
        config_base,                    // (in) config file basename (string)
        number_processes,               // (in) number of multi-processes (int)
        genome_index_base_list,         // (in) genome index base list (vector)
        output_parent_genome_directory, // (in) output_parent_genome_directory (string)
        output_genome_directory_list,   // (in) output genome directory list (vector)
        bowtie_option_cmd,              // (in) bowtie option command (string)
        read_option_cmd,                // (in) bowtie read files and option  (string)
        samtools_cmd);                  // (in) samtools command (string)

    // [Step 2] done
    std::cout << "Done!" << std::endl;
   
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
