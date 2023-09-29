//=============================================================================
// gereme_config.cpp
//   : functions to get values of config file
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified:
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include "utils.h"
#include "parse_config.h"
#include "sigma_config.h"


//=============================================================================
// config keys and variables to string 
//=============================================================================
// [Program_Info]
const std::string SigmaConfig::program_info_str 
    = "[Program_Info]";
const std::string SigmaConfig::bowtie_directory_str 
    = program_info_str + "Bowtie_Directory";
const std::string SigmaConfig::samtools_directory_str 
    = program_info_str + "Samtools_Directory";

// [Data_Info]
const std::string SigmaConfig::data_info_str 
    = "[Data_Info]";
const std::string SigmaConfig::reference_genome_directory_str 
    = data_info_str + "Reference_Genome_Directory";
const std::string SigmaConfig::paired_end_reads_1_str 
    = data_info_str + "Paired_End_Reads_1";
const std::string SigmaConfig::paired_end_reads_2_str 
    = data_info_str + "Paired_End_Reads_2";
const std::string SigmaConfig::single_end_reads_str 
    = data_info_str + "Single_End_Reads";

// [Bowtie_Search]
const std::string SigmaConfig::bowtie_search_str 
    = "[Bowtie_Search]";
const std::string SigmaConfig::bowtie2_preset_str 
    = bowtie_search_str + "Bowtie2_Preset";
const std::string SigmaConfig::maximum_mismatch_count_str 
    = bowtie_search_str + "Maximum_Mismatch_Count";
const std::string SigmaConfig::minimum_fragment_length_str 
    = bowtie_search_str + "Minimum_Fragment_Length";
const std::string SigmaConfig::maximum_fragment_length_str 
    = bowtie_search_str + "Maximum_Fragment_Length";
const std::string SigmaConfig::bowtie_threads_number_str 
    = bowtie_search_str + "Bowtie_Threads_Number";

// [Model_Probability]
const std::string SigmaConfig::model_probability_str 
    = "[Model_Probability]";
const std::string SigmaConfig::mismatch_probability_str 
    = model_probability_str + "Mismatch_Probability";
const std::string SigmaConfig::minimum_relative_abundance_str
    = model_probability_str + "Minimum_Relative_Abundance";

// [Statistics]
const std::string SigmaConfig::statistics_str 
    = "[Statistics]";
const std::string SigmaConfig::bootstrap_iteration_number_str 
    = statistics_str + "Bootstrap_Iteration_Number";


//=============================================================================
// Constructor 
//=============================================================================
SigmaConfig::SigmaConfig(const std::string &config_path) 
  : ParseConfig(config_path)
{}   


//=============================================================================
// Deconstructor 
//=============================================================================
SigmaConfig::~SigmaConfig() 
{}   


//=============================================================================
// Get reads type 
//=============================================================================
void SigmaConfig::getReadsType(bool & paired_end_reads_flag, bool & fasta_reads_flag)
{
    // paired end case
    if (ParseConfig::keyExists(paired_end_reads_1_str) && ParseConfig::keyExists(paired_end_reads_2_str)) 
    {
        // paired_end = true
        paired_end_reads_flag = true;

        // get read file path
        std::string Paired_End_Reads_1 = ParseConfig::getValueOfKey(paired_end_reads_1_str);
        std::string Paired_End_Reads_2 = ParseConfig::getValueOfKey(paired_end_reads_2_str);

        // read files string -> list
        std::vector<std::string> Paired_End_Reads_1_list = Utils::commaStringToVector(Paired_End_Reads_1); 
        std::vector<std::string> Paired_End_Reads_2_list = Utils::commaStringToVector(Paired_End_Reads_2); 

        // check file exists
        for (std::vector<std::string>::iterator it = Paired_End_Reads_1_list.begin() ; it != Paired_End_Reads_1_list.end(); ++it)
            if (!Utils::isFileExist(*it))
                Utils::exitWithError("*** Error: Reads file " + *it + " does not exist. Check the config file!");
        
        for (std::vector<std::string>::iterator it = Paired_End_Reads_2_list.begin() ; it != Paired_End_Reads_2_list.end(); ++it)
            if (!Utils::isFileExist(*it))
                Utils::exitWithError("*** Error: Reads file " + *it + " does not exist. Check the config file!");
        
        // check read file is fasta or fastq
        if (Utils::isFastaFormat(Paired_End_Reads_1) && Utils::isFastaFormat(Paired_End_Reads_2))
            fasta_reads_flag = true;
        else if (Utils::isFastqFormat(Paired_End_Reads_1) && Utils::isFastqFormat(Paired_End_Reads_2))
            fasta_reads_flag = false;
        else
            Utils::exitWithError("*** Error: check read file format in the config file!");
    }
    // single end case
    else if (ParseConfig::keyExists(single_end_reads_str)) 
    {
        // paired_end = false
        paired_end_reads_flag = false;

        // get read file path
        std::string Single_End_Reads = ParseConfig::getValueOfKey(single_end_reads_str);

        // read files string -> list
        std::vector<std::string> Single_End_Reads_list = Utils::commaStringToVector(Single_End_Reads); 

        // check file exists
        for (std::vector<std::string>::iterator it = Single_End_Reads_list.begin() ; it != Single_End_Reads_list.end(); ++it)
            if (!Utils::isFileExist(*it))
                Utils::exitWithError("*** Error: Reads file " + *it + " does not exist. Check the config file!");

        // check read file is fasta or fastq
        if (Utils::isFastaFormat(Single_End_Reads))
            fasta_reads_flag = true;
        else if (Utils::isFastqFormat(Single_End_Reads))
            fasta_reads_flag = false;
        else
            Utils::exitWithError("*** Error: check read file format in the config file!");
    }
}


//=============================================================================
// Get bowtie cmd (bowtie or bowtie2): Default:bowtie2
//=============================================================================
std::string SigmaConfig::getBowtieCmd()
{
    // get key string
    std::string bowtie_cmd = "bowtie2";
    
    return bowtie_cmd;
}


//=============================================================================
// Get bowtie build cmd (bowtie-build or bowtie2-build) 
//=============================================================================
std::string SigmaConfig::getBowtieBuildCmd()
{
    // get key string
    std::string bowtie_build_cmd = "bowtie2-build";
    
    return bowtie_build_cmd;
}


//=============================================================================
// Get bowtie path 
//=============================================================================
std::string SigmaConfig::getBowtiePath()
{
    std::string bowtie_path;
    std::string bowtie_cmd = getBowtieCmd();

    // get key string
    std::string key = bowtie_directory_str;

    // if key is in map
    if (ParseConfig::keyExists(key)) 
    {
        std::string bowtie_directory = ParseConfig::getValueOfKey(key);
        // if provided path doesn't have "/" at last character
        if (*bowtie_directory.rbegin() != Utils::getPathSeparator()) 
            bowtie_path = bowtie_directory + Utils::getPathSeparator() + bowtie_cmd;
        else 
            bowtie_path = bowtie_directory + bowtie_cmd;

        // check executable exist
        if (!Utils::canExec(bowtie_path))
            Utils::exitWithError("*** Error: Cannot find " + bowtie_path + " Check the config file!");
    }
    else 
    {
        // bowtie path is just bowtie command
        bowtie_path = bowtie_cmd;

        // check the program exists
        std::string system_command = "which " + bowtie_path + " > /dev/null 2>&1";
        // system call
        const int res = system(system_command.c_str());
        if (res!=0)
            Utils::exitWithError("*** Error: Cannot find bowtie2 executable!");
    }

    return bowtie_path;
}


//=============================================================================
// Get bowtie-build path 
//=============================================================================
std::string SigmaConfig::getBowtieBuildPath()
{
    // variables
    std::string bowtie_build_path;
    std::string bowtie_build_cmd = getBowtieBuildCmd();

    // get key string
    std::string key = bowtie_directory_str;

    // if key is in map
    if (ParseConfig::keyExists(key)) 
    {
        std::string bowtie_directory = ParseConfig::getValueOfKey(key);

        // if provided path doesn't have "/" at last character
        if (*bowtie_directory.rbegin() != Utils::getPathSeparator()) 
            bowtie_build_path = bowtie_directory + Utils::getPathSeparator() + bowtie_build_cmd;
        else 
            bowtie_build_path = bowtie_directory + bowtie_build_cmd;

        // check executable exist
        if (!Utils::canExec(bowtie_build_path))
            Utils::exitWithError("*** Error: Cannot find " + bowtie_build_path + " Check the config file!");
    }
    else 
    {
        // bowtie build path is just bowtie build command
        bowtie_build_path = bowtie_build_cmd;

        // check the program exists
        std::string system_command = "which " + bowtie_build_path + " > /dev/null 2>&1";
        // system call
        const int res = system(system_command.c_str());
        if (res!=0)
            Utils::exitWithError("*** Error: Cannot find bowtie2 build executable!");
    }

    return bowtie_build_path;
}


//=============================================================================
// Get samtools path 
//=============================================================================
std::string SigmaConfig::getSamtoolsPath()
{
    std::string samtools_path;
    std::string samtools_cmd = "samtools";

    // get key string
    std::string key = samtools_directory_str;

    // if key is in map
    if (ParseConfig::keyExists(key)) 
    {
        std::string samtools_directory = ParseConfig::getValueOfKey(key);

        // if provided path doesn't have "/" at last character
        if (*samtools_directory.rbegin() != Utils::getPathSeparator()) 
            samtools_path = samtools_directory + Utils::getPathSeparator() + samtools_cmd;
        else 
            samtools_path = samtools_directory + samtools_cmd;

        // check executable exist
        if (!Utils::canExec(samtools_path))
            Utils::exitWithError("*** Error: Cannot find " + samtools_path + " Check the config file!");
    }
    else 
    {
        samtools_path = samtools_cmd;

        // check the program exists
        std::string system_command = "which " + samtools_path + " > /dev/null 2>&1";
        // system call
        const int res = system(system_command.c_str());
        if (res!=0)
            Utils::exitWithError("*** Error: Cannot samtools executable!");
    }
    
    return samtools_path;
}


//=============================================================================
// [Data_Info]getReferenceGenomeDirectory 
//=============================================================================
std::string SigmaConfig::getReferenceGenomeDirectory()
{
    // get key string
    std::string key = reference_genome_directory_str;
    std::string reference_genome_directory = "";
    
    // if key is in map
    if (ParseConfig::keyExists(key)) 
    {
        reference_genome_directory = ParseConfig::getValueOfKey(key);
        if (Utils::isDirectory(reference_genome_directory)) {
            return reference_genome_directory;
        }
        else  {
            Utils::exitWithError("*** Error: Reference_Genome_Directory does not exist. Check the config file!");
            return reference_genome_directory;
        }
    }
    // if key is not in map, then default or error
    else 
    {
        Utils::exitWithError("*** Error: Provide Reference_Genome_Directory in the config file.");
        return reference_genome_directory;
    }
}


//=============================================================================
// [Data_Info]getPairedEndReads1
//=============================================================================
std::string SigmaConfig::getPairedEndReads1()
{
    // get key string
    std::string key = paired_end_reads_1_str;
    return ParseConfig::getValueOfKey(key);
}


//=============================================================================
// [Data_Info]getPairedEndReads2
//=============================================================================
std::string SigmaConfig::getPairedEndReads2()
{
    // get key string
    std::string key = paired_end_reads_2_str;
    return ParseConfig::getValueOfKey(key);
}


//=============================================================================
// [Data_Info]getSingleEndReads
//=============================================================================
std::string SigmaConfig::getSingleEndReads()
{
    // get key string
    std::string key = single_end_reads_str;
    return ParseConfig::getValueOfKey(key);
}


//=============================================================================
// [Bowtie_Search]getBowtie2Preset
//=============================================================================
std::string SigmaConfig::getBowtie2Preset()
{
    // bowtie2_preset (default: sensitive)
    std::string bowtie2_preset = "sensitive";

    // get key string
    std::string key = bowtie2_preset_str;

    // if key is in map
    if (ParseConfig::keyExists(key)) 
        bowtie2_preset = ParseConfig::getValueOfKey(key);

    return bowtie2_preset;
}


//=============================================================================
// [Bowtie_Search]getMaximumMismatchCount
//=============================================================================
std::string SigmaConfig::getMaximumMismatchCount()
{
    // maximum_mismatch_count (default: 3)
    std::string maximum_mismatch_count = "3";

    // get key string
    std::string key = maximum_mismatch_count_str;

    // if key is in map
    if (ParseConfig::keyExists(key)) 
        maximum_mismatch_count = ParseConfig::getValueOfKey(key);

    return maximum_mismatch_count;
}


//=============================================================================
// [Bowtie_Search]getMinimumFragmentLength
//=============================================================================
std::string SigmaConfig::getMinimumFragmentLength()
{
    // minimum_fragment_length (default: 0)
    std::string minimum_fragment_length = "0";

    // get key string
    std::string key = minimum_fragment_length_str;

    // if key is in map
    if (ParseConfig::keyExists(key)) 
        minimum_fragment_length = ParseConfig::getValueOfKey(key);

    return minimum_fragment_length;
}


//=============================================================================
// [Bowtie_Search]getMaximumFragmentLength
//=============================================================================
std::string SigmaConfig::getMaximumFragmentLength()
{
    // minimum_fragment_length (default: 500)
    std::string maximum_fragment_length = "500";

	// get key string
	std::string key = maximum_fragment_length_str;

    // if key is in map
    if (ParseConfig::keyExists(key)) 
        maximum_fragment_length = ParseConfig::getValueOfKey(key);

	return maximum_fragment_length;
}


//=============================================================================
// [Bowtie_Search]getBowtieNumberThread
//=============================================================================
std::string SigmaConfig::getBowtieThreadsNumber()
{
    // bowtie threads number
    std::string bowtie_threads_number = "1";

	// get key string
	std::string key = bowtie_threads_number_str;

    // if key is in map
    if (ParseConfig::keyExists(key)) 
        bowtie_threads_number = ParseConfig::getValueOfKey(key);

	return bowtie_threads_number;
}


//=============================================================================
// [Model_Probability]getMismatchProbability
//=============================================================================
double SigmaConfig::getMismatchProbability()
{
    // mismatch probability (default = 0.05)
    double mismatch_probability = 0.05;

    // get key string
    std::string key = mismatch_probability_str;

    // if key is in map
    if (ParseConfig::keyExists(key)) 
        mismatch_probability = Utils::stringToDouble(ParseConfig::getValueOfKey(key));

    return mismatch_probability;
}


//=============================================================================
// [Model_Probability]getMinimumRelativeAbundance
//=============================================================================
double SigmaConfig::getMinimumRelativeAbundance()
{
    // minimum relative abundance (default = 0.01)
    double minimum_relative_abundance = 0.05;

    // get key string
    std::string key = minimum_relative_abundance_str;

    // if key is in map
    if (ParseConfig::keyExists(key)) 
        minimum_relative_abundance = Utils::stringToDouble(ParseConfig::getValueOfKey(key));

    return minimum_relative_abundance;
}


//=============================================================================
// [Statistics]Bootstrap_Iteration_Number
//=============================================================================
unsigned int SigmaConfig::getBootstrapIterationNumber()
{
    // get key string
    std::string key = bootstrap_iteration_number_str;
    std::string value_str = ParseConfig::getValueOfKey(key);
    unsigned int value = Utils::stringToUnsignedInt(value_str);
    return value;
}

