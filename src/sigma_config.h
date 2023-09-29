#ifndef SIGMA_CONFIG_H
#define SIGMA_CONFIG_H

//=============================================================================
// sigma_config.h
//   : sigma_config.h is the header file for handling string and getting
//     value of sigma config file.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified:
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include "utils.h"
#include "parse_config.h"


//=============================================================================
// Class: SigmaConfig
//=============================================================================
class SigmaConfig : public ParseConfig
{
public:
    // Default constructor
    SigmaConfig(const std::string &config_path);

    // Default desstructor
    ~SigmaConfig();

    // Get reads type 
    void getReadsType(bool & paired_end_reads_flag, bool & fasta_reads_flag);

    // [Program_Info] 
    std::string getBowtieCmd();
    std::string getBowtiePath();
    std::string getBowtieBuildCmd();
    std::string getBowtieBuildPath();
    std::string getSamtoolsPath();

    // [Data_Info]
    std::string getReferenceGenomeDirectory();
    std::string getPairedEndReads1();
    std::string getPairedEndReads2();
    std::string getSingleEndReads();

    // [Bowtie_Search]
    std::string getBowtie2Preset();
    std::string getMaximumMismatchCount();
    std::string getMinimumFragmentLength();
    std::string getMaximumFragmentLength();
    std::string getBowtieThreadsNumber();

    // [Model_Probability] 
    double getMismatchProbability();
    double getMinimumRelativeAbundance();

    // [Statistics] 
    unsigned int getBootstrapIterationNumber();


private:
    // config config_path
    static std::string config_path;

    // [Program_Info]
    static const std::string program_info_str;
    static const std::string bowtie_directory_str;
    static const std::string samtools_directory_str;

    // [Data_Info]
    static const std::string data_info_str;
    static const std::string reference_genome_directory_str;
    static const std::string paired_end_reads_1_str;
    static const std::string paired_end_reads_2_str;
    static const std::string single_end_reads_str;

    // [Bowtie_Search]
    static const std::string bowtie_search_str;
    static const std::string bowtie2_preset_str;
    static const std::string maximum_mismatch_count_str;
    static const std::string minimum_fragment_length_str;
    static const std::string maximum_fragment_length_str;
    static const std::string bowtie_threads_number_str;

    // [Model_Probability]
    static const std::string model_probability_str;
    static const std::string mismatch_probability_str;
    static const std::string minimum_relative_abundance_str;

    // [Statistics]
    static const std::string statistics_str;
    static const std::string bootstrap_iteration_number_str;

    // [Program_Info]
    std::string bowtie_selection_def;

};

#endif
