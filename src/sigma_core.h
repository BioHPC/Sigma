#ifndef SIGMA_CORE_H
#define SIGMA_CORE_H

//=============================================================================
// sigma_core.h
//   : This is the header file for methods and functions of SIGMA.
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


// Load config file and get bowtie-build option 
void getBowtieBuildOption(
        const std::string & config_path,                                    // (in)
        std::string & reference_genome_directory,                           // (out)
        std::string & bowtie_build_option_cmd);                             // (out)

// Check genome index is already build in previous 
bool isBowtieIndexExist(
        const std::string & genome_index_base,                              // (in)
        const std::vector<std::string> & genome_fastas_path_sublist);       // (in)

// Load config file and get bowtie option **/
void getBowtieOption(
        const std::string & config_path,                                    // (in)
        std::string & reference_genome_directory,                           // (out)
        std::string & bowtie_option_cmd,                                    // (out)
        std::string & read_option_cmd,                                      // (out)
        std::string & samtools_cmd);                                        // (out)

// Search genome directory list and save to vector **/
void searchGenomeDirectoryList(
        const std::string & working_directory,                              // (in)
        std::string & reference_genome_directory,                           // (in)
        std::vector<std::string> & genome_directory_list,                   // (out)
        std::vector< std::vector<std::string> > & genome_fasta_path_list,   // (out)
        std::vector<std::string> & genome_index_base_list,                  // (out)
        std::string & output_parent_genome_directory,                       // (out)
        std::vector<std::string> & output_genome_directory_list);           // (out)

// Search remaining genome directory list and save to vector **/
void remainGenomeDirectoryList(
        std::string & reference_genome_directory,                           // (in)
        std::vector<int> & remain_genome_index_list);                       // (out)

// get longest read length and number of reads for each file
void getReadsInfo(const std::string & read_file,                            // (in)
                  const bool & fasta_reads_flag,                            // (in)
                  unsigned int & read_length,                               // (out)
                  unsigned int & number_of_reads);                          // (out)

// Count total number of reads
void getReadsLengthAndCount(const std::string & config_path,                // (in)
                            unsigned int & total_number_of_reads,           // (out)
                            unsigned int & longest_read_length);            // (out)

// get alignment rate from bowtie alignment log file
double getAlignmentRate(const std::string & alignment_log_path,             // (in)
                        const std::string & alignment_bam_path);            // (in)

// return common read ID for the paired_end
std::string getReadID(std::string & read_id_raw,                            // (in)
                      bool & paired_end_reads_flag,                         // (in)
                      bool & fasta_reads_flag);                             // (in)

// get values of SAM fields from SAM output line
void getValuesOfSAMLine(const std::string & line,                           // (in)
                        std::string & read_id,                              // (out)
                        unsigned int & read_flag,                           // (out)
                        unsigned int & read_seq_length,                     // (out)
                        unsigned int & read_mismatch_count);                // (out)

// Generate Q matrix
void generateQMatrix(const std::string & config_path,                       // (in)
                     const std::string & working_directory);                // (in)

// Write Sigma output to HTML format **/
void writeHTMLOutput(const std::string & gvector_filepath);                 // (in)

#endif
