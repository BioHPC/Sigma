#ifndef SIGMA_STAT_H
#define SIGMA_STAT_H

//=============================================================================
// sigma_stat.h
//   : This is the header file for methods and functions of SIGMA statistics.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified:
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include "utils.h"
#include "sigma_core.h"
#include "parse_config.h"
#include "sigma_config.h"


// Load Q matrix for bootstrapping **/
void getBootstrapQmatrixData(
        const std::string & qmatrix_filepath,                               // (in)
        std::vector<std::string> & qmatrix_comments_list,                   // (out)
        std::vector<std::string> & qmatrix_data_list);                      // (out)

// Write bootstrap Q matrix **/
void writeBootstrapQmatrix(
        const int & proc_id,                                                // (in)
        const std::string & bootstrap_qmatrix_filepath,                     // (in)
        const std::vector<std::string> & qmatrix_comments_list,             // (in)
        const std::vector<std::string> & qmatrix_data_list);                // (in)

// Get gvector data **/
void getGvectorData(
        const std::string & gvector_filepath,                               // (in)
        std::map<unsigned int, std::string> & gvector_genome_name_map,      // (out)
        std::map<unsigned int, double> & gvector_align_rate_map,            // (out)
        std::map<unsigned int, double> & gvector_pct_chance_map,            // (out)
        std::map<unsigned int, double> & gvector_pct_scaled_map);           // (out)

// Calculate bootstrap gvector statistics **/
void calculateBootstrapGvectorStat(
        const std::string & working_directory,                              // (in)
        const std::string & qmatrix_filename,                               // (in)
        const unsigned int bootstrap_iteration_number);                     // (in)

// Load Q matrix for jackknife **/
void getJackknifeQmatrixData(
        const std::string & qmatrix_filepath,                               // (in)
        std::vector<std::string> & qmatrix_comments_list,                   // (out)
        std::vector<unsigned int> & genome_index_list,                      // (out)
        std::map<unsigned int, std::string> & genome_name_map,              // (out)
        std::map<unsigned int, double> & align_rate_map,                    // (out)
        std::vector<std::string> & read_id_data,                            // (out)
        std::vector< std::vector<unsigned int> > & genome_index_data,       // (out)
        std::vector< std::vector<double> > & qvalue_data);                  // (out)

// Write jackknife Q matrix **/
void writeJackknifeQmatrix(
        const int & proc_id,                                                // (in)
        const int & task_id,                                                // (in)
        const std::string & jackknife_qmatrix_filepath,                     // (in)
        const std::vector<std::string> & qmatrix_comments_list,             // (in)
        const std::vector<unsigned int> & genome_index_list,                // (in)
        const std::map<unsigned int, std::string> & genome_name_map,        // (in)
        const std::map<unsigned int, double> & align_rate_map,              // (in)
        const std::vector<std::string> & read_id_data,                      // (in)
        const std::vector< std::vector<unsigned int> > & genome_index_data, // (in)
        const std::vector< std::vector<double> > & qvalue_data);            // (in)

// Get objective fuction value **/
void getObjFuncValue(
        const unsigned int & longest_read_length,                           // (in)
        const double & mismatch_probability,                                // (in)
        const std::vector< std::vector<unsigned int> > & genome_index_data, // (in)
        const std::vector< std::vector<double> > & qvalue_data,             // (in)
        const std::vector<double> & gvector_pct_chance_list,                // (in)
        const std::vector<double> & gvector_scaled_chance_list,             // (in)
        std::vector<double> & obj_value_list);                              // (out)
        
// Calculate jackknife gvector stat **/
void calculateJackknifeGvectorStat(
        const std::string & working_directory,                              // (in)
        const std::string & qmatrix_filename,                               // (in) 
        const std::string & qmatrix_filepath,                               // (in)
        const unsigned int & longest_read_length,                           // (in)
        const double & mismatch_probability,                                // (in)
        const std::vector<std::string> & qmatrix_comments_list,             // (in)
        const std::vector<unsigned int> & genome_index_list,                // (in)
        const std::map<unsigned int, std::string> & genome_name_map,        // (in)
        const std::map<unsigned int, double> & align_rate_map,              // (in)
        const std::vector<std::string> & read_id_data,                      // (in)
        const std::vector< std::vector<unsigned int> > & genome_index_data, // (in)
        const std::vector< std::vector<double> > & qvalue_data);            // (in)

#endif
