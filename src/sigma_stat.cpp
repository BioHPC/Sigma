//=============================================================================
// sigma_stat.cpp
//   : For methods and functions of SIGMA statistics.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified:
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include <math.h>
#include <numeric>
#include <iomanip>
#include "sigma_stat.h"


//=============================================================================
// Load Q matrix for bootstrapping
//=============================================================================
void getBootstrapQmatrixData(
        const std::string & qmatrix_filepath,               // (in)
        std::vector<std::string> & qmatrix_comments_list,   // (out)
        std::vector<std::string> & qmatrix_data_list)       // (out)
{
    // open Q-matrix -> input_file
    std::ifstream input_file (qmatrix_filepath.c_str());

    // get Q-matrix
    for (std::string line; getline (input_file, line); ) 
    {
        // data line
        if (line [0] == '*')
            qmatrix_data_list.push_back(line);
        // comments line
        else
            qmatrix_comments_list.push_back(line);
    }
}


//=============================================================================
// Load Q matrix  for jackknife 
//=============================================================================
void getJackknifeQmatrixData(
        const std::string & qmatrix_filepath,                           // (in)
        std::vector<std::string> & qmatrix_comments_list,               // (out)
        std::vector<unsigned int> & genome_index_list,                  // (out)
        std::map<unsigned int, std::string> & genome_name_map,          // (out)
        std::map<unsigned int, double> & align_rate_map,                // (out)
        std::vector<std::string> & read_id_data,                        // (out)
        std::vector< std::vector<unsigned int> > & genome_index_data,   // (out)
        std::vector< std::vector<double> > & qvalue_data)               // (out)
{
    // variables
    unsigned int genome_index;
    double qvalue, align_rate;
    std::string genome_name;

    // typedef vector of vector
    typedef std::vector< std::vector<unsigned int> > vector_of_vector_unsigned_int;
    typedef std::vector< std::vector<double> > vector_of_vector_double;

    // open Q-matrix -> input_file
    std::ifstream input_file (qmatrix_filepath.c_str());

    // get Q-matrix
    for (std::string line; getline (input_file, line); ) 
    {
        // line_stream
        std::istringstream line_stream(line);

        // starts with # and +
        if (line [0] == '#' || line [0] == '+') 
        {
            // save to qmatrix_comments_list
            qmatrix_comments_list.push_back(line);
        }
        // starts with @
        else if (line [0] == '@') 
        {
            // save to genome_index_list
            std::vector<std::string> fields_vector;
            for (std::string field; getline(line_stream, field, '\t'); ) { 
                fields_vector.push_back(field);
            }   

            // get each data
            genome_index = Utils::stringToUnsignedInt(fields_vector[1]);
            genome_name  = fields_vector[2];
            align_rate   = Utils::stringToDouble(fields_vector[3]);

            // save each data
            genome_index_list.push_back(genome_index);
			genome_name_map.insert(std::pair<unsigned int, std::string>(genome_index, genome_name));
            align_rate_map.insert(std::pair<unsigned int, double>(genome_index, align_rate));
        }
        // starts with @
        else if (line [0] == '*') 
        {
            // value_type for data
            genome_index_data.push_back(vector_of_vector_unsigned_int::value_type());
            qvalue_data.push_back(vector_of_vector_double::value_type());

            // loop token by tab-delimited and save to vector
            std::vector<std::string> fields_vector;
            for (std::string field; getline(line_stream, field, '\t'); ) {
                fields_vector.push_back(field);
            }

            // push_back read_id to read_id_data
            read_id_data.push_back(fields_vector[1]);

            // load from column 2 (data vector)
            for(unsigned int i = 2; i < fields_vector.size(); i++) {
                // field has pair "genome_index=q_value"
                std::string field=fields_vector[i];
                std::istringstream field_stream(field);

                // get genome_index and push to vector of vector
                std::string genome_index_string;
                getline(field_stream, genome_index_string, '=');
                genome_index = Utils::stringToUnsignedInt(genome_index_string);
                genome_index_data.back().push_back(genome_index);
           
                // get genome_index and push to vector of vector
                std::string qvalue_string;
                getline(field_stream, qvalue_string);
                qvalue = Utils::stringToDouble(qvalue_string);
                qvalue_data.back().push_back(qvalue);
            }
            fields_vector.clear();
        }
    }
}


//=============================================================================
// Write bootstrap Q matrix 
//=============================================================================
void writeBootstrapQmatrix(
        const int & proc_id,                                    // (in)
        const std::string & bootstrap_qmatrix_filepath,         // (in)
        const std::vector<std::string> & qmatrix_comments_list, // (in)
        const std::vector<std::string> & qmatrix_data_list)     // (in)
{
    // if file exists, then remove the file
    Utils::ifFileExistRemove(bootstrap_qmatrix_filepath);

    // new Q matrix -> output_file
    std::ofstream output_file (bootstrap_qmatrix_filepath.c_str());

    // print genome index info
    for (unsigned int i=0; i<qmatrix_comments_list.size(); i++) {
        output_file << qmatrix_comments_list[i] << std::endl;
    }

    // get Q matrix reads count
    unsigned int qmatrix_data_count = qmatrix_data_list.size();

    // print genome data with randomly seletected line
    for (unsigned int i=0; i<qmatrix_data_count; i++) {
        unsigned int random_index = rand() % qmatrix_data_count;
        output_file << qmatrix_data_list[random_index] << std::endl;
    }
}


//=============================================================================
// Write jackknife Q matrix 
//=============================================================================
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
        const std::vector< std::vector<double> > & qvalue_data)             // (in)
{
    // variables
    unsigned int genome_index, jk_genome_index;
    double qvalue, align_rate;
    std::string genome_name;

    // if file exists, then remove the file
    Utils::ifFileExistRemove(jackknife_qmatrix_filepath);

    // get jackknife genome index that will be excluded
    jk_genome_index = genome_index_list[task_id];

    // new Q matrix -> output_file
    std::ofstream output_file (jackknife_qmatrix_filepath.c_str());

    // print comments
    for (unsigned int i=0; i<qmatrix_comments_list.size(); i++) {
        output_file << qmatrix_comments_list[i] << std::endl;
    }

    // print genome index info
    for (unsigned int i=0; i<genome_index_list.size(); i++) {
        genome_index = genome_index_list[i];
        genome_name = genome_name_map.find(genome_index) -> second;
        align_rate = align_rate_map.find(genome_index) -> second;

        // if genome_index < jk_genome_index, then just print
        if (genome_index < jk_genome_index)
            output_file << "@\t" << genome_index << "\t" << genome_name << "\t" << align_rate << "\n";
        // if genome_index > jk_genome_index, then just print
        else if (genome_index > jk_genome_index)
            output_file << "@\t" << genome_index-1 << "\t" << genome_name << "\t" << align_rate << "\n";
    }

    // print q-matrix after removing specifig genome index
    for (unsigned int i=0; i<read_id_data.size(); i++) 
    {
        // print read_id
        if (!(genome_index_data[i].size() == 1 && genome_index_data[i][0] == jk_genome_index))
        {
            output_file << "*\t" << read_id_data[i];

            // loop columns of data
            for (unsigned int j=0; j<genome_index_data[i].size(); j++ ) 
            {
                // get genome index and qvalue of data
                genome_index = genome_index_data[i][j];
                qvalue = qvalue_data[i][j];

                // print without jackknife genome index data
                // if genome_index < jk_genome_index, then just print 
                if (genome_index < jk_genome_index) {
                    output_file << "\t" << genome_index << "=" << qvalue;
                }
                // if genome_index > jk_genome_index, then print genome_index-1
                else if (genome_index > jk_genome_index) {
                    output_file << "\t" << genome_index-1 << "=" << qvalue;
                }
            }
            // end line
            output_file << std::endl;
        }
    }
}


//=============================================================================
// Get gvector data 
//=============================================================================
void getGvectorData(
        const std::string & gvector_filepath,                           // (in)
        std::map<unsigned int, std::string> & gvector_genome_name_map,  // (out)
        std::map<unsigned int, double> & gvector_align_rate_map,        // (out)
        std::map<unsigned int, double> & gvector_pct_chance_map,        // (out)
        std::map<unsigned int, double> & gvector_pct_scaled_map)        // (out)
{
    // variables
    unsigned int genome_index;
    std::string genome_name;
    double align_rate, pct_chance, pct_scaled;

    // open gvector -> input_file
    std::ifstream input_file (gvector_filepath.c_str());

    // get genome list
    for (std::string line; getline (input_file, line); ) 
    {
        // line stringstream
        std::istringstream line_stream(line);

        // loop fields by tab-delimited
        std::vector<std::string> fields_vector;
        for (std::string field; getline(line_stream, field, '\t');) 
        {
            fields_vector.push_back(field);
        }

        // line with @ [GenomeIndex    GenomeName    PercentageAlignmentRate]
        if (line [0] == '@') 
        {
            genome_index = Utils::stringToUnsignedInt(fields_vector[1]);
            genome_name  = fields_vector[2];
            align_rate   = Utils::stringToDouble(fields_vector[3]);
			gvector_genome_name_map.insert(std::pair<unsigned int, std::string>(genome_index, genome_name));
            gvector_align_rate_map.insert(std::pair<unsigned int, double>(genome_index, align_rate));
        }
        // line with * [GenomeIndex    PercentageChance    ScaledPercentageChance]
        else if (line [0] == '*') 
        {
            genome_index = Utils::stringToUnsignedInt(fields_vector[1]);
            pct_chance   = Utils::stringToDouble(fields_vector[2]);
            pct_scaled   = Utils::stringToDouble(fields_vector[3]);
            gvector_pct_chance_map.insert(std::pair<unsigned int, double>(genome_index, pct_chance));
            gvector_pct_scaled_map.insert(std::pair<unsigned int, double>(genome_index, pct_scaled));
        }
    }
}


//=============================================================================
// Calculate bootstrap gvector 
//=============================================================================
void calculateBootstrapGvectorStat(
        const std::string & working_directory,          // (in)
        const std::string & qmatrix_filename,           // (in)
        const unsigned int bootstrap_iteration_number)  // (in)
{
    // variables
    std::string qmatrix_filebase;   // qmatrix filebase (sigma_config.qmatrix.txt -> sigma_config)
    std::string bootstrap_gvector_filename;     // new generated gvector filename
    std::string bootstrap_gvector_filepath;     // new generated gvector filepath
    std::string bootstrap_ipopt_log_filename;   // ipopt log filename
    std::string bootstrap_ipopt_log_filepath;   // ipopt log filepath
    std::string bootstrap_html_filename;        // html filename
    std::string bootstrap_html_filepath;        // html filepath

    unsigned int genome_index;
    double percentage_chance, percentage_scaled;
    std::string genome_name;

    // original gvector 
    qmatrix_filebase = Utils::getFilebase2(qmatrix_filename);
    std::string gvector_filename = qmatrix_filebase + ".gvector.txt";
    std::string gvector_filepath = working_directory + gvector_filename;

    // get gvector_genome_list
    std::map<unsigned int, std::string> gvector_genome_name_map;  // (out)
    std::map<unsigned int, double> gvector_align_rate_map, gvector_pct_chance_map, gvector_pct_scaled_map;
    if (Utils::isFileExist(gvector_filepath)) {
        getGvectorData(gvector_filepath,        // (in)
                       gvector_genome_name_map, // (out)
                       gvector_align_rate_map,  // (out)
                       gvector_pct_chance_map,  // (out)
                       gvector_pct_scaled_map); // (out)
    } 
    else {
        Utils::exitWithError("*** Error: gvector file does not exist in working directory!");
    }

    // bootstrap gvector all data
    std::map<unsigned int, std::vector<double> > bootstrap_pct_chance_data;
    std::map<unsigned int, std::vector<double> > bootstrap_pct_scaled_data;
    
    // loop bootstrap iteration
    for (unsigned int i=0; i<bootstrap_iteration_number; i++) 
    {
        // bootstrap gvector filename and filepath
        bootstrap_gvector_filename = qmatrix_filebase + ".bootstrap." + Utils::intToString(i) + ".gvector.txt";
        bootstrap_gvector_filepath = working_directory + bootstrap_gvector_filename;
        
        // load to bootstrap_pct_chance_data
        if (Utils::isFileExist(bootstrap_gvector_filepath)) 
        {
            // open bootstrap gvector -> input_file
            std::ifstream input_file (bootstrap_gvector_filepath.c_str());

            // get Q-matrix
            for (std::string line; getline (input_file, line); ) 
            {
                // data line
                if (line [0] == '*') 
                {
                    // line stringstream
                    std::istringstream line_stream(line);
            
                    // loop token by tab-delimited and save to vector
                    std::vector<std::string> fields_vector;
                    for (std::string field; getline(line_stream, field, '\t'); ) {
                        fields_vector.push_back(field);
                    }   

                    // get genome index and percentage chance from each bootstrap data
                    genome_index = Utils::stringToUnsignedInt(fields_vector[1]); // key: genome index
                    percentage_chance = Utils::stringToDouble(fields_vector[2]); // val: percentage chance
                    percentage_scaled = Utils::stringToDouble(fields_vector[3]); // val: scaled chance

                    // save to bootstrap_pct_chance_data
                    bootstrap_pct_chance_data[genome_index].push_back(percentage_chance);
                    bootstrap_pct_scaled_data[genome_index].push_back(percentage_scaled);
                }
            }
        }
    }

    // delete bootstrap gvector.txt and ipopt.txt html 
    for (unsigned int i=0; i<bootstrap_iteration_number; i++) 
    {
        // bootstrap gvector filename and filepath
        bootstrap_gvector_filename = qmatrix_filebase + ".bootstrap." + Utils::intToString(i) + ".gvector.txt";
        bootstrap_gvector_filepath = working_directory + bootstrap_gvector_filename;

        // ipopt log filename and filepath
        bootstrap_ipopt_log_filename = qmatrix_filebase + ".bootstrap." + Utils::intToString(i) + ".ipopt.txt";
        bootstrap_ipopt_log_filepath = working_directory + bootstrap_ipopt_log_filename;

        // ipopt log filename and filepath
        bootstrap_html_filename = qmatrix_filebase + ".bootstrap." + Utils::intToString(i) + ".html";
        bootstrap_html_filepath = working_directory + bootstrap_html_filename;

        // if file exists, then remove the file
        Utils::ifFileExistRemove(bootstrap_gvector_filepath);
        Utils::ifFileExistRemove(bootstrap_ipopt_log_filepath);
        Utils::ifFileExistRemove(bootstrap_html_filepath);
    }

    // bootstrap stat data
    std::map<unsigned int, std::vector<double> > bootstrap_stat_data;
    
    // iterate bootstrap_pct_chance_data
    typedef std::map<unsigned int, std::vector<double> >::iterator it_type;
    for(it_type iterator = bootstrap_pct_chance_data.begin(); iterator != bootstrap_pct_chance_data.end(); iterator++) 
    {
        // get key: genome index
        genome_index = iterator->first;

        // get value: percentage chance list
        std::vector<double> fields_vector = iterator->second;

        // size n
        unsigned int n = fields_vector.size();

        // sum, mean, and stdev
        double sum = std::accumulate(fields_vector.begin(), fields_vector.end(), 0.0);
        double mean = sum / n;
        double sq_sum = std::inner_product(fields_vector.begin(), fields_vector.end(), fields_vector.begin(), 0.0);
        double stdev = sqrt(sq_sum / n  - mean * mean);

        // 95% upper confidence bount 
        double upper_confidence_bound = mean + 1.96*( stdev / double(sqrt(n)) );
        // 95% lower confidence bount 
        double lower_confidence_bound = mean - 1.96*( stdev / double(sqrt(n)) );

        // save to bootstrap_stat_data
        bootstrap_stat_data[genome_index].push_back(mean);
        bootstrap_stat_data[genome_index].push_back(stdev);
        bootstrap_stat_data[genome_index].push_back(upper_confidence_bound);
        bootstrap_stat_data[genome_index].push_back(lower_confidence_bound);
    }

    // iterate bootstrap_pct_scaled_data
    typedef std::map<unsigned int, std::vector<double> >::iterator it_type;
    for(it_type iterator = bootstrap_pct_scaled_data.begin(); iterator != bootstrap_pct_scaled_data.end(); iterator++) 
    {
        // get key: genome index
        genome_index = iterator->first;

        // get value: percentage chance list
        std::vector<double> fields_vector = iterator->second;

        // size n
        unsigned int n = fields_vector.size();

        // sum, mean, and stdev
        double sum = std::accumulate(fields_vector.begin(), fields_vector.end(), 0.0);
        double mean = sum / n;
        double sq_sum = std::inner_product(fields_vector.begin(), fields_vector.end(), fields_vector.begin(), 0.0);
        double stdev = sqrt(sq_sum / n  - mean * mean);

        // 95% upper confidence bount 
        double upper_confidence_bound = mean + 1.96*( stdev / double(sqrt(n)) );
        // 95% lower confidence bount 
        double lower_confidence_bound = mean - 1.96*( stdev / double(sqrt(n)) );

        // save to bootstrap_stat_data
        bootstrap_stat_data[genome_index].push_back(mean);
        bootstrap_stat_data[genome_index].push_back(stdev);
        bootstrap_stat_data[genome_index].push_back(upper_confidence_bound);
        bootstrap_stat_data[genome_index].push_back(lower_confidence_bound);
    }

    // output write
    std::string bootstrap_pct_chance_path = working_directory + qmatrix_filebase + ".stat_bootstrap_percentage_scaled_estimation_data.txt";
    std::string bootstrap_pct_scaled_path = working_directory + qmatrix_filebase + ".stat_bootstrap_relative_abundance_estimation_data.txt";
    std::string bootstrap_gvector_stat_path = working_directory + qmatrix_filebase + ".stat_bootstrap.txt";

    // if file exists, then remove the file
    Utils::ifFileExistRemove(bootstrap_pct_chance_path);
    Utils::ifFileExistRemove(bootstrap_pct_scaled_path);
    Utils::ifFileExistRemove(bootstrap_gvector_stat_path);

    // open output files
    std::ofstream bootstrap_pct_chance_out (bootstrap_pct_chance_path.c_str());
    std::ofstream bootstrap_pct_scaled_out (bootstrap_pct_scaled_path.c_str());
    std::ofstream bootstrap_gvector_stat_out (bootstrap_gvector_stat_path.c_str());

    // print comments for bootstrap_pct_chance_out
    bootstrap_pct_chance_out << "#\t@\tGenomeIndex\tGenomeName\tPercentageScaledEstimation" << std::endl;
    bootstrap_pct_chance_out << "#\t*\tGenomeIndex";
    for (unsigned int i=0; i<bootstrap_iteration_number; i++) 
    {
        unsigned int run_num = i + 1;
        std::string run_str = "Run" + Utils::unsignedIntToString(run_num) + "_PercentageScaledEstimation";
        bootstrap_pct_chance_out << "\t" << run_str;
    }
    bootstrap_pct_chance_out << std::endl;

    // print comments for bootstrap_pct_scaled_out
    bootstrap_pct_scaled_out << "#\t@\tGenomeIndex\tGenomeName\tRelativeAbundanceEstimation" << std::endl;
    bootstrap_pct_scaled_out << "#\t*\tGenome_Index";
    for (unsigned int i=0; i<bootstrap_iteration_number; i++) 
    {
        unsigned int run_num = i + 1;
        std::string run_str = "Run" + Utils::unsignedIntToString(run_num) + "_RelativeAbundanceEstimation";
        bootstrap_pct_scaled_out << "\t" << run_str;
    }
    bootstrap_pct_scaled_out << std::endl;

    // print comments for bootstrap_gvector_stat_out
    bootstrap_gvector_stat_out << "#\t@\tGenomeIndex\tGenomeName\tPercentageScaledEstimation\tRelativeAbundanceEstimation" << std::endl;
    bootstrap_gvector_stat_out << "#\t*\tGenomeIndex\tPct_Average\tPct_STDEV\tPct_UpperConfidenceBound\tPct_LowerConfidenceBound";
    bootstrap_gvector_stat_out << "\tRel_Average\tRel_STDEV\tRel_UpperConfidenceBound\tRel_LowerConfidenceBound" << std::endl;

    // print gvector info
    for (std::map<unsigned int, std::string>::iterator iter = gvector_genome_name_map.begin(); iter != gvector_genome_name_map.end(); iter++) 
    {
        // get key: genome index
        genome_index = iter->first;
        // get value: genome name
        genome_name = iter->second;
        // percentage chance
        percentage_chance = gvector_pct_chance_map.find(genome_index) -> second;
        percentage_scaled = gvector_pct_scaled_map.find(genome_index) -> second;

        bootstrap_pct_chance_out << "@\t" << genome_index << "\t" << genome_name << "\t" << percentage_chance << std::endl;
        bootstrap_pct_scaled_out << "@\t" << genome_index << "\t" << genome_name << "\t" << percentage_scaled << std::endl;
        bootstrap_gvector_stat_out << "@\t" << genome_index << "\t" << genome_name << "\t" << percentage_chance << percentage_scaled << std::endl;
    }

    // print bootstrap data
    for (std::map<unsigned int, std::vector<double> >::iterator iter = bootstrap_pct_chance_data.begin(); iter != bootstrap_pct_chance_data.end(); iter++) 
    {
        // get key: genome index
        genome_index = iter->first;

        // get value: percentage chance list
        std::vector<double> fields_vector = iter->second;
        
        // print
        bootstrap_pct_chance_out << "*\t" << genome_index;
        for (unsigned int i=0; i<fields_vector.size(); i++) {
            bootstrap_pct_chance_out << "\t" << fields_vector[i];
        }
        bootstrap_pct_chance_out << std::endl;
    }

    // print bootstrap data
    for (std::map<unsigned int, std::vector<double> >::iterator iter = bootstrap_pct_scaled_data.begin(); iter != bootstrap_pct_scaled_data.end(); iter++) 
    {
        // get key: genome index
        genome_index = iter->first;

        // get value: percentage chance list
        std::vector<double> fields_vector = iter->second;
        
        // print
        bootstrap_pct_scaled_out << "*\t" << genome_index;
        for (unsigned int i=0; i<fields_vector.size(); i++) {
            bootstrap_pct_scaled_out << "\t" << fields_vector[i];
        }
        bootstrap_pct_scaled_out << std::endl;
    }

    // print bootstrap data
    for (std::map<unsigned int, std::vector<double> >::iterator iter = bootstrap_stat_data.begin(); iter != bootstrap_stat_data.end(); iter++) 
    {
        // get key: genome index
        genome_index = iter->first;

        // get value: percentage chance list
        std::vector<double> fields_vector = iter->second;
        
        // print
        bootstrap_gvector_stat_out << "*\t" << genome_index;
        for (unsigned int i=0; i<fields_vector.size(); i++) {
            bootstrap_gvector_stat_out << "\t" << fields_vector[i];
        }
        bootstrap_gvector_stat_out << std::endl;
    }
}


//=============================================================================
// Get objective function value 
//=============================================================================
void getObjFuncValue(
        const unsigned int & longest_read_length,                           // (in)
        const double & mismatch_probability,                                // (in)
        const std::vector< std::vector<unsigned int> > & genome_index_data, // (in)
        const std::vector< std::vector<double> > & qvalue_data,             // (in)
        const std::vector<double> & gvector_pct_chance_list,                // (in)
        const std::vector<double> & gvector_scaled_chance_list,             // (in)
        std::vector<double> & obj_value_list)                               // (out)
{
    // variables
    unsigned int column_number, genome_index;
    double qvalue, obj_value_pct_chance, obj_value_scaled_chance;
    double obj_value_pct_chance_one, obj_value_scaled_chance_one;
    double no_match_qvalue;

    // initialize
    obj_value_pct_chance = 0.0;
    obj_value_scaled_chance = 0.0;

    // to skip log(0) for no match
    // Ex) 100bp nomatch (100 mismatches) --> (0.05)^100 * (0.95)^0 --> scaled --> (0.05)^100 * (0.95)^(0-100)
    no_match_qvalue = pow(mismatch_probability/(1.0-mismatch_probability),longest_read_length);

    // loop number of reads
    for (unsigned int i=0; i<qvalue_data.size(); i++) 
    {
        // save partial objective value for each row
        obj_value_pct_chance_one = 0.0;
        obj_value_scaled_chance_one = 0.0;

        // column number is the number of genomes that matched to the read
        column_number = qvalue_data[i].size();
        
        // loop columns of the row
        for (unsigned int j=0; j<column_number; j++) 
        {
            // get genome index and qvalue from data
            genome_index = genome_index_data[i][j];
            qvalue = qvalue_data[i][j];

            // one row of object value inside log is sum{Q[i][j]*X[j]} for each j 
            obj_value_pct_chance_one += qvalue*((double) gvector_pct_chance_list[genome_index])*0.01 ;
            obj_value_scaled_chance_one += qvalue*((double) gvector_scaled_chance_list[genome_index])*0.01;
        }

        // skip log(0) for no match
        if (obj_value_pct_chance_one < no_match_qvalue)
            obj_value_pct_chance_one = no_match_qvalue;

        if (obj_value_scaled_chance_one < no_match_qvalue)
            obj_value_scaled_chance_one = no_match_qvalue;

        // object value = sum{log[obj_value_one]} for each i
        obj_value_pct_chance += log(obj_value_pct_chance_one);
        obj_value_scaled_chance += log(obj_value_scaled_chance_one);
    }

    // save final objective value
    obj_value_list.push_back(obj_value_pct_chance);
    obj_value_list.push_back(obj_value_scaled_chance);
}
        

//=============================================================================
// Calculate jackknife gvector stat 
//=============================================================================
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
        const std::vector< std::vector<double> > & qvalue_data)             // (in)
{
    // variables
    unsigned int genome_index, jk_genome_index;
    double pct_chance, scaled_chance, jk_pct_chance, jk_scaled_chance, pct_chance_log_ratio, scaled_chance_log_ratio;
    double pct_chance_obj_value, jk_pct_chance_obj_value;
    double scaled_chance_obj_value, jk_scaled_chance_obj_value;
    std::string jk_gvector_filename, jk_gvector_filepath, genome_name, jk_ipopt_log_filename, jk_ipopt_log_filepath;
    std::string jk_html_filename, jk_html_filepath;

    // original gvector 
    std::string qmatrix_filebase = Utils::getFilebase2(qmatrix_filename);
    std::string gvector_filename = qmatrix_filebase + ".gvector.txt";
    std::string gvector_filepath = working_directory + gvector_filename;

    // get gvector_genome_list
    std::map<unsigned int, std::string> gvector_genome_name_map;  // (out)
    std::map<unsigned int, double> gvector_align_rate_map, gvector_pct_chance_map, gvector_pct_scaled_map;

    // jackknife objective value list initialize
    std::vector<double> gvector_pct_chance_list, gvector_scaled_chance_list, obj_value_list;
    std::vector<double> jk_pct_chance_obj_value_list, jk_scaled_chance_obj_value_list;

    // log ratio test value list
    std::vector<double> pct_chance_log_ratio_list;
    std::vector<double> scaled_chance_log_ratio_list;

    // if g-vector file exist, then get g-vector data
    if (Utils::isFileExist(gvector_filepath)) 
    {
        // get g-vector data
        getGvectorData(gvector_filepath,        // (in)
                       gvector_genome_name_map, // (out)
                       gvector_align_rate_map,  // (out)
                       gvector_pct_chance_map,  // (out)
                       gvector_pct_scaled_map); // (out)
    } 
    else {
        Utils::exitWithError("*** Error: gvector file " + gvector_filepath + " does not exist in working directory!");
    }

    // get gvector_pct_chance_list and gvector_scaled_chance_list
    for (unsigned int i=0; i<genome_index_list.size(); i++) 
    {
        genome_index = genome_index_list[i];

        // check key exist
        if (gvector_pct_chance_map.find(genome_index) != gvector_pct_chance_map.end() ) 
            pct_chance = gvector_pct_chance_map.find(genome_index) -> second;
        else
            pct_chance = 0.0;

        // check key exist
        if (gvector_pct_scaled_map.find(genome_index) != gvector_pct_scaled_map.end() ) 
            scaled_chance = gvector_pct_scaled_map.find(genome_index) -> second;
        else
            scaled_chance = 0.0;
      
        // push to the list
        gvector_pct_chance_list.push_back(pct_chance);
        gvector_scaled_chance_list.push_back(scaled_chance);
    }

    // get objective function value
    getObjFuncValue(longest_read_length,        // (in)
                    mismatch_probability,       // (in)
                    genome_index_data,          // (in) 
                    qvalue_data,                // (in)
                    gvector_pct_chance_list,    // (in)
                    gvector_scaled_chance_list, // (in)
                    obj_value_list);            // (out)

    // get pct_chance_obj_value and scaled_chance_obj_value by original gvector
    pct_chance_obj_value = obj_value_list[0];
    scaled_chance_obj_value = obj_value_list[1];

    // loop genome index to get new g-vector 
    for (unsigned int i=0; i<genome_index_list.size(); i++) 
    {
        // get jackknife genome index
        jk_genome_index = genome_index_list[i];

        // bootstrap gvector filename and filepath
        jk_gvector_filename = qmatrix_filebase + ".jackknife." + Utils::intToString(jk_genome_index) + ".gvector.txt";
        jk_gvector_filepath = working_directory + jk_gvector_filename;

        // get gvector_genome_list
        std::map<unsigned int, std::string> jk_gvector_genome_name_map;  // (out)
        std::map<unsigned int, double> jk_gvector_align_rate_map, jk_gvector_pct_chance_map, jk_gvector_pct_scaled_map;

        // if jackknife g-vector file exists, then get jackknife g-vector data
        if (Utils::isFileExist(jk_gvector_filepath)) 
        {
            // get g-vector data
            getGvectorData(jk_gvector_filepath,        // (in)
                           jk_gvector_genome_name_map, // (out)
                           jk_gvector_align_rate_map,  // (out)
                           jk_gvector_pct_chance_map,  // (out)
                           jk_gvector_pct_scaled_map); // (out)
        } 
        else {
            Utils::exitWithError("*** Error: gvector file " + jk_gvector_filepath + " does not exist in working directory!");
        }

        // list initialize
        std::vector<double> jk_gvector_pct_chance_list, jk_gvector_scaled_chance_list;

        // get jk_gvector_pct_chance_list and jk_gvector_scaled_chance_list
        for (unsigned int j=0; j<genome_index_list.size(); j++) 
        {
            // get revised jackknife genome index
            genome_index = genome_index_list[j];
            

            // if genome_index == jk_genome_index, then push_back 0.0 to the list
            if (genome_index == jk_genome_index) {
                jk_gvector_pct_chance_list.push_back(0.0);
                jk_gvector_scaled_chance_list.push_back(0.0);
            }
            // if not, then get the chance value by genome name
            else {
                // get genome name by genome index
                genome_name = genome_name_map.find(genome_index) -> second;

                // find revised jackknife genome index (key) by genome name (value)
                std::map<unsigned int, std::string>::const_iterator it;
                unsigned int target_genome_index = -1;

                // to find the key using a value, will require a linear search of the map
                for (it=jk_gvector_genome_name_map.begin(); it!=jk_gvector_genome_name_map.end(); ++it)
                {
                    if (it->second == genome_name) {
                        target_genome_index = it->first;
                        break;
                    }
                }

                // get percentage chance and scaled chance from revised jackknife g-vector
                jk_pct_chance = jk_gvector_pct_chance_map.find(target_genome_index) -> second;
                jk_scaled_chance = jk_gvector_pct_scaled_map.find(target_genome_index) -> second;

                // push back to list
                jk_gvector_pct_chance_list.push_back(jk_pct_chance);
                jk_gvector_scaled_chance_list.push_back(jk_scaled_chance);
            }
        }

        // jk_obj_value_list
        std::vector<double> jk_obj_value_list;

        // get objective function value
        getObjFuncValue(longest_read_length,           // (in)
                        mismatch_probability,          // (in)
                        genome_index_data,             // (in) 
                        qvalue_data,                   // (in)
                        jk_gvector_pct_chance_list,    // (in)
                        jk_gvector_scaled_chance_list, // (in)
                        jk_obj_value_list);            // (out)

        // get each jackknife objective value (0:percentage chance, 1:scaled chance)
        jk_pct_chance_obj_value_list.push_back(jk_obj_value_list[0]);
        jk_scaled_chance_obj_value_list.push_back(jk_obj_value_list[1]);

        // memory clear
        jk_gvector_genome_name_map.clear();
        jk_gvector_align_rate_map.clear();
        jk_gvector_pct_chance_map.clear();
        jk_gvector_pct_scaled_map.clear();
        jk_gvector_pct_chance_list.clear();
        jk_gvector_scaled_chance_list.clear();
    }

    // delete jackknife gvector.txt and ipopt.txt
    for (unsigned int i=0; i<genome_index_list.size(); i++) 
    {
        // get jackknife genome index
        jk_genome_index = genome_index_list[i];

        // jackknife gvector filename and filepath
        jk_gvector_filename = qmatrix_filebase + ".jackknife." + Utils::intToString(jk_genome_index) + ".gvector.txt";
        jk_gvector_filepath = working_directory + jk_gvector_filename;

        // ipopt log filename and filepath
        jk_ipopt_log_filename = qmatrix_filebase + ".jackknife." + Utils::intToString(jk_genome_index) + ".ipopt.txt";
        jk_ipopt_log_filepath = working_directory + jk_ipopt_log_filename;

        // html filename and filepath
        jk_html_filename = qmatrix_filebase + ".jackknife." + Utils::intToString(jk_genome_index) + ".html";
        jk_html_filepath = working_directory + jk_html_filename;

        // if file exists, then remove the file
        Utils::ifFileExistRemove(jk_gvector_filepath);
        Utils::ifFileExistRemove(jk_ipopt_log_filepath);
        Utils::ifFileExistRemove(jk_html_filepath);
    }

    // likelihood ratio test value list
    for (unsigned int i=0; i<genome_index_list.size(); i++) 
    {
        // get log ratio test value (-2*log(L0/L1))
        pct_chance_log_ratio = fabs(-2.0*(jk_pct_chance_obj_value_list[i] - pct_chance_obj_value));
        scaled_chance_log_ratio = fabs(-2.0*(jk_scaled_chance_obj_value_list[i] - scaled_chance_obj_value));

        // save to list
        pct_chance_log_ratio_list.push_back(pct_chance_log_ratio);
        scaled_chance_log_ratio_list.push_back(scaled_chance_log_ratio);
    }

    // output write
    std::string jk_pct_chance_path = working_directory + qmatrix_filebase + ".stat_jackknife_percentage_scaled_estimation.txt";
    std::string jk_scaled_chance_path = working_directory + qmatrix_filebase + ".stat_jackknife_relative_abundance_estimation.txt";

    // if file exists, then remove the file
    Utils::ifFileExistRemove(jk_pct_chance_path);
    Utils::ifFileExistRemove(jk_scaled_chance_path);

    // open output file
    std::ofstream jk_pct_chance_out(jk_pct_chance_path.c_str());
    std::ofstream jk_scaled_chance_out(jk_scaled_chance_path.c_str());
    std::streamsize jk_pct_chance_precision = jk_pct_chance_out.precision(20);
    std::streamsize jk_scaled_chance_precision = jk_pct_chance_out.precision(20);

    // print comments 
    jk_pct_chance_out << "#\t@\tGenomeIndex\tGenomeName\tPercentageScaledEstimation" << std::endl;
    jk_pct_chance_out << "#\t*\tGenomeIndex\tObjectiveFunctionValue[log]\tGoutObjectiveFunctionValue[log]\tLikelihoodRatio" << std::endl;

    // print comments 
    jk_scaled_chance_out << "#\t@\tGenomeIndex\tGenomeName\tRelativeAbundanceEstimation" << std::endl;
    jk_scaled_chance_out << "#\t*\tGenomeIndex\tObjectiveFunctionValue[log]\tGoutObjectiveFunctionValue[log]\tLikelihoodRatio" << std::endl;

    // print gvector info
    for (unsigned int i=0; i<genome_index_list.size(); i++) 
    {
        // get revised jackknife genome index
        genome_index = genome_index_list[i];
            
        // get value: genome name
        genome_name = genome_name_map.find(genome_index) -> second;

        // percentage chance
        pct_chance = gvector_pct_chance_map.find(genome_index) -> second;
        scaled_chance = gvector_pct_scaled_map.find(genome_index) -> second;

        // print to out
        jk_pct_chance_out << "@\t" << genome_index << "\t" << genome_name << "\t" << pct_chance << std::endl;
        jk_scaled_chance_out << "@\t" << genome_index << "\t" << genome_name << "\t" << scaled_chance << std::endl;
    }

    // print likelihood ratio test
    for (unsigned int i=0; i<genome_index_list.size(); i++) 
    {
        // get revised jackknife genome index
        genome_index = genome_index_list[i];
            
        // get pct_chance_log_ratio_list and scaled_chance_log_ratio_list
        jk_pct_chance_obj_value = jk_pct_chance_obj_value_list[genome_index];
        jk_scaled_chance_obj_value = jk_scaled_chance_obj_value_list[genome_index];
        pct_chance_log_ratio = pct_chance_log_ratio_list[genome_index];
        scaled_chance_log_ratio = scaled_chance_log_ratio_list[genome_index];
        
        // write
        jk_pct_chance_out << "*\t" << genome_index << "\t" << pct_chance_obj_value 
                          << "\t" << jk_pct_chance_obj_value << "\t" << pct_chance_log_ratio << std::endl;
        jk_scaled_chance_out << "*\t" << genome_index << "\t" << scaled_chance_obj_value 
                             << "\t" << jk_scaled_chance_obj_value << "\t" << scaled_chance_log_ratio << std::endl;
    }

    // close
    jk_pct_chance_out.precision(jk_pct_chance_precision);
    jk_scaled_chance_out.precision(jk_scaled_chance_precision);
    jk_pct_chance_out.close();
    jk_scaled_chance_out.close();

}
