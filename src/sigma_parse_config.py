#!/usr/bin/env python

"""
sigma_parse_config.py

Created by Tae-Hyuk (Ted) Ahn on 02/22/2013.
Copyright (c) 2012 Tae-Hyuk Ahn (ORNL). Allrights reserved.
"""

## Import Python package modules
import sys

## Import sigma library
import utils

## global variables
config_path_str       = "config_path"
working_directory_str = "working_directory"

## string class
class ConfigString:
  
    # [Program_Info]
    program_info_str = "[Program_Info]"
    bowtie_directory_str = program_info_str + "Bowtie_Directory"
    samtools_directory_str = program_info_str + "Samtools_Directory"

    # [Data_Info]
    data_info_str = "[Data_Info]"
    reference_genome_directory_str = data_info_str + "Reference_Genome_Directory"
    paired_end_reads_1_str = data_info_str + "Paired_End_Reads_1"
    paired_end_reads_2_str = data_info_str + "Paired_End_Reads_2"
    single_end_reads_str = data_info_str + "Single_End_Reads"

    # [Bowtie_Search]
    bowtie_search_str = "[Bowtie_Search]"
    maximum_mismatch_count_str = bowtie_search_str + "Maximum_Mismatch_Count"
    minimum_fragment_length_str = bowtie_search_str + "Minimum_Fragment_Length"
    maximum_fragment_length_str = bowtie_search_str + "Maximum_Fragment_Length"
    bowtie_threads_number_str = bowtie_search_str + "Bowtie_Threads_Number"

    # [Model_Probability]
    model_probability_str = "[Model_Probability]"
    mismatch_probability_str = model_probability_str + "Mismatch_Probability"
    minimum_relative_abundance = model_probability_str + "Minimum_Relative_Abundance"

    # [Statistics]
    statistics_str = "[Statistics]"
    bootstrap_iteration_number_str = statistics_str + "Bootstrap_Iteration_Number"

    # [Genome_Reconstruction]
    genome_reconstruction_str = "[Genome_Reconstruction]"
    reconstruction_selection_str = genome_reconstruction_str + "Reconstruction_Selection"
    reconstruction_genome_name_str = genome_reconstruction_str + "Reconstruction_Genome_Name"
    reconstruction_cutoff_abundance_str = genome_reconstruction_str + "Reconstruction_Cutoff_Abundance"
    minumum_coverage_length_str = genome_reconstruction_str + "Minumum_Coverage_Length"
    minimum_average_coverage_depth_str = genome_reconstruction_str + "Minimum_Average_Coverage_Depth"

    # [Variants_Calling]
    reads_filtering_str = "[Variants_Calling]"
    filtering_genome_name_str = reads_filtering_str + "Filtering_Genome_Name"

## Parse config
def parse_config(option_map):

    # get config file
    config_path = option_map[config_path_str]

    # Save config values to dictionary
    config_map = {}    # initialize dictionay

    # Call parseConfigKeyValues
    config_map = parse_config_key_value(config_path)

    # return 
    return config_map


## parse config keys and values
def parse_config_key_value (filepath) :

    # get config file
    config_file = open(filepath)

    # initialize config map
    config_map = {}

    # read lines
    for line in config_file.readlines():
        line = line.strip()

        # do not consider comment lines
        if line.startswith("#"):
            continue
        # do not consider blank line
        elif line == "":
            continue
        # consider section name
        elif line.startswith("[") and line.endswith("]"):
            section_name = line
        # consider config key and value
        else:
            line_list = line.split("=")
            if (len(line_list) == 2):
                key = section_name.strip() + line_list[0].strip()
                val = line_list[1].strip()
                config_map[key] = val
            else:
                sys.stderr.write("\n** Check config file!that should be separted with # only once\n")
                utils.die("** Program exit!")

    return config_map
