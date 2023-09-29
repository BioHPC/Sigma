#!/usr/bin/env python

"""
sigma_core.py

sigma_core.py include common modules for sigma

Created by Tae-Hyuk (Ted) Ahn on 03/01/2013.
Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). Allrights reserved.
"""

## Import standard modules
import sys, warnings, os, re
import time
from subprocess import Popen, PIPE, check_call, STDOUT
import random
from math import sqrt

## Import sigma library
import utils
import sigma_parse_config


## global variables
config_path_str       = "config_path"
working_directory_str = "working_directory"


## Import classes
ConfigString = sigma_parse_config.ConfigString


## Search genome fasta files
def search_genome_fasta_path(config_map):

    # retrieve options
    reference_genome_directory = utils.set_directory_path(config_map[ConfigString.reference_genome_directory_str])

    # genome data structure
    genome_directory_list = []
    genome_fasta_path_list = []
    fasta_count = 0 

    # loop the genome directory
    dirs = [d for d in os.listdir(reference_genome_directory) if os.path.isdir(os.path.join(reference_genome_directory, d))]
    for genome_enum, genome_name in enumerate(dirs):

        # get genome directory
        genome_directory = reference_genome_directory + genome_name
        genome_directory_list.append(genome_directory)
        genome_fasta_path_list.append([])

        # loop fasta files
        for filename in os.listdir(genome_directory):

            # check fasta file
            if utils.check_fasta_format(filename):
                fasta_path = utils.set_directory_path(genome_directory) + filename
                genome_fasta_path_list[genome_enum].append(fasta_path)
                fasta_count += 1

    # check no orgaism list
    if len(genome_directory_list) == 0:
        sys.stderr.write("\n** Cannot search genome directory. Check the reference genome directory.\n")
        utils.die("** Program exit!")

    # check no fasta file
    if fasta_count == 0:
        sys.stderr.write("\n** Cannot find genome fasta file. Check the reference genome directory.\n")
        utils.die("** Program exit!")

    return (genome_directory_list, genome_fasta_path_list)


## Get genome index base with genome_directory_list
def get_genome_index_base(genome_directory_list):

    # genome_index_base_list
    genome_index_base_list = []    # initialize dictionary

    # for loop genome_index_base_list
    for genome_directory in genome_directory_list:

        # retrieve genome_directory_child_path
        genome_name = os.path.basename(genome_directory)

        # set bowtie_index_count as 0
        bowtie_index_count = 0

        # get genome_index_base
        genome_index_base = utils.set_directory_path(genome_directory) + genome_name

        # loop genome_index_path
        for filename in os.listdir(genome_directory):
            if utils.check_bowtie_index_format(filename):
                bowtie_index_count += 1

        # if bowtie index exist, then save to list
        if bowtie_index_count > 0:
            genome_index_base_list.append(genome_index_base)

    # check genome index == 0
    if len(genome_index_base_list) == 0:
        sys.stderr.write("\n** Cannot open bowtie2 index files. Check the Reference_Genome_Directory in config file.\n")
        utils.die("** Program exit!")

    return genome_index_base_list


## get output genome directory list
def get_output_genome_directory_list(option_map, config_map, genome_directory_list):

    # output genome directory list
    output_genome_directory_list = []

    # get working directory
    working_directory = option_map[working_directory_str]

    # get reference_genome_directory
    reference_genome_directory = utils.set_directory_path(config_map[ConfigString.reference_genome_directory_str])
    parent_genome_directory_name = "sigma_alignments_output"

    # get parent genome directory
    parent_genome_directory = working_directory + parent_genome_directory_name

    # mkdir parent_genome_directory
    if not os.path.isdir(parent_genome_directory):
        check_call(["mkdir", parent_genome_directory], stdout = PIPE, stderr = sys.stderr)

    # for loop genome_index_base_list
    for genome_directory in genome_directory_list:

        # retrieve genome_directory_child_path
        genome_name = os.path.basename(genome_directory)
        child_genome_directory = utils.set_directory_path(parent_genome_directory) + genome_name
        output_genome_directory_list.append(child_genome_directory)

        # mkdir child_genome_directory
        if not os.path.isdir(child_genome_directory):
            check_call(["mkdir", child_genome_directory], stdout = PIPE, stderr = sys.stderr)

    return output_genome_directory_list


## get output genome directory list (without genome directory)
def get_output_genome_directory_list_wo(option_map, config_map):

    # output genome directory list
    output_genome_directory_list = []

    # get working directory
    working_directory = option_map[working_directory_str]

    # get reference_genome_directory
    reference_genome_directory = utils.set_directory_path(config_map[ConfigString.reference_genome_directory_str])
    parent_genome_directory_name = os.path.basename(os.path.normpath(reference_genome_directory))

    # get parent genome directory
    parent_genome_directory = working_directory + parent_genome_directory_name

    # mkdir parent_genome_directory
    if not os.path.isdir(parent_genome_directory):
        sys.stderr.write("\n** Cannot open output directory. Check the Reference_Genome_Directory in config file exists in the working directory.\n")
        utils.die("** Program exit!")

    # loop the genome directory
    dirs = [d for d in os.listdir(parent_genome_directory) if os.path.isdir(os.path.join(parent_genome_directory, d))]
    for genome_enum, genome_name in enumerate(dirs):

        # get genome directory
        child_genome_directory = utils.set_directory_path(parent_genome_directory) + genome_name

        # loop fasta files
        for filename in os.listdir(child_genome_directory):

            # check bam file exists
            if utils.check_bam_format(filename):

                # append the genome directory into the list only if bam file exists
                output_genome_directory_list.append(child_genome_directory)

    return output_genome_directory_list


## Get fasta or fastq file(s) list in working dir
def get_fastaq_file_list(working_directory):

    # define sipros file extension 
    file_list = []

    # working directory
    if os.path.exists(working_directory):
        for file_name in os.listdir(working_directory):

            # check the fata format
            if utils.check_fasta_format(file_name) or utils.check_fastq_format(file_name):
                file_path_name = working_directory + file_name
                file_list.append(file_path_name)

        if len(file_list) == 0:
            sys.stderr.write("\n** Cannot open fasta format file(s).\n")
            utils.die("** Program exit!")
        file_list = sorted(file_list)

    else:
        sys.stderr.write("\n** Cannot open working_directory %s.\n" % working_directory)
        utils.die("** Program exit!")

    return file_list


## get reads type (paired_end or unpaired_end), (fasta or fastq)
def get_reads_type(config_map):

    # get reads
    paired_end_reads_flag = False
    fasta_reads_flag = False

    # paired_reads_1 and paired_reads_2
    if (ConfigString.paired_end_reads_1_str in config_map) and (ConfigString.paired_end_reads_2_str in config_map):
        paired_end_reads_1 = config_map[ConfigString.paired_end_reads_1_str]
        paired_end_reads_2 = config_map[ConfigString.paired_end_reads_2_str]
        paired_end_reads_flag = True

        # check fasta or fastq
        if utils.check_fasta_format(paired_end_reads_1):
            fasta_reads_flag = True
        else:
            fasta_reads_flag = False

    # unpaired_reads
    elif (ConfigString.unpaired_reads_str in config_map):
        single_end_reads = config_map[ConfigString.single_end_reads_str]
        paired_end_reads_flag = False

        # check fasta or fastq
        if utils.check_fasta_format(single_end_reads):
            fasta_reads_flag = True
        else:
            fasta_reads_flag = False
    else:
        sys.stderr.write("\n** Cannot open reads files. Check the config file.\n")
        utils.die("** Program exit!")

    return (paired_end_reads_flag, fasta_reads_flag)


## check aligned reads
def check_aligned_reads(filename):

    check_align = False
    if utils.check_path_exist(filename):
        check_align = True
        return check_align
    else:
        sys.stderr.write("\n** Cannot open %s.\n" %(filename))
        utils.die("** Program exit!")
    

## get_samout_path
def get_samout_path(config_path, genome_output_directory):

    # get sam output filename
    config_file = os.path.basename(config_path)
    (config_base, config_ext) = os.path.splitext(config_file)
    genome_name = os.path.basename(genome_output_directory)
    samout_path = utils.set_directory_path(genome_output_directory) + genome_name + ".align.sam"

    return samout_path


## get_bamout_path
def get_bamout_path(config_path, genome_output_directory):

    # get sam output filename
    config_file = os.path.basename(config_path)
    (config_base, config_ext) = os.path.splitext(config_file)
    genome_name = os.path.basename(genome_output_directory)
    bamout_path = utils.set_directory_path(genome_output_directory) + genome_name + ".align.bam"

    return bamout_path


## get_filtered_samout_path
def get_filtered_samout_path(config_path, genome_output_directory):

    # get sam output filename
    config_file = os.path.basename(config_path)
    (config_base, config_ext) = os.path.splitext(config_file)
    genome_name = os.path.basename(genome_output_directory)
    filtered_samout_path = utils.set_directory_path(genome_output_directory) + genome_name + ".filtered.sam"

    return filtered_samout_path


## get_filtered_bamout_path
def get_filtered_bamout_path(config_path, genome_output_directory):

    # get sam output filename
    config_file = os.path.basename(config_path)
    (config_base, config_ext) = os.path.splitext(config_file)
    genome_name = os.path.basename(genome_output_directory)
    filtered_bamout_path = utils.set_directory_path(genome_output_directory) + genome_name + ".filtered.bam"

    return filtered_bamout_path


## get_samlog_path
def get_samlog_path(config_path, genome_output_directory):

    # get sam log filename
    config_file = os.path.basename(config_path)
    (config_base, config_ext) = os.path.splitext(config_file)
    genome_name = os.path.basename(genome_output_directory)
    samlog_path = utils.set_directory_path(genome_output_directory) + genome_name + ".align.log"

    return samlog_path

## get_bamlog_path
def get_bamlog_path(config_path, genome_output_directory):

    # get sam log filename
    config_file = os.path.basename(config_path)
    (config_base, config_ext) = os.path.splitext(config_file)
    genome_name = os.path.basename(genome_output_directory)
    bamlog_path = utils.set_directory_path(genome_output_directory) + genome_name + ".align.log"

    return bamlog_path


## get read id for (paired vs unpaied) and (fastq vs fasta)
def get_read_id(read_id, paired_end_flag, fasta_reads_flag):

    # check paired end
    if (paired_end_flag):

        # check fasta read format
        if (fasta_reads_flag):
            read_id_list = read_id.split('.')

            # check list length > 1
            if len(read_id_list) > 1:
                read_id = '.'.join(read_id_list[0:-1])
            else:
                read_id = read_id_list[0]

    return read_id


## get overall alignment rate
def get_alignment_rate(samlog_path, bowtie_selection):

    # initialize 
    alignment_rate = 0.0;

    # for bowtie
    if (bowtie_selection == "1"):
        # find the string
        target_string = "reads with at least one reported alignment"

        # open bowtie log file and read lines
        try:
            for read_line in open(samlog_path):
                read_line = read_line.strip()
                if read_line.find(target_string) >= 0:
                   # now value is xx.xx%
                   alignment_rate = read_line[read_line.find("(")+1:read_line.find(")")]
                   # delete %
                   alignment_rate = alignment_rate[0:-1]
        except:
            alignment_rate = 0.0
            
    # for bowtie
    elif (bowtie_selection == "2"):
        # find the string
        target_string = "overall alignment rate"
        number_total_reads_string = "reads; of these:"
        number_noalign_reads_string = "aligned concordantly 0 times"
        number_total_reads = 1;
        number_noalign_reads = 0;

        # open bowtie log file and read lines
        with open(samlog_path, 'r') as f:
            for read_line in f:
                read_line = read_line.strip()
                if read_line.find(number_total_reads_string) >= 0:
                   number_total_reads = int(read_line.split(" ")[0])
                elif read_line.find(number_noalign_reads_string) >= 0:
                   number_noalign_reads = int(read_line.split(" ")[0])
            # calculate alignment rate
            alignment_rate = 100.0*(float(number_total_reads - number_noalign_reads) / float(number_total_reads))
            

    return alignment_rate 


# Get output_genome_fasta_sequence_list
def get_output_genome_fasta_sequence_list(genome_fasta_id_list, \
        genome_fasta_sequence_list, sequence_all_count_list, diff_genome_fasta_sequence_list, \
        minumum_coverage_length, minimum_average_coverage_depth):

    # loop fasta files
    for genome_fasta_enum, genome_fasta_id in enumerate(genome_fasta_id_list):

        # get diff genome fasta sequence length
        diff_genome_fasta_sequence_size = len(diff_genome_fasta_sequence_list[genome_fasta_enum])

        # initialize
        previous_letter_updated_tag = False

        # loop sequence letters
        for sequence_letter_position in range(0, diff_genome_fasta_sequence_size):

            # get original sequence letter (lower case) from genome_fasta_sequence_list
            original_sequence_letter = genome_fasta_sequence_list[genome_fasta_enum][sequence_letter_position].lower()
            # get diff sequence letter (upper case) from diff_genome_fasta_sequence_list
            diff_sequence_letter = diff_genome_fasta_sequence_list[genome_fasta_enum][sequence_letter_position].upper()

            # first letter
            if (sequence_letter_position == 0):
                # if the letter is not updated
                if (diff_sequence_letter == 'X'):
                    diff_genome_fasta_sequence_list[genome_fasta_enum][sequence_letter_position] = original_sequence_letter
                # if the letter is updated
                else:
                    updated_sequence_start_position = sequence_letter_position
                    previous_letter_updated_tag = True
            # last letter
            elif (sequence_letter_position == (diff_genome_fasta_sequence_size-1) ):
                # if the letter is not updated
                if (diff_sequence_letter == 'X'):
                    diff_genome_fasta_sequence_list[genome_fasta_enum][sequence_letter_position] = original_sequence_letter
                # if the letter is updated
                else:
                    if (previous_letter_updated_tag is False):
                        updated_sequence_start_position = sequence_letter_position 
                    updated_sequence_end_position = sequence_letter_position + 1

                    # get coverage length and depth
                    coverage_length = updated_sequence_end_position - updated_sequence_start_position - 1
                    coverage_count_sum = sum(sequence_all_count_list[genome_fasta_enum][updated_sequence_start_position:updated_sequence_start_position])
                    average_coverage_depth = coverage_count_sum / float(coverage_length)

                    # if the coverage do not greater than or equal to the threshold, update with original sequences
                    if (coverage_length < minumum_coverage_length) or (average_coverage_depth < minimum_average_coverage_depth):
                        diff_genome_fasta_sequence_list[genome_fasta_enum][updated_sequence_start_position:updated_sequence_end_position] \
                            = [x.lower() for x in genome_fasta_sequence_list[genome_fasta_enum][updated_sequence_start_position:updated_sequence_end_position]]
            # inside
            else:
                # if the letter is not updated
                if (diff_sequence_letter == 'X'):
                    diff_genome_fasta_sequence_list[genome_fasta_enum][sequence_letter_position] = original_sequence_letter
                    if (previous_letter_updated_tag is True):
                        updated_sequence_end_position = sequence_letter_position

                        # get coverage length and depth
                        coverage_length = updated_sequence_end_position - updated_sequence_start_position - 1
                        coverage_count_sum = sum(sequence_all_count_list[genome_fasta_enum][updated_sequence_start_position:updated_sequence_end_position])
                        average_coverage_depth = coverage_count_sum / float(coverage_length)

                        # if the coverage do not greater than or equal to the threshold, update with original sequences
                        if (coverage_length < minumum_coverage_length) or (average_coverage_depth < minimum_average_coverage_depth):
                            diff_genome_fasta_sequence_list[genome_fasta_enum][updated_sequence_start_position:updated_sequence_end_position] \
                                = [x.lower() for x in genome_fasta_sequence_list[genome_fasta_enum][updated_sequence_start_position:updated_sequence_end_position]]
                    previous_letter_updated_tag = False
                # if the letter is updated
                else:
                    if (previous_letter_updated_tag is False):
                        updated_sequence_start_position = sequence_letter_position
                    previous_letter_updated_tag = True
    return (diff_genome_fasta_sequence_list)


## Get diff genome fasta sequence list
def get_diff_genome_fasta_sequence_list(genome_fasta_id_list, genome_fasta_sequence_list, \
        sequence_a_tvalue_list, sequence_c_tvalue_list, sequence_g_tvalue_list, sequence_t_tvalue_list, \
        sequence_all_count_list):

    # initialize 
    diff_genome_fasta_sequence_list = []
    diff_sequence_letter_list = []
    diff_sequence_letter_position_list = []

    # loop fasta files
    for genome_fasta_enum, genome_fasta_id in enumerate(genome_fasta_id_list):

        # for each fasta
        diff_genome_fasta_sequence_list.append([])
        diff_sequence_letter_list.append([])
        diff_sequence_letter_position_list.append([])

        # get genome fasta sequence length
        genome_fasta_sequence_size = len(genome_fasta_sequence_list[genome_fasta_enum])

        # loop sequence letters
        for sequence_letter_position in range(0, genome_fasta_sequence_size):

            # initialize letter_tvalue_list
            letter_tvalue_list = []

            # check sequence_all_count_list is 0
            if sequence_all_count_list[genome_fasta_enum][sequence_letter_position] == 0:
                diff_genome_fasta_sequence_list[genome_fasta_enum].append('X')
            # check sequence_all_count_list is not 0
            else:
                # get original letter (upper case) from genome_fasta_sequence_list
                original_sequence_letter = genome_fasta_sequence_list[genome_fasta_enum][sequence_letter_position].upper()
            
                # letter_tvalue_list [tvalue for a, tvalue for c, tvalue for g, tvalue for t]
                letter_tvalue_list.append(sequence_a_tvalue_list[genome_fasta_enum][sequence_letter_position])
                letter_tvalue_list.append(sequence_c_tvalue_list[genome_fasta_enum][sequence_letter_position])
                letter_tvalue_list.append(sequence_g_tvalue_list[genome_fasta_enum][sequence_letter_position])
                letter_tvalue_list.append(sequence_t_tvalue_list[genome_fasta_enum][sequence_letter_position])

                # get output letter by tvalue
                # case 1) find only one maximum tvalue
                max_letter_tvalue = max(letter_tvalue_list)
                max_letter_index_list = [i for i, j in enumerate(letter_tvalue_list) if j == max_letter_tvalue]
                if len(max_letter_index_list) == 1:
                    if max_letter_index_list[0] == 0:
                        diff_sequence_letter = 'A'
                    elif max_letter_index_list[0] == 1:
                        diff_sequence_letter = 'C'
                    elif max_letter_index_list[0] == 2:
                        diff_sequence_letter = 'G'
                    elif max_letter_index_list[0] == 3:
                        diff_sequence_letter = 'T'
                else:
                    diff_sequence_letter = original_sequence_letter

                # update diff_genome_fasta_sequence_list with diff_sequence_letter
                diff_genome_fasta_sequence_list[genome_fasta_enum].append(diff_sequence_letter)

                # if diff_sequence_letter is different to original_sequence_letter
                if diff_sequence_letter != original_sequence_letter:
                    diff_sequence_letter_list[genome_fasta_enum].append(diff_sequence_letter)
                    diff_sequence_letter_position_list[genome_fasta_enum].append(sequence_letter_position)

    return (diff_genome_fasta_sequence_list, diff_sequence_letter_list, diff_sequence_letter_position_list)


## Class for ParseSamout
class ParseSamout:
    # init
    def __init__(self, read_line_list):
        self.read_line_list = read_line_list
    # get read id
    def get_read_id(self, paired_end_reads_flag, fasta_reads_flag):
        read_id = get_read_id(self.read_line_list[0], paired_end_reads_flag, fasta_reads_flag)
        return read_id
    # get read flag
    def get_read_flag(self):
        return  int(self.read_line_list[1])
    # get reference name
    def get_reference_name(self):
        return  self.read_line_list[2]
    # get reference name
    def get_reference_name(self):
        return  self.read_line_list[2]
    # get reference match position (leftmost)
    def get_reference_position(self):
        return  int(self.read_line_list[3])
    # get read sequence
    def get_read_sequence(self):
        return  self.read_line_list[9]
    # get start and end position
    def get_start_end_position(self):
        read_sequence_size = len(self.get_read_sequence())
        start_position = self.get_reference_position() - 1
        end_position = start_position + read_sequence_size
        return  (start_position, end_position)


# initialize sequence tvalue list
def initialize_sequence_count(genome_fasta_id_list, genome_fasta_sequence_list):

    # initialize sequence count list
    sequence_a_tvalue_list = []
    sequence_c_tvalue_list = []
    sequence_g_tvalue_list = []
    sequence_t_tvalue_list = []
    sequence_all_count_list = []

    # loop fasta files
    for genome_fasta_enum, genome_fasta_id in enumerate(genome_fasta_id_list):
        sequence_a_tvalue_list.append([])
        sequence_c_tvalue_list.append([])
        sequence_g_tvalue_list.append([])
        sequence_t_tvalue_list.append([])
        sequence_all_count_list.append([])
        genome_fasta_sequence_size = len(genome_fasta_sequence_list[genome_fasta_enum])

        # loop sequence letters
        for sequence_letter_position in range(0, genome_fasta_sequence_size):
            sequence_a_tvalue_list[genome_fasta_enum].append(0.0)
            sequence_c_tvalue_list[genome_fasta_enum].append(0.0)
            sequence_g_tvalue_list[genome_fasta_enum].append(0.0)
            sequence_t_tvalue_list[genome_fasta_enum].append(0.0)
            sequence_all_count_list[genome_fasta_enum].append(0)

    return (sequence_a_tvalue_list, sequence_c_tvalue_list, sequence_g_tvalue_list, \
            sequence_t_tvalue_list, sequence_all_count_list)


# get refrence name list and genome fasta sequences to list
def get_genome_fasta_id_sequence(genome_fasta_path_sublist, genome_name):

    # initialize
    genome_fasta_id_list = []
    fasta_sequence_list = []
    fasta_sequence_sublist = []
    genome_fasta_id_description_map = {}

    # loop fasta files for the genome
    for fasta_path in genome_fasta_path_sublist:

        # check file exist
        try:
            # open fasta file
            fasta_file = open(fasta_path, 'r')
        except IOError:
            sys.stderr.write("** The file, %s, does not exist!\n" %(fasta_path))
            utils.die("** Program exit!")

        # read line
        for fasta_line in fasta_file:
            fasta_line = fasta_line.strip()

            # if fasta id
            if fasta_line.startswith('>'):
                fasta_id_line_list = fasta_line.split(' ')

                # for genome fasta id
                genome_fasta_id = fasta_id_line_list[0][1:]
                genome_fasta_id_list.append(genome_fasta_id)

                # for genome fasta id description
                genome_fasta_description = 'Description N/A'
                if len(fasta_id_line_list) > 1:
                    genome_fasta_description = ' '.join(fasta_id_line_list[1:])
                genome_fasta_id_description_map[genome_fasta_id] = genome_fasta_description

                # for genome fasta sequence
                if fasta_sequence_sublist != []:
                    fasta_sequence = ''.join(fasta_sequence_sublist)
                    fasta_sequence_list.append(fasta_sequence)
                    fasta_sequence_sublist = []
            # fasta sequence
            else:
                fasta_sequence_sublist.append(fasta_line)
            
        # for genome fasta sequence
        fasta_sequence = ''.join(fasta_sequence_sublist)
        fasta_sequence_list.append(fasta_sequence)
        fasta_sequence_sublist = []
        fasta_file.close()

    return (genome_fasta_id_list, genome_fasta_id_description_map, fasta_sequence_list)

# get target genome list for reonctruction
def get_target_genome(gvector_path, reconstruction_selection, reconstruction_cutoff_abundance, reconstruction_genome_name):

    # genome_name_map(key:genome_index, val:genome_name)
    genome_name_map = {}
    # genome_chance_map
    genome_chance_map = {}
    # target_genome_index_list
    target_genome_index_list = []
    # to check the reconstruction_genome_name
    reconstruction_genome_name_exist_tag = False

    # check g-vector file
    if not utils.check_path_exist(gvector_path):
        utils.die("** Error: No such file or directory: " + gvector_path)
        
    # get target genome list (g_vector percentage >= cutoff_reconstruction_percentage)
    with open(gvector_path, 'r') as f:
		for gvector_line in f:
			gvector_line = gvector_line.strip()
			gvector_columns = gvector_line.split('\t')

			# consider first column is '@': @	Genome_Index	Genome_Name	Alignment_Rate
			if gvector_columns[0] == '@':
				genome_index = int(gvector_columns[1])
				genome_name = gvector_columns[2]
				genome_name_map[genome_index] = genome_name
			# consider first column is '*': *	Genome_Index	Percentage_Chance
			elif gvector_columns[0] == '*':
				genome_index = int(gvector_columns[1])
				# ScaledPercentageChance is considered now.
				# genome_chance = float(gvector_columns[2])
				genome_chance = float(gvector_columns[3])
				genome_chance_map[genome_index] = genome_chance

				# reconstruct all genomes >= cut-off
				if reconstruction_selection == 1:
					# if genome_chance >= Reconstruction_Cutoff_Abundance, then save
					if genome_chance >= reconstruction_cutoff_abundance:
						target_genome_index_list.append(genome_index)
				# reconstruct one target genome
				else:
					genome_name = genome_name_map[genome_index]
					if genome_name == reconstruction_genome_name:
						reconstruction_genome_name_exist_tag = True
						target_genome_index_list.append(genome_index)

    if not reconstruction_genome_name_exist_tag:
        utils.die("** Check the config file! Reconstruction_Genome_Name value does not match to the SIGMA gvector results!")

    if len(target_genome_index_list) == 0:
        utils.die("** No target genomes exist for reconstruction!")

    return (genome_name_map, genome_chance_map, target_genome_index_list)


# get target genome list for filtering
def get_filtering_target_genome(gvector_path, filtering_genome_name):

    # genome_name_map(key:genome_index, val:genome_name)
    genome_name_map = {}
    # genome_chance_map
    genome_chance_map = {}
    # target_genome_index_list
    target_genome_index_list = []
    # to check the filtering_genome_name
    filtering_genome_name_exist_tag = False

    # check g-vector file
    if not utils.check_path_exist(gvector_path):
        utils.die("** Error: No such file or directory: " + gvector_path)
        
    # get target genome list
    with open(gvector_path, 'r') as f:
		for gvector_line in f:
			gvector_line = gvector_line.strip()
			gvector_columns = gvector_line.split('\t')

			# consider first column is '@': @	Genome_Index	Genome_Name	Alignment_Rate
			if gvector_columns[0] == '@':
				genome_index = int(gvector_columns[1])
				genome_name = gvector_columns[2]
				genome_name_map[genome_index] = genome_name
			# consider first column is '*': *	Genome_Index	Percentage_Chance
			elif gvector_columns[0] == '*':
				genome_index = int(gvector_columns[1])
				# if PercentageChance is considered
				genome_chance = float(gvector_columns[2])
				# if ScaledPercentageChance is considered
				#genome_chance = float(gvector_columns[3])
				genome_chance_map[genome_index] = genome_chance
				genome_name = genome_name_map[genome_index]
				if genome_name == filtering_genome_name:
					filtering_genome_name_exist_tag = True
					target_genome_index_list.append(genome_index)

    if not filtering_genome_name_exist_tag:
        utils.die("** Check the config file! Filtering_Genome_Name value does not match to the SIGMA gvector results!")

    if len(target_genome_index_list) == 0:
        utils.die("** No target genomes exist for filtering!")

    return (genome_name_map, genome_chance_map, target_genome_index_list)


# get genome_directory_map and output_genome_directory_map
def get_genome_directory_map(genome_directory_list, output_genome_directory_list):

    # genome_directory_map: key->genome_name, val->genome_directory
    genome_directory_map = {}
    # output_genome_directory_map: key->genome_name, val->output_genome_directory
    output_genome_directory_map = {}

    # genome_directory_map
    for genome_directory in genome_directory_list:
        genome_name = os.path.basename(genome_directory)
        genome_directory_map[genome_name] = genome_directory

    # output_genome_directory_map
    for output_genome_directory in output_genome_directory_list:
        genome_name = os.path.basename(output_genome_directory)
        output_genome_directory_map[genome_name] = output_genome_directory

    return (genome_directory_map, output_genome_directory_map)


# get T(i,j) for j genome
def get_tvalue(qmatrix_path, genome_chance_map, target_genome_index):

    # initialize tvalue_map
    tvalue_map = {}

    # open and read qmatrix
    for qmatrix_line in open(qmatrix_path):
        qmatrix_line = qmatrix_line.strip()
        line_columns = qmatrix_line.split('\t')
        line_columns_size = len(line_columns)

        # read data
        if line_columns[0] == '*':
            # get read_id
            read_id = line_columns[1]

            # loop for sum_qvalue_times_genome_chance
            sum_qvalue_times_genome_chance = 0.0
            for x in range(2, line_columns_size):
                qvalue_set = line_columns[x]
                qvalue_set_list = qvalue_set.split('=')
                genome_index = int(qvalue_set_list[0])
                qvalue = float(qvalue_set_list[1])
                genome_chance = 0.0
                if genome_index in genome_chance_map:
                    genome_chance = genome_chance_map[genome_index]
                sum_qvalue_times_genome_chance += qvalue*genome_chance

            # for target_qvalue_times_genome_chance
            for x in range(2, line_columns_size):
                qvalue_set = line_columns[x]
                qvalue_set_list = qvalue_set.split('=')
                genome_index = int(qvalue_set_list[0])
                qvalue = float(qvalue_set_list[1])
                genome_chance = 0.0
                if genome_index in genome_chance_map:
                    genome_chance = genome_chance_map[genome_index]
            
                # check genome_index to target_genome_index
                if genome_index == target_genome_index:
                    target_qvalue_times_genome_chance = qvalue*genome_chance
                    # Method 1) calculate tvalue = Pr(g_j|e_i) = target_qvalue_times_genome_chance / sum_qvalue_times_genome_chance
                    # tvalue = target_qvalue_times_genome_chance / sum_qvalue_times_genome_chance
                    # Method 2) calculate tvalue = Pr(e_i,g_i) = Pr(e_i|g_j)*Pr(g_j) = target_qvalue_times_genome_chance
                    tvalue = target_qvalue_times_genome_chance
                    tvalue_map[read_id] = tvalue

    return tvalue_map


# get target read IDs based on T-value
def reads_by_tvalue(qmatrix_path, genome_chance_map, target_genome_index):

    # initialize taget reads list
    target_reads_list = set()

    # open and read qmatrix
    for qmatrix_line in open(qmatrix_path):
        qmatrix_line = qmatrix_line.strip()
        line_columns = qmatrix_line.split('\t')
        line_columns_size = len(line_columns)

        # read data
        if line_columns[0] == '*':
            # get read_id
            read_id = line_columns[1]

            # initialize list
            genome_index_list = []
            qvalue_list = []
            tvalue_list = []

            # for target_qvalue_times_genome_chance
            for x in range(2, line_columns_size):
                # get genome index and tvalue
                qvalue_set = line_columns[x]
                qvalue_set_list = qvalue_set.split('=')
                genome_index = int(qvalue_set_list[0])
                qvalue = float(qvalue_set_list[1])
                genome_chance = 0.0
                if genome_index in genome_chance_map:
                    genome_chance = genome_chance_map[genome_index]
                tvalue = qvalue*genome_chance

                # append to list
                genome_index_list.append(genome_index)
                qvalue_list.append(qvalue)
                tvalue_list.append(tvalue)
            
            # find max index and value
            max_tvalue = max(tvalue_list)
            max_index = tvalue_list.index(max_tvalue)
            max_genome_index = genome_index_list[max_index]

            # check max_genome_index to target_genome_index
            if max_genome_index == target_genome_index:
                target_reads_list.add(read_id)

            # clear list
            del genome_index_list[:]
            del qvalue_list[:]
            del tvalue_list[:]

    return target_reads_list


## get_output_genome_path
def get_output_genome_fasta_path_sublist(output_genome_directory, genome_fasta_path_sublist):

    # initialize list
    output_genome_fasta_path_sublist = []

    # for loop 
    for genome_fasta_path in genome_fasta_path_sublist:
        genome_fasta_name = os.path.basename(genome_fasta_path)
        (genome_fasta_name_base, genome_fasta_name_ext) = os.path.splitext(genome_fasta_name)
        output_genome_fasta_name = genome_fasta_name_base + ".full" + genome_fasta_name_ext
        output_genome_fasta_path = utils.set_directory_path(output_genome_directory) + output_genome_fasta_name

        # append to list
        output_genome_fasta_path_sublist.append(output_genome_fasta_path)
    
    return output_genome_fasta_path_sublist

## get_output_diff_sequence_path_sublist
def get_output_diff_sequence_path_sublist(output_genome_directory, genome_fasta_path_sublist):

    # initialize list
    output_diff_sequence_path_sublist = []

    # for loop 
    for genome_fasta_path in genome_fasta_path_sublist:
        genome_fasta_name = os.path.basename(genome_fasta_path)
        (genome_fasta_name_base, genome_fasta_name_ext) = os.path.splitext(genome_fasta_name)
        output_diff_sequence_name = genome_fasta_name_base + ".SNPs.txt"
        output_diff_sequence_path = utils.set_directory_path(output_genome_directory) + output_diff_sequence_name

        # append to list
        output_diff_sequence_path_sublist.append(output_diff_sequence_path)
    
    return output_diff_sequence_path_sublist


## open Q-matrix and count reads
def get_qmatrix_reads_count(qmatrix_path):

    # check path exist
    if utils.check_path_exist(qmatrix_path):

        # initialize
        qmatrix_reads_count = 0

        # open and read qmatrix
        for qmatrix_line in open(qmatrix_path):
            qmatrix_line = qmatrix_line.strip()
            line_columns = qmatrix_line.split('\t')

            # read data
            if line_columns[0] == '*':
                qmatrix_reads_count += 1

        return qmatrix_reads_count

    else:
        sys.stderr.write("\n** Cannot open %s.\n" % qmatrix_path)
        utils.die("** Program exit!")


## get gvector genome list
def get_gvector_genome_list(gvector_path):

    # check file exist
    if utils.check_path_exist(gvector_path):

        # initialize
        gvector_genome_list = []
            
        # open file
        gvector_file = open(gvector_path, 'r')
        for gvector_line in gvector_file:
            gvector_line = gvector_line.strip()
            gvector_columns = gvector_line.split('\t')

            # consider first column is '@': @   Genome_Index    Genome_Name Alignment_Rate
            if gvector_columns[0] == '@':
                gvector_genome_list.append(gvector_line)

        return gvector_genome_list

    else:
        sys.stderr.write("\n** Cannot open %s.\n" % qmatrix_path)
        utils.die("** Program exit!")


## open Q-matrix and get data
def get_qmatrix_data(qmatrix_path):

    # check path exist
    if utils.check_path_exist(qmatrix_path):

        # initialize
        qmatrix_comments_list = []
        qmatrix_data_list = []

        # open and read qmatrix
        for qmatrix_line in open(qmatrix_path):
            qmatrix_line = qmatrix_line.strip()
            line_columns = qmatrix_line.split('\t')

            # read data
            if line_columns[0] == '*':
                qmatrix_data_list.append(qmatrix_line)
            else:
                qmatrix_comments_list.append(qmatrix_line)

        return (qmatrix_comments_list, qmatrix_data_list)

    else:
        sys.stderr.write("\n** Cannot open %s.\n" % qmatrix_path)
        utils.die("** Program exit!")


## write Q-matrix bootstrapping results
def write_bootstrap_qmatrix(bootstrap_qmatrix_path, qmatrix_comments_list, qmatrix_data_list):

    # check path exist
    if utils.check_path_exist(bootstrap_qmatrix_path):
        check_call(["rm", "-f", bootstrap_qmatrix_path], stdout = PIPE, stderr = sys.stderr)
    bootstrap_qmatrix_out = open(bootstrap_qmatrix_path, 'wb')

    # loop for qmatrix_comments_list
    for qmatrix_comments_list_item in qmatrix_comments_list:
        bootstrap_qmatrix_out.write(str(qmatrix_comments_list_item) + '\n')

    # get qmatrix reads count
    qmatrix_data_count = len(qmatrix_data_list)

    # loop for qmatrix_data_list
    for qmatrix_data_index in xrange(0, qmatrix_data_count):
        random_index = random.randint(0, qmatrix_data_count - 1)
        bootstrap_qmatrix_out.write(str(qmatrix_data_list[random_index]) + '\n')


## get_bootstrap_gvector
def get_bootstrap_gvector(bootstrap_gvector_path):

    # check path exist
    if utils.check_path_exist(bootstrap_gvector_path):

        # initialize
        bootstrap_gvector_map = {}

        # open and read bootstrap_gvector_path
        for gvector_line in open(bootstrap_gvector_path):
            gvector_line = gvector_line.strip()
            gvector_columns = gvector_line.split('\t')

            # read data
            if gvector_columns[0] == '*':
                genome_index = int(gvector_columns[1])
                # ScaledPercentageChance is considered now.
                # genome_chance = float(gvector_columns[2])
                genome_chance = float(gvector_columns[3])
                bootstrap_gvector_map[genome_index] = genome_chance

        return bootstrap_gvector_map

    else:
        sys.stderr.write("\n** Cannot open %s!\n" % bootstrap_gvector_path)
        utils.die("** Program exit!")


## calculate_bootstrap_gvector_stat
def calculate_bootstrap_gvector_stat(bootstrap_gvector_list):

    # initialize
    bootstrap_gvector_stat_list = []

    # get size
    for genome_index, genome_chance_list in enumerate(bootstrap_gvector_list):

        bootstrap_gvector_stat_list.append([])

        # size n
        n = len(genome_chance_list)
        # average
        average_percentage_chance = sum(genome_chance_list) / float(n)
        # STD
        stdev_percentage_chance = sqrt( sum([(x-average_percentage_chance)**2 for x in genome_chance_list]) / float(n-1) )
        # 95% upper confidence bount 
        upper_confidence_bound = average_percentage_chance + 1.96*( stdev_percentage_chance / float(sqrt(n)) )
        # 95% lower confidence bount 
        lower_confidence_bound = average_percentage_chance - 1.96*( stdev_percentage_chance / float(sqrt(n)) )

        # append to list
        bootstrap_gvector_stat_list[genome_index].append(average_percentage_chance)
        bootstrap_gvector_stat_list[genome_index].append(stdev_percentage_chance)
        bootstrap_gvector_stat_list[genome_index].append(upper_confidence_bound)
        bootstrap_gvector_stat_list[genome_index].append(lower_confidence_bound)

    return bootstrap_gvector_stat_list


## check bowtie indexes are already built
def check_bowtie_index_built(genome_index_base,
                             genome_fasta_path_sublist):

    # get bowtie index filepath (basename.1.bt2)
    genome_index_bt1_path = genome_index_base + ".1.bt2"
    # get fasta filepath (first fasta file)
    genome_fasta_1_path = genome_fasta_path_sublist[0]

    # check if the bowtie index file exist
    if utils.check_path_exist(genome_index_bt1_path):

        # get file creation time
        index_creation_time = os.path.getctime(genome_index_bt1_path)
        fasta_creation_time = os.path.getctime(genome_fasta_1_path)
       
        # if index_creation_time is newer than fasta_creation_time
        if index_creation_time >= fasta_creation_time:
             return True    
        else:
             return False

    # if bowtie indexes do not exist, return false
    else:
        return False


## count number of reads
def count_number_of_reads(reads_filepath, fasta_reads_flag):

    # number of reads
    number_of_reads = 0

    # check path exist
    if utils.check_path_exist(reads_filepath):

        # check fasta or fastq format
        if fasta_reads_flag:   # fasta format
            for read_line in open(reads_filepath,):
                read_line = read_line.strip()
                # fasta read id starts with >
                if read_line.startswith('>'):
                    number_of_reads += 1
        else:    # fastq format
            for read_line in open(reads_filepath):
                read_line = read_line.strip()
                # fastq read id starts with @
                if read_line.startswith('@'):
                    number_of_reads += 1
    else:
        sys.stderr.write("\n** Cannot open %s.\n" % reads_filepath)
        utils.die("** Program exit!")

    return number_of_reads


