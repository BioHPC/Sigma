#!/usr/bin/env python

"""
utils.py

utils.py includes useful common classes, definitions,
and funcitons. 

Created by Tae-Hyuk (Ted) Ahn on 03/01/2013.
Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). Allrights reserved.
"""

## Import standard modules
import sys, warnings, os, re
from datetime import datetime, date, time
import math


## Class Usage
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


## Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)


## Returns the current time in a nice format
def curr_time():
    curr_time = datetime.now()
    return curr_time.strftime("%c")


## Format time as a pretty string
def format_time(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds)


## Check path exist
def check_path_exist(input_path):
    if os.path.exists(input_path):
        return True
    else:
        return False


## Check fasta format file
def check_fasta_format(filename):

    # file extension
    file_format = ['.fasta', '.fa', '.fas', '.fna', '.ffn']
 
    # check fasta format file
    filename_wo_ext, file_ext = os.path.splitext(filename)
    if file_ext in file_format:
        return True
    else:
        return False


## Check fastq format file
def check_fastq_format(filename):

    # fastq file extension
    file_format = ['.fastq', '.fq', '.faq']
 
    # check fastq format file
    filename_wo_ext, file_ext = os.path.splitext(filename)
    if file_ext in file_format:
        return True
    else:
        return False


## Check bowtie2 index format file
def check_bowtie_index_format(filename):

    # bowtie2 index file extension
    file_format = ['.bt2']
 
    # check bowtie2 format
    filename_wo_ext, file_ext = os.path.splitext(filename)
    if file_ext in file_format:
        return True
    else:
        return False


## Check bam file format file
def check_bam_format(filename):

    # bam file extension
    file_format = ['.bam']
 
    # check bam format
    filename_wo_ext, file_ext = os.path.splitext(filename)
    if file_ext in file_format:
        return True
    else:
        return False



## List to comma separated string
def list_to_comma_string(input_list):

    if len(input_list) > 1:
        converted_str = ','.join(input_list) 
    else:
        converted_str = ''.join(input_list)

    return converted_str


## Comma string to list
def comma_string_to_list(input_string):

    converted_list = []
    input_string = input_string.strip()
    converted_list = input_string.split(',')

    return converted_list


## check integer
def check_int(input_value):
    check_bool = True
    try: 
        value_check = int(input_value)
    except: 
        check_bool = False

    return check_bool


## set directory path (append "/" at last if directory path doesn't have it)
def set_directory_path(directory_path):
    if os.path.isdir(directory_path):
        directory_path = os.path.normpath(directory_path)
        directory_path += '/'
    else:
        sys.stderr.write("\n** Cannot open %s.\n" % directory_path)
        die("** Program exit!")

    return directory_path

        
## get base path with file path and working directory
def get_base_path(file_path, directory_path):

    # Get base path
    filename = os.path.basename(file_path)
    base_out = os.path.splitext(filename)[0]

    # If base common prefix ends with '_' or '.', then remove
    base_out = base_out[:-1] if (base_out[-1] in ('_', '.')) else base_out

    # set directory path
    directory_path = set_directory_path(directory_path)

    # base_out
    base_out = directory_path + base_out

    return base_out


## check file size is 0 or not
def check_file_size_not_zero(filename):

    if check_path_exist(filename):
        # use os.stat to check file size
        file_stat = os.stat(filename)
  
        # st_size ==0 -> false
        if file_stat.st_size == 0:
            return False
        else:
            return True
    else:
        sys.stderr.write("\n** Cannot open %s.\n" %(filename))
        die("** Program exit!")
    

## check pid alive
def check_pid_alive(pid):        
    """ Check For the existence of a unix pid. """
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True
