#!/usr/bin/python

"""
vcf_filter.py

filter vcf file


Created by Tae-Hyuk (Ted) Ahn on 04/01/2013.
Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). Allrights reserved.
"""


import sys, warnings, os, re
from datetime import datetime, date, time
from subprocess import Popen, PIPE, check_call, STDOUT
import getopt

## Version control
version = "0.0.1 (Alpha)"

## Help message
help_message = '''

  [Usage]
    vcf_filter.py [options] -i <input vcf file> -o <output vcf file>

  [Options]
    -q <int or float> : filter out below quality (phred scaled score)
                        Ex) -q 20 is phred-20 (1% error) quality
    -f <float>        : filter out below allele frequency (0.0-1.0)
                        Ex) -f 0.5 filter out a case that has 4 alleles out of 10 reads depth
    -m                : report homozygous (FQ value is negative from samtools)
    -h/--help
    -v/--version

'''

## Class Usage
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg 


## Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)


## string class
class OptionString:

     input_filename_str = "input_filename"
     output_filename_str = "output_filename"
     qual_str = "qual"
     homo_str = "homo"
     freq_str = "freq"


## Parse options
def parse_option(argv, version, help_message):

    # option map dictionay
    option_map = {}

    # get arguments
    try:
        opts, args = getopt.getopt(argv[1:], "hvVi:o:q:f:m",
                                   ["help","version"])

    # Error handling of options
    except getopt.error, msg:
        raise Usage(msg)

    # get program name
    program_name =  sys.argv[0].split("/")[-1]

    # Basic options
    for option, value in opts:
        # help
        if option in ("-h", "--help"):
            raise Usage(help_message)
        # version
        if option in ("-v", "-V", "--version"):
            print "\n%s V%s\n" % (program_name, version)
            sys.exit(0)
        # -i
        if option in ("-i"):
            option_map[OptionString.input_filename_str] = value
        # -o
        if option in ("-o"):
            option_map[OptionString.output_filename_str] = value
        # -q: quality
        if option in ("-q"):
            option_map[OptionString.qual_str] = value
        # -f: frequency
        if option in ("-f"):
            option_map[OptionString.freq_str] = value
        # -m: homozygous
        if option in ("-m"):
            option_map[OptionString.homo_str] = 1

    # check you got two arguments
    if not OptionString.input_filename_str in option_map:
        raise Usage(help_message)
    if not OptionString.output_filename_str in option_map:
        raise Usage(help_message)

    return (option_map)


def work(option_map):

    # input and output
    input_filename  = option_map[OptionString.input_filename_str]
    output_filename = option_map[OptionString.output_filename_str]

    # qual
    QUAL_thread = 0.0
    if OptionString.qual_str in option_map:
       QUAL_thread = float(option_map[OptionString.qual_str])

    # freq
    freq_thread = 0.0
    if OptionString.freq_str in option_map:
       freq_thread = float(option_map[OptionString.freq_str])

    # homo
    HOMO_bool = False
    if OptionString.homo_str in option_map:
       HOMO_bool = True

    # open output 
    output = open(output_filename, 'wb')

    # open the fasta file and write
    with open(input_filename, 'r') as f:
        for line in f:
            line = line.strip()
            if (line.startswith("#")):
                output.write(line + "\n") 
            else:
                line_list = line.split("\t")

                #print line_list

                # QUAL field
                QUAL = float(line_list[5])

                # FQ field
                FQ_obj = re.search(r'FQ=(\W\d*)|FQ=(.*);', line_list[7])
                # FQ handle
                FQ = 0.0
                if FQ_obj.group(1):
                    FQ = float(FQ_obj.group(1))
                elif FQ_obj.group(2):
                    FQ = float(FQ_obj.group(2))
                #print "FQ=", FQ
        
                # DP field
                DP_obj = re.search(r'DP=(.*?);', line_list[7], re.M|re.I)
                # DP field
                DP = 0.0
                if DP_obj.group():
                    DP = float(DP_obj.group(1))
                #print "DP=", DP
        
                # DP4 field
                DP4_obj = re.search(r'DP4=(.*?);', line_list[7], re.M|re.I)
                DP4 = 0.0
                # DP4 field
                if DP4_obj.group():
                    DP4_list = DP4_obj.group(1).split(',')
                    DP4 = float(DP4_list[2]) + float(DP4_list[3])
                #print "DP4=", DP4

                # calculate frequency
                freq = (DP4/DP)
                  
                # check homo
                if (HOMO_bool) and not(OptionString.qual_str in option_map) and not(OptionString.freq_str in option_map):
                    if FQ < 0:
                        output.write(line + "\n") 
                elif not(HOMO_bool) and (OptionString.qual_str in option_map) and not(OptionString.freq_str in option_map):
                    if QUAL >= QUAL_thread:
                        output.write(line + "\n") 
                elif not(HOMO_bool) and not(OptionString.qual_str in option_map) and (OptionString.freq_str in option_map):
                    if (freq >= freq_thread):
                        output.write(line + "\n") 
                elif (HOMO_bool) and (OptionString.qual_str in option_map) and not(OptionString.freq_str in option_map):
                    if (FQ < 0) and (QUAL >= QUAL_thread):
                        output.write(line + "\n") 
                elif (HOMO_bool) and not(OptionString.qual_str in option_map) and (OptionString.freq_str in option_map):
                    if (FQ < 0) and (freq >= freq_thread):
                        output.write(line + "\n") 
                elif not(HOMO_bool) and (OptionString.qual_str in option_map) and (OptionString.freq_str in option_map):
                    if (QUAL >= QUAL_thread) and (freq >= freq_thread):
                        output.write(line + "\n") 
                elif (HOMO_bool) and (OptionString.qual_str in option_map) and (OptionString.freq_str in option_map):
                    if (FQ < 0) and (QUAL >= QUAL_thread) and (freq >= freq_thread):
                        output.write(line + "\n") 
                        
def main(argv=None):

    # try to get arguments and error handling
    try:
        if argv is None:
            argv = sys.argv
        try:
            # get program name
            program_name =  os.path.basename(sys.argv[0])

            # parse option
            (option_map) = parse_option(argv, version, help_message)

            # display work start and time record
            start_time = datetime.now()
            sys.stderr.write("\n*********************************************************\n")
            sys.stderr.write("Beginning %s run (V%s)\n" % ( program_name, version))

            work(option_map)

            # time record, calculate elapsed time, and display work end
            finish_time = datetime.now()
            duration = finish_time - start_time
            sys.stderr.write("Ending %s run\n" % (program_name))
            sys.stderr.write("Total Elapsed Time =  %s [seconds]\n" % (duration))
            sys.stderr.write("*********************************************************\n\n")


        # Error handling
        except Usage, err:
            sys.stderr.write("%s: %s\n" %(os.path.basename(sys.argv[0]), str(err.msg)))
            return 2


    # Error handling
    except Usage, err:
        sys.stderr.write("%s: %s\n" %(os.path.basename(sys.argv[0]), str(err.msg)))
        sys.stderr.write("for help use -h/--help")
        return 2


## If this program runs as standalone, then exit.
if __name__ == "__main__":
    sys.exit(main())

