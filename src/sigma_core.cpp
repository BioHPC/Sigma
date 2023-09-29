//=============================================================================
// sigma_core.cpp
//   : For methods and functions of SIGMA.
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 10/01/2016
// Modified:
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include <math.h>
#include <numeric>
#include <iomanip>
#include "sigma_core.h"


//=============================================================================
// Load config file and get bowtie-build option 
//=============================================================================
void getBowtieBuildOption(
        const std::string & config_path,            // (in) config file path
        std::string & reference_genome_directory,   // (out) reference genome dir
        std::string & bowtie_build_option_cmd)      // (out) bowtie option command
{
    // load config file
    SigmaConfig config(config_path);

    // get variables
    std::string bowtie_build_path = config.getBowtieBuildPath();
    reference_genome_directory = config.getReferenceGenomeDirectory();

    //
    // bowtie-build option
    // Usage: 
    //   bowtie2-build [options]* <reference_in> <ebwt_base>
    // Options:
    //   -f : The reference inputs are fasta files 
    //   ... we don't neeed addtional options for indexing
    //

    // bowtie_option
    bowtie_build_option_cmd = bowtie_build_path + " -f";
}


//=============================================================================
// Check genome index is already build in previous
//=============================================================================
bool isBowtieIndexExist(
        const std::string & genome_index_base,                          // (in) genome index base
        const std::vector<std::string> & genome_fastas_path_sublist)    // (in) genome fasta path sublist
{

    // string varialbe
    std::string genome_index_filepath, genome_fasta_filepath;

    // for bowtie2
    genome_index_filepath = genome_index_base + ".1.bt2";

    // if bowtie index file exist
    if (Utils::isFileExist(genome_index_filepath)) {

        // modifed time
        unsigned int index_time, fasta_time, max_fasta_time;

        // create a file attribute structure
        struct stat index_attrib, fasta_attrib; 

        // get the attributes of index
        stat(genome_index_filepath.c_str(), &index_attrib); 

        // get modified time for index
        index_time = index_attrib.st_mtime;

        // check several fasta files, and get max modifed time
        max_fasta_time = 0;
        for (unsigned int i=0; i<genome_fastas_path_sublist.size(); i++) {

            genome_fasta_filepath = genome_fastas_path_sublist.at(i);

            // get the attributes of fasta
            stat(genome_fasta_filepath.c_str(), &fasta_attrib);

            // get modified time for index
            fasta_time = fasta_attrib.st_mtime;
            if (fasta_time > max_fasta_time)
                max_fasta_time = fasta_time;
        }


        // compare modifed time
        if (index_time > max_fasta_time) 
            return true;
        else
            return false;
    }
    else
        return false;
}


//=============================================================================
// Load config file and get bowtie option 
//=============================================================================
void getBowtieOption(
        const std::string & config_path,            // (in) config file path
        std::string & reference_genome_directory,   // (out) reference genome dir
        std::string & bowtie_option_cmd,            // (out) bowtie option command
        std::string & read_option_cmd,              // (out) bowtie read files and option 
        std::string & samtools_cmd)                 // (out) samtools command
{
    // load config file
    SigmaConfig config(config_path);

    // get variables
    reference_genome_directory = config.getReferenceGenomeDirectory();
    std::string bowtie_path = config.getBowtiePath();
    std::string samtools_path = config.getSamtoolsPath();
    std::string bowtie2_preset  = config.getBowtie2Preset();
    std::string maximum_mismatch_count  = config.getMaximumMismatchCount();
    std::string minimum_fragment_length = config.getMinimumFragmentLength();
    std::string maximum_fragment_length = config.getMaximumFragmentLength();
    std::string bowtie_threads_number = config.getBowtieThreadsNumber();
    
    // for bowtie read option command
    bool paired_end_reads_flag, fasta_reads_flag;
    config.getReadsType(paired_end_reads_flag, fasta_reads_flag);

    // 
    // bowie1 options:
    //  -n: Maximum number of mismatches permitted in the "seed" (default: 2)
    //  -l: The seed length (default: 28)
    //  -v: Alignments may have no more than v mismatches (0~3). Quality values are ignored. 
    //  -e: Maximum permitted total of quality values at all mismatched read positions (default: 70)
    //  --best: Make Bowtie guarantee that reported singleton alignments are "best" in terms of stratum
    //  -p %s : Launch NTHREADS parallel search threads (default: 1)
    //  -t : Print the wall-clock time 
    //  -I %s : The minimum fragment length for valid paired-end alignments (default: 0)
    //  -X %s : The maximum fragment length for valid paired-end alignments (default: 500)
    //
    // 
    // bowie2 options:
    //  -q:fastq, -f:fasta
    //  -1 <paired-end-1> -2 <paired-end-2> or -U <single-end>
    //  -p %s : Launch NTHREADS parallel search threads (default: 1)
    //  -t : Print the wall-clock time 
    //  --no-unal : Suppress SAM records for reads that failed to align.
    //  -I %s : The minimum fragment length for valid paired-end alignments (default: 0)
    //  -X %s : The maximum fragment length for valid paired-end alignments (default: 500)
    //  --no-discordant : A discordant alignment is an alignment where both mates align 
    //                    uniquely, but that does not satisfy the paired-end constraints.
    //                    --no-discordant only report concordant alignments.
    //  --no-mixed : By default, when bowtie2 cannot find a concordant or discordant 
    //               alignment for a pair, it then tries to find alignments for the individual 
    //               mates. This option disables that behavior.
    //  --ignore-quals : Ignore qualities of mismatches
    //  -mp 6,6 : Sets the maximum (MX) and minimum (MN) mismatch penalties to 6
    //  -np 6 : Sets penalty for positions where the read, reference, or both, contain an 
    //          ambiguous character such as N. Default: 1.
    //  --score-min <func> : Sets a function governing the minimum alignment score needed for 
    //                       an alignment to be considered "valid"
    //                       L,0,-0.6 sets the minimum-score function f to f(x) = 0 + -0.6 * x, 
    //                       where x is the read length. For example, if we set L,-18.0,0.0,
    //                       then it allows "3 (18/6)" mismatches without considering length of read.
    //  --gbar <int> : Disallow gaps within <int> positions of the beginning or end of the read. 
    //                 Default: 4. We set as 1000 for not allowing gaps.
    //

    // get bowtie path
    bowtie_option_cmd = bowtie_path;

    bowtie_option_cmd += " -t --no-unal";

    if (bowtie2_preset == "only-mismatch") {
        // calculate min_score based on number of mismatches
        int min_score = (Utils::stringToInt(maximum_mismatch_count))*(-6);
        bowtie_option_cmd += " --no-discordant --no-mixed --ignore-quals --mp 6,6 --np 6 --score-min L," + Utils::intToString(min_score) + ",0.0 --gbar 1000";
    }
    else if (bowtie2_preset == "very-fast") {
        bowtie_option_cmd += " --very-fast";
    }
    else if (bowtie2_preset == "fast") {
        bowtie_option_cmd += " --fast";
    }
    else if (bowtie2_preset == "very-sensitive") {
        bowtie_option_cmd += " --very-sensitive";
    }
    else {
        bowtie_option_cmd += " --sensitive";
    }

    // for fasta or fastq
    if (fasta_reads_flag) {    // fasta
        // for bowtie option by config file
        bowtie_option_cmd += " -f -p " + bowtie_threads_number + " -I " + minimum_fragment_length + " -X " + maximum_fragment_length;
    }
    else {    // fastq
        bowtie_option_cmd += " -q -p " + bowtie_threads_number + " -I " + minimum_fragment_length + " -X " + maximum_fragment_length;
    }

    // for read command
    if (paired_end_reads_flag) {    // for paired-end
        std::string paired_end_reads_1 = config.getPairedEndReads1();
        std::string paired_end_reads_2 = config.getPairedEndReads2();
        read_option_cmd = "-1 " + paired_end_reads_1 + " -2 " + paired_end_reads_2;
    }
    else {    // for single-end
        std::string single_end_reads = config.getSingleEndReads();
        // bowtie2 uss -U option for single-end read
        read_option_cmd = "-U " + single_end_reads;
    }

    // for samtools command
    if (paired_end_reads_flag) {    // for paired-end
        // -h: Include the header in the output. 
        // -b: bam format output
        // -S: input is SAM format
        // -F INT: Skip alignments with bits present in INT [0]
        //    0x0004  the query sequence itself is unmapped +
        //    0x0008  the mate is unmapped
        // sam format
        // samtools_cmd = samtools_path + " view -h -S -F 0x000C";
        // bam format
        samtools_cmd = samtools_path + " view -b -h -S -F 0x000C";
    }
    else {
        // 0x0004  the query sequence itself is unmapped
        // sam format
        // samtools_cmd = samtools_path + " view -h -S -F 0x0004";
        // bam format
        samtools_cmd = samtools_path + " view -b -h -S -F 0x0004";
    }

}


//=============================================================================
// Search genome directory list and save to vector 
//=============================================================================
void searchGenomeDirectoryList(
        const std::string & working_directory,                              // (in), working directory
        std::string & reference_genome_directory,                           // (in), reference genome directory
        std::vector<std::string> & genome_directory_list,                   // (out) genome directory list (vector)
        std::vector< std::vector<std::string> > & genome_fasta_path_list,   // (out) genome fasta path list (vec of vec)
        std::vector<std::string> & genome_index_base_list,                  // (out) genome index base list (vector)
        std::string & output_parent_genome_directory,                       // (out) output_parent_genome_directory (string)
        std::vector<std::string> & output_genome_directory_list)            // (out) output genome directory list (vector)
{
    // variables
    std::string genome_name, genome_directory, genome_index_base, output_genome_directory, 
                filename, fasta_path;

    // reference genome directory check (/A/B/C/ -> /A/B/C)
    if (*reference_genome_directory.rbegin() == Utils::getPathSeparator()) {
        std::string temp_string = reference_genome_directory.substr(0, reference_genome_directory.size()-1);
        reference_genome_directory= temp_string;
    }
    
    // dirent and stat 
    DIR *directory;
    struct dirent *directory_entry;


    // opendir
    directory = opendir(reference_genome_directory.c_str());

    // default output_parent_genome_directory_name
    std::string output_parent_genome_directory_name = "sigma_alignments_output";
    output_parent_genome_directory = working_directory + output_parent_genome_directory_name;

    // for output parent genome directory
    //size_t sep_pos = reference_genome_directory.rfind(Utils::getPathSeparator(), reference_genome_directory.length());
    //if (sep_pos != std::string::npos) 
    //    output_parent_genome_directory_name = reference_genome_directory.substr(sep_pos+1, reference_genome_directory.length()-sep_pos);
    output_parent_genome_directory = working_directory + output_parent_genome_directory_name;

    // typedef vector of vector
    typedef std::vector< std::vector<std::string> > vector_of_vector_string;
    
    // loop subdirectories and files
    while ((directory_entry = readdir(directory)) != NULL) {
        // get directory name (genome name)
        genome_name = directory_entry->d_name;

        // handle others
        genome_directory = reference_genome_directory + Utils::getPathSeparator() + genome_name;
        genome_index_base = genome_directory + Utils::getPathSeparator() + genome_name;
        output_genome_directory = output_parent_genome_directory + Utils::getPathSeparator() + genome_name;

        // only consider directories
        if (genome_name[0] == '.')
            continue;

        // only consider directories
        if (Utils::isDirectory(genome_directory)) {
            // push back to vector
            genome_directory_list.push_back( genome_directory );
            genome_index_base_list.push_back( genome_index_base );
            output_genome_directory_list.push_back( output_genome_directory );
            genome_fasta_path_list.push_back(vector_of_vector_string::value_type());

            // dirent and stat for sub directory
            DIR *sub_directory;
            struct dirent *sub_directory_entry;

            // opendir
            sub_directory = opendir(genome_directory.c_str());

            // loop sub_direcotoreis and files
            while ((sub_directory_entry = readdir(sub_directory)) != NULL) {

                // get filename
                filename = sub_directory_entry->d_name;

                // only consider directories
                if (genome_name[0] == '.')
                    continue;

                // only consider fasta format file
                if (Utils::isFastaFormat(filename)) {
                    // push back gnome_fasta_path to vector of vector
                    fasta_path = genome_directory + Utils::getPathSeparator() + filename;
                    genome_fasta_path_list.back().push_back(fasta_path);
                }
            }
            // close
            closedir(sub_directory);
        }
    }
    // close
    closedir(directory);

    // check validation (genome directory)
    if (genome_directory_list.size() == 0)
        Utils::exitWithError("*** Error: genome directory does not exist. Check reference directory in the config file.\n");

}


//=============================================================================
// Search remain genome directory list and save to vector 
//=============================================================================
void remainGenomeDirectoryList(
        std::string & reference_genome_directory,       // (in), reference genome directory
        std::vector<int> & remain_genome_index_list)    // (out) remain genome index list (vector)
{
    // variables
    int genome_num;
    std::string genome_name, genome_directory, genome_index_base, system_run_log_path;

    // reference genome directory check (/A/B/C/ -> /A/B/C)
    if (*reference_genome_directory.rbegin() == Utils::getPathSeparator()) {
        std::string temp_string = reference_genome_directory.substr(0, reference_genome_directory.size()-1);
        reference_genome_directory= temp_string;
    }
    
    // dirent and stat 
    DIR *directory;
    struct dirent *directory_entry;

    // opendir
    directory = opendir(reference_genome_directory.c_str());

    // loop subdirectories and files
    genome_num = 0;
    while ((directory_entry = readdir(directory)) != NULL) {
        // get directory name (genome name)
        genome_name = directory_entry->d_name;

        // handle others
        genome_directory = reference_genome_directory + Utils::getPathSeparator() + genome_name;
        genome_index_base = genome_directory + Utils::getPathSeparator() + genome_name;
        system_run_log_path = genome_index_base + ".align.log";

        // only consider directories
        if (genome_name[0] == '.')
            continue;

        // only consider directories
        if (Utils::isDirectory(genome_directory)) {
            // only consider directories that don't have system run log path file
            if (!Utils::isFileExist(system_run_log_path)) {
                remain_genome_index_list.push_back( genome_num );
            }
            genome_num++;
        }
    }
    // close
    closedir(directory);
}


//=============================================================================
// get longest read length and number of reads for each file
//=============================================================================
void getReadsInfo(const std::string & filepath,     // (in)
                  const bool & fasta_reads_flag,    // (in)
                  unsigned int & max_read_length,   // (out)
                  unsigned int & number_of_reads)   // (out)
{
    // initialize
    unsigned int read_length = 0;
    max_read_length = 0;
    number_of_reads = 0;

    // open read_file -> input_file
    std::ifstream input_file (filepath.c_str());

    if(input_file == NULL)
        Utils::exitWithError("*** Error: Unable to open file: "+filepath);

    // fasta file format
    if (fasta_reads_flag) 
    {    
        // get line
        for (std::string line; getline (input_file, line); ) { 
            line = line.erase(line.find_last_not_of(" \n\r\t")+1);
            if (line[0] == '>') {
                number_of_reads += 1;
                read_length = 0;
            }
            else {
                read_length += line.length();
                if (read_length > max_read_length)
                    max_read_length = read_length;
            }
        }
    }
    // fastq file format
    else
    {
        // get line
        std::string line;
        while(!input_file.eof())
        {
            for(unsigned int i = 0; i < 4; i++)   
            {
                getline (input_file,line);
                line = line.erase(line.find_last_not_of(" \n\r\t")+1);

                // check squence (second line of four lines)
                if (i == 1)
                    read_length = line.length();
                    if (read_length > max_read_length)
                        max_read_length = read_length;
            }
            number_of_reads += 1;
        }
    }

    // close input_file
    input_file.close();
}


//=============================================================================
// get longest read length and total number of reads
//=============================================================================
void getReadsLengthAndCount(const std::string & config_path,        // (in)
                            unsigned int & total_number_of_reads,   // (out)
                            unsigned int & longest_read_length)     // (out)
{
    // load config file
    SigmaConfig config(config_path);

    // for bowtie read option command
    bool paired_end_reads_flag, fasta_reads_flag;
    config.getReadsType(paired_end_reads_flag, fasta_reads_flag);

    total_number_of_reads = 0;
    longest_read_length = 0;
    unsigned int read_length, number_of_reads;

    // for read command
    if (paired_end_reads_flag) // for paired-end
    {    
        // get reads (comma seprator string)
        std::string paired_end_reads_1 = config.getPairedEndReads1();
        std::string paired_end_reads_2 = config.getPairedEndReads2();

        // read files string -> list
        std::vector<std::string> paired_end_reads_1_list = Utils::commaStringToVector(paired_end_reads_1); 
        std::vector<std::string> paired_end_reads_2_list = Utils::commaStringToVector(paired_end_reads_2); 

        // for read1 files
        for (std::vector<std::string>::iterator it = paired_end_reads_1_list.begin() ; it != paired_end_reads_1_list.end(); ++it)
        {
            std::string filepath = *it;

            // count
            getReadsInfo(filepath,          // (in)
                         fasta_reads_flag,  // (in)
                         read_length,       // (out)
                         number_of_reads);  // (out)

            // update total number of reads
            total_number_of_reads += number_of_reads;

            // update longest read length
            if (read_length > longest_read_length)
                longest_read_length = read_length;
        }

        // for read2 files
        for (std::vector<std::string>::iterator it = paired_end_reads_2_list.begin() ; it != paired_end_reads_2_list.end(); ++it)
        {
            std::string filepath = *it;

            // count
            getReadsInfo(filepath,          // (in)
                         fasta_reads_flag,  // (in)
                         read_length,       // (out)
                         number_of_reads);  // (out)

            // update total number of reads
            total_number_of_reads += number_of_reads;

            // update longest read length
            if (read_length > longest_read_length)
                longest_read_length = read_length;
        }
    }
    else // for single-end
    {    
        // get reads (comma seprator string)
        std::string single_end_reads = config.getSingleEndReads();

        // read files string -> list
        std::vector<std::string> single_end_reads_list = Utils::commaStringToVector(single_end_reads);

        // for read files
        for (std::vector<std::string>::iterator it = single_end_reads_list.begin() ; it != single_end_reads_list.end(); ++it)
        {
            std::string filepath = *it;

            // count
            getReadsInfo(filepath,          // (in)
                         fasta_reads_flag,  // (in)
                         read_length,       // (out)
                         number_of_reads);  // (out)

            // update total number of reads
            total_number_of_reads += number_of_reads;

            // update longest read length
            if (read_length > longest_read_length)
                longest_read_length = read_length;
        }
    }
}


//=============================================================================
// get alignment rate from bowtie alignment log file
//=============================================================================
double getAlignmentRate(const std::string & alignment_log_path, 
                        const std::string & alignment_bam_path)
{
    double alignment_rate = 0.0;

    // if bowtie2 alignment log exists
    if (Utils::isFileExist(alignment_log_path))
    {
        // get the string
        std::string target_string = "overall alignment rate";

        // open read_file -> input_file
        std::ifstream input_file (alignment_log_path.c_str());

        // get line
        for (std::string line; getline (input_file, line); ) { 
            line = line.erase(line.find_last_not_of(" \n\r\t")+1);
            // if line contains the target_string
            if (line.find(target_string) != std::string::npos) {
                // save fields to vector
                std::istringstream line_stream(line);
                std::vector<std::string> line_vector;
                for (std::string field; getline(line_stream, field, ' '); ) {
                    line_vector.push_back(field);
                }
                std::string alignment_rate_string = line_vector[0].substr(0, line_vector[0].size()-1);
                alignment_rate = Utils::stringToDouble(alignment_rate_string);
            }
        }
    }
    // if bowtie2 alignment log does not exist (provided by another alignment tool)
    // currently not support this! 
    else
    {
        Utils::exitWithError("*** Error: Alignment log file " + alignment_log_path + " does not exist!\n");
    }
        
    return alignment_rate;
}


//=============================================================================
// std::pair<std::string, double> MyPair
//=============================================================================
typedef std::pair<std::string, double> stringDoublePair;

struct CompareByKey {
    bool operator() (const stringDoublePair& a, const stringDoublePair& b) const {
    return a.first > b.first;
    };
};

struct CompareByValue {
    bool operator() (const stringDoublePair& a, const stringDoublePair& b) const {
    return a.second > b.second;
    };
};


//=============================================================================
// return common read ID for the paired_end
//=============================================================================
std::string getReadID(std::string & read_id_raw, 
                      bool & paired_end_reads_flag, 
                      bool & fasta_reads_flag) 
{
    std::string read_id = "";

    // only if paired-end and fasta format, then remove .1 and .2
    // fastq format (/1 and /2) extensioin will be removed in the alignment
    // single-end will keep the read ID
    if (paired_end_reads_flag && fasta_reads_flag)
    {
        size_t pos = read_id_raw.find_last_of(".");
        if (pos != std::string::npos)
            read_id.assign(read_id_raw.begin(),read_id_raw.begin()+pos);
        else
            Utils::exitWithError("*** Error: FASTA paired-end read ID" + read_id_raw + " does not satisfy suggested format (.1 and .2)!\n");
    }
    else 
    {    
        read_id = read_id_raw;
    }    


    return read_id;
}


//=============================================================================
// get values of SAM fields from SAM output line
//=============================================================================
void getValuesOfSAMLine(const std::string & line,           // (in)
                        std::string & read_id,              // (out)
                        unsigned int & read_flag,           // (out)
                        unsigned int & read_seq_length,     // (out)
                        unsigned int & read_mismatch_count) // (out)
{
    // line stringstream
    std::istringstream line_stream(line);

    // loop fields by tab-delimited and save to vector
    std::vector<std::string> fields_vector;
    for (std::string field; getline(line_stream, field, '\t'); ) {
        fields_vector.push_back(field);
    }

    // get read ID
    // read_id = getReadID(fields_vector[0], paired_end_reads_flag, fasta_reads_flag);
    read_id = fields_vector[0];

    // get bitwise FLAG
    read_flag = Utils::stringToUnsignedInt(fields_vector[1]);

    // get read sequence and length
    std::string read_seq = fields_vector[9];
    read_seq_length = read_seq.length();

    // get read mismatch count
    read_mismatch_count = 0;
    for (unsigned int j=11; j<fields_vector.size(); j++) 
    {
        // check the field contains 'NM'
        if (fields_vector[j].find("NM") != std::string::npos) 
        {
            std::vector<std::string> elems;
            std::stringstream ss(fields_vector[j]);
            std::string item;
            while (std::getline(ss, item, ':'))
                elems.push_back(item);
            read_mismatch_count = Utils::stringToUnsignedInt(elems.back());
        }
    }
}


//=============================================================================
// generate Q matrix
//=============================================================================
void generateQMatrix(const std::string & config_path,       // (in), config_path
                     const std::string & working_directory) // (in), working_directory
{

    // [Step 1] prepare
    std::cout << "  ** Check alignment results and prepare model: Running -> " << std::flush;

    // load config file
    SigmaConfig config(config_path);

    // get variables
    std::string reference_genome_directory = config.getReferenceGenomeDirectory();
    std::string bowtie_path = config.getBowtiePath();
    std::string samtools_path = config.getSamtoolsPath();
    double mismatch_probability = config.getMismatchProbability();
    double minimum_relative_abundance = config.getMinimumRelativeAbundance();

    // qmatrix output
    std::string qmatrix_filename = "sigma_out.qmatrix.txt";
    std::string qmatrix_filepath = working_directory + qmatrix_filename;
    std::ofstream qmatrix_out (qmatrix_filepath.c_str());

    // get paired_end_reads_flag, fasta_reads_flag
    bool paired_end_reads_flag, fasta_reads_flag;
    config.getReadsType(paired_end_reads_flag, fasta_reads_flag);

    // initialize
    std::string output_parent_genome_directory;
    std::vector<std::string> genome_directory_list, genome_index_base_list, output_genome_directory_list;
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

    // get number of reads and longest read length
    unsigned int total_number_of_reads, longest_read_length;
    getReadsLengthAndCount(config_path,             // (in)
                           total_number_of_reads,   // (out)
                           longest_read_length);    // (out)
    // if paired-end, then longest_read_length should be double
    longest_read_length = longest_read_length * 2;

    // get alignment_rate_pair
    std::string output_genome_directory, genome_name;
    std::string alignment_log_path, alignment_bam_path, alignment_sam_path;
    double alignment_rate;

    std::vector< stringDoublePair > alignment_rate_vector;

    for(std::vector<int>::size_type i = 0; i != output_genome_directory_list.size(); i++) 
    {
        // get genome name
        output_genome_directory = output_genome_directory_list[i];
        genome_name = Utils::getFilename(output_genome_directory);

        // get alignment rate and insert to map
        alignment_log_path = output_genome_directory + Utils::getPathSeparator() + genome_name + ".align.log";
        alignment_bam_path = output_genome_directory + Utils::getPathSeparator() + genome_name + ".align.bam";
        alignment_rate = getAlignmentRate(alignment_log_path, alignment_bam_path);
        alignment_rate_vector.push_back(std::make_pair(output_genome_directory, alignment_rate));
    }

    // sort by value
    std::sort(alignment_rate_vector.begin(), alignment_rate_vector.end(), CompareByValue());

    // only consider genomes if genome alignment_rate > minimum_relative_abundance
    std::vector<std::string> aligned_genome_list;
    typedef std::vector< stringDoublePair >::iterator it_type;
    for(it_type iterator = alignment_rate_vector.begin(); iterator != alignment_rate_vector.end(); iterator++) 
    {
        output_genome_directory = iterator->first;
        alignment_rate = iterator->second;

        // only consider genomes if genome alignment_rate > minimum_relative_abundance
        if (alignment_rate >= minimum_relative_abundance) {
            aligned_genome_list.push_back(output_genome_directory);
        }
    }

    // [Step 1] Done
    std::cout << "Done!" << std::endl;

    // [Step 2] Build model
    std::cout << "  ** Build a probabilistic qmatrix model: Running -> " << std::flush;


    // initizlize map
    std::string system_cmd;
    typedef std::map<unsigned int, std::vector<unsigned int> > map_of_vector_unsigned_int;
    typedef std::map<unsigned int, std::vector<double> > map_of_vector_double;
    map_of_vector_unsigned_int aligned_genome_index_map;    // to save genome index for each read
    map_of_vector_double aligned_qvalue_map;                // to save qvalue for each read
    std::map<std::string, unsigned int> reads_index_map;    // to save unique read IDs

    // initialize variables
    unsigned int unique_read_index = 0; // for reads_index_map value
    std::string read_id, read2_id;      // read id for paired-ends
    unsigned int read_flag, read_seq_length, read_mismatch_count;       // for read
    unsigned int read2_flag, read2_seq_length, read2_mismatch_count;    // for the mate-pair read
    unsigned int read_index = 0;        // for read index
    int read_match_count_scale;         // number of matches is scaled 
    double qvalue_mismatch, qvalue_match, qvalue;   // qvalue

    // loop genomes that have at least minimum_relative_abundance rate
    // then save data (genome index, qvalue) to structure
    for (unsigned int i=0; i<aligned_genome_list.size(); i++) 
    {
        // get output genome directory and genome name
        output_genome_directory = aligned_genome_list[i];
        genome_name = Utils::getFilename(output_genome_directory);

        // get bam and sam filepath
        alignment_bam_path = output_genome_directory + Utils::getPathSeparator() + genome_name + ".align.bam";
        alignment_sam_path = output_genome_directory + Utils::getPathSeparator() + genome_name + ".align.sam";

        // convert BAM to SAM
        if (!Utils::isFileExist(alignment_sam_path))
        {
            // system command
            system_cmd = samtools_path + " view -h " + alignment_bam_path + " -o " + alignment_sam_path;

            // system call
            const int res = system(system_cmd.c_str());
            if (res!=0)
                Utils::exitWithError("*** Error: Failed to convert BAM to SAM. Command: " + system_cmd);
        }

        // open SAM file -> input_file
        std::ifstream input_file (alignment_sam_path.c_str());
     
        // get line
        for (std::string line; getline (input_file, line); ) 
        {
            // if not header lines
            if (line[0] != '@') 
            {
                // get values of SAM fields
                getValuesOfSAMLine(line,                    // (in)
                                   read_id,                 // (out)
                                   read_flag,               // (out)
                                   read_seq_length,         // (out)
                                   read_mismatch_count);    // (out)

                // if paired-end, then load second line
                if (paired_end_reads_flag) {
                    // get next line (mate pair)
                    getline (input_file, line);
                    getValuesOfSAMLine(line,                    // (in)
                                       read2_id,                 // (out)
                                       read2_flag,               // (out)
                                       read2_seq_length,         // (out)
                                       read2_mismatch_count);    // (out)
                    // update
                    read_id = getReadID(read2_id, paired_end_reads_flag, fasta_reads_flag);
                    read_seq_length += read2_seq_length;
                    read_mismatch_count += read2_mismatch_count;
                }

                // qvalue = (0.95)^(# of matches - longest)*(0.5)^(# of mismatches)
                read_match_count_scale = read_seq_length - read_mismatch_count - longest_read_length;
                qvalue_match = pow(1.0 - mismatch_probability, read_match_count_scale);
                qvalue_mismatch = pow(mismatch_probability, read_mismatch_count);
                qvalue = qvalue_match*qvalue_mismatch;

                // update reads_index_map by key->read ID, val->unique_read_index
                // if key exists
                if (reads_index_map.find(read_id) != reads_index_map.end()) {
                    read_index = reads_index_map.find(read_id) -> second;
                }
                else {
                    reads_index_map.insert(std::pair<std::string, unsigned int>(read_id, unique_read_index));
                    read_index = unique_read_index;
                    unique_read_index += 1;
                }

                // insert map 
                aligned_genome_index_map[read_index].push_back(i);
                aligned_qvalue_map[read_index].push_back(qvalue);
            }
        }
       

        // close input_file
        input_file.close();
        
        // delete SAM file
        if (Utils::isFileExist(alignment_sam_path))
        {
            // system command
            system_cmd = "rm -f " + alignment_sam_path;

            // system call
            const int res = system(system_cmd.c_str());
            if (res!=0)
                Utils::exitWithError("*** Error: Failed to delete SAM file. Command: " + system_cmd);
        }
    }

    // count number of mapped reads
    unsigned int number_of_mapped_reads = reads_index_map.size();
    // if paired-end, then multiply 2
    if (paired_end_reads_flag) 
        number_of_mapped_reads = number_of_mapped_reads * 2;
    unsigned int number_of_unmapped_reads = total_number_of_reads - number_of_mapped_reads;

    // precision
    std::streamsize qmatrix_precision = qmatrix_out.precision(15);

    // write comments
    qmatrix_out << "#\t+\tMatrixName\tTotalNumberReads\tNumberMappedReads\tNumberUnmappedReads" << std::endl;
    qmatrix_out << "#\t@\tGenomeIndex\tGenomeName\tPercentageAlignmentRate" << std::endl;
    qmatrix_out << "#\t*\tReadID\tGenomeIndex=QValue" << std::endl;

    // print qmatrix name
    qmatrix_out << "+\t" << qmatrix_filename << "\t" 
                << Utils::unsignedIntToString(total_number_of_reads) << "\t"
                << Utils::unsignedIntToString(number_of_mapped_reads) << "\t"
                << Utils::unsignedIntToString(number_of_unmapped_reads) << std::endl;

    // only consider genomes if genome alignment_rate > minimum_relative_abundance
    unsigned int genome_enum = 0;
    for(it_type iterator = alignment_rate_vector.begin(); iterator != alignment_rate_vector.end(); iterator++) 
    {
        output_genome_directory = iterator->first;
        alignment_rate = iterator->second;
        genome_name = Utils::getFilename(output_genome_directory);

        // only consider genomes if genome alignment_rate > minimum_relative_abundance
        if (alignment_rate >= minimum_relative_abundance) {
            qmatrix_out << "@\t" << Utils::unsignedIntToString(genome_enum) << "\t"
                        << genome_name << "\t" 
                        << Utils::doubleToString(alignment_rate) << std::endl;
            genome_enum++;
        }
    }


    // write qvalue for each aligned read
    unsigned int genome_index;
    for(std::map<std::string, unsigned int>::iterator iter=reads_index_map.begin(); iter!=reads_index_map.end(); iter++)
    {
        // get key: read_id
        read_id = iter->first;
        // get value: unique_read_index
        unique_read_index = iter->second;

        // write read_id
        qmatrix_out << "*\t" << read_id;
        
        // loop columns of data
        for (unsigned int j=0; j<aligned_genome_index_map[unique_read_index].size(); j++ ) 
        {
            // get qvalue
            genome_index = aligned_genome_index_map[unique_read_index][j];
            qvalue = aligned_qvalue_map[unique_read_index][j];
            qmatrix_out << "\t" << genome_index << "=" << qvalue;
        }
        // end line
        qmatrix_out << std::endl;
    }

    // close
    qmatrix_out.precision(qmatrix_precision);
    qmatrix_out.close();

    // clear
    aligned_genome_index_map.clear();
    aligned_qvalue_map.clear();
    reads_index_map.clear();

    // [Step 2] Done
    std::cout << "Done!" << std::endl;
    
}


//=============================================================================
// Write Sigma output to HTML format 
//=============================================================================
void writeHTMLOutput(const std::string & gvector_filepath)              // (in)
{
    // variables
    unsigned int genome_index;
    std::string genome_name, gvector_filebase, html_output_filename;
    double align_rate, pct_chance, pct_scaled;

    // map
    std::vector<std::string> info_list;
    std::map<unsigned int, std::string> gvector_genome_name_map;
    std::map<unsigned int, double> gvector_align_rate_map, gvector_pct_chance_map, gvector_pct_scaled_map;

    // open gvector -> input_file
    std::ifstream input_file (gvector_filepath.c_str());

    // for HTML output
    gvector_filebase = Utils::getFilebase2(gvector_filepath);
    html_output_filename = gvector_filebase + ".html";
    std::ofstream html_output (html_output_filename.c_str());

    // get gvector data
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

        // line with + 
        if (line [0] == '+') 
        {
            for (unsigned int i=0; i<fields_vector.size(); i++) {
                info_list.push_back(fields_vector[i]);
            }
        }
        // line with @ [GenomeIndex    GenomeName    PercentageAlignmentRate]
        else if (line [0] == '@') 
        {
            genome_index = Utils::stringToUnsignedInt(fields_vector[1]);
            genome_name  = fields_vector[2];
            align_rate   = Utils::stringToDouble(fields_vector[3]);

            // save to map
			gvector_genome_name_map.insert(std::pair<unsigned int, std::string>(genome_index, genome_name));
            gvector_align_rate_map.insert(std::pair<unsigned int, double>(genome_index, align_rate));
        }
        // line with * [GenomeIndex    PercentageChance    ScaledPercentageChance]
        else if (line [0] == '*') 
        {
            genome_index = Utils::stringToUnsignedInt(fields_vector[1]);
            pct_chance   = Utils::stringToDouble(fields_vector[2]);
            pct_scaled   = Utils::stringToDouble(fields_vector[3]);

            // save to map
            gvector_pct_chance_map.insert(std::pair<unsigned int, double>(genome_index, pct_chance));
            gvector_pct_scaled_map.insert(std::pair<unsigned int, double>(genome_index, pct_scaled));
        }
    }

    const std::string text1(
        "<html>\n"
        "  <head>\n"
        "    <!--Load the AJAX API-->\n"
        "    <script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>\n"
        "    <script type=\"text/javascript\">\n"
        "\n"
        "      // Load the Visualization API and the packages.\n"
        "      google.load('visualization', '1.0', {packages:['table']}); \n"
        "      google.load('visualization', '1.0', {packages:['corechart']});\n"
        "\n"
        "      // Set a callback to run when the Google Visualization API is loaded.\n"
        "      google.setOnLoadCallback(drawInfoTable);\n"
        "      google.setOnLoadCallback(drawAlignmentTable);\n"
        "      google.setOnLoadCallback(drawSigmaTableChart);\n"
        "      google.setOnLoadCallback(drawSigmaTableChart2);\n"
        "\n"
        "      // Callback that creates and populates a data table, instantiates the table,\n"
        "      // passes in the data and draws it.\n"
        "      function drawInfoTable() \n"
        "      {\n"
        "        // Create the data table\n"
        "        var data = new google.visualization.DataTable();\n"
        "        data.addColumn('string', 'Gvector Name');\n"
        "        data.addColumn('number', 'Total Number of Reads');\n"
        "        data.addColumn('number', 'Number of Mapped Reads');\n"
        "        data.addColumn('number', 'Number of Unmapped Reads');\n"
        "        data.addRows(1);\n"
    );
    html_output << text1 
                << "        data.setCell(0, 0, '" << info_list[1] << "');" << std::endl
                << "        data.setCell(0, 1, " << info_list[2] << ");" << std::endl
                << "        data.setCell(0, 2, " << info_list[3] << ");" << std::endl
                << "        data.setCell(0, 3, " << info_list[4] << ");" << std::endl << std::endl;
    const std::string text2(
        "        // Instantiate and draw the table\n"
        "        var table = new google.visualization.Table(document.getElementById('drawInfoTable_div'));\n"
        "        table.draw(data);\n"
        "      }\n"
        "\n"
        "      // Callback that creates and populates a data table, instantiates the table,\n"
        "      // passes in the data and draws it.\n"
        "      function drawAlignmentTable()\n"
        "      {\n"
        "        // Create the data table\n"
        "        var data = new google.visualization.DataTable();\n"
        "        data.addColumn('string', 'Genome Index');\n"
        "        data.addColumn('string', 'Genome Name');\n"
        "        data.addColumn('number', 'Percentage Alignment Rate');\n"
    );
    html_output << text2 
                << "        data.addRows(" << gvector_genome_name_map.size() << ");" << std::endl;
    for (unsigned int i=0; i<gvector_genome_name_map.size(); i++) {
        genome_index = i;
        genome_name = gvector_genome_name_map.find(genome_index) -> second;
        align_rate = gvector_align_rate_map.find(genome_index) -> second;
        html_output << "        data.setCell(" << i << ", 0, '" << genome_index << "');" << std::endl;
        html_output << "        data.setCell(" << i << ", 1, '" << genome_name << "');" << std::endl;
        html_output << "        data.setCell(" << i << ", 2, " << align_rate << ");" << std::endl;
    }
    const std::string text3(
        "\n"
        "        // % Format\n"
        "        var formatter = new google.visualization.NumberFormat({prefix: '%'});\n"
        "        formatter.format(data, 2); // Apply formatter to second column\n"
        "\n"
        "        // Instantiate and draw the table\n"
        "        var table = new google.visualization.Table(document.getElementById('drawAlignmentTable_div'));\n"
        "        table.draw(data);\n"
        "      }\n"
        "\n"
        "      // Callback that creates and populates a data table, instantiates the table,\n"
        "      // passes in the data and draws it.\n"
        "      function drawSigmaTableChart()\n"
        "      {\n"
        "        // Create the data table\n"
        "        var data = new google.visualization.DataTable();\n"
        "        data.addColumn('string', 'Genome Index : Genome Name');\n"
        "        data.addColumn('number', 'Relative Abundance Estimation');\n"
        "        data.addColumn('number', 'Percentage Scaled Estimation');\n"
    );
    html_output << text3 
                << "        data.addRows(" << gvector_pct_chance_map.size() << ");" << std::endl;
    unsigned int i=0;
    typedef std::map<unsigned int, double>::iterator it_type;
    for(it_type iter = gvector_pct_chance_map.begin(); iter != gvector_pct_chance_map.end(); iter++) {
        genome_index = iter->first;
        genome_name = gvector_genome_name_map.find(genome_index) -> second;
        pct_chance = gvector_pct_chance_map.find(genome_index) -> second;
        pct_scaled = gvector_pct_scaled_map.find(genome_index) -> second;
        html_output << "        data.setCell(" << i << ", 0, '" << genome_index << ":" << genome_name << "');" << std::endl;
        html_output << "        data.setCell(" << i << ", 1, " << pct_scaled << ");" << std::endl;
        html_output << "        data.setCell(" << i << ", 2, " << pct_chance << ");" << std::endl;
        i++;
    }
    const std::string text4(
        "\n"
        "        // % Format\n"
        "        var formatter = new google.visualization.NumberFormat({prefix: '%'});\n"
        "        formatter.format(data, 2); // Apply formatter to second column\n"
        "\n"
        "        // Instantiate and draw the table\n"
        "        var table = new google.visualization.Table(document.getElementById('drawSigmaTable_div'));\n"
        "        table.draw(data);\n"
        "\n"
        "        // Var view for pie chart\n"
        "        var view = new google.visualization.DataView(data);\n"
        "        view.setColumns([0, 1]);\n"
        "\n"
        "        // Set chart options\n"
        "        var options = {'title':'Genome Chance by Sigma',\n"
        "                       'titleTextStyle':{fontSize: 20},\n"
        "                       'width':900,\n"
        "                       'height':700};\n"
        "\n"
        "        // Instantiate and draw our chart, passing in some options.\n"
        "        var chart = new google.visualization.PieChart(document.getElementById('drawSigmaChart_div'));\n"
        "        chart.draw(view, options);\n"
        "      }\n"
        "\n"
        "    </script>\n"
        "  </head>\n"
        "\n"
        "  <body style=\"font-family: Arial;border: 0 none;\">\n"
        "    <p>\n"
        "    <H1>&nbsp;&nbsp;&nbsp;&nbsp;Sigma Output</H1><BR>\n"
        "    </p>\n"
        "\n"
        "    <!--Div that will hold the tables and pie chart-->\n"
        "    <p>\n"
        "    <H3>&nbsp;&nbsp;&nbsp;&nbsp;Info</H3>\n"
        "    <div id='drawInfoTable_div' style=\"width: 800px;\"></div>\n"
        "    <BR><BR>\n"
        "    </p>\n"
        "\n"
        "    <p>\n"
        "    <H3>&nbsp;&nbsp;&nbsp;&nbsp;Alignments</H3>\n"
        "    <div id='drawAlignmentTable_div' style=\"width: 800px;\"></div>\n"
        "    <BR><BR>\n"
        "    </p>\n"
        "\n"
        "    <p>\n"
        "    <H3>&nbsp;&nbsp;&nbsp;&nbsp;Sigma Results</H3>\n"
        "    <div id='drawSigmaTable_div' style=\"width: 800px;\"></div>\n"
        "    <div id='drawSigmaChart_div'></div>\n"
        "    </p>\n"
        "\n"
        "  </body>\n"
        "</html>\n"
    );
    html_output << text4;

    // close
    html_output.close();
}

