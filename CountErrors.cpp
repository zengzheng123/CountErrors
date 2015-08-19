//
//  CountErrors.cpp
//  CountErrors
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 1/10/14.
//  Copyright (c) 2014 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <algorithm>
#include "api/BamReader.h"

//#define _DEBUG

//#define _PARSING_DEBUG
//#define _FASTA_DEBUG

using namespace std;
using namespace BamTools;


const string VERSION = "CountErrors 1.1.1";

string input_fasta_file;
string input_bam_file;
string input_target_file;
string output_file;
int mapping_quality_threshold = 30;
int base_quality_threshold = 20;
int quality_scale = 33;
bool no_duplicate = false;
bool collapse_strand = false;
bool collapse_end = false;
bool collapse_context = false;


void printUsage(string msg = "")
{
    cout << endl;
    cout << VERSION << endl;
    cout << "Usage: " << endl;
    cout << "[REQUIRED ARGUMENTS]" << endl;
    cout << "\t--fasta                 <string>                        Input reference sequence file" << endl;
    cout << "\t--bam                   <string>                        Input bam file" << endl;
    cout << "\t--target                <string>                        Input target interval bed file" << endl;
    cout << "\t--output                <string>                        Output file" << endl;
    cout << endl;
    cout << "[OPTIONAL ARGUMENTS]" << endl;
    cout << "\t--maq                   <int>                           Mapping quality threshold, Default [30]" << endl;
    cout << "\t--baq                   <int>                           Base quality threshold, Default [20]" << endl;
    cout << "\t--nodup                                                 Do not use reads that are marked as duplicate" << endl;
    cout << "\t--collapse_strand                                       Collapse the strand when counting errors" << endl;
    cout << "\t--collapse_end                                          Collapse the read end when counting errors" << endl;
    cout << "\t--collapse_context                                      Collapse the context bases(left and right bases of the reference) when counting errors" << endl;
    cout << "\t--help                                                  Print command line usage" << endl;
    cout << endl;
    if(!msg.empty())
        cerr << msg << endl;
    exit(1);
}


static struct option long_options[] =
{
    {"fasta",                   required_argument,      0,     'f'},
    {"bam",                     required_argument,      0,     'b'},
    {"target",                  required_argument,      0,     't'},
    {"output",                  required_argument,      0,     'o'},
    {"maq",                     required_argument,      0,     'Q'},
    {"baq",                     required_argument,      0,     'q'},
    {"nodup",                   no_argument,            0,     'd'},
    {"collapse_strand",         no_argument,            0,     '1'},
    {"collapse_end",            no_argument,            0,     '2'},
    {"collapse_context",        no_argument,            0,     '3'},
    {"help",                    no_argument,            0,     'h'},
    {0, 0, 0, 0}
};


void parseOption(int argc, const char* argv[])
{
    if(argc == 1)
        printUsage();
    int next_option;
    int option_index = 0;
    do
    {
        next_option = getopt_long(argc, const_cast<char**>(argv), "f:b:t:o:Q:q:d123h", long_options, &option_index);
        switch(next_option)
        {
            case 'f':
                input_fasta_file = optarg;
                break;
            case 'b':
                input_bam_file = optarg;
                break;
            case 't':
                input_target_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'Q':
                mapping_quality_threshold = atoi(optarg);
                break;
            case 'q':
                base_quality_threshold = atoi(optarg);
                break;
            case 'd':
                no_duplicate = true;
                break;
            case '1':
                collapse_strand = true;
                break;
            case '2':
                collapse_end = true;
                break;
            case '3':
                collapse_context = true;
                break;
            case 'h':
                printUsage();
                break;
            case -1:
                break; //parsed all options
            default:
                printUsage(string("[ERROR] Argument error: ") + argv[optind - 1]);
        }
    } while(next_option != -1);
    if(input_fasta_file.empty())
        printUsage("[ERROR] Please specify input fasta file");
    if(input_bam_file.empty())
        printUsage("[ERROR] Please specify input bam file");
    if(input_target_file.empty())
        printUsage("[ERROR] please specify input target interval bed file");
    if(output_file.empty())
        printUsage("[ERROR] Please specify output file1");
    base_quality_threshold += quality_scale;
#ifdef _DEBUG
    cerr << "[DEBUG] Parsing options complete." << endl;
#endif
}


void split(const string& line, char delim, vector<std::string>& parsed_item, bool ignore_empty_item = false)
{
    stringstream my_ss(line);
    string item;
    while (getline(my_ss, item, delim))
    {
#ifdef _PARSING_DEBUG
        cerr << "[DEBUG] Parsed item: " << item << endl;
#endif
        if (ignore_empty_item && item.empty())
            continue;    // use to skip empty item
        parsed_item.push_back(item);
    }
}


class tri_nucleo_entry
{
public:
    tri_nucleo_entry() {}
    
    tri_nucleo_entry(char _ref, char _alt, char _left, char _right, char _strand, int _end): ref(_ref), alt(_alt), left(_left), right(_right), strand(_strand), end(_end)
    {}
    
    ~tri_nucleo_entry() {}
    
    bool operator<(const tri_nucleo_entry& other) const
    {
        if(ref != other.ref)
            return ref < other.ref;
        if(alt != other.alt)
            return alt < other.alt;
        if(left != other.left)
            return left < other.left;
        if(right != other.right)
            return right < other.right;
        if(end != other.end)
        	return end < other.end;
        return strand < other.strand;
    }
    
    friend std::ostream& operator<< (ostream& out_fs, const tri_nucleo_entry& my_entry)
    {
    	out_fs << my_entry.ref << "\t" << my_entry.alt;
        if(my_entry.left != '*')
            out_fs << "\t" << my_entry.left;
        if(my_entry.right != '*')
            out_fs << "\t" << my_entry.right;
        if(my_entry.end != -1)
            out_fs<< "\tR" << my_entry.end;
        if(my_entry.strand != '*')
            out_fs  << "\t" << my_entry.strand;
    	return out_fs;
    }
    
    char ref;
    char alt;
    char left;
    char right;
    char strand;
    int end;
};


void create_hash_entry(map<tri_nucleo_entry, int>& tri_nucleo_count)
{
    char nucleo[] = {'A', 'C', 'G', 'T'};
    size_t nucleo_size = 4;
    char strand[] = {'+', '-'};
    size_t strand_size = 2;
    int end[] = {1, 2};
    size_t end_size = 2;
    for(size_t ref_index = 0; ref_index < nucleo_size; ref_index ++)
        for(size_t alt_index = 0; alt_index < nucleo_size; alt_index ++)
            for(size_t left_index = 0; left_index < nucleo_size; left_index ++)
                for(size_t right_index = 0; right_index < nucleo_size; right_index ++)
                    for(size_t strand_index = 0; strand_index < strand_size; strand_index ++)
                        for(size_t end_index = 0; end_index < end_size; end_index ++)
                        {
                            tri_nucleo_entry new_entry(nucleo[ref_index], nucleo[alt_index], nucleo[left_index], nucleo[right_index], strand[strand_index], end[end_index]);
                            tri_nucleo_count.insert(make_pair(new_entry, 0));
                        }
}


void output_reference_sequence(string output_fastafile, map<string, string>& ref_seq, vector<string>& orignal_header)
{
    size_t base_per_line = 80;
    ofstream out_fs(output_fastafile.c_str());
    for(size_t i = 0; i < orignal_header.size(); i++)
    {
        map<string, string>::iterator it = ref_seq.find(orignal_header[i]);
        size_t chrom_len = it->second.length();
#ifdef _FASTA_DEBUG
        cerr << "[DEBUG] output reference sequence: " << it->first << ": " << chrom_len << endl;
#endif
        size_t current_len = 0;
        out_fs << ">" << it->first << endl;
        while(current_len < chrom_len)
        {
            size_t output_len = min(base_per_line, chrom_len - current_len);
            out_fs << it->second.substr(current_len, output_len) << endl;
            current_len += output_len;
        }
    }
    out_fs.close();
}


void load_reference_sequence_speedup(string fasta_filename, map<string, string>& ref_seq, vector<string>* orignal_header = NULL)
{
    cout << "Loading reference sequence: " << fasta_filename << endl;
    ifstream ref_fs(fasta_filename.c_str());
    if(!ref_fs)
    {
        cerr << "[ERROR] fail to open reference fasta file: " << fasta_filename << endl;
        exit(1);
    }
    string fasta_index_filename = fasta_filename + ".fai";
    ifstream index_fs(fasta_index_filename.c_str());
    if(!index_fs)
    {
        cerr << "[ERROR] fail to open reference fasta index file: " << fasta_index_filename << endl;
        exit(1);
    }
    string line;
    while(getline(index_fs, line))
    {
        vector<string> index_items;
        split(line, '\t', index_items);
        
        string chrom_name = index_items[0];
        int chrom_len = atoi(index_items[1].c_str());
        long long int chrom_offset = atoll(index_items[2].c_str());
        int base_per_line = atoi(index_items[3].c_str());
        int byte_per_line = atoi(index_items[4].c_str());
        int byte_len = chrom_len + (chrom_len / base_per_line) * (byte_per_line - base_per_line);
        if(orignal_header != NULL)
            orignal_header->push_back(chrom_name);
#ifdef _FASTA_DEBUG
        cerr << "[DEBUG] index entry: chrom_name:" << chrom_name << "\tchrom_len:" << chrom_len << "\tchrom_offset:" << chrom_offset << "\tbase:" << base_per_line << "\tbyte" << byte_per_line << endl;
#endif
        string new_seq;
        ref_seq.insert(make_pair(chrom_name, new_seq));
        ref_seq[chrom_name].resize(chrom_len);
#ifdef _FASTA_DEBUG
        cerr << "[DEBUG] ref_seq size has been resized to " << ref_seq[chrom_name].length() << endl;
#endif
        ref_fs.seekg(chrom_offset);
        char* seq_buff = new char[byte_len];
        ref_fs.read(seq_buff, byte_len);
#ifdef _FASTA_DEBUG
        if(ref_fs.gcount() != byte_len)
            cerr << "[DEBUG] gcount does not equal to bytelen " << ref_fs.gcount() << " != " << byte_len << endl;
        cerr << "[DEBUG] last char in the buffer is: " << seq_buff[byte_len - 1] << endl;
#endif
        string::iterator it_target = ref_seq[chrom_name].begin();
        char* it_source = seq_buff;
        for(int i = 0; i < byte_len; i++)
        {
            if(!isspace(*it_source))
            {
                *it_target = toupper(*it_source);
                it_target ++;
            }
            it_source ++;
        }
        delete[] seq_buff;
#ifdef _FASTA_DEBUG
        cout << "[DEBUG] reference sequence loaded for: " << chrom_name << endl;
#endif
    }
    ref_fs.close();
    index_fs.close();
}


void print_count(map<tri_nucleo_entry, int>& tri_nucleo_count)
{
#ifdef _DEBUG
    cerr << "[DEBUG]: printing count to output file: " << output_file << endl;
#endif
    ofstream output_fs(output_file.c_str());
    output_fs << "REF\tALT";
    if(!collapse_context)
        output_fs << "\tLEFT\tRIGHT";
    if(!collapse_end)
        output_fs << "\tREAD";
    if(!collapse_strand)
        output_fs << "\tSTRAND";
    output_fs << "\tCOUNT" << endl;
    for(map<tri_nucleo_entry, int>::iterator it = tri_nucleo_count.begin(); it != tri_nucleo_count.end(); it++)
        if(it->second > 0)
            output_fs << it->first << "\t" << it->second << endl;
    output_fs.close();
}


void CountErrors()
{
    map<string, string> reference_sequence;
    BamReader my_bam_reader;
    if(!my_bam_reader.Open(input_bam_file))
    {
        cerr << "[ERROR] fail to open input bam file: " << input_bam_file << endl;
        exit(1);
    }
    string input_bam_index_file1 = input_bam_file.substr(0, input_bam_file.length() - 3) + "bai";
    string input_bam_index_file2 = input_bam_file + ".bai";
    if(!my_bam_reader.OpenIndex(input_bam_index_file1))
    {
        if(!my_bam_reader.OpenIndex(input_bam_index_file2))
        {
            cerr << "[ERROR] fail to open input bam index file: " << input_bam_index_file1 << ", or " << input_bam_index_file2 << endl;
            exit(1);
        }
    }
    ifstream target_fs(input_target_file.c_str());
    if(!target_fs)
    {
        cerr << "[ERROR] fail to open input target file: " << input_target_file << endl;
        exit(1);
    }
    load_reference_sequence_speedup(input_fasta_file, reference_sequence);
    map<tri_nucleo_entry, int> tri_nucleo_count;
    create_hash_entry(tri_nucleo_count);
    cout << "Processing bam file: " << input_bam_file << endl;
    string line;
    while(getline(target_fs, line))
    {
        if(line.substr(0, 5) == "track" || line.substr(0, 7) == "browser")
            continue;
        vector<string> target_item;
        split(line, '\t', target_item);
        string target_chrom = target_item[0];
        int target_start = atoi(target_item[1].c_str());
        int target_end = atoi(target_item[2].c_str());
        int target_chrom_id = my_bam_reader.GetReferenceID(target_chrom);
        cout << "Processing target region: " << target_chrom << " " << target_start << " " << target_end << endl;
#ifdef _DEBUG
        cerr << "[DEBUG] Processing target region: " << target_chrom << " " << target_start << " " << target_end << endl;
#endif
        my_bam_reader.SetRegion(target_chrom_id, target_start, target_chrom_id, target_end);
        BamAlignment my_bam_alignment;
        while(my_bam_reader.GetNextAlignment(my_bam_alignment))
        {
            if(my_bam_alignment.MapQuality < mapping_quality_threshold || (no_duplicate && my_bam_alignment.IsDuplicate()) )
                continue;
#ifdef _DEBUG
            cerr << "[DEBUG] Processing bam entry: " << my_bam_alignment.Name << " " << my_bam_alignment.AlignmentFlag << " " << my_bam_reader.GetReferenceData()[my_bam_alignment.RefID].RefName \
            << " " << my_bam_alignment.Position << " ";
            for(size_t i = 0; i < my_bam_alignment.CigarData.size(); i++)
                cerr <<  my_bam_alignment.CigarData[i].Length <<  my_bam_alignment.CigarData[i].Type;
            cerr << " " << my_bam_alignment.MapQuality << " " << my_bam_alignment.QueryBases << " " << my_bam_alignment.Qualities << endl;
#endif
            int cur_end;
            if(collapse_end)
                cur_end = -1;
            else
                cur_end = my_bam_alignment.IsFirstMate() ? 1 : 2;
            char cur_strand;
            if(collapse_strand)
                cur_strand = '*';
            else
                cur_strand = my_bam_alignment.IsReverseStrand() ? '-' : '+';
            int cur_base = 0;
            int cur_pos = my_bam_alignment.Position;
            for(size_t i = 0; i < my_bam_alignment.CigarData.size(); i++)
            {
                CigarOp& cur_cigar =my_bam_alignment.CigarData[i];
                if(cur_cigar.Type == 'M')
                {
                    for(int j = 0; j < cur_cigar.Length; j++)
                    {
                        if(cur_pos >= target_start && cur_pos < target_end && my_bam_alignment.Qualities[cur_base] >= base_quality_threshold) //target region is half open region
                        {
                            char cur_ref = reference_sequence[target_chrom][cur_pos];
                            char cur_alt = my_bam_alignment.QueryBases[cur_base];
                            char cur_left;
                            char cur_right;
                            if(collapse_context)
                            {
                                cur_left = '*';
                                cur_right = '*';
                            }
                            else
                            {
                                cur_left = reference_sequence[target_chrom][cur_pos - 1];
                                cur_right = reference_sequence[target_chrom][cur_pos + 1];
                            }
                            tri_nucleo_entry new_entry(cur_ref, cur_alt, cur_left, cur_right, cur_strand, cur_end);
#ifdef _DEBUG
                            cerr << "[DEBUG] " << new_entry << "\tcur_base=" << cur_base << "\tcur_pos=" << cur_pos << "\tcur_qual=" << my_bam_alignment.Qualities[cur_base] << endl;
#endif
                            tri_nucleo_count[new_entry] ++;
                        }
                        cur_pos ++;
                        cur_base ++;
                    }
                }
                else if(cur_cigar.Type == 'N' || cur_cigar.Type == 'D')
                {
                    cur_pos += cur_cigar.Length;
                }
                else if(cur_cigar.Type == 'I')
                {
                    cur_base += cur_cigar.Length;
                }
                else if(cur_cigar.Type == 'S')
                {
                    cur_base += cur_cigar.Length;
                }
                else if(cur_cigar.Type == 'H')
                {
                    // do nothing
                }
            }
        }
    }
    print_count(tri_nucleo_count);
    target_fs.close();
    my_bam_reader.Close();
}

int main(int argc, const char * argv[])
{
    parseOption(argc, argv);
    CountErrors();
    return 0;
}


























