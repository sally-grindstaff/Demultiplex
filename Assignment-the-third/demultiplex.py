#!/usr/bin/env python

import argparse
import gzip
import Bioinfo
import re
import csv

def get_args():
    parser = argparse.ArgumentParser(description='Takes four FASTQ files (read 1, read 2, index 1, index 2) made from a paired-end, dual-matched indexing library prep as input. \
    Creates read 1 and read 2 FASTQ files for each correctly-matched pair of indexes. Also creates a  pair of read 1 and read 2 FASTQ files containing index-swapped reads and a \
    pair of read 1 and read 2 FASTQ files containing reads with low-quality or unknown indexes. Outputs a file containing statistics about how many read-pairs are successfully matched, \
    how many read-pairs come from each index/sample, how many read-pairs display index swapping, and how many read-pairs have unknown or low-quality indexes.')
    parser.add_argument('-r1', '--read1', help='Filepath to read 1 input file. Must be a zipped FASTQ file.', required=True)
    parser.add_argument('-r2', '--read2', help='Filepath to read 2 input file. Must be a zipped FASTQ file.', required=True)
    parser.add_argument('-i1', '--index1', help='Filepath to index 1 input file. Must be a zipped FASTQ file.', required=True)
    parser.add_argument('-i2', '--index2', help='Filepath to index 2 input file. Must be a zipped FASTQ file.', required=True)
    parser.add_argument('-b', '--barcodes', help='Name of file containing all intended barcode (aka index) sequences. Must be tab-separated text file with the following columns in the \
    following order: sample, group, treatment, index, index sequence. First row should have column headers, and all other rows should have index/sample information. There should be no \'/\' anywhere in this file', required=True)
    parser.add_argument('-o', '--output_dir', help='Path to directory in which output files should be written. Must already exist. User must include forward slash at the end of the directory name.', required=True)
    return parser.parse_args()

args = get_args()

#Might be nice/necessary to check that the input files are really FASTQs.
#Maybe should use try/except so that both zipped and unzipped FASTQs can be used as input.

barcode_dict = {}
#keys are index sequences
#values are lists where first item is barcode name, second item is sample name, third item is group name, and fourth item is treatment name.

#parse barcodes file to populate barcode_dict:
with open(args.barcodes, 'r') as fh_b:
    is_first_line = True
    for line in fh_b:
        if is_first_line == True:
            is_first_line = False
            continue
        line = line.strip()
        temp_list = line.split('\t')
        barcode_dict[temp_list[4]] = [temp_list[3],temp_list[0],temp_list[1],temp_list[2]]

#open zipped input FASTQs and assign them to variable names. Not using with open because I want to keep these files open until I've looped through all of them.
r1 = gzip.open(args.read1, 'rt')
r2 = gzip.open(args.read2, 'rt')
i1 = gzip.open(args.index1, 'rt')
i2 = gzip.open(args.index2, 'rt')

##use these for small unit tests where input files are not zipped
# r1 = open(args.read1, 'rt')
# r2 = open(args.read2, 'rt')
# i1 = open(args.index1, 'rt')
# i2 = open(args.index2, 'rt')

files_dict = {}
#keys are either barcode sequences or bucket name (swapped, unknown) and values are lists containing opened R1 and R2 file handles

# loop through barcode dict and open r1 and r2 files for each barcode for writing
for barcode_seq in barcode_dict:
    files_dict[barcode_seq] = [open(args.output_dir+f'{barcode_dict[barcode_seq][1]}_{barcode_dict[barcode_seq][2]}_{barcode_dict[barcode_seq][3]}_R1.fq', 'w'),open(args.output_dir+f'{barcode_dict[barcode_seq][1]}_{barcode_dict[barcode_seq][2]}_{barcode_dict[barcode_seq][3]}_R2.fq', 'w')]

# open files for writing for swapped and unknown indexes
files_dict['swapped'] = [open(args.output_dir+'swapped.R1.fq', 'w'), open(args.output_dir+'swapped.R2.fq', 'w')]
files_dict['unknown'] = [open(args.output_dir+'unknown.R1.fq', 'w'), open(args.output_dir+'unknown.R2.fq', 'w')]

def rev_comp(seq: str) -> str:
    '''Takes a sequence in the 5'->3' direction and returns its reverse complement in the 5'->3' direction'''
    temp_list = []
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N'}
    for base in seq:
        temp_list.insert(0, complement_dict[base])
    return ''.join(temp_list)

def check_qual(scorestring: str) -> bool:
    '''Takes a scorestring and returns True if every converted Phred+33 score is greater than or equal to 30, otherwise returns False.'''
    passes = True
    for score in scorestring:
        if Bioinfo.convert_phred(score) < 30:
            passes = False
            break
    return passes

def write_to_fastq(bucket: str):
    '''Takes a bucket/category name ('unknown,' 'swapped', or barcode sequence) and outputs one record to the R1.fastq file for that category and one record to the R2.fastq file for that category. Returns None.'''
    files_dict[bucket][0].write(header+'\n')
    files_dict[bucket][0].write(record_dict[r1][1]+'\n')
    files_dict[bucket][0].write(record_dict[r1][2]+'\n')
    files_dict[bucket][0].write(record_dict[r1][3]+'\n')
    files_dict[bucket][1].write(header+'\n')
    files_dict[bucket][1].write(record_dict[r2][1]+'\n')
    files_dict[bucket][1].write(record_dict[r2][2]+'\n')
    files_dict[bucket][1].write(record_dict[r2][3]+'\n')
    
swapped_counter = 0 #initialize counters for total swapped, unknown, and matched read pairs
unknown_counter = 0
matched_counter = 0
swapped_dict = {} #counts each combination of swapped indexes. key is barcode1-barcode2 and value is number of times that combo was seen. Note that ATG-CCC is the same as CCC-ATG, so only one of these combinations will appear in this dict. Therefore, must try both combinations when incrementing counter.
matched_dict = {} #counts each combination of matched indexes key is barcode1-barcode2 (where barcode1 = barcode2) and value is number of times that combo was seen.
outer_count = 0
for barcode1 in barcode_dict:
    matched_dict[barcode1] = 0 #initialize matched dict
    outer_count += 1
    inner_count = 0
    for barcode2 in barcode_dict: #initialize swapped dict
        inner_count += 1
        if inner_count > outer_count: #this makes it so that we are not counting the same combination twice just because it appears in a different order
            swapped_dict[barcode1+'-'+barcode2] = 0

while True:
    #initialize dictionary to store one record per file at a time; this gets initialized and repopulated after dealing with one record per file
    record_dict = {r1:[0,0,0,0],r2:[0,0,0,0],i1:[0,0,0,0],i2:[0,0,0,0]} #keys are file handles, values are lists where each item is one line
    #populate dictionary. Each value is a list where the first item is the first line of the record, the second item is the secod line of the record, etc.
    for key in record_dict:
        record_dict[key][0] = key.readline().strip()
        record_dict[key][1] = key.readline().strip()
        record_dict[key][2] = key.readline().strip()
        record_dict[key][3] = key.readline().strip()
    if record_dict[r1][0] == '': #ends loop when the end of the file is reached
        break
    record_dict[i2][1] = rev_comp(record_dict[i2][1]) #change i2 to its rev comp for easier comparison and ease of modifying the header
    header = re.sub(r'(.*)\s.*', r'\1', record_dict[r1][0])+'-'+record_dict[i1][1]+'-'+record_dict[i2][1] #modify header to header-i1-i2 (i2 is already rev comp'd). Take out anything after the space because that will be unique to the read, and I want both reads in the same pair to have the same header for now.
    #check quality scores, Ns, and whether indexes are in barcodes dict
    if record_dict[i1][1] in barcode_dict and record_dict[i2][1] in barcode_dict and check_qual(record_dict[i1][3]) and check_qual(record_dict[i2][3]):
        #check that index 1 matches rev comp index 2, assign to correct matched bucket based on barcode sequence if True
        if record_dict[i1][1] == record_dict[i2][1]:
            bucket = record_dict[i1][1]
            matched_counter += 1
            matched_dict[record_dict[i1][1]] += 1
        else: #if index 1 doesn't match index 2, assign to swapped bucket
            bucket = 'swapped'
            swapped_counter += 1
            if record_dict[i1][1]+'-'+record_dict[i2][1] in swapped_dict: #if statement is necessary here because the keys in the swapped counter dictionary have the barcodes in a certain order, but i1-i2 and i2-i1 are equivalent in this case because I just want to count how many times the (unordered) combination occurs
                swapped_dict[record_dict[i1][1]+'-'+record_dict[i2][1]] += 1
            else:
                swapped_dict[record_dict[i2][1]+'-'+record_dict[i1][1]] += 1
    else: #if quality conditions aren't met, assign to unknown bucket
        bucket = 'unknown'
        unknown_counter += 1
    write_to_fastq(bucket) #write the records for read 1 and read 2 to the appropriate file

for file in files_dict:
    files_dict[file][0].close()
    files_dict[file][1].close()

r1.close()
r2.close()
i1.close()
i2.close()

total_read_pairs = matched_counter + unknown_counter + swapped_counter

with open(args.output_dir+'stats.txt', 'w') as stats:
    stats.write('Number (and percentage) of read pairs with unknown or low-quality indexes: '+str(unknown_counter)+' ('+str(100*unknown_counter/total_read_pairs)+'%)\n')
    stats.write('Number (and percentage) of read pairs with index-swapping: '+str(swapped_counter)+' ('+str(swapped_counter*100/total_read_pairs)+'%)\n')
    stats.write('Number (and percentage) of read pairs with correctly-matched indexes: '+str(matched_counter)+' ('+str(matched_counter*100/total_read_pairs)+'%)\n')
with open(args.output_dir+'stats.tsv', 'w') as stats2:
    stats2_tsv = csv.writer(stats2, delimiter = '\t')
    stats2_tsv.writerow(['#Index 1', 'Index 2', 'Number of read pairs', 'Percentage out of all read pairs'])
    for key in matched_dict:
        stats2_tsv.writerow([key,key,str(matched_dict[key]),str(100*matched_dict[key]/total_read_pairs)])
    for key in swapped_dict:
        barcode1 = re.sub(r'(.+)-.*',r'\1',key)
        barcode2 = re.sub(r'.*-(.*)',r'\1',key)
        stats2_tsv.writerow([barcode1,barcode2,str(swapped_dict[key]),str(100*swapped_dict[key]/total_read_pairs)])
