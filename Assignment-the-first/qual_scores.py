#!/usr/bin/env python

import Bioinfo
from matplotlib import pyplot as plt
import gzip

read1_file = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
read2_file = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz'
index1_file = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz'
index2_file = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz'

# test_file = '/projects/bgmp/sgrindst/bioinformatics/Bi622/Demultiplex/TEST-input_FASTQ/test.R1.fq'
# test2_file = '/projects/bgmp/sgrindst/bioinformatics/Bi622/Demultiplex/TEST-input_FASTQ/test.R2.fq'

def init_list(lst: list, read_length: int, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    for i in range(read_length):
        lst.append(value)
    return lst

def populate_list(file, read_length):
    """Takes a FASTQ file and returns a list in which each item is the total sum of Q scores across all reads for a base position. The index of the base position in the read matches the index of the sum in the list. Assumes all reads are the same length."""
    qscore_list = []
    qscore_list = init_list(qscore_list, read_length)
    with gzip.open(file, 'rt') as fh:
        num_lines = 0
        for line in fh:
            line = line.strip()
            num_lines += 1
            if num_lines % 4 == 0:
                i = 0
                for score in line:
                    qscore_list[i] += Bioinfo.convert_phred(score)
                    i += 1
        i = 0
        for score in qscore_list:
            qscore_list[i] = score/(num_lines/4)
            i +=1
    return qscore_list

read1_list = populate_list(read1_file, 101)
read2_list = populate_list(read2_file, 101)
index1_list = populate_list(index1_file, 8)
index2_list = populate_list(index2_file, 8)
# test_list = populate_list(test_file,9)
# test2_list = populate_list(test2_file,4)

def make_hist(some_list:list, name:str):
    plt.bar(range(len(some_list)), some_list)
    plt.xlabel('# Base Pair')
    plt.ylabel('Mean quality score')
    plt.title(f'Quality Score Distribution in {name}')
    plt.savefig(f'{name}_qscore_distribution.png')
    plt.close()

make_hist(read1_list,'1294_S1_L008_R1_001')
make_hist(read2_list,'1294_S1_L008_R4_001')
make_hist(index1_list,'1294_S1_L008_R2_001')
make_hist(index2_list,'1294_S1_L008_R3_001')
# make_hist(test_list, 'test')
# make_hist(test2_list, 'test2')
