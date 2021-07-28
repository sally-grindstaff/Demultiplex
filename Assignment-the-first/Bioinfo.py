#!/usr/bin/env python

def validate_base_seq(seq: str,RNAflag: bool =False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNAbases = set('ATGCatcg')
    RNAbases = set('AUGCaucg')
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

if __name__ == '__main__':
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score for quality scores encoded with Phred+33"""
    dec: int = ord(letter)
    phred: int = dec - 33
    return phred

if __name__ == '__main__':
    assert convert_phred('I') == 40, 'Wrong phred score'
    print('Converted character into phred score correctly')

def gc_content(DNA: str) -> int:
    '''Takes DNA (or RNA) sequence and returns GC content of the sequence in decimal format as a fraction of 1.'''
    GC_count = 0
    if validate_base_seq(DNA):
        DNA = DNA.upper()
        for letter in DNA:
            if letter == 'G' or letter == 'C':
                GC_count += 1
    else:
        print('Not a valid sequence')
    return GC_count/len(DNA)

if __name__ == '__main__':
    assert gc_content('ATGCGA') == 0.5, 'gc_content does not find correct GC content'
    print('Calculated gc content correctly')

def qual_score(phred_score: str) -> int:
    '''Calculates average quality score of a read, given a quality score string in Phred+33 notation'''
    sum: int = 0
    for letter in phred_score:
        sum += convert_phred(letter)
    return sum/len(phred_score)

if __name__ == '__main__':
    qual_score('JJII') == 40.5, 'Incorrect average phred quality score'
    print('Correct average phred quality score')

def fasta_fixer(fasta: str):
    '''Takes FASTA filename as input and outputs a FASTA file that has exactly 2 lines per record (one header line and one sequence line).'''
    output = fasta[0:-2] + 'fixed.fa'
    seq = ''
    with open(fasta, 'r') as fh, open(output, 'x') as out:
        for line in fh:
            line = line.strip()
            if line[0] == '>' and seq != '':
                out.write(header+'\n')
                out.write(seq+'\n')
                seq = ''
            if line[0] == '>':
                header = line
            else:
                seq += line
        out.write(header+'\n')
        out.write(seq+'\n')


DNA_bases = 'ATCGNatcgn'
RNA_bases = 'AUCGNaucgn'
    