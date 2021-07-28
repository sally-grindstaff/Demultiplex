# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 |
| 1294_S1_L008_R2_001.fastq.gz | index1 |
| 1294_S1_L008_R3_001.fastq.gz | index2 |
| 1294_S1_L008_R4_001.fastq.gz | read2 |

2. Per-base NT distribution
    1. 
    (1294_S1_L008_R1_001_qscore_distribution.png)
    (1294_S1_L008_R2_001_qscore_distribution.png)
    (1294_S1_L008_R3_001_qscore_distribution.png)
    (1294_S1_L008_R4_001_qscore_distribution.png)
    2. ```
        The cutoff I chose for average quality score for reads from files R1 and R4 (read 1 and read 2) is Q20. This is because the probability that a base with a Q-score of 20 is correct is 99%, which is acceptable because the reads from this experiment (mRNA-seq) will be aligned to a reference in later steps, so we don't need to be too stringent with the cutoff; we would like to have as much data as possible while still aligning things correctly, and 1 base being wrong out of every 100 seems acceptable for the purpose of alignment.
        The cutoff I chose for index reads (files R2 and R3) is Q30 for each individual base. In other words, if any one quality score in the index read has a quality score lower than 30, the entire read pair is thrown out. This is because every single base matters to correctly map the read pair to the sample it came from, and if any one base is wrong, we may end up mapping the read pair to the wrong sample. Therefore, applying the quality score cutoff to each base rather than to the average across the read makes sense, because it allows us to have confidence that every single base in the index is correct. I chose a cutoff of Q30 for the index sequences, because this is the negative log of Illumina's average error rate (0.1%).
       ```
    3. ```
       zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep 'N' | wc -l
       zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep 'N' | wc -l
       3976613 (R2) + 3328051 (R3)= 7,304,664 total indexes with Ns
       ```
    
## Part 2
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
