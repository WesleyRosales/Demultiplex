# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz |Read 1 |
| 1294_S1_L008_R2_001.fastq.gz |Index 1 |
| 1294_S1_L008_R3_001.fastq.gz |Index 2 |
| 1294_S1_L008_R4_001.fastq.gz |Read 2  |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.<br>
    ![Index 1 histogram](./index1_hist.png)<br>
    ![Index 2 histogram](./index2_hist.png)<br>
    ![Read 1 histogram](./read1_hist.png)<br>
    ![Read 2 histogram](./read2_hist.png)<br>
    ## Note: If images do not show, they are in the Assignment_the_first folder
    2. 20 is a good quality score cutoff because a quality score of 20 corresponds to a 99% accuracy. It is important to maintain high accuracy for both the indexes and the reads because a miscalled base can have dramatic effects on the analysis downstream.
    3. 1294_S1_L008_R2_001.fastq.gz has 3976613 indexes with "N" </br>
       command: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n "2~4p" | grep "N" | wc -l<br>
       1294_S1_L008_R3_001.fastq.gz has 3328051 indexes with "N"<br>
       command: zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n "2~4p" | grep "N" | wc -l
    
## Part 2
1. Define the problem<br>
Sequencing of pooled DNA needs to be demultiplexed into groups discerned by dual-matches barcodes, index hops, and poor quality indices
2. Describe output<br>
2 fastq files per barcode (1 for read 1, 1 for read 2)<br>
2 fastq files for index hopped reads (1 for read 1, 1 for read 2)<br>
2 fastq files for undetermined or poor quality barcodes (1 for read 1, 1 for read 2)<br>
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
#
    1. open barcodes.txt file and save barcodes in a dictionary with barcode name as key and barcode sequence as value
    2. Open each file as read
    3. open all necessary files to write (24 barcode files, 2 index swapped files, 2 poor quality files) as append to continuously add to them
    NOTE: fastq are matched to each other so record 827 of read1 is related to record 827 of read 2, index 1, and index 2 fastq files.
    4. IN A WHILE LOOP:
        for each record, save each line (except "+") as a variables.

        reverse complement index 2 sequences

        check for "" in any variables to indicate end of file

        if index 1 == rev comp index 2 AND index 1 is barcodes dictionary:
            make a fastq record for read 1 and read 2 with header having "index1-index1" appended to it
            add that record to file to be written in proper barcode file

        if index 1 != rev comp index 2 AND index 1 and rev comp index 2 are in barcodes dictionary:
            make a fastq record for read 1 and read 2 with header having "index1-rev comp index2" appended to it
            add that record to appropriate swapped file to be written

        if index 1 or index 2 have qscore < 20 OR they have "N" in sequence:
            make a fastq record for read 1 and read 2 with header having "index1-rev comp index2" appended to it
            add that record to appropriate poor quality file to be written
    5. close all files

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

Function 1: reverse complement<br>
Doc string: reverse complements an input DNA string; any "N" values stay as "N"<br>
function header: rev_comp(DNA)<br>
tests: "ATCG" -> "CGAT"; "ATNN" -> "NNTA"<br>
return revcomp_string

