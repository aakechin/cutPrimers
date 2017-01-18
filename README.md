# cutPrimers
curPrimers is a tool for trimming primer sequences from amplicon based NGS reads

## Requirements
cutPrimers works on the Python3+ and requires the following packages:
* **Biopython** - you can install it with: `sudo apt-get install python3-biopython`
* **regex** - you can install it with: `sudo apt-get install python3-regex`
* **argparse** - you can install it with: `sudo apt-get install python3-argparse`

## Installation
cutPrimers does not require any installations

## Use
You can see all parameters with 
```
python3 cutPrimers.py -h
```

## Example of use
As an example you can use files from directory examples. First untar them and try to trim them with the following commands:
```
tar -xf example/*.tar.gz
python3 cutPrimers.py -r1 example/1_S1_L001_R1_001.fastq -r2 example/1_S1_L001_R2_001.fastq -pr15 example/primers_R1_5.fa -pr25 example/primers_R2_5.fa -pr13 example/primers_R1_3.fa -pr23 example/primers_R2_3.fa -tr1 example/1_r1_trimmed.fastq -tr2 example/1_r2_trimmed.fastq -utr1 example/1_r1_untrimmed.fastq -utr2 example/1_r2_untrimmed.fastq -t 2
```

## Parameters
```
-h, --help - show this help message and exit
  --readsFile_r1, -r1 - file with R1 reads of one sample
  --readsFile_r2, -r2 - file with R2 reads of one sample
  --primersFileR1_5, -pr15 - fasta-file with sequences of primers on the 5'-end of R1 reads
  --primersFileR2_5, -pr25 - fasta-file with sequences of primers on the 5'-end of R2 reads. Do not use this parameter if you have single-end reads
  --primersFileR1_3, -pr13 - fasta-file with sequences of primers on the 3'-end of R1 reads. It is not required. But if it is determined, -pr23 is necessary
  --primersFileR2_3, -pr23 - fasta-file with sequences of primers on the 3'-end of R2 reads
  --trimmedReadsR1, -tr1 - name of file for trimmed R1 reads
  --trimmedReadsR2, -tr2 - name of file for trimmed R2 reads
  --untrimmedReadsR1, -utr1 - name of file for untrimmed R1 reads
  --untrimmedReadsR2, -utr2 - name of file for untrimmed R2 reads
  --primersStatistics, -stat - name of file for statistics of errors in primers. This works only for paired-end reads with primers at 3'- and 5'-ends
  --error-number, -err - number of errors (substitutions, insertions, deletions) that allowed during searching primer sequence in a read sequence. Default: 5
  --primer-location-buffer, -plb - Buffer of primer location in the read from the end. Default: 10
  --threads, -t - number of threads
```
## Citation
Manuscript is prepared. Now you can cite it by link.
