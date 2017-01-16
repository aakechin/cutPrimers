# cutPrimers
curPrimers is a tool for trimming primer sequences from amplicon based NGS reads

## Requirements
cutPrimers works on the Python3+ and requires the following packages:
* Biopython v1.65 or higher (older version has not been tested) - you can install with:
```sudo apt-get install python3-biopython
* regex - you can install with: ```sudo apt-get install python3-regex```
* argparse - you can install with: ```sudo apt-get install python3-argparse```

## Installation
cutPrimers does not require any installations

## Use
You can see all parameters with ```python3 cutPrimers.py -h```

## Parameters
-h, --help            show this help message and exit
  --readsFile_r1 READSFILE1, -r1 READSFILE1
                        file with R1 reads of one sample
  --readsFile_r2 READSFILE2, -r2 READSFILE2
                        file with R2 reads of one sample
  --primersFileR1_5 PRIMERSFILER1_5, -pr15 PRIMERSFILER1_5
                        fasta-file with sequences of primers on the 5'-end of
                        R1 reads
  --primersFileR2_5 PRIMERSFILER2_5, -pr25 PRIMERSFILER2_5
                        fasta-file with sequences of primers on the 5'-end of
                        R2 reads. Do not use this parameter if you have
                        single-end reads
  --primersFileR1_3 PRIMERSFILER1_3, -pr13 PRIMERSFILER1_3
                        fasta-file with sequences of primers on the 3'-end of
                        R1 reads. It is not required. But if it is determined,
                        -pr23 is necessary
  --primersFileR2_3 PRIMERSFILER2_3, -pr23 PRIMERSFILER2_3
                        fasta-file with sequences of primers on the 3'-end of
                        R2 reads
  --trimmedReadsR1 TRIMMEDREADSR1, -tr1 TRIMMEDREADSR1
                        name of file for trimmed R1 reads
  --trimmedReadsR2 TRIMMEDREADSR2, -tr2 TRIMMEDREADSR2
                        name of file for trimmed R2 reads
  --untrimmedReadsR1 UNTRIMMEDREADSR1, -utr1 UNTRIMMEDREADSR1
                        name of file for untrimmed R1 reads
  --untrimmedReadsR2 UNTRIMMEDREADSR2, -utr2 UNTRIMMEDREADSR2
                        name of file for untrimmed R2 reads
  --primersStatistics PRIMERSSTATISTICS, -stat PRIMERSSTATISTICS
                        name of file for statistics of errors in primers. This
                        works only for paired-end reads with primers at 3'-
                        and 5'-ends
  --error-number ERRNUMBER, -err ERRNUMBER
                        number of errors (substitutions, insertions,
                        deletions) that allowed during searching primer
                        sequence in a read sequence. Default: 5
  --primer-location-buffer PRIMERLOCBUF, -plb PRIMERLOCBUF
                        Buffer of primer location in the read from the end.
                        Default: 10
  --threads THREADS, -t THREADS
                        number of threads
