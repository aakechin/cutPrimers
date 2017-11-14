# cutPrimers
curPrimers is a tool for trimming primer sequences from amplicon based NGS reads

## Requirements
### Linux
cutPrimers works on the Python3+ and requires the following packages:
* **Biopython** - you can install it with: `sudo apt-get install python3-biopython` or download it from http://biopython.org/wiki/Download and install it locally with `python3 setup.py install --user`
* **regex** - you can install it with: `sudo apt-get install python3-regex`  or download it from https://pypi.python.org/pypi/regex/ and install it locally with `python3 setup.py install --user`
* **argparse** - you can install it with: `sudo pip3 install argparse` or download it from https://pypi.python.org/pypi/argparse and install it locally with `python3 setup.py install --user`

### Windows
For use on windows download and install python3.6 from www.python.org/downloads/ (**Attention! Remember to check "Add Python 3.6 to PATH" in the bottom of the installation window!**). After installation, restart your computer.

After that, install the followong packages with respective commands in command line (to run command line, search in Start menu "cmd" and run "cmd.exe"):
* **Biopython** - with: `pip install biopython`. If you do not have Visual Studio C++ already installed, pip will show an error. In that case, download and install it from landinghub.visualstudio.com/visual-cpp-build-tools/
* **regex** - you can install it with: `pip install regex`
* **argparse** - you can install it with: `pip install argparse`

### Mac OS
For use on Mac OS download and install python3.6 from www.python.org/downloads/

After that install the followong packages with respective commands in command line:
* **Biopython** - with: `pip install biopython`
* **regex** - you can install it with: `pip install regex`
* **argparse** - you can install it with: `pip install argparse`

## Installation
cutPrimers does not require any installations

## Use
You can see all parameters with 
```
python3 cutPrimers.py -h
```

## Example of use
As an example you can use files from directory "examples". Trim them with the following commands:
```
python3 cutPrimers.py -r1 example/1_S1_L001_R1_001.fastq.gz -r2 example/1_S1_L001_R2_001.fastq.gz -pr15 example/primers_R1_5.fa -pr25 example/primers_R2_5.fa -pr13 example/primers_R1_3.fa -pr23 example/primers_R2_3.fa -tr1 example/1_r1_trimmed.fastq.gz -tr2 example/1_r2_trimmed.fastq.gz -utr1 example/1_r1_untrimmed.fastq.gz -utr2 example/1_r2_untrimmed.fastq.gz -t 2
```
As a result you should get four files with the following sizes: 4.3 Mb, 2.3 Mb, 3.7 Mb and 2.3 Mb for files with trimmed R1 reads, untrimmed R1 reads, trimmed R2 reads and untrimmed R2 reads, respectively. All files will be gzipped. cutPrimers can remove primer sequences both from gzipped and native FASTQ-files.

Before using cutPrimers we recommend to use Trimmomatic for removing adaptor sequences from 3'-ends of reads. In this case, you will get more reliable results. You can use cutPrimers with Trimmomatic with the following commands:
```
mkdir example_trimmed
java -jar trimmomatic-0.32.jar PE example/${i}_S${i}_L001_R1_001.fastq.gz example/${i}_S${i}_L001_R2_001.fastq.gz example_trimmed/patient_${i}.r1.ad_trimmed.fastq.gz example_trimmed/patient_${i}.r1.ad_unpaired.fastq.gz example_trimmed/patient_${i}.r2.ad_trimmed.fastq.gz example_trimmed/patient_${i}.r2.ad_unpaired.fastq.gz ILLUMINACLIP:examples/adapters.fa:3:10:7
python3 cutPrimers.py -r1 example_trimmed/patient_${i}.r1.ad_trimmed.fastq.gz -r2 example_trimmed/patient_${i}.r2.ad_trimmed.fastq.gz -pr15 example/primers_R1_5.fa -pr13 example/primers_R1_3.fa -pr25 example/primers_R2_5.fa -pr23 example/primers_R2_3.fa -tr1 example_trimmed/patient_${i}.r1.ad_trimmed.trimmed.fastq.gz -tr2 example_trimmed/patient_${i}.r2.ad_trimmed.trimmed.fastq.gz -utr1 example_trimmed/patient_${i}.r1.ad_trimmed.untrimmed.fastq.gz -utr2 example_trimmed/patient_${i}.r2.ad_trimmed.untrimmed.fastq.gz -t 12 -primer3
java -jar trimmomatic-0.32.jar PE example_trimmed/patient_${i}.r1.ad_trimmed.trimmed.fastq.gz example_trimmed/patient_${i}.r2.ad_trimmed.trimmed.fastq.gz example_trimmed/patient_${i}.r1.ad_trimmed.trimmed.qual_trimmed.fastq.gz example_trimmed/patient_${i}.r1.ad_trimmed.trimmed.qual_unpaired.fastq.gz example_trimmed/patient_${i}.r2.ad_trimmed.trimmed.qual_trimmed.fastq.gz example_trimmed/patient_${i}.r2.ad_trimmed.trimmed.qual_unpaired.fastq.gz LEADING:10 TRAILING:10 MINLEN:10
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
  --min-primer3-length MINPRIMER3LEN, -primer3len MINPRIMER3LEN - minimal length of primer on the 3'-end to trim. Use this parameter, if you are ready to trim only part of primer sequence of the 3'-end of read
  --primer3-absent, -primer3 - if primer at the 3'-end may be absent, use this parameter
  --identify-dimers IDIMER, -idimer IDIMER - use this parameter if you want to get statistics of homo- and heterodimer formation. Choose file to which statistics of primer-dimers will be written. This parameter may slightly decrease the speed of analysis
  --threads, -t - number of threads
```
## Citation
**cutPrimers: A New Tool for Accurate Cutting of Primers from Reads of Targeted Next Generation Sequencing**. Kechin A, Boyarskikh U, Kel A, Filipenko M, 2017, Journal of Computational Biology, 2017 Jul 17. doi: 10.1089/cmb.2017.0096 (https://www.ncbi.nlm.nih.gov/pubmed/28715235)
