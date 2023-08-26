# CRISPRi-dual-sgRNA-screens
## Fastq to counts
This repository contains scripts for alignment of sequencing data from dual-sgRNA CRISPR screening data. These scripts were adapted from: https://github.com/mhorlbeck/ScreenProcessing as well as https://github.com/abearab/CRISPRi-dual-sgRNA-screens .

For libraries with IBC, demultiplexing on only the i5 index using the i7 index (IBC) as a read is performed as detailed in : https://gist.github.com/sumeetg23/a064a36801d2763e94da2e191699fb9f. For alignment of data with IBC, use the python script `dualguide_UMI_fastqgz_to_counts.py`. Example command: 

```bash
python module1/dualguide_UMI_fastqgz_to_counts.py \
  library.csv \
  IBC_sequences.csv output_directory \
  SampleName1* SampleName2*
```

For alignment of data without IBC, use the python script `dualguide_fastqgz_to_counts.py`. Example command:

```bash
python module1/dualguide_fastqgz_to_counts.py \
  library.csv  \
  output_directory \
  SampleName1* SampleName2*
```

In these examples, there are a few inputs to the command:

1) library.csv is the library file containing the following columns: 'sgID_AB','sgID_A','protospacer_A','sgID_B','protospacer_B'. The 'protospacer_A' and 'protospacer_B' columns contain the 20bp nucleotide sequence of the sgRNA in position A and the sgRNA in position B in the library, respectively.

2) IBC_sequences.csv is the IBC sequence file containing only a single column named 'UMI' containing the 8bp nucleotide sequence of the IBC.

3) "output_directory" is the name of the output directory where you will find the output count files.

4) SampleName1* and SampleName2* are fastq files following standard naming convention (eg, SampleName1_S1_L001_R1_001.fastq.gz). Unix wildcards can be used to select multiple files at once. The script will search for all *.fastq.gz, *.fastq, and *.fa(/fasta/fna) files with the given wildcard name and detect the read 1 and read 2 files. 

The script will output two files of counts for each sample: *.all.aligned.counts.txt and *.AB.match.counts.txt. 

*.all.aligned.counts.txt contains the number of counts for every sgRNA A - sgRNA B - IBC combination aligned in the reads (including both pairs denoted in the library table and recombined pairs). The first column denotes the name of the name of the sgRNA A - sgRNA B - IBC combination where the names are separated by "++" delimiters. The second column denotes the number of reads aligning to this element.

*.AB.match.counts.txt contains the number of counts for sgRNA A - sgRNA B - IBC combinations aligned including only pairs denoted in the library table (excluding recombined pairs). The first column denotes the name of the name of the sgRNA A - sgRNA B - IBC combination where the names are separated by "++" delimiters. The second column denotes the number of reads aligning to this element.

These scripts will also output statistics detailing the (i) the number of reads that mapped to an sgRNA/IBC in the library by position and (ii) the number of reads with mapped sgRNAs that do not match the library and thus represent recombined reads to standard output. 

## Merge counts for multiple samples  
This part will merge outputs into a table in which rows are each sgRNA-UMI combination and columns are sample names. 

```bash
python module1/count_files_to_counts_matrix.py alignments results
```
### Output 
Raw counts matrix (`counts.txt`) 
```
columns: sample IDs
rows: sgID_AB counts
```

## Calculate the enrichment of each guide between screen arms
Sample sheet example: 
```
Index,Treat,Rep
i05,T0,2
i06,T0,1
i07,DMSO,1
i08,DMSO,2
i09,DMSO,3
i10,RM6,1
i11,RM6,2
i12,RM6,3 
```
From this, we can then calculate the enrichment of each guide between screen arms in using module 2. which employs DESeq2 and the Wald test between treatment and control input.

Command line script example: 
```bash
Rscript module2/PhenotypeScores.R \
  -i results/counts.txt \
  -s results/samplesheet.txt \
  -t RM6 -c DMSO
```

### Outputs
These outputs will be written in the same folder as the input count matrix.
1. Normalized counts (`counts_DESeq2_normalized.txt`): 
```
columns: sample IDs
rows: sgID_AB counts (DESeq2 normalized)
```
2. Phenotype scores (`{phenotype}{cond1}_vs_{cond2}.txt`)

```
    target  baseMean        log2FoldChange  lfcSE   stat    pvalue  padj
MAP1LC3B_-_87425948.23-P1P2|MAP1LC3B_-_87425925.23-P1P2++MAP1LC3B_-_87425948.23-P1P2|MAP1LC3B_-_87425925.23-P1P2        MAP1LC3B       598.454899564438        0.0345856896391499      0.14241765524649        0.502626099621637       0.777778849921918     0.88940640596369
```







