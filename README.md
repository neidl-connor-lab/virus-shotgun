# virus-shotgun
v2.0 by Jackie T.

## Requirements
Computing cluster access to [bowtie2](https://doi.org/10.1038%2Fnmeth.1923), [samtools](http://www.htslib.org/download/), and [R](https://www.r-project.org/) modules are required. Local installations have not been tested.

| Package   | Version |
| :-------- | :------ |
| bowtie2   | 2.4.2   |
| htslib    | 1.18    |
| samtools  | 1.18    |
| R         | 4.0.2   |

You will also need to install [kraken2](https://ccb.jhu.edu/software/kraken2/) either locally or globally.

## Quick start

### 1. Download repo 

Cloning this repo will create a directory holding the pipeline script and assorted files. Move to the directory where you'd like to save the pipeline; here, we're using `workflows/`.

```
cd workflows/
git clone https://github.com/neidl-connor-lab/virus-shotgun.git
cd virus-shotgun
```

### 2. Run `setup.sh` to set up [LoFreq](https://doi.org/10.1093%2Fnar%2Fgks918), [Kraken2](https://doi.org/10.1186/s13059-019-1891-0), and make an index

Create an alignment index from a reference sequence FASTA file, and write it to `pipeline/indices/`. Most viruses have an official NCBI [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/). This script will also set up LoFreq and Kraken2 in your `pipeline/` directory. 

Run `setup.sh` with the `-h` flag to view the full list of options.

```
./setup.sh -h
```

| Flag | Argument                                  |
| :--- | :---------------------------------------- |
| `-P` | SCC project                               |
| `-N` | job name                                  |
| `-f` | genome reference FASTA file               |
| `-b` | bowtie2 index ID                          |

```
usage: qsub -P PROJECT -N JOBNAME ./setup.sh -f FASTA -b BOWTIE

arguments:
  -f virus genome FASTA file
  -b bowtie2 index path and prefix
  -h show this message and exit
```

Here is an example where we download the [SARS-CoV-2 RefSeq](https://www.ncbi.nlm.nih.gov/nuccore/1798174254) and then run `setup.sh`. If this is your first time running `setup.sh` in the directory, it will unpack LoFreq and set up the [standard Kraken2 database](https://benlangmead.github.io/aws-indexes/k2) as well!

```
# download and decompress the reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz
gunzip GCA_009858895.3_ASM985889v3_genomic.fna.gz

# move the reference into your pipeline directory for safekeeping
mv GCA_009858895.3_ASM985889v3_genomic.fna pipeline/

# submit the indexing job
qsub -P project \
     -N test-index \
     setup.sh \
     -f pipeline/GCA_009858895.3_ASM985889v3_genomic.fna \
     -b sarscov2

# wait until the job is done
qstat -u $(whoami)

# check out the pipeline directory!
ls pipeline/*
```

If you would like to make another index, just run `setup.sh` again! It will be much faster the second time since it will skip the Lofreq and Kraken2 setup.

### 3. Run pipeline

> You _must_ run this script from cloned project directory used in step 2.

View pipeline options and required arguments by running `pipeline.sh` with the `-h` flag.

```
$ pipeline.sh -h
usage: qsub -P PROJECT -N JOBNAME ./pipeline.sh -i INDEX -f FASTA -o ODIR -s SAMPLE -x R1 [-y R2] [-t THLD]
Please submit the job from the pipeline directory!

arguments:
  -i bowtie2 index path and prefix
  -f reference FASTA
  -o output directory
  -s sample ID
  -x FASTQ file; R1 file if paired reads
  -y [OPTIONAL] R2 FASTQ file if paired reads
  -t [OPTIONAL] minimum aligned read depth (default: 10)
  -h print this message and exit
```

The help message indicates the required arguments and how to pass them:

| Flag | Argument                                  |
| :--- | :---------------------------------------- |
| `-P` | SCC project                               |
| `-N` | job name                                  |
| `-i` | path to index created in step 2           |
| `-f` | reference FASTA file                      |
| `-o` | output directory                          |
| `-x` | path to R1 or unpaired FASTQ file         |
| `-y` | path to R2 FASTQ file (paired-read only)  |
| `-t` | minimum aligned read depth threshold      |

Here is an example where we're using the index we created in step 2. The job output will be written to a file named `log-test.qlog`. Fill in your own project allocation and FASTQ files!

```
qsub -P project \
     -N test-pipeline \
     pipeline.sh \
     -i pipeline/indices/sarscov2 \
     -f pipeline/GCA_009858895.3_ASM985889v3_genomic.fna \
     -o data \
     -s test \
     -x test-r1.fq.gz \
     -y test-r2.fq.gz \
     -t 50
```

## Pipeline steps

## 1. Check for contaminants

Raw reads are passed to Kraken2 for metagenomic classification using the standard Kraken database constructed when running `setup.sh` for the first time.

| Flag                | Meaning                |
| :------------------ | :--------------------- |
| `--threads`         | parallelize this job   |
| `--db`              | path to database       |
| `--output`          | raw output filename    |
| `--report`          | report filename        |
| `--use-names`       | use taxon names        |
| `--paired`          | paired reads           |
| `*.fq.gz`           | input file(s)          |

```
# paired
kraken2 --threads 4 --db pipeline/kraken2db --output - --report odir/sample/metagenomics.tsv --use-names --paired paired-r1.fq.gz paired-r2.fq.gz

# unpaired
kraken2 --threads 4 --db pipeline/kraken2db --output - --report odir/sample/metagenomics.tsv --use-names unpaired.fq.gz
```

### 2. Align to reference

Raw FASTQ files are aligned to the previously-constructed reference using [Bowtie2](https://doi.org/10.1038%2Fnmeth.1923). All output files are placed in a subdirectory named with the sample ID. 

| Flag        | Meaning                |
| :---------- | :--------------------- |
| `--threads` | parallelize this job   |
| `-x`        | index path and prefix  |
| `-1`        | (paired) R1 FASTQ file |
| `-2`        | (paired) R2 FASTQ file |
| `-U`        | (unpaired) FASTQ file  |
| `*.sam`     | uncompressed alignment |
| `*.log`     | bowtie2 output stats   |

```
# paired
bowtie2 --threads 4 -x 'pipeline/bowtie/index' -1 'paired-r1.fq.gz' -2 'paired-r2.fq.gz' 1> 'odir/sample/alignment.sam' 2> 'odir/sample/bowtie2.log'

# unpaired
bowtie2 --threads 4 -x 'pipeline/bowtie/index' -U 'unpaired.fq.gz' 1> 'odir/sample/alignment.sam' 2> 'odir/sample/bowtie2.log'
```

The uncompressed SAM output is then compressed to BAM format.

```
# compress
samtools view --threads 4 -b -h 'odir/sample/alignment.sam' > 'odir/sample/alignment-raw.bam'
```

### 3. Process alignment

These steps prep the alignment file for coverage, SNV, and consensus calling. We use samtools to sort the BAM, LoFreq to score insertions and deletions, and then samtools again to index the alignment.

| Flag            | Meaning                      |
| :-------------- | :--------------------------- |
| `--threads`     | parallelize this job         |
| `*-clipped.bam` | alignment from step 2        |
| `*-sorted.bam`  | sorted alignment             |
| `--dindel`      | algorithm for scoring indels |
| `--ref`         | reference FASTA file         |
| `alignment.bam` | final alignment file         |

```
samtools sort --threads 4 'odir/sample/alignment-clipped.bam' > 'odir/sample/alignment-sorted.bam'

lofreq indelqual --dindel --ref 'reference.fa' 'odir/sample/alignment-sorted.bam' > 'odir/sample/alignment.bam'

samtools index 'odir/sample/id.bam'
```

### 4. Calculate coverage

The `samtools depth` command calculates the aligned read depth for each nucleotide of the genome used for alignment. The output is an easy-to-analyze TSV table.

| Flag            | Meaning                 |
| :-------------- | :---------------------- |
| `--threads`     | parallelize this job    |
| `-a`            | include all nucleotides |
| `-H`            | include a file header   |
| `alignment.bam` | final alignment file    |
| `coverage.tsv`  | coverage table          |

```
samtools depth --threads 4 -a -H 'odir/sample/alignment.bam' > 'odir/sample/coverage.tsv'
```

### 5. Assemble consensus

The `samtools consensus` command assembles a consensus by examining the reads aligned to each nucleotide in the reference sequence and calling the most frequent allele.

| Flag            | Meaning                      |
| :-------------- | :--------------------------- |
| `--threads`     | parallelize this job         |
| `--use-qual`    | use quality scores           |
| `--min-depth`   | minimum aligned read depth   |
| `--call-fract`  | minimum SNV frequency        |
| `--output`      | output FASTA file            |
| `alignment.bam` | final aligment file          |

```
samtools consensus --threads 4 --use-qual --min-depth 10 --call-fract 0.5 --output 'odir/sample/consensus.fa' 'odir/sample/alignment.bam'
```

### 6. Quantify SNVs

Use LoFreq to make a detailed table of consensus and sub-consensus SNVs.

| Flag            | Meaning                    |
| :-------------- | :------------------------- |
| `--pp-threads`  | parallelize this job       |
| `--call-indels` | include indels in output   |
| `--min-cov`     | minimum aligned read depth |
| `--ref`         | reference FASTA file       |
| `alignment.bam` | final alignment file       |
| `snvs.vcf`      | output VCF                 |

```
lofreq call-parallel --pp-threads 4 --call-indels --min-cov 10 --ref 'reference.fa' 'odir/sample/alignment.bam' > 'odir/sample/snvs.vcf'
```

Use R to format the VCF file into a more human-readable CSV format.

| Flag      | Meaning    |
| :-------- | :--------- |
| `--vcf`   | LoFreq VCF |
| `--ofile` | CSV output |

```
Rscript pipeline/format.r --vcf 'odir/sample/snvs.vcf' --ofile 'odir/sample/snvs.csv'
```
