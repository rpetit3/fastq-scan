[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/fastq-scan/README.html)
[![Docker Repository on Quay.io](https://quay.io/repository/biocontainers/fastq-scan/status "Docker Repository on Quay.io")](https://quay.io/repository/biocontainers/fastq-scan)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/fastq-scan/badges/downloads.svg)](https://anaconda.org/bioconda/fastq-scan)

*fastq-scan* reads a FASTQ from STDIN and outputs summary statistics (read lengths, per-read qualities, per-base qualities) in JSON format.

# fastq-scan
I wanted a quick method to output simple summary statistics of an input FASTQ in JSON format. There are (likely better) alternatives including [FastQC](https://github.com/s-andrews/FastQC) and [seqtk's fqchk](https://github.com/lh3/seqtk), but they didn't ouput JSON which I wanted. After a brief Google search, I stumbled upon the Biostars question: ["How To Efficiently Parse A Huge Fastq File?"](https://www.biostars.org/p/10353/#10358). From this question I used code from Pierre Lindenbaum's C++ solution as the base for this program.

# Installation
### Bioconda
*fastq-scan* is availble on [BioConda](https://bioconda.github.io/recipes/fastq-scan/README.html).
```
conda install -c bioconda fastq-scan
```

### From Source
```
git clone git@github.com:rpetit3/fastq-scan.git
cd fastq-scan
make
```
This will compile the program using g++, and I'll let you handle where to put it. I've testing this on gcc version 7.3.0, but would guess since its pretty basic most gcc versions should be fine. You can then run `make test` to make sure everything is working as expected.

```
make test
./fastq-scan -h
Usage: cat FASTQ | fastq-scan [options]
Version: 0.3

Optional arguments:
    -g INT   Genome size for calculating estimated sequencing coverage. (Default 1)
    -p INT   ASCII offset for input quality scores, can be 33 or 64. (Default 33)
    -v       Print version information and exit
    -h       Show this message and exit

cat example.fq | ./fastq-scan -g 150000
{
    "qc_stats": {
        "total_bp":7500,
        "coverage":0.05,
        "read_total":75,
        "read_min":100,
        "read_mean":100,
        "read_std":0,
        "read_median":100,
        "read_max":100,
        "read_25th":100,
        "read_75th":100,
        "qual_mean":34.0267,
        "qual_std":0.711306,
        "qual_median":34,
        "qual_25th":34,
        "qual_75th":34
    },
    "read_lengths": {

        "100":75
    },
    "per_base_quality": {
        "1":30.7467,        "2":31.5467,        "3":31.5467,        "4":35.44,        "5":34.24,
        "6":34.12,        "7":34.7067,        "8":34.24,        "9":36.9333,        "10":37.0667,
        "11":35.88,        "12":36.0667,        "13":36.72,        "14":38.2667,        "15":37.48,
        "16":38.2133,        "17":36.7467,        "18":37.8267,        "19":36.3333,        "20":37.2933,
        "21":37.9867,        "22":37.1067,        "23":37.4133,        "24":38.2667,        "25":36.6133,
        "26":36.2,        "27":36.3067,        "28":35.8533,        "29":36.5067,        "30":37.72,
        "31":37.3333,        "32":36.0133,        "33":37.4933,        "34":36.1067,        "35":36.76,
        "36":34.8533,        "37":36.3733,        "38":35.1867,        "39":36.0133,        "40":35.3067,
        "41":35.6,        "42":36.7867,        "43":35.52,        "44":37.3333,        "45":36.6533,
        "46":36.8,        "47":35.9867,        "48":35.4533,        "49":35.2,        "50":37.2533,
        "51":35.04,        "52":36,        "53":35.28,        "54":36.16,        "55":35.2,
        "56":33.6133,        "57":36.0533,        "58":34.4533,        "59":35.88,        "60":35.3733,
        "61":35.6933,        "62":34.8267,        "63":35.1067,        "64":35.2933,        "65":32.2667,
        "66":34.4267,        "67":33.9333,        "68":33.6667,        "69":32.6133,        "70":33.4267,
        "71":32.8267,        "72":32.96,        "73":33.5467,        "74":33.1067,        "75":31.8667,
        "76":30.72,        "77":30.6133,        "78":30.2133,        "79":31.7467,        "80":33.8933,
        "81":32.72,        "82":33.1733,        "83":31.5867,        "84":32.6933,        "85":32.0667,
        "86":32.2933,        "87":30.7467,        "88":30.6933,        "89":32.48,        "90":31.08,
        "91":31.6133,        "92":31.72,        "93":30.3867,        "94":30.7067,        "95":29.9733,
        "96":31.96,        "97":32.44,        "98":30.2267,        "99":31.2533,        "100":30.2267
    }
}

```

# Example Usage
*fastq-scan* reads from STDIN, so pretty much any FASTQ output can be piped into fastq-scan. There are a few things to be aware of. I've assumed that all FASTQ entries are the four line variant, which should be a safe assumption in 2018. Also, I have a PHRED offset (33 vs 64) guesser function. By default it is set to PHRED33, it could produce errors if there are not any PHRED33 or PHRED64 specific characters in the quality scores.

### Usage
```
./fastq-scan
Usage: cat FASTQ | fastq-scan [options]
Version: 0.3

Optional arguments:
    -g INT   Genome size for calculating estimated sequencing coverage. (Default 1)
    -p INT   ASCII offset for input quality scores, can be 33 or 64. (Default 33)
    -v       Print version information and exit
    -h       Show this message and exit
```

#### *-g* Genome Size
This is an optional parameter that you can use to estimate the sequencing coverage as calulated by (TOTAL_BP / GENOME_SIZE). By default, the genome size is set to 1. If a genome size is not given the total coverage will not be calculated and instead be set to `0.00` in the final JSON output.

#### *-p* ASCII Offset for Quality Scores
This optional parameter can be used to explicitely state the ASCII offset for the input quality scores. It defaults to PHRED+33 scores. Only `33` or `64` are valid inputs for this parameter.

#### *-v* Version
```
./fastq-scan -v
fastq-scan 0.3
```

### *example.fq*
An example FASTQ file, aptly named *example-q33.fq* (also *example-q64.fq*), has been included to demonstrate usage of *fastq-scan*. For those interested this is a small set of simulated reads from the [*Lotus japonicus* (NC_002694)](https://www.ncbi.nlm.nih.gov/nuccore/NC_002694.1) chloroplast genome. The reads were simulated using [ART (vMountRainier)](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm).
```
head example-q33.fq
@NC_002694.1-75
TGTATACAATAAGAATCCATTTATTGACAAATTTCATTCGAAAATTATGAAACATAAATTTTTTTTTATTGGATCAAGAATTCCAATTTTTTAAGTATAA
+
@CCDFFF*HDGHFGAIGJJJGIIHIJCHJJ*EJJHG@?GJFJGI-IBFGJ.JHGB.DDFDCEFA?)FC=JCDFF9DAA;E?>DDDDDF5?BCDDDDD><(
@NC_002694.1-74
TTTTTTGGACTTGAAGAAAAAAAAATACAACTTTGCTGACAATTATTTGTTTGGTCAGAAGAGTCCTCCAAATATTCTGATCTTATATTGATTATCATTT
+
CCCFFFFFHGF=>GJEJIJHJHBJ,)IDJGIJFIDFHBA@IGCIIBJGJIJI)JG0JJCJDIGHDIJI(9JDD?=<BHD<DCEDBBDAAD7DDDA;D?ED
@NC_002694.1-73
AAAAAGTGAAATATTCAGTTAATGAATGCCGAATCTCCGCTCTTATTCTATGAACATTTCATAATCCTATAAATTATCTTTATAAATTATCTTATAGAAT
```

### Execution
```
cat example-q33.fq | ./fastq-scan -g 150000
```


### Example Output
After the FASTQ file is read, simple summary statistics are output in JSON format. Below is the summary statistics of *example.fq*
```
{
    "qc_stats": {
        "total_bp":7500,
        "coverage":0.05,
        "read_total":75,
        "read_min":100,
        "read_mean":100,
        "read_std":0,
        "read_median":100,
        "read_max":100,
        "read_25th":100,
        "read_75th":100,
        "qual_mean":34.0267,
        "qual_std":0.711306,
        "qual_median":34,
        "qual_25th":34,
        "qual_75th":34
    },
    "read_lengths": {

        "100":75
    },
    "per_base_quality": {
        "1":30.7467,        "2":31.5467,        "3":31.5467,        "4":35.44,        "5":34.24,
        "6":34.12,        "7":34.7067,        "8":34.24,        "9":36.9333,        "10":37.0667,
        "11":35.88,        "12":36.0667,        "13":36.72,        "14":38.2667,        "15":37.48,
        "16":38.2133,        "17":36.7467,        "18":37.8267,        "19":36.3333,        "20":37.2933,
        "21":37.9867,        "22":37.1067,        "23":37.4133,        "24":38.2667,        "25":36.6133,
        "26":36.2,        "27":36.3067,        "28":35.8533,        "29":36.5067,        "30":37.72,
        "31":37.3333,        "32":36.0133,        "33":37.4933,        "34":36.1067,        "35":36.76,
        "36":34.8533,        "37":36.3733,        "38":35.1867,        "39":36.0133,        "40":35.3067,
        "41":35.6,        "42":36.7867,        "43":35.52,        "44":37.3333,        "45":36.6533,
        "46":36.8,        "47":35.9867,        "48":35.4533,        "49":35.2,        "50":37.2533,
        "51":35.04,        "52":36,        "53":35.28,        "54":36.16,        "55":35.2,
        "56":33.6133,        "57":36.0533,        "58":34.4533,        "59":35.88,        "60":35.3733,
        "61":35.6933,        "62":34.8267,        "63":35.1067,        "64":35.2933,        "65":32.2667,
        "66":34.4267,        "67":33.9333,        "68":33.6667,        "69":32.6133,        "70":33.4267,
        "71":32.8267,        "72":32.96,        "73":33.5467,        "74":33.1067,        "75":31.8667,
        "76":30.72,        "77":30.6133,        "78":30.2133,        "79":31.7467,        "80":33.8933,
        "81":32.72,        "82":33.1733,        "83":31.5867,        "84":32.6933,        "85":32.0667,
        "86":32.2933,        "87":30.7467,        "88":30.6933,        "89":32.48,        "90":31.08,
        "91":31.6133,        "92":31.72,        "93":30.3867,        "94":30.7067,        "95":29.9733,
        "96":31.96,        "97":32.44,        "98":30.2267,        "99":31.2533,        "100":30.2267
    }
}
```

# Naming
Originally this was named *fastq-stats*, but after a quick Google search (which I didn't do originally!) I found another [*fastq-stats*](https://github.com/ExpressionAnalysis/ea-utils/blob/master/clipper/fastq-stats.cpp) available from [ea-utils](http://expressionanalysis.github.io/ea-utils/). So I decided to rename it to *fastq-scan*, since its similar to the [*Scan*](
https://tvtropes.org/pmwiki/pmwiki.php/Main/EnemyScan) ability found in some video games/movies/tv etc... Basically it 'scans' a FASTQ and provides the user with otherwise hidden information about the FASTQ reads.
