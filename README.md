# PDRi

## Quick Start
```shell
java -jar PDRi.jar reference.fasta assembly.fasta
```

## Introduction

**PDR** (**P**airwise **D**istance **R**econstuction) a genome assembly evaluation metric. It derives from a common concern in genetic studies, and takes completeness, contiguity, and correctness into consideration. PDRi is a implementation of it by integral.

## Installation

PDRi only needs Java(1.8 or above) and BWA to run. Please go to [release](https://github.com/taoistly/PDR/releases) to download PDRi.jar.

## Usage

```shell
java -jar PDRi.jar [Options] <reference genome> <assembly>

Options:
  --threads INT  Threads to use [default: CPU core]
  -k INT         Bin size [default: 1000]
  -d PATH        Temporary folder for intermediate results
  -a TEXT        Executable path of aligner (BWA or minimap2) [default: bwa]
  -e INT         Maximum error for two alignment segment to be jointed
                 [default: 0]
  -m INT         Minimum chromosome length (in bp) to summarize and report.
                 This doesn't effect result. [default: 1% genome]
  -h, --help     Show this message and exit
```

## FAQ

#### Q: How to specify BWA?

By default, PDRi tries to invoke BWA in PATH. If there is no BWA in PATH or you want to use specific BWA version, please provide executable BWA path by option `-a`, e.g. `java -jar PDRi.jar -a /usr/bin/bwa reference.fasta assembly.fasta `

##### Q: PDRi runs too slowly on large genome, how to accelerate it?

For large genome, bwa index building may cost a few hours. PDRi also supports Minimap2 as aligner. It is much faster than BWA without explicit index building, but its alignment precision is slightly lower. For rough test, you may want to use Minimap2 to replace BWA.  
 
## Citation
To be published

## Feedback
If you encounter any problem or have feedback, please feel free to contact me (x86@u.nus.edu). I will try to reply ASAP.