# PDRi

## Quick Start
```shell
java -jar PDRi.jar reference.fasta assembly.fasta
```

## Introduction

**PDR** (**P**airwise **D**istance **R**econstuction) is a genome assembly evaluation metric. It derives from a common concern in genetic studies, and takes completeness, contiguity, and correctness into consideration. PDRi is a implementation of it by integral.

## Installation

PDRi only needs Java(1.8 or above) and BWA to run. Please go to [release](https://github.com/taoistly/PDR/releases) to download PDRi.jar.

## Usage

```shell
java -jar PDRi.jar [Options] <reference genome> <assembly>

Options:
  --threads INT  Threads to use [default: CPU core]
  -k INT         Block size [default: 1000]
  -d PATH        Temporary folder for intermediate files [default: PDRTmp]
  -a TEXT        Executable path of aligner (BWA or minimap2) [default: bwa]
  -e INT         Maximum offset for two alignment segment to be jointed
                 [default: 0]
  -m INT         Minimum chromosome length (in bp) to summarize and report
                 alignment statistics. This doesn't change PDR result.
                 [default: 1% genome]
  -h, --help     Show this message and exit

```

## Output
PDRi will display a progress bar with estimated time during its execution. Finally, it will print a result like:
```
====== Finished ======
Genome payload: 2.923725537E9
PDR Total:      7.8534562466405898E18
PDR Ratio:      0.9187294255357177
[Success] elapsed time: 82m4s
```
`Genome payload` and `PDR total` are intermediate results for debugging purposes.

`PDR Ratio` is the only useful output. For example, a PDR Ratio = 0.9187294255357177 roughly means: 
> 91.873% pairs from the reference have correct distances in your assembly.

More generally (but less precisely), this roughly means:
> 91.873% information from the reference are also found in your assembly.

## FAQ

#### Q: How to specify BWA?

By default, PDRi tries to invoke BWA in PATH. If there is no BWA in PATH or you want to use specific BWA version, please provide executable BWA path by option `-a`, e.g. `java -jar PDRi.jar -a /usr/bin/bwa reference.fasta assembly.fasta `

#### Q: PDRi runs too slowly on large genome, how to accelerate it?

For large genome, BWA index building may cost a few hours. PDRi also supports Minimap2 as aligner. It is much faster than BWA without explicit index building, but its alignment precision is slightly lower. For rough test, you may want to use Minimap2 to replace BWA.  
 
## Citation
> Xie, L., & Wong, L. (2021). PDR: a new genome assembly evaluation metric based on genetics concerns. Bioinformatics, 37(3), 289-295.

## Feedback
If you encounter any problem or have feedback, please feel free to contact me (x86@u.nus.edu) or use github issue. I will reply ASAP.
