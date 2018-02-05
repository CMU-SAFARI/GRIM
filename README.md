# GRIM-Filter 

GRIM-Filter is an algorithm optimized to exploit 3D-stacked memory systems that integrate computation within a logic layer stacked under memory layers, to perform processing-in-memory (PIM). GRIM-Filter quickly filters seed locations by 1) introducing a new representation of coarse-grained segments of the reference genome, and 2) using massively-parallel in-memory operations to identify read presence within each coarse-grained segment.

Our code baseline is taken from [mrFAST_v2.6.1.0](http://mrfast.sourceforge.net/), which is described in detail in the following publications:
* [C. Alkan et al., Personalized copy number and segmental duplication maps using next-generation sequencing, Nature Genetics](https://www.ncbi.nlm.nih.gov/pubmed/19718026) 
* [H. Xin, Accelerating read mapping with FastHASH, BMC Genomics](https://www.ncbi.nlm.nih.gov/pubmed/23369189)

While we use mrFAST as a baseline, GRIM-Filter can be adapted to run with any other read mapper. 

The algorithm of GRIM-Filter is described at: [J.S. Kim et al., GRIM-Filter: Fast Seed Location Filtering in DNA Read Mapping using Processing-in-Memory Technologies, To appear in BMC Genomics](https://arxiv.org/pdf/1711.01177.pdf)

## Prerequisites 

In order to run GRIM-Filter, have the following files: 
* Human Genome FASTA file (e.g., [Human_g1k_v37 Genome](ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/.../human_g1k_v37.fasta.gz))
* Read Sequence data sets (FASTA file) 

## Getting Started

To build mrFAST with GRIM-Filter, simply do: 

```
$ make 
```

To build the hash table used by mrFAST, run the following command:
```
./mrfast --index <Genome FASTA File>
```
There is more information on the parameters for hash table generation in the [mrFAST User Manual](refgen.fasta).

To build the bitvectors that are referenced by GRIM-Filter, run the following
command: 
```
./mrfast --index <Genome FASTA File> -t 0 -k <Number of Bins> -b <Token Size> -f <Number of Tokens the Bitvector can Count (1)>
```

This will generate a <Genome FASTA Filename>.bv file in the same directory as
your Genome FASTA File. 

You can then use the bitvectors by running mrfast with the following command: 
```
./mrfast --search <Genome FASTA File> -b <Token Size> -t 1 -e <error Tolerance (%)> -k <Number of Bins> -q 1 --seq <Read Sequences FASTA File>
```

## Contributors 

* **Jeremie S. Kim** (Carnegie Mellon University) 

