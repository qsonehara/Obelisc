# Obelisc
`obelisc` (**Ob**s**e**rvational **li**nkage **sc**an) is an identical-by-descent (IBD) mapping tool based on the SNP streak approach.

<p><img src="https://github.com/qsonehara/Obelisc/blob/images/Obelisc_overview.png" width=800px></p>

## Overview

`obelisc` is a command line tool to perform two effective nonparametric linkage analyses (i.e., SNP streak-based IBD mapping and runs of homozygosity (ROH) mapping). `obelisc` also supports an intuitive visualization of the analytical results.

`obelisc` requires no predefined inheritance model with parameters and is applicable to both dominant and recessive traits.



## Requirements

- python 3.x
- numpy
- pandas
- matplotlib
- pysnptools
- argparse
- Cython (if you want to build LOCH_MappingTools from cython code)



## Installation

`obelisc` can be installed using `pip`.

```
pip3 install git+https://github.com/qsonehara/Obelisc
```



## Example usage

When you install `obelisc`, you get a simple command-line interface. You can perform SNP streak-based IBD mapping and ROH mapping by running the command `obelisc` as follows:

```
obelisc example
```

Here, `example` is the prefix of a PLINK binary fileset (i.e., example.fam, example.bim, and example.bed) and a text file specifying case individuals (i.e., example.case).

You can run `obelisc` with the example dataset provided as follows:

```
obelisc hapmap3_r2_b37_fwd.consensus.qc.poly.JPT.example
```

You can specify the filename of another case-specifying text file, the size of the scanning window, the minimum size of output stretches, the filename of mapping output, and the filename of the diagram of results as follows;

```
obelisc example -c case_file.case -w 1500 -s 1500 -o example_ibdmapping -f example_out.pdf
```

#### Options

| Option name                    | Descriptions                                                 | Default                            |
| ------------------------------ | ------------------------------------------------------------ | ---------------------------------- |
| `--case`, `-c`                 | the case file name                                           | {the first argument}.case          |
| `--ibd-win-kb`, `-w`           | the sliding window size during IBD mapping (kb) (regarding SNP streak mapping) | 1500                               |
| `--ibd-str-kb`, `-s`           | the minimum size of outputted IBD stretches (kb) (regarding SNP streak mapping) | 1500                               |
| `--ibd-inconsistent-snp`, `-i` | the maximum number of markers in a scanning window inconsistent with the SNP streak condition  (regarding SNP streak mapping) | 1                                  |
| `--ibd-min-snp`, `-m`          | the minimum number of markers a scanning window requires  (regarding SNP streak mapping) | 25                                 |
| `--ibd-win-gap`, `-g`          | the maximum length a pair of nearest markers can be apart (kb) (regarding SNP streak mapping) | 1000                               |
| `--ibd-out`, `-o`              | the prefix of the output text file for SNP streak mapping results | generated with option values above |
| `--roh-win-kb`, `-W`           | the sliding window size during ROH mapping (kb)              | 1500                               |
| `--roh-str-kb`, `-S`           | the minimum size of outputted ROH stretches (kb)             | 1500                               |
| `--roh-inconsistent-snp`, `-I` | the maximum number of heterozygous markers in a scanning window | 1                                  |
| `--roh-min-snp`, `-M`          | the minimum number of markers a scanning window requires  (regarding ROH mapping) | 25                                 |
| `--roh-win-gap`, `-G`          | the maximum length a pair of nearest markers can be apart (kb) (regarding ROH mapping) | 1000                               |
| `--roh-out`, `-O`              | the prefix of the output text file for ROH mapping results   | generated with option values above |
| `--point`, `-p`                | Show the statuses of all of the IBD markers in the cases.    | false                              |
| `--fig`, `-f`                  | the prefix of the output PDF file for both mapping results   | generated with option values above |



## Contact

Kyuto Sonehara ([qsonehara@sg.med.osaka-u.ac.jp](mailto:qsonehara@sg.med.osaka-u.ac.jp))
