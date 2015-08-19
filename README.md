# CountErrors 1.1.1

### Install:
Run "make" in the program directory to compile

### "Usage:
```
[REQUIRED ARGUMENTS]

--fasta                 <string>                        Input reference sequence file
--bam                   <string>                        Input bam file
--target                <string>                        Input target interval bed file
--output                <string>                        Output file

[OPTIONAL ARGUMENTS]

--maq                   <int>                           Mapping quality threshold, Default [30]
--baq                   <int>                           Base quality threshold, Default [20]
--nodup                                                 Do not use reads that are marked as duplicate
--collapse_strand                                       Collapse the strand when counting errors
--collapse_end                                          Collapse the read end when counting errors
--collapse_context                                      Collapse the context bases(left and right bases of the reference) when counting errors
--help                                                  Print command line usage
```
