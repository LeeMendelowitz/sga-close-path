SGA - String Graph Assembler
==========================================

SGA is a de novo genome assembler based on the concept of string graphs. The major goal of SGA is to be very memory efficient, which is achieved by using a compressed representation of DNA sequence reads.

This repository has modified SGA v0.9.35 to add the `sga close-path` feature for 
removing unsupported edges from the string graph using distance estimates generated
by paired reads.

The source code updates for `sga close-path` can be found in these directories:

 - [Executable](https://github.com/LeeMendelowitz/sga-close-path/blob/close_path/src/SGA/close-path.cpp)
 - [Graph Search Code](https://github.com/LeeMendelowitz/sga-close-path/tree/close_path/src/GraphSearch)


--------
For installation and usage instructions see src/README

For running examples see src/examples and the [sga wiki](https://github.com/jts/sga/wiki)

For questions or support contact jared.simpson --at-- sanger.ac.uk
