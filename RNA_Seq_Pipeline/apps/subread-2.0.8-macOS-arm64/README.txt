Installation
--------------
You may use Subread package by directly downloading a binary release suitable for your operating system (no compilation is required), or by installing it to your computer from the source.

Here we describe how to install it from the source. You will need to download the source release of the package from http://subread.sourceforge.net/. The source release includes a keyword 'source' in the name of the tar ball. Uncompress the tar ball, enter the 'src' directory and issue one of the following commands to build it.

*** for Linux OS, use command:

make -f Makefile.Linux

*** For Windows, use command:

make -f Makefile.Windows

When building Subread on a Windows computer, the ``gcc" and ``make" programs will be used. The paths to these two programs need to be included in your environment variable PATH so that these two programs can be accessed.

An alternative way to build Subread on a Windows computer is to use Mingw-w64. This is often a more reliable approach. After installing Mingw-w64, you can just use the following command to build Subread. Note that the path to the two mingw_64 executables ``mingw32-make.exe" and ``gcc.exe" should also be included in your environment variable PATH.

mingw32-make -f Makefile.Windows


*** For Mac OS, use command:

make -f Makefile.MacOS

*** For FreeBSD OS, use command:

gmake -f Makefile.FreeBSD

If the build is successful, a new directory called 'bin' will be created under the home directory of the package (ie. one level up from 'src' directory). The 'bin directory contains all the generated executables. To enable easy access to these executables, you may copy the executables to a system directory such as '/usr/bin' or add the path to the executables to your search paths that are usually stored in your environment variable PATH.

Content
--------------
annotation		Directory including NCBI RefSeq gene annotations for multiple versions of human and mouse genomes.
bin			Directory including executables after compilation (or directly available from a binary release). 
doc			Directory including the Subread Users Guide.
LICENSE			The license agreement for using this package.
README.txt		This file.
src			Directory including source code (binary releases do not have this directory).
test			Directory including test data and scripts.

Citation
--------------
If you use Subread or Subjunc aligners, please cite:

Liao Y, Smyth GK and Shi W. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10):e108, 2013

If you use the featureCounts program, please cite:

Liao Y, Smyth GK and Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014

If you use Subread and featureCounts to map and quantify RNA-seq data, please cite:

Liao Y, Smyth GK and Shi W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Research, 47(8):e47, 2019

Mailing lists
--------------
Please post your questions/suggestions to Bioconductor support site (https://support.bioconductor.org/) or Subread google group (https://groups.google.com/forum/#!forum/subread).
