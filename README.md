# NSC: Negative SkyCube

NSC is an index structure designed to speed up multi-dimensional skyline query answering. It consists of associating to a tuple, the subspaces where it is dominated. Those subspaces are encoded in a way that minimize memory consumption.

NSC has been presented at **CIKM 2016** conference through the paper [Computing and summarizing the negative Skycube](https://dl.acm.org/citation.cfm?doid=2983323.2983759).

Its efficiency wrt building time, memory usage and skyline query answering has been shown in that paper.

Afterward we studied its incremental maintenance. 

### Repository organization

This repository includes NSC sources as well as state of the art methods to build and answer skyline queries.

* dataset/: real datasets used for experimentations

* src/: source code with implemention of 

  * NSC
  * NSC with counters
  * BSkyTree (authors version)
  * Naive method
  * Compressed Skycube (CSC) (our own version)

### Building and answering Skycube

To begin, please clone this repository. This software requires C++ compiler.

Next, compile the source code with [Makefile](https://github.com/karimalami7/NSC/blob/master/src/makefile).

Then, identify the input dataset and methods to be run with arguments to be set in [run](https://github.com/karimalami7/NSC/blob/master/src/run.sh) file.

Arguments:

1. Data type:

* INDE: Independent data.

* ANTI: Anti-correlated data.

* CORR: Correlated data.

* other: when using a real dataset.

2. In case of real dataset, input the path of the file. In case of generated dataset, define the number of distinct values of each dimension.

3. Size of the input data.

4. Number of dimensions. 

5. Number of parallel threads.

6. List of methods to be run: (NSC, NSCwM, TREE, NAIF, CSC).

Outputs are (i) building time (for methods that build a structure), (ii) skyline query execution time for one random query, (iii) Skycube query execution time.

