# NSC: Negative SkyCube

**NSC** is an index structure designed to speed up multi-dimensional skyline query answering. It consists of associating to a tuple, the subspaces where it is dominated. Those subspaces are encoded in a way that minimize memory consumption.

The **NSC** structure has been presented at **CIKM 2016** conference through the paper [Computing and summarizing the negative Skycube](https://dl.acm.org/citation.cfm?doid=2983323.2983759). Additionnal materials, including incremental maintenance of **NSC**, have been published in **Information Systems** through the paper [The negative skycube](https://www.sciencedirect.com/science/article/abs/pii/S0306437919304958). 

If you consider **NSC** in your work, please cite:

    @article{alami2020negative,
      title={The negative skycube},
      author={Alami, Karim and Hanusse, Nicolas and Kamnang-Wanko, Patrick and Maabout, Sofian},
      journal={Information Systems},
      volume={88},
      pages={101443},
      year={2020},
      publisher={Elsevier}
    }

**NSC**'s efficiency wrt building time, memory usage and skyline query answering has been shown in that paper.

Afterward we studied its incremental maintenance. 

### Repository organization

This repository includes **NSC** sources as well as state of the art methods to build and answer skyline queries.

* dataset/: real datasets used for experimentations

* src/: source code with implemention of 

  * NSC
  * NSC with counters
  * BSkyTree (authors version)
  * Naive method
  * Compressed Skycube (CSC) (our own version)

### Building and answering Skycube

To begin, please clone this repository. This software requires C++ compiler.

Next, change current directory to ./src and compile the source code with [Makefile](https://github.com/karimalami7/NSC/blob/master/src/makefile).

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

6. List of methods to be run:

| Methods  | Arguments  |
|---|---|
| NSC  |  NSC |
| NSC with counters  |  NSCwM |
| BSkyTree | TREE|
| Naive method  | NAIF |
| Compressed SkyCube | CSC |

Outputs are (i) building time (for methods that build a structure), (ii) skyline query execution time for one random query, (iii) Skycube query execution time.

## Maintenance of NSC

Experiment maintenace of **NSC** through NSCwM argument.

The interactive menu allows to: 

1. Delete one tuple from the initial data.
2. Insert one tuple into the initial data.
3. Query NSC to retrieve the skyline of a subspace or the whole skycube.
4. Delete multiple tuples from the initial data.
5. Insert multiple tuples into the initial data.
6. Compute impact of every tuple in the **Topmost**.

