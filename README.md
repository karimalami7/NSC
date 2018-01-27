# NSC: Negative SkyCube

NSC is an index structure designed to speed up multi-dimensional skyline query answering. It consists of associating to a tuple, the subspaces where it is dominated. Those subspaces are encoded in a way that minimize memory consumption.

NSC has been accepted in the conference **CIKM 2016** through the paper [Computing and summarizing the negative Skycube](https://dl.acm.org/citation.cfm?doid=2983323.2983759).

Its efficiency wrt building time, memory usage and skyline query answering has been shown in this paper.

Afterward we studied its incremental maintenance. Few state of the art methods worked on managing the skyline upon deletions and insertions.

### Repository organization

This repository includes NSC sources as well as state of the art methods to build and answer skyline queries.

...

### Building and answering Skycube

...
