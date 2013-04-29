# Multidimensional Binpacking Constraint
This directory contains an implementation of a simple constraint
decomposition for the `multibin_packing` constraint.

The decomposition is based on the search for maximal cliques in a conflict graph.
The maximal cliques are used for consistency checking and for posting `alldifferent` constraints on subsets of assignment variables. In addition, in case all the bins have the same capacity with respect to the same dimensions, it is possible to post symmetry breaking constraints.

The conflict graph *G=(V,E)* has an edge for every pair of items that does not fit in the same bin for at least one dimension. A maximal clique in *G* corresponds to a subset of items that must be placed into different bins. The maximum clique of *G* gives a lower bound on the number of bins necessary to place every item.

## Istructions
In order to compile this example you need the following external libraries:

1. Gecode v4.0.0 (www.gecode.org)
2. Cliquer v1.21 (http://users.tkk.fi/pat/cliquer.html)

You must specify in the **Makefile** the directories where you have installed these two
libraries, and then you just have to type `make`. To test that every thinks works, just type:

`$ ./example instance6-18-6-20_1.txt`

The output should be something like:

`2,4,5,4,4,0,0,3,0,2,3,1,5,2,5,1,1,3,`

`Nodes 94 Memory 72960 Time 9.161`

The first list of integers is the certificate of the solution (the first item is packed into bin number 2; the second item in bin number 4; and so on ...).
The other numbers report the number of nodes, the memory consumed (in byte), and run time (in seconds).

## References
We wrote a short abstract about the decomposition approach impelmented in this example:

[1] S. Gualandi and M. Lombardi _A Simple and Effective Decomposition for the Multidimensional Constraints_. Technical Report, April, 2013.

The report is available from the authors by request.