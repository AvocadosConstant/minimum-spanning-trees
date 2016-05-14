# Minimum Spanning Tree Algorithm Analysis
###### CS 375 (Design and Analysis of Algorithms) Final Project

***
### Problem to solve: 
Given a connected undirected graph, what is its minimum spanning tree?

### Algorithms to implement:  
Prim's algorithm and Kruskal's algorithm. Possibly Borůvka's algorithm and the reverse-delete algorithm time permitted.

### Questions to investigate:
Which algorithms are most efficient when? Do some perform better or worse depending on factors such as: Number of vertices, number of edges, maximum and minimum degree, planarity, simple graphs vs multi graphs, other graph invariants, etc.

### Implementation:
We will implement the algorithms and run them on a variety of different graphs with different properties. We will log each algorithms efficiency, effectiveness, and correctness and compare them to each other.

***
### Description:
Our code is written in C++. We have written four different implementations of minimum spanning tree algorithms: Prim's Algorithm (using adjacency lists and using adjacency matrices), and Kruskal's algorithm (using adjacency lists and adjacency matrices).


### Instructions:
Inside src/ run make to compile the code. To run the program, execute the following: 

``` ./MinimumSpanningTree <input-file> <output-file>```

With our current configuration, the <input-file> parameter is necessary, but unused. Currently, the program will generate a couple hundred random graphs on 20 vertices with random sparcity. The details on the graphs along with the runtimes of the algorithms of the graph are written line by line into the output file with the following comma delimited format:

```<number-of-vertices>, <number of edges>, <runtime of prims with adj list> <runtime of prim's with matrix>, <runtime of kruskal's with adj list> <runtime of kruskal's with matrix> ``` 

### Directories:
```input/``` contains several test graphs to run with the program. Currently the code ignores the input file and generates graphs instead.

```pres/``` contains resources for our presentation.

```src/``` contains our sourcecode ```mst.cpp``` and makefile.
