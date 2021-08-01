# vertexCover-with-kernelization
two implementations to VertexCover problem - recursive and Linear Programming

## Introduction
Kernelization is an approach for pre-processing algorithms. It is usually achieved by
replacing the algorithm's original input with a smaller one, called a kernel.
Parameterized complexity is used to analyze the result of the kernelization
algorithms.

A **vertex cover** of a graph G is a set S of vertices of G such that every edge of G has at least one member of S as an endpoint.
The vertex cover problem is an NP-complete problem.

## How to use
#### kernelization_algorithm(G, k):
  This function gets graph G and integer k and returns a smaller equals instance of the same graph or NULL if the is no solution.
  ![alt text](https://github.com/ofekbenda/vertexCover-with-kernelization/tree/master/imgs_for_readme/kernelization1.png?raw=true)

#### kernelization_algorithm_opt(G, k):
  Improvement of kernelization_algorithm(G, k) by union reduction 1 and 2.
    ![alt text](https://github.com/ofekbenda/vertexCover-with-kernelization/tree/master/imgs_for_readme/kernelization_improvement.png?raw=true)
    
#### recursive_algorithem(G, k):
  This function gets graph G and integer k and returns a solution. The correctness of this algorithm is based on the fact that for each edge (u,v) in G, either
u is in the vertex cover or v is in the vertex cover.
![alt text](https://github.com/ofekbenda/vertexCover-with-kernelization/tree/master/imgs_for_readme/recursive.png?raw=true)

#### recursive_algorithem_opt(G, k):
  Improvement of recursive_algorithem(G, k).
  For a vertex v, any vertex cover of G must contain either v or all of its neighbors N(v).
  Since we pick the vertex with the maximum degree, we guarantee to cover each iteration's maximum number of edges.  

#### LPVC(G):
  kernelization algorithm based on LP. For a graph G, we obtain the following linear programming instance:
  ![alt text](https://github.com/ofekbenda/vertexCover-with-kernelization/tree/master/imgs_for_readme/LP.png?raw=true)
  
#### LP_parm(G, k):
  Given a graph G and an integer k, we ask for a vertex cover of G of size at most k, but instead of seeking an FPT
algorithm parameterized by k, the parameter now is k − vc∗(G).
    ![alt text](https://github.com/ofekbenda/vertexCover-with-kernelization/tree/master/imgs_for_readme/VC_above_LP.png?raw=true)
    ![alt text](https://github.com/ofekbenda/vertexCover-with-kernelization/tree/master/imgs_for_readme/vc4.png?raw=true)
    ![alt text](https://github.com/ofekbenda/vertexCover-with-kernelization/tree/master/imgs_for_readme/vc5.png?raw=true)

    ![alt text](https://github.com/ofekbenda/vertexCover-with-kernelization/tree/master/imgs_for_readme/vc3.8.png?raw=true)
    
    
  **main contain examples how to use the functions**

