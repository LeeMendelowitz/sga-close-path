9/29/2012
Working directory for source code to perform graph search in string graph


METHOD FOR FINDING PATHS from A to B, with orientation given for A and B:
    - paths must be between d_min and d_max in length
 
1. Perform BFS from A, accepting nodes which start <= d_max/2.
2. Perfrom BFS from B, accepting nodes which start <= d_max/2.
3. Find nodes that are both A and B nodes. If there are none, then there is no path. From these nodes, 
walk in both directions, using only edges that are either A edges or B edges. Place all nodes encountered into a new set.
This is the set of nodes for the subgraph.
4. Construct the subgraph, considering all edges between selected nodes. Note that some edges may be included in this step that were not
generated in the BFS. (These would be edges to nodes that appeared earlier in the BFS, and therefore were not taken during BFS. These
edges create cycles).
5. Use the search tree algorithm or something similar to it to get all walks.


