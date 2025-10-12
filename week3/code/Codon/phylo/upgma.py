import numpy as np

from .tree import Tree, TreeNode

MAX_FLOAT = np.finfo(np.float64).max


# @cython.boundscheck(False)
# @cython.wraparound(False)
def upgma(distances) -> Tree:
    distances = np.asarray(distances, dtype=np.float64)

    n0 = int(distances.shape[0])
    n1 = int(distances.shape[1])

    if n0 != n1 \
        or not np.allclose(distances.T, distances):
            raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if (distances < 0.0).any():
        raise ValueError("Distances must be positive")


    # Keep track on clustered indices
    nodes = np.array(
        [TreeNode(index=i) for i in range(n0)]
    )
    # Indicates whether an index in the distance matrix has already been
    # clustered and the repsective rows and columns can be ignored
    is_clustered_v = np.full(
        n0, False, dtype=bool
    )
    # Number of indices in the current node (cardinality)
    # (required for proportional averaging)
    cluster_size_v = np.ones(
        n0, dtype=np.float64
    )
    # Distance of each node from leaf nodes,
    # used for calculation of distance to child nodes
    node_heights = np.zeros(
        n0, dtype=np.float64
    )


    # Cluster indices
    distances_v = distances.astype(np.float64, copy=True)
    # Exit loop via 'break'
    while True:

        # Find minimum distance
        dist_min = MAX_FLOAT
        i_min = -1
        j_min = -1

        nv0 = int(distances_v.shape[0])

        for i in range(nv0):
            if is_clustered_v[i]:
                continue
            for j in range(i):
                if is_clustered_v[j]:
                    continue
                dist = distances_v[i,j]
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j
        
        if i_min == -1 or j_min == -1:
            # No distance found -> all leaf nodes are clustered
            # -> exit loop
            break
        
        # Cluster the nodes with minimum distance
        # replacing the node at position i_min
        # leaving the node at position j_min empty
        # (is_clustered_v -> True)
        height = dist_min/2

        child_list: list[TreeNode] = [nodes[i_min], nodes[j_min]]
        dist_list: list[float] = [float(height-node_heights[i_min]), float(height-node_heights[j_min])]
        
        nodes[i_min] = TreeNode(
            children = child_list,
            distances = dist_list
        )
        node_heights[i_min] = height
        # Mark position j_min as clustered
        nodes[j_min] = None
        is_clustered_v[j_min] = True
        # Calculate arithmetic mean distances of child nodes
        # as distances for new node and update matrix
        for k in range(nv0):
            if not is_clustered_v[k] and k != i_min:
                mean = (
                    (
                          distances_v[i_min,k] * cluster_size_v[i_min]
                        + distances_v[j_min,k] * cluster_size_v[j_min]
                    ) / (cluster_size_v[i_min] + cluster_size_v[j_min])
                )
                distances_v[i_min,k] = mean
                distances_v[k,i_min] = mean
        # Updating cluster size of new node
        cluster_size_v[i_min] = cluster_size_v[i_min] + cluster_size_v[j_min]
    

    # As each higher level node is always created on position i_min
    # and i is always higher than j in minimum distance calculation,
    # the root node must be at the last index
    return Tree(nodes[len(nodes)-1])