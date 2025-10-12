import numpy as np
from .tree import Tree, TreeNode
from typing import Optional, List

MAX_FLOAT = np.finfo(np.float64).max

def neighbor_joining(distances) -> Tree:

    # distances = np.asarray(distances, dtype=np.float64)

    n0 = int(distances.shape[0])
    n1 = int(distances.shape[1])

    if n0 != n1 \
        or not np.allclose(distances.T, distances):
            raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if n0 < 4:
        raise ValueError("At least 4 nodes are required")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")


    # Keep track on clustered indices
    nodes = np.array(
        [TreeNode(index=int(i)) for i in range(n0)]
    )
    # Indicates whether an index in the distance matrix has already been
    # clustered and the repsective rows and columns can be ignored
    is_clustered_v = np.full(
        distances.shape[0], False, dtype=bool
    )
    n_rem_nodes = \
        len(distances) - np.count_nonzero(np.asarray(is_clustered_v))
    # The divergence of of a 'taxum'
    # describes the relative evolution rate
    divergence_v = np.zeros(
        n0, dtype=np.float64
    )
    # Triangular matrix for storing the divergence corrected distances
    corr_distances_v = np.zeros(
        (n0,) * 2, dtype=np.float64
    )
    distances_v = distances.astype(np.float64, copy=True)


    nv0 = int(distances_v.shape[0])
    n_corr_v0 = int(corr_distances_v.shape[0])


    while True:
        # Calculate divergence
        for i in range(nv0):
            if is_clustered_v[i]:
                continue
            dist_sum: float = 0.0
            for k in range(nv0):
                if is_clustered_v[k]:
                    continue
                dist_sum += distances_v[i,k]
            divergence_v[i] = dist_sum
        
        # Calculate corrected distance matrix
        for i in range(nv0):
            if is_clustered_v[i]:
                    continue
            for j in range(i):
                if is_clustered_v[j]:
                    continue
                corr_distances_v[i,j] = \
                    (n_rem_nodes - 2) * distances_v[i,j] \
                    - divergence_v[i] - divergence_v[j]

        # Find minimum corrected distance
        dist_min: float = MAX_FLOAT
        i_min: int = -1
        j_min: int = -1
        for i in range(n_corr_v0):
            if is_clustered_v[i]:
                    continue
            for j in range(i):
                if is_clustered_v[j]:
                    continue
                dist = corr_distances_v[i,j]
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j
        
        # Check if all nodes have been clustered
        if i_min == -1 or j_min == -1:
            # No distance found -> all leaf nodes are clustered
            # -> exit loop
            break
        
        # Cluster the nodes with minimum distance
        # replacing the node at position i_min
        # leaving the node at position j_min empty
        # (is_clustered_v -> True)
        node_dist_i = 0.5 * (
            distances_v[i_min,j_min]
            + 1/(n_rem_nodes-2) * (divergence_v[i_min] - divergence_v[j_min])
        )
        node_dist_j = 0.5 * (
            distances_v[i_min,j_min]
            + 1/(n_rem_nodes-2) * (divergence_v[j_min] - divergence_v[i_min])
        )

        if n_rem_nodes > 3:
            # Clustering is not finished
            # -> Create a node with two children
            child_list: List[TreeNode] = [nodes[i_min], nodes[j_min]]
            dist_list: List[float] = [float(node_dist_i), float(node_dist_j)]

            nodes[i_min] = TreeNode(
                children=child_list, 
                distances=dist_list
            )

            # Mark position j_min as clustered
            nodes[j_min] = None
            is_clustered_v[j_min] = True
        else:
            # Clustering is finished
            # Combine ast three nodes into root node
            # Find the index of the remaining one of the three nodes
            # (other than i_min and j_min)
            is_clustered_v[i_min] = True
            is_clustered_v[j_min] = True
            # The index of the remaining one
            k = np.where(~np.asarray(is_clustered_v, dtype=bool))[0][0]
            node_dist_k = 0.5 * (
                distances_v[i_min,k] + distances_v[j_min,k]
                - distances_v[i_min,j_min]
            )

            root_children: List[TreeNode] = [nodes[i_min], nodes[j_min], nodes[k]]
            root_dists: List[float] = [float(node_dist_i), float(node_dist_j), float(node_dist_k)]
            
            root = TreeNode(
                children=root_children, 
                distances=root_dists
            )

            # Clustering is finished -> put into tree and return
            return Tree(root)
        
        # Update distance matrix
        # Calculate distances of new node to all other nodes
        for k in range(nv0):
            if not is_clustered_v[k] and k != i_min:
                dist = 0.5 * (
                    distances_v[i_min,k] + distances_v[j_min,k]
                    - distances_v[i_min,j_min]
                )
                distances_v[i_min,k] = dist
                distances_v[k,i_min] = dist

        # Update the amount of remaining nodes
        n_rem_nodes = \
            len(distances) - np.count_nonzero(np.asarray(is_clustered_v))