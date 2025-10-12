from unittest import test
from phylo.tree import TreeNode
from phylo.tree import Tree
import phylo.upgma as upgma
from phylo.nj import neighbor_joining
from phylo.util import data_dir, join_path
import numpy as np
import time

def distances() -> np.ndarray:
    path = join_path(data_dir(), "distances.txt")
    rows: list[list[float]] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.replace(",", " ").split()
            rows.append([float(x) for x in parts])
    return np.asarray(rows, dtype=float)


def upgma_newick() -> str:
    # Newick notation of the tree created from 'distances.txt',
    # created via DendroUPGMA
    with open(join_path(data_dir(), "newick_upgma.txt"), "r") as file:
        newick: str = file.read().strip()
    return newick

def tree(distances):
    return upgma.upgma(distances)

d = np.asarray(distances(), dtype=np.int64)
t = tree(d)
newick = upgma_newick()

start = time.time()

@test
def test_upgma(tree: Tree, upgma_newick):
    """
    Compare the results of `upgma()` with DendroUPGMA.
    """
    ref_tree = Tree.from_newick(upgma_newick)
    # Cannot apply direct tree equality assertion because the distance
    # might not be exactly equal due to floating point rounding errors
    for i in range(Tree.__len__(tree)):
        for j in range(Tree.__len__(tree)):
            # Check for equal distances and equal topologies
            assert abs(tree.get_distance(i, j) - ref_tree.get_distance(i, j)) <= 1e-3
            assert tree.get_distance(i, j, topological=True) == ref_tree.get_distance(
                i, j, topological=True
            )
test_upgma(t, newick)

@test
def test_neighbor_joining():
    """
    Compare the results of `neighbor_join()` with a known tree.
    """
    dist = np.array([
        [ 0.0,  5.0,  4.0,  7.0,  6.0,  8.0],
        [ 5.0,  0.0,  7.0, 10.0,  9.0, 11.0],
        [ 4.0,  7.0,  0.0,  7.0,  6.0,  8.0],
        [ 7.0, 10.0,  7.0,  0.0,  5.0,  9.0],
        [ 6.0,  9.0,  6.0,  5.0,  0.0,  8.0],
        [ 8.0, 11.0,  8.0,  9.0,  8.0,  0.0],
    ])  # fmt: skip

    ref_tree = Tree(
        TreeNode(
            [
                TreeNode(
                    [
                        TreeNode(
                            children=[
                                TreeNode(index=0),
                                TreeNode(index=1),
                            ],
                            distances=[1.0, 4.0],
                        ),
                        TreeNode(index=2),
                    ],
                    distances=[1.0, 2.0],
                ),
                TreeNode(
                    children=[
                        TreeNode(index=3),
                        TreeNode(index=4),
                    ],
                    distances=[3.0, 2.0],
                ),
                TreeNode(index=5),
            ],
            distances=[1.0, 1.0, 5.0],
        )
    )

    test_tree = neighbor_joining(dist)

    assert test_tree == ref_tree
test_neighbor_joining()

@test
def test_distances(tree):
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert leaf.distance_to(tree.root) == dist
    # Example topological distances
    assert tree.get_distance(0, 19, True) == 9.0
    assert tree.get_distance(4, 2, True) == 10.0
test_distances(t)

end = time.time()

elapsed_ms = (end - start) * 1000

if __name__ == "__main__":
    print(f'codon\t{elapsed_ms}ms')