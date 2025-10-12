import numpy as np
import copy
from .copyable import Copyable
from typing import Optional, List


@extend
class set:
    def __hash__(self):
        MAX = int.MAX
        MASK = 2 * MAX + 1
        n = len(self)
        h = 1927868237 * (n + 1)
        h &= MASK
        for x in self:
            hx = hash(x)
            h ^= (hx ^ (hx << 16) ^ 89869747)  * 3644798167
            h &= MASK
        h = h * 69069 + 907133923
        h &= MASK
        if h > MAX:
            h -= MASK + 1
        if h == -1:
            h = 590923713
        return h

class TreeNode:
    _index: int
    _distance: float
    _is_root: bool
    _parent: Optional[TreeNode]
    _children: List[TreeNode]

    @overload
    def __init__(
            self, 
            children: List[TreeNode]=[], 
            distances: Optional[List[float]]=None, 
            index: Optional[int]=None):
        
        self._is_root: bool = False
        self._distance: float = 0.0
        self._parent: Optional[TreeNode] = None
        self._children = []

        if index is None:
            # Node is intermediate -> has children
            if children is None or distances is None:
                raise TypeError(
                    "Either reference index (for terminal node) or "
                    "child nodes including the distance "
                    "(for intermediate node) must be set"
                )
            for item in children:
                if not isinstance(item, TreeNode):
                    raise TypeError(
                        f"Expected 'TreeNode', but got '{type(item).__name__}'"
                    )
            for item in distances:
                if not isinstance(item, float):
                    raise TypeError(
                        f"Expected 'float', "
                        f"but got '{type(item).__name__}'"
                    )
            if len(children) == 0:
                raise TreeError(
                    "Intermediate nodes must at least contain one child node"
                )
            if len(children) != len(distances):
                raise ValueError(
                    "The number of children must equal the number of distances"
                )
            for i in range(len(children)):
                for j in range(len(children)):
                    if i != j and children[i] is children[j]:
                        raise TreeError(
                            "Two child nodes cannot be the same object"
                        )
            self._index = -1
            self._children = list(children)
            for child, distance in zip(children, distances):
                child._set_parent(self, distance)
        elif index < 0:
            raise ValueError("Index cannot be negative")
        else:
            # Node is terminal -> has no children
            if children is not None or distances is not None:
                raise TypeError(
                    "Reference index and child nodes are mutually exclusive"
                )
            self._index = index
            self._children = []
    
    @overload
    def __init__(self, index: int):
        if index < 0:
            raise ValueError("Index cannot be negative")
        self._index = index
        self._is_root = False
        self._parent = None
        self._distance = 0.0
        self._children = []

    def _set_parent(self, parent: Optional[TreeNode], distance: float):
        if parent is None:
            raise ValueError("Parent should not be None")
        
        if self._parent is not None or self._is_root:
            raise TreeError("Node already has a parent")
        self._parent = parent
        self._distance = distance
    
    def copy(self):
        if self.is_leaf():
            return TreeNode(index=self._index)
        else:
            distances = [child.distance for child in self._children]
            children_clones = [child.copy() for child in self._children]
            return TreeNode(children_clones, distances)

    @property
    def index(self):
        return None if self._index == -1 else self._index
    
    @property
    def children(self) -> List[TreeNode]:
        return self._children
    
    @property
    def parent(self):
        return self._parent
    
    @property
    def distance(self):
        return None if self._parent is None else self._distance

    def is_leaf(self):
        return False if self._index == -1 else True
    
    def is_root(self):
        return bool(self._is_root)
    
    def as_root(self):
        if self._parent is not None:
            raise TreeError("Node has parent, cannot be a root node")
        self._is_root = True
    
    def distance_to(self, node, topological=False):
        # Sum distances until LCA has been reached
        distance = 0.0
        lca = self.lowest_common_ancestor(node)
        if lca is None:
            raise TreeError("The nodes do not have a common ancestor")
        current_node = self
        while current_node is not lca:
            p = current_node._parent
            if p is None:
                raise TreeError("Broken parent chain on self-side before reaching LCA")
            if topological:
                distance += 1.0
            else:
                distance += current_node._distance
            current_node = p
        current_node = node
        while current_node is not lca:
            p = current_node._parent
            if p is None:
                raise TreeError("Broken parent chain on self-side before reaching LCA")
            if topological:
                distance += 1.0
            else:
                distance += current_node._distance
            current_node = p
        return distance
    
    def lowest_common_ancestor(self, node):
        # Create two paths from the leaves to root
        self_path = _create_path_to_root(self)
        other_path = _create_path_to_root(node)
        lca = None
        # Reverse Iteration through path (beginning from root)
        # until the paths diverge
        for i in range(-1, -min(len(self_path), len(other_path))-1, -1):
            if self_path[i] is other_path[i]:
                # Same node -> common ancestor
                lca = self_path[i]
            else:
                # Different node -> Not common ancestor
                # -> return last common ancewstor found
                break
        return lca
    
    def get_indices(self):
        return np.array(
            [leaf._index for leaf in self.get_leaves()], dtype=np.int64
        )

    def get_leaves(self):
        leaf_list: List[Optional[TreeNode]] = []
        _get_leaves(self, leaf_list)
        return leaf_list
    
    def get_leaf_count(self):
        return _get_leaf_count(self)
    
    def to_newick(self, labels: Optional[List[str]]=None, include_distance: bool=True, 
                  round_distance: Optional[int]=None):
        if self.is_leaf():
            if labels is not None:
                for label in labels:
                    label = str(labels[self._index])
                    # Characters that are part of the Newick syntax
                    # are illegal
                    illegal_chars = [",",":",";","(",")"]
                    for char in illegal_chars:
                        if char in label:
                            raise ValueError(
                                f'Label {label} contains '
                                f'illegal character {char}'
                            )
            else:
                label = str(self._index)
            if include_distance:
                if round_distance is None:
                    return f'{label}:{self._distance}'
                else:
                    precision = round(self._distance, round_distance)
                    return f'{label}:{precision}'
            else:
                return f'{label}'
        else:
            # Build string in a recursive way
            child_strings = [child.to_newick(
                labels, include_distance, round_distance
            ) for child in self._children]
            joined_child_strings = ",".join(child_strings)
            if include_distance:
                if round_distance is None:
                    return f'({joined_child_strings}):{self._distance}'
                else:
                    precision = round(self._distance, round_distance)
                    return (
                        f'({joined_child_strings}):{precision}'
                    )
            else:
                return f'({joined_child_strings})'
    
    @staticmethod
    def from_newick(newick, labels=None):
        subnewick_start_i = -1
        subnewick_stop_i  = -1
        level = 0
        
        # Ignore any whitespace
        newick = "".join(newick.split())

        # Find brackets belonging to sub-newick
        # e.g. (A:0.1,B:0.2):0.5
        #      ^           ^
        for i in range(len(newick)):
            char = newick[i]
            if char == "(":
                subnewick_start_i = i
                break
            if char == ")":
                raise InvalidFileError("Bracket closed before it was opened")
        for i in reversed(range(len(newick))):
            char = newick[i]
            if char == ")":
                subnewick_stop_i = i+1
                break
            if char == "(":
                raise InvalidFileError("Bracket was opened but not closed")
        
        if subnewick_start_i == -1 and subnewick_stop_i == -1:
            # No brackets -> no sub-newwick -> Leaf node
            label_and_distance = newick
            try:
                label, distance_str = label_and_distance.split(":")
                distance = float(distance_str)
            except ValueError:
                distance = 0.0
                label = label_and_distance
            index = int(label) if labels is None else labels.index(label)
            return TreeNode(index=index), distance
        
        else:
            # Intermediate node
            if subnewick_stop_i == len(newick):
                # Node with neither distance nor label
                distance = 0.0
            else:
                label_and_distance = newick[subnewick_stop_i:]
                try:
                    _label_unused, distance_str = label_and_distance.split(":")
                    distance = float(distance_str)
                except ValueError:
                    # No colon -> No distance is provided
                    distance = 0.0
                # Label of intermediate nodes is discarded 
                distance = float(distance)
            
            subnewick = newick[subnewick_start_i+1 : subnewick_stop_i-1]
            if len(subnewick) == 0:
                raise InvalidFileError(
                    "Intermediate node must at least have one child"
                )
            # Parse childs
            # Split subnewick at ',' if ',' is at current level
            # (not in a subsubnewick)
            comma_pos = []
            for i, char in enumerate(subnewick):
                if char == "(":
                    level += 1
                elif char == ")":
                    level -= 1
                elif char == ",":
                    if level == 0:
                        comma_pos.append(i)
                if level < 0:
                    raise InvalidFileError(
                        "Bracket closed before it was opened"
                    )
        
            children = []
            distances = []
            # Recursive tree construction
            for i, pos in enumerate(comma_pos):
                if i == 0:
                    # (A,B),(C,D),(E,F)
                    # -----
                    child, dist = TreeNode.from_newick(
                        subnewick[:pos], labels=labels
                    )
                else:
                    # (A,B),(C,D),(E,F)
                    #       -----
                    prev_pos = comma_pos[i-1]
                    child, dist = TreeNode.from_newick(
                        subnewick[prev_pos+1 : pos], labels=labels
                    )
                children.append(child)
                distances.append(dist)
            # Node after last comma
            # (A,B),(C,D),(E,F)
            #             -----
            if len(comma_pos) != 0:
                child, dist = TreeNode.from_newick(
                    subnewick[comma_pos[-1]+1:], labels=labels
                )
            else:
                # Single child node:
                child, dist = TreeNode.from_newick(
                    subnewick, labels=labels
                )
            children.append(child)
            distances.append(dist)
            return TreeNode(children, distances), distance

    def __str__(self):
        return self.to_newick()
    
    def __eq__(self, item):
        if not isinstance(item, TreeNode):
            return False
        node = item
        if self._distance != item._distance:
            return False
        if self._index != -1:
            return self._index == item._index
        else:
            if set(self._children) != set(node._children):
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        # Order of children is not important -> set
        children_set = set(self._children) \
                       if self._children is not None else None
        return hash((self._index, children_set, self._distance))

class Tree(Copyable):
    _leaves: List[Optional[TreeNode]]
    _root: Optional[TreeNode]
    leaves_unsorted: List[Optional[TreeNode]]
    leaf_count: int

    def __init__(self, root):
        if root is None:
            raise ValueError("Root cannot be none")
        
        root.as_root()
        self._root = root
        
        leaves_unsorted: List[Optional[TreeNode]] = self._root.get_leaves()
        leaf_count: int = len(leaves_unsorted)
        indices = np.array(
            [leaf.index for leaf in leaves_unsorted]
        )
        self._leaves: List[Optional[TreeNode]] = [None] * leaf_count
        
        for i in range(len(indices)):
            index = indices[i]
            if index >= leaf_count or index < 0:
                raise TreeError("The tree's indices are out of range")
            self._leaves[index] = leaves_unsorted[i]
    
    def __copy_create__(self):
        return Tree(self._root.copy())
    
    @property
    def root(self):
        return self._root
    
    @property
    def leaves(self):
        return copy.copy(self._leaves)

    def get_distance(self, index1, index2, topological=False):
        return self._leaves[index1].distance_to(
            self._leaves[index2], topological
        )
    
    def to_newick(self, labels: Optional[List[str]]=None, include_distance: bool=True, 
                  round_distance: Optional[int]=None):
        return self._root.to_newick(
            labels, include_distance, round_distance
        ) + ";"
    
    @staticmethod
    def from_newick(newick: str, labels=None):
            
        newick = newick.strip()
        if len(newick) == 0:
            raise InvalidFileError("Newick string is empty")

        # Remove terminal colon as required by 'TreeNode.from_newick()'
        if newick[-1] == ";":
            newick = newick[:-1]
        root, distance = TreeNode.from_newick(newick, labels)
        return Tree(root)

    def __str__(self):
        return self.to_newick()
    
    def __len__(self) -> int:
        return len(self._leaves)
    
    def __eq__(self, item):
        if not isinstance(item, Tree):
            return False
        return self._root == item._root
    
    def __hash__(self):
        return hash(self._root)

def _get_leaves(node: Optional[TreeNode], leaf_list: List[Optional[TreeNode]]):
    # cdef TreeNode child
    if node._index == -1:
        # Intermediate node -> Recursive calls
        for child in node._children:
            _get_leaves(child, leaf_list)
    else:
        # Node itself is leaf node -> add node -> terminate
        leaf_list.append(node)


def _get_leaf_count(node: TreeNode) -> int:
    count = 0
    if node._index == -1:
        # Intermediate node -> Recursive calls
        for child in node._children:
            count += _get_leaf_count(child)
        return count
    else:
        # Leaf node -> return count of itself = 1
        return 1


def _create_path_to_root(node: Optional[TreeNode]):
    """
    Create a list of nodes representing the path from this node to the
    specified node
    """
    path = []
    current_node = node
    while current_node is not None:
        path.append(current_node)
        current_node = current_node._parent
    return path



def as_binary(tree_or_node):
    if isinstance(tree_or_node, Tree):
        node, _ = _as_binary(tree_or_node.root)
        return Tree(node)
    elif isinstance(tree_or_node, TreeNode):
        node, _ = _as_binary(tree_or_node)
        return _as_binary(node)
    else:
        raise TypeError(
            f'Expected Tree or TreeNode, not {type(tree_or_node).__name__}'
        )

def _as_binary(node):
    """
    The actual logic wrapped by :func:`as_binary()`.
    
    Parameters
    ----------
    node : TreeNode
        The node to be converted.
    
    Returns
    -------
    binary_node: TreeNode
        The converted node.
    distance : float
        The distance of the converted node to its parent
    """
    children = node.children
    if children is None:
        # Leaf node
        return TreeNode(index=node.index), node.distance
    elif len(children) == 1:
        # Intermediate node with one child
        # -> Omit node and directly connect its child to its parent
        # The distances are added
        #
        #      |--            |--   
        #      |              |   
        # --|--|--   ->   ----|--  
        #      |              |   
        #      |--            |-- 
        #
        child, distance = _as_binary(node.children[0])
        if node.is_root():
            # Child is new root -> No distance to parent
            return child, None
        else:
            return child, node.distance + distance
    elif len(children) > 2:
        # Intermediate node with more than two childs
        # -> Create a new node having two childs:
        #    - One of the childs of the original node
        #    - The original node with one child less (distance = 0)
        # Repeat until all children are put into binary nodes
        #
        #   |--          |--
        #   |          --|  |--
        # --|--   ->     |--|
        #   |               |--
        #   |--
        #
        # The remaining children
        rem_children, distances = [list(tup) for tup in zip(
            *[_as_binary(child) for child in children]
        )]
        current_div_node = None
        while len(rem_children) > 0:
            if current_div_node is None:
                # The bottom-most node is created
                #-> Gets two of the remaining childs
                current_div_node = TreeNode(
                    rem_children[:2],
                    distances[:2]
                )
                # Pop the two utilized remaining childs from the list
                rem_children.pop(0)
                rem_children.pop(0)
                distances.pop(0)
                distances.pop(0)
            else:
                # A node is created that gets one remaining child
                # and the intermediate node from the last step
                current_div_node = TreeNode(
                    (current_div_node, rem_children[0]),
                    (0, distances[0.0]) 
                )
                # Pop the utilized remaining child from the list
                rem_children.pop(0)
                distances.pop(0)
        return current_div_node, node.distance
    else:
        # Intermediate node with exactly two childs
        # -> Keep node unchanged
        binary_children, distances = [list(tup) for tup in zip(
            *[_as_binary(child) for child in children]
        )]
        return TreeNode(binary_children, distances), node.distance



class TreeError(Static[Exception]):
    message: str
    def __init__(self, message: str): 
        super().__init__('TreeError', message)
        self.message = message
        
class InvalidFileError(Static[Exception]):
    message: str
    def __init__(self, message: str):
        super().__init__('InvalidFileError', message)
        self.message = message
        
