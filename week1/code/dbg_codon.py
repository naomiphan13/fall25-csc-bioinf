import copy

def reverse_complement(key):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rc = []
    for base in reversed(key):
        rc.append(complement.get(base, base)) 
    return ''.join(rc)

class Node:
    def __init__(self, kmer):
        self._children: set[str] = set()  
        self._count = 0
        self.kmer = kmer
        self.visited = False
        self.depth = 0
        
        self.max_depth_child = ""

    def add_child(self, kmer):
        if kmer:                 
            self._children.add(kmer)

    def increase(self):
        self._count += 1

    def reset(self):
        self.visited = False
        self.depth = 0
        self.max_depth_child = ""

    def get_count(self):
        return self._count

    def get_children(self):
        return list(self._children)

    def remove_children(self, target):
        self._children = self._children - target

class DBG:
    def __init__(self, k, data_list):
        self.k = k
        self.nodes = {}
        self.nodes["__seed__"] = Node("__seed__")
        del self.nodes["__seed__"]

        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list):
        assert len(data_list) > 0
        assert len(data_list[0]) > 0
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list):
        for data in data_list:
            for original in data:
                rc = reverse_complement(original)
                for i in range(len(original) - self.k - 1):
                    kmer1 = original[i: i + self.k]
                    kmer2 = original[i + 1: i + 1 + self.k]
                    self._add_arc(kmer1, kmer2)

                    rc_kmer1 = rc[i: i + self.k]
                    rc_kmer2 = rc[i + 1: i + 1 + self.k]
                    self._add_arc(rc_kmer1, rc_kmer2)

    def _add_node(self, kmer:str):
        if not kmer:
            return
        if kmer not in self.nodes:
            self.nodes[kmer] = Node(kmer)
        self.nodes[kmer].increase()

    def _add_arc(self, kmer1, kmer2):
        if not kmer1 or not kmer2:
            return
        self._add_node(kmer1)
        self._add_node(kmer2)
        self.nodes[kmer1].add_child(kmer2)

    def _get_count(self, child):
        return self.nodes[child].get_count()

    def _get_sorted_children(self, kmer):
        children = self.nodes[kmer].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, start_kmer: str) -> int:
        if start_kmer not in self.nodes:
            return 0
        if self.nodes[start_kmer].visited:
            return self.nodes[start_kmer].depth

        stack: list[tuple[str, int]] = [(start_kmer, 0)]
        self.nodes[start_kmer].visited = True

        while stack:
            kmer, idx = stack[-1]
            node = self.nodes[kmer]

            children = [c for c in self._get_sorted_children(kmer) if c in self.nodes]

            if idx < len(children):
                child = children[idx]
                stack[-1] = (kmer, idx + 1)
                child_node = self.nodes[child]
                if not child_node.visited:
                    child_node.visited = True
                    stack.append((child, 0))
            else:
                max_depth = 0
                max_child = ""
                for c in children:
                    d = self.nodes[c].depth
                    if d > max_depth:
                        max_depth, max_child = d, c
                node.depth = max_depth + 1
                node.max_depth_child = max_child
                stack.pop()

        return self.nodes[start_kmer].depth

    def _reset(self):
        for kmer in self.nodes:
            self.nodes[kmer].reset()

    def _get_longest_path(self):
        max_depth, start = 0, ""
        for kmer in self.nodes:
            depth = self._get_depth(kmer)
            if depth > max_depth:
                max_depth, start = depth, kmer

        if start is None:
            return []

        path = []
        curr = start

        while curr in self.nodes and curr != "":
            path.append(curr)
            curr = self.nodes[curr].max_depth_child
        return path

    def _delete_path(self, path):
        for kmer in path:
            if kmer in self.nodes:
                del self.nodes[kmer]
        path_set = set([k for k in path if k != ""])
        for kmer in list(self.nodes.keys()):
            self.nodes[kmer].remove_children(path_set)

    def _concat_path(self, path):
        if len(path) < 1:
            return None
        concat = copy.copy(path[0])
        for i in range(1, len(path)):
            concat += path[i][-1]
        return concat

    def get_longest_contig(self):
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig
