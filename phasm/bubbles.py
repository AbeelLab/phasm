"""
Identify superbubbles in an assembly graph.
"""

import random
from enum import Enum
from typing import Iterator, Tuple, Hashable, NamedTuple, Set
from collections import deque

import networkx

from phasm.assembly_graph import AssemblyGraph
from phasm.rmq import RangeMinimumQuery, RangeMaximumQuery

Node = Hashable
Bubble = Tuple[Node, Node]


class CandidateType(Enum):
    ENTRANCE = 1
    EXIT = 2


Candidate = NamedTuple('Candidate', [
    ('v', Node), ('type', CandidateType)])


def partition_graph(g: AssemblyGraph) -> Iterator[Tuple[AssemblyGraph, bool]]:
    """This function partitions a directed graph into a set of subgraphs. It
    is partioned in such way that the set of super bubbles of `g` is the same
    as the union of the super bubble sets of all subgraphs returned by this
    function.

    This function yields each partitioned subgraph, together with a flag if
    it is acyclic or not.
    """

    singleton_nodes = []

    for connected_component in networkx.strongly_connected_components(g):
        if len(connected_component) == 1:
            singleton_nodes.extend(connected_component)
            continue

        subgraph = g.subgraph(connected_component)
        for u, v in g.in_edges_iter(connected_component):
            if u not in subgraph:
                subgraph.add_edge('r_', v)

        for u, v in g.out_edges_iter(connected_component):
            if v not in subgraph:
                subgraph.add_edge(u, "re_")

        yield subgraph, False

    # Build subgraph with only singleton strongly connected components
    subgraph = g.subgraph(singleton_nodes)
    for u, v in g.in_edges_iter(singleton_nodes):
        if u not in subgraph:
            subgraph.add_edge('r_', v)

    start_nodes = [n for n in subgraph.nodes_iter()
                   if subgraph.in_degree(n) == 0]
    for n in start_nodes:
        if n == 'r_':
            continue
        subgraph.add_edge('r_', n)

    for u, v in g.out_edges(singleton_nodes):
        if v not in subgraph:
            subgraph.add_edge(u, 're_')

    yield subgraph, True


def is_back_edge(t: networkx.DiGraph, e: Tuple[Node, Node]):
    u, v = e[:2]
    if not t.has_edge(u, v):
        if u in t and v in networkx.ancestors(t, u):
            return True

    return False


def _is_duplicated(g, node):
    return node in g and 'original' in g.node[node]


def _duplicated_id(node):
    return str(node) + '*'


def graph_to_dag(g: AssemblyGraph) -> Tuple[AssemblyGraph, networkx.DiGraph]:
    """Converts a general strongly connected directed graph to a
    directed acyclic graph, by duplicating certain nodes and adding edges in
    such way that no cycle is created.

    This DAG can be used to identify superbubbles in the original graph, for
    details please refer to the paper by Sung et al. [SUNG2015]_.

    This function returns a tuple with the newly created DAG and the DFS tree
    which can be used to identify back edges.

    .. [SUNG2015] Wing-Kin Sung, Sadakane, K., Shibuya, T., Belorkar, A.,
                  & Pyrogova, I. (2015). An O(m log m)-Time Algorithm for
                  Detecting Superbubbles. IEEE/ACM Transactions on
                  Computational Biology and Bioinformatics, 12(4), 770–777.
                  http://doi.org/10.1109/TCBB.2014.2385696
    """

    dag = networkx.DiGraph()

    # Start with creating the DFS tree, which can be used to identify "back
    # edges"
    start_node = 'r_'
    if start_node not in g:
        # This random chosen root should work because the graph is expected to
        # be strongly connected
        start_node = random.choice([n for n in g.nodes_iter() if n != 're_'])

    dfs_tree = networkx.dfs_tree(g, start_node)

    # Add nodes, create two nodes for each node in the original graph
    # Duplicated nodes are marked with an *
    non_artificial_nodes = (n for n in g if n not in {'r_', 're_'})
    for n in non_artificial_nodes:
        dag.add_node(n)
        dag.add_node(_duplicated_id(n), original=n)

    dag.add_nodes_from(['r_', 're_'])

    # Add edges from original graph in such way that it does not create cycles,
    # here we can use the duplicated nodes created above.
    for u, v in g.edges_iter():
        if u == 'r_':
            dag.add_edge(u, v)
        elif v == 're_':
            dag.add_edge(str(u) + '*', v)
        else:
            if is_back_edge(dfs_tree, (u, v)):
                dag.add_edge(u, str(v) + '*')
            else:
                dag.add_edge(u, v)
                dag.add_edge(str(u) + '*', str(v) + '*')

    if 'r_' not in g:
        start_nodes = (n for n in dag.nodes_iter() if dag.in_degree(n) == 0)
        for n in start_nodes:
            if n == 'r_' or n == 're_':
                continue

            dag.add_edge('r_', n)

    if 're_' not in g:
        end_nodes = (n for n in dag.nodes_iter() if dag.out_degree(n) == 0)
        for n in end_nodes:
            if n == 'r_' or n == 're_':
                continue

            dag.add_edge(n, 're_')

    return dag, dfs_tree


class SuperBubbleFinderDAG:
    """Identifies and reports all superbubbles in a directed acyclic graph.

    Based on the O(|V| + |E|)-algorithm by Brankovic et al. [BRANKOVIC2015]_.

    .. [BRANKOVIC2015] Brankovic, L., Iliopoulos, C. S., Kundu, R.,
                       Mohamed, M., Pissis, S. P., & Vayani, F. (2015).
                       Linear-Time Superbubble Identification Algorithm for
                       Genome Assembly. arXiv. Retrieved from
                       http://arxiv.org/abs/1505.04019
    """

    def __init__(self, g: AssemblyGraph):
        self.g = g
        self.topo_sorted = networkx.topological_sort(g)
        self.ordering = {
            v: i for i, v in enumerate(self.topo_sorted)
        }

        # Prepare arrays for minimum/maximum query in a range
        # (Range minimum/maximum query problem)
        self.out_parent = [
            min(self.ordering[p]
                for p in g.predecessors_iter(self.topo_sorted[i]))
            if g.in_degree(self.topo_sorted[i]) else 2**32 - 1
            for i in range(len(self.topo_sorted))
        ]
        self.out_parent_rmq = RangeMinimumQuery(self.out_parent)

        self.out_child = [
            max(self.ordering[p] for p in g.neighbors_iter(
                self.topo_sorted[i]))
            if g.out_degree(self.topo_sorted[i]) > 0 else -1
            for i in range(len(self.topo_sorted))
        ]
        self.out_child_rmq = RangeMaximumQuery(self.out_child)

        self.previous_entrance = {}  # type: Dict[Node, Node]
        self.alternative_entrance = {}  # type: Dict[Node, Node]

        # Build candidates
        self.candidates = []  # type: List[Candidate]
        self.v_to_candidate_index = {}  # type: Dict[Node, int]

        # Keep track of the last seen entrance in the candidates list, used to
        # fill the previous_entrance dictionary
        prev_entr = None

        for v in self.topo_sorted:
            self.alternative_entrance[v] = None
            self.previous_entrance[v] = prev_entr

            if self._is_exit(v):
                self.v_to_candidate_index[v] = len(self.candidates)
                self.candidates.append(Candidate(v, CandidateType.EXIT))

            if self._is_entrance(v):
                # Only keep entrance index if both exit and entrance
                self.v_to_candidate_index[v] = len(self.candidates)
                self.candidates.append(Candidate(v, CandidateType.ENTRANCE))
                prev_entr = v

    def _is_entrance(self, v: Node):
        for n in self.g.neighbors_iter(v):
            if self.g.in_degree(n) == 1:
                return True

        return False

    def _is_exit(self, v: Node):
        for n in self.g.predecessors_iter(v):
            if self.g.out_degree(n) == 1:
                return True

        return False

    def __iter__(self):
        while self.candidates:
            if self.candidates[-1].type == CandidateType.ENTRANCE:
                self.candidates.pop()
            else:
                yield from self._check_candidates(0, len(self.candidates)-1)

    def _check_candidates(self, start_index: int, end_index: int):
        start = self.candidates[start_index].v
        exit = self.candidates[end_index].v

        if (not start or not exit or self.ordering[start] >=
                self.ordering[exit]):
            self.candidates.pop()
            return

        possible_start = self.previous_entrance[exit]
        valid = None
        while self.ordering[possible_start] >= self.ordering[start]:
            valid = self._validate_superbubble(possible_start, exit)
            if not valid:
                break

            if (valid == possible_start or valid ==
                    self.alternative_entrance[possible_start]):
                break

            self.alternative_entrance[possible_start] = valid
            possible_start = valid

        self.candidates.pop()
        if possible_start == valid:
            # At this point we have found a correct superbubble entrance-exit
            # pair
            yield (possible_start, exit)

            # Check for nested superbubbles
            while self.candidates[-1].v != possible_start:
                if self.candidates[-1].type == CandidateType.EXIT:
                    yield from self._check_candidates(
                        self.v_to_candidate_index[possible_start]+1,
                        len(self.candidates)-1)
                else:
                    self.candidates.pop()

    def _validate_superbubble(self, start_v: Node, exit_v: Node):
        start = self.ordering[start_v]
        exit = self.ordering[exit_v]

        outchild_idx = self.out_child_rmq.query(start, exit)
        outchild_ord = self.out_child[outchild_idx]

        outparent_idx = self.out_parent_rmq.query(start+1, exit+1)
        outparent_ord = self.out_parent[outparent_idx]

        if outchild_ord != exit:
            return None

        if outparent_ord == start:
            return start_v

        outparent_v = self.topo_sorted[outparent_ord]
        if outparent_v in self.v_to_candidate_index:
            outparent_can = self.candidates[
                self.v_to_candidate_index[outparent_v]]
            if outparent_can.type == CandidateType.ENTRANCE:
                return outparent_v

        return self.previous_entrance[outparent_v]


def find_superbubbles(g: AssemblyGraph) -> Iterator[Bubble]:
    """Find superbubbles in a general directed graph.

    This algorithm first partitions the graph into several subgraphs,
    converts each subgraph to a DAG if necessary, and then finds superbubbles
    in this constructed DAG. The partitioning and DAG construction happens in
    such a way that superbubbles in the original graph can be easily found.

    For more details please refer to the paper by Sung et al. [SUNG2015]_.

    .. [SUNG2015] Wing-Kin Sung, Sadakane, K., Shibuya, T., Belorkar, A., &
                  Pyrogova, I. (2015). An O(m log m)-Time Algorithm for
                  Detecting Superbubbles. IEEE/ACM Transactions on
                  Computational Biology and Bioinformatics, 12(4), 770–777.
                  http://doi.org/10.1109/TCBB.2014.2385696
    """

    for partition, acyclic in partition_graph(g):
        if acyclic:
            superbubble_finder = SuperBubbleFinderDAG(partition)
            yield from (b for b in superbubble_finder if 'r_' not in b
                        and 're_' not in b)
        else:
            dag, dfs_tree = graph_to_dag(partition)
            superbubble_finder = SuperBubbleFinderDAG(dag)

            superbubbles = set(iter(superbubble_finder))

            # The graph_to_dag algorithm duplicates nodes, so we have to filter
            # some reported superbubbles out
            for s, t in superbubbles:
                # Ignore any superbubble involving artificial roots and sinks
                if s == 'r_' or t == 're_':
                    continue

                if _is_duplicated(g, s) and _is_duplicated(g, t):
                    orig_s = g.node[s]['original']
                    orig_t = g.node[t]['original']

                    if (orig_s, orig_t) in superbubbles:
                        # (s, t) and (s*, t*) are both reported as
                        # superbubbles, so this is a valid bubble
                        yield orig_s, orig_t
                elif not _is_duplicated(g, s) and not _is_duplicated(g, t):
                    if (_duplicated_id(s), _duplicated_id(t)) in superbubbles:
                        # (s, t) and (s*, t*) are both reported as
                        # superbubbles, so this is a valid bubble
                        yield s, t
                elif not _is_duplicated(g, s) and _is_duplicated(g, t):
                    orig_t = g.node[t]['original']

                    # A (s, t*) superbubble is always a valid bubble
                    yield s, orig_t


def superbubble_nodes(g: AssemblyGraph, source: Node,
                      sink: Node) -> Set[Node]:
    """Find all nodes inside a superbubble."""

    queue = deque([source])
    visited = set([source, sink])

    while queue:
        current = queue.popleft()

        for neighbour in g.neighbors_iter(current):
            if neighbour not in visited:
                queue.append(neighbour)
                visited.add(neighbour)

    return visited
