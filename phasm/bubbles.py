"""
Identify superbubbles in an assembly graph.
"""

import random
from typing import Iterator, Tuple

import networkx

from phasm.assembly_graph import AssemblyGraph


def partition_graph(g: AssemblyGraph) -> Iterator[Tuple[AssemblyGraph, bool]]:
    """This function partitions a directed graph into a set of directed acyclic
    subgraphs. It is partioned in such way that the set of super bubbles of `g`
    is the same as the union of the super bubble sets of all subgraphs returned
    by this function.

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

    start_nodes = (n for n in subgraph.nodes_iter()
                   if subgraph.in_degree(n) == 0)
    for n in start_nodes:
        if n == 'r_':
            continue
        subgraph.add_edge('r_', n)

    for u, v in g.out_edges_iter(singleton_nodes):
        if v not in subgraph:
            subgraph.add_edge(u, 're_')

    yield subgraph, True


def is_back_edge(t: networkx.DiGraph, e: Tuple[str, str]):
    u, v = e[:2]
    if not t.has_edge(u, v):
        if u in t and v in networkx.ancestors(t, u):
            return True

    return False


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
    non_artificial_nodes = (n for n in g if n not in {'r_', 're_'})
    for n in non_artificial_nodes:
        dag.add_nodes_from([n, str(n) + '*'])

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
