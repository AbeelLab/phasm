import logging

import networkx
import pytest  # noqa

from phasm.bubbles import partition_graph, graph_to_dag, SuperBubbleFinderDAG

logging.basicConfig(level=logging.DEBUG)


def get_test_graph():
    g = networkx.DiGraph()

    g.add_edges_from([
        ("1", "2"),
        ("1", "6"),
        ("2", "3"),
        ("2", "7"),
        ("3", "4"),
        ("3", "5"),
        ("4", "5"),
        ("5", "2"),
        ("5", "6"),
        ("6", "7"),
        ("7", "8")
    ])

    return g


def get_test_graph2():
    g = networkx.DiGraph()

    g.add_edges_from([
        ("1", "2"),
        ("1", "3"),
        ("2", "3"),
        ("3", "4"),
        ("3", "5"),
        ("3", "11"),
        ("4", "8"),
        ("5", "6"),
        ("5", "9"),
        ("6", "7"),
        ("6", "10"),
        ("7", "8"),
        ("8", "13"),
        ("8", "14"),
        ("9", "10"),
        ("10", "7"),
        ("11", "12"),
        ("12", "8"),
        ("13", "14"),
        ("13", "15"),
        ("15", "14")
    ])

    return g


def test_partition():
    logger = logging.getLogger('partition')
    g = get_test_graph()

    i = 0
    for partition, acyclic in partition_graph(g):
        if i == 0:
            expected_nodes = {'5', '4', '3', '2', 'r_', 're_'}
            expected_edges = {
                ("5", "2"),
                ("5", "re_"),
                ("4", "5"),
                ("2", "3"),
                ("2", "re_"),
                ("3", "4"),
                ("3", "5"),
                ("r_", "2")
            }

            assert set(partition.nodes()) == expected_nodes
            assert set(partition.edges()) == expected_edges
        elif i == 1:
            expected_nodes = {'8', '7', '6', '1', 'r_', 're_'}
            expected_edges = {
                ("7", "8"),
                ("6", "7"),
                ("1", "6"),
                ("1", "re_"),
                ("r_", "7"),
                ("r_", "6"),
                ("r_", "1")
            }

            assert set(partition.nodes()) == expected_nodes
            assert set(partition.edges()) == expected_edges
        else:
            raise ValueError("Unexpected partition, got more than 2")

        logger.debug("%s", partition.nodes())
        logger.debug("%s", partition.edges())

        assert networkx.is_directed_acyclic_graph(partition) == acyclic

        i += 1


def test_graph_to_dag():
    logger = logging.getLogger('graph_to_dag')
    g = get_test_graph()

    for partition, acyclic in partition_graph(g):
        if acyclic:
            continue

        dag, dfs_tree = graph_to_dag(partition)

        logger.debug("%s", dag.nodes())
        logger.debug("%s", dag.edges())
        logger.debug("DFS nodes: %s", dfs_tree.nodes())
        logger.debug("DFS edges: %s", dfs_tree.edges())

        assert networkx.is_directed_acyclic_graph(dag)


def test_find_superbubble_finder_dag():
    g = get_test_graph2()
    superbubbles = SuperBubbleFinderDAG(g)

    superbubble_set = set(iter(superbubbles))
    assert superbubble_set == {
        ('8', '14'), ('3', '8'), ('11', '12'), ('5', '7'), ('1', '3')}


if __name__ == '__main__':
    test_find_superbubble_finder_dag()
