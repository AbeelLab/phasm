import logging

import networkx
import pytest  # noqa

from phasm.bubbles import (partition_graph, graph_to_dag, SuperBubbleFinderDAG,
                           find_superbubbles, BubbleError)

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


@pytest.fixture
def test_graph():
    return get_test_graph()


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


@pytest.fixture
def test_graph2():
    return get_test_graph2()


def get_test_graph3():
    g = networkx.DiGraph()

    g.add_edges_from([
        (1, 2),
        (1, 3),
        (2, 5),
        (2, 4),
        (3, 9),
        (5, 6),
        (6, 4),
        (6, 7),
        (7, 8),
        (8, 10),
        (8, 11),
        (9, 8),
        (9, 11),
        (10, 12),
        (11, 12)
    ])

    return g


@pytest.fixture
def test_graph3():
    return get_test_graph3()


def get_test_graph4():
    g = networkx.DiGraph()

    g.add_edges_from([
        ('r_', 1),
        (1, 2),
        (1, 7),
        (2, 3),
        (2, 4),
        (2, 5),
        (2, 6),
        (3, 11),
        (4, 8),
        (5, 8),
        (6, 8),
        (7, 9),
        (7, 10),
        (8, 11),
        (9, 12),
        (10, 12),
        (11, 13),
        (12, 13),
        (13, 14),
        (13, 15),
        (14, 16),
        (15, 17),
        (15, 18),
        (16, 20),
        (17, 19),
        (17, 24),
        (18, 24),
        (19, 16),
        (19, 21),
        (20, 23),
        (21, 20),
        (21, 22),
        (22, 23),
        (23, 're_'),
        (24, 're_')
    ])

    return g


@pytest.fixture(scope="module")
def test_graph4():
    return get_test_graph4()


def test_partition(test_graph):
    logger = logging.getLogger('partition')
    g = test_graph

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
                ("r_", "1"),
                ("8", "re_")
            }

            assert set(partition.nodes()) == expected_nodes
            assert set(partition.edges()) == expected_edges
        else:
            raise ValueError("Unexpected partition, got more than 2")

        logger.debug("%s", partition.nodes())
        logger.debug("%s", partition.edges())

        assert networkx.is_directed_acyclic_graph(partition) == acyclic

        i += 1


def test_graph_to_dag(test_graph):
    logger = logging.getLogger('graph_to_dag')
    g = test_graph

    for partition, acyclic in partition_graph(g):
        if acyclic:
            continue

        dag, dfs_tree = graph_to_dag(partition)

        logger.debug("%s", dag.nodes())
        logger.debug("%s", dag.edges())
        logger.debug("DFS nodes: %s", dfs_tree.nodes())
        logger.debug("DFS edges: %s", dfs_tree.edges())

        assert networkx.is_directed_acyclic_graph(dag)


def test_find_superbubble_finder_dag(test_graph2, test_graph3, test_graph4):
    g = test_graph2
    superbubbles = SuperBubbleFinderDAG(g)

    superbubble_set = set(iter(superbubbles))
    assert superbubble_set == {
        ('8', '14'), ('3', '8'), ('11', '12'), ('5', '7'), ('1', '3')}

    g = test_graph3
    with pytest.raises(BubbleError):
        superbubble_finder = SuperBubbleFinderDAG(g)
        superbubble_set = set(iter(superbubble_finder))

    # Add artificial sink vertex
    sink_nodes = [n for n in g.nodes_iter() if g.out_degree(n) == 0]
    for n in sink_nodes:
        g.add_edge(n, '__sink__')

    superbubble_finder = SuperBubbleFinderDAG(g)
    superbubble_set = set(iter(superbubble_finder))

    assert superbubble_set == {(1, '__sink__'), (5, 6), (3, 9)}

    superbubbles = SuperBubbleFinderDAG(test_graph4)
    superbubble_set = set(iter(superbubbles))

    assert superbubble_set == {(1, 13), (2, 11), (7, 12)}


def test_find_superbubbles(test_graph, test_graph2, test_graph4):
    g = test_graph
    logger = logging.getLogger('find_superbubbles')

    superbubbles = set(iter(find_superbubbles(g)))
    logger.debug('Found superbubbles %s', superbubbles)
    assert superbubbles == {("3", "5"), ("7", "8")}

    superbubbles = set(iter(find_superbubbles(test_graph4)))
    assert superbubbles == {(1, 13), (2, 11), (7, 12)}

    superbubbles = set(iter(find_superbubbles(test_graph2,
                                              report_nested=False)))
    assert superbubbles == {('1', '3'), ('3', '8'), ('8', '14')}


if __name__ == '__main__':
    g1 = get_test_graph2()
    g2 = get_test_graph3()
    g3 = get_test_graph4()
    test_find_superbubble_finder_dag(g1, g2, g3)
