import pytest  # noqa

from phasm.assembly_graph import (AssemblyGraph, remove_tips, node_path_edges,
                                  remove_transitive_edges, clean_graph)


def test_tip_removal():
    g = AssemblyGraph()
    g.add_edges_from(node_path_edges(['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7',
                                      'v8', 'v9', 'v10', 'v11', 'v12', 'v13',
                                      'v14', 'v15']))

    g.add_edges_from([
        ('v5', 'vt1'),
        ('vt1', 'vt2'),
        ('vt2', 'vt3')
    ])

    num_incoming, num_outgoing = remove_tips(g, 3)
    num_isolated_nodes = clean_graph(g)

    edges = set(g.edges())
    assert ('v5', 'vt1') not in edges
    assert ('vt1', 'vt2') not in edges
    assert ('vt2', 'vt3') not in edges

    assert ('v1', 'v2') in edges
    assert ('v8', 'v9') in edges
    assert ('v14', 'v15') in edges

    assert num_incoming == 0
    assert num_outgoing == 3
    assert num_isolated_nodes == 3

    g = AssemblyGraph()
    g.add_edges_from(node_path_edges(['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7',
                                      'v8', 'v9', 'v10', 'v11', 'v12', 'v13',
                                      'v14', 'v15']))

    g.add_edges_from([
        ('vt1', 'vt2'),
        ('vt2', 'vt3'),
        ('vt3', 'v6')
    ])

    num_incoming, num_outgoing = remove_tips(g, 3)
    num_isolated_nodes = clean_graph(g)

    edges = set(g.edges())
    assert ('vt1', 'vt2') not in edges
    assert ('vt2', 'vt3') not in edges
    assert ('vt3', 'v6') not in edges

    assert ('v1', 'v2') in edges
    assert ('v8', 'v9') in edges
    assert ('v14', 'v15') in edges

    assert num_incoming == 3
    assert num_outgoing == 0
    assert num_isolated_nodes == 3


def test_transitive_reduction():
    g = AssemblyGraph()
    g.add_edges_from(node_path_edges(['v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7',
                                      'v8', 'v9', 'v10', 'v11', 'v12']))

    for u, v, data in g.edges_iter(data=True):
        data['weight'] = 1

    g.add_edges_from([
        ('v2', 'v4', {'weight': 3}),
        ('v8', 'v10', {'weight': 3}),
        ('v8', 'v11', {'weight': 5})
    ])

    edges_to_remove = remove_transitive_edges(g)
    g.remove_edges_from(edges_to_remove)

    edges = set(g.edges())

    assert ('v2', 'v4') not in edges
    assert ('v8', 'v11') not in edges

    assert len(edges_to_remove) == 3
