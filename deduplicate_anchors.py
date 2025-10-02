"""
Clustering script for collapsing anchors that represent the same underlying sequence variation.

Script can be used on a list of anchors or more creatively on extendors (anchor-targets) or even targets themselves.
The output will give a list of sets that represent the sequences that have the same underlying sequence when tiled by
a given (maximum) shift distance.

Original author: Tavor Baharov
Edits and Maintenance provided by: Conor Messer
"""

import networkx as nx
import numpy as np
import scipy
import time
import sys
import os
import pandas as pd

import logging

logger = logging.getLogger(__name__)


def generate_adjacency_table(anchor_list, max_shift_dist):
    """
    Build boolean table that labels each pair of input sequences as identical within the given shift distance

    :param anchor_list: list of sequences
    :param max_shift_dist: int, maximum shift distance to compare sequences
    :return: sparse (n x n) matrix, where True denotes input i and input j are identical up to given shift distance
    """
    n = len(anchor_list)
    adj_table = scipy.sparse.csr_matrix((n, n), dtype=bool)
    for shift in range(1, max_shift_dist + 1):
        # todo only needed at max_shift_dist? Is it possible to align at shift=1 but not shift=5???
        seq_shift_l = np.asarray([a[shift:] for a in anchor_list])
        seq_shift_r = np.asarray([a[:-shift] for a in anchor_list])

        ### use hash table
        table_l = _generate_hash_table(seq_shift_l)
        table_r = _generate_hash_table(seq_shift_r)
        adj_table_update = _create_adj_table(table_l, table_r, n)

        adj_table = adj_table + adj_table_update

    logger.debug(f"adjacency sparse table size: {_sparse_size(adj_table) / 1000000} Mb")

    return adj_table


def _generate_hash_table(anchors):
    """
    Generates a dict that collects the indices of all equivalent sequences given

    :param anchors: list of str
    :return: dict of {str, anchor: list of int, indices} recording indices of all sequence pileups
    """
    hash_table = {}
    for i, anch in enumerate(anchors):
        if anch in hash_table.keys():
            prev_ids = hash_table[anch]
        else:
            prev_ids = []

        prev_ids.append(i)

        hash_table[anch] = prev_ids

    return hash_table


def _create_adj_table(table_1, table_2, n):
    """
    Generate adjacency table from two given hash tables

    TODO: need to do a little more research on most efficient sparse matrix representations.
    lil_matrix seems to be a good one to build incrementally, but should weigh memory vs. compute trade-offs of differing schemes

    :param table_1: dict of {anchors: list of anchor indices}, giving pileups of first (left shifted) sequences
    :param table_2: dict of {anchors: list of anchor indices}, giving pileups of second (right shifted) sequences
    :param n: total number of sequences
    :return: scipy.sparse matrix annotating overlapping sequences as determined by given dictionaries
    """
    # using csr_table tuple creation
    rows = []
    cols = []
    for anch in table_1.keys():
        for i in table_1[anch]:
            if anch in table_2:
                for j in table_2[anch]:
                    rows.append(i)
                    cols.append(j)
    data = [True] * len(rows)
    update_adj_table = scipy.sparse.csr_matrix(
        (data, (rows, cols)), shape=(n, n), dtype=bool
    )

    return update_adj_table


def generate_clusters(adjacency_table):
    """
    Generate clusters (connected components in a graph) from a given adjacency table

    :param adjacency_table: scipy or numpy boolean matrix, where [i,j]=True means that i and j are equivalent
    :return: list of sets of int, denoting the clusters of indices given by the adjacency table
    """
    ### from this similarity matrix, generate clusters
    # todo any speed-ups?
    adj_graph = nx.to_networkx_graph(adjacency_table)
    assemblies = list(nx.connected_components(adj_graph))

    logger.debug(f"Graph size: {_graph_size(adj_graph) / 1000000} Mb")
    logger.debug(f"Assemblies size: {sys.getsizeof(assemblies) / 1000000} Mb")

    ### order clusters by size
    assemblies.sort(key=len, reverse=True)
    return assemblies


def cluster_anchors(anchor_list, max_shift_dist=5):
    """
    Cluster given anchors by exact matches within given shift distance.

    :param anchor_list: list of sequences
    :param max_shift_dist: int, maximum shift distance to compare sequences
    :return: list of lists, where each inner list is a cluster of sequences
    """
    start = time.time()

    # get length of input sequences
    k = len(anchor_list[0])
    assert max_shift_dist <= k
    assert np.asarray(
        [len(anchor) == k for anchor in anchor_list]
    ).all()  # assert all sequences are same length
    logger.info(f"Number of input sequences: {len(anchor_list)}")
    logger.info(f"kmer size: {k}")
    logger.info(f"Max shift distance: {max_shift_dist}")
    logger.debug(f"Anchor list size: {sys.getsizeof(anchor_list)}")

    # generate similarity matrix,
    # i,j=1 indicates that read i is equivalent to j within shift distance max_shift_dist
    adj_table = generate_adjacency_table(anchor_list, max_shift_dist)
    adj_mat_time = time.time()
    logger.info(f"Time taken for adjacency matrix construction: {adj_mat_time - start}")
    logger.debug(f"Adjacency table size: {_sparse_size(adj_table)}")

    assemblies = generate_clusters(adj_table)
    networkx_time = time.time()
    logger.info(
        f"Time taken for networkx graph construction: {networkx_time - adj_mat_time}"
    )
    logger.info(f"Number of anchor overlaps: {adj_table.sum()}")

    # create output sets of sequences
    output = [{anchor_list[i] for i in list(cc)} for cc in assemblies]
    logger.info(f"Number of clusters: {len(output)}")

    return output


def _sparse_size(matrix):
    return (
        sys.getsizeof(matrix.data)
        + sys.getsizeof(matrix.indices)
        + sys.getsizeof(matrix.indptr)
    )


def _graph_size(graph):
    edge_mem = sum([sys.getsizeof(e) for e in graph.edges])
    node_mem = sum([sys.getsizeof(n) for n in graph.nodes])
    return edge_mem + node_mem


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to input file containing sequences",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to output file containing clusters",
    )
    parser.add_argument("-s", "--maxShiftDist", default=5, type=int)
    parser.add_argument("-d", "--log_debug", action="store_false")
    parser.add_argument("--temp_dir", type=str, default="Unused")

    args = parser.parse_args()

    level = logging.DEBUG if args.log_debug else logging.INFO
    logging.basicConfig(stream=sys.stdout, level=level)

    assert os.path.exists(args.input)
    anchors = list(pd.read_csv(args.input, header=None)[0])

    clusters = cluster_anchors(anchors, args.maxShiftDist)

    with open(args.output, "w") as f:
        for i, cluster in enumerate(clusters):
            for anchor in cluster:
                f.write(f"{i}\t{anchor}\n")
