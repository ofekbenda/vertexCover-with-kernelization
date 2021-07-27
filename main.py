import json
import random
import timeit

import networkx as nx
import numpy as np
from scipy.optimize import linprog
from matplotlib import pyplot as plt


# def kernelization_algorithm_1(G, k):
#     graph_dict = nx.to_dict_of_lists(G, None)
#
#     # step 1: If there exists an isolated vertex, remove it
#     for key in list(graph_dict.keys()):
#         print(graph_dict.keys())
#         if not graph_dict[key]:
#             G.remove_node(key)
#             print(f"G: {G.nodes}")
#             kernelization_algorithm_1(G, k)
#
#     # step 2: If there exists a vertex with at least k+1 neighbors, remove it and decrement k by 1.
#     for key in list(graph_dict.keys()):
#         if len(graph_dict[key]) >= k + 1:
#             G.remove_node(key)
#             k -= 1
#             kernelization_algorithm_1(G, k-1)
#
#     # step 3: If G contains more than k^2 edges, return NO
#     if k < 0 or ((len(G.edges) > k**2) or (len(G.nodes) > 2 * k**2)):
#         return False
#     return True if len(graph_dict) == 0 else G

photo_arr = []

def kernelization_algorithm(G, k):
    changed = True
    while(changed):
        changed = False
        to_remove = []

        # step 1: If there exists an isolated vertex, remove it
        for node in G.nodes:
            if G.degree(node) == 0:
                to_remove.append(node)
                changed = True

        for node in to_remove:
            G.remove_node(node)

        to_remove = []
        # step 2: If there exists a vertex with at least k+1 neighbors, remove it and decrement k by 1.
        for node in G.nodes:
            if G.degree(node) >= k + 1:
                to_remove.append(node)
                changed = True

        for node in to_remove:
            G.remove_node(node)
            k -= 1

    # step 3: If G contains more than k^2 edges, return NO
    if k < 0 or ((len(G.edges) > k**2) or (len(G.nodes) > 2 * k**2)):
        return None
    return G if len(G.edges) == 0 else G


def kernelization_algorithm_opt(G, k):
    changed = True
    while(changed):
        changed = False
        to_remove = []
        # step 1: If there exists an isolated vertex, remove it +
        # step 2: If there exists a vertex with at least k+1 neighbors, remove it and decrement k by 1.
        reduce_k = 0
        for node in G.nodes:
            if G.degree(node) == 0:
                to_remove.append(node)
                changed = True

            if G.degree(node) >= k + 1:
                to_remove.append(node)
                reduce_k += 1
                changed = True

        for node in to_remove:
            G.remove_node(node)

        k -= reduce_k

    # step 3: If G contains more than k^2 edges, return NO
    if k < 0 or ((len(G.edges) > k**2) or (len(G.nodes) > 2 * k**2)):
        return None
    return G
    # return True if len(G.nodes) == 0 else G


def recursive_algorithem(G, k):
    if G is None:
        return G

    if k <= 0 and len(G.edges) > 0:
        return None

    if len(G.edges) == 0 and k < 0:
        return None

    if len(G.edges) == 0 and k >= 0:
        return []

    u, v = list(G.edges)[0]
    nodes_without_u = list(G.nodes)
    nodes_without_v = list(G.nodes)

    nodes_without_u.remove(u)
    nodes_without_v.remove(v)

    S1 = recursive_algorithem(G.subgraph(nodes_without_u), k - 1)
    S2 = recursive_algorithem(G.subgraph(nodes_without_v), k - 1)

    if S1 is not None and (S2 is None or (len(S1) <= len(S2))):
        S1.append(u)
        return S1
    if S2 is not None and (S1 is None or (len(S2) < len(S1))):
        S2.append(v)
        return S2

    return None

# def prev_rec(G,k):
#     if G is None:
#         return False
#
#     if k == 0 and len(G.edges) > 0:
#         return False
#
#     if len(G.edges) == 0:
#         return True
#
#     u, v = list(G.edges)[0]
#     nodes_without_u = list(G.nodes)
#     nodes_without_v = list(G.nodes)
#
#     nodes_without_u.remove(u)
#     nodes_without_v.remove(v)
#
#     if recursive_algorithem(G.subgraph(nodes_without_u), k - 1):
#         return True
#     if recursive_algorithem(G.subgraph(nodes_without_v), k - 1):
#         return True
#
#     return False


def recursive_algorithem_opt(G, k):
    if G is None:
        return None

    if k <= 0 and len(G.edges) > 0:
        return None

    if len(G.edges) == 0 and k < 0:
        return None

    if len(G.edges) == 0 and 0 <= k:
        return []

    # step 1: find a vertex v ∈ V (G) of maximum degree in G
    v = None
    max_deg = -1
    for n in G.nodes:
        if G.degree(n) > max_deg:
            max_deg = G.degree(n)
            v = n

    # If deg(v) = 1, then G contains only isolated vertices or an signal edges. The instance has a trivial solution.
    if G.degree(v) == 1:
        return [v]

    # If deg(v) != 1, recursive call:
    # add n to S - vertex cover group or add all n's neighbors
    v_neighbors = list(G.neighbors(v))


    nodes_without_v = list(G.nodes)
    nodes_without_v_neighbors = list(G.nodes)

    nodes_without_v.remove(v)
    for neighbor in v_neighbors:
        nodes_without_v_neighbors.remove(neighbor)

    S1 = recursive_algorithem_opt(G.subgraph(nodes_without_v), k-1)
    S2 = recursive_algorithem_opt(G.subgraph(nodes_without_v_neighbors), k - len(v_neighbors))

    if S1 is not None and (S2 is None or (len(S1) <= len(S2))):
        S1.append(v)
        return S1
    if S2 is not None and (S1 is None or (len(S2) < len(S1))):
        S2 += v_neighbors
        return S2
    else:
        return None


# VC kernelization algorithm based on LP
def LPVC(G):
    n = len(G.nodes)
    m = len(G.edges)
    A = []
    nodes_list = np.array(list(G.nodes))

    # Set the coefficients of the linear objective function
    c = [1] * n

    # Set the inequality constraints matrix
    # for each uv in E(G) -> x_u + x_v >= 1
    for e in G.edges:
        node_u = np.where(nodes_list == e[0])[0][0]
        node_v = np.where(nodes_list == e[1])[0][0]
        curr_array = [0] * n
        curr_array[node_u] = -1
        curr_array[node_v] = -1
        A = A + [curr_array]

    # for each v in V(G) -> 0 <= x_v <= 1
    bounds = [(0, 1)] * n

    # Set the inequality constraints vec
    b = [-1] * m

    return linprog(c, A, b, None, None, bounds, "simplex")


# Reduction 5 - Invoke the algorithm of lemma 3.7
def LP_VC5(G, k):
    min = k
    node_index = 0
    min_lp = []
    nodes_list = np.array(list(G.nodes))
    sol = [0] * len(nodes_list)

    for v in G.nodes:
        G1 = nx.Graph()
        G1.add_nodes_from(G.nodes)
        G1.add_edges_from(G.edges)
        G1.remove_node(v)
        lp = LPVC(G1)
        if round(lp.fun, 2) <= k-1 and round(lp.fun, 2) < min:
            node_index = np.where(nodes_list == v)[0][0]
            min = round(lp.fun, 2)
            min_lp = lp
    for i in range(0, len(nodes_list)):
        if i == node_index:
            sol[i] = 1.0
        elif i < node_index:
            sol[i] = min_lp.x[i]
        else:
            sol[i] = min_lp.x[i-1]
    return sol


# Reduction 4
def LP_red4(G, lp_x):

    nodes_list = np.array(list(G.nodes))
    v0 = []
    v_half = []
    v1 = []

    # Find the VC from lp using reduction 4
    for v in G.nodes:
        node_v = np.where(nodes_list == v)[0][0]
        if lp_x[node_v] < 0.5:
            v0.append(v)
        if lp_x[node_v] == 0.5:
            v_half.append(v)
        if lp_x[node_v] > 0.5:
            v1.append(v)
    return [v0, v_half, v1]


# Theorem 3.8
def LP_38(G):
    if len(G.nodes) > 0:
        v = np.array(list(G.nodes))[0]

        G1 = nx.Graph()
        G2 = nx.Graph()
        G1.add_nodes_from(G.nodes)
        G2.add_nodes_from(G.nodes)
        G1.add_edges_from(G.edges)
        G2.add_edges_from(G.edges)

        # pick v to the vertex cover
        G1.remove_node(v)
        res_v = [v] + LP_parm(G1, k - 1)

        # pick N[v] to the vertex cover
        G2.remove_nodes_from(list(G.neighbors(v)))
        res_Nv = [v] + LP_parm(G2, k - len(list(G.neighbors(v))))

        if len(res_v) < len(res_Nv):
            return res_v
        else:
            return res_Nv


# VC parameterized algorithm above LP
def LP_parm(G, k):
    # if G is None:
    #     return None

    n = len(G.nodes)

    if len(G.edges) == 0:
        return []

    lp = LPVC(G)

    if round(lp.fun, 2) > k:
        return None

    # step 1 - reduction 4
    res_red4 = LP_red4(G, lp.x)
    v0 = res_red4[0]
    v_half = res_red4[1]
    v1 = res_red4[2]

    # step 2: if x is all-1/2-solution - reduction VC.5.
    if len(v_half) == n:
        lp_x = LP_VC5(G, k)
        res_red4 = LP_red4(G, lp_x)
        v0 = res_red4[0]
        v1 = res_red4[2]

    # step 3: reduction 4
    for v in v0:
        G.remove_node(v)
    for v in v1:
        G.remove_node(v)

    # step 5: if x was all-1/2-solution - apply 3.8
    if len(v_half) == n and len(G.edges) > 0:
        return v1 + LP_38(G)

    return v1 + LP_parm(G, k - len(v1))


if __name__ =='__main__':
    # spacial graphs
    # G = nx.cycle_graph(20)
    # k = 10

    # G = nx.complete_bipartite_graph(10,10)
    # k = 10

    G = nx.balanced_tree(2, 4)
    k = 10

    nx.draw_networkx(G, with_labels = True)
    plt.show()
    # Recursive
    start = timeit.default_timer()
    res1 = recursive_algorithem(G, k)
    end = timeit.default_timer()
    print("recursive without kernelization:")
    print(f"runtime: {end - start}")
    print(f"result: {res1}")

    # Recursive with opt
    start = timeit.default_timer()
    res2 = recursive_algorithem_opt(G, k)
    end = timeit.default_timer()
    print("recursive opt without kernelization:")
    print(f"runtime: {end - start}")
    print(f"result: {res2}")

    # LP
    start = timeit.default_timer()
    res3 = LP_parm(G, k)
    end = timeit.default_timer()
    print("LP without kernelization:")
    print(f"runtime: {end - start}")
    print(f"result: {res3}")
    # --------------------------------------------------

    # g = nx.DiGraph()
    # g.add_nodes_from([int(key) for key in data.keys()])
    #
    # for k, v in data.items():
    #     g.add_edges_from(([(int(k), t) for t in v]))
    #
    # print(g.adj)
    #
    # k = int(len(g.edges)/2)
    # start = timeit.default_timer()
    # print(f"lp: {LP_parm(g, k)}")
    # end = timeit.default_timer()
    # print(f"runtime: {end-start}")


    # ------------------------  TEST REC VS REC_OPT VS LP WITHOUT KERNELIZATION ------------------
    #  runtime_rec_without_kernel = []
    # runtime_rec_opt_without_kernel = []
    # runtime_LP_without_kernel = []

    # without kernelization: recursive vs. LP
    # print("WITHOUT KERNELIZATION:")
    # k = 10
    # for n in range(10, 2000, 50):
    #     m = random.randint(int(n/4), n)
    #     G = nx.gnm_random_graph(n, m)

        # k = int(m/2)
        # print(" **** ")
        # print("number of nodes: ", G.number_of_nodes())
        # print("number of edges: ", G.number_of_edges())
        # print(" **** ")
        # G1 = kernelization_algorithm(G, k)
        # Recursive
        # start = timeit.default_timer()
        # recursive_algorithem(G, k)
        # end = timeit.default_timer()
        # runtime_rec_without_kernel.append(end-start)

        # Recursive with opt
        # start = timeit.default_timer()
        # recursive_algorithem_opt(G, k)
        # end = timeit.default_timer()
        # runtime_rec_opt_without_kernel.append(end - start)

        # # LP
        # start = timeit.default_timer()
        # LP_parm(G, k)
        # end = timeit.default_timer()
        # runtime_LP_without_kernel.append(end - start)

    # print(runtime_rec_without_kernel)
    # print(runtime_rec_opt_without_kernel)
    # print(runtime_LP_without_kernel)

    # plt.plot(range(10, 2000, 50), runtime_rec_without_kernel, label="Recursive algorithm")
    # plt.plot(range(10, 2000, 50), runtime_rec_opt_without_kernel, label="Recursive improvement algorithm")
    # plt.plot(range(10, 2000, 50), runtime_LP_without_kernel, label="LP algorithm")
    # plt.title("Without Kernelization")
    # plt.xlabel("Number of nodes")
    # plt.ylabel("Runtime")
    # plt.legend()
    # plt.show()


    #----------------------- TEST KERNELIZATION -----------------------------
    # G = nx.gnm_random_graph(100, 70)
    # k = 20
    # Recursive with opt without kernelization
    # start = timeit.default_timer()
    # res = recursive_algorithem_opt(G, k)
    # end = timeit.default_timer()
    # print("Recursive with opt without kernelization:")
    # print(f"runtime: {end - start}")
    # print(f"result: {res}")
    #
    # # Recursive with opt with kernelization
    # start = timeit.default_timer()
    # G1 = kernelization_algorithm(G, k)
    # res1 = recursive_algorithem_opt(G1, k)
    # end = timeit.default_timer()
    # print("Recursive with opt with kernelization:")
    # print(f"runtime: {end - start}")
    # print(f"result: {res1}")

    # lp without kernelization
    # start = timeit.default_timer()
    # res = recursive_algorithem_opt(G, k)
    # end = timeit.default_timer()
    # print("LP without kernelization:")
    # print(f"runtime: {end - start}")
    # print(f"result: {res}")
    #
    # # lp with kernelization
    # start = timeit.default_timer()
    # G1 = kernelization_algorithm(G, k)
    # res1 = recursive_algorithem_opt(G1, k)
    # end = timeit.default_timer()
    # print("LP with kernelization:")
    # print(f"runtime: {end - start}")
    # print(f"result: {res1}")



