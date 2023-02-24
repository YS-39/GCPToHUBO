# from numpy import *
import numpy as np
from scipy import special
from sympy import Symbol


def x(vertex1: int, vertex2: int):
    vertex = ((vertex1, vertex2),)

    return vertex


def x_hubo(vertex: int, digit: int):

    return 'x_{%d, %d}' % (vertex, digit)


def asc_binary_list(maxcolor: int):
    binary_list = []
    NB = int(np.ceil(np.log2(maxcolor)))
    for i in range(2 ** NB):
        binary_list.append(np.binary_repr(i, width=NB))

    return binary_list, NB


def desc_binary_list(maxcolor: int):
    binary_list = []
    NB = int(np.ceil(np.log2(maxcolor)))
    for i in reversed(range(2 ** NB)):
        binary_list.append(np.binary_repr(i, width=NB))

    return binary_list, NB


def state_rep(vertex: int, maxcolor: int, binary_desc=False):
    if binary_desc:
        binary_list, NB = desc_binary_list(maxcolor)
    else:
        binary_list, NB = asc_binary_list(maxcolor)
    sigma = []
    sigma_str = str()
    sigma_list = []
    for i in range(2 ** NB):
        for r in range(NB):
            x = Symbol('x_{%d,%d}' % (vertex, r + 1))
            if binary_list[i][r] == '0':
                sigma.append('(' + str(1 - int(binary_list[i][r]) + (2 * int(binary_list[i][r]) - 1) * x) + ')')
            else:
                sigma.append(str(1 - int(binary_list[i][r]) + (2 * int(binary_list[i][r]) - 1) * x))
    j = 0
    for k in range(len(sigma)):
        if j == NB - 1:
            sigma_str += sigma[k]
            sigma_list.append(sigma_str)
            sigma_str = str()
            j = 0
        else:
            sigma_str += sigma[k] + '*'
            j += 1

    return sigma_list


def GCPToHUBO(vertex_comb: tuple, num_vertex: int, maxcolor: int, binary_desc=False):
    vertex_list = []
    objfunc = str()
    costfunc = str()
    penaltyfunc = str()
    NB = int(np.ceil(np.log2(maxcolor)))
    for v in range(1, num_vertex + 1):
        vertex_list.append(state_rep(v, maxcolor, binary_desc))
    for c in range(len(vertex_comb)):
        if c == len(vertex_comb) - 1:
            for i in range(maxcolor):
                if i == maxcolor - 1:
                    costfunc += vertex_list[vertex_comb[c][0] - 1][i] + '*' + vertex_list[vertex_comb[c][1] - 1][i]
                else:
                    costfunc += vertex_list[vertex_comb[c][0] - 1][i] + '*' + vertex_list[vertex_comb[c][1] - 1][
                        i] + ' + '
        else:
            for i in range(maxcolor):
                costfunc += vertex_list[vertex_comb[c][0] - 1][i] + '*' + vertex_list[vertex_comb[c][1] - 1][i] + ' + '

    if not 2 ** NB == maxcolor:
        for u in range(num_vertex):
            if u == num_vertex - 1:
                for j in range(maxcolor, 2 ** NB):
                    if j == 2 ** NB - 1:
                        penaltyfunc += vertex_list[u][j]
                    else:
                        penaltyfunc += vertex_list[u][j] + ' + '
            else:
                for j in range(maxcolor, 2 ** NB):
                    penaltyfunc += vertex_list[u][j] + ' + '

    objfunc = costfunc + ' + ' + penaltyfunc

    return objfunc, costfunc, penaltyfunc


def GCP_Qubit_Count(num_edge: int, num_vertex: int, maxcolor: int):
    qubo_qubit = int(
        num_vertex * maxcolor + np.ceil(np.log2(num_edge * maxcolor + num_vertex * (maxcolor - 1) ** 2)) + 1)
    hubo_qubit = int(num_vertex * np.ceil(np.log2(maxcolor)) + np.ceil(np.log2(num_edge)) + 1)
    hubo_or_qubit = int(num_vertex * (np.ceil(np.log2(maxcolor)) + 1) + np.ceil(np.log2(num_edge)) + 1)

    return qubo_qubit, hubo_qubit, hubo_or_qubit


def GCP_CNOT_Count(num_edge: int, num_vertex: int, maxcolor: int):
    NB = int(np.ceil(np.log2(maxcolor)))
    hubo_asc_list = []
    hubo_or_list = []
    qubo = 2 * num_vertex * maxcolor * (
            np.ceil(np.log2(num_edge * maxcolor + num_vertex * (
                    maxcolor - 1) ** 2)) + 1) + 6 * 2 * num_edge * maxcolor + num_vertex * maxcolor * (
                   maxcolor - 1) * (1 / 2) * (
                   np.ceil(np.log2(num_edge * maxcolor + num_vertex * (maxcolor - 1) ** 2)) + 1)
    for k in range(1, 2 * NB + 1):
        if k == 1:
            G_kCR = 2 * (num_edge * special.binom(2 * NB, k) - special.binom(NB, k) * num_vertex * 5) * (
                    np.ceil(np.log2(num_edge)) + 1)
            hubo_asc_list.append(G_kCR)
        elif 1 < k <= NB:
            G_kCR = 6 * (k - 1) * (num_edge * special.binom(2 * NB, k) - special.binom(NB,
                                                                                       k) * num_vertex * 5) * (
                            np.ceil(np.log2(num_edge)) + 1)
            hubo_asc_list.append(G_kCR)
        else:
            G_kCR = 6 * (k - 1) * num_edge * special.binom(2 * NB, k) * (np.ceil(np.log2(num_edge)) + 1)
            hubo_asc_list.append(G_kCR)
    hubo_asc = sum(hubo_asc_list)

    hubo_pf = 6 * (NB - 1) * num_vertex * (2 ** NB - maxcolor) * (np.ceil(np.log2(num_edge)) + 1) + 6 * (
            2 * NB - 1) * num_edge * maxcolor * (np.ceil(np.log2(num_edge)) + 1)

    for k in range(1, NB + 2):
        if k == 1:
            G_kCR = 2 * ((2 ** k) * num_edge * special.binom(NB + 1, NB + 1 - k) - special.binom(NB + 1,
                                                                                                 k) * num_vertex * 5) * (
                            np.ceil(np.log2(num_edge)) + 1)
            hubo_or_list.append(G_kCR)
        else:
            G_kCR = 6 * (k - 1) * ((2 ** k) * num_edge * special.binom(NB + 1, NB + 1 - k) - special.binom(
                NB + 1, k) * num_vertex * 5) * (np.ceil(np.log2(num_edge)) + 1)
            hubo_or_list.append(G_kCR)
    hubo_or = sum(hubo_or_list)

    return int(qubo), int(hubo_asc), int(hubo_pf), int(hubo_or)


def GCP_TGate_Count(num_edge: int, num_vertex: int, maxcolor: int):
    NB = int(np.ceil(np.log2(maxcolor)))
    hubo_asc_list = []
    hubo_or_list = []
    qubo = 8 * 2 * num_edge * maxcolor + num_vertex * maxcolor * (maxcolor - 1) * (1 / 2) * (
            np.ceil(np.log2(num_edge * maxcolor + num_vertex * (maxcolor - 1) ** 2)) + 1)
    for k in range(2, 2 * NB + 1):
        if 2 <= k <= NB:
            G_kCR = 8 * (k - 1) * (num_edge * special.binom(2 * NB, k) - special.binom(NB,
                                                                                       k) * num_vertex * 5) * (
                            np.ceil(np.log2(num_edge)) + 1)
            hubo_asc_list.append(G_kCR)
        else:
            G_kCR = 8 * (k - 1) * num_edge * special.binom(2 * NB, k) * (np.ceil(np.log2(num_edge)) + 1)
            hubo_asc_list.append(G_kCR)
    hubo_asc = sum(hubo_asc_list)

    hubo_pf = 8 * (NB - 1) * num_vertex * (2 ** NB - maxcolor) * (np.ceil(np.log2(num_edge)) + 1) + 8 * (
            2 * NB - 1) * num_edge * maxcolor * (np.ceil(np.log2(num_edge)) + 1)

    for k in range(2, NB + 2):
        G_kCR = 8 * (k - 1) * ((2 ** k) * num_edge * special.binom(NB + 1, NB + 1 - k) - special.binom(
            NB + 1, k) * num_vertex * 5) * (np.ceil(np.log2(num_edge)) + 1)
        hubo_or_list.append(G_kCR)
    hubo_or = sum(hubo_or_list)

    return int(qubo), int(hubo_asc), int(hubo_pf), int(hubo_or)


def GCP_TotalTGate_Count(num_edge: int, num_vertex: int, maxcolor: int):
    NB = int(np.ceil(np.log2(maxcolor)))
    hubo_asc_list = []
    hubo_or_list = []
    qubo = 8 * 2 * num_edge * maxcolor + num_vertex * maxcolor * (maxcolor - 1) * (1 / 2) * (
            np.ceil(np.log2(num_edge * maxcolor + num_vertex * (maxcolor - 1) ** 2)) + 1)
    for k in range(2, 2 * NB + 1):
        if 2 <= k <= NB:
            G_kCR = 8 * (k - 1) * (num_edge * special.binom(2 * NB, k) - special.binom(NB,
                                                                                       k) * num_vertex * 5) * (
                            np.ceil(np.log2(num_edge)) + 1)
            hubo_asc_list.append(G_kCR)
        else:
            G_kCR = 8 * (k - 1) * num_edge * special.binom(2 * NB, k) * (np.ceil(np.log2(num_edge)) + 1)
            hubo_asc_list.append(G_kCR)
    hubo_asc = sum(hubo_asc_list)

    hubo_pf = 8 * (NB - 1) * num_vertex * (2 ** NB - maxcolor) * (np.ceil(np.log2(num_edge)) + 1) + 8 * (
            2 * NB - 1) * num_edge * maxcolor * (np.ceil(np.log2(num_edge)) + 1)

    for k in range(2, NB + 2):
        G_kCR = 8 * (k - 1) * ((2 ** k) * num_edge * special.binom(NB + 1, NB + 1 - k) - special.binom(
            NB + 1, k) * num_vertex * 5) * (np.ceil(np.log2(num_edge)) + 1)
        hubo_or_list.append(G_kCR)
    hubo_or = sum(hubo_or_list)

    total_qubo = int(qubo) + 2 * int(qubo) * np.sqrt(2 ** (num_vertex * maxcolor))
    total_hubo_asc = int(hubo_asc) + 2 * int(hubo_asc) * np.sqrt(2 ** (num_vertex * NB))
    total_hubo_pf = int(hubo_pf) + 2 * int(hubo_pf) * np.sqrt(2 ** (num_vertex * NB))
    total_hubo_or = int(hubo_or) + 2 * int(hubo_or) * np.sqrt(2 ** (num_vertex * (NB + 1)))

    return int(total_qubo), int(total_hubo_asc), int(total_hubo_pf), int(total_hubo_or)