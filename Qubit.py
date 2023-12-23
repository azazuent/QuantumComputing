import random
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from itertools import product

KET0 = np.array([[1],
                 [0]], dtype=complex)  # Вектор КЕТ0 - (1, 0)

KET1 = np.array([[0],
                 [1]], dtype=complex)  # КЕТ1

KET_PLUS = np.array([[1],
                     [1]], dtype=complex) / np.sqrt(2)

KET_MINUS = np.array([[1],
                      [-1]], dtype=complex) / np.sqrt(2)

H = np.array([
    [1, 1],
    [1, -1]
], dtype=complex) / np.sqrt(2)  # Оператор Адамара, как указан

X = np.array([
    [0, 1],
    [1, 0]
], dtype=complex)  # Оператор Паули Х, как указан

Y = np.array([
    [0, -1j],
    [1j, 0]
], dtype=complex)  # Оператор Паули Y, как указан

Z = np.array([
    [1, 0],
    [0, -1]
], dtype=complex)  # Оператор Паули Z, как указан

Uf = np.array([
    [0, 1, 0, 0],
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
], dtype=complex)

I = np.array([
    [1, 0],
    [0, 1]
], dtype=complex)

CNOT = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]
], dtype=complex)

NCNOT = np.array([
    [0, 1, 0, 0],
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
], dtype=complex)

CZ = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, -1]], dtype=complex)

CH = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    np.array([0, 0, 1, 1]) / np.sqrt(2),
    np.array([0, 0, 1, -1]) / np.sqrt(2)], dtype=complex)


def R(k):
    return np.array([
        [1, 0],
        [0, np.e ** ((2 * np.pi * 1j) / (2 ** k))]
    ], dtype=complex)


def q_rng() -> bool:
    qubit = Qubit()
    qubit.h()
    return qubit.measure()


def tensordot_arb(*q_entities):
    res_q_entity = 1
    for i in range(len(q_entities)):
        res_q_entity = sp.kron(res_q_entity, q_entities[i])
    return res_q_entity


def get_p_of_state(entangled, pos=0, relation_to_0=False):
    n_qubits = round(np.log2(entangled.size))

    projector = 1
    for i in range(n_qubits):
        operator = I
        if i == pos:
            target = KET0.copy() if relation_to_0 else KET1.copy()
            operator = target @ target.transpose()
        projector = sp.kron(projector, operator)

    isolated_qubit = projector @ entangled @ entangled.conjugate().transpose()
    p = np.trace(isolated_qubit)

    return round(np.real(p), 5)


# def measure_system(system, left_bound=0, right_bound=None):
#     size = round(np.log2(system.size))
#
#     if right_bound is None:
#         right_bound = size - 1
#
#     probs = []
#
#     size_of_measured = right_bound - left_bound + 1
#
#     system_dagger = system.conjugate().transpose()
#     prefix = tensordot_arb(*[I for _ in range(left_bound)])
#     postfix = tensordot_arb(*[I for _ in range(size - right_bound - 1)])
#
#     proj0 = KET0 @ KET0.transpose()
#     proj1 = KET1 @ KET1.transpose()
#
#     for i in range(2 ** size_of_measured):
#         projector = tensordot_arb(*[proj0 if c == '0' else proj1
#                                     for c in format(i, f"0{size_of_measured}b")[::-1]])
#         projector = tensordot_arb(prefix, projector, postfix)
#         probs.append(np.trace(projector @ system @ system_dagger))
#     probs = np.real(probs)
#     plt.bar([i for i in range(2 ** size_of_measured)], probs)
#     plt.show()
#     return probs


def measure_system(system, left_bound=0, right_bound=None):
    size = round(np.log2(system.size))

    if right_bound is None:
        right_bound = size - 1

    size_left = left_bound
    size_of_measured = right_bound - left_bound + 1
    size_right = size - size_of_measured - left_bound

    probs = np.zeros(2 ** size_of_measured)

    combsL = list(map("".join, product("01", repeat=size_left)))
    combsM = list(map("".join, product("01", repeat=size_of_measured)))
    combsR = list(map("".join, product("01", repeat=size_right)))

    for cM in combsM:
        for cL in combsL:
            for cR in combsR:
                probs[int(cM[::-1], 2)] += abs(system[int(cL + cM + cR, 2)]) * abs(system[int(cL + cM + cR, 2)])
    probs = np.round(np.real(probs), 5)
    plt.bar([i for i in range(2 ** size_of_measured)], probs)
    plt.show()
    return np.real(probs)


def shoot(probs, shots=1024):
    stats = {i: 0 for i in range(len(probs))}
    positions = [i for i in range(len(probs))]
    for i in range(shots):
        index = random.choices(positions, probs)[0]
        stats[index] += 1
    stats = {key: stats[key] for key in stats if stats[key] != 0}
    return stats


def get_p_of_all_states(system, relation_to_0=False):
    n_qubits = round(np.log2(system.size))
    return [get_p_of_state(system, i, relation_to_0) for i in range(n_qubits)]


def turn(angle: float, axe):
    if axe == 'x':
        return sp.linalg.expm(-1j * angle * X / 2)
    elif axe == 'y':
        return sp.linalg.expm(-1j * angle * Y / 2)
    elif axe == 'z':
        return sp.linalg.expm(-1j * angle * Z / 2)
    else:
        raise Exception("Invalid axe value")


def gen_c_operator(control, target, operation, reverse=False):
    proj_basis = KET1
    if reverse:
        proj_basis = KET0

    size = abs(control - target) + 1

    c_operator = 1
    c_operator_complement = 1

    for i in range(size):

        complement_operator = I
        operator = I

        if i == control:
            operator = proj_basis @ proj_basis.transpose()
            complement_operator = (X @ proj_basis) @ (X @ proj_basis).transpose()

        if i == target:
            operator = operation

        c_operator = tensordot_arb(c_operator, operator)
        c_operator_complement = tensordot_arb(c_operator_complement, complement_operator)

    c_operator += c_operator_complement

    return c_operator


def gen_QFT(size):
    operator = tensordot_arb(*[I for _ in range(size)])
    for i in range(size):
        prefix = [I for _ in range(i)]
        postfix = [I for _ in range(size - i - 1)]
        operator = operator @ tensordot_arb(*prefix, H, *postfix)

        for j in range(2, size - i + 1):
            i_operator = gen_c_operator(j - 1, 0, R(j))
            i_operator = tensordot_arb(*prefix, i_operator, *postfix[j - 1:])
            operator = operator @ i_operator

    return operator


def gen_swap(size):
    operator = gen_c_operator(0, size - 1, X)
    operator = operator @ gen_c_operator(size - 1, 0, X)
    operator = operator @ gen_c_operator(0, size - 1, X)

    return operator


def n_mod_21(n: int):
    if n == 2:
        operator = tensordot_arb(KET0 @ KET0.transpose(), I, I, I, I) + \
                   tensordot_arb(KET1 @ KET1.transpose(), I, I, gen_swap(2))

        operator = operator @ \
                   (tensordot_arb(KET0 @ KET0.transpose(), I, I, I, I) +
                    tensordot_arb(KET1 @ KET1.transpose(), gen_swap(2), I, I))

        operator = operator @ tensordot_arb(I, I, gen_c_operator(2, 0, X)) @ gen_c_operator(4, 0, X)
        operator = operator @ tensordot_arb(I, I, I, gen_swap(2))
        operator = operator @ tensordot_arb(gen_swap(4), I)
        operator = operator @ tensordot_arb(I, I, gen_swap(2), I)
        operator = operator @ tensordot_arb(I, gen_swap(2), I, I)

        return operator
    elif n == 4:
        return gen_swap(5) @ tensordot_arb(I, I, gen_swap(3))
    elif n == 5:
        return tensordot_arb(X, I, X, I, X) @ tensordot_arb(I, I, gen_swap(3)) @ gen_swap(5)
    elif n == 8:
        return tensordot_arb(gen_swap(4), I)
    elif n == 13:
        return tensordot_arb(I, I, X, I, X)
    elif n == 16:
        return tensordot_arb(I, I, gen_swap(3)) @ gen_swap(5)
    elif n == 17:
        return tensordot_arb(X, I, X, I, X) @ tensordot_arb(I, I, gen_swap(3)) @ tensordot_arb(gen_swap(3), I, I)
    elif n == 20:
        return tensordot_arb(X, I, X, I, X)
    else:
        raise Exception(f"Modulus operator for {n} has not been implemented yet")


class Qubit:

    def __init__(self):  # Конструктор, инициализируем значением КЕТ0
        self.state = KET0.copy()

    def x(self):  # Для применения оператора Х, перемножаем
        self.state = X @ self.state  # текущее состояние с матрицей X

    def h(self):  # Для применения Адамара, перемножаем
        self.state = H @ self.state  # текущее состояние с матрицей H

    def z(self):
        self.state = Z @ self.state

    def ry(self, angle: float, axe):
        self.state = turn(angle, axe) @ self.state

    def measure(self):  # Функция измерения:
        pr1 = get_p_of_state(self.state)
        sample = np.random.random() <= pr1
        return bool(1 if sample else 0)

    def reset(self):
        self.state = KET0.copy()
