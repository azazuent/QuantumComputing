import random

import numpy as np
import scipy as sp

KET0 = np.array([[1],
                 [0]], dtype=complex)  # Вектор КЕТ0 - (1, 0)

KET1 = np.array([[0],
                 [1]], dtype=complex)  # КЕТ1

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
        [0, np.e**((2 * np.pi * 1j)/(2**k))]
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
    p = sum([isolated_qubit[i][i] for i in range(entangled.size)])

    return round(np.real(p), 5)


def measure_system(system):
    probs = get_p_of_all_states(system)
    measurement = [1 if random.random() < prob else 0 for prob in probs]
    return measurement


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
            i_operator = gen_c_operator(j-1, 0, R(j))
            i_operator = tensordot_arb(*prefix, i_operator, *postfix[j-1:])
            operator = operator @ i_operator

    return operator


def gen_swap(size):
    operator = gen_c_operator(0, size - 1, X)
    operator = operator @ gen_c_operator(size - 1, 0, X)
    operator = operator @ gen_c_operator(0, size - 1, X)

    return operator

def mod2x21():
    op1 = tensordot_arb(KET0 @ KET0.transpose(), I, I, I, I)
    op1 += tensordot_arb(KET1 @ KET1.transpose(), I, I, gen_swap(2))

    op2 = tensordot_arb(KET0 @ KET0.transpose(), I, I, I, I)
    op2 += tensordot_arb(KET1 @ KET1.transpose(), gen_swap(2), I, I)

    op3 = tensordot_arb(I, I, gen_c_operator(2, 0, X))
    op4 = gen_c_operator(4, 0, X)

    op5 = tensordot_arb(I, I, I, gen_swap(2))

    op6 = tensordot_arb(gen_swap(4), I)

    op7 = tensordot_arb(I, I, gen_swap(2), I)
    op8 = tensordot_arb(I, gen_swap(2), I, I)

    return op1 @ op2 @ op3 @ op4 @ op5 @ op6 @ op7 @ op8

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
