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


def entangle(*q_entities):
    res_q_entity = q_entities[0]
    for i in range(1, len(q_entities)):
        res_q_entity = sp.kron(res_q_entity, q_entities[i])
    return res_q_entity


def get_p_of_1(entangled, pos=0, relation=1):
    n_qubits = round(np.log2(entangled.size))

    operators = [I for _ in range(n_qubits)]
    target = KET0.copy() if relation == 0 else KET1.copy()
    operators[pos] = target @ target.transpose()

    projector = operators[0]
    for i in range(1, n_qubits):
        projector = sp.kron(projector, operators[i])

    isolated_qubit = projector @ entangled @ entangled.conjugate().transpose()
    p = sum([isolated_qubit[i][i] for i in range(entangled.size)])

    return round(np.real(p), 5)


def turn(angle: float, axe):
    if axe == 'x':
        return sp.linalg.expm(-1j * angle * X / 2)
    elif axe == 'y':
        return sp.linalg.expm(-1j * angle * Y / 2)
    elif axe == 'z':
        return sp.linalg.expm(-1j * angle * Z / 2)
    else:
        raise Exception("Invalid axe value")


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
        pr1 = get_p_of_1(self.state)
        sample = not np.random.random() <= pr1
        return bool(0 if sample else 1)

    def reset(self):
        self.state = KET0.copy()
