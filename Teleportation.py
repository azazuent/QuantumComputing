from Qubit import *
import numpy as np

qubit_1 = KET0
qubit_2 = KET0
qubit_3 = KET0

qubit_1 = turn(np.pi / 6, 'y') @ qubit_1
print(get_p_of_state(qubit_1))

entangled = entangle(qubit_1, qubit_2, qubit_3)

op1 = entangle(I, H, I)
op2 = entangle(I, CNOT)
op3 = entangle(CNOT, I)
op4 = entangle(H, I, I)
op5 = entangle(I, CNOT)
op6 = np.array([[1, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, -1, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, -1]], dtype=complex)

entangled = op1 @ entangled
entangled = op2 @ entangled
entangled = op3 @ entangled
entangled = op4 @ entangled
entangled = op5 @ entangled
entangled = op6 @ entangled

print(get_p_of_state(entangled, 2))
