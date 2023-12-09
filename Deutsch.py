from Qubit import *
import numpy as np


def deutsch(oracle):
    n_qubits = round(np.log2(np.sqrt(oracle.size)))

    qubits = [KET0 for _ in range(n_qubits)]
    qubits[-1] = KET1

    whole_hadamard = tensordot_arb(*[H for _ in range(n_qubits)])

    entangled = whole_hadamard @ tensordot_arb(*qubits)

    entangled = oracle @ entangled

    entangled = whole_hadamard @ entangled

    return measure_system(entangled)


print("f(x) = 0: ", deutsch(tensordot_arb(I, I)))
print("f(x) = 1: ", deutsch(tensordot_arb(I, X)))
print("f(x) = x: ", deutsch(CNOT))
print("f(x) = !x: ", deutsch(tensordot_arb(X, I) @ CNOT @ tensordot_arb(X, I)))

print("f(x) = x mod 2 (3 qubits): ", deutsch(gen_c_operator(0, 8, X)))
