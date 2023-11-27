from Qubit import *


def deutsch(oracle):

    qubit1 = Qubit()
    qubit2 = Qubit()

    qubit2.x()

    qubit1.h()
    qubit2.h()

    entangled = entangle(qubit1.state, qubit2.state)

    entangled = oracle @ entangled

    hi = entangle(H, I)

    entangled = hi @ entangled

    return get_p_of_state(entangled, 0)


print(deutsch(entangle(I, I)))
print(deutsch(entangle(I, X)))
print(deutsch(CNOT))
print(deutsch(entangle(X, I) @ CNOT @ entangle(X, I)))
