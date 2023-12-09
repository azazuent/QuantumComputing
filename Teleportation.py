from Qubit import *
import numpy as np


def q_teleport(q_entity, bell_state=0):

    if not 0 <= bell_state <= 3:
        raise Exception("Bell state indexes are in range [0, 3]")

    op1 = tensordot_arb(I, H, I)

    if bell_state == 1 or bell_state == 3:
        op1 = tensordot_arb(I, X, I) @ op1

    if bell_state == 2 or bell_state == 3:
        op2 = tensordot_arb(I, NCNOT)
    else:
        op2 = tensordot_arb(I, CNOT)

    op5 = op2

    op3 = tensordot_arb(CNOT, I)

    op4 = tensordot_arb(H, I, I)

    if bell_state == 1:
        op4 = tensordot_arb(H, I, Z)
    elif bell_state == 3:
        op4 = tensordot_arb(Z, Z, I) @ op4

    op6 = gen_c_operator(0, 2, Z)

    q_entity = op1 @ q_entity
    q_entity = op2 @ q_entity
    q_entity = op3 @ q_entity
    q_entity = op4 @ q_entity
    q_entity = op5 @ q_entity
    q_entity = op6 @ q_entity

    return q_entity


qubit_1 = KET0
qubit_2 = KET0
qubit_3 = KET0

qubit_1 = turn(np.pi / 6, 'y') @ qubit_1

entangled = tensordot_arb(qubit_1, qubit_2, qubit_3)

print("Initial probabilities: ", get_p_of_all_states(entangled))

print("Fi+ teleport: ", get_p_of_all_states(q_teleport(entangled.copy(), bell_state=0)))

print("Fi- teleport: ", get_p_of_all_states(q_teleport(entangled.copy(), bell_state=1)))

print("Psi+ teleport: ", get_p_of_all_states(q_teleport(entangled.copy(), bell_state=2)))

print("Psi- teleport: ", get_p_of_all_states(q_teleport(entangled.copy(), bell_state=3)))

