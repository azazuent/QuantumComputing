from Qubit import *


def simon(oracle):

    n_qubits = round(np.log2(oracle.shape[0]))

    half_h = [H for _ in range(n_qubits // 2)]
    half_i = [I for _ in range(n_qubits // 2)]

    simon_op = tensordot_arb(*half_h, *half_i)

    qubits = tensordot_arb(*[KET0 for _ in range(n_qubits)])

    qubits = simon_op @ qubits
    qubits = oracle @ qubits
    qubits = simon_op @ qubits

    return measure_system(qubits)


oracle = tensordot_arb(gen_c_operator(0, 3, X), I, I)
oracle = tensordot_arb(I, I, gen_c_operator(0, 3, X)) @ oracle

print("Simon`s algorithm result:\n", simon(oracle))
