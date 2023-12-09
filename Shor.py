from math import gcd, sqrt, ceil, log
from random import randint
from Qubit import *


def calc_base(n):
    if n == 1:
        return 1
    for i in range(2, ceil(sqrt(n)) + 1):
        possible_base = log(n, i)
        if possible_base.is_integer():
            return i
    return None


def find_order(x):
    n_qubits = ceil(log(x, 2)) + 1
    system = tensordot_arb(*[KET0 for _ in range(n_qubits)], KET1, *[KET0 for _ in range(4)])
    system = tensordot_arb(*[H for _ in range(n_qubits)], *[I for _ in range(5)]) @ system
    print(get_p_of_all_states(system))
    for i in range(n_qubits):
        mod = tensordot_arb(*[I for _ in range(i)],
                            KET1 @ KET1.transpose(),
                            *[I for _ in range(n_qubits - i - 1)],
                            mod2x21())
        mod += tensordot_arb(*[I for _ in range(i)],
                             KET0 @ KET0.transpose(),
                             *[I for _ in range(n_qubits - i + 5 - 1)])
        system = mod @ system
        print(get_p_of_all_states(system))
    system = tensordot_arb(gen_QFT(n_qubits), *[I for _ in range(5)]) @ system

    for i in range(n_qubits // 2):
        system = tensordot_arb(*[I for _ in range(i)],
                               gen_swap(n_qubits - 2 * i),
                               *[I for _ in range(i)],
                               *[I for _ in range(5)]) @ system

    print(get_p_of_all_states(system)[:n_qubits])


def shor(n):
    if not (n & 0b1):
        return 2

    if calc_base(n):
        return calc_base(n)

    x = randint(2, n - 1)
    if gcd(x, n) > 1:
        return gcd(x, n), n // gcd(x, n)

    r = 6
    # r = find_order(n)

    guess1 = x ** (r // 2) - 1
    guess2 = x ** (r // 2) + 1
    if 1 in [gcd(guess1, n), gcd(guess2, n)]:
        print(x)
    return gcd(n, guess1), gcd(n, guess2)


find_order(21)
