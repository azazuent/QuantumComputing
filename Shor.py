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


def find_order(x, n):

    n_qubits = 6
    system = tensordot_arb(*[KET_PLUS for _ in range(n_qubits)], KET1, *[KET0 for _ in range(4)])

    mod_op = n_mod_21(x)
    print("Calculating modulus operators")
    for i in range(n_qubits):

        mod = tensordot_arb(*[I for _ in range(i)],
                            KET1 @ KET1.transpose(),
                            *[I for _ in range(n_qubits - i - 1)],
                            mod_op)
        mod += tensordot_arb(*[I for _ in range(i)],
                             KET0 @ KET0.transpose(),
                             *[I for _ in range(n_qubits - i + 5 - 1)])
        mod_op = mod_op @ mod_op

        system = mod @ system

    print(f"Generating QFT for {n_qubits} qubits")
    system = tensordot_arb(gen_QFT(n_qubits), *[I for _ in range(5)]) @ system
    for i in range(n_qubits // 2):
        system = tensordot_arb(*[I for _ in range(i)],
                               gen_swap(n_qubits - 2 * i),
                               *[I for _ in range(i)],
                               *[I for _ in range(5)]) @ system
    print("Isolating first register (might take a while)")
    measurement = measure_system(system, 0, n_qubits-1)
    print("Taking probes")
    shoot(measurement, 4000)


def shor(n):

    if not (n & 0b1):
        return 2

    if calc_base(n):
        return calc_base(n)

    if not n == 21:
        raise Exception("This implementation doesn't factor this number yet")

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


find_order(16, 21)
