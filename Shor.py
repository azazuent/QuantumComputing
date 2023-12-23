from math import gcd, sqrt, ceil, log
from Qubit import *


def calc_base(n):
    if n == 1:
        return 1
    for i in range(2, ceil(sqrt(n)) + 1):
        possible_base = log(n, i)
        if possible_base.is_integer():
            return i, round(possible_base)
    return None


def find_order(x):

    n_qubits = 7
    system = tensordot_arb(*[KET_PLUS for _ in range(n_qubits)], KET1, *[KET0 for _ in range(4)])

    mod_op = n_mod_21(x)
    print("Calculating modulus operators...")
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

    print(f"Generating QFT for {n_qubits} qubits...")
    system = tensordot_arb(gen_QFT(n_qubits), *[I for _ in range(5)]) @ system
    for i in range(n_qubits // 2):
        system = tensordot_arb(*[I for _ in range(i)],
                               gen_swap(n_qubits - 2 * i),
                               *[I for _ in range(i)],
                               *[I for _ in range(5)]) @ system

    print("Measuring...")
    probs = measure_system(system, 0, n_qubits-1)
    measurement = shoot(probs)
    measurement = sorted(measurement.items(), key=lambda i: i[1], reverse=True)
    for i in measurement:
        if i[0] != 0:
            possible_r = 2 ** n_qubits // i[0]
            if x ** possible_r % 21 == 1:
                return possible_r


def shor(n):

    if not (n & 0b1):
        print("Provided number is even")
        return 2, n // 2

    if calc_base(n):
        base, power = calc_base(n)
        print(f"Provided number is {base} ^ {power}")
        return calc_base(n)[0]

    if not n == 21:
        raise Exception("This implementation doesn't factor this number yet")

    pick_pool = [i for i in range(2, 21)]
    while True:
        try:
            x = random.choice(pick_pool)
            #x = 13
            pick_pool.remove(x)
            print(f"Randomly picked number {x}")

            gcd_value = gcd(x, n)
            print(f"gcd({x}, {n}) = {gcd_value}")

            if gcd_value > 1:
                return gcd_value, n // gcd_value

            print("Provided number will be factorised with quantum computation")

            r = find_order(x)
            print(f"Period for {x}^n mod 21 is {r}")
            if r & 0b1:
                raise Exception("Period value is odd, let's pick another number")

            guess1 = x ** (r // 2) - 1
            guess2 = x ** (r // 2) + 1
            if 21 in [gcd(guess1, n), gcd(guess2, n)]:
                raise Exception("One of the factors is 21, let's pick another number")
            return gcd(n, guess1), gcd(n, guess2)
        except Exception as e:
            print(str(e))

try:
    n = int(input("Enter the number you need factorised: "))
    print(f"Factors of {n} are: {', '.join([str(i) for i in sorted(shor(n))])}")
except Exception as e:
    print(str(e))
