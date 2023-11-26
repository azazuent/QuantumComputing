from Qubit import *


def sample_bit(q: Qubit):
    q.h()
    result = q.measure()
    q.reset()
    return result


def prepare_message_qubit(message, basis, q: Qubit):
    if message:
        q.x()
    if basis:
        q.h()


def measure_message_qubit(basis: bool, q: Qubit):
    if basis:
        q.h()
    result = q.measure()
    q.reset()
    return result


def convert_to_hex(bits):
    return hex(int("".join(["1" if bit else "0" for bit in bits]), 2))


def send_single_bit_with_bb84(
 your_qubit: Qubit,
 eve_qubit: Qubit
 ):
    [your_message, your_basis] = [
        sample_bit(your_qubit) for _ in range(2)
        ]

    eve_basis = sample_bit(eve_qubit)

    if your_message:
        your_qubit.x()
    if your_basis:
        your_qubit.h()

    # ОТПРАВКА КУБИТА...

    eve_result = measure_message_qubit(eve_basis, your_qubit)
    return (your_message, your_basis), (eve_result, eve_basis)


def simulate_bb84(n_bits: int):
    your_device = Qubit()
    eve_device = Qubit()
    key = []
    n_rounds = 0
    while len(key) < n_bits:
        n_rounds += 1
        ((your_message, your_basis), (eve_result, eve_basis)) = send_single_bit_with_bb84(your_device, eve_device)
        if your_basis == eve_basis:
            assert your_message == eve_result
            key.append(your_message)
    print(f"Потребовалось {n_rounds} раундов, чтобы сгенерировать + {n_bits} -битовый ключ.")
    return key


def apply_one_time_pad(message, key):
    return [message_bit ^ key_bit for (message_bit, key_bit) in zip(message, key)]

if __name__ == "__main__":
    print("Генерирование 96-битового ключа путем симулирования BB84...")
    key = simulate_bb84(96)
    print(f"Получен ключ {convert_to_hex(key)}.")
    message = [
    1, 1, 0, 1, 1, 0, 0, 0,
    0, 0, 1, 1, 1, 1, 0, 1,
    1, 1, 0, 1, 1, 1, 0, 0,
    1, 0, 0, 1, 0, 1, 1, 0,
    1, 1, 0, 1, 1, 0, 0, 0,
    0, 0, 1, 1, 1, 1, 0, 1,
    1, 1, 0, 1, 1, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 1,
    1, 1, 0, 1, 1, 0, 0, 0,
    0, 0, 1, 1, 1, 1, 0, 1,
    1, 1, 0, 1, 1, 1, 0, 0,
    1, 0, 1, 1, 1, 0, 1, 1
    ]
    print(f"Использование ключа для отправки секретного сообщения: {convert_to_hex(message)}.")
    encrypted_message = apply_one_time_pad(message, key)
    print(f"Зашифрованное сообщение: {convert_to_hex(encrypted_message)}.")
    decrypted_message = apply_one_time_pad(encrypted_message, key)
    print(f"Ева расшифровала, получив: {convert_to_hex(decrypted_message)}.")

