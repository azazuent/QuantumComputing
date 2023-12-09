from Qubit import *
from numpy import pi

games = 10000
wins = 0
for j in range(games):
    a = q_rng()
    b = q_rng()
    A = E = 0
    if (a and b) == (A != E):
        wins += 1
print("Классическая стратегия - доля побед: " + str(100*wins/games) + "%")

wins = 0
for j in range(games):
    a = q_rng()
    b = q_rng()
    A = E = 0
    entangled = tensordot_arb(KET0, KET0)
    if not a:
        entangled = tensordot_arb(turn(pi / 2, 'y'), I) @ entangled
    else:
        entangled = tensordot_arb(turn(0, 'y'), I) @ entangled
    if not b:
        entangled = tensordot_arb(I, turn(np.pi / 4, 'y')) @ entangled
    else:
        entangled = tensordot_arb(I, turn(3 * np.pi / 4, 'y')) @ entangled
    prob = np.power(entangled, 2)
    r = np.random.random()
    if r < prob[0]:
        A = E = 0
        r -= prob[0]
    elif r < prob[1]:
        A = 0
        E = 1
        r -= prob[1]
    elif r < prob[2]:
        A = 1
        E = 0
    else:
        A = E = 1
    if (a and b) == (A != E):
        wins += 1
print("Квантовая стратегия - доля побед: " + str(100*wins/games) + "%")
