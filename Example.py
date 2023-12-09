from Qubit import *


q1 = Qubit()
q2 = Qubit()
q3 = Qubit()
q4 = Qubit()

entangled = tensordot_arb(q1.state, q2.state, q3.state, q4.state)
entangled = tensordot_arb(H, I, I, I) @ entangled
entangled = tensordot_arb(CH, I, I) @ entangled
entangled = tensordot_arb(I, CH, I) @ entangled
entangled = tensordot_arb(I, I, CH) @ entangled

results = [get_p_of_state(entangled, i) for i in range(4)]
print(*results)
