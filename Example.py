from Qubit import *


q1 = Qubit()
q2 = Qubit()
q3 = Qubit()
q4 = Qubit()

entangled = entangle(q1.state, q2.state, q3.state, q4.state)
entangled = entangle(H, I, I, I) @ entangled
entangled = entangle(CH, I, I) @ entangled
entangled = entangle(I, CH, I) @ entangled
entangled = entangle(I, I, CH) @ entangled

results = [get_p_of_1(entangled, i) for i in range(4)]
print(*results)
