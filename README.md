# Python CHP Stabilizer Simulator

A simple reference python implementation of Scott Aaronson and Daniel Gottesman's CHP simulator
as defined in [their 2004 paper "Improved Simulation of Stabilizer Circuits"](https://arxiv.org/abs/quant-ph/0406196).
This simulator is capable of simulating quantum stabilizer circuits in polynomial time and space.
Specifically, it uses `O(q^2*m + q*c)` time and `O(q^2)` space where
`q` is the number of qubits,
`m` is the number of measurements,
and `c` is the number of Hadamard/CNOT/Phase gates.

# Usage

Here is an example of simulating [a circuit](https://algassert.com/quirk#circuit=%7B%22cols%22%3A%5B%5B1%2C1%2C%22H%22%5D%2C%5B%22X%22%2C1%2C%22%E2%80%A2%22%5D%2C%5B1%2C%22X%22%2C%22%E2%80%A2%22%5D%2C%5B%22Z%5E%C2%BD%22%2C%22Z%5E%C2%BD%22%5D%2C%5B%22H%22%2C%22H%22%2C%22H%22%5D%2C%5B%22Measure%22%2C%22Measure%22%2C%22Measure%22%5D%2C%5B%22Chance3%22%5D%5D%7D):

```python
from stabilizer_sim import StabilizerSim
sim = StabilizerSim(num_qubits=3)

# Desired circuit:
# 0: -------X-------S---H---M---
#           |
# 1: -------|---X---S---H---M---
#           |   |
# 2: ---H---@---@-------H---M---

sim.hadamard(2)
sim.cnot(2, 0)
sim.cnot(2, 1)
sim.phase(0)
sim.phase(1)
sim.hadamard(0)
sim.hadamard(1)
sim.hadamard(2)

# Show internal simulator state.
print(sim, '\n')
# prints:
#   -Y..
#   -.Y.
#   +..X
#   ----
#   +X.X
#   +.XX
#   +YYZ

# Perform measurements
v0 = sim.measure(0)
v1 = sim.measure(1)
v2 = sim.measure(2)
print(v0)
print(v1)
print(v2)
# prints [note: one of four possible results for this circuit]:
#   True (random)
#   False (random)
#   False (determined)

# Check pattern the outputs should satisfy.
assert not v0.determined
assert not v1.determined
assert v2.determined
assert bool(v0) ^ bool(v1) ^ bool(v2)
```

# Installation

At the moment there is no special installation method such as a pypi package.
Just copy-paste the `stabilizer_sim` folder of this repository into your project.
The only runtime dependency is `numpy`.
