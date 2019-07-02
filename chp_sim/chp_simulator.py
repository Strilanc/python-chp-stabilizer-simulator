import random
from typing import Union, Any

import numpy as np


class ChpSimulator:
    """The bare minimum needed for the CHP simulation.

    Reference:
        "Improved Simulation of Stabilizer Circuits"
        Scott Aaronson and Daniel Gottesman
        https://arxiv.org/abs/quant-ph/0406196
    """

    def __init__(self, num_qubits):
        self._n = num_qubits
        self._table = np.eye(2 * num_qubits + 1, dtype=np.bool)
        self._x = self._table[:, :self._n]
        self._z = self._table[:, self._n:-1]
        self._r = self._table[:, -1]

    def cnot(self, control: int, target: int) -> None:
        """Applies a CNOT gate between two qubits.

        Args:
            control: The control qubit of the CNOT.
            target: The target qubit of the CNOT.
        """
        self._r[:] ^= self._x[:, control] & self._z[:, target] & (
                self._x[:, target] ^ self._z[:, control] ^ True)
        self._x[:, target] ^= self._x[:, control]
        self._z[:, control] ^= self._z[:, target]

    def hadamard(self, qubit: int) -> None:
        """Applies a Hadamard gate to a qubit.

        Args:
            qubit: The qubit to apply the H gate to.
        """
        self._r[:] ^= self._x[:, qubit] & self._z[:, qubit]
        # Perform a XOR-swap
        self._x[:, qubit] ^= self._z[:, qubit]
        self._z[:, qubit] ^= self._x[:, qubit]
        self._x[:, qubit] ^= self._z[:, qubit]

    def phase(self, qubit: int) -> None:
        """Applies an S gate to a qubit.

        Args:
            qubit: The qubit to apply the S gate to.
        """
        self._r[:] ^= self._x[:, qubit] & self._z[:, qubit]
        self._z[:, qubit] ^= self._x[:, qubit]

    def measure(self,
                qubit: int,
                *,
                bias: Union[float, int, bool] = 0.5) -> 'MeasureResult':
        """Computational basis (Z basis) measurement.

        Args:
            qubit: The index of the qubit to measure.
            bias: When the measurement result is random, this is the probability
                of getting a True result value instead of False. Mostly useful
                for making reproducible unit tests.

        Returns:
            A MeasurementResult instance whose `value` attribute is the outcome
            of the measurement and whose `determined` attribute indicates
            whether the outcome was deterministic or random.
        """
        n = self._n
        for p in range(n):
            if self._x[p+n, qubit]:
                return self._measure_random(qubit, p, bias)
        return self._measure_determined(qubit)

    def _measure_random(self,
                        a: int,
                        p: int,
                        bias: Union[float, int, bool]) -> 'MeasureResult':
        n = self._n
        assert self._x[p+n, a]
        self._table[p, :] = self._table[p + n, :]
        self._table[p + n, :] = 0
        self._z[p + n, a] = 1
        self._r[p + n] = random.random() < bias

        for i in range(2*n):
            if self._x[i, a] and i != p and i != p + n:
                self._row_mult(i, p)
        return MeasureResult(value=self._r[p + n], determined=False)

    def _measure_determined(self, a: int) -> 'MeasureResult':
        n = self._n
        self._table[-1, :] = 0
        for i in range(n):
            if self._x[i, a]:
                self._row_mult(-1, i + n)
        return MeasureResult(value=self._r[-1], determined=True)

    def _row_product_sign(self, i: int, k: int) -> bool:
        """Determines the sign of two rows' Pauli Products."""
        pauli_phases = sum(
            pauli_product_phase(self._x[i, j],
                                self._z[i, j],
                                self._x[k, j],
                                self._z[k, j])
            for j in range(self._n)
        )
        assert not pauli_phases & 1, (
            "Expected commuting rows but got {}, {} from \n{}".format(
                i, k, self))
        p = (pauli_phases >> 1) & 1
        return bool(self._r[i] ^ self._r[k] ^ p)

    def _row_mult(self, i: int, k: int) -> None:
        """Multiplies row k's Paulis into row i's Paulis."""
        self._r[i] = self._row_product_sign(i, k)
        self._x[i, :self._n] ^= self._x[k, :self._n]
        self._z[i, :self._n] ^= self._z[k, :self._n]

    def __str__(self):
        """Represents the state as a list of Pauli products.

        Each Pauli product is what the X or Z observable of a qubit at the
        current time was equal to at time zero (after accounting for
        measurements).
        """

        def _cell(row: int, col: int) -> str:
            k = int(self._x[row, col]) + 2 * int(self._z[row, col])
            return ['.', 'X', 'Z', 'Y'][k]

        def _row(row: int) -> str:
            result = '-' if self._r[row] else '+'
            for col in range(self._n):
                result += str(_cell(row, col))
            return result

        z_obs = [_row(row) for row in range(self._n)]
        sep = ['-' * (self._n + 1)]
        x_obs = [_row(row) for row in range(self._n, 2 * self._n)]
        return '\n'.join(z_obs + sep + x_obs)

    def _repr_pretty_(self, p: Any, cycle: bool) -> None:
        p.text(str(self))


def pauli_product_phase(x1: bool, z1: bool, x2: bool, z2: bool) -> int:
    """Determines the power of i in the product of two Paulis.

    For example, X*Y = iZ and so this method would return +1 for X and Y.

    The input Paulis are encoded into the following form:

        x z | Pauli
        ----+-------
        0 0 | I
        1 0 | X
        1 1 | Y
        0 1 | Z
    """
    # Analyze by case over first gate.

    if x1 and z1:  # Y gate.
        # No phase for YI = Y
        # -1 phase for YX = -iZ
        # No phase for YY = I
        # +1 phase for YZ = +iX
        return int(z2) - int(x2)

    if x1:  # X gate.
        # No phase for XI = X
        # No phase for XX = I
        # +1 phase for XY = iZ
        # -1 phase for XZ = -iY
        return z2 and 2*int(x2) - 1

    if z1:  # Z gate.
        # No phase for ZI = Z
        # +1 phase for ZX = -iY
        # -1 phase for ZY = iX
        # No phase for ZZ = I
        return x2 and 1 - 2*int(z2)

    # Identity gate.
    return 0


class MeasureResult:
    """A measurement's output and whether it was random or not."""

    def __init__(self, value: bool, determined: bool):
        self.value = bool(value)
        self.determined = bool(determined)

    def __bool__(self):
        return self.value

    def __eq__(self, other):
        if isinstance(other, (bool, int)):
            return self.value == other
        if isinstance(other, MeasureResult):
            return self.value == other.value and self.determined == other.determined
        return NotImplemented

    def __str__(self):
        return '{} ({})'.format(self.value,
                                ['random', 'determined'][self.determined])

    def __repr__(self):
        return 'MeasureResult(value={!r}, determined={!r})'.format(
            self.value,
            self.determined)
