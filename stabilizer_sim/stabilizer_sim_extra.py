from typing import Union

from .stabilizer_sim_core import StabilizerSimCore, MeasureResult


class StabilizerSim(StabilizerSimCore):
    """CHP stabilizer simulator with some utility methods for compound gates.

    Reference:
        "Improved Simulation of Stabilizer Circuits"
        Scott Aaronson and Daniel Gottesman
        https://arxiv.org/abs/quant-ph/0406196
    """

    def x(self, a: int) -> None:
        """Pauli X operation."""
        self.hadamard(a)
        self.phase(a)
        self.phase(a)
        self.hadamard(a)

    def y(self, a: int) -> None:
        """Pauli Y operation."""
        self.phase(a)
        self.phase(a)
        self.hadamard(a)
        self.phase(a)
        self.phase(a)
        self.hadamard(a)

    def z(self, a: int) -> None:
        """Pauli Z operation."""
        self.phase(a)
        self.phase(a)

    def sqrt_x(self, a: int):
        """+90 degree rotation around X axis."""
        self.hadamard(a)
        self.phase(a)
        self.hadamard(a)

    def sqrt_x_dag(self, a: int) -> None:
        """-90 degree rotation around X axis."""
        self.hadamard(a)
        self.phase(a)
        self.phase(a)
        self.phase(a)
        self.hadamard(a)

    def sqrt_z(self, a: int):
        """+90 degree rotation around Z axis."""
        self.phase(a)

    def sqrt_z_dag(self, a: int) -> None:
        """-90 degree rotation around Z axis."""
        self.phase(a)
        self.phase(a)
        self.phase(a)

    def h_xz(self, a: int) -> None:
        """180 degree rotation around X+Z."""
        self.hadamard(a)

    def h_yz(self, a: int) -> None:
        """180 degree rotation around Y+Z."""
        self.hadamard(a)
        self.phase(a)
        self.hadamard(a)
        self.phase(a)
        self.phase(a)

    def h_xy(self, a: int) -> None:
        """180 degree rotation around X+Y."""
        self.hadamard(a)
        self.phase(a)
        self.phase(a)
        self.hadamard(a)
        self.phase(a)

    def measure_x(self, a: int, *, bias: Union[bool, int, float]=0.5
                  ) -> MeasureResult:
        """X basis measurement."""
        self.hadamard(a)
        v = self.measure(a, bias=bias)
        self.hadamard(a)
        return v

    def measure_y(self, a: int, *, bias: Union[bool, int, float]=0.5
                  ) -> MeasureResult:
        """Y basis measurement."""
        self.h_yz(a)
        v = self.measure(a, bias=bias)
        self.h_yz(a)
        return v

    def measure_z(self, a: int, *, bias: Union[bool, int, float]=0.5
                  ) -> MeasureResult:
        """Z basis measurement."""
        return self.measure(a, bias=bias)

    def measure_x_and_reset(self, a: int, *, bias: Union[bool, int, float]=0.5
                            ) -> MeasureResult:
        """X basis measurement followed by a reset."""
        self.hadamard(a)
        return self.measure_z_and_reset(a, bias=bias)

    def measure_y_and_reset(self, a: int, *, bias: Union[bool, int, float]=0.5
                            ) -> MeasureResult:
        """Y basis measurement followed by a reset."""
        self.h_yz(a)
        return self.measure_z_and_reset(a, bias=bias)

    def measure_z_and_reset(self, a: int, *, bias: Union[bool, int, float]=0.5
                            ) -> MeasureResult:
        """Z basis measurement followed by a reset."""
        v = self.measure(a, bias=bias)
        if v.value:
            self.x(a)
        return v

    def xnot(self, a: int, b: int) -> None:
        """An X gate controlled by an X-axis control."""
        self.hadamard(a)
        self.cnot(a, b)
        self.hadamard(a)

    def cz(self, a: int, b: int) -> None:
        """A Z gate controlled by a Z-axis control."""
        self.hadamard(b)
        self.cnot(a, b)
        self.hadamard(b)
