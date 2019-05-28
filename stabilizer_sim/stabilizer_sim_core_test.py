import cirq

from .stabilizer_sim_core import pauli_product_phase, StabilizerSimCore
from . import MeasureResult


def test_pauli_product_phase():
    q = cirq.LineQubit(0)
    mapping = {
        0: cirq.I(q),
        1: cirq.X(q),
        2: cirq.Z(q),
        3: cirq.Y(q)
    }
    for k1, v1 in mapping.items():
        for k2, v2 in mapping.items():
            x1 = bool(k1 & 1)
            z1 = bool(k1 & 2)
            x2 = bool(k2 & 1)
            z2 = bool(k2 & 2)
            p = pauli_product_phase(x1, z1, x2, z2)
            if k1 == 0 or k2 == 0:
                assert p == 0
            else:
                assert [1, 1j, None, -1j][p] == (v1 * v2).coefficient


def test_identity():
    s = StabilizerSimCore(num_qubits=1)
    assert s.measure(0) == MeasureResult(value=False, determined=True)


def test_bit_flip():
    s = StabilizerSimCore(num_qubits=1)
    s.hadamard(0)
    s.phase(0)
    s.phase(0)
    s.hadamard(0)
    assert s.measure(0) == MeasureResult(value=True, determined=True)


def test_identity_2():
    s = StabilizerSimCore(num_qubits=2)
    assert s.measure(0) == MeasureResult(value=False, determined=True)
    assert s.measure(1) == MeasureResult(value=False, determined=True)


def test_bit_flip_2():
    s = StabilizerSimCore(num_qubits=2)
    s.hadamard(0)
    s.phase(0)
    s.phase(0)
    s.hadamard(0)
    assert s.measure(0) == MeasureResult(value=True, determined=True)
    assert s.measure(1) == MeasureResult(value=False, determined=True)


def test_epr():
    s = StabilizerSimCore(num_qubits=2)
    s.hadamard(0)
    s.cnot(0, 1)
    v1 = s.measure(0)
    assert not v1.determined
    v2 = s.measure(1)
    assert v2.determined
    assert v1.value == v2.value


def test_phase_kickback_consume_s_state():
    s = StabilizerSimCore(num_qubits=2)
    s.hadamard(1)
    s.phase(1)
    s.hadamard(0)
    s.cnot(0, 1)
    v1 = s.measure(1)
    assert not v1.determined
    if v1:
        s.phase(0)
        s.phase(0)
    s.phase(0)
    s.hadamard(0)
    assert s.measure(0) == MeasureResult(value=True, determined=True)


def test_phase_kickback_preserve_s_state():
    s = StabilizerSimCore(num_qubits=2)

    # Prepare S state.
    s.hadamard(1)
    s.phase(1)

    # Prepare test input.
    s.hadamard(0)

    # Kickback.
    s.cnot(0, 1)
    s.hadamard(1)
    s.cnot(0, 1)
    s.hadamard(1)

    # Check.
    s.phase(0)
    s.hadamard(0)
    assert s.measure(0) == MeasureResult(value=True, determined=True)
    s.phase(1)
    s.hadamard(1)
    assert s.measure(1) == MeasureResult(value=True, determined=True)


def test_kickback_vs_stabilizer():
    sim = StabilizerSimCore(num_qubits=3)
    sim.hadamard(2)
    sim.cnot(2, 0)
    sim.cnot(2, 1)
    sim.phase(0)
    sim.phase(1)
    sim.hadamard(0)
    sim.hadamard(1)
    sim.hadamard(2)
    assert str(sim).strip() == """
-YII
-IYI
+IIX
----
+XIX
+IXX
+YYZ
    """.strip().replace('I', '.')
    v0 = sim.measure(0, bias=0)
    assert str(sim).strip() == """
+XIX
-IYI
+IIX
----
+ZII
+IXX
+ZYY
    """.strip().replace('I', '.')
    v1 = sim.measure(1, bias=0)
    assert str(sim).strip() == """
+XIX
+IXX
+IIX
----
+ZII
+IZI
-ZZZ
    """.strip().replace('I', '.')
    v2 = sim.measure(2, bias=0)
    assert str(sim).strip() == """
+XIX
+IXX
+IIX
----
+ZII
+IZI
-ZZZ
    """.strip().replace('I', '.')
    assert v0 == MeasureResult(value=False, determined=False)
    assert v1 == MeasureResult(value=False, determined=False)
    assert v2 == MeasureResult(value=True, determined=True)


def test_s_state_distillation_low_depth():
    for _ in range(100):
        sim = StabilizerSimCore(num_qubits=9)

        stabilizers = [
            (0, 1, 2, 3),
            (0, 1, 4, 5),
            (0, 2, 4, 6),
            (1, 2, 4, 7)
        ]
        checks = [
            {'s': [0], 'q': stabilizers[0]},
            {'s': [1], 'q': stabilizers[1]},
            {'s': [2], 'q': stabilizers[2]},
        ]

        stabilizer_measurements = []
        anc = 8
        for stabilizer in stabilizers:
            sim.hadamard(anc)
            for k in stabilizer:
                sim.cnot(anc, k)
            sim.hadamard(anc)
            v = sim.measure(anc)
            assert not v.determined
            if v.value:
                sim.hadamard(anc)
                sim.phase(anc)
                sim.phase(anc)
                sim.hadamard(anc)
            stabilizer_measurements.append(v)

        qubit_measurements = []
        for k in range(7):
            sim.phase(k)
            sim.hadamard(k)
            qubit_measurements.append(sim.measure(k))

        if sum([int(e.value) for e in stabilizer_measurements + qubit_measurements]) & 1:
            sim.phase(7)
            sim.phase(7)

        sim.phase(7)
        sim.hadamard(7)
        r = sim.measure(7)

        assert r == MeasureResult(value=False, determined=True)

        for c in checks:
            rvs = [int(stabilizer_measurements[k].value) for k in c['s']]
            rms = [int(qubit_measurements[k].value) for k in c['q']]
            assert sum(rvs + rms) & 1 == 0


def test_s_state_distillation_low_space():
    for _ in range(100):
        sim = StabilizerSimCore(num_qubits=5)

        phasors = [
            (0,),
            (1,),
            (2,),
            (0, 1, 2),
            (0, 1, 3),
            (0, 2, 3),
            (1, 2, 3),
        ]

        anc = 4
        for phasor in phasors:
            sim.hadamard(anc)
            for k in phasor:
                sim.cnot(anc, k)
            sim.hadamard(anc)
            sim.phase(anc)
            sim.hadamard(anc)
            v = sim.measure(anc)
            assert not v.determined
            if v.value:
                for k in phasor + (anc,):
                    sim.hadamard(k)
                    sim.phase(k)
                    sim.phase(k)
                    sim.hadamard(k)

        for k in range(3):
            assert sim.measure(k) == MeasureResult(value=False, determined=True)
        sim.phase(3)
        sim.hadamard(3)
        assert sim.measure(3) == MeasureResult(value=True, determined=True)
