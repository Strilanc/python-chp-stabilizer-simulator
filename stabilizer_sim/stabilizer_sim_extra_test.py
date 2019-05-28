import collections
import itertools
from typing import Set

from . import StabilizerSim, MeasureResult


def test_s_state_distillation_low_depth():
    for _ in range(100):
        sim = StabilizerSim(num_qubits=9)

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
            for k in stabilizer:
                sim.xnot(anc, k)
            v = sim.measure_z_and_reset(anc)
            assert not v.determined
            stabilizer_measurements.append(v)

        qubit_measurements = []
        for k in range(7):
            sim.phase(k)
            sim.hadamard(k)
            qubit_measurements.append(sim.measure(k))

        if sum([int(e.value) for e in stabilizer_measurements + qubit_measurements]) & 1:
            sim.z(7)

        for c in checks:
            rvs = [int(stabilizer_measurements[k].value) for k in c['s']]
            rms = [int(qubit_measurements[k].value) for k in c['q']]
            assert sum(rvs + rms) & 1 == 0

        assert sim.measure_y(7) == MeasureResult(value=True, determined=True)


def test_s_state_distillation_low_space():
    for _ in range(100):
        sim = StabilizerSim(num_qubits=5)

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
            for k in phasor:
                sim.xnot(anc, k)
            sim.phase(anc)
            v = sim.measure_x_and_reset(anc)
            assert not v.determined
            if v.value:
                for k in phasor:
                    sim.x(k)

        for k in range(3):
            assert sim.measure(k) == MeasureResult(value=False, determined=True)
        assert sim.measure_y(3) == MeasureResult(value=False, determined=True)


def test_count_s_state_distillation_failure_cases():
    def distill(errors: Set[int]) -> str:
        sim = StabilizerSim(num_qubits=5)

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
        for e, phasor in enumerate(phasors):
            for k in phasor:
                sim.xnot(anc, k)
            sim.phase(anc)

            if e in errors:
                sim.z(anc)

            v = sim.measure_x_and_reset(anc)
            assert not v.determined
            if v.value:
                for k in phasor:
                    sim.x(k)

        result = sim.measure_y(3)
        checks = [sim.measure(k) for k in range(3)]
        assert result.determined
        assert all(e.determined for e in checks)
        good_result = result.value == False
        checks_passed = not any(e.value for e in checks)
        if checks_passed:
            if good_result:
                return 'good'
            else:
                return 'ERROR'
        else:
            if good_result:
                return 'victim'
            else:
                return 'caught'

    def classify(errs) -> collections.Counter:
        result = collections.Counter()
        for err in errs:
            result[distill(err)] += 1
        return result

    nones = list(itertools.combinations(range(7), 0))
    singles = list(itertools.combinations(range(7), 1))
    doubles = list(itertools.combinations(range(7), 2))
    triples = list(itertools.combinations(range(7), 3))

    assert classify(nones) == {'good': 1}
    assert classify(singles) == {'caught': 3, 'victim': 4}
    assert classify(doubles) == {'caught': 12, 'victim': 9}
    assert classify(triples) == {'caught': 12, 'victim': 16, 'ERROR': 7}
