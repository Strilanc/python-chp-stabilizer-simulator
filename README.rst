Python CHP Stabilizer Simulator
-------------------------------

A simple reference python implementation of Scott Aaronson and Daniel Gottesman's CHP simulator
as defined in
`their 2004 paper "Improved Simulation of Stabilizer Circuits" <https://arxiv.org/abs/quant-ph/0406196>`__.
This simulator is capable of simulating quantum stabilizer circuits in polynomial time and space.
Specifically, it uses ``O(q^2*m + q*c)`` time and ``O(q^2)`` space where
``q`` is the number of qubits,
``m`` is the number of measurements,
and ``c`` is the number of Hadamard/CNOT/Phase gates.

Installation
------------

The ``chp_sim`` package is available on pypi and can be installed using ``pip``:

.. code-block:: bash

    python -m pip install chp_sim

Alternatively, you can just copy paste the ``chp_sim`` directory of the github
repository into your project.
The only runtime dependency is ``numpy``.

Usage
-----

Here is an example of simulating
`a circuit <https://algassert.com/quirk#circuit=%7B%22cols%22%3A%5B%5B1%2C1%2C%22H%22%5D%2C%5B%22X%22%2C1%2C%22%E2%80%A2%22%5D%2C%5B1%2C%22X%22%2C%22%E2%80%A2%22%5D%2C%5B%22Z%5E%C2%BD%22%2C%22Z%5E%C2%BD%22%5D%2C%5B%22H%22%2C%22H%22%2C%22H%22%5D%2C%5B%22Measure%22%2C%22Measure%22%2C%22Measure%22%5D%2C%5B%22Chance3%22%5D%5D%7D>`__:

.. code-block:: python

    import chp_sim
    sim = chp_sim.ChpSimulator(num_qubits=3)

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


Packaging
---------

(Notes to self on how to release a new version.)

1. Edit the source code as needed and run tests.

    .. code-block:: bash

        pytest

2. Build the wheel.

    .. code-block:: bash

        python3 setup.py -q bdist_wheel
        ls dist

3. Upload to test pypi.

    .. code-block:: bash

        twine upload dist/*.whl --repository-url=https://test.pypi.org/legacy/ --username="${TEST_TWINE_USERNAME}" --password="${TEST_TWINE_PASSWORD}"

4. Verify the test package works.

    .. code-block:: bash

        mkvirtualenv test --python=/usr/bin/python3
        pip install numpy
        pip install chp_sim --index-url=https://test.pypi.org/simple/
        python -c "import chp_sim; print(chp_sim.__version__); print(chp_sim.ChpSimulator(4))"

5. Upload to prod pypi.

    .. code-block:: bash

        twine upload dist/*.whl --username="${PROD_TWINE_USERNAME}" --password="${PROD_TWINE_PASSWORD}"

6. Verify the prod package works.

    .. code-block:: bash

        mkvirtualenv test --python=/usr/bin/python3
        pip install chp_sim
        python -c "import chp_sim; print(chp_sim.__version__); print(chp_sim.ChpSimulator(4))"
