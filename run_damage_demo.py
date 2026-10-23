
import numpy as np
from reentrant_fea.ductile_kernel import State, integrate

E = 69000.0
nu = 0.3

plastic_table = [
    (200.0, 0.0),
    (220.0, 0.01),
    (250.0, 0.10),
]

fracture_locus = [
    (-1.0, 0.047),
    (0.0, 0.047),
    (1.0 / 3.0, 0.047),
    (1.0, 0.047),
]

state = State.virgin()

for step in range(1, 101):
    strain = np.array([
        0.001 * step, 0, 0, 0, 0, 0
    ])

    effective, damaged, candidate = integrate(
        strain,
        state,
        E,
        nu,
        plastic_table,
        fracture_locus,
        characteristic_length=0.3,
        failure_displacement=0.001,
    )

    state = candidate

    print(
        f"Step={step:03d} | "
        f"Effective stress={effective[0]:.2f} MPa | "
        f"Damaged stress={damaged[0]:.2f} MPa | "
        f"State={state}"
    )
