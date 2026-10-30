# User-supplied plastic hardening data

The user explicitly confirmed on **2026-10-27** that density is
`2.67e-9` Mg/mm³, elastic constants are `E = 69000` MPa and `nu = 0.3`, and the
supplied plastic data is **already true stress versus equivalent true plastic
strain**. [The material confirmation](user_material_confirmation.json) records
these fixed inputs and binds the supplied table, existing default source,
precision audit and current configuration with SHA-256 hashes. No further
engineering-to-true conversion is appropriate. Exploratory slopes read from a
printed coupon figure do not override this confirmation or establish that these
inputs are wrong.

`user_plastic_table_2026-10-26.inp` preserves the complete `*Plastic` block
supplied by the user in this conversation: 178 points, ordered as true stress
(MPa) and equivalent true plastic strain (dimensionless), following Abaqus input
conventions. Its independent experimental origin, a paper page containing those
numeric rows, and a raw tensile-test archive have not been established by this
audit. It must not be labelled newly recovered paper experimental data.

The audit compares **every row and both columns** against the existing
`PLASTIC_DATA` in `reference_setup_numeric_source.py`, using exact decimal
numbers and no tolerance. Counts agree at 178, but the numbers are **not exactly
identical**:

| Quantity | Exact matching values | Differing values | Largest absolute difference |
|---|---:|---:|---:|
| True stress | 0 | 178 | 0.000499 MPa (row 133) |
| Plastic strain | 8 | 170 | 0.00000005 (rows 50, 128, 155) |

Every one of the 356 supplied numbers matches six-significant-digit formatting
of the existing source's binary64 values. This explains the small differences;
it is a formatting-consistency observation, not exact numeric identity or
independent experimental validation. The existing higher-precision default table
and previous run inputs remain unchanged. No solver was rerun for this audit.

`user_plastic_table_audit.json` includes all 178 differing rows, source-file and
numeric SHA-256 hashes, exact comparison details and the six-digit check. Numeric
hashes use ASCII `stress,plastic_strain` rows, one LF each, fixed-point decimal
notation, insignificant trailing fractional zeros removed and zero normalized
to `0`. The supplied numeric SHA-256 is
`3dda8eda3646153ce0d0738c99eaa684d912e3c6a3b8a868a5c2596fc31ef024`;
the existing source's numeric SHA-256 is
`b4fa8c2dfe2e342ad788afcc4f11e2c8c09ea272e5d1e3a3d47f88fc2afb6982`.

Reproduce the audit from the repository root:

```sh
/workspace/fea-venv/bin/python scripts/audit_plastic_table.py \
  validation/materials/user_plastic_table_2026-10-26.inp \
  --output validation/materials/user_plastic_table_audit.json
```

Both model generators now accept this block directly through `--plastic`, as
well as the existing literal Python table. For example, to **generate** a native
deck with the supplied rounded values for a future separate run:

```sh
/workspace/fea-venv/bin/python run_radioss.py \
  --config configs/dong2019_brick_damage.json \
  --plastic validation/materials/user_plastic_table_2026-10-26.inp \
  --output results/opencourant/user_plastic_thin \
  --duration-seconds .01 --drop-ccx-viscosity
```

This does not launch the solver. Duration, AMS, energy and contact qualification
requirements still apply to any subsequent explicit simulation. The input
reader rejects duplicate or missing plastic curves, unsupported rate/field
options, extra columns and invalid hardening values; it executes no Python or
Abaqus source. It performs no engineering-to-true conversion, so data must
already use the documented true-stress/plastic-strain interpretation.

A hardening table does not supply the missing fracture-locus coordinates or
establish lattice calibration and agreement with the paper.
