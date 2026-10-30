"""Read the original plastic table without importing/executing Abaqus."""
import ast
import re
from pathlib import Path

def from_reference(path):
    s=Path(path).read_text(encoding='utf-8-sig')
    m=re.search(r'\bPLASTIC_DATA\s*=\s*\(',s)
    if not m: raise ValueError('PLASTIC_DATA table not found')
    start=s.index('(',m.start())
    depth=0
    end=None
    for i in range(start,len(s)):
        if s[i]=='(': depth+=1
        elif s[i]==')':
            depth-=1
            if depth==0:
                end=i+1
                break
    if end is None: raise ValueError('Unclosed PLASTIC_DATA')
    values=ast.literal_eval(s[start:end])
    if any(len(p)!=2 for p in values): raise ValueError('Invalid plastic table')
    if any(values[i+1][1]<=values[i][1] for i in range(len(values)-1)):
        raise ValueError('Plastic strains must increase')
    return values


def from_abaqus_plastic(path):
    """Read one rate-independent *Plastic block as true stress/plastic strain.

    Only two-column isotropic data are supported. Temperature, field variables,
    multiple tables and rate-dependent curves need explicit separate support;
    they must not be silently collapsed into a different constitutive law.
    No Python or Abaqus code is executed.
    """
    lines = Path(path).read_text(encoding='utf-8-sig').splitlines()
    headers = [i for i, line in enumerate(lines)
               if re.match(r'^\s*\*Plastic(?:\s|,|$)', line, re.I)]
    if len(headers) != 1:
        raise ValueError('Require exactly one *Plastic block; missing or duplicate curves are unsupported')
    options = [field.strip().lower() for field in lines[headers[0]].strip().split(',')[1:]]
    if any(option not in ('hardening=isotropic',) for option in options):
        raise ValueError('Only a rate-independent isotropic *Plastic block is supported')
    values = []
    for line in lines[headers[0]+1:]:
        line = line.strip()
        if not line or line.startswith('**'):
            continue
        if line.startswith('*'):
            break
        fields = [field.strip() for field in line.split(',')]
        if len(fields) != 2:
            raise ValueError('Each *Plastic row requires true stress and plastic strain only')
        try:
            row = tuple(float(field.replace('D', 'E').replace('d', 'e')) for field in fields)
        except ValueError as error:
            raise ValueError('Nonnumeric *Plastic row') from error
        values.append(row)
    validate_plastic_table(values)
    return tuple(values)


def load_plastic(path):
    """Dispatch to literal Python table or Abaqus data without executing either."""
    path = Path(path)
    text = path.read_text(encoding='utf-8-sig')
    if path.suffix.lower() == '.inp' or re.search(r'^\s*\*Plastic(?:\s|,|$)', text, re.I | re.M):
        return from_abaqus_plastic(path)
    return from_reference(path)


def validate_plastic_table(plastic):
    """Require a valid true-stress/equivalent-plastic-strain hardening table."""
    import numpy as np
    table = np.asarray(plastic, dtype=float)
    if (table.ndim != 2 or table.shape[1] != 2 or len(table) < 2
            or not np.isfinite(table).all() or table[0, 1] != 0
            or np.any(table[:, 0] <= 0) or np.any(np.diff(table[:, 1]) <= 0)
            or np.any(np.diff(table[:, 0]) < 0)):
        raise ValueError('Plastic table requires positive nondecreasing true stress and increasing plastic strain starting at zero; softening is not ductile damage')
    return table


def validate_material(config, plastic):
    """Validate rate-independent isotropic J2 input (true stress, plastic strain)."""
    import numpy as np
    E, nu = float(config['young']), float(config['poisson'])
    if not np.isfinite([E, nu]).all() or E <= 0 or not -1 < nu < .5:
        raise ValueError('Invalid isotropic elasticity')
    table = validate_plastic_table(plastic)
    strain = float(config['compression_strain'])
    if not np.isfinite(strain) or not 0 < strain < 1:
        raise ValueError('compression_strain must be between zero and one')
    return table
