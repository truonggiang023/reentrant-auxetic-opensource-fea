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


def validate_material(config, plastic):
    """Validate rate-independent isotropic J2 input (true stress, plastic strain)."""
    import numpy as np
    E, nu = float(config['young']), float(config['poisson'])
    if not np.isfinite([E, nu]).all() or E <= 0 or not -1 < nu < .5:
        raise ValueError('Invalid isotropic elasticity')
    table = np.asarray(plastic, dtype=float)
    if (table.ndim != 2 or table.shape[1] != 2 or len(table) < 2
            or not np.isfinite(table).all() or table[0, 1] != 0
            or np.any(table[:, 0] <= 0) or np.any(np.diff(table[:, 1]) <= 0)
            or np.any(np.diff(table[:, 0]) < 0)):
        raise ValueError('Plastic table requires positive nondecreasing true stress and increasing plastic strain starting at zero; softening is not ductile damage')
    strain = float(config['compression_strain'])
    if not np.isfinite(strain) or not 0 < strain < 1:
        raise ValueError('compression_strain must be between zero and one')
    return table
