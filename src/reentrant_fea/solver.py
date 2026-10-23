"""CalculiX nonlinear elastic-plastic baseline (no false damage/contact claims)."""
from pathlib import Path
from .plastic import validate_material
import subprocess

def write_deck(path,nodes,elements,top,bottom,c,plastic):
    """Rigid-grip idealization: every top node shares U2, every bottom node fixed.
    Not a deformable steel platen or bonded face-to-face contact implementation.
    """
    validate_material(c, plastic)
    path=Path(path)
    if len(plastic)<2: raise ValueError('At least two plastic points required')
    lines=['*HEADING','Independent ReEntrant WEDGE6 baseline','*NODE']
    lines.extend('%d, %.10g, %.10g, %.10g'%(i+1,*p) for i,p in enumerate(nodes))
    lines.append('*ELEMENT,TYPE=C3D6,ELSET=CORE')
    lines.extend('%d, %s'%(i+1,', '.join(str(int(v)+1) for v in el)) for i,el in enumerate(elements))
    lines += ['*NSET,NSET=TOP']
    lines.extend(str(int(n)+1)+',' for n in top)
    lines += ['*NSET,NSET=BOTTOM']
    lines.extend(str(int(n)+1)+',' for n in bottom)
    lines += ['*MATERIAL,NAME=ALSI10MG','*ELASTIC','%.10g, %.10g'%(c['young'],c['poisson']), '*PLASTIC']
    lines.extend('%.10g, %.10g'%(stress,strain) for stress,strain in plastic)
    lines += ['*SOLID SECTION,ELSET=CORE,MATERIAL=ALSI10MG', '',
              '*BOUNDARY','BOTTOM,1,3,0.0','TOP,1,1,0.0','TOP,3,3,0.0',
              '*STEP,NLGEOM=YES,INC=100000','*STATIC','0.001,1.0,1e-9,0.02',
              '*BOUNDARY','TOP,2,2,%.10g'%(-c['compression_strain']*(nodes[:,1].max()-nodes[:,1].min())),
              '*NODE FILE','U,RF','*EL FILE','S,E,PEEQ','*END STEP']
    path.write_text('\n'.join(lines)+'\n',encoding='utf8')
    return path

def run_ccx(deck,executable='ccx'):
    deck=Path(deck).resolve()
    return subprocess.run([executable,'-i',deck.stem],cwd=deck.parent,check=True)
