"""Experimental CalculiX deck with rigid bonded grips and optional self-contact.

No native ductile softening/deletion is provided by this deck. Enabling damage
without a verified solver extension raises an error rather than faking failure.
"""
from pathlib import Path
from .plastic import validate_material
import numpy as np

def write_advanced_deck(path, nodes, elements, top, bottom, config, plastic,
                        *, self_contact=False, ductile_damage=False):
    validate_material(config, plastic)
    if ductile_damage:
        raise NotImplementedError('Coupled ductile damage and element deletion require a verified constitutive solver extension; CalculiX *PLASTIC alone cannot provide them.')
    if self_contact:
        raise NotImplementedError('Same-surface self contact is not verified for CalculiX 2.20')
    nodes=np.asarray(nodes,dtype=float)
    elements=np.asarray(elements,dtype=int)
    top=np.asarray(top,dtype=int)
    bottom=np.asarray(bottom,dtype=int)
    if len(set(top).intersection(bottom)):
        raise ValueError('TOP and BOTTOM node sets overlap')
    if len(top)<1 or len(bottom)<1:
        raise ValueError('Missing end grip nodes')
    if not np.isfinite(nodes).all(): raise ValueError('Nonfinite mesh coordinates')
    if len(plastic)<2: raise ValueError('At least two plastic data points needed')
    ymin,ymax=float(nodes[:,1].min()),float(nodes[:,1].max())
    nt,nb=len(nodes)+1,len(nodes)+2
    rt,rb=len(nodes)+3,len(nodes)+4
    # Reference nodes are independent of core, controlled by rigid-body kinematics.
    lines=['*HEADING','ReEntrant experimental nonlinear bonded-grip model',
           '*NODE']
    lines += [f'{i+1}, {x:.10g}, {y:.10g}, {z:.10g}' for i,(x,y,z) in enumerate(nodes)]
    lines += [f'{nt}, {np.mean(nodes[top,0]):.10g}, {ymax:.10g}, {np.mean(nodes[top,2]):.10g}',
              f'{nb}, {np.mean(nodes[bottom,0]):.10g}, {ymin:.10g}, {np.mean(nodes[bottom,2]):.10g}',
              f'{rt}, 0., 0., 0.', f'{rb}, 0., 0., 0.',
              '*ELEMENT,TYPE=C3D6,ELSET=CORE']
    lines += [f'{i+1}, '+', '.join(str(int(j)+1) for j in el) for i,el in enumerate(elements)]
    def nset(name, ids):
        lines.append(f'*NSET,NSET={name}')
        ids=list(ids)
        for k in range(0,len(ids),12): lines.append(', '.join(str(int(v)+1) for v in ids[k:k+12]))
    nset('TOP_GRIP',top)
    nset('BOTTOM_GRIP',bottom)
    lines += ['*NSET,NSET=TOP_RP',str(nt),'*NSET,NSET=BOTTOM_RP',str(nb),
              '*MATERIAL,NAME=ALSI10MG','*ELASTIC',f"{config['young']}, {config['poisson']}",
              '*PLASTIC']
    lines += [f'{stress:.12g}, {strain:.12g}' for stress,strain in plastic]
    lines += ['*SOLID SECTION,ELSET=CORE,MATERIAL=ALSI10MG',
              '*RIGID BODY,NSET=TOP_GRIP,REF NODE='+str(nt)+',ROT NODE='+str(rt),
              '*RIGID BODY,NSET=BOTTOM_GRIP,REF NODE='+str(nb)+',ROT NODE='+str(rb)]
    lines += ['*BOUNDARY','BOTTOM_RP,1,3,0.0',f'{rt},1,3,0.',f'{rb},1,3,0.','TOP_RP,1,1,0.0',
              'TOP_RP,3,3,0.0',
              '*STEP,NLGEOM=YES,INC=100000','*STATIC',
              '0.0001,1.0,1e-6,0.01',
              '*BOUNDARY',f"TOP_RP,2,2,{-config['compression_strain']*(ymax-ymin):.10g}",
              '*NODE FILE','U,RF','*EL FILE','S,E,PEEQ',
              '*NODE PRINT,NSET=TOP_RP','U,RF',
              '*END STEP']
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text('\n'.join(lines)+'\n',encoding='utf8')
    return path
