"""Create a minimal CalculiX C3D6 elastic-plastic mesh baseline from Abaqus INP.
Not a replacement for Abaqus ductile damage, self-contact, or tied steel platens.
"""
import argparse
from pathlib import Path

def sections(path):
    section=None; lines=[]
    for raw in Path(path).read_text(errors='replace').splitlines():
        s=raw.strip()
        if not s or s.startswith('**'): continue
        if s.startswith('*'):
            if section is not None: yield section,lines
            section=s.upper();lines=[]
        elif section is not None: lines.append(s)
    if section is not None: yield section,lines

def read_mesh(path):
    nodes={}; elements={}
    for header,rows in sections(path):
        if header.startswith('*NODE'):
            for row in rows:
                p=[v.strip() for v in row.split(',') if v.strip()]
                if len(p)>=4: nodes[int(p[0])]=tuple(map(float,p[1:4]))
        elif header.startswith('*ELEMENT') and 'TYPE=C3D6' in header.replace(' ',''):
            for row in rows:
                p=[v.strip() for v in row.split(',') if v.strip()]
                if len(p)==7: elements[int(p[0])]=tuple(map(int,p[1:]))
    if not nodes or not elements: raise ValueError('Need *NODE and *ELEMENT, TYPE=C3D6 in mesh INP')
    return nodes,elements

def read_plastic(path):
    for header,rows in sections(path):
        if header.startswith('*PLASTIC'):
            return [r for r in rows if len(r.split(','))>=2]
    raise ValueError('No *PLASTIC table in input file')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True)
    ap.add_argument('--output',default='core_baseline.inp')
    ap.add_argument('--plastic-inp',help='Abaqus-generated INP with exact *PLASTIC table')
    ap.add_argument('--height',type=float,help='Expected specimen height, for checking only')
    ap.add_argument('--displacement',type=float,required=True,help='Top prescribed negative U2 magnitude, mm')
    ap.add_argument('--tolerance',type=float,default=1e-4)
    args=ap.parse_args()
    nodes,elems=read_mesh(args.input)
    ys=[p[1] for p in nodes.values()]; y0=min(ys);y1=max(ys)
    if args.height is not None and abs(y1-y0-args.height)>args.tolerance:
        raise ValueError('Height mismatch: mesh %.6f vs expected %.6f'%(y1-y0,args.height))
    bottom=[i for i,p in nodes.items() if abs(p[1]-y0)<=args.tolerance]
    top=[i for i,p in nodes.items() if abs(p[1]-y1)<=args.tolerance]
    if not bottom or not top: raise ValueError('Missing top or bottom nodes')
    lines=['** ELASTIC-PLASTIC BASELINE ONLY: NO DAMAGE / NO CONTACT / NO PLATENS',
           '*NODE']
    lines += ['%d, %.12g, %.12g, %.12g'%(i,*p) for i,p in sorted(nodes.items())]
    lines += ['*ELEMENT, TYPE=C3D6, ELSET=CORE']
    lines += ['%d, %s'%(i,', '.join(map(str,c))) for i,c in sorted(elems.items())]
    for name,ids in [('BOTTOM',bottom),('TOP',top)]:
        lines += ['*NSET, NSET='+name]
        lines += [', '.join(map(str,ids[k:k+12])) for k in range(0,len(ids),12)]
    lines += ['*MATERIAL, NAME=ALSI10MG','*ELASTIC','69000., 0.30','*DENSITY','2.67E-9']
    if args.plastic_inp:
        lines += ['*PLASTIC']+read_plastic(args.plastic_inp)
    else:
        lines += ['** WARNING: NO PLASTIC TABLE SUPPLIED; ELASTIC ONLY']
    lines += ['*SOLID SECTION, ELSET=CORE, MATERIAL=ALSI10MG',
              '*STEP, NLGEOM=YES, INC=1000','*STATIC','0.001, 1.0, 1e-8, 0.02',
              '** Simplified rigid-face displacement constraints, NOT platen tie',
              '*BOUNDARY','BOTTOM, 1, 3, 0.0','TOP, 1, 1, 0.0',
              'TOP, 2, 2, %.12g'%(-abs(args.displacement)),
              'TOP, 3, 3, 0.0','*NODE FILE','U, RF','*EL FILE','S, E, PE',
              '*END STEP']
    Path(args.output).write_text('\n'.join(lines)+'\n')
    print('Wrote',args.output,'nodes',len(nodes),'wedge elements',len(elems),
          'bottom nodes',len(bottom),'top nodes',len(top))
    print('WARNING: BASELINE ONLY; does not include damage, element deletion, contact or platen mechanics')

if __name__=='__main__': main()
