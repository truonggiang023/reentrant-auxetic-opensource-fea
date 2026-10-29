"""Streaming ASCII FRD result fields, including contact and state variables."""
from pathlib import Path


def fields(path):
    time = None
    name = None
    data = {}
    last = None
    with Path(path).open() as stream:
        for line in stream:
            if line.startswith('  100CL'):
                time = float(line.split()[2])
            elif line.startswith(' -4'):
                name=line.split()[1];data={};last=None
            elif name and line.startswith(' -1'):
                last=int(line[3:13]);values=line[13:].rstrip('\n')
                data[last]=[float(values[i:i+12]) for i in range(0,len(values),12) if values[i:i+12].strip()]
            elif name and line.startswith(' -2') and last is not None:
                values=line[13:].rstrip('\n')
                data[last].extend(float(values[i:i+12]) for i in range(0,len(values),12) if values[i:i+12].strip())
            elif name and line.startswith(' -3'):
                yield time,name,data
                name=None


def contact_envelope(path, committed=True):
    """Largest reported nodal contact overclosure and pressure (not Gauss data)."""
    penetration=pressure=0.;frames=0
    source=committed_fields(path) if committed else fields(path)
    for time,name,data in source:
        if name!='CONTACT':continue
        frames+=1
        for row in data.values():
            penetration=max(penetration,-row[0]);pressure=max(pressure,row[3])
    if not frames:raise ValueError('No CONTACT output in FRD')
    return {'contact_frames':frames,'max_reported_penetration_mm':penetration,
            'max_reported_pressure_MPa':pressure}


def committed_fields(path, status=None):
    """Exclude unconverged 'best solution' frames written after a solver failure.

    FRD timestamps have limited precision. Match them to accepted TOTAL times
    in the status file, not simply the latest available FRD dataset.
    """
    import bisect
    status=Path(status) if status else Path(path).with_suffix('.sta')
    times=[]
    for line in status.read_text().splitlines():
        p=line.split()
        if len(p)==7 and p[0].isdigit() and p[2].isdigit():times.append(float(p[4]))
    times=sorted(set(times))
    for time,name,data in fields(path):
        i=bisect.bisect_left(times,time)
        nearby=times[max(0,i-1):i+1]
        if any(abs(t-time)<=max(1e-10,abs(t)*6e-6) for t in nearby):
            yield time,name,data


def projected_overlap(nodes, elements, displacement, tolerance=1e-8):
    """Overlap of front-plane element projections, for plane-strain diagnostics.

    This is a geometric postprocessor, not a contact force calculation. It can
    flag missed wall interactions; it cannot prove complete 3-D contact validity.
    """
    import numpy as np
    from shapely.geometry import Polygon
    from shapely.strtree import STRtree
    zmin=nodes[:,2].min();polys=[];ids=[];invalid=[]
    for eid,el in enumerate(elements,1):
        face=el[:len(el)//2]
        if not np.all(abs(nodes[face,2]-zmin)<1e-7):continue
        xy=np.array([nodes[v,:2]+np.array(displacement[int(v)+1][:2]) for v in face])
        area=.5*np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1))
        poly=Polygon(xy)
        if area<=0 or not poly.is_valid:
            invalid.append(eid);continue
        polys.append(poly);ids.append(eid)
    tree=STRtree(polys);pairs=0;maximum=0.;total=0.
    for i,poly in enumerate(polys):
        for j in tree.query(poly):
            if j<=i:continue
            area=poly.intersection(polys[int(j)]).area
            if area>tolerance:pairs+=1;maximum=max(maximum,area);total+=area
    return {'front_plane_elements':len(polys),'invalid_projected_elements':invalid,
            'overlap_pairs':pairs,'max_overlap_mm2':maximum,'sum_pair_overlap_mm2':total}


def read_core_deck(path):
    """Read this project's single-step static paper deck for postprocessing."""
    import numpy as np
    text=Path(path).read_text();nodes={};elements=[];mode='';applied=None;period=None;steps=0
    if '*DYNAMIC' in text:raise ValueError('Dynamic strain mapping requires its actual amplitude history')
    for line in text.splitlines():
        if line.startswith('*STEP'):steps+=1
        if line=='*NODE':mode='nodes';continue
        if line.startswith('*ELEMENT'):
            mode='core' if 'ELSET=CORE' in line else '';continue
        if line.startswith('*STATIC'):mode='static';continue
        if line.startswith('*USER MATERIAL') and int(line.split('CONSTANTS=')[1])>21:
            raise ValueError('Unsafe oversized USER MATERIAL table; rerun through the sidecar interface')
        if line.startswith('*'):mode='';continue
        if line.startswith(('TOP_PLATEN_NODES,2,2,', 'RP_TOP_PLATEN,2,2,')):applied=-float(line.split(',')[-1])
        values=line.split(',')
        if mode=='nodes':nodes[int(values[0])]=np.array(list(map(float,values[1:4])))
        elif mode=='core':elements.append(list(map(int,values[1:])))
        elif mode=='static':period=float(values[1]);mode=''
    if steps!=1 or not period or not applied or not elements:raise ValueError('Unsupported or incomplete static paper deck')
    if sorted(nodes)!=list(range(1,len(nodes)+1)):raise ValueError('Expected contiguous generated node IDs')
    xyz=np.array([nodes[i] for i in range(1,len(nodes)+1)])
    e=np.array(elements,dtype=int)-1;core=xyz[np.unique(e)]
    height=float(np.ptp(core[:,1]));width=float(np.ptp(core[:,0]));depth=float(np.ptp(core[:,2]))
    return xyz,e,{'height_mm':height,'width_mm':width,'depth_mm':depth,
                  'period':period,'end_nominal_strain':applied/height,
                  'core_nodes':len(np.unique(e)),'core_elements':len(e),
                  'native_material':'*USER MATERIAL' in text}
