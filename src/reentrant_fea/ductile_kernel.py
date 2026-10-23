"""Standalone small-strain 3D J2 + scalar ductile damage constitutive kernel.

NOT wired into CalculiX. Deletion flag is a request to a host FE solver, not
physical removal by this module. Voigt order [xx,yy,zz,xy,yz,xz], engineering
shear strains; stresses in MPa and lengths in mm.
"""
from dataclasses import dataclass, replace
import numpy as np

@dataclass
class State:
    plastic_strain: np.ndarray
    peeq: float = 0.0
    omega: float = 0.0
    damage: float = 0.0
    post_u: float = 0.0
    deleted: bool = False

    @classmethod
    def virgin(cls):
        return cls(np.zeros(6))


def elastic_matrix(E, nu):
    if E <= 0 or not (-1.0 < nu < 0.5):
        raise ValueError('Invalid isotropic elasticity')
    lam = E*nu/((1+nu)*(1-2*nu))
    mu = E/(2*(1+nu))
    C = np.zeros((6,6))
    C[:3,:3] = lam
    for i in range(3): C[i,i] += 2*mu
    for i in range(3,6): C[i,i] = mu
    return C


def deviator(s):
    out = np.array(s, dtype=float, copy=True)
    out[:3] -= np.mean(out[:3])
    return out


def j2_equivalent(s):
    d = deviator(s)
    return np.sqrt(max(0., 1.5*(np.dot(d[:3],d[:3])+2*np.dot(d[3:],d[3:]))))


def interpolate_locus(eta, locus):
    x = np.asarray(locus,dtype=float)
    if x.ndim != 2 or x.shape[1] != 2 or len(x) < 1 or not np.isfinite(x).all() or np.any(x[:,1]<=0):
        raise ValueError('Invalid fracture locus')
    x = x[np.argsort(x[:,0])]
    if len(x)>1 and np.any(np.diff(x[:,0])<=0): raise ValueError('Duplicate triaxiality')
    return float(np.interp(eta,x[:,0],x[:,1]))


def flow_stress(peeq, hardening):
    table = np.asarray(hardening,dtype=float)
    if table.ndim != 2 or table.shape[1]!=2 or len(table)<2 or not np.isfinite(table).all():
        raise ValueError('Invalid hardening table')
    if np.any(np.diff(table[:,1])<=0) or np.any(np.diff(table[:,0])<0) or table[0,1]<0 or table[0,0]<=0:
        raise ValueError('Hardening requires positive increasing strain and nondecreasing stress')
    return float(np.interp(peeq,table[:,1],table[:,0]))


def integrate(strain, previous, E, nu, hardening, locus, characteristic_length, failure_displacement,
              residual_stiffness=1e-6):
    """Total small strain update. Return (effective stress, damaged stress, trial state).

    The returned state must be COMMITTED only after global equilibrium converges.
    Stress degradation is coupled at the material-point level; host FEM must use
    a consistent/verified tangent and handle element deletion explicitly.
    """
    if characteristic_length<=0 or failure_displacement<=0 or not 0<residual_stiffness<=1:
        raise ValueError('Invalid length, displacement or residual stiffness')
    if previous.deleted:
        return np.zeros(6), np.zeros(6), replace(previous,plastic_strain=previous.plastic_strain.copy())
    strain = np.asarray(strain,dtype=float)
    if strain.shape != (6,) or not np.isfinite(strain).all(): raise ValueError('Strain must be finite Voigt-6')
    C=elastic_matrix(E,nu)
    G=E/(2*(1+nu))
    trial=C@(strain-previous.plastic_strain)
    qtrial=j2_equivalent(trial)
    yield0=flow_stress(previous.peeq,hardening)
    dp=0.0
    if qtrial > yield0+1e-10:
        lo,hi=0.,qtrial/(3*G)
        for _ in range(75):
            mid=(lo+hi)/2
            f=qtrial-3*G*mid-flow_stress(previous.peeq+mid,hardening)
            if f>0: lo=mid
            else: hi=mid
        dp=(lo+hi)/2
    if dp:
        dev=deviator(trial)
        # Engineering plastic shear strain = 2 * tensor plastic shear strain.
        dep=np.zeros(6)
        dep[:3]=1.5*dp*dev[:3]/qtrial
        dep[3:]=3*dp*dev[3:]/qtrial
        ep=previous.plastic_strain+dep
    else:
        ep=previous.plastic_strain.copy()
    effective=C@(strain-ep)
    q=j2_equivalent(effective)
    eta=float(np.mean(effective[:3])/q) if q>1e-12 else 0.0
    omega_increment=dp/interpolate_locus(eta,locus)
    old=previous.omega
    omega=old+omega_increment
    post_fraction=0. if omega_increment==0 else np.clip((omega-1.)/omega_increment,0.,1.)
    post_u=previous.post_u+characteristic_length*dp*post_fraction
    damage=float(np.clip(post_u/failure_displacement,0.,1.))
    deleted=damage>=1.
    state=State(ep,previous.peeq+dp,omega,damage,post_u,deleted)
    stress=(max(residual_stiffness,1-damage))*effective
    return effective, stress, state
