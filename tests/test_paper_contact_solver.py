import numpy as np
import pytest
from reentrant_fea.paper_contact_solver import write_paper_contact_deck


def test_paper_deck_has_sliding_platens_with_outer_face_reference_points(tmp_path):
    # Two triangular prisms form a rectangular slab, with Y-normal end faces.
    nodes=np.array([
        [0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],
        [0.,0.,2.],[1.,0.,2.],[1.,1.,2.],[0.,1.,2.]])
    els=np.array([[0,1,2,4,5,6],[0,2,3,4,6,7]])
    cfg={"young":69000.,"poisson":.3,"compression_strain":.65}
    path=tmp_path/"candidate.inp"
    info=write_paper_contact_deck(path,nodes,els,cfg,[(225.,0.),(300.,.02)])
    s=path.read_text()
    assert "*ELEMENT,TYPE=C3D6,ELSET=CORE" in s
    assert "*ELEMENT,TYPE=C3D8,ELSET=PLATENS" in s
    assert "CORE_TOP,TOP_PLATEN_CONTACT" in s
    assert "CORE_BOTTOM,BOTTOM_PLATEN_CONTACT" in s
    assert "RP_TOP_PLATEN,2,2,-0.65" in s
    assert "CORE_PLANE_STRAIN,3,3,0." in s
    assert "*RIGID BODY,NSET=TOP_PLATEN_NODES,REF NODE=26,ROT NODE=28" in s
    assert "*NODE PRINT,NSET=RP_TOP_PLATEN\nU,RF" in s
    assert info["nominal_area"] == 2.
    assert "*DAMAGE" not in s  # no false damage claim
    assert info["top_faces"]>0 and info["bottom_faces"]>0


def test_optional_self_contact_is_explicit(tmp_path):
    nodes=np.array([
        [0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],
        [0.,0.,2.],[1.,0.,2.],[1.,1.,2.],[0.,1.,2.]])
    els=np.array([[0,1,2,4,5,6],[0,2,3,4,6,7]])
    config={"young":69000.,"poisson":.3,"compression_strain":.65}
    out=tmp_path/"a.inp"
    info=write_paper_contact_deck(out,nodes,els,config,
                                 [(225.,0.),(300.,.02)],self_contact=True)
    assert "CORE_LATERAL,CORE_LATERAL" not in out.read_text()
    assert info["wall_contact_pairs"] == 0  # A single solid slab has no wall gap.


def test_frictionless_omits_friction_card(tmp_path):
    nodes=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                    [0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
    els=np.array([[0,1,2,4,5,6],[0,2,3,4,6,7]])
    cfg={'young':69000,'poisson':.3,'compression_strain':.12,
         'friction':0.,'stick_slope':0.}
    deck=tmp_path/'frictionless.inp'
    write_paper_contact_deck(deck,nodes,els,cfg,[(225,0),(300,.2)])
    assert '*FRICTION' not in deck.read_text()
