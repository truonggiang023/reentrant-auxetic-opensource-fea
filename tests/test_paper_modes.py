import numpy as np
import pytest
from reentrant_fea.paper_modes import write_paper_deck


def _mesh():
    nodes = np.array([
        [0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
        [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.],
    ])
    elems = np.array([[0, 1, 2, 4, 5, 6], [0, 2, 3, 4, 6, 7]])
    return nodes, elems, np.array([2, 3, 6, 7]), np.array([0, 1, 4, 5])


def test_paper_mode_has_three_snapshots_contact_and_out_of_plane(tmp_path):
    nodes, elems, top, bottom = _mesh()
    config = {"young": 69000., "poisson": 0.3, "mesh_size": 0.1}
    p = write_paper_deck(
        tmp_path / "reentrant.inp", nodes, elems, top, bottom,
        config, [(224.67, 0.0), (351.15, 0.03814)]
    )
    deck = p.read_text()
    assert deck.count("*STEP,NLGEOM=YES") == 3
    assert deck.count("*CONTACT PAIR") == 2
    assert "*SOLID SECTION,ELSET=TOP_PLATEN,MATERIAL=STEEL" in deck
    assert "CORE_OUT_OF_PLANE,3,3,0." in deck
    assert "TOP_PLATEN_NODES,2,2,-0.12" in deck
    assert "TOP_PLATEN_NODES,2,2,-0.35" in deck
    assert "TOP_PLATEN_NODES,2,2,-0.65" in deck


def test_damage_is_never_silently_faked(tmp_path):
    nodes, elems, top, bottom = _mesh()
    with pytest.raises(NotImplementedError):
        write_paper_deck(
            tmp_path / "reentrant.inp", nodes, elems, top, bottom,
            {"young": 69000., "poisson": 0.3},
            [(224.67, 0.0), (351.15, 0.03814)], damage=True
        )
