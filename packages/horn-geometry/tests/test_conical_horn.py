"""Validate exported conical STEP volume and boundary areas using gmsh/OCC."""
import gmsh
import numpy as np
import pytest
from horn_geometry.generator import create_conical_horn


def test_conical_export_geometry(tmp_path):
    r,R,length=.025,.1,.2
    step=create_conical_horn(r,R,length,tmp_path/'horn.step')
    gmsh.initialize()
    try:
        gmsh.model.occ.importShapes(str(step));gmsh.model.occ.synchronize()
        volumes=gmsh.model.getEntities(3)
        assert len(volumes)==1
        volume=gmsh.model.occ.getMass(*volumes[0])
        assert volume==pytest.approx(np.pi*length*(R*R+R*r+r*r)/3,rel=1e-3)
        planes=[]
        for dim,tag in gmsh.model.getEntities(2):
            x,y,z=gmsh.model.occ.getCenterOfMass(dim,tag)
            if abs(z)<1e-7 or abs(z-length)<1e-7:
                planes.append((z,gmsh.model.occ.getMass(dim,tag)))
        planes.sort()
        assert len(planes)==2
        assert planes[0][1]==pytest.approx(np.pi*r*r,rel=1e-3)
        assert planes[1][1]==pytest.approx(np.pi*R*R,rel=1e-3)
    finally:
        gmsh.finalize()
