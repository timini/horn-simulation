"""CAD eligibility checks for the axisymmetric modal aperture approximation."""
import gmsh
import numpy as np


class ModalGeometryError(ValueError):
    """The CAD volume is not eligible for the axisymmetric aperture model."""


def verify_axisymmetric_disk_horn(volumes, inlet, outlet):
    """Reject annular/offset/noncircular ports and non-axisymmetric interiors.

    Boolean overlap under two noncommensurate rotations checks the complete
    volume, including internal walls, rather than inferring shape from area.
    STEP splines are compared within 0.1 micrometre radial CAD tolerance;
    the overlap tolerance accounts for the corresponding displaced shell.
    """
    if len(volumes)!=1:
        raise ModalGeometryError('Modal aperture requires one connected axisymmetric volume')
    radii=[]
    for surfaces in (inlet,outlet):
        if len(surfaces)!=1:
            raise ModalGeometryError('Modal aperture requires single circular disk ports')
        curves=gmsh.model.getBoundary([(2,surfaces[0])],oriented=False)
        if len(curves)!=1:
            raise ModalGeometryError('Modal aperture requires circular disk ports without holes')
        center=gmsh.model.occ.getCenterOfMass(2,surfaces[0])
        area=gmsh.model.occ.getMass(2,surfaces[0])
        # STEP lofts can encode an exact circle as a rational B-spline.
        lower,upper=gmsh.model.getParametrizationBounds(*curves[0])
        points=np.asarray(gmsh.model.getValue(*curves[0],np.linspace(lower[0],upper[0],257))).reshape(-1,3)
        radius=np.sqrt(area/np.pi)
        radii.append(radius)
        if not np.allclose(np.linalg.norm(points[:,:2],axis=1),radius,rtol=1e-6,atol=1e-7):
            raise ModalGeometryError('Modal aperture requires centered circular disk ports')
        if np.linalg.norm(center[:2])>max(1e-7,radius*1e-6):
            raise ModalGeometryError('Modal aperture ports must be centered on the z axis')
    mass=gmsh.model.occ.getMass(*volumes[0])
    if not np.isfinite(mass) or mass<=0:raise ModalGeometryError('Invalid acoustic volume')
    for angle in (.713,1.231):
        first=gmsh.model.occ.copy(volumes);second=gmsh.model.occ.copy(volumes)
        gmsh.model.occ.rotate(second,0,0,0,0,0,1,angle)
        intersection,_=gmsh.model.occ.intersect(first,second,removeObject=True,removeTool=True)
        common=sum(gmsh.model.occ.getMass(*entity) for entity in intersection if entity[0]==3)
        gmsh.model.occ.remove(intersection,recursive=True)
        if not np.isfinite(common) or abs(common/mass-1)>max(1e-6,4e-7/min(radii)):
            raise ModalGeometryError('Modal aperture requires a rotationally invariant acoustic volume')
    gmsh.model.occ.synchronize()
