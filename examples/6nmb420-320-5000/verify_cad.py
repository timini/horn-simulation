"""Check a physical STEP export in millimetres using the CAD kernel."""
import json,sys
from pathlib import Path
import gmsh
import numpy as np
source=Path(sys.argv[1]);destination=Path(sys.argv[2])
g= json.loads(Path(sys.argv[3]).read_text())
gmsh.initialize();gmsh.option.setNumber('General.Terminal',0)
try:
 gmsh.option.setString('Geometry.OCCTargetUnit','MM')
 gmsh.model.occ.importShapes(str(source));gmsh.model.occ.synchronize()
 volumes=gmsh.model.getEntities(3);assert len(volumes)==1
 faces=gmsh.model.getBoundary(volumes,oriented=False)
 ports={}
 for dim,tag in faces:
  box=gmsh.model.getBoundingBox(dim,tag)
  if abs(box[5]-box[2])<1e-4:
   z=gmsh.model.occ.getCenterOfMass(dim,tag)[2]
   ports['throat' if abs(z)<1e-4 else 'mouth']={'z_mm':z,'diameter_mm':float(2*np.sqrt(gmsh.model.occ.getMass(dim,tag)/np.pi))}
 assert set(ports)=={'throat','mouth'}
 np.testing.assert_allclose(ports['mouth']['z_mm'],1000*g['length'],rtol=1e-8,atol=1e-5)
 for name in ports:np.testing.assert_allclose(ports[name]['diameter_mm'],2000*g[name+'_radius'],rtol=1e-6,atol=1e-5)
 volume=gmsh.model.occ.getMass(*volumes[0])*1e-6
 destination.write_text(json.dumps({'passed':True,'import_unit':'mm','ports':ports,'acoustic_volume_litres':volume,'scope':'Physical STEP scale and circular port dimensions; not manufacturing qualification'},indent=2)+'\n')
finally:gmsh.finalize()
