"""Generate original Gmsh geometry and public-schema inputs for a pinned solver."""
from pathlib import Path
import json
import numpy as np


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def case_definitions(protocol):
    for geometry in protocol['geometries']:
        for size in protocol['mesh_sizes_m']:
            yield f"{geometry['name']}-{size:g}", geometry, size


def frequencies(protocol):
    return np.geomspace(protocol['frequency_min_hz'], protocol['frequency_max_hz'],
                        protocol['frequency_count'])


def prepare_cases(root, protocol):
    import gmsh
    for name, geometry, size in case_definitions(protocol):
        out = root / name
        out.mkdir()
        radius = geometry['throat_radius_m']
        mouth = geometry['mouth_radius_m']
        length = geometry['length_m']
        gmsh.initialize()
        try:
            gmsh.model.add(name)
            if radius == mouth:
                volume = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, length, radius)
            else:
                volume = gmsh.model.occ.addCone(0, 0, 0, 0, 0, length, radius, mouth)
            gmsh.model.occ.synchronize()
            gmsh.model.addPhysicalGroup(3, [volume], 1)
            gmsh.model.setPhysicalName(3, 1, 'air')
            tags = {2: [], 3: [], 4: []}
            for _, tag in gmsh.model.getEntities(2):
                box = gmsh.model.getBoundingBox(2, tag)
                if np.allclose([box[2], box[5]], 0, atol=1e-6):
                    group = 2
                elif np.allclose([box[2], box[5]], length, atol=1e-6):
                    group = 3
                else:
                    group = 4
                tags[group].append(tag)
            for tag, entities in tags.items():
                gmsh.model.addPhysicalGroup(2, entities, tag)
                gmsh.model.setPhysicalName(2, tag, {2: 'inlet', 3: 'mouth', 4: 'walls'}[tag])
            gmsh.write(str(out / 'horn.step'))
            gmsh.option.setNumber('Mesh.MeshSizeMin', size)
            gmsh.option.setNumber('Mesh.MeshSizeMax', size)
            gmsh.option.setNumber('Mesh.MshFileVersion', 4.1)
            gmsh.option.setNumber('Mesh.Binary', 0)
            gmsh.model.mesh.generate(3)
            gmsh.write(str(out / 'horn.msh'))
        finally:
            gmsh.finalize()
        mesh_id, region_id = 'mesh:horn', 'region:horn'
        boundaries = [dict(
            id=n, name=n, kind=k, region_id=region_id, parameters={},
            group=dict(mesh_id=mesh_id, dimension=2, name=n, tag=t))
            for n, k, t in [('inlet', 'moving', 2),
                            ('mouth', 'plane_wave_tube_termination', 3), ('walls', 'rigid', 4)]]
        driver = {**protocol['driver'], 'motion_axis': [0, 0, 1], 'motion_profile': 'rigid_translation'}
        system = dict(
            id='system:reference', name='Single-horn numerical reference', model_version=1,
            metadata={}, interfaces=[],
            meshes=[dict(id=mesh_id, name='horn', file='horn.msh', purpose='fem_volume',
                         scale_to_m=1., translation_m=[0, 0, 0])],
            regions=[dict(id=region_id, name='horn', kind='bounded_air', mesh_ids=[mesh_id],
                          density_kg_per_m3=protocol['air_density_kg_m3'],
                          sound_speed_m_per_s=protocol['sound_speed_m_s'],
                          volume_groups=[dict(mesh_id=mesh_id, dimension=3, name='air', tag=1)],
                          loss_model=dict(bulk_loss_factor=0.))],
            components=[dict(id='driver', name='Synthetic characterized piston',
                             kind='electrodynamic_transducer', boundary_ids=['inlet'], parameters=driver)],
            excitation_ports=[dict(id='voltage', name='Native reference voltage',
                                   kind='voltage', component_id='driver')], boundaries=boundaries)
        project = dict(schema_version=9, physical_system=system, symmetry='off',
                       stitch_exterior_meshes=False, imported_meshes=[], observation_planes=[],
                       component_channel_by_id={'driver': 'main'},
                       project_preferences=dict(freq_min_hz=protocol['frequency_min_hz'],
                           freq_max_hz=protocol['frequency_max_hz'], freq_count=protocol['frequency_count'],
                           polar_angle_step_deg=5., polar_observation_distance_m=1.,
                           spherical_sampling_enabled=False))
        write_json(out / 'project.blab.json', project)
        write_json(out / 'request.json', dict(schema_version=1,
            frequencies_hz=frequencies(protocol).tolist(), include_project_observations=False,
            retain=['fem_nodal_pressure']))
        write_json(out / 'definition.json', dict(geometry=geometry, mesh_size_m=size,
            mass_convention='synthetic dry diaphragm, no rear load', driver=protocol['driver']))
