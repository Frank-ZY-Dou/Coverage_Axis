import trimesh
import numpy as np

def sample_vertices_from_mesh(input_obj_path, output_obj_path, num_samples=3000):
    mesh = trimesh.load_mesh(input_obj_path)

    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError("Loaded mesh is not a valid trimesh object.")

    sampled_indices = np.random.choice(len(mesh.vertices), num_samples, replace=False)
    sampled_vertices = mesh.vertices[sampled_indices]
    sampled_normals = mesh.vertex_normals[sampled_indices]

    new_mesh = trimesh.Trimesh(vertices=sampled_vertices, faces=[], vertex_normals=sampled_normals)

    new_mesh.export(output_obj_path)
    print(f"Sampled mesh saved to {output_obj_path}")


input_obj = "./input/01Ants-12_mesh.obj"
output_obj = "./input/01Ants-12_mesh_ori_pc.obj"
sample_vertices_from_mesh(input_obj, output_obj)
