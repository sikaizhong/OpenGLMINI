#include<vector.h>
#include<mesh_io.h>
#include <iostream>
#include <sstream>
namespace CGL {
	void load_obj_file(std::string filename, Mesh* mesh) {
		auto& verts = mesh->m_verts;
		auto& norms = mesh->m_normals;
		auto& tex_coord = mesh->m_tex_coord;
		auto& facet_vrt = mesh->facet_vrt;
		auto& facet_tex = mesh->facet_tex;
		auto& facet_nrm = mesh->facet_nrm;
		std::ifstream in;
		in.open(filename, std::ifstream::in);
		if (in.fail()) return;
		std::string line;
		while (!in.eof()) {
			std::getline(in, line);
			std::istringstream iss(line.c_str());
			char trash;
			if (!line.compare(0, 2, "v ")) {
				iss >> trash;
				vec3 v;
				for (int i = 0; i < 3; i++) iss >> v[i];
				verts.push_back(v);
			}
			else if (!line.compare(0, 3, "vn ")) {
				iss >> trash >> trash;
				vec3 n;
				for (int i = 0; i < 3; i++) iss >> n[i];
				norms.push_back(normalize(n));
				//norms.push_back(n);

			}
			else if (!line.compare(0, 3, "vt ")) {
				iss >> trash >> trash;
				ivec2 uv;
				for (int i = 0; i < 2; i++) iss >> uv[i];
				tex_coord.push_back({ uv.x, 1 - uv.y });
			}
			else if (!line.compare(0, 2, "f ")) {
				int f, t, n;
				iss >> trash;
				int cnt = 0;
				while (iss >> f >> trash >> t >> trash >> n) {
					facet_vrt.push_back(--f);
					facet_tex.push_back(--t);
					facet_nrm.push_back(--n);
					cnt++;
				}
				if (3 != cnt) {
					std::cerr << "Error: the obj file is supposed to be triangulated" << std::endl;
					return;
				}
			}
		}

	}



}