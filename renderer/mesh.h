#pragma once 
#include<string>
#include<vector>
#include<vector.h>
namespace CGL {
	class Mesh {
	public:
		std::vector<vec3> m_verts;
		std::vector<vec3> m_normals;
		std::vector<ivec2> m_tex_coord;
		std::vector<unsigned> m_facets;
		std::vector<int> facet_vrt{};
		std::vector<int> facet_tex{};  // per-triangle indices in the above arrays
		std::vector<int> facet_nrm{};
	};



}