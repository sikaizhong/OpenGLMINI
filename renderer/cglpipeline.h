#pragma once
#include<vector.h>
#include<color.h>
#include<matrix.h>
namespace CGL {
	struct VertexData {
		vec3 m_local_pos;
		vec3 m_world_pos;
		vec3 m_world_normal;
		vec4 m_color;
		vec2 m_tex;
		vec4 m_clipping_pos;
		vec3 m_screen_pos;
	};
	struct FragmentData
	{
	public:
		vec3 m_pos;  //World space position
		vec3 m_nor;  //World space normal
		vec2 m_tex;	//World space texture coordinate
		vec4 m_color;// object color
		ivec2 m_spos;//Screen space position
		//mat3 m_tbn;  //Tangent, bitangent, normal matrix;
	};




	inline void basic_vs(VertexData& vertex, mat4& model, mat4& view, mat4& projection) {
		mat4 scale;
		scale(0, 0) = scale(1, 1) = scale(2, 2) = 1.0;
		float data[4] = { vertex.m_world_pos[0], vertex.m_world_pos[1], vertex.m_world_pos[2],1.0f };
		vertex.m_clipping_pos = projection * view * model * scale *
			vec4(data);
		vertex.m_screen_pos[0] = vertex.m_clipping_pos[0];
		vertex.m_screen_pos[1] = vertex.m_clipping_pos[1];
		vertex.m_screen_pos[2] = vertex.m_clipping_pos[2];
	}
	inline Color basic_fs(FragmentData&fragment) {
		return Color(fragment.m_color[0], fragment.m_color[1],
			fragment.m_color[2], fragment.m_color[3]);
	}
}