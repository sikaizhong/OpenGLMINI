#pragma once
#define NOMINMAX
#include<color.h>
#include<matrix.h>
#include<cmath>
#include<vector.h>
namespace CGL {
	class VertexData {
	public:
		vec3 m_local_pos;
		vec3 m_world_pos;
		vec3 m_world_normal;
		vec4 m_cpos;
		vec3 m_screen_pos;

		vec3 m_clipping_normal;

		vec4 m_color;
		vec2 m_tex;
		
		VertexData() = default;
		VertexData(const VertexData&) = default;

		static VertexData lerp(
			const VertexData& v0,
			const VertexData& v1,
			float frac)
		{
			//Linear interpolation
			VertexData result;
			
			result.m_cpos = (1.0f - frac) * v0.m_cpos + frac * v1.m_cpos;
			result.m_world_pos = (1.0f - frac) * v0.m_world_pos + frac * v1.m_world_pos;


			return result;
		}








	};
	struct FragmentData
	{
	public:
		vec3 m_pos;  //World space position
		vec3 m_nor;  //World space normal
		vec3 m_clipping_normal;
		vec2 m_tex;	//World space texture coordinate
		vec4 m_color;// object color
		ivec2 m_spos;//Screen space position
		//mat3 m_tbn;  //Tangent, bitangent, normal matrix;
	};


	inline vec3 reflect(const vec3& L, const vec3& N) {
		return 2 * dot(L, N) * N - L;
	}

	inline void basic_vs(VertexData& vertex, mat4& model, mat4& view, mat4& projection) {
		mat4 scale;
		scale(0, 0) = scale(1, 1) = scale(2, 2) = 1.8;
		float data[4] = { vertex.m_world_pos[0], vertex.m_world_pos[1], vertex.m_world_pos[2],1.0f };
		vertex.m_cpos = projection * view * model * scale *
			vec4(data);
	
		vertex.m_screen_pos[0] = vertex.m_cpos[0];
		vertex.m_screen_pos[1] = vertex.m_cpos[1];
		vertex.m_screen_pos[2] = vertex.m_cpos[2];
	}
	inline Color basic_fs(FragmentData&fragment) {
		return Color(fragment.m_color[0], fragment.m_color[1],
			fragment.m_color[2], fragment.m_color[3]);
	}

	inline void Phong_vs(VertexData& vertex, mat4& model, mat4& view, mat4& projection) {
		mat4 scale;
		scale(0, 0) = scale(1, 1) = scale(2, 2) =0.8f;
		float data[4] = { vertex.m_world_pos[0], vertex.m_world_pos[1], vertex.m_world_pos[2],1.0f };
		vertex.m_cpos = //projection *
			view * scale *
			vec4(data);
		vertex.m_cpos = projection * vertex.m_cpos; 
		for (int x = 0; x < 4; x++) {
			//vertex.m_cpos[x] = vertex.m_cpos[x] / vertex.m_cpos[3];
		}

		float data_n[4] = { vertex.m_world_normal[0], vertex.m_world_normal[1], vertex.m_world_normal[2],1.0f};

		vertex.m_clipping_normal = vec3((projection * view * model * scale * vec4(data_n)).data());
	}
	inline Color Phong_fs(FragmentData& fragment) {
		// ambient
		vec3 lightColor = vec3(1.0f, 1.0f, 1.0f);
		vec3 lightPos = vec3(1.2f, 1.0f, 2.0f);
		vec3 viewPos = vec3(0, 0, 1);
		vec3 objectColor = vec3(200, 200, 200);
		float ambientStrength = 0.2;
		vec3 ambient = ambientStrength * lightColor;
		vec3 Normal = fragment.m_nor;
		vec3 FragPos = fragment.m_pos;
		// diffuse 
		vec3 norm = normalize(Normal);
		vec3 lightDir = normalize(lightPos - FragPos);
		float diff = std::max(dot(norm, lightDir), 0.0f);
		vec3 diffuse = diff * lightColor;
		// specular
		float specularStrength = 0.5;
		vec3 viewDir = normalize(viewPos - FragPos);
		vec3 reflectDir = reflect(-lightDir, norm);
		float spec = pow(std::max(dot(viewDir, reflectDir), 0.0f), 32);
		vec3 specular = specularStrength * spec * lightColor;
		Color res;
		vec3 m = (ambient + diffuse + specular);

		res.r = std::min(255.0f,(m.x) * objectColor.x);
		res.g = std::min(255.0f, (m.y) * objectColor.y);
		res.b = std::min(255.0f, (m.z) * objectColor.z);
		return res;

	}
}