#pragma once 
#include<memory>
#include<shader.h>
namespace CGL {
	class Program final {
	public:
		void compile(std::shared_ptr<VertexShader>& vs, std::shared_ptr<FragmentShader>& fs) {
			m_vertex_shader = vs;
			m_fragment_shader = fs;
		}
	public:
		std::shared_ptr<VertexShader> m_vertex_shader;
		std::shared_ptr<FragmentShader> m_fragment_shader;
	};
}