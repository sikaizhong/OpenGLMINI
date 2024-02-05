#pragma once
#include<vector.h>
#include<matrix.h>
#include<color.h>
#include<framebuffer.h>
#include<cglpipeline.h>
#include<program.h>
#include<mesh.h>
namespace CGL {
	class Rasterizer final {
	public:
		Rasterizer(int width, int height):m_width(width), m_height(height) {
			m_cur_framebuffer = std::make_shared<CFramebuffer>(width, height);
			init_default_shader();
		}
		void run();
		void init_default_shader();
		void cgl_clear_color(float r, float g, float b, float a) {
			color.b = r * 255;
			color.g = g * 255;
			color.r = b * 255;
		}
		void cgl_clear(int state = 0) {
			std::fill(m_cur_framebuffer->CDepthBuffer.begin(), m_cur_framebuffer->CDepthBuffer.end(), -10000000000.0f);
			for (int x = 0; x < m_cur_framebuffer->CColorBuffer.size() / 3; x++) {
				m_cur_framebuffer->CColorBuffer[3 * x + 0] = color.r;
				m_cur_framebuffer->CColorBuffer[3 * x + 1] = color.g;
				m_cur_framebuffer->CColorBuffer[3 * x + 2] = color.b;
			}
		}
		std::vector<VertexData> homogeneous_clipping(VertexData& v0, VertexData& v1, VertexData& v2, float near, float far);

		std::vector<VertexData> homogeneous_clipping_fast_aux(VertexData& v0, VertexData& v1, VertexData& v2, float near, float far);
		std::vector<VertexData> homogeneous_clipping_SutherlandHodgeman_aux(VertexData& v0, VertexData& v1, VertexData& v2, float near, float far);

		std::vector<VertexData> Rasterizer::clipingSutherlandHodgemanAux(
			const std::vector<VertexData>& polygon,
			const int& axis,
			const int& side);

		bool depth_test(int x, int y, float new_depth, int screen_width, int screen_height);


		void set_up_vertex_array(std::vector<float>& vertices, std::vector<unsigned>& triangles);
		void set_up_vertex_array(std::vector<float>& vertices,std::vector<float>&color, std::vector<unsigned>& triangles);
		void set_up_vertex_array(std::vector<vec3>& vertices, std::vector<int>& triangles);
		void set_up_vertex_array(Mesh* msh);

		void draw_by_index(int p_count = 0);
		void draw_triangle_barycentric_aux(VertexData& v0, VertexData& v1, VertexData& v2, int screen_width, int screen_height);

		void set_viewport(int x, int y, int height, int width);
		void cgl_look_at(const vec3& eye, const vec3& center, const vec3& up);
		void cgl_perspective(float fovy, float aspect, float near, float far);

		void cgl_orth(float left, float right, float bottom, float top, float near, float far);


		std::shared_ptr<CFramebuffer> m_cur_framebuffer;
		std::shared_ptr<Program> m_current_program;
		mat4 m_md_mat;
		mat4 m_v_mat;
		mat4 m_prj_mat;
		mat4 m_vp_mat;
	protected:
		std::vector<VertexData> m_vertex_buffer;
		std::vector<unsigned> m_index_buffer;
		std::vector<FragmentData> m_frag_buffer;
	private:
		int m_init_render_state = 0b00000001;// nothing is enabled;
		int m_width;
		int m_height;
		Color color;   
	};




}