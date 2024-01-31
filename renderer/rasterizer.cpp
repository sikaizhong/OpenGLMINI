#include<rasterizer_state.h>
#include<rasterizer.h>
namespace CGL {
	std::vector<float> vertex = {
	 -0.5f, -0.5f, 0.0f,
	 0.5f, -0.5f, 0.0f,
	 0.0f,  0.5f, 0.0f
	};
	std::vector<unsigned> index = { 0,1,2 };
	inline bool is_point_inside_triangle(const vec3& A, const vec3& B, const vec3& C, int x, int y) {
		// Calculate barycentric coordinates
		float detT = (B.y - C.y) * (A.x - C.x) + (C.x - B.x) * (A.y - C.y);
		float alpha = ((B.y - C.y) * (x - C.x) + (C.x - B.x) * (y - C.y)) / detT;
		float beta = ((C.y - A.y) * (x - C.x) + (A.x - C.x) * (y - C.y)) / detT;
		float gamma = 1.0f - alpha - beta;
		// Check if point is inside the triangle
		return alpha >= 0 && beta >= 0 && (alpha + beta < 1);
	}

	inline bool is_point_inside_triangle(const vec3& A, const vec3& B, const vec3& C, int x, int y,
		float&alpha,float&beta,float&gamma) {
		// Calculate barycentric coordinates
		float detT = (B.y - C.y) * (A.x - C.x) + (C.x - B.x) * (A.y - C.y);
		alpha = ((B.y - C.y) * (x - C.x) + (C.x - B.x) * (y - C.y)) / detT;
		beta = ((C.y - A.y) * (x - C.x) + (A.x - C.x) * (y - C.y)) / detT;
		gamma = 1.0f - alpha - beta;
		// Check if point is inside the triangle
		return alpha >= 0 && beta >= 0 && (alpha + beta < 1);
		// do not use this one;
		//return alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 && gamma <= 1;
		/*float dot00 = (C.x - A.x) * (C.x - A.x) + (C.y - A.y) * (C.y - A.y);
		float dot01= (C.x - A.x) * (B.x - A.x) + (C.y - A.y) * (B.y - A.y);
		float dot02= (C.x - A.x) * (x - A.x) + (C.y - A.y) * (y - A.y);
		float dot11 = (B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y);
		float dot12= (B.x - A.x) * (x - A.x) + (B.y - A.y) * (y - A.y);
		float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
		alpha = (dot11 * dot02 - dot01 * dot12) * invDenom;
		beta = (dot00 * dot12 - dot01 * dot02) * invDenom;
		gamma = 1 - alpha - beta;
		return (alpha >= 0) && (beta >= 0) && (alpha + beta < 1);*/


	}

	void Rasterizer::set_up_vertex_array(std::vector<float>& vertices, std::vector<unsigned>& triangles) {
		m_vertex_buffer.clear();
		m_vertex_buffer.resize(vertices.size() / 3);
		for (int x = 0; x < vertices.size() / 3; x++) {
			m_vertex_buffer[x].m_world_pos[0] = vertices[3 * x + 0];
			m_vertex_buffer[x].m_world_pos[1] = vertices[3 * x + 1];
			m_vertex_buffer[x].m_world_pos[2] = vertices[3 * x + 2];
		}
		m_index_buffer.clear();
		m_index_buffer.resize(triangles.size());
		std::copy(triangles.begin(), triangles.end(), m_index_buffer.begin());
	}

	void Rasterizer::draw_by_index(int p_count) {
		m_frag_buffer.resize(m_width * m_height);
		for (auto& v : m_vertex_buffer)
			m_current_program->m_vertex_shader->execute(v, m_md_mat, m_v_mat, m_prj_mat);

		// send it to the screen coordinates;
		for (auto& v : m_vertex_buffer) {
			auto& data = m_vp_mat * v.m_clipping_pos;
			v.m_screen_pos[0] = data[0];
			v.m_screen_pos[1] = data[1];
			v.m_screen_pos[2] = data[2];
		}
		for (int x = 0; x < m_index_buffer.size() / 3; x++) {
			draw_triangle_barycentric_aux(m_vertex_buffer[m_index_buffer[3 * x + 0]],
				m_vertex_buffer[m_index_buffer[3 * x + 1]], m_vertex_buffer[m_index_buffer[3 * x + 2]],
				m_width, m_height);
		}
	}
	// get necessary info  for fragment;
	void Rasterizer::draw_triangle_barycentric_aux(VertexData& v0, VertexData& v1,
		VertexData& v2, int screen_width, int screen_height) {
		// get bounding box;
		float aabb_x0 = std::min(std::min(v0.m_screen_pos.x, v1.m_screen_pos.x), v2.m_screen_pos.x);
		float aabb_x1 = std::max(std::max(v0.m_screen_pos.x, v1.m_screen_pos.x), v2.m_screen_pos.x);
		float aabb_y0 = std::min(std::min(v0.m_screen_pos.y, v1.m_screen_pos.y), v2.m_screen_pos.y);
		float aabb_y1 = std::max(std::max(v0.m_screen_pos.y, v1.m_screen_pos.y), v2.m_screen_pos.y);

		// make sure it stays in the screen
		int aabb_x_min = std::floor(aabb_x0);
		int aabb_x_max = std::ceil(aabb_x1);
		int aabb_y_min = std::floor(aabb_y0);
		int aabb_y_max = std::ceil(aabb_y1);

		// This is clipping;
		aabb_x_max = std::min(aabb_x_max, screen_width - 1);
		aabb_y_max = std::min(aabb_y_max, screen_height - 1);

		aabb_x_min = std::max(aabb_x_min, 0);
		aabb_y_min = std::max(aabb_y_min, 0);
		for (int x = aabb_x_min; x <= aabb_x_max; x++) {
			for (int y = aabb_y_min; y <= aabb_y_max; y++) {
				float alpha, beta, gamma;
				if (is_point_inside_triangle(v0.m_screen_pos, v1.m_screen_pos, v2.m_screen_pos,
					x, y,alpha,beta,gamma)) {
					int idx = y * m_height + x;
					m_frag_buffer[idx].m_color[0] = alpha * 255;
					m_frag_buffer[idx].m_color[1] = beta * 255;
					m_frag_buffer[idx].m_color[2] = gamma * 255;
					auto color = m_current_program->m_fragment_shader->execute(m_frag_buffer[idx]);
					m_cur_framebuffer->CColorBuffer[3 * idx + 0] = color.r;
					m_cur_framebuffer->CColorBuffer[3 * idx + 1] = color.g;
					m_cur_framebuffer->CColorBuffer[3 * idx + 2] = color.b;
				}
			}
		}
	}

	void Rasterizer::init_default_shader() {
		m_current_program = std::make_shared<Program>();
		std::shared_ptr<VertexShader> vs = std::make_shared<VertexShader>();
		std::shared_ptr<FragmentShader> fs = std::make_shared<FragmentShader>();
		vs->execute = basic_vs;
		fs->execute = basic_fs;
		m_current_program->compile(vs, fs);
	}


	void Rasterizer::set_viewport(int x, int y, int w, int h) {
		 //ndc space -> screen space
		float hw  =  w * 0.5f;
		float hh  =  h * 0.5f;
		float hd = 255.f/2.0f;
		m_vp_mat[0][0] = hw;    m_vp_mat[0][1] = 0.0f;  m_vp_mat[0][2] = 0.0f; m_vp_mat[0][3] = x+hw;
		m_vp_mat[1][0] = 0.0f;	m_vp_mat[1][1] = hh;    m_vp_mat[1][2] = 0.0f; m_vp_mat[1][3] = y+hh;
		m_vp_mat[2][0] = 0.0f;  m_vp_mat[2][1] = 0.0f;  m_vp_mat[2][2] = hd;   m_vp_mat[2][3] = hd;
		m_vp_mat[3][0] = 0.0f;  m_vp_mat[3][1] = 0.0f;  m_vp_mat[3][2] = 0.0f; m_vp_mat[3][3] = 0.0f;
	}
	void Rasterizer::run() {
		set_viewport(0, 0, m_width, m_height);
		cgl_clear_color(0.2f, 0.4f, 0.6f, 1.0f);
		cgl_clear();
		set_up_vertex_array(vertex, index);
		draw_by_index();
	}


}