#include<rasterizer_state.h>
#include<rasterizer.h>
#include<glm/glm.hpp>
#include<glm/gtc/matrix_transform.hpp>
#undef GLM_FORCE_RADIANS
namespace CGL {
	void convert(glm::mat4& gm, mat4& res) {
		for (int x = 0; x < 4; x++) {
			for (int y = 0; y < 4; y++) {
				res[x][y] = gm[y][x];
			}
		}
	}
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
		float& alpha, float& beta, float& gamma) {
		// Calculate barycentric coordinates
		float detT = (B.y - C.y) * (A.x - C.x) + (C.x - B.x) * (A.y - C.y);
		alpha = ((B.y - C.y) * (x - C.x) + (C.x - B.x) * (y - C.y)) / detT;
		beta = ((C.y - A.y) * (x - C.x) + (A.x - C.x) * (y - C.y)) / detT;
		gamma = 1.0f - alpha - beta;
		// Check if point is inside the triangle
		return alpha >= 0 && beta >= 0 && (alpha + beta < 1);
		// do not use this one;
		//return alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 && gamma <= 1;
	}

	void Rasterizer::set_viewport(int x, int y, int w, int h) {
		// [-1,1]*[-1,1]*[-1,1] to the screen space;
		//ndc space -> screen space
		float hw = w * 0.5f;
		float hh = h * 0.5f;
		float hd = 2.0f / 2.0f;
		m_vp_mat[0][0] = hw;    m_vp_mat[0][1] = 0.0f;  m_vp_mat[0][2] = 0.0f; m_vp_mat[0][3] = x + hw;
		m_vp_mat[1][0] = 0.0f;	m_vp_mat[1][1] = hh;    m_vp_mat[1][2] = 0.0f; m_vp_mat[1][3] = y + hh;
		m_vp_mat[2][0] = 0.0f;  m_vp_mat[2][1] = 0.0f;  m_vp_mat[2][2] = hd;   m_vp_mat[2][3] = 0.0f;
		m_vp_mat[3][0] = 0.0f;  m_vp_mat[3][1] = 0.0f;  m_vp_mat[3][2] = 0.0f; m_vp_mat[3][3] = 1.0f;
	}
	void Rasterizer::cgl_look_at(const vec3& eye, const vec3& center, const vec3& up) {
		// translate the camera to the origin: T;
		// rotate the camera such that the 3 axis are along with the canonical axis:R;
		// camera gazes at z directions; if we put camera at (0,0,positive), gaze at (0,0,0)
		// it basic transforms z into z-positive and keep other unchanged;
		// this will generate a negative z values at the most time which is not wanted 
		// due to positive z near and far;
		// we could flip it 	
		// view=R*T;	
		vec3 z_axis = normalize((eye - center));//n
		vec3 x_axis = normalize(cross(up, z_axis));//u
		vec3 y_axis = normalize(cross(z_axis, x_axis));//v
		// unflipped
		m_v_mat[0][0] = x_axis.x; m_v_mat[0][1] = x_axis.y;  m_v_mat[0][2] = x_axis.z; m_v_mat[0][3] = -dot(x_axis, eye);
		m_v_mat[1][0] = y_axis.x; m_v_mat[1][1] = y_axis.y;  m_v_mat[1][2] = y_axis.z; m_v_mat[1][3] = -dot(y_axis, eye);
		m_v_mat[2][0] = z_axis.x; m_v_mat[2][1] = z_axis.y;  m_v_mat[2][2] = z_axis.z; m_v_mat[2][3] = -dot(z_axis, eye);
		m_v_mat[3][0] = 0;        m_v_mat[3][1] = 0;         m_v_mat[3][2] = 0;        m_v_mat[3][3] = 1.0f;
	
		//auto gm=glm::lookAt(glm::vec3(eye[0], eye[1], eye[2]),
		//	glm::vec3(center[0], center[1], center[2]), glm::vec3(up[0], up[1], up[2]));
		//mat4 res;
		//convert(gm, m_v_mat);
		//res;


	}

	void Rasterizer::cgl_perspective(float fovy, float aspect, float near, float far) {
		//Setup perspective matrix (camera space -> homogeneous space)
		float rFovy = fovy * 3.14159265358979323846 / 180.0f;
		const float tanHalfFovy = std::tan(rFovy * 0.5f);
		float f_n = far - near;
		m_prj_mat[0][0] = 1.0f / (aspect * tanHalfFovy); m_prj_mat[0][1] = 0.0f;			 m_prj_mat[0][2] = 0.0f;					    m_prj_mat[0][3] = 0.0f;
		m_prj_mat[1][0] = 0.0f;						    m_prj_mat[1][1] = 1.0f / tanHalfFovy;m_prj_mat[1][2] = 0.0f;				        m_prj_mat[1][3] = 0.0f;
		m_prj_mat[2][0] = 0.0f;						    m_prj_mat[2][1] = 0.0f;			     m_prj_mat[2][2] = (far + near) / f_n;	        m_prj_mat[2][3] = 2.0f * near * far / f_n;
		m_prj_mat[3][0] = 0.0f;						    m_prj_mat[3][1] = 0.0f;				 m_prj_mat[3][2] = -1.0f;	                    m_prj_mat[3][3] = 0.0f;
		//vec4 data = vec4(0.5, 0.5, 0.01, 1);
		//auto res = m_prj_mat * data;
		//res;
		//auto gm = glm::perspective(rFovy, aspect, near, far);
		//mat4 res;
		//convert(gm, res);
		//res;
	
	}   

	



	void Rasterizer::cgl_orth(float left, float right, float bottom, float top, float near, float far) {
		//Setup orthogonal matrix (camera space -> homogeneous space)
		// goal: make the the point in the frustum [l,r]*[b,t]*[-n,-f] to [-1,1]*[-1,1]*[-1,1];
		// n and f are positive;
		// translation: move to center (right+left)/2, (top+bot)/2, -(near+far)/2; 
		// scaling: 1/(right-left), 1/(top-bottom), 1/(far-near);
		// scaling+translation;
		// after model transformation: I put it at the origin and make is diagonal be 1;
		// after view transformation: z values is translated regarding to the camera pos;
		// it is always negative if we put camera quite far away;
		// we can flip z zxis in view matrix or here;
		m_prj_mat[0][0] = 2.0f / (right - left); m_prj_mat[0][1] = 0.0f;                  m_prj_mat[0][2] = 0.0f;                         m_prj_mat[0][3] = -(right + left) / (right - left);
		m_prj_mat[1][0] = 0.0f;				     m_prj_mat[1][1] = 2.0f / (top - bottom); m_prj_mat[1][2] = 0.0f;                         m_prj_mat[1][3] = -(top + bottom) / (top - bottom);
		m_prj_mat[2][0] = 0.0f;                  m_prj_mat[2][1] = 0.0f;                  m_prj_mat[2][2] = 2.0f / (far - near);          m_prj_mat[2][3] = (far + near) / (far - near);
		m_prj_mat[3][0] = 0.0f;                  m_prj_mat[3][1] = 0.0f;                  m_prj_mat[3][2] = 0;                            m_prj_mat[3][3] = 1.0f;
	}

	void Rasterizer::run() {
		set_viewport(0, 0, m_width, m_height);
		cgl_clear_color(0.2f, 0.4f, 0.6f, 1.0f);
		cgl_clear();
		vec3 eye(0, 0, 2.0f);
		vec3 center(0, 0, 0.0f);
		vec3 up(0, 1, 0);
		cgl_look_at(eye, center, up);
		cgl_perspective(60, 1.0, 0.01, 200.0f);

		//cgl_orth(-1.0f, 1.0f, -1.0f, 1.0f, 0.01f, 200.0f);
		//set_up_vertex_array(vertex, index);
		draw_by_index();
	}







	bool Rasterizer::depth_test(int x, int y, float new_depth, int screen_width, int screen_height) {
		if (m_init_render_state & CGL_DEPTH_TEST) {
			int index = y * screen_width + x;
			if (m_cur_framebuffer->CDepthBuffer[index] < new_depth) {
				m_cur_framebuffer->CDepthBuffer[index] = new_depth;
				return true;
			}
			return false;
		}
		return true;
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

	void Rasterizer::set_up_vertex_array(std::vector<vec3>& vertices, std::vector<int>& triangles) {
		m_vertex_buffer.clear();
		m_vertex_buffer.resize(vertices.size());
		for (int x = 0; x < vertices.size(); x++) {
			m_vertex_buffer[x].m_world_pos[0] = vertices[x][0];
			m_vertex_buffer[x].m_world_pos[1] = vertices[x][1];
			m_vertex_buffer[x].m_world_pos[2] = vertices[x][2];
		}
		m_index_buffer.clear();
		m_index_buffer.resize(triangles.size());
		std::copy(triangles.begin(), triangles.end(), m_index_buffer.begin());
	}

	void Rasterizer::set_up_vertex_array(Mesh* msh) {
		auto& vertices = msh->m_verts;
		auto& normals = msh->m_normals;
		auto& triangles = msh->facet_vrt;
		m_vertex_buffer.clear();
		m_vertex_buffer.resize(vertices.size());
		for (int x = 0; x < vertices.size(); x++) { 
			m_vertex_buffer[x].m_world_pos[0] = vertices[x][0];
			m_vertex_buffer[x].m_world_pos[1] = vertices[x][1];
			m_vertex_buffer[x].m_world_pos[2] = vertices[x][2];
		}
		for (int x = 0; x < vertices.size(); x++) {
			m_vertex_buffer[x].m_world_normal[0] = normals[x][0];
			m_vertex_buffer[x].m_world_normal[1] = normals[x][1];
			m_vertex_buffer[x].m_world_normal[2] = normals[x][2];
		}


		m_index_buffer.clear();
		m_index_buffer.resize(triangles.size());
		std::copy(triangles.begin(), triangles.end(), m_index_buffer.begin());
	}



	void Rasterizer::draw_by_index(int p_count) {
		m_frag_buffer.resize(m_width * m_height);
		for (auto& v : m_vertex_buffer)
			m_current_program->m_vertex_shader->execute(v, m_md_mat, m_v_mat, m_prj_mat);
		//std::ofstream out("z.txt");
		//for (auto& v : m_vertex_buffer)
		//	out << v.m_cpos[2] << std::endl;
		//exit(0);
		//getchar();
		

		// clipping 
		// really has to get the real vertex data, not just index;
		std::vector<VertexData> vertex_back_up;
		for (int x = 0; x < m_index_buffer.size() / 3; x++) {
			auto vs = homogeneous_clipping_SutherlandHodgeman_aux(m_vertex_buffer[m_index_buffer[3 * x + 0]],
				m_vertex_buffer[m_index_buffer[3 * x + 1]], m_vertex_buffer[m_index_buffer[3 * x + 2]],
				0.01, 200);
			for (int m = 0; m < int(vs.size()) - 2; m++) {
					vertex_back_up.push_back(vs[0]);
					vertex_back_up.push_back(vs[m + 1]);
					vertex_back_up.push_back(vs[m + 2]);
			}
		}
		// get ndc 
		// send it to the screen coordinates;
		for (auto& v : vertex_back_up) {
			//v.m_cpos[0] = v.m_cpos[0] / v.m_cpos[3];
			//v.m_cpos[1] = v.m_cpos[1] / v.m_cpos[3];
			//v.m_cpos[2] = v.m_cpos[2] / v.m_cpos[3];
			//v.m_cpos[3] = v.m_cpos[3] / v.m_cpos[3];

			auto& data = m_vp_mat * v.m_cpos;
			v.m_screen_pos[0] = data[0]  /data[3];
			v.m_screen_pos[1] = data[1] / data[3];
			v.m_screen_pos[2] = data[2] / data[3];
		}
		
		for (int x = 0; x < vertex_back_up.size() / 3; x++) {
				draw_triangle_barycentric_aux(vertex_back_up[3*x+0],
					vertex_back_up[3 * x + 1], vertex_back_up[3 * x + 2],
				m_width, m_height);
		}


		//for (int x = 0; x < m_index_buffer.size() / 3; x++) {
		//	draw_triangle_barycentric_aux(m_vertex_buffer[m_index_buffer[3 * x + 0]],
		//		m_vertex_buffer[m_index_buffer[3 * x + 1]], m_vertex_buffer[m_index_buffer[3 * x + 2]],
		//		m_width, m_height);
		//}
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
					x, y, alpha, beta, gamma)) {
					float z = alpha * v0.m_screen_pos.z + beta * v1.m_screen_pos.z + gamma * v2.m_screen_pos.z;
					if (depth_test(x, y, z, m_width, m_height)) {
						int idx = y * m_height + x;

						m_frag_buffer[idx].m_pos[0] = alpha * v0.m_world_pos[0] + beta * v1.m_world_pos[0] + gamma * v2.m_world_pos[0];
						m_frag_buffer[idx].m_pos[1] = alpha * v0.m_world_pos[1] + beta * v1.m_world_pos[1] + gamma * v2.m_world_pos[1];
						m_frag_buffer[idx].m_pos[2] = alpha * v0.m_world_pos[2] + beta * v1.m_world_pos[2] + gamma * v2.m_world_pos[2];
						//alpha = 1;
						//beta = gamma = 0;
						vec3 v10 = v1.m_world_pos - v0.m_world_pos;
						vec3 v20 = v2.m_world_pos - v0.m_world_pos;
						vec3 norm = cross(v10, v20);
						m_frag_buffer[idx].m_nor[0] = norm[0];
						m_frag_buffer[idx].m_nor[1] = norm[1];
						m_frag_buffer[idx].m_nor[2] = norm[2];



						//m_frag_buffer[idx].m_nor[0] = alpha * v0.m_world_normal[0] + beta * v1.m_world_normal[0] + gamma * v2.m_world_pos[0];
						//m_frag_buffer[idx].m_nor[1] = alpha * v0.m_world_normal[1] + beta * v1.m_world_normal[1] + gamma * v2.m_world_pos[1];
						//m_frag_buffer[idx].m_nor[2] = alpha * v0.m_world_normal[2] + beta * v1.m_world_normal[2] + gamma * v2.m_world_pos[2];


						//m_frag_buffer[idx].m_color[0] = alpha * 255;
						//m_frag_buffer[idx].m_color[1] = beta * 255;
						//m_frag_buffer[idx].m_color[2] = gamma * 255;
						auto color = m_current_program->m_fragment_shader->execute(m_frag_buffer[idx]);
						m_cur_framebuffer->CColorBuffer[3 * idx + 0] = color.r;
						m_cur_framebuffer->CColorBuffer[3 * idx + 1] = color.g;
						m_cur_framebuffer->CColorBuffer[3 * idx + 2] = color.b;
					}
				}
			}
		}
	}

	std::vector<VertexData> Rasterizer::homogeneous_clipping(VertexData& v0, VertexData& v1, VertexData& v2, float near, float far) {
		return homogeneous_clipping_fast_aux(v0, v1, v1, near, far);
	}


	std::vector<VertexData> Rasterizer::homogeneous_clipping_fast_aux(VertexData& v0, VertexData& v1, VertexData& v2, float near, float far) {
		std::vector<VertexData> res;
		return res;
	}

	std::vector<VertexData> Rasterizer::homogeneous_clipping_SutherlandHodgeman_aux(VertexData& v0, VertexData& v1, VertexData& v2, float near, float far) {
		//return { v0,v1,v2 };
		std::vector<VertexData> res;
		auto is_point_inside_frustum = [](const vec4& p, float& near, float& far)->bool {
			return (p.x <= p.w && p.x >= -p.w) && (p.y <= p.w && p.y >= -p.w)
				&& (p.z <= p.w && p.z >= -p.w) 
				&& (p.w <= far && p.w >= near);
			};
		//Totally inside
		if (is_point_inside_frustum(v0.m_cpos, near, far) &&
			is_point_inside_frustum(v1.m_cpos, near, far) &&
			is_point_inside_frustum(v2.m_cpos, near, far))
		{
			return { v0,v1,v2 };
			//return {};
		}
		//return {};
		// totally outside;
		//Totally outside
		if (v0.m_cpos.w < near && v1.m_cpos.w < near && v2.m_cpos.w < near)
			return{};
		if (v0.m_cpos.w > far && v1.m_cpos.w > far && v2.m_cpos.w > far)
			return{};
		if (v0.m_cpos.x > v0.m_cpos.w && v1.m_cpos.x > v1.m_cpos.w && v2.m_cpos.x > v2.m_cpos.w)
			return{};
		if (v0.m_cpos.x < -v0.m_cpos.w && v1.m_cpos.x < -v1.m_cpos.w && v2.m_cpos.x < -v2.m_cpos.w)
			return{};
		if (v0.m_cpos.y > v0.m_cpos.w && v1.m_cpos.y > v1.m_cpos.w && v2.m_cpos.y > v2.m_cpos.w)
			return{};
		if (v0.m_cpos.y < -v0.m_cpos.w && v1.m_cpos.y < -v1.m_cpos.w && v2.m_cpos.y < -v2.m_cpos.w)
			return{};
		if (v0.m_cpos.z > v0.m_cpos.w && v1.m_cpos.z > v1.m_cpos.w && v2.m_cpos.z > v2.m_cpos.w)
			return{};
		if (v0.m_cpos.z < -v0.m_cpos.w && v1.m_cpos.z < -v1.m_cpos.w && v2.m_cpos.z < -v2.m_cpos.w)
			return{};
		//return {v0,v1,v2};

		//return {};
		std::vector<VertexData> insideVertices;
		std::vector<VertexData> tmp = { v0, v1, v2 };
		enum Axis { X = 0, Y = 1, Z = 2 };

		//w=x plane & w=-x plane
		{
			insideVertices = clipingSutherlandHodgemanAux(tmp, Axis::X, +1);
			tmp = insideVertices;
			
			insideVertices = clipingSutherlandHodgemanAux(tmp, Axis::X, -1);
			tmp = insideVertices;
			
			//return insideVertices;
		}
		
		//w=y plane & w=-y plane
		{
			insideVertices = clipingSutherlandHodgemanAux(tmp, Axis::Y, +1);
			tmp = insideVertices;
			
			insideVertices = clipingSutherlandHodgemanAux(tmp, Axis::Y, -1);
			tmp = insideVertices;
		}

		//w=z plane & w=-z plane
		{
			insideVertices = clipingSutherlandHodgemanAux(tmp, Axis::Z, +1);
			tmp = insideVertices;

			insideVertices = clipingSutherlandHodgemanAux(tmp, Axis::Z, -1);
			tmp = insideVertices;
		}
		
		//w=1e-5 plane
		{
			insideVertices = {};
			int numVerts = tmp.size();
			constexpr float wClippingPlane = 1e-5;
			for (int i = 0; i < numVerts; ++i)
			{
				const auto& begVert = tmp[(i - 1 + numVerts) % numVerts];
				const auto& endVert = tmp[i];
				float begIsInside = (begVert.m_cpos.w < wClippingPlane) ? -1 : 1;
				float endIsInside = (endVert.m_cpos.w < wClippingPlane) ? -1 : 1;
				//One of them is outside
				if (begIsInside * endIsInside < 0)
				{
					// t = (w_clipping_plane-w1)/((w1-w2)
					float t = (wClippingPlane - begVert.m_cpos.w) / (begVert.m_cpos.w - endVert.m_cpos.w);
					auto intersectedVert = VertexData::lerp(begVert, endVert, t);
					insideVertices.push_back(intersectedVert);
				}
				//If current vertices is inside
				if (endIsInside > 0)
				{
					insideVertices.push_back(endVert);
				}
			}
		}

		return insideVertices;


	}

	std::vector<VertexData> Rasterizer::clipingSutherlandHodgemanAux(
		const std::vector<VertexData>& polygon,
		const int& axis,
		const int& side)
	{
		std::vector<VertexData> insidePolygon;

		int numVerts = polygon.size();
		for (int i = 0; i < numVerts; ++i)
		{
			const auto& begVert = polygon[(i - 1 + numVerts) % numVerts];
			const auto& endVert = polygon[i];
			char begIsInside = ((side * (begVert.m_cpos[axis]) <= begVert.m_cpos.w) ? 1 : -1);
			char endIsInside = ((side * (endVert.m_cpos[axis]) <= endVert.m_cpos.w) ? 1 : -1);
			//One of them is outside
			if (begIsInside * endIsInside < 0)
			{
				// t = (w1 - y1)/((w1-y1)-(w2-y2))
				float t = (begVert.m_cpos.w - side * begVert.m_cpos[axis])
					/ ((begVert.m_cpos.w - side * begVert.m_cpos[axis]) - (endVert.m_cpos.w - side * endVert.m_cpos[axis]));
				auto intersectedVert = VertexData::lerp(begVert, endVert, t);
				insidePolygon.push_back(intersectedVert);
			}
			//If current vertices is inside
			if (endIsInside > 0)
			{
				insidePolygon.push_back(endVert);
			}
		}
		return insidePolygon;
	}


	void Rasterizer::init_default_shader() {
		m_current_program = std::make_shared<Program>();
		std::shared_ptr<VertexShader> vs = std::make_shared<VertexShader>();
		std::shared_ptr<FragmentShader> fs = std::make_shared<FragmentShader>();
		//vs->execute = basic_vs;
		//fs->execute = basic_fs;
		vs->execute = Phong_vs;
		fs->execute = Phong_fs;
		m_current_program->compile(vs, fs);
	}


	
	


}