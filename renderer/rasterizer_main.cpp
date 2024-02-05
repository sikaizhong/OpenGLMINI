#include<iostream>
#include<mainwindow.h>
#include<application.h>
#include<global_context.h>
#include<rasterizer.h>
#include<mesh_io.h>
using namespace std;


int main() {
	CGL::Application app;
	app.initialize();
	CGL::Mesh msh;
	std::string filename = std::string("../../model/african_head/african_head.obj");
	//std::string filename = std::string("../../model/diablo3_pose/diablo3_pose.obj");

	CGL::load_obj_file(filename, &msh);
	auto& raster = CGL::g_runtime_global_context.m_render_system;
	raster->set_up_vertex_array(&msh);
	//raster->set_up_vertex_array(msh.m_verts, msh.facet_vrt);
	//raster->
	app.tick();
	app.shut_down();
}