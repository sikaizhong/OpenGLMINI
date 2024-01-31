#include<global_context.h>
#include<mainwindow.h>
#include<application.h>
#include<rasterizer.h>
namespace CGL {

	void Application::initialize() {
		g_runtime_global_context.start_systems(std::string(""));
	}
	bool Application::tick() {
		g_runtime_global_context.m_render_system->run();
		g_runtime_global_context.m_window_system->run();
		return true;
	}
	void Application::shut_down() {
	}
}