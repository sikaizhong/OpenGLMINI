#include<global_context.h>
#include<rasterizer.h>
#include<mainwindow.h>
namespace CGL {
	RuntimeGlobalContext g_runtime_global_context;
	void RuntimeGlobalContext::start_systems(const std::string& config_file_path) {
		m_window_system = std::make_shared<WindowSystem>(1000, 1000);
		m_render_system = std::make_shared<Rasterizer>(1000,1000);

	}
	void RuntimeGlobalContext::shut_down_systems() {

	}
}