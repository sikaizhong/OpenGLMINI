#pragma once 
#include <string>
#include <memory>
namespace CGL {
	class WindowSystem;
	class Rasterizer;
	class RuntimeGlobalContext {
	public:
		void start_systems(const std::string& config_file_path);
		void shut_down_systems();

	public:
		std::shared_ptr<WindowSystem> m_window_system;
		std::shared_ptr<Rasterizer>   m_render_system;

	};
	extern RuntimeGlobalContext g_runtime_global_context;

}