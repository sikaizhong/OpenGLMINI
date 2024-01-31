#pragma once 
namespace CGL {

	class Application {

	public:
		void initialize();
		bool tick();
		void shut_down();
	private:
		bool m_run;
	};


}