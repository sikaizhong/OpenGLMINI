#pragma once
#include <windows.h>
#include <gdiplus.h>
#pragma comment(lib, "gdiplus.lib")
# pragma comment(lib, "secur32.lib")
# pragma comment(lib, "winmm.lib")
# pragma comment(lib, "dmoguids.lib")
# pragma comment(lib, "wmcodecdspuuid.lib")
# pragma comment(lib, "msdmo.lib")
# pragma comment(lib, "Strmiids.lib")

#include<framebuffer.h>

namespace CGL {
	class WindowSystem final {
	public:
		WindowSystem(int width, int height);
		~WindowSystem() {};
		void run();
	private:
		static LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);
		static void set_up_bitmap(CFramebuffer* framebuffer);
		void calculate_fps();
	
	private:
		HWND m_hwnd;
		HDC m_hdc;
		ULONG_PTR   m_gdi_plus_token;


		int m_width;
		int m_height;

		static Gdiplus::Bitmap* m_pbitmap;
	};


}