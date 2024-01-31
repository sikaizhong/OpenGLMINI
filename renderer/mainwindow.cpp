#define NOMINMAX
#include<mainwindow.h>
#include<global_context.h>
#include<rasterizer.h>
namespace CGL {
	Gdiplus::Bitmap* WindowSystem::m_pbitmap = nullptr;
	WindowSystem::WindowSystem(int width, int height) : m_width(width), m_height(height),m_hdc(nullptr) {
		Gdiplus::GdiplusStartupInput gdiplusStartupInput;
		GdiplusStartup(&m_gdi_plus_token, &gdiplusStartupInput, NULL);
		// Create a window
		WNDCLASS wc = {};
		wc.lpfnWndProc = WindowProc;
		wc.hInstance = GetModuleHandle(NULL);
		wc.lpszClassName = "GDIPlusWindowClass";
		RegisterClass(&wc);
		m_hwnd = CreateWindowEx(
			0,                              // Optional window styles
			"GDIPlusWindowClass",          // Window class
			"OpenGLMINI",          // Window title
			WS_OVERLAPPEDWINDOW,            // Window style
			// Size and position
			CW_USEDEFAULT, CW_USEDEFAULT, m_width, m_height,
			NULL,       // Parent window    
			NULL,       // Menu
			GetModuleHandle(NULL),  // Instance handle
			NULL        // Additional application data
		);
		if (m_hwnd == NULL) {
			exit(0);
		}
		// Show the window
		ShowWindow(m_hwnd, SW_SHOWNORMAL);
		UpdateWindow(m_hwnd);
	}

	void WindowSystem::run() {
		MSG msg = {};
		while (GetMessage(&msg, NULL, 0, 0)) {
			TranslateMessage(&msg);
			DispatchMessage(&msg);
			auto& raster = g_runtime_global_context.m_render_system;
			set_up_bitmap(raster->m_cur_framebuffer.get());
			InvalidateRect(m_hwnd, NULL, TRUE);// force to call paint;
			calculate_fps();
		}
	}

	LRESULT CALLBACK WindowSystem::WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam) {
		switch (uMsg) {
		case WM_CREATE: {
			HMENU hMenu = CreateMenu();
			HMENU hFileMenu = CreateMenu();
			//HMENU hLoadMenu = CreateMenu();
			//AppendMenu(hFileMenu, MF_STRING, IDM_EXIT, "Exit");
			//AppendMenu(hFileMenu, MF_STRING, IDM_LOAD, "Load");
		    //AppendMenu(hMenu, MF_POPUP, reinterpret_cast<UINT_PTR>(hFileMenu), "File");
			SetMenu(hwnd, hMenu);
		} break;
		case WM_PAINT: {
			PAINTSTRUCT ps;
			HDC hdc = BeginPaint(hwnd, &ps);
			if (m_pbitmap != nullptr) {
				// Draw the image using GDI+
				Gdiplus::Graphics graphics(hdc);
				graphics.DrawImage(m_pbitmap, 0, 0, m_pbitmap->GetWidth(), m_pbitmap->GetHeight());
			}
			EndPaint(hwnd, &ps);
			return 0;
		}
		case WM_KEYDOWN:
			// Handle key down event
			if (wParam == VK_F1) {
				//render(raster);
				//raster.cgl_enable(CGL_DEPTH_TEST);
				//set_up_bitmap(&raster.framebuffer);
				//InvalidateRect(hwnd, NULL, TRUE);// force to call paint;
			}
			if (wParam == VK_F2) {
				//raster.cgl_disable(CGL_DEPTH_TEST);
			}
			if (wParam == VK_ESCAPE) {
				DestroyWindow(hwnd);
			}
			break;
		case WIM_CLOSE:
			DestroyWindow(hwnd);
			return 0;
		case WM_DESTROY:
			PostQuitMessage(0);
			return 0;
		default:
			return DefWindowProc(hwnd, uMsg, wParam, lParam);
		}
	}
	void WindowSystem::set_up_bitmap(CFramebuffer* framebuffer) {
		delete m_pbitmap;
		//BYTE
		// Create a Bitmap from raw data
		int stride = framebuffer->m_width * 3; // Assuming 24-bit RGB (3 bytes per pixel)
		Gdiplus::PixelFormat pixelFormat = PixelFormat24bppRGB; // Assuming 24-bit RGB
		m_pbitmap = new Gdiplus::Bitmap(framebuffer->m_width, framebuffer->m_height, stride, pixelFormat, framebuffer->CColorBuffer.data());
		m_pbitmap->RotateFlip(Gdiplus::Rotate180FlipX);
	}

	void WindowSystem::calculate_fps() {
		static float  fps = 0;
		static int    frameCount = 0;
		static float  currentTime = 0.0f;
		static float  lastTime = 0.0f;
		frameCount++;
		currentTime = timeGetTime() * 0.001f;
		std::string title = "OpenGLMINI - FPS: " + std::to_string(fps);
		if (currentTime - lastTime > 1.0f) {
			fps = (float)frameCount / (currentTime - lastTime);
			lastTime = currentTime;
			frameCount = 0;
		}
		SetWindowTextA(m_hwnd, title.c_str());
	}

}