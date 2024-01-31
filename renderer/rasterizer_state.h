#pragma once 
#include<memory>
namespace CGL {
#define CGL_DEPTH_TEST                     0b00000001
#define CGL_BACK_FACE_CULLING              0b00000010
	struct RenderState {
		static const int depth_test = 0b00000001;
		static const int face_culling = 0b00000010;
		static const int STATE_C = 0b00000100;
		static const int STATE_D = 0b00001000;
		static const int STATE_E = 0b00010000;
		static const int STATE_F = 0b00100000;
		static const int STATE_G = 0b01000000;
		static const int STATE_H = 0b10000000;
	};
}