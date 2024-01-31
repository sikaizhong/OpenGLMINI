#pragma once
#include <limits>
#include<vector>
namespace CGL {
	class CFramebuffer {
	public:
		CFramebuffer(int width, int height) {
			m_height= height;
			m_width= width;
			CColorBuffer.resize(width * height * 3);// support rgba;
			CDepthBuffer.resize(width * height);// support rgba;
			std::fill(CDepthBuffer.begin(), CDepthBuffer.end(), -1000000000.0f);
			std::fill(CColorBuffer.begin(), CColorBuffer.end(), 0);
		}

		void clear() { CColorBuffer.clear(); CDepthBuffer.clear(); }
	public:
		int m_height;
		int m_width;
		std::vector<unsigned char > CColorBuffer;
		std::vector<float> CDepthBuffer;
		std::vector<float> CMultiSampleBuffer;

	};
}