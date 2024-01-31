#pragma once 
#include<memory>
#include<functional>
#include<matrix.h>
#include<vector.h>
#include<color.h>
namespace CGL {
	class VertexData;
	class IShader {
	public:
	};

	class VertexShader :public IShader {

	public:
		std::function<void(VertexData&, mat4&, mat4&, mat4&)> execute;
	};

	class FragmentShader : public IShader {
	public:
		std::function<Color(FragmentData&)> execute;
	};
}