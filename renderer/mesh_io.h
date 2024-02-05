#pragma once 
#include <vector>
#include <string>
#include <fstream>
#include<mesh.h>
namespace CGL {
	void load_obj_file(std::string path, Mesh* mesh);
}