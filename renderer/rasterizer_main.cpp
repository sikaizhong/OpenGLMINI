#include<iostream>
#include<mainwindow.h>
#include<application.h>
using namespace std;


int main() {
	CGL::Application app;
	app.initialize();
	app.tick();
	app.shut_down();
}