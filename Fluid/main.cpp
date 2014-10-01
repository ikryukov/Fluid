//
//  main.cpp
//  Fluid
//
//  Created by Ilya Kryukov on 01.10.14.
//  Copyright (c) 2014 Ilya Kryukov. All rights reserved.
//

#include <iostream>

#include "IRender.h"

#include "GL/glew.h"

#include "GLFW/glfw3.h"

#include "glm/glm.hpp"

#include "IRender.h"
#include "OGLRender.h"

int main(int argc, const char * argv[]) {
	GLFWwindow* window;
	
	/* Initialize the library */
	if (!glfwInit())
		return -1;
	
	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}
	float ratio;
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	ratio = width / (float) height;
	
	glViewport(0, 0, width, height);
	
	OGLRender render;
	render.init(width, height);
	
	glfwSwapInterval(1);
	/* Make the window's context current */
	glfwMakeContextCurrent(window);
	
	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		render.render();
		
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	
	glfwDestroyWindow(window);
	glfwTerminate();
    return 0;
}
