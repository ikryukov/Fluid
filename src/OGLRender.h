//
//  OGLRender.h
//  Fluid
//
//  Created by Ilya Kryukov on 02.10.14.
//  Copyright (c) 2014 Ilya Kryukov. All rights reserved.
//

#ifndef Fluid_OGLRender_h
#define Fluid_OGLRender_h

#include "IRender.h"

#include "GL/glew.h"

#include "Solver.h"

class OGLRender : public IRender {
	
private:
	int m_width;
	int m_height;
	
public:
	void init(int width, int height);
	void render();
	void render(Solver* fsolver);

};

#endif
