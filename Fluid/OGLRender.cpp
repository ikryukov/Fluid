//
//  OGLRender.cpp
//  Fluid
//
//  Created by Ilya Kryukov on 02.10.14.
//  Copyright (c) 2014 Ilya Kryukov. All rights reserved.
//

#include "OGLRender.h"

#include "glm/glm.hpp"
#include "glm/ext.hpp"

using namespace glm;

void OGLRender::init(int width, int height)
{
	m_width = width;
	m_height = height;
	glViewport(0, 0, m_width, m_height);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (float)width/(float) height, 0.1, 1000.0);
	glMatrixMode(GL_MODELVIEW);
}

void OGLRender::render()
{
	glColor3f(1.0f, 0.85f, 0.35f);
	glBegin(GL_TRIANGLES);
	{
		glVertex3f(  0.0,  0.6, 0.0);
		glVertex3f( -0.2, -0.3, 0.0);
		glVertex3f(  0.2, -0.3 ,0.0);
	}
	glEnd();
}

void OGLRender::render(Solver* fsolver)
{
	glDepthMask(GL_TRUE);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
	int gridSizeX = fsolver->getGridSizeX(),
	gridSizeY = fsolver->getGridSizeY(),
	gridSizeZ = fsolver->getGridSizeZ();
	float dx = fsolver->getVoxelSize();
	const vector<Density> &density = fsolver->getDensity();
	const vector<float> &curlMagnitude = fsolver->getCurlMagnitude();
	const vector<int> solid = fsolver->getSolid();
	
	static double angleX = 0.0f, angleY = 0.0f;
	double distZ = -100.0;
	int displayMode = 0;
	int displayVelocities = 0;
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (float)m_width/(float) m_height, 0.1, 1000.0);
	gluLookAt(0,0,distZ,0,0,0,0,1,0);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	//glRotated(angleY, 1, 0, 0);
	//glRotated(angleX, 0, 1, 0);
	
	angleX += 0.5;
	angleY += 0.5;
	
	if (angleX > 360) angleX = 0.0;
	if (angleY > 360) angleY = 0.0;
	
	glTranslatef(-gridSizeX/2.0f, -gridSizeY/2.0f, -gridSizeZ/2.0f);
	//glTranslatef(-gridSizeX/2.0f, 0.0f, 0.0f);
	
	glCullFace(GL_NONE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	
	glBegin(GL_QUADS);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glVertex3i(0,   0, 0);
	glVertex3i(0,   gridSizeY, 0);
	glVertex3i(gridSizeX, gridSizeY, 0);
	glVertex3i(gridSizeX, 0, 0);
	
	glVertex3i(0,   0, gridSizeZ);
	glVertex3i(0,   gridSizeY, gridSizeZ);
	glVertex3i(gridSizeX, gridSizeY, gridSizeZ);
	glVertex3i(gridSizeX, 0, gridSizeZ);
	
	glVertex3i(0, 0,   0);
	glVertex3i(0, 0,   gridSizeZ);
	glVertex3i(0, gridSizeY, gridSizeZ);
	glVertex3i(0, gridSizeY, 0);
	
	glVertex3i(gridSizeX, 0,   0);
	glVertex3i(gridSizeX, 0,   gridSizeZ);
	glVertex3i(gridSizeX, gridSizeY, gridSizeZ);
	glVertex3i(gridSizeX, gridSizeY, 0);
	
	glVertex3i(0,   0, 0);
	glVertex3i(0,   0, gridSizeZ);
	glVertex3i(gridSizeX, 0, gridSizeZ);
	glVertex3i(gridSizeX, 0, 0);
	
	glVertex3i(0,   gridSizeY, 0);
	glVertex3i(0,   gridSizeY, gridSizeZ);
	glVertex3i(gridSizeX, gridSizeY, gridSizeZ);
	glVertex3i(gridSizeX, gridSizeY, 0);
	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_QUADS);
	glColor4f(.5f, .5f, 1.0f, 1.0f);
	for (int z=0, pos = 0; z<gridSizeZ; ++z) {
		for (int y=0; y<gridSizeY; ++y) {
			for (int x=0; x<gridSizeX; ++x, ++pos) {
				if (!solid[pos])
					continue;
				glVertex3i(x,   y, z);
				glVertex3i(x,   y+1, z);
				glVertex3i(x+1, y+1, z);
				glVertex3i(x+1, y, z);
				
				glVertex3i(x,   y, z+1);
				glVertex3i(x,   y+1, z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y, z+1);
				
				glVertex3i(x, y,   z);
				glVertex3i(x, y,   z+1);
				glVertex3i(x, y+1, z+1);
				glVertex3i(x, y+1, z);
				
				glVertex3i(x+1, y,   z);
				glVertex3i(x+1, y,   z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y+1, z);
				
				glVertex3i(x,   y, z);
				glVertex3i(x,   y, z+1);
				glVertex3i(x+1, y, z+1);
				glVertex3i(x+1, y, z);
				
				glVertex3i(x,   y+1, z);
				glVertex3i(x,   y+1, z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y+1, z);
				
			}
		}
	}
	glEnd();
	glDepthMask(GL_FALSE);
	glBegin(GL_QUADS);
	for (int z=0, pos = 0; z<gridSizeZ; ++z) {
		for (int y=0; y<gridSizeY; ++y) {
			for (int x=0; x<gridSizeX; ++x, ++pos) {
				if (solid[pos])
					continue;
				if (displayMode == 0) {
					float f = .5f * density[pos].density;
					glColor4f(1.0f, 1.0f, 1.0f, f);
				} else {
					float f = .01f * curlMagnitude[pos];
					glColor4f(0.0f, 0.0f, 1.0f, f);
				}
				glVertex3i(x,   y, z);
				glVertex3i(x,   y+1, z);
				glVertex3i(x+1, y+1, z);
				glVertex3i(x+1, y, z);
				
				glVertex3i(x,   y, z+1);
				glVertex3i(x,   y+1, z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y, z+1);
				
				glVertex3i(x, y,   z);
				glVertex3i(x, y,   z+1);
				glVertex3i(x, y+1, z+1);
				glVertex3i(x, y+1, z);
				
				glVertex3i(x+1, y,   z);
				glVertex3i(x+1, y,   z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y+1, z);
				
				glVertex3i(x,   y, z);
				glVertex3i(x,   y, z+1);
				glVertex3i(x+1, y, z+1);
				glVertex3i(x+1, y, z);
				
				glVertex3i(x,   y+1, z);
				glVertex3i(x,   y+1, z+1);
				glVertex3i(x+1, y+1, z+1);
				glVertex3i(x+1, y+1, z);
				
			}
		}
	}
	
	glEnd();
}