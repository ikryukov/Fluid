//
//  Solver.h
//  Fluid
//
//  Created by Ilya Kryukov on 02.10.14.
//  Copyright (c) 2014 Ilya Kryukov. All rights reserved.
//

#ifndef __Fluid__Solver__
#define __Fluid__Solver__

#include <vector>
#include <glm/vec3.hpp>

using namespace glm;

class Solver {
public:
	Solver(int gridX, int gridY, int gridZ, float dx);
	~Solver();
	void simulate(float dt);
	
	inline int getGridSizeX() const { return m_gridX; }
	inline int getGridSizeY() const { return m_gridY; }
	inline int getGridSizeZ() const { return m_gridZ; }
	inline float getVoxelSize() const { return m_dx; }

private:
	int m_gridX, m_gridY, m_gridZ;
	float m_dx, m_invDx;

	// Utils
	int getIdx(int x, int y, int z);

};

#endif /* defined(__Fluid__Solver__) */
