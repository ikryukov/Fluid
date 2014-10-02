//
//  Solver.cpp
//  Fluid
//
//  Created by Ilya Kryukov on 02.10.14.
//  Copyright (c) 2014 Ilya Kryukov. All rights reserved.
//

#include "Solver.h"

Solver::Solver(int gridX, int gridY, int gridZ, float dx)
: m_gridX(gridX), m_gridY(gridY), m_gridZ(gridZ), m_dx(dx), m_time(0.0f)
{
	int numVelocity = (m_gridX + 1) * (m_gridY + 1) * (m_gridZ + 1);
	m_numVoxels = m_gridX * m_gridY * m_gridZ;
	m_Slice = gridX * gridY;
	m_VelocitySlice = (gridX + 1) * (gridY + 1);
	m_invDx = 1.0f / m_dx;
	
	m_Density0.resize(m_numVoxels);
	m_Density1.resize(m_numVoxels);
	
	m_Curl.resize(m_numVoxels);
	m_CurlMagnitude.resize(m_numVoxels);
	m_vcForce.resize(m_numVoxels);
	m_Divergence.resize(m_numVoxels);
	m_Pressure.resize(m_numVoxels);
	m_Solid.resize(m_numVoxels);
	
	for (int x = 0; x < m_gridX; ++x)
	{
		for (int y = 0; y < m_gridY; ++y)
		{
			for (int z = 0; z < m_gridZ; ++z)
			{
				m_Solid[z + (y + (x * m_gridY))* m_gridZ] = 0;
			}
		}
	}
	
	for (int i = 0; i < 3; ++i)
	{
		m_F[i].resize(numVelocity);
		m_u0[i].resize(numVelocity);
		m_u1[i].resize(numVelocity);
	}
}

Solver::~Solver()
{
	
}