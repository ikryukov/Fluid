//
//  Solver.cpp
//  Fluid
//
//  Created by Ilya Kryukov on 02.10.14.
//  Copyright (c) 2014 Ilya Kryukov. All rights reserved.
//

#include "Solver.h"
#include <iostream>

Solver::Solver(int gridX, int gridY, int gridZ, float dx)
: m_gridX(gridX), m_gridY(gridY), m_gridZ(gridZ), m_dx(dx)
{
	m_invDx = 1.0f / m_dx;
}

Solver::~Solver()
{
	
}

void Solver::simulate(float dt)
{
	std::cout << "solver::simulate " << std::endl;
}

int Solver::getIdx(int x, int y, int z)
{
	return z + (y + (x * m_gridY)) * m_gridZ;
}