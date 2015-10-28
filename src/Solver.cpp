//
//  Solver.cpp
//  Fluid
//
//  Created by Ilya Kryukov on 02.10.14.
//  Copyright (c) 2014 Ilya Kryukov. All rights reserved.
//

#include "Solver.h"
#include <algorithm>
#include <iostream>
#include <cmath>

Solver::Solver(int gridX, int gridY, int gridZ, float h)
	: m_gridX(gridX), m_gridY(gridY), m_gridZ(gridZ), m_h(h)
	, m_kCfl(1.0f)
    , Gravity(vec3(0, -9.8f, 0))
{
	for (float x = 0.0f; x < 10.0f; x += 1.0f)
	{
		for (float y = 0.0f; y < 4.0f; y += 1.0f)
		{
			for (float z = 0.0f; z < 4.0f; z += 1.0f)
			{
				Marker p;
				p.pos = vec3(x, y, z);
				m_markers.push_back(p);
			}
		}
	}
}

Solver::~Solver()
{
	
}

void Solver::simulate()
{
	std::cout << "solver::simulate" << std::endl;
	float dt = calculateCFL();
	std::cout << "time step dt = " << dt << std::endl;;
	dynamicGridUpdate();
    advanceVelocityField();
    applyExternalForces(dt);
}

int Solver::getIdx(int x, int y, int z)
{
	return z + (y + (x * m_gridY)) * m_gridZ;
}

float Solver::calculateCFL()
{
	float dt = 0.0f;
	float uMax = 1.0f;
	//for (int i = 0; i < m_cells.size(); ++i)
	for (std::map<int, Cell>::iterator it = m_mapCells.begin(); it != m_mapCells.end(); ++it)
	{
		uMax = std::max(uMax, it->second.u.x);
	}
	dt = m_kCfl * m_h / uMax;
	return dt;
}

void Solver::dynamicGridUpdate()
{
	for (std::map<int, Cell>::iterator it = m_mapCells.begin(); it != m_mapCells.end(); ++it)
	{
		it->second.mLayer = -1;
	}
	// update cells that currently have fluid in them
	for (int i = 0; i < m_markers.size(); ++i)
	{
		Marker& p = m_markers[i];

		ivec3 cellIdx(p.pos.x / m_h, p.pos.y / m_h, p.pos.z / m_h);
		Cell* pCell = NULL;
		bool isExist = getCell(cellIdx, &pCell);
		if (!isExist)
		{
			if (isInsideSimulationBounds(cellIdx))
			{
				Cell c;
				c.mType = Cell::tFluid;
				c.mLayer = 0;
				addCellToGrid(c, cellIdx);
			}
		}
		else if (pCell->mType != Cell::tSolid)
		{
			pCell->mType = Cell::tFluid;
			pCell->mLayer = 0;
		}
	}
	// create a buffer zone around the fluid
	for (int i = 1; i <= std::max(2, (int) m_kCfl); ++i)
	{
		for (std::map<int, Cell>::iterator it = m_mapCells.begin(); it != m_mapCells.end(); ++it)
		{
			Cell& C = it->second;
			if ((C.mType == Cell::tFluid || C.mType == Cell::tAir) && (C.mLayer == i - 1))
			{
				//for each of the six neighbors of C, N
				int dx[6] = { 1, -1, 0, 0, 0, 0 };
				int dy[6] = { 0, 0, 1, -1, 0, 0 };
				int dz[6] = { 0, 0, 0, 0, 1, -1 };
				for (int k = 0; k < 6; ++k)
				{
					Cell* pN = NULL;
					ivec3 nIdx(C.idx.x + dx[k], C.idx.y + dy[k], C.idx.z + dz[k]);
					bool isExist = getCell(nIdx, &pN);
					if (isExist)
					{
						if (pN->mLayer == -1 && pN->mType != Cell::tSolid)
						{
							pN->mType = Cell::tAir;
							pN->mLayer = i;
						}
                    }
                    else
                    {
                        Cell N;
                        N.mLayer = i;
                        if (isInsideSimulationBounds(nIdx))
                        {
                            N.mType = Cell::tAir;
                        }
                        else
                        {
                            N.mType = Cell::tSolid;
                        }
                        addCellToGrid(N, nIdx);
                    }
					
				}
			}
		}
	}
	std::map<int, Cell>::iterator it = m_mapCells.begin();
	while (it != m_mapCells.end())
	{
		if (it->second.mLayer == -1)
		{
			m_mapCells.erase(it);
			it = m_mapCells.begin();
		}
		else
		{
			++it;
		}		
	}
}

void Solver::advanceVelocityField()
{
    for (std::map<int, Cell>::iterator it = m_mapCells.begin(); it != m_mapCells.end(); ++it)
    {
        vec3 pos = it->second.pos;
        vec3 particle = traceParticle(pos, -0.5f);
        vec3 newU = getVelocity(particle);
        it->second.unew = newU;
    }
    for (std::map<int, Cell>::iterator it = m_mapCells.begin(); it != m_mapCells.end(); ++it)
    {
        it->second.u = it->second.unew;
    }
}

void Solver::applyExternalForces(float dt)
{
    for (std::map<int, Cell>::iterator it = m_mapCells.begin(); it != m_mapCells.end(); ++it)
    {
        if (it->second.mType == Cell::tFluid)
            it->second.u += dt * Gravity;
    }
}

int Solver::hash(ivec3 v)
{
	int r = 541 * v.x + 79 * v.y + 31 * v.z;
	return r;
}

bool Solver::getCell(ivec3 cellIdx, Cell** pCell)
{
	int hk = hash(cellIdx);
	if (m_mapCells.count(hk) == 0)
	{
		return false;
	}
	else 
	{
		*pCell = &m_mapCells[hk];
		return true;
	}
}

Cell& Solver::cell(int i, int j, int k)
{
    ivec3 idx(i, j, k);
    int hk = hash(idx);
    // TODO: check existance
    return m_mapCells[hk];
}

bool Solver::isInsideSimulationBounds(ivec3 cellIdx)
{
	return cellIdx.x >= 0 && cellIdx.x < m_gridX &&
		cellIdx.y >= 0 && cellIdx.y < m_gridY &&
		cellIdx.z >= 0 && cellIdx.z < m_gridZ;
}

void Solver::addCellToGrid(Cell& c, ivec3 cellIdx)
{
	int hk = hash(cellIdx);
	c.idx = cellIdx;
    c.pos = vec3(cellIdx.x * m_h, cellIdx.y * m_h, cellIdx.z * m_h);
	m_mapCells[hk] = c;
}

vec3 Solver::getVelocity(vec3 p)
{
    vec3 v;
    v.x = getInterpolatedValue(p.x / m_h, p.y / m_h - 0.5f, p.z / m_h - 0.5f, 0);
    v.y = getInterpolatedValue(p.x / m_h - 0.5f, p.y / m_h, p.z / m_h - 0.5f, 1);
    v.z = getInterpolatedValue(p.x / m_h - 0.5f, p.y / m_h - 0.5f, p.z / m_h, 2);
    return v;
}

float Solver::getInterpolatedValue(float x, float y, float z, int index)
{
    int i = floor(x);
    int j = floor(y);
    int k = floor(z);
    return (i+1-x) * (j+1-y) * (k+1-z) * cell(i, j, k).u[index] +
    (x-i) * (j+1-y) * (k+1-z) * cell(i+1, j, k).u[index] + (i+1-x) * (y-j) * (k+1-z) * cell(i, j+1, k).u[index] + (x-i) * (y-j) * (k+1-z) * cell(i+1, j+1, k).u[index] + (i+1-x) * (j+1-y) * (z-k) * cell(i, j, k+1).u[index] + (x-i) * (j+1-y) * (z-k) * cell(i+1, j, k+1).u[index] + (i+1-x) * (y-j) * (z-k) * cell(i, j+1, k+1).u[index] + (x-i) * (y-j) * (z-k) * cell(i+1, j+1, k+1).u[index];
}

vec3 Solver::traceParticle(vec3 p, float t)
{
    vec3 v = getVelocity(p);
    v = getVelocity(p + 0.5f * t * v);
    return p + t * v;
}