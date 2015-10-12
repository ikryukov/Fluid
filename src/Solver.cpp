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

Solver::Solver(int gridX, int gridY, int gridZ, float h)
	: m_gridX(gridX), m_gridY(gridY), m_gridZ(gridZ), m_h(h)
	, m_kCfl(1.0f)
{
	for (float x = 0.0f; x < 10.0f; x += 1.0f)
	{
		for (float y = 0.0f; y < 5.0f; y += 1.0f)
		{
			for (float z = 0.0f; z < 10.0f; z += 1.0f)
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
					ivec3 nIdx(C.idx.x + dx[i], C.idx.y + dy[i], C.idx.z + dz[i]);
					bool isExist = getCell(nIdx, &pN);
					if (isExist)
					{
						if (pN->mLayer == -1 && pN->mType != Cell::tSolid)
						{
							pN->mType = Cell::tAir;
							pN->mLayer = i;
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

bool Solver::isInsideSimulationBounds(ivec3 cellIdx)
{
	return cellIdx.x >= 0 && cellIdx.x < m_gridX &&
		cellIdx.y >= 0 && cellIdx.y < m_gridY &&
		cellIdx.z >= 0 && cellIdx.z < m_gridZ;
}

void Solver::addCellToGrid(Cell c, ivec3 cellIdx)
{
	int hk = hash(cellIdx);
	c.idx = cellIdx;
	m_mapCells[hk] = c;
}