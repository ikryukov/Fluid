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
#include <map>
#include <glm/vec3.hpp>

using namespace glm;

struct Marker
{
	vec3 pos;
};

struct Cell
{
	enum Type
	{
		tUnknown = 0,
		tAir,
		tFluid,
		tSolid,
		tLast
	};
	vec3 u;
    vec3 unew;
    vec3 pos;
	ivec3 idx;
	float p;
	Type mType;
	int mLayer;
};

class Solver {
public:
	Solver(int gridX, int gridY, int gridZ, float h);
	~Solver();
	void simulate();
	
	inline int getGridSizeX() const { return m_gridX; }
	inline int getGridSizeY() const { return m_gridY; }
	inline int getGridSizeZ() const { return m_gridZ; }
	
	float calculateCFL();
	void dynamicGridUpdate();
    void advanceVelocityField();

	// Spatial grid
	std::map<int, Cell> m_mapCells;
	int m_gridX, m_gridY, m_gridZ;

private:
	
	std::vector<Marker> m_markers;
	float m_h; // width of a grid cell
	float m_kCfl;

	std::vector<Cell*> m_fluidCells;
	std::vector<Cell*> m_airCells;

	// Utils
	bool getCell(ivec3 cellIdx, Cell** pCell); // true if cell exists
    Cell& cell(int i, int j, int k);
	bool isInsideSimulationBounds(ivec3 cellIdx);
	int getIdx(int x, int y, int z);
	void addCellToGrid(Cell& c, ivec3 cellIdx);
	int hash(ivec3 v);
    vec3 traceParticle(vec3 p, float t);
    vec3 getVelocity(vec3 p);
    float getInterpolatedValue(float x, float y, float z, int index);

};

#endif /* defined(__Fluid__Solver__) */
