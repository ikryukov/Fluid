//
//  Solver.h
//  Fluid
//
//  Created by Ilya Kryukov on 02.10.14.
//  Copyright (c) 2014 Ilya Kryukov. All rights reserved.
//

#ifndef __Fluid__Solver__
#define __Fluid__Solver__

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <glm/vec3.hpp>

using namespace boost::numeric::ublas;
using namespace glm;

typedef coordinate_matrix<float> SparseMatrix;
typedef vector<float> Vector;

struct Density {
	float density;
	float temp;
	
	inline Density() : density(0), temp(273) {
	}
	
	inline Density(float density, float temp)
	: density(density), temp(temp) {
	}
	
	inline Density operator+(const Density &v) const {
		return Density(density + v.density, temp + v.temp);
	}
	
	inline Density operator*(float f) const {
		return Density(density * f, temp * f);
	}
};

inline std::ostream& operator<<(std::ostream &os, const Density &d) {
	os << "(" << d.density << ", " << d.temp << ")";
	return os;
}

class Solver {
public:
	Solver(int gridX, int gridY, int gridZ, float dx);
	~Solver();
	void simulate(float dt);
private:
	int m_gridX, m_gridY, m_gridZ;
	float m_dx, m_invDx;
	float m_time;
	int m_numVoxels, m_Slice, m_VelocitySlice;
	
	vector<Density> m_Density0, m_Density1;
	vector<vec3> m_Curl, m_vcForce;
	vector<float> m_CurlMagnitude;
	vector<float> m_u0[3];
	vector<float> m_u1[3];
	vector<float> m_F[3];
	vector<float> m_Divergence;
	vector<float> m_Pressure;
	vector<int> m_Solid;
	
	// Pressure matrix
	vector<float> m_ADiag, m_APlusX, m_APlusY, m_APlusZ;
	
	void constructPressureMatrix();
	
	void project(float dt);
	void advectVelocity(float dt);
	void advectDensity(float dt);
	void calcForces();
	void integrate(float dt);
	
	// Utils
	int getIdx(int x, int y, int z);

};

#endif /* defined(__Fluid__Solver__) */
