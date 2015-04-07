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
#include "aabb.h"


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
	
	inline int getGridSizeX() const { return m_gridX; }
	inline int getGridSizeY() const { return m_gridY; }
	inline int getGridSizeZ() const { return m_gridZ; }
	inline float getVoxelSize() const { return m_dx; }
	inline const vector<Density> &getDensity() const { return m_Density0; }
	inline const vector<float> &getCurlMagnitude() const { return m_CurlMagnitude; }
	/*inline const Vector &getDivergence() const { return m_divergence; }
	 inline const Vector &getPressure() const { return m_pressure; }
	 inline const Vector &getVelocity(int i) const { return m_u0[i]; }*/
	inline const vector<int> getSolid() const { return m_Solid; }
	
	/**
	 * Return an interpolated velocity value for the given position
	 */
	inline vec3 getVelocity(vec3 p) const {
		/* Re-scale to [0, gridX] x [0,gridY] x [0,gridZ] */
		p.x /= m_dx; p.y /= m_dx; p.z /= m_dx;
		
		return vec3(
					   getVelocityComponent(p.x, p.y - .5f, p.z - .5f, 0),
					   getVelocityComponent(p.x - .5f, p.y, p.z - .5f, 1),
					   getVelocityComponent(p.x - .5f, p.y - .5f, p.z, 2)
					   );
	}
	
	inline Density getDensity(float x, float y, float z) const {
		return getDensityCatmullRom(x, y, z);
	}
	
private:
	int m_gridX, m_gridY, m_gridZ;
	float m_dx, m_invDx;
	float m_time;
	int m_numVoxels, m_Slice, m_VelocitySlice;
	AABB m_aabb;
	
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
	
	// Temporary vectors for the PCG iteration
	vector<float> m_cgP, m_cgW, m_cgZ, m_cgR, m_cgQ;
	// Incomplete cholesky preconditioner
	vector<float> m_precond;
	
	void constructPressureMatrix();
	void computePrecondition();
	
	void project(float dt);
	void advectVelocity(float dt);
	void advectDensity(float dt);
	void calcForces();
	void integrate(float dt);
	
	float PCG(const vector<float> &b, vector<float> &x, int max_its, float tol);
	void axpy_prod_fast(const vector<float> &x, vector<float> &y) const;
	void precond_solve(const vector<float> &b, vector<float> &x);

	/**
	 * Perform a particle tracing step using the RK2 integrator
	 */
	vec3 traceParticle(const vec3 &p, float dt) const {
		return clip(p, getVelocity(p + getVelocity(p) * (dt * .5f)) * dt);
	}
	
	/**
	 * Return one of the velocity components by performing tri-linear
	 * interpolation over the staggered grid
	 */
	float getVelocityComponent(float x, float y, float z, int c) const {
		const int i = (int) x, j = (int) y, k = (int) z;
		const int pos=i+j*(m_gridX+1)+k*m_VelocitySlice;
		
		/* Return zero for positions outside of the grid */
		if (i < 0 || j < 0 || k < 0 || i >= m_gridX ||
			j >= m_gridY || k >= m_gridZ)
			return 0;
		
		const float alpha = x-i,
		beta = y-j,
		gamma = z-k,
		A1 = m_u0[c][pos],
		B1 = m_u0[c][pos+1],
		C1 = m_u0[c][pos+m_gridX+1],
		D1 = m_u0[c][pos+m_gridX+2];
		const float A2 = m_u0[c][pos+m_VelocitySlice],
		B2 = m_u0[c][pos+m_VelocitySlice+1],
		C2 = m_u0[c][pos+m_VelocitySlice+m_gridX+1],
		D2 = m_u0[c][pos+m_VelocitySlice+m_gridX+2];
		
		return (1-gamma) * ((1-alpha) * (1-beta) * A1 + alpha * (1-beta) * B1
							+ (1-alpha) * beta * C1 + alpha*beta*D1)
		+ gamma * ((1-alpha) * (1-beta) * A2 + alpha * (1-beta) * B2
				   + (1-alpha) * beta * C2 + alpha*beta*D2);
	}
	
	inline Density getDensityCatmullRom(float x, float y, float z) const {
		x = (x / m_dx) - .5f;
		y = (y / m_dx) - .5f;
		z = (z / m_dx) - .5f;
		
		const int i = (int) x, j = (int) y, k = (int) z;
		
		if (i < 0 || j < 0 || k < 0 || i > m_gridX-1 ||
			j > m_gridY-1 || k > m_gridZ-1)
			return Density();
		
		const float alpha = x-i,
		beta = y-j,
		gamma = z-k;
		
		const Density
		A = (z-1>= 0) ?     getDensityCatmullRomY(x, y, z-1, alpha, beta) : Density(),
		B =                 getDensityCatmullRomY(x, y, z,   alpha, beta),
		C = (z+1<m_gridZ) ? getDensityCatmullRomY(x, y, z+1, alpha, beta) : Density(),
		D = (z+2<m_gridZ) ? getDensityCatmullRomY(x, y, z+2, alpha, beta) : Density();
		
		const float gamma2 = gamma*gamma, gamma3 = gamma2*gamma;
		Density d = A*(-0.5f*gamma + gamma2 - 0.5f*gamma3) +
		B*(1.0f - gamma2*(5.0f/2.0f) + gamma3*(3.0f/2.0f)) +
		C*(0.5f*gamma + 2*gamma2 - gamma3*(3.0f/2.0f)) +
		D*(-0.5f*gamma2 + 0.5f*gamma3);
		if (d.density < std::min(B.density, C.density) || d.temp < std::min(B.temp, C.temp)
			|| d.density > std::max(B.density, C.density) || d.temp > std::max(B.temp, C.temp)) {
			return B*(1.0f - gamma) + C*gamma;
		}
		return d;
	}
	
	Density getDensityCatmullRomY(int x, int y, int z, float alpha, float beta) const {
		const Density
		A = (y-1>=0) ?      getDensityCatmullRomX(x, y-1, z, alpha) : Density(),
		B =                 getDensityCatmullRomX(x, y,   z, alpha),
		C = (y+1<m_gridY) ? getDensityCatmullRomX(x, y+1, z, alpha) : Density(),
		D = (y+2<m_gridY) ? getDensityCatmullRomX(x, y+2, z, alpha) : Density();
		
		const float beta2 = beta*beta, beta3 = beta2*beta;
		Density d = A*(-0.5f*beta + beta2 - 0.5f*beta3) +
		B*(1.0f - beta2*(5.0f/2.0f) + beta3*(3.0f/2.0f)) +
		C*(0.5f*beta + 2*beta2 - beta3*(3.0f/2.0f)) +
		D*(-0.5f*beta2 + 0.5f*beta3);
		if (d.density < std::min(B.density, C.density) || d.temp < std::min(B.temp, C.temp)
			|| d.density > std::max(B.density, C.density) || d.temp > std::max(B.temp, C.temp)) {
			return B*(1.0f - beta) + C*beta;
		}
		return d;
	}
	
	inline Density getDensityCatmullRomX(int x, int y, int z, float t) const {
		const int pos=x+y*m_gridX+z*m_Slice;
		const Density A = (x-1 >= 0) ? m_Density0[pos-1] : Density(),
		B = m_Density0[pos],
		C = (x+1<m_gridX) ? m_Density0[pos+1] : Density(),
		D = (x+2<m_gridX) ? m_Density0[pos+2] : Density();
		const float t2 = t*t, t3 = t2*t;
		Density d = A*(-0.5f*t + t2 - 0.5f*t3) +
		B*(1.0f - t2*(5.0f/2.0f) + t3*(3.0f/2.0f)) +
		C*(0.5f*t + 2*t2 - t3*(3.0f/2.0f)) +
		D*(-0.5f*t2 + 0.5f*t3);
		
		/* Switch to trilinear interpolation in the case of an overshoot */
		if (d.density < std::min(B.density, C.density) || d.temp < std::min(B.temp, C.temp)
			|| d.density > std::max(B.density, C.density) || d.temp > std::max(B.temp, C.temp)) {
			return B*(1.0f - t) + C*t;
		}
		return d;
	}
	
	// Utils
	int getIdx(int x, int y, int z);
	
	inline vec3 clip(vec3 p, vec3 v) const {
		Ray3 r(p, v);
		float nearT, farT;
		if (m_aabb.rayIntersect(r, nearT, farT)) {
			if (nearT < 0)
				nearT = farT;
			if (farT < 0)
				std::cout << "Internal error while clipping to AABB!" << std::endl;
			if (farT > 1)
				return p+v;
			std::cout << "Clipping to AABB" << std::endl;
			
			return p + (v*nearT);
		}
		return p+v;
	}


};

#endif /* defined(__Fluid__Solver__) */
