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
: m_gridX(gridX), m_gridY(gridY), m_gridZ(gridZ), m_dx(dx), m_time(0.0f)
{
	m_aabb = AABB(vec3(0,0,0), vec3(gridX,gridY,gridZ));
	int numVelocity = (m_gridX + 1) * (m_gridY + 1) * (m_gridZ + 1);
	m_numVoxels = m_gridX * m_gridY * m_gridZ;
	m_Slice = gridX * gridY;
	m_VelocitySlice = (gridX + 1) * (gridY + 1);
	m_invDx = 1.0f / m_dx;
	
	m_precond.resize(m_numVoxels);
	
	m_Density0.resize(m_numVoxels);
	m_Density1.resize(m_numVoxels);
	
	m_Curl.resize(m_numVoxels);
	m_CurlMagnitude.resize(m_numVoxels);
	m_vcForce.resize(m_numVoxels);
	m_Divergence.resize(m_numVoxels);
	m_Pressure.resize(m_numVoxels);
	m_Solid.resize(m_numVoxels);
	
	m_ADiag.resize(m_numVoxels);
	m_APlusX.resize(m_numVoxels);
	m_APlusY.resize(m_numVoxels);
	m_APlusZ.resize(m_numVoxels);
	
	m_cgW.resize(m_numVoxels);
	m_cgP.resize(m_numVoxels);
	m_cgZ.resize(m_numVoxels);
	m_cgR.resize(m_numVoxels);
	m_cgQ.resize(m_numVoxels);
	
	for (int x = 0; x < m_gridX; ++x)
	{
		for (int y = 0; y < m_gridY; ++y)
		{
			for (int z = 0; z < m_gridZ; ++z)
			{
				int pos = getIdx(x, y, z);
				m_Solid[pos] = 0;
			}
		}
	}
	
	for (int i = 0; i < 3; ++i)
	{
		m_F[i].resize(numVelocity);
		m_u0[i].resize(numVelocity);
		m_u1[i].resize(numVelocity);
	}
	
	constructPressureMatrix();
	// TODO: calc precondition
	computePrecondition();
}

Solver::~Solver()
{
	
}

void Solver::constructPressureMatrix()
{
	std::cout << "FluidSolver: Constructing pressure matrix" << std::endl;
	for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				bool fluid = !m_Solid[pos];
				bool fluidRight = (x != m_gridX-1) && !m_Solid[pos+1];
				bool fluidBelow = (y != m_gridY-1) && !m_Solid[pos+m_gridX];
				bool fluidBehind = (z != m_gridZ-1) && !m_Solid[pos+m_Slice];
				
				if (fluid && fluidRight) {
					m_ADiag[pos] += 1;
					m_ADiag[pos+1] += 1;
					m_APlusX[pos] = -1;
				}
				
				if (fluid && fluidBelow) {
					m_ADiag[pos] += 1;
					m_ADiag[pos+m_gridX] += 1;
					m_APlusY[pos] = -1;
				}
				
				if (fluid && fluidBehind) {
					m_ADiag[pos] += 1;
					m_ADiag[pos+m_Slice] += 1;
					m_APlusZ[pos] = -1;
				}
			}
		}
	}
}

void Solver::computePrecondition()
{
	std::cout << "FluidSolver: Computing the MIC preconditioner" << std::endl;
	const float rho = .25f, tau = 0.97f;
	for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (m_Solid[pos])
					continue;
				float termLeft = 0.0f, termAbove = 0.0f, termFront = 0.0f,
				termLeft2 = 0.0f, termAbove2 = 0.0f, termFront2 = 0.0f;
				if (x > 0 && !m_Solid[pos-1]) {
					termLeft  = m_APlusX[pos-1] * m_precond[pos-1];
					termLeft2 = m_APlusX[pos-1] * (m_APlusY[pos-1] + m_APlusZ[pos-1]) * m_precond[pos-1] * m_precond[pos-1];
				}
				
				if (y > 0 && !m_Solid[pos-m_gridX]) {
					termAbove  = m_APlusY[pos-m_gridX] * m_precond[pos-m_gridX];
					termAbove2 = m_APlusY[pos-m_gridX] * (m_APlusX[pos-m_gridX] + m_APlusZ[pos-m_gridX]) * m_precond[pos-m_gridX] * m_precond[pos-m_gridX];
				}
				
				if (z > 0 && !m_Solid[pos-m_Slice]) {
					termFront  = m_APlusZ[pos-m_Slice] * m_precond[pos-m_Slice];
					termFront2 = m_APlusZ[pos-m_Slice] * (m_APlusX[pos-m_Slice] + m_APlusY[pos-m_Slice]) * m_precond[pos-m_Slice] * m_precond[pos-m_Slice];
				}
				
				float e = m_ADiag[pos] - termLeft*termLeft - termAbove*termAbove
				- termFront*termFront - tau * (termLeft2 + termAbove2 + termFront2);
				if (e < rho * m_ADiag[pos])
					e = m_ADiag[pos];
				m_precond[pos] = 1.0f / std::sqrt(e);
			}
		}
	}
}

void Solver::simulate(float dt)
{
	std::cout << "solver::simulate " << m_time << std::endl;
	//TODO: calc CFL condition
	project(dt);
	advectVelocity(dt);
	advectDensity(dt);
	calcForces();
	
	float targetTemp = 273+8.4, rateDensity = 10.0f;
	
	int rad = 4;
	for (int i=-rad; i<=rad; ++i) {
		for (int j=-rad; j<=rad; ++j) {
			int x=m_gridX-1, y = m_gridY-1 - 6+i, z = m_gridZ/2+j;
			float tmp = rad*rad - i*i-j*j;
			if (tmp < 0)
				continue;
			int dist = -std::max((int) std::sqrt(tmp)-1, 0);
			m_Density0[dist + x + y * m_gridX + z * m_Slice].density = dt*rateDensity;
			m_Density0[dist + x + y * m_gridX + z * m_Slice].temp = targetTemp;
			m_F[0][dist + x + y * (m_gridX+1) + z * m_VelocitySlice] = -3.8f;
		}
	}
	
	int z;
	{
		for (z=0; z<m_gridZ; ++z) {
			int velIdx = z * m_VelocitySlice;
			for (int y=0; y<m_gridY; ++y, ++velIdx) {
				for (int x=0; x<m_gridX; ++x, ++velIdx) {
					m_u0[0][velIdx] += dt*m_F[0][velIdx];
					m_u0[1][velIdx] += dt*m_F[1][velIdx];
					m_u0[2][velIdx] += dt*m_F[2][velIdx];
				}
			}
		}
	}
	
	integrate(dt);
	m_time += dt;
}

void Solver::project(float dt)
{
	/* Calculate the divergence */
#pragma omp parallel
	{
#pragma omp for schedule(dynamic, CHUNKSIZE)
		for (int z=0; z<m_gridZ; ++z) {
			int velIdx = z * m_VelocitySlice, pos = z * m_Slice;
			for (int y=0; y<m_gridY; ++y, ++velIdx) {
				for (int x=0; x<m_gridX; ++x, ++pos, ++velIdx) {
					if (m_Solid[pos]) {
						m_Divergence[pos] = 0;
						continue;
					}
					m_Divergence[pos] =
					(m_u0[0][velIdx+1]           - m_u0[0][velIdx]
					 + m_u0[1][velIdx+m_gridX+1]   - m_u0[1][velIdx]
					 + m_u0[2][velIdx+m_VelocitySlice]  - m_u0[2][velIdx])
					* m_invDx;
				}
			}
		}
	}
	
	PCG(m_Divergence * (m_dx*m_dx), m_Pressure, 100, 1e-4);
	
	/* Apply the computed gradients */
#pragma omp parallel
	{
#pragma omp for schedule(dynamic, CHUNKSIZE)
		for (int z=0; z<m_gridZ; ++z) {
			int pos = z * m_Slice;
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x, ++pos) {
					if (m_Solid[pos])
						continue;
					int velIdx = x + y * (m_gridX+1) + z * m_VelocitySlice;
					if (x < m_gridX-1 && !m_Solid[pos + 1])
						m_u1[0][velIdx+1] = m_u0[0][velIdx+1] +
						(m_Pressure[pos+1] - m_Pressure[pos]) * m_invDx;
					if (y < m_gridY-1 && !m_Solid[pos + m_gridX])
						m_u1[1][velIdx+m_gridX+1] = m_u0[1][velIdx+m_gridX+1] +
						(m_Pressure[pos+m_gridX] - m_Pressure[pos]) * m_invDx;
					if (z < m_gridZ-1 && !m_Solid[pos + m_Slice])
						m_u1[2][velIdx+m_VelocitySlice] = m_u0[2][velIdx+m_VelocitySlice] +
						(m_Pressure[pos+m_Slice] - m_Pressure[pos]) * m_invDx;
				}
			}
		}
	}
	
	m_u0[0].swap(m_u1[0]);
	m_u0[1].swap(m_u1[1]);
	m_u0[2].swap(m_u1[2]);
}

void Solver::advectVelocity(float dt)
{
	int z;
	
	/* Advect the velocity field */
	{
		for (z=0; z<m_gridZ; ++z) {
			int pos = z*m_Slice;
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x, ++pos) {
					int velIdx = x + y * (m_gridX+1) + z * m_VelocitySlice;
					if (m_Solid[pos])
						continue;
					
					/* Advect X velocities */
					if (x < m_gridX-1 && !m_Solid[pos + 1]) {
						vec3 p((x+1.0f)*m_dx, (y+.5f)*m_dx, (z+0.5f)*m_dx);
						vec3 p2 = traceParticle(p, -dt);
						m_u1[0][velIdx+1] = getVelocity(p2).x;
					}
					
					/* Advect Y velocities */
					if (y < m_gridY-1 && !m_Solid[pos + m_gridX]) {
						vec3 p((x+0.5f)*m_dx, (y+1.0f)*m_dx, (z+0.5f)*m_dx);
						vec3 p2 = traceParticle(p, -dt);
						m_u1[1][velIdx+m_gridX+1] = getVelocity(p2).y;
					}
					
					/* Advect Z velocities */
					if (z < m_gridZ-1 && !m_Solid[pos + m_Slice]) {
						vec3 p((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+1.0f)*m_dx);
						vec3 p2 = traceParticle(p, -dt);
						m_u1[2][velIdx+m_VelocitySlice] = getVelocity(p2).z;
					}
				}
			}
		}
	}
	
	m_u0[0].swap(m_u1[0]);
	m_u0[1].swap(m_u1[1]);
	m_u0[2].swap(m_u1[2]);
}

void Solver::advectDensity(float dt)
{
	/* Advect the density field */
	int z;
	
	{
		for (z=0; z<m_gridZ; ++z) {
			int pos = z*m_Slice;
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x, ++pos) {
					if (pos >= m_numVoxels) {
						break;
					}
					if (m_Solid[pos])
						continue;
					vec3 p((x+.5f)*m_dx, (y+.5f)*m_dx, (z+.5f)*m_dx);
					vec3 p2 = traceParticle(p, -dt);
					m_Density1[pos] = getDensity(p2.x, p2.y, p2.z);
				}
			}
		}
	}
	m_Density0.swap(m_Density1);
}

void Solver::calcForces()
{
	int z;
	std::fill(m_F[0].begin(), m_F[0].end(), 0);
	std::fill(m_F[1].begin(), m_F[1].end(), 0);
	std::fill(m_F[2].begin(), m_F[2].end(), 0);
	
	/* Calculate the magnitude of the curl */
	{
		for (z=0; z<m_gridZ; ++z) {
			int pos = z * m_Slice;
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x, ++pos) {
					if (m_Solid[pos])
						continue;
					vec3 velTop     = getVelocity(vec3((x+0.5f)*m_dx, (y+0.0f)*m_dx, (z+0.5f)*m_dx));
					vec3 velBot     = getVelocity(vec3((x+0.5f)*m_dx, (y+1.0f)*m_dx, (z+0.5f)*m_dx));
					vec3 velLeft    = getVelocity(vec3((x+0.0f)*m_dx, (y+0.5f)*m_dx, (z+0.5f)*m_dx));
					vec3 velRight   = getVelocity(vec3((x+1.0f)*m_dx, (y+0.5f)*m_dx, (z+0.5f)*m_dx));
					vec3 velFront   = getVelocity(vec3((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+0.0f)*m_dx));
					vec3 velBack    = getVelocity(vec3((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+1.0f)*m_dx));
					
					vec3 dudx = (velRight - velLeft) / m_dx;
					vec3 dudy = (velBot - velTop) / m_dx;
					vec3 dudz = (velBack - velFront) / m_dx;
					vec3 curl(dudy.z - dudz.y, dudz.x - dudx.z, dudx.y-dudy.x);
					m_Curl[pos] = curl;
					m_CurlMagnitude[pos] = curl.length();
				}
			}
		}
	}
	
	/* Add the vorticity confinement force */
	{
		for (z=0; z<m_gridZ; ++z) {
			int pos = z * m_Slice;
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x, ++pos) {
					if (y==0 || x==0 || z == 0 || x==m_gridX-1
						|| y==m_gridY-1 || z == m_gridZ-1 || m_Solid[pos]) {
						m_Curl[pos] = vec3(0,0,0);
						continue;
					}
					vec3 N(
							  m_CurlMagnitude[pos+1]- m_CurlMagnitude[pos-1],
							  m_CurlMagnitude[pos+m_gridX]- m_CurlMagnitude[pos-m_gridX],
							  m_CurlMagnitude[pos+m_Slice]- m_CurlMagnitude[pos-m_Slice]
							  );
					float length = N.length();
					if (length < Epsilon)
						continue;
					m_vcForce[pos] = cross(N / length, m_Curl[pos]) * (m_dx * 0.8f);
				}
			}
		}
	}
	
	/* Add the buoyancy force */
	float ambient = 273.0f;
	/*	int nCells = 0;
	 for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
	 for (int x=0; x<m_gridX; ++x, ++pos) {
	 if (m_solid[pos])
	 continue;
	 ambient += m_d0[pos].temp;
	 if (m_d0[pos].temp != 273)
	 ++nCells;
	 }
		}
	 }
	 ambient /= nCells;
	 cout << ambient << endl;
	 */
	const float a = 0.0625f*0.5f, b=0.025f;
	{
		for (z=0; z<m_gridZ; ++z) {
			int pos = z * m_Slice;
			int velIdx = z * m_VelocitySlice;
			for (int y=0; y<m_gridY; ++y, ++velIdx) {
				for (int x=0; x<m_gridX; ++x, ++velIdx, ++pos) {
					if (m_Solid[pos])
						continue;
					
					const Density top   = getDensity((x+.5f)*m_dx, y*m_dx, (z+.5f)*m_dx);
					
					if (x != 0 && !m_Solid[pos-1])
						m_F[0][velIdx] += (m_vcForce[pos].x + m_vcForce[pos-1].x)/2.0f;
					
					if (y!= 0 && !m_Solid[pos-m_gridX]) {
						m_F[1][velIdx] -= -a*top.density + b*(top.temp-ambient);
						m_F[1][velIdx] += (m_vcForce[pos].y + m_vcForce[pos-m_gridX].y)/2.0f;
					}
					
					if (z != 0 && !m_Solid[pos-m_Slice])
						m_F[2][velIdx] += (m_vcForce[pos].z + m_vcForce[pos-m_Slice].z)/2.0f;
				}
			}
		}
	}

}

void Solver::integrate(float dt)
{
	
}

float Solver::PCG(const vector<float> &b, vector<float> &x, int max_its, float tol)
{
	float beta=0, lastRho=0, rho=0, tolSquared = tol*tol;
	int k=0;
	
	axpy_prod_fast(x, m_cgR);
	m_cgR = b - m_cgR;
	precond_solve(m_cgR, m_cgZ);
	rho = inner_prod(m_cgR, m_cgZ);
	
	while (k<max_its && rho > tolSquared) {
		if (k == 0) {
			noalias(m_cgP) = m_cgZ;
		} else {
			beta = rho / lastRho;
			m_cgP = m_cgZ + m_cgP * beta;
		}
		
		axpy_prod_fast(m_cgP, m_cgW);
		const float alpha = rho / inner_prod(m_cgP, m_cgW);
		noalias(x) += alpha * m_cgP;
		noalias(m_cgR) -= alpha * m_cgW;
		
		precond_solve(m_cgR, m_cgZ);
		lastRho = rho;
		rho = inner_prod(m_cgR, m_cgZ);
		k++;
	}
	//cout << "  + PCG: Residual after " << k << " iterations : " << std::sqrt(rho) << endl;
	
	return std::sqrt(rho);
}

void Solver::axpy_prod_fast(const vector<float> &x, vector<float> &y) const {
	int i;
	for (i = 0; i < m_Slice; ++i) {
		float result = m_ADiag[i] * x[i]
		+ m_APlusX[i] * x[i+1]
		+ m_APlusY[i] * x[i+m_gridX]
		+ m_APlusZ[i] * x[i+m_Slice];
		if (i-1 >= 0)
			result += m_APlusX[i-1] * x[i-1];
		if (i-m_gridX >= 0)
			result += m_APlusY[i-m_gridX] * x[i-m_gridX];
		if (i-m_Slice >= 0)
			result += m_APlusZ[i-m_Slice] * x[i-m_Slice];
		
		y[i] = result;
	}
	for (i=m_numVoxels-m_Slice; i<m_numVoxels; ++i) {
		float result = m_ADiag[i] * x[i]
		+ m_APlusX[i-1] * x[i-1]
		+ m_APlusY[i-m_gridX] * x[i-m_gridX]
		+ m_APlusZ[i-m_Slice] * x[i-m_Slice];
		
		if (i+1 < m_numVoxels)
			result += m_APlusX[i] * x[i+1];
		if (i+m_gridX < m_numVoxels)
			result += m_APlusY[i] * x[i+m_gridX];
		if (i+m_Slice < m_numVoxels)
			result += m_APlusZ[i] * x[i+m_Slice];
		
		y[i] = result;
	}
	
	int end = m_numVoxels-m_Slice;

		for (i=m_Slice; i<end; ++i) {
			y[i] = m_ADiag[i] * x[i]
			+ m_APlusX[i] * x[i+1]
			+ m_APlusY[i] * x[i+m_gridX]
			+ m_APlusZ[i] * x[i+m_Slice]
			+ m_APlusX[i-1] * x[i-1]
			+ m_APlusY[i-m_gridX] * x[i-m_gridX]
			+ m_APlusZ[i-m_Slice] * x[i-m_Slice];
		}
}

void Solver::precond_solve(const vector<float> &b, vector<float> &x)
{
	/* Solve the lower triangular system Lq = b */
	for (int i=0; i<m_numVoxels; ++i) {
		if (m_Solid[i])
			continue;
		float temp = b[i];
		if (i > 0)
			temp -= m_APlusX[i-1] * m_precond[i-1] * m_cgQ[i-1];
		if (i > m_gridX)
			temp -= m_APlusY[i-m_gridX] * m_precond[i-m_gridX] * m_cgQ[i-m_gridX];
		if (i > m_Slice)
			temp -= m_APlusZ[i-m_Slice] * m_precond[i-m_Slice] * m_cgQ[i-m_Slice];
		m_cgQ[i] = temp * m_precond[i];
	}
	
	/* Solve the upper triangular system L^Tx = q */
	for (int i=m_numVoxels-1; i>=0; --i) {
		if (m_Solid[i])
			continue;
		float temp = m_cgQ[i];
		if (i+1 < m_numVoxels)
			temp -= m_APlusX[i] * m_precond[i] * x[i+1];
		if (i+m_gridX < m_numVoxels)
			temp -= m_APlusY[i] * m_precond[i] * x[i+m_gridX];
		if (i+m_Slice < m_numVoxels)
			temp -= m_APlusZ[i] * m_precond[i] * x[i+m_Slice];
		x[i] = temp * m_precond[i];
	}
}

int Solver::getIdx(int x, int y, int z)
{
	return z + (y + (x * m_gridY)) * m_gridZ;
}