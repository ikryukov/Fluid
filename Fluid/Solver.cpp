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
}

Solver::~Solver()
{
	
}

void Solver::constructPressureMatrix()
{
	for (int x = 0; x < m_gridX; ++x) {
		for (int y = 0; y < m_gridY; ++y) {
			for (int z = 0; z < m_gridZ; ++z) {
				
				int pos_center = getIdx(x, y, z);
				int pos_right  = getIdx(x + 1, y, z);
				int pos_below  = getIdx(x, y + 1, z);
				int pos_behind = getIdx(x, y, z + 1);
				
				bool fluid_center = !m_Solid[pos_center];
				bool fluid_right = (x != m_gridX - 1) && !m_Solid[pos_right];
				bool fluid_below = (y != m_gridY - 1) && !m_Solid[pos_below];
				bool fluid_behind = (z != m_gridZ - 1) && !m_Solid[pos_behind];
				
				if (fluid_center && fluid_right) {
					m_ADiag[pos_center] += 1;
					m_ADiag[pos_right] += 1;
					m_APlusX[pos_center] = -1;
				}
				
				if (fluid_center && fluid_below) {
					m_ADiag[pos_center] += 1;
					m_ADiag[pos_below] += 1;
					m_APlusY[pos_center] = -1;
				}
				
				if (fluid_center && fluid_behind) {
					m_ADiag[pos_center] += 1;
					m_ADiag[pos_behind] += 1;
					m_APlusZ[pos_center] = -1;
				}
			}
		}
	}
}


void Solver::simulate(float dt)
{
	//TODO: calc CFL condition
	project(dt);
	advectVelocity(dt);
	advectDensity(dt);
	calcForces();
	integrate(dt);
	m_time += dt;
}

void Solver::project(float dt)
{
	// Calculate divergence
	for (int x = 0 ; x < m_gridX; ++x) {
		int velIdx = x * m_VelocitySlice;
		int pos = x * m_Slice;
		for (int y = 0; y < m_gridY; ++y) {
			for (int z = 0; z < m_gridZ; ++z) {
				if (m_Solid[pos]) {
					m_Divergence[pos] = 0;
					continue;
				}
				m_Divergence[pos] = (m_u0[0][velIdx + 1] - m_u0[0][velIdx]
									 + m_u0[1][velIdx + m_gridZ + 1] - m_u0[1][velIdx]
									 + m_u0[2][velIdx + m_VelocitySlice] - m_u0[2][velIdx])
				* m_invDx;
				++pos;
				++velIdx;
			}
			++velIdx;
		}
	}
	PCG(m_Divergence * (m_dx * m_dx), m_Pressure, 100, 1e-4);
	
	// Apply the computed gradient
	for (int x = 0; x < m_gridX; ++x)
	{
		int pos = x * m_Slice;
		for (int y = 0; y < m_gridY; ++y)
		{
			for (int z = 0; z < m_gridZ; ++z)
			{
				if (m_Solid[pos])
				{
					continue;
				}
				int velIdx = z + y * (m_gridZ + 1) + x * m_VelocitySlice;
				if (z < m_gridZ - 1 && !m_Solid[pos + 1])
				{
					m_u1[0][velIdx + 1] = m_u0[0][velIdx + 1] + (m_Pressure[pos + 1] - m_Pressure[pos]) * m_invDx;
				}
				if (y < m_gridY - 1 && m_Solid[pos + m_gridZ])
				{
					m_u1[1][velIdx + m_gridZ + 1] = m_u0[1][velIdx + m_gridZ + 1] + (m_Pressure[pos + m_gridZ] - m_Pressure[pos]) * m_invDx;
				}
				if (x < m_gridX - 1 && !m_Solid[pos + m_Slice])
				{
					m_u1[2][velIdx + m_VelocitySlice] = m_u0[2][velIdx + m_VelocitySlice] + (m_Pressure[pos + m_Slice] - m_Pressure[pos]) * m_invDx;
				}
				++pos;
			}
		}
	}
	
	m_u0[0].swap(m_u1[0]);
	m_u0[1].swap(m_u1[1]);
	m_u0[2].swap(m_u1[2]);
}

void Solver::advectVelocity(float dt)
{
	
}

void Solver::advectDensity(float dt)
{
	
}

void Solver::calcForces()
{
	
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