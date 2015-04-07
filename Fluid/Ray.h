//
//  Ray.h
//  Fluid
//
//  Created by Ilya Kryukov on 27.03.15.
//  Copyright (c) 2015 Ilya Kryukov. All rights reserved.
//

#ifndef Fluid_Ray_h
#define Fluid_Ray_h

#include "glm/glm.hpp"

const float Epsilon = 1e-7;

using namespace glm;

/** \brief Simple three-dimensional ray class with
 minimum / maximum extent information */
class Ray3 {
public:
	/// Ray origin
	vec3 o;
	/// Minimum range for intersection tests
	mutable float mint;
	/// Ray direction
	vec3 d;
	/// Maximum range for intersection tests
	mutable float maxt;
	
	/// Construct a new ray
	Ray3() : mint(Epsilon), maxt(std::numeric_limits<float>::max()) {
	}
	
	/// Construct a new ray
	Ray3(vec3 _o, vec3 _d)
	: o(_o), mint(Epsilon),  d(_d), maxt(std::numeric_limits<float>::max()) {
	}
	
	/// Return 3d coordinates of a point on the ray
	vec3 operator() (float t) const { return o + d*t; }
	
	/// Return a string representation of this ray
	std::string toString() const {
		std::ostringstream oss;
		//oss << "Ray3[orig=" << o.toString() << ", dest=" << d.toString() << "]";
		return oss.str();
	}
};

#endif
