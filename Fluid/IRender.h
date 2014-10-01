//
//  IRender.h
//  Fluid
//
//  Created by Ilya Kryukov on 01.10.14.
//  Copyright (c) 2014 Ilya Kryukov. All rights reserved.
//

#ifndef Fluid_IRender_h
#define Fluid_IRender_h

class IRender {

public:
	virtual void init(int width, int height) = 0;
	virtual void render() = 0;
};

#endif
