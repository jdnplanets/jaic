/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
Random number generator class.
*/

 // Rnd.h

#pragma once
#include <random>

/*object for sampling random numbers*/
class Rnd {
    public:
        //constructor: set initial random seed and distribution limits
        Rnd(): mt_gen{std::random_device()()}, rnd_dist{0,1.0} {}
        float operator() () {return rnd_dist(mt_gen);}
    
    protected:
        std::mt19937 mt_gen;	    //random number generator
        std::uniform_real_distribution<float> rnd_dist;  //uniform distribution
    };
    
    extern Rnd rnd;		//tell the compiler that an object of type Rnd called rnd is defined somewhere

    