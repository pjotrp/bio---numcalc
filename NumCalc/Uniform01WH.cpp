/*
 * File Uniform01WH.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday September 24 2004
 */

#include "Uniform01WH.h"

#include <cmath>

//** Class constructor: *******************************************************/
Uniform01WH::Uniform01WH(long seed) {
	setSeeds(seed);
}

//** Class destructor: *******************************************************/
Uniform01WH::~Uniform01WH() {}

//** Other methodes: *********************************************************/
void Uniform01WH::setSeed(long seed) {
	setSeeds(seed);
}

void Uniform01WH::setSeeds(long seed1, long seed2, long seed3) {
	ix = seed1;
	iy = seed2;
	iz = seed3;
}

double Uniform01WH::drawNumber() const {
	double i;
	ix = (171 * ix) % 30269;
	iy = (172 * iy) % 30307;
	iz = (170 * iz) % 30323;
	return modf((double) ix / 30269. + (double) iy / 30307. + (double) iz / 30323, &i);
}
