/*
 * File Uniform01QD.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday September 24 2004
 */

#include "Uniform01QD.h"

#include <cmath>

//** Class constructor: *******************************************************/
Uniform01QD::Uniform01QD(long seed) {
	setSeed(seed);
}

//** Class destructor: *******************************************************/
Uniform01QD::~Uniform01QD() {}

//** Other methodes: *********************************************************/
void Uniform01QD::setSeed(long seed) {
	_seed = seed;
}

double Uniform01QD::drawNumber() const {
	_seed = _seed * 1103515245 + 12345;
	if (_seed < 0) _seed = -_seed;
	return _seed / 2147483648.0;
}
