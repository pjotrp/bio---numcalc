/*
 * File Uniform01K.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday September 24 2004
 */

#include "Uniform01K.h"

const long Uniform01K::MAXNUMBER = 1000000000;
const long Uniform01K::ZERO = 0;
const long Uniform01K::MODSEED = 256434901;

//** Class constructor: *******************************************************/
Uniform01K::Uniform01K(long seed) {
	setSeed(seed);
}

//** Class destructor: *******************************************************/
Uniform01K::~Uniform01K() {}

//** Other methodes: *********************************************************/
void Uniform01K::setSeed(long seed) {
	long tmp1, tmp2;
	unsigned short i, ii, j;
	tmp1 = MODSEED - (seed < 0 ? -seed : seed);
	tmp1 %= MAXNUMBER;
	_tab[55] = tmp1;
	tmp2 = 1;

	for (i=1 ; i<=54 ; i++) {
		ii = (21 * i) % 55;
		_tab[ii] = tmp2;
		tmp2 = tmp1 - tmp2;

		if (tmp2 < ZERO) tmp2 += MAXNUMBER;

		tmp1 = _tab[ii];
	}

	for (j=1 ; j<=4 ; j++)
		for (i=1 ; i<=55 ; i++) {
			_tab[i] -= _tab[1 + (i + 30) % 55];
			if (_tab[i] < ZERO) _tab[i] += MAXNUMBER;
		}
	_it1 = 0;
	_it2 = 31;
}

double Uniform01K::drawNumber() const {
	if (++_it1 == 56) _it1 = 1;
	if (++_it2 == 56) _it2 = 1;

	long tmp = _tab[_it1] - _tab[_it2];
	if (tmp < ZERO) tmp += MAXNUMBER;
	_tab[_it2] = tmp;

	return 1. * tmp / MAXNUMBER;
}
