//
// File Uniform01WH.cpp
// Author : Sylvain Gaillard
// Last modification : Friday September 24 2004
//

/*
Copyright or � or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
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
