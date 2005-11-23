//
// File: NumTools.h
// Created by: Julien Dutheil
// Created on: Mon Nov 10 12:06:55 2003
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

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

#include "NumTools.h"

/******************************************************************************/

double NumTools::abs(double a) { return abs<double>(a); }
float  NumTools::abs(float  a) { return abs<float >(a); }
long   NumTools::abs(long   a) { return abs<long  >(a); }
int    NumTools::abs(int    a) { return abs<int   >(a); }
short  NumTools::abs(short  a) { return abs<short >(a); }

/******************************************************************************/

double NumTools::sign(double a) { return sign<double>(a); }
float  NumTools::sign(float  a) { return sign<float >(a); }
long   NumTools::sign(long   a) { return sign<long  >(a); }
int    NumTools::sign(int    a) { return sign<int   >(a); }
short  NumTools::sign(short  a) { return sign<short >(a); }

/******************************************************************************/

double NumTools::max(double a, double b) { return max<double>(a, b); }
float  NumTools::max(float  a, float  b) { return max<float >(a, b); }
long   NumTools::max(long   a, long   b) { return max<long  >(a, b); }
int    NumTools::max(int    a, int    b) { return max<int   >(a, b); }
short  NumTools::max(short  a, short  b) { return max<short >(a, b); }

/******************************************************************************/

double NumTools::min(double a, double b) { return min<double>(a, b); }
float  NumTools::min(float  a, float  b) { return min<float >(a, b); }
long   NumTools::min(long   a, long   b) { return min<long  >(a, b); }
int    NumTools::min(int    a, int    b) { return min<int   >(a, b); }
short  NumTools::min(short  a, short  b) { return min<short >(a, b); }

/******************************************************************************/

double NumTools::sign(double a, double b) { return sign<double>(a, b); }
float  NumTools::sign(float  a, float  b) { return sign<float >(a, b); }
long   NumTools::sign(long   a, long   b) { return sign<long  >(a, b); }
int    NumTools::sign(int    a, int    b) { return sign<int   >(a, b); }
short  NumTools::sign(short  a, short  b) { return sign<short >(a, b); }

/******************************************************************************/

double NumTools::sqr(double a) { return sqr<double>(a); }
float  NumTools::sqr(float  a) { return sqr<float >(a); }
long   NumTools::sqr(long   a) { return sqr<long  >(a); }
int    NumTools::sqr(int    a) { return sqr<int   >(a); }
short  NumTools::sqr(short  a) { return sqr<short >(a); }

/******************************************************************************/

