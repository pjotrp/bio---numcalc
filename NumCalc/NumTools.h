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

#ifndef _NUMTOOLS_H_
#define _NUMTOOLS_H_

#include "Functions.h"

/**
 * @brief Some utilitary function for numerical calculus.
 */
class NumTools
{
public:

  /**
   * @brief Get the magnitude of a value.
   *
   * This template function may work with any type for which the operators
   * < and - are defined.
   *
   * @param a The value for which the magnitude must be returned.
   * @return The magnitude of the value a.
   */ 
  template<class T> static T abs(T a) { return a < 0 ? -a : a; }

  /**
   * @brief Get the sign of a value.
   *
   * This template function may work with any type for which the operators
   * < and == are defined.
   *
   * @param a The value for which the sign must be returned.
   * @return -1 if a < 0, 0 if a = 0, 1 else.
   */ 
  template<class T> static T sign(T a) { return a < 0 ? -1 : (a == 0 ? 0 : 1); }

  /**
   * @brief Get the max between 2 values.
   *
   * This template function may work with any type for which the operator
   * > is defined.
   *
   * @param a, b The two values to compare.
   * @return a if a > b, b else.
   */ 
  template<class T> static T max(T a, T b) { return a > b ? a : b; }

  /**
   * @brief Get the min between 2 values.
   *
   * This template function may work with any type for which the operator
   * < is defined.
   *
   * @param a, b The two values to compare.
   * @return a if a < b, b else.
   */ 
  template<class T> static T min(T a, T b) { return a < b ? a : b; }

  /**
   * @brief Get the magnitude of a times the sign of b.
   *
   * @param a The value whose magnitude must be used.
   * @param b The value whose sign must be used.
   * @return abs<T>(a) * sign<T>(b).
   */	 
  template<class T> static T sign(T a, T b) { return abs<T>(a) * sign<T>(b); }

  /**
   * @brief Get the square of a number.
   *
   * @param a The value.
   * @return @f$ a^2 @f$.
   */ 
  template<class T> static T sqr(T a) { return a * a; }

  /**
   * @brief Get a number to a given power.
   *
   * @param a The value.
   * @param b The power
   * @return @f$ a^b @f$.
   */ 
  template<class T> static T pow(T a, T b) { return exp(b * log(a)); }

  /**************************************************************************/

  template<class T> static void swap(T & a, T & b)
  {
	  T swap = a;
	  a = b;
	  b = swap;	
  }

  template<class T> static void shift(T & a, T & b, T c)
  {
  	a = b; b = c;
  }

  template<class T> static void shift(T & a, T & b, T & c, T d)
  {
  	a = b; b = c; c = d;
  }

  /**************************************************************************/

  template<class T> static T fact(T n) { return (n == 0) ? 1 : n * fact(n - 1); }

  /**************************************************************************/
  
  /**
   * @Brief Find one root of the given function.
   *
   * @param f The function to analyse.
   * @param param The name of the parameter to solve.
   * @param a Lower bound of initial interval.
   * @param b Upper bound of initial interval.
   * @param tolerance The final precision requested.
   * @return The value of the parameter for which the function is zero.
   * @throw Exception If something bad happened or if the initial interval do not contains a root.
   */
  static double uniRoot(Function & f, const string & param, double a, double b, double tolerance) throw (Exception);
  
  /**************************************************************************/

};

#endif	//_NUMTOOLS_H_

