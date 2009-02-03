//
// File: BasicVectorTools.h
// Created by: Sylvain Gaillard
// Created on: Tue Jan 14:58:50 CET 2009
//

/*
Copyright or © or Copr. CNRS, (January 13, 2009)

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

#ifndef _NUMCALCAPPLICATIONTOOLS_H_
#define _NUMCALCAPPLICATIONTOOLS_H_

#include <Utils/StringTokenizer.h>
#include "VectorTools.h"
using namespace bpp;

namespace bpp {
  class NumCalcApplicationTools {
    public:
      NumCalcApplicationTools();
      virtual ~NumCalcApplicationTools();

    public:
      /**
       * @brief Build a vector of integers as describe by a string
       *
       * Build a vector of integers following a description like:
       * "2, 5, 7-10, 4" &rarr; [2, 5, 7, 8, 9, 10, 4]
       *
       * @param s The string to parse.
       * @param delim Delimiter between elements.
       * @param seqdelim Delimiter between min and max for a sequence.
       * @return A vector containing the integers
       */
      static vector<int> seqFromString(const string & s, const string & delim = ",", const string & seqdelim = "-") {
        vector<int> seq;
        StringTokenizer * st = new StringTokenizer(s, delim, true);
        while (st->hasMoreToken()) {
          StringTokenizer * st2 = new StringTokenizer(st->nextToken(), seqdelim, true);
          if (st2->numberOfRemainingTokens() > 1) {
            vector<int> tmp = VectorTools::seq(TextTools::toInt(st2->getToken(0)), TextTools::toInt(st2->getToken(1)), 1);
            VectorTools::append(seq, tmp);
          } else {
            seq.push_back(TextTools::toInt(st2->getToken(0))); }  
        }  
        return seq;
      }

  };
}

#endif  //_NUMCALCAPPLICATIONTOOLS_H_