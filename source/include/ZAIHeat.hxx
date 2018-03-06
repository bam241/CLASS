#ifndef _ZAIHeat_
#define _ZAIHeat_

/*!
 \file
 \brief Header file for ZAIHeat classes.


 @author BaM
 @author BLG
 @version 2.0
 */

#include <map>

#include <iostream>
#include "TObject.h"
#include "ZAI.hxx"

using namespace std;

class IsotopicVector;

//-----------------------------------------------------------------------------//
//! Defines the decay heat of a ZAI
/*!
The aims of this class is to handle decay heat
Activity-to-Heat factors are from (92SDOe + M.e.Brendan, dec 98).
Values are in CLASS_PATH/data/HeatTox.dat

@author BLG
@author BaM
@version 1.0
*/
//________________________________________________________________________

class ZAIHeat {
 public:
  /*!
   \name Constructor/Desctructor
   */
  //@{

  ZAIHeat();  //!< Normal Constructor.

  ~ZAIHeat();  //!< Normal Destructor.
  //@}

  /*!
   \name Fucntions returning decay Heat [W]
   */
  //@{
  double GetHeat(ZAI zai) const;  //!< get with ZAI
  double GetHeat(const int Z, const int A, const int I) const {
    return GetHeat(ZAI(Z, A, I));
  }  //!< Get with Z, A

  double GetHeat(const IsotopicVector IV) const;  // return Heat of IV [W]
  //@}
 private:
  map<ZAI, double> fZAIHeat;  //! ZAI Heat list
};

#endif
