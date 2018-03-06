#ifndef _ZAITox_
#define _ZAITox_

/*!
 \file
 \brief Header file for ZAITox classes.


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
//! Defines the Radiotoxicity of a ZAI
/*!
The aims of this class is to handle radiotoxicity
Activity-to-Sievert factors are from (??).
Values are in CLASS_PATH/data/HeatTox.dat

@author BLG
@author BaM
@version 1.0
*/
//________________________________________________________________________

class ZAITox {
 public:
  /*!
   \name Constructor/Desctructor
   */
  //@{

  ZAITox();  //!< Normal Constructor.

  ~ZAITox();  //!< Normal Destructor.
  //@}

  /*!
   \name Fucntions returning radiotoxicity [Sv]
   */
  //@{
  double GetRadioTox(ZAI zai) const;  //!< get with ZAI
  double GetRadioTox(const int Z, const int A, const int I) const {
    return GetRadioTox(ZAI(Z, A, I));
  }  //!< Get with Z, A

  double GetRadioTox(const IsotopicVector IV) const;  // return Heat of IV [W]
  //@}
 private:
  map<ZAI, double> fZAITox;  //! ZAI Radiotox list
};

#endif
