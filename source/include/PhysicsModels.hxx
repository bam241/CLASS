
#ifndef _PhysicsModels_
#define _PhysicsModels_

/*!
 \file
 \brief Header file for XS_INTERPOLATOR class.


 @author BaM
 @author BLG
 @version 1.0
 */
#include "EquivalenceModel.hxx"
#include "EvolutionData.hxx"
#include "IrradiationModel.hxx"
#include "XSModel.hxx"

using namespace std;
typedef long long int cSecond;

//-----------------------------------------------------------------------------//
//! Container object of XSModel, EquivalenceModel and IrradiationModel

/*!
 Define a contener of all physics models used for a specific couple
(reactor,fuel).

These class aim is basicaly to store 3 differents physics model :

 The 2 following are data base related (for one Reactor and one fuel type ) :
User can either define his own (see manual) or uses the provided ones  :
\li XSModel : Mean cross section prediction (Closest, MLP )
\li EquivalenceModel : Fissile content prediction ( Linear,Quadratique, MLP ,
Baker & Ross, ...)

 This one is bateman solvers related :
\li IrradiationModel : can be Runge Kutta 4 or Matrix


 @author BaM
 @author BLG
 @version 1.0
 */
//________________________________________________________________________

class PhysicsModels : public CLASSObject {
 public:
  /*!
   \name Constructor/Desctructor
   */
  //@{

  PhysicsModels();  //!< Default Constructor

  //{
  /// XS, EM, IM Contructor
  /*!
   \param XS : The XSModel (Mean cross section predictor)
   \param EM : The EquivalenceModel (Fissile content predictor)
   \param IM : The IrradiationModel (Bateman solver)
   */
  PhysicsModels(XSModel* XS, EquivalenceModel* EM, IrradiationModel* IM);
  //}

  //{
  /// CLASSLogger Contructor
  /*!
   \param log : The CLASSLogger
   \param XS : The XSModel (Mean cross section predictor)
   \param EM : The EquivalenceModel (Fissile content predictor)
   \param IM : The IrradiationModel (Bateman solver)
   */
  PhysicsModels(CLASSLogger* log, XSModel* XS, EquivalenceModel* EM,
                IrradiationModel* IM);
  //}

  ~PhysicsModels() { ; }
  //@}

  //{

  //{
  /// GenerateEvolutionData
  /*!
   Call the 3 Physics models to compute the depletion calculation for the right
   fresh fuel with
   the right mean cross sections
   \param IV : The fresh fuel composition
   \param cycletime : The irradiation time [s]
   \param Power : The thermal (as always in CLASS) Power [W]
   */
  EvolutionData GenerateEvolutionData(IsotopicVector IV, double cycletime,
                                      double Power);
  //}

  XSModel* GetXSM() { return fXSM; }  ///< return the mean cross section predictor

  EquivalenceModel* GetEQM() { return fEQM; }  ///< returnthe EQM

  IrradiationModel* GetIM() { return fIM; }  ///< return the Bateman solver

  PhysicsModels* GetPhysicsModels() { return this; }  ///< return this
  
  void SetXSM(XSModel* XSM) { fXSM = XSM; }
  void SetEQM(EquivalenceModel* EQM) { fEQM = EQM;}
  void SetIM(IrradiationModel* IM)  { fIM = IM; } 

 private:
  XSModel* fXSM;  ///< The XSModel (mean cross sections prediction)
  EquivalenceModel*
      fEQM;  ///< The EquivalenceModel (fresh fissile content prediction)
  IrradiationModel*
      fIM;  ///< The IrradiationModel (The Bateman's solver)
};

#endif
