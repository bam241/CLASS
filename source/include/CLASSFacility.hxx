
#ifndef _CLASSFACILITY_
#define _CLASSFACILITY_

#include <fstream>
#include <string>

#include "CLASSObject.hxx"
#include "IsotopicVector.hxx"

#include "TNamed.h"

using namespace std;

class Scenario;
//-----------------------------------------------------------------------------//
//! Defines the common properties of all facilities

/*!
 Define a CLASS Facility.
 The aim of these class is to gather all the commom properties of the
 facilities.

 @author BaM
 */
//________________________________________________________________________

class CLASSFacility : public CLASSObject {
 public:
  /*!
   \name Constructor/Desctructor
   */
  //@{

  //{
  /// Normal Constructor.
  /*!
  Make a new Facility
  \param type identification type of the facility :
  \li 4 Reactor,
  \li 8 Pool,
  \li 16 FabricationPlant.
   */

  CLASSFacility(int type = 0);
  //}

  //{
  /// Special Constructor.
  /*!
   Make a new Facility
   \param log : used for the log.
   \param type identification type of the facility :
   \li 4 Reactor,
   \li 8 Pool,
   \li 16 FabricationPlant.

   */
  CLASSFacility(CLASSLogger* log, int type = 0);
  //}

  //{
  /// Special Constructor.
  /*!
   Make a new Facility
   \param log : used for the log.
   \param cycletime duration of the cycle [s],
   \param type identification type of the facility :
   \li 4 Reactor,
   \li 8 Pool,
   \li 16 FabricationPlant.

   */
  CLASSFacility(CLASSLogger* log, cSecond cycletime, int type = 0);
  //}

  //{
  /// Special Constructor.
  /*!
   Make a new Facility
   \param log : used for the log.
   \param creationtime creation date of the Facility [s],
   \param lifetime operating duration [s],
   \param type identification type of the facility :
   \li 4 Reactor,
   \li 8 Pool,
   \li 16 FabricationPlant.

   */
  CLASSFacility(CLASSLogger* log, cSecond creationtime, cSecond lifetime,
                int type = 0);
  //}

  //{
  /// Special Constructor.
  /*!
   Make a new Facility
   \param log : used for the log.
   \param creationtime creation date of the Facility [s],
   \param lifetime operating duration [s],
   \param cycletime duration of the cycle [s],
   \param type identification type of the facility :
   \li 4 Reactor,
   \li 8 Pool,
   \li 16 FabricationPlant.

   */
  CLASSFacility(CLASSLogger* log, cSecond startingtime, cSecond lifetime,
                cSecond cycletime, int type = 0);
  //}

  //********* Get Method *********//
  /*!
   \name Get Function
   */
  //@{

  int GetId() const { return fId; }  //!< Return the Facility Parc'Is
  IsotopicVector GetInsideIV() const {
    return fInsideIV;
  }  //!< Return the IV contained in the Facility

  int GetFacilityType() const {
    return fFacilityType;
  }  //!< Return the Facility Type id

  cSecond GetInternalTime() const {
    return fInternalTime;
  }  //!< Return Creation Time

  cSecond GetCycleTime() const {
    return fCycleTime;
  }  //!< Return the cycle time of the Facility
  cSecond GetCreationTime() const {
    return fCreationTime;
  }  //!< Return the creation time of the Facility
  cSecond GetLifeTime() const {
    return fLifeTime;
  }  //!< Return the life time of the Facility
#ifndef __ROOTCLING__
  Scenario* GetParc() { return fParc; }  //!< return the pointer to the Park
#endif

  IsotopicVector GetCumulativeIVIn() {
    return fCumulativeIVIn;
  }  //!< return the cumulative sum of all incoming IV
  IsotopicVector GetCumulativeIVOut() {
    return fCumulativeIVOut;
  }  //!< return the cumulative sum of all outcoming IV
  //@}

  //********* Set Method *********//
  /*!
   \name Set Function
   */
  //@{
  void SetId(int id) { fId = id; }  //!< Set The Facility Parc'Id
#ifndef __ROOTCLING__
  void SetParc(Scenario* parc) {
    fParc = parc;
  }  //!< Set the Pointer to the Parc
#endif
  void SetFacilityType(int type) {
    fFacilityType = type;
  }  //!< Set the facility type :
     /// \li 2 reactor Studown
     /// \li 4 start/End of reactor cycle,
     /// \li 8 end of Cooling,
     /// \li 16 fuel Fabrication
  using CLASSObject::SetName;
  using CLASSObject::GetName;

  void SetInsideIV(IsotopicVector const& isotopicvector) {
    fInsideIV = isotopicvector;
  }  //!< Set the IV inside the Facility

  void SetCreationTime(cSecond CTtime) {
    fCreationTime = CTtime;
  }  //!< Set the creation Time
  void SetLifeTime(cSecond Ltime) {
    fLifeTime = Ltime;
  }  //!< Set the life time of the facility
  void SetShutDownTime(cSecond SDTime) {
    fLifeTime = SDTime - fCreationTime;
  }  //!< Set the shutdown time of the facility

  void SetInCycleTime(cSecond ICtime) {
    fInCycleTime = ICtime;
    fIsStarted = true;
  }  //!< Set the In cycle time
  void SetInternalTime(cSecond INtime) {
    fInternalTime = INtime;
  }  //!< Set the Internal Time
  virtual void SetCycleTime(cSecond Ctime) {
    fCycleTime = Ctime;
  }  //!< Set the cycle time (Cycle of the loading Plan)

  //@}

  /*!
   \name Evolution Method
   */
  //@{
  virtual void ApplyZAIThreshold(
      int z = 90);  //!< Put all nuclei below the threshold in -2 -2 -2 ZAI...
  void AddCumulativeIVIn(IsotopicVector const& IV) {
    fCumulativeIVIn += IV;
  }  //!< Add the Input IsotopicVector in the cumulative IV IN
  void AddCumulativeIVOut(IsotopicVector const& IV) {
    fCumulativeIVOut += IV;
  }  //!< Add the Input IsotopicVector in the cumulative IV OUT
  virtual void Evolution(cSecond t) = 0;  //!< Performs the Evolution to time t
  virtual void Dump() {}  //!< Write Modification (IV In/Out, filling the TF...)

  //@}
 protected:
  bool fIsStarted;       ///< True if Running, False Otherwise
  bool fIsShutDown;      ///< True if the facility is stoped, False Otherwise
  bool fIsAtEndOfCycle;  ///< True if Reaching the end of a Facility cycle

  cSecond fInternalTime;  ///< Internal Clock [s]
  cSecond
      fInCycleTime;    ///< Time spent since the beginning of the last Cycle [s]
  cSecond fCycleTime;  ///< Cycle duration Time [s]

  IsotopicVector
      fInsideIV;  ///< All IV in the Facility (fuel for reactor, total for all others...)
  IsotopicVector
      fCumulativeIVIn;  ///< All IV in the Facility (fuel for reactor, total for all others...)
  IsotopicVector
      fCumulativeIVOut;  ///< All IV in the Facility (fuel for reactor, total for all others...)

  //********* Internal Parameter *********//
 private:
  int fId;            //!< Identity of the Facility inside the Scenario
  int fFacilityType;  ///< Type of facility :
                      /// \li 4 reactor,
                      /// \li 8 Pool,
                      /// \li 16 FabricationPlant.

#ifndef __ROOTCLING__
  Scenario* fParc;  //!< Pointer to the main Scenario
#endif
  cSecond fCreationTime;  ///< CLASS Universal Time of Creation [s]
  cSecond
      fLifeTime;  ///< Time of life Of the Reactor (operating's duration) [s]

  ClassDef(CLASSFacility, 1);
};

#endif
