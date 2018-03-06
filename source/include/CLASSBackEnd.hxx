#ifndef _CLASSBACKEND_
#define _CLASSBACKEND_

#include <fstream>
#include <map>
#include <string>

#include "CLASSFacility.hxx"
#include "DecayDataBank.hxx"
#include "IsotopicVector.hxx"

#include "TNamed.h"

//________________________________________________________________________
//! Class defining the common properties of all back end fuel cycle facilities
/*!
 Define a CLASS Facility.
 The aim of these class is to gather all the commom properties of the
 facilities which are involve in the BackEnd Fuel cycle.

 @author BaM
*/
//________________________________________________________________________

class CLASSBackEnd : public CLASSFacility {
 public:
  //{
  /// Default Constructor.
  /*!Create an empty CLASSBackEnd
   \param type
   \li -2 :SeparationPlant
   \li -1 :Storage
   \li 8 :Pool
   */
  CLASSBackEnd(int type = 0);
  //}

  //{
  /// CLASSLogger Constructor.
  /*!
   Create an empty CLASSBackEnd loading a CLASSLogger
   \param log : used for the log.
   \param type
   \li -2 :SeparationPlant
   \li -1 :Storage
   \li 8 :Pool     */
  CLASSBackEnd(CLASSLogger* log, int type = 0);
  //}

  //{
  /// Cycle time Constructor.
  /*!
   Create an empty CLASSBackEnd loading a CLASSLogger
   \param log : used for the log.
   \param cycletime Cycle time of the CLASSBackend (e.g. Cooling time for the
   pool) in [s],
   \param type -2 :SeparationPlant -1 : Storage ; 8 :Pool
   */
  CLASSBackEnd(CLASSLogger* log, cSecond cycletime, int type = 0);
  //}
  //@}

  //********* Get Method *********//
  /*!
   \name Get Function
   */
  //@{

  std::vector<IsotopicVector> GetIVArray()
      const;  //!< Return the IsotopicVector Array
  std::vector<cSecond> GetIVArrayArrivalTime()
      const;  //!<Vector of arrival time of each IV in the CLASSBackEnd

  int GetIVNumber() const;  //!< Return the number of Isotopic Vector present in
                            //!the CLASSBackEnd object
  bool GetStorageType()
      const;  //!< Return the storageType : True if it is a Storage
  IsotopicVector GetIV(int i) const;

#ifndef __ROOTCLING__
  DecayDataBank*
  GetDecayDataBank();  //!< Return the pointer to the decay DataBank
  CLASSBackEnd* GetOutBackEndFacility()
      const;  //!< Return the pointer to the OUtBackEndFacility
  virtual std::map<cSecond, int>
  GetTheBackEndTimePath();  //!< Get the full path
#endif

  //@}

  //********* Set Method *********//
  /*!
   \name Set Function
   */
  //@{
  void SetIsStorageType(bool val = true);  //!< Set the fIsStorage bool
  virtual void SetIVArray(std::vector<IsotopicVector> const&
                              ivarray);  //!< Set The isotopicVector Array
  void SetIVArrayArrivalTime(
      std::vector<cSecond> const&
          IVArrayArrivalTime);  //!< Set Arrival Time in Back end

#ifndef __ROOTCLING__
  void SetDecayDataBank(DecayDataBank* decayDB);  //!< Set the Decay DataBank
  virtual void SetOutBackEndFacility(
      CLASSBackEnd* befacility);  //!< Set an out Facility
#endif

  using CLASSFacility::SetName;

  //@}

  /*!
   \name BackEndFacility specific Method
   */
  //@{
  virtual void ApplyZAIThreshold(
      int z = 90);  //!< Put all nuclei below the threshold in -2 -2 -2 ZAI...
  virtual void AddIV(
      IsotopicVector isotopicvector);  //!< Add an Isotopicvector to the IVArray

  void ClearIVArray();  //!< Empty the IVArray removing all fuel stored

  //@}
  virtual void Evolution(cSecond t) {}  //!< Performs the Evolution until time t

  void UpdateInsideIV();

 protected:
  IsotopicVector GetDecay(
      IsotopicVector const& isotopicvector,
      cSecond t);  //!< Get IsotopicVector decay at time t [s]

  //********* Internal Parameter *********//
  std::vector<IsotopicVector>
      fIVArray;  ///< Vector containning all the fuel stored.
  std::vector<cSecond>
      fIVArrayArrivalTime;  ///< Vector containning the arrival time of each fuel in [s]

#ifndef __ROOTCLING__
  CLASSBackEnd* fOutBackEndFacility;  //!< Facility getting the fuel at the end
                                      //!of the cycle
#endif

 private:
  bool fIsStorageType;  //!< True if there is no out CLASSBackEnd facility (like
                        //!a Storage...)

#ifndef __ROOTCLING__
  DecayDataBank* fDecayDataBase;  //!< Pointer to the DecayDataBank
#endif

  ClassDef(CLASSBackEnd, 2);
};

///// inline functions ///////////////////////////////////////////////////////
inline std::vector<IsotopicVector> CLASSBackEnd::GetIVArray() const {
  return fIVArray;
}
inline std::vector<cSecond> CLASSBackEnd::GetIVArrayArrivalTime() const {
  return fIVArrayArrivalTime;
}
inline int CLASSBackEnd::GetIVNumber() const { return fIVArray.size(); }
inline bool CLASSBackEnd::GetStorageType() const { return fIsStorageType; }
inline IsotopicVector CLASSBackEnd::GetIV(int i) const {
  if (i < (int)fIVArray.size()) {
    return fIVArray[i];
  } else {
    return IsotopicVector();
  }
}

#ifndef __ROOTCLING__
inline DecayDataBank* CLASSBackEnd::GetDecayDataBank() {
  return fDecayDataBase;
}
inline CLASSBackEnd* CLASSBackEnd::GetOutBackEndFacility() const {
  return fOutBackEndFacility;
}
#endif

inline void CLASSBackEnd::SetIsStorageType(bool val) { fIsStorageType = val; }
inline void CLASSBackEnd::SetIVArrayArrivalTime(
    std::vector<cSecond> const& IVArrayArrivalTime) {
  fIVArrayArrivalTime = IVArrayArrivalTime;
}
inline void CLASSBackEnd::SetIVArray(
    std::vector<IsotopicVector> const& ivarray) {
  fIVArray = ivarray;
}

#ifndef __ROOTCLING__
inline void CLASSBackEnd::SetDecayDataBank(DecayDataBank* decayDB) {
  fDecayDataBase = decayDB;
}
inline void CLASSBackEnd::SetOutBackEndFacility(CLASSBackEnd* befacility) {
  fOutBackEndFacility = befacility;
  fIsStorageType = false;
}
#endif

#endif
