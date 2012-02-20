// These classes are not fuly implemented...
// this is not a big job to do it...so you can

#ifndef _TMKcode_HXX_
#define _TMKcode_HXX_

#include <iostream>
#include <vector>
using namespace std;
#include <libValErr/ValErr.hxx>

//________________________________________________________________________
//@{
// 		A cycle of a KCODE.
//
//	@author WEC
//	@version 1.0
//@}
class TMKcodeCycle
{

 public :
	virtual unsigned Read(istream& s, bool verbose=1);	//@- read the cycle in "m" file 
	virtual unsigned Write(ostream& s);					//@- write a TMKcodeCycle 
	virtual unsigned Add(TMKcodeCycle *Other);			//@- sum 2 cycles
	virtual void operator *= (const float& factor);		//@- multply each attribute by factor
	virtual ~TMKcodeCycle(){};     //@- Virtual destructor to avoid GCC 4.0 warnings

	float fKEffCollision;    //@- collision keff estimator
	float fKEffAbsorption;   //@- absorption keff estimator
	float fKEffTL;           //@- track length keff estimator
	float fLifeCollision;    //@- collision prompt life time 
	float fLifeAbsorption;   //@- absorption prompt life time

};

//________________________________________________________________________
//@{
// 		A cycle of a KCODE.
//
//	@author WEC
//	@version 1.0
//@}

class TMKcodeCycle19 : public TMKcodeCycle 
{

 public :
	unsigned Read(istream& s, bool verbose=1);		//@- read the cycle in "m" file 
	unsigned Write(ostream& s);						//@- write a TMKcodeCycle19 
	virtual unsigned Add(TMKcodeCycle *Other);		//@- sum 2 cycles
	virtual void operator *= (const float& factor);	//@- multply each attribute by factor
	virtual	~TMKcodeCycle19(){};    //@- Virtual destructor to avoid GCC 4.0 warnings

	ValErr_t fKEffCollisionMoy;		//@- average over preceding cycles of fKEffCollision
	ValErr_t fKEffAbsorptionMoy;	//@- average over preceding cycles of fKEffAbsorption
	ValErr_t fKEffTLMoy;			//@- average over preceding cycles of fKEffTL
	ValErr_t fKEffAll;				//@- average over preceding cycles of mixed Keff estimator
	ValErr_t fKEffAllSkip;			//@- average over skip cycles of mixed Keff estimator
	ValErr_t fLifeAll;				//@- mixed prompt life time estimator
	float fHistories;				//@- Neutron number for the cycle
	float fKEffAllFM;				//@- FOM of fKEffAllAverage

};
//________________________________________________________________________
//@{
// 		All KCODE cycles.
//
//	@author WEC
//	@version 1.0
//@}
class TMKcode 
{
 public :
	//@{ 
	//	Normal Constructor.
	// @param ActiveCycles: active cycle in a Kcode
	// @param InactiveCycles: inactive cycle in a Kcode
	// @param Taille: size 5 or 19 (records per cycle)
	//@}
	TMKcode(unsigned long ActiveCycles,unsigned long InactiveCycles,unsigned long Taille);
	virtual ~TMKcode();		//@- Destructor
	unsigned Read(istream& s, bool verbose=1);	//@- Read KCODE cycles from "m" file
	unsigned Write(ostream& s);					//@- write KCODE cycles 
  
	TMKcodeCycle* GetCycle(unsigned long n){return fKcodeCycleTab[n];}						//@- return cycle n (size=5)
	TMKcodeCycle19* GetCycle19(unsigned long n){return (TMKcodeCycle19*)fKcodeCycleTab[n];}	//@- return cycle n (size=19)
  
	//@{
	//	Add 2 Kcodes.
	//
	// The sum is made cycle per cycle. It returns 1 even if number of cycles is different (keeps only the minimum number).
	//@}
	unsigned Add(TMKcode *Other); 
	unsigned long GetTotalCycle() {return fTotalCycles;}//@- returns the number of effective cycles
	unsigned long GetTaille() {return fTaille;}			//@- returns the size of cycles (5 or 19)
	TMKcode& operator *= (const float& factor);			//@- multply each attribute of each cycle by factor
	
 protected :
	void Init(){/* not implemented */}
	void Free(){/* not implemented */}

	unsigned long fTotalCycles;				//@- total number of effective cycles
	unsigned long fInactiveCycles;			//@- number of inactive cycles
	unsigned long fTaille;					//@- size(records) of each cycle (5 or 19)
	vector<TMKcodeCycle*> fKcodeCycleTab;	//@- recorded cycle vector
};											


#endif
