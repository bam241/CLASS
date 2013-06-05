#ifndef __DataBank_HXX__
#define __DataBank_HXX__

/*!
 \file
 \brief Header file for DataBank class. 
 The aim of this Class is to store the evolution Database of the all decay nuclei.
 
 @author BaM, Marc
 @version 0.
 */

#include <map>
#include <string>
#include <vector>

using namespace std;
typedef long long int cSecond;

class IsotopicVector;
class ZAI;
class EvolutionData;
class LogFile;



///< A ZAIIDataBase defined a database which contain the evolution of faction of all product, subproduct (or sub(sub...sub)product) for a nucleus. 
/*!
 The aim of this class is to handle the evolution all Information of all Nuclueus product, subproduct (or sub(sub...sub)product) as the fonction of the time.
 
 @author BaM
 @version 1.0
 */

template <class T> 
class DataBank {

public :
//********* Constructor/Destructor Method *********//
	///< Normal Constructor.
	DataBank();
	
	DataBank(LogFile* Log, string DB_index_file );
	
	///< Normal Destructor.
	~DataBank();

//********* Get Method *********//
	LogFile*			GetLog()			{ return fLog; }			//!< Return the Pointer to Log
	map<T ,EvolutionData >	GetDataBank()	const	{ return fDataBank; }		//!< Return the DataBank
	string 				GetDataBaseIndex()	const	{ return fDataBaseIndex; }		//!< Return the index Name
	string				GetFuelType()		const	{ return fFuelType; }			//!< Return the fuel type of the DB
	vector<double>			GetFuelParameter()	const	{ return fFuelParameter; }		//!< Return the Fuel parameter of the DB
	pair<double,double>		GetBurnUpRange()	const	{ return fBurnUpRange;}			//!< Return the BurnUp range of the DB
	bool 				IsDefine(const T& key)	const;						//!< True the key is define, false unstead

	map<double, EvolutionData>	GetDistancesTo(IsotopicVector isotopicvector, double t = 0) const;	//! Return a map containing the distance of each EvolutionData in the DataBase to the set IV at the t time
	EvolutionData	GetClosest(IsotopicVector isotopicvector, double t = 0) const;	//! Return the closest

//********* Set Method *********//
	void SetDataBank(map<T ,EvolutionData > mymap)	{ fDataBank = mymap; } 

	void SetDataBaseIndex(string database) { fDataBaseIndex = database; }
	EvolutionData GenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power); //!< Genration of a New EvolutionData From the one already present
	void SetUpdateReferenceDBatEachStep(bool val) {fUpdateReferenceDBatEachStep = val;}

//********* Modification Method *********//
	IsotopicVector	Evolution(const T &key, double dt);	///< Return the Product IsotopicVector evolution from zai during a dt time
	void		ReadDataBase();				///< ...
	void		CalculateDistanceParameter();///< Calculate automaticly the weight for each ZAI in the distance calculation from the mean XS of the DataBank
	void		SetDistanceParameter(IsotopicVector DistanceParameter);///< Define mannually the weight for each ZAI in the distance calculation 
	void		SetDistanceType(int DistanceType);///< Define the way to decide if two isotopic vectors are close. 0 is for the standard norme, 1 for each ZAI weighted with its XS, 2 for each ZAI weighted with coefficient given by the user


//********* Printing Method *********//
	void Print() const;
	
protected :
	
	map<T, EvolutionData>	fDataBank;
 	string				fDataBaseIndex;
 	LogFile*			fLog;

	bool		fUpdateReferenceDBatEachStep;


 	string 				fFuelType;
 	pair<double,double>		fBurnUpRange;
 	vector<double>			fFuelParameter;

	int 	fDistanceType;		///< 0 is for the standard norme,
					///< 1 for each ZAI weighted with its XS,
					///< 2 for each ZAI weighted with coefficient given by the user
	
	T	fDistanceParameter;	///< weight for each ZAI in the distance calculation


};



#endif
