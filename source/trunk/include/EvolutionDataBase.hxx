#ifndef __EVOLUTIONDATABASE_HXX__
#define __EVOLUTIONDATABASE_HXX__

/*!
 \file
 \brief Header file for EvolutionDataBase class. 
 The aim of this Class is to store the evolution Database of the all decay nuclei.
 
 @author BaM
 @version 0.
 */

#include <map>
#include <string>
#include <vector>

using namespace std;

class IsotopicVector;
class ZAI;
class EvolutiveProduct;
class LogFile;



///< A ZAIIDataBase defined a database which contain the evolution of faction of all product, subproduct (or sub(sub...sub)product) for a nucleus. 
/*!
 The aim of this class is to handle the evolution all information of all Nuclueus product, subproduct (or sub(sub...sub)product) as the fonction of the time.
 
 @author BaM
 @version 1.0
 */

template <class T> 
class EvolutionDataBase {

public :
//********* Constructor/Destructor Method *********//
	///< Normal Constructor.
	EvolutionDataBase(LogFile* Log, string DB_index_file );
	
	///< Normal Destructor.
	~EvolutionDataBase();

//********* Get Method *********//
	LogFile*			GetLog()			{ return fLog; }			//!< Return the Pointer to Log
	map<T ,EvolutiveProduct >	GetEvolutionDataBase()	const	{ return fEvolutionDataBase; }		//!< Return the EvolutionDataBase
	string 				GetDataBaseIndex()	const	{ return fDataBaseIndex; }		//!< Return the index Name
	string				GetFuelType()		const	{ return fFuelType; }			//!< Return the fuel type of the DB
	vector<double>			GetPFuelParameter()	const	{ return fFuelParameter; }		//!< Return the Fuel parameter of the DB
	pair<double,double>		GetBurnUpRange()	const	{ return fBurnUpRange;}			//!< Return the BurnUp range of the DB
	bool 				IsDefine(const T& key)	const;						//!< True the key is define, false unstead

	map<double, EvolutiveProduct>	GetDistancesTo(IsotopicVector isotopicvector, double t = 0) const;	//! Return a map containing the distance of each EvolutiveProduct in the DataBase to the set IV at the t time

//********* Set Method *********//

	void SetDataBaseIndex(string database) { fDataBaseIndex = database; }
	EvolutiveProduct GenerateDB(IsotopicVector isotopicvector, double cycletime, double Power); //!< Genration of a New EvolutiveProduct From the one already present


//********* Modification Method *********//
	IsotopicVector	Evolution(const T &key, double dt);	///< Return the Product IsotopicVector evolution from zai during a dt time
	void		ReadDataBase();				///< ...


//********* Printing Method *********//
	void Print() const;
	
protected :
	
	map<T, EvolutiveProduct>	fEvolutionDataBase;
 	string				fDataBaseIndex;
 	LogFile*			fLog;


 	string 				fFuelType;
 	pair<double,double>		fBurnUpRange;
 	vector<double>			fFuelParameter;


};



#endif
