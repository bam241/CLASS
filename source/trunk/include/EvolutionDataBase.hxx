#ifndef __EVOLUTIONDATABASE_HXX__
#define __EVOLUTIONDATABASE_HXX__

/*!
 \file
 \brief Header file for EvolutionDataBase class. 
 */


#include <map>
#include <string>


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


class EvolutionDataBase 
{
	
public :
//********* Constructor/Destructor Method *********//

	EvolutionDataBase(LogFile* Log, string DB_index_file = "Default_Index.dat" );
	~EvolutionDataBase();

//********* Get Method *********//

	map<ZAI ,EvolutiveProduct* >	GetEvolutionDataBase() const {return fEvolutionDataBase;}
	string 	GetDataBaseIndex() const {return fDataBaseIndex;}
	bool 	IsDefine(const ZAI& zai) const;
	
//********* Set Method *********//

	void SetDataBaseIndex(string database) {fDataBaseIndex = database;}
	
//********* Modification Method *********//
	IsotopicVector	DecayProduction(const ZAI &zai, double dt); 
					///< Return the Product IsotopicVector evolution from zai during a dt time
	bool AddEvolutiveProduct(const ZAI& zai);



//********* Printing Method *********//
	void Print() const;
	
protected :
	
	map<ZAI ,EvolutiveProduct* >	fEvolutionDataBase;
 	string				fDataBaseIndex;
 	LogFile*			fLog;
 	

};


#endif
