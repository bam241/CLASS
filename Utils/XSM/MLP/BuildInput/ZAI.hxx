#ifndef MN_ZAI_H_
#define MN_ZAI_H_

#include <map>
using namespace std;

class ZAI {
	public :
		int Z;
		int A;
		int I;

	
		ZAI(int z, int a, int i=0) { Z=z;A=a; I=i;}
	
		bool operator <(const ZAI& zai) const	{return (Z != zai.Z)?  
									(Z < zai.Z) : ( (A != zai.A)?
												 (A < zai.A) : (I < zai.I) );}
};

class IsotopicVector {
	public :
	
	IsotopicVector() {}	///< Normal Constructor.
 	~IsotopicVector() {};	 ///< Normal Destructor.
	double	GetZAIIsotopicQuantity(const ZAI& zai) const
	{
		map<ZAI ,double> IsotopicQuantity = IVquantity;
		
		map<ZAI ,double>::iterator it;
		it = IsotopicQuantity.find(zai);
		
		if ( it != IsotopicQuantity.end() )
			return it->second;
		else
			return 0;
	}
	map<ZAI,double> IVquantity;
	double MASS;
};

#endif //MN_ZAI_H_
