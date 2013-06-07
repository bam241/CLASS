using namespace std;


#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
string dtoa(double num)
{
	ostringstream os(ostringstream::out);
	os<<setprecision(3)<<num;
	return os.str();
}

int main(int argc, char** argv)
{

	CLASS gCLASS;
	gCLASS.SetStockManagement(true);
	cSecond year = 3600*8766; // 3600*24*365.25;//Year duration in seconds
	gCLASS.SetTimeStep(year/4);
	
	
	cout << "DB Definition... \t" << endl;
	cout << "\t -decay \t" ;
	DataBank<ZAI>* DecayDB = new DataBank<ZAI>(gCLASS.GetLog(), "/home/cramal/CLASS/DataBase/Decay.idx");
	gCLASS.SetDecayDataBase(DecayDB);
	cout << "\t...OK!" << endl;
	
	cout << "\t -MOX \t" ;
	DataBank<IsotopicVector>* MOX = new DataBank<IsotopicVector>(gCLASS.GetLog(), "/home/cramal/CLASS/DataBase/DB_40_45/DB.idx");
	cout << "Calculation of Distance Parameters : " << endl;
	MOX->SetDistanceType(1);//Use coefficient calculated by CLASS to estimate the distance between an Isotopic Vector and an EvolutionData from the MOX DataBank
	cout << "\t...OK!" << endl;

	cout << "\t -UOx \t" ;
	EvolutionData DB_REP_UOX = EvolutionData(gCLASS.GetLog(), "/home/cramal/CLASS/DataBase/REPUOX.dat");
	cout << "\t...OK!" << endl;
	cout << "DB Definition...\t \t \t...Done!" << endl;

	cout << "Storage Definition... \t"<<endl  ;
	Storage *Stock = new Storage(gCLASS.GetLog());
	Stock->SetName("Stock");
	Storage *ReUsable = new Storage(gCLASS.GetLog());
	ReUsable->SetName("ReUsable");
	cout << "\t...OK!" << endl;
	
	cout << "Pool Definition... \t" <<endl ;
	Pool *Cooling_UOX = new Pool(gCLASS.GetLog(),Stock, gCLASS.GetAbsoluteTime(), (double)year*5 );
	Cooling_UOX->SetName("Cooling Pool for UOX");
	Pool *Cooling_MOX = new Pool(gCLASS.GetLog(),Stock, gCLASS.GetAbsoluteTime(), (double)year*5 );
	Cooling_MOX->SetName("Cooling Pool for MOX");
	cout << "\t...OK!" << endl;

	cout << "FabricationPlant Definition... \t" <<endl ;
	FabricationPlant *FP_MOX = new FabricationPlant(gCLASS.GetLog(),Stock, ReUsable);
	FP_MOX->SetName("Fabrication Plant for MOX");
	cout << "\t...OK!" << endl;


	double UOX_Cycle = year*4;
	double MOX_Cycle = year*4;
	
	cout << "Reactor Definition... \t"<<endl ;
	int nUOX=4;
	Reactor *REP_UOX[nUOX];
	int nMOX=1;
	Reactor *REP_MOX[nMOX];

	for(int i=0;i<nUOX;i++){

	REP_UOX[i] = new Reactor(gCLASS.GetLog(),DB_REP_UOX, Cooling_UOX,0, UOX_Cycle*(10+3*i), 3e9,100,40,0.8 ) ;
	string Name;
	Name="UOX "+dtoa(i);
	REP_UOX[i]->SetName(Name.c_str());
	}

	for(int i=0;i<nMOX;i++){
	string Name;
	REP_MOX[i] = new Reactor(gCLASS.GetLog(), MOX, FP_MOX, Cooling_MOX, UOX_Cycle*(10+3*i) , MOX_Cycle*(10+3*(nMOX-i)), 3e9, 100, 40,0.8);
	Name="MOX "+dtoa(i);
	REP_MOX[i]->SetName(Name.c_str());
	}

	cout << "\t...OK!" << endl;
	cout << endl << endl;

	cout << "Addind Phase : " << endl;
	cout << "\t -FabricationPlant" ;
	gCLASS.AddFabricationPlant(FP_MOX);
	cout << "\t...OK!" << endl;

	cout << "\t -Pool \t\t" ;
	gCLASS.AddPool(Cooling_UOX);
	gCLASS.AddPool(Cooling_MOX);
	cout << "\t...OK!" << endl;

	cout << "\t -Reactor \t\t" ;
		
		
	for(int i = 0; i< nUOX; i++ ) gCLASS.AddReactor(REP_UOX[i]);

	for(int i = 0; i< nMOX; i++ ) gCLASS.AddReactor(REP_MOX[i]);
	
	cout << "...OK!" << endl;

	cout << "\t -Storage \t\t" ;
	gCLASS.AddStorage(Stock);

	gCLASS.AddStorage(ReUsable);

	cout << "...OK!" << endl;

	cout << "Addind Phase...\t \t \t...Done!" << endl;
	cout << endl << endl;
	


	cout << "Beginning the Evolution" << endl;
	gCLASS.Evolution((double)year*100);

}
//==========================================================================================
// Compilation
//==========================================================================================
/*

\rm CLASS* ; g++ -o CLASS_exec Example_CLASS.cpp -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -Wunused-result 


*/
