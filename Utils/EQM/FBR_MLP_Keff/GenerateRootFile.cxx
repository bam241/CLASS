/**********************************************************/
//	Make the input file for the MLPs training
//
//	This programs reads a set of .dat files which are the
//	results of a depletion calculation (see manual and 
//  looks for XS_CLOSEST). From the reading it fills a 
//  TTree (Data) and write it in a file named
//	TrainingInput.root . 
//	The file TrainingInput.cxx is the list of MLP outputs 
// (cross sections)
//
//@author BaM, BaL
/**********************************************************/
#include "include/GenerateRootFile.hxx"
#include <TH1F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include "include/StringLine.hxx"
#include <TString.h>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <sstream>

using namespace std;
//--------------------------------------------------------------------------------------------------
/*************************
 		MAIN
*************************/
int main(int argc, char ** argv){
	
	if(argc!=2)
	{
		cout << "Usage : BuildInputTree Path" << endl;
		cout << " Where Path is the path to the folder containing the Evolution Datas" << endl;
		cout << " i.e the (.dat) files" << endl;
		exit(0);
	}	

	fEvolutionDataFolder = argv[1];

	InitMass();	//Load nuclei masses
	CheckJob();	// looks fot the .dat files in the fEvolutionDataFolder

	FilePath = "DB_TMP/";
	DataPath = FilePath + "Data/";
	string Command = "mkdir -p " + DataPath;
	system(Command.c_str());

	cout << "Reading .dat files ..." << endl;
	for(int i = 0; i < (int)JobName.size(); i++)
	{
		ReadAndFill(JobName[i]);
		if (i%100 == 0)
			cout << "\r" << i << " .dat files read" <<flush;
	}
	cout << "Filling the TTree ..." << endl;
	DumpInputNeuron("TrainingInput.root");
	cout << "Training input generated in file : TrainingInput.root " << endl;
	cout << "Names of MLP outputs in file : TrainingInput.cxx " << endl;
	system("rm -r DB_TMP");

}
//--------------------------------------------------------------------------------------------------
// Function definitions
//--------------------------------------------------------------------------------------------------
//Convert int to string
string itoa(int num)
{
	ostringstream os(ostringstream::out);
	os << num;
	return os.str();
}
//--------------------------------------------------------------------------------------------------
void DumpInputNeuron(string filename)
{
	TFile*   fOutFile = new TFile(filename.c_str(),"RECREATE");
	TTree*   fOutT = new TTree("Data", "Data");

/**********************INITIALISATIONNN********************/
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
	//@@@Change the input value according to your fresh fuel compo 
	//-> here MOX FUEL
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
	double U5   			 = 0;
	double U8   			 = 0;
	double Pu8 				 = 0;
	double Pu9 				 = 0;
	double Pu10 			 = 0;
	double Pu11				 = 0;
	double Pu12 	    	 = 0;
	double Am1 	    		 = 0;

	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
	double k_eff 			 = 0;
	

/**********************BRANCHING**************************************************/
/**********************Fresh fuel**************************************************/
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
	//@@@Change the input value according to your fresh fuel compo 
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
	fOutT->Branch(	"U5"			,&U5			,"U5/D"				);
	fOutT->Branch(	"U8"			,&U8			,"U8/D"				);
	fOutT->Branch(	"Pu8"			,&Pu8			,"Pu8/D"			);
	fOutT->Branch(	"Pu9"			,&Pu9			,"Pu9/D"			);
	fOutT->Branch(	"Pu10"			,&Pu10			,"Pu10/D"			);
	fOutT->Branch(	"Pu11"			,&Pu11			,"Pu11/D"			);
	fOutT->Branch(	"Pu12"			,&Pu12			,"Pu12/D"			);
	fOutT->Branch(	"Am1"			,&Am1			,"Am1/D"			);
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
	fOutT->Branch(	"k_eff"			,&k_eff			,"k_eff/D"			);

/**********************FILLING THE TTREE**************************************************/

	//Fill containing all the output of the networks to train

	 int NumOfBase=fActinideCompoInit.size();
	for(int b=0;b<NumOfBase;b++) 
	{ 
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
	//@@@Change the input value according to your fresh fuel compo 
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
															//   (Z , A ,I)
		U5 			=  fActinideCompoInit[b].GetZAIIsotopicQuantity(ZAI(92,235,0));
		U8 			=  fActinideCompoInit[b].GetZAIIsotopicQuantity(ZAI(92,238,0));
		Pu8  	  	=  fActinideCompoInit[b].GetZAIIsotopicQuantity(ZAI(94,238,0));
		Pu9  	  	=  fActinideCompoInit[b].GetZAIIsotopicQuantity(ZAI(94,239,0));
		Pu10 	  	=  fActinideCompoInit[b].GetZAIIsotopicQuantity(ZAI(94,240,0));
		Pu11 	  	=  fActinideCompoInit[b].GetZAIIsotopicQuantity(ZAI(94,241,0));
		Pu12 	  	=  fActinideCompoInit[b].GetZAIIsotopicQuantity(ZAI(94,242,0));
		Am1 	  	=  fActinideCompoInit[b].GetZAIIsotopicQuantity(ZAI(95,241,0));
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
	//@@@Change the TStep according the step you want
	//@@@ I have assumed all your depletion calculations have 
	//@@@ the same time binning
	//@@@ hoping is what you have done ....		
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
		int Tstep = 0; //at begining of cycle if 0

	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//				
	
		k_eff=fkeff[b][Tstep];

		fOutT->Fill();		
	}

	fOutFile->Write();
	delete fOutT;
	fOutFile-> Close();
	delete fOutFile;

}
//--------------------------------------------------------------------------------------------------
void InitMass()
{
//Set tghe mass of the nuceli

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
//@@@Change : ADD THE MASS OF THE NUCLEI PRESENT IN YOUR FRESH FUEL
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
	
	ZAI U238  = ZAI(92,238,0);
	ZAImass.insert(pair<ZAI,double>(U238, 238050788.247e-6));
	
	ZAI U234  = ZAI(92,234,0);
	ZAImass.insert(pair<ZAI,double>(U234, 234041000.000e-6));

	ZAI U235  = ZAI(92,235,0);
	ZAImass.insert(pair<ZAI,double>(U235, 235043929.918e-6));
	
	ZAI Pu238 = ZAI(94,238,0);
	ZAImass.insert(pair<ZAI,double>(Pu238,238049559.894e-6));
	
	ZAI Pu239 = ZAI(94,239,0);
	ZAImass.insert(pair<ZAI,double>(Pu239,239052163.381e-6));
	
	ZAI Pu240 = ZAI(94,240,0);
	ZAImass.insert(pair<ZAI,double>(Pu240,240053813.545e-6));
	
	ZAI Pu241 = ZAI(94,241,0);
	ZAImass.insert(pair<ZAI,double>(Pu241,241056851.456e-6));
	
	ZAI Pu242 = ZAI(94,242,0);
	ZAImass.insert(pair<ZAI,double>(Pu242,242058742.611e-6));
	
	ZAI Am241 = ZAI(95,241,0);
	ZAImass.insert(pair<ZAI,double>(Am241,241056829.144e-6));
	
	ZAI O16 = ZAI(8,16,0);
	ZAImass.insert(pair<ZAI,double>(O16, 0.));
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
	
}
//--------------------------------------------------------------------------------------------------
void CheckJob()
{	//LOAD THE LIST OF EvolutionData
	cout << "Scanning " << fEvolutionDataFolder << " for .dat files ..." << endl;
	cout << "Please wait ..."<< endl;

	string Command = "find "+ fEvolutionDataFolder + " -name \"*.dat\" > JOB.tmp";
	system(Command.c_str());
	
	ifstream JOB("JOB.tmp");
	if (JOB.is_open())
	{
		while (!JOB.eof())
		{
			string tmp;
			getline(JOB,tmp);
			JobName.push_back(tmp);
		}
	} JOB.close(); JobName.pop_back();
	
	// Remove temporary files...
	Command = "\\rm -f JOB.tmp";
	system(Command.c_str());
	random_shuffle(JobName.begin(), JobName.end());
	cout << "Scan complete" <<endl; 
}

//--------------------------------------------------------------------------------------------------
void ReadAndFill(string jobname)
{	//Read a .dat file and fill XS maps and the fuel initial composition
//Read a .dat file and fill XS maps and the fuel initial composition

	vector<double>	vT;
	vector<int>	Z;
	vector<int>	A;
	vector<int>	I;
	vector<double>	Q;
	
	ifstream DecayDB(jobname.c_str());							// Open the File
	if(!DecayDB)
	{
		cout << "!!Warning!! !!!EvolutiveProduct!!! \n Can't open \"" << jobname << "\"\n" << endl;
	}
	
	string line;
	int start = 0;
	
	getline(DecayDB, line);
	
	/******Getting Time vecotr ....******/
	if( StringLine::NextWord(line, start, ' ') != "time")
	{
		cout << "!!Bad Trouble!! !!!EvolutiveProduct!!! Bad Database file : " <<  jobname << endl;
		exit (1);
	}
	
	while(start < (int)line.size())
		vT.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));

	fNOfTimeStep=int(vT.size());


	start = 0;
	getline(DecayDB, line);
	string tmp = StringLine::NextWord(line, start, ' ');
	vector<double> vKeff;
	if ( tmp == "keff"  || tmp == "Keff" )
	{
		while(start < (int)line.size())
		vKeff.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
	}
	fkeff.push_back(vKeff);


	/****Getting Inventories***/	
	getline(DecayDB, line);
	do
	{	
		start = 0;
		int z;
		string tmp2 = StringLine::NextWord(line, start, ' ');
		if (tmp2 == "Inv") {
			z = atoi(StringLine::NextWord(line, start, ' ').c_str());
		}
		else z = atoi(tmp2.c_str());
		int a = atoi(StringLine::NextWord(line, start, ' ').c_str());
		int i = atoi(StringLine::NextWord(line, start, ' ').c_str());
		
		if(a!=0 && z!=0)
		{
			
			ZAI zaitmp(z, a, i);
			Z.push_back(z);
			A.push_back(a);
			I.push_back(i);
			if(!fIsAllNucleiAlreadyFill)
			{	
				fAllNuclei.push_back(zaitmp);
				fTime=vT;	
			}	

			long double q = atof(StringLine::NextWord(line, start, ' ').c_str());
			Q.push_back(q);	

		}

		getline(DecayDB, line);
		start = 0;
		tmp2 = StringLine::NextWord(line, start, ' ');
		
		if(line == "" || line == "CrossSection" || tmp2 == "XSFis" || tmp2 == "XSCap" || tmp2 == "XSn2n") break;
	}while (!DecayDB.eof() );

	if(fAllNuclei.size()!=0)
	 fIsAllNucleiAlreadyFill=true;
	
	DecayDB.close();
	
	double N = 0;
	for(int i=0; i < (int)Z.size()-2; i++)
	{
		if(  Z[i]>89 ) 
			N += Q[i];
	}
	

	IsotopicVector CompoBasei; 

	for(int i=0; i < (int)Z.size()-2; i++)
	{
		if(Z[i]>89)
		{
			ZAI zai = ZAI(Z[i], A[i], I[i]);
			CompoBasei.IVquantity.insert(pair<ZAI,double>(zai,Q[i]/N));
		}	
	}
	
	fActinideCompoInit.push_back(CompoBasei);

GoodJobName.push_back(jobname);


}
/*-------------------------------------------------------------------------------------------------
COMPILATION :

g++ -o GenerateRootFile GenerateRootFile.cxx `root-config --cflags` `root-config --libs`

*/