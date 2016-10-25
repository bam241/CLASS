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
#include "Gene.hxx"
#include <TH1F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include "StringLine.hxx"
#include <TString.h>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <sstream>
#include <numeric>
#include <functional>
#include <algorithm>    
#include <ctime>

using namespace std;

const ZAIMass cZAIMass; //atomic masses

string ElNames[110]={"  ","H","He","Li","Be",
						"B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
						"K","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge",
						"As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd",
						"Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
						"Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W",
						"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra",
						"Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
						"Lr","Rf","Db","Sg","Bh","Hs","Mt"};

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
//Give a name to a cross section 
//!!!!!!!!!!!!!!!!!!!!!!!!!!//
// YOU MUST KEEP THE FORMAT :
// XS_Z_A_I_fis	(fission cross section)
// XS_Z_A_I_cap (n,gamma cross section)
// XS_Z_A_I_n2n	(n,2n cross section)
//!!!!!!!!!!!!!!!!!!!!!!!!!!//
string NameXS(ZAI act,string xs)
{
	stringstream Name;
	Name<<"XS_"<<act.Z<<"_"<<act.A<<"_"<<act.I<<"_"<<xs;

	return Name.str();
}

//--------------------------------------------------------------------------------------------------
void  FillMapName()
{
	cout<<endl;
	cout<< "============================================="<<endl;
	cout<< "Looking for fresh fuel composition @ t=0 ..."<<endl;
	cout<< "-> It will be the input variables of the MLPs"<<endl;
	cout<<endl;

	vector<ZAI> ZAI_T0 =  fActinideCompoInit[0].GetNonZeroZAIList();

	for(int zai = 0 ; zai < ZAI_T0.size() ; zai++ )
	{	stringstream ssname;
		int Z = ZAI_T0[zai].Z;
		int A = ZAI_T0[zai].A;
		int I = ZAI_T0[zai].I;
		string sI="";
		if(I == 1)
			sI = "m";
		else if(I == 2) 
			sI = "2m";
		else if(I == 3) 
			sI = "3m";

		ssname<<A<<ElNames[Z]<<sI;

		cout<< "Add this nuclei with this name [y/n] ?  (if you don't know say yes to all)  "<<endl;
		cout<< "[Z\tA\tI\tName]"<<endl;
		cout<<"\033[36m"<<Z<<"\t"<<A<<"\t"<<I<<"\t"<<ssname.str()<<"\033[0m"<<endl;

		if (UserSayYes())
			fMapName.insert(pair<ZAI,string> ( ZAI(Z,A,I) , ssname.str()) );

	}

	bool UserWantToAdd = true;
	while(UserWantToAdd)
	{
		cout<< "Do you want to add additional nuclei  [y/n] ? (if you don't know say no) "<<endl;
		if(UserSayYes())
		{	
			int Z = 0;
			int A = 0;
			int I = 0;
			string Name;
			cout << "Z -> "; cin >> Z ; cout <<" A -> "; cin >> A ; cout <<" I -> "; cin >> I ;cout <<" Name -> "; cin >> Name;
			fMapName.insert(pair<ZAI,string> ( ZAI(Z,A,I) ,Name) );
		}

		else
			UserWantToAdd = false;
	}	

	cout<< "============================================="<<endl;


}
//--------------------------------------------------------------------------------------------------
bool UserSayYes()
{
		bool AnswerIsNotGiven = true;
		bool isYES = false;
		while(AnswerIsNotGiven)
		{
			string answer;
			cin>>answer;
	
			if(answer == "y" || answer == "yes" || answer == "Yes" || answer == "Y")
			{
				isYES = true;
				AnswerIsNotGiven = false;
			}
			else if (answer == "n" || answer == "no" || answer == "No" || answer == "N")
			{
				isYES = false;
				AnswerIsNotGiven = false;
			}

			else{
				cout << "Yes OR No ???!"<<endl;
			}
		}
	return isYES;	
}

//--------------------------------------------------------------------------------------------------
void CreateInfoFile()
{

	/**************************************************/
	// GETTING USER INFO
	/**************************************************/
    double MeanHMMass = std::accumulate(fHMMass.begin(), fHMMass.end(), 0) / fHMMass.size();
    string ReactorType, FuelType, Author,Mail,XSBase,HLCut,EnergyDisc,FPYBase,SABase,Geom,AddInfo,DepCode ;
    double Power = 0;
    ReadInfo(StringLine::ReplaceAll(JobName[0],".dat",".info" ), ReactorType,FuelType,Power);

	cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-++-+-+-+"<endl;
	cout<<"Data_Base_Info.nfo FILE GENERATOR"<<endl;
	cout<<"+-+-+-+-+-+-+-+-+-+-+-+-+-++-+-+-+"<endl;
	cout<<endl;

	cout<<"-> Author name : "<<endl;
	cin>> Author;
	cout<<endl;

	cout<<"-> email adress : "<<endl;
	cin>> Mail;
	cout<<endl;

	cout<<"-> Depletion code used : "<<endl;
	cin>> DepCode;
	cout<<endl;

	cout<<"-> Cross section data base (e.g ENSDF7.1) : "<<endl;
	cin>> XSBase;
	cout<<endl;

	cout<<"-> Fission Yield data base (e.g ENSDF7.1) : "<<endl;
	cin>> FPYBase;
	cout<<endl;

	cout<<"-> S(alpha,beta) data base (e.g ENSDF7.1) : "<<endl;
	cin>> SABase;
	cout<<endl;

	cout<<"-> Geometry simulated  (e.g Cubic Assembly with mirror boundary) : "<<endl;
	cin>> Geom;
	cout<<endl;

	cout<<"-> Half life cut [s] (if any) : "<<endl;
	cin>> HLCut;
	cout<<endl;

	cout<<"-> Multi group treatment (yes/no if yes number of groups) : "<<endl;
	cin>> EnergyDisc;
	cout<<endl;	

	cout<<"-> Additional informations : "<<endl;
	cin>> AddInfo;
	cout<<endl;	


	cout<<"-> Reactor type (e.g PWR, FBR,...) : "<<endl;
	cout<<"Found in a .info file :"<<endl;
	cout<<"\033[36m"<<ReactorType<<"\033[0m"<<endl;
	cout<< "Is that corect ? [y/n] "
	if(!UserSayYes())
	{	cout<<"\t So what it is ?"<<endl;
		cin>> ReactorType;
	}	
	cout<<endl;

	cout<<"-> Fuel type (e.g UOX, MOX,...) : "<<endl;
	cout<<"Found in a .info file :"<<endl;
	cout<<"\033[36m"<<FuelType<<"\033[0m"<<endl;
	cout<< "Is that corect ? [y/n] "
	if(!UserSayYes())
	{	cout<<"\t So what it is ?"<<endl;
		cin>> FuelType;
	}	
	cout<<endl;

	cout<<"-> Simulated heavy metal mass (tons) : "<<endl;
	cout<< "\t Calculated from your evolution datas the AVERAGE heavy metal mass is : "
	cout<<"\033[36m"<<MeanHMMass<<"\033[0m"<<" tons"<<endl;
	cout<< "Is that corect ? [y/n] "
	if(!UserSayYes())
	{	cout<<"\t So what it is ?"<<endl;
		cin>> MeanHMMass;
	}
	cout<<endl;

	cout<<"-> Simulated thermal power (W) : "<<endl;
	cout<<"Found in a .info file :"<<endl;
	cout<<"\033[36m"<<Power<<"\033[0m"<<endl;
	cout<< "Is that corect ? [y/n] "
	if(!UserSayYes())
	{	cout<<"\t So what it is ?"<<endl;
		cin>> Power;
	}	
	cout<<endl;

	/**************************************************/
	// BUILDING FILE
	/**************************************************/
	ofstream InfoFile("Data_Base_Info.nfo");

	InfoFile << "To be used with : XSM_MLP.cxx"<<endl;
	InfoFile << endl;
	InfoFile << "Reactor Type :"<<endl;
	InfoFile << "K_REACTOR "<< ReactorType <<endl;
	InfoFile << endl;
	InfoFile << "Fuel Type :"<<endl;
	InfoFile << "K_FUEL "<< FuelType <<endl;
	InfoFile << endl;
	InfoFile << "Heavy Metal [t] :"<<endl;
	InfoFile << "K_MASS  "<< MeanHMMass <<endl;
	InfoFile << endl;
	InfoFile << "Thermal Power [W] :"<<endl;
	InfoFile << "K_POWER  "<< Power <<endl;
	InfoFile << "Irradiation time steps [s] :"<<endl;
	InfoFile << "K_TIMESTEP";
		for( int t = 0 ; t < fNOfTimeStep ; t++ )
			InfoFile <<" "<<fTime[t];
	InfoFile << endl<<endl;
	InfoFile << "Z A I Name (input MLP) :"<<endl;
	for(map<ZAI,string> it = fMapName.begin() ; it != fMapName.end() ; it++ )
			InfoFile <<"K_ZAINAME "<<it->first.Z <<" " <<it->first.A <<" " <<it->first.I<<" " << it->second <<endl ;
	InfoFile <<endl;
	InfoFile << "Fuel range (Z A I min max) :"<<endl;
	for(map<ZAI,string> it = fMapName.begin() ; it != fMapName.end() ; it++ )
			InfoFile <<"K_ZAIL "<<it->first.Z <<" " <<it->first.A <<" " <<it->first.I<<" " << std::min_element(GetAllCompoOf(it->first))  << " "  <<std::max_element(GetAllCompoOf(it->first)) <<endl ;
  
  	InfoFile << endl;
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t ); 
    InfoFile <<" Date : "<< now->tm_mday<< '/'<<  (now->tm_mon + 1) << '/'  <<(now->tm_year + 1900)<< endl;

    string , FuelType, Author,Mail,XSBase,HLCut,EnergyDisc,FPYBase,SABase,Geom,AddInfo, DepCode ;



}
//--------------------------------------------------------------------------------------------------
vector<double> GetAllCompoOf(ZAI zai)
{
	vector<double> AllCompoOfZAI;

		for( int b = 0 ; b < fActinideCompoInit.size() ; b++ )
			AllCompoOfZAI.push_back(fActinideCompoInit[b].GetZAIIsotopicQuantity(zai)); 

return AllCompoOfZAI;
}

//--------------------------------------------------------------------------------------------------
void EvolutionData::ReadInfo(string InfoDBFile,string &ReactorType,string &FuelType,double &Power,)
{	
	ifstream InfoDB(InfoDBFile.c_str());				// Open the File
	if(!InfoDB)
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Can't open \"" << InfoDBFile << "\"\n" << endl;
	}
	else
	{
		int start = 0;
		string line;
		getline(InfoDB, line);
		if ( tlc(StringLine::NextWord(line, start, ' ')) ==  "reactor")
			ReactorType = StringLine::NextWord(line, start, ' ');
		
		start = 0;
		getline(InfoDB, line);
		if (tlc(StringLine::NextWord(line, start, ' ')) ==  "fueltype")
			FuelType = StringLine::NextWord(line, start, ' ');
	
		start = 0;
		getline(InfoDB, line);
		if ( tlc(StringLine::NextWord(line, start, ' ')) ==  "cycletime")
			fCycleTime = atof(StringLine::NextWord(line, start, ' ').c_str());
	
		getline(InfoDB, line); // Assembly HM Mass DONT TRUST THIS ONE CALCULATED WITH A instead of real atomic mass
	
		start = 0;
		getline(InfoDB, line);
		if ( tlc(StringLine::NextWord(line, start, ' ')) ==  "constantpower")
			Power = atof(StringLine::NextWord(line, start, ' ').c_str());
		InfoDB.close();
	}

}


//--------------------------------------------------------------------------------------------------
void DumpInputNeuron(string filename)
{
	TFile*   fOutFile = new TFile(filename.c_str(),"RECREATE");
	TTree*   fOutT = new TTree("Data", "Data");


/**********************INITIALISATIONNN********************/
	FillMapName();

	////////////////////////////////////////////////////////
	// INIT FRESH FUEL COMPOSITION and TIME
	////////////////////////////////////////////////////////

	double *FreshCompo = new double[fMapName.size()]; 

	for(int i = 0 ; i < fMapName.size() ; i++ )
		FreshCompo[i] = 0;

	////////////////////////////////////////////////////////
	double Time 		     = 0;	

/**********************init map********************/
	map < ZAI,vector<double> > mAllXS;
	map < ZAI, vector<double>  > mAllInventories;
	for(int act=0;act<int(fAllNuclei.size());act++ )	
	{	 
		vector<double>  InitVect;
		for(int Tstep=0 ;Tstep<fNOfTimeStep;Tstep++)
		{
			InitVect.push_back(0);
		}
		mAllInventories.insert( pair<ZAI,vector<double> >(fAllNuclei[act],InitVect) );
	}

	for(int act=0;act<int(fAllNuclei.size());act++ )	
	{	
		vector< double>  InitVect;
		for(int xs=0;xs<3;xs++)
		{
			
			InitVect.push_back(0);
				
		}

		mAllXS.insert(pair<ZAI,vector< double> >(fAllNuclei[act], InitVect) );
	}


/**********************BRANCHING**************************************************/
/**********************Fresh fuel**************************************************/

	map<ZAI,string>::iterator it;
	int index = 0;
	for(it = fMapName.begin() ; it != fMapName.end() ; it++ )
	{	string Name = it->second + "/D"; 
		fOutT->Branch( it->second.c_str() , &FreshCompo[index] , Name.c_str());
		index++;
	}

	////////////////////////////////////////////////////////
	fOutT->Branch(	"Time"			,&Time			,"Time/D"			);


/**********************cross section**************************************************/

	string XSType[3]={"fis","cap","n2n"}; //!!!!!!DO NOT TOUCH THIS

	for(int act=0;act<int(fAllNuclei.size());act++ )
	{	
		for(int xs=0;xs<3;xs++)
		{
				string NamedXS= NameXS(fAllNuclei[act],XSType[xs]);
				string NameXSBis=NamedXS+"/D";

		if(	xs==0 && fXSFis[0][fAllNuclei[act]].size()!=0)	//Check if the reaction we want to branch exists
			 if(fXSFis[0][fAllNuclei[act]][0]!=0) 
			 		fOutT->Branch(	NamedXS.c_str() ,  &mAllXS[fAllNuclei[act]][xs],   NameXSBis.c_str()	);			

		if(xs==1 && fXSCap[0][fAllNuclei[act]].size()!=0)
			 if(fXSCap[0][fAllNuclei[act]][0]!=0) 
			 		fOutT->Branch(	NamedXS.c_str() ,  &mAllXS[fAllNuclei[act]][xs],   NameXSBis.c_str()	);				

		if(xs==2 && fXSN2N[0][fAllNuclei[act]].size()!=0)
			 if(fXSN2N[0][fAllNuclei[act]][0]!=0) 
					fOutT->Branch(	NamedXS.c_str() ,  &mAllXS[fAllNuclei[act]][xs],   NameXSBis.c_str()	);				


		}
	}
/**********************FILLING THE TTREE**************************************************/

	//Fill containing all the output of the networks to train
	 ofstream  InputNetwork("TrainingInput.cxx");

	 int NumOfBase=fActinideCompoInit.size();
	for(int b=0;b<NumOfBase;b++) 
	{ 
		///////////////////////////////////////////////////////
		int index =0;
		for(it = fMapName.begin() ; it != fMapName.end() ; it++ )
		{	
			int Z = it->first.Z;
			int A = it->first.Z;
			int I = it->first.Z;
			FreshCompo[index] = fActinideCompoInit[b].GetZAIIsotopicQuantity(ZAI(Z,A,I));

			index++;
		}	
		///////////////////////////////////////////////////////
	
			for(int Tstep=0 ;Tstep<fNOfTimeStep;Tstep++ )	
			{	
	 			Time=fTime[Tstep];
				for(int act=0;act<int(fAllNuclei.size());act++)
				{
	
					if(fXSFis[b][fAllNuclei[act]].size()!=0) // reaction may not be present
					{	
						if(fXSFis[b][fAllNuclei[act]][Tstep] !=0 )
						{	mAllXS[fAllNuclei[act]][0]   = 	fXSFis[b][fAllNuclei[act]][Tstep];
							if(b==0 && Tstep==0)
								InputNetwork<<"OUTPUT.push_back(\""<<NameXS(fAllNuclei[act],XSType[0])<<"\");"<<endl;
						}
					}		
					if(fXSCap[b][fAllNuclei[act]].size()!=0)
				 	{	
				 		if(fXSCap[b][fAllNuclei[act]][Tstep] !=0) 
				 		{	mAllXS[fAllNuclei[act]][1]   = 	fXSCap[b][fAllNuclei[act]][Tstep];
							if(b==0 && Tstep==0)		 			
				 				InputNetwork<<"OUTPUT.push_back(\""<<NameXS(fAllNuclei[act],XSType[1])<<"\");"<<endl;
						}
				 	}	
				 	if( fXSN2N[b][fAllNuclei[act]].size()!=0)
				 	{	
				 		if(fXSN2N[b][fAllNuclei[act]][Tstep] !=0)
				 		{	mAllXS[fAllNuclei[act]][2]   = 	fXSN2N[b][fAllNuclei[act]][Tstep];
							if(b==0 && Tstep==0)
				 				InputNetwork<<"OUTPUT.push_back(\""<<NameXS(fAllNuclei[act],XSType[2])<<"\");"<<endl;
						}
				 	}
				 	
				 } 
				 fOutT->Fill();	
			}		
	}

	fOutFile->Write();
	delete fOutT;
	fOutFile-> Close();
	delete fOutFile;

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

	//XS  
	 map<ZAI, vector <double> > mapFistmp;
	 map<ZAI, vector <double> > mapCaptmp;
	 map<ZAI, vector <double> > mapN2Ntmp;
	do
	{
		
		start = 0;
		int z;
		string tmp2 = StringLine::NextWord(line, start, ' ');
		if (tmp2 == "XSFis") 
		{
			z = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int a = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int i = atoi(StringLine::NextWord(line, start, ' ').c_str());
		
			if(a!=0 && z!=0)
			{
				
				ZAI zaitmp(z, a, i);

				vector<double> XSTime;
				for(int i = 0; i < (int)vT.size(); i++)
				{	long double q = atof(StringLine::NextWord(line, start, ' ').c_str());
					XSTime.push_back(q);
				}
				mapFistmp.insert( pair<ZAI,vector<double> >(zaitmp,XSTime));		
			}
		}
		else if (tmp2 == "XSCap") 
		{
			z = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int a = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int i = atoi(StringLine::NextWord(line, start, ' ').c_str());
		
			if(a!=0 && z!=0)
			{
				
				ZAI zaitmp(z, a, i);
				vector<double> XSTime;
				for(int i = 0; i < (int)vT.size(); i++)
				{	long double q = atof(StringLine::NextWord(line, start, ' ').c_str());
					XSTime.push_back(q);
				}
				mapCaptmp.insert( pair<ZAI,vector<double> >(zaitmp,XSTime));	
				
			}
		}	
		else if (tmp2 == "XSn2n") 
		{
			z = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int a = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int i = atoi(StringLine::NextWord(line, start, ' ').c_str());
		
			if(a!=0 && z!=0)
			{		
				ZAI zaitmp(z, a, i);
				vector<double> XSTime;
				for(int i = 0; i < (int)vT.size(); i++)
				{	long double q = atof(StringLine::NextWord(line, start, ' ').c_str());
					XSTime.push_back(q);
				}
				mapN2Ntmp.insert( pair<ZAI,vector<double> >(zaitmp,XSTime));					
			}
		}
		getline(DecayDB, line);
		start = 0;
		tmp2 = StringLine::NextWord(line, start, ' ');
		
		if(line == "") break;
	}while (!DecayDB.eof() );
	
	fXSFis.push_back(mapFistmp);
	fXSCap.push_back(mapCaptmp);
	fXSN2N.push_back(mapN2Ntmp);
	
	DecayDB.close();
	
	double N = 0;
	for(int i=0; i < (int)Z.size()-2; i++)
	{
		if(  Z[i]>89 ) 
			N += Q[i];
	}
	

	IsotopicVector CompoBasei; 
	IsotopicVector CompoBaseiUnormalize;

	for(int i=0; i < (int)Z.size()-2; i++)
	{
		if(Z[i]>89)
		{
			ZAI zai = ZAI(Z[i], A[i], I[i]);
			CompoBasei.IVquantity.insert(pair<ZAI,double>(zai,Q[i]/N));
			CompoBaseiUnormalize.IVquantity.insert(pair<ZAI,double>(zai,Q[i]));
		}	
	}
	
	fActinideCompoInit.push_back(CompoBasei);
	fHMMass.push_back(cZAIMass(CompoBaseiUnormalize));

GoodJobName.push_back(jobname);

}
/*-------------------------------------------------------------------------------------------------
COMPILATION :

g++ -o Gene Gene.cxx `root-config --cflags` `root-config --libs`

*/