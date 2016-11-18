//	DESCRIPTION
// This program convert MURE output tp EvolutionData forma
//
//@author BaM
#include "BinaryFormat2.hxx"
#include "StringLine.hxx"
#include "ZAI.hxx"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using std::string;
using namespace std;

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
//------------------------------------- METHODS ------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------

string itoa(int num)
{
	ostringstream os(ostringstream::out);
	os << "" << num;
	return os.str();
}

int ReadCommentVersion(ifstream& in, string filename)
{
	if(!in.good())
	{
		cerr << "ERROR in MureRead:: Cannot find ASCII file " << filename <<  " to Read" << endl;
		return -1;
	}
	//
	//Go throught the comments at the beginning of the file and Read the version number of Data Format
	//
	string TheComment;
	int MureDataVersion;
	bool bComment=true;
	while (bComment)
	{
		in>>TheComment;
		if (TheComment=="%")
		{
			getline(in,TheComment); // the "in>>" command keep the cursor at the end of line
		}
		else if (TheComment=="V")
		{
			in>>MureDataVersion;
			bComment=false;
		}
		else
		{
			cout << "No comments or format version at the beginning of data files" << endl;
			bComment=false;
			MureDataVersion=0;
		}
	}
	return MureDataVersion;
}

double GetPropOf( vector <ZAI> vZAI , ZAI zai)
{

	for(int nuc = 0 ; nuc < vZAI.size() ; nuc++ )
	{
		if( vZAI[nuc].Z() == zai.Z() &&  vZAI[nuc].A() == zai.A() && vZAI[nuc].I() == zai.I() )
			return vZAI[nuc].Prop();
	}

return 0.0;
}

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
//------------------------------------- MAIN ---------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------



int main(int argc, char** argv)
{

	// =========================================================================================
	//  VARIABLES
	// =========================================================================================
	string sPWD = getenv("PWD");

	string s_tmp = ""; int i_tmp = 0; double tmp = 0;
	string DataASCIIFile = "DATA_";	 int NDataASCIIFile = 0;
	string DataBinFile   = "BDATA_"; int NDataBinFile   = 0;
	string DataFile		 = "";		 int NDataFile	    = 0;
	string LastDataFile  = "";

	string OutDataFile	 = "";
	string OutDataFileInfo	 = "";
	string OutLOG	 = "";

	bool IsThereASCCIFile = false;
	bool IsThereBinFile   = false;

	vector <int>    vCellNumber;
	vector <string> vCellComment;

	int NumberOfCells = 1;

	vector < double > vTime;
	vector < double > vKeff;
	vector < double > vFlux;
	vector < vector < double > > vFlux1;

	vector <ZAI >                    zai1;
	vector <vector <ZAI > >          zai2;
	vector <vector <vector<ZAI > > > zai3;

	vector <ZAI >                    zai_FS_1;
	vector <vector <ZAI > >          zai_FS_2;
	vector <vector <vector<ZAI > > > zai_FS_3;


	vector <ZAI >                    zai_SC_1;
	vector <vector <ZAI > >          zai_SC_2;
	vector <vector <vector<ZAI > > > zai_SC_3;


	vector <ZAI >                    zai_SN_1;
	vector <vector <ZAI > >          zai_SN_2;
	vector <vector <vector<ZAI > > > zai_SN_3;


	vector <bool >                   HasToBePrint;

	// =========================================================================================
	//  WELCOME !
	// =========================================================================================

	// MURE output path
	if(argc < 4)
	{
		cout << "MURE2CLASS usage: " << endl;
		cout << "Arguments should be :" << endl;
		cout << "\t1 Path," << endl;
		cout << "\t2 ReactorType," << endl;
		cout << "\t3 FuelType," << endl;
		cout << "Optional : "<<endl ;
		cout << "\t4 WantedCell [Default :0]" <<  endl;
		cout << "\t5 Step to Skip [Default :0]" << endl;
		exit(1);
	}
	string DBPath 			= argv[1];
	string ReactorType		= argv[2];
	string FuelType			= argv[3];
	int WantedCell = 0;
	int StepToSkip = 0;
	if(argc == 5)
		WantedCell = atoi(argv[4]);
	if(argc == 6)
		StepToSkip = atoi(argv[5]);

	// Name of the output file
	if (DBPath.back() == '/')
		DBPath.pop_back();

 	size_t found = DBPath.find_last_of("/");
  	string OutName = DBPath.substr(found+1);
	OutDataFile     = "./" + OutName + ".dat";
	OutDataFileInfo = "./" + OutName + ".info";
	OutLOG          = "./" + OutName + ".log";

	ofstream OutputLog(OutLOG.c_str());

	// Check how many DATAxxx and BDATAxxx files there are
	string Command = "find " + DBPath + " -name \"" + DataASCIIFile + "*\" > ASCII.tmp";
	system(Command.c_str());
	Command = "find " + DBPath + " -name \"" + DataBinFile   + "*\" > BIN.tmp";
	system(Command.c_str());

	ifstream ASCIITMP("ASCII.tmp");
	if (ASCIITMP.is_open())
	{
		while (!ASCIITMP.eof())
		{
			getline(ASCIITMP,s_tmp);
			i_tmp++;
		}
	} ASCIITMP.close(); i_tmp--; NDataASCIIFile = i_tmp;

	i_tmp=0;
	ifstream BINTMP("BIN.tmp");
	if (BINTMP.is_open())
	{
		while (!BINTMP.eof())
		{
			getline(BINTMP,s_tmp);
			i_tmp++;
		}
	} BINTMP.close(); i_tmp--; NDataBinFile = i_tmp;

	// Remove temporary files...
	Command = "\\rm -f ASCII.tmp BIN.tmp";
	system(Command.c_str());

	// Say if there is Binary or ASCII file
	if      (NDataASCIIFile>0) {IsThereASCCIFile=true; DataFile = DataASCIIFile; NDataFile = NDataASCIIFile;}
	else if (NDataBinFile>0)   {IsThereBinFile  =true; DataFile = DataBinFile;   NDataFile = NDataBinFile;}
	else {cout << endl << "There is no MURE DATA file in the given path... EXIT!" << endl << endl; exit(1);}

	// =========================================================================================
	//  READ (B)DATA files...
	// =========================================================================================

	OutputLog << endl;
	OutputLog << "===================================================" << endl;
	OutputLog << "---------------------------------------------------" << endl;
	if (IsThereBinFile) OutputLog << endl << "Binary file detected... " << endl << endl;
	else cout << endl << "ASCII file detected... " << endl << endl;
	//	/*sleep(1)*/;

	for (int t=0; t<NDataFile; t++)
	{
		string SFX = "";
		if      (t<=9)				SFX = "00" + itoa(t);
		else if (t>=10 && t<=99)	SFX = "0" + itoa(t);
		else if (t>=100 && t<=999) 	SFX = "" + itoa(t);
		else {cout << endl << "there is more than 1000 DATA files, not yet implemented... EXIT!" << endl << endl; exit(1);}
		string s_File = DBPath + "/" + DataFile + SFX;
		if (t%10==0) {string ff = DataFile + SFX; cout << "Reading file " << ff << endl; /*sleep(1)*/;}
		LastDataFile = DataFile + SFX;

		// READ BINARY FILE
		if (IsThereBinFile)
		{
			ifstream in(s_File.c_str(),ios::binary);
			if(!in.good()){cout << "Cannot find Binary file " << s_File << endl; exit(0);}
			//
			//Read The version number of writing Mure evolving Data Format
			//
			char TextV;
			int MureDataVersion=0;
			in.read((char*)&TextV, sizeof(char));

			if(TextV=='V') in.read((char*)&MureDataVersion, sizeof(int));
			else {in.close(); in.clear(); in.open(s_File.c_str(),ios::binary);}

			string TheComment;

			// Read the File Header
			FileHeader FH;
			FH.read(in);
			NumberOfCells = FH.NCells;

			// --------------------------- TIME ----------------------
			vTime.push_back(FH.Time);
			vKeff.push_back(FH.K);

			CellHeader CH;
			NucleusRecord NR;
			bool STOP2READ=false;

			for(int c=0 ; c<NumberOfCells; c++)
			{
				if(STOP2READ)break;
				CH.read(in);
				vFlux.push_back(CH.Flux);

				// read the Cell Comment
				int StringSize; in.read((char*)&StringSize, sizeof(int)); TheComment.resize(StringSize);
				char tmp[StringSize+1];	in.read(tmp, (StringSize+1)*sizeof(char)); tmp[StringSize]='\0';
				TheComment=tmp;

				bool SkipCell=false;
				if(CH.CellNumber<0 && TheComment!="WasteStorage") SkipCell=true;
				if(TheComment=="WasteStorage") STOP2READ=true;

				if (t==0)
				{
					vCellNumber.push_back(CH.CellNumber);
					vCellComment.push_back(TheComment);
				}

				//read the spatial variables
				int NSpatialVariables; in.read((char*)&NSpatialVariables, sizeof(int));

				vector<double> SpatialVariables(NSpatialVariables);
				vector<string> SpatialVariableNames;

				for(int j=0; j<NSpatialVariables; j++) in.read((char*)&SpatialVariables[j],sizeof(double));

				//read the spatial variable names
				for(int j=0; j<NSpatialVariables; j++)
				{
					int StringSize; in.read((char*)&StringSize,sizeof(int));
					s_tmp=""; s_tmp.resize(StringSize); char tmp[StringSize+1];
					in.read(tmp, (StringSize+1)*sizeof(char)); tmp[StringSize]='\0';
					s_tmp=tmp; SpatialVariableNames.push_back(s_tmp);
				}

				for(int i=0; i<CH.NNucleusRecords; i++)
				{
					// Read the Nucleus Record
					NR.read(in);

					if(!SkipCell)
					{
						ZAI zai(NR.Z, NR.A, NR.I, NR.Proportion);
						zai1.push_back(zai);
						//cout<<zai.Z()<<" "<<zai.A()<<" "<<zai.I()<<" "<<zai.Prop()<<endl;
						//if (t==0 && NR.Proportion>=2e+25) cout << NR.Z << "  " << NR.A << "  " << NR.I << "  " << NR.Proportion << endl;
					}

					ReactionRecord RR;
					for (int k=0; k<NR.NReactionRecords; k++)
					{
						RR.read(in);
						if(RR.Sigma==0.0) RR.Sigma=1e-10;

						if(!SkipCell)
						{
							ZAI zai(NR.Z, NR.A, NR.I, RR.Sigma);
							if (RR.Code==18)  zai_FS_1.push_back(zai);
							if (RR.Code==102) zai_SC_1.push_back(zai);
							if (RR.Code==16)  zai_SN_1.push_back(zai);

						}





					}
				}
				if(!SkipCell)
				{
					zai2.push_back(zai1);
					zai1.clear();
					zai_FS_2.push_back(zai_FS_1);
					zai_FS_1.clear();
					zai_SC_2.push_back(zai_SC_1);
					zai_SC_1.clear();
					zai_SN_2.push_back(zai_SN_1);
					zai_SN_1.clear();
				}
			}
			in.close();
			zai3.push_back(zai2); zai2.clear();
			zai_FS_3.push_back(zai_FS_2); zai_FS_2.clear();
			zai_SC_3.push_back(zai_SC_2); zai_SC_2.clear();
			zai_SN_3.push_back(zai_SN_2); zai_SN_2.clear();

			vFlux1.push_back(vFlux); vFlux.clear();

		}
		// READ ASCII FILE
		else if (IsThereASCCIFile)
		{
			ifstream in(s_File.c_str());
			int TranspCode=1;
			if(!in.good()) {if(!in.good()){OutputLog << "Cannot find ASCII file " << s_File << endl; exit(0);}}
			//
			//Read The version number of writing Mure evolving Data Format
			//
			int MureDataVersion=ReadCommentVersion(in, s_File);
			if(MureDataVersion==0) {in.close(); in.clear(); in.open(s_File.c_str());}
			else if(MureDataVersion<0) {OutputLog << "Old version of MURE... Check!!! exit(1)"; exit(0);}

			// Read out the TIME the KEFF and the KEFF ERROR
			double Time,Keff,Keff_Err;
			in>>Time>>Keff>>Keff_Err;
			in>>NumberOfCells;

			// --------------------------- TIME ----------------------
			vTime.push_back(Time);
			vKeff.push_back(Keff);

			bool STOP2READ=false;

			for(int c=0 ; c<NumberOfCells; c++)
			{
				if(STOP2READ)break;
				int CellNumber;
				double Volume,Flux,FluxErr;
				in>>CellNumber>>Volume>>Flux>>FluxErr;

				vFlux.push_back(Flux);

				int NNucleusRecords;
				in>>NNucleusRecords;

				// read the Cell Comment and Spatial Variables
				string TheComment;
				getline(in,TheComment); // the "in>>" command keep the cursor at the end of line
				getline(in,TheComment); // so two "getline" commands are requested

				bool SkipCell=false;
				if(CellNumber<0 && TheComment!="WasteStorage") SkipCell=true;
				if(TheComment=="WasteStorage") STOP2READ=true;

				if (t==0)
				{
					vCellNumber.push_back(CellNumber);
					vCellComment.push_back(TheComment);
				}

				int NSpatialVariables;
				in>>NSpatialVariables;

				vector<double> SpatialVariables;
				vector<string> SpatialVariableNames;

				for(int j=0; j<NSpatialVariables; j++)
				{
					string TheString;
					double val;
					in>>TheString>>val;
					SpatialVariableNames.push_back(TheString);
					SpatialVariables.push_back(val);
				}

				for(int i=0; i<NNucleusRecords; i++)
				{
					// Read the Nucleus Record
					int Z,A,I;
					double Mass,Proportion;
					in>>Z>>A>>I>>Mass>>Proportion;
					cout << Z << endl;
					if(!SkipCell)
					{
						ZAI zai(Z, A, I, Proportion);
						zai1.push_back(zai);
					}

					int NReactionRecords;
					in>>NReactionRecords;
					for (int k=0; k<NReactionRecords; k++)
					{
						int Code;
						double Sigma,SigmaErr;
						in>>Code>>Sigma>>SigmaErr;
						if(!SkipCell)
						{
							ZAI zai(Z, A, I, Proportion);
							if (Code==18)  zai_FS_1.push_back(zai);
							if (Code==102) zai_SC_1.push_back(zai);
							if (Code==16)  zai_SN_1.push_back(zai);

						}

						if(Sigma==0.0) Sigma=1e-10;
					}




				}
				if(!SkipCell)
				{
					zai2.push_back(zai1);
					zai1.clear();
					zai_FS_2.push_back(zai_FS_1);
					zai_FS_1.clear();
					zai_SC_2.push_back(zai_SC_1);
					zai_SC_1.clear();
					zai_SN_2.push_back(zai_SN_1);
					zai_SN_1.clear();
				}
			}
			in.close();
			zai3.push_back(zai2); zai2.clear();
			zai_FS_3.push_back(zai_FS_2); zai_FS_2.clear();
			zai_SC_3.push_back(zai_SC_2); zai_SC_2.clear();
			zai_SN_3.push_back(zai_SN_2); zai_SN_2.clear();
			vFlux1.push_back(vFlux); vFlux.clear();


		}
		// else ERROR
		else {OutputLog << endl << "DATA files are neither binary nor ASCII?? ... EXIT!" << endl; exit(1);}
	}
	OutputLog << endl << "Last file read : " << LastDataFile << endl << endl;
	OutputLog << "All (B)DATA files have been read..." << endl << endl;
	OutputLog << "---------------------------------------------------" << endl;
	OutputLog << "===================================================" << endl;

	// =========================================================================================
	//  Manage several Cells
	// =========================================================================================

	// Time bins
	int Nt = zai3.size();
	// Number of cells
	int Nc = zai3[0].size();

	// Number of nuclides
	int Ni = zai3[0][0].size();

	double CycleTime = (vTime[vTime.size()-1] - vTime[0]) / 3600. / 24 / 365.4;

	// Manage the cells...
	bool SumOfCell = false;

	if (vCellNumber.size()>=2)
	{
		OutputLog << endl;
		OutputLog << "===================================================" << endl;
		OutputLog << "-------------- WARNING ----------------------------" << endl;
		OutputLog << "===================================================" << endl << endl;
		OutputLog << "THERE IS MORE THAN ONE CELL... Cells are : " << endl << endl; /*sleep(1)*/;
		for(int i=0; i<vCellNumber.size(); i++)
		{
			OutputLog << "index : " << i << " - Cell number  : " << vCellNumber[i] << endl;
			OutputLog << "Cell comment : " << vCellComment[i] << endl << endl;
		}
		if(WantedCell==0) SumOfCell = true;

		if (!SumOfCell && WantedCell>Nc) {OutputLog << "The index you choose is not defined... EXIT!" << endl; exit(1);}

		zai2.clear();
		zai_FS_2.clear();
		zai_SC_2.clear();
		zai_SN_2.clear();
		vFlux.clear();
		for (int t=0; t<Nt; t++)
		{
			for (int i=0; i<Ni; i++)
			{
				double SUM = 0;
				for (int c=0; c<Nc; c++) SUM += zai3[t][c][i].Prop();
				ZAI zai(zai3[t][0][i].Z(),zai3[t][0][i].A(),zai3[t][0][i].I(),SUM);
				zai1.push_back(zai);
			}
			for (int i=0; i<(int)zai_FS_3[0][0].size(); i++)
			{
				double SUM = 0;
				for (int c=0; c<Nc; c++) SUM += zai_FS_3[t][c][i].Prop()/Nc;
				ZAI zai(zai_FS_3[t][0][i].Z(),zai_FS_3[t][0][i].A(),zai_FS_3[t][0][i].I(),SUM);
				zai_FS_1.push_back(zai);
			}
			for (int i=0; i<(int)zai_SC_3[0][0].size(); i++)
			{
				double SUM = 0;
				for (int c=0; c<Nc; c++) SUM += zai_SC_3[t][c][i].Prop()/Nc;
				ZAI zai(zai_SC_3[t][0][i].Z(),zai_SC_3[t][0][i].A(),zai_SC_3[t][0][i].I(),SUM);
				zai_SC_1.push_back(zai);
			}
			for (int i=0; i<(int)zai_SN_3[0][0].size(); i++)
			{
				double SUM = 0;
				for (int c=0; c<Nc; c++) SUM += zai_SN_3[t][c][i].Prop()/Nc;
				ZAI zai(zai_SN_3[t][0][i].Z(),zai_SN_3[t][0][i].A(),zai_SN_3[t][0][i].I(),SUM);
				zai_SN_1.push_back(zai);
			}

			double SUM = 0;
			for (int j=0; j<(int)vFlux1[t].size(); j++) SUM += vFlux1[t][j];

			vFlux.push_back(SUM);

			zai2.push_back(zai1);
			zai1.clear();
			zai_FS_2.push_back(zai_FS_1);
			zai_FS_1.clear();
			zai_SC_2.push_back(zai_SC_1);
			zai_SC_1.clear();
			zai_SN_2.push_back(zai_SN_1);
			zai_SN_1.clear();

		}
		OutputLog << endl;
		OutputLog << "---------------------------------------------------" << endl;
		OutputLog << "-------------- Sum of cells done ------------------" << endl;
		OutputLog << "---------------------------------------------------" << endl << endl;
	}
	else WantedCell=0;


	// =========================================================================================
	//  Manage the cutoff and calculate total N and M at t=0
	// =========================================================================================

	double NTotal = 0;
	double MTotalFissile = 0;
	if (SumOfCell)
	{
		for (int i=0; i<Ni; i++)
		{
			double Ni = zai2[0][i].Prop();
			NTotal += Ni;
			if (zai2[0][i].Z()>=90) MTotalFissile += Ni*zai2[0][i].A()/6.023e23;
		}
	}
	else
	{
		for (int i=0; i<Ni; i++)
		{
			double Ni = zai3[0][WantedCell][i].Prop();
			NTotal += Ni;
			if (zai3[0][WantedCell][i].Z()>=90) MTotalFissile += Ni*zai3[0][WantedCell][i].A()/6.023e23;
		}
	}

	// =========================================================================================
	//  Write the DATABASE and convert Number to Mass - Manage Normalization Factor
	// =========================================================================================

	ofstream Output(OutDataFile.c_str());
	Output.precision(16);
	Output << "time";
	for(int t=StepToSkip; t<vTime.size(); t++) Output << " " << vTime.at(t)-vTime.at(StepToSkip);
	Output << endl;

	Output << "keff";
	for(int t=StepToSkip; t<vKeff.size(); t++) Output << " " << vKeff.at(t);
	Output << endl;

	Output << "flux";
	if (SumOfCell)
	{
		for(int t=StepToSkip; t<vTime.size(); t++) Output << " " << vFlux.at(t);
		Output << endl;
	}
	else
	{
		for(int t=StepToSkip; t<vTime.size(); t++) Output << " " << vFlux1[t][WantedCell];
		Output << endl;
	}
	int NPrinted=0;
	{
		if (SumOfCell)
		{
			for(int i=0; i < zai2[0].size(); i++)
			{
				Output << "Inv " << zai2[0][i].Z() << " " << zai2[0][i].A() << " " << zai2[0][i].I() << " ";
				for (int t=StepToSkip; t<vTime.size(); t++)
				{
					double Val = zai2[t][i].Prop() ;
					Output << Val << " ";
					if (t==StepToSkip) NPrinted++;

				}
				Output << endl;
			}
		}
		else
		{
			for(int i=0; i < zai3[0][0].size(); i++)
			{	
				Output << "Inv " << zai3[0][0][i].Z() << " " << zai3[0][0][i].A() << " " << zai3[0][0][i].I() << " ";
				for (int t=StepToSkip; t<vTime.size(); t++)
				{
					double Val = zai3[t][WantedCell][i].Prop() ;
					Output << Val << " ";
					if (t==StepToSkip) NPrinted++;
				}
				Output << endl;
			}
		}
	}

	if (SumOfCell)
	{

		for(int i=0; i< (int)zai_FS_2[0].size(); i++)
		{

			Output << "XSFis " << zai_FS_2[0][i].Z() << " " << zai_FS_2[0][i].A() << " " << zai_FS_2[0][i].I() << " ";
			for (int t=StepToSkip; t<vTime.size(); t++)
			{
				double Val = zai_FS_2[t][i].Prop();
				Output << Val << " ";

			}
			Output << endl;
		}
	}
	else
	{
		for(int i=0; i< (int)zai_FS_3[0][0].size(); i++)
		{

			Output << "XSFis " << zai_FS_3[0][0][i].Z() << " " << zai_FS_3[0][0][i].A() << " " << zai_FS_3[0][0][i].I() << " ";
			for (int t=StepToSkip; t<vTime.size(); t++)
			{
				double Val = zai_FS_3[t][WantedCell][i].Prop();
				Output << Val << " ";
			}
			Output << endl;
		}
	}

	if (SumOfCell)
	{
		for(int i=0; i< (int)zai_SC_2[0].size(); i++)
		{
			Output << "XSCap " << zai_SC_2[0][i].Z() << " " << zai_SC_2[0][i].A() << " " << zai_SC_2[0][i].I() << " ";
			for (int t=StepToSkip; t<vTime.size(); t++)
			{
				double Val = zai_SC_2[t][i].Prop();
				Output << Val << " ";

			}
			Output << endl;
		}
	}
	else
	{
		for(int i=0; i< (int)zai_SC_3[0][0].size(); i++)
		{
			Output << "XSCap " << zai_SC_3[0][0][i].Z() << " " << zai_SC_3[0][0][i].A() << " " << zai_SC_3[0][0][i].I() << " ";
			for (int t=StepToSkip; t<vTime.size(); t++)
			{
				double Val = zai_SC_3[t][WantedCell][i].Prop();
				Output << Val << " ";
			}
			Output << endl;
		}
	}

	if (SumOfCell)
	{
		for(int i=0; i<  (int)zai_SN_2[0].size(); i++)
		{
			Output << "XSn2n " << zai_SN_2[0][i].Z() << " " << zai_SN_2[0][i].A() << " " << zai_SN_2[0][i].I() << " ";
			for (int t=StepToSkip; t<vTime.size(); t++)
			{
				double Val = zai_SN_2[t][i].Prop();
				Output << Val << " ";

			}
			Output << endl;
		}
	}
	else
	{
		for(int i=0; i<  (int)zai_SN_3[0][0].size(); i++)
		{
			Output << "XSn2n " << zai_SN_3[0][0][i].Z() << " " << zai_SN_3[0][0][i].A() << " " << zai_SN_3[0][0][i].I() << " ";
			for (int t=StepToSkip; t<vTime.size(); t++)
			{
				double Val = zai_SN_3[t][WantedCell][i].Prop();
				Output << Val << " ";
			}
			Output << endl;
		}
	}
	Output.close();


	double Power = 0 ;//Calculated Power @ t = 0 (as we supposed it is constant)
	for(int c = 0 ; c < NumberOfCells ; c++)
	{
		for(int nuc = 0 ; nuc < zai_FS_3[0][c].size() ; nuc++)		//ef * Nf * sigma_f *phi
		{

			double EnergyPerFisison_f = 1.9679e6*zai_FS_3[0][c][nuc].A()-2.601e8; //eV
			double XS_f = zai_FS_3[0][c][nuc].Prop() * 1e-24;
			double N_f  = GetPropOf( zai3[0][c], zai_FS_3[0][c][nuc] );

			Power +=   EnergyPerFisison_f * N_f * XS_f * vFlux1[0][c]; //eV.s-1

			cout << "Z " <<zai_FS_3[0][c][nuc].Z() << " A " << zai_FS_3[0][c][nuc].A() <<" I " << zai_FS_3[0][c][nuc].I() <<endl;
			cout << "\t "<< EnergyPerFisison_f <<" eV " <<endl;
			cout << "\t "<< XS_f << "cm^2"<<endl;
			cout << "\t "<< N_f << " at"<<endl;

		}
	}
	Power *= 1.60218e-19 ; //W

	ofstream OutputInfo(OutDataFileInfo.c_str());
	OutputInfo << "Reactor " << ReactorType << endl;
	OutputInfo << "Fueltype " << FuelType << endl;
	OutputInfo << "CycleTime " << CycleTime << endl;
	OutputInfo << "AssemblyHeavyMetalMass " << MTotalFissile << " g" << endl;
	OutputInfo << "ConstantPower " << Power << " W" << endl;

	OutputInfo << "Nnuclei " << NPrinted << endl;

	OutputInfo.close();

	// =========================================================================================
	//  BYE!
	// =========================================================================================

	OutputLog << "===================================================" << endl;
	OutputLog << "---------------------------------------------------" << endl;

	OutputLog << endl << "The database " << OutDataFile << " has been generated..." << endl;
	OutputLog << "The database information " << OutDataFileInfo << " has been generated..." << endl << endl;
	OutputLog << NPrinted << " nuclides have been written" << endl << endl;

	OutputLog << "---------------------------------------------------" << endl;
	OutputLog << "===================================================" << endl;
	
}


/*
 g++ -std=c++11 -o MURE2CLASS MURE2CLASS.cxx
 */
