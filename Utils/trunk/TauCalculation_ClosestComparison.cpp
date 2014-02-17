using namespace std;


#include <cmath>
#include <math.h>
#include <vector>
//#include <TMath.h>
#include <iostream>
#include <TGraphErrors.h>
#include <TFile.h>
#include <ctime>
#include <TH1F.h>
#include <TH2F.h>
#include "CLASSHeaders.hxx"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <math.h>
using namespace std;

string itoa(int num)
{
	ostringstream os(ostringstream::out);
	os<<""<<num;
	return os.str();
}

int CurveColor(int graph_num)
{
	int ColorTable[]={1,kBlue,kRed,kGreen+1,kMagenta,kCyan+2,kOrange-3,kRed+2,kBlue-2,
		kSpring+9,kGreen+3,kAzure+8,kMagenta+2,kYellow+2,kBlue-9,kOrange+2};
	return ColorTable[graph_num%10];

}

int main(int argc, char** argv)
{




	map<ZAI, double> fZAIMass;





	fZAIMass.insert( pair< ZAI,double >( ZAI(92,238,0), 238050788.247e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(92,235,0), 235043929.918e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(94,238,0), 238049559.894e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(94,239,0), 239052163.381e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(94,240,0), 240053813.545e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(94,241,0), 241056851.456e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(94,242,0), 242058742.611e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(95,241,0), 241056829.144e-6 ) );



	LogFile  *logfile = new LogFile("CLASS_log");
	DataBank<IsotopicVector> MOX = DataBank<IsotopicVector>(logfile, "../DB.idx", true, false);
	DataBank<IsotopicVector> MOX_new = DataBank<IsotopicVector>();


	TFile *outfile = new TFile ("CLASSresult.root", "RECREATE");

	map<IsotopicVector ,EvolutionData > evolutiondataBase = MOX.GetDataBank();

	map<IsotopicVector ,EvolutionData >::iterator it;

	int loop = 0;
	int start = time(NULL);

	TH1F* Tau_Phi[10];
	TH1F* Tau_N[10];
	TH1F* Tau_N_Total[10];

	TH2F* Tau_N_Indiv[9][10];
	TH2F* Tau_N_Full[10];

	string NucleiName[9] = {"U238", "U235", "Pu238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Np237"};
	string SpeciesName[10] = {"U", "Pu", "Am", "Np", "Cm", "Th", "Pa", "Actinide", "PF", "Other"};
	TH2F* Tau_N_Indiv_Species[10][10];


	string name;
	string TitleName;

	double TauPhiTot[10];
	double TauNTot[10];
	double TauN_Nico[10] ;
	double TauN_MAX[10] ;



	for( int i = 0; i < 10; i++)
	{

		TauPhiTot[i] = 0;
		TauNTot[i] = 0;

		int bin_number = 140;
		double x[bin_number];
		double dE = exp( 1./bin_number * log( 1./1e-7) );
		for(int k=0 ; k<=bin_number; k++)
			x[k]=pow(dE,k)*1e-7;

		for (int j =0; j < 9; j++)
		{
			name = "TauN_" + NucleiName[j] + "_Step_" + itoa(i);
			TitleName = "#tau_{N} " + NucleiName[j] + " Step " + itoa(i);
			Tau_N_Indiv[j][i] = new TH2F(name.c_str(), TitleName.c_str(),bin_number, x, bin_number, x);
			Tau_N_Indiv[j][i]->SetTitle(TitleName.c_str());
			Tau_N_Indiv[j][i]->SetLineColor(CurveColor(j));
			Tau_N_Indiv[j][i]->SetMarkerColor(CurveColor(j));
			Tau_N_Indiv[j][i]->SetMarkerStyle(10);
			Tau_N_Indiv[j][i]->SetLineColor(CurveColor(j));
			Tau_N_Indiv[j][i]->SetMarkerColor(CurveColor(j));

			Tau_N_Indiv[j][i]->SetXTitle("#tau_{N_{i}}");
			Tau_N_Indiv[j][i]->SetYTitle("#frac{N_{i}}{N_{tot}}");
			Tau_N_Indiv[j][i]->GetXaxis()->CenterTitle();
			Tau_N_Indiv[j][i]->GetYaxis()->CenterTitle();
			Tau_N_Indiv[j][i]->GetYaxis()->SetTitleOffset(1.25);
		}
		for (int j =0; j < 10; j++)
		{
			name = "TauN_" + SpeciesName[j] + "_Step_" + itoa(i);
			TitleName = "#tau_{N} " + SpeciesName[j] + " Step " + itoa(i);
			Tau_N_Indiv_Species[j][i] = new TH2F(name.c_str(), TitleName.c_str(),bin_number, x, bin_number, x);
			Tau_N_Indiv_Species[j][i]->SetTitle(TitleName.c_str());
			Tau_N_Indiv_Species[j][i]->SetLineColor(CurveColor(j+9));
			Tau_N_Indiv_Species[j][i]->SetMarkerColor(CurveColor(j+9));
			Tau_N_Indiv_Species[j][i]->SetMarkerStyle(10);
			Tau_N_Indiv_Species[j][i]->SetLineColor(CurveColor(j+9));
			Tau_N_Indiv_Species[j][i]->SetMarkerColor(CurveColor(j+9));

			Tau_N_Indiv_Species[j][i]->SetXTitle("#tau_{N_{i}}");
			Tau_N_Indiv_Species[j][i]->SetYTitle("#frac{N_{i}}{N_{tot}}");
			Tau_N_Indiv_Species[j][i]->GetXaxis()->CenterTitle();
			Tau_N_Indiv_Species[j][i]->GetYaxis()->CenterTitle();
			Tau_N_Indiv_Species[j][i]->GetYaxis()->SetTitleOffset(1.25);
		}
		name = "TauN_Full_Step_" + itoa(i);
		TitleName = "#tau_{N} Step " + itoa(i);

		Tau_N_Full[i] = new TH2F(name.c_str(), name.c_str(),bin_number, x, bin_number, x);
		Tau_N_Full[i]->SetTitle(TitleName.c_str());
		Tau_N_Full[i]->SetLineColor(8+1);
		Tau_N_Full[i]->SetMarkerColor(8+1);
		Tau_N_Full[i]->SetMarkerStyle(10);
		Tau_N_Full[i]->SetLineColor(8+1);
		Tau_N_Full[i]->SetMarkerColor(8+1);

		Tau_N_Full[i]->SetXTitle("#tau_{N}");
		Tau_N_Full[i]->SetYTitle("#frac{N}{N_{tot}}");
		Tau_N_Full[i]->GetXaxis()->CenterTitle();
		Tau_N_Full[i]->GetYaxis()->CenterTitle();
		Tau_N_Full[i]->GetYaxis()->SetTitleOffset(1.25);
	}

	for( it = evolutiondataBase.begin(); it != evolutiondataBase.end(); it++ )
	{

		map<ZAI ,TGraph* > QuantityRef =  (*it).second.GetEvolutionData();
		map<ZAI ,TGraph* > XSFissionRef =  (*it).second.GetFissionXS();
		map<ZAI ,TGraph* > XSCaptureRef =  (*it).second.GetCaptureXS();
		map<ZAI ,TGraph* > XSn2nRef =  (*it).second.Getn2nXS();



		map<IsotopicVector ,EvolutionData > evolutiondataBase_new = evolutiondataBase;
		map<IsotopicVector ,EvolutionData >::iterator it_bis = evolutiondataBase_new.begin();
		for (int k = 0; k < loop; k++)
			it_bis++;

		evolutiondataBase_new.erase(  it_bis   );
		MOX_new.SetDataBank(evolutiondataBase_new);

		MOX_new.LoadFPYield("SpontaneousFPyield.dat","FPyield_0.025.dat");

		MOX_new.SetFissionEnergy(90, 232, 0, 1.962e8);

		MOX_new.SetFissionEnergy(92, 233, 0, 1.99e8);

		MOX_new.SetFissionEnergy(92, 235,0, 201920000);
		MOX_new.SetFissionEnergy(92, 238, 0, 205520000);
		MOX_new.SetFissionEnergy(94, 239, 0, 209990000);
		MOX_new.SetFissionEnergy(94, 241, 0, 213600000);
		MOX_new.SetShortestHalfLife(3600);


		MOX.LoadFPYield("SpontaneousFPyield.dat","FPyield_0.025.dat");

		MOX.SetFissionEnergy(90, 232, 0, 1.962e8);

		MOX.SetFissionEnergy(92, 233, 0, 1.99e8);

		MOX.SetFissionEnergy(92, 235,0, 201920000);
		MOX.SetFissionEnergy(92, 238, 0, 205520000);
		MOX.SetFissionEnergy(94, 239, 0, 209990000);
		MOX.SetFissionEnergy(94, 241, 0, 213600000);
		MOX.SetShortestHalfLife(3600);

		EvolutionData EVOBuild = MOX_new.GenerateEvolutionData( (*it).second.GetIsotopicVectorAt(0.).GetActinidesComposition() ,it->second.GetFinalTime() , it->second.GetPower());


		double Na = 6.02214129e23;	//N Avogadro
		double M_ref = 0;
		double M = 0;
		double Power_ref =  (*it).second.GetPower() *(1+9.72932285999612498e-03);
		{
			map<ZAI, double >::iterator it2 ;


			IsotopicVector IVtmp = (*it).second.GetIsotopicVectorAt(0.).GetActinidesComposition().GetActinidesComposition() + EVOBuild.GetIsotopicVectorAt(0.).GetActinidesComposition();
			map<ZAI, double >isotopicquantity = IVtmp.GetIsotopicQuantity();

			for( it2 = isotopicquantity.begin(); it2 != isotopicquantity.end(); it2++ )
			{
				M_ref += (*it).second.GetIsotopicVectorAt(0.).GetActinidesComposition().GetActinidesComposition().GetZAIIsotopicQuantity( (*it2).first )*cZAIMass.fZAIMass.find( (*it2).first )->second/Na*1e-6;
				M += EVOBuild.GetIsotopicVectorAt(0.).GetActinidesComposition().GetZAIIsotopicQuantity( (*it2).first )*cZAIMass.fZAIMass.find( (*it2).first )->second/Na*1e-6;
			}
		}


		map<ZAI ,TGraph* > QuantityBuild =  EVOBuild.GetEvolutionData();
		map<ZAI ,TGraph* > XSFissionBuild =  EVOBuild.GetFissionXS();
		map<ZAI ,TGraph* > XSCaptureBuild =  EVOBuild.GetCaptureXS();
		map<ZAI ,TGraph* > XSn2nBuild =  EVOBuild.Getn2nXS();


		double CycleTime = it->second.GetFinalTime();
		double TauN[10] ;

		double TauPhi[10];
		double NiTot[10] ;
		double PhiTot[10];
		for (int i= 0; i < 10; i++)
			TauN[i] = TauPhi[i] = NiTot[i] = PhiTot[i] = 0;


		vector<double> N_tot[10];
		vector<double> DN_tot[10];

		double N[9][10];
		vector<double> N_species[10][10];
		vector<double> DN_species[10][10];
		
		double DN[9][10];
		for (int i= 0; i < 10; i++)
			for( int j =0; j < 9; j++)
				N[j][i] = DN[j][i] = 0;


		map<ZAI ,TGraph* >::iterator it_Q_Ref;

		for(it_Q_Ref = QuantityRef.begin(); it_Q_Ref != QuantityRef.end(); it_Q_Ref++)
		{

			map<ZAI ,TGraph* >::iterator it_Fiss_Ref = XSFissionRef.find( (*it_Q_Ref).first );
			map<ZAI ,TGraph* >::iterator it_Cap_Ref = XSCaptureRef.find( (*it_Q_Ref).first );
			map<ZAI ,TGraph* >::iterator it_n2n_Ref = XSn2nRef.find( (*it_Q_Ref).first );


			map<ZAI ,TGraph* >::iterator it_Q_Build = QuantityBuild.find( (*it_Q_Ref).first );
			map<ZAI ,TGraph* >::iterator it_Fiss_Build = XSFissionBuild.find( (*it_Q_Ref).first );
			map<ZAI ,TGraph* >::iterator it_Cap_Build = XSCaptureBuild.find( (*it_Q_Ref).first );
			map<ZAI ,TGraph* >::iterator it_n2n_Build = XSn2nBuild.find( (*it_Q_Ref).first );

			if( it_Fiss_Ref != XSFissionRef.end() && it_n2n_Ref != XSn2nRef.end() && it_Cap_Ref != XSCaptureRef.end()
			   && it_Q_Build != QuantityBuild.end() && it_Fiss_Build != XSFissionBuild.end() && it_Cap_Build != XSCaptureBuild.end() && it_n2n_Build != XSn2nBuild.end() )
				for (int i= 0; i < 10; i++)
				{
					N_tot[i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
					DN_tot[i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));


					if((*it_Q_Ref).first.A() == 238 && (*it_Q_Ref).first.Z() == 92)
					{
						N[0][i] = it_Q_Ref->second->Eval(CycleTime/10*i);
						DN[0][i] = abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i);
					}


					// U238
					if((*it_Q_Ref).first.A() == 238 && (*it_Q_Ref).first.Z() == 92)
					{
						N[0][i] = it_Q_Ref->second->Eval(CycleTime/10*i);
						DN[0][i] = abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i);
					}

					// U235
					else if((*it_Q_Ref).first.A() == 235 && (*it_Q_Ref).first.Z() == 92)
					{
						N[1][i] = it_Q_Ref->second->Eval(CycleTime/10*i);
						DN[1][i] = abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i);

					}
					else if((*it_Q_Ref).first.A() == 238 && (*it_Q_Ref).first.Z() == 94)
					{
						N[2][i] = it_Q_Ref->second->Eval(CycleTime/10*i);
						DN[2][i] = abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i);
					}

					// Pu239
					else if((*it_Q_Ref).first.A() == 239 && (*it_Q_Ref).first.Z() == 94)
					{
						N[3][i] = it_Q_Ref->second->Eval(CycleTime/10*i);
						DN[3][i] = abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i);
					}

					// Pu240
					else if((*it_Q_Ref).first.A() == 240 && (*it_Q_Ref).first.Z() == 94)
					{
						N[4][i] = it_Q_Ref->second->Eval(CycleTime/10*i);
						DN[4][i] = abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i);
					}

					// Pu241
					else if((*it_Q_Ref).first.A() == 241 && (*it_Q_Ref).first.Z() == 94)
					{
						N[5][i] = it_Q_Ref->second->Eval(CycleTime/10*i);
						DN[5][i] = abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i);
					}
					// Pu242
					else if((*it_Q_Ref).first.A() == 242 && (*it_Q_Ref).first.Z() == 94)
					{
						N[6][i] = it_Q_Ref->second->Eval(CycleTime/10*i);
						DN[6][i] = abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i);
					}

					// Am241
					else if((*it_Q_Ref).first.A() == 241 && (*it_Q_Ref).first.Z() == 95)
					{
						N[7][i] = it_Q_Ref->second->Eval(CycleTime/10*i);
						DN[7][i] = abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i);
					}
					//NP237
					else if((*it_Q_Ref).first.A() == 237 && (*it_Q_Ref).first.Z() == 93)
					{
						N[8][i] = it_Q_Ref->second->Eval(CycleTime/10*i);
						DN[8][i] = abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i);
					}


					//U
					if((*it_Q_Ref).first.Z() == 92)
					{
						N_species[0][i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
						DN_species[0][i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));
					}
					//Pu
					else if((*it_Q_Ref).first.Z() == 94)
					{
						N_species[1][i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
						DN_species[1][i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));
					}
					//Am
					else if((*it_Q_Ref).first.Z() == 95)
					{
						N_species[2][i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
						DN_species[2][i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));
					}
					//Np
					else if((*it_Q_Ref).first.Z() == 93)
					{
						N_species[3][i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
						DN_species[3][i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));
					}
					//Cm
					else if((*it_Q_Ref).first.Z() == 96)
					{
						N_species[4][i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
						DN_species[4][i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));
					}
					//Th
					else if((*it_Q_Ref).first.Z() == 90)
					{
						N_species[5][i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
						DN_species[5][i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));
					}
					//Pa
					else if((*it_Q_Ref).first.Z() == 91)
					{
						N_species[6][i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
						DN_species[6][i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));
					}
					//Actinides
					if((*it_Q_Ref).first.Z() >= 90)
					{
						N_species[7][i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
						DN_species[7][i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));
					}
					//PF
					else if((*it_Q_Ref).first.Z() > 20 && (*it_Q_Ref).first.Z() < 75)
					{
						N_species[8][i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
						DN_species[8][i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));
					}
					//Other
					else
					{
						N_species[9][i].push_back(it_Q_Ref->second->Eval(CycleTime/10*i));
						DN_species[9][i].push_back(abs( it_Q_Ref->second->Eval(CycleTime/10*i) - it_Q_Build->second->Eval(CycleTime/10*i) )/it_Q_Ref->second->Eval(CycleTime/10*i));
					}

					

					NiTot[i] += it_Q_Ref->second->Eval(CycleTime/10*i);


				}



		}

		for(int i=0; i < 10; i++)
		{

			for (int j =0; j < 9; j++)
				Tau_N_Indiv[j][i]->Fill(DN[j][i], N[j][i]/NiTot[i] );


			for (int j =0; j < (int)N_tot[i].size(); j++)
				Tau_N_Full[i]->Fill(DN_tot[i][j],N_tot[i][j]/NiTot[i] );


			for (int j =0; j < 10; j++)
				for (int k =0; k < (int)DN_species[j][i].size(); k++)
					Tau_N_Indiv_Species[j][i]->Fill(DN_species[j][i][k], N_species[j][i][k]/NiTot[i] );



		}

		loop++;
		//		if (loop == 30) break;
		int now  = time(NULL);
		int duree=difftime(now, start);
		int reste=(100.*duree/((float)loop*100.0/(float) evolutiondataBase.size()))-duree;
		if(loop%10 == 0)
			cout << "Step : " << loop << " || Still need : "  << reste/60 << " min " << reste%60 << " sec                                                 \r" << flush;

		MOX_new.Clear();
		EVOBuild.DeleteEvolutionData();

	}

	double T[10];
	for(int i=0; i < 10; i++)
	{

		for (int j =0; j < 9; j++)
			Tau_N_Indiv[j][i]->Write();
		for (int j =0; j < 10; j++)
			Tau_N_Indiv_Species[j][i]->Write();


		Tau_N_Full[i]->Write();
	}

	outfile->Close();


}
//==========================================================================================
// Compilation
//==========================================================================================
/*
 
 \rm CLASS* ; g++ -Wno-deprecated -o CLASS_exec TauCalculation.cpp -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp
 
 */
