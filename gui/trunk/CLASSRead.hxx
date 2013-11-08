#ifndef _CLASSRead
#define _CLASSRead_

#include <TROOT.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>

#include <vector>
#include <iostream>

#include "ZAI.hxx"
#include "CLASSPlotElement.hxx"

using namespace std;

	//________________________________________________________________________
	//
	//		CLASSRead
	//@{
	//	CLASSRead is a class that allows to retrieve results in Binary or Ascii Data
	// 	files written by CLASS. By default binary files are read.
	//
	// @author FMS
	// @version 0.1
	//
	//      Modifications: @author RC
	//       - The module XSGUI: has been developped in the lattice code DRAGON to format the result in the same format as ASCII file  DATA_
	//         The consequences in CLASSGui are basicly introduction of new nuclei corresponding to chemical element (same Z different A).
	//          artificial Z,A number has been assigned to this nuclei:
	//                           Z_nat=Z+200
	//                           A_nat=A_moy+400 (where A_avg correspond to the most abondant isotope)
	//@}
	//________________________________________________________________________
class CLASSRead
{
public :
	
	CLASSRead(TString filemname);
	~CLASSRead();		//@- destructor
	
	vector< vector< TString > >	GetReactorName()	{return fReactorName;}
	vector< vector< TString > >	GetPoolName()		{return fPoolName;}
	vector< vector< TString > >	GetFabricationName()	{return fFabricationName;}
	vector< vector< TString > >	GetStockName()		{return fStockName;}
	vector< ZAI >		GetZAIvector()		{return fZAIvector;}
	vector< TTree* >	GetData()		{return fData;}
	
	string GetBranchInName(CLASSPlotElement toplot);
	string GetLegendOutName(CLASSPlotElement toplot);
	string GetTittleOutName(CLASSPlotElement toplot);

	
	void ReadName();
	void ReadZAI();
	void PlotTTree(vector<CLASSPlotElement> toplot, string opt = "L");
	void Plot(vector<CLASSPlotElement> toplot, string opt = "L");
	void PlotPower(vector<CLASSPlotElement> toplot, string opt = "L");
	void PlotTTreePower(vector<CLASSPlotElement> toplot, string opt = "L");
	void AddFile(TString filemname);
	TCanvas* fCNuclei;
	TCanvas* fCPower;
	
	
private :
		
	
	vector< vector< TString > > fReactorName;
	vector< vector< TString > > fPoolName;
	vector< vector< TString > > fFabricationName;
	vector< vector< TString > > fStockName;
	vector< ZAI >	fZAIvector;
	
	
	vector< TTree* > fData;
	TFile*	fFileIn;

	TGraph** fGraph;
	TLatex** fLegend;
	int fNumberGraphIterator;

	TGraph**  fGraphPower;
	TLatex** fLegendPower;
	int fNumberGraphPowerIterator;
	
};



#endif
