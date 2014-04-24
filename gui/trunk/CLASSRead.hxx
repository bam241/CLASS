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
	// Thank to Francisco, Guerin, Louard, Shi for their participation to the XML ouput implementation
	// @author BLG
	// @author BAM
	// @version 1.0
	// @version 0.1
	//
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

	void ConvertxmlTTreeMass(vector<CLASSPlotElement> toplot, string filename);
	void ConvertXmlMass(vector<CLASSPlotElement> toplot, string opt);           
	void ConvertxmlTTreePower(vector<CLASSPlotElement> toplotPower, string opt);
	void ConvertXmlPower(vector<CLASSPlotElement> toplotPower, string opt );

	void Write(string filename, string fileformat="ASCII");
	void ASCIIWrite(string filename = "ASCII");


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
	TGraph* fGraphSumOfSelected;

	TLatex** fLegend;
	TLatex* fLegendSumOfSelected;
	int fNumberGraphIterator;

	TGraph**  fGraphPower;
	TLatex** fLegendPower;
	int fNumberGraphPowerIterator;


	double Xmin;
	double Xmax;
	double Ymin;
	double Ymax;

	
};



#endif
