#ifndef _CLASSWIN_
#define _CLASSWIN_




#include <time.h>
#include <vector>
#include <iostream>


#include <TROOT.h>
#include <TStyle.h>
#include <TGFrame.h>
#include "TGLayout.h"
#include <TGTab.h>
#include <TGButton.h>
#include <TGMsgBox.h>
#include <TGComboBox.h>
#include <TGTextBuffer.h>
#include <TGLabel.h>
#include <TGTextEntry.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLegend.h>
#include <TH2.h>
#include "TColor.h"
#include "THStack.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TGClient.h>
#include <TGTextView.h>
#include <TGSlider.h>
#include <TImage.h>

#include "CLASSRead.hxx"

using namespace std; 
//________________________________________________________________________
//
//		Windows and Widgets classes.
//@{
//		MainWin is the main window with Tab (1 for inventory, 1 for cross-section
//		and 1 for fluxes). SubWin is used for small windows that appears when
//		some buttons are pressed (e.g. Save As button, Exec Macro button, ...)
// @author PTO
// @version 1.0
// @modified PG,RC,BL
//@}

//Here is the widget code: it correspond to signal emitted when a widget is activated 
enum CommandId{
	M_BUTTON_PLOT=20,
	M_BUTTON_SAVE,
	M_BUTTON_MACRO,
	M_BUTTON_QUIT,
	M_BUT_OK,
	M_BUT_CANCEL,
	M_CHECK_PLOTALL};
			  
			  
class MainWin :  public  TGMainFrame   
{
public:

	MainWin(CLASSRead *DATA); 
	MainWin(){}
	~MainWin();		//@- destructor
	
	
	void Start();
	bool ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2); //@- widget signal handler method
	void CloseWindow();			//@- destroy the main window
	void Plot();				//@- general Plot method


	void FillItemTab(int i);
	
	void FillNucTab();			//@- fill the Inventory Tab foil
	void FillTotalTab();
	void FillReactorTab();			//@- fill the Cross-section Tab foil
	void FillStockTab();			//@- fill the Flux Tab foil
	void FillPoolTab();			//@- fill the Reaction Rate Tab foil
	void FillFabricationTab();		//@- fill the Breeding Ratio Tab foil
	void FillNucFoil(int NCheck,int Ncol,int Nline);

	string fSaveFileName;
	
private:
	
	
	CLASSRead* fDATA;

	
	
	int		fGraphLineWidth;
	double		fGraphMarkerSize;

	int	fNumberOfParc;
	int	fNumberOfTOT;
	int *	fNumberOfReactor;
	int *	fNumberOfStock;
	int *	fNumberOfPool;
	int *	fNumberOfFab;

	vector<TGCheckButton*>	fSaveCheckButton;		//@- save each of these check boxes
	
	//
	//widgets
	//
	FontStruct_t fLabelFontB;
	FontStruct_t fLabelFontS;
	
	//main tab
	TGCompositeFrame *	fGeneF0;		//@- Common widget frame
	TGCompositeFrame *	fPlotSaveQuitFrame;
	
	TGCompositeFrame **	fParcTabFoil;
	TGTab *			fParcTab; //Parc table
	
	TGCompositeFrame ***	fFacilitiesTabFoil;
	
	TGTab **		fFacilitiesTab;// tab
	TGTab ***		fItemTab;
	
	TGCheckButton ***	fCheckArrayTotal;
	TGCheckButton ***	fCheckArrayReactor;
	TGCheckButton ***	fCheckArrayStock;
	TGCheckButton ***	fCheckArrayPool;
	TGCheckButton ***	fCheckArrayFab;
	
	TGCompositeFrame ***	fTabFoilReactor;
	TGCompositeFrame ***	fTabFoilStock;
	TGCompositeFrame ***	fTabFoilPool;
	TGCompositeFrame ***	fTabFoilFab;
	
	
	
	TGTextButton *	fButtonPlot;		//@- Press buttons to Plot,Save,...
	TGTextButton *	fButtonSave;
	TGTextButton *	fButtonQuit;
	int		fMainWidth;		//@- Width of fGeneF0
	
	//Inventory tab
	TGTab *			fTabNuc;		//@- Tab of Inventory nuclei
	TGCompositeFrame **	fTabFoilNuc;		//@- its foils: AM,FP0,...,FPN, gaz, misc
	TGComboBox *		fTimeScaleNuc;		//@- time scale combo box
	
	
	//@- specific frames
	TGCompositeFrame *	fMiscGrpFNuc;
	TGCompositeFrame *	fXGprFNuc;
	TGCompositeFrame *	fYGprFNuc;
	TGCompositeFrame *	fMiscHzFrame1;
	TGCompositeFrame *	fMiscHzFrame2;
	TGCompositeFrame *	fAxisHzFrame;
	
	TGCheckButton **	fCheckArrayNuc;	//@- array of check box for inventory nuclei tab
	int			fNselectedNucleus;	//@- number of check box selected in the previous array

	
	
};




class SubWin: public  TGTransientFrame
{
public :
	SubWin(string event,const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
               UInt_t options = kVerticalFrame);   	//@- Normal constructor
	~SubWin();
	
	
	Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2); //@- widget signal handler method
	void CloseWindow(){delete this;}		//@- Close windows

	
private:

	TGLayoutHints *		fL0;			//@- position hints in the sub win
	TGCompositeFrame *	fS0;			//@- main frame of the sub win
	TGHorizontalFrame *	fSH1;			//@- misc frame
	TGHorizontalFrame *	fSH2;			//@- misc frame
	
	TGTextButton *		fButtonOK;		//@- OK an Cancel button
	TGTextButton *		fButtonCan;		//@- OK an Cancel button
	
	TGLabel	*		fLmsg;			//@- a label for the write text widget
	TGLabel	*		fNucLabel[11];
	TGTextView *		fTGVNuclei;
	TGTextEntry *		TEName;		//@- the write text widget in save as and exec macro
 	TGTextEntry *		TEExtract[11];		//@- the extract nuclei
	
	MainWin	*		fParent; 		//@- parent of the sub win (i.e. MainWin)
	   
  };


#endif