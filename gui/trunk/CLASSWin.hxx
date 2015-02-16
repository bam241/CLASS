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

// Thank to Francisco, Guerin, Louard, Shi for their participation to the XML ouput implementation
// @author BLG
// @author BAM
// @version 1.0

//@}

//Here is the widget code: it correspond to signal emitted when a widget is activated 
enum CommandId{
	M_BUTTON_PLOT=20,
	M_BUTTON_SAVE,
	M_BUTTON_MACRO,
	M_BUTTON_QUIT,
	M_BUT_OK,
	M_BUT_CANCEL,
	M_CHECK_PLOTALL,
	M_RADIO_XML_SAVE,
	M_RADIO_ASCII_SAVE,
	M_RADIO_DECAY_HEAT,
	M_RADIO_RADIOTOX,
	M_CB_SCENAR_Time,
	M_CHECK_LINEAR_Tox,
	M_TE_toxfirst,
	M_TE_toxlast,
	M_TE_toxnstep,
	M_CHECK_AM_NUC,
	M_CHECK_FP_NUC,
	M_CHECK_BY_MOTHER,
	M_CHECK_INSIDE,
	M_CHECK_CUMIN,
	M_CHECK_CUMOUT,
	M_BUTTON_MOTHER_MORE};
			  
			  
class MainWin :  public  TGMainFrame   
{
public:

	MainWin(CLASSRead *DATA,vector<string> VFileName); 
	MainWin(){}
	~MainWin();		//@- destructor
	
	
	void Start(vector<string> VFileName);
	bool ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2); //@- widget signal handler method
	void CloseWindow();			//@- destroy the main window
	void Plot();				//@- general Plot method
        void Conversionxml();                   //@- general Conversionxml method


	void FillItemTab(int i);
	
	void FillNucTab();			//@- fill the Inventory Tab foil
	void FillTotalTab();
	void FillReactorTab();			//@- fill the Cross-section Tab foil
	void FillStockTab();			//@- fill the Flux Tab foil
	void FillPoolTab();			//@- fill the Reaction Rate Tab foil
	void FillFabricationTab();		//@- fill the Breeding Ratio Tab foil
	void FillNucFoil(int NCheck,int Ncol,int Nline);

	string fSaveFileName;
	string fSaveFileFormat;
	
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
	
	//Configuration Tab

	TGTab *			fPlotConfigTab; //Plot configuratio tab
	TGCompositeFrame **	fPlotConfigFoil;
	TGCompositeFrame	*fInventoryFrame;
	
	
		//Radio Or Decay
	TGCompositeFrame *	fDecayOrRadioFrame;
	TGRadioButton		*fButtonRadiotox,*fButtonHeat;
	
		//By mother
	bool			fMotherIsVisible,fIsByMother,fIsLinear;
	TGPictureButton 	*fByMotherMore;  //@- hide or show by mother
	int			fTimeStep;
	TGCheckButton		*fByMotherButton ;
	TGComboBox		*fScenarTimeCBox;
	TGCompositeFrame	*fTimeParametersFrame,*fByMotherFrame,*fScenarTimeFrame;
	TGTextEntry		*TEtoxfirst,*TEtoxlast,*TEtoxnstep;// the write text widget in time evolution parameters
	TGCheckButton		*fCheckLinear;

	double			fToxTimeFirst;	// first time steps for radiotoxicity calculations
	double			fToxTimeLast;	// last time steps for radiotoxicity calculations
	int			fToxNstep;	// number of time steps for radiotoxicity calculations

	
	//misc frame
	TGCompositeFrame	*fMiscFrame,*fMiscHzFrame;
	TGCheckButton		*fCheckAMNuc,*fCheckFPNuc,*fCheckSumOfSelected;
	
	//Factories Arrays
	TGCheckButton ***	fCheckArrayTotal;
	TGCheckButton ***	fCheckArrayReactor;
	TGCheckButton ***	fCheckArrayStock;
	TGCheckButton ***	fCheckArrayPool;
	TGCheckButton ***	fCheckArrayFab;
	//Factories Foils
	TGCompositeFrame ***	fTabFoilReactor;
	TGCompositeFrame ***	fTabFoilStock;
	TGCompositeFrame ***	fTabFoilPool;
	TGCompositeFrame ***	fTabFoilFab;
	

	TGCheckButton ** fCheckIVPlot ;

	
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
	void SaveAs();

	
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

	TGRadioButton *fRadioASCIISave;
	TGRadioButton *fRadioXMLSave;
	   
  };


#endif
