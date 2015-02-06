//
// This program train and test a MLP using TMVA from a training data in a form of a TTRee
//
//@author Root_tmva_team modified by BaL 
#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>
#include <sstream>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVARegGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#endif

using namespace TMVA;
   
int main(int argc, char const *argv[])
{
   if(argc!=1)
   {
      std::cout << "TrainMLP Usage :" << std::endl;
      std::cout << "\tTrainMLP YourTrainingData.root" << std::endl;
   }   
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegression" << std::endl;


   // Create a new root output file
   TString outfileName( "TMVA_MOX_Equivalence.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );//UPDATE


   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/ 
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in 
   // front of the "Silent" argument in the option string

   TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile, 
                                               "!V:!Silent:Color:DrawProgressBar" );

   // If you wish to modify default settings 
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the TMVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   factory->AddVariable( "BU"           , "Burn Up"     , "GWj_tMLi"        , 'F' );
   factory->AddVariable( "U5_enrichment","U5 enrichment", "FractionIsotopic", 'F' );
   factory->AddVariable( "Pu8"          , "Pu 238"      , "FractionIsotopic", 'F' );
   factory->AddVariable( "Pu9"          , "Pu 239"      , "FractionIsotopic", 'F' );
   factory->AddVariable( "Pu10"         , "Pu 240"      , "FractionIsotopic", 'F' );
   factory->AddVariable( "Pu11"         , "Pu 241"      , "FractionIsotopic", 'F' );
   factory->AddVariable( "Pu12"         , "Pu 242"      , "FractionIsotopic", 'F' );
   factory->AddVariable( "Am1"          , "Am 241"      , "FractionIsotopic", 'F' );


   // You can add so-called "Spectator variables", which are not used in the MVA training, 
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the 
   // input variables, the response values of all trained MVAs, and the spectator variables
   //factory->AddSpectator( "spec1:=var1*2",  "Spectator 1", "units", 'F' );
   //factory->AddSpectator( "spec2:=var1*3",  "Spectator 2", "units", 'F' );

   // Add the variable carrying the regression target
   factory->AddTarget("teneur");

   // It is also possible to declare additional targets for multi-dimensional regression, ie:
   // -- factory->AddTarget( "fvalue2" );
   // BUT: this is currently ONLY implemented for MLP

   // Read training and test data (see TMVAClassification for reading ASCII files)
   // load the signal and background event samples from ROOT trees
   TFile *input(0);
   TString fname = argv[1]; //Training Data input file
   if (!gSystem->AccessPathName( fname )) 
      input = TFile::Open( fname ); // check if file in local directory exists 
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVARegression           : Using input file: " << input->GetName() << std::endl;

   // --- Register the regression tree

   TTree *regTree = (TTree*)input->Get("Data");

   // global event weights per tree (see below for setting event-wise weights)
   Double_t regWeight  = 1.0;   

   // You can add an arbitrary number of regression trees
   factory->AddRegressionTree( regTree, regWeight );

   // This would set individual event weights (the variables defined in the 
   // expression need to exist in the original TTree)
  // factory->SetWeightExpression( "var1", "Regression" );

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycut = ""; // for example: TCut mycut = "abs(var1)<0.5 && abs(var2-0.5)<1";

   // tell the factory to use all remaining events in the trees after training for testing:
   factory->PrepareTrainingAndTestTree( mycut, 
                                        "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:

   // ---- Book MVA methods
   //
   std::stringstream Name;
   Name<<"MLP_Equivalence"<<0;
   // Neural network (MLP)                                                                                 
   factory->BookMethod( TMVA::Types::kMLP, Name.str().c_str(), "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+10:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );


   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();    

   // --------------------------------------------------------------
   
   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;      

   delete factory;

   // Launch the GUI for the root macros
if (!gROOT->IsBatch()) TMVARegGui( outfileName );
}
/*

g++ -o Train_MLP `root-config --cflags` Train_MLP.cxx `root-config --glibs` -lTMVA -I$ROOTSYS/tmva/test/
*/