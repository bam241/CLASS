/***********************************/
// Train one MLP "INDICE" from the
//file ../BuildInput/TrainingInput.root
//
//@author Root_tmva_Team modified by BaL
/***********************************/

#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
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

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#endif

using namespace TMVA;
std::vector<std::string>OUTPUT;
void LOAD_OUTPUT() 
{

   #include  "../BuildInput/TrainingInput.cxx"

}
void Train_XS_Time(int INDICE) 
{

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegression" << std::endl;

// --------------------------------------------------------------------------------------------------
// --- Here the preparation phase begins

   // Create a new root    OUTPUT file
   std::stringstream OutTMVA;
   OutTMVA<<"Training_output_"<<INDICE<<".root";
   TString outfileName( OutTMVA.str().c_str() );
   TFile*    OUTPUTFile = TFile::Open( outfileName, "UPDATE" );//RECREATE

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory will
   // then run the performance analysis for you.
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/ 
   //
   // The second argument is the    OUTPUT file for the training results
   // All TMVA    OUTPUT can be suppressed by removing the "!" (not) in 
   // front of the "Silent" argument in the option string


   TMVA::Factory *factory = new TMVA::Factory( "TMVARegression",    OUTPUTFile, 
                                               "!V:!Silent:Color:DrawProgressBar" );

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@Change : change the  each first arguments with the name used in your training sample
 //   factory->AddVariable( "TheNameOfYourNucleusInTheTrainingFile"           , "whatever"    , "whatever", 'F' ); 
   factory->AddVariable( "U5"           , "U 235"    , "FractionIsotopic", 'F' ); 
   factory->AddVariable( "U8"           , "U 238"    , "FractionIsotopic", 'F' ); 
   factory->AddVariable( "Pu8"          , "Pu 238"   , "FractionIsotopic", 'F' );
   factory->AddVariable( "Pu9"          , "Pu 239"   , "FractionIsotopic", 'F' );
   factory->AddVariable( "Pu10"         , "Pu 240"   , "FractionIsotopic", 'F' );
   factory->AddVariable( "Pu11"         , "Pu 241"   , "FractionIsotopic", 'F' );
   factory->AddVariable( "Pu12"         , "Pu 242"   , "FractionIsotopic", 'F' );
   factory->AddVariable( "Am1"          , "Am 241"   , "FractionIsotopic", 'F' );
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   factory->AddVariable( "Time"         , "Time"     , "seconds"         , 'F' );

   // Add the variable carrying the regression target
   factory->AddTarget(   OUTPUT[INDICE].c_str() ); //The name of the MLP output

   // It is also possible to declare additional targets for multi-dimensional regression, ie:
   // -- factory->AddTarget( "fvalue2" );
   // BUT: this is currently ONLY implemented for MLP

   // Read training and test data (see TMVAClassification for reading ASCII files)
   // load the signal and background event samples from ROOT trees
   TFile *input(0);

   TString fname = "../BuildInput/TrainingInput.root";
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

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycut = ""; // for example: TCut mycut = "abs(var1)<0.5 && abs(var2-0.5)<1";

   // tell the factory to use all remaining events in the trees after training for testing:
   factory->PrepareTrainingAndTestTree( mycut, 
                                        "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // If no numbers of events are given, half of the events in the tree are used 
   // for training, and the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );  

   // ---- Book MVA methods
   //
   // please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   std::stringstream Name;
   Name<<   OUTPUT[INDICE];
   // Neural network (MLP)                                                                                    
      factory->BookMethod( TMVA::Types::kMLP, Name.str().c_str(), "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N,N:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();    

   // --------------------------------------------------------------
   
   // Save the    OUTPUT
      OUTPUTFile->Close();

   std::cout << "==> Wrote root file: " <<    OUTPUTFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;      

   delete factory;

}
int  main(int argc, char const *argv[])
{
   LOAD_OUTPUT();
   if(argc != 2 || argv[1] == "-h")
     {
      std::cout << "Usage : TrainXS i"<< std::endl;
      std::cout << "With i the cross section index " << std::endl;
      std::cout << "File ../BuildInput/TrainingInput.cxx indicates\n indice ranging from 0 to "<< OUTPUT.size()-1 << std::endl;
      exit(0);
     } 

   Train_XS_Time(atoi(argv[1])) ;

   return 0;
}
 /*
 COMPILATION :

 g++ -o Train_XS  `root-config --cflags` Train_XS.cxx `root-config --glibs` -lTMVA
 
 */
