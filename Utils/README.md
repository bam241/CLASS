# Utils folder
## EQM
This folder contains programs to help you with the generation of equivalence models (EQM). An EQM describes how the FabricationPlant build a fuel.

### FBR\_MLP\_Keff
#### Generate
The programs in this folder generate weigths files to be used with :  
**\$CLASS_PATH/Model/Equivalence/EQM\_FBR\_MLP\_Keff.cxx**  
This model is related to Fast Breeder Reactor. The idea is that the fissile content is such as the keff of the reactor is equal to a user defined value at time T (often Begining of Cycle, Middle or End Of Cycle). Where T depends on the weigths models. For instance, the weights in \$CLASS\_PATH/DATA\_BASES/FBR\_Na/MOX/EQModel/MLP\_K\_EFF\_BOC aims to predict the keff at begining of cycle for a FBR-Na loaded with MOX fuel. Thus, programs in this folder build Artificial Neural Network (ANN) weights files to estimate keff at time T according your depletion (or just neutron transport if T=BOC) calculations results.

##### Usage:
* Generate a folder containing results of many depletions calculations of FBR reactor (each calculation should defer in terms of fresh fuel composition). It has to be formated as **EvolutionData**s (see USEGUIDE.pdf in ../documentation for format). If you are using MURE or SMURE depletion code you can use the [MURE2CLASS](#mure2class) utility to convert MURE outputs in EvolutionData (.dat files)
* Compile the file .cxx with :

```bash
g++ -o Generate_FBR_Keff Generate_FBR_Keff.cxx `root-config --cflags` `root-config --libs`
```
* Execute with

```bash
./Generate_FBR_Keff PATH_TO_YOUR_EVOLUTION_DATAS
```
* Follow program instructions

#### Test
The programs in this folder test the regression performances of the model you created thanks to \$CLASS_PATH/Utils/EQM/FBR\_MLP\_Keff/Generate/Generate\_FBR.cxx.
##### Usage:


### MLP\_Kinf
#### Generate
The programs in this folder generate weigths files to be used with :  
**\$CLASS\_PATH/source/Model/Equivalence/EQM\_MLP\_Kinf.cxx**  
This model is related to non breeder reactor. See USEGUIDE.pdf in ../documentation for additional informations.

##### Usage:
* Generate a folder containing results of many depletions calculations of a non-breeder reactor (each calculation should defer in terms of fresh fuel composition). It has to be formated as **EvolutionData**s (see \$CLASS\_PATH/documentation/Manuel/USEGUIDE.pdf for format). If you are using MURE or SMURE depletion code you can use the [MURE2CLASS](#mure2class) utility to convert MURE outputs in EvolutionData (.dat files) 
* Compile the file .cxx with      
```
g++ -o Generate_MLP_Kinf Generate_MLP_Kinf.cxx `root-config --cflags` `root-config --libs`
```
* Execute with   
```
./Generate_MLP_Kinf PATH_TO_YOUR_EVOLUTION_DATAS
```
* Follow program instructions


#### Test
The programs in this folder test the regression performances of the model you created thanks to \$CLASS\_PATH/Utils/EQM/MLP\_Kinf/Generate/Generate\_MLP\_Kinf.cxx.
##### Usage:


### PWR\_MOX\_MLP
#### Generate
The programs in this folder generate weigths files to be used with :  
**\$CLASS\_PATH/source/Model/Equivalence/EQM\_PWR\_MLP\_MOX.cxx**   
This model is made for PWR loaded with MOX (U,Pu)O2 fuel. See USEGUIDE.pdf in ../documentation for additional informations.
#####Usage: 

* Generate a folder containing results of many depletions calculations of a non-breeder reactor (each calculation should defer in terms of fresh fuel composition). It has to be formated as **EvolutionData**s (see USEGUIDE.pdf for format). If you are using MURE or SMURE depletion code you can use the [MURE2CLASS](#mure2class) utility to convert MURE outputs in EvolutionData (.dat files)   

#### Test
The programs in this folder test the regression performances of the model you created thanks to \$CLASS\_PATH/Utils/EQM/PWR\_MOX\_MLP/Generate/Generate_PWR.cxx.
#####Usage:


## MURE2CLASS
Convert a (S)MURE evolution result to CLASS format EvolutionData (.dat file)
##### Usage:
* Compile the file .cxx with      
```  g++ -std=c++11 -o MURE2CLASS MURE2CLASS.cxx```

* Execute with   
```
./MURE2CLASS PATH_TO_YOUR_MURE_RESULT ReactorType FuelType
```
with ReactorType and FuelType is any strings you want.

## XSM / MLP
These utils allow you to build mean cross section predictors based on artificial neural networks (MLP). The weights and .nfo file produced by this utils are planned to be used with $CLASS\_PATH/source/Model/XSM/XSM\_MLP.cxx
### Generate 
#####Usage: 

* Generate a folder containing results of many depletions calculations (each calculation should defer in terms of fresh fuel composition). It has to be formated as **EvolutionData**s (see USEGUIDE.pdf for format). If you are using MURE or SMURE depletion code you can use the [MURE2CLASS](#mure2class) utility to convert MURE outputs in EvolutionData (.dat files) 
 
* Compile file Generate_XSM.cxx with : 
``` g++ -o Generate_XSM Generate_XSM.cxx `root-config --cflags` `root-config --libs` ``` 

* Execute with 
``` ./Generate_MLP_Kinf PATH_TO_YOUR_EVOLUTION_DATAS ```

* Follow instructions.


### Test
The program in this folder aims to test regression performance of your mean cross section predictors.
Once you have trained your ANN with Generate_MLP_Kinf, you can test their performances using indications in file EvaluateTrainingCommands.dat.

## cgui
Graphical user interface (web) to define a CLASS input. Usefull if you have poor knowledge in C++. Note that this user friendly interface does not have all the possibilities given by the C++ one. More informations can be found [here](cgui/README.md) (in french).
