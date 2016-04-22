#ifndef _BINARYFORMAT2_
#define _BINARYFORMAT2_

#include <fstream>
using namespace std; 

/*!
 \file
 \brief Header file for binary output files management.
 
 This file is to be included in any units which manipulate (read/void write) MURE
 binary output files. These structures allow easy way to handle 
 binary output of MURE evolution data (BDATA_* files).
 
 @author JW
 @author PTO
 @version 1.0
*/

//! Header of MURE output binary file
struct FileHeader 
{
	FileHeader() {Version=1;}	//!< Constructor
	short Version;				//!< Binary gile version
	float Time;					//!< Time of MCNP step at which this file is printed
	float K;					//!< keff
	float Kerr;					//!< keff error
	int NCells;					//!< Number of evolving cells
	void write(ofstream &out);	//!< Write the header into a stream
	void read(ifstream &in);	//!< Read the header from a stream
};

//! Header of an evolving cell in a binary file
struct CellHeader 
{
	int CellNumber;		//!< MCNP number of this cell
	float Volume;		//!< Volume of the cell
	float Flux;			//!< Flux in this cell
	float FluxErr;		//!< Flux error 
	int NNucleusRecords;		//!< Number of nuclei records in this cell
	void write(ofstream &out);	//!< Write the header into a stream
	void read(ifstream &in);	//!< Read the header from a stream
};

//! Record of a nucleus in a binary file	
struct NucleusRecord 
{
	short Z;	//!< Proton number of the nucleus
	short A;	//!< Nucleon number of the nucleus
	short I;	//!< Isomeric state of the nucleus
	float Mass; //!< Atomic mass
	float Proportion;			//!< Number of this nuclei in cell
	short NReactionRecords;		//!< Number of reaction records of this nucleus
	void write(ofstream &out);	//!< Write the record into a stream
	void read(ifstream &in);	//!< Read the record from a stream
};

//! Record of a reaction in a binary file
struct ReactionRecord 
{
	short Code;		//!< Reaction code
	float Sigma;	//!< Reaction cross-section
	float SigmaErr;	//!< Cross-section error
	void write(ofstream &out);	//!< Write the record into a stream
	void read(ifstream &in);	//!< Read the record from a stream
};
//--------------------------------------------------------
//	Implementation of methods
//-------------------------------------------------------
void FileHeader::write(ofstream &out)
{
	out.write((char*)&Version,sizeof(Version));
	out.write((char*)&Time,sizeof(Time));
	out.write((char*)&K,sizeof(K));
	out.write((char*)&Kerr,sizeof(Kerr));
	out.write((char*)&NCells,sizeof(NCells));
}

void FileHeader::read(ifstream &in)
{
	in.read((char*)&Version,sizeof(Version));
	in.read((char*)&Time,sizeof(Time));
	in.read((char*)&K,sizeof(K));
	in.read((char*)&Kerr,sizeof(Kerr));
	in.read((char*)&NCells,sizeof(NCells));
}
	
void CellHeader::write(ofstream &out)
{
	out.write((char*)&CellNumber,sizeof(CellNumber));
	out.write((char*)&Volume,sizeof(Volume));
	out.write((char*)&Flux,sizeof(Flux));
	out.write((char*)&FluxErr,sizeof(FluxErr));
	out.write((char*)&NNucleusRecords,sizeof(NNucleusRecords));
}

void CellHeader::read(ifstream &in)
{
	in.read((char*)&CellNumber,sizeof(CellNumber));
	in.read((char*)&Volume,sizeof(Volume));
	in.read((char*)&Flux,sizeof(Flux));
	in.read((char*)&FluxErr,sizeof(FluxErr));
	in.read((char*)&NNucleusRecords,sizeof(NNucleusRecords));
}
	
void NucleusRecord::write(ofstream &out)
{
	out.write((char*)&Z,sizeof(Z));
	out.write((char*)&A,sizeof(A));
	out.write((char*)&I,sizeof(I));
	out.write((char*)&Mass,sizeof(Mass));
	out.write((char*)&Proportion,sizeof(Proportion));
	out.write((char*)&NReactionRecords,sizeof(NReactionRecords));
}

void NucleusRecord::read(ifstream &in)
{
	in.read((char*)&Z,sizeof(Z));
	in.read((char*)&A,sizeof(A));
	in.read((char*)&I,sizeof(I));
	in.read((char*)&Mass,sizeof(Mass));
	in.read((char*)&Proportion,sizeof(Proportion));
	in.read((char*)&NReactionRecords,sizeof(NReactionRecords));
}
	
void ReactionRecord::write(ofstream &out)
{
	out.write((char*)&Code,sizeof(Code));
	out.write((char*)&Sigma,sizeof(Sigma));
	out.write((char*)&SigmaErr,sizeof(SigmaErr));
}

void ReactionRecord::read(ifstream &in)
{
	in.read((char*)&Code,sizeof(Code));
	in.read((char*)&Sigma,sizeof(Sigma));
	in.read((char*)&SigmaErr,sizeof(SigmaErr));
}
	

#endif
