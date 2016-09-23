#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import os

import html.htmlHelper as htmlH
import html.model
import html.backend

import cgitb; cgitb.enable(display=1,logdir=".",format="text")
import cgi

import json


def strHeader ( models ) :
	xsmInc = set([ "#include \"XS/XSM_MLP.hxx\"" for m in models if m["type"] == "Equivalence_Model" ])
	eqmInc = set([ "#include \"Equivalence/"+m["em"]+".hxx\"" for m in models if m["type"] == "Equivalence_Model" ])
	imInc  = set([ "#include \"Irradiation/"+m["im"]["model"]+".hxx\"" for m in models if m["type"] == "Equivalence_Model" ])

	s = ""
	s += "\n".join(xsmInc)
	s += "\n".join(eqmInc)
	s += "\n".join(imInc)

	return s


def strLogger ( log ) :
	"""
		Returns as a string the definition of a logger in CLASS
	"""
	s = """
	//Definition of the Log file : CLASS messages output
	int Std_output_level  = %(termLevel)s;
	int File_output_level = %(fileLevel)s;
	CLASSLogger *Logger   = new CLASSLogger("%(file)s",Std_output_level,File_output_level);"""

	data = {\
		'termLevel' : log["termLevel"],\
		'fileLevel' : log["fileLevel"],\
		'file'      : log["file"]\
	}

	return s % data

def strScenario ( config ) :
	"""
		Returns as a string the configuration of the gCLASS scenario
	"""
	s = """
	Scenario *gCLASS=new Scenario( %(timeStart)s*year , Logger );
	gCLASS->SetStockManagement( %(SetStockManagement)s );
	gCLASS->SetTimeStep( %(timeStep)s*year );
	gCLASS->SetOutputFileName( "%(outputFile)s" );
	gCLASS->SetZAIThreshold( %(SetZAIThreshold)s );
	gCLASS->SetSoberTerminalOutput();
	"""

	data = {\
		'timeStart'          : config["timeStart"],\
		'SetStockManagement' : config["SetStockManagement"],\
		'timeStep'           : config["timeStep"],\
		'outputFile'         : config["outputFile"],\
		'SetZAIThreshold'    : config["SetZAIThreshold"]\
	}

	return s % data


def strStock ( unity ) :
	"""
		Returns as a string the code of one stock
	"""

	initStock = [ "\t\t{id}->AddToStock( ZAI({z},{a},{i}) * {q} );".format(id=unity['id'],z=zaiQ[0],a=zaiQ[1],i=zaiQ[2],q=zaiQ[3]) for zaiQ in json.loads("["+unity["iv"]+"]") ]
	
	data = {\
		'id'   : unity["id"],\
		'init' : "\n".join( initStock ),\
	}

	s = """
	Storage * %(id)s = new Storage ( gCLASS->GetLog() );
		%(id)s->SetName( "%(id)s" );
%(init)s
	gCLASS->Add( %(id)s );"""

	return s % data

def strStocks ( backend ) :
	"""
		Returns as a string the code of all stocks (loop on `strStock()`)
	"""
	stocks = [ strStock(unity) for unity in backend if unity["type"] == "Stock" ]

	return "\n".join(stocks)

def strPool ( unity ) :
	"""
		Returns as a string the code for one pool
	"""
	data = {\
		'id'              : unity["id"],\
		'tcooling'        : unity["tcooling"],\
		'backendFacility' : unity["backendFacility"],\
	}

	s = """
	Pool * %(id)s = new Pool( gCLASS->GetLog() , %(backendFacility)s , %(tcooling)s*year );
		%(id)s->SetName( "%(id)s" );
	gCLASS->Add( %(id)s );"""

	return s % data

def strSPAddZAIEff ( sp , ivEff ) :
	"""
		Returns as a string the code to add an IsotopicVector to a SP
	"""
	ivId = ivEff['id']
	adds = [ "{id}.Add( {z},{a},{i} , {e} );".format(id=ivId,z=zaiE[0],a=zaiE[1],i=zaiE[2],e=zaiE[3]) for zaiE in json.loads("["+ivEff['vector']+"]") ]

	data = {\
		'sp'      : sp,\
		'ivId'    : ivId,\
		'add'     : "\n\t\t\t".join(adds),\
		'backend' : ivEff['out'],\
	}

	s = """{
			IsotpicVector %(ivId)s;
			%(add)s
			%(sp)s->SetBackEndDestination( %(backend)s , %(ivId)s , 0*year );
		}"""

	return s % data



def strSP ( unity ) :
	"""
		Returns as a string the code for one separation plant
	"""
	ivSep = [ strSPAddZAIEff( unity['id'] , ivEff ) for ivEff in unity["ivlist"] ]

	data = {\
		'id'     : unity["id"],\
		'ivList' : " ".join(ivSep),\
	}

	s = """
	SeprationPlant * %(id)s = new SeprationPlant( gCLASS->GetLog() );
		%(ivList)s
	gCLASS->AddSeparationPlant( %(id)s );"""

	return s % data


def strBackendFac ( backend ) :
	"""
		Returns as a string the code for all separation plants and pools
	"""
	separationPlants = [ strSP(unity) for unity in backend if unity["type"] == "SeparationPlant" ]
	pools = [ strPool(unity) for unity in backend if unity["type"] == "Pool" ]

	s  = "\n".join(separationPlants)
	s += "\n"
	s += "\n".join(pools)
	return s


def strModels ( models ) :
	return " "


def strReactors ( models ) :
	return " "



form = cgi.FieldStorage()

print( """%(contentType)s

	<!doctype html>
	<html>
	<head>
		<meta charset="utf-8" />
		<title>Run</title>
	</head>
	<body>
	Run for your life !
	</body>
	</html>
""" % { 'contentType': "Content-type: text/html;" })

f = open("./form.txt","w")
f.write( "plop\n" )

f.write( str(form) )
f.write( "\n----\n" )

scenario = json.loads( str( form.getvalue("data") ) )
f.write( str(scenario) )

config  = scenario["config"]
backend = scenario["backend"]
models  = scenario["models"]

f.write( "\n---- config \n"  )
f.write( str(config) )

f.write( "\n---- backend \n" )
f.write( str(backend) )

f.write( "\n---- models \n"  )
f.write( str(models) )

f.write( "\n" )

data = { \
	'modelsHeader': strHeader( models ),\
	'logger'      : strLogger( config["log"] ),\
	'scenario'    : strScenario( config ),\
	'stocks'      : strStocks( backend ),\
	'backendFac'  : strBackendFac( backend ),\
	'models'      : strModels( models ),\
	'reactors'    : strReactors( models ),\
	'timeEnd'     : config["timeEnd"],\
	'inputFile'   : config["inputFile"],\
}

inputCLASS = """
#include "CLASSHeaders.hxx"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>

%(modelsHeader)s

using namespace std;

int main(int argc, char** argv)
{
	cSecond year = 3600*24.*365.25; //< seconds in one year


/* ----- LOG MANAGEMENT --------------------------------------------------- */
%(logger)s

/* ----- CONFIGURATION PATH ----------------------------------------------- */
	string CLASS_PATH = getenv("CLASS_PATH");
	if ( CLASS_PATH == "" )
	{
		cout << "Please setenv CLASS_PATH to your CLASS installation folder in your .bashs or .tcshrc" << endl;
		exit(1);
	}
	string PATH_TO_DATA = CLASS_PATH + "/DATA_BASES/";

/* ----- Define the Scenario ---------------------------------------------- */
%(scenario)s

// Decay data base -----------------------------------------------------------
	// The decay data base is taken from the file Decay.idx
	// This decay data base will be used for all the decay calculations in this Scenario
	DecayDataBank* DecayDB = new DecayDataBank( gCLASS->GetLog() , PATH_TO_DATA + "DECAY/ALL/Decay.idx" );
	gCLASS->SetDecayDataBase(DecayDB);

// Backend Facilities --------------------------------------------------------
	// Stock -----------------------------------------------------------------
	%(stocks)s

	// SeparationPlant/Pool --------------------------------------------------
	%(backendFac)s


// Reactor data base ---------------------------------------------------------
	%(models)s

// Reactors ------------------------------------------------------------------
	%(reactors)s

/* ----- Start simulation ------------------------------------------------- */
	gCLASS->Evolution( (double)year*%(timeEnd)s );

	delete gCLASS;
	return 0;
}

//============================================================================
// Compilation
//============================================================================
/*
 rm CLASS* ; g++ -o CLASS_Exec %(inputFile)s -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -lTMVA -fopenmp -lgomp -Wunused-result
*/

"""

f.write( " ---- " )
f.write( inputCLASS % data )

f.close()
