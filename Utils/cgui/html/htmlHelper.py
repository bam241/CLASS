#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import cgi
import os, sys
import fnmatch, glob, re

def contentType ( type ) :
	return """Content-type: """ + type + """; charset=utf-8"""

def printHTMLfile ( file ) :
	"""
		print the content of the HTML file
	"""
	f = open(file,'rb')
	print f.read()
	f.close()

def strHTMLfile ( file ) :
	"""
		return the content of HTML file as a string
	"""
	with open( file , 'rb' ) as f :
		data = f.read()
	return data

def header (title) :
	"""
		return as a string the content of the HTML header file with the correct title page
	"""
	return strHTMLfile("html/header.html")%{'title':title}

def footer () :
	"""
		return as a string the content of the HTML footer file
	"""
	return strHTMLfile("html/footer.html")

def findEQModels () :
	"""
		find xml models in $CLASS_PATH/**/EQModel
	"""
	EQModels = []
	for root, dirnames, filenames in os.walk(os.environ['CLASS_PATH']+"/DATA_BASES") :
		for dir in fnmatch.filter(dirnames, 'EQModel'):
			EQModels.append(os.path.join(root, dir))

	matches = []
	for path in EQModels :
		for root, dirnames, filenames in os.walk( path ) :
			for file in filenames :
				if fnmatch.fnmatch(file, '*.xml') :
					matches.append(os.path.join(root,file))

	return matches

def findXSModels () :
	"""
		find nfo models in $CLASS_PATH/**/EQModel
	"""
	XSModels = []
	for root, dirnames, filenames in os.walk(os.environ['CLASS_PATH']+"/DATA_BASES") :
		for dir in fnmatch.filter(dirnames, 'XSModel'):
			XSModels.append(os.path.join(root, dir))

	matches = []
	for path in XSModels :
		for root, dirnames, filenames in os.walk( path ) :
			for file in filenames :
				if fnmatch.fnmatch(file, '*.nfo') :
					matches.append(os.path.join(root,""))

	return matches

def findMLPWeight () :
	"""
		find all xml files in $CLASS_PATH/DATA_BASES in a dir EQModel
	"""
	# all dir `EQModel`
	EQModels = []
	for root, dirnames, filenames in os.walk(os.environ['CLASS_PATH']+"/DATA_BASES") :
		for dir in fnmatch.filter(dirnames, 'EQModel'):
			EQModels.append(os.path.join(root, dir))

	matches = []
	for path in EQModels :
		for root, dirnames, filenames in os.walk( path ) :
			for file in filenames :
				if fnmatch.fnmatch(file, '*.xml') :
					matches.append(os.path.join(root,file))

	return matches


def findFixedFuel () :
	"""
		find models in $CLASS_PATH/**/FixedFuel
	"""
	matches = []
	for root, dirnames, filenames in os.walk(os.environ['CLASS_PATH']+"/DATA_BASES") :
		for dir in fnmatch.filter(dirnames, 'FixedFuel'):
			matches.append(os.path.join(root, dir))

	return matches

def listEQModels () :
	"""
		list class in `source/Model/Equivalence`
	"""
	r = []

	for file in glob.glob( os.environ['CLASS_PATH']+"/source/Model/Equivalence/*.hxx" ) :
		for line in open(file) :
			if "public EquivalenceModel" in line :
				r.append(re.sub(":public.*$","",line.replace("class","").replace('\n','').replace(' ','')))

	return r


def listIModels () :
	"""
		list irradiation models
	"""
	return [ re.sub(".hxx","",os.path.basename(eqm)) for eqm in glob.glob( os.environ['CLASS_PATH']+"/source/Model/Irradiation/*.hxx" ) ]

def listXSModels () :
	"""
		list XS Model
	"""
	return [ re.sub(os.environ['CLASS_PATH']+"/DATA_BASES/","",xsm) for xsm in findXSModels() ]

def listMLPWeight ():
	"""
		list MLP Weight xml files with the correct path
	"""
	return [ re.sub(os.environ['CLASS_PATH']+"/DATA_BASES/","",mlp) for mlp in findMLPWeight() ]


def listEvoDataPATH () :
	"""
		return a list of EvoData with the absolute PATH
	"""
	return [ glob.glob( ed+"/*.dat") for ed in findFixedFuel() ][0]

def listEvoData () :
	"""
		return a list of EvoData with the relative PATH
	"""
	return [ re.sub(os.environ['CLASS_PATH']+"/DATA_BASES/","",ed) for ed in listEvoDataPATH() ]

def listStorageManagement () :
	"""
		return a list of different storage management
	"""
	return ["FiFo","LiFo","Mix","Rand"]


def labelAddUnity ( list , model , name ) :
	"""
		return HTML code for a label button to add unities
	"""
	return """<button class="label label-success btn btn-success btn-right" data-ng-click=\"""" + list +""".push( """ + model + """ )\" >Add """ + name + """</button> <span class="badge btn-right" >{{ """ + list + """.length }}</span>"""

def printSectionTitle ( title , label ) :
	print ("""&nbsp;
	<section class="panel panel-primary" >
		<header class="panel-heading panel-title" >
			""" + title + """
			""" + label + """
		</header>
	</section>
		""")

def sectionTitle ( title , label ) :
	return """&nbsp;
	<section class="panel panel-primary" >
		<header class="panel-heading panel-title" >
			%(title)s
			%(label)s
		</header>
	</section>
		"""%{'title':title , 'label':label}