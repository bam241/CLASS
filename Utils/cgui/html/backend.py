#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import htmlHelper as htmlH

backendType = ['Pool','SeparationPlant','Stock']

emptyBackendFacility = """{id:'',type:'',ivlist:[],iv:'',tcooling:5}"""

def strBackend ( unity ) :
	"""
		return the formulaire of an backend facility as a string
	"""
	return htmlH.strHTMLfile( "html/backend/" + unity + ".html" )

def backendList () :
	"""
		browse backends and stringify the switch and formulaire for each one
	"""
	r = ""
	for unity in backendType :
		r += """<!-- %(unity)s -->
			<div class="panel-body" ng-switch-when="%(unity)s" >
				%(form)s
			</div>
			""" % {'unity':unity,'form':strBackend(unity)}

	return r

def form () :
	backend = htmlH.strHTMLfile("html/backend/backend.html").format( backendList() )
	return backend
