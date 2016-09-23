#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Launcher of Python Web server on 8080 port
It does not need any permission on a lot of configuration. It can run on local
or on a distant server for computing.
"""

port=8080

import os , sys

python_version = sys.version_info.major

# import switch between Python's version
if python_version == 3 :
	import http.server
else :
	import BaseHTTPServer
	import CGIHTTPServer

import cgitb; cgitb.enable(format="html")

# save pid of server into a file (`server.pid`)
pid = str(os.getpid())
pid_file = open("server.pid","w")
pid_file.write( pid )
pid_file.close()

# prepare sever
server_address = ("", port)
cgi_dir = ["/"]

if python_version == 3 :
	handler = http.server.CGIHTTPRequestHandler
	handler.cgi_directories = cgi_dir
	httpd = http.server.HTTPServer(server_address, handler)
else :
	handler = CGIHTTPServer.CGIHTTPRequestHandler
	handler.cgi_directories = cgi_dir
	httpd = CGIHTTPServer.BaseHTTPServer.HTTPServer(server_address, handler)


# launch server
print("""Server HTTP on localhost -- port: """ + str(server_address[1]) + """\n\tpid: """ + pid)
try:
	httpd.serve_forever()
except KeyboardInterrupt:
	print("""Server killed""")
	os.remove("server.pid")
	httpd.socket.close()
