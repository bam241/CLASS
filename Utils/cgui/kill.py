#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Kill Python Web server with the pid
"""

import os

# display a beautiful page to notify the server will be killed soon
import html.htmlHelper as htmlH

data = { \
'contentType': htmlH.contentType("text/html") ,\
'header'     : htmlH.header("Kill") ,\
'footer'     : htmlH.footer() \
}

print ("""%(contentType)s

%(header)s

<div class="shining" style="height:100%%; width:100%%" >
	<div class="row">
		<div class="col-sm-6" style="margin: 25%% 2em; background:rgba(255,255,255,0.7); padding: 1em; border-radius:1em;" >
			<h4>Server killed ;(</h4>
		</div>
	</div>
</div>

%(footer)s
""" % data )

# read PID file
pid_file = open("server.pid","r")
pid = pid_file.readline()
pid_file.close()

# kill this pid process
os.kill(int(pid),5)
os.remove("server.pid")


# R.I.P. Python Web Server ;(
# ...
# I will be back !
