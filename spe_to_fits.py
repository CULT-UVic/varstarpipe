#!/usr/bin/python

import sys
import piUtils
import glob
import pyfits

#author: ctooley@uvic.ca
#questions? use google or email me
#this is a simple script meant to ease the use of piUtils.py - you need to have
#piUtils.py in the same dir as this script.

if len( sys.argv ) <= 1 :
	print "Usage: " + sys.argv[0] + " [filename.spe [filename.spe [...]]]"
	print "   OR"
	print "       " + sys.argv[0] + " *.spe"
	sys.exit(0)


for arg in sys.argv[1:] :
	for speFile in glob.glob( './' + arg ) :
		name = speFile.replace( '.spe', '' )
		name = speFile.replace( '.SPE', '' )
		print speFile
		piUtils.speToFits(speFile, fitsfile=name, clobber=True)
		
