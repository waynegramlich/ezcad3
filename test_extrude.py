#!/usr/bin/env python

from EZCAD3 import *

class Tube(Part):

    def __init__(self, up, name):
	""" *Tube*: Initialize *Base* object (i.e. *self*) with a parent of *up*. """

	# Verify argument types:
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)

	# Initialize the superclass:
	Part.__init__(self, up, name)


    def construct(self):
	""" *Tube*: Construct the *Base* object (i.e. *self*). """

	# Some constants:
	zero = L()
	degrees0 = Angle()

	# Perform the extrusion:
	xy = L(inch=2.0)
	height = L(inch=1.0)
	width = height
	start = P( -xy, zero, zero)
	end =   P(zero, zero, zero)
	material = Material("aluminum", "")
	color = Color("green")
	thickness = L(inch="1/8")
	
	self.rectangular_tube_extrude("rectangular tube", material, color,
	  width, height, thickness, start, zero, end, zero, degrees0)
	  
if __name__== "__main__":
    ezcad = EZCAD3(0)
    tube= Tube(None, "extrusion")
    tube.process(ezcad)
