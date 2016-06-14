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

	# Set *debug* to *True* to trace:
	debug = False
	#debug = True
	if debug:
	    print("=>Base.construct(*)")

	# Use *tube* instead of *self*:
	tube = self

	# Perform the extrusion:
	zero = L()
	thickness = L(inch="1/8")
	dx = L(inch=1.0)
	dy = dx

	x1 = -dx/2
	x2 = x1 + thickness
	x3 = -x2
	x4 = -x1	

	y1 = x1
	y2 = x2
	y3 = x3
	y4 = x4

	r = zero

	outer_contour = Contour("outer tube contour")
	outer_contour.bend_append("outer SW", P(x1, y1, zero), r)
	outer_contour.bend_append("outer NW", P(x1, y4, zero), r)
	outer_contour.bend_append("outer NE", P(x4, y4, zero), r)
	outer_contour.bend_append("outer SE", P(x4, y1, zero), r)

	inner_contour = Contour("inner tube contour")
	inner_contour.bend_append("inner SW", P(x2, y2, zero), r)
	inner_contour.bend_append("inner NW", P(x2, y3, zero), r)
	inner_contour.bend_append("inner NE", P(x3, y3, zero), r)
	inner_contour.bend_append("inner SE", P(x3, y2, zero), r)
	
	contours = [outer_contour, inner_contour]

	xy = L(inch=2.0)
	start = P(-xy,  xy, zero)
	end =   P( xy, -xy, zero)
	degrees0 = Angle()
	material = Material("aluminum", "")
	color = Color("green")

	tube.extrude("rectangular tube", material, color, contours, start, end, degrees0)
	  
if __name__== "__main__":
    ezcad = EZCAD3(0)
    tube= Tube(None, "extrusion")
    tube.process(ezcad)
