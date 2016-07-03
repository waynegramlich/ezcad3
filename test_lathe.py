#!/usr/bin/env python

from EZCAD3 import *

class Lathe(Part):

    def __init__(self, up, name):
	""" *Lathe*: Initialize *Base* object (i.e. *self*) with a
	    parent of *up*.
	"""

	# Verify argument types:
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)

	# Initialize the *Part* superclass:
	Part.__init__(self, up, name)


    def construct(self):
	""" *Lathe*: Construct the *Base* object (i.e. *self*).
	"""

	# Use *lathe* instead of *self*:
	lathe = self

	# Define the *material* and *color*:
	material = Material("plastic", "acrylic")
	color = Color("green")

	# Define some constants:
	zero = L()
	one = L(inch=1.0)
	two = L(inch=2.0)

	# Define some points:
	start = P(zero, two, zero)
	peak =  P(one, one, zero)
	end =   P(zero, zero, zero)

	# Construct the exterior contour:
	r1 = L(inch="1/32")
	contour = Contour("profile")
	contour.bend_append("start", start, r1)
	contour.bend_append("peak",  peak,  r1)
	contour.bend_append("end",   end,   r1)

	# Lathe it out:
	lathe.tool_prefer("Laser_007")
	lathe.lathe("lathe", material, color, start, end, contour, 32)

def test_lathe():
    ezcad = EZCAD3(0, directory="/tmp")
    lathe= Lathe(None, "lathe")
    lathe.process(ezcad)

if __name__== "__main__":
    test_lathe()
