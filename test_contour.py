#!/usr/bin/env python

from EZCAD3 import *

class Plate(Part):

    def __init__(self, up, name):
	""" *Plate*: Initialize *Base* object (i.e. *self*) with a
	    parent of *up*.
	"""

	# Verify argument types:
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)

	# Initialize the superclass:
	Part.__init__(self, up, name)


    def construct(self):
	""" *Plate*: Construct the *Base* object (i.e. *self*).
	"""

	# Set *debug* to *True* to trace:
	debug = False
	#debug = True
	if debug:
	    print("=>Base.construct(*)")

	# Use *plate* instead of *self*:
	plate = self

	# 
	plate.dx_l = dx = L(inch=4.0)
	plate.dy_l = dy = L(inch=4.0)
	plate.dz_l = dz = L(inch="1/2")
	plate.color_c = color = Color("blue")
	#print("dx={0:i} dy={1:i} dz={2:i}".format(dx, dy, dz))

	# Start with a solid block of the right dimensions:
	zero = L()

	# Define the extra material to be add to the next block operation:
	extra = L(inch="1/8")
	plate.extra_xyz(extra, extra, zero)

	# Compute the *plate* boundaries:
	x1 = -dx/2
	x2 = -dx/4
	x3 =  dx/4
	x4 =  dx/2

	y1 = -dy/2
	y2 = -dy/4
	y3 =  dy/4
	y4 =  dy/2

	z1 = -dz
	z2 =  zero

	# Compute the two corners and the *material*:
	corner1 = P(x1, y1, z1)
	corner2 = P(x4, y4, z2)
	material = Material("plastic", "acrylic")
	plate.block("Plate block", material, color, corner1, corner2, "t", "")

	# Construct the exterior contour:
	r1 = L(inch="2/16")
	contour = Contour("exterior_contour")

	clockwise = False
	if clockwise:
	    contour.bend_append("SW corner",     P(x1,   y1,   z2), r1)
	    contour.bend_append("NW corner",     P(x1,   y4,   z2), r1)
	    contour.bend_append("NE corner",     P(x4,   y4,   z2), r1)
	    contour.bend_append("E corner",      P(x4,   zero, z2), r1)
	    contour.bend_append("O corner", 	 P(zero, zero, z2), r1)
	    contour.bend_append("S corner",      P(zero, y1,   z2), r1)
	else:
	    contour.bend_append("SW corner",     P(x1,   y1,   z2), r1)
	    contour.bend_append("S corner",      P(zero, y1,   z2), r1)
	    contour.bend_append("O corner",      P(zero, zero, z2), r1)
	    contour.bend_append("E corner",      P(x4,   zero, z2), r1)
	    contour.bend_append("NE corner",     P(x4,   y4,   z2), r1)
	    contour.bend_append("NW corner",     P(x1,   y4,   z2), r1)

	plate.tool_prefer("Laser_007")
	plate.contour("exterior contour", contour, P(x1, y1, z2), P(x1, y1, z1), extra/2, "")

	# Drill some holes:
	no_36 = L(inch=0.1065)
	hole1_end =   P(x2, y2, z1)
	hole1_start = P(x2, y2, z2)
	plate.countersink_hole("hole1", no_36, no_36*2, hole1_start, hole1_end, "t")

	hole2_end =   P(x2, y3, z1)
	hole2_start = P(x2, y3, z2)
	plate.countersink_hole("hole2", no_36, no_36*2, hole2_start, hole2_end, "t")

	hole3_end =   P(x3, y3, z1)
	hole3_start = P(x3, y3, z2)
	plate.countersink_hole("hole3", no_36, no_36*2, hole3_start, hole3_end, "t")

	if debug:
	    print("<=Base.construct(*)")

def test_contour():
    ezcad = EZCAD3(0, directory="/tmp")
    plate= Plate(None, "plate")
    plate.process(ezcad)

if __name__== "__main__":
    test_contour()
