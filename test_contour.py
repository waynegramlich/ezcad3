#!/usr/bin/env python

from EZCAD3 import *

class Plate(Part):

    def __init__(self, up):
	""" *Plate*: Initialize *Base* object (i.e. *self*) with a
	    parent of *up*.
	"""

	# Verify argument types:
	assert up == None or isinstance(up, Part)

	# Initialize the superclass:
	Part.__init__(self, up)


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
	plate.dz_l = dz = L(inch="1/8")
	plate.color_c = color = Color("blue")
	#print("dx={0:i} dy={1:i} dz={2:i}".format(dx, dy, dz))

	# Start with a solid block of the right dimensions:
	zero = L()


	# Define the extra material to be add to the next block operation:
	extra = L(inch="1/8")
	plate.extra_xyz(extra, extra, zero)

	# Compute the *plate* boundaries:
	x1 = -dx/2
	y1 = -dy/2
	z1 = -dz
	x2 =  dx/2
	y2 =  dy/2
	z2 =  zero

	# Compute the two corners and the *material*:
	corner1 = P(x1, y1, z1)
	corner2 = P(x2, y2, z2)
	material = Material("plastic", "acrylic")
	plate.block(comment = "Initial block of material for plate",
	  corner1 = corner1,
	  corner2 = corner2,
	  material = material,
	  color = color,
	  top = "t")

	# Mount up the block:
	plate.vice_position("Mount Block", plate.t, plate.tn, plate.tw)
	plate.tooling_plate("Tooling Plate", "2r 2c")
	plate.tooling_plate_mount("Tooling Plate Mount")

	# Construct the exterior contour:
	r1 = L(inch="1/4")
	contour = Contour("exterior_contour")

	clockwise = True
	if clockwise:
	    contour.bend_append("SW corner",     P(x1,   y1,   z2), r1)
	    contour.bend_append("NW corner",     P(x1,   y2,   z2), r1)
	    contour.bend_append("NE corner",     P(x2,   y2,   z2), r1)
	    contour.bend_append("E corner",      P(x2,   zero, z2), r1)
	    contour.bend_append("O corner", P(zero, zero, z2), r1)
	    contour.bend_append("S corner",      P(zero, y1,   z2), r1)
	else:
	    contour.bend_append("SW corner",     P(x1,   y1,   z2), r1)
	    contour.bend_append("S corner",      P(zero, y1,   z2), r1)
	    contour.bend_append("O corner", P(zero, zero, z2), r1)
	    contour.bend_append("E corner",      P(x2,   zero, z2), r1)
	    contour.bend_append("NE corner",     P(x2,   y2,   z2), r1)
	    contour.bend_append("NW corner",     P(x1,   y2,   z2), r1)

	plate.contour("exterior contour", contour, P(x1, y1, z2), P(x1, y1, z1), extra/2, "")

	if debug:
	    print("<=Base.construct(*)")

if __name__== "__main__":
    ezcad = EZCAD3(0)
    plate= Plate(None)
    plate.process(ezcad)
