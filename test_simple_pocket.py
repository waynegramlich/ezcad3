#!/usr/bin/env python

from EZCAD3 import *

class Base(Part):

    def __init__(self, up, name):
	""" *Base*: Initialize *Base* object (i.e. *self*) with a parent of *up*. """

	# Verify argument types:
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)

	# Initialize the superclass:
	Part.__init__(self, up, name)


    def construct(self):
	""" *Base*: Construct the *Base* object (i.e. *self*). """

	# Set *debug* to *True* to trace:
	debug = False
	debug = True
	if debug:
	    print("=>Base.construct(*)")

	# Use *base* instead of *self*:
	base = self

	# Grab some values from *box*:
	base.dx_l = dx = L(inch=4.0)
	base.dy_l = dy = L(inch=4.0)
	base.dz_l = dz = L(inch=0.25)
	base.thickness_l = thickness = L(inch=.25)
	base.wall_thickness_l = wall_thickness = L(inch=.25)
	base.material_m = material = Material("plastic")
	base.color_c = color = Color("green")
	print("dx={0:i} dy={1:i} dz={2:i}".format(dx, dy, dz))

	# Start with a solid block of the right dimensions:
	zero = L()
	extra = L(mm=10.0)
	degrees0 = Angle(deg=0.0)
	base.extra_xyz(extra, extra, zero)
	corner1 = P(-dx/2, -dy/2, zero)
	corner2 = P( dx/2,  dy/2, -dz)
	base.block("Initial block of material", material, color, corner1, corner2, "t", "")

	# Mount up the block:
	base.vice_position("Mount Block", base.t, base.tn, base.tw)
	base.tooling_plate("Tooling Plate", "2r 2c")
	base.tooling_plate_mount("Tooling Plate Mount")
	base.cnc_fence()

	# Pocket out the body of the box:
	corner1 = P(-dx/2 + thickness, -dy/2 + thickness, -dz+L(inch="1/8"))
	corner2 = P( dx/2 - thickness,  dy/2 - thickness, zero)
	radius = L(inch=1.0)
	pocket_t = "t"
	flags = ""
	print("pocket bottom_corner={0:i} top_corner={1:i}".format(corner1, corner2))
	base.simple_pocket("Box Pocket", corner1, corner2, radius, pocket_t, degrees0, flags)

	if debug:
	    print("<=Base.construct(*)")

def test_simple_pocket():
    ezcad = EZCAD3(0, directory="/tmp")
    base = Base(None, "Base")
    base.process(ezcad)
    
def test_bounding_box():
    Bounding_Box.unit_test()

if __name__== "__main__":
    test_simple_pocket()
