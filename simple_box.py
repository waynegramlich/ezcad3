#!/usr/bin/env python

from EZCAD3 import *   # The EZCAD (revision 3) classes:

class Simple_Box(Part):

    def __init__(self, up, name):
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)
	Part.__init__(self, up, name)
	self.base_ = Simple_Box_Base(self, "Simple_Box_Base")
	self.cover_ = Simple_Box_Cover(self, "Simple_Box_Cover")
	self.dx_l = L(mm=100.0)
	self.dy_l = L(mm=100.0)
	self.dz_l = L(mm=25.00)
	self.wall_thickness_l = L(mm=10.0)
	self.material_m = Material("plastic", "ABS")

    def construct(self):
	pass

class Simple_Box_Base(Part):

    def __init__(self, up, name):
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)
	Part.__init__(self, up, name)

    def construct(self):
	# Grab some values from *box*:
	box = self.up
	dx = box.dx_l
	dy = box.dy_l
	dz = box.dz_l
	wall_thickness = box.wall_thickness_l
	material = box.material_m

	# Start with a solid block of the right dimensions:
	zero = L()
	color = Color("green")
	#extra = L(mm=6.0)
	#self.extra(extra, extra, zero)
	self.block("Initial block of material",	material, color,
	  P(-dx/2, -dy/2, zero), P( dx/2,  dy/2, -dz), "t", "")

	# Mount up the block:
	self.vice_position("Mount Block", self.t, self.tn, self.tw)
	self.tooling_plate("Tooling Plate", "2r 2c")
	self.tooling_plate_mount("Tooling Plate Mount")
	#self.cnc_flush()

	# Pocket out the body of the box:
	bottom_corner = P(-dx/2 + wall_thickness, -dy/2 + wall_thickness, -dz + wall_thickness)
	top_corner =    P( dx/2 - wall_thickness,  dy/2 - wall_thickness, zero)
	self.simple_pocket("Box Pocket", bottom_corner, top_corner, L(inch="1/8"), "t", Angle(), "")

class Simple_Box_Cover(Part):

    def __init__(self, up, name):
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)
	Part.__init__(self, up, name)

    def construct(self):
	# Grab some values from *parent* and *base*:
	box = self.up
	dx = box.dx_l
	dy = box.dy_l
	dz = box.dz_l
	material = box.material_m
	wall_thickness = box.wall_thickness_l
	base = box.base_
	base_height = base.t.z - base.b.z
	color = Color("blue")
	zero = L()

	# Compute local values:
	self.lip_thickness_l = lip_thickness = wall_thickness / 2
	self.gap_l = gap = L(mm = 0.1)

	# Do the top part of the cover:
	self.block("Cover Top", material, color,
           base.tsw, base.tne + P(zero, zero, wall_thickness), "t", "")

	# Do the lip part of the cover:
	self.block("Cover Lip", material, color,
	  base.tsw + P(wall_thickness + gap, wall_thickness + gap, -lip_thickness),
	  base.tne - P(wall_thickness + gap, wall_thickness + gap, zero), "t", "t")

if __name__== "__main__":
    ezcad = EZCAD3(0)               # Using EZCAD 3.0
    simple_box = Simple_Box(None, "Simple_Box")   # Initialize top-level sub-assembly
    simple_box.process(ezcad)       # Process the design
