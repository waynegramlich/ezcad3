#!/usr/bin/env python

from EZCAD3 import *   # The EZCAD (revision 3) classes:

class Simple_Box(Part):

    def __init__(self, up):
	assert isinstance(up, Part) or up == None
	Part.__init__(self, self)
	self.base_ = Simple_Box_Base()
	self.cover_ = Simple_Box_Cover()

    def construct(self):
	pass

class Simple_Box_Base(Part):

    def __init__(self, up):
	assert isinstance(up, Part) or up == None
	assert isinstance(up, Part) or up == None
	Part.__init__(self, self)

    def construct(self):
	# Grab some values from *box*:
	box = self.up
	self.box.dx_l = dx = L(mm=100.0)
	self.box.dy_l = dy = L(mm=100.0)
	self.box.dz_l = dz = L(mm=25.00)
	self.wall_thickness_l = wall_thickness = L(mm=10.0)
	self.material_m = Material("plastic")
	self.color_c = Color("green")

	# Start with a solid block of the right dimensions:
	zero = L()
	extra = L(mm=6.0)
	self.extra(extra, extra, zero)
	self.block(comment = "Initial block of material",
	  corner1 = P(-dx/2, -dy/2, zero),
	  corner2 = P( dx/2,  dy/2, -dz),
	  material = material,
	  color = color,
	  top = "t")

	# Mount up the block:
	self.vice_position("Mount Block", self.t, self.tn, self.tw)
	self.tooling_plate("Tooling Plate", "2r 2c")
	self.tooling_plate_mount("Tooling Plate Mount")
	self.cnc_flush()

	# Pocket out the body of the box:
	bottom_corner = P(-dx/2 + thickness, -dy/2 + thickness, dz+L(mm=3.0))
	top_corner =    P( dx/2 - thickness,  dy/2 - thickness, zero)
	
	self.simple_pocket(comment = "Box Pocket",
	  bottom_corner1 = bottom_corner,
	  top_corner = top_corner,
	  radius = L(inch="1/8"),
	  pocket_top ="t")

class Simple_Box_Cover(Part):

    def __init__(self, up):
	assert isinstance(up, Part) or up == None
	Part.__init__(self, up, place)

    def construct(self):
	# Grab some values from *parent* and *base*:
	box = self.up
	dx = box.dx_l
	dy = box.dy_l
	dz = box.dz_l
	material = box.material_m
	wall_thickness = box.wall_thickness_l
	base = box.base_
	base_height = base.height_l
	zero = L()

	# Compute local values:
	self.lip_thickness_l = lip_thickness = wall_thickness / 2
	self.gap_l = gap = L(mm = 0.1)

	# Do the top part of the cover:
	self.block(comment = "Cover Top",
	  material = material,
	  color = Color("green"),
	  corner1 = base.tsw,
	  corner2 = base.tne + P(z = wall_thickness))

	# Do the lip part of the cover:
	self.block(comment = "Cover Lip",
	  corner1 = base.tsw + \
	  P(wall_thickness + gap, wall_thickness + gap, -lip_thickness),
	  corner2 = base.tne - \
	  P(wall_thickness + gap, wall_thickness + gap, zero),
          welds = "t")

if __name__== "__main__":
    ezcad = EZCAD3(0)               # Using EZCAD 3.0
    simple_box = Simple_Box(None)   # Initialize top-level sub-assembly
    simple_box.process(ezcad)       # Process the design
