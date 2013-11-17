#!/usr/bin/env python

from EZCAD3 import *   # The EZCAD (revision 3) classes:

class Simple_Box(Part):

    def __init__(self, up, dx=L(mm=100.0), dy=L(mm=50.0),
      dz=L(25.0), wall_thickness=L(mm=5.0),
      material=Material("plastic", "ABS")):
	# Initialize the *Part*:
	Part.__init__(self, up)

	# Remember the initialization values:
	self.dx_l = dx
	self.dy_l = dy
	self.dz_l = dz
	self.wall_thickness_l = wall_thickness
	self.material_m = material

	# Instantiate the sub-*Part*'s:
	self.base_ = Simple_Box_Base(self)
	self.cover_ = Simple_Box_Cover(self)

    def construct(self):
	pass

class Simple_Box_Base(Part):

    def __init__(self, up, place = True):
	Part.__init__(self, up, place)

    def construct(self):
	# Grab some values from *box*:
	box = self.up
	dx = box.dx_l
	dy = box.dy_l
	dz = box.dz_l
	wall_thickness = box.wall_thickness_l
	material = box.material_m

	# Add another member variable:
	self.height_l = height = dz - wall_thickness
	zero = L()

	# Start with a solid block of the right dimensions:
	height = dz - wall_thickness
	self.block(comment = "Initial block of material",
	  material = material,
	  color = Color("blue"),
	  corner1 = P(-dx/2, -dy/2, zero),
	  corner2 = P( dx/2,  dy/2, height))

	# Pocket out the body of the box:
	self.simple_pocket(comment = "Box Pocket",
	  corner1 = self.bsw + \
            P(wall_thickness, wall_thickness, wall_thickness),
	  corner2 = self.tne - \
            P(wall_thickness, wall_thickness, zero),
          pocket_top = "t")

class Simple_Box_Cover(Part):

    def __init__(self, up, place = True):
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
