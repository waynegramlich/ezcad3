#!/usr/bin/env python

from EZCAD3 import *

class Block(Part):

    def __init__(self, up, name):
	""" *Block*: Initialize *Base* object (i.e. *self*) with a
	    parent of *up*.
	"""

	# Verify argument types:
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)

	# Initialize the superclass:
	Part.__init__(self, up, name)


    def construct(self):
	""" *Block*: Construct the *Plate* object (i.e. *self*). """

	# Use *block* instead of *self*:
	block = self

	# Some constants:
	zero = L()

	# Specify *dx*, *dy*, and *dz*:
	block.dx_l = dx = L(inch=4.0)
	block.dy_l = dy = L(inch=4.0)
	block.dz_l = dz = L(inch="1/2")
	block.color_c = color = Color("blue")
	block.material_m = material = Material("aluminum", "")

	# Define the extra material to be add to the next block operation:
	extra = L(inch="1/8")
	block.extra_xyz(extra, extra, zero)

	# Compute the *block* boundaries:
	x1 = -dx
	x2 = x1 + dx/4
	x4 =  zero
	x3 = x4 - dx/4

	y1 = -dy
	y2 = y1 + dy/4
	y4 =  zero
	y3 = y4 - dy/4

	z1 = -dz
	z2 =  zero
	z3 =  dz

	# Do the bottom block:
	corner1 = P(x1, y1, z1)
	corner2 = P(x4, y4, z2)
	block.block("bottom block", material, color, corner1, corner2, "t", "t")

	# Weld n the top block:
	corner1 = P(x2, y2, z2)
	corner2 = P(x3, y3, z3)
	block.block("block", material, color, corner1, corner2, "t", "b")

def test_block():
    ezcad = EZCAD3(0, directory="/tmp")
    block = Block(None, "block")
    block.process(ezcad)


if __name__== "__main__":
    test_block()
