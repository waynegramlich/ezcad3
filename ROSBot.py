#!/usr/bin/env python

# Common metric screw lengths:
#
#     http://www.bjg-design.com/designbook/mscrlng1.htm
#
#
# Excerpt from McMaster Carr.
# Class 12.9 Black Alloy Steel (Metric)- DIN 912.
# Screw
# Length	Screw Size
# 	M1.6  M2   M3   M4   M5   M6   M8   M10  M12  M14  M16  M20  M24  M30
#    3	  x    -    -    -    -    -    -    -    -    -    -    -    -    -  
#    4	  x    -    -    -    -    -    -    -    -    -    -    -    -    -  
#    5	  x    x    -    -    -    -    -    -    -    -    -    -    -    -  
#    6	  x    x    x    x    -    -    -    -    -    -    -    -    -    -  
#    8	  -    x    x    x    -    -    -    -    -    -    -    -    -    -  
#   10	  -    -    x    x    x    x    x    -    -    -    -    -    -    -  
#   12	  -    x    x    x    x    x    -    -    -    -    -    -    -    -  
#   14	  -    -    x    x    x    x    -    -    -    -    -    -    -    -  
#   16	  -    -    x    x    x    x    x    x    -    -    -    -    -    -  
#   18	  -    -    -    x    x    x    x    -    -    -    -    -    -    -  
#   20	  -    -    x    x    x    x    x    x    -    -    -    -    -    -  
#   22	  -    -    -    -    x    x    x    -    -    -    -    -    -    -  
#   25	  -    -    x    x    x    x    x    x    x    -    -    -    -    -  
#   30	  -    -    x    x    x    x    x    x    x    x    x    -    -    -  
#   35	  -    -    x    x    x    x    x    x    x    x    x    -    -    -  
#   40	  -    -    -    x    x    x    x    x    x    x    x    -    -    -  
#   45	  -    -    -    -    -    x    x    x    x    -    x    -    -    -  
#   50	  -    -    -    -    x    x    x    x    x    x    x    x    -    -  
#   55	  -    -    -    -    x    x    x    x    x    -    -    -    -    -  
#   60	  -    -    -    -    x    x    x    x    x    x    x    x    x    -  
#   65	  -    -    -    -    -    x    x    x    x    x    x    x    x    -  
#   70	  -    -    -    -    -    x    x    x    x    -    -    x    x    -  
#   75	  -    -    -    -    -    -    x    -    x    -    -    -    -    -  
#   80	  -    -    -    -    -    -    x    x    x    x    -    x    x    x  
#   90	  -    -    -    -    -    -    x    x    x    x    -    x    x    x  
#  100	  -    -    -    -    -    -    x    x    x    x    -    x    x    x  
#  110	  -    -    -    -    -    -    x    x    x    x    -    x    x    x  
#  120	  -    -    -    -    -    -    x    x    x    x    -    x    x    x  
#  150	  -    -    -    -    -    -    x    x    x    x    -    x    x    x 


import math
from EZCAD3 import *

class Base(Part):
    def __init__(self, up):
	""" *Base*: initialize. """
	Part.__init__(self, up)
	self.dx_l = L(inch=12)
	self.dy_l = L(inch=12)
	self.z0_l = L()
	self.dz_l = L(mm=5.0)

    def configure(self, dx = None, dy = None, z0 = None, dz = None):
	""" *Base*: configure. """
	# Check argument types:
	none_type = type(None)
	assert type(dx) == none_type or isinstance(dx, L)
	assert type(dy) == none_type or isinstance(dx, L)
	assert type(z0) == none_type or isinstance(z0, L)
	assert type(dz) == none_type or isinstance(dz, L)

	if isinstance(dx, L):
	    self.dx_l = dx
	if isinstance(dy, L):
	    self.dy_l = dy
	if isinstance(z0, L):
	    self.z0_l = z0
	if isinstance(dz, L):
	    self.dz_l = dz

    def construct(self):
	""" *Base*: consturct. """
	dx = self.dx_l
	dy = self.dy_l
	dz = self.dz_l

	x0 = -(dx / 2.0)
	x1 =   dx / 2.0

	y0 = -(dy / 2.0)
	y1 =   dy / 2.0

	z0 = self.z0_l
	z1 = z0 + dz / 2
	z2 = z1 + dz / 2

	#print("Base.construct(): x0={0} x1={1} y0={2} y1={3}". \
	#  format(x0, x1, y0, y1))

	z = L()
	outer_contour = Contour(name = "Base")
	outer_contour.bend_append(P(x0, y0, z0), z)
	outer_contour.bend_append(P(x0, y1, z0), z)
	outer_contour.bend_append(P(x1, y1, z0), z)
	outer_contour.bend_append(P(x1, y0, z0), z)
	
	x_pitch = L(inch=1)
	y_pitch = L(inch=1)
	column_count = int(dx / x_pitch)
	column_span = x_pitch * (column_count - 1)
	x_start = -column_span / 2.0
	row_count = int(dy / y_pitch)
	row_span = y_pitch * (row_count - 1)
	y_start = -row_span / 2.0

	#print("x_start={0} column_span={1} column_count={2}".
	#  format(x_start, column_span, column_count))
	#print("y_start={0} row_span={1} row_count={2}".
	#  format(y_start, row_span, row_count))

	hole_dx = L(mm=3)
	hole_dy = L(mm=3)

	inner_contours = []
	for row in range(row_count):
	    y = y_start + y_pitch * row
	    for column in range(column_count):
		x = x_start + x_pitch * column

		# Hole center is at (x, y):
		x0 = x - hole_dx / 2
		x1 = x + hole_dx / 2
		y0 = y - hole_dy / 2
		y1 = y + hole_dy / 2

		# Construct the *hole_contour*:
		hole_contour = Contour(
		  name = "Base Hole[{0},{1}]".format(row, column))
		hole_contour.bend_append(P(x0, y0, z0), z)
		hole_contour.bend_append(P(x0, y1, z0), z)
		hole_contour.bend_append(P(x1, y1, z0), z)
		hole_contour.bend_append(P(x1, y0, z0), z)

		inner_contours.append(hole_contour)

	self.extrude(comment = "Base",
	  material = Material("wood", "plywood"),
	  color = Color("beige"),
	  outer_contour = outer_contour,
	  inner_contours = inner_contours,
	  start = P(z, z, z0),
	  end = P(z, z, z0 + dz))

	self.dxf_write(center = P(0, 0, z1))

class Castor_Assembly(Part):
    def __init__(self, up):
	""" *Castor_Assembly*: Initialization """
	self.thin_l = L(mm=3)
	self.thick_l = L(mm=5)
	Part.__init__(self, up)
	self.bottom_side_ = Castor_Bottom(self)
	self.corner1_p = P()
	self.corner2_p = P()
	self.east_side_ = Castor_East_West(self)
	self.north_side_ = Castor_North_South(self)
	self.plate_ = Castor_Plate(self)
	self.screw_bne_ = Fastener(self, "CA.side_bne")
	self.screw_bnw_ = Fastener(self, "CA.side_bnw")
	self.screw_bse_ = Fastener(self, "CA.side_bse")
	self.screw_bsw_ = Fastener(self, "CA.side_bsw")
	self.screw_tne_ = Fastener(self, "CA.side_tne")
	self.screw_tnw_ = Fastener(self, "CA.side_tnw")
	self.screw_tse_ = Fastener(self, "CA.side_tse")
	self.screw_tsw_ = Fastener(self, "CA.side_tsw")
	self.screw_bottom_ne_ = Fastener(self, "CA.bottom_ne")
	self.screw_bottom_nw_ = Fastener(self, "CA.bottom_ne")
	self.screw_bottom_se_ = Fastener(self, "CA.bottom_se")
	self.screw_bottom_sw_ = Fastener(self, "CA.bottom_sw")
	self.screw_top_e_ = Fastener(self, "CA.top_e")
	self.screw_top_w_ = Fastener(self, "CA.top_w")
	self.south_side_ = Castor_North_South(self)
	self.west_side_ = Castor_East_West(self)
	self.wheel_ = Castor_Wheel(self)

    def configure(self, corner1 = P(), corner2 = P()):
	""" Castor_Assembly: Configuration. """

	# Do argument type checking:
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)

	# Verify that corners are ordered properly:
	assert corner1.x < corner2.x
	assert corner1.y < corner2.y
	assert corner1.z < corner2.z

	# Load up *corner1* and *corner2*:
	self.corner1_p = corner1
	self.corner2_p = corner2

    def construct(self):
	""" Castor_Assembly: Construction. """

	# Grab some values from *castor* (== *self*):
	bottom_side = self.bottom_side_
	corner1 = self.corner1_p
	corner2 = self.corner2_p
	east_side = self.east_side_
	north_side = self.north_side_
	plate = self.plate_
	screw_bottom_ne = self.screw_bottom_ne_
	screw_bottom_nw = self.screw_bottom_nw_
	screw_bottom_se = self.screw_bottom_se_
	screw_bottom_sw = self.screw_bottom_sw_
	screw_bne = self.screw_bne_
	screw_bnw = self.screw_bnw_
	screw_bse = self.screw_bse_
	screw_bsw = self.screw_bsw_
	screw_tne = self.screw_tne_
	screw_tnw = self.screw_tnw_
	screw_tse = self.screw_tse_
	screw_tsw = self.screw_tsw_
	screw_top_e = self.screw_top_e_
	screw_top_w = self.screw_top_w_
	south_side = self.south_side_
	thick = self.thick_l
	thin = self.thin_l
	west_side = self.west_side_
	wheel = self.wheel_

	screw_length = L(mm=10)

	# Compute some X locations:
	## Start from east working towards X center:
	self.x0_l = x0 = corner1.x		# East surface of *east_side*
	self.x1_l = x1 = x0 + thin		# West surface of *east_side*
	self.x2_l = x2 = x0 + screw_length	# End of *east_side* screws
	self.x3_l = x3 = x2 + 2 * thick		# East tongue of *bottom_side*
	## Now start from west woriking towards center:
	self.x8_l = x8 = corner2.x		# West surface of *west_side*
	self.x7_l = x7 = x8 - thin		# East surface of *west_side*
	self.x6_l = x6 = x7 - screw_length	# End of *west_side* screws
	self.x5_l = x5 = x6 - 2 * thick		# West tongue of *bottom_side*
	self.x4_l = x4 = (x0 + x8)/2		# X Box center

	# Compute some Y locations:
	## Start from south working towards Y center:
	self.y0_l  = y0 = corner1.y		# South surface of *east_side*
	self.y1_l  = y1 = y0 + thin		# South surface of *south_side*
	self.y2_l  = y2 = y1 + thick/2		# Center of *south_side*
	self.y3_l  = y3 = y2 + thick/2		# North surface of *south_side*
	self.y4_l  = y4 = y3 + thick		# South tongue of *bottom_side*
	## Now start from north working towards Y center:
	self.y10_l = y10 = corner2.y		# North surface of *east_side*
	self.y9_l  = y9  = y10 - thin		# North surface of *north_side*
	self.y8_l  = y8  = y9 - thick/2		# Center of *north_side*
	self.y7_l  = y7  = y8 - thick/2		# South surface of *north_side*
	self.y6_l  = y6  = y7 - thick		# North tongue of *bottom_side*
	self.y5_l  = y5  = (y0 + y10)/2		# Y Box Center

	# Compute some Z locations:
	## Start from bottom working towards Z center:
	self.z0_l = z0 = corner1.z		# floor surface
	self.z1_l = z1 = z0 + wheel.diameter_l/2 # Wheel axel
	self.z2_l = z2 = z1 + wheel.diameter_l/2 # Bottom surface of *east_side*
	self.z3_l = z3 = z2 + thin		# Bottom of *bottom_side*
	self.z4_l = z4 = z3 + thick/2		# Center of *bottom_side*
	self.z5_l = z5 = z4 + thick/2		# Top surface of *bottom_side*
	self.z6_l = z6 = z5 + thick		# Bottom surface of *south_side*
	## Now from top working towards Z center:
	self.z10_l = z10 = corner2.z		# Top surface of caster asm.
	self.z9_l = z9 = z10 - screw_length + thick # End of top screws
	self.z8_l = z8 = z9 - thick		# Top tongue of *south_side*
	self.z7_l = z7 = (z2 + z8)/2		# Z Box Center

	# Configure everything:
	bottom_side.configure("Bottom Side", z3, z5)
	east_side.configure("West Side", x7, x8)
	west_side.configure("East Side", x0, x1)
	south_side.configure("South Side", y1, y3)
	north_side.configure("North Side", y7, y9)
	plate.configure(P(x4, y5, z3))
	wheel.configure(P(x4, y5, z3 - wheel.diameter_l/2))

	# Install the 8 screws that hold the castor sides together:
	screw_flags = "#4-40"
	screw_bne.configure(comment = "Castor BNE",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x8, y8, z4),
	  end = P(x6, y8, z4),
	  flags = screw_flags)
	screw_bne.drill(part = east_side, select = "close")
	screw_bne.nut_ledge(part = north_side, flags="N")

	screw_bnw.configure(comment = "Castor BNW",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x0, y8, z4),
	  end = P(x2, y8, z4),
	  flags = screw_flags)
	screw_bnw.drill(part = west_side, select = "close")
	screw_bnw.nut_ledge(part = north_side, flags="N")
	
	screw_bse.configure(comment = "Castor BSE",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x8, y2, z4),
	  end = P(x6, y2, z4),
	  flags = screw_flags)
	screw_bse.drill(part = east_side, select = "close")
	screw_bse.nut_ledge(part = south_side, flags="S")

	screw_bsw.configure(comment = "Castor BSW",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x0, y2, z4),
	  end = P(x2, y2, z4),
	  flags = screw_flags)
	screw_bsw.drill(part = west_side, select = "close")
	screw_bsw.nut_ledge(part = south_side, flags="S")

	screw_tne.configure(comment = "Castor TNE",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x8, y8, z9),
	  end = P(x6, y8, z9),
	  flags = screw_flags)
	screw_tne.drill(part = east_side, select = "close")
	screw_tne.nut_ledge(part = north_side, flags="N")

	screw_tnw.configure(comment = "Castor TNW",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x0, y8, z9),
	  end = P(x2, y8, z9),
	  flags = screw_flags)
	screw_tnw.drill(part = west_side, select = "close")
	screw_tnw.nut_ledge(part = north_side, flags="N")
	
	screw_tse.configure(comment = "Castor TSE",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x8, y2, z9),
	  end = P(x6, y2, z9),
	  flags = screw_flags)
	screw_tse.drill(part = east_side, select = "close")
	screw_tse.nut_ledge(part = south_side, flags="S")

	screw_tsw.configure(comment = "Castor TSW",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x0, y2, z9),
	  end = P(x2, y2, z9),
	  flags = screw_flags)
	screw_tsw.drill(part = west_side, select = "close")
	screw_tsw.nut_ledge(part = south_side, flags="S")
	
	# Drill 4 holes for the creeper castor:
	pitch_dx = L(inch=1)
	pitch_dy = L(inch="1-5/8")
	screw_bottom_ne.configure(comment = "Caster Bottom NE Screw",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x4 + pitch_dx/2, y5 + pitch_dy/2, z5),
	  end =   P(x4 + pitch_dx/2, y5 + pitch_dy/2, z2),
	  flags = screw_flags)
	screw_bottom_ne.drill(part = bottom_side, select = "close")
	screw_bottom_ne.drill(part = plate, select = "close")

	screw_bottom_nw.configure(comment = "Caster Bottom NW Screw",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x4 - pitch_dx/2, y5 + pitch_dy/2, z5),
	  end =   P(x4 - pitch_dx/2, y5 + pitch_dy/2, z2),
	  flags = screw_flags)
	screw_bottom_nw.drill(part = bottom_side, select = "close")
	screw_bottom_nw.drill(part = plate, select = "close")

	screw_bottom_se.configure(comment = "Caster Bottom SE Screw",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x4 + pitch_dx/2, y5 - pitch_dy/2, z5),
	  end =   P(x4 + pitch_dx/2, y5 - pitch_dy/2, z2),
	  flags = screw_flags)
	screw_bottom_se.drill(part = bottom_side, select = "close")
	screw_bottom_se.drill(part = plate, select = "close")

	screw_bottom_sw.configure(comment = "Caster Bottom SW Screw",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P(x4 - pitch_dx/2, y5 - pitch_dy/2, z5),
	  end =   P(x4 - pitch_dx/2, y5 - pitch_dy/2, z2),
	  flags = screw_flags)
	screw_bottom_sw.drill(part = bottom_side, select = "close")
	screw_bottom_sw.drill(part = plate, select = "close")

	# Install the top mount screws:
	screw_top_e.configure(comment = "Caster Top North Screw",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P((x7 + x8)/2, y5, z10 + thick * 3),
	  end =   P((x7 + x8)/2, y5, z9),
	  flags = screw_flags)
	screw_top_e.nut_ledge(part = east_side, flags="E")
	  
	screw_top_w.configure(comment = "Caster Top North Screw",
	  color = Color("black"),
	  material = Material("Steel", "Stainless"),
	  start = P((x0 + x1)/2, y5, z10 + thick * 3),
	  end =   P((x0 + x1)/2, y5, z9),
	  flags = screw_flags)
	screw_top_w.nut_ledge(part = west_side, flags="E")

	#north_side.invisible_set()

class Castor_Bottom(Part):
    def __init__(self, up):
	""" *Castor_Bottom*: Initialization. """
	Part.__init__(self, up)
	self.name_s = "<no name>"
	self.z2_l = L()
	self.z3_l = L()

    def configure(self, name, z2, z3):
	""" *Castor_Bottom: Configuration. """
	assert isinstance(name, str)
	assert isinstance(z2, L)
	assert isinstance(z3, L)
	self.name_s = name
	self.z2_l = z2
	self.z3_l = z3


    def construct(self):
	""" *Castor_Bottom*: Construction. """

	# Grab some values from *self*:
	assembly = self.up
	z2 = self.z2_l
	z3 = self.z3_l

	# Grab some valus from *assembly*:
	x0 = assembly.x0_l
	x1 = assembly.x1_l
	x3 = assembly.x3_l
	x4 = assembly.x4_l
	x5 = assembly.x5_l
	x7 = assembly.x7_l
	x8 = assembly.x8_l

	y1 = assembly.y1_l
	y3 = assembly.y3_l
	y4 = assembly.y4_l
	y5 = assembly.y5_l
	y6 = assembly.y6_l
	y7 = assembly.y7_l
	y9 = assembly.y9_l

	# Do the initial block:
	self.block(comment = "Bottom Side",
	  color = Color("green"),
	  material = Material("wood", "plywood"),
	  corner1 = P(x1, y3, z2),
	  corner2 = P(x7, y7, z3),
	  top = "t")
	#print("bottom:c1={0} c2={1}".format(P(x1, y1, z2), P(x7, y7, z3)))

	# Do the West tongue:
	self.block(comment = "Bottom Side West Tongue",
	  corner1 = P(x0, y4, z2),
	  corner2 = P(x1, y6, z3),
	  welds = "e",
	  top = "t")

	# Do the East tongue:
	self.block(comment = "Bottom Side West Tongue",
	  corner1 = P(x7, y4, z2),
	  corner2 = P(x8, y6, z3),
	  welds = "w",
	  top = "t")

	# Do the South tongue:
	self.block(comment = "Bottom Side South Tongue",
	  corner1 = P(x3, y1, z2),
	  corner2 = P(x5, y3, z3),
	  welds = "n",
	  top = "t")

	# Do the North tongue:
	self.block(comment = "Bottom Side North Tongue",
	  corner1 = P(x3, y7, z2),
	  corner2 = P(x5, y9, z3),
	  welds = "s",
	  top = "t")

	self.dxf_write(center = P(x4, y5, (z2 + z3)/2),
	  plane_normal = P(0, 0, L(mm=1)))

class Castor_East_West(Part):
    def __init__(self, up):
	""" *Castor_East_West*: Initialization. """
	Part.__init__(self, up)
	self.name_s = "<no name>"
	self.x0_l = L()
	self.x1_l = L()

    def configure(self, name, x0, x1):
	""" *Castor_East_West*: Initialization. """
	assert isinstance(name, str)
	assert isinstance(x0, L)
	assert isinstance(x1, L)
	self.name_s = name
	self.x0_l = x0
	self.x1_l = x1

    def construct(self):
	""" *Castor_East_West*: Construction. """
	assembly = self.up
	name = self.name_s
	x0 = self.x0_l
	x1 = self.x1_l

	y0 = assembly.y0_l
	y1 = assembly.y1_l
	y3 = assembly.y3_l
	y4 = assembly.y4_l
	y5 = assembly.y5_l
	y6 = assembly.y6_l
	y7 = assembly.y7_l
	y9 = assembly.y9_l
	y10 = assembly.y10_l

	z2 = assembly.z2_l
	z3 = assembly.z3_l
	z5 = assembly.z5_l
	z6 = assembly.z6_l
	z7 = assembly.z7_l
	z8 = assembly.z8_l
	z10 = assembly.z10_l

	self.block(comment = name,
	  color = Color("red"),
	  material = Material("wood", "plywood"),
	  corner1 = P(x0, y0, z2),
	  corner2 = P(x1, y10, z10),
	  top = "e")

	# Open bottom tongue slot:
	self.simple_pocket(comment = "Bottom Tongue Slot",
	  bottom_corner = P(x0, y4, z3),
	  top_corner =    P(x1, y6, z5),
	  pocket_top = "e")

	# Remove the south tongue slot:
	self.simple_pocket(comment = "{0} South Tongue Slot".format(name),
	  bottom_corner = P(x0, y1, z6),
	  top_corner =    P(x1, y3, z8),
	  pocket_top = "e")

	# Remvoe the north tongue:
	self.simple_pocket(comment = "{0} North Tongue Slot".format(name),
	  bottom_corner = P(x0, y7, z6),
	  top_corner =    P(x1, y9, z8),
	  pocket_top = "e")
	  
	self.dxf_write(center = P((x0 + x1)/2, y4, z7),
	  plane_normal=P(L(mm=1), 0, 0))


class Castor_North_South(Part):
    def __init__(self, up):
	""" *Castor_North_South*: Initialization. """
	Part.__init__(self, up)
	self.name = "<no name>"
	self.y0_l = L()
	self.y1_l = L()

    def configure(self, name, y0, y1):
	""" *Castor_North_South*: Initialization. """
	assert isinstance(name, str)
	assert isinstance(y0, L)
	assert isinstance(y1, L)
	self.name_s = name
	self.y0_l = y0
	self.y1_l = y1

    def construct(self):
	""" *Castor_North_South*: Construction. """

	# Grab some values from *self*:
	name = self.name_s
	y0 = self.y0_l
	y1 = self.y1_l
	assembly = self.up

	# Grab some values from *assembly*:
	x0 = assembly.x0_l
	x1 = assembly.x1_l
	x3 = assembly.x3_l
	x4 = assembly.x4_l
	x5 = assembly.x5_l
	x7 = assembly.x7_l
	x8 = assembly.x8_l

	z2 = assembly.z2_l
	z3 = assembly.z3_l
	z5 = assembly.z5_l
	z6 = assembly.z6_l
	z7 = assembly.z7_l
	z8 = assembly.z8_l
	z10 = assembly.z10_l

	# The initial chunk of plywood fits between the east and west sides:
	self.block(comment = name,
	  color = Color("blue"),
	  material = Material("wood", "plywood"),
	  corner1 = P(x1, y0, z2),
	  corner2 = P(x7, y1, z10),
	  top = "n")

	# Cut open bottom tongue slot:
	self.simple_pocket(comment = "Bottom Tongue Slot",
	  bottom_corner = P(x3, y0, z3),
	  top_corner =    P(x5, y1, z5),
	  pocket_top = "n")

	# Add in an west tongue:
	self.block(comment = "{0} West Tongue",
	  corner1 = P(x0, y0, z6),
	  corner2 = P(x1, y1, z8),
	  welds = "e",
	  top = "n")

	# Add in an east tongue:
	self.block(comment = "{0} East Tongue",
	  corner1 = P(x7, y0, z6),
	  corner2 = P(x8, y1, z8),
	  welds = "w",
	  top = "n")
	  
	self.dxf_write(center = P(x4, (y0+y1)/2, z7),
	  plane_normal = P(0, L(mm=1), 0))

class Castor_Plate(Part):
    def __init__(self, up):
	""" *Castor_Plate*: Initalization. """
	Part.__init__(self, up)
	self.center_p = P()

    def configure(self, center = P()):
	""" *Castor_Plate*: Configuration. """
	assert isinstance(center, P)
	self.center_p = center

    def construct(self):
	""" *Castor_Plate*: Construction. """
	center = self.center_p
	x = center.x
	y = center.y
	z = center.z
	
	plate_dx = L(inch="1-9/16")
	plate_dy = L(inch="2-1/4")
	plate_dz = L(inch="5/32")
	self.block(comment = "Castor Plate",
	  color = Color("yellow"),
	  material = Material("plastic", "unknown"),
	  corner1 = P(x - plate_dx/2, y - plate_dy/2, z - plate_dz),
	  corner2 = P(x + plate_dx/2, y + plate_dy/2, z),
	  top = "t")

class Castor_Wheel(Part):
    def __init__(self, up):
	""" *Castor_Wheel*: Initialization. """
	Part.__init__(self, up)
	self.center_p = P()

    def configure(self, center = P()):
	""" *Castor_Wheel*: Configuration. """
	assert isinstance(center, P)
	self.center_p = center

    def construct(self):
	""" *Castor_Wheel*: Construction. """
	center = self.center_p
	self.diameter_l = diameter = L(inch="2-3/4")
	dy = L(inch=.5)
	self.cylinder(comment = "Castor Wheel",
	  color = Color("orange"),
	  material = Material("Plastic", "Nylon"),
	  diameter = diameter,
	  start = P(center.x, center.y - dy/2, center.z),
	  end =   P(center.x, center.y + dy/2, center.z))

class Dual_Slot_Encoder(Part):
    def __init__(self, up):
	""" Dual Slot Encoder: initialize """
	Part.__init__(self, up)
	self.dual_slot_encoder_pcb_ = Dual_Slot_Encoder_PCB(self)
	self.slot_interrupter1_ = slot_interrupter1 = Slot_Interrupter(self)
	self.slot_interrupter2_ = slot_interrupter2 = Slot_Interrupter(self)
	self.disk_offset_dy_l = L(mm=3)
	self.center_p = P()

    def configure(self, center = None):
	""" Dual Slot Encoder: configure """
	assert isinstance(center, P)
	self.center_p = center
	
    def construct(self):
	""" Dual Slot Encoder: construct """
	disk_offset_dy = self.disk_offset_dy_l
	dual_slot_encoder_pcb = self.dual_slot_encoder_pcb_
	motor_assembly = self.up
	slot_interrupter1 = self.slot_interrupter1_
	slot_interrupter2 = self.slot_interrupter2_
	center = self.center_p

	x = center.x
	y = center.y
	z = center.z

	dual_slot_encoder_pcb.configure(
	  center = P(center.x, center.y, center.z))
	slot_interrupter1.configure(center =
	  P(x + slot_interrupter1.dx_l/2, y + disk_offset_dy, z))
	slot_interrupter2.configure(center =
	  P(x - slot_interrupter2.dx_l/2, y + disk_offset_dy, z))

	#dual_slot_encoder_pcb.invisible_set()

class Dual_Slot_Encoder_PCB(Part):
    """ Dual_Slot_Encoder_PCB: initialize """
    def __init__(self, up):
	Part.__init__(self, up)
	self.diameter_l = L(mm=3.175)
	self.dx_l = L(mm=50)
	self.dy_l = L(mm=25)
	self.dz_l = L(mm=1)
	self.holes_pitch_dx_l = L(mm=40)
	self.holes_pitch_dy_l = L(mm=15)
	self.center_p = P()

    def configure(self, center = None):
	""" Dual_Slot_Encoder_PCB: configure """
	assert isinstance(center, P)
	self.center_p = center

    def construct(self):
	""" Dual_Slot_Encoder_PCB: construct """
	diameter = self.diameter_l
	dx = self.dx_l
	dy = self.dy_l
	dz = self.dz_l
	center = self.center_p

	x0 = center.x - dx/2
	x1 = center.x + dx/2
	y0 = center.y - dy/2
	y1 = center.y + dy/2
	z0 = center.z
	z1 = center.z + dz

	# Create the PCB:
	z = L()
	self.block(comment = "Dual_Slot_Encoders_PCB",
	  material = Material("fiberglass", "FR4"),
	  color = Color("green"),
	  corner1 = P(center.x - dx/2, y0, z0),
	  corner2 = P(center.x + dx/2, y1, z1),
	  top = "t")

	# Drill the mounting holes:
	hx0 = x0 + L(mm=5)
	hx1 = x1 - L(mm=5)
	hy0 = y0 + L(mm=5)
	hy1 = y1 - L(mm=5)
	holes = [(hx0, hy0, "SW"),
	         (hx0, hy1, "NW"),
		 (hx1, hy0, "SE"),
		 (hx1, hy1, "NE")]
	for hole in holes:
	    x = hole[0]
	    y = hole[1]
	    text = hole[2]
	    #print("[{0}]: x={1} y={2}".format(text, x, y))
	    self.hole(comment = "{0} Mounting Hole".format(text), 
	      diameter = diameter,
	      start = P(x, y, z1),
	      end =   P(x, y, z0),
	      flags = "t")

class Encoder_Disk(Part):
    def __init__(self, up):
	""" *Encoder_Disk*: Initialization. """
	Part.__init__(self, up)
	self.diameter_l = diameter = L(mm=47) #L(mm=50)
	radius = diameter / 2
	self.actual_radius_l = actual_radius = radius - L(mm=0.75)
	self.slot_inner_radius_l = actual_radius - L(mm=5.00)
	self.slot_outer_radius_l = actual_radius - L(mm=0.50)
	self.slots_count_i = 36
	self.dy_l = L(mm=1.00)
	self.y_center_l = L(0)

    def configure(self, center = P()):
	""" *Encoder_Disk*: Configuration. """
	assert isinstance(center, P)
	self.center_p = center

    def construct(self):	
	""" *Encoder_Disk*: Construction. """
	# Grab some values from *self*:
	center = self.center_p
	diameter = self.diameter_l
	dy = self.dy_l
	actual_radius = self.actual_radius_l
	slots_count = self.slots_count_i
	slot_inner_radius = self.slot_inner_radius_l
	slot_outer_radius = self.slot_outer_radius_l

	pi = Angle.PI
        slot_length = slot_outer_radius - slot_inner_radius
	slot_median_radius = slot_inner_radius + .75 * slot_length
	slot_median_circumference = 2 * pi * slot_median_radius

	# Basic concept behind an encoder disk is that it consists of a
	# bunch of beam interrupter teeth interspersed with gaps between
	# the teeth.  The gaps are called slots and the non-gaps are
	# called teeth.  Ideally, the slot width is the same as the
	# no tooth width.  The way we think about the disk is that it
	# has 4 places of interest
        # 
	# * Slot center
	# * Slot to tooth edge
	# * Tooth center
	# * Tooth to Slot edge.
	#
	# Ideally these locations are equidistant apart.
	#
	# The encoder consists of 2 IR (Infrared) slot interrupters
	# with a rather narrow IR beam (.3mm) that can be occluded
	# by an encoder disk tooth.  What we want to ensure is that
	# one IR beam is pointing at an edge and the other beam
	# is pointing at either a tooth or a slot center.  This will
	# ensure that we get a nice quadrature signal out.
	#
	# Now we can go through the math.  Let R be the radius from the disk
	# center to the IR beam.  The circumference of the IR beam circle
	# is 2*PI*R.  Let S be the number of slots (which happens to equal
	# the number of teeth.)  2*PI*R/S is the approximate slot center
	# to slot center distance.  We divide this by 4 to get a slot
	# center to edge distance.  D = 2*PI*R/(4*S) = PI*R/(2*S).  Let
	# B be the distance between the beams.  We need the following
	# to be true:
	#
	# N*D ~= B, where N is odd (= 1, 3, 5, 7, 9, ...)
	#
	# In our particular case, R and B are fixed, so we can only
	# adjust N and S.  This means:
	#
	#    N*D = B				(1)
	#    N*PI*R/(2*S) = B			(2)
	#    N*PI*R = 2*S*B			(3)
	#    (N*PI*R)/(2*B) = S			(4)

	r = slot_median_radius = slot_inner_radius + .75 * slot_length
	b = distance_between_slot_interrupter_beams = L(mm=3.50)
        # = 3	# s = 30
	n = 5	# s = 50
	#n = 7	# s = ??
	pi = math.pi
	fractional_slots = (n * pi * r) / (2 * b)
	self.slot_count_i = slots_count = int(round(fractional_slots))
	print("b={0} r={1} n={2} frac_slots={3} slots_count={4}".
	  format(b, r, n, fractional_slots, slots_count))

	# Grab *x*, *y*, and *z* from *center*
	x = center.x
	y = center.y
	z = center.z

	# Compute some y coordinates:
	y0 = y - dy/2
	y1 = y
	y2 = y + dy/2

	zero = L()
	two_pi = Angle.TWO_PI
	slot_angle = Angle(rad=two_pi) / slots_count
	half_slot_angle = slot_angle / 2
	#print("slot_angle={0} half_slot_angle={1}".
	#  format(slot_angle, half_slot_angle))

	slot_contours = []
	outer_contour = Contour(name = "Encoder Disk")
	for slot_index in range(slots_count):
	    # Compute the angles for this iteration:
	    angle0 = slot_angle * slot_index
	    angle1 = angle0 + half_slot_angle

	    # Compute X/Y coordinates for the slot corners:
	    x0 = x + actual_radius.cosine(angle0)
	    z0 = z + actual_radius.sine(angle0)
	    x3 = x + actual_radius.cosine(angle1)
	    z3 = z + actual_radius.sine(angle1)

	    # Append (x,y) for the outer contour:
	    outer_contour.bend_append(P(x0, y0, z0), zero,
	      name="[{0}].xz0".format(slot_index))
	    outer_contour.bend_append(P(x3, y0, z3), zero,
	      name="[{0}].xz3".format(slot_index))

	    x0 = x + slot_outer_radius.cosine(angle0)
	    z0 = z + slot_outer_radius.sine(angle0)
	    x1 = x + slot_inner_radius.cosine(angle0)
	    z1 = z + slot_inner_radius.sine(angle0)
	    x2 = x + slot_inner_radius.cosine(angle1)
	    z2 = z + slot_inner_radius.sine(angle1)
	    x3 = x + slot_outer_radius.cosine(angle1)
	    z3 = z + slot_outer_radius.sine(angle1)

	    #print("[{0}]:angle0={1} angle1={2}".
	    #  format(slot_index, angle0, angle1))
	    #print("[{0}]:x0={1} z0={2}".format(slot_index, x0, z0))
	    #print("[{0}]:x1={1} z1={2}".format(slot_index, x1, z1))
	    #print("[{0}]:x2={1} z2={2}".format(slot_index, x2, z2))
	    #print("[{0}]:x3={1} z3={2}".format(slot_index, x3, z3))

	    # Append a *slot_contour* to *slot_contours*:
	    slot_contour = Contour(name = "Encoder_Disk_Slot[{0}]".
	      format(slot_index))
	    slot_contour.bend_append(P(x0, y0, z0), zero,
	      name="[{0}].0".format(slot_index))
	    slot_contour.bend_append(P(x1, y0, z1), zero,
	      name="[{0}].1".format(slot_index))
	    slot_contour.bend_append(P(x2, y0, z2), zero,
	      name="[{0}].2".format(slot_index))
	    slot_contour.bend_append(P(x3, y0, z3), zero,
	      name="[{0}].2".format(slot_index))
	    slot_contours.append(slot_contour)

	self.extrude(comment = "Encoder Disk",
	  material = Material("Plastic", "Derlin"),
	  color = Color("yellow"),
	  outer_contour = outer_contour,
	  inner_contours = slot_contours,
	  start = P(x, y0, z),
	  end   = P(x, y2, z))

	#diameter = L(inch=0.0890),    # #2-56 close fit
	self.hole(comment = "Screw_Hole",
	  diameter = L(mm=1.75),       # Much tighter fit
	  start = P(x, y2, z),
	  end =   P(x, y0, z),
	  flags = "t")

	self.dxf_write(center = P(x, y1, z), plane_normal = P(0, L(mm=1), 0))

class Encoder_Mount(Part):
    def __init__(self, up):
	""" Encoder_Mount: initialize """
	self.thick_l = up.thick_l
	self.thin_l = up.thin_l
	Part.__init__(self, up)
	self.dx_l = L(mm=60)
	self.dy_l = L(mm=35)
	self.dz_l = self.thick_l
	self.tongue_dx_l = L(mm=50)
	self.center_p = P()

    def configure(self, dx = None, dy = None, dz = None,
      tongue_dx = None, center = None):
	""" Encoder_Mount: configure """
	none_type = type(None)
	assert type(dx) == none_type or isinstance(dx, L)
	assert type(dy) == none_type or isinstance(dy, L)
	assert type(dz) == none_type or isinstance(dz, L)
	assert type(tongue_dx) == none_type or isinstance(tongue_dx, L)
	assert isinstance(center, P)

	if isinstance(dx, L):
	    self.dx_l = dx
	if isinstance(dy, L):
	    self.dy_l = dy
	if isinstance(dz, L):
	    self.dz_l = dz
	if isinstance(tongue_dx, L):
	    self.tongue_dx_l = tongue_dx
	self.center_p = center

    def construct(self):
	""" Encoder_Mount: construct """
	dx = self.dx_l
	dy = self.dy_l
	dz = self.dz_l
	tongue_dx = self.tongue_dx_l
	center = self.center_p
	thick = self.thick_l
	thin = self.thin_l

	# Go get some related *Part*'s:
	motor_assembly = self.up
	dual_slot_encoder = motor_assembly.dual_slot_encoder_
	dual_slot_encoder_pcb = dual_slot_encoder.dual_slot_encoder_pcb_

	# Grab some values from the Encoder PCB:
	pcb_y_center = dual_slot_encoder_pcb.center_p.y
	pcb_dx = dual_slot_encoder_pcb.dx_l
	pcb_dy = dual_slot_encoder_pcb.dy_l
	holes_pitch_dx = dual_slot_encoder_pcb.holes_pitch_dx_l
	holes_pitch_dy = dual_slot_encoder_pcb.holes_pitch_dy_l

	x = center.x
	y = center.y
	z = center.z

	# Compute some X locations in ascending order:
	x0 = x - dx/2			# West mount edge
	x1 = x - tongue_dx/2		# West tongue edge
	x2 = x - pcb_dx/2 - L(mm=15.0)	# West west_wire_hole edge
	x3 = x - pcb_dx/2 - L(mm=5.0)	# East west_wire_hole edge
	x4 = x - pcb_dx/2 + L(mm=2.5)	
	x5 = x - pcb_dx/2 + L(mm=8.5)
	x6 = x + pcb_dx/2 - L(mm=8.5)
	x7 = x + pcb_dx/2 - L(mm=2.5)
	x8 = x + pcb_dx/2 + L(mm=5.0)	# West east_wire_hole edge
	x9 = x + pcb_dx/2 + L(mm=15.0)	# East west_wire_hole edge
	x10 = x + tongue_dx/2		# East tongue edge
	x11 = x + dx/2			# East mount edge

	# Compute some Y locations in ascending order:
	y0 = y - dy/2
	y1 = y - dy/2 + L(mm=3)
	y2 = pcb_y_center - pcb_dy/2 + L(mm=2.0)
	y3 = pcb_y_center - pcb_dy/2 + L(mm=8.5)
	y4 = pcb_y_center + pcb_dy/2 - L(mm=8.5)
	y5 = pcb_y_center + pcb_dy/2 - L(mm=2.0)
	y6 = y + dy/2 - L(mm=3)
	y7 = y + dy/2
	#print("Encoder_Mount:construct: y0={0} y_center={1} y7={2} dy={3}".
	#  format(y0, y_center, y7, dy))

	# Compute some Z locations in ascending order:
	z0 = z - dz/2		# Bottom surface
	z1 = z + dz/2		# Top surface

	# Compute the outer contour:
	zero = L()
	outer_contour = Contour(name = "Encoder_Mount")
	# Corner: x0, y0:
	outer_contour.bend_append(P(x1, y0, z0), zero)
	outer_contour.bend_append(P(x1, y1, z0), zero)
	outer_contour.bend_append(P(x0, y1, z0), zero)

	# Corner: x0, y7:
	outer_contour.bend_append(P(x0, y6, z0), zero)
	outer_contour.bend_append(P(x1, y6, z0), zero)
	outer_contour.bend_append(P(x1, y7, z0), zero)

	# Corner: x11, x11:
	outer_contour.bend_append(P(x10, y7, z0), zero)
	outer_contour.bend_append(P(x10, y6, z0), zero)
	outer_contour.bend_append(P(x11, y6, z0), zero)

	# Corner: x11, y0:
	outer_contour.bend_append(P(x11, y1, z0), zero)
	outer_contour.bend_append(P(x10, y1, z0), zero)
	outer_contour.bend_append(P(x10, y0, z0), zero)

	# Compute the PCB cut out inner contour:
	pcb_cutout = Contour(name = "PCB_Cutout")
	curve_radius = L(mm=0.5)
	# Corner (x4, y2):
	pcb_cutout.bend_append(P(x5, y2, z0), curve_radius)
	pcb_cutout.bend_append(P(x5, y3, z0), curve_radius)
	pcb_cutout.bend_append(P(x4, y3, z0), curve_radius)
	# Corner (x4, y3):
	pcb_cutout.bend_append(P(x4, y4, z0), curve_radius)
	pcb_cutout.bend_append(P(x5, y4, z0), curve_radius)
	pcb_cutout.bend_append(P(x5, y5, z0), curve_radius)
	# Corner (x5, y3):
	pcb_cutout.bend_append(P(x6, y5, z0), curve_radius)
	pcb_cutout.bend_append(P(x6, y4, z0), curve_radius)
	pcb_cutout.bend_append(P(x7, y4, z0), curve_radius)
	# Corner (x5, y2):
	pcb_cutout.bend_append(P(x7, y3, z0), curve_radius)
	pcb_cutout.bend_append(P(x6, y3, z0), curve_radius)
	pcb_cutout.bend_append(P(x6, y2, z0), curve_radius)

	# West Wire hole:
	wire_hole_west = Contour(name = "West Wire Hole")
	wire_hole_west.bend_append(P(x2, y3, z0), curve_radius)
	wire_hole_west.bend_append(P(x2, y4, z0), curve_radius)
	wire_hole_west.bend_append(P(x3, y4, z0), curve_radius)
	wire_hole_west.bend_append(P(x3, y3, z0), curve_radius)

	# East Wire hole:
	wire_hole_east = Contour(name = "East Wire Hole")
	wire_hole_east.bend_append(P(x8, y3, z0), curve_radius)
	wire_hole_east.bend_append(P(x8, y4, z0), curve_radius)
	wire_hole_east.bend_append(P(x9, y4, z0), curve_radius)
	wire_hole_east.bend_append(P(x9, y3, z0), curve_radius)

	# Do the extrusion:
	self.extrude(comment = "PCB Mount Block",
	  material = Material("Wood", "Plywood"),
	  color = Color("purple"),
	  outer_contour = outer_contour,
	  inner_contours = [pcb_cutout, wire_hole_west, wire_hole_east],
	  start = P(x, y, z0),
	  end =   P(x, y, z1))

class GM3_Motor(Part):
    def __init__(self, up):
	""" *GM3_Motor*: Initialization. """
	Part.__init__(self, up)
	self.shaft_center = P()

    def configure(self, shaft_center = P()):
	""" *GM3_Motor*: Configuration. """
	none_type = type(None)
	assert isinstance(shaft_center, P)
	self.shaft_center_p = shaft_center

    def construct(self):
	""" *GM3_Motor*: Construction. """
	zero = L()

	# Grab some values from *self*:
	shaft_center = self.shaft_center_p

	# Diameters:
	self.back_holes_diameter_l =    back_holes_diameter =    L(mm=2.90)
	self.encoder_hub_diameter_l =   encoder_hub_diameter =   L(mm=9.73)
	self.encoder_shaft_diameter_l = encoder_shaft_diameter = L(mm=6.83)
	self.front_hole_diameter_l =    front_hole_diameter =    L(mm=2.76)
	self.pin_diameter_l =           pin_diameter =           L(mm=4.48)
	self.screw_hole_diameter_l =    screw_hole_diameter =    L(mm=2) # Guess
	self.wheel_hub_diameter_l =     wheel_hub_diameter =     L(mm=9.20)
	self.wheel_shaft_diameter_l =   wheel_shaft_diameter =   L(mm=6.90)

	# X dimensions:
	self.dx_l =            dx =            L(mm=69.38)
	self.front_dx_l =      front_dx =      L(mm=10.80)
	self.back_dx_l =       back_dx =       front_dx - L(mm=65.51)
	self.back_holes_dx_l = back_holes_dx = front_dx - L(mm=30.81)
	self.front_tip_dx_l =  front_tip_dx =  back_dx + dx
	self.front_hole_dx_l = front_hole_dx = front_dx + L(mm=2.62)
	self.pin_dx_l =        pin_dx =        front_dx - L(mm=22.23)
	self.strap_hole_dx_l = strap_hole_dx =  -L(mm=30.00)
	self.strap_pocket_dx_l = strap_pocket_dx = L(mm=4.00)

	# Y dimensions:
	self.wheel_shaft_dy_l =   wheel_shaft_dy =   L(mm=9.40)  # spec=9.20
	self.encoder_shaft_dy_l = encoder_shaft_dy = L(mm=10.20) # spec=8.93
	self.encoder_hub_dy_l =   encoder_hub_dy =   L(mm=2.00)  # Guess
	self.pin_dy_l =		  pin_dy =           L(mm=1.00)  # Guess
	self.dy_l =               dy =               L(mm=18.30) # spec=18.64

	# Z dimensions:
	self.dz_l =               dz =                     L(mm=22.23)
	self.back_holes_pitch_dz_l = back_holes_pitch_dz = L(mm=17.44)
	self.strap_pocket_dz_l = strap_pocket_dz =         L(mm=8.00)

	# Grab the center values:
	x = shaft_center.x
	y = shaft_center.y
	z = shaft_center.z

	# Various Y position is ascending values:
	y0 = y - dy/2 - wheel_shaft_dy	# Wheel shaft edge
	y1 = y - dy/2 - pin_dy		# Alignment pin edge
	y2 = y - dy/2			# Wheel side of motor body
	y3 = y				# Motor center
	y4 = y3 + dy/2			# Encoder side of motor body
	y5 = y4 + encoder_hub_dy	# Encoder hub edge
	y6 = y4 + encoder_shaft_dy	# Edge of encoder shaft

	# Contstuct the motor block:
	self.block(comment = "GM3 Motor Block",
	  material = Material("plastic", "ABS"),
	  color = Color("red"),
	  corner1 = P(x + back_dx,  y2, z - dz/2),
	  corner2 = P(x + front_dx, y4, z + dz/2),
	  top = "t")

	# Add various shafts, hubs, and pins:
	self.cylinder(comment = "wheel shaft",
	  diameter = wheel_shaft_diameter,
	  start = P(x, y0, z),
	  end =   P(x, y3, z))	# Weld
	self.cylinder(comment = "encoder shaft",
	  diameter = wheel_hub_diameter,
	  start = P(x, y6, z),
	  end =   P(x, y3, z))	# Weld
	self.cylinder(comment = "encoder hub",
	  diameter = encoder_hub_diameter,
	  start = P(x, y5, z),
	  end =   P(x, y3, z))	# Weld
	self.cylinder(comment = "align pin",
	  diameter = pin_diameter,
	  start = P(x + pin_dx, y1, z),
	  end =   P(x + pin_dx, y2, z)) # Weld

	# Drill some holes:
	self.hole(comment = "Mounting Hole 1",
	  diameter = back_holes_diameter,
	  start = P(x + back_holes_dx, y4, z + back_holes_pitch_dz/2),
	  end =   P(x + back_holes_dx, y2, z + back_holes_pitch_dz/2),
	  flags = "t")
	self.hole(comment = "Mounting Hole 2",
	  diameter = back_holes_diameter,
	  start = P(x + back_holes_dx, y4, z - back_holes_pitch_dz/2),
	  end =   P(x + back_holes_dx, y2, z - back_holes_pitch_dz/2),
	  flags = "t")
	self.hole(comment = "Shaft Center Screw Hole",
	  diameter = screw_hole_diameter,
	  start = P(x, y0, z),
	  end =   P(x, y5, z),
	  flags = "t")

class GM3_Wheel(Part):
    def __init__(self, up):
	""" *GM3_Wheel*: Initialization. """
	Part.__init__(self, up)
	self.center_p = P()
	self.diameter_l = L(mm=69)
	self.width_l = L(mm=7.62)

    def configure(self, center = P()):
	""" *GM3_Wheel*: Configuration. """
	assert isinstance(center, P)
	self.center_p = center

    def construct(self):
	""" *GM3_Wheel*: Construction. """

	# Grab some values from *self*:
	center = self.center_p
	diameter = self.diameter_l
	width = self.width_l

	# Extract *x*, *y*, *z* from *center*:
	x = center.x
	y = center.y
	z = center.z

	# Do the wheel:
	self.cylinder(comment = "Wheel",
	  material = Material("Plastic", "ABS"),
	  color = Color("blue"),
	  diameter = diameter,
	  start = P(x, y - width/2, z),
	  end =   P(x, y + width/2, z))

	# Do the screw hole for the wheel:
	self.hole(comment = "Screw Hole",
	  diameter = L(mm=8),
	  start = P(x, y + width/2, z),
	  end =   P(x, y - width/2, z),
	  flags = "t")

class Front_Back_Motor_Side(Part):
    def __init__(self, up):
	""" *Front_Back_Motor_Side*: initialize. """
	self.thick_l = up.thick_l
	self.thin_l = up.thin_l
	Part.__init__(self, up)
	self.dx_l = L(mm=5)
	self.dy_l = L(mm=50)
	self.dz_l = L(mm=100)
	self.inside_dy_l = L(mm=40)
	self.center_p = P()
	self.z_tongue_bottom = L(mm=5)
	self.z_tongue_top = L(mm=95)

    def configure(self, dx = None, dy = None, dz = None, inside_dy = None,
      z_tongue_bottom = None, z_tongue_top = None, center = None):
	""" *Front_Back_Motor_Side*: configure. """
	# Check argument types:
	none_type = type(None)
	assert type(dx) == none_type or isinstance(dx, L)
	assert type(dy) == none_type or isinstance(dy, L)
	assert type(dz) == none_type or isinstance(dz, L)
	assert type(inside_dy) == none_type or isinstance(inside_dy, L)
	assert isinstance(center, P)
	assert type(z_tongue_bottom) == none_type or \
	  isinstance(z_tongue_bottom, L)
	assert type(z_tongue_top) == none_type or isinstance(z_tongue_top, L)

	# Load up *self*:
	if isinstance(dx, L):
	    self.dx_l = dx
	if isinstance(dy, L):
	    self.dy_l = dy
	if isinstance(dz, L):
	    self.dz_l = dz
	if isinstance(inside_dy, L):
	    self.inside_dy_l = inside_dy
	self.center_p = center
	if isinstance(z_tongue_bottom, L):
	    self.z_tongue_bottom_l = z_tongue_bottom
	if isinstance(z_tongue_top, L):
	    self.z_tongue_top_l = z_tongue_top

    def construct(self):
	""" *Front_Back_Motor_Side*: construct. """
	# Grab some values from {self}.
	dx = self.dx_l
	dy = self.dy_l
	dz = self.dz_l
	inside_dy = self.inside_dy_l
	center = self.center_p
	z_tongue_bottom = self.z_tongue_bottom_l
	z_tongue_top = self.z_tongue_top_l
	zero = L()

	x_center = center.x
	y_center = center.y
	z_center = center.z

	# Some X locations:
	x0 = x_center - dx/2
	x1 = x_center
	x2 = x_center + dx/2

	# Some Y locations:
	y0 = y_center - dy/2		# South tongue Y
	y1 = y_center - inside_dy / 2	# South inside Y
	y2 = y_center			# Screw mount
	y3 = y_center + inside_dy / 2	# North inside Y
	y4 = y_center + dy/2		# North tongue Y
	#print("y_center={0} inside_dy={1} dy={2}".
	#  format(y_center, inside_dy, dy))
	#print("y0={0} y1={1} y2={2} y3={3} y4={4}".format(y0, y1, y2, y3, y4))

	# Some Z locations:
	z0 = z_center - dz/2	# Bottom
	z1 = z_tongue_bottom	# Bottom of tongue
	z2 = z_tongue_top	# Top of tongue
	z5 = z_center + dz/2	# Top
	z4 = z5 - L(mm=10)	# Upper mount screw
	z3 = z4 - L(mm=10)	# Upper mount screw

	#print("FBMS:x0={0} x1={1} x2={2}".format(x0, x1, x2))
	#print("FBMS:z0={0} z5={1}".format(z0, z5))

	# Do the outer contour:
	outer_contour = Contour(name = "Front_Back_Motor_Side")
	outer_contour.bend_append(P(0, y1, z0), zero)
	if True:
	    # Install inner tongue:
	    outer_contour.bend_append(P(0, y1, z1), zero)
	    outer_contour.bend_append(P(0, y0, z1), zero)
	    outer_contour.bend_append(P(0, y0, z2), zero)
	    outer_contour.bend_append(P(0, y1, z2), zero)
	outer_contour.bend_append(P(0, y1, z5), zero)
	outer_contour.bend_append(P(0, y3, z5), zero)
	if True:
	    # Install outer tongue:
	    outer_contour.bend_append(P(0, y3, z2), zero)
	    outer_contour.bend_append(P(0, y4, z2), zero)
	    outer_contour.bend_append(P(0, y4, z1), zero)
	    outer_contour.bend_append(P(0, y3, z1), zero)
	outer_contour.bend_append(P(0, y3, z0), zero)

	self.extrude(comment = "Front_Back_Motor_Side",
	  material = Material("Wood", "Plywood"),
	  color = Color("lime"),
	  outer_contour = outer_contour,
	  inner_contours = [],
	  start = P(x0, 0, 0),
	  end = P(x2, 0, 0))

class Inner_Outer_Motor_Side(Part):
    def __init__(self, up):
	self.thick_l = up.thick_l
	self.thin_l = up.thin_l
	Part.__init__(self, up)
	self.bottom_z_l =          L(mm=-5)
	self.dx_l =                L(mm=80)
	self.dy_l =                self.thin_l
	self.encoder_tongue_dz_l = L(mm=10)
	self.top_z_l =             L(mm=10)
	self.center_p =            P()
	self.side_tongue_dx_l =    L(mm=5)
	self.back_center_x_l =     L(mm=-20)
	self.front_center_x_l =    L(mm=20)
	self.shaft_z_l = L()
	self.side_tongue_bottom_z_l = L()
	self.side_tongue_top_z_l = L()

    def configure(self, bottom_z = None, dx = None,
      side_tongue_dx = None, side_tongue_bottom_z = None,
      side_tongue_top_z = None, back_center_x = None, front_center_x = None,
      top_z = None, center = None, shaft_z = None):
	""" Inner_Outer_Motor_Side: configure """
	# Check argument types:
	assert isinstance(bottom_z, L)
	assert isinstance(dx, L)
	assert isinstance(side_tongue_dx, L)
	assert isinstance(back_center_x, L)
	assert isinstance(front_center_x, L)
	assert isinstance(top_z, L)
	assert isinstance(center, P)
	assert isinstance(shaft_z, L)
    
	if isinstance(bottom_z, L):
	    self.bottom_z_l = bottom_z
	if isinstance(dx, L):
	    self.dx_l = dx
	if isinstance(top_z, L):
	    self.top_z_l = top_z
	self.center_p = center
	if isinstance(side_tongue_dx, L):
	    self.side_tongue_dx_l = side_tongue_dx
	if isinstance(back_center_x, L):
	    self.back_center_x_l = back_center_x
	if isinstance(front_center_x, L):
	    self.front_center_x_l = front_center_x
	if isinstance(side_tongue_bottom_z, L):
	    self.side_tongue_bottom_z_l = side_tongue_bottom_z
	if isinstance(side_tongue_bottom_z, L):
	    self.side_tongue_top_z_l = side_tongue_top_z
	if isinstance(shaft_z, L):
	    self.shaft_z_l = shaft_z

    def construct(self):
	""" Inner_Outer_Motor_Side: construct """

	# Grab some values from *self*:
	back_center_x = self.back_center_x_l
	bottom_z = self.bottom_z_l
	dx = self.dx_l
	dy = self.dy_l
	encoder_tongue_dz = self.encoder_tongue_dz_l
	front_center_x = self.front_center_x_l
	side_tongue_bottom_z = self.side_tongue_bottom_z_l
	side_tongue_dx = self.side_tongue_dx_l
	side_tongue_top_z = self.side_tongue_top_z_l
	shaft_z = self.shaft_z_l
	top_z = self.top_z_l
	center = self.center_p

	x_center = center.x
	y_center = center.y
	z_center = center.z

	#print("IOMSl.construct: name={0} x={1}:{2} y={3}:{4} z={5}:{6}".
	#  format(self, -dx/2, dx/2, y_center-dy/2, y_center+dy/2,
	#  bottom_z, top_z))

	# Grab some *Part*'s from the *motor_assembly*:
	center = self.center_p
	motor_assembly = self.up
	gm3_motor = motor_assembly.gm3_motor_
	encoder_mount = motor_assembly.encoder_mount_
	dual_slot_encoder = motor_assembly.dual_slot_encoder_
	dual_slot_encoder_pcb = dual_slot_encoder.dual_slot_encoder_pcb_

	x_center = center.x
	y_center = center.y
	z_center = center.z

	# Compute X positions:
	pcb_dx = dual_slot_encoder_pcb.dx_l
	encoder_tongue_dx = encoder_mount.tongue_dx_l
	x0 = x_center - dx / 2
	x1 = back_center_x - side_tongue_dx/2
	x2 = back_center_x + side_tongue_dx/2
	x3 = x_center - encoder_tongue_dx / 2
	x4 = x_center - pcb_dx / 2 - L(mm=1)
	x5 = x_center + pcb_dx / 2 + L(mm=1)
	x6 = x_center + encoder_tongue_dx / 2
	x7 = front_center_x - side_tongue_dx/2
	x8 = front_center_x + side_tongue_dx/2
	x9 = x_center + dx / 2

	# Cable slot X positions:
        cable_width = L(inch=0.550)
	#print("cable_width={0:i}".format(cable_width))
	x20 = x_center - L(mm=50.0) - cable_width/2
	x21 = x20 + cable_width
	x22 = x_center - L(mm=25.0) - cable_width/2
	x23 = x22 + cable_width
	x24 = x_center - L(mm=00.0) - cable_width/2
	x25 = x24 + cable_width
	x26 = x_center + L(mm=25.0) - cable_width/2
	x27 = x26 + cable_width
	x28 = x_center + L(mm=50.0) - cable_width/2
	x29 = x28 + cable_width

	y0 = y_center - dy / 2
	y1 = y_center + dy / 2
	#print("IOMS:y0={0} y_center={1} y1={2} dy={3}".
	#  format(y0, y_center, y1, dy))

	encoder_mount_z_center = encoder_mount.center_p.z
	encoder_mount_dz = encoder_mount.dz_l

	z_center = (bottom_z + top_z)/2
	z0 = bottom_z
	z1 = side_tongue_bottom_z - L(mm=0.2)
	z2 = z_center - encoder_tongue_dz/2
	z3 = encoder_mount_z_center - encoder_mount_dz/2 - L(mm=5)
	z4 = encoder_mount_z_center - encoder_mount_dz/2
	z5 = encoder_mount_z_center + encoder_mount_dz/2
	z6 = encoder_mount_z_center + encoder_mount_dz/2 + L(mm=5)
	z7 = z_center + encoder_tongue_dz/2
	z8 = side_tongue_top_z + L(mm=0.2)
	z9 = top_z - L(inch = 0.033)
	z10 = top_z
	#print("IOMS:z0={0} z1={1} z8={2} z10={3}".format(z0, z1, z8, z10))
	#print("IOMS:z3={0} z4={1} em_z_cent={2} z5={3} z6={4} em_dz={5}".
	#  format(z3, z4, encoder_mount_z_center, z5, z6, encoder_mount_dz))

	# Do the outer contour:
	zero = L()
	outer_contour = Contour(name = "Inner_Outer_Motor_Side")
	outer_contour.bend_append(P(x0, 0, z0), zero)
	outer_contour.bend_append(P(x0, 0, z10), zero)

	# First cable slot:
	outer_contour.bend_append(P(x20, 0, z10), zero)
	outer_contour.bend_append(P(x20, 0, z9), zero)
	outer_contour.bend_append(P(x21, 0, z9), zero)
	outer_contour.bend_append(P(x21, 0, z10), zero)

	# Second cable slot:
	outer_contour.bend_append(P(x22, 0, z10), zero)
	outer_contour.bend_append(P(x22, 0, z9), zero)
	outer_contour.bend_append(P(x23, 0, z9), zero)
	outer_contour.bend_append(P(x23, 0, z10), zero)

	# Third cable slot:
	outer_contour.bend_append(P(x24, 0, z10), zero)
	outer_contour.bend_append(P(x24, 0, z9), zero)
	outer_contour.bend_append(P(x25, 0, z9), zero)
	outer_contour.bend_append(P(x25, 0, z10), zero)

	# Fourth cable slot:
	outer_contour.bend_append(P(x26, 0, z10), zero)
	outer_contour.bend_append(P(x26, 0, z9), zero)
	outer_contour.bend_append(P(x27, 0, z9), zero)
	outer_contour.bend_append(P(x27, 0, z10), zero)

	# Fifth cable slot:
	outer_contour.bend_append(P(x28, 0, z10), zero)
	outer_contour.bend_append(P(x28, 0, z9), zero)
	outer_contour.bend_append(P(x29, 0, z9), zero)
	outer_contour.bend_append(P(x29, 0, z10), zero)

	outer_contour.bend_append(P(x9, 0, z10), zero)
	outer_contour.bend_append(P(x9, 0, z0), zero)

	# Do the *encoder_pocket* for the encoder mount:
	encoder_pocket = Contour(name = "IOMS:Encoder Pocket")
	encoder_pocket.bend_append(P(x3, y0, z4), zero)
	if True:
	    # Add in the bottom PCB gap:
	    encoder_pocket.bend_append(P(x4, y0, z4), zero)
	    encoder_pocket.bend_append(P(x4, y0, z3), zero)
	    encoder_pocket.bend_append(P(x5, y0, z3), zero)
	    encoder_pocket.bend_append(P(x5, y0, z4), zero)
	encoder_pocket.bend_append(P(x6, y0, z4), zero)
	encoder_pocket.bend_append(P(x6, y0, z5), zero)
	if True:
	    # Add in the top PCB gap:
	    encoder_pocket.bend_append(P(x5, y0, z5), zero)
	    encoder_pocket.bend_append(P(x5, y0, z6), zero)
	    encoder_pocket.bend_append(P(x4, y0, z6), zero)
	    encoder_pocket.bend_append(P(x4, y0, z5), zero)
	encoder_pocket.bend_append(P(x3, y0, z5), zero)
	#print("encoder_pocket={0}".format(encoder_pocket))

	# Do the *back_pocket*:
	back_pocket = Contour(name = "IOMS:Back Pocket")
	back_pocket.bend_append(P(x1, 0, z1), zero)
	back_pocket.bend_append(P(x2, 0, z1), zero)
	back_pocket.bend_append(P(x2, 0, z8), zero)
	back_pocket.bend_append(P(x1, 0, z8), zero)
	#print("back_pocket={0}".format(back_pocket))

	# Do the *front_pocket*:
	front_pocket = Contour(name = "IOMS:Front_Pocket")
	front_pocket.bend_append(P(x7, 0, z1), zero)
	front_pocket.bend_append(P(x8, 0, z1), zero)
	front_pocket.bend_append(P(x8, 0, z8), zero)
	front_pocket.bend_append(P(x7, 0, z8), zero)
	#print("front_pocket={0}".format(front_pocket))

	# Do the extrusion:
	self.extrude(comment = "Inner_Outer_Motor_Side",
	  material = Material("Wood", "Plywood"),
	  color = Color("cyan"),
	  outer_contour = outer_contour,
	  inner_contours = [back_pocket, front_pocket, encoder_pocket],
	  start = P(x_center, y0, shaft_z),
	  end =   P(x_center, y1, shaft_z))
	
	# Do the encoder hub hole:
	self.hole(comment = "Hub Hole 1",
	  diameter = gm3_motor.encoder_hub_diameter_l + L(mm=1),
	  start = P(x_center, y1, shaft_z),
	  end =   P(x_center, y0, shaft_z),
	  flags = "t")

	# Do the motor mount holes:
	pin_diameter = gm3_motor.pin_diameter_l
	pin_dx = gm3_motor.pin_dx_l
	back_holes_diameter = gm3_motor.back_holes_diameter_l + L(mm=1)
	back_holes_dx = gm3_motor.back_holes_dx_l
	back_holes_pitch_dz = gm3_motor.back_holes_pitch_dz_l
	strap_pocket_dx = gm3_motor.strap_pocket_dx_l
	strap_pocket_dz = gm3_motor.strap_pocket_dz_l
	strap_hole_dx = gm3_motor.strap_hole_dx_l
	self.hole(comment = "Mount Hole 1",
	  diameter = back_holes_diameter,
	  start = P(x_center + back_holes_dx, y1,
	  shaft_z + back_holes_pitch_dz/2),
	  end =   P(x_center + back_holes_dx, y0,
	  shaft_z + back_holes_pitch_dz/2),
	  flags = "t")
	self.hole(comment = "Mount Hole 2",
	  diameter = back_holes_diameter,
	  start = P(x_center + back_holes_dx, y1,
	  shaft_z - back_holes_pitch_dz/2),
	  end =   P(x_center + back_holes_dx, y0,
	  shaft_z - back_holes_pitch_dz/2),
	  flags = "t")
	self.hole(comment = "Alignment Pin",
          diameter = pin_diameter,
	  start = P(x_center + pin_dx, y1, shaft_z),
	  end =   P(x_center + pin_dx, y0, shaft_z),
	  flags = "t")
	self.simple_pocket(comment = "Strap Hole Pocket",
	  bottom_corner =
	  P(x_center + strap_hole_dx - strap_pocket_dx/2, y0 - L(mm=1),
	  shaft_z - strap_pocket_dz/2),
	  top_corner =
	  P(x_center + strap_hole_dx + strap_pocket_dx/2, y1 + L(mm=1),
	  shaft_z + strap_pocket_dx/2),
	  pocket_top = "w")


class Right_Motor_Assembly(Part):
    def __init__(self, up):
	""" Right_Motor_Assembly: initialize """
	self.thin_l = up.thin_l
	self.thick_l = up.thick_l
	Part.__init__(self, up)
	self.back_motor_side_ = Front_Back_Motor_Side(self)
	self.back_screw_ = Fastener(self, "RA:b")
	self.dual_slot_encoder_ = Dual_Slot_Encoder(self)
	self.encoder_disk_ = Encoder_Disk(self)
	self.encoder_mount_ = Encoder_Mount(self)
	self.encoder_ne_screw_ = Fastener(self, "RA:ene")
	self.encoder_nw_screw_ = Fastener(self, "RA:enw")
	self.encoder_se_screw_ = Fastener(self, "RA:ese")
	self.encoder_sw_screw_ = Fastener(self, "RA:esw")
	self.front_motor_side_ = Front_Back_Motor_Side(self)
	self.front_screw_ = Fastener(self, "RA:f")
	self.gm3_motor_ = GM3_Motor(self)
	self.gm3_wheel_ = GM3_Wheel(self)
	self.inner_be_screw_ = Fastener(self, "RA:ibe")
	self.inner_bw_screw_ = Fastener(self, "RA:ibw")
	self.inner_motor_side_ = Inner_Outer_Motor_Side(self)
	self.inner_te_screw_ = Fastener(self, "RA:ite")
	self.inner_tw_screw_ = Fastener(self, "RA:itw")
	self.outer_be_screw_ = Fastener(self, "RA:obe")
	self.outer_bw_screw_ = Fastener(self, "RA:obw")
	self.outer_motor_side_ = Inner_Outer_Motor_Side(self)
	self.outer_te_screw_ = Fastener(self, "RA:ote")
	self.outer_tw_screw_ = Fastener(self, "RA:otw")
	self.top_back_screw_ = Fastener(self, "RA:tb")
	self.top_center_p = P()
	self.top_front_screw_ = Fastener(self, "RA:tf")
	self.top_screw_pitch_l = L()
	self.z0_l = L()
	
    def configure(self, top_center = P(), top_screw_pitch = L(), z0 = L()):
	""" *Right_Motor_Assembly*: Configuration. """
	assert isinstance(top_center, P)
	assert isinstance(top_screw_pitch, L)
	assert isinstance(z0, L)
	self.top_center_p = top_center
	self.top_screw_pitch_l = top_screw_pitch
	self.z0_l = z0
	
    def construct(self):
	""" *Right_Motor_Assembly*: construct """

	# Grab the various *Part*' from *self*:
	back_motor_side = self.back_motor_side_
	back_screw = self.back_screw_
	dual_slot_encoder = self.dual_slot_encoder_
	dual_slot_encoder_pcb = dual_slot_encoder.dual_slot_encoder_pcb_
	encoder_disk = self.encoder_disk_
	encoder_mount = self.encoder_mount_
	encoder_se_screw = self.encoder_se_screw_
	encoder_sw_screw = self.encoder_sw_screw_
	encoder_ne_screw = self.encoder_ne_screw_
	encoder_nw_screw = self.encoder_nw_screw_
	front_motor_side = self.front_motor_side_
	front_screw = self.front_screw_
	gm3_motor = self.gm3_motor_
	gm3_wheel = self.gm3_wheel_
	inner_be_screw = self.inner_be_screw_
	inner_bw_screw = self.inner_bw_screw_
	inner_motor_side = self.inner_motor_side_
	inner_te_screw = self.inner_te_screw_
	inner_tw_screw = self.inner_tw_screw_
	outer_be_screw = self.outer_be_screw_
	outer_bw_screw = self.outer_bw_screw_
	outer_motor_side = self.outer_motor_side_
	outer_te_screw = self.outer_te_screw_
	outer_tw_screw = self.outer_tw_screw_
	thick = self.thick_l
	thin = self.thin_l
	top_back_screw = self.top_back_screw_
	top_center = self.top_center_p
	top_front_screw = self.top_front_screw_
	top_screw_pitch = self.top_screw_pitch_l
	z0 = self.z0_l

	slot_interrupter1 = dual_slot_encoder.slot_interrupter1_

	# Compute interesting X locations in ascending order:
	zero = L()
	top_screw_pitch = L(inch=5)
	dx = top_screw_pitch + 2 * thin + thick
	x0 = top_center.x - dx/2	# West edge of motor assembly
	x1 = x0 + thin			# West edge of back motor side
	x2 = x1 + thick/2		# Center of back motor side
	x3 = x2 + thick/2		# East edge of back motor side
	x4 = x3 + thick			# West tongue edge
	x5 = x4 + L(mm=2)		# End of screw mounts
	x6 = top_center.x - dual_slot_encoder_pcb.dx/2 # West end of encoder pcb
	x7 = top_center.x		# Center
        x8 = top_center.x + dual_slot_encoder_pcb.dx/2 # East end of encoder pcb
	x14 = top_center.x + dx/2	# East edge of motor assembly
	x13 = x14 - thin		# East edge of front motor side
	x12 = x13 - thick/2		# Center of front motor side
	x11 = x12 - thick/2		# West edge of front motor side
	x10 = x11 - thick		# East tongue edge
	x9 = x10 - thin			# End of screw mounts
	assert x12 - x2 == top_screw_pitch, \
	  "top_screw_pitch missed {0} != {1}".format(x12 - x2, top_screw_pitch)

	# Compute interesting Y locations in ascending order:
	one = L(mm=1)
	y11 = top_center.y				# Center of assembly
	y4 = y11 - L(mm=20)				# Center inner_side
	y5 = y4 + thin/2				# North inner side
	y3 = y4 - thin/2				# South inner side
	# The motor used to be on the south side of *inner_side*.  Now it
	# is on the north side.  Thus *y1*, *y2*, *y3* are now collectively
	# greater than *y3*, *y4*, and *y5*.
	y1 = y5 + gm3_motor.dy_l/2			# Motor center
	y2 = y1 + gm3_motor.dy_l/2                      # Disk side Motor edge
	y0 = y5 - gm3_motor.wheel_shaft_dy_l		# Wheel center

	y16 = y11 + L(mm=20) + thin/2			# North outer side
	y15 = y16 - thin/2				# Outer side center
	y14 = y15 - thin/2				# South outer side
	y13 = y16 - L(mm=10)				# End mount screws

	y7 = y3 + L(mm=10)				# End of mount screws
	y8 = y2 + gm3_motor.encoder_shaft_dy_l		# Encoder shaft edge
	y9 = y8 + encoder_disk.dy_l/2			# Encoder disk center

	# Compute the PCB values:
	y10 = y9 - dual_slot_encoder.disk_offset_dy_l	# PCB center
	y6 = y10 - dual_slot_encoder_pcb.dy/2		# East PCB edge
	y12 = y10 + dual_slot_encoder_pcb.dy/2		# West PCB edge

	#print("TLA:y0={0} y1={1} y2={2} y3={3} y4={4} y5={5} y6={6} y7={7}".
	#  format(y0, y1, y2, y3, y4, y5, y6, y7))
	#print("TLA:y8={0} y9={1} y10={2} y11={3} y12={4} y13={5} y14={6}".
	#  format(y8, y9, y10, y11, y12, y13, y14))
	#print("TLA:y15={0} y16={1}".format(y15, y16))
	
	# Compute interesting Z locations in acending order:
	slot_dz = slot_interrupter1.dz_l.absolute()
	slot_gap_dz = slot_interrupter1.slot_dz_l.absolute()
	slot_base_dz = slot_dz - slot_gap_dz
	# z0 comes from self.z0_l			# The floor
	z3 = z0 + gm3_wheel.diameter_l/2		# Shaft center-line
	z1 = z3 - gm3_motor.dz / 2 - L(mm=5)		# Outer side edge
	z2 = z1 + L(mm=10)				# Bottom tongue edge
	# Note *z3* defined before *z1* and *z2*:
	z4 = z3 + encoder_disk.diameter_l / 2		# Top encoder disk edge
	z5 = z4 + L(mm=0.5)				# Top of disk gap
	z6 = z5 + slot_base_dz				# PCB bottom
	z7 = z6 + dual_slot_encoder_pcb.dz_l/2		# PCB middle
	z8 = z7 + encoder_mount.dz_l/2			# Encoder mount center
	z9 = z8 + encoder_mount.dz_l/2			# Encoder mount top
	z12 = top_center.z				# Top edge of assembly
	z11 = z12 - L(mm=7)				# End of screw mount
	z10 = z11 - L(mm=10)				# Top tongue edge
	z13 = z12 + L(mm=5)				# Top of base

	# Configure the various *Part*'s:
	dual_slot_encoder.configure(center = P(x7, y10, z6))
	encoder_disk.configure(center = P(x7, y9, z3))
	gm3_motor.configure(shaft_center = P(x7, y1, z3))
	gm3_wheel.configure(center = P(x7, y0, z3))
	inner_motor_side.configure(dx = dx,
	  side_tongue_dx = x3 - x1, back_center_x = x2, front_center_x = x12,
	  center = P(x7, y4, 0), bottom_z = z1, top_z = z12,
	  side_tongue_bottom_z = z2, side_tongue_top_z = z10, shaft_z = z3)
	outer_motor_side.configure(dx = dx,
	  side_tongue_dx = x3 - x1, back_center_x = x2, front_center_x = x12,
	  center = P(x7, y15, 0), bottom_z = z1, top_z = z12,
	  side_tongue_bottom_z = z2, side_tongue_top_z = z10, shaft_z = z3)
	#print("RMA:z2={0} z10={1}".format(z2, z10))

	#print("RMA:x0={0} x11={1} x13={2} x14={3}".format(x0, x11, x13, x14))
	#print("RMA:z1={0} z12={1}".format(z1, z12))
	back_motor_side.configure(dx = x3 - x1,
	  dy = y16 - y3, inside_dy = y14 - y5,
	  dz = z12 - z1, z_tongue_bottom = z2, z_tongue_top = z10,
	  center = P((x1 + x3)/2, (y16 + y3)/2, (z1 + z12)/2))
	front_motor_side.configure(dx = x13 - x11,
	  dy = y16 - y3, inside_dy = y14 - y5,
	  dz = z12 - z1, z_tongue_bottom = z2, z_tongue_top = z10,
	  center = P((x13 + x11)/2, (y16 + y3)/2, (z1 + z12)/2))

	#inner_outer_dx = inner_motor_side.dx_l
	#inner_outer_pitch_dy = y15 - y4
	#inner_outer_dy = inner_motor_side.dy_l/2 + \
	#  inner_outer_pitch_dy +  outer_motor_side.dy_l/2
	#inner_outer_center_y = (y4 + y15)/2
	encoder_mount.configure(dx = x11 - x3, tongue_dx = x10 - x4,
	  dy = y16 - y3, center = P(x7, y11, z8), dz = thick)

	# Configure the *inner_motor_side* fasteners:
	screw_flags = "M3x.05"
	screw_flags = "#4-40"
	inner_be_screw.configure(comment = "Inner BE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x12, y3, z1 + L(mm=5)),
	  end = P(x12, y7, z1 + L(mm=5)),
	  flags = screw_flags)
	#print("RMA:y2={0} y5={1}".format(y2, y5))
	inner_be_screw.drill(part = inner_motor_side, select = "close")
	inner_be_screw.nut_ledge(part = front_motor_side, flags="E")

	inner_bw_screw.configure(comment = "Inner BW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x2, y3, z1 + L(mm=5)),
	  end = P(x2, y7, z1 + L(mm=5)),
	  flags = screw_flags)
	inner_bw_screw.drill(part = inner_motor_side, select = "close")
	inner_bw_screw.nut_ledge(part = back_motor_side, flags="W")

	inner_te_screw.configure(comment = "Inner TE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x12, y3, z12 - L(mm=5)),
	  end = P(x12, y7, z12 - L(mm=5)),
	  flags = screw_flags)
	inner_te_screw.drill(part = inner_motor_side, select = "close")
	inner_te_screw.nut_ledge(part = front_motor_side, flags="E")

	inner_tw_screw.configure(comment = "Inner TW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x2, y3, z12 - L(mm=5)),
	  end = P(x2, y7, z12 - L(mm=5)),
	  flags = screw_flags)
	inner_tw_screw.drill(part = inner_motor_side, select = "close")
	inner_tw_screw.nut_ledge(part = back_motor_side, flags="W")

	# Configure the *outer_motor_side* fasteners:
	outer_be_screw.configure(comment = "Outer BE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x12, y16, z1 + L(mm=5)),
	  end = P(x12, y13, z1 + L(mm=5)),
	  flags = screw_flags)
	outer_be_screw.drill(part = outer_motor_side, select = "close")
	outer_be_screw.nut_ledge(part = front_motor_side, flags="E")

	outer_bw_screw.configure(comment = "Outer BW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x2, y16, z1 + L(mm=5)),
	  end = P(x2, y13, z1 + L(mm=5)),
	  flags = screw_flags)
	outer_bw_screw.drill(part = outer_motor_side, select = "close")
	outer_bw_screw.nut_ledge(part = back_motor_side, flags="W")

	outer_te_screw.configure(comment = "Outer TE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x12, y16, z12 - L(mm=5)),
	  end = P(x12, y13, z12 - L(mm=5)),
	  flags = screw_flags)
	outer_te_screw.drill(part = outer_motor_side, select = "close")
	outer_te_screw.nut_ledge(part = front_motor_side, flags="E")

	outer_tw_screw.configure(comment = "Outer TW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x2, y16, z12 - L(mm=5)),
	  end = P(x2, y13, z12 - L(mm=5)),
	  flags = screw_flags)
	outer_tw_screw.drill(part = outer_motor_side, select = "close")
	outer_tw_screw.nut_ledge(part = back_motor_side, flags="W")

	# Configure the *outer_motor_side* fasteners:
	back_screw.configure(comment = "Back Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x1, y11, z8),
	  end = P(x5, y11, z8),
	  flags = screw_flags)
	back_screw.drill(part = back_motor_side, select = "close")
	back_screw.nut_ledge(part = encoder_mount, flags="T")

	front_screw.configure(comment = "Front Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x13, y11, z8),
	  end = P(x9, y11, z8),
	  flags = screw_flags)
	front_screw.drill(part = front_motor_side, select = "close")
	front_screw.nut_ledge(part = encoder_mount, flags="T")

	# Configure the top screw fasteners:
	top_back_screw.configure(comment = "Top Back Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x2, y11, z13 + L(mm=10)),
	  end = P(x2, y11, z11),
	  flags = screw_flags)
	top_back_screw.nut_ledge(part = back_motor_side, flags="W")

	top_front_screw.configure(comment = "Top Front Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x12, y11, z13 + L(mm=10)),
	  end = P(x12, y11, z11),
	  flags = screw_flags)
	top_front_screw.nut_ledge(part = front_motor_side, flags="E")

	encoder_se_screw.configure(comment = "Encoder SE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x8 - L(mm=5), y6 + L(mm=5), z9),
	  end =   P(x8 - L(mm=5), y6 + L(mm=5), z6),
	  flags = screw_flags)
	encoder_se_screw.drill(part = encoder_mount, select = "close")

	encoder_sw_screw.configure(comment = "Encoder SW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x6 + L(mm=5), y6 + L(mm=5), z9),
	  end =   P(x6 + L(mm=5), y6 + L(mm=5), z6),
	  flags = screw_flags)
	encoder_sw_screw.drill(part = encoder_mount, select = "close")

	encoder_ne_screw.configure(comment = "Encoder NE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x8 - L(mm=5), y12 - L(mm=5), z9),
	  end =   P(x8 - L(mm=5), y12 - L(mm=5), z6),
	  flags = screw_flags)
	encoder_ne_screw.drill(part = encoder_mount, select = "close")

	encoder_nw_screw.configure(comment = "Encoder NW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x6 + L(mm=5), y12 - L(mm=5), z9),
	  end =   P(x6 + L(mm=5), y12 - L(mm=5), z6),
	  flags = screw_flags)
	encoder_nw_screw.drill(part = encoder_mount, select = "close")

	# Deal with visibility here:
	gm3_wheel.invisible_set()
	#encoder_mount.invisible_set()
	#outer_motor_side.invisible_set()
	#inner_motor_side.invisible_set()
	#back_motor_side.invisible_set()
	#front_motor_side.invisible_set()

	# Deal .dxf files here:
	encoder_mount.dxf_write(center = P(0, y11, z8))
	inner_motor_side.dxf_write(
	  center = P(0, y4, z8), plane_normal = P(0, L(mm=1), 0))
	front_motor_side.dxf_write(
	  center = P(x12, y11, z8), plane_normal = P(L(mm=1), 0, 0))

class ROSBot(Part):
    def __init__(self, up):
	""" *ROSBot*: Initialization. """
	self.thick_l = L(mm=5)
	self.thin_l = L(mm=3)
	Part.__init__(self, up)
	self.base_ = base = Base(self)
	self.right_motor_assembly_ = Right_Motor_Assembly(self)
	self.castor_assembly_ = Castor_Assembly(self)

    def construct(self):
	""" *ROSBot*: Construction. """
	right_motor_assembly = self.right_motor_assembly_
	castor_assembly = self.castor_assembly_
	thick = self.thick_l
	thin = self.thin_l
	base = self.base_

	# Grab some values from *base*:
	base_z0 = base.z0_l
	base_dz = base.dz_l

	z0 = base_z0 - L(mm=100)
	dx = L(inch=2) + thin
	dy = L(inch=3) + thin
	castor_x = L(inch=-4)
	castor_assembly.configure(
	  corner1 = P(castor_x - dx/2, -dy/2, z0),
	  corner2 = P(castor_x + dx/2,  dy/2, base_z0 - base_dz))

	top_center = P(0, 0, base_z0 - thick)
	top_screw_pitch = L(inch=5)
	right_motor_assembly.configure(
	  top_center = P(L(inch=2.5), -L(inch=5), base_z0 - thick),
	  top_screw_pitch = L(inch=5),
	  z0 = z0)

	#base.invisible_set()

class Slot_Interrupter(Part):
    def __init__(self, up):
	""" Slot_Interrupter: initialize """
	Part.__init__(self, up)
	self.dx_l = L(mm=2.6)
	self.dy_l = L(mm=4.5)
	self.dz_l = L(mm=-4.5)
	self.slot_dy_l = L(mm=2)
	self.slot_dz_l = L(mm=-2.1)
	self.center_p = P()

    def configure(self, center = None):
	""" Slot_Interrupter: configure """
	# Check argument types:
	none_type = type(None)
	assert isinstance(center, P)

	# Change any configuration values:
	self.center_p = center

    def construct(self):
	""" Slot_Interrupter: construct """
	zero = L()

	dx = self.dx_l
	dy = self.dy_l
	dz = self.dz_l
	slot_dy = self.slot_dy_l
	slot_dz = self.slot_dz_l
	center = self.center_p

	x_center = center.x
	y_center = center.y
	z_center = center.z

	x0 = x_center - dx/2
	x1 = x_center + dx/2

	y0 = y_center - dy/2
	y1 = y_center - slot_dy/2
	y2 = y_center + slot_dy/2
	y3 = y_center + dy/2

	z0 = z_center
	z1 = z_center + dz - slot_dz
	z2 = z_center + dz
	#print("Slot_Interrupter:z0/z1/z2={0}/{1}/{2}".format(z0, z1, z2))

	z = L()
	self.block(comment = "Slot Encoder Block",
	  material = Material("plastic", "plastic"),
	  color = Color("blue"),
	  corner1 = P(x0, y0, z0),
	  corner2 = P(x1, y3, z2),
	  top = "t")

	pocket_top = "t"
	if dz < zero:
	    pocket_top = "b"

	self.simple_pocket(comment = "Remove Slot",
	  bottom_corner = P(x0, y1, z1),
	  top_corner = P(x1, y2, z2),
	  pocket_top = pocket_top)

if __name__ == "__main__":
	ezcad = EZCAD3(0, adjust = L(mm = -0.15))
	test = ROSBot(None)
	test.process(ezcad)
