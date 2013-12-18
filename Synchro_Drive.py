#!/usr/bin/env python

from EZCAD3 import *

# Assemblies:

class Synchro_Drive(Part):

    def __init__(self, up):
	Part.__init__(self, up)
	self.wheel_assembly_ = Wheel_Assembly(self)
	self.motor_assembly_ = Motor_Assembly(self)

    def construct(self):
	pass

class Motor_Assembly(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.motor_base_ = Motor_Base(self)

    def construct(self):
	motor_base = self.motor_base_
	synchro_drive = self.up
	wheel_assembly = synchro_drive.wheel_assembly_
	turn_table = wheel_assembly.turn_table_
	gear_box = wheel_assembly.gear_box_

	self.base_sz_l = base_sz = gear_box.nz_l + turn_table.height_l

class Wheel_Assembly(Part):
    def __init__(self, up):
	Part.__init__(self, up)

	# List all of the sub-*Part*'s that make up a wheel assembly:
	#  teeth_count = 19, belt_class = "XL",
        self.twist_timing_pulley_ = Timing_Pulley(self,
	  teeth_count = 44, belt_class = "MXL",
	  belt_width = L(inch = "1/8") + L(inch = .005),
	  lip_extra = L(mm = 1.5), lip_width = L(mm = 5),
	  bearing_diameter = L(inch = "7/8"), bearing_width = L(inch = "9/32"),
	  shaft_diameter = L(inch = "1/2"))
	self.bearing_ = Bearing(self)
	self.bevel_gear_ = \
	  Bevel_Gear(self, part_name = "A 1M 4-Y16016", color = Color("red"))
	self.gear_box_ = Gear_Box(self)
	self.gear_box_cover_ = Gear_Box_Cover(self)
	self.gear_box_shelf_ = Gear_Box_Shelf(self)
	self.gear_box_top_ = Gear_Box_Top(self)
	self.horizontal_shaft_ = Horizontal_Shaft(self)
	self.non_wheel_side_gear_box_side_ = \
	  Gear_Box_Side(self, False)
	self.shim_ = Shim(self)
	self.turn_table_ = Turn_Table(self)
	self.vertical_shaft_ = Vertical_Shaft(self)
	self.wheel_ = Wheel(self)
	self.wheel_side_gear_box_side_ = \
	  Gear_Box_Side(self, True)
 
    def construct(self):
	""" *Wheel_Assembly*:  """
	# The origin of this assembly is at the logical intersection
	# of the vertical shaft and the horizontal shaft.  The horizontal
	# axis is aligned with the Y axis; this causes the wheel to be
	# oriented in the roll in the X direction.

	# *synchro_drive is the up-level parent of *self*:
	synchro_drive = self.up

	# Grab some values from *bearing*:
	bearing = self.bearing_
	bearing_width = bearing.width_l

	# Grab some values from *bevel_gear*:
	bevel_gear = self.bevel_gear_
	bevel_gear_outside_diameter = bevel_gear.outside_diameter_l
	bevel_gear_major_diameter = bevel_gear.major_diameter_l

	# Grep some some values from *gear_box_side":
	gear_box_side = self.wheel_side_gear_box_side_
	gear_box_side_bearing_lip_dy = \
	  gear_box_side.bearing_lip_dy_l

	# Grab some values from *motor_base* (from *motor_assembly*):
	motor_assembly = synchro_drive.motor_assembly_
	motor_base = motor_assembly.motor_base_
	motor_base_dz = motor_base.dz_l

	# Grab some values from *shim*:
	shim = self.shim_
	shim_dz = shim.dz_l

	# Grab some values from *turn_table*:
	turn_table = self.turn_table_
	turn_table.no_automatic_place()
	turn_table_height = turn_table.height_l

	# Grab some values from *vertical_shaft*:
	vertical_shaft = self.vertical_shaft_
	vertical_shaft_diameter = vertical_shaft.diameter_l

	zero = L()
	y0 = zero
	y1 = y0 + bevel_gear_major_diameter
	y2 = y1 + shim_dz
	y3 = y2 + bearing_width / 2
	y4 = y2 + bearing_width
	y5 = y4 + gear_box_side_bearing_lip_dy

	# Define some constants:
	gear_box = self.gear_box_
	self.gear_box_shelf_sz_l = \
	  gear_box_shelf_sz = bevel_gear.nz_l + shim_dz

	# Update the timing pulley:
	self.twist_timing_pulley_.update(
	  lip_width = turn_table_height + L(mm = 1.000) + motor_base_dz,
	  shaft_diameter = vertical_shaft_diameter + L(inch = "1/4"))

	# Place various items into the wheel assembly:
	one = L(mm = 1.0)
	x_axis = P(one, zero, zero)
	degrees_90 = Angle(deg = 90)

	# Place the 4 bearings:
	self.place(bearing, name = "Wheel Side Bearing", 
	  axis = x_axis, rotate = degrees_90, translate = P(zero, -y3, zero))
	self.place(bearing, name = "Non Wheel Side Bearing", 
	  axis = x_axis, rotate = degrees_90, translate = P(zero,  y3, zero))
	self.place(bearing, name = "Lower Bearing", translate =
	  P(zero, zero, gear_box_shelf_sz + bearing_width / 2))

	# Place the horizontal bevel gear.  (The vertical one auto places):
	self.place(bevel_gear, name = "Horizontal_Bevel_Gear", 
	  axis = x_axis, rotate = degrees_90)

	# Place the shims on next to the bevel gears:
	self.place(shim, name = "Vertical Shaft Shim",
	  translate = P(zero, zero, y1))
	self.place(shim, name = "Horzontal Shaft Shim",
	  axis = x_axis, rotate = degrees_90,
	  translate = P(zero, -y1, zero))

	# Drop in the turn table:
	self.place(self.turn_table_, name = "XTurn_Table",
	  translate = P(zero, zero, gear_box.nz_l))

	# Drop in the two pullies:
	self.place(self.twist_timing_pulley_, name = "Twist Turn Pulley",
	  translate = P(zero, zero, gear_box.nz_l))

# Various parts:

class Bearing(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	# R6ZZ Shield Ball Bearing from VBX.Com:
	self.bore_l = bore = L(inch = "3/8")
	self.diameter_l = diameter = L(inch = "7/8")
	self.width_l = width = L(inch = "9/32")

	self.no_automatic_place()

	zero = L()
	self.cylinder(comment = "Bearing Cylinder",
	  material = Material("steel", ""),
	  color = Color("blue"),
	  diameter = diameter, 
	  start = P(zero, zero, -width / 2),
	  end = P(zero, zero, width / 2),
	  sides = 30)
	self.hole(comment = "Bore Hole",
	  diameter = bore,
	  start = P(zero, zero, -width / 2),
	  end = P(zero, zero, width / 2),
	  flags = "t")

class Bevel_Gear(Part):
    def __init__(self, up, part_name, color):
	Part.__init__(self, up)

	# Check argument types:
	assert isinstance(part_name, str)
	assert isinstance(color, Color)

	# Lookup up the *bevel_gear*:
	bevel_gear_table = Bevel_Gear_Table()
	table  = bevel_gear_table.table
	assert part_name in table
	bevel_gear_data = table[part_name]

	# Load up *self*:
	self.part_name_s = part_name
	self.bevel_gear_data_o = bevel_gear_data
	self.color = color
	self.pitch_i = bevel_gear_data.pitch
	self.teeth_i = bevel_gear_data.teeth
	self.pitch_diameter_l = bevel_gear_data.pitch_diameter
	self.outside_diameter_l = bevel_gear_data.outside_diameter
	self.bore_l = bevel_gear_data.bore
	self.height_l = bevel_gear_data.height
	self.hub_diameter_l = bevel_gear_data.hub_diameter
	self.hub_heigh_l = bevel_gear_data.hub_height
	self.major_diameter_l = bevel_gear_data.major_diameter

    def construct(self):
	# Grab some values from *bevel_gear*:
	bevel_gear_data = self.bevel_gear_data_o
	material = bevel_gear_data.material
	outside_diameter = bevel_gear_data.outside_diameter
	bore = bevel_gear_data.bore
	height = bevel_gear_data.height
	hub_diameter = bevel_gear_data.hub_diameter
	hub_height = bevel_gear_data.hub_height
	major_diameter = bevel_gear_data.major_diameter
	color = self.color
	zero = L()

	self.nz_l = major_diameter

	self.cylinder(comment = "Bevel Gear Teeth",
	  material = material,
	  color = Color("red"),
	  diameter = outside_diameter, 
	  start = P(zero, zero, major_diameter - height),
	  end = P(zero, zero, major_diameter - hub_height),
	  sides = 30)
	self.cylinder(comment = "Bevel Gear Hub",
	  diameter = hub_diameter, 
	  start = P(zero, zero, major_diameter - hub_height),
	  end = P(zero, zero, major_diameter),
	  sides = 30)
	self.hole(comment = "Bore Hole",
	  diameter = bore,
	  start = P(zero, zero, major_diameter),
	  end = P(zero, zero, major_diameter - height),
	  flags = "t")
	
class Bevel_Gear_Data:
    # From www.sdp-si.com:
    #
    # Mat   Teeth  PD   OD  Bore  Height  Hub Diam  Hub Height  Maj Diam
    #
    # 48 Pitch:
    # Actel   15  .313  .350 .125  7/32    1/4       1/8         5/16
    # Nylon   18  .375  .404 .125  9/32    21/64     1/4         13/32
    # Nylon   24  .500  .529 .187  3/8     3/8       1/4         17/3
    # 32 Pitch:
    # Nylon   24  .750  .800 .187  13/32   1/2       3/16        11/16
    # Nylon   16  .500  .550 .187  11/32   13/32     3/16        1/2
    # 24 Pitch:
    # Nylon   24 1.000 1.045 .253 9/16     5/8       3/16        7/8
    # Nylon   30 1.250 1.310 .250 37/64    5/8       5/16        1-1/32
    # Nylon   36 1.500 1.560 .312 39/64    11/16     5/16        1-3/16
    # 16 Pitch:
    # Nylon   16 1.000 1.090 .375 3/4      3/4       7/16        1-1/16
    # ~Nylon  32 2.000 2.088 .500 7/8      1-1/4     3/8         1-9/16

    # Part number: 1M 4-ppttt  Where pp is pitch and ttt is teeth
    #
    #  Part Number   1-24  25-99  100-249  250-999  1000+
    # A 1M 4-Y48015: 2.85   1.73    1.31     1.11   0.87
    # A 1M 4-Y48018: 2.85   1.73    1.31     1.11   0.87
    # A 1M 4-Y48024: 2.85   1.73    1.31     1.11   0.86
    #
    # A 1M 4-Y32024: 2.96   1.83    1.39     1.19   0.94
    # A 1M 4-Y32016: 2.92   1.80    1.36     1.16   0.92
    #
    # A 1M 4-Y24024: 3.39   2.31    1.87     1.67   1.45
    # A 1M 4-Y24030: 3.41   2.34    1.89     1.70   1.47
    # A 1M 4-Y24036: 3.46   2.37    1.92     1.72   1.50
    #
    # A 1M 4-Y16016: 3.41   2.34    1.89     1.70   1.47
    # A 1M 4-Y16032: 5.63   4.28    3.67     3.43   3.17

    def __init__(self, part_name = None, material = None, pitch = None,
      teeth = None, pitch_diameter = None, outside_diameter = None,
      bore = None, height = None, hub_diameter = None, hub_height = None,
      major_diameter = None, cost1_24 = None, cost25_99 = None,
      cost100_249 = None, cost250_999 = None, cost_1000_plus = None):

        # Check argument types:
	assert isinstance(part_name, str)
	assert isinstance(material, Material)
	assert isinstance(pitch, int)
	assert isinstance(teeth, int)
	assert isinstance(pitch_diameter, L)
	assert isinstance(outside_diameter, L)
	assert isinstance(bore, L)
	assert isinstance(height, L)
	assert isinstance(hub_diameter, L)
	assert isinstance(hub_height, L)
	assert isinstance(major_diameter, L)
      	assert isinstance(cost1_24, float)
	assert isinstance(cost25_99, float)
	assert isinstance(cost100_249, float)
	assert isinstance(cost250_999, float)
	assert isinstance(cost_1000_plus, float)

	# Load up *self*:
	self.part_name = part_name
	self.material = material
	self.pitch = pitch
	self.teeth = teeth
	self.pitch_diameter = pitch_diameter
	self.outside_diameter = outside_diameter
	self.bore = bore
	self.height = height
	self.hub_diameter = hub_diameter
	self.hub_height = hub_height
	self.major_diameter = major_diameter
	self.cost1_24 = cost1_24 
	self.cost25_99 = cost25_99
	self.cost100_249 = cost100_249
	self.cost250_999 = cost250_999
	self.cost_1000_plus = cost_1000_plus

class Gear_Box(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	# Grab some values *wheel_assembly*:
	wheel_assembly = self.up
	bearing = wheel_assembly.bearing_
	bevel_gear = wheel_assembly.bevel_gear_
	shim = wheel_assembly.shim_
	gear_box_side = wheel_assembly.wheel_side_gear_box_side_

	# Figure out *y5* the distance from the center to the outside Y wall:
	zero = L()
	y0 = zero
	y1 = y0 + bevel_gear.major_diameter_l
	y2 = y1 + shim.dz_l
	y3 = y2 + bearing.width_l / 2
	y4 = y2 + bearing.width_l
	y5 = y4 + gear_box_side.bearing_lip_dy_l

	# Specify *dx*, *dy*, *nz*, and *sz*:
	dx = L(mm = 40.00)
	dy  = y5 * 2
	self.nz_l = nz = L(mm = 75.00)
	self.sz_l = sz = L(mm = -10.00)

	# Now specify the virtual box for the gear box:
	self.virtual_box(comment = "Virtual box for gear box",
	  corner1 = P(-dx / 2, -dy / 2, sz),
	  corner2 = P( dx / 2,  dy / 2, nz))

class Gear_Box_Cover(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	self.dx_l = dx = L(mm = 2.50)

	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_
	gear_box_nz = gear_box.nz_l
	gear_box_sz = gear_box.sz_l

	self.block(comment = "Bevel Gear Box Cover Block",
	  material = Material("plastic", "ABS"),
	  color = Color("chartreuse"),
	  corner1 = P(-gear_box.dx / 2,
	    -gear_box.dy / 2, gear_box_sz),
	  corner2 = P(-gear_box.dx / 2 + dx,
	     gear_box.dy / 2, gear_box_nz))

class Gear_Box_Shelf(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	# Grab some values from *wheel_assembly:
	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_
	gear_box_shelf_sz = wheel_assembly.gear_box_shelf_sz_l

	# Grab some values form *bearing*:
	bearing = wheel_assembly.bearing_
	bearing_diameter = bearing.diameter_l
	bearing_width = bearing.width_l

	# Grab some values from *gear_box_cover*:
	gear_box_cover = wheel_assembly.gear_box_cover_
	gear_box_cover_dx = gear_box_cover.dx_l

	# Grab some values from *gear_box_side*:
	gear_box_side = wheel_assembly.wheel_side_gear_box_side_
	gear_box_side_thin_dy = gear_box_side.thin_dy_l

	# Grab some values from *twist_timing_pulley*:
	twist_timing_pulley = wheel_assembly.twist_timing_pulley_
	twist_timing_pulley_shaft_diameter = \
	  twist_timing_pulley.shaft_diameter_l

	# Remember some values into *self*:
	self.bearing_lip_dz_l = bearing_lip_dz = L(mm = 2.00)
	self.dz_l = dz = bearing_width + bearing_lip_dz
	
	dx2 = gear_box.dx / 2 - gear_box_cover_dx
	dy2 = gear_box.dy / 2 - gear_box_side_thin_dy
	z0 = gear_box_shelf_sz
	z1 = z0 + bearing_width
	z2 = z0 + dz

	# Generate the shelf block:
	self.block(comment = "Shelf Block",
	  material = Material("plastic", "ABS"),
	  color = Color("gold"),
	  corner1 = P(-dx2, -dy2, z0),
	  corner2 = P( dx2,  dy2, z2))

	# Cut a hole for the shaft:
	zero = L()
	self.hole(comment = "Shaft Hole",
	  diameter = twist_timing_pulley_shaft_diameter,
	  start = P(zero, zero, z0),
	  end = P(zero, zero, z2),
	  flags = "t")

	# Cut a bearing hole:
	self.hole(comment = "Bearing Hole",
	  diameter = bearing_diameter,
	  start = P(zero, zero, z0),
	  end = P(zero, zero, z1),
	  flags = "f")

class Gear_Box_Side(Part):
    def __init__(self, up, wheel_side):
	# Deal with *wheel_side* argument:
	assert isinstance(wheel_side, bool)
	if wheel_side:
	    name_prefix = "Wheel_Side"
	else:
	    name_prefix = "Non_Wheel_Side"
	self.wheel_side_b = wheel_side

	Part.__init__(self, up,
	  name = name_prefix + "_Gear_Box_Side")

    def construct(self):
	""" *Gear_Box_Side*: """

	wheel_side = self.wheel_side_b

	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	bearing = wheel_assembly.bearing_
	gear_box_top = wheel_assembly.gear_box_top_
	
	gear_box_top_base_dz = gear_box_top.base_dz_l

	# Grab the basic bevel gear box dimensions from *wheel_assembly*:
	gear_box = wheel_assembly.gear_box_
	gear_box_nz = gear_box.nz_l
	gear_box_sz = gear_box.sz_l
	gear_box_shelf_sz = wheel_assembly.gear_box_shelf_sz_l
	
	# Grab values from *bearing*:
	bearing_diameter = bearing.diameter_l
	bearing_width = bearing.width_l

	# Grab values from *gear_box_cover*:
	gear_box_cover = wheel_assembly.gear_box_cover_
	gear_box_cover_dx = gear_box_cover.dx_l

	# Grab values from *gear_box_shelf*:
	gear_box_shelf = wheel_assembly.gear_box_shelf_
	gear_box_shelf_dz = gear_box_shelf.dz_l

	horizontal_shaft = wheel_assembly.horizontal_shaft_
	horizontal_shaft_diameter = horizontal_shaft.diameter_l

	# Figure out the various dimensions of the shelf:
	self.shelf_nub_dy_l = shelf_nub_dy = L(mm = 5.00)
	self.shelf_nub_dz_l = shelf_nub_dz = L(mm = 5.00)
	self.shelf_nub_sz_l = shelf_nub_sz = \
	  gear_box_shelf_sz + gear_box_shelf_dz
	self.shelf_nub_nz_l = shelf_nub_nz = shelf_nub_sz + shelf_nub_dz

	self.bearing_lip_dy_l = bearing_lip_dy = L(mm = 2.00)
	self.thin_dy_l = thin_dy = bearing_width + bearing_lip_dy
	self.thick_dy_l = thick_dy = thin_dy + shelf_nub_dy
		
	# Identify various points on the Y axis:
	y1 = gear_box.dy / 2
	y2 = y1 - bearing_lip_dy
	y3 = y1 - thin_dy
	y4 = y1 - thick_dy
	if wheel_side:
	    y1 = -y1
	    y2 = -y2
	    y3 = -y3
	    y4 = -y4

	# Do the main block:
	corner1 = P(-gear_box.dx / 2 + gear_box_cover_dx,
	  y1, gear_box_sz)
	corner2 = P( gear_box.dx / 2 - gear_box_cover_dx,
	  y3, gear_box_nz - gear_box_top_base_dz)
	self.block(comment = self._name + " Main Block",
	  material = Material("plastic", "abs"),
	  color = Color("orange"),
	  corner1 = corner1,
	  corner2 = corner2)

	# Do the shelf block:
	corner1 = \
	  P(-gear_box.dx / 2 + gear_box_cover_dx, y3, shelf_nub_sz)
	corner2 = \
	  P( gear_box.dx / 2 - gear_box_cover_dx, y4, shelf_nub_nz)
	self.block(comment = self._name + " Shelf Block",
	  material = Material("plastic", "abs"),
	  color = Color("orange"),
	  corner1 = corner1,
	  corner2 = corner2)

	# Do the bearing hole:
	zero = L()
	self.hole(comment = "Bearing hole",
	  diameter = bearing_diameter,
	  start = P(zero, y3, zero),
	  end = P(zero, y2, zero),
	  flags = "f")

	# Do a full shaft hole for the *wheel_side*:
	if wheel_side:
	    diameter = (horizontal_shaft_diameter + bearing_diameter) / 2
	    self.hole(comment = "Shaft hole",
	      diameter = diameter,
	      start = P(zero, y3, zero),
	      end = P(zero, y1, zero),
	      flags = "t")

class Bevel_Gear_Table:
    def __init__(self):
	bevel_gears = [
	  Bevel_Gear_Data(
	    part_name = "A 1M 4-Y48015",
	    material = Material("plastic", "actel"),
	    pitch = 48,
	    teeth = 15,
	    pitch_diameter = L(inch = 0.313),
	    outside_diameter = L(inch = 0.350),
	    bore = L(inch = 0.125),
	    height = L(inch = "7/32"),
	    hub_diameter = L(inch = "1/4"),
	    hub_height = L(inch = "1/8"),
	    major_diameter = L(inch = "5/16"),
	    cost1_24 = 2.85,
	    cost25_99 = 1.73,
	    cost100_249 = 1.31,
	    cost250_999 = 1.11,
	    cost_1000_plus = 0.87),
	  Bevel_Gear_Data(
	    part_name = "A 1M 4-Y16016",
	    material = Material("plastic", "nylon"),
	    pitch = 16,
	    teeth = 16,
	    pitch_diameter = L(inch = 1.000),
	    outside_diameter = L(inch = 1.090),
	    bore = L(inch = 0.375),
	    height = L(inch = "3/4"),
	    hub_diameter = L(inch = "3/4"),
	    hub_height = L(inch = "7/16"),
	    major_diameter = L(inch = "1-1/16"),
	    cost1_24 = 3.41,
	    cost25_99 = 2.34,
	    cost100_249 = 1.89,
	    cost250_999 = 1.70,
	    cost_1000_plus = 1.47) ]

	self.table = table = {}
	for bevel_gear in bevel_gears:
	    part_name = bevel_gear.part_name
	    assert not part_name in table
	    table[part_name] = bevel_gear

class Gear_Box_Top(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	self.base_dz_l = base_dz = L(mm = 4.00)
	self.mount_dx_l = mount_dx = L(mm = 2.00)
	self.mount_dy_l = mount_dy = L(mm = 2.00)
	self.mount_dz_l = mount_dz = L(mm = 6.00)

	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_
	gear_box_nz = gear_box.nz_l

	# Grab some values from *gear_box_cover*:
	gear_box_cover = wheel_assembly.gear_box_cover_
	gear_box_cover_dx = gear_box_cover.dx_l

	# Grab some values from *gear_box_side*:
	gear_box_side = wheel_assembly.wheel_side_gear_box_side_
	gear_box_side_thin_dy = gear_box_side.thin_dy_l

	# Grab some values from *turn_table*:
	turn_table = wheel_assembly.turn_table_
	turn_table_size = turn_table.size_l
	turn_table_bottom_hole_pitch1 = turn_table.bottom_hole_pitch1_l
	turn_table_bottom_hole_diameter1 = turn_table.bottom_hole_diameter1_l

	# Grab some values from *twist_timing_pulley*:
	twist_timing_pulley = wheel_assembly.twist_timing_pulley_
	twist_timing_pulley_holes_count = twist_timing_pulley.holes_count_i
	twist_timing_pulley_hole_diameter = twist_timing_pulley.hole_diameter_l
	twist_timing_pulley_holes_radius = twist_timing_pulley.holes_radius_l
	twist_timing_pulley_shaft_diameter = \
	  twist_timing_pulley.shaft_diameter_l

	z3 = gear_box_nz
	z2 = z3 - base_dz / 2
	z1 = z3 - base_dz
	z0 = z1 - mount_dz

	# Get the basic base in place:
	self.block(comment = "Bevel Gear Box Top Block",
	  material = Material("plastic", "ABS"),
	  color = Color("azure"),
	  corner1 = P(-turn_table_size / 2, -turn_table_size / 2, z1),
	  corner2 = P( turn_table_size / 2,  turn_table_size / 2, z3))

	# Add some stuff to mount to:
	mount_dy2 = gear_box.dy / 2 - gear_box_side_thin_dy
	self.block(comment = "Bevel Gear Box Top Mount Block",
	  corner1 =
	  P(-gear_box.dx / 2 + gear_box_cover_dx, -mount_dy2, z0),
	  corner2 =
	  P( gear_box.dx / 2 - gear_box_cover_dx,  mount_dy2, z2))

	# Do the twist timing pulley holes:
	angle_delta = Angle(360) / twist_timing_pulley_holes_count
	for index in range(twist_timing_pulley_holes_count):
	    angle = index * angle_delta
	    x = twist_timing_pulley_holes_radius.cosine(angle)
	    y = twist_timing_pulley_holes_radius.sine(angle)
	    self.hole(comment = "Hole {0}".format(index),
	      diameter = twist_timing_pulley_hole_diameter,
	      start = P(x, y, z3),
	      end = P(x, y, z1),
	      flags = "t")

	# Do the shaft hole:
	zero = L()
	self.hole(comment = "Shaft Hole",
	  diameter = twist_timing_pulley_shaft_diameter,
	  start = P(zero, zero, z3),
	  end = P(zero, zero, z0),
	  flags = "t")

	# Do the turn table mounting holes:
	for x_sign in [-1, 1]:
	    x = x_sign * turn_table_bottom_hole_pitch1 / 2
	    for y_sign in [-1, 1]:
		y = y_sign * turn_table_bottom_hole_pitch1 / 2
		self.hole(comment =
		  "Turn Table Hole [{0}, {1}]".format(x_sign, y_sign),
		  diameter = turn_table_bottom_hole_diameter1,
		  start = P(x, y, z0),
		  end = P(x, y, z3),
		  flags = "t")

	# Pocket out some ..
	self.simple_pocket(comment = "Top Pocket",
	  corner1 =
	    P(-gear_box.dx / 2 + mount_dx + gear_box_cover_dx,
	    -mount_dy2 + mount_dy, z0),
	  corner2 =
	    P( gear_box.dx / 2 - mount_dx - gear_box_cover_dx,
	    mount_dy2 - mount_dy, z1))

class Shim(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	self.no_automatic_place()

	wheel_assembly = self.up
	bearing = wheel_assembly.bearing_
	bearing_diameter = bearing.diameter_l
	bearing_bore = bearing.bore_l

	self.bore_l = bore = bearing_bore
	self.diameter_l = diameter = (bearing_diameter + bore) / 2
	self.dz_l = dz = L(inch = .010)
	
	zero = L()
	self.cylinder(comment = "Shim Cylinder",
	  material = Material("steel", "stainless"),
	  color = Color("coral"),
	  diameter = diameter,
	  start = P(zero, zero, zero),
	  end = P(zero, zero, dz))

class Timing_Pulley(Part):
    def __init__(self, up, belt_class = "MXL", teeth_count = -1,
      belt_width = L(), lip_extra = L(), lip_width = L(),
      holes_count = 4, hole_diameter = L(inch = 0.0760),
      bearing_diameter = L(), bearing_width = L(), shaft_diameter = L()):
	Part.__init__(self, up)

	# Check argument types:
	zero = L()
      	assert isinstance(belt_class, str)
	assert isinstance(teeth_count, int)
	assert isinstance(belt_width, L) and belt_width > zero
	assert isinstance(lip_extra, L) and lip_extra > zero
	assert isinstance(lip_width, L) and lip_width > zero
	assert isinstance(holes_count, int) and holes_count >= 0
	assert isinstance(hole_diameter, L) and hole_diameter > zero
	assert isinstance(bearing_diameter, L)
	assert isinstance(bearing_width, L)
	assert isinstance(shaft_diameter, L)

	# Remember the arguments;
	self.belt_class_s = belt_class
	self.teeth_count_i = teeth_count
	self.belt_width_l = belt_width
	self.lip_extra_l = lip_extra
	self.lip_width_l = lip_width
	self.holes_count_i = holes_count
	self.hole_diameter_l = hole_diameter
	self.bearing_diameter_l = bearing_diameter
	self.bearing_width_l = bearing_width
	self.shaft_diameter_l = shaft_diameter

    def construct(self):
	""" *Timing_Pulley* construct method. """
	self.no_automatic_place()

	# Grab some value out of *self*:
	belt_class = self.belt_class_s
	teeth_count = self.teeth_count_i
	belt_width = self.belt_width_l
	lip_width = self.lip_width_l
	lip_extra = self.lip_extra_l
	holes_count = self.holes_count_i
	hole_diameter = self.hole_diameter_l
	bearing_diameter = self.bearing_diameter_l
	bearing_width = self.bearing_width_l
	shaft_diameter = self.shaft_diameter_l

	# The following URL has some great diagrams for timing belt
	# dimensions:
	#
	#    http://econobelt.com/timing_belts.htm
	#
	#         |<---------P---------->|
	#      -->|    |<--W             |                    |
	#         |                      |                    V
	#         +----+                 +----+             -----
	#        /      \               /      \
	#       /        \             /        \             H
	#      /          \           /          \              
	# ----+            +---------+            +-------  -----
	#    /              \                                 ^
	#   /<------A------->\                                |
	#
	#   P = Belt tooth pitch
	#   H = Belt tooth height
	#   A = Belt point angle
	#   W = Width width
	#   R = Transition radius
	#
	#    Belt Class       P       W       H       A       R
	#       MXL		.080"	.030"	.020"	40 deg  .005"
	#       XL		.200"	.051"	.050"	50 deg  .015"
	#       L		.375"	.128"	.075"	40 deg  .020"
	#
	# We need to figure out how much of the belt is on the
	# bottom and how much is transition:
	#
	#           A
	#          /|\
	#         / | \
	#        /  |  \ R
	#       /   |   \
	#      /    |    \
	#     /     |     \
	#    B------C------D
	#
	# Triagnle ABD is an isosceles triangle.
	# Angle <BAD is the point angle (i.e. 40 or 50 degrees).
	# Segment AC is the belt tooth height (i.e. .020", .050", or .075")
	#
	# From Mathworld:
	#
	#  http://mathworld.wolfram.com/IsoscelesTriangle.html
	#
	#     h = R * cos(1/2 * theta)			(1)
	#     x = R * sin(1/2 * theta)			(2)
	# 
	# where:
	#
	#     theta = <BAD
	#     h = segment AC
	#     x = CD
	#
	# Solving equations (1) and (2) for R:
	#
	#     R = h / cos(1/2 * theta)			(3)
	#     R = x / sin(1/2 * theta)			(4)
	#
	# Equating R from (3) and (4):
	#
	#     h / cos(1/2 * theta) = x / sin(1/2 * theta)	(5)
	#
	# Now solve for x:
	#
	#     x = h * sin(1/2 * theta) / cos(1/2 * theta)	(6)
	#     x = h * tan(1/2 * theta)			(7)
	#
	# The amount of the belt is under the slope is 2*x:
	#
	#    Belt Class      H       A	  x	   2x
	#       MXL		.020"	40 deg	.00727"  .01454"
	#       XL		.050"	50 deg  .02331"  .04662"
	#       L		.075"	40 deg  .02729"  .05458"
	#
	#    Belt Class       P       W       H      2x  	P-W-2x
	#       MXL		.080"	.030"	.020"  .014"	.036"
	#       XL		.200"	.051"	.050"  .047"	.102"
	#       L		.375"	.128"	.075"  .055"    .192"
	#
	# P-W-2x is the amount of belt on the bottom and 2x is the
	# the amount of transitionl belt:

	# Construct the *belt_info* table:
	belt_info = {}
	belt_info["MXL"] = (.080, .030, .020, .014, .036)
	belt_info["XL"]  = (.200, .051, .050, .047, .102)
	belt_info["L"] =   (.375, .128, .075, .055, .192)

	# Read out the *belt_record*:
	belt_record = belt_info[belt_class]
	p = belt_record[0]
	w = belt_record[1]
	h = belt_record[2]
	x2 = belt_record[3]
	pw2x = belt_record[4]

	pitch = L(inch = p)
	height = L(inch = h)

	# Compute the ratios of upper, transition and lower belt:
	w_ratio = w / p
	x2_ratio = x2 / p
	pw2x_ratio = pw2x / p

	# Compute the *tooth_angle* which is the repetition angle:
	tooth_angle = Angle(deg = 360) / teeth_count

	#   |<-- P/2 -->|
	#
	#   +-----------+
	#   |          /
	#   |         /
	#   |        /
	#   |       /
	#   |      /
	#   | A/2 /
	#   |    /
	#   |   / R
	#   |  /
	#   | /
	#   |/
	#   +
	#
	#   P/2 = R * sin()		(1)
	#   R = P/(2*sin(a))		(2)

	# Compute the tooth radis (R in equation (2) .):
	tooth_radius = pitch / (2 * (tooth_angle / 2).sine())
	#print("tooth_angle={0:d}, tooth_angle.sine()={1}".
	#  format(tooth_angle, (tooth_angle / 2).sine()))
	#print("pitch={0}, tooth_radius={1} circum={2} p*tr={3}".
	#  format(pitch, tooth_radius,
	#  tooth_radius * 2 * 3.1415629, pitch * teeth_count))

	# Compute the various Z heights:
	zero = L()
	z0 = zero			# Bottom
	z1 = z0 + lip_width / 2		# Extend into lip to force a weld
	z2 = z0 + lip_width		# Lip top and gear bottom
	z3 = z2 + belt_width		# Gear top

	r1 = tooth_radius
	r2 = r1 - height
	self.holes_radius_l = holes_radius = (r2 + bearing_diameter / 2) / 2

	self.cylinder(comment = "Pulley Lip",
	  material = Material("plastic", "ABS"),
	  color = Color("crimson"),
	  diameter = tooth_radius * 2 + lip_extra,
	  start = P(zero, zero, z0),
	  end = P(zero, zero, z2),
	  sides = teeth_count)
	
	quick = False
	quick = True
	if quick:
	    self.cylinder(comment = "Pulley Gear",
	    diameter = tooth_radius * 2,
	    start = P(zero, zero, z1),
	    end = P(zero, zero, z3),
	    sides = teeth_count)
	else:
	    # Compute *outer_path*:
	    outer_path = []
	    for index in range(teeth_count):
		# Compute the 4 angles:
		angle0 = tooth_angle * (index + 0.00)
		angle1 = tooth_angle * (index + w_ratio)
		angle2 = tooth_angle * (index + w_ratio + x2_ratio / 2)
		angle3 = \
		  tooth_angle * (index + w_ratio + x2_ratio / 2 + pw2x_ratio)

		# Add everything to the *outer_path*:
		outer_path.append(
		  Bend(P(r2.cosine(angle0), r2.sine(angle0), zero), zero))
		outer_path.append(
		  Bend(P(r2.cosine(angle1), r2.sine(angle1), zero), zero))
		outer_path.append(
		  Bend(P(r1.cosine(angle2), r1.sine(angle2), zero), zero))
		outer_path.append(
		  Bend(P(r1.cosine(angle3), r1.sine(angle3), zero), zero))

	    # Now extrude the shape:
	    self.extrude(comment = "Gear",
	      outer_path = outer_path,
	      start = P(zero, zero, z1),
	      end = P(zero, zero, z3))

	# Put in the 
	if holes_count > 0 and hole_diameter > zero:
	    angle_delta = Angle(deg = 360) / holes_count
	    for index in range(holes_count):
		angle = angle_delta * index
		x = holes_radius.cosine(angle)
		y = holes_radius.sine(angle)
		self.hole(comment = "Screw hole {0}".format(index),
		  diameter = hole_diameter,
		  start = P(x, y, z3),
		  end = P(x, y, z0),
		  flags = "t")

	if shaft_diameter > zero:
	    self.hole(comment = "Shaft Hole",
	      diameter = shaft_diameter,
	      start = P(zero, zero, z3),
	      end = P(zero, zero, zero),
	      flags = "t")

	if bearing_diameter > zero and bearing_width > zero:
	    self.hole(comment = "Bearing Hole",
	      diameter = bearing_diameter,
	      start = P(zero, zero, z3),
	      end = P(zero, zero, z3 - bearing_width),
	      flags = "f")

    def update(self, belt_class = None, teeth_count = None,
      belt_width = None, lip_extra = None, lip_width = None,
      holes_count = None, hole_diameter = None,
      bearing_diameter = None, bearing_width = None, shaft_diameter = None):

	# Check argument types:
	none_type = type(None)
	zero = L()
      	assert type(belt_class) == none_type or \
	  isinstance(belt_class, str)
	assert type(teeth_count) == none_type or \
	  isinstance(teeth_count, int)
	assert type(belt_width) == none_type or \
	  isinstance(belt_width, L) and belt_width > zero
	assert type(lip_extra) == none_type or \
	  isinstance(lip_extra, L) and lip_extra > zero
	assert type(lip_width) == none_type or \
	  isinstance(lip_width, L) and lip_width > zero
	assert type(holes_count == none_type) or \
	  isinstance(holes_count, int) and holes_count >= 0
	assert type(hole_diameter) == none_type or \
	  isinstance(hole_diameter, L) and hole_diameter > zero
	assert type(bearing_diameter) == none_type or \
	  isinstance(bearing_diameter, L)
	assert type(bearing_width) == none_type or \
	  isinstance(bearing_width, L)
	assert type(shaft_diameter) == none_type or \
	  isinstance(shaft_diameter, L)

	# Remember the arguments;
	if type(belt_class) != none_type:
	    self.belt_class_s = belt_class
	if type(teeth_count) != none_type:
	    self.teeth_count_i = teeth_count
	if type(belt_width) != none_type:
	    self.belt_width_l = belt_width
	if type(lip_extra) != none_type:
	    self.lip_extra_l = lip_extra
	if type(lip_width) != none_type:
	    self.lip_width_l = lip_width
	if type(holes_count) != none_type:
	    self.holes_count_i = holes_count
	if type(hole_diameter) != none_type:
	    self.hole_diameter_l = hole_diameter
	if type(bearing_diameter) != none_type:
	    self.bearing_diameter_l = bearing_diameter
	if type(bearing_width) != none_type:
	    self.bearing_width_l = bearing_width
	if type(shaft_diameter) != none_type:
	    self.shaft_diameter_l = shaft_diameter

class Motor_Base(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	motor_assembly = self.up
	synchro_drive = motor_assembly.up
	wheel_assembly = synchro_drive.wheel_assembly_

	turn_table = wheel_assembly.turn_table_
	turn_table_size = turn_table.size_l
	turn_table_bore = turn_table.bore_l

	self.dx_l = dx = turn_table_size
	self.dy_l = dy = turn_table_size
	self.dz_l = dz = L(mm = 3.0)

	z0 = motor_assembly.base_sz_l
	z1 = z0 + dz

	self.block(comment = "Motor Base",
	  material = Material("plastic", "ABS"),
	  color = Color("lavender"),
	  corner1 = P(-dx/2, -dy/2, z0),
	  corner2 = P( dx/2,  dy/2, z1))

	zero = L()
	self.hole(comment = "Shaft Hole",
	  diameter = turn_table_bore,
	  start = P(zero, zero, z1),
	  end = P(zero, zero, z0),
	  flags = "t")

class Turn_Table(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	zero = L()
	self.turn_table_upper_sheet_ = Turn_Table_Sheet(self, is_upper = True)
	self.turn_table_bottom_sheet_ = Turn_Table_Sheet(self, is_upper = False)

    def construct(self):
	# Provide some 
	self.size_l = size = L(inch = 3.00)
	self.bore_l = bore = L(inch = 1.31)
	self.top_hole_pitch1_l = top_hole_pitch1 = L(inch = 2.56)
	self.top_hole_pitch2_l = top_hole_pitch2 = L(inch = 2.13)
	self.bottom_hole_pitch1_l = bottom_hole_pitch1 = L(inch = 2.56)
	self.top_hole_diameter1_l = top_hole_diameter1 = L(inch = 0.16)
	self.top_hole_diameter2_l = top_hole_diameter2 = L(inch = 0.094)
	self.bottom_hole_diameter1_l = bottom_hole_diameter1 = L(inch = 0.16)
	self.sheet_thickness_l = sheet_thickness = L(inch = 0.03)
	self.height_l = height = L(inch = 0.26)

class Turn_Table_Sheet(Part):
    def __init__(self, up, is_upper = False):
	name = "Bottom_Turn_Table_Sheet"
	if is_upper:
	    name = "Top_Turn_Table_Sheet"
	Part.__init__(self, up, name = name)
	self.is_upper = is_upper
 
    def construct(self):
	
	#print("Turn_Table_Sheet:size={0} sheet_thickness={1}".
	#  format(size, sheet_thickness))

	# Grab some values from *turn_table*:
	turn_table = self.up
	turn_table_size = turn_table.size_l
	turn_table_bore = turn_table.bore_l
	turn_table_top_hole_pitch1 = turn_table.top_hole_pitch1_l
	turn_table_top_hole_pitch2 = turn_table.top_hole_pitch2_l
	turn_table_bottom_hole_pitch1 = turn_table.bottom_hole_pitch1_l
	turn_table_top_hole_diameter1 = turn_table.top_hole_diameter1_l
	turn_table_top_hole_diameter2 = turn_table.top_hole_diameter2_l
	turn_table_bottom_hole_diameter1 = turn_table.bottom_hole_diameter1_l
	turn_table_sheet_thickness = turn_table.sheet_thickness_l
	turn_table_height = turn_table.height_l

	# Make the original sheet:
	zero = L()
	is_upper = self.is_upper
	comment = "Bottom Turn Table Sheet"
	z0 = zero
	if is_upper:
	    comment = "Top Turn Table Sheet"
	    z0 = turn_table_height - turn_table_sheet_thickness
	z1 = z0 + turn_table_sheet_thickness
	#print("is_upper={0} comment={1} z0={2} z1={3}".
	#  format(is_upper, comment, z0, z1))

	self.block(comment = comment,
	  material = Material("Steel", "Galvanized"),
	  color = Color("purple"),
	  corner1 = P(-turn_table_size / 2, -turn_table_size / 2, z0),
	  corner2 = P( turn_table_size / 2,  turn_table_size / 2, z1))

	# Bore the hole through the center:
	self.hole(comment = "Bore Hole",
	  diameter = turn_table_bore,
	  start = P(zero, zero, zero),
	  end =   P(zero, zero, turn_table_sheet_thickness),
	  flags = "t")

	# Do the the holes:
	for x_sign in [-1, 1]:
	    x1 = x_sign * turn_table_top_hole_pitch1 / 2
	    x2 = x_sign * turn_table_top_hole_pitch2 / 2
	    x3 = x_sign * turn_table_bottom_hole_pitch1 / 2
	    for y_sign in [-1, 1]:
		y1 = y_sign * turn_table_top_hole_pitch1 / 2
		y2 = y_sign * turn_table_top_hole_pitch2 / 2
		y3 = y_sign * turn_table_bottom_hole_pitch1 / 2
		if is_upper:
		    self.hole(comment = "Top Hole1 [{0}, {1}]".
		      format(x_sign, y_sign),
		      diameter = turn_table_top_hole_diameter1,
		      start = P(x1, y1, z1),
		      end = P(x1, y1, z0),
		      flags = "t")
		    self.hole(comment = "Top Hole2 [{0}, {1}]".
		      format(x_sign, y_sign),
		      diameter = turn_table_top_hole_diameter2,
		      start = P(x2, y2, z1),
		      end = P(x2, y2, z0),
		      flags = "t")
		else:
		    self.hole(comment = "Bottom Hole1 [{0}, {1}]".
		      format(x_sign, y_sign),
		      diameter = turn_table_bottom_hole_diameter1,
		      start = P(x3, y3, z1),
		      end = P(x3, y3, z0),
		      flags = "t")

class Wheel(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	# Define some constants:
	self.diameter_l = diameter = L(inch = "4-7/8")
	self.width_l = width = L(inch = 0.8)
	self.radius_l = radius = diameter / 2

	zero = L()
	self.cylinder(comment = "Basic Wheel",
	  material = Material(),
	  color = Color("lime"),
	  diameter = diameter, 
	  start = P(zero, -radius - width / 2, zero),
	  end = P(zero, -radius + width / 2, zero),
	  sides = 30)

class Horizontal_Shaft(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_

	# Grab *wheel*:
	wheel = wheel_assembly.wheel_

	# Some values from *bevel_gear*:
	bevel_gear = wheel_assembly.bevel_gear_
	bevel_gear_data = bevel_gear.bevel_gear_data_o
	bevel_gear_bore = bevel_gear_data.bore

	# Grab some values from *gear_box_side*:
	gear_box_side = wheel_assembly.wheel_side_gear_box_side_
	gear_box_side_bearing_lip_dy = \
	  gear_box_side.bearing_lip_dy_l

	# Load up *self*:
	self.diameter_l = diameter = bevel_gear_bore

	# Build the shaft:
	zero = L()
	self.cylinder(comment = "Horizontal Shaft",
	  material = Material("Steel", ""),
	  color = Color("cyan"),
	  diameter = diameter,
	  start = P(zero, wheel.w.x, zero),
          end = P(zero,
	    gear_box.dy / 2 - gear_box_side_bearing_lip_dy, zero))

class Vertical_Shaft(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	bevel_gear = self.up.bevel_gear_
	bevel_gear_data = bevel_gear.bevel_gear_data_o
	bevel_gear_bore = bevel_gear_data.bore
	bevel_gear_major_diameter = bevel_gear_data.major_diameter
	bevel_gear_height = bevel_gear_data.height

	self.diameter_l = diameter = bevel_gear_bore

	zero = L()
	self.cylinder(comment = "Horizontal Shaft",
	  material = Material("Steel", ""),
	  color = Color("cyan"),
	  diameter = diameter,
	  start = P(zero, zero, bevel_gear_major_diameter - bevel_gear_height),
          end = P(zero, zero, L(mm = 100.0)))

if __name__ == "__main__":
    ezcad = EZCAD3(0)
    test = Synchro_Drive(None)
    test.process(ezcad)

