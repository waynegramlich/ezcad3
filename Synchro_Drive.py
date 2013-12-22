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
	synchro_drive = self.up
	wheel_assembly = synchro_drive.wheel_assembly_
	turn_table = wheel_assembly.turn_table_
	gear_box = wheel_assembly.gear_box_

	self.base_bz_l = base_bz = gear_box.t.z + turn_table.height_l

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

	self.bottom_bearing_ = Bearing(self)
	self.gear_box_ = Gear_Box(self)
	self.gear_box_bottom_ = Gear_Box_Bottom(self)
	self.gear_box_bottom_screws_ = Gear_Box_Screws(self, label="Bottom")
	self.gear_box_shelf_ = Gear_Box_Shelf(self)
	self.gear_box_shelf_screws_ = Gear_Box_Screws(self, label="Shelf")
	self.gear_box_top_ = Gear_Box_Top(self)
	self.gear_box_top_screws_ = Gear_Box_Screws(self, label="Top")
	self.horizontal_shaft_ = Horizontal_Shaft(self)
	self.north_bearing_ = Bearing(self)
	self.north_gear_box_side_ = Gear_Box_Side(self, False)
	self.south_bearing_ = Bearing(self)
	self.south_bevel_gear_ = \
	  Bevel_Gear(self, part_name = "A 1M 4-Y16016", color = Color("red"))
	self.south_gear_box_side_ = Gear_Box_Side(self, True)
	self.south_shim_ = Shim(self)
	self.top_bearing_ = Bearing(self)
	self.top_bevel_gear_ = \
	  Bevel_Gear(self, part_name = "A 1M 4-Y16016", color = Color("red"))
	self.top_shim_ = Shim(self)
	self.turn_table_ = Turn_Table(self)
	self.turn_table_bottom_screws_ = Turn_Table_Bottom_Screws(self)
	self.vertical_shaft_ = Vertical_Shaft(self)
	self.west_gear_box_cover_ = Gear_Box_Cover(self)
	self.wheel_ = Wheel(self)
 
    def construct(self):
	""" *Wheel_Assembly*: Construct method. """

	# The origin of this assembly is at the logical intersection
	# of the vertical shaft and the horizontal shaft.  The horizontal
	# axis is aligned with the Y axis; this causes the wheel to be
	# oriented to roll in the direction of the X axis.

	# Grab some values from *wheel_assembly* (i.e. *self*):
	wheel_assembly = self
	gear_box = wheel_assembly.gear_box_
	gear_box_cover = wheel_assembly.west_gear_box_cover_
	gear_box_bottom = self.gear_box_bottom_
	gear_box_bottom_screws = self.gear_box_bottom_screws_
	gear_box_side = self.south_gear_box_side_
	gear_box_shelf = self.gear_box_shelf_
	gear_box_shelf_screws = self.gear_box_shelf_screws_
	gear_box_top = self.gear_box_top_
	gear_box_top_screws = self.gear_box_top_screws_
	south_bevel_gear = wheel_assembly.south_bevel_gear_
	south_shim = wheel_assembly.south_shim_
	north_bearing = wheel_assembly.north_bearing_
	top_bearing = wheel_assembly.top_bearing_
	synchro_drive = wheel_assembly.up
	turn_table = self.turn_table_
	twist_timing_pulley = wheel_assembly.twist_timing_pulley_
	bottom_bearing = wheel_assembly.bottom_bearing_
	top_bevel_gear = wheel_assembly.top_bevel_gear_
	vertical_shaft = wheel_assembly.vertical_shaft_
	top_shim = wheel_assembly.top_shim_
	south_bearing = wheel_assembly.south_bearing_

	# Grab some values from *synchro_drive*:
	motor_assembly = synchro_drive.motor_assembly_
	motor_base = motor_assembly.motor_base_

	# Compute various locations along the Y axis:
	zero = L()
	y0 = zero
	y1 = y0 + south_bevel_gear.major_diameter_l
	y2 = y1 + south_shim.width_l
	y3 = y2 + south_bearing.width_l / 2
	y4 = y2 + south_bearing.width_l
	y5 = y4 + gear_box_side.bearing_lip_dy_l

	# The pockets that screws go to are all the same for the top, shelf,
	# and bottom:
	self.pocket_wall_dx_l = pocket_wall_dx = L(mm = 3.00)
	self.pocket_wall_dy_l = pocket_wall_dy = L(mm = 3.00)
	self.pocket_dx_l = \
	  gear_box.dx - 2 * gear_box_cover.dx -  2 * pocket_wall_dx
	self.pocket_dy_l = \
	  gear_box.dy - 2 * gear_box_side.dy - 2 * pocket_wall_dy
	self.pocket_dz_l = L(mm = 6.00)

	# Update the timing pulley:
	self.twist_timing_pulley_.update(
	  lip_width = turn_table.height_l + L(mm = 1.000) + motor_base.dz,
	  shaft_diameter = vertical_shaft.diameter_l + L(inch = "1/4"))

	# Place various items into the wheel assembly:
	one = L(mm = 1.0)
	x_axis = P(one, zero, zero)
	degrees_90 = Angle(deg = 90)

	# Place the 4 bearings:
	south_bearing.place(axis = x_axis,
	  rotate = degrees_90, translate = P(zero, -y3, zero))
	north_bearing.place(axis = x_axis,
	  rotate = degrees_90, translate = P(zero,  y3, zero))
	bottom_bearing.place(translate = P(zero, zero, gear_box_shelf.t.z -
	  gear_box_shelf.bearing_lip_dz_l - bottom_bearing.width_l / 2))
	top_bearing.place(translate = P(zero, zero,
	  gear_box.t.z + twist_timing_pulley.dz - bottom_bearing.width_l / 2))

	# Place the horizontal bevel gear.  (The vertical bevel gear is
	# is already placed relative to the origin.):
	south_bevel_gear.place(axis = x_axis, rotate = degrees_90)

	# Place the shims on next to the bevel gears:
	top_shim.place(translate = P(zero, zero, y1))
	south_shim.place(axis = x_axis,
	  rotate = degrees_90, translate = P(zero, -y1, zero))

	# Drop in the turn table:
	turn_table.place(translate = P(zero, zero, gear_box.t.z))

	# Drop in the two pullies:
	twist_timing_pulley.place(translate = P(zero, zero, gear_box.t.z))

	# Configure the screws:
	gear_box_bottom_screws.configure(
	  z = gear_box_bottom.screw_z_l,
	  gear_box_part = gear_box_bottom)
	gear_box_shelf_screws.configure(
	  z = gear_box_shelf.screw_z_l,
	  gear_box_part = gear_box_shelf)
	gear_box_top_screws.configure(
	  z = gear_box_top.screw_z_l,
	  gear_box_part = gear_box_top)

# Various parts:

class Bearing(Part):
    def construct(self):
	# R6ZZ Shield Ball Bearing from VBX.Com:
	self.bore_l = bore = L(inch = "3/8")
	self.diameter_l = diameter = L(inch = "7/8")
	self.width_l = width = L(inch = "9/32")

	zero = L()
	self.cylinder(comment = "Bearing Cylinder",
	  material = Material("steel", ""),
	  color = Color("blue"),
	  diameter = diameter, 
	  start = P(zero, zero, -width / 2),
	  end =   P(zero, zero,  width / 2),
	  sides = 30)
	self.hole(comment = "Bore Hole",
	  diameter = bore,
	  start = P(zero, zero, -width / 2),
	  end =   P(zero, zero,  width / 2),
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

class Gear_Box(Part):
    def construct(self):
	# Grab some values *wheel_assembly*:
	wheel_assembly = self.up
	bearing = wheel_assembly.south_bearing_
	bevel_gear = wheel_assembly.top_bevel_gear_
	south_shim = wheel_assembly.south_shim_
	top_shim = wheel_assembly.top_shim_
	gear_box_side = wheel_assembly.south_gear_box_side_

	# Figure out *y5* the distance from the center to the outside Y wall:
	zero = L()
	y0 = zero
	y1 = y0 + bevel_gear.major_diameter_l
	y2 = y1 + south_shim.width_l
	y3 = y2 + bearing.width_l / 2
	y4 = y2 + bearing.width_l
	y5 = y4 + gear_box_side.bearing_lip_dy_l

	# Specify *dx*, *dy*, *tz*, and *bz*:
	dx = L(mm = 40.00)
	dy  = y5 * 2
	tz = L(mm = 75.00)
	bz = L(mm = -25.00)

	# Now specify the virtual box for the gear box:
	self.virtual_box(comment = "Virtual box for gear box",
	  corner1 = P(-dx / 2, -dy / 2, bz),
	  corner2 = P( dx / 2,  dy / 2, tz))

class Gear_Box_Bottom(Part):
    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_
	gear_box_cover = wheel_assembly.west_gear_box_cover_
	gear_box_side = wheel_assembly.south_gear_box_side_

	# Specify the *base_dz* height and pocket dx/dy/dz:
	base_dz = L(mm = 5.00)
	pocket_dx = wheel_assembly.pocket_dx_l
	pocket_dy = wheel_assembly.pocket_dy_l
	pocket_dz = wheel_assembly.pocket_dz_l

	z0 = gear_box.b.z
	z1 = z0 + base_dz
	z2 = z1 + pocket_dz / 2
	z3 = z2 + pocket_dz / 2
	self.screw_z_l = z2

	# Get the basic base in place:
	self.block(comment = "Gear Box Bottom Block",
	  material = Material("plastic", "ABS"),
	  color = Color("grey"),
	  corner1 = P(gear_box.w.x + gear_box_cover.dx,
	    gear_box.s.y + gear_box_side.dy, z0),
	  corner2 = P(gear_box.e.x - gear_box_cover.dx,
	    gear_box.n.y - gear_box_side.dy, z3))

	# Pocket out some the mounting area:
	self.simple_pocket(comment = "Gear Box Bottom Pocket",
	  bottom_corner = P(-pocket_dx / 2, -pocket_dy / 2, z1),
	  top_corner =    P( pocket_dx / 2,  pocket_dy / 2, z3))

class Gear_Box_Cover(Part):
    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_

	# Construct the cover out of a block of material:
	dx = L(mm = 2.50)
	self.block(comment = "Bevel Gear Box Cover Block",
	  material = Material("plastic", "ABS"),
	  color = Color("chartreuse"),
	  corner1 = P(gear_box.w.x,      gear_box.s.y, gear_box.b.z),
	  corner2 = P(gear_box.w.x + dx, gear_box.n.y, gear_box.t.z))

	# Put some cover mount holes here:

class Gear_Box_Screws(Part):
    def __init__(self, up, label):
	""" *Gear_Box_Screws*: """

	# Initialize the *Part* super-class:
	Part.__init__(self, up)
	assert isinstance(label, str)

	# There are two screws on each of the for sides {East,West,North,South}.
	# For each side there is a screw one side and one on the other.
	# Thus, *en_screw* is on the East side with towards the North:
	self.en_screw_ = Fastener(self)
	self.es_screw_ = Fastener(self)
	self.wn_screw_ = Fastener(self)
	self.ws_screw_ = Fastener(self)
	self.ne_screw_ = Fastener(self)
	self.nw_screw_ = Fastener(self)
	self.se_screw_ = Fastener(self)
	self.sw_screw_ = Fastener(self)
	self.z_l = L()
	self._label = label
	self._gear_box_part = None

    def configure(self, z = L(), gear_box_part = None):
	""" *Gear_Box_Screws*: """

	# Check argument types:
	assert isinstance(z, L)
	assert isinstance(gear_box_part, Part)
	
	# Load up *self*:
	self.z_l = z
	self._gear_box_part = gear_box_part

    def construct(self):
	""" *Gear_Box_Screws*: """

	#print("=>Gear_Box_Screws.construct():label='{0}'".format(self._label))

	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_
	west_gear_box_cover = wheel_assembly.west_gear_box_cover_
	#east_gear_box_cover = wheel_assembly.east_gear_box_cover_
	east_gear_box_cover = wheel_assembly.west_gear_box_cover_
	south_gear_box_side = wheel_assembly.south_gear_box_side_
	north_gear_box_side = wheel_assembly.north_gear_box_side_
	gear_box_top = wheel_assembly.gear_box_top_

	# Grab the pocket size from *wheel_assembly*:
	pocket_dx = wheel_assembly.pocket_dx_l
	pocket_dy = wheel_assembly.pocket_dy_l

	# Compute various locations along the X axis:
	x0 = gear_box.w.x
	x1 = x0 + west_gear_box_cover.dx
	x2 = x1 + wheel_assembly.pocket_wall_dy_l
	x3 = x2 + pocket_dx * 0.20
	x4 = x3 + pocket_dx * 0.60
	x5 = x4 + pocket_dx * 0.20
	x6 = x5 + wheel_assembly.pocket_wall_dy_l
	x7 = x6 + east_gear_box_cover.dx

	# Compute various locations along the Y Axis:
	y0 = gear_box.s.y
	y1 = y0 + south_gear_box_side.dy
	y2 = y1 + wheel_assembly.pocket_wall_dy_l
	y3 = y2 + pocket_dy * 0.10
	y4 = y3 + pocket_dy * 0.80
	y5 = y4 + pocket_dy * 0.10
	y6 = y5 + wheel_assembly.pocket_wall_dy_l
	y7 = y6 + north_gear_box_side.dy

	# *z* is the height for the screws:
	z = self.z_l

	# Grab the 8 screws from *self*:
	en_screw = self.en_screw_
	es_screw = self.es_screw_
	wn_screw = self.wn_screw_
	ws_screw = self.ws_screw_
	ne_screw = self.ne_screw_
	nw_screw = self.nw_screw_
	se_screw = self.se_screw_
	sw_screw = self.sw_screw_

	# Configure all 8 screws:
	en_screw.configure(
	  comment = "East Side North Screw",
	  diameter_pitch = "#4-40",
	  start = P(x7, y4, z),
	  end = P(x5, y4, z))
	es_screw.configure(
	  comment = "East Side South Screw",
	  diameter_pitch = "#4-40",
	  start = P(x7, y3, z),
	  end = P(x5, y3, z))
	wn_screw.configure(
	  comment = "West Side North Screw",
	  diameter_pitch = "#4-40",
	  start = P(x0, y4, z),
	  end = P(x2, y4, z))
	ws_screw.configure(
	  comment = "West Side South Screw",
	  diameter_pitch = "#4-40",
	  start = P(x0, y3, z),
	  end = P(x2, y3, z))
	ne_screw.configure(
	  comment = "North Side East Screw",
	  diameter_pitch = "#4-40",
	  start = P(x4, y7, z),
	  end = P(x4, y5, z))
	nw_screw.configure(
	  comment = "North Side West Screw",
	  diameter_pitch = "#4-40",
	  start = P(x3, y7, z),
	  end = P(x3, y5, z))
	se_screw.configure(
	  comment = "South Side East Screw",
	  diameter_pitch = "#4-40",
	  start = P(x4, y0, z),
	  end = P(x4, y2, z))
	sw_screw.configure(
	  comment = "South Side West Screw",
	  diameter_pitch = "#4-40",
	  start = P(x3, y0, z),
	  end = P(x3, y2, z))

	# Drill 8 holes for in *gear_box_part*:
	gear_box_part = self._gear_box_part
	en_screw.drill(gear_box_part)
	es_screw.drill(gear_box_part)
	wn_screw.drill(gear_box_part)
	ws_screw.drill(gear_box_part)
	ne_screw.drill(gear_box_part)
	nw_screw.drill(gear_box_part)
	se_screw.drill(gear_box_part)
	sw_screw.drill(gear_box_part)

	# Drill the sides:
	ne_screw.drill(north_gear_box_side)
	nw_screw.drill(north_gear_box_side)
	se_screw.drill(south_gear_box_side)
	sw_screw.drill(south_gear_box_side)

	# Drill the covers:
	#en_screw.drill(east_gear_box_cover)
	#es_screw.drill(east_gear_box_cover)
	wn_screw.drill(west_gear_box_cover)
	ws_screw.drill(west_gear_box_cover)

	#print("<=Gear_Box_Screws.construct():label='{0}'".format(self._label))

class Gear_Box_Shelf(Part):
    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	bearing = wheel_assembly.south_bearing_
	gear_box = wheel_assembly.gear_box_
	gear_box_cover = wheel_assembly.west_gear_box_cover_
	gear_box_side = wheel_assembly.south_gear_box_side_
	twist_timing_pulley = wheel_assembly.twist_timing_pulley_
	top_bevel_gear = wheel_assembly.top_bevel_gear_
	top_shim = wheel_assembly.top_shim_

	# Remember some values into *self*:
	self.bearing_lip_dz_l = bearing_lip_dz = L(mm = 2.00)
	pocket_dx = wheel_assembly.pocket_dx_l
	pocket_dy = wheel_assembly.pocket_dy_l
	pocket_dz = wheel_assembly.pocket_dz_l

	# Compute various Z coordinates:
	z0 = top_bevel_gear.major_diameter_l + top_shim.dz - pocket_dz
	z1 = z0 + pocket_dz / 2
	z2 = z1 + pocket_dz / 2
	z3 = z2 + bearing.width_l
	z4 = z3 + bearing_lip_dz
	self.screw_z_l = z1

	dx = gear_box.dx - 2 * gear_box_cover.dx
	dy = gear_box.dy - 2 * gear_box_side.dy

	# Generate the shelf block:
	self.block(comment = "Shelf Block",
	  material = Material("plastic", "ABS"),
	  color = Color("gold"),
	  corner1 = P(-dx / 2, -dy / 2, z0),	# should z0
	  corner2 = P( dx / 2,  dy / 2, z4))

	# Cut a hole for the shaft:
	zero = L()
	self.hole(comment = "Shaft Hole",
	  diameter = twist_timing_pulley.shaft_diameter_l,
	  start = P(zero, zero, z0),
	  end = P(zero, zero, z4),
	  flags = "t")

	# Cut a bearing hole:
	self.hole(comment = "Bearing Hole",
	  diameter = bearing.diameter_l,
	  start = P(zero, zero, z0),
	  end = P(zero, zero, z3),
	  flags = "f")

	# Cut out the pocket:
	self.simple_pocket(comment = "Mount Pocket",
	  bottom_corner = P(-pocket_dx / 2, -pocket_dy / 2, z0),
	  top_corner = P( pocket_dx / 2,  pocket_dy / 2, z1))

class Gear_Box_Side(Part):
    def __init__(self, up, south_side):
	# Deal with *south_side* argument:
	assert isinstance(south_side, bool)
	if south_side:
	    name_prefix = "South_Side"
	else:
	    name_prefix = "North_Side"
	self.south_side_b = south_side

	Part.__init__(self, up,
	  name = name_prefix + "_Gear_Box_Side")

    def construct(self):
	""" *Gear_Box_Side*: """

	south_side = self.south_side_b

	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	bearing = wheel_assembly.south_bearing_
	gear_box = wheel_assembly.gear_box_
	gear_box_bottom = wheel_assembly.gear_box_bottom_
	gear_box_cover = wheel_assembly.west_gear_box_cover_
	gear_box_shelf = wheel_assembly.gear_box_shelf_
	gear_box_top = wheel_assembly.gear_box_top_
	horizontal_shaft = wheel_assembly.horizontal_shaft_
	
	self.bearing_lip_dy_l = bearing_lip_dy = L(mm = 2.00)
		
	# Identify various points on the Y axis:
	y1 = gear_box.dy / 2
	y2 = y1 - bearing_lip_dy
	y3 = y2 - bearing.width_l
	if south_side:
	    y1 = -y1
	    y2 = -y2
	    y3 = -y3

	# Do the main block:
	self.block(comment = self._name + " Main Block",
	  material = Material("plastic", "abs"),
	  color = Color("orange"),
	  corner1 =
	    P(gear_box.w.x + gear_box_cover.dx, y1, gear_box.b.z),
	  corner2 = P(gear_box.e.x - gear_box_cover.dx, y3,
	    gear_box.t.z - gear_box_top.base_dz_l))

	# Do the bearing hole:
	zero = L()
	self.hole(comment = "Bearing hole",
	  diameter = bearing.diameter_l,
	  start = P(zero, y3, zero),
	  end = P(zero, y2, zero),
	  flags = "f")

	# Do a full shaft hole for the *south_side*:
	if south_side:
	    diameter = (horizontal_shaft.diameter_l + bearing.diameter_l) / 2
	    self.hole(comment = "Shaft hole",
	      diameter = diameter,
	      start = P(zero, y3, zero),
	      end = P(zero, y1, zero),
	      flags = "t")

class Gear_Box_Top(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_
	gear_box_cover = wheel_assembly.west_gear_box_cover_
	gear_box_side = wheel_assembly.south_gear_box_side_
	turn_table = wheel_assembly.turn_table_
	twist_timing_pulley = wheel_assembly.twist_timing_pulley_

	# Grab the pocket size from *wheel_assembly*:
	pocket_dx = wheel_assembly.pocket_dx_l
	pocket_dy = wheel_assembly.pocket_dy_l
	pocket_dz = wheel_assembly.pocket_dz_l
	self.base_dz_l = base_dz = L(mm = 5.00)

	# Compute various locations along the Z Axis:
	z0 = gear_box.t.z - self.dz
	z1 = z0 + pocket_dz / 2
	z2 = z1 + pocket_dz / 2
	z3 = z2 + base_dz / 2
	z4 = z3 + base_dz / 2
	self.screw_z_l = z1

	# Get the basic base in place:
	self.block(comment = "Bevel Gear Box Top Block",
	  material = Material("plastic", "ABS"),
	  color = Color("azure"),
	  corner1 = P(-turn_table.size_l / 2, -turn_table.size_l / 2, z2),
	  corner2 = P( turn_table.size_l / 2,  turn_table.size_l / 2, z4))

	# Add some stuff to mount to (merge it into the block):
	self.block(comment = "Bevel Gear Box Top Mount Block",
	  corner1 = P(gear_box.w.x + gear_box_cover.dx,
	    gear_box.s.y + gear_box_side.dy, z0),
	  corner2 = P(gear_box.e.x - gear_box_cover.dx,
	    gear_box.n.y - gear_box_side.dy, z3))

	# Do the twist timing pulley holes:
	angle_delta = Angle(360) / twist_timing_pulley.holes_count_i
	for index in range(twist_timing_pulley.holes_count_i):
	    angle = index * angle_delta
	    x = twist_timing_pulley.holes_radius_l.cosine(angle)
	    y = twist_timing_pulley.holes_radius_l.sine(angle)
	    self.hole(comment = "Hole {0}".format(index),
	      diameter = twist_timing_pulley.hole_diameter_l,
	      start = P(x, y, z4),
	      end = P(x, y, z2),
	      flags = "t")

	# Do the shaft hole:
	zero = L()
	self.hole(comment = "Shaft Hole",
	  diameter = twist_timing_pulley.shaft_diameter_l,
	  start = P(zero, zero, z4),
	  end = P(zero, zero, z0),
	  flags = "t")

	# Do the turn table mounting holes:
	for x_sign in [-1, 1]:
	    x = x_sign * turn_table.bottom_hole_pitch1_l / 2
	    for y_sign in [-1, 1]:
		y = y_sign * turn_table.bottom_hole_pitch1_l / 2
		self.hole(comment =
		  "Turn Table Hole [{0}, {1}]".format(x_sign, y_sign),
		  diameter = turn_table.bottom_hole_diameter1_l,
		  start = P(x, y, z0),
		  end = P(x, y, z4),
		  flags = "t")

	# Perform the top pocket:
	self.simple_pocket(comment = "Top Pocket",
	  top_corner = P(-pocket_dx / 2, -pocket_dy / 2, z0),
	  bottom_corner = P( pocket_dx / 2,  pocket_dy / 2, z2))

class Motor_Base(Part):
    def construct(self):
	#print("=>Motor_Base.construct()")
	motor_assembly = self.up
	synchro_drive = motor_assembly.up
	wheel_assembly = synchro_drive.wheel_assembly_

	turn_table = wheel_assembly.turn_table_
	turn_table_size = turn_table.size_l
	turn_table_bore = turn_table.bore_l

	self.dx_l = dx = turn_table_size
	self.dy_l = dy = turn_table_size
	self.dz_l = dz = L(mm = 3.0)

	z0 = motor_assembly.base_bz_l
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
	  end =   P(zero, zero, z0),
	  flags = "t")
	#print("<=Motor_Base.construct()")

class Shim(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):

	wheel_assembly = self.up
	bearing = wheel_assembly.south_bearing_
	bearing_diameter = bearing.diameter_l
	bearing_bore = bearing.bore_l

	self.bore_l = bore = bearing_bore
	self.diameter_l = diameter = (bearing_diameter + bore) / 2
	self.width_l = width = L(inch = .010)
	
	zero = L()
	self.cylinder(comment = "Shim Cylinder",
	  material = Material("steel", "stainless"),
	  color = Color("coral"),
	  diameter = diameter,
	  start = P(zero, zero, zero),
	  end = P(zero, zero, width))

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
	#quick = True
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

class Turn_Table(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	zero = L()
	self.turn_table_bottom_sheet_ = Turn_Table_Sheet(self, is_upper = False)
	self.turn_table_top_sheet_ = Turn_Table_Sheet(self, is_upper = True)

    def construct(self):
	# Provide some constants to work with:
	self.size_l = size = L(inch = 3.00)
	self.bore_l = bore = L(inch = 1.31)
	self.bottom_hole_pitch1_l = bottom_hole_pitch1 = L(inch = 2.56)
	self.bottom_hole_pitch2_l = bottom_hole_pitch2 = L(inch = 2.13)
	self.top_hole_pitch1_l = top_hole_pitch1 = L(inch = 2.56)
	self.bottom_hole_diameter1_l = bottom_hole_diameter1 = L(inch = 0.16)
	self.bottom_hole_diameter2_l = bottom_hole_diameter2 = L(inch = 0.094)
	self.top_hole_diameter1_l = top_hole_diameter1 = L(inch = 0.16)
	self.sheet_thickness_l = sheet_thickness = L(inch = 0.03)
	self.height_l = height = L(inch = 0.26)

class Turn_Table_Bottom_Screws(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.turn_table_bottom_screw0_ = Fastener(self)
	self.turn_table_bottom_screw1_ = Fastener(self)
	self.turn_table_bottom_screw2_ = Fastener(self)
	self.turn_table_bottom_screw3_ = Fastener(self)

    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	gear_box_top = wheel_assembly.gear_box_top_
	turn_table = wheel_assembly.turn_table_
	turn_table_bottom_sheet = turn_table.turn_table_bottom_sheet_

	# Grab the screws from *self*:
	screw0 = self.turn_table_bottom_screw0_
	screw1 = self.turn_table_bottom_screw1_
	screw2 = self.turn_table_bottom_screw2_
	screw3 = self.turn_table_bottom_screw3_

	# Compute various distances along Z axis:
	z0 = gear_box_top.b.z
	z1 = z0 + wheel_assembly.pocket_dz_l
	z2 = z1 + gear_box_top.base_dz_l
	z3 = z2 + turn_table.sheet_thickness_l

	# Do the the holes:
	for x_sign in [-1, 1]:
	    x = x_sign * turn_table.bottom_hole_pitch2_l / 2
	    for y_sign in [-1, 1]:
		y = y_sign * turn_table.bottom_hole_pitch2_l / 2

		# Select the *screw* based on *x_sign* and *y_sign*:
		if x_sign == -1 and y_sign == -1:
		    screw = screw0
		elif x_sign == -1 and y_sign == 1:
		    screw = screw1
		elif x_sign == 1 and y_sign == -1:
		    screw = screw2
		elif x_sign == 1 and y_sign == 1:
		    screw = screw3

		screw.configure(
		  comment = "[{0},{1}] Turn_Table_Screw".format(x_sign, y_sign),
		  diameter_pitch = "#2-56",
		  start = P(x, y, z1),
		  end = P(x, y, z3))

	# Now drill some holes in *gear_box_top*:
	screw0.drill(gear_box_top)
	screw1.drill(gear_box_top)
	screw2.drill(gear_box_top)
	screw3.drill(gear_box_top)

class Turn_Table_Sheet(Part):
    def __init__(self, up, is_upper = False):
	name = "Bottom_Turn_Table_Sheet"
	if is_upper:
	    name = "Top_Turn_Table_Sheet"
	Part.__init__(self, up, name = name)
	self.is_upper_b = is_upper
 
    def construct(self):
	#print("Turn_Table_Sheet:size={0} sheet_thickness={1}".
	#  format(size, sheet_thickness))

	# Grab some values from *turn_table*:
	turn_table = self.up
	turn_table_size = turn_table.size_l
	turn_table_bore = turn_table.bore_l
	turn_table_bottom_hole_pitch1 = turn_table.bottom_hole_pitch1_l
	turn_table_bottom_hole_pitch2 = turn_table.bottom_hole_pitch2_l
	turn_table_top_hole_pitch1 = turn_table.top_hole_pitch1_l
	turn_table_bottom_hole_diameter1 = turn_table.bottom_hole_diameter1_l
	turn_table_bottom_hole_diameter2 = turn_table.bottom_hole_diameter2_l
	turn_table_top_hole_diameter1 = turn_table.top_hole_diameter1_l
	turn_table_sheet_thickness = turn_table.sheet_thickness_l
	turn_table_height = turn_table.height_l

	# Make the original sheet:
	zero = L()
	is_upper = self.is_upper_b
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
	    x1 = x_sign * turn_table_bottom_hole_pitch1 / 2
	    x2 = x_sign * turn_table_bottom_hole_pitch2 / 2
	    x3 = x_sign * turn_table_top_hole_pitch1 / 2
	    for y_sign in [-1, 1]:
		y1 = y_sign * turn_table_bottom_hole_pitch1 / 2
		y2 = y_sign * turn_table_bottom_hole_pitch2 / 2
		y3 = y_sign * turn_table_top_hole_pitch1 / 2
		if is_upper:
		    self.hole(comment = "Top Hole1 [{0}, {1}]".
		      format(x_sign, y_sign),
		      diameter = turn_table_top_hole_diameter1,
		      start = P(x3, y3, z1),
		      end = P(x3, y3, z0),
		      flags = "t")
		else:
		    self.hole(comment = "Bottom Hole1 [{0}, {1}]".
		      format(x_sign, y_sign),
		      diameter = turn_table_bottom_hole_diameter1,
		      start = P(x1, y1, z1),
		      end = P(x1, y1, z0),
		      flags = "t")
		    self.hole(comment = "Bottom Hole2 [{0}, {1}]".
		      format(x_sign, y_sign),
		      diameter = turn_table_bottom_hole_diameter2,
		      start = P(x2, y2, z1),
		      end = P(x2, y2, z0),
		      flags = "t")

class Wheel(Part):
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
    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_
	gear_box_side = wheel_assembly.south_gear_box_side_
	wheel = wheel_assembly.wheel_

	# Some values from *bevel_gear*:
	bevel_gear = wheel_assembly.top_bevel_gear_
	bevel_gear_data = bevel_gear.bevel_gear_data_o
	bevel_gear_bore = bevel_gear_data.bore

	# Load up *self*:
	self.diameter_l = diameter = bevel_gear_bore

	# Build the shaft:
	zero = L()
	self.cylinder(comment = "Horizontal Shaft",
	  material = Material("Steel", ""),
	  color = Color("cyan"),
	  diameter = diameter,
	  start = P(zero, wheel.w.x, zero),
          end =   P(zero, gear_box.n.y - gear_box_side.bearing_lip_dy_l, zero))

class Vertical_Shaft(Part):
    def construct(self):
	wheel_assembly = self.up
	bevel_gear = wheel_assembly.top_bevel_gear_
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

