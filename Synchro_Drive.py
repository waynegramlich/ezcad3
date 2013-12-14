#!/usr/bin/env python

from EZCAD3 import *

class Synchro_Drive(Part):

    def __init__(self, up):
	Part.__init__(self, up)
	self.wheel_assembly_ = Synchro_Drive_Wheel_Assembly(self)

    def construct(self):
	pass

class Synchro_Drive_Wheel_Assembly(Part):
    def __init__(self, up):
	Part.__init__(self, up)

	# List all of the sub-*Part*'s that make up a wheel assembly:
        self.timing_pully_ = Timing_Pully(self, teeth_count = 30,
	  tooth_diameter = L(mm = 1.0), pitch = L(inch = 0.080),
	  width = L(inch = "1/4"))
	self.bearing_ = Bearing(self)
	self.bevel_gear_ = \
	  Bevel_Gear(self, part_name = "A 1M 4-Y16016", color = Color("red"))
	self.horizontal_shaft_ = Horizontal_Shaft(self)
	self.wheel_side_bevel_gear_box_side_ = \
	  Bevel_Gear_Box_Side(self, True)
	self.non_wheel_side_bevel_gear_box_side_ = \
	  Bevel_Gear_Box_Side(self, False)
	self.bevel_gear_box_shelf_ = Bevel_Gear_Box_Shelf(self)
	self.turn_table_ = Turn_Table(self)
	self.vertical_shaft_ = Vertical_Shaft(self)
	self.wheel_ = Wheel(self)
 
    def construct(self):
	# The origin of this assembly is at the logical intersection
	# of the vertical shaft and the horizontal shaft.  The horizontal
	# axis is aligned with the Y axis; this causes the wheel to be
	# oriented in the roll in the X direction.

	turn_table = self.turn_table_
	turn_table.no_automatic_place()

	bearing = self.bearing_
	bearing_width = bearing.width_l

	bevel_gear = self.bevel_gear_
	bevel_gear_width = bevel_gear.width_l
	bevel_gear_outside_diameter = bevel_gear.outside_diameter_l

	bevel_gear_box_side = self.wheel_side_bevel_gear_box_side_
	bevel_gear_box_side_bearing_lip_dy = \
	  bevel_gear_box_side.bearing_lip_dy_l

	# Define some constants:
	self.bevel_gear_box_dx_l = L(mm = 30.00)
	self.bevel_gear_box_dy_l = bevel_gear_box_dy = L(mm = 80.00)
	self.bevel_gear_box_nz_l = bevel_gear_box_nz = L(mm = 75.00)
	self.bevel_gear_box_sz_l = L(mm = -10.00)
	self.bevel_gear_box_shelf_sz_l = bevel_gear.nz_l

	# Place various items into the wheel assembly:
	zero = L()
	one = L(mm = 1.0)
	self.place(self.bearing_, name = "Wheel Side Bearing", 
	  axis = P(L(mm = 1.0), zero, zero), rotate = Angle(deg = 90),
	  translate = P(zero, -bevel_gear_box_dy / 2 +
	    bevel_gear_box_side_bearing_lip_dy + bearing_width / 2, zero))
	self.place(self.bearing_, name = "Non Wheel Side Bearing", 
	  axis = P(L(mm = 1.0), zero, zero), rotate = Angle(deg = 90),
	  translate = P(zero, bevel_gear_box_dy / 2 -
	    bevel_gear_box_side_bearing_lip_dy - bearing_width / 2, zero))
	self.place(bevel_gear, name = "Horizontal_Bevel_Gear", 
	  axis = P(L(mm = 1.0), zero, zero), rotate = Angle(deg = 90))
	self.place(self.turn_table_, name = "XTurn_Table",
	  translate = P(zero, zero, bevel_gear_box_nz))
	self.place(self.timing_pully_, name = "Turn Pully",
	  translate = P(zero, zero, L(mm = 110.0)))

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

class Bevel_Gear_Box_Shelf(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	wheel_assembly = self.up
	bevel_gear_box_shelf_sz = wheel_assembly.bevel_gear_box_shelf_sz_l
	bevel_gear_box_dx = wheel_assembly.bevel_gear_box_dx_l
	bevel_gear_box_dy = wheel_assembly.bevel_gear_box_dy_l

	bearing = wheel_assembly.bearing_
	bearing_diameter = bearing.diameter_l
	bearing_width = bearing.width_l

	bevel_gear_box_side = wheel_assembly.wheel_side_bevel_gear_box_side_
	bevel_gear_box_side_thin_dy = bevel_gear_box_side.thin_dy_l

	self.bearing_lip_dz_l = bearing_lip_dz = L(mm = 2.00)
	self.dz_l = dz = bearing_width + bearing_lip_dz
	
	corner1 = P(-bevel_gear_box_dx / 2,
	  -bevel_gear_box_dy / 2 + bevel_gear_box_side_thin_dy,
	  bevel_gear_box_shelf_sz)
	corner2 = P(bevel_gear_box_dx / 2,
	  bevel_gear_box_dy / 2 - bevel_gear_box_side_thin_dy,
	  bevel_gear_box_shelf_sz + dz)
	self.block(comment = "Shelf Block",
	  material = Material("plastic", "ABS"),
	  color = Color("gold"),
	  corner1 = corner1,
	  corner2 = corner2)

class Bevel_Gear_Box_Side(Part):
    def __init__(self, up, wheel_side):
	# Deal with *wheel_side* argument:
	assert isinstance(wheel_side, bool)
	if wheel_side:
	    name_prefix = "Wheel_Side"
	else:
	    name_prefix = "Non_Wheel_Side"
	self.wheel_side_b = wheel_side

	Part.__init__(self, up,
	  name = name_prefix + "_Bevel_Gear_Box_Side")

    def construct(self):
	""" *Bevel_Gear_Box_Side*: """

	wheel_side = self.wheel_side_b

	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	bearing = wheel_assembly.bearing_

	# Grab the basic bevel gear box dimensions from *wheel_assembly*:
	bevel_gear_box_dx = wheel_assembly.bevel_gear_box_dx_l
	bevel_gear_box_dy = wheel_assembly.bevel_gear_box_dy_l
	bevel_gear_box_nz = wheel_assembly.bevel_gear_box_nz_l
	bevel_gear_box_sz = wheel_assembly.bevel_gear_box_sz_l
	bevel_gear_box_shelf_sz = wheel_assembly.bevel_gear_box_shelf_sz_l
	
	# Grab values from *bearing*:
	bearing_diameter = bearing.diameter_l
	bearing_width = bearing.width_l

	# Grab values from *bevel_gear_box_shelf*:
	bevel_gear_box_shelf = wheel_assembly.bevel_gear_box_shelf_
	bevel_gear_box_shelf_dz = bevel_gear_box_shelf.dz_l

	horizontal_shaft = wheel_assembly.horizontal_shaft_
	horizontal_shaft_diameter = horizontal_shaft.diameter_l

	# Figure out the various dimensions of the shelf:
	self.shelf_nub_dy_l = shelf_nub_dy = L(mm = 5.00)
	self.shelf_nub_dz_l = shelf_nub_dz = L(mm = 5.00)
	self.shelf_nub_sz_l = shelf_nub_sz = \
	  bevel_gear_box_shelf_sz + bevel_gear_box_shelf_dz
	self.shelf_nub_nz_l = shelf_nub_nz = shelf_nub_sz + shelf_nub_dz

	self.bearing_lip_dy_l = bearing_lip_dy = L(mm = 2.00)
	self.thin_dy_l = thin_dy = bearing_width + bearing_lip_dy
	self.thick_dy_l = thick_dy = thin_dy + shelf_nub_dy
		
	# Identify various points on the Y axis:
	y1 = bevel_gear_box_dy / 2
	y2 = y1 - bearing_lip_dy
	y3 = y1 - thin_dy
	y4 = y1 - thick_dy
	if wheel_side:
	    y1 = -y1
	    y2 = -y2
	    y3 = -y3
	    y4 = -y4

	corner1 = \
	  P(-bevel_gear_box_dx / 2, y1, bevel_gear_box_sz)
	corner2 = \
	  P( bevel_gear_box_dx / 2, y3, bevel_gear_box_nz)
	self.block(comment = self._name + " Main Block",
	  material = Material("plastic", "abs"),
	  color = Color("orange"),
	  corner1 = corner1,
	  corner2 = corner2)

	corner1 = \
	  P(-bevel_gear_box_dx / 2, y3, shelf_nub_sz)
	corner2 = \
	  P( bevel_gear_box_dx / 2, y4, shelf_nub_nz)
	self.block(comment = self._name + " Shelf Block",
	  material = Material("plastic", "abs"),
	  color = Color("orange"),
	  corner1 = corner1,
	  corner2 = corner2)

	zero = L()
	self.hole(comment = "Bearing hole",
	  diameter = bearing_diameter,
	  start = P(zero, y3, zero),
	  end = P(zero, y2, zero),
	  flags = "f")

	if wheel_side:
	    diameter = (horizontal_shaft_diameter + bearing_diameter) / 2
	    self.hole(comment = "Shaft hole",
	      diameter = diameter,
	      start = P(zero, y3, zero),
	      end = P(zero, y1, zero),
	      flags = "t")

class Timing_Pully(Part):
    def __init__(self, up, teeth_count = -1,
      pitch = L(), width = L(), tooth_diameter = L()):
	Part.__init__(self, up)

	# Check argument types:
	zero = L()
      	assert isinstance(teeth_count, int)
	assert isinstance(pitch, L) and pitch > zero
	assert isinstance(tooth_diameter, L) and tooth_diameter > zero

	# Remember the arguments;
	self.teeth_count = teeth_count
	self.tooth_diameter = tooth_diameter
	self.pitch = pitch
	self.width = width

    def construct(self):
	self.no_automatic_place()

	teeth_count = self.teeth_count
	tooth_diameter = self.tooth_diameter
	pitch = self.pitch
	width = self.width

	# The angle between each tooth
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

	tooth_radius = pitch / (2 * (tooth_angle / 2).sine())
	#print("tooth_angle={0:d}, tooth_angle.sine()={1}".
	#  format(tooth_angle, (tooth_angle / 2).sine()))
	#print("pitch={0}, tooth_radius={1} circum={2} p*tr={3}".
	#  format(pitch, tooth_radius,
	#  tooth_radius * 2 * 3.1415629, pitch * teeth_count))

	zero = L()
	self.cylinder(comment = "Pully Body",
	  material = Material("plastic", "ABS"),
	  color = Color("crimson"),
	  diameter = tooth_radius * 2,
	  start = P(zero, zero, zero),
	  end = P(zero, zero, width),
	  sides = teeth_count)

	for index in range(teeth_count):
	    angle = tooth_angle * index
	    x = tooth_radius.cosine(angle)
	    y = tooth_radius.sine(angle)
	    self.hole(comment = "Tooth {0}".format(index),
	      diameter = tooth_diameter,
	      start = P(x, y, zero),
	      end = P(x, y, width),
	      flags = "t")

class Turn_Table(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	zero = L()
	self.turn_table_sheet_ = Turn_Table_Sheet(self)

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

	zero = L()
	self.place(self.turn_table_sheet_, name = "Top Sheet",
	  translate = P(zero, zero, height - sheet_thickness))

class Turn_Table_Sheet(Part):
    def __init__(self, up):
	Part.__init__(self, up)
 
    def construct(self):
	turn_table = self.up
	size = turn_table.size_l
	sheet_thickness = turn_table.sheet_thickness_l
	bore = turn_table.bore_l
	
	#print("Turn_Table_Sheet:size={0} sheet_thickness={1}".
	#  format(size, sheet_thickness))

	zero = L()
	self.block(comment = "Turn Table Sheet",
	  material = Material("Steel", "Galvanized"),
	  color = Color("purple"),
	  corner1 = P(-size / 2, -size / 2, zero),
	  corner2 = P( size / 2,  size / 2, sheet_thickness))
	self.hole(comment = "Bore Hole",
	  diameter = bore,
	  start = P(zero, zero, zero),
	  end =   P(zero, zero, sheet_thickness),
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
	wheel_assembly = self.up

	wheel = wheel_assembly.wheel_

	bevel_gear = wheel_assembly.bevel_gear_
	bevel_gear_data = bevel_gear.bevel_gear_data_o
	bore = bevel_gear_data.bore

	bevel_gear_box_dy = wheel_assembly.bevel_gear_box_dy_l

	bevel_gear_box_side = wheel_assembly.wheel_side_bevel_gear_box_side_
	bevel_gear_box_side_bearing_lip_dy = \
	  bevel_gear_box_side.bearing_lip_dy_l


	zero = L()
	self.cylinder(comment = "Horizontal Shaft",
	  material = Material("Steel", ""),
	  color = Color("cyan"),
	  diameter = bore,
	  start = P(zero, wheel.w.x, zero),
          end = P(zero,
	    bevel_gear_box_dy / 2 - bevel_gear_box_side_bearing_lip_dy, zero))

class Vertical_Shaft(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	bevel_gear = self.up.bevel_gear_
	bevel_gear_data = bevel_gear.bevel_gear_data_o
	bore = bevel_gear_data.bore
	major_diameter = bevel_gear_data.major_diameter
	height = bevel_gear_data.height

	zero = L()
	self.cylinder(comment = "Horizontal Shaft",
	  material = Material("Steel", ""),
	  color = Color("cyan"),
	  diameter = bore,
	  start = P(zero, zero, major_diameter - height),
          end = P(zero, zero, L(mm = 100.0)))

if __name__ == "__main__":
    ezcad = EZCAD3(0)
    test = Synchro_Drive(None)
    test.process(ezcad)

