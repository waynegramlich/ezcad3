#!/usr/bin/env python

from EZCAD3 import *

# Assemblies:

class Base(Part):
    def __init__(self, up):
	Part.__init__(self,up)
	self.left_front_synchro_ = Synchro_Drive(self)
	self.right_front_synchro_ = Synchro_Drive(self)

    def construct(self):
	left_front_synchro = self.left_front_synchro_
	right_front_synchro = self.right_front_synchro_

	one = L(mm = 1.00)
	zero = L()
	self.twist_axis_dy_l = twist_axis_dy = L(inch = 16.00)
	z_axis = P(zero, zero, one)

	front_right = P(zero, -twist_axis_dy / 2, zero)
	front_left = P(zero, twist_axis_dy / 2, zero)
	left_front_synchro.place(translate = front_right)
	right_front_synchro.place(translate = front_left,
	  axis = z_axis, rotate = Angle(180))

class Strut_Assembly(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	#self.strut_ = Strut(self)
	self.strut_top_ = Strut_Top(self)

    def construct(self):
	# Grab some relevant *Part*'s:
	synchro_drive = self.up
	motor_assembly = synchro_drive.motor_assembly_
	motor_base = motor_assembly.motor_base_

	# Grab coordnates for point J:
	j_x = motor_base.j_x_l
	j_y = motor_base.j_y_l

class Strut(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	# Grab some relevant *Part*'s:
	strut_assembly = self.up
	synchro_drive = strut_assembly.up
	strut_top = strut_assembly.strut_top_
	motor_assembly = synchro_drive.motor_assembly_

	y0 = motor_assembly.n.y
	y1 = strut_top.n.y
	x0 = L(inch = "1/4")
	x1 = -L(inch = "1/4")

	z0 = motor_assembly.t.z
	z1 = z0 + L(inch = 0.125)

	corner1 = P(x0, y0, z0)
	corner2 = P(x1, y1, z1)
	print("Strut: c1={0} c2={1}".format(corner1, corner2))
	self.block(comment = "Strut",
	  color = Color("violet"),
	  corner1 = corner1,
	  corner2 = corner2,
	  top = "t")

class Strut_Top(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	# Grab some *Part*'s from the *Part* tree:
	strut_assembly = self.up
	synchro = strut_assembly.up
	base = synchro.up
	motor_assembly = synchro.motor_assembly_
	motor_base = motor_assembly.motor_base_
	
	# Grab distance between the two twist axes from *base*:
	twist_axis_dy = base.twist_axis_dy_l

	# Point J and/or M in *motor_base* of the motor assembly
	# is the closest that strut assembly can get to its partner
	# on the other side.  Points Q and R are mirror image of
	# points J and M across the X axis line that goes through
	# the pair origin O.
	#
	#                  +-------+
	#                 /        |
	#                /     o   |
	#               /          |
	#              /          /
	#             /          /
	#            /          /
	#           /          /
	#          M          J Q----------R
	#          |          | | Right    |
	#          |          | | Strut    |
	#          | Left     | | Assembly |
	#          | Motor    | T----------S
	#          | Assembly |  ||      ||
	#          +----------+  ||      ||
	#           ||      ||   ||      ||
	#           ||      || O ||      ||
	#           ||      ||   ||      ||
	#           ||      ||  +----------+
	#           ||      ||  | Right    |
	#          R----------Q | Motor    |
	#          | Left     | | Assembly |
	#          | Strut    | |          |
	#          | Assembly | |          |
	#          S----------T J          M
	#                      /          /
	#                     /          /
	#                    /          /
	#                   /          /
	#                  |          /
	#                  |   o     /
	#                  |        /
	#                  +-------+
	#
	# +Y
	#  ^
	#  |
	#  O----> +X
	#
	# where:
	#    "O"  is the origin of the *Synchro_Pair*.
	#    "o"  is the twist axis for each motor assembly.  Note that
	#         the twist axis aligns in Y with "O".
	#    "||" Is a strut pair (upper and lower.)
	#    "J"  is an inside point for the motor base assembly.
	#    "M"  is an outside point for the motor base assembly,
	#         where Jy = My on the same motor assembly.
	#    "Q"  Mirror of "J"
	#    "R"  Mirror of "M"
	
	# Grab X/Y coordinates of J and M from *motor_base*:
	j_x = motor_base.j_x_l
	j_y = motor_base.j_y_l
	m_x = motor_base.m_x_l
	m_y = motor_base.m_y_l

	# Length of the top strut in Y:
	self.dy_l = dy = L(mm = 50.00)

	# Compute the X/Y locations of Q, R, S, and T.  Remember, the
	# origin coordinates are still around the sychro twist axis.
	self.q_x_l = q_x = j_x
	self.q_y_l = q_y = twist_axis_dy - j_y
	self.r_x_l = r_x = m_x
	self.r_y_l = r_y = twist_axis_dy - m_y
	self.s_x_l = s_x = r_x
	self.s_y_l = s_y = r_y - dy
	self.t_x_l = t_x = q_x
	self.t_y_l = t_y = q_y - dy

	z0 = motor_base.b.z
	z3 = motor_base.t.z

	self.block(comment = "Strut Top Block",
	  color = Color("orange"),
	  corner1 = P(q_x, q_y, z0),
	  corner2 = P(s_x, s_y, z3),
	  top = "t")

class Synchro_Drive(Part):

    def __init__(self, up):
	Part.__init__(self, up)
	self.wheel_assembly_ = Wheel_Assembly(self)
	self.motor_assembly_ = Motor_Assembly(self)
	self.strut_assembly_ = Strut_Assembly(self)

    def construct(self):
	zero = L()
	one = L(mm = 1.00)
	motor_assembly = self.motor_assembly_
	motor_assembly.place(
	  axis = P(zero, zero, one))
	#  rotate = Angle(deg = 130))

class Motor_Assembly(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.motor_base_ = Motor_Base(self)
	self.motor_base_cover_ = Motor_Base_Cover(self)
	self.motor_base_screws_ = Motor_Base_Screws(self)
	self.drive_belt_ = Belt(self)
	self.drive_magnet_ = Magnet(self)
	self.drive_motor_ = Motor(self)
	self.drive_motor_pulley_ = Pulley(self)
	self.drive_motor_pulley_screws_ = Drive_Motor_Pulley_Screws(self)
	self.drive_motor_pulley_top_ = Pulley_Top(self)
	self.end_screws_ = End_Screws(self)
	self.pcb_ = PCB(self)
	self.twist_belt_ = Belt(self)
	self.twist_magnet_ = Magnet(self)
	self.twist_motor_ = Motor(self)
	self.twist_motor_pulley_ = Pulley(self)
	self.twist_motor_pulley_screws_ = Twist_Motor_Pulley_Screws(self)
	self.twist_motor_pulley_top_ = Pulley_Top(self)

    def construct(self):
	# Grab some values from *self*:
	motor_assembly = self
	drive_belt = motor_assembly.drive_belt_
	drive_magnet = motor_assembly.drive_magnet_
	drive_motor = motor_assembly.drive_motor_
	drive_motor_pulley = motor_assembly.drive_motor_pulley_
	drive_motor_pulley_top = motor_assembly.drive_motor_pulley_top_
	#drive_motor_pulley_top = motor_assembly.drive_motor_pulley_top_
	motor_base = motor_assembly.motor_base_
	motor_base_cover = motor_assembly.motor_base_cover_
	pcb = motor_assembly.pcb_
	synchro_drive = self.up
	twist_magnet = motor_assembly.twist_magnet_
	twist_motor = motor_assembly.twist_motor_
	twist_motor_pulley = motor_assembly.twist_motor_pulley_
	twist_motor_pulley_top = motor_assembly.twist_motor_pulley_top_
	twist_belt = motor_assembly.twist_belt_

	# Define some constants:
	zero = L()
	sqrt2 = math.sqrt(2)

	motor_base_cover.invisible_set()

	# Grab some values from *wheel_assembly:
	wheel_assembly = synchro_drive.wheel_assembly_
	gear_box = wheel_assembly.gear_box_
	turn_table = wheel_assembly.turn_table_
	twist_wheel_pulley = wheel_assembly.twist_wheel_pulley_
	drive_wheel_pulley = wheel_assembly.drive_wheel_pulley_
	wheel = wheel_assembly.wheel_

	# If the wheel radius is R and its witdh is W, (R+W/2)*sqrt(2)
	# is a pretty close approximation to the radius of the *wheel_sweep*,
	# which is the circle that entirly encloses the wheel location:
	wheel_sweep = (wheel.diameter_l / 2 + wheel.width_l / 2) * sqrt2
	self.motor_shaft_radius_l = motor_shaft_radius = \
	  wheel_sweep + twist_motor.spur_gear_body_diameter_l / 2 - \
	  twist_motor.shaft_offset_l

	# The timing belts are basically specified by the number of teeth
	# on the belt.  Furthermore, the number of teeth seems to be a
	# multiple of 5.  We start by computing the desired distance between
	# the wheel vertical shaft and the twist motor shaft.  Then
	# that is reduced to a *teeth_radius*.  Next, we compute the minimum
	# of *belt_teeth*.  Next, we round up to a multiple of 5.  Finally,
	# we recompute the final updated *motor_shaft_radius*:
	pulley_teeth = twist_motor_pulley.teeth_count_i
	pulley_pitch = twist_motor_pulley.pitch_l
	teeth_radius = int(math.ceil(motor_shaft_radius / pulley_pitch))
	#print("motor_shaft_radius={0} teeth_radius={1} pulley_teeth={2}".
	#  format(motor_shaft_radius, teeth_radius, pulley_teeth))
	belt_teeth = pulley_teeth + 2 * teeth_radius
	while belt_teeth % 5 != 0:
	    belt_teeth += 1	# Round up
	motor_shaft_radius = \
	  pulley_pitch * float((belt_teeth - pulley_teeth) / 2)
	#print("belt_teeth={0} motor_shaft_radius={1}".
	# format(belt_teeth, motor_shaft_radius))

	# Normally the position of the twist and drive motor shafts
	# is computed using polar coordinates from the main shaft.
	# All of these various positions of everything are written
	# into a file called pcb_holes.txt which is shown below:
	#
	# pcb corner1: pcb_x=150.0 pcb_y=150.0 x=5.0 y=70.0
	# pcb corner2: pcb_x=50.0 pcb_y=50.0 x=105.0 y=170.0
	# drive_shaft: pcb_x=107.890135052 pcb_y=121.189510914
	#   x=47.1098649476 y=98.810489086
	# twist_shaft: pcb_x=77.4986166491 pcb_y=140.925993734
	#   x=77.5013833509 y=79.0740062661
	# J #4-40:: pcb_x=147.0 pcb_y=145.0 x=8.0 y=75.0
	# K #4-40:: pcb_x=147.0 pcb_y=100.0 x=8.0 y=120.0
	# L #4-40:: pcb_x=53.0 pcb_y=100.0 x=102.0 y=120.0
	# M #4-40:: pcb_x=53.0 pcb_y=147.0 x=102.0 y=73.0
	# N #6-32: pcb_x=143.65 pcb_y=90.65 x=11.35 y=129.35
	# O #6-32: pcb_x=143.65 pcb_y=56.35 x=11.35 y=163.65
	# P #6-32: pcb_x=56.35 pcb_y=90.65 x=98.65 y=129.35
	# Q #6-32: pcb_x=56.35 pcb_y=56.35 x=98.65 y=163.65
	#
	# After the PCB has been manufactured using these numbers,
	# everything has to be made relative to the M corner on
	# the Motor_Base.  The M corner is (PCBc2.x-3, PCBc1.y+3).

	motor_base_m_x = motor_base.m_x_l
	motor_base_m_y = motor_base.m_y_l
	pcb_corner1_x = motor_base_m_x - L(mm=102) + L(mm=5) - L(mm=3)
	pcb_corner1_y = motor_base_m_y - L(mm=73) + L(mm=70) + L(mm=3)

	# Offsets from PCB corner1 to twist and drive motor shaft:
	twist_motor_dx = L(mm=77.5013833509 - 5)
	twist_motor_dy = L(mm=79.0740062661 - 70)
	drive_motor_dx = L(mm=47.1098649476 - 5)
	drive_motor_dy = L(mm=98.8104890861 - 70)

	self.twist_motor_x_l = twist_motor_x = pcb_corner1_x + twist_motor_dx
	self.twist_motor_y_l = twist_motor_y = pcb_corner1_y + twist_motor_dy
	self.drive_motor_x_l = drive_motor_x = pcb_corner1_x + drive_motor_dx
	self.drive_motor_y_l = drive_motor_y = pcb_corner1_y + drive_motor_dy
	  
	#print("twist_motor=({0},{1}".format(twist_motor_x, twist_motor_y))
	#print("drive_motor=({0},{1}".format(drive_motor_x, drive_motor_y))

	#self.twist_motor_angle_a = twist_motor_angle = Angle(deg = 66)
	#self.drive_motor_angle_a = drive_motor_angle = Angle(deg = 48)
	self.twist_motor_angle_a = twist_motor_angle = Angle(deg = 48)
	self.drive_motor_angle_a = drive_motor_angle = Angle(deg = 66)

	# Use this code *before* PCB is manufactured:
	# Compute the twist/drive motor shaft locations:
	#self.twist_motor_x_l = twist_motor_x = \
	#  motor_shaft_radius.cosine(twist_motor_angle)
	#self.twist_motor_y_l = twist_motor_y = \
	#  motor_shaft_radius.sine(twist_motor_angle)
	#self.drive_motor_x_l = drive_motor_x = \
	#  motor_shaft_radius.cosine(drive_motor_angle)
	#self.drive_motor_y_l = drive_motor_y = \

	self.drive_motor_shaft_x_l = drive_motor_shaft_x = drive_motor_x
	self.drive_motor_shaft_y_l = drive_motor_shaft_y = \
	  drive_motor_y - drive_motor.shaft_offset_l
	self.twist_motor_shaft_x_l = twist_motor_shaft_x = twist_motor_x
	self.twist_motor_shaft_y_l = twist_motor_shaft_y = \
	  twist_motor_y - twist_motor.shaft_offset_l

	# Configure *drive_belt*:
	drive_wheel_pulley_lip_bz = \
	  drive_wheel_pulley.b.z + drive_wheel_pulley.lip_height_l
	drive_wheel_pulley_lip_tz = drive_wheel_pulley.t.z
	drive_wheel_pulley_lip_cz = \
	  (drive_wheel_pulley_lip_bz + drive_wheel_pulley_lip_tz) / 2
	drive_belt_height = L(inch = "1/8")
	drive_belt.configure(
	  r1 = drive_wheel_pulley.tooth_radius_l,
	  x1 = zero,
	  y1 = zero,
	  z1 = drive_wheel_pulley_lip_cz - drive_belt_height / 2,
	  r2 = drive_wheel_pulley.tooth_radius_l,
	  x2 = drive_motor_shaft_x,
	  y2 = drive_motor_shaft_y,
	  z2 = drive_wheel_pulley_lip_cz + drive_belt_height / 2,
	  width = drive_wheel_pulley.tooth_height_l)

	# Configure *drive_magnet*:
	drive_magnet_dz = L(inch = "1/8")
	drive_magnet.configure(
	  diameter = L(inch = "1/4"),
	  dz = drive_magnet_dz,
	  x = drive_motor_shaft_x,
	  y = drive_motor_shaft_y,
	  z_bottom = pcb.magnet_z_top_l - drive_magnet_dz)

	# Configure *drive_motor_pulley*:
	drive_motor_pulley_lip_height = \
	  (drive_wheel_pulley.t.z - drive_wheel_pulley.belt_width_l - \
	  drive_wheel_pulley.belt_width_extra_l) - twist_motor_pulley.b.z
	drive_motor_pulley.configure(
	  belt_class = drive_wheel_pulley.belt_class_s,
	  belt_width = drive_wheel_pulley.belt_width_l,
	  belt_width_extra = drive_wheel_pulley.belt_width_extra_l,
	  holes_radius = L(mm = 10.00),
	  lip_extra = drive_wheel_pulley.lip_extra_l,
	  lip_height = drive_motor_pulley_lip_height,
	  name = "Drive_Motor_Pulley",
	  set_screw_dz = drive_motor_pulley_lip_height / 2,
	  shaft_diameter = drive_motor.shaft_diameter_l,
	  teeth_count = drive_wheel_pulley.teeth_count_i,
	  x = drive_motor_shaft_x,
	  y = drive_motor_shaft_y,
	  z_bottom = twist_motor_pulley.b.z)

	# Configure *drive_motor_pulley_top*:
	drive_motor_pulley_top_lip_dz = L(mm = 6.00)
	drive_motor_pulley_top_z_bottom = drive_motor_pulley.t.z
	drive_motor_pulley_top_hub_dz = pcb.magnet_z_top_l - \
	  drive_motor_pulley_top_z_bottom - drive_motor_pulley_top_lip_dz
	drive_motor_pulley_top.configure(
	  hub_diameter = drive_motor.shaft_diameter_l + L(mm = 6),
	  hub_dz = drive_motor_pulley_top_hub_dz,
	  lip_diameter = drive_wheel_pulley.dx,
	  lip_dz = drive_motor_pulley_top_lip_dz,
	  magnet_diameter = drive_magnet.dx,
	  magnet_dz = drive_magnet.dz,
	  shaft_diameter = drive_motor.shaft_diameter_l,
	  x = drive_motor_shaft_x,
	  y = drive_motor_shaft_y,
	  z_bottom = drive_motor_pulley.t.z)

	# Do *pcb* configuration:
	pcb_x1 = motor_base.m_x_l - L(100)
	pcb_y1 = motor_base.m_y_l
	pcb_dz = L(mm = 1.6)
	#pcb_x1 = L(mm = 5)
	#pcb_y1 = L(mm = 70)
	pcb.configure(x1 = pcb_x1, y1 = pcb_y1,
	  x2 = pcb_x1 + L(mm = 100), y2 = pcb_y1 + L(mm = 100),
	  dz = pcb_dz, z_bottom = motor_base.t.z - pcb_dz)
	#print("pcb.bsw={0} pcb.tne={1}".format(pcb.bsw, pcb.tne))

	# Do the *twist_belt* configuration:
	twist_wheel_pulley_lip_bz = \
	  twist_wheel_pulley.b.z + twist_wheel_pulley.lip_height_l
	twist_wheel_pulley_lip_tz = twist_wheel_pulley.t.z
	twist_wheel_pulley_lip_cz = \
	  (twist_wheel_pulley_lip_bz + twist_wheel_pulley_lip_tz) / 2
	twist_belt_height = L(inch = "1/8")
	twist_belt.configure(
	  r1 = twist_wheel_pulley.tooth_radius_l,
	  x1 = zero,
	  y1 = zero,
	  z1 = twist_wheel_pulley_lip_cz - twist_belt_height / 2,
	  r2 = twist_wheel_pulley.tooth_radius_l,
	  x2 = twist_motor_shaft_x,
	  y2 = twist_motor_shaft_y,
	  z2 = twist_wheel_pulley_lip_cz + twist_belt_height / 2,
	  width = twist_wheel_pulley.tooth_height_l)

	# Configure *twist_magnet*:
	twist_magnet_dz =  L(inch = "1/8")
	twist_magnet.configure(
	  diameter = L(inch = "1/4"),
	  dz = twist_magnet_dz,
	  x = twist_motor_shaft_x,
	  y = twist_motor_shaft_y,
	  z_bottom = pcb.magnet_z_top_l - twist_magnet_dz)

	# Do the *twist_motor* configuration:
	twist_motor.configure(
	  x = twist_motor_x,
	  y = twist_motor_y,
	  z = motor_base.b.z,
	  spur_gear_body_dz = L(mm = 27.00))
	drive_motor.configure(
	  x = drive_motor_x,
	  y = drive_motor_y,
	  z = motor_base.b.z,
	  spur_gear_body_dz = L(mm = 18.50))

	# Do the *twist_motor_pulley* configuration:
	twist_motor_pulley.configure(
	  belt_class = twist_wheel_pulley.belt_class_s,
	  belt_width = twist_wheel_pulley.belt_width_l,
	  belt_width_extra = twist_wheel_pulley.belt_width_extra_l,
	  lip_extra = twist_wheel_pulley.lip_extra_l,
	  lip_height = L(mm = 2.00),
	  name = "Twist_Motor_Pulley",
	  holes_radius = L(mm = 10.00),
	  shaft_diameter = twist_motor.shaft_diameter_l,
	  teeth_count = twist_wheel_pulley.teeth_count_i,
	  x = twist_motor_shaft_x,
	  y = twist_motor_shaft_y,
	  z_bottom = motor_base.b.z + motor_base.floor_dz_l + L(mm = 0.50))

	# Configure *twist_motor_pulley_top*:
	twist_motor_pulley_top_lip_dz = L(mm = 8.00)
	twist_motor_pulley_top_z_bottom = twist_motor_pulley.t.z
	twist_motor_pulley_top_hub_dz = pcb.magnet_z_top_l - \
	  twist_motor_pulley_top_z_bottom - twist_motor_pulley_top_lip_dz
	twist_motor_pulley_top.configure(
	  hub_diameter = twist_motor.shaft_diameter_l + L(mm = 6),
	  hub_dz = twist_motor_pulley_top_hub_dz,
	  lip_diameter = twist_wheel_pulley.dx,
	  lip_dz = twist_motor_pulley_top_lip_dz,
	  magnet_diameter = twist_magnet.dx,
	  magnet_dz = twist_magnet.dz,
	  set_screw_dz = L(mm = 4.00),
	  shaft_diameter = twist_motor.shaft_diameter_l,
	  x = twist_motor_shaft_x,
	  y = twist_motor_shaft_y,
	  z_bottom = twist_motor_pulley.t.z)

	zero = L()
	self.base_bz_l = base_bz = gear_box.t.z + turn_table.height_l


class Wheel_Assembly(Part):
    def __init__(self, up):
	Part.__init__(self, up)

	# List all of the sub-*Part*'s that make up a wheel assembly:

	self.bottom_bearing_ = Bearing(self)
	self.drive_wheel_pulley_ = Pulley(self)
	self.drive_wheel_pulley_top_ = Pulley_Top(self)
	self.drive_wheel_pulley_screws_ = Drive_Wheel_Pulley_Screws(self)
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
	self.twist_wheel_pulley_ = Pulley(self)
	self.twist_wheel_pulley_top_ = Pulley_Top(self)
	self.twist_wheel_pulley_screws_ = Twist_Wheel_Pulley_Screws(self)
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
	bottom_bearing = wheel_assembly.bottom_bearing_
	drive_wheel_pulley = wheel_assembly.drive_wheel_pulley_
	drive_wheel_pulley_top = wheel_assembly.drive_wheel_pulley_top_
	gear_box = wheel_assembly.gear_box_
	gear_box_cover = wheel_assembly.west_gear_box_cover_
	gear_box_bottom = wheel_assembly.gear_box_bottom_
	gear_box_bottom_screws = wheel_assembly.gear_box_bottom_screws_
	gear_box_side = wheel_assembly.south_gear_box_side_
	gear_box_shelf = wheel_assembly.gear_box_shelf_
	gear_box_shelf_screws = wheel_assembly.gear_box_shelf_screws_
	gear_box_top = wheel_assembly.gear_box_top_
	gear_box_top_screws = wheel_assembly.gear_box_top_screws_
	north_bearing = wheel_assembly.north_bearing_
	south_bearing = wheel_assembly.south_bearing_
	south_bevel_gear = wheel_assembly.south_bevel_gear_
	south_shim = wheel_assembly.south_shim_
	top_bearing = wheel_assembly.top_bearing_
	top_bevel_gear = wheel_assembly.top_bevel_gear_
	top_shim = wheel_assembly.top_shim_
	turn_table = wheel_assembly.turn_table_
	twist_wheel_pulley = wheel_assembly.twist_wheel_pulley_
	twist_wheel_pulley_top = wheel_assembly.twist_wheel_pulley_top_
	vertical_shaft = wheel_assembly.vertical_shaft_

	synchro_drive = wheel_assembly.up
	motor_assembly = synchro_drive.motor_assembly_
	motor_base = motor_assembly.motor_base_
	twist_motor_pulley = motor_assembly.twist_motor_pulley_

	belt_width_extra = L(inch=.012)

	lip_height = (twist_motor_pulley.t.z - \
	  twist_motor_pulley.belt_width_l - \
	  twist_motor_pulley.belt_width_extra_l) - turn_table.b.z
	twist_wheel_pulley.configure(
	  bearing_diameter = top_bearing.diameter_l,
	  bearing_sides = top_bearing.sides_i,
	  bearing_width = top_bearing.width_l,
	  belt_class = "MXL",
	  belt_width = L(inch = "1/8"),
	  belt_width_extra = belt_width_extra,
	  lip_extra = L(mm = 0.5),
	  lip_height = lip_height,
	  name = "Twist_Wheel_Pulley",
	  shaft_diameter = L(inch = "1/2"),
	  teeth_count = 46,
	  z_bottom = gear_box.t.z)

	twist_wheel_pulley_top.configure(
	  lip_diameter = twist_wheel_pulley.lip_diameter_l,
	  lip_dz = L(mm = 2.00),
	  shaft_diameter = (top_bearing.diameter_l + top_bearing.bore_l) / 2,
	  z_bottom = twist_wheel_pulley.t.z)

	drive_wheel_pulley.configure(
	  belt_class = "MXL",
	  belt_width = L(inch = "1/8"),
	  belt_width_extra = belt_width_extra,
	  holes_radius = L(mm = 10.00),
	  lip_extra = L(mm = 0.50),
	  lip_height = L(mm = 2.00),
	  name = "Drive_Wheel_Pulley",
	  shaft_diameter = vertical_shaft.diameter_l,
	  teeth_count = 46,
	  z_bottom = twist_wheel_pulley_top.t.z + L(mm = 1.00))

	#  hub_diameter = vertical_shaft.dx + L(mm = 10),
	drive_wheel_pulley_top.configure(
	  hub_diameter = drive_wheel_pulley.dx,
	  hub_dz = L(mm = 6.00),
	  lip_diameter = drive_wheel_pulley.dx,
	  lip_dz = L(mm = 2.00),
	  set_screw_dz = L(mm = 4.00),
	  shaft_diameter = vertical_shaft.dx,
	  z_bottom = drive_wheel_pulley.t.z)

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
	self.pocket_wall_dx_l = pocket_wall_dx = L(mm = 4.00)
	self.pocket_wall_dy_l = pocket_wall_dy = L(mm = 4.00)
	self.pocket_dx_l = \
	  gear_box.dx - 2 * gear_box_cover.dx -  2 * pocket_wall_dx
	self.pocket_dy_l = \
	  gear_box.dy - 2 * gear_box_side.dy - 2 * pocket_wall_dy
	self.pocket_dz_l = L(mm = 8.00)

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
	p = P(zero, zero, twist_wheel_pulley.t.z - top_bearing.width_l / 2)
	top_bearing.place(translate = p)

	# Place the horizontal bevel gear.  (The vertical bevel gear is
	# is already placed relative to the origin.):
	south_bevel_gear.place(axis = x_axis, rotate = degrees_90)

	# Place the shims on next to the bevel gears:
	top_shim.place(translate = P(zero, zero, y1))
	south_shim.place(axis = x_axis,
	  rotate = degrees_90, translate = P(zero, -y1, zero))

	# Drop in the turn table:
	turn_table.place(translate = P(zero, zero, gear_box.t.z))

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
	self.sides_i = sides = 60

	zero = L()
	self.cylinder(comment = "Bearing Cylinder",
	  material = Material("steel", ""),
	  color = Color("blue"),
	  diameter = diameter, 
	  start = P(zero, zero, -width / 2),
	  end =   P(zero, zero,  width / 2),
	  sides = sides)
	self.hole(comment = "Bore Hole",
	  diameter = bore,
	  start = P(zero, zero, -width / 2),
	  end =   P(zero, zero,  width / 2),
	  flags = "t",
	  sides = sides)

class Belt(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	zero = L()
	self.x1_l = zero
	self.y1_l = zero
	self.x2_l = zero
	self.y2_l = zero

    def configure(self, r1 = None, x1 = None, y1 = None, z1 = None,
      r2 = None, x2 = None, y2 = None, z2 = None, width = None):
	# Check argument_types:
	none_type = type(None)
	assert type(x1) == none_type or isinstance(x1, L)
	assert type(y1) == none_type or isinstance(y1, L)
	assert type(z1) == none_type or isinstance(z1, L)
	assert type(r1) == none_type or isinstance(r1, L)
	assert type(x2) == none_type or isinstance(x2, L)
	assert type(y2) == none_type or isinstance(y2, L)
	assert type(z2) == none_type or isinstance(z2, L)
	assert type(r2) == none_type or isinstance(r2, L)
	assert type(width) == none_type or isinstance(width, L)
	
	# Load up configuration values:
	if isinstance(r1, L):
	    self.r1_l = r1
	if isinstance(x1, L):
	    self.x1_l = x1
	if isinstance(y1, L):
	    self.y1_l = y1
	if isinstance(z1, L):
	    self.z1_l = z1
	if isinstance(r2, L):
	    self.r2_l = r2
	if isinstance(x2, L):
	    self.x2_l = x2
	if isinstance(y2, L):
	    self.y2_l = y2
	if isinstance(z2, L):
	    self.z2_l = z2
	if isinstance(width, L):
	    self.width_l = width

    def construct(self):
	# Grab values from *self*:
	r1 = self.r1_l
	x1 = self.x1_l
	y1 = self.y1_l
	z1 = self.z1_l
	r2 = self.r2_l
	x2 = self.x2_l
	y2 = self.y2_l
	z2 = self.z2_l
	width = self.width_l
	zero = L()

	# Figure out direction vector from (*x1*, *y1*) to (*x2*, *y2*).
	dx = x2 - x1
	dy = y2 - y1
	bearing = dy.arc_tangent2(dx)

	# Compute the 4 corners of the belt (*x1a*, *y1a*), (*x1b*, *y1b*),
	# (*x2a*, *y2a*), and (*x2b*, *y2b*):
	sqrt2_r1 = math.sqrt(2) * r1
	angle1a = bearing + Angle(deg = 135)
	x1a = x1 + sqrt2_r1.cosine(angle1a)
	y1a = y1 + sqrt2_r1.sine(angle1a)
	angle1b = bearing - Angle(deg = 135)
	x1b = x1 + sqrt2_r1.cosine(angle1b)
	y1b = y1 + sqrt2_r1.sine(angle1b)
	sqrt2_r2 = math.sqrt(2) * r2
	angle2a = bearing + Angle(deg = 45)
	x2a = x2 + sqrt2_r2.cosine(angle2a)
	y2a = y2 + sqrt2_r2.sine(angle2a)
	angle2b = bearing - Angle(deg = 45)
	x2b = x2 + sqrt2_r2.cosine(angle2b)
	y2b = y2 + sqrt2_r2.sine(angle2b)

	outer_contour = Contour()
	outer_contour.bend_append(P(x2a, y2a, z1), r2)
	outer_contour.bend_append(P(x2b, y2b, z1), r2)
	outer_contour.bend_append(P(x1b, y1b, z1), r1)
	outer_contour.bend_append(P(x1a, y1a, z1), r1)

	# Compute center of belt:
	center_z1 = (P(x1, y1, z1) + P(x2, y2, z1)) / 2
	center_z2 = (P(x1, y1, z2) + P(x2, y2, z2)) / 2

	inner_contour = outer_contour.adjust(
	  delta = -width,
	  start = center_z1,
	  end = center_z2)

	self.extrude(
	  comment = "Belt ({0},{1},{2}) - ({3},{4},{5})".
	   format(r1, x1, y1, r2, x2, y2),
	  material = Material("plastic", "ABS"),
	  color = Color("brown"),
	  outer_contour = outer_contour,
	  inner_contours = [inner_contour],
	  start = center_z1,
	  end = center_z2)

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

class Drive_Motor_Pulley_Screws(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.screw0_ = Fastener(self, comment="screw0")
	self.screw1_ = Fastener(self, comment="screw1")
	self.screw2_ = Fastener(self, comment="screw2")
	self.screw3_ = Fastener(self, comment="screw3")

    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	drive_motor_pulley = wheel_assembly.drive_motor_pulley_
	drive_motor_pulley_top = wheel_assembly.drive_motor_pulley_top_

	# Create the 6 screws:
	screw0 = self.screw0_
	screw1 = self.screw1_
	screw2 = self.screw2_
	screw3 = self.screw3_
	screws = [screw0, screw1, screw2, screw3]

	# Compute some Z axis locations:
	z0 = drive_motor_pulley.b.z
	z1 = drive_motor_pulley_top.b.z + drive_motor_pulley_top.lip_dz_l

	pulley_x = drive_motor_pulley.x_l
	pulley_y = drive_motor_pulley.y_l

	holes_count = len(screws)
	angle_delta = Angle(deg = 360) / holes_count
	holes_radius = drive_motor_pulley.holes_radius_l
	#print("Drive_Motor_Pulley_Screws:holes_radius={0}".
	#  format(holes_radius))
	angle_adjust = Angle(deg = 45)
	for index in range(holes_count):
	    screw = screws[index]
	    angle = index * angle_delta + angle_adjust
	    x = pulley_x + holes_radius.cosine(angle)
	    y = pulley_y + holes_radius.sine(angle)
	    #print("holes_radius={0} x={1} y={2}".format(holes_radius, x, y))
	    screw.configure(comment = "Pulley Hole {0}".format(index),
	      flags = "#0-80:hi:fh",
	      start = P(x, y, z0),
	      end = P(x, y, z1),
	      sides_angle = Angle(deg = 30))
	    screw.drill(drive_motor_pulley, select = "close")
	    screw.drill(drive_motor_pulley_top, select = "close")

	#print("Pulley_Screw.dz={0:.3i}in".format(z3 - z0))

class Drive_Wheel_Pulley_Screws(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.screw0_ = Fastener(self, comment="screw0")
	self.screw1_ = Fastener(self, comment="screw1")
	self.screw2_ = Fastener(self, comment="screw2")
	self.screw3_ = Fastener(self, comment="screw3")

    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	drive_wheel_pulley = wheel_assembly.drive_wheel_pulley_
	drive_wheel_pulley_top = wheel_assembly.drive_wheel_pulley_top_

	# Create the 6 screws:
	screw0 = self.screw0_
	screw1 = self.screw1_
	screw2 = self.screw2_
	screw3 = self.screw3_
	screws = [screw0, screw1, screw2, screw3]

	# Compute some Z axis locations:
	z0 = drive_wheel_pulley.b.z
	z1 = drive_wheel_pulley_top.t.z

	holes_count = len(screws)
	angle_delta = Angle(deg = 360) / holes_count
	holes_radius = drive_wheel_pulley.holes_radius_l
	angle_adjust = Angle(deg = 45)
	for index in range(holes_count):
	    screw = screws[index]
	    angle = index * angle_delta + angle_adjust
	    x = holes_radius.cosine(angle)
	    y = holes_radius.sine(angle)
	    #print("holes_radius={0} x={1} y={2}".format(holes_radius, x, y))
	    screw.configure(comment = "Pulley Hole {0}".format(index),
	      flags = "#0-80:hi:fh",
	      start = P(x, y, z0),
	      end = P(x, y, z1),
	      sides_angle = Angle(deg = 30))
	    screw.drill(drive_wheel_pulley, select = "close")
	    screw.drill(drive_wheel_pulley_top, select = "close")

	#print("Pulley_Screw.dz={0:.3i}in".format(z3 - z0))

class End_Screws(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.screw0_ = Fastener(self, comment="screw0")
	self.screw1_ = Fastener(self, comment="screw1")
	self.screw2_ = Fastener(self, comment="screw2")
	self.screw3_ = Fastener(self, comment="screw3")

    def construct(self):
	# Grab some values from *motor_assembly*:
	motor_assembly = self.up
	motor_base = motor_assembly.motor_base_

	# Grab some value from *motor_base*:
	c_x = motor_base.c_x_l
	c_y = motor_base.c_y_l
	d_x = motor_base.d_x_l
	d_y = motor_base.d_y_l
	floor_dz = motor_base.floor_dz_l
	wall_width = motor_base.wall_width_l

	# Create the 6 screws:
	screw0 = self.screw0_
	screw1 = self.screw1_
	screw2 = self.screw2_
	screw3 = self.screw3_

	# Compute some Z axis locations:
	x0 = d_x + wall_width
	x3 = c_x - wall_width
	dx = x3 - x0
	x1 = x0 + .20 * dx
	x2 = x0 + .80 * dx

	y0 = motor_base.d_y_l
	y1 = y0 + wall_width

	z0 = motor_base.b.z
	z1 = z0 + floor_dz
	z4 = motor_base.t.z
	dz = z4 - z1
	z2 = z1 + .25 * dz
	z3 = z1 + .75 * dz

	screw_records = [ \
	  (screw0, x1, z2), \
	  (screw1, x2, z2), \
	  (screw2, x1, z3), \
	  (screw3, x2, z3) ]

	for index in range(len(screw_records)):
	    screw_record = screw_records[index]
	    screw = screw_record[0]
	    x = screw_record[1]
	    z = screw_record[2]
	    
	    #print("End_Screws:start={0} end={1}".format(start, end))
	    screw.configure(comment = "End Screw {0}".format(index),
	      flags = "#6-32:hi",
	      start = P(x, y0, z),
	      end = P(x, y1, z),
	      sides_angle = Angle(deg = 30))
	    screw.drill(motor_base, select = "close")

	#print("Pulley_Screw.dz={0:.3i}in".format(z3 - z0))

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
	    gear_box.n.y - gear_box_side.dy, z3),
	  top = "t")

	# Pocket out some the mounting area:
	self.simple_pocket(comment = "Gear Box Bottom Pocket",
	  bottom_corner = P(-pocket_dx / 2, -pocket_dy / 2, z1),
	  top_corner =    P( pocket_dx / 2,  pocket_dy / 2, z3))

class Gear_Box_Cover(Part):
    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_
	gear_box_top = wheel_assembly.gear_box_top_

	# Construct the cover out of a block of material:
	dx = L(mm = 2.50)
	self.block(comment = "Bevel Gear Box Cover Block",
	  material = Material("plastic", "ABS"),
	  color = Color("chartreuse"),
	  corner1 = P(gear_box.w.x,      gear_box.s.y, gear_box.b.z),
	  corner2 = P(gear_box.w.x + dx,
	    gear_box.n.y, gear_box.t.z - gear_box_top.base_dz_l),
	  top = "w")

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
	self.en_screw_ = Fastener(self, comment="en_screw")
	self.es_screw_ = Fastener(self, comment="es_screw")
	self.wn_screw_ = Fastener(self, comment="wn_screw")
	self.ws_screw_ = Fastener(self, comment="ws_screw")
	self.ne_screw_ = Fastener(self, comment="ne_screw")
	self.nw_screw_ = Fastener(self, comment="nw_screw")
	self.se_screw_ = Fastener(self, comment="se_screw")
	self.sw_screw_ = Fastener(self, comment="sw_screw")
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
	sides_angle = Angle(deg = 30)
	en_screw.configure(
	  comment = "East Side North Screw",
	  flags = "#4-40:hi",
	  start = P(x7, y4, z),
	  end = P(x5, y4, z),
	  sides_angle = sides_angle)
	es_screw.configure(
	  comment = "East Side South Screw",
	  flags = "#4-40:hi",
	  start = P(x7, y3, z),
	  end = P(x5, y3, z),
	  sides_angle = sides_angle)
	wn_screw.configure(
	  comment = "West Side North Screw",
	  flags = "#4-40:hi",
	  start = P(x0, y4, z),
	  end = P(x2, y4, z),
	  sides_angle = sides_angle)
	ws_screw.configure(
	  comment = "West Side South Screw",
	  flags = "#4-40:hi",
	  start = P(x0, y3, z),
	  end = P(x2, y3, z),
	  sides_angle = sides_angle)
	ne_screw.configure(
	  comment = "North Side East Screw",
	  flags = "#4-40:hi",
	  start = P(x4, y7, z),
	  end = P(x4, y5, z))
	nw_screw.configure(
	  comment = "North Side West Screw",
	  flags = "#4-40:hi",
	  start = P(x3, y7, z),
	  end = P(x3, y5, z))
	se_screw.configure(
	  comment = "South Side East Screw",
	  flags = "#4-40:hi",
	  start = P(x4, y0, z),
	  end = P(x4, y2, z))
	sw_screw.configure(
	  comment = "South Side West Screw",
	  flags = "#4-40:hi",
	  start = P(x3, y0, z),
	  end = P(x3, y2, z))

	# Drill 8 holes for in *gear_box_part*:
	gear_box_part = self._gear_box_part
	en_screw.drill(gear_box_part, select = "close")
	es_screw.drill(gear_box_part, select = "close")
	wn_screw.drill(gear_box_part, select = "close")
	ws_screw.drill(gear_box_part, select = "close")
	ne_screw.drill(gear_box_part, select = "close")
	nw_screw.drill(gear_box_part, select = "close")
	se_screw.drill(gear_box_part, select = "close")
	sw_screw.drill(gear_box_part, select = "close")

	# Drill the sides:
	ne_screw.drill(north_gear_box_side, select = "close")
	nw_screw.drill(north_gear_box_side, select = "close")
	se_screw.drill(south_gear_box_side, select = "close")
	sw_screw.drill(south_gear_box_side, select = "close")

	# Drill the covers:
	#en_screw.drill(east_gear_box_cover)
	#es_screw.drill(east_gear_box_cover)
	wn_screw.drill(west_gear_box_cover, select = "close")
	ws_screw.drill(west_gear_box_cover, select = "close")

	#print("<=Gear_Box_Screws.construct():label='{0}'".format(self._label))

class Gear_Box_Shelf(Part):
    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	bearing = wheel_assembly.south_bearing_
	gear_box = wheel_assembly.gear_box_
	gear_box_cover = wheel_assembly.west_gear_box_cover_
	gear_box_side = wheel_assembly.south_gear_box_side_
	twist_wheel_pulley = wheel_assembly.twist_wheel_pulley_
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
	z3 = z2 + top_shim.dz
	z4 = z3 + bearing.width_l
	z5 = z4 + bearing_lip_dz
	self.screw_z_l = z1

	dx = gear_box.dx - 2 * gear_box_cover.dx
	dy = gear_box.dy - 2 * gear_box_side.dy

	# Generate the shelf block:
	self.block(comment = "Shelf Block",
	  material = Material("plastic", "ABS"),
	  color = Color("gold"),
	  corner1 = P(-dx / 2, -dy / 2, z0),	# should z0
	  corner2 = P( dx / 2,  dy / 2, z5),
	  top = "b")

	# Cut a hole for the shaft:
	zero = L()
	self.hole(comment = "Shaft Hole",
	  diameter = twist_wheel_pulley.shaft_diameter_l,
	  start = P(zero, zero, z0),
	  end = P(zero, zero, z5),
	  flags = "t")

	# Cut a bearing hole:
	self.hole(comment = "Bearing Hole",
	  diameter = bearing.diameter_l,
	  start = P(zero, zero, z0),
	  end = P(zero, zero, z4),
	  flags = "f",
	  sides = 60)

	# Cut out the pocket:
	self.simple_pocket(comment = "Mount Pocket",
	  bottom_corner = P(-pocket_dx / 2, -pocket_dy / 2, z0),
	  top_corner = P( pocket_dx / 2,  pocket_dy / 2, z2))

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
	top = "s"
	if south_side:
	    y1 = -y1
	    y2 = -y2
	    y3 = -y3
	    top = "n"

	# Do the main block:
	self.block(comment = self._name + " Main Block",
	  material = Material("plastic", "abs"),
	  color = Color("orange"),
	  corner1 =
	    P(gear_box.w.x + gear_box_cover.dx, y1, gear_box.b.z),
	  corner2 = P(gear_box.e.x - gear_box_cover.dx, y3,
	    gear_box.t.z - gear_box_top.base_dz_l),
	  top = top)

	# Do the bearing hole:
	zero = L()
	self.hole(comment = "Bearing hole",
	  diameter = bearing.diameter_l,
	  start = P(zero, y3, zero),
	  end = P(zero, y2, zero),
	  flags = "f",
	  sides = 60)

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
	twist_wheel_pulley = wheel_assembly.twist_wheel_pulley_

	# Grab the pocket size from *wheel_assembly*:
	pocket_dx = wheel_assembly.pocket_dx_l
	pocket_dy = wheel_assembly.pocket_dy_l
	pocket_dz = wheel_assembly.pocket_dz_l
	self.base_dz_l = base_dz = L(mm = 5.00)

	# Compute various locations along the Z Axis:
	z0 = gear_box.t.z - base_dz - pocket_dz
	z1 = z0 + pocket_dz / 2
	z2 = z1 + pocket_dz / 2
	z3 = z2 + base_dz / 2
	z4 = z3 + base_dz / 2
	self.screw_z_l = z1
	#print("Gear_Box_Top.construct: z0={0} z1={1} z2={2} z3={3} z4={4}".
	#  format(z0, z1, z2, z3, z4))

	# Get the basic base in place:
	self.block(comment = "Bevel Gear Box Top Block",
	  material = Material("plastic", "ABS"),
	  color = Color("azure"),
	  corner1 = P(-turn_table.size_l / 2, -turn_table.size_l / 2, z2),
	  corner2 = P( turn_table.size_l / 2,  turn_table.size_l / 2, z4),
	  top = "b")

	# Add some stuff to mount to (merge it into the block):
	self.block(comment = "Bevel Gear Box Top Mount Block",
	  corner1 = P(gear_box.w.x + gear_box_cover.dx,
	    gear_box.s.y + gear_box_side.dy, z0),
	  corner2 = P(gear_box.e.x - gear_box_cover.dx,
	    gear_box.n.y - gear_box_side.dy, z3),
	  top = "b")

	# Do the shaft hole:
	zero = L()
	self.hole(comment = "Shaft Hole",
	  diameter = twist_wheel_pulley.shaft_diameter_l,
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
	  bottom_corner = P( pocket_dx / 2,  pocket_dy / 2, z2),
	  pocket_top = "b")

# The D42DAI from http://www.kjmagnetics.com/ is a 1/4" by 1/8"
# thick nickel plated diametrically magnetized magnet for $.38/each.
# The D82DIA is 1/2" x 1/8" for $1.13/each.

class Magnet(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	zero = L()
	Part.diameter_l = L(inch = "1/4")
	Part.x_l = zero
	Part.y_l = zero
	Part.z_bottom_l = zero
	Part.dz_l = L(inch = "1/8")

    def configure(self, diameter = None,
      x = None, y = None, z_bottom = None, dz = None):
	# Perform argument checks:
	none_type = type(None)
	assert type(diameter) == none_type or isinstance(diameter, L)
	assert type(x) == none_type or isinstance(x, L)
	assert type(y) == none_type or isinstance(y, L)
	assert type(z_bottom) == none_type or isinstance(z_bottom, L)
	assert type(dz) == none_type or isinstance(dz, L)

	# Load up *self*:
	if isinstance(diameter, L):
	    self.diameter_l = diameter
	if isinstance(x, L):
	    self.x_l = x
	if isinstance(y, L):
	    self.y_l = y
	if isinstance(z_bottom, L):
	    self.z_bottom_l = z_bottom
	if isinstance(dz, L):
	    self.dz_l = dz

    def construct(self):
	# Grab some values from *self*:
	x = self.x_l
	y = self.y_l
	z_bottom = self.z_bottom_l 
	dz = self.dz_l

	# Represent the magnet using a cylinder:
	self.cylinder(comment = "Magnet",
	  material = Material("iron", "magnetic"),
	  color = Color("cyan"),
	  diameter = self.diameter_l,
	  start = P(x, y, z_bottom),
	  end = P(x, y, z_bottom + dz))

class Motor(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	zero = L()
	self.shaft_diameter_l = zero
	self.spur_gear_body_dz_l = zero
	self.x_l = zero
	self.y_l = zero
	self.z_l = zero

    def configure(self, x = None, y = None, z = None, spur_gear_body_dz = None):
	# Check argument types:
	none_type = type(None)
	assert type(x) == none_type or isinstance(x, L)
	assert type(y) == none_type or isinstance(x, L)
	assert type(z) == none_type or isinstance(x, L)
	assert type(spur_gear_body_dz) == none_type or \
	  isinstance(spur_gear_body_dz, L)

	# Load up *self*:
	if isinstance(x, L):
	    self.x_l = x
	if isinstance(y, L):
	    self.y_l = y
	if isinstance(z, L):
	    self.z_l = z
	if isinstance(spur_gear_body_dz, L):
	    self.spur_gear_body_dz_l = spur_gear_body_dz

    def construct(self):
	x = self.x_l
	y = self.y_l
	z = self.z_l
	spur_gear_body_dz = self.spur_gear_body_dz_l

	zero = L()

	self.spur_gear_hub_dz_l = spur_gear_hub_dz = L(mm = 4.60)

	weld = L(mm = 1.00)
	z0 = z - L(mm = 26.00) - spur_gear_body_dz	# Bottom
	z1 = z0 + L(mm = 3.00)				# Motor shaft hub
	z2 = z1 + L(mm = 3.00)				# Motor back
	z3 = z2 + L(mm = 20.00)				# Motor body
	z4 = z3 + spur_gear_body_dz			# Spur gear body
	z5 = z4 + spur_gear_hub_dz			# Spur gear hub body
	z6 = z5 + L(mm = 4.00)				# Round shaft
	z7 = z6 + L(mm = 12.00)				# Shaft with flat

	self.shaft_diameter_l = d1 = L(mm = 6.00)	# Motor shaft diameter
	d2 = L(mm = 9.00)		# Motor hub diameter
	self.spur_gear_hub_diameter_l = d3 = L(mm = 12.00) # Spur gear hub diam.
	d4 = L(mm = 33.00)		# Motor body diameter
	self.spur_gear_body_diameter_l = d5 = L(mm = 37.00)# Spur gear body dia.

	r0 = zero
	self.shaft_offset_l = r1 = L(mm = 7.00)		# Shaft offset
	self.mount_hole_radius_l = r2 = L(mm = 31.00) / 2 # Mount hole radius

	self.mount_angle0_a = a0 = Angle(0)
	self.mount_angle1_a = a1 = a0 + Angle(deg = 120)
	self.mount_angle2_a = a2 = a1 + Angle(deg = 60)
	self.mount_angle3_a = a3 = a2 + Angle(deg = 60)

	self.cylinder(comment = "Motor shaft hub",
	  material = Material("metal", "steel"),
	  color = Color("magenta"),
	  diameter = d1,
	  start = P(x, y, z0),
	  end = P(x, y, z1 + weld))
	self.cylinder(comment = "Motor body",
	  diameter = d4,
	  start = P(x, y, z1),
	  end = P(x, y, z3 + weld),
	  sides = 30)
	self.cylinder(comment = "Spur gear body",
	  diameter = d5,
	  start = P(x, y, z3),
	  end = P(x, y, z4),
	  sides = 30)
	self.cylinder(comment = "Spur gear hub",
	  diameter = d3,
	  start = P(x, y - r1, z4 - weld),
	  end = P(x, y - r1, z5))
	self.cylinder(comment = "Shaft",
	  diameter = r1,
	  start = P(x, y - r1, z5 - weld),
	  end = P(x, y - r1, z7))

class Motor_Base(Part):
    def construct(self):

	#print("=>Motor_Base.construct()")
	motor_assembly = self.up
	synchro_drive = motor_assembly.up
	drive_motor = motor_assembly.drive_motor_
	pcb = motor_assembly.pcb_
	twist_motor = motor_assembly.twist_motor_
	wheel_assembly = synchro_drive.wheel_assembly_
	turn_table = wheel_assembly.turn_table_
	vertical_shaft = wheel_assembly.vertical_shaft_
	wheel = wheel_assembly.wheel_

	turn_table_size = turn_table.size_l
	turn_table_bore = turn_table.bore_l
	turn_table_bottom_hole_pitch1 = turn_table.bottom_hole_pitch1_l

	x_adjust = L(mm=-3)
	y_adjust = L(mm=-4)

	self.floor_dz_l = floor_dz = twist_motor.spur_gear_hub_dz_l
	#self.dz_l = dz = \
	#  vertical_shaft.t.z + L(mm = 2) - motor_assembly.base_bz_l
	self.dz_l = dz = L(mm = 33.00)
	self.ny_l = ny = L(mm = 123.00) + y_adjust
	self.ex_l = ex = L(mm = 105.00) + x_adjust
	self.wall_width_l = wall_width = L(mm = 6.00)

	self.z0_l = z0 = motor_assembly.base_bz_l
	self.z1_l = z1 = z0 + L(mm = 0.20)
	self.z2_l = z2 = z1 + L(mm = 0.20)
	self.z3_l = z3 = z0 + floor_dz / 2
	self.z4_l = z4 = z3 + floor_dz / 2
	self.z5_l = z5 = z0 + dz - pcb.dz
	self.z6_l = z6 = z0 + dz
	#print("Motor_Base: z0={0} z3={1} z4={2} z6={3}".format(z0, z3, z4, z6))

	# Below is the crude (not to scale) ASCII art for the locations
	# of the base contour (letters A-M are the corners).  O indicates
	# a screw hole and S indicates the shaft hole.
	#
	#                U---V                     Y---Z
	#                |   |                     |   |
	#                |   |                     |   |
	#                |   |                     |   |
	#                |   |                     |   |
	#                |   |                     |   |
	#                |   |                     |   |
	#                |   |                     |   |
	#                K...W---------------------X...L
	#                | O                         O |
	#                |                             |
	#                |                             |
	#                |                             |
	#                |                             |
	#                J O                         O M
	#               /                             /
	#              /                             /
	#             ~	                            ~
	#            /                             /
	#   HH------I O ----------------+         /
	#   |       |                   .        /
	#   |   O   |               O   .       /
	#   |       |                   .      /
	#   G-------H                   .     /
	#   | O                         .    /
	#   |                           .   /
	#   |             S             .  /
	#   |                           . /
	#   | O                        O./
	#   F-------E           B-------A
	#   |       |           |       |
	#   |   O   |           |   O   |
	#   |       | O       O |       |
	#   EE------D-----------C------BB

	z = zero = L()
	shaft_offset = twist_motor.shaft_offset_l
	ww = wall_width
	ww2 = ww / 2
	ttx2 = turn_table.dx / 2
	tty2 = turn_table.dy / 2
	clearance = cl =L(mm = 20)
	jx = kx = L(mm = 5) + x_adjust
	#jy = motor_assembly.twist_motor_y_l - shaft_offset
	jy = L(mm = 75) + y_adjust
	#my = motor_assembly.drive_motor_y_l - shaft_offset
	my = L(mm = 70) + y_adjust
	#print("my={0} pcb.s.y={1} pcb.n.y={2}".format(my, pcb.s.y, pcb.n.y))

	# Define locations of each contour corner and screw (if needed):
	self.a_x_l = a_x = ttx2 + ww
	self.a_y_l = a_y = -tty2 + cl
	self.screw_a_x_l = a_x - ww2
	self.screw_a_y_l = a_y + ww2

	self.b_x_l = b_x = ttx2 - cl
	self.b_y_l = b_y = -tty2 + cl

	self.c_x_l = c_x = ttx2 - cl
	self.c_y_l = c_y = -tty2 - ww
	self.screw_c_x_l = c_x - ww2
	self.screw_c_y_l = c_y + ww2

	self.d_x_l = d_x = -ttx2 + cl
	self.d_y_l = d_y = -tty2 - ww
	self.screw_d_x_l = d_x + ww2
	self.screw_d_y_l = d_y + ww2

	self.e_x_l = e_x = -ttx2 + cl
	self.e_y_l = e_y = -tty2 + cl

	self.f_x_l = f_x = -ttx2 - ww
	self.f_y_l = f_y = -tty2 + cl
	self.screw_f_x_l = f_x + ww2
	self.screw_f_y_l = f_y + ww2

	self.g_x_l = g_x = -ttx2 - ww
	self.g_y_l = g_y = tty2 - cl
	self.screw_g_x_l = g_x + ww2
	self.screw_g_y_l = g_y - ww2

	self.h_x_l = h_x = -ttx2 + cl
	self.h_y_l = h_y = tty2 - cl

	self.i_x_l = i_x = -ttx2 + cl
	self.i_y_l = i_y = tty2 + ww
	self.screw_i_x_l = i_x + ww2
	self.screw_i_y_l = i_y

	self.j_x_l = j_x = jx
	self.j_y_l = j_y = jy
	self.screw_j_x_l = j_x + ww2
	self.screw_j_y_l = j_y

	self.k_x_l = k_x = kx
	self.k_y_l = k_y = ny
	self.screw_k_x_l = k_x + ww2
	self.screw_k_y_l = k_y - ww2

	self.l_x_l = l_x = ex
	self.l_y_l = l_y = ny
	self.screw_l_x_l = l_x - ww2
	self.screw_l_y_l = l_y - ww2

	yyy = L(mm = 50.00)
	strut_width = L(inch = 0.500)

	self.m_x_l = m_x = ex
	self.m_y_l = m_y = my
	self.screw_m_x_l = m_x - ww2
	self.screw_m_y_l = m_y + ww2

	pcb_n_y = m_y + L(mm = 100)


	self.u_x_l = u_x = k_x
	self.u_y_l = u_y = pcb_n_y

	self.v_x_l = v_x = u_x + strut_width
	self.v_x_l = v_y = pcb_n_y

	self.w_x_l = w_x = v_x
	self.w_y_l = w_y = k_y

	self.x_x_l = x_x = l_x - strut_width
	self.x_y_l = x_y = w_y

	self.y_x_l = y_x = x_x
	self.y_y_l = y_y = pcb_n_y

	self.z_x_l = z_x = l_x
	self.z_y_l = z_y = pcb_n_y

	self.hh_x = hh_x = g_x
	self.hh_y = hh_y = i_y

	self.ee_x = ee_x = f_x
	self.ee_y = ee_y = d_y

	self.bb_x = bb_x = a_x
	self.bb_y = bb_y = c_y

	self.attach_radius_l = attach_radius = L(inch = 0.24)
	self.turn_radius_l = turn_radius = (a_x - b_x) / 2 - L(inch = .001)

	# Make sure the bottom is made out of a single *bottom_contour*
	# extrusion:
	r = L(mm = 10)
	lip_y = L(mm = 4.5) + L(mm=4)
	helper_contour = Contour()
	helper_contour.bend_append(P(m_x + L(mm=3), m_y, z), r)		# M
	helper_contour.bend_append(P(z_x + L(mm=3), z_y + lip_y, z), r)	# Z
	helper_contour.bend_append(P(y_x - L(mm=30), y_y + lip_y, z), r) # Y
	helper_contour.bend_append(P(x_x, x_y, z), z)			# X
	helper_contour.bend_append(P(w_x, w_y, z), z)			# W
	helper_contour.bend_append(P(v_x + L(mm=30), v_y + lip_y, z), r) # V
	helper_contour.bend_append(P(u_x - L(mm=30), u_y + lip_y, z), r) # U
	helper_contour.bend_append(P(k_x, k_y, z), z)			# K
	helper_contour.bend_append(P(j_x, j_y, z), z)			# J
	#helper_contour.bend_append(P(i_x, i_y, z), z)			# I
	helper_contour.bend_append(P(hh_x, hh_y + L(mm=20), z), r)	# HH
	helper_contour.bend_append(P(ee_x, ee_y - lip_y, z), z)		# EE
	helper_contour.bend_append(P(bb_x + L(mm=20), bb_y - lip_y, z), r) # BB
	#helper_contour.bend_append(P(a_x, a_y, z), z)			# A

	self.extrude(comment = "Motor Bottom Helper Extrusion",
	  material = Material("plastic", "ABS"),
	  color = Color("lavender"),
	  outer_contour = helper_contour,
	  start = P(zero, zero, z0),
	  end = P(zero, zero, z2))

	# Make sure the bottom is made out of a single *bottom_contour*
	# extrusion:
	bottom_contour = Contour()
	bottom_contour.bend_append(P(m_x, m_y, z), z)			# M
	bottom_contour.bend_append(P(z_x, z_y, z), attach_radius)	# Z
	bottom_contour.bend_append(P(y_x, y_y, z), attach_radius)	# Y
	bottom_contour.bend_append(P(x_x, x_y, z), z)			# X
	bottom_contour.bend_append(P(w_x, w_y, z), z)			# W
	bottom_contour.bend_append(P(v_x, v_y, z), attach_radius)	# V
	bottom_contour.bend_append(P(u_x, u_y, z), attach_radius)	# U
	bottom_contour.bend_append(P(j_x, j_y, z), z)			# J
	bottom_contour.bend_append(P(i_x, i_y, z), z)			# I
	bottom_contour.bend_append(P(hh_x, hh_y, z), turn_radius)	# HH
	bottom_contour.bend_append(P(ee_x, ee_y, z), turn_radius)	# EE
	bottom_contour.bend_append(P(bb_x, bb_y, z), turn_radius)	# BB
	bottom_contour.bend_append(P(a_x, a_y, z), z)			# A

	self.extrude(comment = "Motor Bottom Contour Extrusion",
	  material = Material("plastic", "ABS"),
	  color = Color("lavender"),
	  outer_contour = bottom_contour,
	  start = P(zero, zero, z1),
	  end = P(zero, zero, z4))

	# Now construct the wall out of a wall extrusion:
	self.outer_countour_o = outer_contour = Contour()
	outer_contour.bend_append(P(m_x, m_y, z), z)		# M
	outer_contour.bend_append(P(l_x, l_y, z), z)		# L
	outer_contour.bend_append(P(k_x, k_y, z), z)		# K
	outer_contour.bend_append(P(j_x, j_y, z), z)		# J
	outer_contour.bend_append(P(i_x, i_y, z), z)		# I
	outer_contour.bend_append(P(h_x, h_y, z), turn_radius)	# H
	outer_contour.bend_append(P(g_x, g_y, z), z)		# G
	outer_contour.bend_append(P(f_x, f_y, z), z)		# F
	outer_contour.bend_append(P(e_x, e_y, z), turn_radius)	# E
	outer_contour.bend_append(P(d_x, d_y, z), z)		# D
	outer_contour.bend_append(P(c_x, c_y, z), z)		# C
	outer_contour.bend_append(P(b_x, b_y, z), turn_radius)	# B
	outer_contour.bend_append(P(a_x, a_y, z), z)		# A

	# *start* and *end* specify the extrude direction axis:
	start = P(zero, zero, z0)
	end = P(zero, zero, z3)

	# Make an *inner_contour* that is *wall_width* inside *outer_contour*:
	inner_contour = outer_contour.adjust(
	  delta = -wall_width,
	  start = start,
	  end = end)

	#print("outer_contour={0}".format(outer_contour))
	#print("inner_contour={0}".format(inner_contour))

	# Now put a wall around and force it to "weld" into the base:
	self.extrude(comment = "Motor Base Wall",
	  outer_contour = outer_contour,
	  inner_contours = [inner_contour],
	  start = P(zero, zero, z3),
	  end = P(zero, zero, z6))

	self.attach_dx_l = attach_dx = L(inch = "1/2")
	pcb_n_y = m_y + L(mm = 100)

	# Now put on some attach panels:
	# FIXME: *pcb* bounding box is broken:
	r = attach_radius
	east_attach_contour = Contour()
	east_attach_contour.bend_append(P(k_x,             k_y - ww2, z), z)
	east_attach_contour.bend_append(P(k_x + attach_dx, k_y - ww2, z), z)
	east_attach_contour.bend_append(P(k_x + attach_dx, pcb_n_y,   z), r)
	east_attach_contour.bend_append(P(k_x,             pcb_n_y,   z), r)
	self.extrude(comment = "East Attach",
	  outer_contour = east_attach_contour,
	  start = P(z, z, z3),
	  end =   P(z, z, z6))

	west_attach_contour = Contour()
	west_attach_contour.bend_append(P(l_x,             l_y - ww2, z), z)
	west_attach_contour.bend_append(P(l_x - attach_dx, l_y - ww2, z), z)
	west_attach_contour.bend_append(P(l_x - attach_dx, pcb_n_y,   z), r)
	west_attach_contour.bend_append(P(l_x,             pcb_n_y,   z), r)
	self.extrude(comment = "West Attach",
	  outer_contour = west_attach_contour,
	  start = P(z, z, z3),
	  end =   P(z, z, z6))

	#self.block(comment = "West Attach",
	#  corner1 = P(l_x - attach_dx, l_y - ww2, z3),
	#  corner2 = P(l_x,             pcb_n_y,   z5),
	#  top = "w")

	# Define the screw locatons for the attach panels:
	self.screw_n_x_l = k_x + attach_dx / 2
	self.screw_n_y_l = k_y + attach_dx / 2
	self.screw_o_x_l = k_x + attach_dx / 2
	self.screw_o_y_l = pcb_n_y - attach_dx / 2
	self.screw_p_x_l = l_x - attach_dx / 2
	self.screw_p_y_l = k_y + attach_dx / 2
	self.screw_q_x_l = l_x - attach_dx / 2
	self.screw_q_y_l = pcb_n_y - attach_dx / 2

	# Make a hole for the shaft/wheel pulley/drive pulley to poke through:
	# Leave some extra *clearance* so that the belts can be tightened.
	self.hole(comment = "Twist/Drive Shaft Hole",
	  diameter = turn_table_bore + L(mm = 16.00),
	  start = P(zero, zero, z4),
	  end =   P(zero, zero, z0),
	  flags = "t",
	  sides = 60)

	# Now put in the holes for the twist motor:
	hub_diameter = twist_motor.spur_gear_hub_diameter_l
	twist_hub_diameter = hub_diameter + L(mm = 0.25)
	twist_motor_x = motor_assembly.twist_motor_x_l
	twist_motor_y = motor_assembly.twist_motor_y_l
	twist_shaft_y = twist_motor_y - twist_motor.shaft_offset_l
	self.hole(comment = "Twist Motor Shaft Hole",
	  diameter = twist_hub_diameter,
	  start = P(twist_motor_x, twist_shaft_y, z4),
	  end =   P(twist_motor_x, twist_shaft_y, z0),
	  flags = "t",
	  sides = 30)

	# Now put in the hole for the drive motor:
	drive_hub_diameter = twist_hub_diameter
	drive_motor_x = motor_assembly.drive_motor_x_l
	drive_motor_y = motor_assembly.drive_motor_y_l
	drive_shaft_y = drive_motor_y - drive_motor.shaft_offset_l
	self.hole(comment = "Drive Motor Shaft Hole",
	  diameter = drive_hub_diameter,
	  start = P(drive_motor_x, drive_shaft_y, z4),
	  end =   P(drive_motor_x, drive_shaft_y, z0),
	  flags = "t",
	  sides = 30)

	#print("<=Motor_Base.construct()")

	# Now put on some attach panels:
	# FIXME: *pcb* bounding box is broken:
	# Install the the turn table mouting holes.  These are oversized
	# so that the two timing belts can be tightened:
	hole_pitch = turn_table_bottom_hole_pitch1
	for dx in (-hole_pitch/2, hole_pitch/2):
	    for dy in (-hole_pitch/2, hole_pitch/2):
		self.hole(comment = "Turntable_hole ({0},{1})".format(dx, dy),
		  diameter = L(inch = 0.500),
		  start = P(dx, dy, z4),
		  end = P(dx, dy, z0),
		  flags = "t")

	# Make room for the PCB to fit in:
	# FIXME: Should use the PCB bounding box which is currently broken!!!:
	extra = L(mm = .01)
	self.simple_pocket(comment = "PCB pocket",
	  bottom_corner = P(j_x - extra, m_y, z5),
	  top_corner = P(m_x + extra, z_y + extra, z6))

class Motor_Base_Cover(Part):
    def construct(self):
	#print("=>Motor_Base.construct()")
	motor_assembly = self.up
	motor_base = motor_assembly.motor_base_

	zero = L()
	z0 = motor_base.t.z
	z1 = z0 + L(mm = 3.0)

	outer_contour = motor_base.outer_countour_o
	self.extrude(
	  comment="Motor_Base_Cover",
	  material = Material("plastic", "ABS"),
	  color = Color("yellow"),
	  outer_contour = outer_contour,
	  start = P(zero, zero, z0),
	  end = P(zero, zero, z1))

class Motor_Base_Screws(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.screw_a_ = Fastener(self, comment="screw_a")
	self.screw_c_ = Fastener(self, comment="screw_c")
	self.screw_d_ = Fastener(self, comment="screw_d")
	self.screw_f_ = Fastener(self, comment="screw_f")
	self.screw_g_ = Fastener(self, comment="screw_g")
	self.screw_i_ = Fastener(self, comment="screw_i")
	self.screw_j_ = Fastener(self, comment="screw_j")
	self.screw_k_ = Fastener(self, comment="screw_k")
	self.screw_l_ = Fastener(self, comment="screw_l")
	self.screw_m_ = Fastener(self, comment="screw_m")
	self.screw_n_ = Fastener(self, comment="screw_n")
	self.screw_o_ = Fastener(self, comment="screw_o")
	self.screw_p_ = Fastener(self, comment="screw_p")
	self.screw_q_ = Fastener(self, comment="screw_q")
	self.drive_screw_a_ = Fastener(self, comment="drive_screw_a")
	self.drive_screw_b_ = Fastener(self, comment="drive_screw_b")
	self.twist_screw_a_ = Fastener(self, comment="twist_screw_a")
	self.twist_screw_b_ = Fastener(self, comment="twist_screw_b")

    def construct(self):
	ma = motor_assembly = self.up
	mb = motor_assembly.motor_base_
	pcb = motor_assembly.pcb_
	drive_motor = motor_assembly.drive_motor_
	twist_motor = motor_assembly.twist_motor_
	mbc = motor_assembly.motor_base_cover_

	self._pcb_holes_file = open("pcb_holes.txt", "w")
	self._x1 = x1 = pcb.x1_l
	self._y1 = y1 = pcb.y1_l

	self.x_y_log("pcb corner1", x1, y1)
	self.x_y_log("pcb corner2", pcb.x2_l, pcb.y2_l)
	self.x_y_log("drive_shaft",
	  ma.drive_motor_shaft_x_l, ma.drive_motor_shaft_y_l)
	self.x_y_log("twist_shaft",
	  ma.twist_motor_shaft_x_l, ma.twist_motor_shaft_y_l)

	screw_records = ( \
	  (self.screw_a_, "A", mb.screw_a_x_l, mb.screw_a_y_l), \
	  (self.screw_c_, "C", mb.screw_c_x_l, mb.screw_c_y_l), \
	  (self.screw_d_, "D", mb.screw_d_x_l, mb.screw_d_y_l), \
	  (self.screw_f_, "F", mb.screw_f_x_l, mb.screw_f_y_l), \
	  (self.screw_g_, "G", mb.screw_g_x_l, mb.screw_g_y_l), \
	  (self.screw_i_, "I", mb.screw_i_x_l, mb.screw_i_y_l), \
	  (self.screw_j_, "J", mb.screw_j_x_l, mb.screw_j_y_l), \
	  (self.screw_k_, "K", mb.screw_k_x_l, mb.screw_k_y_l), \
	  (self.screw_l_, "L", mb.screw_l_x_l, mb.screw_l_y_l), \
	  (self.screw_m_, "M", mb.screw_m_x_l, mb.screw_m_y_l)  )
	for screw_record in screw_records:
	    screw = screw_record[0]
	    screw_letter = screw_record[1]
	    screw_x = screw_record[2]
	    screw_y = screw_record[3]
	    screw.configure(
	      comment = "Screw {0}".format(screw_letter),
	      start = P(screw_x, screw_y, mbc.t.z),
	      end = P(screw_x, screw_y, mb.b.z),
	      flags = "#4-40")
	    screw.drill(
	      part = mb,
	      select = "close")
	    screw.drill(
	      part = mbc,
	      select = "close")
	    if "JKLM".find(screw_letter) >= 0:
		screw.drill(
		  part = pcb,
		  select = "close")
		self.x_y_log("{0} #4-40:".format(screw_letter),
		  screw_x, screw_y)

	screw_records = ( \
	  (self.screw_n_, "N", mb.screw_n_x_l, mb.screw_n_y_l), \
	  (self.screw_o_, "O", mb.screw_o_x_l, mb.screw_o_y_l), \
	  (self.screw_p_, "P", mb.screw_p_x_l, mb.screw_p_y_l), \
	  (self.screw_q_, "Q", mb.screw_q_x_l, mb.screw_q_y_l)  )
	for screw_record in screw_records:
	    screw = screw_record[0]
	    screw_letter = screw_record[1]
	    screw_x = screw_record[2]
	    screw_y = screw_record[3]
	    screw.configure(
	      comment = "Screw {0}".format(screw_letter),
	      start = P(screw_x, screw_y, mb.t.z),
	      end = P(screw_x, screw_y, mb.b.z),
	      flags = "#6-32")
	    screw.drill(
	      part = mb,
	      select = "close")
	    screw.drill(
	      part = pcb,
	      select = "close")
	    self.x_y_log("{0} #6-32".format(screw_letter), screw_x, screw_y)

	self._pcb_holes_file.close();

	# Place all the motor mount screws:
	motor_records = ( \
	  (twist_motor, "Twist", self.twist_screw_a_, self.twist_screw_b_),
	  (drive_motor, "Drive", self.drive_screw_a_, self.drive_screw_b_) )
	for motor_record in motor_records:
	    motor = motor_record[0]
	    motor_label = motor_record[1]
	    screw_a = motor_record[2]
	    screw_b = motor_record[3]
	    
	    # There are 4 possible mount angles -- *angle0*, *angle1*,
	    # *angle2*, and *angle3*.  We only do *angle0* and *angel2*.
	    # For one of the motors *angle1* interferes with the motor
	    # base "north" wall.  For, both motors, angle3 is too close
	    # to the motor hub:
	    angle0 = motor.mount_angle0_a
	    angle2 = motor.mount_angle2_a

	    motor_x = motor.x_l
	    motor_y = motor.y_l
	    mount_angle0 = motor.mount_angle0_a
	    mount_angle2 = motor.mount_angle2_a
	    mount_radius = motor.mount_hole_radius_l

	    #print("{0} motor=({1},{2})".format(motor_label, motor_x, motor_y))

	    screw_records = ( \
	      ("A", mount_angle0, screw_a),
	      ("B", mount_angle2, screw_b) )

	    # Do each screw:
	    for screw_record in screw_records:
		screw_label = screw_record[0]
		angle = screw_record[1]
		screw = screw_record[2]
		assert isinstance(screw, Fastener)

		x = motor_x + mount_radius.cosine(angle)
		y = motor_y + mount_radius.sine(angle)
		z = motor_z = mb.b.z - L(mm = 10)

		#print("{0} screw {1}=({2},{3})".
		#  format(motor_label, screw_label, x, y))

		screw.configure(
		  comment = "{0} Screw {1}".format(motor_label, screw_label),
		  start = P(x, y, mb.z3_l),
		  end = P(x, y, motor_z),
		  flags = "M3x.05")

		# Make a recess hole for each screw head:
		washer_diameter = L(mm = 7.00)
		mb.hole(comment = "{0} Screw {1} Recess".
		  format(motor_label, screw_label),
		  diameter = washer_diameter + L(mm = 1.0),
		  start = P(x, y, mb.z4_l),
		  end = P(x, y, mb.z3_l),
		  flags = "f")

		# Now drill the holes:
		screw.drill(motor, select = "thread")
		screw.drill(mb, select = "close")

    def x_y_log(self, label, x, y):
	x1 = self._x1
	y1 = self._y1
	self._pcb_holes_file.write("{0}: pcb_x={1} pcb_y={2} x={3} y={4}\n".
	  format(label, L(mm = 150) - (x - x1), L(mm = 150) - (y - y1), x, y))

class PCB(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	zero = L()
	self.dz_l = zero
	self.x1_l = zero
	self.y1_l = zero
	self.x2_l = zero
	self.y2_l = zero
	self.z_bottom_l = zero

    def configure(self, dz = None, x1 = None, y1 = None, x2 = None, y2 = None,
      z_bottom = None):
	# Check argument types:
	none_type = type(None)
	assert type(x1) == none_type or isinstance(x1, L)
	assert type(y1) == none_type or isinstance(y1, L)
	assert type(x2) == none_type or isinstance(x2, L)
	assert type(y2) == none_type or isinstance(y2, L)
	assert type(dz) == none_type or isinstance(dz, L)
	assert type(z_bottom) == none_type or isinstance(z_bottom, L)

	# Load up *self*:
	if isinstance(dz, L):
	    self.dz_l = dz
	if isinstance(x1, L):
	    self.x1_l = x1
	if isinstance(y1, L):
	    self.y1_l = y1
	if isinstance(x2, L):
	    self.x2_l = x2
	if isinstance(y2, L):
	    self.y2_l = y2
	if isinstance(z_bottom, L):
	    self.z_bottom_l = z_bottom

    def construct(self):
	z0 = self.z_bottom_l
	z1 = z0 + self.dz_l
	#print("PCB.construct: z0={0} z1={1}".format(z0, z1))

	x1 = self.x1_l
	y1 = self.y1_l
	x2 = self.x2_l
	y2 = self.y2_l
	center_z0 = P((x1 + x2) / 2, (y1 + y2) / 2, z0)
	center_z1 = P((x1 + x2) / 2, (y1 + y2) / 2, z1)

	self.as50xx_z_bottom_l = as50xx_z_bottom = z0 - L(mm = 0.90)
	self.magnet_z_top_l = as50xx_z_bottom - L(mm = 1.00)

	z = zero = L()
	outer_contour = Contour()
	outer_contour.bend_append(P(x2, y1, z), z)
	outer_contour.bend_append(P(x2, y2, z), z)
	outer_contour.bend_append(P(x1, y2, z), z)
	outer_contour.bend_append(P(x1, y1, z), z)

	inner_contour = outer_contour.adjust(
	  delta = -L(mm = 10),
	  maximum_radius = z,
	  start = center_z0,
	  end = center_z1)

	self.extrude(comment = "PCB block",
	  material = Material("pcb", "fr4"),
	  color = Color("light_green"),
	  outer_contour = outer_contour,
	  inner_contours = [inner_contour],
	  start = center_z0,
	  end = center_z1)

class Pulley(Part):
    def __init__(self, up):
	# Initialize *Part* super class:
	Part.__init__(self, up)

	# Initialize configuration varaibles for *self*:
	zero = L()
	self.belt_class_s = "MXL"
	self.belt_width_l = L(inch = "1/8")
	self.belt_width_extra_l = zero
	self.bearing_diameter_l = zero
	self.bearing_sides_i = 16
	self.bearing_width_l = zero
	self.holes_radius_l = zero
	self.lip_extra_l = L(inch = "1/8")
	self.lip_height_l = L(inch = 10)
	self.lib_name_s = "No_Name"
	self.pitch_l = zero
	self.set_screw_dz_l = zero
	self.shaft_diameter_l = zero
	self.teeth_count_i = 20
	self.z_bottom_l = zero
	self.x_l = zero
	self.y_l = zero

    def configure(self,
      bearing_diameter = None, bearing_sides = None, bearing_width = None,
      belt_class = None, belt_width = None, belt_width_extra = None,
      holes_radius = None, lip_extra = None, lip_height = None, name = None,
      set_screw_dz = None, shaft_diameter = None, teeth_count = None,
      x = None, y = None, z_bottom = None):

	# Check argument types:
	zero = L()
	none_type = type(None)
	assert type(bearing_diameter) == none_type or \
	  isinstance(bearing_diameter, L)
	assert type(bearing_sides) == none_type or \
	  isinstance(bearing_sides, int)
	assert type(bearing_width) == none_type or isinstance(bearing_width, L)
      	assert type(belt_class) == none_type or  isinstance(belt_class, str)
	assert type(belt_width) == none_type or \
	  isinstance(belt_width, L) and belt_width >= zero
	assert type(belt_width_extra) == none_type or \
	  isinstance(belt_width_extra, L) and belt_width_extra >= zero
	assert type(holes_radius) == none_type or \
	  isinstance(holes_radius, L)
	assert type(lip_extra) == none_type or \
	  isinstance(lip_extra, L) and lip_extra >= zero
	assert type(lip_height) == none_type or isinstance(lip_height, L)
	assert type(name) == none_type or isinstance(name, str)
	assert type(set_screw_dz) == none_type or isinstance(set_screw_dz, L)
	assert type(shaft_diameter) == none_type or \
	  isinstance(shaft_diameter, L)
	assert type(teeth_count) == none_type or isinstance(teeth_count, int)
	assert type(x) == none_type or isinstance(x, L)
	assert type(y) == none_type or isinstance(y, L)
	assert type(z_bottom) == none_type or isinstance(z_bottom, L)

	# Remember the arguments;
	if isinstance(bearing_diameter, L):
	    self.bearing_diameter_l = bearing_diameter
	if isinstance(bearing_sides, int):
	    self.bearing_sides_i = bearing_sides
	if isinstance(bearing_width, L):
	    self.bearing_width_l = bearing_width
	if isinstance(belt_class, str):
	    self.belt_class_s = belt_class
	if isinstance(belt_width, L):
	    self.belt_width_l = belt_width
	if isinstance(belt_width_extra, L):
	    self.belt_width_extra_l = belt_width_extra
	if isinstance(holes_radius, L):
	    self.holes_radius_l = holes_radius
	if isinstance(lip_extra, L):
	    self.lip_extra_l = lip_extra
	if isinstance(lip_height, L):
	    self.lip_height_l = lip_height
	if isinstance(name, str):
	    self.name_s = name
	if isinstance(set_screw_dz, L):
	    self.set_screw_dz_l = set_screw_dz
	if isinstance(shaft_diameter, L):
	    self.shaft_diameter_l = shaft_diameter
	if isinstance(teeth_count, int):
	    self.teeth_count_i = teeth_count
	if isinstance(x, L):
	    self.x_l = x
	if isinstance(y, L):
	    self.y_l = y
	if isinstance(z_bottom, L):
	    self.z_bottom_l = z_bottom

    def construct(self):
	""" *Pulley* construct method. """

	# Grab some value out of *self*:
	bearing_diameter = self.bearing_diameter_l
	bearing_sides = self.bearing_sides_i
	bearing_width = self.bearing_width_l
	belt_class = self.belt_class_s
	belt_width = self.belt_width_l
	belt_width_extra = self.belt_width_extra_l
	lip_height = self.lip_height_l
	lip_extra = self.lip_extra_l
	name = self.name_s
	shaft_diameter = self.shaft_diameter_l
	set_screw_dz = self.set_screw_dz_l
	teeth_count = self.teeth_count_i
	z_bottom = self.z_bottom_l

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
	#       MXL         .080"   .030"   .020"  40 deg   .005"
	#       XL          .200"   .051"   .050"  50 deg   .015"
	#       L           .375"   .128"   .075"  40 deg   .020"
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

	self.pitch_l = pitch = L(inch = p)
	self.tooth_height_l = height = L(inch = h)

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
	self.tooth_radius_l = tooth_radius = \
	  pitch / (2 * (tooth_angle / 2).sine())
	#print("tooth_angle={0:d}, tooth_angle.sine()={1}".
	#  format(tooth_angle, (tooth_angle / 2).sine()))
	#print("pitch={0}, tooth_radius={1} circum={2} p*tr={3}".
	#  format(pitch, tooth_radius,
	#  tooth_radius * 2 * 3.1415629, pitch * teeth_count))

	# Compute the various Z heights:
	zero = L()
	z0 = z_bottom			# Bottom
	z1 = z0 + lip_height / 2	# Extend into lip to force a weld
	z2 = z1 + lip_height / 2	# Lip top and gear bottom
	z3 = z2 + belt_width_extra / 2 	# Gear top
	z4 = z3 + belt_width / 2
	z5 = z4 + belt_width / 2
	z6 = z5 + belt_width_extra / 2	# Gear top
	#print(("Pulley.construct(): " +
	#  "z0={0} z1={1} z2={2} z3={3} z4={4} z5={5} z6={6}").
	#  format(z0, z1, z2, z3, z4, z5, z6))

	shaft_radius = shaft_diameter / 2
	bearing_radius = bearing_diameter / 2

	# Compute various locations along the radius
	r0 = shaft_radius
	r1 = bearing_radius
	r2 = (tooth_radius - height + bearing_radius) / 2
	r3 = tooth_radius - 2 * height
	r4 = tooth_radius - height
	r5 = tooth_radius
	r6 = tooth_radius + lip_extra
	#print("r0={0} r1={1} r2={2} r3={3} r4={4} r5={5} r6={6}".
	#  format(r0, r1, r2, r3, r4, r5, r6))

	if bearing_radius > zero:
	    self.holes_radius_l = holes_radius =  r2
	self.tooth_diameter_l = tooth_diameter = 2 * tooth_radius
	self.lip_diameter_l = lip_diameter = 2 * r6

	x = self.x_l
	y = self.y_l
	self.cylinder(comment = "{0} Pulley Lip".format(name),
	  material = Material("plastic", "ABS"),
	  color = Color("crimson"),
	  diameter = lip_diameter,
	  start = P(x, y, z0),
	  end = P(x, y, z2),
	  sides = teeth_count)
	#print("Pulley.construct():after lip:dx={0} dy={1} dz={2} x={3} y={4}".
	#  format(self.dx, self.dy, self.dz, x, y))
	
	quick = False
	#quick = True
	if quick:
	    self.cylinder(comment = "{0} Pulley Gear".format(name),
	    diameter = tooth_radius * 2,
	    start = P(zero, zero, z1),
	    end = P(zero, zero, z6),
	    sides = teeth_count)
	else:
	    # Compute *outer_path*:
	    outer_contour = Contour()
	    for index in range(teeth_count):
		# Compute the 4 angles:
		angle0 = tooth_angle * (index + 0.00)
		angle1 = tooth_angle * (index + w_ratio)
		angle2 = tooth_angle * (index + w_ratio + x2_ratio / 2)
		angle3 = \
		  tooth_angle * (index + w_ratio + x2_ratio / 2 + pw2x_ratio)

		# Add everything to the *outer_path*:
		outer_contour.bend_append(
		  P(x + r4.cosine(angle0), y + r4.sine(angle0), zero), zero)
		outer_contour.bend_append(
		  P(x + r4.cosine(angle1), y + r4.sine(angle1), zero), zero)
		outer_contour.bend_append(
		  P(x + r3.cosine(angle2), y + r3.sine(angle2), zero), zero)
		outer_contour.bend_append(
		  P(x + r3.cosine(angle3), y + r3.sine(angle3), zero), zero)

	    # Now extrude the shape:
	    self.extrude(comment = "{0} Gear".format(name),
	      outer_contour = outer_contour,
	      start = P(x, y, z1),
	      end = P(x, y, z6))
	    #print("Pulley.construct():after extrd:dx={0} dy={1} dz={2}".
	    #  format(self.dx, self.dy, self.dz))

	if shaft_diameter > zero:
	    self.hole(comment = "{0} Shaft Hole".format(name),
	      diameter = shaft_diameter,
	      start = P(x, y, z6),
	      end = P(x, y, z0),
	      flags = "t")

	#print("bearing_diameter={0} bearing_width={1}".
	#  format(bearing_diameter, bearing_width))
	if bearing_diameter > zero and bearing_width > zero:
	    self.hole(comment = "{0} Bearing Hole".format(name),
	      diameter = bearing_diameter,
	      start = P(x, y, z6),
	      end = P(x, y, z6 - bearing_width),
	      flags = "f",
	      sides = bearing_sides)

	# Do the set screw hole if requested:
	if set_screw_dz > zero:
	    #print("{0} Set_Screw: x={1} y={2} z0={3} set_screw_dz={4}".
	    #  format(name, x, y, z0, set_screw_dz))
	    self.hole(comment = "{0} Set Screw Hole".format(name),
	      diameter = L(inch = .0890),	# #43 drill for #4-40 tap
	      start = P(x, y + lip_diameter, z0 + set_screw_dz),
	      end = P(x, y, z0 + set_screw_dz),
	      flags = "t")

class Pulley_Top(Part):
    def __init__(self, up):
	# Initialize *Part* super class:
	Part.__init__(self, up)

	# Initialize configuration varaibles for *self*:
	zero = L()
	one = L(mm = 1)
	self.hub_diameter_l = zero
	self.hub_dz_l = zero
	self.lip_diameter_l = zero
	self.lip_dz_l = zero
	self.magnet_diameter_l = zero
	self.magnet_dz_l = zero
	self.set_screw_dz_l = zero
	self.shaft_diameter_l = -one
	self.x_l = zero
	self.y_l = zero
	self.z_bottom_l = zero

    def configure(self, hub_diameter = None, hub_dz = None, ip_diameter = None,
      lip_diameter = None, magnet_diameter = None, magnet_dz = None,
      lip_dz = None, shaft_diameter = None, set_screw_dz = None,
      x = None, y = None, z_bottom = None):

	# Check argument types:
	zero = L()
	none_type = type(None)
	assert type(hub_diameter) == none_type or isinstance(hub_diameter, L)
	assert type(hub_dz) == none_type or isinstance(hub_dz, L)
	assert type(lip_diameter) == none_type or isinstance(lip_diameter, L)
	assert type(lip_dz) == none_type or isinstance(lip_dz, L)
	assert type(magnet_diameter) == none_type or \
	  isinstance(magnet_diameter, L)
	assert type(magnet_dz) == none_type or isinstance(magnet_dz, L)
	assert type(set_screw_dz) == none_type or isinstance(set_screw_dz, L)
	assert type(shaft_diameter) == none_type or \
	  isinstance(shaft_diameter, L)
	assert type(x) == none_type or isinstance(x, L)
	assert type(y) == none_type or isinstance(y, L)
	assert type(z_bottom) == none_type or isinstance(z_bottom, L)

	# Remember the arguments:
	if isinstance(hub_diameter, L):
	    self.hub_diameter_l = hub_diameter
	if isinstance(hub_dz, L):
	    self.hub_dz_l = hub_dz
	if isinstance(lip_diameter, L):
	    self.lip_diameter_l = lip_diameter
	if isinstance(lip_dz, L):
	    self.lip_dz_l = lip_dz
	if isinstance(magnet_diameter, L):
	    self.magnet_diameter_l = magnet_diameter
	if isinstance(magnet_dz, L):
	    self.magnet_dz_l = magnet_dz
	if isinstance(set_screw_dz, L):
	    self.set_screw_dz_l = set_screw_dz
	if isinstance(shaft_diameter, L):
	    self.shaft_diameter_l = shaft_diameter
	if isinstance(x, L):
	    self.x_l = x
	if isinstance(y, L):
	    self.y_l = y
	if isinstance(z_bottom, L):
	    self.z_bottom_l = z_bottom

    def construct(self):
	""" *Pulley_Top* construct method. """

	z0 = self.z_bottom_l				# Bottom
	z1 = z0 + self.lip_dz_l / 2			# Half way thru lip
	z2 = z1 + self.lip_dz_l / 2			# Lip top
	z3 = z2 + self.hub_dz_l / 2			# Set screw hight
	z4 = z3 + self.hub_dz_l / 2 - self.magnet_dz_l	# Magnet bottom
	z5 = z4 + self.magnet_dz_l			# Top

	# Grab some value out of *self*:
	zero = L()
	x = self.x_l
	y = self.y_l
	#print("Pulley_Top.construct():x={0} y={1} z0={2} z1={3} z2={4} z3={5}".
	#  format(x, y, z0, z1, z2, z3))

	# Do the lip:
	lip_diameter = self.lip_diameter_l
	self.cylinder(comment = "Lip Cylinder",
	  material = Material("plastic", "ABS"),
	  color = Color("pink"),
	  diameter = lip_diameter,
	  start = P(x, y, z0),
	  end = P(x, y, z2),
	  sides = 60)

	# Do the hub (weld it into the lip using *z1*):
	hub_diameter = self.hub_diameter_l
	if hub_diameter > zero:
	    self.cylinder(comment = "Hub Cylinder",
	    diameter = hub_diameter,
	    start = P(x, y, z1),
	    end = P(x, y, z5),
	    sides = 30)

	# Do the shaft hole:
	shaft_diameter = self.shaft_diameter_l
	if isinstance(shaft_diameter, L) and shaft_diameter > zero:
	    self.hole(comment = "Shaft Hole",
	      diameter = shaft_diameter,
	      start = P(x, y, z5),
	      end = P(x, y, z0),
	      flags = "t",
	      sides = 60)

	# Do the set screw hole:
	set_screw_dz = self.set_screw_dz_l
	if set_screw_dz > zero:
	    self.hole(comment = "Set Screw Hole",
	      diameter = L(inch = .0890),	# #43 drill for #4-40 tap
	      start = P(x, y + lip_diameter, z0 + set_screw_dz),
	      end = P(x, y, z0 + set_screw_dz),
	      flags = "t")

	# Do magnet hole:
	magnet_diameter = self.magnet_diameter_l
	if magnet_diameter > zero:
	    self.hole(comment = "Magnet hole",
	      diameter = magnet_diameter,
	      start = P(x, y, z5),
	      end = P(x, y, z4),
	      flags = "f")

class Shim(Part):
    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):

	wheel_assembly = self.up
	bearing = wheel_assembly.south_bearing_
	bearing_diameter = bearing.diameter_l
	bearing_bore = bearing.bore_l

	self.bore_l = bore = bearing_bore = L(inch = 0.388)
	self.diameter_l = diameter = L(inch = 0.625)
	self.width_l = width = L(inch = 0.055)
	
	zero = L()
	self.cylinder(comment = "Shim Cylinder",
	  material = Material("steel", "stainless"),
	  color = Color("coral"),
	  diameter = diameter,
	  start = P(zero, zero, zero),
	  end = P(zero, zero, width))

class Test(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.pulley_top_ = Pulley_Top(self)
	#self.test_block_ = Test_Block(self)
	#self.test_block_screws_ = Test_Block_Screws(self)
	#self.bearing_ = Bearing(self)

    def construct(self):
	pulley_top = self.pulley_top_
	pulley_top.configure(
	  hub_diameter = L(mm = 20),
	  hub_dz = L(mm = 8),
	  magnet_diameter = L(mm = 8),
	  magnet_dz = L(mm = 2.50),
	  set_screw = True,
	  shaft_diameter = L(mm = 6),
	  lip_diameter = L(mm = 25),
	  lip_dz = L(mm = 2))
	#bearing = self.bearing_
	#test_block = self.test_block_
	#zero = L()
	#bearing.place(translate =
	#  P(zero, zero, test_block.bearing_lip_dz_l + bearing.width_l / 2))

class Test_Block(Part):
    def construct(self):
	test = self.up
	bearing = test.bearing_

	dx = L(inch = 2.000)
	dy = L(inch = 2.000)

	self.pocket_wall_dx_l = pocket_wall_dx = L(mm = 6.00)
	self.pocket_wall_dy_l = pocket_wall_dy = L(mm = 6.00)
	self.pocket_dx_l = pocket_dx = dx - 2 * pocket_wall_dx
	self.pocket_dy_l = pocket_dy = dy - 2 * pocket_wall_dy
	self.pocket_dz_l = pocket_dz = L(mm = 14.00)
	self.bearing_lip_dz_l = bearing_lip_dz = L(mm = 2.00)

	x0 = -dx / 2
	x1 = x0 + pocket_wall_dx
	x2 = x1 + pocket_dx / 2
	x3 = x2 + pocket_dx / 2
	x4 = x3 + pocket_wall_dx

	y0 = -dy / 2
	y1 = y0 + pocket_wall_dy
	y2 = y1 + pocket_dy / 2
	y3 = y2 + pocket_dy / 2
	y4 = y3 + pocket_wall_dy

	z0 = L()
	z1 = z0 + bearing_lip_dz
	z2 = z1 + bearing.width_l
	z3 = z2 + pocket_dz

	self.block(comment = "Test Block",
	  material = Material("plastic", "ABS"),
	  color = Color("cyan"),
	  corner1 = P(x0, y0, z0),
	  corner2 = P(x4, y4, z3),
	  top = "t")
	self.simple_pocket(comment = "Pocket",
	  bottom_corner = P(x1, y1, z2),
	  top_corner =    P(x3, y3, z3),
	  pocket_top = "t")
	self.hole(comment = "bearing hole",
	  diameter = bearing.diameter_l,
	  start = P(x2, y2, z3),
	  end = P(x2, y2, z1),
	  sides = bearing.sides_i,
	  flags = "f")
	self.hole(comment = "shaft hole",
	  diameter = (bearing.diameter_l + bearing.bore_l)/2,
	  start = P(x2, y2, z3),
	  end = P(x2, y2, z0),
	  flags = "t")

class Test_Block_Screws(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.screw_2_56_bottom_ = Fastener(self)
	self.screw_4_40_bottom_ = Fastener(self)
	self.screw_6_32_bottom_ = Fastener(self)
	self.screw_10_24_bottom_ = Fastener(self)
	self.screw_2_56_side_ = Fastener(self)
	self.screw_4_40_side_ = Fastener(self)
	self.screw_6_32_side_ = Fastener(self)
	self.screw_10_24_side_ = Fastener(self)

    def construct(self):
        test = self.up
	test_block = test.test_block_
	bearing = test.bearing_

	pocket_wall_dx = test_block.pocket_wall_dx_l
	pocket_wall_dy = test_block.pocket_wall_dy_l
	pocket_dx = test_block.pocket_dx_l
	pocket_dy = test_block.pocket_dy_l

	x0 = test_block.w.x
	x1 = x0 + pocket_wall_dx
	x2 = x1 + pocket_dx * 0.15
	x3 = x2 + pocket_dx * 0.35
	x4 = x3 + pocket_dx * 0.35
	x5 = x4 + pocket_dx * 0.15
	x6 = x5 + pocket_wall_dx
	#print("x0={0} x1={1} x2={2} x3={3} x4={4} x5={5} x6={6}".
	#  format(x0, x1, x2, x3, x4, x5, x6))

	y0 = test_block.s.y
	y1 = y0 + pocket_wall_dy
	y2 = y1 + pocket_dy * 0.15
	y3 = y2 + pocket_dy * 0.35
	y4 = y3 + pocket_dy * 0.35
	y5 = y4 + pocket_dy * 0.15
	y6 = y5 + pocket_wall_dy
	#print("y0={0} y1={1} y2={2} y3={3} y4={4} y5={5} y6={6}".
	#  format(y0, y1, y2, y3, y4, y5, y6))

	z0 = test_block.b.z
	z1 = z0 + test_block.bearing_lip_dz_l
	z2 = z1 + bearing.width_l
	z3 = z2 + test_block.pocket_dz_l / 2
	z4 = z3 + test_block.pocket_dz_l / 2
	#print("z0={0} z1={1} z2={2} z3={3} z4={4}".format(z0, z1, z2, z3, z4))

	screw_2_56_bottom  = self.screw_2_56_bottom_
	screw_4_40_bottom  = self.screw_4_40_bottom_
	screw_6_32_bottom  = self.screw_6_32_bottom_
	screw_10_24_bottom = self.screw_10_24_bottom_
	screw_2_56_side  = self.screw_2_56_side_
	screw_4_40_side  = self.screw_4_40_side_
	screw_6_32_side  = self.screw_6_32_side_
	screw_10_24_side = self.screw_10_24_side_

	screw_2_56_bottom.configure(
	  comment = "2-56 bottom",
	  flags = "#2-56:hi",
	  start = P(x2, y2, z0),
	  end = P(x2, y2, z2))
	screw_4_40_bottom.configure(
	  comment = "4-40 bottom",
	  flags = "#4-40:hi",
	  start = P(x2, y4, z0),
	  end = P(x2, y4, z2))
	screw_6_32_bottom.configure(
	  comment = "6-32 bottom",
	  flags = "#6-32:hi",
	  start = P(x4, y2, z0),
	  end = P(x4, y2, z2))
	screw_10_24_bottom.configure(
	  comment = "10-24 bottom",
	  flags = "#10-24:hi",
	  start = P(x4, y4, z0),
	  end = P(x4, y4, z2))

	zero = L()
	degrees30 = Angle(deg = 30)
	screw_2_56_side.configure(
	  comment = "2-56 side",
	  flags = "#2-56:hi",
	  start = P(x0, y3, z3),
	  end = P(x1, y3, z3))
	screw_4_40_side.configure(
	  comment = "4-40 side",
	  flags = "#4-40:hi",
	  start = P(x3, y0, z3),
	  end = P(x3, y1, z3),
	  sides_angle = degrees30)
	screw_6_32_side.configure(
	  comment = "6-32 side",
	  flags = "#6-32:hi",
	  start = P(x6, y3, z3),
	  end = P(x5, y3, z3))
	screw_10_24_side.configure(
	  comment = "10-24 side",
	  flags = "#10-24:hi",
	  start = P(x3, y6, z3),
	  end = P(x3, y5, z3),
	  sides_angle = degrees30)

	screw_2_56_bottom.drill(test_block, select = "close")
	screw_4_40_bottom.drill(test_block, select = "close")
	screw_6_32_bottom.drill(test_block, select = "close")
	screw_10_24_bottom.drill(test_block, select = "close")
	screw_2_56_side.drill(test_block, select = "close")
	screw_4_40_side.drill(test_block, select = "close")
	screw_6_32_side.drill(test_block, select = "close")
	screw_10_24_side.drill(test_block, select = "close")

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
	self.turn_table_bottom_screw0_ = Fastener(self, comment="screw0")
	self.turn_table_bottom_screw1_ = Fastener(self, comment="screw1")
	self.turn_table_bottom_screw2_ = Fastener(self, comment="screw2")
	self.turn_table_bottom_screw3_ = Fastener(self, comment="screw3")

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
		  flags = "#2-56:hi",
		  start = P(x, y, z3),
		  end = P(x, y, z1))

	# Now drill some holes in *gear_box_top*:
	screw0.drill(gear_box_top, select = "close")
	screw1.drill(gear_box_top, select = "close")
	screw2.drill(gear_box_top, select = "close")
	screw3.drill(gear_box_top, select = "close")

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
	  corner2 = P( turn_table_size / 2,  turn_table_size / 2, z1),
	  top = "t")

	# Bore the hole through the center:
	self.hole(comment = "Bore Hole",
	  diameter = turn_table_bore,
	  start = P(zero, zero, z1),
	  end =   P(zero, zero, z0),
	  flags = "t",
	  sides = 60)

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

class Twist_Motor_Pulley_Screws(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.screw0_ = Fastener(self, comment="screw0")
	self.screw1_ = Fastener(self, comment="screw1")
	self.screw2_ = Fastener(self, comment="screw2")
	self.screw3_ = Fastener(self, comment="screw3")

    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	twist_motor_pulley = wheel_assembly.twist_motor_pulley_
	twist_motor_pulley_top = wheel_assembly.twist_motor_pulley_top_

	# Create the 6 screws:
	screw0 = self.screw0_
	screw1 = self.screw1_
	screw2 = self.screw2_
	screw3 = self.screw3_
	screws = [screw0, screw1, screw2, screw3]

	# Compute some Z axis locations:
	z0 = twist_motor_pulley.b.z
	z1 = twist_motor_pulley_top.b.z + twist_motor_pulley_top.lip_dz_l

	pulley_x = twist_motor_pulley.x_l
	pulley_y = twist_motor_pulley.y_l

	holes_count = len(screws)
	angle_delta = Angle(deg = 360) / holes_count
	holes_radius = twist_motor_pulley.holes_radius_l
	#print("Twist_Motor_Pulley_Screws:holes_radius={0}".
	#  format(holes_radius))
	angle_adjust = Angle(deg = 45)
	for index in range(holes_count):
	    screw = screws[index]
	    angle = index * angle_delta + angle_adjust
	    x = pulley_x + holes_radius.cosine(angle)
	    y = pulley_y + holes_radius.sine(angle)
	    #print("holes_radius={0} x={1} y={2}".format(holes_radius, x, y))
	    screw.configure(comment = "Pulley Hole {0}".format(index),
	      flags = "#0-80:hi:fh",
	      start = P(x, y, z0),
	      end = P(x, y, z1),
	      sides_angle = Angle(deg = 30))
	    screw.drill(twist_motor_pulley, select = "close")
	    screw.drill(twist_motor_pulley_top, select = "close")

	#print("Pulley_Screw.dz={0:.3i}in".format(z3 - z0))

class Twist_Wheel_Pulley_Screws(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.screw0_ = Fastener(self, comment="screw0")
	self.screw1_ = Fastener(self, comment="screw1")
	self.screw2_ = Fastener(self, comment="screw2")
	self.screw3_ = Fastener(self, comment="screw3")

    def construct(self):
	# Grab some values from *wheel_assembly*:
	wheel_assembly = self.up
	gear_box = wheel_assembly.gear_box_
	gear_box_top = wheel_assembly.gear_box_top_
	twist_wheel_pulley = wheel_assembly.twist_wheel_pulley_
	twist_wheel_pulley_top = wheel_assembly.twist_wheel_pulley_top_

	# Create the 6 screws:
	screw0 = self.screw0_
	screw1 = self.screw1_
	screw2 = self.screw2_
	screw3 = self.screw3_
	screws = [screw0, screw1, screw2, screw3]

	# Compute some Z axis locations:
	z0 = gear_box.t.z - gear_box_top.base_dz_l
	z1 = z0 + gear_box_top.base_dz_l
	z2 = z1 + twist_wheel_pulley.dz
	z3 = z2 + twist_wheel_pulley_top.dz
	#print("Twist_Wheel_Pulley_Screws:z0={0} z1={1} z2={2}".
	#  format(z0, z1, z2))

	holes_count = len(screws)
	angle_delta = Angle(deg = 360) / holes_count
	holes_radius = twist_wheel_pulley.holes_radius_l
	angle_adjust = Angle(deg = 45)
	for index in range(holes_count):
	    screw = screws[index]
	    angle = index * angle_delta + angle_adjust
	    x = holes_radius.cosine(angle)
	    y = holes_radius.sine(angle)
	    #print("holes_radius={0} x={1} y={2}".format(holes_radius, x, y))
	    screw.configure(comment = "Pulley Hole {0}".format(index),
	      flags = "#0-80:hi:fh",
	      start = P(x, y, z3),
	      end = P(x, y, z0),
	      sides_angle = Angle(deg = 30))
	    screw.drill(gear_box_top, select = "close")
	    screw.drill(twist_wheel_pulley, select = "close")
	    screw.drill(twist_wheel_pulley_top, select = "close")

	#print("Pulley_Screw.dz={0:.3i}in".format(z3 - z0))

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

	#print("Horizontal_Shaft.dy={0:i}in".format(self.dy))

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
          end = P(zero, zero, L(mm = 108.5)))

	#print("Vertical_Shaft.dz={0:i}in".format(self.dz))

class Bearing_Starter(Part):
    def construct(self):
	zero = L()
	self.cylinder(comment = "Basic Cylinder",
	  material = Material("plastic", "ABS"),
	  color = Color("cyan"),
	  diameter = L(inch = 1.200),
	  start = P(),
	  end = P(zero, zero, L(inch = 0.500)),
	  sides = 60)
	self.hole(comment = "Bearing Hole",
	  diameter = L(inch = "7/8") + L(inch = 0.025),
	  start = self.t,
	  end = P(zero, zero, self.t.z - L(inch = "1/8")),
	  flags = "f",
	  sides = 60)

class Bearing_Ender(Part):
    def construct(self):
	zero = L()
	self.cylinder(comment = "Basic Cylinder",
	  material = Material("plastic", "ABS"),
	  color = Color("cyan"),
	  diameter = L(inch = "7/8"),
	  start = P(),
	  end = P(zero, zero, L(inch = "3/8")),
	  sides = 60)
	self.hole(comment = "Shaft_Hole",
	  diameter = L(inch = "3/8"),
	  start = self.t,
	  end = self.b,
	  flags = "t",
	  sides = 60)
	self.hole(comment = "Bearing Lip",
	  diameter = L(inch = "3/4"),
	  start = self.t,
	  end = P(zero, zero, self.t.z - L(inch = "1/16")),
	  flags = "f",
	  sides = 60)

class Bevel_Gear_Jig(Part):
    def construct(self):
	# Mat   Teeth  PD   OD  Bore  Height  Hub Diam  Hub Height  Maj Diam
	# 16 Pitch:
	# Nylon   16 1.000 1.090 .375 3/4      3/4       7/16        1-1/16
	zero = L()
	dx = L(inch = 1.500)
	dy = L(inch = 1.500)
	dz = L(inch = "7/16")

	pocket_dx = L(inch = "1/4")
	pocket_dy = L(inch = "1/4")

	hub_height = L(inch = "7/16")
	hub_diameter = L(inch = "3/4") + L(inch = .010)
	outside_diameter = L(inch = 1.090) + L(inch = .010)
	height = L(inch = "3/4")

	z0 = zero
	z1 = z0 + hub_height / 2
	z2 = z1 + hub_height / 2
	z3 = height

	self.block(comment = "Jig Block",
	  material = Material("plastic", "ABS"),
	  color = Color("cyan"),
	  corner1 = P(-dx/2, -dy/2, z0),
	  corner2 = P( dx/2,  dy/2, z3),
	  top = "t")
	self.simple_pocket(comment = "SW Pocket",
	  bottom_corner = P(-dx/2, -dy/2, z0),
	  top_corner = P(-dx/2 + pocket_dx, -dy/2 + pocket_dy, z3))
	self.simple_pocket(comment = "SE Pocket",
	  bottom_corner = P( dx/2, -dy/2, z0),
	  top_corner = P(dx/2 - pocket_dx,  -dy/2 + pocket_dy, z3))
	self.hole(comment = "Teeth Hole",
	  diameter = outside_diameter,
	  start = P(zero, zero, z3),
	  end = P(zero, zero, z2),
	  sides = 60)
	self.hole(comment = "Hub Hole",
	  diameter = hub_diameter, 
	  start = P(zero, zero, z2),
	  end = P(zero, zero, z0),
	  sides = 60)
	self.hole(comment = "Set Screw Hole",
	  diameter = L(inch = .0890),
	  start = P(zero, dy/2, z1),
	  end = P(zero, zero, z1),
	  flags = "t")
	
class Counter_Sink(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.block_ = Counter_Sink_Block(self)
	self.screw_ = Fastener(self)

    def construct(self):
	block = self.block_
	screw = self.screw_
	screw.configure(
	  comment = "screw",
	  material = Material("steel", "stainless"),
	  color = Color("black"),
	  flags = "#10-24:hi:fh",
	  start = block.t,
	  end = block.b)
	screw.drill(block, select = "close")

class Counter_Sink_Block(Part):
    def construct(self):
	zero = L()
	ten = L(mm = 10)
	self.block(comment = "Block",
	  material = Material("plastic", "ABS"),
	  color = Color("cyan"),
	  corner1 = P(-ten, -ten, -ten),
	  corner2 = P( ten,  ten,  ten),
	  top = "t")

if __name__ == "__main__":
    ezcad = EZCAD3(0, adjust = L(mm = -0.10))
    #ezcad = EZCAD3(0)
    #test = Synchro_Drive(None)
    test = Base(None)
    #test = Test(None)
    #test = Bearing_Starter(None)
    #test = Bearing_Ender(None)
    #test = Bevel_Gear_Jig(None)
    #test = Counter_Sink(None)
    test.process(ezcad)








