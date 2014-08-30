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
	self.dx_l = L(cm=25)
	self.dy_l = L(cm=25)
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
	outer_contour = Contour()
	outer_contour.bend_append(P(x0, y0, z0), z)
	outer_contour.bend_append(P(x0, y1, z0), z)
	outer_contour.bend_append(P(x1, y1, z0), z)
	outer_contour.bend_append(P(x1, y0, z0), z)
	
	x_pitch = L(cm=1)
	y_pitch = L(cm=1)
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
		hole_contour = Contour()
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

class Dual_Slot_Encoder(Part):
    def __init__(self, up):
	""" Dual Slot Encoder: initialize """
	Part.__init__(self, up)
	self.dual_slot_encoder_pcb_ = Dual_Slot_Encoder_PCB(self)
	self.slot_interrupter1_ = slot_interrupter1 = Slot_Interrupter(self)
	self.slot_interrupter2_ = slot_interrupter2 = Slot_Interrupter(self)
	self.disk_offset_dy_l = L(mm=3)
	self.y_center_l = L()
	self.z_center_l = L()

    def configure(self, y_center = None, z_center = None):
	""" Dual Slot Encoder: configure """
	none_type = type(None)
	assert type(y_center) == none_type or isinstance(y_center, L)
	assert type(z_center) == none_type or isinstance(z_center, L)

	# Load up any specified values:    
	if isinstance(y_center, L):
	    self.y_center_l = y_center
	if isinstance(z_center, L):
	    self.z_center_l = z_center
	
    def construct(self):
	""" Dual Slot Encoder: construct """
	disk_offset_dy = self.disk_offset_dy_l
	dual_slot_encoder_pcb = self.dual_slot_encoder_pcb_
	motor_assembly = self.up
	slot_interrupter1 = self.slot_interrupter1_
	slot_interrupter2 = self.slot_interrupter2_
	y_center = self.y_center_l
	z_center = self.z_center_l

	dual_slot_encoder_pcb.configure(
	  y_center = y_center, z_center = z_center)
	slot_interrupter1.configure(x_center =  slot_interrupter1.dx_l/2,
	  y_center = y_center + disk_offset_dy, z_center = z_center)
	slot_interrupter2.configure(x_center = -slot_interrupter2.dx_l/2,
	  y_center = y_center + disk_offset_dy, z_center = z_center)

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
	self.y_center_l = L(0)
	self.z_center_l = L(0)

    def configure(self, y_center = None, z_center = None):
	""" Dual_Slot_Encoder_PCB: configure """
	none_type = type(None)
	assert type(y_center) == none_type or isinstance(y_center, L)
	assert type(z_center) == none_type or isinstance(z_center, L)

	# Load up any specified values:
	if isinstance(y_center, L):
	    self.y_center_l = y_center
	if isinstance(z_center, L):
	    self.z_center_l = z_center

    def construct(self):
	""" Dual_Slot_Encoder_PCB: construct """
	diameter = self.diameter_l
	dx = self.dx_l
	dy = self.dy_l
	dz = self.dz_l
	y_center = self.y_center_l
	z_center = self.z_center_l

	x0 = -dx/2
	x1 =  dx/2
	y0 = y_center - dy/2
	y1 = y_center + dy/2
	z0 = z_center
	z1 = z_center + dz

	# Create the PCB:
	z = L()
	self.block(comment = "Dual_Slot_Encoders_PCB",
	  material = Material("fiberglass", "FR4"),
	  color = Color("green"),
	  corner1 = P(-dx/2, y0, z0),
	  corner2 = P( dx/2, y1, z1),
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
	Part.__init__(self, up)
	self.diameter_l = diameter = L(mm=50)
	radius = diameter / 2
	self.slot_inner_radius_l = radius - L(mm=4.00)
	self.slot_outer_radius_l = radius - L(mm=0.50)
	self.slots_count_i = 36
	self.dy_l = L(mm=1.00)
	self.y_center_l = L(0)

    def configure(self, y_center = None):
	none_type = type(None)
	assert type(y_center) == none_type or isinstance(y_center, L)
	if isinstance(y_center, L):
	    self.y_center_l = y_center

    def construct(self):	
	diameter = self.diameter_l
	radius = diameter/2
	slot_inner_radius = self.slot_inner_radius_l
	slot_outer_radius = self.slot_outer_radius_l
	slots_count = self.slots_count_i
	dy = self.dy_l
	y_center = self.y_center_l

	y0 = y_center - dy/2
	y1 = y_center
	y2 = y_center + dy/2

	zero = L()
	two_pi = Angle.TWO_PI
	slot_angle = Angle(rad=two_pi) / slots_count
	half_slot_angle = slot_angle / 2
	#print("slot_angle={0} half_slot_angle={1}".
	#  format(slot_angle, half_slot_angle))

	slot_contours = []
	outer_contour = Contour()
	for slot_index in range(slots_count):
	    # Compute the angles for this iteration:
	    angle0 = slot_angle * slot_index
	    angle1 = angle0 + half_slot_angle

	    # Append (x,y) for the outer contour:
	    x0 = radius.cosine(angle0)
	    z0 = radius.sine(angle0)
	    x1 = radius.cosine(angle1)
	    z1 = radius.sine(angle1)
	    outer_contour.bend_append(P(x0, y0, z0), zero,
	      name="[{0}].xz0".format(slot_index))
	    outer_contour.bend_append(P(x1, y0, z1), zero,
	      name="[{0}].xz1".format(slot_index))

	    # Compute X/Y coordinates for the slot corners:
	    x0 = slot_inner_radius.cosine(angle0)
	    z0 = slot_inner_radius.sine(angle0)
	    x1 = slot_inner_radius.cosine(angle1)
	    z1 = slot_inner_radius.sine(angle1)
	    x2 = slot_outer_radius.cosine(angle1)
	    z2 = slot_outer_radius.sine(angle1)
	    x3 = slot_outer_radius.cosine(angle0)
	    z3 = slot_outer_radius.sine(angle0)
	
	    #print("[{0}]:angle0={1} angle1={2}".
	    #  format(slot_index, angle0, angle1))
	    #print("[{0}]:x0={1} z0={2}".format(slot_index, x0, z0))
	    #print("[{0}]:x1={1} z1={2}".format(slot_index, x1, z1))
	    #print("[{0}]:x2={1} z2={2}".format(slot_index, x2, z2))
	    #print("[{0}]:x3={1} z3={2}".format(slot_index, x3, z3))

	    # Append a *slot_contour* to *slot_contours*:
	    slot_contour = Contour()
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
	  material = Material("Plastic", "Derilin"),
	  color = Color("yellow"),
	  outer_contour = outer_contour,
	  inner_contours = slot_contours,
	  start = P(0, y0, 0),
	  end   = P(0, y2, 0))

	self.hole(comment = "Screw_Hole",
	  diameter = L(mm=3),
	  start = P(0, y2, 0),
	  end =   P(0, y0, 0),
	  flags = "t")

	self.dxf_write(center = P(0, y1, 0), plane_normal = P(0, L(mm=1), 0))

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
	self.y_center_l = L()
	self.z_center_l = L()

    def configure(self, dx = None, dy = None, dz = None,
      tongue_dx = None, y_center = None, z_center = None):
	""" Encoder_Mount: configure """
	none_type = type(None)
	assert type(dx) == none_type or isinstance(dx, L)
	assert type(dy) == none_type or isinstance(dy, L)
	assert type(dz) == none_type or isinstance(dz, L)
	assert type(tongue_dx) == none_type or isinstance(tongue_dx, L)
	assert type(y_center) == none_type or isinstance(y_center, L)
	assert type(z_center) == none_type or isinstance(z_center, L)

	if isinstance(dx, L):
	    self.dx_l = dx
	if isinstance(dy, L):
	    self.dy_l = dy
	if isinstance(dz, L):
	    self.dz_l = dz
	if isinstance(tongue_dx, L):
	    self.tongue_dx_l = tongue_dx
	if isinstance(y_center, L):
	    self.y_center_l = y_center
	if isinstance(z_center, L):
	    self.z_center_l = z_center

    def construct(self):
	""" Encoder_Mount: construct """
	dx = self.dx_l
	dy = self.dy_l
	dz = self.dz_l
	tongue_dx = self.tongue_dx_l
	y_center = self.y_center_l
	z_center = self.z_center_l
	thick = self.thick_l
	thin = self.thin_l

	# Go get some related *Part*'s:
	motor_assembly = self.up
	dual_slot_encoder = motor_assembly.dual_slot_encoder_
	dual_slot_encoder_pcb = dual_slot_encoder.dual_slot_encoder_pcb_

	# Grab some values from the Encoder PCB:
	pcb_y_center = dual_slot_encoder_pcb.y_center_l
	pcb_dx = dual_slot_encoder_pcb.dx_l
	pcb_dy = dual_slot_encoder_pcb.dy_l
	holes_pitch_dx = dual_slot_encoder_pcb.holes_pitch_dx_l
	holes_pitch_dy = dual_slot_encoder_pcb.holes_pitch_dy_l

	# Compute some X locations in ascending order:
	x0 = -dx/2			# West mount edge
	x1 = -tongue_dx/2		# West tongue edge
	x2 = -pcb_dx/2 + L(mm=2.5)
	x3 = -pcb_dx/2 + L(mm=8.5)
	x4 =  pcb_dx/2 - L(mm=8.5)
	x5 =  pcb_dx/2 - L(mm=2.5)
	x6 =  tongue_dx/2		# East tongue edge
	x7 =  dx/2			# East mount edge

	# Compute some Y locations in ascending order:
	y0 = y_center - dy/2
	y1 = y_center - dy/2 + L(mm=3)
	y2 = pcb_y_center - pcb_dy/2 + L(mm=2.5)
	y3 = pcb_y_center - pcb_dy/2 + L(mm=8.5)
	y4 = pcb_y_center + pcb_dy/2 - L(mm=8.5)
	y5 = pcb_y_center + pcb_dy/2 - L(mm=2.5)
	y6 = y_center + dy/2 - L(mm=3)
	y7 = y_center + dy/2
	print("Encoder_Mount:construct: y0={0} y_center={1} y7={2} dy={3}".
	  format(y0, y_center, y7, dy))

	# Compute some Z locations in ascending order:
	z0 = z_center - dz/2		# Bottom surface
	z1 = z_center + dz/2		# Top surface

	# Compute the outer contour:
	zero = L()
	outer_contour = Contour()
	# Corner: x0, y0:
	outer_contour.bend_append(P(x1, y0, z0), zero)
	outer_contour.bend_append(P(x1, y1, z0), zero)
	outer_contour.bend_append(P(x0, y1, z0), zero)

	# Corner: x0, y7:
	outer_contour.bend_append(P(x0, y6, z0), zero)
	outer_contour.bend_append(P(x1, y6, z0), zero)
	outer_contour.bend_append(P(x1, y7, z0), zero)

	# Corner: x7, x7:
	outer_contour.bend_append(P(x6, y7, z0), zero)
	outer_contour.bend_append(P(x6, y6, z0), zero)
	outer_contour.bend_append(P(x7, y6, z0), zero)

	# Corner: x7, y0:
	outer_contour.bend_append(P(x7, y1, z0), zero)
	outer_contour.bend_append(P(x6, y1, z0), zero)
	outer_contour.bend_append(P(x6, y0, z0), zero)

	# Compute the PCB cut out inner contour:
	pcb_cutout = Contour()
	curve_radius = L(mm=0.5)
	# Corner (x2, y2):
	pcb_cutout.bend_append(P(x3, y2, z0), curve_radius)
	pcb_cutout.bend_append(P(x3, y3, z0), curve_radius)
	pcb_cutout.bend_append(P(x2, y3, z0), curve_radius)
	# Corner (x2, y3):
	pcb_cutout.bend_append(P(x2, y4, z0), curve_radius)
	pcb_cutout.bend_append(P(x3, y4, z0), curve_radius)
	pcb_cutout.bend_append(P(x3, y5, z0), curve_radius)
	# Corner (x3, y3):
	pcb_cutout.bend_append(P(x4, y5, z0), curve_radius)
	pcb_cutout.bend_append(P(x4, y4, z0), curve_radius)
	pcb_cutout.bend_append(P(x5, y4, z0), curve_radius)
	# Corner (x3, y2):
	pcb_cutout.bend_append(P(x5, y3, z0), curve_radius)
	pcb_cutout.bend_append(P(x4, y3, z0), curve_radius)
	pcb_cutout.bend_append(P(x4, y2, z0), curve_radius)

	# Do the extrusion:
	self.extrude(comment = "PCB Mount Block",
	  material = Material("Wood", "Plywood"),
	  color = Color("purple"),
	  outer_contour = outer_contour,
	  inner_contours = [pcb_cutout],
	  start = P(0, y_center, z0),
	  end =   P(0, y_center, z1))

class GM3_Motor(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.y_center_l = L()
	self.y_direction_f = 1.0

    def configure(self, y_center = None, y_direction = None):
	none_type = type(None)
	assert type(y_center) == none_type or isinstance(y_center, L)
	assert type(y_direction) == none_type or isinstance(y_direction, float)
	if isinstance(y_center, L):
	    self.y_center_l = y_center
	if isinstance(y_direction, float):
	    self.y_direction_l = y_direction

    def construct(self):
	zero = L()

	# Grab some values from *self*:
	y_center = self.y_center_l
	y_direction = self.y_direction_f

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
	self.wheel_shaft_dy_l =   wheel_shaft_dy =   L(mm=9.20)
	self.encoder_shaft_dy_l = encoder_shaft_dy = L(mm=8.93)
	self.encoder_hub_dy_l =   encoder_hub_dy =   L(mm=2.00) # Guess
	self.pin_dy_l =		  pin_dy =           L(mm=1.00) # Guess
	self.dy_l =               dy =               L(mm=18.64)

	# Z dimensions:
	self.dz_l =                  dz =                  L(mm=22.23)
	self.back_holes_pitch_dz_l = back_holes_pitch_dz = L(mm=17.44)
	self.strap_pocket_dz_l = strap_pocket_dz =         L(mm=8.00)

	# Various Y position is ascending values:
	y0 = y_center - dy/2 - wheel_shaft_dy	# Wheel shaft edge
	y1 = y_center - dy/2 - pin_dy		# Alignment pin edge
	y2 = y_center - dy/2			# Wheel side of motor body
	y3 = y_center				# Motor center
	y4 = y3 + dy/2				# Encoder side of motor body
	y5 = y4 + encoder_hub_dy		# Encoder hub edge
	y6 = y4 + encoder_shaft_dy		# Edge of encoder shaft

	# Flip *y0* through *y5* if *y_direction* is negative:
	if y_direction < 0.0:
	    y0 = -y0
	    y1 = -y1
	    y2 = -y2
	    y3 = -y3
	    y4 = -y4
	    y5 = -y5
	    y6 = -y6

	# Contstuct the motor block:
	self.block(comment = "GM3 Motor Block",
	  material = Material("plastic", "ABS"),
	  color = Color("red"),
	  corner1 = P(back_dx,  y2, -dz/2),
	  corner2 = P(front_dx, y4,  dz/2),
	  top = "t")

	# Add various shafts, hubs, and pins:
	self.cylinder(comment = "wheel shaft",
	  diameter = wheel_shaft_diameter,
	  start = P(0, y0, 0),
	  end =   P(0, y3, 0))	# Weld
	self.cylinder(comment = "encoder shaft",
	  diameter = wheel_hub_diameter,
	  start = P(0, y6, 0),
	  end =   P(0, y3, 0))	# Weld
	self.cylinder(comment = "encoder hub",
	  diameter = encoder_hub_diameter,
	  start = P(0, y5, zero),
	  end =   P(0, y3, 0))	# Weld
	self.cylinder(comment = "align pin",
	  diameter = pin_diameter,
	  start = P(pin_dx, y1, 0),
	  end =   P(pin_dx, y2, 0)) # Weld

	# Drill some holes:
	self.hole(comment = "Mounting Hole 1",
	  diameter = back_holes_diameter,
	  start = P(back_holes_dx, y4,  back_holes_pitch_dz/2),
	  end =   P(back_holes_dx, y2,  back_holes_pitch_dz/2),
	  flags = "t")
	self.hole(comment = "Mounting Hole 2",
	  diameter = back_holes_diameter,
	  start = P(back_holes_dx, y4, -back_holes_pitch_dz/2),
	  end =   P(back_holes_dx, y2, -back_holes_pitch_dz/2),
	  flags = "t")
	self.hole(comment = "Shaft Center Screw Hole",
	  diameter = screw_hole_diameter,
	  start = P(0, y0, 0),
	  end =   P(0, y5, 0),
	  flags = "t")

class GM3_Wheel(Part):
    def __init__(self, up):
	Part.__init__(self, up)
	self.diameter_l = L(mm=69)
	self.width_l = L(mm=7.62)

    def construct(self):
	diameter = self.diameter_l
	width = self.width_l

	z = L()
	self.cylinder(comment = "Wheel",
	  material = Material("Plastic", "ABS"),
	  color = Color("blue"),
	  diameter = diameter,
	  start = P(z, -width/2, z),
	  end =   P(z,  width/2, z))

	self.hole(comment = "Screw Hole",
	  diameter = L(mm=8),
	  start = P(z,  width/2, z),
	  end =   P(z, -width/2, z),
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
	self.x_center_l = L()
	self.y_center_l = L()
	self.z_center_l = L()
	self.z_tongue_bottom = L(mm=5)
	self.z_tongue_top = L(mm=95)

    def configure(self, dx = None, dy = None, dz = None,
      inside_dy = None, x_center = None, y_center = None,
      z_tongue_bottom = None, z_tongue_top = None, z_center = None):
	""" *Front_Back_Motor_Side*: configure. """
	# Check argument types:
	none_type = type(None)
	assert type(dx) == none_type or isinstance(dx, L)
	assert type(dy) == none_type or isinstance(dy, L)
	assert type(dz) == none_type or isinstance(dz, L)
	assert type(inside_dy) == none_type or isinstance(inside_dy, L)
	assert type(x_center) == none_type or isinstance(x_center, L)
	assert type(y_center) == none_type or isinstance(y_center, L)
	assert type(z_center) == none_type or isinstance(z_center, L)
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
	if isinstance(x_center, L):
	    self.x_center_l = x_center
	if isinstance(y_center, L):
	    self.y_center_l = y_center
	if isinstance(z_center, L):
	    self.z_center_l = z_center
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
	x_center = self.x_center_l
	y_center = self.y_center_l
	z_center = self.z_center_l
	z_tongue_bottom = self.z_tongue_bottom_l
	z_tongue_top = self.z_tongue_top_l
	zero = L()

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
	outer_contour = Contour()
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
	self.y_center_l =          L(mm=0)
	self.side_tongue_dx_l =    L(mm=5)
	self.back_center_x_l =     L(mm=-20)
	self.front_center_x_l =    L(mm=20)
	self.side_tongue_bottom_z_l = L()
	self.side_tongue_top_z_l = L()

    def configure(self, bottom_z = None, dx = None, dy = None,
      side_tongue_dx = None, side_tongue_bottom_z = None,
      side_tongue_top_z = None, back_center_x = None, front_center_x = None,
      encoder_tongue_dz = None, top_z = None, y_center = None):
	""" Inner_Outer_Motor_Side: configure """
	# Check argument types:
	none_type = type(None)
	assert type(bottom_z) == none_type or isinstance(bottom_z, L)
	assert type(dx) == none_type or isinstance(dx, L)
	assert type(side_tongue_dx) == none_type or \
	  isinstance(side_tongue_dx, L)
	assert type(back_center_x) == none_type or \
	  isinstance(back_center_x, L)
	assert type(front_center_x) == none_type or \
	  isinstance(front_center_x, L)
	assert type(dy) == none_type or isinstance(dy, L)
	assert type(top_z) == none_type or isinstance(top_z, L)
	assert type(y_center) == none_type or isinstance(y_center, L)
	assert type(encoder_tongue_dz) == none_type or \
	  isisntance(encoder_tongue_dz, L)
    
	if isinstance(bottom_z, L):
	    self.bottom_z_l = bottom_z
	if isinstance(dx, L):
	    self.dx_l = dx
	if isinstance(dy, L):
	    self.dy_l = dy
	if isinstance(top_z, L):
	    self.top_z_l = top_z
	if isinstance(y_center, L):
	    self.y_center_l = y_center
	if isinstance(encoder_tongue_dz, L):
	    self.encoder_tongue_dz_l = encoder_tongue_dz
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
	top_z = self.top_z_l
	y_center = self.y_center_l

	#print("IOMSl.construct: name={0} x={1}:{2} y={3}:{4} z={5}:{6}".
	#  format(self, -dx/2, dx/2, y_center-dy/2, y_center+dy/2,
	#  bottom_z, top_z))

	# Grab some *Part*'s from the *motor_assembly*:
	motor_assembly = self.up
	gm3_motor = motor_assembly.gm3_motor_
	encoder_mount = motor_assembly.encoder_mount_
	dual_slot_encoder = motor_assembly.dual_slot_encoder_
	dual_slot_encoder_pcb = dual_slot_encoder.dual_slot_encoder_pcb_

	pcb_dx = dual_slot_encoder_pcb.dx_l
	encoder_tongue_dx = encoder_mount.tongue_dx_l
	x0 = -dx / 2
	x1 = back_center_x - side_tongue_dx/2
	x2 = back_center_x + side_tongue_dx/2
	x3 = -encoder_tongue_dx / 2
	x4 = -pcb_dx / 2
	x5 =  pcb_dx / 2
	x6 =  encoder_tongue_dx / 2
	x7 = front_center_x - side_tongue_dx/2
	x8 = front_center_x + side_tongue_dx/2
	x9 =  dx / 2

	y0 = y_center - dy / 2
	y1 = y_center + dy / 2
	print("IOMS:y0={0} y_center={1} y1={2} dy={3}".
	  format(y0, y_center, y1, dy))

	encoder_mount_z_center = encoder_mount.z_center_l
	encoder_mount_dz = encoder_mount.dz_l

	z_center = (bottom_z + top_z)/2
	z0 = bottom_z
	z1 = side_tongue_bottom_z
	z2 = z_center - encoder_tongue_dz/2
	z3 = encoder_mount_z_center - encoder_mount_dz/2 - L(mm=5)
	z4 = encoder_mount_z_center - encoder_mount_dz/2
	z5 = encoder_mount_z_center + encoder_mount_dz/2
	z6 = encoder_mount_z_center + encoder_mount_dz/2 + L(mm=5)
	z7 = z_center + encoder_tongue_dz/2
	z8 = side_tongue_top_z
	z9 = top_z
	#print("IOMS:z0={0} z1={1} z8={2} z9={3}".format(z0, z1, z8, z9))
	#print("IOMS:z3={0} z4={1} em_z_cent={2} z5={3} z6={4} em_dz={5}".
	#  format(z3, z4, encoder_mount_z_center, z5, z6, encoder_mount_dz))

	# Do the outer contour:
	zero = L()
	outer_contour = Contour()
	outer_contour.bend_append(P(x0, 0, z0), zero)
	outer_contour.bend_append(P(x0, 0, z9), zero)
	outer_contour.bend_append(P(x9, 0, z9), zero)
	outer_contour.bend_append(P(x9, 0, z0), zero)

	# Do the *encoder_pocket* for the encoder mount:
	encoder_pocket = Contour()
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
	back_pocket = Contour()
	back_pocket.bend_append(P(x1, 0, z1), zero)
	back_pocket.bend_append(P(x2, 0, z1), zero)
	back_pocket.bend_append(P(x2, 0, z8), zero)
	back_pocket.bend_append(P(x1, 0, z8), zero)
	#print("back_pocket={0}".format(back_pocket))

	# Do the *front_pocket*:
	front_pocket = Contour()
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
	  inner_contours = [encoder_pocket, back_pocket, front_pocket],
	  start = P(0, y0, 0),
	  end =   P(0, y1, 0))
	
	# Do the encoder hub hole:
	self.hole(comment = "Hub Hole 1",
	  diameter = gm3_motor.encoder_hub_diameter_l + L(mm=1),
	  start = P(0, y1, 0),
	  end =   P(0, y0, 0),
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
	  start = P(back_holes_dx, y1,  back_holes_pitch_dz/2),
	  end =   P(back_holes_dx, y0,  back_holes_pitch_dz/2),
	  flags = "t")
	self.hole(comment = "Mount Hole 2",
	  diameter = back_holes_diameter,
	  start = P(back_holes_dx, y1, -back_holes_pitch_dz/2),
	  end =   P(back_holes_dx, y0, -back_holes_pitch_dz/2),
	  flags = "t")
	self.hole(comment = "Alignment Pin",
          diameter = pin_diameter,
	  start = P(pin_dx, y1, 0),
	  end =   P(pin_dx, y0, 0),
	  flags = "t")
	self.simple_pocket(comment = "Strap Hole Pocket",
	  bottom_corner =
	  P(strap_hole_dx - strap_pocket_dx/2, y0, -strap_pocket_dz/2),
	  top_corner =
	  P(strap_hole_dx + strap_pocket_dx/2, y1, strap_pocket_dx/2),
	  pocket_top = "w")
	#self.hole(comment = "Strap Hole",
	#  diameter = strap_hole_diameter,
	#  start = P(strap_hole_dx, y1, 0),
	#  end = P(strap_hole_dx, y0, 0),
	#  flags = "t")

class Right_Motor_Assembly(Part):
    def __init__(self, up):
	""" Right_Motor_Assembly: initialize """
	self.thin_l = up.thin_l
	self.thick_l = up.thick_l
	Part.__init__(self, up)
	self.back_motor_side_ = Front_Back_Motor_Side(self)
	self.back_screw_ = Fastener(self)
	self.dual_slot_encoder_ = Dual_Slot_Encoder(self)
	self.encoder_disk_ = Encoder_Disk(self)
	self.encoder_mount_ = Encoder_Mount(self)
	self.front_motor_side_ = Front_Back_Motor_Side(self)
	self.front_screw_ = Fastener(self)
	self.gm3_motor_ = GM3_Motor(self)
	self.gm3_wheel_ = GM3_Wheel(self)
	self.inner_motor_side_ = Inner_Outer_Motor_Side(self)
	self.outer_motor_side_ = Inner_Outer_Motor_Side(self)
	self.inner_be_screw_ = Fastener(self)
	self.inner_bw_screw_ = Fastener(self)
	self.inner_te_screw_ = Fastener(self)
	self.inner_tw_screw_ = Fastener(self)
	self.outer_be_screw_ = Fastener(self)
	self.outer_bw_screw_ = Fastener(self)
	self.outer_te_screw_ = Fastener(self)
	self.outer_tw_screw_ = Fastener(self)
	self.top_back_screw_ = Fastener(self)
	self.top_front_screw_ = Fastener(self)
	self.encoder_se_screw_ = Fastener(self)
	self.encoder_sw_screw_ = Fastener(self)
	self.encoder_ne_screw_ = Fastener(self)
	self.encoder_nw_screw_ = Fastener(self)
	
    def construct(self):
	""" Right_Motor_Assembly: construct """

	# Grab the various *Part*' from *self*:
	back_motor_side = self.back_motor_side_
	back_screw = self.back_screw_
	dual_slot_encoder = self.dual_slot_encoder_
	dual_slot_encoder_pcb = dual_slot_encoder.dual_slot_encoder_pcb_
	encoder_disk = self.encoder_disk_
	encoder_mount = self.encoder_mount_
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
	slot_interrupter1 = dual_slot_encoder.slot_interrupter1_
	top_back_screw = self.top_back_screw_
	top_front_screw = self.top_front_screw_
	thick = self.thick_l
	thin = self.thin_l
	encoder_se_screw = self.encoder_se_screw_
	encoder_sw_screw = self.encoder_sw_screw_
	encoder_ne_screw = self.encoder_ne_screw_
	encoder_nw_screw = self.encoder_nw_screw_

	# Compute interesting X locations in ascending order:
	zero = L()
	dx = L(mm=130)
	x0 = -dx/2			# West edge of motor assembly
	x1 = x0 + L(mm=3)		# West edge of back motor side
	x2 = x1 + thick/2		# Center of back motor side
	x3 = x2 + thick/2		# East edge of back motor side
	x4 = x3 + L(mm=5)		# West tongue edge
	x5 = x4 + L(mm=2)		# End of screw mounts
	x6 = -dual_slot_encoder_pcb.dx/2 # West end of encoder pcb
	x7 = zero			# Center
        x8 = dual_slot_encoder_pcb.dx/2	# East end of encoder pcb
	x14 = dx/2			# East edge of motor assembly
	x13 = x14 - L(mm=3)		# East edge of front motor side
	x12 = x13 - thick/2		# Center of front motor side
	x11 = x12 - thick/2		# West edge of front motor side
	x10 = x11 - L(mm=5)		# East tongue edge
	x9 = x10 - L(mm=3)		# End of screw mounts

	# Compute interesting Y locations in ascending order:
	one = L(mm=1)
	y0 = zero					# Wheel center
	y1 = y0 + gm3_motor.wheel_shaft_dy_l + gm3_motor.dy_l/2 # Motor center
	y2 = y1 + gm3_motor.dy_l/2                      # Disk side Motor edge
	y3 = y0 + gm3_motor.wheel_shaft_dy_l - thin	# South inner side
	y4 = y3 + thin/2				# Center inner side
	y5 = y4 + thin/2				# North inner side
	#y6 done below:
	y7 = y3 + L(mm=10)				# End of mount screws
	y8 = y2 + gm3_motor.encoder_shaft_dy_l		# Encoder shaft edge
	y9 = y8 + encoder_disk.dy_l/2			# Encoder disk center
	y10 = y9 - dual_slot_encoder.disk_offset_dy_l	# PCB center

	y6 = y10 - dual_slot_encoder_pcb.dy/2		# East PCB edge
	y12 = y10 + dual_slot_encoder_pcb.dy/2		# West PCB edge

	y11 = y4 + L(mm=20)				# Center of assembly
	y16 = y11 + L(mm=20) + thin/2			# North outer side
	y15 = y16 - thin/2				# Outer side center
	y14 = y15 - thin/2				# South outer side
	y13 = y16 - L(mm=10)				# End mount screws

	print("TLA:y0={0} y1={1} y2={2} y3={3} y4={4} y5={5} y6={6} y7={7}".
	  format(y0, y1, y2, y3, y4, y5, y6, y7))
	print("TLA:y8={0} y9={1} y10={2} y11={3} y12={4} y13={5} y14={6}".
	  format(y8, y9, y10, y11, y12, y13, y14))
	print("TLA:y15={0} y16={1}".format(y15, y16))
	
	# Compute interesting Z locations in acending order:
	slot_dz = slot_interrupter1.dz_l.absolute()
	slot_gap_dz = slot_interrupter1.slot_dz_l.absolute()
	slot_base_dz = slot_dz - slot_gap_dz
	z0 = -gm3_motor.dz / 2 - L(mm=5)		# Outer side edge
	z1 = z0 + L(mm=10)				# Bottom tongue edge
	z2 = zero					# Shaft center-line
	z3 = z2 + encoder_disk.diameter_l / 2		# Top encoder disk edge
	z4 = z3 + L(mm=0.5)				# Top of disk gap
	z5 = z4 + slot_base_dz				# PCB bottom
	z6 = z5 + dual_slot_encoder_pcb.dz_l/2		# PCB middle
	z7 = z6 + encoder_mount.dz_l/2			# Encoder mount center
	z8 = z7 + encoder_mount.dz_l/2			# Encoder mount top
	z9 = z4 + L(mm=15)				# Top tongue edge
	z10 = z9 + L(mm=3)				# End of screw mount
	z11 = z10 + L(mm=7)				# Top edge of assembly
	z12 = z11 + L(mm=5)				# Top of mount base

	# Configure the various *Part*'s:
	gm3_motor.configure(y_center = y1)
	inner_motor_side.configure(dx = dx,
	  side_tongue_dx = x3 - x1, back_center_x = x2, front_center_x = x12,
	  y_center = y4, bottom_z = z0, top_z = z11,
	  side_tongue_bottom_z = z1, side_tongue_top_z = z9)
	outer_motor_side.configure(dx = dx,
	  side_tongue_dx = x3 - x1, back_center_x = x2, front_center_x = x12,
	  y_center = y15, bottom_z = z0, top_z = z11,
	  side_tongue_bottom_z = z1, side_tongue_top_z = z9)
	#print("RMA:z1={0} z9={1}".format(z1, z9))
	encoder_disk.configure(y_center = y9)
	dual_slot_encoder.configure(y_center = y10, z_center = z5)

	#print("RMA:x0={0} x11={1} x13={2} x14={3}".format(x0, x11, x13, x14))
	#print("RMA:z0={0} z11={1}".format(z0, z11))
	back_motor_side.configure(dx = x3 - x1, x_center = (x1 + x3)/2,
	  dy = y16 - y3, inside_dy = y14 - y5, y_center = (y16 + y3)/2,
	  dz = z11 - z0, z_tongue_bottom = z1, z_tongue_top = z9,
	  z_center = (z0 + z11)/2)
	front_motor_side.configure(dx = x13 - x11, x_center = (x13 + x11)/2,
	  dy = y16 - y3, inside_dy = y14 - y5, y_center = (y16 + y3)/2,
	  dz = z11 - z0, z_tongue_bottom = z1, z_tongue_top = z9,
	  z_center = (z0 + z11)/2)

	#inner_outer_dx = inner_motor_side.dx_l
	#inner_outer_pitch_dy = y15 - y4
	#inner_outer_dy = inner_motor_side.dy_l/2 + \
	#  inner_outer_pitch_dy +  outer_motor_side.dy_l/2
	#inner_outer_center_y = (y4 + y15)/2
	encoder_mount.configure(dx = x11 - x3, tongue_dx = x10 - x4,
	  dy = y16 - y3, y_center = y11, z_center = z7, dz = thick)

	# Configure the *inner_motor_side* fasteners:
	screw_flags = "M3x.05"
	screw_flags = "#4-40"
	inner_be_screw.configure(comment = "Inner BE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x12, y3, z0 + L(mm=5)),
	  end = P(x12, y7, z0 + L(mm=5)),
	  flags = screw_flags)
	#print("RMA:y2={0} y5={1}".format(y2, y5))
	inner_be_screw.drill(part = inner_motor_side, select = "close")
	inner_be_screw.nut_ledge(part = front_motor_side, flags="E")

	inner_bw_screw.configure(comment = "Inner BW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x2, y3, z0 + L(mm=5)),
	  end = P(x2, y7, z0 + L(mm=5)),
	  flags = screw_flags)
	inner_bw_screw.drill(part = inner_motor_side, select = "close")
	inner_bw_screw.nut_ledge(part = back_motor_side, flags="W")

	inner_te_screw.configure(comment = "Inner TE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x12, y3, z11 - L(mm=5)),
	  end = P(x12, y7, z11 - L(mm=5)),
	  flags = screw_flags)
	inner_te_screw.drill(part = inner_motor_side, select = "close")
	inner_te_screw.nut_ledge(part = front_motor_side, flags="E")

	inner_tw_screw.configure(comment = "Inner TW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x2, y3, z11 - L(mm=5)),
	  end = P(x2, y7, z11 - L(mm=5)),
	  flags = screw_flags)
	inner_tw_screw.drill(part = inner_motor_side, select = "close")
	inner_tw_screw.nut_ledge(part = back_motor_side, flags="W")

	# Configure the *outer_motor_side* fasteners:
	outer_be_screw.configure(comment = "Outer BE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x12, y16, z0 + L(mm=5)),
	  end = P(x12, y13, z0 + L(mm=5)),
	  flags = screw_flags)
	outer_be_screw.drill(part = outer_motor_side, select = "close")
	outer_be_screw.nut_ledge(part = front_motor_side, flags="E")

	outer_bw_screw.configure(comment = "Outer BW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x2, y16, z0 + L(mm=5)),
	  end = P(x2, y13, z0 + L(mm=5)),
	  flags = screw_flags)
	outer_bw_screw.drill(part = outer_motor_side, select = "close")
	outer_bw_screw.nut_ledge(part = back_motor_side, flags="W")

	outer_te_screw.configure(comment = "Outer TE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x12, y16, z11 - L(mm=5)),
	  end = P(x12, y13, z11 - L(mm=5)),
	  flags = screw_flags)
	outer_te_screw.drill(part = outer_motor_side, select = "close")
	outer_te_screw.nut_ledge(part = front_motor_side, flags="E")

	outer_tw_screw.configure(comment = "Outer TW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x2, y16, z11 - L(mm=5)),
	  end = P(x2, y13, z11 - L(mm=5)),
	  flags = screw_flags)
	outer_tw_screw.drill(part = outer_motor_side, select = "close")
	outer_tw_screw.nut_ledge(part = back_motor_side, flags="W")

	# Configure the *outer_motor_side* fasteners:
	back_screw.configure(comment = "Back Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x1, y11, z7),
	  end = P(x5, y11, z7),
	  flags = screw_flags)
	back_screw.drill(part = back_motor_side, select = "close")
	back_screw.nut_ledge(part = encoder_mount, flags="T")

	front_screw.configure(comment = "Front Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x13, y11, z7),
	  end = P(x9, y11, z7),
	  flags = screw_flags)
	front_screw.drill(part = front_motor_side, select = "close")
	front_screw.nut_ledge(part = encoder_mount, flags="T")

	# Configure the top screw fasteners:
	top_back_screw.configure(comment = "Top Back Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x2, y11, z12),
	  end = P(x2, y11, z10),
	  flags = screw_flags)
	top_back_screw.nut_ledge(part = back_motor_side, flags="W")

	top_front_screw.configure(comment = "Top Front Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x12, y11, z12),
	  end = P(x12, y11, z10),
	  flags = screw_flags)
	top_front_screw.nut_ledge(part = front_motor_side, flags="E")

	encoder_se_screw.configure(comment = "Encoder SE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x8 - L(mm=5), y6 + L(mm=5), z8),
	  end =   P(x8 - L(mm=5), y6 + L(mm=5), z5),
	  flags = screw_flags)
	encoder_se_screw.drill(part = encoder_mount, select = "close")

	encoder_sw_screw.configure(comment = "Encoder SW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x6 + L(mm=5), y6 + L(mm=5), z8),
	  end =   P(x6 + L(mm=5), y6 + L(mm=5), z5),
	  flags = screw_flags)
	encoder_sw_screw.drill(part = encoder_mount, select = "close")

	encoder_ne_screw.configure(comment = "Encoder NE Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x8 - L(mm=5), y12 - L(mm=5), z8),
	  end =   P(x8 - L(mm=5), y12 - L(mm=5), z5),
	  flags = screw_flags)
	encoder_ne_screw.drill(part = encoder_mount, select = "close")

	encoder_nw_screw.configure(comment = "Encoder NW Screw",
	  material = Material("Steel", "x"),
	  color = Color("black"),
	  start = P(x6 + L(mm=5), y12 - L(mm=5), z8),
	  end =   P(x6 + L(mm=5), y12 - L(mm=5), z5),
	  flags = screw_flags)
	encoder_nw_screw.drill(part = encoder_mount, select = "close")

	# Deal with visibility here:
	#gm3_wheel.invisible_set()
	#encoder_mount.invisible_set()
	#outer_motor_side.invisible_set()
	#inner_motor_side.invisible_set()
	#back_motor_side.invisible_set()
	front_motor_side.invisible_set()

	# Deal .dxf files here:
	encoder_mount.dxf_write(center = P(0, y11, z7))
	inner_motor_side.dxf_write(
	  center = P(0, y4, z7), plane_normal = P(0, L(mm=1), 0))
	front_motor_side.dxf_write(
	  center = P(x12, y11, z7), plane_normal = P(L(mm=1), 0, 0))

class ROSBot(Part):
    def __init__(self, up):
	self.thick_l = L(mm=5)
	self.thin_l = L(mm=3)
	Part.__init__(self, up)
	self.base_ = base = Base(self)
	self.right_motor_assembly_ = Right_Motor_Assembly(self)

    def construct(self):
	right_motor_assembly = self.right_motor_assembly_
	z = L()
	right_motor_assembly.place(translate = P(z, z, L(mm=50)))

class Slot_Interrupter(Part):
    def __init__(self, up):
	""" Slot_Interrupter: initialize """
	Part.__init__(self, up)
	self.dx_l = L(mm=2.6)
	self.dy_l = L(mm=4.5)
	self.dz_l = L(mm=-4.5)
	self.slot_dy_l = L(mm=2)
	self.slot_dz_l = L(mm=-2.1)
	self.x_center_l = L(mm=0.0)
	self.y_center_l = L(mm=0.0)
	self.z_center_l = L(mm=0.0)

    def configure(self, dx = None, dy = None, dz = None,
      slot_dy = None, slot_dz = None,
      x_center = None, y_center = None, z_center = None):
	""" Slot_Interrupter: configure """
	# Check argument types:
	none_type = type(None)
	assert type(dx) == none_type or isinstance(dx, L)
	assert type(dy) == none_type or isinstance(dy, L)
	assert type(dz) == none_type or isinstance(dz, L)
	assert type(slot_dy) == none_type or isinstance(slot_dy, L)
	assert type(slot_dz) == none_type or isinstance(slot_dz, L)
	assert type(x_center) == none_type or isinstance(x_center, L)
	assert type(y_center) == none_type or isinstance(y_center, L)
	assert type(z_center) == none_type or isinstance(z_center, L)

	# Note that *dy* and *slot_dz* can be negative:

	# Change any configuration values:
	if isinstance(dx, L):
	    self.dx_l = dx
	if isinstance(dx, L):
	    self.dy_l = dy
	if isinstance(dz, L):
	    self.dz_l = dz
	if isinstance(slot_dy, L):
	    self.slot_dy_l = slot_dy
	if isinstance(slot_dz, L):
	    self.slot_dz_l = slot_dz
	if isinstance(x_center, L):
	    self.x_center_l = x_center
	if isinstance(y_center, L):
	    self.y_center_l = y_center
	if isinstance(z_center, L):
	    self.z_center_l = z_center

    def construct(self):
	""" Slot_Interrupter: construct """
	zero = L()

	dx = self.dx_l
	dy = self.dy_l
	dz = self.dz_l
	slot_dy = self.slot_dy_l
	slot_dz = self.slot_dz_l

	x_center = self.x_center_l
	y_center = self.y_center_l
	z_center = self.z_center_l

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
	ezcad = EZCAD3(0, adjust = L(mm = -0.10))
	test = ROSBot(None)
	test.process(ezcad)
