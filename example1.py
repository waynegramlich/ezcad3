#!/usr/bin/env python

from EZCAD3 import *

class Camera(Part):

    def __init__(self, up):
	Part.__init__(self, up)
	self.board_ = Camera_Circuit_Board(self)
	self.sensor_ = Camera_Sensor(self)

    def construct(self):
	pass

class Camera_Circuit_Board(Part):

    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	self.length_l = length = L(mm=40.30)
	self.width_l = width = L(mm=25.00)
	self.height_l = height = L(mm=1.40)
	zero = L()

	sensor = self.up.sensor_

	self.block(
	  comment="pcb",
	  material=Material("plastic", "abs"),
	  color=Color("green", alpha = 0.5),
	  corner1=P(-length/2, -width/2, sensor.bz - height),
	  corner2=P( length/2,  width/2, sensor.bz))
	#  center=P(-length.half(), -width.half()), rotate=Angle.deg(90.0))

	print("length={0:m}".format(length))
	#print("Camera_Circuit_Board.box={0:m}\n".format(self.box))

class Camera_Sensor(Part):

    def __init__(self, up):
	Part.__init__(self, up)

    def construct(self):
	self.length_l = length = L(mm=8.90)
	self.width_l = width = L(mm=8.90)
	self.height_l = height =L(mm=0.15)
	zero = L()

	self.block(
	  comment="MT9V022M",
	  material=Material("silicon", ""),
	  color=Color("black", alpha = 0.5),
	  corner1=P(-length/2, -width/2, -height),
	  corner2=P( length/2,  width/2, zero))

	self.hole(comment="sample_hole",
	  diameter = L(mm=3.0), start=self.t, end=self.b)

def main():

    ezcad = EZCAD3(0)
    print("Instantiate camera")
    camera = Camera(None)
    print("Process camera")    
    ezcad.process(camera)
    #print("camera.box={0:m}".format(camera.box))

main()

