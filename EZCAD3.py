#!/usr/bin/python

## @package EZCAD3
#
# More details here...

import glob
import hashlib
import math
import numpy
import os
import os.path
import re
import subprocess
import xml.etree.ElementTree as ET

def int_compare(left, right):
    """ *int*: Return -1, 0, or 1 depending upon whether *left* sorts before,
	at, or after *right*.
    """

    assert isinstance(left, int)
    assert isinstance(right, int)
    if left < right:
	return -1
    elif left == right:
	return 0
    return 1

#  The L, P, and Angle classes are listed first because Python does not
# have a good mechanism for dealing with forward references to class names.

## @brief *L* corresponds to a length.  This class solves the issue
#  mixed units.  Internally, all lengths are turned into millimeters.
class L:
    """ A *L* represents a length. """

    ## @brief Initialize *self* (an *L* object) to the sum of
    #         *mm* + *cm* + *inch* + *ft*.
    #  @param self is the *L* object to initialize.
    #  @param mm is the number of millimeters to add in.
    #  @param cm is the number of centiimeters to add in.
    #  @param inch is the number of inches to add together in.
    #  @param ft is the number of feet to add together in.
    #
    # <I>__init__</I>() will initialize *self (an *L* object) to contain
    # the sum of *mm*, *cm*, *inch*, and *ft*.  Usually, only one of these
    # values is specfied.  The argment values must be either an *int* or
    # *float*.  The *inch* parameter can also be specified as a string
    # of the form "W-N/D" where is a whole number, N is a numerator, and
    # D is a denumerator.  In this case, the result is computed as
    # float(W)+float(N)/float(D).
    def __init__(self, mm=0.0, cm=0.0, inch=0.0, ft=0.0):
	""" *L*: Initialize *L* with sum of *mm* + *cm* + *inch* + *ft*. """

	# Check argment types:
	assert isinstance(mm, float) or isinstance(mm, int)
	assert isinstance(cm, float) or isinstance(cm, int)
	assert isinstance(inch, float) or isinstance(inch, int) or \
	       isinstance(inch, str)
	assert isinstance(ft, float) or isinstance(ft, int)

	# Deal with *inch* argument when it is a string:
	if isinstance(inch, str):
	    # We have a string, parse it:
            whole_fraction = inch.split("-")
	    if len(whole_fraction) == 2:
		whole = float(whole_fraction[0])
		fraction = whole_fraction[1]
	    elif len(whole_fraction) == 1:
		whole = 0.0
		fraction = inch
	    else:
		assert False, "Poorly formed inches string '{0}'".format(inches)
	    numerator_denominator = fraction.split("/")
	    assert len(numerator_denominator) == 2, \
	      "Bad fraction '{0}'".format(fraction)
	    numerator = float(numerator_denominator[0])
	    denominator = float(numerator_denominator[1])
	    inch = whole + numerator / denominator
		
	# Load up *self*:
	scale = 1.0
	self._mm = scale * (mm + cm * 10.0 + inch * 25.4 + ft * (12.0 * 25.4))

    ## @brief Returns *self* + *length*.
    #  @param self is the first *L* object to sum.
    #  @param length is the second *L* object to sum.
    #  @returns sum of *self* + *length*.
    #
    # <I>__add__</I>() returns the *self* + *length*.
    def __add__(self, length):
	""" *L*: Return *self* + *length*. """

	# Check argument types:
	assert isinstance(length, L)

	# Return result:
	return L(self._mm + length._mm)

    ## @brief Returns *self* / *number*.
    #  @param *self* is the *L* object to divide.
    #  @param *number* the value to divide into *self*
    #  @returns *self* / *number*
    #
    # <I>__div__</I>() returns *self* / *number*.  *number* must be either
    # be an *int*, *float*, or *L*.  For a divisor that is an *int* or *float*,
    # the returned result is of type *L*.  For a divisor that is an *L* type,
    # the returned quotient result is of type *float*.
    def __div__(self, divisor):
	""" *L*: Return {self} / {scalar}. """

	# Check argument types:
	assert isinstance(divisor, float) or \
	  isinstance(divisor, int) or isinstance(divisor, L)

	if isinstance(divisor, L):
	    result = self._mm / divisor._mm
	else:
	    result = L(mm = self._mm / divisor)

	return result

    ## @brief Returns *True* if *self* equals *length*.
    #  @param self is the first argument of the equality test.
    #  @param length is the second argument of the equality test.
    #  @returns *True* if *self* equals *length*.
    #
    # <I>__eq__</I>() returns *True* if *self* equals *length* and *False*
    # otherwise.
    def __eq__(self, length):
	""" *L*: Return *self* <= *length*. """

	# Check arguments types:
	assert isinstance(length, L)

	# Return result:
	return self._mm == length._mm

    ## @brief Returns *True* if *self* is greater than or equal to *length*.
    #  @param self is the first argument of the test.
    #  @param length is the second argument of the test.
    #  @returns *True* if *self* is greater than or equal to *length*.
    #
    # <I>__ge__</I>() returns *True* if *self* is greater than or equal
    # to *length* and *False* otherwise.
    def __ge__(self, length):
	""" *L*: Return *self* >= *length* """

	# Check argument types:
	assert isinstance(length, L)

	# Return result:
	return self._mm >= length._mm

    ## @brief Returns *True* if *self* is greater than to *length*.
    #  @param self is the first argument of the test.
    #  @param length is the second argument of the test.
    #  @returns *True* if *self* is greater than to *length*.
    #
    # <I>__gt__</I>() returns *True* if *self* is greater than *length*
    # and *False* otherwise.
    def __gt__(self, length):
	""" *L*: Return *self* > *length* """

	# Check argument types:
	assert isinstance(length, L)

	# Return result:
	return self._mm > length._mm

    ## @brief Returns *True* if *self* is less than or equal to *length*.
    #  @param self is the first argument of the test.
    #  @param length is the second argument of the test.
    #  @returns *True* if *self* is less than or equal to *length*.
    #
    # <I>__le__</I>() returns *True* if *self* is less than or equal
    # to *length* and *False* otherwise.
    def __le__(self, length):
	""" *L*: Return *self* <= *length*. """

	# Check argument types:
	assert isinstance(length, L)

	# Return result:
	return self._mm <= length._mm

    ## @brief Returns *True* if *self* is less than to *length*.
    #  @param self is the first argument of the test.
    #  @param length is the second argument of the test.
    #  @returns *True* if *self* is less than *length*.
    #
    # <I>__lt__</I>() returns *True* if *self* is less than *length*
    # and *False* otherwise.
    def __lt__(self, length):
	""" *L*: Return *self* < *length*. """

	# Check argument types:
	assert isinstance(length, L)

	# Return result:
	return self._mm < length._mm

    ## @brief Returns *self* as a *format*'ed string.
    #  @param self is the *L* object to format.
    #  @param format is the format control string.
    #  @returns formatted string.
    #
    # <I>__format__</I>() will return a formatted version of string.
    # The last letter of *format* specifies the unit formats -- "m"
    # is millimeters, "c" is centimeters, "i" is inches, and "f" is
    # is feet.  After the units suffix is stripped off, the remaining
    # suffix is treated as a *float* suffix (e.g. ".3" means 3 decimal
    # places after the decimal point.)
    def __format__(self, format):
	""" *L*: Return *self* as a *format*'ed string. """

	# Check argument types:
	assert isinstance(format, str)

	# Deal with units suffix:
	value = self._mm
	if format.endswith("c"):
	    # Centimeters:
	    value /= 10.0
	    format = format[:-1]
	elif format.endswith("f"):
	    # Feet:
            value /= (25.4 * 12)
	    format = format[:-1]
	elif format.endswith("i"):
            # Inches:
	    value /=  25.4
	    format = format[:-1]
	elif format.endswith("m"):
	    # We are already in millimeters:
	    format = format[:-1]

	# Now do the final format:
	if len(format) == 0:
	    result = "{0}".format(value)
	else:
	    result = ("{0:" + format + "f}").format(value)
	return result

    ## @brief Ruturn *self* &times; *number*.
    #  @param *self* is the first argument of the product.
    #  @param *number* is the second argument of the product.
    #  @returns *self* &times; *number*.
    #
    # <I>__mul__</I>() returns *self* &times; *number*.  *number* must
    # be either a *float* or an *int*.
    def __mul__(self, number):
	""" *L*: Multiply *self* by *number*. """

	# Check argument types:
	assert isinstance(number, float) or isinstance(number, int)

	# Return result:
	return L(self._mm * number)

    ## @brief Returns *True* if *self* is not equal to *length*.
    #  @param self is the first argument of the inequality test.
    #  @param length is the second argument of the inequality test.
    #  @returns *True* if *self* is not equal to *length*.
    #
    # <I>__eq__</I>() returns *True* if *self* is not equal to *length*
    # and *False* otherwise.
    def __ne__(self, length):
	""" L: Return {True} if {self} is not equal to {length}. """

	# Check argument types:
	assert isinstance(length, L)

	# Return result:
	return self._mm != length._mm

    ## @brief Return the negative of *self*.
    #  @param self is the *L* object to negate.
    #  @returns the negative of *self*.
    #
    # <I>__neg__</I>() returns the negative of *self*.
    def __neg__(self):
	""" *L*: Return -*self*. """

	# Return result:
	return L(-self._mm)

    ## @brief Return *self* &times; *number*.
    #  @param *self* is the first argument of the product.
    #  @param *number* is the second argument of the product.
    #  @returns *self* &times; *number*.
    #
    # <I>__mul__</I>() returns *self* &times; *number*.  *number* must
    # be either a *float* or an *int*.
    def __rmul__(self, number):
	""" *L*: Multiply *self* by *number*. """

	# Check argument types:
	assert isinstance(number, float) or isinstance(number, int)

	# Return result:
	return L(self._mm * number)

    ## @brief Return *self* converted to a string.
    #  @param self is the *L* object to convert.
    #  @returns *self* as a string.
    #
    # <I>__str__</I>() returns *self* converted to a string in units of
    # millimeters.
    def __str__(self):
	""" *L*: Return *self* as a formated string. """

	# Return result:
	return str(self._mm)

    ## @brief Returns *self* - *length*.
    #  @param self is the first *L* object to subtract from.
    #  @param length is the second *L* to subtract.
    #  @returns *self* - *length*.
    #
    # <I>__add__</I>() returns *self* - *length*.
    def __sub__(self, length):
	""" *L*: Return *self* - *length* . """

	# Check argument types:
	assert isinstance(length, L)

	# Return result:
	return L(self._mm - length._mm)

    ## @brief Returns the absolute value of *self*.
    #  @param self is the *L* object to convert to an absolut
    #  @returns |*self*|
    #
    # *absolute*() returns the absolute value of *self*.
    def absolute(self):
	""" *L*: Return |*self*|. """

	# Perform computation:
	result = self
	mm = self._mm
	if mm < 0.0:
	    result = L(-mm)
	return result

    ## @brief Returns the arc tangent of *self* / *dx* with resulting in the
    #         angle in the "correct" quadrant.
    #  @param *self* is the "Y" value of the arc tangent.
    #  @param *dx* is the "X" value fo the arc tangent.
    #  @returns arc tangent of *self* / *dx*.
    #
    # *arc_tangent2*() returns the arc tangent of *self* / *dx* with
    # resulting in the angle in the "correct" quadrant.  The result is
    # an *Angle*.
    def arc_tangent2(self, dx):
	""" L: Return the arctangent of {self} over {dx} being careful
	    that the returned angle is in the correct quadrant. """

	# Check argument type:
	assert isinstance(dx, L)

	# Return result:
	return Angle(rad = math.atan2(self._mm, dx._mm))

    def average(self, length):
	""" *L*: Return the average of the two *L* objects (i.e. *self* and *length)
	    as an *L* object. """
	return L((self._mm + length._mm) / 2)

    def centimeters(self):
	""" L: Return {self} as a scalar measured in centimeters. """

	return self._mm / 10.0

    def distance(self, dy):
	""" L: Return the length of the line between (*self*,*dy) and the origin.  """
	dx_mm = self._mm
	dy_mm = dy._mm
	distance_mm = math.sqrt(dx_mm * dx_mm + dy_mm * dy_mm)
	return L(mm=distance_mm)

    def cosine(self, angle):
	""" L: Return {self} * cos(angle). """

	assert isinstance(angle, Angle)
	return L(self._mm * angle.cosine())

    def inches(self):
	""" L: Return {self} as a scalar in units of inches. """

	return self._mm / 25.4

    def maximum(self, length):
	""" L: Return the maximum of {self} and {length}. """

	assert isinstance(self, L)
	result = self
	if length > result:
	    result = length
	return result

    def minimum(self, length):
	""" L: Return the minimum of {self} and {length}. """

	assert isinstance(self, L)
	result = self
	if length < result:
	    result = length
	return result

    def meters(self):
	""" L: Return {self} as a scalar in units of meters. """

	return self._mm / 1000.0

    def millimeters(self):
	""" L: Return a {self} as a scalar in units of millimeters. """

	return self._mm

    def sine(self, angle):
	""" L: Return {self} * sin(angle). """

	return L(self._mm * angle.sine())

class P:
    """ {P} represents a point in 3-space. """

    def __init__(self, x = 0, y = 0, z = 0):
	""" P: Intialize {self} to contain {part}, {x}, {y}, {z}. """

	# Deal with default arguments:
	assert isinstance(x, L) or isinstance(x, int) or isinstance(x, float)
	assert isinstance(y, L) or isinstance(y, int) or isinstance(y, float)
	assert isinstance(z, L) or isinstance(z, int) or isinstance(z, float)
	if not isinstance(x, L):
	    x = float(x)
	    assert x == 0.0, "Only 0 does not need to be an L() type"
	    x = L(mm = x)
	if not isinstance(y, L):
	    y = float(y)
	    assert y == 0.0, "Only 0 does not need to be an L() type"
	    y = L(mm = y)
	if not isinstance(z, L):
	    z = float(z)
	    assert z == 0.0, "Only 0 does not need to be an L() type"
	    z = L(mm = z)

	# Load up *self*:
	self.x = x
	self.y = y
	self.z = z

    def __add__(self, point):
	""" P: Add {point} to {self}. """

	assert isinstance(point, P)
	return P(self.x + point.x, self.y + point.y, self.z + point.z)

    def __div__(self, number):
	""" P: Return the result of dividing {self} by {number}. """

	number = float(number)
	assert isinstance(number, int) or isinstance(number, float)
	return P(self.x / number, self.y / number, self.z / number)

    def __eq__(self, point):
	""" P: Return {True} if {self} is equal to {point}. """

	assert isinstance(point, P)
	return self.x == point.x and self.y == point.y and self.z == point.z

    def __format__(self, format):
	""" *P*: Return *self* formatted as a string. """

	assert isinstance(format, str)
	if format != "":
            format = ":" + format
	format_string = \
	    "[{0" + format + "}, {1" + format + "}, {2" + format + "}]"

	return format_string.format(self.x, self.y, self.z)

    def __mul__(self, scalar):
	""" P: Return the result of muliplying {self} by {scalar}. """

	return P(self.x * scalar, self.y * scalar, self.z * scalar)

    def __rmul__(self, scalar):
	""" P: Return the result of muliplying {self} by {scalar}. """

	return P(self.x * scalar, self.y * scalar, self.z * scalar)

    def __ne__(self, point):
	""" P: Return {True} if {self} is not equal to {point}. """

	assert isinstance(point, P)
	return self.x != point.x or self.y != point.y or self.z != point.z

    def __neg__(self):
	""" P: Return the negative of {self}. """

	return P(-self.x, -self.y, -self.z)

    def __str__(self):
	""" P: Return {self} as a formatted string. """

	return "({0}, {1}, {2})".format(self.x, self.y, self.z)

    def __sub__(self, point):
	""" P: Subtract {point} from {self}. """

	assert isinstance(point, P)
	return P(self.x - point.x, self.y - point.y, self.z - point.z)

    def angle_between(self, point):
	""" P dimensions: Return the angle between {self} and {point}. """

	# Verify argument types:
	assert isinstance(point, P)

	# The math:
	#
	# a . b = ||a|| ||b|| cos( <AB )		(1)
	# (a . b) / (||a|| ||b||) = cos( <AB )		(2)
	# acos( (a . b) / (||a|| ||b||) ] = <AB		(3)
	# <AB = acos( (a . b) / (||a|| ||b||) ]		(4)

	# Grab the X/Y/Z coordinates in millimeters:
	x1 = self.x._mm
	y1 = self.y._mm
	z1 = self.z._mm
	x2 = point.x._mm
	y2 = point.y._mm
	z2 = point.z._mm

	# Compute the dot product, lengths and numerator:
	dot_product = x1 * x2 + y1 * y2 + z1 * z2
	length1 = math.sqrt(x1 * x1 + y1 * y1 + z1 * z1)
	length2 = math.sqrt(x2 * x2 + y2 * y2 + z2 * z2)
	numerator = length1 * length2

        result = Angle()
	#print("dot={0} num={1} xxx={2}".format(dot_product, numerator, xxx))
	if numerator > 0.0:
	    result = Angle(rad = math.acos(dot_product / numerator))
	else:
	    #assert False, "What happened here?"
	    pass
	return result

    def cross_product(self, point):
	""" *P*: Return the cross product of *self* with *point. """

	# Check argument types:
	assert isinstance(point, P)

	ux = self.x._mm
	uy = self.y._mm
	uz = self.z._mm
	vx = point.x._mm
	vy = point.y._mm
	vz = point.z._mm

	return P(L(mm = uy * vz - uz * vy),
	  L(mm = uz * vx - ux * vz), L(mm = ux * vy - uy * vx))

    def distance(self, point):
	assert isinstance(point, P)
	dx = (self.x - point.x)._mm
	dy = (self.y - point.y)._mm
	dz = (self.z - point.z)._mm
	length = L(mm = math.sqrt(dx * dx + dy * dy + dz * dz))
	return length

    def dot_product(self, point):
	""" *Point*: Return dot product of *self* . *point* . """
	assert isinstance(point, P)
	return L(mm = self.x._mm * point.x._mm +
	  self.y._mm * point.y._mm + self.z._mm * point.z._mm)

    def half(self):
	""" P dimensions: Return {self} / 2. """

	result = self / 2.0
	return result	

    def length(self):
	""" P dimensions: Return the length of self. """

	x = self.x._mm
	y = self.y._mm
	z = self.z._mm
	return L(mm=math.sqrt(x * x + y * y + z * z))

    def matrix_create(self):
	""" Matrix public: Create a matrix that corresponds to {self}. """

	result = Matrix(   [[
	  self.x._mm,
	  self.y._mm,
	  self.z._mm,
	  1.0              ]] )

	return result

    def normalize(self):
	""" *P*: Return *self* normalized to have a length of 1. """

	x = self.x._mm
	y = self.y._mm
	z = self.z._mm
	length = math.sqrt(x * x + y * y + z * z)
	return P(L(mm = x / length), L(mm = y / length), L(mm = z / length))
	

    def points(self, dx, dy, dz):
	""" Part construct: Return a list of points centered around
	    {self} that are separated by {dx} in X, {dy} in Y and {dx}
	    in Z.  If all of {dx}, {dy}, {dz} are non-zero, 8 {P}'s are
	    returned.  If one of {dx}, {dy}, and {dz} is zero, 4 {P}'s
	    are returned.  If two of {dx}, {dy}, and {dz} are zero, only
	    2 {P}'s are returned. """

	# Extract some values from {part}:
	part = self.part
	x = self.x
	y = self.y
	z = self.z

	# Construct {x_list}, {y_list}, {z_list} to have either 1 or 2 value
	# depending upon whether {dx}, {dy}, {dz} is zero, repsectively:
	zero = L.inch(0)
	if dx == zero:
	    x_list = (x)
	else:
            half_dx = dx / 2
	    x_list = (x - half_dx, x + half_dx)

	if dy == zero:
	    y_list = (y)
	else:
            half_dy = dy / 2
	    y_list = (y - half_dy, y + half_dy)

	if dz == zero:
	    z_list = (z)
	else:
            half_dz = dz / 2
	    z_list = (z - half_dz, z + half_dz)

	# Now iterate over {x_list}, {y_list}, and {z_list}
	# to generate {result}:
	result = []
	for x in x_list:
	    for y in y_list:
		for z in z_list:
		    point = P(part, x, y, z)
		    result.append(point)

	return result

    def twice(self):
	""" P dimensions: Return {self} * 2. """

	return self * 2.0

    def x_adjust(self, x):
	""" P dimensions: Return copy of {self} with {x} added to the
	    x field. """

	return P(self.part, self.x + x, self.y, self.z)

    def xy_adjust(self, x, y):
	""" P dimensions: Return copy of {self} with {x} added to the
	    x field and {y} added to the y field. """

	return P(self.part, self.x + x, self.y + y, self.z)

    def xyz_adjust(self, x, y, z):
	""" P dimensions: Return copy of {self} with {x} added to the
	    x field, {y} added to the y field, and {z} added to the z field. """

	return P(self.part, self.x + x, self.y + y, self.z + z)

    def xz_adjust(self, x, z):
	""" P dimensions: Return copy of {self} with {x} added to the
	    x field and {z} added to the z field. """

	return P(self.part, self.x + x, self.y, self.z + z)

    def y_adjust(self, y):
	""" P dimensions: Return copy of {self} with {y} added to the
	    y field. """

	return P(self.part, self.x, self.y + y, self.z)

    def yz_adjust(self, y, z):
	""" P dimensions: Return copy of {self} with {y} added to the
	    y field and {z} added to the z field. """

	return P(self.part, self.x, self.y + y, self.z + z)

    def z_adjust(self, z):
	""" P dimensions: Return copy of {self} with {z} added to the
	    z field. """

	return P(self.part, self.x, self.y, self.z + z)


## The *Angle* class represents an angle.
#
# An *Angle* object represents an angle.  *Angle*'s can be specified
# in either degrees or radians.
class Angle:
    """ An {Angle} represents an angle. """

    ## *PI* respresents the irrational number &pi;.
    PI = math.pi

    ## *TWO_PI* represents 2&pi;.
    TWO_PI = 2 * PI

    ## Initialize *Angle* to be *scalar_radians*.
    #  @param self *Angle* object to initialize.
    #  @param scalar_radians is the value of the angle in radians.
    #
    # <I>__init__</I>() will initialize *self* (an *Angle* object) to
    # *scalar_radians*.
    def __init__(self, deg = 0.0, rad = 0.0):
	""" *Angle*: Initialize self to contain *degrees* + *radians*. """

	# Check argument types:
	assert isinstance(deg, float) or isinstance(deg, int)
	assert isinstance(rad, float) or isinstance(rad, int)

	# Initliaze the final value:
	radians = float(rad)
	if deg != 0.0:
	    radians += float(deg) * math.pi / 180.0
	self.radians = radians
	#print("Angle.__init__(deg={0}, rad={1})=>{2}".
	#  format(deg, rad, radians))

    ## @brief Return the sum of two *Angle*'s.
    #  @param self is the first *Angle* object.
    #  @param angle is the second *Angle* object.
    #  @returns the sum of *self* and *angle* as a normalized *Angle* object.
    #
    # <I>__add__</I>() will add *self* to *angle* and return the result as an
    # *Angle* object.
    def __add__(self, angle):
	""" Angle: Return *angle* added to *self*. """

	assert isinstance(angle, Angle)
	return Angle(rad = self.radians + angle.radians)

    ## @brief Divides an *Angle* by *divisor*.
    #  @param self is the *Angle* object to divide.
    #  @param divisor the amount to divide the *Angle* bye.
    #  @returns an *Angle* the corresponds to *self* / *divisor.
    #
    # <I>__div__</I>() will divide *self* by *divisor*.  *divisor* must be
    # either a *float* or an *int*.
    def __div__(self, divisor):
	""" Angle: Return *self* divided by *divisor*. """

	result = None
	if isinstance(divisor, float) or isinstance(divisor, int):
	    result = Angle(rad = self.radians / float(divisor))
	elif isinstance(divisor, Angle):
	    result = self.radians/divisor.radians
	else:
	    assert False, "divisor is neither an Angle nor an integer/float"
	return result

    ## @brief Return *True* if angles are equal.
    #  @param self is the first *Angle* to compare.
    #  @param angle is the second *Angle* to compare.
    #  @returns *True* if *self* equals *angle*.
    #
    # <I>__eq__</I>() will return *True* if *self* equals *angle* and *False*
    # otherwise.  Due to the fact that floating point arithmetic is
    # involved, to *Angle*'s could be very close to one another and
    # still not exactly match.
    def __eq__(self, angle):
	""" Angle: Return *True* if *self* is equal to *angle*. """

	assert isinstance(angle, Angle)
	return self.radians == angle.radians

    ## @brief Return a formatted version of *self* using *format* to
    #         control formatting.
    #  @param self is the *Angle* to format.
    #  @param format is the the formatting control string.
    #  @returns a formatted version of *self* as a *str*.
    #
    # <I>__format__</I>() will return a formatted version of *self* using
    # *format* to control formatting.  Right now *format* is ignoreed and
    # the value returned is *self* converted to degrees.
    def __format__(self, format):
	""" Angle: Format *self* into a string and return it. """

	# Check argument types:
	assert isinstance(format, str)

	#print("Angle.__format__(*, '{0}')".format(format))

	value = self.radians
	if format.endswith("r"):
	    # Do radians:
	    format = format[:-1]
	elif format.endswith("d"):
	    # Do degrees:
            format = format[:-1]
            value = value * 180.0 / math.pi
        else:
	    # Assume degrees output by default:
            value *= 180.0 / math.pi

	if len(format) == 0:
	    result = "{0}".format(value)
	else:
	    result = ("{0:" + format + "f}").format(value)
	#print("Angle.__format__()=>'{0}'".format(result))
	return result

    ## @brief Return *True* if *self* is greater than or equal to *angle*.
    #  @param self is the first *Angle* to compare.
    #  @param angle is the second *Angle* to compare.
    #  @returns *True* if *self* is greater than or equal to *angle*.
    #
    # <I>__ge__</I>() will return *True* if *self* is greater than or 
    # equal to *angle* and *False* otherwise.
    def __ge__(self, angle):
	""" Angle: Return *True* if *self* is greater than or equal
	    to *angle*. """

	assert isinstance(angle, Angle)
	return self.radians >= angle.radians

    ## @brief Return *True* if *self* is greater than *angle*.
    #  @param self is the first *Angle* to compare.
    #  @param angle is the second *Angle* to compare.
    #  @returns *True* if *self* is greater than *angle*.
    #
    # <I>__gt__</I>() will return *True* if *self* is greater than *angle*
    # and *False* otherwise.
    def __gt__(self, angle):
	""" Angle: Return *True* if *self* is greater than *angle*. """

	assert isinstance(angle, Angle)
	return self.radians > angle.radians

    ## @brief Return *True* if *self* is less than *angle*.
    #  @param self is the first *Angle* to compare.
    #  @param angle is the second *Angle* to compare.
    #  @returns *True* if *self* is less than *angle*.
    #
    # <I>__gt__</I>() will return *True* if *self* is less than *angle*
    # and *False* otherwise.
    def __le__(self, angle):
	""" Angle: Return *True* if *self* is less than or equal to *angle*. """

	assert isinstance(angle, Angle)
	return self.radians <= angle.radians

    ## @brief Return *True* if *self* is less than or equal to *angle*.
    #  @param self is the first *Angle* to compare.
    #  @param angle is the second *Angle* to compare.
    #  @returns *True* if *self* is less than or equal to *angle*.
    #
    # <I>__ge__</I>() will return *True* if *self* is less than or 
    # equal to *angle* and *False* otherwise.
    def __lt__(self, angle):
	""" Angle: Return *True* if *self* is less than *angle*. """

	assert isinstance(angle, Angle)
	return self.radians < angle.radians

    ## @brief Return *self* &times; *multiplier*.
    #  @param self is the *Angle* object to multiply.
    #  @param multiplier is the number to mulitply *self* by.
    #  @returns *self* &times *multiplier*.
    #
    # <I>__mul__</I>() returns *self* (an *Angle* object) &times; *multiplier*.
    def __mul__(self, multiplier):
	""" Angle: Return *self* divided by *multiplier*. """

	assert isinstance(multiplier, float) or isinstance(multiplier, int)
	return Angle(rad = self.radians * multiplier)

    ## @brief Return *True* if angles are equal.
    #  @param self is the first *Angle* to compare.
    #  @param angle is the second *Angle* to compare.
    #  @returns *True* if *self* equals *angle*.
    #
    # <I>__ne__</I>() will return *True* if *self* is not equal to *angle*
    # and *False* otherwise.
    def __ne__(self, angle):
	""" Angle: Return {True} if {self} is not equal to {angle}. """

	assert isinstance(angle, Angle)
	return self.radians != angle.radians

    ## @brief Returns -*self*.
    #  @param self is the *Angle* object to negate.
    #
    # <I>__neg__</I>() will return -*self*.
    def __neg__(self):
	""" Angle: Return negative of {self}. """

	return Angle(rad = -self.radians)

    ## @brief Return *self* &times; *multiplier*.
    #  @param self is the *Angle* object to multiply.
    #  @param multiplier is the number to mulitply *self* by.
    #  @returns *self* &times *multiplier*.
    #
    # <I>__rmul__</I>() returns *self* (an *Angle* object) &times; *multiplier*.
    def __rmul__(self, multiplier):
	""" Angle: Return {self} multiplied by {multiplier}. """

	assert isinstance(multiplier, float) or isinstance(multiplier, int)
	return Angle(rad = self.radians * multiplier)

    ## @brief Returns *self* converted to degrees.
    #  @param self is the *Angle* object to convert to degrees.
    #  @returns *self* converted to degrees.
    #
    # <I>__str__</I>() Returns *self* converted to degrees as a *str* object.
    def __str__(self):
	""" Angle: Return {self} as a string. """

	return str(self.radians * 180.0 / Angle.pi)

    ## @brief Returns *self* - *angle*.
    #  @param self is the *Angle* object to subtract from.
    #  @param angle is the *Angle* to subtract.
    #  @returns *self* - *angle*.
    #
    # <I>__sub__</I>() returns *self* - *angle*.
    def __sub__(self, angle):
	""" Angle: Return {angle} subtracted from {self}. """

	assert isinstance(angle, Angle)
	return Angle(rad = self.radians - angle.radians)

    def absolute(self):
	""" *Angle*: Return absolution value of *Angle* object(i.e. *self*).
	"""

	result = self
	if self.radians < 0.0:
	    result = -self
	return result

    ## @brief Returns the cosine of *self*.
    #  @param self is the *Angle* to compute the cosine of
    #  @returns the cosine of *self*.
    #
    #  *cosine*() returns the cosine of *self* (an *Angle* object.)
    def cosine(self):
	""" Angle: Return the cosine of {self}. """

	return math.cos(self.radians)

    ## @brief Returns *scalar_degrees* as an *Angle* object.
    #  @param scalar_degrees is the number to convert into an *Angle*.
    #  @returns *scalar_degrees* as an *Angle*.
    #
    # *deg*() will return *scalar_degrees* as an *Angle* object.
    @staticmethod
    def deg(scalar_degrees):
	""" Angle: Convert {scalar_degrees} into an {Angle} and return it. """

	assert isinstance(scalar_degrees, float) or \
	  isinstance(scalar_degrees, int)
	return Angle(rad = scalar_degrees * math.pi / 180.0)

    ## @brief Returns *self* converted degrees as a *float*.
    #  @param self is the *Angle* to convert to degrees.
    #  @returns *self* converted to degrees.
    #
    # *degrees*() will return *self* converted to degrees
    def degrees(self):
	""" Angle: Convert {self} back into degrees and return it. """

	return self.radians * 180.0 / math.pi

    ## @brief Returns *self*/2.
    #  @param self is the *Angle* object to divide by 2.
    #  @returns *self/2.
    #
    # *half*() returns the *self*/2.
    def half(self):
	""" Angle: Return half of {self}. """

	return Angle(self.radians / 2.0)

    ## @brief Return a normalized value between -&pi; and &pi;.
    #  @param self is the *Angle* to normalize.
    #  @returns a normalized value for *self*.
    #
    # *normalize*() will return a normalized angle between -&pi; and &pi;.
    def normalize(self):
	""" Angle: Return a normalized angle. """

	scalar_radians = self.scalar_radians
	while (scalar_radians > PI):
            scalar_radians -= TWO_PI
	while (scalar_radians < PI):
	    scalar_radians += TWO_PI
	return Angle(scalar_radaians)

    ## @brief Returns *scalar_radians* as an *Angle* object.
    #  @param scalar_radians is the number to convert into an *Angle*.
    #  @returns *scalar_radians* as an *Angle*.
    #
    # *deg*() will return *scalar_radians* as an *Angle* object.
    @staticmethod
    def rad(scalar_radians):
	""" Angle: Convert *scalar_radians* radians into an *Angle* and
	    return it. """

	return Angle(scalar_radians)

    ## @brief Returns *self* converted to radians.
    #  @param self is the *Angle* object to convert to radians.
    #  @returns *self* converted to radians.
    #
    # *radians*() returns *self* converted to radians.
    def radians(self):
	""" Angle: Convert *self* into radians and return it. """

	return self.radians

    ## @brief Returns the sine of *self*.
    #  @param self is the *Angle* to compute the sine of
    #  @returns the sine of *self*.
    #
    #  *sine*() returns the sine of *self* (an *Angle* object.)
    def sine(self):
	""" Angle: Return the sine of {self}. """

	return math.sin(self.radians)

    ## @brief Returns the tangent of *self*.
    #  @param self is the *Angle* to compute the tangent of
    #  @returns the tangent of *self*.
    #
    #  *tangent*() returns the tangent of *self* (an *Angle* object.)
    def tangent(self):
	""" Angle: Return the tangent of {self}. """

	return math.tan(self.radians)

    ## @brief Returns *self* times 2.
    #  @param self is the *Angle* to double.
    #  @returns twice the value of *self*.
    #
    # *twice*() returns twice the value of *self*.
    def twice(self):
	""" Angle: Return twice of {self}. """

	return Angle(self.length * 2.0)

class Bend:

    def __init__(self, point, radius, plane = None, name = ""):
	assert isinstance(point, P)
	assert isinstance(radius, L)

	if not isinstance(plane, L):
	    plane = L()
	self.point = point
	self.point_tmp = point
	self.radius = radius
	self._px = 0.0
	self._py = 0.0
	self._name = name
	self._plane = plane

    def __format__(self, format):
	return "Bend(point={0}, radius={1} name='{2}'". \
	  format(self.point, self.radius, self._name)

	#return ("P={0} radius={1} P=({2},{3}) Cx=({4},{5})" + \
	#  " Before=({6},{7}) After=({8},{9}) bf={10:.2f} af={11:.2f}"). \
	#  format(self.point, self.radius,
	#  self._px, self._py, self._center_x, self._center_y,
	#  self._before_tangent_x, self._before_tangent_y,
	#  self._after_tangent_x,  self._after_tangent_y,
	#  self._before_fraction, self._after_fraction)

    def compute(self, before_bend, after_bend, trace = -1000000):
	""" *Bend*: compute bend information. """
	assert isinstance(before_bend, Bend)
	assert isinstance(after_bend, Bend)
	assert isinstance(trace, int)

	if trace >= 0:
	    print("{0}=>Bend.compute({1})".format(trace * ' ', self))

	#print("Bend.compute(): before_bend={0}".format(before_bend))
	#print("Bend.compute(): after_bend={0}".format(after_bend))

        # Below is an ASCII art picture of a *Bend*.  B represents the
	# *Bend* point (i.e. *self*.*point*).  *D* and *E* are the
	# adjecent *Bend* objects and their associated *point* objects.
	# The ultimate goal here is to compute to compute C, which is the
	# center of a circle of radius R (i.e. *self*.*radius*) that
	# touches to line segments BD and BE tangentially.
	#
	#           C
	#           |
	#    D      |      E
	#     \     |     /
	#      \    |    /
	#       \   |   /
	#        \  |  /
	#         \ | /
	#          \|/
	#           B
	#
	#
	# The center of the circle must be on segment BC which is
	# the line segement that bisects the angle <DBE.  Thus, the
	# angles <DBC and <CBE are equal.  We will compute <CBE:
	#
	#        BD . BE = |BD| * |BE| * cos(<DBE)                (1)
	#
	#        cos(<DBE) = (BD . BE) / (|BD| * |BE|)            (2)
	#
	#        <DBE = acos( (BD . BE) / (|BD| * |BE|) )         (3)
        #
	#        <CBE = <DBE / 2                                  (4)
	#
	#        <CBE = acos( (BD . BE) / (|BD| * |BE|) ) / 2     (5)
	#
	#
	# The crude picture above is redrawn with the line BE horizontal
	# and segment BC going up at an angle.  Even though the ASCII
	# art does not show it, angle <CBE is the same for both the
	# diagram above and below.  A circle of radius R with touch
	# to segment BE at location T.  Thus, angle <BTC is 90 degrees.
	# The length of segment |CT| is R.  The length of segment |BC|
	# is defined as L.
	#
	#          C
	#         /|
	#        / |
	#     L /  |R
	#      /   |
	#     /    |
	#    B-----T-----E
	#
	# Using trigonometery:
	#
	#        R = L * sin(<CBE)                            (1)
	#
	# Solving for L:
	#
	#        L = R / sin(<CBE)                            (2)
	#
	# Now we do this:

	# Extract the radius:
	radius = self.radius._mm

	# Extract the X/Y coordinates for B (i.e. *self*):
	b = self
	bx = b._px
	by = b._py

	if trace >= 0:
	    print("{0}Bend.compute(): bx={1} by={2}".
	      format(trace * ' ', bx, by))

	# Extract the X/Y coordinates for D (i.e. *bend1*.*\_point*):
	d = before_bend
	dx = d._px
	dy = d._py

	# Extract the X/Y coordinates for E (i.e. *bend2*.*\_point*):
	e = after_bend
	ex = e._px
	ey = e._py

	#print("before:{0} at:{1} after:{2}".format(d.point, b.point, e.point))

	# Compute the direction vector DB:
	dbx = dx - bx
	dby = dy - by

	# Compute the direction vector EB:
	ebx = ex - bx
	eby = ey - by

	#print("Point:{0}: db=({1},{2}) eb=({3},{4})".
	#  format(self.point, dbx, dby, ebx, eby))

	# Compute the length of |DB|:
	db_length = math.sqrt(dbx * dbx + dby * dby)

	# Compute the length of |EB|:
	eb_length = math.sqrt(ebx * ebx + eby * eby)

	#print("Point:{0}: db_length=({1}) eb_length=({2})".
	# format(self.point, db_length, eb_length))

	# Compute normalized DB direction vector -- nDB:
	ndbx = dbx / db_length
	ndby = dby / db_length

	# Compute normalized EB direction vector -- nEB:
	nebx = ebx / eb_length
	neby = eby / eb_length

	#print("Point:{0}: ndb=({1},{2}) neb=({3},{4})".
	# format(self.point, ndbx, ndby, nebx, neby))

	# Compute the dot product of nDB . nEB:
	dot_product = ndbx * nebx + ndby * neby
	
	# Compute <DBE:
	angle_dbe = math.acos(dot_product)

	# Compute <CBE:
	angle_cbe = angle_dbe / 2.0

	r2d = 180.0 / Angle.PI
	#print("<dbe={0} <cbe={1}".format(angle_dbe * r2d, angle_cbe * r2d))

	# Now compute length:
	bc_length = radius / math.sin(angle_cbe)

	#print ("bc_length={0}".format(bc_length))

	# We need a direction vector.  Compute F which is at the midpoint
	# between *ndb* and *neb*:
	fx = (ndbx + nebx) / 2.0
	fy = (ndby + neby) / 2.0

	# Compute the length of |FB|:
	f_length = math.sqrt(fx * fx + fy * fy)

	#print("f=({0},{1}) f_length={2}".format(fx, fy, f_length))

	# Compute the normalized vector nFB:
	nfx = fx / f_length
	nfy = fy / f_length

	#print("nfb=({0},{1})".format(nfx, nfy))

	# Compute C:
	cx = bx + nfx * bc_length
	cy = by + nfy * bc_length

	#print("c=({0},{1})".format(cx, cy))

	# Load the center into *self*:
	self._center_x = cx
	self._center_y = cy

	# Now we compute the tangent points::
	tangent_length = bc_length * math.cos(angle_cbe)
	self._before_tangent_x = bx + tangent_length * ndbx
	self._before_tangent_y = by + tangent_length * ndby
	self._after_tangent_x = bx + tangent_length * nebx
	self._after_tangent_y = by + tangent_length * neby

	# The fraction is how much of the segment length is used
	# for bend radius.  The two fractions for a given segment must
	# sum to be less than or equal to 1.0:
	self._before_fraction = tangent_length / db_length
	self._after_fraction = tangent_length / eb_length

	if trace >= 0:
	    print("{0}<=Bend.compute({1})".format(trace * ' ', self))


    def copy(self):
	""" *Bend* """
	result = Bend(self.point, self.radius)
	result._px = self._px
	result._py = self._py
	return result

class XBlock:

    def __init__(self, part, tool, vice_x, vice_y):
	""" *Block*: Initialize the *Block* object (i.e. *self*) with *part*,
            *tool*, *vice_x*, and *vice_y*. """

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(tool, Tool)
	assert isinstance(vice_x, L)
	assert isinstance(vice_y, L)

	# Generate a UID:
	shop = part._shop_get()
	uid = shop._blocks_uid_get()

	# Initialzie the *Block* object (i.e. *self*)
	self._feed = Speed(in_per_sec=0.0)
	self._operations = []	# list[Code_Operation]
	self._part = part
	self._priority = part._priority_get()
	self._program_number = 0
	self._spindle = Hertz()
	self._tool = tool
	self._text = ""
	self._uid = uid
	self._vice_x = vice_x
	self._vice_y = vice_y

    def _comment_get(self):
	""" *Block*: Get the comment field of the *Block* object (i.e. *self*) to *comment*. """

	return self._comment

    def _comment_set(self, comment):
	""" *Block*: Set the comment field of the *Block* object (i.e. *self*) to *comment*. """

	# Verify argument types:
	assert isinstance(comment, str)
	
	# Load *comment* into the *Block* object (i.e. *self*):
	self._comment = comment

    def _program_number_get(self):
	""" *Block*: Return the program number field of the *Block* object (i.e. *self*.) """

	return self._program_number

    def _part_get(self):
	""" *Block*: Return the part field of the *Block* object (i.e. *self*.) """

	return self._part

    def _spindle_get(self):
	""" *Block*: Return the spindle field of the *Block* object (i.e. *self*.) """

	return self._spindle

    def _spindle_set(self, spindle):
	""" *Block*: Set the spindle speed field of the *Block* object (i.e. *self*) to *spindle*.
	"""

	# Verify argument types:
	assert isinstance(spindle, Hertz)
	
	# Load *spindle* into the *Block* object (i.e. *self*):
	self._spindle = spindle

    def _text_get(self):
	""" *Block*: Return the text field of the *Block* object (i.e. *self*.) """

	
	# Load *text* into the *Block* object (i.e. *self*):
	return self._text

    def _text_set(self, text):
	""" *Block*: Set the text field of the *Block* object (i.e. *self*) to *text*. """

	# Verify argument types:
	assert isinstance(text, str)
	
	# Load *text* into the *Block* object (i.e. *self*):
	self._text = text

    def _tool_get(self):
	""" *Block*: Return the tool field of the *Block* object (i.e. *self*.) """
	
	return self._tool

    def _uid_get(self):
	""" *Block*: Return the unique identifier field of the *Block* object (i.e. *self*.) """
	
	return self._uid

    def _vice_x_get(self):
	""" *Block*: Return the vice X field of the *Block* object (i.e. *self*.) """
	
	return self._vice_x

    def _vice_y_get(self):
	""" *Block*: Return the vice Y field of the *Block* object (i.e. *self*.) """
	
	return self._vice_y

    def drill_append(self, diameter, f, x, y, z_start, z_stop):
	""" *Block*: ... """

	# Verify argument types:
	assert isinstance(diameter, L)
	assert isinstance(f, L)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(z_start, L)
	assert isinstance(z_stop, L)

	operation_drill = Code_Drill(diameter, f, x, y, z_start, z_stop)
	block._operations.append(operation_drill)

class Hertz:
    def __init__(self, frequency=0.0, rpm=0.0, rps=0.0):
	""" *Hertz*: Initialize the *Hertz* object (i.e. *self*) to be the sum of
	    *frequency*, *rpm* (Rotations Per Minute) and *rps* (Rotations Per Second). """

	# Verify argument types:
	assert isinstance(frequency, float)
	assert isinstance(rpm, float)
	assert isinstance(rps, float)

	# Load up self:
	self._frequency = frequency + rpm / 60.0 + rps

    def __format__(self, format):
	""" *Hertz*: Return the *Hertz* object (i.e. *self*) formatted by *format* which
	    must be either "rps" or "rpm". """

	frequency = self._frequency
	result = ""
	if format == "" or format == "frequency" or format == "rps":
	    result = "{0}".format(frequency)
	if format == "rpm":
	    result = "{0}".format(60.0 * frequency)
	else:
	    result = "Unknown Hertz format '{0}'".format(format)
        return result

    def __eq__(hertz1, hertz2):
	""" *Hertz*: Return true if the *hertz1* is equal to *hertz2. """

	assert isinstance(hertz2, Hertz)
	return hertz1._frequency > hertz2._frequency

    def __gt__(hertz1, hertz2):
	""" *Hertz*: Return true if the *hertz1* is greater than *hertz2. """

	assert isinstance(hertz2, Hertz)
	return hertz1._frequency > hertz2._frequency

    def __ge__(hertz1, hertz2):
	""" *Hertz*: Return true if the *hertz1* is greater than or equal to *hertz2. """

	assert isinstance(hertz2, Hertz)
	return hertz1._frequency >= hertz2._frequency

    def __lt__(hertz1, hertz2):
	""" *Hertz*: Return true if the *hertz1* is less than *hertz2. """

	assert isinstance(hertz2, Hertz)
	return hertz1._frequency < hertz2._frequency

    def __le__(hertz1, hertz2):
	""" *Hertz*: Return true if the *hertz1* is less than or equal to *hertz2. """

	assert isinstance(hertz2, Hertz)
	return hertz1._frequency <= hertz2._frequency

    def __ne__(hertz1, hertz2):
	""" *Hertz*: Return true if the *hertz1* is not equal to *hertz2. """

	assert isinstance(hertz2, Hertz)
	return hertz1._frequency != hertz2._frequency

    def frequency(self):
	""" *Hertz*: Return the frequence of the * """
	return self._frequency

class Material:
    def __init__(self, generic = "plastic", specific = "ABS"):
	""" *Material*: Initialize *Material* object (i.e. *self*)
	    to contain *generic* and *specific*.
	"""

	# Verify argument types:
	assert isinstance(generic, str)
	assert isinstance(specific, str)

	# Load up *self*
	self._generic = generic.lower()
	self._specific = specific.lower()

    def __format__(self, format):
	""" *Material*: Return *Maerial* object (i.e. *self*) formatted
	    as a string.
	"""

	return "[{0}, {1}]".format(self._generic, self._specific)

    def _find(self, material_name):
	""" *Material*: Return True if *material_name* matches *self*.
	"""

	# Verify argument types:
	assert isinstance(material_name, str)
	return self._generic == material_name

    def _needs_coolant(self):
	""" *Material*: Return *True* if the *Material* object (i.e. *self*) should
	    be machined with coolant on and *False* otherwise.
	"""

	return self._generic != "plastic"

    def _generic_get(self):
	""" *Material*: Return generic name for *Material* object (i.e. *self*).
	"""

	return self._generic

class Bounding_Box:

    def __init__(self, points = []):
	""" *Bounding_Box*: Initialize the *Bounding_Box* object (i.e. *self*)
	    with each point in *points*.
	"""

	# Verify argument types:
	assert isinstance(points, list) or isinstance(points, tuple)
	for point in points:
	    assert isinstance(point, P)

	# Intialize the *Bounding_Box* object (i.e. *self*):
	zero = L()
	self._is_empty = True	# *True* if no points have been added in yet.
	self._east = zero	# Highest X coordinate
	self._west = zero	# Lowest X coordinate
	self._north = zero	# Highest Y coordinate
	self._south = zero	# Lowest Y cooridinate
	self._top = zero	# Highest Z coordinate
	self._bottom = zero	# Lowest Z coordinate

	# Merge each *point* in *points* into the *Bounding_Box* object (i.e. *self*):
	for point in points:
	    self.point_expand(point)

    def __eq__(self, bounding_box2):
	""" *Bounding_Box*: Return *True* if the *Bounding_Box* object
	    (i.e. *self*) matches *bounding_box2* and *False* otherwise.
	"""
	
	# Verify argument types:
	assert isinstance(bounding_box2, Bounding_Box)

	# Use *bounding_box1* instead of *self*.
	bounding_box1 = self

	# Test for equality:
	result = False
	if bounding_box1._is_empty:
	    # *bounding_box1* is empty:
	    if bounding_box2._is_empty:
		# Both *bounding_box1* and *bounding_box2* are empty; so we return *True*:
		result = True
	    else:
		# *bounding_box1* is empty, but *bounding_box2* is non-empty; clearly not equal:
		result = False
	else:
	    # *bounding_box1* is non-empty:
	    if bounding_box2._is_empty:
		# *bounding_box1* is non-empty, but *bounding_box2* is empty; clearly not equal:
		result = False
	    else:
		# Both *bounding_box1* and *bounding_box2* are non-empty; compare all 6 pairs
		# bounding box coordinates:
		result = True
		if bounding_box1._east   != bounding_box2._east:
		    result = False
		if bounding_box1._west   != bounding_box2._west:
		    result = False
		if bounding_box1._north  != bounding_box2._north:
                    result = False
		if bounding_box1._south  != bounding_box2._south:
		    result = False
		if bounding_box1._top    != bounding_box2._top:
                    result = False
		if bounding_box1._bottom != bounding_box2._bottom:
		    result = False
	return result

    def __format__(self, format_text):
	""" *Bounding_Box*: Return the values of the *Bounding_Box* object
	    (i.e. *self*) as a string.
	"""

	# Verify argument types:
	assert isinstance(format_text, str)

	if self._is_empty:
	    result = "[None]"
	else:
	    result = "[bsw={0} tne={1}]".format(self.bsw_get(), self.tne_get())
	return result

    def __ne__(self, bounding_box2):
	""" *Bounding_Box*: Return *False* if the *Bounding_Box* object
	    (i.e. *self*) matches *bounding_box2* and *True* otherwise.
	"""
	
	# Verify argument types:
	assert isinstance(bounding_box2, Bounding_Box)

	return not self.__eq__(bounding_box2)

    def bounding_box_expand(self, bounding_box):
	""" *Bounding_Box*: Expand the *Bounding_Box_Object* (i.e. *self*)
	    to contain *bounding_box*.
	"""	

	# Verify argument types:
	assert isinstance(bounding_box, Bounding_Box)

	# Since both boxes are aligned in the X/Y/Z axes, we can get
	# by with just expanding based on the two corners:
	if not extra_bounding_box.is_empty:
	    self.point_expand(extra_bounding_box.tne_get())
	    self.point_expand(extra_bounding_box.bsw_get())

    def b_get(self):
	""" *Bounding_Box*: Return the center of the bottom rectangle of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east.average(self._west), self._north.average(self._south), self._bottom)
	
    def be_get(self):
	""" *Bounding_Box*: Return the center of the bottom/east edge of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east, self._north.average(self._south), self._bottom)

    def bn_get(self):
	""" *Bounding_Box*: Return the center of the bottom/north edge of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east.average(self._west), self._north, self._bottom)

    def bne_get(self):
	""" *Bounding_Box*: Return the center of the bottom/north/east corner of the
	    *Bounding_Box* object (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east, self._north, self._bottom)

    def bnw_get(self):
	""" *Bounding_Box*: Return the center of the bottom/north/west corner of the
	    *Bounding_Box* object (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._west, self._north, self._bottom)

    def bs_get(self):
	""" *Bounding_Box*: Return the center of the bottom/south edge of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east.average(self._west), self._south, self._bottom)

    def bse_get(self):
	""" *Bounding_Box*: Return the center of the bottom/south/eest corner of the
	    *Bounding_Box* object (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east, self._south, self._bottom)

    def bsw_get(self):
	""" *Bounding_Box*: Return the bottom/south/west corner corner of the
	    *Bounding_Box* object (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._west, self._south, self._bottom)

    def bw_get(self):
	""" *Bounding_Box*: Return the center of the bottom/west edge of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._west, self._north.average(self._south), self._bottom)

    def c_get(self):
	""" *Bounding_Box*: Return the center point of the *Bounding_Box* object (i.e. *self*.) """

	return P(self._east.average(self._west),
	  self._north.average(self._south), self._bottom.average(self._top))

    def e_get(self):
	""" *Bounding_Box*: Return the center of east the *Bounding_Box* object (i.e. *self*)
	    rectangle as a point (i.e. *P*) object. """

	return P(self._east, self._north.average(self._south), self._top.average(self._bottom))

    def n_get(self):
	""" *Bounding_Box*: Return the center of north the *Bounding_Box* object (i.e. *self*)
	    rectangle as a point (i.e. *P*) object. """

	return P(self._east.average(self._west), self._north, self._top.average(self._bottom))

    def ne_get(self):
	""" *Bounding_Box*: Return the center of north/east the *Bounding_Box* (i.e. *self*)
	    edge as a point (i.e. *P*) object. """

	return P(self._east, self._north, self._top.average(self._bottom))

    def nw_get(self):
	""" *Bounding_Box*: Return the center of north/east the *Bounding_Box* (i.e. *self*)
	    edge as a point (i.e. *P*) object. """

	return P(self._west, self._north, self._top.average(self._bottom))

    def s_get(self):
	""" *Bounding_Box*: Return the center of south the *Bounding_Box* object (i.e. *self*)
	    rectangle as a point (i.e. *P*) object. """

	return P(self._east.average(self._west), self._south, self._top.average(self._bottom))

    def se_get(self):
	""" *Bounding_Box*: Return the center of south/east the *Bounding_Box* object (i.e. *self*)
	    edge as a point (i.e. *P*) object. """

	return P(self._east, self._south, self._top.average(self._bottom))

    def sw_get(self):
	""" *Bounding_Box*: Return the center of south/west the *Bounding_Box* object (i.e. *self*)
	    edge as a point (i.e. *P*) object. """

	return P(self._west, self._south, self._top.average(self._bottom))

    def w_get(self):
	""" *Bounding_Box*: Return the center of west the *Bounding_Box* object (i.e. *self*)
	    rectangle as a point (i.e. *P*) object. """

	return P(self._west, self._north.average(self._south), self._top.average(self._bottom))

    def t_get(self):
	""" *Bounding_Box*: Return the center of the top rectangle of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east.average(self._west), self._north.average(self._south), self._top)
	
    def te_get(self):
	""" *Bounding_Box*: Return the center of the top/east edge of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east, self._north.average(self._south), self._top)

    def tn_get(self):
	""" *Bounding_Box*: Return the center of the top/north edge of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east.average(self._west), self._north, self._top)

    def tne_get(self):
	""" *Bounding_Box*: Return the center of the top/north/east corner of the
	    *Bounding_Box* object (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east, self._north, self._top)

    def tnw_get(self):
	""" *Bounding_Box*: Return the center of the top/north/west corner of the
	    *Bounding_Box* object (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._west, self._north, self._top)

    def ts_get(self):
	""" *Bounding_Box*: Return the center of the top/south edge of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east.average(self._west), self._south, self._top)

    def tse_get(self):
	""" *Bounding_Box*: Return the center of the top/south/eest corner of the
	    *Bounding_Box* object (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._east, self._south, self._top)

    def tsw_get(self):
	""" *Bounding_Box*: Return the top/south/west corner corner of the
	    *Bounding_Box* object (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._west, self._south, self._top)

    def tw_get(self):
	""" *Bounding_Box*: Return the center of the top/west edge of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object. """

	return P(self._west, self._north.average(self._south), self._top)


    def matrix_apply(self, matrix):
	""" *Bounding_Box*: Return a new *Bounding_Box* object that bounds
	    the current *Bounding_Box* object (i.e. self) after each corner
	    is relocated using *matrix*.
	"""

	# Verify argument types:
	assert isinstance(matrix, Matrix)

	# Contstruct the other 8 points that make up the bounding box:
	tne = self.tne_get()
	tnw = self.tnw_get()
	tse = self.tse_get()
	tsw = self.tsw_get()
	bne = self.bne_get()
	bnw = self.bnw_get()
	bse = self.bse_get()
	bsw = self.bsw_get()

	# Now construct the new *result* bounding box by applying *matrix*
	# to all 8 corners of the original bounding box and return it:
	result = Bounding_Box(
	  [matrix.point_multiply(tne),
	  matrix.point_multiply(tnw),
	  matrix.point_multiply(tse),
	  matrix.point_multiply(tsw),
	  matrix.point_multiply(bne),
	  matrix.point_multiply(bnw),
	  matrix.point_multiply(bse),
	  matrix.point_multiply(bsw)])
	return result

    def point_expand(self, point):
	""" *Bounding_Box*: Expand the *Bounding_Box* object (i.e. *self*)
	    to enclose *point*.
	"""

	# Verify argument types:
	assert isinstance(point, P)

	# Extract *x*/*y*/*z* from *point*:
	x = point.x
	y = point.y
	z = point.z

	# Deal with an empty bounding box first:
	if self._is_empty:
	    # The *Bounding_Box* object (i.e. *self*) so mark it as non_empty and make
            # exactly enclose *point*:
	    self._east = x
	    self._west = x
	    self._north = y
	    self._south = y
	    self._top = z
	    self._bottom = z
	else:
	    # Update the minimum and maximum *x*:
	    self._east =   x.maximum(self._east)
	    self._west =   x.minimum(self._west)
	    self._north =  y.maximum(self._north)
	    self._south =  y.minimum(self._south)
	    self._top =    z.maximum(self._top)
	    self._bottom = z.minimum(self._bottom)

	# Sure that we mark everything as not empty:
	self._is_empty = False

    @staticmethod
    def unit_test():
	zero = L()
	x1 = L(-1.0)
	x2 = L( 1.0)
	y1 = L(-2.0)
	y2 = L( 2.0)
	z1 = L(-3.0)
	z2 = L( 3.0)
	
	point1 = P(x1, y1, z1)
	point2 = P(x2, y2, z2)

	bounding_box = Bounding_Box()
	print("bounding_box 1: {0}".format(bounding_box))
	bounding_box.point_expand(point1)
	print("bounding_box 2: {0}".format(bounding_box))
	bounding_box.point_expand(point2)
	print("bounding_box 3: {0}".format(bounding_box))

	b = bounding_box.b_get()
	be = bounding_box.be_get()
	bn = bounding_box.bn_get()
	bne = bounding_box.bne_get()
	bnw = bounding_box.bnw_get()
	bs = bounding_box.bs_get()
	bse = bounding_box.bse_get()
	bsw = bounding_box.bsw_get()
	bw = bounding_box.bw_get()

	c = bounding_box.c_get()
	e = bounding_box.e_get()
	n = bounding_box.n_get()
	ne = bounding_box.ne_get()
	nw = bounding_box.nw_get()
	s = bounding_box.s_get()
	se = bounding_box.se_get()
	sw = bounding_box.sw_get()
	w = bounding_box.w_get()

	t = bounding_box.t_get()
	te = bounding_box.te_get()
	tn = bounding_box.tn_get()
	tne = bounding_box.tne_get()
	tnw = bounding_box.tnw_get()
	ts = bounding_box.ts_get()
	tse = bounding_box.tse_get()
	tsw = bounding_box.tsw_get()
	tw = bounding_box.tw_get()

	assert b == P(zero, zero, z1), "{0}".format(b)
	assert be == P(x2, zero, z1), "{0}".format(be)
	assert bn == P(zero, y2, z1), "{0}".format(bn)
	assert bne == P(x2, y2, z1), "{0}".format(bne)
	assert bnw == P(x1, y2, z1), "{0}".format(bnw)
	assert bs == P(zero, y1, z1), "{0}".format(bs)
	assert bse == P(x2, y1, z1), "{0}".format(bse)
	assert bsw == P(x1, y1, z1), "{0}".format(bsw)
	assert bw == P(x1, zero, z1), "{0}".format(bw)

	assert c == P(zero, zero, zero), "{0}".format(c)
	assert e == P(x2, zero, zero), "{0}".format(e)
	assert n == P(zero, y2, zero), "{0}".format(n)
	assert ne == P(x2, y2, zero), "{0}".format(ne)
	assert nw == P(x1, y2, zero), "{0}".format(nw)
	assert s == P(zero, y1, zero), "{0}".format(s)
	assert se == P(x2, y1, zero), "{0}".format(se)
	assert sw == P(x1, y1, zero), "{0}".format(sw)
	assert w == P(x1, zero, zero), "{0}".format(w)

	print("t={0}".format(t))
	assert t == P(zero, zero, z2), "{0}".format(t)
	assert te == P(x2, zero, z2), "{0}".format(te)
	assert tn == P(zero, y2, z2), "{0}".format(tn)
	assert tne == P(x2, y2, z2), "{0}".format(tne)
	assert tnw == P(x1, y2, z2), "{0}".format(tnw)
	assert ts == P(zero, y1, z2), "{0}".format(ts)
	assert tse == P(x2, y1, z2), "{0}".format(tse)
	assert tsw == P(x1, y1, z2), "{0}".format(tsw)
	assert tw == P(x1, zero, z2), "{0}".format(tw)

##FIXME: Is this used any more?!!!
#
#class Code_Block:
#    """ *Code_Block*: 
#    """
#
#    def __init__(self, part, tool, vice_x, vice_y):
#	""" *Code_Block*: Initialize the *Code_Block* object (i.e. *self*)
#	    to contain *part*, *tool*, *vice_x*, and *vice_y*.
#	"""
#
#	# Verify argument types:
#	assert isinstance(part, Part)
#	assert isinstance(tool, Tool)
#	assert isinstance(vice_x, L)
#	assert isinstance(vice_y, L)
#
#	self.code = Code()		# Parent *Code* ooject:
#	self.comment = ""		# Comment for block
#	self.speed = Speed()		# Nominal feed rate for this block
#	self.operations = []		# Operations for block
#	self.part = part		# *Part* that owns this *Code_Block*
#	self.priority = 0		# Priority block group
#	self.program_number = 0		# Program number for block
#	self.spindle = Hertz()		# Nominal speed rate for this block
#	self.text = ""			# Final text of block
#	#self.tool = Tool()		# *Tool* associated with block
#	self.uid = 0			# Uniquie id for Block (for debugging)
#	self.vice_x = L()
#	self.vice_y = L()
#
#    def _text_get(self):
#	""" *Code_Block*: Teturn the text field of the *Code_Block* object (i.e. *self*.)
#	"""
#
#	
#	# Load *text* into the *Block* object (i.e. *self*):
#	return self._text
#
#    def _text_set(self, text):
#	""" *Code_Block*: Set the text field of the *Code_Block* object (i.e. *self*) to *text*.
#	"""
#
#	# Verify argument types:
#	assert isinstance(text, str)
#	
#	# Load *text* into the *Code_Block* object (i.e. *self*):
#	self._text = text


# *Speed* class:

class Speed:
    """ *Speed* is a class that represents a speed.
    """

    def __init__(self, mm_per_sec = 0.0, in_per_sec = 0.0, cm_per_sec = 0.0,
      ft_per_min = 0.0, m_per_sec = 0.0):
	""" *Speed*: Initialize a *Speed* object using zero, one or more
	    of *mm_per_sec*, *in_per_sec*, o *cm_per_sec*, *ft_per_sec*,
	    or *m_per_sec*.
	"""

	self._mm_per_sec =			\
	  mm_per_sec +				\
	  in_per_sec * 25.4 +			\
	  cm_per_sec * 10.0 +			\
	  ft_per_min * 60.0 * 12.0 * 25.4 +	\
	  m_per_sec * 60.0 / 1000.0

    def __eq__(speed1, speed2):
	""" *Speed*: Return true if the *speed1* is equal to *speed2. """

	return speed1._mm_per_sec > speed2._mm_per_sec

    def __gt__(speed1, speed2):
	""" *Speed*: Return true if the *speed1* is greater than *speed2. """

	return speed1._mm_per_sec > speed2._mm_per_sec

    def __ge__(speed1, speed2):
	""" *Speed*: Return true if the *speed1* is greater than or equal to *speed2. """

	return speed1._mm_per_sec >= speed2._mm_per_sec

    def __lt__(speed1, speed2):
	""" *Speed*: Return true if the *speed1* is less than *speed2. """

	return speed1._mm_per_sec < speed2._mm_per_sec

    def __le__(speed1, speed2):
	""" *Speed*: Return true if the *speed1* is less than or equal to *speed2. """

	return speed1._mm_per_sec <= speed2._mm_per_sec

    def __ne__(speed1, speed2):
	""" *Speed*: Return true if the *speed1* is not equal to *speed2. """

	return speed1._mm_per_sec != speed2._mm_per_sec

    def __div__(self, divisor):
	""" *Speed*: Return the *Speed* object (i.e. *self) divided by the *divisor* and
	    return a *Hertz* object.
	"""

	# Verify argument types:
	if isinstance(divisor, L):
	    inches_per_second = self.inches_per_second()
	    inches = divisor.inches()
	    rps = inches_per_second / inches
	    hertz =  Hertz(rps=rps)
	    return hertz
	elif isinstance(divisor, int) or isinstance(divisor, float):
	    return Speed(mm_per_sec=self._mm_per_sec / divisor)
	else:
	    assert False, "Bad Speed divsor type"

    def __format__(self, format_text):
	""" *Speed*: Return the *Speed* object (i.e. *self*) as a string.
	"""

	# Figure out what *scale* to use:
	scale = 1.0
	if format_text.endswith('F'):
	    # Feet per minute:
	    scale = 60.0 * 25.4 * 12.0
	    format_text = format_text[:-1]
	elif format_text.endswith('M'):
	    # Meters per minute:
	    scale = 60.0 / 1000.0
	    format_text = format_text[:-1]
	elif format_text.endswith('c'):
	    # Centimeters per second:
	    scale = 10.0
	    format_text = format_text[:-1]
	elif format_text.endswith('i'):
	    # Inches per second:
	    scale =  25.4
	    format_text = format_text[:-1]
	elif format_text.endswith('m'):
	    # Millimeters per second:
	    scale = 1.0
	    format_text = format_text[:-1]
	
	# Format the value as a string and return it:
	value = self._mm_per_sec / scale
	speed_format_text = "{0:" + format_text + "f}"
	#print("base={0} format_text={1} scaled={0} text='{1}'".
	#  format(self._mm_per_sec, format_text, value, speed_format_text))
	return speed_format_text.format(value)

    def millimeters_per_second(self):
	""" *Speed*: Return the speed a *float* in millimeters per second.
	"""

	return self._mm_per_sec

    def centimeters_per_second(self):
	""" *Speed*: Return the speed as a *float* in centimeters per second.
	"""

	return self._mm_per_sec / 10.0

    def inches_per_second(self):
	""" *Speed*: Return the speed as a *float* in inches per second.
	"""

	return self._mm_per_sec / 25.4

    def feet_per_minute(self):
	""" *Speed*: Return the speed as a *float* in feet per minute.
	"""

	return self._mm_per_sec * 25.4 * 12.0 * 60.0

    def meters_per_minute(self):
	""" *Speed*: Return the speed as a *float* in feet per minute.
	"""

	return self._mm_per_sec / 1000.0 * 60.0

class Speed_Range:
    """ *Speed_Range* is a class that represents a range of speeds.
    """

    def __init__(self, low_speed, high_speed):
	""" *Speed_Range*: Initialize a *Speed_Range* object (i.e. *self*)
	    to contain *low_speed* and *high_speed*.
	"""

	# Verify argument types:
	assert isinstance(low_speed, Speed)
	assert isinstance(high_speed, Speed)
	assert low_speed.millimeters_per_second() < high_speed.millimeters_per_second()

	# Load up *self*:
	self._low_speed = low_speed
	self._high_speed = high_speed

    def __format__(self, format_text):
	""" *Speed_Range*: Format a *Speed_Range* object (i.e. *self*) into
	    a text string and return it.
	"""

	# Verify argument types:
	assert isinstance(format_text, str)

	# Return a string of the form "low-high":
	range_format_text = "{0:" + format_text + "}-{1:" + format_text + "}"
	#print("Speed_Range.__format__(*, '{0}'): '{1}'".
	#  format(format_text, range_format_text))
	return range_format_text.format(self._low_speed, self._high_speed)

    def _high_speed_get(self):
	""" *Speed_Range*: Return the high speed for the *Speed_Range* object (i.e. *self*). """

	return self._high_speed

    def _low_speed_get(self):
	""" *Speed_Range*: Return the low speed for the *Speed_Range* object (i.e. *self*). """

	return self._low_speed

#class Code_Time:
#
#    def __init__(self, seconds = 0.0):
#	""" *Code_Time*: Initialize a *Code_Time* object. """
#
#	self.seconds = seconds	# Everything is converted to secconds.

class Color:
    COLORS = {
      "alice_blue": 0xf0f8ff,
      "antique_white": 0xfaebd7,
      "aqua": 0x00ffff,
      "aquamarine": 0x7fffd4,
      "azure": 0xf0ffff,
      "beige": 0xf5f5dc,
      "bisque": 0xffe4c4,
      "black": 0x000000,
      "blanched_almond": 0xffebcd,
      "blue": 0x0000ff,
      "blue_violet": 0x8a2be2,
      "brown": 0xa52a2a,
      "burlywood": 0xdeb887,
      "cadet_blue": 0x5f9ea0,
      "chartreuse": 0x7fff00,
      "chocolate": 0xd2691e,
      "coral": 0xf08080,
      "corn_flower_blue": 0x6495ed,
      "corn_silk": 0xfff8dc,
      "crimson": 0xdc143c,
      "cyan": 0x00ffff,
      "dark_blue": 0x00008b,
      "dark_cyan": 0x008b8b,
      "dark_goldenrod": 0xb8860b,
      "dark_gray": 0xa9a9a9,
      "dark_green": 0x006400,
      "dark_grey": 0xa9a9a9,
      "dark_khaki": 0xbdb76b,
      "dark_magenta": 0x8b008b,
      "dark_olive_green": 0x556b2f,
      "dark_orange": 0xff8c00,
      "dark_orchid": 0x9932cc,
      "dark_red": 0x8b0000,
      "dark_salmon": 0xe9967a,
      "dark_sea_green": 0x8fbc8f,
      "dark_slate_blue": 0x483d8b,
      "dark_slate_gray": 0x2f4f4f,
      "dark_slate_grey": 0x2f4f4f,
      "dark_turquoise": 0x40e0d0,
      "dark_violet": 0x9f00d3,
      "deep_pink": 0xff1493,
      "deep_sky_blue": 0x00bfff,
      "dim_gray": 0x696969,
      "dim_grey": 0x696969,
      "dodger_blue": 0x1e90ff,
      "fire_brick": 0xb22222,
      "floral_white": 0xfffaf0,
      "forest_green": 0x228b22,
      "fuchsia": 0xff00ff,
      "gainsboro": 0xdcdcdc,
      "ghost_white": 0xf8f8ff,
      "gold": 0xffd700,
      "goldenrod": 0xdaa520,
      "gray": 0x808080,
      "green": 0x008000,
      "green_yellow": 0xadff2f,
      "grey": 0x808080,
      "honey_dew": 0xf0fff0,
      "hot_pink": 0xff1493,
      "indian_red": 0xcd5c5c,
      "indigo": 0x4b0082,
      "ivory": 0xfffff0,
      "khaki": 0xf0e68c,
      "lavender": 0xe6e6fa,
      "lavender_blush": 0xfff0f5,
      "lawn_green": 0x7cfc00,
      "lemon_chiffon": 0xfffacd,
      "light_blue": 0xadd8e6,
      "light_coral": 0xf08080,
      "light_cyan": 0xe0ffff,
      "light_goldenrod_yellow": 0xfafad2,
      "light_gray": 0xd3d3d3,
      "light_green": 0x90ee90,
      "light_grey": 0xd3d3d3,
      "light_pink": 0xffb6c1,
      "light_salmon": 0xffa07a,
      "light_sea_green": 0x20b2aa,
      "light_sky_blue": 0x87cefa,
      "light_slate_gray": 0x778899,
      "light_slate_grey": 0x778899,
      "light_steel_blue": 0xb0c4de,
      "light_yellow": 0xffffe0,
      "lime": 0x00ff00,
      "lime_green": 0x2e8b57,
      "linen": 0xfaf0e6,
      "magenta": 0xff00ff,
      "maroon": 0x800000,
      "medium_aquamarine": 0x66cdaa,
      "medium_blue": 0x0000cd,
      "medium_orchid": 0xba55d3,
      "medium_purple": 0x9370db,
      "medium_sea_green": 0x3cb371,
      "medium_slate_blue": 0x66cdaa,
      "medium_spring_green": 0x00fa9a,
      "medium_turquoise": 0x48d1cc,
      "medium_violet_red": 0xc71585,
      "mid_night_blue": 0x191970,
      "mint_cream": 0xf5fffa,
      "misty_rose": 0xffe4e1,
      "moccasin": 0xffe4b5,
      "navajo_white": 0xffdead,
      "navy": 0x000080,
      "old_lace": 0xfdf5e6,
      "olive": 0x808000,
      "olive_drab": 0x6b8e23,
      "orange": 0xffa500,
      "orange_red": 0xff4500,
      "orchid": 0xda70d6,
      "pale_goldenrod": 0xeee8aa,
      "pale_green": 0x98fb98,
      "pale_turquoise": 0xafeeee,
      "pale_violet_red": 0xdb7093,
      "papaya_whip": 0xffefd5,
      "peach_puff": 0xffdab9,
      "peru": 0xcd8f3f,
      "pink": 0xffc0cb,
      "plum": 0xdda0dd,
      "powder_blue": 0xb0e0e6,
      "purple": 0x800080,
      "red": 0xff0000,
      "rosy_brown": 0xbc8f8f,
      "royal_blue": 0x4169e1,
      "saddle_brown": 0x8b2be2,
      "salmon": 0xfa8072,
      "sandy_brown": 0xf4a460,
      "sea_green": 0x2e8b57,
      "sea_shell": 0xfff5ee,
      "sienna": 0xa0522d,
      "silver": 0xc0c0c0,
      "sky_blue": 0x87ceeb,
      "slate_blue": 0x6a5acd,
      "slate_gray": 0x708090,
      "slate_grey": 0x708090,
      "snow": 0xfffafa,
      "spring_green": 0x00ff7f,
      "steel_blue": 0x4682b4,
      "tan": 0xd2b48c,
      "teal": 0x008080,
      "thistle": 0xd8bfd8,
      "tomato": 0xff6347,
      "turquoise": 0x40e0d0,
      "violet": 0xee82ee,
      "wheat": 0xf5deb3,
      "white": 0xffffff,
      "white_smoke": 0xf5f5f5,
      "yellow": 0xffff00,
      "yellow_green": 0x9acd32,
    }

    def __init__(self,
      name = "", red = -1.0, green = -1.0, blue = -1.0, alpha = 1.0):
	""" *Color*: Initialize *self*"""

	# Deal with SVG color:
	colors = Color.COLORS
	if name in colors:
	    rgb = colors[name]
	    red_byte = (rgb >> 16) & 0xff
	    green_byte = (rgb >> 8) & 0xff
	    blue_byte = rgb & 0xff
	    #print("Color bytes=({0}, {1}, {2})". \
	    #  format(red_byte, green_byte, blue_byte))
            red = float(red_byte) / 255.0
	    green = float(green_byte) / 255.0
	    blue = float(blue_byte) / 255.0

	# Keep red/green/blue/alpha between 0.0 and 1.0:
	if red < 0.0:
	    red = 0.0
	elif red > 1.0:
	    red = 1.0
	if green < 0.0:
	    green = 0.0
	elif green > 1.0:
            green = 1.0
	if blue < 0.0:
            blue = 0.0
	elif blue > 1.0:
	    blue = 1.0
	if alpha < 0.0:
	    alpha = 0.0
	elif alpha > 1.0:
	    alpha = 1.0

	# Load up *self*:
	self.red = red
	self.green = green
	self.blue = blue
	self.alpha = alpha

	#print("Color({0}, {1}, {2}, {3})".
	#  format(self.red, self.green, self.blue, self.alpha))

    def __format__(self, format):
	""" *Color*: Return formatted version of *self*. """

	return "[r={0:.2f}, g={1:.2f}, b={2:.2f}, a={3:.2f}]".format(
	  self.red, self.green, self.blue, self.alpha)

# These classes are unused???

class Indexed_Point:
    def __init__(self, x, y, label, index):
	""" *Indexed_Point*: """
	self.x = x
	self.y = y
	self.label = label
	self.index = index

class Indexed_Points:
    def __init__(self):
	self._indexed_points = []
	self._table = {}
	self._paths = []

    def lookup(self, x, y, label):
	# Check argument types:
	assert isinstance(x, float)
	assert isinstance(y, float)
	assert isinstance(label, str)

	# Is *point* alread in *table*?:
	key = (x, y)
	table = self._table
	if key in table:
	    # Yes; lookup previously stored index:
	    indexed_point = table[key]
	else:
	    # No; create new *index* and store it away:
	    indexed_points = self._indexed_points
	    index = len(indexed_points)
	    indexed_point = Indexed_Point(x, y, label, index)
	    indexed_points.append(indexed_point)
	    table[key] = indexed_point
	return indexed_point

    def polygon_append(self, lines, indent):
	# Check argument types:
	assert isinstance(lines, list)
	assert isinstance(indent, int) and indent >= 0

	# Output the .scad contents:
	spaces = " " * indent
	lines.append("{0}polygon(".format(spaces))
	lines.append("{0} points = [".format(spaces))
	self.points_append(lines, spaces)
	lines.append("{0} ], paths = [".format(spaces))
	self.paths_append(lines, spaces)
	lines.append("{0} ]);".format(spaces))

    def path_create(self):
	path = []
	self._paths.append(path)
	return path

    def points_append(self, lines, spaces):
	# Check argument types:
	assert isinstance(lines, list)
	assert isinstance(spaces, str)

	indexed_points = self._indexed_points
	size = len(indexed_points)
	for index in range(size):
	    indexed_point = indexed_points[index]
	    suffix = ","
	    if index >= size - 1:
		suffix = ""
	    lines.append("{0}  [{1:.3f}, {2:.3f}]{3}".
	      format(spaces, indexed_point.x, indexed_point.y, suffix))

    def paths_append(self, lines, spaces):
	paths = self._paths
	paths_size = len(paths)
	for paths_index in range(paths_size):
	    path = paths[paths_index]
	    path_size = len(path)
	    lines.append("{0}  [".format(spaces))
	    for path_index in range(path_size):
		indexed_point = path[path_index]
		path_suffix = ","
		if path_index >= path_size - 1:
		    path_suffix = ""
		lines.append("{0}    {1}{2}".
		  format(spaces, indexed_point.index, path_suffix))

	    paths_suffix = ","
	    if paths_index >= paths_size + 1:
		paths_suffix = ""
	    lines.append("{0}  ]{1}".format(spaces, paths_suffix))

class Contour:

    def __init__(self, name = None):
	""" *Contour*: Initialize *self* with *bends* and *extrude_axis*. """

	# Load *self*
	self._bends = []
	self._name = "No Name"
	if isinstance(name, str):
	    self._name = name

    def __format__(self, format):
	""" *Contour*: """
	result = "{"
	bends = self._bends
	for index in range(len(bends)):
	    bend = bends[index]
	    result += "[{0}]: {1}, {2} ". \
              format(index, bend.point, bend.radius)
	return result + "}"

    def adjust(self, delta = L(), start = P(), end = P(),
      maximum_radius = L(mm = 1.e10), minimum_radius = L(mm = 0)):
	""" *Contour*: """

	# Check argument types:
	assert isinstance(delta, L)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(maximum_radius, L)
	assert isinstance(minimum_radius, L)

	trace = False
	#trace = True
	if trace:
	    print("=>Contour.adjust():len(bends)={0}".format(len(self._bends)))

	extrude_axis = end - start

	# The process is to loop around *bend* and adjust segment in or out
	# by *delta* depending upon when the the sign of *delta* is negative
	# or positive.  The following steps occur:
	#
	# 1. Take a pair of Bends in sequence to get the end-points of
	#    a line.  Compute the direction vector and a perpendicular
	#    vector.  Move each end-point by *delta* along the perpendicular.
	#
	# 2. Find the intersections of each of the adjusted line segments.
	#    This is done via the determanent method describe at
	#
	#        http://mathworld.wolfram.com/Line-LineIntersection.html
	#
	# 3. Adjust the bend radius up or down appropriately.

	self.project(axis = extrude_axis)

	adjusted_contour = Contour("Adjusted {0}".format(self._name))
	bends = self._bends
	bends_size = len(bends)
	for bends_index in range(bends_size):
	    # The intersection computation is easier to transcribe if we
	    # use the same notation as used in mathwolrd.com description.
	    # The first line goes through points (x1, y1) and (x2, y2) and
	    # the second line goes through points (x3, y3) and (x4, y4).
	    # We will use (*ox1*, *oy1*), (*ox2*, *oy2*), (*ox3*, *oy3*)
	    # and (*ox4*, *oy4*) to represent the [o]iginal points.
	    # We will use (*ax1*, *ay1*), (*ax2*, *ay2*), (*ax3*, *ay3*)
	    # and (*ax4*, *ay4*) to represent the [a]djusted points.
            # Note that using this notation, (*ox2*, *oy2*) is the same
	    # as (*ox3*, *oy3*), but (*ax2*, *ay2*) is not the same as
	    # (*ax3*, *ay3).

	    # *bend1* occurs *bend2* occurs before *bend3*:
	    bend1 = bends[(bends_index - 1) % bends_size]
	    bend2 = bends[bends_index]
	    bend3 = bends[(bends_index + 1) % bends_size]
	    if trace:
		print("bends_index={0}".format(bends_index))
	    
	    # Get the [o]riginal X/Y values:
	    ox1 = bend1._px
	    oy1 = bend1._py
	    ox2 = bend2._px
	    oy2 = bend2._py
	    ox3 = bend2._px
	    oy3 = bend2._py
	    ox4 = bend3._px
	    oy4 = bend3._py
	    if trace:
		print(
		  "[{0}]:o1=({1},{2}) o2=({3},{4}) o3=({5},{6}) o4=({7},{8})".
		  format(bends_index, ox1, oy1, ox2, oy2, ox3, oy3, ox4, oy4))

	    # Compute vector (*odx12*, *ody12*) from (*ox1*, *oy1*) to
	    # (*ox2*, *oy2*):
	    odx12 = ox2 - ox1
	    ody12 = oy2 - oy1

	    # Compute the vector (*odx34*, *ody34*) from (*ox3*, *oy3*) to
	    # (*ox4*, *oy4*):
	    odx34 = ox4 - ox3
	    ody34 = oy4 - oy3
	    if trace:
		print("[{0}]:od12=({1},{2}) od34=({3},{4})".
		  format(bends_index, odx12, ody12, odx34, ody34))

	    # Compute the lengths of (*dx12*, *dy12*) and (*dx34*, *dy34*):
	    length12 = math.sqrt(odx12 * odx12 + ody12 * ody12)
	    length34 = math.sqrt(odx34 * odx34 + ody34 * ody34)
	    if trace:
		print("[{0}]:len12={1} len34={2}".
		  format(bends_index, length12, length34))

	    # Due to constraint propagation, sometimes we get some points
	    # that land on top of one another.  Just fake it for now:
	    if length12 <= 0.0:
		length12 = 0.1
	    if length34 <= 0.0:
		length34 = 0.1

	    # Compute perpendicular normal (*ndx12*, *ndy12*) from
	    # (*dx12*, *dy12*) and normal (*ndx34*, *ndy34*) from
	    # (*dx34*, *dy34*):
	    ndx12 =  ody12 / length12
	    ndy12 = -odx12 / length12
	    ndx34 =  ody34 / length34
	    ndy34 = -odx34 / length34
	    if trace:
		print("[{0}]:n12=({1},{2}) n34=({3},{4})".
		  format(bends_index, ndx12, ndy12, ndx34, ndy34))

	    # Compute [a]djusted bend points (*ax1*, *ay1*), (*ax2*, *ay2*),
	    # (*ax3*, *ay3*), and (*ax4*, *ay4*):
	    d = delta._mm
	    ax1 = ox1 + d * ndx12
	    ay1 = oy1 + d * ndy12
	    ax2 = ox2 + d * ndx12
	    ay2 = oy2 + d * ndy12
	    ax3 = ox3 + d * ndx34
	    ay3 = oy3 + d * ndy34
	    ax4 = ox4 + d * ndx34
	    ay4 = oy4 + d * ndy34
	    if trace:
		print(
		  "[{0}]:a1=({1},{2}) a2=({3},{4}) a3=({5},{6}) a4=({7},{8})".
		  format(bends_index, ax1, ay1, ax2, ay2, ax3, ay3, ax4, ay4))

	    # Now find the intersection of (*x1a*, *y1a*) - (*x2a*, *y2a*)
	    # and (*x2b*, *y2b*) - (*x3b*, *y3b*).  The computed formulas
	    # are taken from:
	    #
	    #    http://mathworld.wolfram.com/Line-LineIntersection.html
	    #
	    #
	    #        | |x1 y1|        |        | |x1 y1|        |
	    #        | |     |  x1-x2 |        | |     |  y1-y2 |
	    #        | |x2 y2|        |        | |x2 y2|        |
	    #        |                |        |                |
	    #        | |x3 y3|        |        | |x3 y3|        |
	    #        | |     |  x3-x4 |        | |     |  y3-y4 |
	    #        | |y4 y4|        |        | |y4 y4|        |
	    #    x = ------------------    y = ------------------
	    #        |  x1-x2   y1-y2 |        |  x1-x2   y1-y2 |
	    #        |  x3-x4   y3-y4 |        |  x3-x4   y3-y4 |
	    #
	    
	    # Compute *det12* and *det34*:
	    det12 = ax1 * ay2 - ax2 * ay1	# = det2(x1, y1, x2, y2)
	    det34 = ax3 * ay4 - ax4 * ay3	# = det2(x3, y3, x4, y4)
	    # Compute x numerator = det2(det12, x1 - x2, det34, x3 - x4):
	    x_numerator = det12 * (ax3 - ax4) - det34 * (ax1 - ax2)
	    # Compute y numerator = det2(det12, y1 - y2, det34, y3 - y4):
	    y_numerator = det12 * (ay3 - ay4) - det34 * (ay1 - ay2)
	    # Compute denomonator = det2(x1 - x2, y1 - y2, x3 - x4, y3 - y4):
	    denominator = (ax1 - ax2) * (ay3 - ay4) - (ax3 - ax4) * (ay1 - ay2)
	    # Sometimes constraint propagation causes a *denominator* of 0.
	    # Set it to non-zero to fake it:
	    if denominator == 0.0:
		denominator = 0.01
	    # Compute final [i]ntersection location (*ix*, *iy):
	    ix = x_numerator / denominator
	    iy = y_numerator / denominator
	    if trace:
		print("[{0}]:det12={1} det34={2}".
		  format(bends_index, det12, det34))
		print("[{0}]:x_num={1} y_num={2} denom={3}".
		  format(bends_index, x_numerator, y_numerator, denominator))
		print("[{0}]:ix={1} iy={2}".format(bends_index, ix, iy))

	    # Now we get to adjust the radius.  When *delta* < 0 we are
	    # making the contour smaller.  The direction of the turn is
	    # what we use to tell the difference:
	    bearing12 = math.atan2(ody12, odx12)
	    bearing34 = math.atan2(ody34, ody34)
	    delta_bearing = bearing34 - bearing12
	    pi = math.pi
	    if delta_bearing < pi:
		delta_bearing += 2.0 * pi
	    elif delta_bearing > pi:
		delta_bearing -= 2.0 * pi

	    new_radius = bend2.radius._mm
	    if delta_bearing < 0.0:
		# Turn clockwise:
		new_radius += delta._mm
	    else:
		# Turn counter-clockwise:
		new_radius -= delta._mm

	    if new_radius < minimum_radius._mm:
		new_radius = minimum_radius._mm
	    elif new_radius > maximum_radius._mm:
		new_radius = maximum_radius._mm

	    zero = L()
	    adjusted_bend = adjusted_contour.bend_append(
	      P(L(mm = ix), L(mm = iy), zero), L(mm = new_radius),
	      plane = bend2._plane)
	    #print("Adjusted_Bend[px={0},py={1} ix={2} iy={3}]".
	    #  format(adjusted_bend._px, adjusted_bend._py, ix, iy))

	# Reproject the contour into a compatible plane:
	adjusted_contour.unproject(axis = extrude_axis)

	if trace:
	    print("<=Contour.adjust():len(bends)={0}".format(len(self._bends)))
	return adjusted_contour

    def bend_append(self, point = P(), radius = L(), plane = None, name = ""):
	""" *Contour*: """

	# Check argument types:
	assert isinstance(point, P)
	assert isinstance(radius, L)
	assert isinstance(name, str)

	# Create and append the bend:
	bend = Bend(point, radius, name = name, plane = plane)
	self._bends.append(bend)
	return bend

    def bends_compute(self, axis = None, trace = -1000000):
	""" *Contour*: """

	if trace >= 0:
	    print("{0}=>Contour.bends_compute(axis={1}, self={2})".
	      format(trace * ' ', axis, self))

	# Check argument types:
	assert isinstance(axis, P)

	# Project all the 3D points down to 2D:
	self.project(axis)

	# Find the bend centers:
	bends = self._bends
	bends_size = len(bends)
	for index in range(bends_size):
	    # Extract three *Bend*'s in sequence
	    before_bend = bends[(index - 1) % bends_size]
	    bend = bends[index]
	    after_bend = bends[(index + 1) % bends_size]

	    # Compute the bend radius center point:
	    #print("bends_compute:before={0} at={1} after={2}".
	    #  format(before_bend, bend, after_bend))
	    bend.compute(before_bend, after_bend, trace = trace + 1)
	if trace >= 0:
	    print("{0}<=Contour.bends_compute(axis={1}, self={2})".
	      format(trace * ' ', axis, self))

    def bounding_box_compute(self, start, end):
	""" *Contour*: Compute bounding box for *self* extruded in
	    *extrude_axis* direction. """

	# Check argument _types:
	big = L(mm=987654321.0)

	# Initialize the X/Y/Z minimum/maximum values:
	ex = -big
	wx = big
	ny = -big
	sy = big
	tz = -big
	bz = big

	# Grap *dx*, *dy*, and *dz* from *extrude_axis*:
	start_x = start.x
	start_y = start.y
	start_z = start.z
	end_x = end.x
	end_y = end.y
	end_z = end.z

	for bend in self._bends:
	    # Extract *x*, *y*, *z* from *point*:
	    point = bend.point_tmp
	    x = point.x
	    y = point.y
	    z = point.z

	    # Adjust the bounding box bounaries:
	    ex = ex.maximum(start_x + x).maximum(end_x + x)
	    wx = wx.minimum(start_x + x).minimum(end_x + x)
	    ny = ny.maximum(start_y + y).maximum(end_y + y)
	    sy = sy.minimum(start_y + y).minimum(end_y + y)
	    tz = tz.maximum(start_z + z).maximum(end_z + z)
	    bz = tz.minimum(start_z + z).minimum(end_z + z)
	    
	# Return the 8 points that are the bounding box of *self*:
	bounding_box = Bounding_Box(ex, wx, ny, sy, tz, bz)
	return bounding_box

    def path_append(self, extrude_axis = None, indexed_points = None,
      maximum_angle = Angle(deg=16.0), trace = -1000000):
	""" *Contour*: Append a path to *indexed_points*. """
	assert isinstance(indexed_points, Indexed_Points)
	assert isinstance(maximum_angle, Angle)
	assert isinstance(extrude_axis, P)
	assert isinstance(trace, int)

	if trace >= 0:
	    print("{0}=>Contour.path_append(extrude_axis={1}, self={2})".
	          format(trace * ' ', extrude_axis, self))

	# Compute the bend information:
	self.bends_compute(extrude_axis, trace = trace + 1)

	pi = math.pi
	r2d = 180.0 / pi 

	path = indexed_points.path_create()
	bends = self._bends
	bends_size = len(bends)
	for index in range(bends_size):
	    before_bend = bends[(index - 1) % bends_size]
	    bend = bends[index]
	    #bend_after = bends[(index + 1) % bends_size]
	    before_total_fraction = \
	      before_bend._after_fraction + bend._before_fraction
	    if before_total_fraction > 1.00001:
		#assert False,
		#  "We have a bogus edge {0} + {1} = {2} > 1.00001". \
		#  format(before_bend._after_fraction, bend._before_fraction,
		#  before_total_fraction)
		pass

	    #print("Bend[{0}]:{1}".format(index, bend))

	    radius = bend.radius._mm
            after_tangent_x = bend._after_tangent_x
            after_tangent_y = bend._after_tangent_y
            before_tangent_x = bend._before_tangent_x
            before_tangent_y = bend._before_tangent_y
	    center_x = bend._center_x
	    center_y = bend._center_y

	    if before_total_fraction < 1.0:
		indexed_point = indexed_points.lookup(
		  before_tangent_x, before_tangent_y,
		  format("Bend[{0}]:before tangent".format(index)))
		path.append(indexed_point)

	    before_dy = before_tangent_y - center_y
	    before_dx = before_tangent_x - center_x
	    before_angle = math.atan2(before_dy, before_dx)

	    after_dy = after_tangent_y - center_y
	    after_dx = after_tangent_x - center_x
	    after_angle = math.atan2(after_dy, after_dx)

	    #print("before_angle={0:.4f} after_angle={1:.4f}".
	    #  format(before_angle * r2d, after_angle * r2d))

	    # Compute and normalize *delta_angle* from *before_angle* to
	    # *after_angle*.  *delta_angle* could wind up being positiver
	    # or negative:
	    delta_angle = after_angle - before_angle
	    while delta_angle > pi:
		delta_angle -= 2.0 * pi
            while delta_angle < -pi:
		delta_angle += 2.0 * pi

	    # We want to divide *delta_angle* by an integer *step_count*
	    # such that the resulting *step_angle* is as close to
	    # *maximum_angle* as possible.  We start by simply dividing
	    # computing *fractional* by divided *delta_angle* by
	    # *maximum_angle*:
	    max_angle = maximum_angle.radians	
	    fractional = delta_angle / max_angle
	    #print(
	    #  "delta_angle={0:.4f} max_angle={1:.4f} fractional={2:.4f}".
	    #  format(delta_angle * r2d, max_angle * r2d, fractional))

	    # Compute the *step_count* which is the number of chunks
	    # we divide *delta_angle* by to get the *step_angle*.  We
	    # want *step_angle* to be as close to *maximum_angle* without
	    # going over.  Remember *step_angle* can be positive or negative:
	    step_count = int(abs(fractional)) + 1
	    step_angle = delta_angle / float(step_count)
	    #print("step_count={0} step_angle={1:.4f}".
	    #  format(step_count, step_angle * r2d))

	    # Now output the step-wise approximation for the arc:
	    for step_index in range(1, step_count):
		angle = before_angle + float(step_index) * step_angle
		x = center_x + radius * math.cos(angle)
		y = center_y + radius * math.sin(angle)
		indexed_point = indexed_points.lookup(x, y,
		  "Bend[{0}]:Step[{1}]".format(index, step_index))
		path.append(indexed_point)

	    if radius > 0.0:
		indexed_point = indexed_points.lookup(
		  after_tangent_x, after_tangent_y,
		  "Bend[{0}]:after tangent".format(index))
		path.append(indexed_point)

	if trace >= 0:
	    print("{0}<=Contour.path_append(extrude_axis={1})".
	      format(trace * ' ', extrude_axis))

    def project(self, axis = None):
	""" *Contour*: """

	#print("=>Contour.project(axis={0}".format(axis))

	# Check argument types:
	assert isinstance(axis, P)
	
	# Extract X/Y/Z values from *axis*:
	axis_x = axis.x._mm
	axis_y = axis.y._mm
	axis_z = axis.z._mm

	# Figure out which axis to use:
	axis_x_is_zero = axis_x == 0.0
	axis_y_is_zero = axis_y == 0.0
	axis_z_is_zero = axis_z == 0.0

	# Figure out which axis to use:
	use_x_axis = not axis_x_is_zero and axis_y_is_zero and axis_z_is_zero
	use_y_axis = axis_x_is_zero and not axis_y_is_zero and axis_z_is_zero
	use_z_axis = axis_x_is_zero and axis_y_is_zero and not axis_z_is_zero

	# Remember these values:
	self.use_x_axis = use_x_axis
	self.use_y_axis = use_y_axis
	self.use_z_axis = use_z_axis

	# For now, constrain *axis* to be in X, Y, or Z direction:
	if not (use_x_axis or use_y_axis or use_z_axis):
	    print("Contour.project: axis must align with X, Y, or Z axis: {0}".
	      format(axis))
	    use_z_axis = True
	    assert EZCAD3.ezcad._update_count != 4

	# For each *bend* in bends, perform the project from 3D down to 2D:
	plane = None
	for bend in self._bends:
	    # Grab the X/Y/Z coordinates for *point*:
	    point = bend.point
	    px = point.x._mm
	    py = point.y._mm
	    pz = point.z._mm
	    #print("Contour.project[{0}]: px/py/pz = {1}/{2}/{3}".
	    #  format(self._name, px, py, pz))

	    if use_x_axis:
		# Axis is aligned with X:
		# This is not obvious.  We put the Z coordinates in the
		# X axis in order to ensure that when we rotate around
		# the Y axis by 90 degrees the Z cordinates are up and
		# down in Z:
		bend._px = pz
		bend._py = py
		if plane == None:
		    plane = px
		else:
		    assert plane == px, \
		      "All points must be in same X plane (plane={0} px={1})". \
		      format(plane, px)
	    elif use_y_axis:
		# Axis is aligned with Y:
		bend._px = px
		bend._py = pz
		if plane == None:
                    plane = py
		else:
		    assert plane == py, \
		      "All points must be in same Y plane"
	    elif use_z_axis:
		# Axis is aligned with Y:
		bend._px = px
		bend._py = py
		if plane == None:
                    plane = pz
		else:
		    assert plane == pz, \
		      "All points must be in same Y plane"
            else:
		assert False, "FIXME: Allow arbitray axis direction"
            bend._plane = plane

	# Keep track of *plane*
	self.plane = plane

	#print("<=Contour.project(axis={0}".format(axis))

    def unproject(self, axis = None):
	""" *Contour*: """

	#print("=>Contour.project(axis={0}".format(axis))

	# Check argument types:
	assert isinstance(axis, P)
	
	# Extract X/Y/Z values from *axis*:
	axis_x = axis.x._mm
	axis_y = axis.y._mm
	axis_z = axis.z._mm

	# Figure out which axis to use:
	axis_x_is_zero = axis_x == 0.0
	axis_y_is_zero = axis_y == 0.0
	axis_z_is_zero = axis_z == 0.0

	# Figure out which axis to use:
	use_x_axis = not axis_x_is_zero and axis_y_is_zero and axis_z_is_zero
	use_y_axis = axis_x_is_zero and not axis_y_is_zero and axis_z_is_zero
	use_z_axis = axis_x_is_zero and axis_y_is_zero and not axis_z_is_zero

	# Remember these values:
	self.use_x_axis = use_x_axis
	self.use_y_axis = use_y_axis
	self.use_z_axis = use_z_axis

	# For now, constrain *axis* to be in X, Y, or Z direction:
	if not (use_x_axis or use_y_axis or use_z_axis):
	    print("Contour.project: axis must align with X, Y, or Z axis: {0}".
	    format(axis))
	    use_z_axis = True
	    #assert EZCAD3.ezcad._update_count != 4

	# For each *bend* in bends, perform the project from 2D back up 3D:
	for bend in self._bends:
	    # Grab the X/Y/Z coordinates for *point*:
            point = bend.point
	    point_x = point.x
	    point_y = point.y
	    plane = bend._plane

	    #print("Contour.unproject[{0}]: point_x/point_y = {1}/{2}".
	    #  format(self._name, point_x, point_y))

	    if use_x_axis:
		# Axis is aligned with X:
		# This is not obvious.  We put the Z coordinates in the
		# X axis in order to ensure that when we rotate around
		# the Y axis by 90 degrees the Z cordinates are up and
		# down in Z:
		point = P(plane, point_y, point_x)
	    elif use_y_axis:
		# Axis is aligned with Y:
		point = P(point_x, plane, point_y)
	    elif use_z_axis:
		# Axis is aligned with Y:
		point = P(point_x, point_y, plane)
            else:
		assert False, "FIXME: Allow arbitray axis direction"
	    bend.point = point

	#print("<=Contour.project(axis={0}".format(axis))

class EZCAD3:
    """ EZCAD3 is the top level engine that executes the design. """

    DIMENSIONS_MODE = 0
    CNC_MODE = 1
    STL_MODE = 2
    VISUALIZATION_MODE = 3

    def __init__(self, minor = None, adjust = L(), directory="."):
	""" *EZCAD*: Initialize the contents of {self} to contain
	    *major* and *minor* version numbers. """

	# FIXME: *adjust* should be moved over to the *Part* library!!!
	# It is used for increasing the size of things in STL mode.

	#print "EZCAD.__init__() called"

	# Check argument types:
	none_type = type(None)
	assert minor == 0
	assert isinstance(adjust, L)
	assert isinstance(directory, str)

	# Load up {self}:
	self._adjust = adjust
	self._directory = directory
	self._mode = EZCAD3.DIMENSIONS_MODE
	self._major = 3
	self._minor = minor
	self._parts_stack = []
	self._shop = Shop("Wayne's Shop")
	self._update_count = 0
	self._xml_indent = 0
	self._xml_stream = None
	self._bases_table = {}

	self._bounding_box_dispatch = {
	  "b":   Bounding_Box.b_get,
	  "be":  Bounding_Box.be_get,
	  "bn":  Bounding_Box.bn_get,
	  "bne": Bounding_Box.bne_get,
	  "bnw": Bounding_Box.bnw_get,
	  "bs":  Bounding_Box.bs_get,
	  "bse": Bounding_Box.bse_get,
	  "bsw": Bounding_Box.bsw_get,
	  "bw":  Bounding_Box.bw_get,
	  "c":   Bounding_Box.c_get,
	  "e":   Bounding_Box.e_get,
	  "n":   Bounding_Box.n_get,
	  "ne":  Bounding_Box.ne_get,
	  "nw":  Bounding_Box.nw_get,
	  "s":   Bounding_Box.s_get,
	  "se":  Bounding_Box.se_get,
	  "sw":  Bounding_Box.sw_get,
	  "t":   Bounding_Box.t_get,
	  "te":  Bounding_Box.te_get,
	  "tn":  Bounding_Box.tn_get,
	  "tne": Bounding_Box.tne_get,
	  "tnw": Bounding_Box.tnw_get,
	  "ts":  Bounding_Box.ts_get,
	  "tse": Bounding_Box.tse_get,
	  "tsw": Bounding_Box.tsw_get,
	  "tw":  Bounding_Box.tw_get,
	  "w":   Bounding_Box.w_get,
	}


	EZCAD3.ezcad = self

    def _directory_get(self):
	""" *EZCAD3*: Return the directory to read/write files from/into from the *EZCAD3* object
	    (i.e. *self*).
        """

	return self._directory

    def process(self, part):
	""" *EZCAD3*: Perform all of the processing starting at *part*.
	"""

	assert isinstance(part, Part)
	part.process(self)

class Operation:
    """ *Operation* is a class that represents a manufacturing operation.
    """

    # These constants are used to sort the order of operations between
    # operations:
    ORDER_NONE =			0
    ORDER_DOWEL_PIN =			1
    ORDER_END_MILL_EXTERIOR =		2
    ORDER_MILL_DRILL_EXTERIOR =		3
    ORDER_MILL_DRILL_CHAMFER =		4
    ORDER_MILL_DRILL_COUNTERSINK =	5
    ORDER_MILL_DRILL =			6
    ORDER_END_MILL_DRILL =		7
    ORDER_END_MILL_ROUND_POCKET =	8
    ORDER_END_MILL_SIMPLE_POCKET =	9
    ORDER_MILL_DRILL_POCKET_CHAMFER =	10
    ORDER_MILL_DOVE_TAIL_CHAMFER = 	11
    ORDER_DOUBLE_ANGLE_V_GROOVE =	12
    ORDER_DOUBLE_ANGLE_CHAMFER =	13
    ORDER_VERTICAL_LATHE =		14
    ORDER_LAST =			15

    # These constants identify what order to do operations in:
    KIND_CONTOUR = 0
    KIND_DOWEL_PIN = 1
    KIND_DRILL = 2
    KIND_ROUND_POCKET = 3
    KIND_SIMPLE_EXTERIOR = 4
    KIND_SIMPLE_POCKET = 5
    KIND_VERTICAL_LATHE = 6

    POCKET_KIND_FLAT = 0
    POCKET_KIND_THROUGH = 1

    def __init__(self,
      name, kind, part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed):
	""" *Operation*: Initialize an *Operation* object to contain
	    *name*, *kind*, *part*, *comment*, *sub_priority*, *tool*,
	    *order*, *follows*.
	"""

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(kind, int)
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(order, int) and Operation.ORDER_NONE < order < Operation.ORDER_LAST
	assert follows is None or isinstance(follows, Operation)
	assert isinstance(feed_speed, Speed)
	assert isinstance(spindle_speed, Hertz)

	# Load up *self*:
	self._name = name
	self._part = part
	self._feed_speed = Speed()
	self._spindle_speed = Hertz()
	self._comment = comment
	self._follows = follows
	self._index = -1
	self._position = 0 # part.position_count
	self._order = order
	self._priority = part._priority_get()
	self._sub_priority = sub_priority
	self._tool = tool
	self._vice_x = part._vice_x_get()
	self._vice_y = part._vice_y_get()

    def _compare(self, operation2):
	""" *Operation*: Return -1, 0, or 1 depending upon whether *operation1*
	    (i.e. *self*) should sort before, at, or after *operation2*.
	"""

	# Verify argument types:
	assert isinstance(operation2, Operation)

	# Use *operation1* instead of *self*:
	operation1 = self

	# Sort by *priority* first:
	result = int_compare(operation1._priority, operation2._priority)
	if result != 0:
	    return result

	# Sort by *sub_priority* second:
	result = int_compare(operation1._sub_priority, operation2.sub_priority)
	if result != 0:
	    return result

	# Sort by *tool* third:
	result = operation1._tool.compare(operation2._tool)
	if result != 0:
	    return result

	# Sort by *order* fourth:
	result = int_compare(operation1._order, operation2._order)
	if result != 0:
	    return result

	# Sort by {kind} fifth:
	result = int_compare(opertaion1._kind, operation2._kind)
	if result != 0:
	    return result

	# Sort by {vice_x} sixth:
	result = in_compare(operation1._vice_x,	operation2._vice_x)
	if result != 0:
	    return result

	# Sort by {vice_y} seventh:
	result = in_compare(operation1._vice_y,	operation2._vice_y)
	if result != 0:
	    return result

	# If we reach this point, we need to compare at the sub-class
	# level.  Since this routine is called from the sub-class
	# *compare* method, the desired comparison occurs next.
	assert type(operation1) == type(operation2)

	return result

    def _comment_get(self):
	""" *Operation*: Return the comment field of the *Operation* object (i.e. *self*). """

	return self._comment

    def _feed_speed_get(self):
	""" *Operation*: Return the feed_speed field of the *Operation* object (i.e. *self*). """

	return self._feed_speed

    def _name_get(self):
	""" *Operation*: Return the name field of the *Operation* object (i.e. *self*). """

	return self._name

    def _part_get(self):
	""" *Operation*: Return the part field of the *Operation* object (i.e. *self*). """

	return self._part

    def _position_get(self):
	""" *Operation*: Return the position field of the *Operation* object (i.e. *self*). """

	return self._position

    def _priority_get(self):
	""" *Operation*: Return the priority field of the *Operation* object (i.e. *self*). """

	return self._priority

    def _spindle_speed_get(self):
	""" *Operation*: Return the spindle_speed field of the *Operation* object (i.e. *self*). """

	return self._spindle_speed

    def _tool_get(self):
	""" *Operation*: Return the tool field of the *Operation* object (i.e. *self*). """

	return self._tool

    def _vice_x_get(self):
	""" *Operation*: Return vice Z coordinate for the *Operation* object (i.e. *self*). """

	return self._vice_x

    def _vice_y_get(self):
	""" *Operation*: Return vice Z coordinate for the *Operation* object (i.e. *self*). """

	return self._vice_y

class Operation_Contour(Operation):
    """ *Operation_Contour* is a class that represents a contour operation.
    """

    # An {Operation_Contour} corresponds to a shape to be milled out.
    # Specifically a contour is an N-sided polygon in the X/Y plane
    # that is repesented by N points P1, ..., Pn.  (The polygon has
    # no intenral holes.) Associated with each point is a radius, Ri,
    # that specifies the rounding radius for each corner.  Ri is set
    # to 0 for no rounding.
    #
    # When we perform an exterior mill on a contour, it traces out a
    # path with rounded corners.  A good way of thinking about it is
    # to think of rolling a circle or radius r that has the same radius
    # as the mill bit and track the center point of the circle.  As it
    # is going along the straight parts of the contour, the circle center
    # rolls in a straight line as well; this line is offset from the
    # contours by r.  When the circle rolls around a corner with a radius
    # of zero, the circle center traces out an arc of radius r.  When
    # the circle rolls around a corner with a radius of Ri, the circle
    # center traces out an arc of radisu Ri + r.  Thus, the circle center
    # traces an arc at every corner, where each corner has an radius of
    # Ri + r.  (Corners with no rounding have Ri = 0.)
    #
    # Now to make things more complicated, we can adjust the mill path
    # both inward and outward.  There is a variable called {contour_offset},
    # which we will call o, and another variable called {finish_offset}
    # which we will call f.  The contour offset can be positive or negativej
    # and specifies an amount to grow or shrink the contour outline.
    # A positive contour offset (o > 0) causes the contour to expand
    # out by an amount o.  And the negative contour offset (o < 0) and
    # causes the contour to reduce by an amount o.  The finish offset
    # is like the contour offset, but it is always positive (f >=0.)
    # When f > 0, the contour is milled first using climb milling
    # at path that is an additional distince f from the final outline.  The
    # final pass is done with no distance additional f distance using climb
    # milling.  Since f and o behave similarly, there is no real need to
    # treat them differently; they can simply be added together to figure out
    # the desired exterior contour path.
    #
    # When we expand the contour out by o+f (o+f > 0), all of the straight
    # line paths are offset to the outside by o+f and all of the arcs have
    # a radius of Ri+o+f.  That is easy.  When o+f < 0, we have to be more
    # careful.  The straight line offset to the inside by o+f, but arc radius
    # never goes less than Ri.  Technically, the arc radius at each corner
    # is max(Ri, Ri+o+f).  Since max(Ri, Ri+o+f) works for all values of o+f,
    # that is what we use.  When we roll the circle with a radius of r around
    # the expanded for the expanded/contracted contour, the radius traced by
    # the circle center is max(r+Ri, r+Ri+f+o).
    #
    # Now comes the details of figuring all of this stuff out.  The two
    # line segements attached to a point Pi are Pi-1:Pi and Pi:Pi+1 and
    # have an angle a.  (Modular arithemetic by N is used for the indicies.)
    # Now we need to create a line that bisects this angle.  All arc centers
    # live on this line.  Let's start with the simple case where f+o=0.
    # If Ri=0, the center of the path arc, Ci, is the same as Pi and the
    # arc radius r.  Now if Ri>0, the center, Ci, is moved in such that
    # Ci has a perpendicular distance from both line segments that is Ri.
    # The law of sines can be used on the right triangle which as the
    # three angles a/2, 90-a/2, and 90.  Thus,
    # 
    #        r      |Pi:Ci|
    #     ------- = ------- = |Pi:Ci| = distance between Pi and Ci.
    #     sin a/2   sin 90
    #
    # We know r and a, so we can easily compute |Pi:Ci|.  Using basic
    # vector arithemetic, the X and Y coordinates of Ci = (Cx,Cy) can
    # be computed.
    #
    # Now we the contour offset and finish offset and replace r above
    # with max(Ri, Ri+f+o) and we can easily compute Ci for the cases
    # where f+o != 0.

    def __init__(self, part, comment, sub_priority, tool, order, follows,
      z_start, z_stop, corners, offset, effective_tool_radius, passes):
	""" *Operation_Contour*:
	"""

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(order, int)
	assert isinstance(follows, Operation)
	assert isinstance(z_start, L)
	assert isinstance(z_stop, L)
	assert isinstance(corners, list)
	assert isinstance(offset, L)
	assert isinstance(effective_tool_radius, L)
	assert isinstance(passes, int)
	for corner in corners:
	    assert isinstance(corner, Corner)

	# Initialize super class:
	Operation.__init__(self, "Contour", KIND_CONTOUR,
	  part, comment, sub_priority, tool, order, follows)

	# Load up the rest of *self*:
	self.z_start = z_start
	self.z_stop = z_stop
	self.corners = corners
	self.offset = offset
	self.effective_tool_radius = effective_tool_radius
	self.passes = passes

    def _cnc_generate(self):
	""" *Operation_Contour*: Ggenerate the CNC code for *self*.
	"""

	# Use *contour* instead of *self*:
	contour = self

	# Grab some values from *contour*:
	part = contour.part
	shop = part._shop
	code = shop.code
	tool = operation.tool

	# Record *comment* into *code*:
	comment = operation.comment
	code.line_comment(comment)

	# Grab some values from *tool*:
	s = tool.spindle
	f = tool.feed
	tool_diameter = tool.diameter
	z_feed = f / 2

	# Grap some more values from *contour*:
	z_start = contour.z_start
	z_stop = contour.z_stop
	z_depth = z_start - z_stop
	passes = contour.passes
	corners = contour.corners
	offset = contour.offset
	radius = contour.effective_tool_radius
	depth_per_pass = z_depth / float(passes)

	# Now generate the G-Code.  We visit each corner once, and the
	# first corner twice.  Thus, we need to loop through
	# {size} + 1 times:

	code._z_safe_assert("contour", comment)

	plunge_offset = tool_diameter
	#call d@("plunge_offset temporary set to zero\n")
	zero = L()
	plunge_offset = zero

	for index in range(passes):
	    code.line_comment("Pass {0} of {1}".format(index + 1, passes))

	    # Get cutter down to the correct depth:
	    z = z_start - depth_per_pass * float(index + 1)

	    #call d@(form@("tool=%v% plunge_offset=%i%\n\") %
	    #  f@(tool.name) / f@(plunge_offset))

	    code.contour(corners, plunge_offset, offset, radius, True, z, f, s)

	code.z_safe_retract(z_feed, s)

    def compare(self, contour2):
	""" *Operation_Contour* will return -1, 0, 1 depending upon whether
	    *self* should sort before, at, or after *contour2*.
	"""

	# Verify argument types:
	assert isinstance(contour2, Operation_Contour)

	# Use *contour1* instead of *self*:
	contour1 = self

	# Make sure that we do the super class compare first:
	result = XML.compare(contour1, contour2)
	if result != 0:
	    return result

	# Compare on *z_start* second:
	result = -contour1.z_start.compare(contour2.z_start)
	if result != 0:
	    return result

	# Compare on *z_stop* third:
	result = -contour1.z_stop.compare(contour2.z_stop)
	if result != 0:
	    return result

	# Compare length of *corners* fourth:
	corners1 = contour1.corners
	corners2 = contour2.corners
	corners_size = corners1.size
	result = int_compare(corners_size, corners2.size)
	if result != 0:
	    return result

	# Compare each *corner* in order fifth:
	for index in range(corners_size):
	    corner1 = corners1[index]
	    corner2 = corners2[index]
	    result = corner1.compare(corner2)
	    if result != 0:
		break

	return result

class Operation_Dowel_Pin(Operation):
    """ *Operation_Dowel_Pin* is a class that implements a dowel pin operation.
    """

    def __init__(self,
      part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed,
      diameter, edge_x, edge_y, original_x, original_y, original_z, plunge_x, plunge_y,
      tip_depth, z_stop):
	""" *Operation_Dowel_Pin*: Initialize an *Operation_Dowel_Pin* object (i.e. *self*)
	    to contain *part*, *comment*, *sub_priority*, *tool*, *order*, *follows*,
	    *feed_speed*, *spindle_speed*, *diameter*, *edge_x*, *edge_y*, *original_y*,
	    *original_z*, *plunge_x*, *plunge_y*, *tip_depth*, and *z_stop*.
	"""

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(sub_priority, int)
	assert follows is None or isinstance(follows, Operation)
	assert isinstance(feed_speed, Speed)
	assert isinstance(spindle_speed, Hertz)
	assert isinstance(diameter, L)
	assert isinstance(edge_x, L)
	assert isinstance(edge_y, L)
	assert isinstance(original_x, L)
	assert isinstance(original_y, L)
	assert isinstance(original_z, L)
	assert isinstance(plunge_x, L)
	assert isinstance(plunge_y, L)
	assert isinstance(tip_depth, L)
	assert isinstance(z_stop, L)

	# Initialize super class:
	Operation.__init__(self, "Dowel_Pin", Operation.KIND_DOWEL_PIN,
	  part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed)

	# Load up the rest of *self*:
	self._diameter = diameter
	self._edge_x = edge_x
	self._edge_y = edge_y
	self._original_x = original_x
	self._original_y = original_y
	self._original_z = original_z
	self._plunge_x = plunge_x
	self._plunge_y = plunge_y
	self._tip_depth = tip_depth
	self._z_stop = z_stop

    def _cnc_generate(self):
	""" *Operation_Dowel_Pin*: Generate the CNC G-code for a an
	    *Operation_Dowel_Pin* object (i.e. *self*.)
	"""

	# Use *dowel_pin* instead of *self*.
	dowel_pin = self

	# Grab the *part*, *shop*, *code*, *vice*, and *jaw_width* from *dowel_pin*:
	part = dowel_pin._part_get()
	shop = part._shop_get()
	code = shop._code_get()
	vice = shop._vice_get()
	jaw_width = vice._jaw_width_get()

	    
	# Grab some values out of *dowel_pin*:
	comment = dowel_pin._comment
	diameter = dowel_pin._diameter
	edge_x = dowel_pin._edge_x
	edge_y = dowel_pin._edge_y
	plunge_x = dowel_pin._plunge_x
	plunge_y = dowel_pin._plunge_y
	z_stop = dowel_pin._z_stop
	tip_depth = dowel_pin._tip_depth

	# Compute *radius* and *half_jaw_width*:
	radius = diameter / 2
	half_jaw_width = jaw_width / 2

	# Define the speed and feed for these operations:
	ipm10 = Speed(in_per_sec=10.0)
	rpm0 = Hertz()

	# Output the *dowel_pin* comment supplied by the user:
	code._line_comment(comment)

	# Output some dimenions for G-code debugging purposes:
	code._line_comment(
	  "{0} Initial Dimensions: {1} x {2} x {3}".format(part._name_get(),
	  part._dx_original_get(), part._dy_original_get(), part._dz_original_get()))

	#call d@(form@("part:%v% cnc_gen@Op_Dowel_Pin:%i% %i%\n\") %
	#  f@(part.name) % f@(plunge_x) / f@(plunge_y))

	# Rapid over to the plunge point:
	code._xy_rapid(plunge_x, plunge_y)

	# Now pause to let operator see if Z-safe is at the right height:
	code._command_begin()
	code._unsigned("M", 6)
	code._unsigned("T", 9)
	code._comment("Operator may check that Z-safe is correct")
	code._command_end()

	# Output some information about the *dowel_pin* for G-code debugging:
	code._line_comment(
	  "z_stop={0:i} tip_depth={1:i}".format(z_stop, tip_depth))

	# Move slowly down to *z_stop*:
	code._z_feed(ipm10, rpm0, z_stop, "dowel_pin")

	# Move slowly to (*edge_x*, *edge_y*).  This may cause the material in the
	# vice to slide over:
	code._xy_feed(ipm10, rpm0, edge_x, edge_y)

	# Now pause again, to let the operator move piece up against the
	# the dowel pin (if it is not already up against there):
	code._command_begin()
	code._unsigned("M", 6)
	code._unsigned("T", 2)
	code._comment("Operator should place part against dowel pin")
	code._command_end()

	# Slowy retract away from the part edge back to *plunge_x* and get
	# back up to Z safe:
	code._xy_feed(ipm10, rpm0, plunge_x, plunge_y)
	code._z_safe_retract(ipm10, rpm0)

    def compare(self, dowel_pin2):
	""" *Operation_Dowel_Pin*: Return -1, 0, or 1 depending upon whether
	    *dowel_pin2* (i.e. *self*) should sort before, at, or after
	    *dowel_pin2*.
	"""

	# Verify argument types:
	assert isinstance(dowel_pin, Operation)

	# Use *dowel_pin1* instead of *self*:
	dowel_pin1 = self

	# Perform superclass comparisson first:
	result = Operation.compare(dowel_pin1, dowel_pin2)
	if result != 0:
	    return result
	
	# It would be an error to have two dowel pin operations next to one another
	# in the operations list.  So rather than write some code to perform the
	# compare, we just assert *False*:
	assert False, "Two Dowel_Pin operations detected"

	return result

class Operation_Drill(Operation):
    """ *Operation_Drill* is a class that corresponds to a drill operation.
    """

    # Orginally copied from:
    #
    #     http://www.precisiontwistdrill.com/techhelp/article08281995.asp
    #
    # It is long gone...
    #
    # As manufacturing facilities move from conventional drilling methods to
    # CNC and "high performance", we shouldn't forget drilling basics.  They
    # still apply in the high tech drilling environment. This is a review of
    # some drilling basics:
    #
    #   * Use the shortest drill possible for the specific application.
    #     Longer drills are (1) more costly, (2) break easier and (3)
    #     drill bellmouthed holes.
    #   * Avoid the tendency to over speed and under feed.  Excessive
    #     speed causes (1) premature outer corner drill wear, (2) material
    #     work hardening, (3) long, stringy chips, (4) reduced drill
    #     life and (5) increased cost per hole.
    #   * Emphasizing feed rate: (1) helps break up chips (2) reduces
    #     premature outer corner drill wear, (3) reduces material work
    #     hardening, (4) extends drill life and (5) reduces cost per hole.
    #   * Use split point drills for drilling alloy materials; benefits
    #     include: (1) start at the point of contact (self-centering),
    #     (2) drill with less torque and thrust and (3) break up chips.
    #   * A hole of three drill diameters or deeper should be considered
    #     a deep hole. Therefore, you should peck drill just enough to
    #     prevent chips from packing in the flutes, because chip clogging
    #     is the major cause of drill breakage. You should decrease speeds
    #     and feeds as follows:
    #
    #    +----------------------------------------------------------+
    #    |     Speed and feed Reduction (Based upon hole depth)     |
    #    +-----------------------+-----------------+----------------+
    #    |     Hole Depth        |                 |                |
    #    |   to Diameter Ratio   |                 |                |
    #    |(times drill diameter) | Speed Reduction | Feed Reduction |
    #    +-----------------------+-----------------+----------------+
    #    |           3           |      10%        |      10%       |
    #    |           4           |      20%        |      10%       |
    #    |           5           |      30%        |      20%       |
    #    |           6           |     35-40%      |      20%       |
    #    +-----------------------+-----------------+----------------+
    #
    #   * Parabolic drills should be reduced as follows:
    #
    #    +-------------------------------------------------------------+
    #    | Speed Reduction -- Parabolic Drills (Based upon hole depth) |
    #    +------------------------------+------------------------------+
    #    | Hole Depth to Diameter Ratio |   Speed                      |
    #    | (times drill diameter)       | Reduction                    |
    #    +------------------------------+------------------------------+
    #    |            3                 |      0                       |
    #    |            4                 |      0                       |
    #    |            5                 |      5%                      |
    #    |         6 to 8               |      10%                     |
    #    |         8 to 11              |      20%                     |
    #    |        11 to 14              |      30%                     |
    #    |        14 to 17              |      50%                     |
    #    |        17 to 20              |      50%                     |
    #    +------------------------------+------------------------------+
    #
    #   * When drilling harder materials (i.e. above R/c 35); (1) reduce
    #     speeds and feeds to prevent points from burning and drilling
    #     breakage and (2) use cobalt drills as their higher hardness
    #     and heavy-duty construction are designed from drilling
    #     harder-materials
    #   * Use drills with a black oxide surface treatment when drilling
    #     ferrous materials.  The black oxide treatment holds the coolants
    #     and lubricants to the surface of the drill retarding material
    #     build-up. This treatment also improves toughness.

    def __init__(self, comment, sub_priority, tool, order, follows,
      diameter, hole_kind, x, y, z_start, z_stop, is_countersink):
	""" *Operation_Drill*: Initialize *Operation_Drill* to contain
	    *diameter*, *hole_kind*, *x*, *y*, *z_start*, *z_stop*, and
	    *counter_sink*.
	"""

	# Initialize super class:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(order, int)
	assert isinstance(follows, Operation)
	assert isinstance(diameter, L)
	assert isinstance(hole_kind, Hole_Kind)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(z_start, L)
	assert isinstance(z_stop, L)
	assert isinstance(is_countersink, bool)

	# Initialize the super class:
	Operation.__init__(self, "Drill", KIND_DRILL,
	  part, comment, sub_priority, tool, order, follows)

	# Load up *self*:
	self.diameter = diameter
	self.hole_kind = hole_kind
	self.x = x
	self.y = y
	self.z_start = z_start
	self.z_stop = z_stop
	self.is_countersink = is_countersink

    def _cnc_generate(self, is_last):
	""" *Operation_Drill*: Generate the CNC G-code for an
	    *Operation_Drill* object (i.e. *self*).
	"""
	# Verify argument types:
	assert isinstance(is_last, bool)

	# Use *drill* instead of *self*:
	drill = self

	part = drill.part
	shop = part._shop
	code = shop.code
	comment = operation.comment
	tool = operation.tool
	f = tool.feed
	s = tool.spindle

	diameter = drill.diameter
	x = drill.x
	y = drill.y
	z_start = drill.z_start
	z_stop = drill.z_stop
	radius = diameter / 2
	half_radius = radius / 2

	is_laser = tool.is_laser()

	#call d@(form@("Part=%v% Drill %i% hole at (%i%, %i%) is_laser=%l%\n\")%
	#  f@(part.name) % f@(diameter) % f@(x) % f@(y) / f@(is_laser))

	if is_laser:
	    code.dxf_circle(x, y, diameter)
	else:
	    hole_kind = drill_hole_kind
	    if hole_kind == HOLE_KIND_THROUGH:
		took_kind = tool.kind
		assert tool_kind == TOOL_KIND_DRILL
		tool_drill = tool.drill
		point_angle = tool_drill.point_angle
		tip_depth = tool.tip_depth(point_angle)
		#call line_comment@(code,
		#  form@("z_stop=%i% diameter=%i% tip_depth=%i%") %
		#  f@(z_stop) % f@(diameter) / f@(tip_depth))
		z_stop = z_stop - tip_depth - L(inch=0.040)
	    elif hole_kind == HOLE_KIND_TIP:
		pass
	    elif hole_kind == HOLE_KIND_FLAT:
		assert False, "Flat holes can't be done with point drills"
	    else:
		assert False, "Unknown hole kind"

	    code._z_safe_assert("drill", comment)

	    diameter_divisor = 3.0
	    material = part.material
	    #switch material.named_material
	    #  case plastic
	    #	# Plastic can really clog up the drill flutes, try to short the
	    #	# the drill spirals:
	    #	diameter_divisor := 0.5    

	    depth = z_start - z_stop
	    trip_depth = diameter / diameter_divisor
	    if depth > trip_depth:
		# "Deep" hole:

		# Output G73 G98 F# Q# R# S# X# Y# Z# (comment):
		# Output G83 G98 F# Q# R# S# X# Y# Z# (comment):
		# Output O900 [x] [y] [z_safe] [z_start] [z_stop] [z_step]
		#    [z_back] [feed] [speed] (comment):
		# Output O910 [x] [y] [z_safe] [z_start] [z_stop] [z_step]
		#    [z_back] [feed] [speed] (comment):

		z_rapid = part.z_rapid
		depth = z_rapid - z_stop
		pecks = int(depth / trip_depth) + 1
		# Add just a little to make sure {pecks} * {q} > {depth};
		# The drill will *never* go below the Z value:
		q = (depth / float(pecks)) + L(inch=.005)

		subroutine_code = 900
		if diameter <= L(inch="3/32"):
		    # Use full retract pecking:
		    subroutine_code = 910

		code.begin()
		code.subroutine_call(subroutine_code)
		code.length("[]", x - code.vice_x)
		code.length("[]", y - code.vice_y)
		code.length("[]", part.z_safe)
		code.length("[]", z_start)
		code.length("[]", z_stop)
		code.length("[]", radius)
		code.length("[]", half_radius)
		code.speed("[]", f)
		code.hertz("[]", s)
		code.comment(comment)
		code.end()

		#call begin@(code)
		#call mode_motion@(code, 73)
		#call mode_motion@(code, 83)
		#call mode_canned_cycle_return@(code, 98)
		#call speed@(code, "F", f)
		#call length@(code, "Q", code_length@(q))
		#call length@(code, "R1", z_rapid)
		#call hertz@(code, "S", s)
		#call length@(code, "X", code_length@(x))
		#call length@(code, "Y", code_length@(y))
		#call length@(code, "Z1", code_length@(z_stop))
		#call comment@(code, comment)
		#call end@(code)
	    else:	
		# Regular hole:

		# Output G81 G98 F# R# S# X# Y# Z# (comment):

		code.begin()
		code.mode_motion(81)
		code.mode_canned_cycle_return(98)
		code.speed("F", f)
		code.length("R1", part.z_rapid)
		code.hertz("S", s)
		code.length("X", x)
		code.length("Y", y)
		code.length("Z1", z_stop)
		code.comment(comment)
		code.end()

	    if is_last:
		# Output a G80 to exit canned cycle mode:
		code.begin()
		code.mode_motion(80)
		code.comment("End canned cycle")
		code.end()

		# Forget the G98 or G99 code:
		code.g2 = -1

	    # Do any dwell:
	    cnc_drill_count = part.cnc_drill_count + 1
	    part.cnc_drill_count = cnc_drill_count
	    if not drill.is_countersink and \
	      cnc_drill_count % part.cnc_drill_pause == 0:
		# Perform a Spindle Stop, Pause, Spindle_Start, Pause:
		code.begin()
		code.unsigned("M", 5)
		code.end()

		code.begin()
		code.unsigned("G11", 4)
		code.time("P", Code_Time(seconds=3.0))
		code.end()

		code.begin()
		code.unsigned("M", 3)
		code.hertz("S", s)
		code.end()

		code.begin()
		code.unsigned("G11", 4)
		code.time("P", Code_Time(seconds=8.0))
		code.end()

    def compare(self, drill2):
	""" *Operation_Drill*: Return -1, 0, 1 if *drill1* (i.e. *self*)
	    should sort before, at, or after *drill2*.
	"""

	# Verify argument types:
	assert isinstance(drill2, Operation)

	# Use *drill1* instead of *self*:
	drill1 = self

	# Perform superclass comparison first:
	result = Operation.compare(drill1, drill2)
	if result != 0:
	    return result

	# Perform *diameter* comparison second:
	result = drill1.diameter.compare(drill2.diameter)
	if result != 0:
	    return result

	# The remaining comparisons are for disambiguation only:
	result = drill1.x.compare(drill2.x)
	if result != 0:
	    return result
	result = drill1.y.compare(drill2.y)
	if result != 0:
	    return result
	result = -drill1.z_start.compare(drill2.z_start)
	if result != 0:
	    return result
	result = -drill1.z_stop.compare(drill2.z_stop)
	if result != 0:
	    return result

	return result

class Operation_Round_Pocket(Operation):
    """ *Operation_Round_Pocket* is a class that implements a round pocket
	manufacturing operation.
    """

    def __init_(self, part, comment, sub_priority, tool, order, follows,
      diameter, hole_kind, x, y, z_start, z_stop):
	""" *Operation_Round_Pocket* will intitialize an
	    *Operation_Round_Pocket* object (i.e. *self*) to contain
	    *part*, *comment*, *sub_priority*, *tool*, *order*, *follows*,
	    *diameter*, *hole_kind*, *x*, *y*, *z_start*, and *z_stop*.
	"""

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(order, int)
	assert isinstance(follows, Operation)
	assert isinstance(diameter, L)
	assert isinstance(hole_kind, Hole_Kind)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(z_start, L)
	assert isinstance(z_stop, L)

	# Initialize superclass:
	Operation.__init__(self, "Round_Pocket", KIND_ROUND_POCKET,
	  part, comment, sub_priority, tool, order, follows)

	# Load up *self*:
	self.diameter = diameter
	self.hole_kind = hole_kind
	self.x = x
	self.y = y
	self.z_start = z_start
	self.z_stop = z_stop

    def _cnc_generate(self):
	""" *Operation_Round_Pocket*: Generate the CNC G-code for an
	    *Operation_Round_Pocket* (i.e. *self*).
	"""

	# Extract some values from {round_pocket}:
	diameter = round_pocket.diameter
	x = round_pocket.x
	y = round_pocket.y
	z_start = round_pocket.z_start
	z_stop = round_pocket.z_stop
	assert z_start >= z_stop

	# Compute some values based on {diameter}:
	maximum_depth = diameter / 3.0
	radius = diameter / 2

	# Extract some values from {part} and {operation}.
	shop = part._shop
	code = shop.code
	tool = operation.tool

	# Extract some values from {tool}:
	tool_diameter = tool.diameter
	comment = operation.comment
	f = tool.feed
	s = tool.spindle
	hole_kind = round_pocket.hole_kind

	# Figure out if {tool} is a laser:
	is_laser = is_laser(tool)

	# Compute some values based on {tool_diameter}:
	tool_radius = tool_diameter / 2 
	half_tool_radius = tool_radius / 2

	#call d@(form@("cnc_gen_round_pocket_pock: part=%v% laser=%l%\n\") %
	#  f@(part.name) / f@(is_laser))
	if is_laser:
	    # We just cut a simple circle:
	    code.dxf_circle(x, y, radius - tool_radius)
	else:
	    # We do all the work to mill out the round_pocket pocket;

	    # Deal with through holes:
	    is_through = False
	    if hole_kind == HOLE_KIND_THROUGH:
		is_through = True
		z_stop = z_stop - L(inch=0.025)

	    code.line_comment(
	      "z_start={0} z_stop={1}".format(z_start, z_stop))
	    
	    z_depth = z_start - z_stop
	    passes = int(z_depth / maximum_depth) + 1
	    depth_per_pass = z_depth / float(passes)
	    assert passes < 100

	    #call line_comment@(code,
	    #  form@("z_depth=%i% passes=%i% depth_per_pass=%i%") %
	    #  f@(z_depth) % f@(passes) / f@(depth_per_pass))

	    # Move to position:
	    code.line_comment(comment)
	    code._z_safe_assert("round_pocket_pocket", comment)
	    code.xy_rapid(x, y)

	    z_feed = f / 4.0
	    shave = L(inch=0.005)
	    for depth_pass in range(passes):
		code.line_comment(
		  "{0} round_pocket pocket [pass {1} of {2}]".
		  format(comment, depth_pass + 1, passes))
		
		# Get to proper depth:
		#call line_comment@(code,
		#  read_only_copy@(form@("x=%i% y=%i% x_value=%i% y_value=%i%")%
		#    f@(x) % f@(y) % f@(x_value@(code)) / f@(y_value@(code))))
		code.xy_feed(f, s, x, y)
		z = z_start - depth_per_pass * float(depth_pass + 1)
		code.z_feed(z_feed, s, z, "round_pocket_pocket")

		radius_remove = radius - shave - tool_radius
		if is_through:
		    code.ccw_circle(radius_remove, f, s, x, y)
		else:
		    # We have to mow out all the invening space:
		    radius_passes = int(radius_remove /  half_tool_radius) + 1
		    pass_remove = radius_remove / float(radius_passes)

		    for radius_index in range(radius_passes):
			code.ccw_circle(
			  pass_remove * float(radius_index + 1), f, s, x, y)

	    # Do a "spring pass" to make everybody happy:
	    code.line_comment(code,
	      "{0} round_pocket pocket 'spring' pass".foramt(comment))
	    path_radius = radius - tool_radius
	    half_path_radius = path_radius / 2
	    code.xy_feed(f, s, x, y)
	    code.xy_ccw_feed(f, half_path_radius, s,
	      x + half_path_radius, y + half_path_radius)
	    code.xy_ccw_feed(f, half_path_radius, s, x, y + path_radius)
	    code.ccw_circle(radius - tool_radius, f, s, x, y)
	    code.xy_ccw_feed(f, half_path_radius, s,
	      x - half_path_radius, y + half_path_radius)
	    code.xy_ccw_feed(f, half_path_radius, s, x, y)
	    code.z_safe_retract(z_feed, s)

    def compare(self, pocket2):
	""" *Operation_Pocket*: Return -1, 0, 1 if *pocket1* (i.e. *self*)
	    should sort before, at, or after *pocket2*.
	"""

	# Verify argument types:
	assert isinstance(pocket2, Operation)

	# Use *pocket1* instead of *self*:
	pocket1 = self

	# Perform superclass comparison first:
	result = Operation.compare(pocket1, pocket2)
	if result != 0:
	    return result

	# We want deeper starts to show up later in sort:
	result = -hole1.z_start.compare(hole2.z_start)
	if result != 0:
	    return result

	# We want deeper stops to show up later is sort:
	result = -hole1.z_stop.compare(hole2.z_stop)
	if result != 0:
	    return result

	# The remaining tests are for disambiguation only:
	result = hole1.diameter.compare(hole2.diameter)
	if result != 0:
	    return result
	result = hole1.x.compare(hole2.x)
	if result != 0:
	    return result
	result = hole1.y.compare(hole2.y)
	if result != 0:
	    return result

	return result

class Operation_Simple_Exterior(Operation):
    """ *Operation_Simple_Exterior* is a class that implements a simple
	exterior contour manufacturing step.
    """

    def __init__(self, part, comment, sub_priority, tool, order, follows,
      x1, y1, x2, y2, z_start, z_stop, passes, corner_radius, tool_radius):
	""" *Operation_Simple_Exterior* will initialize an
	    *Operation_Simple_Exterior* object (i.e. *self*) with *part*,
	    *comment*, *sub_priority*, *tool*, *order*, *follows*, *x1*, *y1*,
	    *x2*, *y2*, *z_start*, *z_stop*, *passes*, *corner_radius*,
	    and *tool_radius*.
	"""

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(order, int)
	assert isinstance(follows, Operation)
	assert isinstance(x1, L)
	assert isinstance(y1, L)
	assert isinstance(x2, L)
	assert isinstance(y2, L)
	assert isinstance(z_start, L)
	assert isinstance(z_stop, L)
	assert isinstance(passes, int)
	assert isinstance(corner_radius, L)
	assert isinstance(tool_radius, L)

	# Initialize superclass:
	Operation.__init__(self, "Round_Pocket", KIND_ROUND_POCKET,
	  part, comment, sub_priority, tool, order, follows)

	# Load up rest of *self:
	self.x1 = x1
	self.y1 = y1
	self.x2 = x2
	self.y2 = y2
	self.z_start = z_start
	self.z_stop = z_stop
	self.passes = passes
	self.corner_radius = corner_radius
	self.tool_radius = tool_radius

    def _cnc_generate(self):
	""" *Operation_Simple_Exterior*: Generate the CNC G-code for an
	    *Operation_Simple_Exterior* object (i.e. self).
	"""

	shop = part._shop
	code = shop.code
	comment = operation.comment
	tool = operation.tool
	f = tool.feed
	s = tool.spindle
    
	corner_radius = exterior.corner_radius
	passes = exterior.passes
	tool_radius = exterior.tool_radius
	x1 = exterior.x1
	y1 = exterior.y1
	x2 = exterior.x2
	y2 = exterior.y2
	z_start = exterior.z_start
	z_stop = exterior.z_stop
	#call d@(form@("comment=%v% tool=%v% z_start=%i% z_stop=%i%\n\") %
	#  f@(comment) % f@(tool.name) % f@(z_start) / f@(z_stop))

	# Compute some more coordinates:
	x1mtr = x1 - tool_radius
	x1pcr = x1 + corner_radius
	x2mcr = x2 - corner_radius
	x2ptr = x2 + tool_radius
	y1mtr = y1 - tool_radius
	y1pcr = y1 + corner_radius
	y2mcr = y2 - corner_radius
	y2ptr = y2 + tool_radius

	r = tool_radius + corner_radius		# Total arc radius

	code._z_safe_assert("simple_exterior", comment)
	code.xy_rapid(x1mtr, y2ptr)		# Top Left corner

	# Output {comment}:
	z_delta = (z_start - z_stop) / float(passes)
	z_feed = f / 4.0
	index = 0
	while index < passes:
	    # Put put a reasonable comment:
	    trim_comment = \
	      "%s%: Trim to size [pass %d% of %d%]". \
	      format(comment, index + 1, passes)
	    code.line_comment(trim_comment)

	    # Bring tool down into material:
	    z = z_start - z_delta * float(index + 1)
	    code.z_feed(z_feed, s, z, "simple_exterior")

	    # Now trace out one contour path:
	    code.xy_feed(    f,    s, x2mcr, y2ptr)	# Top Right
	    code.xy_cw_feed( f, r, s, x2ptr, y2mcr)	# Right Top
	    code.xy_feed(    f,    s, x2ptr, y1pcr)	# Right Bottom
	    code.xy_cw_feed( f, r, s, x2mcr, y1mtr)	# Bottom Right
	    code.xy_feed(    f,    s, x1pcr, y1mtr)	# Bottom Left
	    code.xy_cw_feed( f, r, s, x1mtr, y1pcr)	# Left Bottom
	    code.xy_feed(    f,    s, x1mtr, y2mcr)	# Left Top
	    code.xy_cw_feed( f, r, s, x1pcr, y2ptr)	# Top Left
	    code.xy_feed(    f,    s, x1mtr, y2ptr)	# Top Left corner

	    index = index + 1

	# Get us back to Z safe:
	code.z_safe_retract(z_feed, s)

    def compare(self, exterior2):
	""" *Operation_Exterior*: Return -1, 0, 1 if *exterior1* (i.e. *self*)
	    should sort before, at, or after *exterior2*.
	"""

	# Verify argument types:
	assert isinstance(exterior2, Operation)

	# Use *exterior1* instead of *self*:
	exterior1 = self

	# Perform superclass comparison first:
	result = Operation.compare(exterior1, exterior2)
	if result != 0:
	    return result

	# Compare *z_start* second:
	result = -exterior1.z_start.compare(exterior2.z_start)
	if result != 0:
	    return result

	# Compare *z_stop* third:
	result = -exterior1.z_stop.compare(exterior2.z_stop)
	if result != 0:
	    return result

	# The remaining comparisons are for disambiguation only:
	result = exterior1.x1.compare(exterior2.x1)
	if result != 0:
	    return result
	result = exterior1.y1.compare(exterior2.y1)
	if result != 0:
	    return result
	result = exterior1.x2.compare(exterior2.x2)
	if result != 0:
	    return result
	result = exterior1.y2.compare(exterior2.y2)
	if result != 0:
	    return result
	result = exterior1.corner_radius.compare(exterior2.corner_radius)
	if result != 0:
	    return result
	result = exterior1.tool_radius.compare(exterior2.tool_radius)

	return result

class Operation_Simple_Pocket(Operation):
    """ *Operation_Simple_Pocket* is a class the performs a simple pocket
	manufacturing operation.
    """

    def __init__(self,
      part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed,
      x1, y1, x2, y2, z1, z2, corner_radius, tool_radius, pocket_kind):
	""" *Operation_Simple_Pocket*: Initialize an *Operation_Simple_Pocket*
	    object (i.e. *self*) to contain *part*, *comment*, *sub_priority*,
	    *tool*, *order*, *follows*, *x1*, *y1*, *x2*, *y2*, *z1*,
	    *z2*, *corner_radius*, *tool_radius*, and *pocket_kind*.
	"""

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(order, int)
	assert isinstance(follows, Operation) or follows == None
	assert isinstance(feed_speed, Speed)
	assert isinstance(spindle_speed, Hertz)
	assert isinstance(x1, L)
	assert isinstance(y1, L)
	assert isinstance(x2, L)
	assert isinstance(y2, L)
	assert isinstance(z2, L)
	assert isinstance(z1, L)
	assert isinstance(corner_radius, L)
	assert isinstance(tool_radius, L)
	assert isinstance(pocket_kind, int)
	assert z1 < z2, \
	  "z1 ({0:i}) is not less than z2 ({1:i})".format(z1, z2)

	# Initialize superclass:
	operation_kind = Operation.KIND_SIMPLE_POCKET
	Operation.__init__(self, "Simple_Pocket", operation_kind,
	  part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed)

	# Load up the rest of *self*:
	self._x1 = x1
	self._y1 = y1
	self._x2 = x2
	self._y2 = y2
	self._z_start = z2
	self._z_stop = z1
	self._corner_radius = corner_radius
	self._tool_radius = tool_radius
	self._pocket_kind = pocket_kind

    def _corner_radius_get(self):
	""" *Operation_Simple_Pocket: Return the corner radius for the *Operation_Simple_Pocket*
	    object (i.e. *self*).
	"""

	return self._corner_radius

    def _cnc_generate(self):
	""" *Operation_Simple_Pocket*: Generate the CNC G-code for a
	    *Operation_Simple_Pocket* object (i.e. *self*).
	"""

	# Use *pocket* instead of *self*:
	pocket = self

	# Grap some values from *pocket*:
	x1 = pocket._x1
	y1 = pocket._y1
	x2 = pocket._x2
	y2 = pocket._y2
	z_start = pocket._z_start
	z_stop = pocket._z_stop
	tool_radius = pocket._tool_radius
	corner_radius = pocket._corner_radius
	pocket_kind = pocket._pocket_kind
	comment = pocket._comment
	assert z_stop < z_start

	# Extract some values from {part} and {operation}
	part = pocket._part
	shop = part._shop_get()
	code = shop._code_get()
	tool = pocket._tool
	f = tool._feed_speed_get()
	s = tool._spindle_speed_get()

	# Start with {comment}:
	code._line_comment(comment)

	# Output the pocket operations:
	code._line_comment(
	  "[{0:i},{1:i}]:[{2:i},{3:i}][{4:i}:{5:i}] tool_radius={6:i} corner_radius={7:i}".
	  format(x1, y1, x2, y2, z_stop, z_start, tool_radius, corner_radius))

	# Define some constants:
	zero = L()
	z_extra = zero

	# Figure out the minimum pocket span:
	dx = x2 - x1
	dy = y2 - y1
	minimum_span = zero
	if dx > dy:
	    minimum_span = dy
	else:
	    minimum_span = dx
	half_minimum = minimum_span / 2

	# Emit the preparatory G-Code:

	is_laser = tool._is_laser_get()

	# Compute the number of depth passes required:
	maximum_operation_depth = tool_radius / 2
	if is_laser:
	    maximum_operation_depth = tool.maximum_z_depth

	total_cut = z_start - z_stop + z_extra
	print("z_start={0:i} z_stop={1:i} z_extra={2:i}".format(z_start, z_stop, z_extra))
	passes = int(total_cut / maximum_operation_depth)
	while maximum_operation_depth* float(passes) < total_cut:
	    passes = passes + 1
	assert passes > 0
	assert maximum_operation_depth * float(passes) >= total_cut
	z_step = total_cut / float(passes)

	# *tool_radius* is the radius of the selected end-mill,
	# which we will henceforth, call R. We need to be careful on
	# our inner passes because the end-mill can only cut R/sqrt(2)
	# along the diagnoal.  If each pass is separated by exactly
	# 2R, there will be little triangular islands left behind.
	#
	# We now define T=R/sqrt(2).  To avoid the little island
	# problems, each inner pass must overlap the next pass out
	# by at least (R-T).  Diagramatically:
	#
	# Pass 3 |-------|-------|
	#           R        R
	#                        |  T  |   R   |
	# Pass 2               |-------|-------|
	#                          R       R
	#                                      |  T  |   R   |
	# Pass 1                             |-------|-------|
	#                                        R       R
	#
	# where the over lap between two passes is R-T.  In the example
	# above the offset from the center for pass 1 is R, pass 2 is
	# R+T+R, and for pass 3 is R+T+R+T+R.  Of course, the center to
	# edge distance may not be 2R+N*(R+T), for some value of N.
	# The better way to organize this is to offset the end-mill
	# from the edge by R on the pass N, R+T+R on the pass N-1,
	# and R+T+R+T+R on pass N-2, etc.

	r = tool_radius
	# 0.65 < 0.707 = 1/sqrt(2)
	#t :@= smul@(r, 0.65)
	t = r / 2

	code._z_safe_assert("simple_pocket", comment)

	# Compute the total number of rectangular paths needed:
	#call d@(form@("comment=%v% r=%i% t=%i%\n\") %
	#   f@(comment) % f@(r) / f@(t))
	paths = 0
	remaining = half_minimum
	while remaining > zero:
	    #call d@(form@("comment=%v% paths=%d% remaining=%i%\n\") %
	    #  f@(comment) % f@(paths) / f@(remaining))

	    if paths == 0:
		# {r + r} does not work if pocket width equals tool width,
		# so we just use {r + t} every time now:
		#remaining := remaining - (r + r)
		remaining = remaining - (r + t)
	    else:
		remaining = remaining - (r + t)
	    paths += 1
	#call d@(form@("comment=%v% paths=%d% remaining=%i% (final)\n\") %
	#  f@(comment) % f@(paths) / f@(remaining))

	# Generate {passes} deep depth passes over the pocket:
	for depth_pass in range(passes):
	    # Output "(Depth Pass # of #)" comment:
	    code._line_comment("Depth Pass {0} of {1}".
	      format(depth_pass + 1, passes))

	    # Compute the plunge depth for each pass:
	    z_plunge = z_start - z_step * float(depth_pass + 1)

	    if pocket_kind == Operation.POCKET_KIND_THROUGH:
		# We only need to do the exterior path to a depth of {z_plunge}:
		code.simple_pocket_helper(pocket, r, s, f, z_plunge, depth_pass != 0)

	    elif pocket_kind == Operation.POCKET_KIND_FLAT:
		# Generate {paths} rectangular passes over the pocket:
		for path in range(paths):
		    # Provide a context comment:
		    code._line_comment("Rectangular Path {0} of {1}".format(path + 1, paths))

		    # Compute the {offset} for this path:
		    offset = r + (r + t) * float((paths - 1) - path)

		    # If {offset} exceeds {half_minimum} it will cause
		    # {px1} > {px2} or {py1} > {py2}, which is bad.
		    # We solve the problem by not  letting {offset}
		    # exceed {half_minimum}:
		    if offset > half_minimum:
			offset = half_minimum

		    # We need to get the tool to the starting point safely:
		    rapid_move = False
		    if path == 0:
			# This is the first cut at a new depth:
			if depth_pass == 0:
			    # We are at {z_safe}, so we can move rapidly to the
			    # plunge point:
			    rapid_move = True
			else:
			    # We are sitting at the outer edge at the old depth,
			    # and there is no material between where are now
			    # and ({start_x}, {start_y}); a rapid could be too
			    # fast, so we do a linear move:
			    rapid_move = False
		    else:
			# We are at the previous internal path at this depth
			# and need to move out and cut material as we go:
			rapid_move = False

		    code._simple_pocket_helper(pocket, offset, s, f, z_plunge, rapid_move)
	    else:
		assert False, "Unknown pocket kind: {0}".format(pocket_kind)

	# Return the tool to a safe location above the material:
	code._z_safe_retract(f, s)
	code._line_comment("Simple Pocket Done")

    @staticmethod
    def _compare(simple_pocket, simple_pocket2):
	""" *Operation_Simple_Pocket* returns -1, 0, or 1 if *simple_pocket1*
	    (i.e. *simple_pocket1*) should sort before, at or after *simple_pocket*.
	"""

	# Verify argument types:
	assert isinstance(simple_pocket2, Operation)

	# Use *simple_pocket1* instead of *self*:
	simple_pocket1 = self

	# Perform operation level compare first:
	result = Operation.compare(simple_pocket1, simple_pocket2)
	if result != 0:
	    return result

	# Compare *z_start* first:
	result = -pocket1._z_start.compare(pocket2._z_start)
	if result != 0:
	    return result

	# Compare *z_stop* second:
	result = -pocket1._z_stop.compare(pocket2._z_stop)
	if result != 0:
	    return result

	# The remaining comparisons are really only for disambiguation:
	result = pocket1._x1.compare(pocket2._x1)
	if result != 0:
	    return result
	result = pocket1._y1.compare(pocket2._y1)
	if result != 0:
	    return result
	result = pocket1._x2.compare(pocket2._x2)
	if result != 0:
	    return result
	result = pocket1._y2.compare(pocket2._y2)
	if result != 0:
	    return result
	result = pocket1._corner_radius.compare(pocket2._corner_radius)
	if result != 0:
	    return result
	result = pocket1._tool_radius.compare(pocket2._tool_radius)
	if result != 0:
	    return result
	result = pocket1._pocket_kind - pocket2.__pocket_kind
	if result != 0:
	    return result

	return result

    def _x1_get(self):
	""" *Operation_Simple_Pocket*: Return the x1 field of the *Operation_Simple_Pocket*
	    (i.e. *self*). """

	return self._x1

    def _x2_get(self):
	""" *Operation_Simple_Pocket*: Return the x2 field of the *Operation_Simple_Pocket*
	    (i.e. *self*). """

	return self._x2

    def _y1_get(self):
	""" *Operation_Simple_Pocket*: Return the y1 field of the *Operation_Simple_Pocket*
	    (i.e. *self*). """

	return self._y1

    def _y2_get(self):
	""" *Operation_Simple_Pocket*: Return the y2 field of the *Operation_Simple_Pocket*
	    (i.e. *self*). """

	return self._y2

class Operation_Vertical_Lathe(Operation):
    """ *Operation_Vertical_Lathe* is a class that represents a vertical
	lathe manufacturing operation.
    """

    def __init__(self, part, comment, sub_priority, tool, order, follows,
      x, y, inside_diameter, outside_diameter, z_start, z_stop):
	""" *Operation_Veritcal_Lathe*: Initialize an *Operation_Vertical_Lathe*
	    object (i.e. *self*) to contain *part*, *comment*, *sub_priority*,
	    *tool*, *order*, *follows*, *x*, *y*, *inside_diameter*,
	    *outside_diameter*, *z_start*, and *z_stop*.
	"""

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(order, int)
	assert isinstance(follows, Operation)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(inside_diameter, L)
	assert isinstance(outside_diameter, L)
	assert isinstance(z_start, L)
	assert isinstance(z_stop, L)

	# Initialize superclass:
	Operation.__init__(self, "Vertical_Lathe", KIND_VERTICAL_LATHE,
	  part, comment, sub_priority, tool, order, follows)

	# Load up the rest of *self*:
	self.x = x
	self.y = y
	self.inside_diameter = inside_diameter
	self.outside_diameter = outside_diameter
	self.z_start = z_start
	self.z_stop = z_stop

    def _cnc_generate(self):
	""" *Operation_Vertical_Lathe*:
	"""

	shop = part._shop
	code = shop.code
	zero = L()

	# Grap some values from {lathe}:
	x = lathe.x
	y = lathe.y
	inside_diameter = lathe.inside_diameter
	outside_diameter = lathe.outside_diameter
	z_start = lathe.z_start
	z_stop = lathe.z_stop
	z_extra = zero

	# Grab some values from {operation}:
	comment = operation.comment
	tool = operation.tool
	f = tool.feed
	s = tool.spindle

	# Grab some values from {tool}:
	tool_diameter = tool.diameter
	tool_radius = tool_diameter / 2

	# Deal with {outside_diameter} < {inside_diameter}:
	if outside_diameter < inside_diameter + tool_diameter * 2:
	    outside_diameter = inside_diameter + tool_diameter * 2
	#call d@(form@("outside_diameter=%i% inside_diameter=%i%\n\") %
	#  f@(outside_diameter) / f@(inside_diameter))
	inside_radius = inside_diameter / 2
	outside_radius = outside_diameter / 2

	#call d@(form@("outside_diameter=%i% inside_diameter=%i%\n\") %
	#  f@(outside_diameter) / f@(inside_diameter))
	#call d@(form@("outside_radius=%i% inside_radius=%i%\n\") %
	#  f@(outside_radius) / f@(inside_radius))

	# {epsilon} is small enough to make nudge {gap_passes} down when
	# {outside_diameter} = {inside_diameter} + {twice@(tool_diameter)}
	epsilon = L(inch=0.000000001)

	# Compute {gap_step} from {inside_radius}, {outside_radius},
	# and {tool_diameter}:
	gap_total = outside_radius - inside_radius
	gap_passes = int((gap_total - epsilon) / tool_diameter) + 1
	gap_step = zero
	if gap_passes > 1:
	    gap_step = (gap_total - tool_diameter) / float(gap_passes - 1)
	assert gap_passes < 10000

	#call d@(form@("vertical_lathe: x=%i% y=%i% odiam=%i% idiam=%i%\n\") %
	#  f@(x) % f@(y) % f@(outside_diameter) / f@(inside_diameter))
	#call d@(form@(
	#  "vertical_lathe:gap_tot=%i% gap_pass=%d% gap_step=%i%\n\")%
	#  f@(gap_total) % f@(gap_passes) / f@(gap_step))

	# Compute {depth_step} from {z_start}, {z_stop} and {tool_diameter}:
	depth_total = z_start - z_stop
	depth_passes = int(depth_total / (tool_diameter / 4.0)) + 1
	depth_step = depth_total / float(depth_passes)
	assert depth_passes < 10000

	# Start code generation with {comment}:
	code.line_comment(comment)

	# Get the tool placed at the right place to {plunge}:
	plunge_x = x
	plunge_y = y + outside_radius - tool_radius
	code._z_safe_assert("vertical_lathe", comment)
	code.xy_rapid(plunge_x, plunge_y)

	#call d@(form@("vertical_lathe: plunge_x=%i% plunge_y=%i%\n\") %
	#  f@(plunge_x) / f@(plunge_y))
	#call d@(form@(
	#  "vertical_lathe: gapall=%i% gappasses=%d% gapstep=%i%\n\") %
	#  f@(gap_total) % f@(gap_passes) / f@(gap_step))

	#call d@(form@("tool_diameter=%i% tool_radius=%i%\n\") %
	#  f@(tool_diameter) / f@(tool_radius))

	# Now iterate down to the desired depth:
	depth_index = 0
	while depth_index < depth_passes:
	    # Compute {z_depth} from {depth_step} and {depth_index}:
	    z_depth = z_start - depth_step * float(depth_index + 1)

	    # Make sure we are in the right place to start:
	    code.xy_feed(f, s, plunge_x, plunge_y)

	    # Get to proper {z_depth}:
	    code.line_comment("{0}: Depth Pass {1} of {2}".
	      format(comment, depth_index + 1, depth_passes))
	    code.z_feed(f/ 4.0, s, z_depth, "vertical lathe")

	    # Iterate across the gap:
	    gap_index = 0
	    while gap_index < gap_passes:
		# Compute radius of this pass:
		gap_radius = \
		  outside_radius - tool_radius - gap_step * float(gap_index)

		# Peform the circle operation:
		code.line_comment("{0}: Gap pass {1} of {2}".
		  format(comment, gap_index + 1,gap_passes))
		code.ccw_circle(gap_radius, f, s, x, y)

		#call d@(form@("vertical_lathe[%d%,%d%]: gap_radius=%i%\n\") %
		#  f@(depth_index) % f@(gap_index) / f@(gap_radius))

		gap_index += 1
	    depth_index += 1
    
	# Return the tool to a safe location above the material:
	code.z_safe_retract(f, s)

	code.line_comment("Vertical Lathe Done")

    def compare(self, lathe2):
	""" *Operation_Vertical_Lathe* will return -1, 0, or 1 depending
	    upon whether *lathe1* (i.e. *self*) should sort
	    before, at, or after *lathe2*.
	"""

	# Verify argument types:
	assert isinstance(vertical_lathe2, Operation)

	# Use *lathe1* instead of *self*:
	lathe1 = self

	# Do *Operation* level compare first:
	result = Operation.compare(lathe1, lathe2)
	if result != 0:
	    return result

	# Compare against *z_start* second:
	result = lathe1.z_start.compare(lathe2.z_start)
	if result != 0:
	    return result

	# Compare against *z_stop* third:
	result = lathe1.z_stop.compare(lathe2.z_stop)
	if result != 0:
	    return result

	# The remaining comparisons are for disambiguation only:
	result = lathe1.x.compare(lathe2.x)
	if result != 0:
	    return result
	result = lathe1.y.compare(lathe2.y)
	if result != 0:
	    return result
	result = lathe1.inside_diameter.compare(lathe2.inside_diameter)
	if result != 0:
	    return result
	result = lathe1.outside_diameter.compare(lathe2.outside_diameter)

	return result

class Part:
    """ A {Part} specifies either an assembly of parts or a single
	physical part. """

    # Flavors of values that can be stored in a {Part}:
    def __init__(self, up, name = None):
	""" *Part*: Initialize *self* to have a parent of *up*. """

        # Check argument types:
	assert up == None or isinstance(up, Part)
	assert name == None or isinstance(name, str)

	# Some useful abbreviations:
	zero = L()

	ezcad = EZCAD3.ezcad
	self._bounding_box = Bounding_Box()
	self._color = None
	self._dx_original = zero
	self._dxf_x_offset = zero
	self._dxf_y_offset = zero
	self._dy_original = zero
	self._dz_original = zero
	self._dxf_scad_lines = []
	self._ezcad = ezcad
	self._is_part = False	
	self._material = Material("aluminum", "")
	self._name = name
	self._operations = []
	self._plunge_x = zero
	self._plunge_y = zero
	self._position = Matrix()
	self._position_count = 0
	self._priority = 0
	#self._places = {}
	self._reposition = Matrix()
	self._rotate = None
	self._shop = ezcad._shop
	self._scad_difference_lines = []
	self._scad_union_lines = []
	self._signature_hash = None
	self._top_surface_set = False
	self._tool_preferred = ""
	self._translate = None
	self._visible = True
	self._vice_x = L()
	self._vice_y = L()
	self._z_floor = L(inch=0.5)
	self._z_rapid = L(inch=0.5)
	self._z_safe = L(inch=0.5)
	self.up = up

	# Make sure directory exists:
	directory = ezcad._directory_get()
	if not os.path.exists(directory):
	    os.makdirs(directory)

	#print("<=Part.__init__(*, '{0}', *, place={1})".format(name, place))

    def ___bounding_box_check(self, indent):
	box_points_list = self._box_points_list
	box_points_list_size = len(box_points_list)
	assert box_points_list_size > 0
	current_box_points = box_points_list[-1]
	previous_box_points = current_box_points
	if box_points_list_size > 1:
	    previous_box_points = box_points_list[-2]
	print("{0}{1}: {2} {3}".format(indent * " ", self._name,
	  len(previous_box_points), len(previous_box_points)))

	for attribute_name in dir(self):
	    if not attribute_name.startswith("_") and \
	      attribute_name.endswith("_"):
		sub_part = getattr(self, attribute_name)
		assert isinstance(sub_part, Part)
		sub_part._bounding_box_check(indent + 1)

    def xxxxxx_bounding_box_set(self, bounding_box):
	""" *Part*: For the *Part* object (i.e. *self*), set the bounding
	    box to *bounding_box*.
	"""

	# Verify argument types:
	assert isinstance(bounding_box, Bounding_Box)
	assert not bounding_box.is_empty()

	# Grab *tne* and *bsw* from *bounding_box*:
	tne = bounding_box._tne
	bsw = bounding_box._bsw

	# Grab the 6 X/Y/Z coordinates:
	self.tz = tz = tne.z
	self.ny = ny = tne.y
	self.ex = ex = tne.x
	self.bz = bz = bsw.z
	self.sy = sy = bsw.y
	self.wx = wx = bsw.x

	# Compute the center coordinates:
	self.cx = cx = (ex + wx) / 2
	self.cy = cy = (ny + sy) / 2
	self.cz = cz = (tz + bz) / 2

	# Fill in the slices:
	#      tnw--tn---tne     nw----n----ne     bnw--bn---bne
	#      |     |     |     |     |     |     |     |     |
	#      tw----t----te     w-----c-----e     bw----b----be
	#      |     |     |     |     |     |     |     |     |
	#      tsw--ts---tse     sw----s----se     bsw--bs---bse
	#       (Top Slice)      (Center Slice)    (Bottom Slice)

	# Top Slice:
	#      tnw--tn---tne
	#      |     |     |
	#      tw----t----te
	#      |     |     |
	#      tsw--ts---tse
	self.tnw = P(wx, ny, tz)
	self.tn =  P(cx, ny, tz)
	self.tne = P(ex, ny, tz)
	self.tw =  P(wx, cy, tz)
	self.t  =  P(cx, cy, tz)
	self.te  = P(ex, cy, tz)
	self.tsw = P(wx, sy, tz)
	self.ts =  P(cx, sy, tz)
	self.tse = P(ex, sy, tz)

	# Center slice:
	#      nw----n----ne
	#      |     |     |
	#      w-----c-----e
	#      |     |     |
	#      sw----s----se
	#      (Center Slice)
	self.nw = P(wx, ny, cz)
	self.n =  P(cx, ny, cz)
	self.ne = P(ex, ny, cz)
	self.w =  P(wx, cy, cz)
	self.c  = P(cx, cy, cz)
	self.e  = P(ex, cy, cz)
	self.sw = P(wx, sy, cz)
	self.s =  P(cx, sy, cz)
	self.se = P(ex, sy, cz)

	# Bottom slice:
	#      bnw--bn---bne
	#      |     |     |
	#      bw----b----be
	#      |     |     |
	#      bsw--bs---bse
	#      (Bottom Slice)
	self.bnw = P(wx, ny, bz)
	self.bn =  P(cx, ny, bz)
	self.bne = P(ex, ny, bz)
	self.bw =  P(wx, cy, bz)
	self.b  =  P(cx, cy, bz)
	self.be  = P(ex, cy, bz)
	self.bsw = P(wx, sy, bz)
	self.bs =  P(cx, sy, bz)
	self.bse = P(ex, sy, bz)

	# Save the *boundig_box*:
	self._bounding_box = bounding_box

    def xxxxxxx___bounding_box_update(self, bounding_box, comment, place,
      trace = -1000000):
	""" *Part*: Updated *self* with *bounding_box*, *place* and
	    *comment*. """

	if trace >= 0:
	    print("{0}=>Part._bounding_box_update(comment={1})".
	      format(trace * ' ', comment))

	# Check argument types:
	assert isinstance(bounding_box, Bounding_Box)
	assert isinstance(comment, str)
	assert isinstance(place, Place)

	# Grab some values from *bounding_box*:
	ex = bounding_box.ex
	wx = bounding_box.wx
	ny = bounding_box.ny
	sy = bounding_box.sy
	tz = bounding_box.tz
	bz = bounding_box.bz

	# Grab *forward_matrix* from *place*:
	forward_matrix = place._forward_matrix

	# Compute the bounding box points before rotation and translation:
	tne = P(ex, ny, tz)
	tnw = P(wx, ny, tz)
	tse = P(ex, sy, tz)
	tsw = P(wx, sy, tz)
	bne = P(ex, ny, bz)
	bnw = P(wx, ny, bz)
	bse = P(ex, sy, bz)
	bsw = P(wx, sy, bz)

	# Compute bounding box after rotation and translation:
	self._box_point_update(comment + "[TNE]",
	  forward_matrix.point_multiply(tne), trace = trace + 1)
	self._box_point_update(comment + "[TNW]",
	  forward_matrix.point_multiply(tnw), trace = trace + 1)
	self._box_point_update(comment + "[TSE]",
	  forward_matrix.point_multiply(tse), trace = trace + 1)
	self._box_point_update(comment + "[TSW]",
	  forward_matrix.point_multiply(tsw), trace = trace + 1)
	self._box_point_update(comment + "[BNE]",
	  forward_matrix.point_multiply(bne), trace = trace + 1)
	self._box_point_update(comment + "[BNW]",
	  forward_matrix.point_multiply(bnw), trace = trace + 1)
	self._box_point_update(comment + "[BSE]",
	  forward_matrix.point_multiply(bse), trace = trace + 1)
	self._box_point_update(comment + "[BSW]",
	  forward_matrix.point_multiply(bsw), trace = trace + 1)

	if trace >= 0:
	    print("{0}<=Part._bounding_box_update(comment={1})".
	      format(trace * ' ', comment))

    ## @brief Updates *name*'d *point* in *self* for bounding box calcuation.
    #  @param self is the *Part* to update.
    #  @param name is the name of the *P* object.
    #  @param point is the *P* to update.
    #
    # <I>_box_point_update</I>() will update the *P* named *name* in *self*
    # (a *Part* object) to be *point*.  If the value of *point* has changed
    # from the last time it was updated, the bounding box for *self* is
    # recomputed.
    def __box_point_update(self, name, point, trace = -1000000):
	""" *Part*: Insert/update the point named *name* to *point*. """

	# Check argument types:
	assert isinstance(name, str)
	assert isinstance(point, P)

	# Set *trace* to *True* if tracing is desired:
	if trace >= 0:
	    print("{0}=>Part._box_point_update('{1}', {2})".
	      format(trace * ' ', name, point))

	# Get *current_box_points* from *box_points_list*:
	box_points_list = self._box_points_list
	ezcad = self._ezcad
	update_count = ezcad._update_count
	while update_count >= len(box_points_list):
	    box_points_list.append([])
	box_points_list_size = len(box_points_list)
	current_box_points = box_points_list[update_count]

	# Take the current *point* onto the end of *current_box_points*:
	box_points_index = len(current_box_points)
	current_box_points.append(point)

	#FIXME: For debugging only!!!
	point._name = name

	# See if *point* changed from the last time:
	if update_count == 0:
	    # First time, just force a box update:
	    self._box_recompute("Part._box_point_update(): first point")
	    if trace >= 0:
		print("{0}  Part._box_update:({1}):first".
		  format(trace * ' ', point))
	else:
	    # Fetch *previous_point*:
	    previous_box_points = box_points_list[update_count - 1]
	    previous_box_points_size = len(previous_box_points)
	    if box_points_index < previous_box_points_size:
		previous_point = previous_box_points[box_points_index]

		# Force box update if the point changed:
		if previous_point != point:
		    self._box_recompute("Part._box_point_update(): changed")
		    if trace >= 0:
			print("{0}  Part._box_update:({1}):changed from {2}".
			  format(trace * ' ', point, previous_point))
		elif previous_point._name != point._name:
                    #FIXME: This test is for debugging only!!!
		    distance = previous_point.distance(point)
		    if distance >= L(inch = .00001):
			print("Part._box_point_update('{0}',{1}):skew:{0}!={1}".
			  format(name, point,
			  previous_point._name, point._name))
	    else:
		# How did this happen?:
		print("Part._box_point_update('{0}', {1}): skew {2}<{3}: {4}".
		  format(name, point, box_points_index,
		  previous_box_points_size,
		  previous_box_points[previous_box_points_size - 1]._name))

	# We only need to hang on to the *previous_box_points*; remove
	# any previous one to save a little memory and catch accessing
	# any stale information:
	if box_points_index == 0 & update_count >= 2:
	    box_points_list[update_count - 2] = None

	if trace >= 0:
	    print("{0}<=Part._box_point_update('{1}', {2})".
	      format(trace * ' ', name, point))


    ## @brief Recomputes the bounding box for *self* and any parent *Box*'s.
    #  @param *self* the bounding *Box* to recompute.
    #
    # <I>_recompute</I>() will recompute the bounding box for *self* (a bounding
    # *Box*) and enclosing parent bounding *Box*'s.
    def __box_recompute(self, label):
	""" *Part*: (Internal use only) Recompute the bounding box corners. """

	ezcad = self._ezcad
	update_count = ezcad._update_count

	assert isinstance(label, str)
	# For debugging:
	trace = False
	#trace = True
	if trace:
	    print("=>Part._box_recompute({0}, {1})".format(self._name, label))

	# Initialize the bounding box values with bogus big positive/negative
	# values:
	big = L(mm=123456789.0)
	ex = -big
	wx = big
	ny = -big
	sy = big
	tz = -big
	bz = big

	# Sweep through each *point*:
	ezcad = self._ezcad
	update_count = ezcad._update_count
	box_points_list = self._box_points_list
	box_points = []
	if update_count < len(box_points_list):
	    box_points = box_points_list[update_count]
	for point in box_points:
	    # Update *wx* and *ex* with *x* as appropriate.
	    x = point.x
            if x > ex:
		ex = x
	    if x < wx:
		wx = x
	    # Note that we can not assume that *ex* > *wx*
	    # because the initial value for *ex* is -*big* and the
	    # initial value for *wx* is *big*.  Thus, we use two
	    # sequential **if** statements,  rather than a single
	    # **if**-**elif** combination:

	    # Update *ny* and *sy* with *y* as appropriate.
            y = point.y
	    if y > ny:
		ny = y
	    if y < sy:
		sy = y

	    # Update *tz* and *bz* with *z* as appropriate.
            z = point.z
	    if z > tz:
		tz = z
	    if z < bz:
		bz = z

	# Determine if the bounding box changed:
	if self.ex != ex or self.wx != wx or \
	  self.ny != ny or self.sy != sy or \
	  self.tz != tz or self.bz != bz:
	    if update_count > 10 or trace:
		if self.ex != ex:
		    print("{0}:ex:{1} => {2}".format(self._name, self.ex, ex))
		if self.wx != wx:
		    print("{0}:wx:{1} => {2}".format(self._name, self.wx, wx))
		if self.ny != ny:
		    print("{0}:ny:{1} => {2}".format(self._name, self.ny, ny))
		if self.sy != sy:
		    print("{0}:sy:{1} => {2}".format(self._name, self.sy, sy))
		if self.tz != tz:
		    print("{0}:tz:{1} => {2}".format(self._name, self.tz, tz))
		if self.bz != bz:
		    print("{0}:bz:{1} => {2}".format(self._name, self.bz, bz))

	    # Bounding box changed:
	    self.ex = ex
	    self.wx = wx
	    self.ny = ny
	    self.sy = sy
	    self.tz = tz
	    self.bz = bz

	    # Update the *dx*/*dy*/*dz*:
	    self.dx = ex - wx
	    self.dy = ny - sy
	    self.dz = tz - bz

            # Keep track if we have changed:
	    #print("Part._box_recompute('{0}'): changed".format(label))
	    self._box_changed_count += 1

	    # Compute averate X/Y/Z:
	    cx = (ex + wx) / 2.0
	    cy = (ny + sy) / 2.0
	    cz = (tz + bz) / 2.0

	    # There are 27 points to compute -- 9 on the top slice, 9 on
	    # the middle slice, and 9 on the bottom slice:

	    # Top slice:
	    #	TNW	TN	TNE
	    #	TW	T	TE
	    #	TSW	TS	TSE
	    self.t =   P(cx, cy, tz)
	    self.te =  P(ex, cy, tz)
	    self.tn =  P(cx, ny, tz)
	    self.tne = P(ex, ny, tz)
	    self.tnw = P(wx, ny, tz)
	    self.ts =  P(cx, sy, tz)
	    self.tse = P(ex, sy, tz)
	    self.tsw = P(wx, sy, tz)
	    self.tw =  P(wx, cy, tz)

	    # Middle slice:
	    #	NW	N	NE
	    #	W	C	E
	    #	SW	S	SE
	    self.e =   P(ex, cy, cz)
	    self.n =   P(cx, ny, cz)
	    self.ne =  P(ex, ny, cz)
	    self.nw =  P(wx, ny, cz)
	    self.c =   P(cx, cy, cz)
	    self.s =   P(cx, sy, cz)
	    self.se =  P(ex, sy, cz)
	    self.sw =  P(wx, sy, cz)
	    self.w =   P(wx, cy, cz)

	    # Bottom slice:
	    #	BNW	BN	BNE
	    #	BW	B	BE
	    #	BSW	BS	BSE
	    self.b =   P(cx, cy, bz)
	    self.be =  P(ex, cy, bz)
	    self.bn =  P(cx, ny, bz)
	    self.bne = P(ex, ny, bz)
	    self.bnw = P(wx, ny, bz)
	    self.bs =  P(cx, sy, bz)
	    self.bse = P(ex, sy, bz)
	    self.bsw = P(wx, sy, bz)
	    self.bw =  P(wx, cy, bz)

	    # Recursively update parent *Box*'s until top-most bounding
	    # box is reached:
	    up = self.up
	    if type(up) != type(None):
		up._box_recompute("up._box_recompute")

	if trace:
	    print("<=Part._box_recompute({0}, {1})".format(self._name, label))


    def _color_material_update(self, color, material):
	""" *Part*: Update *color* and *material* for *self*. """

	#print("=>Part._color_material_update({0}, {1}): c={2} m={3}".
	#  format(color, material, self._color, self._material))

	# Check argument types:
	none_type = type(None)
	assert type(color) == none_type or isinstance(color, Color)
	assert type(material) == none_type or isinstance(material, Material)

	# Update the *material* and *color* as appropriate:
	if type(self._material) == none_type and isinstance(material, Material):
	    self._material = material
	if type(self._color) == none_type and isinstance(color, Color):
	    self._color = color
	if isinstance(self._material, Material):
	    self._is_part = True

	#print("<=Part._color_material_update({0}, {1}): c={2} m={3}".
	#  format(color, material, self._color, self._material))

    def _cylinder(self, lines = None, indent = 0, is_solid = True,
      comment = None, material = None, color = None, adjust = None,
      diameter = None, start = None, end = None, sides = None,
      sides_angle = None, trace = -1000000, start_diameter = None,
      end_diameter = None):
	""" *Part*: Deal with commonality between holes (i.e. material
	    removal) and cylinders made out of a *material*. """

	# Check argument types:
	none_type = type(None)
	assert isinstance(comment, str)
	assert type(lines) == none_type or isinstance(lines, list)
	assert isinstance(indent, int)
	assert isinstance(is_solid, bool)
	assert type(material) == none_type or isinstance(material, Material)
	assert type(color) == none_type or isinstance(color, Color)
	assert type(diameter) == none_type or isinstance(diameter, L)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(sides, int)
	assert isinstance(adjust, L)
	assert isinstance(sides_angle, Angle)
	assert type(start_diameter) == none_type or \
	  isinstance(start_diameter, L)
	assert type(end_diameter) == none_type or isinstance(end_diameter, L)
	assert isinstance(diameter, L) or \
	  (isinstance(start_diameter, L) and isinstance(end_diameter, L))

	#if trace <= 0 and \
	#   type(start_diameter) != none_type and \
	#  type(end_diameter) != none_type and \
	#  start_diameter != end_diameter:
	#    trace = 0

	if trace >= 0:
	    print(("{0}=>Part._cylinder(comment='{1}', is_solid={2}," + 
	     " diameter={3}, start_diam={4}, end_diam={5}, start={6}, end={7}").
	     format(trace * ' ', comment, is_solid, diameter, start_diameter,
	     end_diameter, start, end))

	# Update the *color* and *material* of this part:
	self._color_material_update(color, material)
	color = self._color
	material = self._material

	# Compute center axis of the *cylinder*, its *length* and its *center*:
	axis = end - start
	length = axis.length()
	half_length = length / 2
	center = (start + end) / 2

	# If appropriate, we rotate the hole before translation:
	rotate_angle = None
	orthogonal_axis = None
	zero = L()
	if axis.x != zero or axis.y != zero:
	    # We have to tilt the cylinder.print("axis={0}".format(axis))
	    z_axis = P(zero, zero, L(mm=1.0))
	    rotate_angle = z_axis.angle_between(axis)
	    orthogonal_axis = z_axis.cross_product(axis)
	    #print("rotate_angle={0:d}".format(rotate_angle))
	    #print("orthogonal_axis={0:m}".format(orthogonal_axis))
	elif axis.x == zero and axis.y == zero and axis.z < zero:
	    orthogonal_axis = P(L(mm=1.0), zero, zero)
	    rotate_angle = Angle(deg=180.0)

	# Update bounding box:
	if is_solid:
	    #print("comment={0} radius={1} half_length={2}". \
	    #  format(comment, radius, half_length))
	    
	    maximum_diameter = zero
	    if isinstance(diameter, L):
		maximum_diameter = diameter
	    if isinstance(start_diameter, L):
		maximum_diameter = maximum_diameter.maximum(start_diameter)
	    if isinstance(end_diameter, L):
		maximum_diameter = maximum_diameter.maximum(end_diameter)
            radius = maximum_diameter / 2

	    bounding_box = Bounding_Box(radius, -radius,
	      radius, -radius, half_length, -half_length)
	    if trace >= 0:
		print("{0}  bounding_box={1}".format(trace * ' ', bounding_box))
		print("{0}  center={1}".format(trace * ' ', center))
	    place = Place(part = None, name = comment,
	      center = None, axis = orthogonal_axis, rotate = rotate_angle,
	      translate = center, trace = trace + 1)
	    self._bounding_box_update(bounding_box,
	      comment, place, trace = trace + 1)

	# Extract some values from {ezcad}:
	if self._ezcad._mode == EZCAD3.CNC_MODE:
	    # Make sure *lines* is a list:
	    assert isinstance(lines, list)
	    assert length > zero, "Cylinder '{0}' has no height".format(comment)

	    # Make sure we have a reasonable number of *sides*:
	    if sides < 0:
		sides = 16

	    # Output a comment to see what is going on:
	    spaces = " " * indent
	    lines.append("{0}// '{1}'".format(spaces, comment))

	    # Ouput the color if appropriate:
	    if is_solid:
		self._scad_color(lines, color, indent = indent)
	    self._scad_transform(lines, indent = indent,
	      center = None, axis = orthogonal_axis, rotate = rotate_angle,
	      translate = center)

	    if sides_angle != Angle():
		lines.append("{0}rotate(a={1}, v=[0, 0, 1])".
		  format(spaces, sides_angle))

	    if type(start_diameter) == none_type:
		start_diameter = diameter
	    if type(end_diameter) == none_type:
		end_diameter = diameter
            r1 = end_diameter / 2
	    if r1 > zero:
		r1 += adjust
	    r2 = start_diameter / 2
	    if r2 > zero:
		r2 += adjust
	    assert r1 + r2 > zero, \
	      "sd={0} ed={1} a={2} r1={3} r2={4}".format(start_diameter,
	      end_diameter, adjust, r1, r2)
	    if r1 == r2:
		command = \
		  "{0}  cylinder(r={1:m}, h={2:m}, center = true, $fn={3});". \
		  format(spaces, r1, length, sides)
	    else:
		command = ("{0}  cylinder(r1={1:m}, r2={2:m}, h={3:m}," +
		  " center = true, $fn={4});").format(spaces,
		  r1, r2, length, sides)
	    if trace >= 0:
		print("{0}Part._cylinder: r1={1} r2={2} command='{3}'".
		  format(trace * ' ', r1, r2, command))

	    # Output the cylinder in a vertical orientation centered on the
	    # origin.  It is processed before either rotation or translation:
            lines.append(command)
            lines.append("")
	
	if trace >= 0:
	    print("{0}<=Part._cylinder()".format(trace * ' '))

    def _dimensions_update(self, ezcad, trace):
	""" *Part*: Update the dimensions of the *Part* object (i.e. *self*)
	    and all of its children *Part*'s.
	"""

	assert isinstance(ezcad, EZCAD3)
	self._ezcad = ezcad

	# Do any requested tracing:
	if trace >= 0:
	    print("{0}=>Part._dimensions_update('{1}')". \
	      format(' ' * trace, self._name))
	    #print("{0}Part._dimensions_update:places={1}". \
	    #  format(' ' * trace, self._places))

	# Start with nothing *changed*:
	changed = 0

	# First record the current values load into *self*.
	sub_parts = []
	before_values = {}
	name = self._name
	for attribute_name in dir(self):
	    attribute = getattr(self, attribute_name)
	    if attribute_name.startswith("__"):
		pass
	    elif attribute_name.endswith("_"):
		assert isinstance(attribute, Part), \
		  "{0}.{1} is not a Part".format(name, attribute_name)
		changed += attribute._dimensions_update(ezcad, trace + 1)
		sub_parts.append(attribute)
	    elif attribute_name.endswith("_l"):
		assert isinstance(attribute, L), \
		  "{0}.{1} is not an L (i.e. length)". \
		  format(name, attribute_name)
		before_values[attribute_name] = attribute
	    elif attribute_name.endswith("_p"):
		assert isinstance(attribute, P), \
		  "{0}.{1} is not a P (i.e. point)". \
		  format(name, attribute_name)
		before_values[attribute_name] = attribute
	    elif attribute_name.endswith("_a"):
		assert isinstance(attribute, Angle), \
		  "{0}.{1} is not an Angle".format(name, attribute_name)
		before_values[attribute_name] = attribute
	    elif attribute_name.endswith("_b"):
		assert isinstance(attribute, bool), \
		  "{0}.{1} is not an bool".format(name, attribute_name)
		before_values[attribute_name] = attribute
	    elif attribute_name.endswith("_f"):
		assert isinstance(attribute, float), \
		  "{0}.{1} is not a float".format(name, attribute_name)
		before_values[attribute_name] = attribute
	    elif attribute_name.endswith("_i"):
		assert isinstance(attribute, int), \
		  "{0}.{1} is not an int".format(name, attribute_name)
		before_values[attribute_name] = attribute
	    else:
		# ignore *attribute_name*:
		pass

	# Reset the bounding box for *self*:
	before_bounding_box = self._bounding_box
	self._bounding_box = Bounding_Box()

	# Peform dimension updating for *self*:
	self.construct()

	# See if anything changed:
	for attribute_name in before_values.keys():
	    before_value = before_values[attribute_name]
	    after_value = getattr(self, attribute_name)
	    assert type(before_value) == type(after_value), \
	      "{0}.{1} before type ({2}) different from after type ({3})". \
	      format(name, attribute_name, type(before_value),
	      type(after_value))
	    if before_value != after_value:
		changed += 1
		if changed > 0:
		    #print("here 2")
		    pass
		if ezcad._update_count > 10 or trace >= 0:
		    print("{0}Part._dimensions_update:{1}.{2} ({3}=>{4})". \
		      format(' ' * trace, name, attribute_name,
		      before_value, after_value))

	# Now update the *bounding_box* for *self*:
	after_bounding_box = self._bounding_box
	for sub_part in sub_parts:
	    after_bounding_box.bounding_box_expand(sub_part._bounding_box)
	if before_bounding_box != after_bounding_box:
	    changed += 1
	    self._bounding_box = after_bounding_box

	    #self._bounding_box_set(after_bounding_box)

	# Update bounding box with placed *Part* bounding boxes:
	# for place in self._places.values():
	#for sub_part in sub_parts:
	#    # Grab some values from *place* and *part*:
	#    #place_part = place._part
	#    #place_name = place._name
	#    place = Place(center = self._center, axis = self._axis,
	#      rotate = self._rotate, translate = self._translate)
	#    forward_matrix = place._forward_matrix
	#
	#    if trace >= 0:
	#	print("{0}Part._dimensions_update:place={1}". \
	#	  format(' ' * trace, place))
	#	#print("{0}Part._dimensions_update:Merge {1:m} into {2:m}". \
	#	#  format(' ' * trace, place_box, box))
	#
	#    self._box_point_update("[TNE]",
	#      forward_matrix.point_multiply(sub_part.tne))
	#    self._box_point_update("[TNW]",
	#      forward_matrix.point_multiply(sub_part.tnw))
	#    self._box_point_update("[TSE]",
	#      forward_matrix.point_multiply(sub_part.tse))
	#    self._box_point_update("[TSW]",
	#      forward_matrix.point_multiply(sub_part.tsw))
	#    self._box_point_update("BNE]",
	#      forward_matrix.point_multiply(sub_part.bne))
	#    self._box_point_update("[BNW]",
	#      forward_matrix.point_multiply(sub_part.bnw))
	#    self._box_point_update("[BSE]",
	#      forward_matrix.point_multiply(sub_part.bse))
	#    self._box_point_update("[BSW]",
	#      forward_matrix.point_multiply(sub_part.bsw))

	# Determine whether *box* has changed:
	#after_box_changed_count = self._box_changed_count
	#if before_box_changed_count != after_box_changed_count:
	#    if ezcad._update_count > 10:
	#	print("{0} bounding box changed".format(name))
	#    changed += 1

	if trace >= 0:
	    print("{0}<=Part._dimensions_update('{1}')=>{2}". \
	      format(' ' * trace, self._name, changed))

	return changed

    def _dx_original_get(self):
	""" *Part*: Return the DX original field of the *Part* object (i.e. *self*)
	"""

	return self._dx_original

    def _dxf_x_offset_get(self):
	""" *Part*: Return the DXF offset X field of the *Part* object (i.e. *self*)
	"""

	return self._dxf_x_offset


    def _dxf_y_offset_get(self):
	""" *Part*: Return the DXF offset Y field of the *Part* object (i.e. *self*)
	"""

	return self._dxf_y_offset

    def _dy_original_get(self):
	""" *Part*: Return the DY original field of the *Part* object (i.e. *self*)
	"""

	return self._dy_original

    def _dz_original_get(self):
	""" *Part*: Return the DZ original field of the *Part* object (i.e. *self*)
	"""

	return self._dz_original

    def _edge_x_get(self):
	""" *Part*: Return the edge X field of the *Part* object (i.e. *self*)
	"""

	return self._edge_x

    def _edge_y_get(self):
	""" *Part*: Return the edge Y field of the *Part* object (i.e. *self*)
	"""

	return self._edge_y

    def _ezcad_get(self):
	""" *Part*: Return the *EZCAD* object from the *Part* object (i.e. *self*)
	"""

	return self._ezcad

    def _flush(self, program_number):
	""" *Part*: Flush out the CNC code for the *Part* object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(program_number, int)
	assert isinstance(self, Part)

	# Use *part* instead of *self*:
	part = self

	#call d@(form@("=>flush@Part(%v%, %d%)\n\") %
	#  f@(part.name) / f@(program_number))
	original_program_number = program_number

	#call show@(part, "before flush")
	shop = part._shop_get()
	assert isinstance(shop, Shop)
	vice = shop._vice_get()
	jaw_width = vice._jaw_width_get()

	# Compute (*plunge_x*, *plunge_y*) which is the vertical axis over
	# which is to the left of the the part or the vice:
	edge_x = part._edge_x_get()
	edge_y = part._edge_y_get()
	vice_x = part._vice_x_get()
	vice_y = part._vice_y_get()
	assert isinstance(vice_x, L)
	plunge_x = vice_x
	if plunge_x > edge_x:
	    plunge_x = edge_x
	assert isinstance(jaw_width, L)
	assert isinstance(plunge_x, L)
	if plunge_x > jaw_width:
	    plunge_x = plunge_x - L(inch=0.7)
	part._plunge_xy_set(plunge_x, edge_y)
    
	code = shop._code_get()
	assert isinstance(code, Code)
	code._z_safe_set(part._z_safe_get())
	code._z_rapid_set(part._z_rapid_get())
	operations = part._operations_get()
	size = len(operations)

	#call show@(part, "before sort", 1t)

	# FIXME move the sort into *operations_regroup*:
	# Sort *operations* to group similar operations together:
	operations.sort(cmp=Operation._compare)

	print("len(operations)={0}".format(len(operations)))
	for operation in operations:
	    print("operation name:{0}".format(operation._name_get()))

	#call show@(part, "after sort", 1t)

	# Do we need to user regroups instead *operations*:
	part._operations_regroup()

	#if regroups != 0
	#	call show@(part, "after_regroup", 0f)

	# FIXME: The code below should replace all the index stuff:
	current_tool = None
	operation_group = None
	operation_groups = []
	for operation in operations:
	    # Grab the *tool* from *operation*:
	    tool = operation._tool_get()
	    assert isinstance(tool, Tool)

	    # Start a new *operation_group* if the *tool* is different:
	    if current_tool != tool:
		# This code is always exectuted the first time through:
		current_tool = tool
		operation_group = []
		operation_groups.append(operation_group)

	    # Tack *operation* onto *operation_group*:
	    assert isinstance(operation_group, list)
	    operation_group.append(operation)
	print("operation_groups=", operation_groups)

	# Open the top-level *part_ngc_stream* file that invokes each tool operation
	# in a separate .ngc file:
	
	ezcad = self._ezcad
	directory = ezcad._directory_get()
	part_ngc_stream_file_name = os.path.join(directory, "O{0}.ngc".format(program_number))
	part_ngc_stream = open(part_ngc_stream_file_name, "w")
	assert part_ngc_stream != None, "Unable to open {0}".format(part_ngc_stream_file_name)

	# Output some heading lines for *part_ngc_stream*:
	part_ngc_stream.write("( Part: {0})".format(part._name_get()))
	part_ngc_stream.write("( Tooling table: )\n")

	# Now vist each *operation_group* in *operation_groups*:
	for index, operation_group in enumerate(operation_groups):
	    # Output a couple of lines *part_ngc_write*:
	    ngc_program_number = program_number + 1 + index
	    tool_number = tool._number_get()
	    tool_name = tool._name_get()
	    part_ngc_stream.write("( T{0} {1} )\n".format(tool_number, tool_name))
	    part_ngc_stream.write("O{0} call\n".format(ngc_program_number))

	    # Sweep through the *operation_group*:
	    part._flush_helper(operation_group, ngc_program_number)

	# Write out the final lines to *part_ngc_stream*:
	part_ngc_stream.write("G53 Y0.0 ( Move the work to the front )\n")
	part_ngc_stream.write("M2\n")
	part_ngc_stream.close()

	# Empty out *operations* :
	del operations[:]

	# Compute the next *program number* to be a the next multiple of 10 and return it:
	program_number = program_number + len(operation_groups) + 1
	program_number = (program_number + 9) / 10 * 10
	return program_number

    def _flush_helper(self, operations, ngc_program_number):
	""" *Part*: Output the G-code for *operations* to a "On.ngc" file where,
	    N is the *ngc_program_number* using the *Part* object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(operations, list) and len(operations) > 0
	for operation in operations:
	    assert isinstance(operation, Operation)
	assert isinstance(ngc_program_number, int)

	# Use *part* instead of *self*:
	part = self
	part_name = part._name

	# Set *debug* to True to trace this routine:
	debug = True
	#debug = False
	if debug:
	    print("=>Part._flush_helper('{0}', *, {1}, *)".
	      format(part_name, ngc_program_number))
	    print("operations=", operations)

	zero = L()

	# Grab some values from *part*:
	shop = part._shop_get()

	# Grab the first *operation*:
	operation = operations[0]
	vice_x = operation._vice_x_get()
	vice_y = operation._vice_y_get()
	tool = operation._tool_get()

	# Get the *feed_speed* and *spindle_speed*:
	feed_speed = tool._feed_speed_get()
	spindle_speed = tool._spindle_speed_get()

	assert tool != None

	# FIXME: The code below should be methodized!!!
	#block = Block(part, tool, vice_x, vice_y)
	#block._spindle_set(s)
	#block._comment_set(operation._comment_get())
	#block._text_set(block._text_get())

	#commands = code._commands
	#blocks = code._blocks
	#print("before len(blocks)={0} len(commands)={1}".format(len(blocks), len(commands)))

	# Figure out if *all_operations_are_drills*:
	all_operations_are_drills = True
	for operation in operations:
	    if not isinstance(operation, Operation_Drill):
		all_operations_are_drills = False
		break

	# FIXME: Move this code to *Operaton_Drill* class!!!
	# Reorder the drill operations to minimize traverses:
	if all_operations_are_drills:
	    for index in range(first_index, last_index + 1):
		current_operation = operations[index]
		assert isinstance(current_operation, Operation_Drill)
		current_drill = current_operation

		current_x = current_drill.x
		current_y = current_drill.y	

		minimum_distance = L(mm=123456789.0)
		match_index = -1
		for search_index in range(index + 1, last_index + 1):
		    search_drill = operations[search_index]
		    assert isinstance(search_drill, Operation_Drill)
		    search_x = search_drill.x
		    search_y = search_drill.y
		    search_distance = \
		      (search_x - current_x).diagonal_2d(search_y - current_y)
                    if search_distance < minimum_distance:
			minimum_distance = search_distance
			match_index = search_index

		# Move the matching operation to be next after {operation}
		# in operations:
		if match_index >= 0:
		    if debug:
			print("minimum_distance[{0}]={1}\n".
			  format(index, minimum_distance))

			# Swap the two drill operations:
			match_operation = operations[match_index]
			operations[match_index] = operations[index + 1]
			operations[index + 1] = match_operation

	# Now do a cnc generate for each *operaton* in *operations*:
	code = None
	for operation in operations:
	    if debug:
		print("operation name:{0}".format(operation._name_get()))
	    # Grab the *spindle_speed*:
	    spindle_speed = operation._spindle_speed_get()

	    # Initialize the *code* object:
	    if code == None:
		# Grab the *code* object and start the code generation:
		code = shop._code_get()
		code._start(part, tool, ngc_program_number, spindle_speed)
		code._dxf_xy_offset_set(part._dxf_x_offset_get(), part._dxf_y_offset_get())
		code._z_rapid_set(part._z_rapid_get())
		code._z_safe_set(part._z_safe_get())
		code._z_set(part._z_safe_get())
		code._vice_xy_set(operation._vice_x_get(), operation._vice_y_get())
		comment = operation._comment_get()
		code._z_safe_assert("flush_helper1", comment)

		# FXIME: Do we really need to do the stuff below?
		# Each block is set to these two values beforehand:
		#code._z_set(part._z_safe_get())
		#code._s_set(s)
		#code._g1_set(0)
		#code._block_append(block)

	    # Perform the CNC generation step for *operation*:
	    operation._cnc_generate()

	if code != None:
	    code._z_safe_retract_actual()
	    code._dxf_xy_offset_set(zero, zero)
	    code._finish()


    ## @brief Formats *self* into a string and returns it.
    #  @param format is the format control string (currently ignored).
    #  @returns a string representation of *self*.
    #
    # *<I>__format__</I>() will return a *self* (a *Box*) object as
    # formatted string.  Currently *format* is ignored, but maybe some
    # time in the future it will contol formaating of the returned string.

    def __format__(self, format):
	""" *Part*: Return formated version of *self*. """

	# Verify argument types:
	assert isinstance(format, str)

	return "[name={0} bb={1}]".format(self._name, self._bounding_box)

	if format != "":
	    format = ":" + format
	format_string = "{0}[{1" + format + "}:{2" + format + "},{3" + \
	  format + "}:{4" + format + "}:{5" + format + "}:{6" + format + "}]"
	#print("format_string='{0}'".format(format_string))

	print("[{0}:{1},{2}:{3},{4}:{5}]".format(
	  self.wx._mm, self.ex._mm,
	  self.sy._mm, self.ny._mm,
	  self.bz._mm, self.tz._mm))

	return format_string.format(self._name,
	  self.wx, self.ex,
	  self.sy, self.ny,
	  self.bz, self.tz)

    def __getattr__(self, name):
	""" *Part*: ..."""

	ezcad = self._ezcad
	bounding_box_dispatch = ezcad._bounding_box_dispatch

	update_count = ezcad._update_count
	first_update = update_count == 0

	if not name.startswith("_") and name.endswith("_"):
	    if first_update:
		pass
	elif name in bounding_box_dispatch:
	    accessor = bounding_box_dispatch[name]
            return accessor(self._bounding_box)
	elif len(name) >= 2 and name[-2] == '_':
	    if name.endswith("_a"):
		if first_update:
		    return Angle()
	    elif name.endswith("_b"):
		if first_update:
		    return False
	    elif name.endswith("_c"):
		if first_update:
		    return Color()
	    elif name.endswith("_f"):
		if first_update:
		    return 0.0
	    elif name.endswith("_i"):
		if first_update:
		    return 0
	    elif name.endswith("_l"):
		if first_update:
		    return L()
	    elif name.endswith("_m"):
		if first_update:
		    return Material()
	    elif name.endswith("_o"):
		if first_update:
		    return None
	    elif name.endswith("_s"):
		if first_update:
		    return ""
	    elif name.endswith("_p"):
		if first_update:
		    return P()
	elif name.endswith("_pl"):
	    if first_update:
		return Place()
	else:
	    raise AttributeError(		
	      "Part instance has no attribute '{0}' count={1}".
	      format(name, update_count))
	raise AttributeError("Part instance has no attribute named '{0}'". \
	  format(name))

    def _manufacture(self, ezcad):
	""" *Part*: Visit a *Part* object (i.e. *self*) and all is
	    sub-*Part*'s and perform any manufacturing steps.
	"""

	# Verify argument types:
	assert isinstance(ezcad, EZCAD3)

	#print("=>Part._manufacture:{0}".format(self._name))

	# Use *part* instead of *self*:
	part = self

	# Make sure *part* is connected ot *ezcad*:
	part._ezcad = ezcad

	# Figure out the mode name:
	mode = ezcad._mode
	mode_name = "?mode?"
	if mode == EZCAD3.VISUALIZATION_MODE:
	    mode_name = "Visual"
	elif mode == EZCAD3.STL_MODE:
	    mode_name = "STL"
	elif mode == EZCAD3.CNC_MODE:
	    mode_name = "CNC"

	# Set *debug* to *True* for debugging:
	debug = False
	#debug = True
	if debug:
	    print("=>Part._manufacture('{0}'):{1}".
	      format(part._name, mode_name))

	# First manufacture any child *Part*'s:
	for attribute_name in dir(part):
	    if not attribute_name.startswith("_") and \
	      attribute_name.endswith("_"):
		child_part = getattr(part, attribute_name)
		assert isinstance(child_part, Part), \
		  "{0}.{1} is not a Part".format(part.name, attribute_name)
		child_part._manufacture(ezcad)

	# Now run construct this *part*
	part.construct()

	# Do the visualization steps:
	if mode == EZCAD3.VISUALIZATION_MODE:
	    ezcad = part._ezcad_get()
	    directory = ezcad._directory_get()
	    wrl_file_name = os.path.join(directory, "{0}.wrl".format(part._name))
	    wrl_file = open(wrl_file_name, "w")
	    part.wrl_write(wrl_file, file_name = wrl_file_name)
	    wrl_file.close()

	# Now generate any CNC files:
	if mode == EZCAD3.CNC_MODE:
	    # Flush out all of the pending CNC operations:

	    # Set *cnc_debug* to *True* to trace CNC:
	    cnc_debug = False
	    cnc_debug = True
	    if cnc_debug:
		print("=>Part._manfacture('{0}'):CNC".format(part._name))

	    shop = ezcad._shop
	    program_base = shop._program_base_get()
	    program_number = part._flush(program_base)

	    # We want the program base number to start with a mulitple of 10.
	    remainder = program_number % 10
	    if remainder != 0:
		program_number += 10 - remainder
	    shop.program_base = program_number
	    if cnc_debug:
		print("<=Part._manfacture('{0}'):CNC".format(part._name))

	# Now generate any .stl files:
	if ezcad._mode == EZCAD3.STL_MODE:
	    # Now manufacture this node:
	    #scad_difference_lines = []
	    #scad_union_lines = []
	    #self._scad_difference_lines = scad_difference_lines
	    #self._scad_union_lines = scad_union_lines
	    #self.construct()
	    #self._scad_difference_lines = None
	    #self._scad_union_lines = None

	    scad_difference_lines = self._scad_difference_lines
	    scad_union_lines = self._scad_union_lines

	    sub_parts = []
	    sub_part_names = []
	    # Find all the *sub_parts* from *self*:
	    for attribute_name in dir(self):
	    	if not attribute_name.startswith("_") and \
	    	  attribute_name.endswith("_"):
	    	    sub_part = getattr(self, attribute_name)
	    	    assert isinstance(sub_part, Part)
		    sub_parts.append(sub_part)
		    sub_part_names.append(sub_part._name)

	    # Remove duplicates from *sub_parts* and *sub_part_names* and sort:
	    sub_parts = list(set(sub_parts))
	    sub_parts.sort(key = lambda sub_part: sub_part._name)
	    sub_part_names = list(set(sub_part_names))
	    sub_part_names.sort()

	    # The *signature* for the Part (i.e. *self*) is the concatenation
	    # of *scad_union_lines* and *scad_difference_lines*:
	    signature = "\n".join(
	      scad_union_lines + scad_difference_lines + sub_part_names)
            signature_hash = hashlib.sha1(signature).hexdigest()
	    self._signature_hash = signature_hash

            # We want to create a *name* the consists of *base_name* followed
	    # by a digit.  We want a unique *name* for each *signature*.

	    # *base_name* is the most specific class name:
	    base_name = self.__class__.__name__

	    # *bases_table* keeps track of each different *base_name*.
	    bases_table = ezcad._bases_table

	    # We deterimine there is a signature table entry for *base_name*:
	    if base_name in bases_table:
		base_signatures_table = bases_table[base_name]
	    else:
		base_signatures_table = {}
		bases_table[base_name] = base_signatures_table

            # Now ensure see whether we have encountered this *signature*
	    # before:
	    if signature in base_signatures_table:
		# Yes we have, so reuse the previously assigned name:
		name = base_signatures_table[signature]
		self._name = name
	    else:
		# No we have not so we create a new one *and* write out
		# the associated *name*.scad file:
		name = "{0}{1}".format(base_name, len(base_signatures_table))
		base_signatures_table[signature] = name
		self._name = name

		# Now we see whether we need to write out the file and
		# run it through openscad:
		directory = ezcad._directory_get()
		stl_file_name = os.path.join(directory, "{0}_{1}.stl".format(name, signature_hash))
		if os.path.isfile(stl_file_name):
		    # Since the .stl file already exists, we must have already
		    # written out the .scad file.  Thus, there is nothing
		    # more to do:
		    pass
		else:
		    # We need to write out the .scad file and generate the
		    # assocatied .stl file:

		    # Deal with the part *places*:
		    lines = []
		    #place_parts = {}

		    # Output the use list:
		    for sub_part_name in sub_part_names:
			lines.append("use <{0}.scad>;".format(sub_part_name))

		    # Write out the module:
		    lines.append("module {0}() {{".format(name))

		    # Get the difference() followed union():
		    lines.append("  difference() {")
		    lines.append("    union() {")

		    # Output the *scan_union_lines*:
		    for union_line in scad_union_lines:
			lines.append(union_line)

		    # Close off union():
		    lines.append("    }")

		    # Output *scad_difference_lines*:
		    for difference_line in scad_difference_lines:
			lines.append(difference_line)

		    # Close off difference():
		    lines.append("  }")

		    # Perform all the placements:
		    for sub_part in sub_parts:
			#print("Part._manufacture.place={0}".format(place))
			self._scad_transform(lines, center = sub_part._center,
			  axis = sub_part._axis, rotate = sub_part._rotate,
			  translate = sub_part._translate);
			lines.append("{0}();".format(sub_part._name))

		    # Close off the module:
		    lines.append("}")
		    lines.append("")

		    # Call the module we just produced:
		    lines.append("{0}();".format(name))
		    lines.append("")

		    # Write out *scad_file*:
		    directory = ezcad._directory_get()
		    scad_file_name = os.path.join(directory, "{0}.scad".format(name))
		    scad_file = open(scad_file_name, "w")
		    scad_file.write("\n".join(lines))
		    scad_file.close()

		    # Delete any previous *.stl file:
		    directory = ezcad._directory_get()
		    glob_pattern = os.path.join(directory, "{0}_*.stl".format(name))
		    previous_stl_files = glob.glob(glob_pattern)
		    for previous_stl_file in previous_stl_files:
			os.remove(previous_stl_file)

		    # Run the command that convert the .scad file into the
		    # associated .stl file:
		    if self._is_part:
			ignore_file = open("/dev/null", "w")
			scad_file_name = os.path.join(directory, "{0}.scad".format(name))
			command = [ "openscad", "-o", stl_file_name,  scad_file_name ]
			#print("command=", command)
			subprocess.call(command, stderr=ignore_file) 
			ignore_file.close()

		    # Write out DXF file:
		    dxf_scad_lines = self._dxf_scad_lines
		    if len(dxf_scad_lines) > 0:
			directory = ezcad._directory_get()
			dxf_scad_file_name = os.path.join(directory, "{0}_Dxf.scad".format(name))
			dxf_scad_file = open(dxf_scad_file_name, "wa")
			dxf_scad_file.write("use <{0}.scad>;\n".format(name))
			for dxf_scad_line in dxf_scad_lines:
			    dxf_scad_file.write(dxf_scad_line)
			    dxf_scad_file.write("\n")
			dxf_scad_file.write("{0}();\n".format(name))
			dxf_scad_file.close()
			self._dxf_scad_lines = []

			dxf_file_name = os.path.join(directory, "{0}.dxf".format(name))
			dxf_scad_file_name = os.path.join(directory, "{0}_Dxf.scad".format(name))
			command = [ "openscad" , "-o", dxf_file_name, dxf_scad_file_name ]
			print("openscad command = '{0}'".format(" ".join(command)))
			ignore_file = open("/dev/null", "w")
			subprocess.call(command, stderr=ignore_file)
			ignore_file.close()
    
	if debug:
	    print("<=Part._manufacture('{0}'):{1}".
	      format(part._name, mode_name))

    def _material_get(self):
	""" *Part*: Return the matrial associated with the *Part* object (i.e. *self*.)"""

	return self._material

    def _name_get(self):
	""" *Part*: Return the name of the *Part* object (i.e. *self*). """

	return self._name

    def _operation_append(self, operation):
	""" Part*: Append *operation* to the operations list in the *Part* object (i.e. *self*). """

	# Verify argument types:
	assert isinstance(operation, Operation)

	# Append *operation* to *operations*:
	self._operations.append(operation)

    def _operation_regroup_check(self,
      first_index1, last_index1, first_index2, last_index2):
	""" *Part*:
	"""

	# Verify argument types:
	isinstance(first_index1, int)
	isinstance(last_index1, int)
	isinstance(first_index2, int)
	isinstance(last_index2, int)

	# Use *part* instead of *self*:
	part = self

	operations = part.operations
	index = first_index2
	while index <= last_index2:
	    operation = operations[index]
	    follows = operation.follows
	    if follows != None:
		follows_index = follows.index
		if first_index1 <= follows_index and \
		  follows_index <= last_index1:
		    return False
	    index += 1
	return True

    def _operations_get(self):
	""" *Part*: Return the operations list for the *Part* object (i.e. *self*).
	"""
        
	return self._operations

    def _operations_index(self):
	""" *Part*: Sweep through each operation in the *Part* object
	    (i.e. *self*) and make sure the index field matches the index.
	"""

	part = self
	operations = part._operations
	size = len(operations)
	for index in range(size):
	    operations[index].index = index

    def _operations_regroup(self):
	""" *Part*: Regroup the *Operation*'s in the *Part* object (i.e. *self*)
	    to try to minimize tool changes.  No *Operation* will moved in
	    front of its *follows* field if the {follows} field is non-null.
	"""

	# Use *part* instead of *self*:
	part = self

	elimination_count = 0
	part._operations_index()
	operations = part._operations
	size = len(operations)

	# Find a group of operations with the same *priority* field value:
	first_index = 0
	while first_index < size:
	    operation = operations[first_index]
	    priority = operation._priority_get()

	    # Now find the last operation in *operations* with the same
	    # *priority*:
	    last_index = first_index
	    while last_index + 1 < size and \
	      operations[last_index + 1]._priority_get() == priority:
		last_index = last_index + 1

	    # Now pass these values down to the next level of helper routine:
	    elimination_count = elimination_count + \
	      part._operations_regroup_helper1(first_index, last_index)

	    first_index = last_index + 1
	return elimination_count

    def _operations_regroup_helper1(self, first_index, last_index):
	""" *Part*: Regroup the *first_index*'th through *last_indxex*'th
	    *Operation*'s of the *Part* object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(first_index, int)
	assert isinstance(last_index, int)

	# Use *part* instead of *self*:
	part = self

	# Find a group of operations that use the same {Tool}:
	elimination_count = 0
	operations = part._operations
	tool1_first_index = first_index
	while tool1_first_index <= last_index:
	    operation = operations[tool1_first_index]
	    tool = operation._tool_get()

	    # Find the last operation in this sequence with the same {tool}:
	    tool1_last_index = tool1_first_index
	    while tool1_last_index + 1 <= last_index and \
	      operations[tool1_last_index + 1].tool == tool:
		tool1_last_index += 1

	    # Pass these values down to the next level of helper routine:
	    elimination_count += part._operations_regroup_helper2(
	      tool1_first_index, tool1_last_index, last_index)

	    tool1_first_index += 1

	# Return the number of tool changes removed.
	return elimination_count

    def _operations_regroup_helper2(self,
      tool1_first_index, tool1_last_index, last_index):
	""" *Part*:
	"""

	# Verify argument types:
 	assert isinstance(tool1_first_index, int)
	assert isinstance(tool1_last_index, int)
	assert isinstance(last_index, int)

	# Use *part* instead of *self*:
	part = self

	elimination_count = 0
	operations = part._operations
	tool = operations[tool1_first_index]._tool_get()

	# Find the next sequence of operations that match {tool}:
	for index in range(tool1_last_index + 1, last_index + 1):
	    operation = operations[index]
	    if operation.tool == tool:
		# We have another instance of {tool}:
		tool2_first_index = index

		# Find the last operation in the sequence that matches {tool}:
		tool2_last_index = tool2_first_index
		while tool2_last_index + 1 <= last_index and \
		  operations[tool2_last_index + 1].tool == tool:
		    tool2_last_index += 1

		# We now have two tool sequences that use the same tool:
		#call d@(form@("*****Tool:%v% [%d% - %d%] and [%d% - %d%]\n\") %
		#  f@(tool.name) % f@(tool1_first_index) %
		# f@(tool1_last_index) %
		#  f@(tool2_first_index) / f@(tool2_last_index))

		if part._operation_regroup_check(part,
		  tool1_last_index + 1, tool2_first_index - 1,
		  tool2_first_index, tool2_last_index):
		    #call d@("**************reorder possible\n\")
		    part._operation_regroup_helper3(tool1_last_index + 1,
		      tool2_first_index, tool2_last_index)
		    elimination_count += 1

		# We can continue looking for another group of operations with
		# the same {tool} to move:
		index = tool2_last_index
	return elimination_count

    def _operation_regroup_helper3(self, to_index, from_first_index, from_last_index):
	""" *Part*:
	"""

	# Verify argument types:
	assert isinstance(to_index, int)
	assert isinstance(from_first_index, int)
	assert isinstance(from_last_index, int)

	# Use *part* instead of *self*:
	part = self

	# FIXME: There is probably a more Pythony wayd of doing this!!!
	operations = part._operations
	for index in range(from_first_index, from_last_index + 1):
	    operation = operations[index]
	    del operations[index]
	    operations.insert(to_index, operation)
	    to_index = to_index + 1
	part._operations_index()

    def _plunge_x_get(self):
	""" *Part*: Return the plunge X of the *Part* (i.e. *self*). """

	return self._plunge_x

    def _plunge_y_get(self):
	""" *Part*: Return the plunge Y of the *Part* (i.e. *self*). """

	return self._plunge_y

    def _plunge_xy_set(self, plunge_x, plunge_y):

	""" *Part*: Set the plunge X/Y coordinate to (*plunge_x*, *plunge_y*) for the
	    *Part* object (i.e. *self*)
	"""

	# Verify argument types:
	assert isinstance(plunge_x, L)
	assert isinstance(plunge_y, L)

	# Set the plunge X/Y coordinates in the *Part* object (i.e. *self*):
	self._plunge_x = plunge_x
	self._plunge_y = plunge_y

    def _priority_get(self):
	""" *Part*: Return the priority of the *Part* (i.e. *self*). """

	return self._priority

    def _tools_dowel_pin_search(self):
	""" *Part*: Find and return a *Tool_Dowel_Pin* object.
	"""

	# Use *part* instead of *self*:
	part = self

	zero = L()
	dowel_pin = part._tools_search(Tool_Dowel_Pin._match, zero, zero, "dowel_pin")
	assert isinstance(dowel_pin, Tool_Dowel_Pin)
	return dowel_pin

    def _tools_end_mill_search(self, maximum_diameter, maximum_z_depth, from_routine):
	""" *Part*: Search for an end mill with a diameter that is less than or equal to
	    *maximum_diameter* and can mill to a depth of *maximum_z_depth*.  (*from_routine*
	    is used for debugging.
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(maximum_diameter, L) and maximum_diameter > zero
	assert isinstance(maximum_z_depth, L)
	assert isinstance(from_routine, str)

	result = \
	  self._tools_search(Tool_End_Mill._match, maximum_diameter, maximum_z_depth, from_routine)

	return result

    def _position_reset(self):
	""" *Part*: Reset the position and reposition matrix for the
	    *Part* object (i.e. *self*).
	"""

	# Make sure that all previous CNC operations are done before
	# proceeding to future CNC operations:
	self.cnc_fence()

	# Grap *position* and *reposition* *Matrix* objects from *self*:
	position = self._position
	reposition = self._reposition

	# Reset *position* and *reposition* to identity matrix:
	position.identity_store()
	position.identity_store()

    def _shop_get(self):
	""" *Part*: Return the *Shop* object associated with the *Part* object (i.e. *self*.) """

	shop = self._shop
	assert isinstance(shop, Shop)
	return shop

    def _tools_search(self, match_routine, parameter1, parameter2, from_routine):
	""" *Part*: Search for the best *Tool* object using *match_routine*,
	    *parameter1*, and *parameter2*.  *from_routine* is used for debugging*:
	"""

	# Verify argument types:
	assert isinstance(parameter1, L)
	assert isinstance(parameter2, L)
	assert isinstance(from_routine, str)

	# Use *part* instead of *self*:
	part = self

	# Set *debug* to *True* if debugging is required:
	debug = False
	#debug = True
	if debug:
	    print("=>Part._tools_search('{0}', {4}, {1}, {2}, '{3}')\n".
	      format(part._name, parameter1, parameter2, from_routine,
	      match_routine))

	# For now grab the material of the first part and assume the
	# rest are compatible.  In reality, what we want is to grab
	# all the materials and only accept tools that work with all
	# materials.  That is too much work for now:
	shop = part._shop_get()
	part_material = part._material_get()

	best_tool = None
	best_priority = -1.0

	# For now, hardwire the mill speeds.  We need to actually get the
	# speeds from the *Mill* object:
	maximum_spindle = Hertz(rpm=5000.0)
	minimum_spindle = Hertz(rpm=150.0)

	# FIXME: This code probably belongs over in *Shop*:
	# Now search through {tools}:
	zero = L()
	delta = L(inch=.000001)
	tools = shop._tools_get()
	if debug:
	    print("  {0} available tools".format(len(tools)))
	for tool in tools:
	    # Keep track of search results for this tool:
	    tool._search_results_clear()

	    # Lookup the *speed_range* for *material* and *tool*:
	    tool_material = tool._material_get()
	    speed_range = shop._surface_speeds_lookup(part_material, tool_material)
	    speed_range_ok = isinstance(speed_range, Speed_Range)

	    # Always log *surface_speeds_ok*:
	    tool._search_results_append(speed_range_ok,
	      "Surface speeds are acceptable for part_material '{0}' and tool material '{1}'".
	      format(part_material, tool_material))

	    if debug:
		print("  Tool:'{0} {1}'".format(tool._name_get(), speed_range_ok))

	    if speed_range_ok:
		# See whether *tool* has a chance of being a match:
		#if debug:
		#    print("    speed_range={0:.0F}\n".format(speed_range))
		priority = match_routine(tool, parameter1, parameter2, from_routine)
		if priority >= 0.0:
		    # Tool is an acceptable match:
		    if debug:
			print("Tool: '{0}' priority:{1}\n".
			  format(tool._name_get(), tool._diameter_get()))

		    tool_preferred = part._tool_preferred
		    if tool_preferred != "":
			#print("tool={0} preferred={1}\n".
			#  format(tool.name, tool_preferred))
			if tool.name == tool_preferred:
			    priority = priority + 100.0

		    # Select the {Surface_Speed} for {tool_material}:
		    tool.priority = priority
		    tool.spindle = Hertz()
		    tool.feed = Speed()
		    tool_material = tool._material_get()
		    surface_speed = None
		    if tool._flutes_count_get() > 0:
			# We have {surface_speed}, compute {low_speed} and
			# {high_speed}:
			low_speed = speed_range._low_speed_get()
			high_speed = speed_range._high_speed_get()

			#print("Tool diameter={0}\n".format(tool.diameter))

			pi = 3.14159265358979323846
			desired_spindle = low_speed / (tool._diameter_get() * pi)

			#print("low_speed={0} diam={1}\n".
			#  format(low_speed, tool.diameter))
			#print("desired_spindle={0}\n".
			#  format(desired_spindle))

			# Deal with spindle speed limits:
			if desired_spindle > maximum_spindle:
			    # Speed is to fast, clamp it to *maximum_spindle*:
			    desired_spindle = maximum_spindle
			assert desired_spindle >= minimum_spindle,				\
			  "Disired spindle ({0:rpm}) is less than minimum spindle({1:rpm})".	\
			  format(desired_spindle, minimum_spindle)

			# Record everything we figured out back in {tool}:
			tool._spindle_speed_set(desired_spindle)
			chip_load = L(inch=0.001)
			flutes_count = tool._flutes_count_get()

			# Temporarily store the *feed_speed*, *spindle_speed*  and *priority*
			# in the *tool*.  These values are only good until the end of this
			# routine.  They will possibly get overridden during a subsequent
			# tool search:
			feed_speed = Speed(mm_per_sec=
			  ((chip_load * float(flutes_count)).millimeters() *
			   desired_spindle.frequency()))
			tool._feed_speed_set(feed_speed)
			tool._priority_set(priority)

		    # Record {tool} as a match in {tools_matach}:
		    if priority > best_priority:
			best_priority = priority
			best_tool = tool

	if best_tool == None:
	    closest_tool = None
	    largest_success_count = 0
	    for tool in tools:
		# Check {tool} to see if it is exceptable:
		success_count = 0
		search_results = tool._search_results_get()
		search_results_size = len(search_results)

		for search_result in search_results:
		    if search_result[0] == '1':
			success_count += 1

		if success_count > largest_success_count:
		    largest_success_count = success_count
		    closest_tool = tool

	    if closest_tool != None:
		closest_tool.search_results_show()

	if debug:
	    best_tool_name = "?"
	    if best_tool != None:
		best_tool_name = best_tool._name_get()
	    print("<=Part._tools_search({0}, *, {1}, {2}, {3}) => {4}\n\n".
	      format(part._name, parameter1, parameter2, from_routine,
	      best_tool_name))

	assert isinstance(best_tool, Tool), "Could not find a tool that worked"

	return best_tool

    def _vice_x_get(self):
	""" *Part*: Return the vice X field of the *Part* object (i.e. *self*)
	"""

	return self._vice_x


    def _vice_y_get(self):
	""" *Part*: Return the vice Y field of the *Part* object (i.e. *self*)
	"""

	return self._vice_y


    def _z_safe_get(self):
	""" *Part*: Return the Z safe height for vertical air plunging for *Part* (i.e. *self*.)
	"""
	
	return self._z_safe

    def _z_rapid_set(self, z_safe):
	""" *Part*: Set the Z safe height for vertical air plunging for *Part*
	    (i.e. *self*) to *z_save*.
	"""
	
	# Verify argument types:
	assert isintance(z_safe, L)

	# Load *z_rapid* into the *Part* object (i.e. *self*):
	self._z_rapid = z_safe

    def _z_rapid_get(self):
	""" *Part*: Return the Z rapid height for horizontal movement above hold-down tooling
	    for *Part* (i.e. *self*.)
	"""
	
	assert isinstance(self._z_rapid, L)
	return self._z_rapid

    def _z_rapid_set(self, z_rapid):
	""" *Part*: Set the Z rapid height for horizontal movement above hold-down tooling
	    for *Part* (i.e. *self*) to *z_rapid*
	"""
	
	# Verify argument types:
	assert isintance(z_rapid, L)

	# Load *z_rapid* into the *Part* object (i.e. *self*):
	self._z_rapid = z_rapid

    # Public methods come here:

    def block(self, comment = "no comment", material = None, color = None,
      corner1 = None, corner2 = None, welds = "", trace = False,
      center = None, axis = None, rotate = None, translate = None, top = None):
	""" {Part} construct: Create a block with corners at {corner1} and
	    {corner2}.  The block is made of {material} and visualized as
	    {color}. """

	if trace:
	    print("=>Part.block(comment={0}".format(comment))

	#print "block_corners('{0}', {1}, {2}, '{3}', '{4}')".format( \
	#  self.name, corner1, corner2, color, material)

	# Check argument status:
	none_type = type(None)
	assert type(comment) == none_type or isinstance(comment, str)
	assert type(corner1) == none_type or isinstance(corner1, P)
	assert type(corner2) == none_type  or isinstance(corner2, P)
	assert type(material) == none_type or isinstance(material, Material)
	assert type(welds) == none_type or isinstance(welds, str)
	assert type(color) == none_type or isinstance(color, Color)
	assert type(center) == none_type or isinstance(center, P)
	assert type(axis) == none_type or isinstance(axis, P)
	assert type(rotate) == none_type or isinstance(rotate, Angle)
	assert type(translate) == none_type or isinstance(translate, P)
	assert isinstance(top, str)

	# Deal with argument defaults:
	self._color_material_update(color, material)
	color = self._color
	material = self._material

	if type(corner1) == none_type:
	    zero = L(0.0)
	    corner1 = P(zero, zero, zero)
	if type(corner2) == none_type:
	    one = L.mm(1.0)
	    corner2 = P(one, one, one)
	assert isinstance(welds, str)

	# Record the *color* and *material*:
	self._material = material
	self._color = color

	# Make sure that the corners are diagonal from bottom south west
	# to top north east:
	x1 = min(corner1.x, corner2.x)
	x2 = max(corner1.x, corner2.x)
	y1 = min(corner1.y, corner2.y)
	y2 = max(corner1.y, corner2.y)
	z1 = min(corner1.z, corner2.z)
	z2 = max(corner1.z, corner2.z)
	#print("Part.box:{0:m}:{1:m},{2:m}:{3:m},{4:m}:{5:m}". \
	#  format(x1, x2, y1, y2, z1, z2))

	bounding_box = self._bounding_box
	bounding_box.point_expand(P(x1, y1, z1))
	bounding_box.point_expand(P(x2, y2, z2))

	#self._dx_origin = x2 - x1
	#self._dy_origin = y2 - y1
	#self._dz_origin = z2 - z1

	#place = Place(part = None, name = comment, center = center,
	#  axis = axis, rotate = rotate, translate = translate)
	#forward_matrix = place._forward_matrix

	#tne = P(x2, y2, z2)
	#bsw = P(x1, y1, z1)
	#tsw = P(x1, y1, z2)
	#bnw = P(x1, y2, z1)
	#tnw = P(x1, y2, z2)
	#bse = P(x2, y1, z1)
	#tse = P(x2, y1, z2)
	#bne = P(x2, y2, z1)

	#bounding_box = self._bounding_box
	#bounding_box.point_expand(tne)
	#bounding_box.point_expand(tnw)
	#bounding_box.point_expand(tse)
	#bounding_box.point_expand(tsw)
	#bounding_box.point_expand(bne)
	#bounding_box.point_expand(bnw)
	#bounding_box.point_expand(bse)
	#bounding_box.point_expand(bsw)

	#print("before box={0:m}".format(box))
	#self._box_point_update(comment + "[TNE]",
	#  forward_matrix.point_multiply(tne))
	#self._box_point_update(comment + "[TNW]",
	#  forward_matrix.point_multiply(tnw))
	#self._box_point_update(comment + "[TSE]",
	#  forward_matrix.point_multiply(tse))
	#self._box_point_update(comment + "[TSW]",
	#  forward_matrix.point_multiply(tsw))
	#self._box_point_update(comment + "[BNE]",
	#  forward_matrix.point_multiply(bne))
	#self._box_point_update(comment + "[BNW]",
	#  forward_matrix.point_multiply(bnw))
	#self._box_point_update(comment + "[BSE]",
	#  forward_matrix.point_multiply(bse))
	#self._box_point_update(comment + "[BSW]",
	#  forward_matrix.point_multiply(bsw))
	#print("after box={0:m}".format(box))


	ezcad = self._ezcad

	self._is_part = True
	#print("Part.block:{0}._is_part = True".format(self._name))
        
	if ezcad._mode == EZCAD3.CNC_MODE:
	    union_lines = self._scad_union_lines

	    #FIXME: Should only check after the dimensions update:
	    #assert x1 < x2, \
	    # "{0}.block '{1}': equal X coordinates: corner1={2} corner2={3}". \
	    # format(self._name, comment, corner1, corner2)
	    #assert y1 < y2, \
	    # "{0}.block '{1}': equal Y coordinates: corner1={2} corner2={3}". \
	    # format(self._name, comment, corner1, corner2)
	    #assert z1 < z2, \
	    # "{0}.block '{1}': equal Z coordinates: corner1={2} corner2={3}". \
	    # format(self._name, comment, corner1, corner2)


	    # Now make the block a little bigger for "welding":
	    adjust = ezcad._adjust * 2
	    weld_extra = adjust.absolute() + L(mm=.01)
	    if welds.find("t") >= 0:
		z2 += weld_extra
	    if welds.find("b") >= 0:
		z1 -= weld_extra
	    if welds.find("n") >= 0:
		y2 += weld_extra
	    if welds.find("s") >= 0:
		y1 -= weld_extra
	    if welds.find("e") >= 0:
		x2 += weld_extra	
		#print("weld e:{0}".format(x2))
	    if welds.find("w") >= 0:
		x1 -= weld_extra
		#print("weld w:{0}".format(x1))

	    # Deal with *adjust*:
	    zero = L()
	    adjust = ezcad._adjust
	    if top == "t" or top == "b":
		#print("block adjust={0}".format(adjust))
		x1 -= adjust
		x2 += adjust
		y1 -= adjust
		y2 += adjust
	    elif top == "n" or top == "s":
		x1 -= adjust
		x2 += adjust
		z1 -= adjust
		z2 += adjust
	    elif top == "e" or top == "w":
		y1 -= adjust
		y2 += adjust
		z1 -= adjust
		z2 += adjust
	    else:
		assert adjust == zero, "Part.Block: top argument not set"

	    #print "c1=({0},{1},{2}) c2=({3},{4},{5})".format( \
	    #  x1, y1, z1, x2, y2, z2)

            # The transforms are done in reverse order:
	    self._scad_transform(union_lines, center = center,
	      axis = axis, rotate = rotate, translate = translate)

	    # Get the lower south west corner positioned:
	    union_lines.append(
	      "      translate([{0:m}, {1:m}, {2:m}])".format(x1, y1, z1))

	    # The color can be output any old time before the cube:
	    union_lines.append(
	      "      color([{0}, {1}, {2}, {3}])".
	      format(color.red, color.green, color.blue, color.alpha))

	    # Finally, we can output the cube:
	    union_lines.append(
	      "        cube([{0:m}, {1:m}, {2:m}]);".format(
	      x2 - x1, y2 - y1, z2 - z1))
	    #print("cube: x2(={0}) - x1(={1}) = {2}".format(x2, x1, x2 - x1))

	if trace:
	    print("<=Part.block(comment={0}".format(comment))

    def cnc_fence(self):
	""" *Part*: Erect a virtual fence that prevents any CNC operations
	    that already exist from moving forward across the fence and
	    and future CNC operations from moving backward across the fence.
	"""

	# This routine will conceptually erect a fence that prevents the
	# CNC planner from moving any operations across the fence.  Any
	# operations requested before the call to this routine will happen
	# before all requested operations after this routine call.

	# Each CNC operation is taged with the value of *part._priority*
	# when it is created.  By incrementing this field, all operations
	# get grouped by priority first, and other characteristics second.
	self._priority += 1

    def construct(self):
	assert False, \
	  "No construct() method defined for part '{0}'".format(self._name)

    def cylinder(self, comment = "NO_COMMENT",
      material = Material(), color = Color(),
      diameter = L(mm = 1.0), start = P(), end = P(z = L(mm = 1.0)),
      start_diameter = None, end_diameter = None,
      sides = -1, sides_angle = Angle(), welds = "", flags = "", adjust = L(),
      trace = -1000000):
	""" *Part*: Place a *diameter* wide cylinder from *start* to *end*. """

	if trace >= 0:
	    print("{0}=>Part.cylinder(diam={1} start={2} end={3})".
	      format(trace * ' ', diameter, start, end))

	# Check argument types:
	none_type = type(None)
	assert isinstance(comment, str)
	assert type(material) == none_type or isinstance(material, Material)
	assert type(color) == none_type or isinstance(color, Color)
	assert isinstance(diameter, L)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(sides, int)
	assert isinstance(sides_angle, Angle)
	assert isinstance(welds, str)
	assert isinstance(flags, str)
	assert type(start_diameter) == none_type or \
	  isinstance(start_diameter, L)
	assert type(end_diameter) == none_type or \
	  isinstance(end_diameter, L)

	ezcad = self._ezcad
	adjust = ezcad._adjust

	union_lines = self._scad_union_lines
	self._cylinder(lines = union_lines, indent = 6, is_solid = True,
	  comment = comment, material = material, color = color,
	  diameter = diameter, start = start, end = end, sides = sides,
	  adjust = adjust, sides_angle = sides_angle, trace = trace + 1,
	  start_diameter = start_diameter, end_diameter = end_diameter)

	if trace >= 0:
	    print("{0}<=Part.cylinder()".format(trace * ' '))

    def dowel_pin(self, comment):
	""" *Part*: Request that a dowel pin be used to align the *Part* object (i.e. *self*)
	    in the vice with a comment of *comment*.
	"""

	# Verify argument types
	assert isinstance(comment, str)

	# Use *part* instead of *self*:
	part = self

	shop = part._shop
	assert isinstance(shop, Shop)
	if shop._cnc_generate_get():
	    # Find the *tool_dowel_pin* to use :
	    tool_dowel_pin = part._tools_dowel_pin_search()
	    assert isinstance(tool_dowel_pin, Tool_Dowel_Pin), "No dowel pin tool found"

	    # Immediately grab the *feed_speed* and *spindle_speed* that were computed
	    # for *tool_dowel_pin*:
	    feed_speed = tool_dowel_pin._feed_speed_get()
	    spindle_speed = tool_dowel_pin._spindle_speed_get()

	    debug = False
	    debug = True
	    if debug:
		print("=>Part.dowel_pin('{0}', {1})".
		  format(part._name, comment))
	    diameter = tool_dowel_pin._diameter_get()

	    # Figure out the new bounding box orientation:
	    new_bounding_box = part._bounding_box.matrix_apply(part._position)

	    # Compute the Z minimum of the bounding box (negative number):
	    bounding_box_z_minimum = new_bounding_box.bsw_get().z

	    # Now compute {z_depth}, the desired maximum tool descent
	    # as a positive number:
	    tip_depth = tool_dowel_pin._tip_depth_get()
	    z_depth = tip_depth - bounding_box_z_minimum
	    if debug:
		print("dowel_pin: part={0} bb={1} tip_depth={2}".
		  format(part._name, new_bounding_box, tip_depth))
		print("dowel_pin: before z_depth={0}".format(z_depth))

	    # Figure out how deep the dowel pin actually can go:
	    maximum_z_depth = tool_dowel_pin._maximum_z_depth_get()
	    if z_depth > maximum_z_depth:
		z_depth = maximum_z_depth
	    if debug:
		print("after z_depth={0}".format(z_depth))

	    if comment == "":
		comment = "Mount {0} x {1} x {2} piece of {3} in vice". \
		  format(part._dx_original, part._dy_original,
		  part._dz_original, _part.material)

	    if debug:
		print("dowel_pin@({0}): ex={1} ey={2} vx={3} vy={4}".
		  format(part._name, part._edge_x, part._edge_y,
		  part._vice_x, part._vice_y))

	    vice = shop._vice_get()
	    jaw_width = vice._jaw_width_get()
	    half_jaw_width = jaw_width / 2
	    radius = diameter / 2
	    edge_x = part._edge_x - radius
	    edge_y = part._edge_y
	    plunge_x = edge_x

	    if debug:
		print("dowel_pin: {0} jaw_width={1}".
		  format(part._name, jaw_width))

	    # Make sure we always plunge outside of the vise jaws:
	    if plunge_x > -half_jaw_width:
		plunge_x = -half_jaw_width

	    # If we are close the jaw edge, move out a little more:
	    if plunge_x > -jaw_width:
		plunge_x = plunge_x - L(inch=0.70)
	    plunge_y = edge_y

	    if debug:
		print("dowel_pin: {0}: plunge_x={1}".
		  format(part._name, plunge_x))

	    dowel_pin_order = Operation.ORDER_DOWEL_PIN
	    operation = Operation_Dowel_Pin(self,
	      comment, 0, tool_dowel_pin, dowel_pin_order, None, feed_speed, spindle_speed,
	      diameter, edge_x, edge_y, part._dx_original, part._dy_original, part._dz_original,
	      plunge_x, plunge_y, tip_depth, -z_depth)
	    self._operation_append(operation)

	    if debug:
		print("dowel_pin({0} edgex={1}% edgey={2}".
		  format(part._name, edge_x, edge_y))
	    if debug:
		print("<=dowel_pin@Part('{0}', '{1}')".
		  format(part._name, comment))

    def dowel_position_set(self, edge_x, edge_y):
	""" *Part*:
	"""

	# This routine will sets the dowel position for {part} to
	# ({edge_x},{edge_y}).

	# Verify argument types:
	assert isinstance(edge_x, L)
	assert isinstance(edge_y, L)

	self._edge_x = edge_x
	self._edge_y = edge_y

    def dxf_write(self, comment = "no_comment", center = None,
      plane_normal = P(0, 0, L(mm=1)), rotate = Angle()):
	""" *Part*: Write out a .dxf file. """
	# Check argument
	none_type = type(None)
	assert isinstance(comment, str)
	assert isinstance(plane_normal, P)
	assert type(center) == none_type or isinstance(center, P)
	assert isinstance(rotate, Angle)

	if not isinstance(center, P):
	    center = self.c

	scad_lines = []
	z_axis = P(0, 0, L(mm=1))
	# OpenSCAD is screwy.  We list the operations in reverse order:
	# The projection is the *last* operation.  See below:
	scad_lines.append("projection(cut = true)")

	# If we want to rotate in the X/Y plane, we do this here.
	# Will occur *after* the rotation and translation below:
	if rotate != Angle(deg=0):
	    scad_lines.append("rotate(a={0:d}, v=[0, 0, 1])".format(rotate))

	# If we need to do a rotation, we do this here.  This will occur
	# *after* the translation (see below):
	if plane_normal != z_axis:
	    rotate_axis = plane_normal.cross_product(z_axis)
	    scad_lines.append("rotate(a={0:d}, v=[{1:m}, {2:m}, {3:m}])".format(
	      Angle(deg=90), rotate_axis.x, rotate_axis.y, rotate_axis.z))

	# If we need to specify the *center*, we list this command last
	# so it occurs first:
	if center != P():
	    scad_lines.append("translate([{0:m}, {1:m}, {2:m}])".
	      format(-center.x, -center.y, -center.z))

	self._dxf_scad_lines = scad_lines

    def extrude(self, comment = "no_comment", material = None, color = None,
      outer_contour = None, inner_contours = None, start = None, end = None,
      center = None, axis = None, rotate = None, translate = None,
      trace = -1000000):
	""" *Part*: """

	if trace >= 0:
	    print("{0}=>Part.extrude({1})".format(trace * ' ', outer_contour))

	# Check argument types:
	assert isinstance(comment, str)
	none_type = type(None)
	if type(inner_contours) == none_type:
	    inner_contours = []
	if type(material) == none_type:
	    material = self._material
	    if type(material) == none_type:
		material = Material("plastic", "abs")
	if type(color) == none_type:
	    color = self._color
	    if type(color) == none_type:
		color = Color()
	assert isinstance(outer_contour, Contour)
	assert isinstance(inner_contours, list)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(trace, int)

	# Do contour adjustments:
	zero = L()
	ezcad = self._ezcad
	adjust = ezcad._adjust

	#print("Part.extrude: adjust = {0}".format(adjust))
	if adjust != zero:
	    outer_contour = outer_contour.adjust(
	      delta = -adjust, start = start, end = end)
	    for index in range(len(inner_contours)):
		inner_contour = inner_contours[index]
		assert isinstance(inner_contour, Contour)
		inner_contour = inner_contour.adjust(
		  delta = -adjust, start = start, end = end)
		inner_contours[index] = inner_contour

	# Record the *color* and *material*:
	self._material = material
	self._color = color

	extrude_axis = end - start
	height = extrude_axis.length()
	place = Place(part = None, name = comment, center = center,
	  axis = axis, rotate = rotate, translate = translate)

	# Figure out what axis is the extrude axis:
	extrude_x = extrude_axis.x.absolute()
	extrude_y = extrude_axis.y.absolute()
	extrude_z = extrude_axis.z.absolute()
	is_x_axis_extrude = False
	is_y_axis_extrude = False
	is_z_axis_extrude = False
	if extrude_x > extrude_y and extrude_x > extrude_z:
	    is_x_axis_extrude = True
	if extrude_y > extrude_x and extrude_y > extrude_z:
	    is_y_axis_extrude = True
	else:
	    is_z_axis_extrude = True

	# Compute bounding box:
	zero = L()
	bounding_box = outer_contour.bounding_box_compute(start, end)
	self._bounding_box_update(bounding_box,
	  comment, place, trace = trace + 1)
	if trace >= 0:
	    print("{0} Part.extrude:bounding_box={0}".
	     format(trace * ' ', bounding_box))

	self._is_part = True
	ezcad = self._ezcad
	assert isinstance(ezcad, EZCAD3)

	if ezcad._mode == EZCAD3.CNC_MODE:
	    if trace >= 0:
		print("{0}Part.extrude:manufacture".format(trace * ' '))
	    # Set up the transform and extrude:
            scad_union_lines = self._scad_union_lines
	    pad = " " * 6

	    # Perform the extrusion in the correct direction:
	    if is_x_axis_extrude:
		# Extrude in X after a 90 degree flip around Y axis:
		scad_union_lines.append("{0}// X axis extrude".format(pad))
		scad_union_lines.append(
		  "{0}translate([{1:m}, 0, 0])".format(pad, end.x))
		scad_union_lines.append(
		  "{0}rotate(a={1}, v=[0, 1, 0])".format(pad, Angle(deg=-90)))
	    elif is_y_axis_extrude:
		# Extrude in Y after a 90 degree flip around X axis:
		scad_union_lines.append("{0}// Y axis extrude".format(pad))
		scad_union_lines.append(
		  "{0}translate([0, {1:m}, 0])".format(pad, end.y))
		scad_union_lines.append(
		  "{0}rotate(a={1}, v=[1, 0, 0])".format(pad, Angle(deg=90)))
	    else:
		# Extrude in Z:
		scad_union_lines.append("{0}// Z axis extrude".format(pad))
		scad_union_lines.append(
		  "{0}translate([0, 0, {1:m}])".format(pad, start.z))

	    # Now do the linear extrude operation:
            scad_union_lines.append(
	      "{0}linear_extrude(height = {1})".format(pad, height))
	    
	    # Construct the polygon information using *indexed_points*:
	    indexed_points = Indexed_Points()
	    outer_contour.path_append(
	      indexed_points = indexed_points,
	      extrude_axis = extrude_axis, trace = trace + 1)
	    for inner_contour in inner_contours:
		inner_contour.path_append(
		  indexed_points = indexed_points,
		  extrude_axis = extrude_axis)
	    indexed_points.polygon_append(scad_union_lines, 6)

	if trace >= 0:
	    print("{0}<=Part.extrude()".format(trace * ' '))
	    
    def hole(self, comment = "NO_COMMENT", diameter = None,
      start = None, end = None, sides = -1, sides_angle = Angle(),
      top = "t", flags = "", start_diameter = None, end_diameter = None,
      trace = -10000):
	""" Part construct: Make a {diameter} hole in {part} with starting
	    at {start_point} and ending at {end_point}.  {comment} will
	    show in any error messages and any generated G-code.  The
	    allowed flag letters in {flags} are:

	      One of 't' (default), 'f', or 'p':
		't'	through hole (i.e. Through)
		'f'	flat hole (i.e. Flat)
		'p'	tip hole (drill tip stops at {end_point} (i.e. tiP)

	      Allowed additional flags:
		'u' upper hole edge should be chamfered (i.e. Upper)
		'l' lower hole edge should be chamfered (i.e. Lower)
		'm' hole is to be milled (i.e. Milled)
	"""

	if trace >= 0:
	    print("{0}=>Part.hole(d={1}, st={2}, end={3}, std={4} ed={5})".
	      format(' ' * trace, diameter, start, end,
	       start_diameter, end_diameter))

	# Check argument types:
	none_type = type(None)
	assert isinstance(comment, str)
	assert isinstance(sides, int)
	assert type(diameter) == none_type or isinstance(diameter, L)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(flags, str)
	assert isinstance(top, str)
	assert isinstance(sides, int)
	assert isinstance(sides_angle, Angle)
	assert type(start_diameter) == none_type or \
	  isinstance(start_diameter, L)
	assert type(end_diameter) == none_type or isinstance(end_diameter, L)
	assert isinstance(diameter, L) or \
	  (isinstance(start_diameter, L) and isinstance(end_diameter, L))

	flat_hole = flags.find("f") >= 0
	through_hole = not flat_hole
	show = flags.find("s") >= 0

	if self._ezcad._mode == EZCAD3.CNC_MODE:
            axis  = start - end
            normalized = axis.normalize()
	    start = start + normalized
	    if through_hole:
		end = end - normalized

	lines = self._scad_difference_lines
	if show:
	    lines = self._scad_union_lines

	ezcad = self._ezcad
	adjust = ezcad._adjust
	self._cylinder(lines = lines, indent = 4, is_solid = False,
	  comment = comment, material = None, color = None, adjust = -adjust,
	  diameter = diameter, start = start, end = end, sides = sides,
	  sides_angle = sides_angle, start_diameter = start_diameter,
	  end_diameter = end_diameter, trace = trace + 1)

	if trace >= 0:
	    print("{0}<=Part.hole(d={1}, st={2}, end={3}, std={4} ed={5})".
	      format(' ' * trace, diameter, start, end,
	       start_diameter, end_diameter))


    def invisible_set(self, invisible = True):
	""" {Part}: Make part visible/invisible. """

	assert isinstance(invisible, bool)
	self._visible = not invisible

    def process(self, ezcad):
	""" *Part*: Perform all the manufacturing processing for a *Part*
	    object (i.e. *self*) and all its sub-*Part*'s.
	"""

	# Verify argument types:
	assert isinstance(ezcad, EZCAD3)

	# Use *part* instead of *self*:
	part = self

        # For debugging, set *debug* to *True*:
	debug = False
	debug = True
	if debug:
	    print("=>Part.process('{0}')".format(part._name))

	# Do the dimensions propogate phase:
	ezcad._mode = EZCAD3.DIMENSIONS_MODE
	ezcad._update_count = 0
	changed = 1
	while changed != 0:
	    # Find all the child *Part*'s:
	    print("Dimensions update {0}".format(ezcad._update_count + 1))
	    changed = part._dimensions_update(ezcad, -1000000)
	    #changed = part._dimensions_update(ezcad, 0)
	    #part._bounding_box_check(0)
	    print("Part.process: {0} dimension(s) changed\n".format(changed))
	    ezcad._update_count += 1

	print("*************Part:{0} bounding_box={1}".format(self._name, self._bounding_box))

	# Now visit *part* and all of its children in CNC mode:
	ezcad._mode = EZCAD3.CNC_MODE
	shop = ezcad._shop
	shop._cnc_generate_set(True)
	part._manufacture(ezcad)
	shop._cnc_generate_set(False)
	ezcad._update_count += 1

	# Now visit *part* and all of its children in STL mode:
	ezcad._mode = EZCAD3.STL_MODE
	part._manufacture(ezcad)
	ezcad._update_count += 1

	# Now visit *part* and all of its children in visualization mode:
	ezcad._mode = EZCAD3.VISUALIZATION_MODE
	part._manufacture(ezcad)
	ezcad._update_count += 1

	if debug:
	    print("<=Part.process('{0}')".format(part._name))

    def reposition(self, vice_y):
	""" *Part*:
	"""

	# This routine will inform {part} that the parts have been
	# repositioned in the vice.  {vice_y} specifies the new
	# location of the vice edge from the origin.

	#call d@(form@("=>reposition@Part(%v%, %i%)\n\") %
	#  f@(part.name) / f@(vice_y))

	# Verify argument types:
	assert isinstance(vice_y, L)

	# Use *part* instead of *self*:
	part = self

	part.cnc_fence()

	part._position_count += 1
	part._vice_y = vice_y
	part._cnc_drill_count = 0

    def rotate(self, nx, ny, nz, angle, vice_y):
	""" *Part*: Rotate the *Part* object (i.e.e *self*) the vector
	    defind (*nx*, *ny*, *nz*) by *angle*, where *vice_y* is distance
	    from the origin to the vice edge.
	"""

	# Verify argument types:
	assert isinstance(nx, L)
	assert isinstance(ny, L)
	assert isinstance(nz, L)
	assert isinstance(angle, Angle)
	assert isinstance(vice_y, L)

	# Use *part* instead of *self*:
	part = self

	# Clear out any pending operations:
	part.cnc_fence()

	# Get a temporary *matrix*:
	matrix = part._shop._matrix_get()

	# Load *matrix* up with rotate coefficients with {angle} negated.
	matrix.rotate_store(nx, ny, nz, angle)

	# Update the *resposition*:
	reposition = part._reposition
	#call multiply_store@(reposition, matrix, reposition)
	reposition.multiply_store(reposition, matrix)

	position = part._position
	matrix.rotate_store(nx, ny, nz, angle)
	position.multiply_store(matrix, position)

	# Verify that we have an identity matrix:
	#matrix.multiply_store(position, reposition)
	#print("pos * repos=\n\{0}\n".format(matrix))

	# Remember that we have repositioned:
	part._position_count += 1

	# Update the position of the vice edge:
	part.vice_y = vice_y

    def translate(self, center):
	""" *Part*: Move the origin of the *Part* object (i.e. *self*) by
	    *center*.
	"""

	# Verify argument types:
	assert isinstance(center, P)

	# Use *part* instead of *self*:
	part = self

	# Mark all previous operations to be done before this reposition:
	part.cnc_fence()

	# Extract some values from {part}:
	matrix = self._shop._matrix_get()
	position = self._position
	reposition = self._reposition

	# Add positive translation to {position}:
	matrix.translate_store(center.x, center.y, center.z)
	position.multiply_store(matrix, position)

	# Add negateive translation to {reposition}:
	# Note that {dx}, {dy}, {dz} are negated, *and* the matrices
	# are multiplied in the opposite order as well:
	matrix.translate_store(-center.x, -center.y, -center.z)
	reposition.multiply_store(reposition, matrix)

	# Remember that we have adjusted the position:
	part._position_count += 1

    def virtual_box(self, comment = "no comment",
      corner1 = None, corner2 = None,
      center = None, axis = None, rotate = None, translate = None):
	""" {Part} vertual_box: Create a virtual box at {corner1} and
	    {corner2}. """

	# Deal with argument defaults:
	none_type = type(None)

	if type(corner1) == none_type:
	    zero = L(0.0)
	    corner1 = P(zero, zero, zero)
	if type(corner2) == none_type:
	    one = L.mm(1.0)
	    corner2 = P(one, one, one)

	# Check argument types:
	assert isinstance(comment, str)
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)
	assert type(center) == none_type or isinstance(center, P)
	assert type(axis) == none_type or isinstance(axis, P)
	assert type(rotate) == none_type or isinstance(rotate, Angle)
	assert type(translate) == none_type or isinstance(translate, P)

	# Make sure that the corners are diagonal from bottom south west
	# to top north east:
	x1 = min(corner1.x, corner2.x)
	x2 = max(corner1.x, corner2.x)
	y1 = min(corner1.y, corner2.y)
	y2 = max(corner1.y, corner2.y)
	z1 = min(corner1.z, corner2.z)
	z2 = max(corner1.z, corner2.z)
	#print("Part.box:{0:m}:{1:m},{2:m}:{3:m},{4:m}:{5:m}". \
	#  format(x1, x2, y1, y2, z1, z2))

	place = Place(part = None, name = comment, center = center,
	  axis = axis, rotate = rotate, translate = translate)
	forward_matrix = place._forward_matrix

	tne = P(x2, y2, z2)
	bsw = P(x1, y1, z1)
	tsw = P(x1, y1, z2)
	bnw = P(x1, y2, z1)
	tnw = P(x1, y2, z2)
	bse = P(x2, y1, z1)
	tse = P(x2, y1, z2)
	bne = P(x2, y2, z1)

	#print("before box={0:m}".format(box))
	self._box_point_update(comment + "[TNE]",
	  forward_matrix.point_multiply(tne))
	self._box_point_update(comment + "[TNW]",
	  forward_matrix.point_multiply(tnw))
	self._box_point_update(comment + "[TSE]",
	  forward_matrix.point_multiply(tse))
	self._box_point_update(comment + "[TSW]",
	  forward_matrix.point_multiply(tsw))
	self._box_point_update(comment + "[BNE]",
	  forward_matrix.point_multiply(bne))
	self._box_point_update(comment + "[BNW]",
	  forward_matrix.point_multiply(bnw))
	self._box_point_update(comment + "[BSE]",
	  forward_matrix.point_multiply(bse))
	self._box_point_update(comment + "[BSW]",
	  forward_matrix.point_multiply(bsw))
	#print("after box={0:m}".format(box))

    def wrl_write(self, wrl_file, indent = 0, parts_table = {}, file_name = ""):
	""" *Part*: Write *self* to *wrl_file*. """
	# Check argument types:
	assert isinstance(wrl_file, file)
	assert isinstance(indent, int)
	assert isinstance(parts_table, dict)
	assert isinstance(file_name, str)

	# Make sure the top level starts with an empty *parts_table*:
	if indent == 0:
	    parts_table = {}
	    wrl_file.write("#VRML V2.0 utf8\n")
	    wrl_file.write("# Generated by EZCAD3\n")
	    wrl_file.write("Viewpoint { description \"Initial view\" position 0 0 150 }\n")
	    wrl_file.write("NavigationInfo { type \"EXAMINE\" }\n")

	# Do some preparation work:
	name = self._name
	spaces = " " * indent

	#print("{0}=>Part.wrl_write({1}, {2}, {3}, {4}):enter".
	#  format(spaces, name, indent, parts_table.keys(), file_name))

	# Figure out whether to generate USE or DEF:
	if self._visible:
	    if name in parts_table:
		wrl_file.write("{0}USE x{1}\n".format(spaces, name))
	    else:
		# Remember that we have defined *self* in the .wrl file:
		parts_table[name] = self
    
		# Decide whether we are a part or an assembly:
		if self._is_part and self._color.alpha > 0.0:
		    # Read in the .stl file that was generated by OpenSCAD:
		    name = self._name
		    ezcad = self._ezcad_get()
		    directory = ezcad._directory_get()
		    stl_file_name = os.path.join(
		      directory, "{0}_{1}.stl".format(name, self._signature_hash))
		    stl_file = open(stl_file_name, "r")
		    stl_lines = stl_file.readlines()
		    stl_file.close()
    
		    # Extract the *triangles* and *vertices* from the read
		    # in content:
		    triangles = []
		    offsets = []
		    vertices = {}
		    size = len(stl_lines)
		    assert stl_lines[0][:5] == "solid"
		    index = 1
		    while index + 4 < size:
			# Extract *vertex1*, *vertex2*, and *vertex3*:
			xlist = stl_lines[index + 2].split()
			vertex1 = \
			  (float(xlist[1]), float(xlist[2]), float(xlist[3]))
			xlist = stl_lines[index + 3].split()
			vertex2 = \
			  (float(xlist[1]), float(xlist[2]), float(xlist[3]))
			xlist = stl_lines[index + 4].split()
			vertex3 = \
			  (float(xlist[1]), float(xlist[2]), float(xlist[3]))
    
			# Get *offset1* for *vertex1*:
			if vertex1 in vertices:
			    offset1 = vertices[vertex1]
			else:
			    offset1 = len(offsets)
			    offsets.append(vertex1)
			    vertices[vertex1] = offset1
    
			# Get *offset2* for *vertex2*:
			if vertex2 in vertices:
			    offset2 = vertices[vertex2]
			else:
			    offset2 = len(offsets)
			    offsets.append(vertex2)
			    vertices[vertex2] = offset2
    
			# Get *offset3* for *vertex3*:
			if vertex3 in vertices:
			    offset3 = vertices[vertex3]
			else:
			    offset3 = len(offsets)
			    offsets.append(vertex3)
			    vertices[vertex3] = offset3
    
			# Create a triangle using the offsets:
			triangles.append( (offset1, offset2, offset3,) )
			index += 7
		    # We are done with *stl_lines*:
		    stl_lines = None
		    spaces = " " * indent
    
		    # Write out "DEF name Shape {":
		    wrl_file.write(
		      "{0}DEF x{1} Shape {{\n".format(spaces, name))
    
		    # Output appearance *color* and material properties::
		    color = self._color
		    assert isinstance(color, Color)
		    wrl_file.write(
		      "{0} appearance Appearance {{\n".format(spaces))
		    wrl_file.write(
		      "{0}  material Material {{\n".format(spaces))
		    wrl_file.write(
		      "{0}   diffuseColor {1} {2} {3}\n".
		      format(spaces, color.red, color.green, color.blue))
		    #wrl_file.write(
		    #  "{0}   transparency {1}\n".format(spaces, color.alpha))
		    wrl_file.write(
		      "{0}  }}\n".format(spaces))
		    wrl_file.write(
		      "{0} }}\n".format(spaces))
	    
		    # Start "geometry IndexedFaceSet {...}":
		    wrl_file.write(
		      "{0} geometry IndexedFaceSet {{\n".format(spaces))
    
		    # Output the *vertices* in a "Coordinate {...}":
		    wrl_file.write(
		      "{0}  coord Coordinate {{\n".format(spaces))
		    wrl_file.write(
		      "{0}   point[\n".format(spaces))
		    for offset in offsets:
			wrl_file.write(
			  "{0}    {1} {2} {3}\n".
			  format(spaces, offset[0], offset[1], offset[2]))
		    wrl_file.write(
		      "{0}   ]\n".format(spaces))
		    wrl_file.write(
		      "{0}  }}\n".format(spaces))
    
		    # Output the *triangles* in a "coordIndex [...]"::
		    wrl_file.write(
		      "{0}  coordIndex [\n".format(spaces))
		    for triangle in triangles:
			# Output each *triangle* (except last one) with
			# trailing comma:
			wrl_file.write(
			  "{0}   {1} {2} {3} -1\n".
			  format(spaces, triangle[0], triangle[1], triangle[2]))
		    wrl_file.write(
		      "{0}  ]\n".format(spaces))
    
		    # Close "geometry IndexdFaceSet{...}":
		    wrl_file.write(
		      "{0} }}\n".format(spaces))
		    # Close "... Shape {...}":
		    wrl_file.write(
		      "{0}}}\n".format(spaces, name))
		else:
		    # Start a named "Group {...}":
		    wrl_file.write(
		      "{0}DEF x{1} Group {{\n".format(spaces, self._name))
		    wrl_file.write(
		      "{0} children [\n".format(spaces))
    
		    # Output each *place* in *places*
		    sub_parts = []
		    for attribute_name in dir(self):
			if not attribute_name.startswith("_") and \
			  attribute_name.endswith("_"):
			    sub_part = getattr(self, attribute_name)
			    assert isinstance(sub_part, Part)
			    sub_parts.append(sub_part)
		    sub_parts = list(set(sub_parts))
		    sub_parts.sort(key = lambda sub_part: sub_part._name)
    
		    for sub_part in sub_parts:
			# Extract some values from *place*:
			center = sub_part._center
			axis = sub_part._axis
			rotate = sub_part._rotate
			translate = sub_part._translate
    
			# Figure out if we have to do a "Transform...":
			none_type = type(None)
			zero_rotate = type(rotate) == none_type or rotate == Angle()
			zero_translate = \
			  type(translate) == none_type or translate == P()
			if zero_rotate and zero_translate:
			    # We have neither a rotation nor a translation;
			    # so we output *part* without a "Transform..."
			    sub_part.wrl_write(wrl_file,
			      indent + 2, parts_table, file_name)
			else:
			    # We have either a rotation and/or a translation;
			    # So we need to wrap *part* in a "Transform ...":
			    wrl_file.write(
			      "{0}  Transform {{\n".format(spaces))
    
			    # If appropriate, write out "rotation"
			    if not zero_rotate:
				# Move rotation center if not (0,0,0):
				if type(center) != none_type and center != P():
				    wrl_file.write(
				      "{0}   center {1} {2} {3}\n".
				      format(spaces, center.x, center.y, center.z))
				# Write out the rotation:
				wrl_file.write(
				      "{0}   rotation {1} {2} {3} {4:r}\n".format(
				      spaces, axis.x, axis.y, axis.z, rotate))
    
			    # If appropriate, write out the "translation ..."
			    if not zero_translate:
				wrl_file.write(
				  "{0}   translation {1} {2} {3}\n".format(spaces,
				  translate.x, translate.y, translate.z))
    
			    # Now we can write out *part* wrapped in
			    # "children [...]":
			    wrl_file.write(
			      "{0}   children [\n".format(spaces))
			    sub_part.wrl_write(wrl_file,
			      indent + 4, parts_table, file_name)
			    wrl_file.write(
			      "{0}   ]\n".format(spaces))
    
			    # Close out "Transform {...}":
			    wrl_file.write(
			      "{0}  }}\n".format(spaces))
    
		    # We are done writing out *places*, so we can close out
		    # "children [ ...]":
		    wrl_file.write(
		      "{0} ]\n".format(spaces))
    
		    # Close out "Group { ... }":
		    wrl_file.write(
		      "{0}}}\n".format(spaces))
    
	    #print("{0}<=Part.wrl_write({1}, {2}, {3}, {4}):leave".
	    #  format(spaces, name, indent, parts_table.keys(), file_name))
    
    # ===================================
    def __str__(self):
	""" Part: Return {self} as a formatted string. """

	return str(self._name)

    def angle(self, angle_path):
	""" Part dimensions: Return the {Angle} associated with {angle_path}
	    starting from {self}. """
	
	return self.value_lookup(angle_path, Part.ANGLES)

    def angle_set(self, angle_name, angle_value):
	""" Part dimensions: Store an {Angle} named {angle_name} into {self}
	    with a value of  {angle_value}.  In addtion, {angle_value} is
	    returned. """

	return self.value_set(angle_name, angle, Part.ANGLES)

    def boundary_trim(self, comment, corner_radius, flags):
	""" Part construction: Using the last selected top surface
	    from {vice_position}() trim the 4 corners with a radius of
	    {corner_radius}. """

	#print "=>Part.boundary_trim({0}, '{1}', {2}, '{3}')". \
	#  format(self._name, comment, corner_radius, flags)

	ezcad = self._ezcad
	xml_stream = ezcad._xml_stream
	if xml_stream != None:
	    # Grab the 6 surfaces:
	    t = self.t
	    b = self.b
	    n = self.n
	    s = self.s
	    e = self.e
	    w = self.w

	    # Grab the 8 corners:
	    tne = self.tne
	    tnw = self.tnw
	    tse = self.tse
	    tsw = self.tsw
	    bne = self.bne
	    bnw = self.bnw
	    bse = self.bse
	    bsw = self.bsw

	    # Compute {x_extra}, {y_extra}, {z_extra}:
	    x_extra = (w.x - self.xw.x).maximum(self.xe.x - e.x)
	    y_extra = (s.y - self.xs.y).maximum(self.xn.y - n.y)
	    z_extra = (b.z - self.xb.z).maximum(self.xt.z - t.z)

	    assert self._top_surface_set, \
	      "Top surface for part '{0}' is not set".format(self._name)
	    top_surface = self._top_surface
	    if top_surface == t:
		# Top surface:
		self.corner(comment + ":TNW", tnw, corner_radius)
		self.corner(comment + ":TNE", tne, corner_radius)
		self.corner(comment + ":TSE", tse, corner_radius)
		self.corner(comment + ":TSW", tsw, corner_radius)
		extra = x_extra.maximum(y_extra)
		self.contour(comment, t, b, extra, flags)
	    elif top_surface == b:
		# Bottom surface:
		self.corner(comment + ":BNW", bnw, corner_radius)
		self.corner(comment + ":BSW", bsw, corner_radius)
		self.corner(comment + ":BSE", bse, corner_radius)
		self.corner(comment + ":BNE", bne, corner_radius)
		extra = x_extra.maximum(y_extra)
		self.contour(comment, b, t, extra, flags)
	    elif top_surface == e:
		# East surface:
		self.corner(comment + ":BNE", bne, corner_radius)
		self.corner(comment + ":BSE", bse, corner_radius)
		self.corner(comment + ":TSE", tse, corner_radius)
		self.corner(comment + ":TNE", tne, corner_radius)
		extra = y_extra.maximum(z_extra)
		self.contour(comment, e, w, extra, flags)
	    elif top_surface == w:
		# West surface:
		self.corner(comment + ":BNW", bnw, corner_radius)
		self.corner(comment + ":BSW", bsw, corner_radius)
		self.corner(comment + ":TSW", tsw, corner_radius)
		self.corner(comment + ":TNW", tnw, corner_radius)
		extra = y_extra.maximum(z_extra)
		self.contour(comment, w, e, extra, flags)
	    elif top_surface == n:
		# North surface:
		self.corner(comment + ":TNE", tne, corner_radius)
		self.corner(comment + ":TNW", tnw, corner_radius)
		self.corner(comment + ":BNW", bnw, corner_radius)
		self.corner(comment + ":BNE", bne, corner_radius)
		extra = x_extra.maximum(z_extra)
		self.contour(comment, n, s, extra, flags)
	    elif top_surface == s:
		# South surface:
		self.corner(comment + ":TSE", tse, corner_radius)
		self.corner(comment + ":TSW", tsw, corner_radius)
		self.corner(comment + ":BSW", bsw, corner_radius)
		self.corner(comment + ":BSE", bse, corner_radius)
		extra = x_extra.maximum(z_extra)
		self.contour(comment, s, n, extra, flags)
	    else:
		assert False, \
		  "Top surface for {0} is {1} which is not $N/$S/$E/$W/$T/$B". \
		  format(self.name, top_surface)

    def block_diagonal(self, diagonal, color, material):
	""" Part construct: Turn {self} into a block centered on the origin
	    with {diagonal} stretching equally across to the two opposite
	    corners.  The block is made of {material} and visualized
	    as {color}. """

	corner2 = diagonal.half()
	corner1 = -corner2
	self.block_corners(corner1, corner2, color, material)

    def chamfers(self, upper_chamfer, lower_chamfer):
	""" {Part} construct: Set the chamfers for contours and simple
	    pockets.  {upper_chamfer} specifies the chamfer of the
	    upper edge and the {lower_chamfer} specifies the chamer
	    of the lower edge.  Zero disables the chamfer for the
	    specified edge. """

	# Check argument types:
	assert isinstance(upper_chamfer, L)
	assert isinstance(lower_chamfer, L)

	# Extract some values from {ezcad}:
	ezcad = self._ezcad
	xml_indent = ezcad._xml_indent
	xml_stream = ezcad._xml_stream

	# Output <Chamfers Upper= Lower= />:
	xml_stream.write('{0}<Chamfers Upper="{1}" Lower="{2}"/>\n'. \
	  format(xml_indent, upper_chamfer, lower_chamfer))

    def construct_mode(self):
	""" Part dimensions: Return {True} if {self} should be constructed. """

	result = self.ezcad.construct_mode()
	return result

    def contour(self, comment, start_point, end_point, extra, flags):
	""" Part construct: Cause the current list of corners associated
	    with {self} to be removed starting at a depth of {start_point}
	    and ending at a depth of {end_point} using an exterior contour.
	    {extra} is the amount of extra material being removed. {flags}
	    can contain the letter 'u' for an upper chamfer, 'l' for a lower
	    chamfer, and 't' for a contour that cuts entirely through {self}.
	    {comment} is used in error messages and any generated G-code. """

	# Verify that each flag in {flags} is OK:
	for flag in flags:
	    assert flag == 'u' or flag == 'l' or flag == 't', \
	      'Part.Contour: Bad flag "{0}" in flags "{1}" for part {2}'. \
	      format(flag, flags, self.name)

	# Extract some values from {ezcad}:
	ezcad = self._ezcad
	xml_indent = ezcad._xml_indent
	xml_stream = ezcad._xml_stream

	# Make sure we are in construct mode:
	#assert self.ezcad.construct_mode(), \
	#  "Part.Contour: Called when Part {0} is not in consturct mode". \
	#  format(self.name)

	if xml_stream != None:
	    # Make sure we have some {corners} for the contour:
	    corners = self._corners
	    corners_size = len(corners)
	    assert corners_size >= 3, \
	      "Part '{0}' only has {1} corners which is insufficient". \
	      format(self._name, corners_size)

	    # Now figure out the desired order for {corners}.  Depending upon
	    # what {top_surface} is being used, we reverse the contour
	    # direction:
	    top_surface = self._top_surface
	    reverse = None
	    if top_surface == self.e:
		reverse = False
	    elif top_surface == self.w:
		reverse = True
	    elif top_surface == self.n:
		reverse = False
	    elif top_surface == self.s:
		reverse = True
	    elif top_surface == self.t:
		reverse = False
	    elif top_surface == self.b:
		reverse = False
	    assert reverse != None, \
	      "Part {0} top surface must be one of T, B, N, S, E, or W". \
	      format(self.name)
	    if reverse:
		corners.reverse()

	    # Now output {corners} and clear it for the next contour:
	    for corner in self._corners:
		xml_stream.write(corner)
	    self._corners = []

	    # <Contour SX= SY= SZ= EX= EY= EZ= Extra= Flags= Comment= />:
	    xml_stream.write('{0}<Contour SX="{1}" SY="{2}" SZ="{3}"'. \
	      format(" " * xml_indent,
	      start_point.x, start_point.y, start_point.z))
	    xml_stream.write(' EX="{0}" EY="{1}" EZ="{2}"'. \
	      format(end_point.x, end_point.y, end_point.z))
	    xml_stream.write(' Extra="{0}" Flags="{1}" Comment="{2}"/>\n'. \
	      format(extra, flags, comment))

    def contour_reverse(self):
	""" """

	self._corners.reverse()

    def corner(self, comment, corner_point, radius):
	""" Part construct: Add a corner with a radius of {radius} to
	    {self} using {corner_point} to specify the corner location.
	    {comment} will show up any error messages or generated G-code. """

	# Extract some values from {ezcad}:
	ezcad = self._ezcad
	xml_stream = ezcad._xml_stream
	if xml_stream != None:
	    corner_xml = \
	      ('{0}<Corner Radius="{1}" CX="{2}" CY="{3}" CZ="{4}"' + \
	      ' Comment="{5}"/>\n').format(" " * ezcad._xml_indent,
	      radius, corner_point.x, corner_point.y, corner_point.z, comment)
	    self._corners.append(corner_xml)

    def done(self):
	""" Part (tree_mode): Mark {self} as done.  {self} is
	    always returned. """

	#print "=>Part.done('%s')" % (self.name)

	# Verify that parts stack is properly synchronized:
	if self.ezcad.part_tree_mode():
	    # Verify that the top of the parts stack matches {self}:
	    popped_part = self.ezcad.part_pop()
	    assert popped_part == self, \
	      "Part '%s' does not match popped part '%s'; mismathed done()?" % \
	      (self.name, popped_part.name)

	if self.construct_mode():
	    # Flush any remaining screw holes for current surface:
	    if self.top_surface_set:
		self.screw_holes(True)

       	    # Look for any screws that have not been processed:
	    screw_levels = self._screw_levels
            for screw_level_table in screw_levels.values():
		for screw_level in screw_level_table.values():
		    assert screw_level.done, \
		      "Screw {0} on {1} of part {2} has not been output". \
		      format(screw_level.screw.name, screw_level.surface, \
		      self.name)

	    # Close out this part in XML stream:
	    ezcad = self.ezcad
	    xml_stream = ezcad.xml_stream
	    ezcad.xml_indent_pop()
	    xml_stream.write("{0}</Part>\n".format(ezcad._xml_indent))

	#print "<=Part.done('%s')" % (self.name)
	return self

    def dxf_place(self, dxf_name, x_offset, y_offset):
	""" Part dimensions: Place {self} into the DXF file named {dxf_name}
	    with an offset of {x_offset} and {y_offset}. """

	if self.construct_mode():
	    ezcad = self.ezcad
	    xml_stream = ezcad.xml_stream
	    xml_stream.write(
	      '{0}<DXF_Place DX="{1}" DY="{2}" DXF_Name="{3}"/>\n'. \
	      format(ezcad._xml_indent, x_offset, y_offset, dxf_name))

    def extra_ewnstb(self, east, west, north, south, top, bottom):
	""" {Part} manufacture: Update the extra bounding box for {self}
	    to be extended by {east} and {west} in the X axis, {north} and
	    {south} in the Y axis, and {top} and {bottom} in the Z axis. """

	#print "extra_ewnstb: e={0} w={1} n={2} s={3} t={4} b={5}". \
	#  format(east, west, north, south, top, bottom)

	# Arugment type checking:
	assert isinstance(east, L)
	assert isinstance(west, L)
	assert isinstance(north, L)
	assert isinstance(south, L)
	assert isinstance(top, L)
	assert isinstance(bottom, L)

	# Look up {bsw} and {tne}, the bounding box corner {P}'s:
	bsw = self.bsw
	tne = self.tne

	# Compute the dimensions of extra bounding box:
	b = bsw.z - bottom
	s = bsw.y - south
	w = bsw.x - west
	t = tne.z + top
	n = tne.y + north
	e = tne.x + east

	# Compute the mid-point values for east-west, north-south and
	# top-bottom:
	ew2 = (e + w)/2
	ns2 = (n + s)/2
	tb2 = (t + b)/2

	# Install 6 bounding box surface {P}'s into {points}:
	self.xb = P(ew2, ns2,   b)
	self.xe = P(  e, ns2, tb2)
	self.xn = P(ew2,   n, tb2)
	self.xs = P(ew2,   s, tb2)
	self.xt = P(ew2, ns2,   t)
	self.xw = P(  w, ns2, tb2)

	# Install 12 bounding box edge {P}'s into {points}:
	self.xbe = P(  e, ns2,   b)
	self.xbn = P(ew2,   n,   b)
	self.xbs = P(ew2,   s,   b)
	self.xbw = P(  w, ns2,   b)
	self.xne = P(  e,   n, tb2)
	self.xnw = P(  w,   n, tb2)
	self.xse = P(  e,   s, tb2)
	self.xsw = P(  w,   s, tb2)
	self.xte = P(  e, ns2,   t)
	self.xtn = P(ew2,   n,   t)
	self.xts = P(ew2,   s,   t)
	self.xtw = P(  w, ns2,   t)

	# Install 8 bounding box corner {P}'s into {points}:
	self.xbne = P(e, n, b)
	self.xbnw = P(w, n, b)
	self.xbse = P(e, s, b)
	self.xbsw = P(w, s, b)
	self.xtne = P(e, n, t)
	self.xtnw = P(w, n, t)
	self.xtse = P(e, s, t)
	self.xtsw = P(w, s, t)

	xbsw = self.xbsw
	xtne = self.xtne
	ezcad = self._ezcad
	xml_stream = ezcad._xml_stream
	if xml_stream != None:
	    xml_stream.write( ('{0}<Extra C1X="{1}" C1Y="{2}" C1Z="{3}"' + \
	      ' C2X="{4}" C2Y="{5}" C2Z="{6}"/>\n'). \
	      format(' ' * ezcad._xml_indent,
	      xbsw.x, xbsw.y, xbsw.z, xtne.x, xtne.y, xtne.z))

    def extra_xyz(self, dx, dy, dz):
	""" {Part}: Add some extra material the block of {self} by {dx},
	    {dy}, and {dz} in the X, Y, and Z dimensions. """

	# Argument type checking:
	assert isinstance(dx, L)
	assert isinstance(dy, L)
	assert isinstance(dz, L)

	# Pass everything on to {extra_ewnstb}:
	half_dx = dx/2
	half_dy = dy/2
	half_dz = dz/2
	self.extra_ewnstb(half_dx, half_dx, half_dy, half_dy, half_dz, half_dz)

    def extrusion(self, color, material, kind, start_point, end_point, \
      a_width, a_thickness, b_width, b_thickness, rotate):
	""" Part dimensions: Create a {kind} extrusion manufactured out of
	    {material} that goes from {start_point} to {end_point} rotated
	    by {rotate}.  {a_width}, {a_thickness}, {b_width} and
	    {b_thickness} specify the extra dimentions of the extrusion.
	    The returned part is visualized using {color}. The extrusion
	    must be aligned with one of the X, Y or Z axes. """

	# Check argument types:
	assert isinstance(color, str)
	assert isinstance(material, str)
	assert isinstance(start_point, P)
	assert isinstance(end_point, P)
	assert isinstance(a_width, L)
	assert isinstance(a_thickness, L)
	assert isinstance(b_width, L)
	assert isinstance(b_thickness, L)
	assert isinstance(rotate, Angle)

	# Define soem useful abreviations:
	cosine = Angle.cosine
	inch = L.inch
	sine = Angle.sine
	zero = inch(0)

	# Make sure we are in the right mode:
	assert not self.part_tree_mode(), \
	  "Part.extrusion: Part '{0}' is in part tree mode".format(self.name)

	# Compute {extrusion_axis}, the axis along with the tube is aligned:
	extrusion_axis = end_point - start_point
	extrusion_axis_x = extrusion_axis.x
	extrusion_axis_y = extrusion_axis.y
	extrusion_axis_z = extrusion_axis.z
	#print "start={0} end={1} extrusion_axis={2}". \
	#  format(start_point, end_point, extrusion_axis)

	# Extract coordinates from {start_point} and {end_point}:
	x1 = start_point.x
	y1 = start_point.y
	z1 = start_point.z
	x2 = end_point.x
	y2 = end_point.y
	z2 = end_point.z

	if kind == 'L':
	    # angLe extrusion:
	    # We have angle material:
	    #
	    #    |    |<----- a_width ------>|
	    #    v                             
	    # ------- +----------------------O -----
	    #         |                      |   ^
	    # ------- +-------------------+  |   |
	    #    ^                        |  |   |
	    #    |                        |  |   |
	    # b_thickness                 |  |   |
	    #                             |  | b_width
	    #                             |  |   |
	    #                             |  |   |
	    #                             |  |   |
	    #                             |  |   v
	    #                             +--+ -----
	    #
	    #             a_thickness --->|  |<---
	    #
	    # O = origin
	    # A goes negative from origin
	    # B goes negative from origin
	    a1 = -a_width
	    a2 = zero
	    b1 = zero
	    b2 = -b_width
	elif kind == 'C':
	    # We have channel material:

	    #       ->|  |<-- b_thickness
	    #
	    # ------- +--+                +--+ 
	    #    ^    |  |                |  | 
	    #    |    |  |                |  | 
	    #    |    |  |                |  | 
	    #    |    |  |                |  | 
	    #    |    |  |                |  | a_thickness
	    # b_width |  |                |  |   |
	    #    |    |  |                |  |   V
	    #    |    |  +----------------+  | -----
	    #    V    |                      |
	    # ------- O----------------------+ -----
	    #                                    ^
	    #         |<----- a_width ------>|   |
	    a1 = zero
	    a2 = a_width
	    b1 = zero
	    b2 = b_wdith
	elif kind == 'I':
	    assert False, "I beam not defined"
	elif kind == 'T':
	    assert False, "T beam not defined"
	elif kind == 'X':
	    assert False, "X beam not defined"
	elif kind == 'Z':
	    assert False, "Z beam not defined"
	else:
	    assert False, "Unrecognized extrusion '{0}'".format(kind)

	# Figure out what axis we are aligned on:
	if extrusion_axis_x != zero and \
	  extrusion_axis_y == zero and extrusion_axis_z == zero:
	    # Extrusion aligned on X:
	    # (a1, b1) (a1, b2) (a2, b1) (a2, b2)
	    a1_b1_distance = L.distance2(a1, b1)
	    a1_b1_angle = a1.arc_tangent2(b1) + rotate
	    a1_b1_y = y1 + a1_b1_distance.cosine(a1_b1_angle)
	    a1_b1_z = z1 + a1_b1_distance.sine(a1_b1_angle)
	    #print("extrusion: a1b1d={0} a1b1a={1} a1b1x={2} a1b1y={3}". \
	    #  format(a1_b1_distance, a1_b1_angle, a1_b1_y, a1_b1_z))

	    a1_b2_distance = L.distance2(a1, b2)
	    a1_b2_angle = a1.arc_tangent2(b2) + rotate
            a1_b2_y = y1 + a1_b2_distance.cosine(a1_b2_angle)
            a1_b2_z = z1 + a1_b2_distance.sine(a1_b2_angle)
	    #print("extrusion: a1b2d={0} a1b2a={1} a1b2x={2} a1b2y={3}". \
	    #  format(a1_b2_distance, a1_b2_angle, a1_b2_y, a1_b2_z))

	    a2_b1_distance = L.distance2(a2, b1)
	    a2_b1_angle = a2.arc_tangent2(b1) + rotate
            a2_b1_y = y1 + a2_b1_distance.cosine(a2_b1_angle)
            a2_b1_z = z1 + a2_b1_distance.sine(a2_b1_angle)
	    #print("extrusion: a2b1d={0} a2b1a={1} a2b1x={2} a2b1y={3}". \
	    #  format(a2_b1_distance, a2_b1_angle, a2_b1_y, a2_b1_z))

	    a2_b2_distance = L.distance2(a2, b2)
	    a2_b2_angle = a2.arc_tangent2(b2) + rotate
            a2_b2_y = y1 + a2_b2_distance.cosine(a2_b2_angle)
            a2_b2_z = z1 + a2_b2_distance.sine(a2_b2_angle)
	    #print("extrusion: a2b2d={0} a2b2a={1} a2b2x={2} a2b2y={3}". \
	    #  format(a2_b2_distance, a2_b2_angle, a2_b2_y, a2_b2_z))

	    e = max(x1, x2)
	    w = min(x1, x2)
	    n = max(max(a1_b1_y, a1_b2_y), max(a2_b1_y, a2_b2_y))
	    s = min(min(a1_b1_y, a1_b2_y), min(a2_b1_y, a2_b2_y))
	    t = max(max(a1_b1_z, a1_b2_z), max(a2_b1_z, a2_b2_z))
	    b = min(min(a1_b1_z, a1_b2_z), min(a2_b1_z, a2_b2_z))
	    #print("extrusion: e={0} w={1} n={2} s={3} t={4} b={5}". \
	    #  format(e, w, n, s, t, b))
	elif extrusion_axis_x == zero and \
	  extrusion_axis_y != zero and extrusion_axis_z == zero:
	    # Extrusion aligned on Y:
	    assert False, "Fix Y alignment"

	    e = max(xa1, xa2)
	    w = min(xa1, xa2)
	    n = max(y1, y2)
	    s = min(y1, y2)
            t = max(zb1, zb2)
            b = min(zb1, zb2)
	elif extrusion_axis_x == zero and \
	  extrusion_axis_y == zero and extrusion_axis_z != zero:
	    # Extrusion aligned on Z:
	    assert False, "Fix Z alignment"

	    e = max(xa1, xa2)
	    w = min(xa1, xa2)
	    n = max(yb1, yb2)
	    s = min(yb1, yb2)
            t = max(z1, z2)
            b = min(z1, z2)
	else:
	    assert not self.construct_mode(), \
 	      "Angle extrusion is not aligned with X, Y, or Z axis"
	    t = z2
            b = z1
	    n = y2
	    s = y1
	    w = x1
	    e = x2

	# Record the angle corners in {self}:
	self.point_xyz_set(kind + "_extrusion_corner_bsw", w, s, b)
	self.point_xyz_set(kind + "_extrusion_corner_tne", e, n, t)

	# Record the material in {self}:
	self._material = material

	if self.dimensions_mode():
            self.bounding_box_update(w, s, b, e, n, t)

	# In consturct mode, output the <Block ...> line:
	if self.construct_mode():
	    ezcad = self.ezcad
	    xml_stream = ezcad.xml_stream
	    xml_stream.write('{0}<Extrusion'.format(ezcad._xml_indent))
	    xml_stream.write(' Kind="{0}"'.format(kind))
	    xml_stream.write(' SX="{0}" SY="{1}" SZ="{2}"'.format(x1, y1, z1))
	    xml_stream.write(' EX="{0}" EY="{1}" EZ="{2}"'.format(x2, y2, z2))
	    xml_stream.write(' A_Width="{0}" A_Thickness="{1}"'. \
	      format(a_width, b_thickness))
	    xml_stream.write(' B_Width="{0}" B_Thickness="{1}"'. \
	      format(b_width, b_thickness))
	    xml_stream.write(' Rotate="{0}"'.format(rotate))
	    xml_stream.write(' Color="{0}" Transparency="{1}" Material="{2}"'. \
	      format(color, self.transparency, material))
	    xml_stream.write(' Comment="{0}"/>\n'.format(self.name))

    def hole_through(self, comment, diameter, start_point, flags, \
      countersink_diameter = L(0.0)):
	""" Part construct: Make a {diameter} hole in {self} with starting
	    at {start_point} and going all the way through {self}.
	    {comment} will show in any error messages and any generated
	    G-code.  The allowed flag letters in {flags} are:

	      Zero, one or more of the following:
		't' through hole (i.e. Through)
		'u' upper hole edge should be chamfered (i.e. Upper)
		'l' lower hole edge should be chamfered (i.e. Lower)
	"""

	ezcad = self._ezcad
	xml_stream = ezcad._xml_stream

	if xml_stream != None:

	    # Define some useful abbreviations:
	    inch = L.inch
	    zero = inch(0)

	    # Extract some values from {ezcad}:
	    ezcad = self._ezcad
	    xml_indent = ezcad._xml_indent

	    # Write out <Hole_Through Diameter= SX= SY= SZ= Flags= Comment= />:
	    xml_stream.write('{0}<Hole_Through Diameter="{1}"'. \
	      format(" " * xml_indent, diameter))
	    xml_stream.write(' Countersink_Diameter="{0}"'. \
	      format(countersink_diameter))
	    xml_stream.write(' SX="{0}" SY="{1}" SZ="{2}"'. \
	      format(start_point.x, start_point.y, start_point.z))
	    xml_stream.write(' Flags="{0}" Comment="{1}"/>\n'. \
	      format(flags, comment))

    def length(self, length_path):
	""" Part dimensions: Return the {Angle} associated with {length_path}
	    starting from {self}. """
	
	return self.value_lookup(length_path, Part.LENGTHS)

    def length_set(self, length_name, length_value):
	""" Part dimensions: Store a {L} named {length_name} into {self}
	    with a value of {length_value}.  In addtion, {length_value} is
	    returned. """

	return self.value_set(length_name, length_value, Part.LENGTHS)

    def part(self, part_path):
	""" Part dimensions: Return the {Part} associated with {part_pat}
	    starting from {self}. """

	path_parts = part_path.split('/')
	part = self
	for path_part in path_parts:
	    if path_part == "..":
		assert part.parent != None, \
		  "Path %s goes past root" % (part_path)
		part = part.parent
	    else:
		parts = part.parts
		assert path_part in parts, \
		  "No part named '{0}' in path '{1}' starting from part {2}". \
		  format(path_part, part_path, self.name)
		part = parts[path_part]
	return part

    def part_tree_mode(self):
	""" Part internal: Return {True} if {self} is in part tree mode. """

	return self.ezcad.part_tree_mode()

    def path_parts(self):
	""" Part internal: Return a list the contains the sequence of
	    {Part}'s from the root {Part} to {self}.  {self} is always the
	    last item on the returned list. """

	parts = []
	part = self
	while part.parent != None:
	    parts.append(part)
	    part = part.parent

	parts.reverse()
	#print "parts=", parts
	return parts

    def place(self,
      center = None, axis = None, rotate = None, translate = None):
	""" *Part*: Place *self* at *translate* rotated by *rotate* around
	    *axis* centered on *center*. """

	#print("=>Part.place({0}, part='{1}',name='{2}' ...)". \
	#  format(self._name, part._name, name))

	# Check argument types:
	none_type = type(None)
	assert type(center) == none_type or isinstance(center, P)
	assert type(axis) == none_type or isinstance(axis, P)
	assert type(rotate) == none_type or isinstance(rotate, Angle)
	assert type(translate) == none_type or isinstance(translate, P)

	# Load up *self8:
	self._center = center
	self._axis = axis
	self._rotate = rotate
	self._translate = translate

	# Create *place* and stuff into *_places*:
	#place = Place(part = part, name = name,
	#  center = center, axis = axis, rotate = rotate, translate = translate)
	#self._places[name] = place

	#print("<=Part.place({0}, part='{1}',name='{2}' ...)". \
	#  format(self._name, part._name, name))

    def xxxno_automatic_place(self):
	""" *Part*: Disable automatic placement of *self*. """

	name = self._name
	places = self.up._places
	if name in places:
	    del places[name]

    #def point(self, point_path):
    #	""" Part dimensions: Return the {P} associated with {point_path}
    #	    starting from {self}. """
    #
    #	return self.value_lookup(point_path, Part.POINTS)

    def point_map(self, place_path):
	""" Part dimensions: Return the {P} that corresponds to
	    the point specified by {place_path} in the {self} reference
	    frame. """

	assert not self.part_tree_mode(), \
	  "Part.point_subtract: Not allowed in part tree mode"

	# Wait until after all part names have been defined:
	if self.dimensions_define_mode():
	    zero = L.inch(0)
	    result = P(self, zero, zero, zero)
	else:
	    # Look up the {Part}/{P}/{Matrix} triple for {place_path1}:
	    part_point_matrix = self.point_subtract_helper(place_path)
	    part = part_point_matrix[0]
	    point = part_point_matrix[1]
	    matrix = part_point_matrix[2]

	    # Perform the final subtraction and refernce to {part2}:
	    result = matrix.point_multiply(point, self)
	    #print "result={0}".format(result)

	return result

    def relative_path(self, to_part):
	""" Part internal: Return the relative path from {self}
	    to {to_part}. """

	#print "Part.relative_path('{0}', '{1}')". \
	#  format(self.name, to_part.name)

	from_path_parts = self.path_parts()
	to_path_parts = to_part.path_parts()

	#print "from_path_parts=", from_path_parts
	#print "to_path_parts=", to_path_parts

	while len(from_path_parts) != 0 and len(to_path_parts) != 0 and \
	  from_path_parts[0] == to_path_parts[0]:
	    from_path_parts.pop(0)
	    to_path_parts.pop(0)
	
	separator = ''
	path = ""
	while len(from_path_parts) != 0:
	    from_path_parts.pop()
	    path = path + separator + ".."
	    separator = '/'

	while len(to_path_parts) != 0:
	    part = to_path_parts.pop()
	    path = path + separator + part.name
	    separator = '/'

	return path

    def rod(self, color, material, start_point, end_point, diameter, sides):
	""" Part dimensions: Create a rod for {self} out of {material} that
	    goes from {start_point} to {end_point} and is {diameter} round.
	    The rod is visualized with {color}.  {sides} specifies how many
	    sides to use for visualization; setting {sides} to zero gives a
	    reasonable visualization.  The rod must be aligned with
	    one of the X, Y or Z axes. """

	self.tube(color, material, start_point, end_point, diameter, \
	  diameter.half() - L.inch(.0001), sides)

    def _scad_color(self, lines, color, indent = 0):
	if isinstance(color, Color):
	    red = color.red
	    green = color.green
	    blue = color.blue
	    alpha = color.alpha
	    spaces = " " * indent
	    if alpha < 1.0:
		lines.append("{0}color([{1},{2},{3},{4}])".
	          format(spaces, red, green, blue, alpha))
	    else:
		lines.append("{0}color([{1},{2},{3}])".
	          format(spaces, red, green, blue))

    def _scad_transform(self, lines, indent = 0,
       center = None, axis = None, rotate = None, translate = None):

	none_type = type(None)
	assert type(lines) == none_type or isinstance(lines, list)
	assert isinstance(indent, int)
	assert type(center) == none_type or isinstance(center, P)
	assert type(axis) == none_type or isinstance(axis, P)
	assert type(rotate) == none_type or isinstance(rotate, Angle)
	assert type(translate) == none_type or isinstance(translate, P)

	ezcad = self._ezcad
	if ezcad._mode == EZCAD3.CNC_MODE:
            assert isinstance(lines, list)
	    spaces = " " *indent

	    # What we want to do is 4 transforms:
	    #
            #  1) Move the object such that its *center* is at the origin:
	    #  2) Rotate the boject around *axis* by *angle*:
	    #  3) Restore the object to where it was before step 1.
	    #  4) Move the object over by *translate*:
	    #
	    # Note that steps 1-3 are the rotate steps and step 4 is the
	    # the translate step.  They are done this way so that they are
	    # independent of one another.  Thus, a rotation does not affect
	    # the translation and vice versa.
	    #
	    # To further complicate things, we have to output these
	    # transforms in reverse order (i.e. 4, 3, 2, 1.)

	    # Perform step 4, the *translate* step:
	    if isinstance(translate, P):
		lines.append(
		  "{0}translate([{1:m}, {2:m}, {3:m}])".format(
		  spaces, translate.x, translate.y, translate.z))

	    # Steps 1-3 only make sence if exists and is *angle* is non-zero:
	    if isinstance(rotate, Angle) and rotate != Angle():
		# Again, we do the transforms in reverse order.
		# Start with step 3:
		if isinstance(center, P):
		    # Rotate around *center*:
		    assert isinstance(center, P)
		    lines.append(
		      "{0}translate([{1:m}, {2:m}, {3:m}])".format(
		      spaces, center.x, center.y, center.z))

		# Step 2: Rotate around *axis* or Z-axis:
		# at the origin:
		if isinstance(axis, P):
		    # Rotate around *axis*:
		    lines.append(
		      "{0}rotate(a={1}, v=[{2:m}, {3:m}, {4:m}])".format(
		      spaces, rotate, axis.x, axis.y, axis.z))
		else:
		    # Rotate around the Z axis:
		    lines.append(
		      "{0}rotate(a={1}, v=[0, 0, 1])".format(spaces, rotate))

		# This is the transform that moves the part to *center*
		# prior to rotation:
		if isinstance(center, P):
		    # Step 1: move to *center*; notice negative signs:
		    lines.append(
		      "{0}translate([{1:m}, {2:m}, {3:m}])".format(spaces,
		      -center.x, -center.y, -center.z))

    def scalar(self, scalar_path):
	""" Part dimensions: Return the {Angle} associated with {scalar_path}
	    starting from {self}. """
	
	return self.value_lookup(scalar_path, Part.SCALARS)

    def scalar_set(self, scalar_name, scalar_value):
	""" Part dimensions: Store a scaler into {self} named {scalar_name}
	    with a value of {scalar_value}.  In addition, {scalar_value}
	    is returned. """

	return self.value_set(scalar_name, scalar_value, Part.SCALARS)

    def screw(self, screw_path):
	""" Part dimensions: Return the {Screw} associated with {screw_path}
	    starting from {self}. """
	
	screw = self.value_lookup(screw_path, Part.SCREWS)
	return screw

    def screw_anchor_set(self, screw_name, anchor_point):
	""" Part dimensions: Set the anchor for the {Screw} named {screw_name}
	    in {self} to {anchor_point}. """

	if self.dimensions_refine_mode():
	    screw = self.screw_find(screw_name)
	    screw.anchor_set(anchor_point)

    def screw_attach(self, screw_path, surface, flags):
	""" Part dimensions: Attach the {Part} specified by {screw_path}
	    to the last {Screw} created with {Part.screw_create}.
	    The screw is attached to {surface} with hole {flags}. """

	if self.dimensions_mode():
	    screw_last = self.screw_last
	    assert screw_last != None, \
	      "Part '{0}' does not have a previously created Screw". \
	      format(self.name)

	    screw_last.attach(screw_path, surface, flags)

    def screw_create(self, screw_name, thread, direction):
	""" Part dimensions: Store a new {Screw} named {screw_name} into
	    {self}.  The new {screw} will contain {thread} and {direction}.
	    In addtion, {screw_value} is returned. """

	assert not self.part_tree_mode(), \
	  "Can not create screw in part {0} before dimensions mode". \
	  format(self.name)

	# Stuff it into {self} in the screws table:
	if self.dimensions_define_mode():
	    screw = Screw(self, screw_name, thread, direction)
	    self.value_set(screw_name, screw, Part.SCREWS)
	else:
	    screw = self.value_lookup(screw_name, Part.SCREWS)
	self.screw_last = screw
	return screw

    def screw_depth_set(self, screw_name, depth):
	""" Part dimensions: Set screw depth of the {Screw} named {screw_depth}
	    to {depth} for {self}. """

	if self.dimensions_refine_mode():
	    screw_level = self.screw_level_find(screw_name)
	    screw_level.depth_set(depth)

    def screw_diameter(self, screw, flags):
	""" Part internal: Return the diamter associated with {screw}
	    using {flags} to select the style of hole. """

	# http://www.electroimpact.com/company/QMS-0003.pdf

	# Create letter and number drill table:
	drill = {}
	drill["80"] = 0.0135
	drill["79"] = 0.0145
	drill["78"] = 0.0160
	drill["77"] = 0.0180
	drill["76"] = 0.0200
	drill["75"] = 0.0210
	drill["74"] = 0.0225
	drill["73"] = 0.0240
	drill["72"] = 0.0250
	drill["71"] = 0.0260
	drill["70"] = 0.0280
	drill["69"] = 0.0292
	drill["68"] = 0.0310
	drill["67"] = 0.0320
	drill["66"] = 0.0330
	drill["65"] = 0.0350
	drill["64"] = 0.0360
	drill["63"] = 0.0370
	drill["62"] = 0.0380
	drill["61"] = 0.0390
	drill["60"] = 0.0400
	drill["59"] = 0.0410
	drill["58"] = 0.0420
	drill["57"] = 0.0430
	drill["56"] = 0.0465
	drill["55"] = 0.0520
	drill["54"] = 0.0550
	drill["53"] = 0.0595
	drill["52"] = 0.0635
	drill["51"] = 0.0670
	drill["50"] = 0.0700
	drill["49"] = 0.0730
	drill["48"] = 0.0760
	drill["47"] = 0.0785
	drill["46"] = 0.0810
	drill["45"] = 0.0820
	drill["44"] = 0.0860
	drill["43"] = 0.0890
	drill["42"] = 0.0935
	drill["41"] = 0.0960
	drill["40"] = 0.0980
	drill["39"] = 0.0995
	drill["38"] = 0.1015
	drill["37"] = 0.1040
	drill["36"] = 0.1065
	drill["35"] = 0.1100
	drill["34"] = 0.1110
	drill["33"] = 0.1130
	drill["32"] = 0.1160
	drill["31"] = 0.1200
	drill["30"] = 0.1285
	drill["29"] = 0.1360
	drill["28"] = 0.1405
	drill["27"] = 0.1440
	drill["26"] = 0.1470
	drill["25"] = 0.1495
	drill["24"] = 0.1520
	drill["23"] = 0.1540
	drill["22"] = 0.1570
	drill["21"] = 0.1590
	drill["20"] = 0.1610
	drill["19"] = 0.1660
	drill["18"] = 0.1695
	drill["17"] = 0.1730
	drill["16"] = 0.1770
	drill["15"] = 0.1800
	drill["14"] = 0.1820
	drill["13"] = 0.1850
	drill["12"] = 0.1890
	drill["11"] = 0.1910
	drill["10"] = 0.1935
	drill["9"] = 0.1960
	drill["8"] = 0.1990
	drill["7"] = 0.2010
	drill["6"] = 0.2040
	drill["5"] = 0.2055
	drill["4"] = 0.2090
	drill["3"] = 0.2130
	drill["2"] = 0.2210
	drill["1"] = 0.2280
	drill["A"] = 0.2340
	drill["B"] = 0.2380
	drill["C"] = 0.2420
	drill["D"] = 0.2460
	drill["E"] = 0.2500
	drill["F"] = 0.2570
	drill["G"] = 0.2610
	drill["H"] = 0.2660
	drill["I"] = 0.2720
	drill["J"] = 0.2770
	drill["K"] = 0.2810
	drill["L"] = 0.2900
	drill["M"] = 0.2950
	drill["N"] = 0.3020
	drill["O"] = 0.3160
	drill["P"] = 0.3230
	drill["Q"] = 0.3320
	drill["R"] = 0.3390
	drill["S"] = 0.3480
	drill["T"] = 0.3580
	drill["U"] = 0.3680
	drill["V"] = 0.3770
	drill["W"] = 0.3860
	drill["X"] = 0.3970
	drill["Y"] = 0.4040
	drill["Z"] = 0.4130

	# Create metric screw table {m}:
	# (75% M, 75% I, 50% M, 50% I, close M, close I, loose M, loose I)
	m = {}
	m["M1.5x0.35"] = (1.15, "56", 1.25, "55", 1.60, "1/16", 1.65, "52")
	m["M1.6x0.35"] = (1.25, "55", 1.35, "54", 1.70, "51", 1.75, "50")
	m["M1.8x0.35"] = (1.45," 53", 1.55, "1/16", 1.90, "49", 2.00, "5/64")
	m["M2x0.45"] = (1.55, "1/16", 1.70, "51", 2.10, "45", 2.20, "44")
	m["M2x0.40"] = (1.60, "52", 1.75, "50", 2.10, "45", 2.20, "44")
	m["M2.2x0.45"] = (1.75, "50", 1.90, "48", 2.30, "3/32", 2.40, "41")
	m["M2.5x0.45"] = (2.05, "46", 2.20, "44", 2.65, "37", 2.75, "7/64")
	m["M3x0.60"] = (2.40, "41", 2.60, "37", 3.15, "1/8", 3.30, "30")
	m["M3x0.50"] = (2.50, "39", 2.70, "36", 3.15, "1/8", 3.30, "30")
	m["M3.5x0.60"] = (2.90, "32", 3.10, "31", 3.70, "27", 3.85, "24")
	m["M4x0.75"] = (3.25, "30", 3.50, "28", 4.20, "19", 4.40, "17")
	m["M4x0.70"] = (3.30, "30", 3.50, "28", 4.20, "19", 4.40, "17")
	m["M4.5x0.75"] = (3.75, "25", 4.00, "22", 4.75, "13", 5.00, "9")
	m["M5x1.0"] = (4.00, "21", 4.40, "11/64", 5.25, "5", 5.50, "7/32")
	m["M5x0.90"] = (4.10, "20", 4.40, "17", 5.25, "5", 5.50, "7/32")
	m["M5x0.80"] = (4.20, "19", 4.50, "16", 5.25, "5", 5.50, "7/32")
	m["M5.5x0.90"] = (4.60, "14", 4.90, "10", 5.80, "1", 6.10, "B")
	m["M6x1.0"] = (5.00, "8", 5.40, "4", 6.30, "E", 6.60, "G")
	m["M6x0.75"] = (5.25, "4", 5.50, "7/32", 6.30, "E", 6.60, "G")
	m["M7x1.0"] = (6.00, "B", 6.40, "E", 7.40, "L", 7.70, "N")
	m["M7x0.75"] = (6.25, "D", 6.50, "F", 7.40, "L", 7.70, "N")
	m["M8x1.25"] = (6.80, "H", 7.20, "J", 8.40, "Q", 8.80, "S")
	m["M8x1.0"] = (7.00, "J", 7.40, "L", 8.40, "Q", 8.80, "S")
	m["M9x1.25"] = (7.80, "N", 8.20, "P", 9.50, "3/8", 9.90, "25/64")
	m["M9x1.0"] = (8.00, "O", 8.40, "21/64", 9.50, "3/8", 9.90, "25/64")
	m["M10x1.50"] = (8.50, "R", 9.00, "T", 10.50, "Z", 11.00, "7/16")
	m["M10x1.25"] = \
	    (8.80, "11/32", 9.20, "23/64", 10.50, "Z", 11.00, "7/16")
	m["M10x1.0"] = \
	    (9.00, "T", 9.40, "U", 10.50, "Z", 11.00, "7/16")
	m["M11x1.50"] = \
	    (9.50, "3/8", 10.00, "X", 11.60, "29/64", 12.10, "15/32")
	m["M12x1.75"] = \
	    (10.30, "13/32", 10.90, "27/64", 12.60, "1/2", 13.20, "33/64")
	m["M12x1.50"] = \
	    (10.50, "Z", 11.00, "7/16", 12.60, "1/2", 13.20, "33/64")
	m["M12x1.25"] = \
	    (10.80, "27/64", 11.20, "7/16", 12.60, "1/2", 13.20, "33/64")
	m["M14x2.0"] = \
	    (12.10, "15/32", 12.70, "1/2", 14.75, "37/64", 15.50, "39/64")
	m["M14x1.50"] = \
	    (12.50, "1/2", 13.00, "33/64", 14.75, "37/64", 15.50, "39/64")
	m["M14x1.25"] = \
	    (12.80, "1/2", 13.20, "33/64", 14.75, "37/64", 15.50, "39/64")
	m["M15x1.50"] = \
	    (13.50, "17/32", 14.00, "35/64", 15.75, "5/8", 16.50, "21/32")
	m["M16x2.0"] = \
	    (14.00, "35/64", 14.75, "37/64", 16.75, "21/32", 17.50, "11/16")
	m["M16x1.50"] = \
	    (14.50, "37/64", 15.00, "19/32", 16.75, "21/32", 17.50, "11/16")
	m["M17x1.50"] = \
	    (15.50, "39/64", 16.00, "5/8", 18.00, "45/64", 18.50, "47/64")
	m["M18x2.50"] = \
	    (15.50, "39/64", 16.50, "41/64", 19.00, "3/4", 20.00, "25/32")
	m["M18x2.0"] = \
	    (16.00, "5/8", 16.75, "21/32", 19.00, "3/4", 20.00, "25/32")
	m["M18x1.50"] = \
	    (16.50, "21/32", 17.00, "43/64", 19.00, "3/4", 20.00, "25/32")
	m["M19x2.50"] = \
	    (16.50, "21/32", 17.50, "11/16", 20.00, "25/32", 21.00, "53/64")
	m["M20x2.50"] = \
	    (17.50, "11/16", 18.50, "23/32", 21.00, "53/64", 22.00, "55/64")
	m["M20x2.0"] = \
	    (18.00, "45/64", 18.50, "47/64", 21.00, "53/64", 22.00, "55/64")
	m["M20x1.50"] = \
	    (18.50, "47/64", 19.00, "3/4", 21.00, "53/64", 22.00, "55/64")

	diameter = None
	if screw.find("M") == 0:
	    # We have metric:
	    #print "metric_screw"
	    if screw in m:
		# We have found the metric screw:
		values = m[screw]


		# Use {flags} to select {drill_name}:
		drill_name = None
		if flags.find("d") >= 0:
		    # threaDed:
		    if self._material.find("steel") == 0:
			# Hard material => 50% tap:
			drill_name = values[3]
		    else:
			# Soft material => 75% tap:
			drill_name = values[1]
		elif flags.find("p") >= 0:
		    # Loose fit:
		    drill_name = values[7]
		else:
		    # Close fit: 
		    drill_name = values[5]

		# Deal with fractional, number and letter drills:
		if drill_name.find('/') >= 0:
		    # We have a fractional drill:
		    numerator_denominator = drill_name.split('/')
		    diameter = float(numerator_denominator[0]) / \
		      float(numerator_denominator[1])
		else:
		    # We have a number or letter drill:
		    diameter = drill[drill_name]
	else:
	    # Try imperial:
            # (75%, 50%, close, loose, countersink)
	    values = None
	    if screw.find("0-") == 0:
		# Number 0 screw:
		if screw == "0-80":
		    values = (.0469, .0520, .0635, .0700, None)
	    elif screw.find("2-") == 0:
		# Number 2 screw:
		if screw == "2-56":
		    values = (.0700, .0730, .0890, .0960, 0.214)
	    elif screw.find("4-") == 0:
		# Number 4 screw:
	 	if screw == "4-40":
		    values = (.0890, .0960, .1160, .1285, 0.272)
	    elif screw.find("6-") == 0:
		# Number 6 screw:
		if screw == "6-32":
		    values = (.1065, .1160, .1440, .1495, 0.324)
	    elif screw.find("8-") == 0:
		# Number 8 screw:
		if screw == "8-32":
		    values = (.1260, .1440, .1695, .1770, 0.376)
	    elif screw.find("10-") == 0:
		# Number 10 screw:
		if screw == "10-32":
		    values = (.1517, .1590, .1960, .2010, 0.428)
		elif screw == "10-24":
		    values = (.1495, .1610, .1960, .2010, 0.428)
	    elif screw.find("12-") == 0:
		# Number 12 screw:
		if screw == "12-24":
		    values = (.1649, .1770, .2210, .2280, None)
	    assert values != None, "Screw '{0}' is not recognized".format(screw)

	    if flags.find("d") >= 0:
		# threaDed:
		if self._material._find("steel"):
		    # Hard material => 50% tap:
		    diameter = values[1]
		else:
		    # Soft material => 75% tap:
		    diameter = values[0]
	    elif flags.find("s") >= 0:
		# Loose fit:
		diameter = values[3]
	    elif flags.find("!") >= 0:
		# Return counersink
		diameter = values[4]
		assert diameter != None, \
		  "Need to fill in table value for countersink diamete"
	    else:
		# Close fit: 
		diameter = values[2]

	#print "flags={0} values={1} diameter={2}". \
	#  format(flags, values, diameter)

	assert diameter != None, "Screw {0} not recognized".format(screw)
	return L(inch=diameter)

    def screw_find(self, screw_name):
	""" Part dimensions: Return the {Screw} associated with {screw_name}
	    in {self}. """

	# Iterate across {screw_levels} in {self}:
	screw_level = self.screw_level_find(screw_name)
	return screw_level.screw

    def screw_hole(self, comment, screw, start_point, end_point, flags):
	""" Part construct: Make a hole for {screw} in {self} with starting
	    at {start_point} and ending at {end_point}.  {comment} will
	    show up in any error messages and any generated G-code.
	    The allowed flag letters in {flags} are:

	      One of 't' (default), 'f', or 'p':
		't' through hole (i.e. Through)
		'f' flat hole (i.e. Flat)
		'p' tip hole (drill tip stops at {end_point} (i.e. tiP)

	      One of 'd', 'c' (default), or 's':
		'd' hole is to be tapped (i.e. threaDed)
		'c' hole is a close fit (i.e. Close fit)
		's' hole is a loose fit (i.e. LooSe fit)

	      Allowed additional flags:	  
		'u' upper hole edge should be chamfered (i.e. Upper)
		'l' lower hole edge should be chamfered (i.e. Lower)
		'm' hole is to be milled (i.e. Milled)
		'i' imperial vs. metric swapping allowed (i.e. Interchangeable)
		'a' counterink hole for a flAthead screw
	"""

	# Check argument types:
	assert isinstance(comment, str)
	assert isinstance(screw, str)
	assert isinstance(start_point, P)
	assert isinstance(end_point, P)
	assert isinstance(flags, str)

	#print "Part.screw_hole: Part={0} screw={1} flags={2}". \
	#  format(self.name, screw, flags)

	screw_diameter = self.screw_diameter(screw, flags)
	countersink_diameter = L.inch(0)
	if flags.find('a') != 0:
	    countersink_diameter = self.screw_diameter(screw, '!')

	self.hole(comment, screw_diameter, \
	  start_point, end_point, flags, countersink_diameter)
	
    def screw_holes(self, errors_suppress = False):
	""" Part construct: Perform all pending screws for the currently
	    mounted surface. """

	#print "=>Part.screw_holes({0})".format(self.name)

	# Some useful abbreviations:
	big = L(inch=123456789.0)
	
	ezcad = self._ezcad
	xml_stream = ezcad._xml_stream
	if xml_stream != None:
	    # Grab the six surfaces from {self}:
	    b = self.b
	    e = self.e
	    n = self.n
	    s = self.s
	    t = self.t
	    w = self.w

	    # Extract the bounding box of {self}:
	    t_z = t.z
	    b_z = b.z
	    n_y = n.y
	    s_y = s.y
	    e_x = e.x
	    w_x = w.x

	    # Figure out {screw_levels_table} based on {top_surface}:
	    screw_levels = self._screw_levels
	    top_surface = self._top_surface
	    screw_levels_table = {}
	    if self._top_surface_set:
		if top_surface == t:
		    screw_levels_table = screw_levels["T"]
		    top_surface_name = "Top"
		    t_z = big
		    b_z = -big
		elif top_surface == b:
		    screw_levels_table = screw_levels["B"]
		    top_surface_name = "Bottom"
		    t_z = big
		    b_z = -big
		elif top_surface == n:
		    screw_levels_table = screw_levels["N"]
		    top_surface_name = "North"
		    n_y = big
		    s_y = -big
		elif top_surface == s:
		    screw_levels_table = screw_levels["S"]
		    top_surface_name = "South"
		    n_y = big
		    s_y = -big
		elif top_surface == e:
		    screw_levels_table = screw_levels["E"]
		    top_surface_name = "East"
		    e_x = big
		    w_x = -big
		elif top_surface == w:
		    screw_levels_table = screw_levels["W"]
		    top_surface_name = "West"
		    e_x = big
		    w_x = -big
		else:
		    assert False, \
		      "Unexpected top surface for part {0} is {1}". \
		      format(self.name, top_surface)

	    if len(screw_levels_table) == 0:
		assert errors_suppress, \
		  "Part '{0}' does not have any attached holes on {1} surface".\
		  format(self.name, top_surface_name)
	    else:
		for screw_level in screw_levels_table.values():
		    #print "screw_level=", screw_level
		    screw = screw_level.screw
		    #print "screw.name=", screw.name
		    trace = False
		    #trace = screw.name.find("skin_west_bottom") >= 0

		    # Grap {anchor_point_mapped} from {screw}:
		    anchor_point_mapped = screw.anchor_point_mapped
		    if trace:
			print "=>Part.screw_holes({0})".format(self.name)
			print "anchor_point_mapped={0}". \
			  format(anchor_point_mapped)

		    remap_matrix = screw_level.reverse_matrix
		    if trace:
			print "screw_level_forward_matrix=\n{0}". \
			  format(screw_level.forward_matrix.mat	)
			print "screw_level_reverse_matrix=\n{0}". \
			  format(screw_level.reverse_matrix.mat)
			#print "remap_matrix=\n{0}".format(remap_matrix.mat)
		    anchor_point_remapped = \
		      remap_matrix.point_multiply(anchor_point_mapped, self)
		    x = anchor_point_remapped.x
		    y = anchor_point_remapped.y
		    z = anchor_point_remapped.z
		    if trace:
			bse = self.point("$BSE")
			tnw = self.point("$TNW")
			print "bse={0}".format(bse)
			print "tnw={0}".format(tnw)
			print "anchor_point_remapped={0}". \
			  format(anchor_point_remapped)
			#print \
			#  "w_x={0} e_x={1} s_y={2} n_y={3} b_z={4} t_z={5}". \
			#  format(w_x, e_x, s_y, n_y, b_z, t_z)

		    # Make sure everything is in on the part:
		    assert w_x <= x and x <= e_x, \
		      ("X (={0}) for screw {1} not between {2}" + \
		      " and {3} (part={4})"). \
		      format(x, screw.name, w_x, e_x, self.name)
		    assert s_y <= y and y <= n_y, \
		      ("Y (={0}) for screw {1} not between {2}" + \
		      " and {3} (part={4})"). \
		      format(y, screw.name, s_y, n_y, self.name)
		    assert b_z <= z and z <= t_z, \
		      ("Z (={0}) for screw {1} not between {2}" + \
		      " and {3} (part={4})"). \
		      format(z, screw.name, b_z, t_z, self.name)

		    # Compute {start_point} and {end_point} based on
		    # {top_surface} and {depth}:
		    depth = screw_level.depth
		    if top_surface == t:
			start_point = self.point_new(x, y, t.z)
			end_point = start_point.z_adjust(-depth)
		    elif top_surface == b:
			start_point = self.point_new(x, y, b.z)
			end_point = start_point.z_adjust(depth)
		    elif top_surface == n:
			start_point = self.point_new(x, n.y, z)
			end_point = start_point.y_adjust(-depth)
		    elif top_surface == s:
			start_point = self.point_new(x, s.y, z)
			end_point = start_point.y_adjust(depth)
		    elif top_surface == e:
			start_point = self.point_new(e.x, y, z)
			end_point = start_point.x_adjust(-depth)
		    elif top_surface == w:
			start_point = self.point_new(w.x, y, z)
			end_point = start_point.x_adjust(depth)
		    else:
			assert False

		    # No do either a through hole or a hole to {depth}:
		    if not screw_level.done:
			thread = screw.thread
			flags = screw_level.flags
			#print "Part.screw_holes: part={0} screw={1}" + \
			#  " thread={2} flags={3}". \
			#  format(self.name, screw_level.screw.name, \
			#  thread, flags)
			if depth <= L():
			    # Drill all the way through:
			    self.screw_through(screw.name,thread, \
			      start_point, flags)
			else:
			    # Drill to the specified {end_point}:
			    self.screw_hole(screw.name, \
			      thread, start_point, end_point, flags)

			# Remember that we did this {screw_level}:
			screw_level.done = True

		if trace:
		    print "<=Part.screw_holes({0})".format(self.name)

    def screw_level_find(self, screw_name):
	""" Part internal: Return the {Screw_Level} associated with {screw_name}
	    in {self}."""

	# Iterate across {screw_levels} in {self}:
	screw_level = None
	for screw_level_table in self.screw_levels.values():
            # Is {screw} name in {screw_level_table}:
            if screw_name in screw_level_table:
		# Got it:
		screw_level = screw_level_table[screw_name]
	assert screw_level != None, "No screw named '{0}' in part '{1}'". \
	  format(screw_name, self.name)

	return screw_level


    def screw_through(self, comment, \
      screw, start_point, flags, countersink_diameter = L(0.0)):
	""" Part construct: Make a hole for {screw} in {self} with starting
	    at {start_point} all the way through {part}.  {comment} will
	    show in any error messages and any generated G-code.
	    The allowed flag letters in {flags} are:

	      One of 'd', 'c' (default), or 's':
		'd' hole is to be tapped (i.e. threaDed)
		'c' hole is a close fit (i.e. Close fit)
		's' hole is a loose fit (i.e. LooSe fit)

	      Allowed additional flags:	  
		'u' upper hole edge should be chamfered (i.e. Upper)
		'l' lower hole edge should be chamfered (i.e. Lower)
		'm' hole is to be milled (i.e. Milled)
		'i' imperial vs. metric swapping allowed (i.e. Interchangeable)
		'a' counterink hole for a flAthead screw

	    """

	screw_diameter = self.screw_diameter(screw, flags)
	if flags.find('a') >= 0:
	    countersink_diameter = self.screw_diameter(screw, '!')

	#print "Part.screw_through: Part={0} screw={1} flags={2} cs={3}". \
	#  format(self.name, screw, flags, countersink_diameter)

	self.hole_through(comment, screw_diameter, \
	  start_point, flags, countersink_diameter)

    def show(self, indent):
	""" Part debug: Show {self} indented with {indent}. """

	print "%sPart: name='%s'" % (indent, self.name)
	for part in self.parts.values():
	    part.show(indent + " ")

    def simple_pocket(self, comment = "no comment",
      bottom_corner = None, top_corner = None, radius = None, pocket_top = "t",
      center = None, axis = None, rotate = None, translate = None):
	""" *Part*: Create a simple rectangular pocket in the *Part* object (i.e. *self*)
	    bounding corners of *bottom_corner* and *top_corner*, a corner radius if *radius*.
	"""

	print("block_corners(name='{0}', top_corner={1:i}, bottomcorner={2:i}, radius={3:i})".
	  format(self._name, top_corner, bottom_corner, radius))

	# Deal with argument defaults:
	if not isinstance(bottom_corner, P):
	    zero = L()
	    bottom_corner = P(zero, zero, zero)
	if not isinstance(top_corner, P):
	    one = L(cm=1.0)
	    top_corner = P(one, one, one)

	# Check argument types:
	assert isinstance(comment, str)
	assert isinstance(bottom_corner, P)
	assert isinstance(top_corner, P)
	assert isinstance(radius, L)
	assert isinstance(pocket_top, str)
	assert center == None or isinstance(center, P)
	assert axis == None or isinstance(axis, P)
	assert rotate == None or isinstance(rotate, Angle)
	assert translate == None or isinstance(translate, P)

	# Make sure that the corners are diagonal from bottom south west
	# to top north east:
	x1 = min(bottom_corner.x, top_corner.x)
	x2 = max(bottom_corner.x, top_corner.x)
	y1 = min(bottom_corner.y, top_corner.y)
	y2 = max(bottom_corner.y, top_corner.y)
	z1 = min(bottom_corner.z, top_corner.z)
	z2 = max(bottom_corner.z, top_corner.z)
	#print("Part.box:{0:m}:{1:m},{2:m}:{3:m},{4:m}:{5:m}". \
	#  format(x1, x2, y1, y2, z1, z2))

	# Deal with *top* argument:
	extra = L(mm = 1.0)
	ezcad = self._ezcad
	adjust = ezcad._adjust
	#print("simple_pocket:B:x1/x2 y1/y2 z1/z2 = {0}/{1} {2}/{3} {4}/{5}".
	#  format(x1, x2, y1, y2, z1, z2))
	if pocket_top == "t":
	    z2 += extra
	    x1 += adjust
	    x2 -= adjust
	    y1 += adjust
	    y2 -= adjust
	elif pocket_top == "b":
	    z1 -= extra
	    x1 += adjust
	    x2 -= adjust
	    y1 += adjust
	    y2 -= adjust
	elif pocket_top == "n":
	    y2 += extra
	    x1 += adjust
	    x2 -= adjust
	    z1 += adjust
	    z2 -= adjust
	elif pocket_top == "s":
	    y1 -= extra
	    x1 += adjust
	    x2 -= adjust
	    z1 += adjust
	    z2 -= adjust
	elif pocket_top == "e":
	    x2 += extra
	    y1 += adjust
	    y2 -= adjust
	    z1 += adjust
	    z2 -= adjust
	elif pocket_top == "w":
	    x1 -= extra
	    y1 += adjust
	    y2 -= adjust
	    z1 += adjust
	    z2 -= adjust
	else:
            assert False, \
	      "pocket_top = '{0}' instead of 't', 'b', 'n', 's', 'e', or 'w'". \
	      format(pocket_top)
	#print("simple_pocket:A:x1/x2 y1/y2 z1/z2 = {0}/{1} {2}/{3} {4}/{5}".
	#  format(x1, x2, y1, y2, z1, z2))

	dz = z1 - z2

	place = Place(part = None, name = comment, center = center,
	  axis = axis, rotate = rotate, translate = translate)
	forward_matrix = place._forward_matrix

	ezcad = self._ezcad
	if ezcad._mode == EZCAD3.CNC_MODE:
	    difference_lines = self._scad_difference_lines

	    assert x1 < x2, "x1={0} should be less than x2={1}".format(x1, x2)
	    assert y1 < y2, "y1={0} should be less than y2={1}".format(y1, y2)
	    assert z1 < z2, "z1={0} should be less than z2={1}".format(z1, z2)

	    #print "c1=({0},{1},{2}) c2=({3},{4},{5})".format( \
	    #  x1, y1, z1, x2, y2, z2)

            # The transforms are done in reverse order:

	    self._scad_transform(difference_lines, center = center,
	      axis = axis, rotate = rotate, translate = translate)

	    # Get the lower south west corner positioned:
	    difference_lines.append(
	      "      translate([{0:m}, {1:m}, {2:m}])".format(x1, y1, z1))

	    # Finally, we can output the cube:
	    difference_lines.append(
	      "        cube([{0:m}, {1:m}, {2:m}]);".format(
	      x2 - x1, y2 - y1, z2 - z1))

	    maximum_diameter = 2 * radius
	    end_mill_tool = self._tools_end_mill_search(maximum_diameter, dz, "simple_pocket")
	    print("end_mill_toll=", end_mill_tool)
	    print("end_mill={0}".format(end_mill_tool._name_get()))
	    end_mill_diameter = end_mill_tool._diameter_get()
	    end_mill_radius = end_mill_diameter / 2
	    end_mill_feed_speed = end_mill_tool._feed_speed_get()
	    end_mill_spindle_speed = end_mill_tool._spindle_speed_get()

	    print("x1={0:i} x2={1:i} y1={2:i} y2={3:i}".format(x1, x2, y1, y2))
	    print("z1={0:i} z2={1:i}".format(z1, z2))
	    operation_order = Operation.ORDER_END_MILL_SIMPLE_POCKET
	    operation_simple_pocket = Operation_Simple_Pocket(self, comment, 0,
	      end_mill_tool, operation_order, None, end_mill_feed_speed, end_mill_spindle_speed,
	      x1, y1, x2, y2, z1, z2, radius, end_mill_radius, Operation.POCKET_KIND_FLAT)
	    self._operation_append(operation_simple_pocket)

    def xxx_simple_pocket(self, comment, \
      corner1_point, corner2_point, radius, flags):
	""" Part construct: Mill a pocket in {self} where {corner1_point}
	    and {corner2_point} specify a diagonal across the pocket.  The
	    radius of the inside corners is {radius}.  {flags} can have
	    the character 't' for a pocket that goes {self} and 'u' to
	    specify that an upper chamfer is requested.  {comment} will
	    show up in error messages and any generated G-code. """

	# Perform argument type checking:
	assert isinstance(comment, str)
	assert isinstance(corner1_point, P)
	assert isinstance(corner2_point, P)
	assert isinstance(radius, L)
	assert isinstance(flags, str)

	# Verify that each flag in {flags} is OK:
	for flag in flags:
	    assert flag == 'u' or flag == 't', \
	      'simple_pocket("{0}"): Bad flag "{1}" in "{2}" for part {3}'. \
	      format(comment, flag, flags, self.name)

	# Extract some values from {ezcad}:
	ezcad = self._ezcad
	xml_indent = ezcad._xml_indent
	xml_stream = ezcad._xml_stream

	if xml_stream != None:
	    # <Simple_Pocket C1X= C1Y= C1Z= C2X= C2Y= C2Z= ...
	    #  ... Radius= Flags= Comment= />:
	    xml_stream.write('{0}<Simple_Pocket'.format(" " * xml_indent))
	    xml_stream.write(' C1X="{0}" C1Y="{1}" C1Z="{2}"'. \
	      format(corner1_point.x, corner1_point.y, corner1_point.z))
	    xml_stream.write(' C2X="{0}" C2Y="{1}" C2Z="{2}"'. \
	      format(corner2_point.x, corner2_point.y, corner2_point.z))
	    xml_stream.write(' Radius="{0}" Flags="{1}" Comment="{2}"/>\n'. \
	      format(radius, flags, comment))

    def manufacture(self):
	""" Part: Override this method to do any manufacturing steps. """

	print("No manufacture method for '{0}'".format(self))

    def string(self, string_path):
	""" Part dimensions: Return the {String} associated with {string_path}
	    starting from {self}. """
	
	return self.value_lookup(scalar_path, Part.STRINGS)

    def string_set(self, string_name, string_value):
	""" P dimensions: Store Set the value of {string_name} in {self}
	    to {string_value}.  In addition, {string_value} is returned. """

	return self.value_set(self, string_name, string_value, Part.STRINGS)

    def tool_prefer(self, tool_name):
	""" Part construct: Cause {tool_name} to be the preferred tool
	    for {self}.  If {tool_name} is an empty string, the preferred
	    tool is cleared. """

	ezcad = self._ezcad
	xml_stream = ezcad._xml_stream
	if xml_stream != None:
	    xml_stream.write( \
	      '{0}<Tool_Prefer Tool_Name="{1}"/>\n'. \
	      format(" " * ezcad._xml_indent, tool_name))

    def tooling_plate(self, comment, flags, trace = False):
	""" Part construct: Cause a grid of tooling plate holes to be
	drilled into {self}.  {flags} are used to adjust the position
	the tooling plate holes.  The grid is 2 row by 2 column
	(i.e. 2 by 2) but {flags} can be modified to specifically
	set the number of rows and/or columns.  For example,
	"3r 4c" establishes a grid that 3 rows by 4 columns.
	The grid holes are initially distributed to just fill
	the {sefl} bounding box in the X and Y dimensions.
	Additional flags can be specified to move the holes
	around a little.  Tooling plate holes are specified by
	row and column 	where a number specifies the row and an
	upper case letter specifies the column.  The row numbers
	start with 1 on the "north" side and move southwards.
	The column letters start with the letter 'A' on the "west"
	and move eastwards.  Lower case letters 'e' (for East),
	'n' (for North), 's' (for South), and 'w' (for west) are used
	to move a tooling plate hole east, north, south, or west
	respectively.  The lower case letter 'x' (for eXclude)
	can be used to remove a tooling hole  A tooling hole
	can be moved by following it by one or more movement letters.
	The following example "5r 5c 2Bnw 2Dne 3Cx 4Bsw 4Dse" creates
	a 5x5 grid of tooling holes, the middle hole is removed "3Cr"
	and 4 of the tooling holes on the inside are moved outward
	from the center. """

	#if self.construct_mode():
	if True:
	    # Parse {flags}:

	    # {adjusts} lists any holes that have been moved:
	    adjusts = {}

	    # {removes} lists any holes that been removed or ignored:
	    removes = {}

	    # Initialize some variables for the parsing process:
	    column = 0
	    columns = 2
	    number = 0
	    row = 0
	    row_column = (0, 0)
	    rows = 2
	    adjusts[row_column] = (0, 0)

	    # Iterate across {flags}:
	    for flag in flags:
		if '0' <= flag and flag <= '9':
		    # We have a decimal digit, add it to {number}:
		    number = number * 10 + ord(flag) - ord('0')
		else:
		    # We have something other than a decimal digit:
		    if 'A' <= flag and flag <= 'Z':
			# We have a column letter
			column = ord(flag) - ord('A')
			row = number - 1
			row_column = (row, column)
			if not (row_column in adjusts):
			    adjusts[row_column] = (0, 0)
		    elif flag == 'c':
			columns = number
			if trace:
			    print "columns=", columns
		    elif flag == 'e':
			adjust = adjusts[row_column]
			adjusts[row_column] = (adjust[0] + 1, adjust[1])
		    elif flag == 'i':
			removes[row_column] = 'i'
		    elif flag == 'n':
			adjust = adjusts[row_column]
			adjusts[row_column] = (adjust[0], adjust[1] + 1)
		    elif flag == 'r':
			rows = number
			if trace:
			    print "rows=", rows
		    elif flag == 's':
			adjust = adjusts[row_column]
			adjusts[row_column] = (adjust[0], adjust[1] - 1)
		    elif flag == 'w':
			adjust = adjusts[row_column]
			adjusts[row_column] = (adjust[0] - 1, adjust[1])
		    elif flag == 'x':
			removes[row_column] = 'x'
		    elif flag == ' ':
		        number = 0
		    else:
			assert False, \
			  "Bad flag letter '{0}' in flags '{1}'". \
			  format(flag, flags)

		    # All non-digit operations reset {number} at the end:
		    number = 0

	    if trace:
		print "removes=", removes
		print "adjusts=", adjusts

	    ezcad = self._ezcad
	    xml_stream = ezcad._xml_stream
	    if xml_stream != None:
		xml_stream.write('{0}<Tooling_Plate Rows="{1}" '. \
		  format(" " * ezcad._xml_indent, rows))
		xml_stream.write('Columns="{0}" Comment="{1}">\n'. \
		  format(columns, comment))
	    
		ezcad.xml_indent_push()
		for row in range(0, rows):
		    for column in range(0, columns):
			row_column = (row, column)
			adjust = (0, 0)
			if row_column in adjusts:
			    adjust = adjusts[row_column]
			remove = ""
			if row_column in removes:
			    remove = removes[row_column]

			if trace:
			    print "[{0},{1}] = {2}{3}". \
			      format(row, column, adjust, remove)

			xml_stream.write('{0}<Tooling_Hole'. \
			  format(" " * ezcad._xml_indent))
			xml_stream.write(' Row="{0}" Column="{1}"'. \
			  format(row, column))
			xml_stream.write(' Adjust_X="{0}" Adjust_Y="{1}"'. \
			  format(adjust[0], adjust[1]))
			xml_stream.write(' Flags="{0}"/>\n'. format(remove))
		    
		ezcad.xml_indent_pop()
		xml_stream.write('{0}</Tooling_Plate>\n'. \
		  format(" " *ezcad._xml_indent))

    def tooling_plate_mount(self, comment):
	""" Part construction: Cause the mounting plate that holds {self}
	    to be mounted. """

	#if self.construct_mode():
	ezcad = self._ezcad
	xml_stream = ezcad._xml_stream
	if xml_stream != None:
	    ezcad._xml_stream.write( \
	      '{0}<Tooling_Plate_Mount Comment="{1}"/>\n'. \
	      format(" " * ezcad._xml_indent, comment))

    def tube(self, color, material, \
      start_point, end_point, diameter, wall_thickness, sides):
	""" Part dimensions: Create a tube for {self} out of {material}
	    goes from {start_point} to {end_point} is {diameter} round
	    and has a wall thickness of {wall_thickness}.  The tube is
	    visualized with {color}.  {sides} specifies how many sides
	    to use for visualization; setting {sides} to zero gives a
	    reasonable visualization of a tube.  The tube must be
	    aligned with one of the X, Y or Z axes. """

	# Make sure we
	#assert not self.part_tree_mode(), \
	#  "Part.tube: Part '{0}' is in part tree mode".format(self.name)

	# Compute {tube_axis}, the axis along with the tube is aligned:
	tube_axis = end_point - start_point

	# Extract coordinates from {start_point} and {end_point}:
	x1 = start_point.x
	y1 = start_point.y
	z1 = start_point.z
	x2 = end_point.x
	y2 = end_point.y
	z2 = end_point.z

	zero = L.inch(0)
	radius = diameter.half()
	if tube_axis.x != zero:
	    # Tube aligned on X:
	    t = z1 + radius
	    b = z1 - radius
	    n = y1 + radius
	    s = y1 - radius
	    if x1 < x2:
		w = x1
		e = x2
	    else:
		w = x2
		e = x1
	elif tube_axis.y != zero:
	    # Tube aligned on Y:
	    t = z1 + radius
	    b = z1 - radius
	    if y1 < y2:
		s = y1
		n = y2
	    else:
		s = y2
		n = y1
	    e = x1 + radius
	    w = x1 - radius
	elif tube_axis.z != zero:
	    # Tube aligned on Z:
	    if z1 < z2:
		b = z1
		t = z2
	    else:
		b = z2
		t = z1
	    n = y1 + radius
	    s = y1 - radius
	    e = x1 + radius
	    w = x1 - radius
	else:
	    assert not self.construct_mode(), \
 	      "Tube is not aligned with X, Y, or Z axis"
	    t = z2
            b = z1
	    n = y2
	    s = y1
	    w = x1
	    e = x2

	# Compute the tube corners:
	bsw = self.point(w, s, b)
	tne = self.point(e, n, t)
	#print "bsw={0} tne={1}".format(bsw, tne)

	# Record the material in {self}:
	self._material = material

	#if self.dimensions_mode():
        #    self.bounding_box_update(w, s, b, e, n, t)

	self._bounding_box_update(w, s, b, e, n, t)

	# In consturct mode, output the <Tube ...> line:
	#if self.construct_mode():
	self._xml_lines.append( ('<Tube SX="{0}" SY="{1}" SZ="{2}"' + \
	  ' EX="{3}" EY="{4}" EZ="{5}" Outer_Diameter="{6}"' + \
	  ' Wall_Thickness="{7}" Sides="{8}" Color="{9}"' + \
	  ' Transparency="{10}" Material="{11}" Comment="{12}"/>') . \
	  format(x1, y1, z1, x2, y2, z2, diameter, wall_thickness, sides, 
	  color, self._transparency, material, self._name))

    def value_lookup(self, path, flavor):
	""" Part internal: Return the value for for {path} starting
	    from {self}.  {flavor} specifies the type to lookup. """

	if self.dimensions_refine_mode() or self.construct_mode():
	    # We are in a mode where {path} is actually looked up:

	    # Split the {part_path} and {point_name} from {point_path}:
	    part_path_value_name = path.split('.')
	    if len(part_path_value_name) == 1:
		part = self
		value_name = path
	    elif len(part_path_value_name) == 2:
		part_path = part_path_value_name[0]
		value_name = part_path_value_name[1]

		# Now iterate across {part_path} starting from {self}:
		part = self
		for part_name in part_path.split('/'):
		    # {part_name} is the next component of {part_path}:
		    if part_name == "..":
			# Move up the part tree towards the root:
			part = part.parent

			# Check for going too far up tree:
			if part.name == "":
		            assert False, \
			      "Path '%s' goes too far up the Part tree" % (path)
		    else:
			# Move down the tree:
			parts = part.parts
			if part_name in parts:
			    part = parts[part_name]
			else:
			    assert False, \
			      "Part path '%s' does not contain part '%s'" % \
			      (path, part.name)
	    else:
		# We have either an empty path, or one with two or more '.':
		assert False, "Path '%s' is not valid" % (point_path)

	    # Use {flavor} to select the {values} dictionary to use:
	    if flavor == Part.SCALARS:
		values = self.scalars
	    elif flavor == Part.ANGLES:
		values = self.angles
	    elif flavor == Part.LENGTHS:
		values = self.lengths
	    elif flavor == Part.POINTS:
		values = self.points
	    elif flavor == Part.SCREWS:
		values = self.screws
	    else:
		assert False
	    
	    # Now lookup {value_name} in {values}:
	    if value_name in values:
		value = values[value_name]
	    else:
		assert False, \
		  "Can not find '{0}' in path '{1}' starting from part {2}". \
		  format(value_name, path, self.name)
	else:
	    # We are in a mode where {path} is ignored, and an empty
	    # value is returned:
	    if flavor == Part.SCALARS:
		value = 0.0
	    elif flavor == Part.ANGLES:
		value = Angle(0.0)
	    elif flavor == Part.LENGTHS:
		value = L(0.0)
	    elif flavor == Part.POINTS:
		zero = L(0.0)
		value = P(self, zero, zero, zero)
	    elif flavor == Part.SCREWS:
		# Create a bogus {Screw}, that will not actually be used:
		value = Screw(self, "", "", "TB")
	    else:
		assert False

	return value

    def value_set(self, value_name, value, flavor):
	""" Part (internal): Store {value} into {self} under the
	    name {value_name} using the {flavor} dictionary.  In all cases,
	    {value} is returned. """

	# Select the {values} dictionary based on {flavor}:
	if flavor == Part.ANGLES:
	    values = self.angles
	    flavor_name = "Angle"
	elif flavor == Part.LENGTHS:
	    values = self.lengths
	    flavor_name = "L"
	elif flavor == Part.POINTS:
	    values = self.points
	    flavor_name = "P"
	elif flavor == Part.SCALARS:
	    values = self.scalars
	    flavor_name = "Scalar"
	elif flavor == Part.SCREWS:
	    values = self.screws
	    flavor_name = "Screw"
	else:
	    assert False

	# Routine behavior depends upon mode:
	if self.dimensions_define_mode():
	    # In define mode, we make sure that no duplicates occur:
	    assert not (value_name in values), \
	      "{0} name '{1}' is a duplicate in Part '{2}'". \
	      format(flavor_name, value_name, self.name)
	    values[value_name] = value
	elif self.dimensions_refine_mode():
	    # In refine mode, we make sure that {value_name} is defined:
	    assert value_name in values, \
	      "{0} name '{1}' is undefined in part '{2}'". \
	      format(flavor_name, angle_name, self.name)

	    # Now check to see if the value is actually changed:
	    previous_value = values[value_name]
	    if value != previous_value:
		values[value_name] = value
		self.dimensions_changed("value_set['%s']" % (value_name))
	elif self.part_tree_mode():
	    # Make sure we are not in part tree mode:
	    assert False, "Setting value in part tree mode"
	# else we are in construct mode and do nothing:

	return value

    def vertical_lathe(self, comment, axis_start_point, axis_end_point, \
      inner_diameter, outer_diameter, flags):
	""" Part construct: Mill out a tube of material from {self} along
	    the axis from {axis_start_point} to {axis_end_point} where
	    {inner_diameter} and {outer_diameter} specify the tube
	    inside and outside diameter.  {flags} can be 'i' to specify
	    that only the inside diameter matters, allowing for an end-mill
	    that extends past outside_diameter.  {comment} is used in
	    error messages and any generated G-code. """

	# Verify that each flag in {flags} is OK:
	for flag in flags:
	    assert flag == 'i', \
	      "verical_lathe('{0}'): Bad flag '{1}' in '{2}' for part {3}". \
	      format(comment, flag, flags, self.name)

	# Extract some values from {ezcad}:
	ezcad = self._ezcad
	xml_indent = ezcad._xml_indent
	xml_stream = ezcad._xml_stream

	if xml_stream != None:
	    # <Contour SX= SY= SZ= EX= EY= EZ= Engagement= Flags= Comment= />:
	    xml_stream.write('{0}<Vertical_Lathe'.format(" " * xml_indent))
	    xml_stream.write(' SX="{0}" SY="{1}" SZ="{2}"'. \
	      format(axis_start_point.x, axis_start_point.y,
	      axis_start_point.z))
	    xml_stream.write(' EX="{0}" EY="{1}" EZ="{2}"'. \
	      format(axis_end_point.x, axis_end_point.y, axis_end_point.z))
	    xml_stream.write(' Inner_Diameter="{0}" Outer_Diameter="{1}"'. \
	      format(inner_diameter, outer_diameter))
	    xml_stream.write(' Flags="{0}" Comment="{1}"/>\n'. \
	      format(flags, comment))

    def vice_position(self, comment, surface_point, north_point, west_point):
	""" Part construct: Cause {self} to be mounted in a vice with
	    {surface_point} as the top surface, {north_point} as the edge
	    mounted against the top vice edge, and {west_point} as the edge
	    left pointing edge of {self}.  {comment} is the attached
	    to any generated G-code. """

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(surface_point, P)
	assert isinstance(north_point, P)
	assert isinstance(west_point, P)

	# Use *part* instead of self:
	part = self

	# Before we get too far, flush any pending screw holes
	# from the previous mounting:
	if part._top_surface_set:
            part.screw_holes(True)

	# Extract some values from {ezcad}:
	ezcad = part._ezcad
	xml_indent = ezcad._xml_indent
	xml_stream = ezcad._xml_stream

	corner1 = part.tne
	corner2 = part.bsw
	cx = (corner1.x + corner2.x) / 2
	cy = (corner1.y + corner2.y) / 2
	cz = (corner1.z + corner2.z) / 2

	# <Vice_Position TX= TY= TZ= NX= NY= NZ= WX= WY= WZ= CX= CY= CZ=/>:
	#if xml_stream != None:
	#    xml_stream.write( \
	#      '{0}<Vice_Position TX="{1}" TY="{2}" TZ="{3}"'. \
	#      format(" " * xml_indent, surface_point.x, surface_point.y,
	#      surface_point.z))
	#    xml_stream.write( \
	#      ' NX="{0}" NY="{1}" NZ="{2}"'.format( \
	#      north_point.x, north_point.y, north_point.z))
	#    xml_stream.write(' WX="{0}" WY="{1}" WZ="{2}"'.format( \
	#      west_point.x, west_point.y, west_point.z))
	#    xml_stream.write(' CX="{0}" CY="{1}" CZ="{2}"'.format(cx, cy, cz))
	#    xml_stream.write(' Comment="{0}"/>\n'.format(comment))

	part._top_surface = surface_point
	part._top_surface_set = True

	if True:
	    # <Vice_Position TX= TY= TZ= NX= NY= NZ= WX= WY= WZ= CX= CY= CZ= />:
	    #debug = True
	    debug = False
	    trace = False

	    # To debug an individual part:
	    #if part.name = "Tine_Tang":
	    #    debug = True

	    if trace or debug:
		print("<Vice_Position Comment={0} ... />\n".format(comment))

	    # We need {zero} and {one} to define the Z axis below:
	    one = L(mm=1.0)

	    # Load up *Points*'s *t*, *n*, *w*, and *c*:
	    t = surface_point
	    n = north_point
	    w = west_point
	    c = P(cx, cy, cz)

	    # Show what we have when debugging:
	    if debug:
		print("vice_pos0:\n\tt={0}\n\tn={1}\n\tw={2}\n\tc={3}\n".
		format(t, n, w, c))

	    # Force the {position} matrix to identity.  This also forces
	    # the {reposition} matrix to identity as well:
	    part._position_reset()
	    position = part._position
	    reposition = part._reposition
	    shop = part._shop_get()
	    matrix = shop._matrix_get()

	    # Normalize everything to be centered around {c}:
	    print("c={0}".format(c))
	    print("-c={0}".format(-c))
	    part.translate(-c)
	    # Note: This updates the *position* and *reposition* matrices:

	    # Update *t*, *n*, *w*, and *c*.
	    t = position.point_multiply(t)
	    n = position.point_multiply(n)
	    w = position.point_multiply(w)
	    c = position.point_multiply(c)

	    # Show what we have when debuging:
	    if debug:
		print("vice_pos1:\n\tt={0}\n\tn={1}\n\tw={2}\n\tc={3}\n".
		  format(t, n, w,c))

	    # Rotate the part so that the surface is pointing up.  Start
	    # by computing {top_angle}, the angle between ({tx},{ty},{tz})
	    # and the positive Z axis normal:
	    zero = L()
	    z_axis = P(zero, zero, one)
	    top_angle = t.angle_between(z_axis)

	    if debug:
		print("top_angle={0}\n".format(top_angle))

	    # See if we need to rotate {part} by {top_angle}:
	    if top_angle.absolute() > Angle(deg=0.0001):
		# We need to rotate {part} by {top_angle}:
		top_axis = zero
		if top_angle.absolute() > Angle(deg=179.9):
		    # We have to entirely flip the board.  Unfortunately,
		    # a cross product will not work on colinear segments,
		    # so we use the ({n.x}, {n.y}, {n.z}) x (0, 0, 1) to
		    # establish the axis of rotation:
		    top_axis = n.cross_product(z_axis)

		    # Force *top_angle* to 180 degrees, to get rid of any
		    # rounding errors:
		    top_angle = Angle(deg=180.0)
		else:
		    # Have to partially flip the board. Use the cross product
		    # ({t.x}, {t.y}, {t.z}) x (0, 0, 1) to compute {top_axis},
		    # the axis about which to rotate the board:
		    top_axis = t.cross_product(z_axis)

		# Normalize *top_axis*:
		top_axis = top_axis.normalize()

		# Show what we have during debugging:
		if debug:
		    print("top_axis={0}".format(top_axis))

		# Now rotate {part} by {top_angle_degrees} around the
		# normal ({top_nx}, {top_ny}, {top_nz}):
		part.rotate(top_axis.x,
		  top_axis.y, top_axis.z, top_angle, zero)

		# Now T=({tx},{ty},{tz}) surface is on top (i.e. parallel to
		# to Z axis and pointing up:

		# Now that we have done the first transformation, we want to
		# work with T, N, and W in their new locations:
		t = position.point_multiply(t)
		n = position.point_multiply(n)
		w = position.point_multiply(w)
		c = position.point_multiply(c)

		# Show what we have when debugging:
		if debug:
		    print("vice_pos2:\n\tt={0}\n\tn={1}\n\tw={2}\n\tc={3}\n".
		      fromat(t, n, w, c))

	    # Now shift {part} down so that {t} is on Z=0 plane:
	    part.translate(P(zero, zero, -t.z))

	    # Like before, move T, N, W, and C to their new homes:
	    t = position.point_multiply(t)
	    n = position.point_multiply(n)
	    w = position.point_multiply(w)
	    c = position.point_multiply(c)

	    # Show what we have during debugging:
	    if debug:
		print("vice_pos3:\n\tt={0}\n\tn={1}\n\tw={2}\n\tc={3}\n".
		  format(t, n, w, c))

	    # Now rotate {part} around positive Z axis so that N =
	    # ({nx},{ny},{nz}) is facing north:

	    # Start by computing NT = N - T.
	    nt = n - t
	    vice_angle = -z_axis.angle_between(nt)

	    # Compute {vice_angle} between positive Y-axis (i.e. north) and NT:
	    if debug:
		print("NT={0} vice_angle={1}\n").format(nt, vice_angle)

	    # Are we so close to 0 degrees, that we should not bother:
	    if vice_angle.absolute() > Angle(deg=0.1):
		# No, perform the rotation:
		part.rotate(zero, zero, one, vice_angle, nt.length())

		# As before, update T, N and W to their new locations:
		t = position.point_multiply(t)
		n = position.point_multiply(n)
		w = position.point_multiply(w)

		# Show what we have during debugging:
		if debug:
		    print("vice_pos4:\n\tt={0}\n\tn={1}\n\tw={2}\n\tc={3}\n".
		      format(t, n, w, c))

	    # Start by computing NT = N - T and WT = W - T:
	    wt = w - t

	    # Compute {dowel_angle} NTW:
	    dowel_angle = nt.angle_between(wt)

	    # Verify that that the angle <NTW is 90 degrees:
	    if dowel_angle < Angle(deg=89.9) or dowel_angle > Angle(deg=90.1):
		# It is not, let somebody know:
		print("nt={0} wt={1} part='{2}' dowel_angle={3}\n".
		  format(nt, wt, part._name, dowel_angle))

	    part.reposition(n.y)

	    #extra1 :@= copy@(part.extra1)
	    #extra2 :@= copy@(part.extra2)
	    #call matrix_apply@(extra1, position)
	    #call matrix_apply@(extra2, position)

	    if debug:
		print("vice_pos:dowel_pos_set({0}, {1}, 0.0)\n".
		  format(part.name, w.x))
	    part.dowel_position_set(w.x, zero)
	    #extra_x :@= minimum@(extra1.x, extra2.x)
	    #call d@(form@("dowel_set(*, %i%, 0.0)\n") / f@(extra_x))
	    #call dowel_position_set@(part, extra_x, zero)
	    part.dowel_pin(comment)

class Fastener(Part):
    """ *Fastener*: """
        
    # Imperial/Metric hole equivalents:
    #
    #    http://www.schsm.org/html/metric_imperial_screw_equivale.html
    #
    # Metric close/free fit:
    #
    #    http://www.csgnetwork.com/screwmetmachtable.html
    # 
    #				close (in)	free (mm=in)
    # M1.6x0.35 ~= #0-80	  0.0635	1.7 = 0.0669*
    # M1.8x0.35 ~= #1-72	  0.0760*	1.9 = 0.0748
    # M2.2x0.45 ~= #2-56	  0.0890*	2.4 = 0.0866
    # M3.0x0.50 ~= #4-40	  0.1160	3.2 = 0.1260*
    # M3.5x0.60 ~= #6-40	  0.1440	3.7 = 0.1457*
    # M4.0x0.70 ~= #8-32	  0.1695*	4.3 = 0.1693
    # M5.0x0.80 ~= #10-32	  0.1960	5.3 = 0.2087*
    # M6.0x1.00 ~= 1/4-20	  0.2570*	6.4 = 0.2520

    def __init__(self, up, comment = None):
	""" *Fastener*: """
	assert isinstance(comment, str)

	Part.__init__(self, up)
	zero = L()
	self.comment_s = "NO_COMMENT"
	self.color = Color("black")
	self.material = Material("steel", "stainless")
	self.start_p = P()
	self.end_p = P()
	self.flags_s = ""
	self.major_diameter_l = zero
	self.pitch_l = zero
	self.thread75_l = zero
	self.thread50_l = zero
	self.close_fit_l = zero
	self.free_fit_l = zero	
	self.flat_head_diameter_l = zero
	self.hex_insert_b = False
	self.nut_height_l = zero
	self.hex_nut_edge_width_l = zero
	self.hex_nut_tip_width_l = zero
	self.sides_angle_a = Angle()
	self.flat_head_point_angle_a = Angle()

	if isinstance(comment, str):
	    self.comment_s = comment

    def configure(self, comment = None, material = None, color = None,
      flags = None, start = None, end = None, sides_angle = None,
      head_washer_diameter = None, tail_washer_diameter = None):
     	""" *Fastener*: """

	# Check argument types:
	none_type = type(None)
	assert type(comment) == none_type or isinstance(comment, str)
	assert type(material) == none_type or isinstance(material, Material)
	assert type(color) == none_type or isinstance(color, Color)
	assert type(flags) == none_type or isinstance(flags, str)
	assert type(start) == none_type or isinstance(start, P)
	assert type(end) == none_type or isinstance(end, P)
	assert type(sides_angle) == none_type or isinstance(sides_angle, Angle)
	
	major_diameter = None

	if isinstance(comment, str):
	    self.comment_s = comment
	if isinstance(material, Material):
	    self.material = material
	if isinstance(color, Color):
	    self.color = color
	if isinstance(sides_angle, Angle):
	    self.sides_angle_a = sides_angle
	if isinstance(flags, str):
	    self.flags_s = flags
	    for flag in flags.split(':'):
		if flag == "#0-80":
		    major_diameter = L(inch = .0600)
		    pitch = L(inch = "1/80")
		    thread75 = L(inch = "3/64")
		    thread50 = L(inch = .0520)
		    close_fit = L(inch = 0.0635)
		    free_fit = L(inch = 0.0700)
		    # Nut dims from page 1568 of Mach. Handbook (26th ed.):
		    nut_height = L(inch = 0.050)
		    hex_nut_edge_width = L(inch = "5/32")
		    hex_nut_tip_width = L(inch = 0.180)
		    flat_head_diameter = L(inch = 0.119)
		elif flag == "#1-72":
		    major_diameter = L(inch = .0730)
		    pitch = L(inch = "1/72")
		    thread75 = L(inch = .0595)
		    thread50 = L(inch = .0635)
		    close_fit = L(inch = 0.0760)
		    free_fit = L(inch = 0.0810)
		    nut_height = L(inch = 0.050)
		    hex_nut_edge_width = L(inch = "5/32")
		    hex_nut_tip_width = L(inch = 0.180)
		    flat_head_diameter = L(inch = 0.146)
		elif flag == "#2-56":
		    major_diameter = L(inch = .0860)
		    pitch = L(inch = "1/56")
		    thread75 = L(inch = .0700)
		    thread50 = L(inch = .0730)
		    close_fit = L(inch = 0.0890)
		    free_fit = L(inch = 0.0960)
		    nut_height = L(inch = 0.066)
		    hex_nut_edge_width = L(inch = "3/16")
		    hex_nut_tip_width = L(inch = 0.217)
		    flat_head_diameter = L(inch = 0.172)
		elif flag == "#4-40":
		    major_diameter = L(inch = .1120)
		    pitch = L(inch = "1/40")
		    thread75 = L(inch = .0890)
		    thread50 = L(inch = .0960)
		    close_fit = L(inch = 0.1160)
		    free_fit = L(inch = 0.1285)
		    nut_height = L(inch = 0.098)
		    hex_nut_edge_width = L(inch = "1/4")
		    hex_nut_tip_width = L(inch = 0.289)
		    flat_head_diameter = L(inch = 0.225)
		elif flag == "#6-32":
		    major_diameter = L(inch = .1380)
		    pitch = L(inch = "1/32")
		    thread75 = L(inch = .1065)
		    thread50 = L(inch = .1160)
		    close_fit = L(inch = 0.1440)
		    free_fit = L(inch = 0.1495)
		    nut_height = L(inch = 0.114)
		    hex_nut_edge_width = L(inch = "5/16")
		    hex_nut_tip_width = L(inch = 0.361)
		    flat_head_diameter = L(inch = 0.279)
		elif flag == "#10-24":
		    major_diameter = L(inch = 0.1900)
		    pitch = L(inch = "1/24")
		    thread75 = L(inch = 0.1495)
		    thread50 = L(inch = 0.1610)
		    close_fit = L(inch = 0.1960)
		    free_fit = L(inch = 0.2010)
		    nut_height = L(inch = 0.130)
		    hex_nut_edge_width = L(inch = "3/8")
		    hex_nut_tip_width = L(inch = 0.433)
		    flat_head_diameter = L(inch = 0.385)
		elif flag == "#10-32":
		    major_diameter = L(inch = 0.1900)
		    pitch = L(inch = "1/32")
		    thread75 = L(inch = 0.1590)
		    thread50 = L(inch = 0.1695)
		    close_fit = L(inch = 0.1960)
		    free_fit = L(inch = 0.2010)
		    nut_height = L(inch = 0.130)
		    hex_nut_edge_width = L(inch = "3/8")
		    hex_nut_tip_width = L(inch = 0.433)
		    flat_head_diameter = L(inch = 0.385)
		elif flag == "M3x.05":
		    major_diameter = L(mm = 3.00)
		    pitch = L(mm = 0.05)
		    thread75 = L(mm = 2.50)
		    thread50 = L(mm = 2.50)
		    close_fit = L(mm = 3.90)
		    free_fit = L(mm = 3.90)
		    nut_height = L(inch = 0.130)
		    hex_nut_edge_width = L(mm = 5.50)
		    hex_nut_tip_width = L(mm = 6.35)
		    flat_head_diameter = L(mm = 6.30)
		elif flag == "fh":
		    #print("Fastener.drill(): flat_head")
		    self.flat_head_point_angle_a = Angle(deg = 81)
		elif flag == "hi":
		    # Hex insert
		    self.hex_insert_b = True
		else:
		    assert False, \
		      "Unrecognized flag ('{0}') in flags('{1}')". \
		      format(flag, flags)

	if isinstance(major_diameter, L):
	    self.major_diameter_l = major_diameter
	    self.pitch_l = pitch
	    self.thread75_l = thread75
	    self.thread50_l = thread50
	    self.close_fit_l = close_fit
	    self.free_fit_l = free_fit
	    self.nut_height_l = nut_height
	    self.hex_nut_edge_width_l = hex_nut_edge_width
	    self.hex_nut_tip_width_l = hex_nut_tip_width
	    self.flat_head_diameter_l = flat_head_diameter

	if isinstance(start, P):
	    self.start_p = start
	if isinstance(end, P):
	    self.end_p = end

    def construct(self):
	""" *Fastener*: """

	assert self.comment_s != "NO_COMMENT", \
	  "Fastener name is not set (not configured!)"
	self.cylinder(comment = self.comment_s,
	  material = self.material,
	  color = self.color,
	  diameter = self.major_diameter_l,
	  start = self.start_p,
	  end = self.end_p)

    def nut_ledge(self, part = None, flags = ""):
        """ *Fastener*: Cut out ledge for a screw and a nut. """

	# Check argument types:
	none_type = type(None)
	assert isinstance(part, Part)
	assert isinstance(flags, str)

	if part._ezcad._mode == EZCAD3.CNC_MODE:
	    # Grab some values from *self*:
	    start = self.start_p
	    end = self.end_p
	    nut_height = self.nut_height_l - L(mm=0.38)
	    hex_nut_edge_width = self.hex_nut_edge_width_l - L(mm=0.28)
	    hex_nut_tip_width = self.hex_nut_tip_width_l
	    close_fit = self.close_fit_l
	    major_diameter = self.major_diameter_l

	    # Parse *flags*:
	    surface_normal = P(0, 0, L(mm=1))
	    for flag in flags.split(":"):
		if flag == 'N' or flag == 'S':
		    surface_normal = P(0, L(mm=1), 0)
		elif flag == 'E' or flag == 'W':
		    surface_normal = P(L(mm=1), 0, 0)
		elif flag == 'T' or flag == 'B':
		    pass
		else:
		    assert False, "Unrecognized flag '{0}'".format(flag)

	    # Compute the various screw axes:
	    screw_axis = (start - end).normalize()
	    edge_axis = screw_axis.cross_product(surface_normal)

	    # Generate the screw cut out block:

	    # Compute *p1* and *p2* which are the two diagonally opposite
	    # points on the screw removal cube:
	    d_edge = (close_fit / 2).millimeters()
	    p1 = end - \
	      edge_axis * d_edge - \
	      surface_normal * (major_diameter*2).millimeters()
	    p2 = start + \
 	      edge_axis * (close_fit/2).millimeters() + \
	      surface_normal * (major_diameter*2).millimeters()

	    # Compute the cube dimensions (*dx*, *dy*, *dz*) and the
	    # cube center point (*x*, *y*, *z*):
	    x = (p1.x + p2.x)/2
	    y = (p1.y + p2.y)/2
	    z = (p1.z + p2.z)/2
	    dx = (p2.x - p1.x).absolute()
	    dy = (p2.y - p1.y).absolute()
	    dz = (p2.z - p1.z).absolute()

	    # Print stuff out for debugging:
	    #print("Fastener:nut_ledge:close_fit={0:m} d_edge={1}".
	    #  format(close_fit, d_edge))
	    #print("Fastener:nut_ledge:start={0} end={1}".format(start, end))
	    #print("Fastener.nut_ledge:screw_axis={0} surf_norm={1} edge={2}".
	    #  format(screw_axis, surface_normal, edge_axis))
	    #print("Fastener:nut_ledge:p1={0} p2={1}".format(p1, p2))
	    #print("Fastener:nut_ledge:x={0} y={1} z={2}".format(x, y, z))
	    #print("Fastener:nut_ledge:dx={0} dy={1} dz={2}".format(dx, dy, dz))

	    # Output the OpenSCAD command to *difference_lines*:
	    difference_lines = part._scad_difference_lines
	    difference_lines.append(
	      "    // '{0}' nut Ledge screw remove".format(self.comment_s))
	    difference_lines.append(
	      "    translate([{0:m}, {1:m}, {2:m}])".format(x, y, z))
	    difference_lines.append(
	      "      cube([{0:m}, {1:m}, {2:m}], center=true);".
	      format(dx, dy, dz))

	    # Generate the hex cut out block:
	    # Compute the distances along the screw.
	    screw1 = nut_height.millimeters() / 2
	    screw2 = screw1 + nut_height.millimeters()
	    edge1 = (hex_nut_edge_width / 2).millimeters()
	    edge2 = -edge1
	    surface1 = (major_diameter * 2).millimeters()
	    surface2 = -surface1

	    # Compute the two cube corners *p1* and *p2*.
	    p1 = end + \
    	      screw_axis * screw1 - \
	      edge_axis * edge1 - \
	      surface_normal * surface1
	    p2 = end + \
    	      screw_axis * screw2 - \
	      edge_axis * edge2 - \
	      surface_normal * surface2

	    # Compute translate for cube center (*x*, *y*, *z*):
	    x = (p1.x + p2.x) / 2
	    y = (p1.y + p2.y) / 2
	    z = (p1.z + p2.z) / 2

	    # Compute dimensions of cube (*dx*, *dy*, *dz*):
	    dx = (p2.x - p1.x).absolute()
	    dy = (p2.y - p1.y).absolute()
	    dz = (p2.z - p1.z).absolute()

	    # Output the nut removal:
	    difference_lines.append(
	      "    // '{0}' Nut Ledge nut remove:".format(self.comment_s))
	    difference_lines.append(
	      "    translate([{0:m}, {1:m}, {2:m}])".format(x, y, z))
	    difference_lines.append(
	      "      cube([{0:m}, {1:m}, {2:m}], center=true);".
	      format(dx, dy, dz))
	    difference_lines.append("")
	

    def drill(self, part = None, select = None, trace = -1000000):
	""" *Fastener*: """

	if trace >= 0:
	    print("{0}=>Fastener.drill({1}, select='{2}', start={3}, end={4})".
	      format(' ' * trace, part, select, self.start_p, self.end_p))

	# Check argument types:
	none_type = type(None)
	assert type(part) == none_type or isinstance(part, Part)
	assert isinstance(select, str)

	# Drill the hole:
	if isinstance(part, Part):
	    start = self.start_p
            end = self.end_p

	    if select == "thread":
		material = self._material
		generic = material._generic
		if generic == "steel" or generic == "iron":
		    diameter = self.thread50_l
		else:
		    diameter = self.thread75_l
	    elif select == "close":
		diameter = self.close_fit_l
	    elif select == "free":
		diameter = self.free_fit_l
	    else:
		assert False, "select='{0}', not 'thread', 'close' or 'free'".\
		  format(select)

	    part.hole(comment = "'{0} Drill'".format(self.comment_s),
	      diameter = diameter,
	      start = start,
	      end = end,
	      flags = "t")

	    zero = L()
	    flat_head_point_angle = self.flat_head_point_angle_a
	    if flat_head_point_angle > Angle():
		# Most flat head screws have an 80 degree point angle.
		# Thus, in the crude diagram below.  The head diameter
		# (HD) is the distance AD which is the sum of AC and CD.
		# The point angle is <DBA which is the sum of <CBD and
		# <CBA = b + b = 80.  Thus, angle <CBA = b = 80/2 = 40 degrees.
		#
		#  --- A
		#   ^  |\
		#   |  |a\
		#   |  |  \
		#   |  |  b\---------------------+
		#  HD  C----B                    |
		#   |  |  b/---------------------+
		#   |  |  /
		#   |  |a/
		#   V  |/
		#  --- D
                #
		# What we need is the distance CB.  Triagle ABC is a right
		# triangle.  We know angle <CBA is 40 degrees, and thus angle
		# <CAB is 90 - 40 = 50 degrees.
		#
		# Let r = the length of the hyptoneuse (i.e. distance AB).
		# Let h = the length along the screw axis (i.e. distance BC).
		# Let w = the half width of the screw head (i.e. distance AC).
		#
		#    r * sin(a) = h              (1)
		#    r * sin(b) = w              (2)
		#    a + b = 90 degrees          (3)
		#
		# We know w and we know b, and we want h:
		#
		#    h/sin(a) = w/sin(b)         (4)
		#    h = w * sin(a)/sin(b)       (5)
		#    h = w * sin(90-b)/sin(b)    (6)

		flat_head_diameter = self.flat_head_diameter_l
		w = flat_head_diameter / 2
		b = flat_head_point_angle / 2
		h = w * ((Angle(deg = 90) - b).sine() / b.sine())
		direction = (start - end)
		if direction.length() <= zero:
		    print("Fastener.drill: part={0}: zero length screw head".
		      format(part))
		else:
		    normalized_direction = direction.normalize()
		    if trace >= 0:
			print("{0}head_radius={1} half_point_angle={2} h={3}".
			  format(' ' * trace, w, b, h))

		    flat_head_start = start - normalized_direction * h._mm
		    flat_head_end = start + normalized_direction * h._mm
		    #print(
                    #    "Part[{0}]:start={1} end={2} normalize={3} fh_end={4}".
		    #    format(part._name,
		    #    start, end, normalized_direction, flat_head_end))
		    part.hole(comment = "Flat Head:" + self.comment_s,
		      start_diameter = 2 * flat_head_diameter,
		      end_diameter = L(),
		      start = flat_head_start,
		      end = flat_head_end,
		      trace = trace + 1)

	    if self.hex_insert_b:
		direction = end - start
		direction_length = direction.length()
		if direction_length <= zero:
		    print("Fastener.drill(): part={0}: zero length screw".
		      format(part))
		else:
		    normalized_direction = direction.normalize()
		    nut_height = self.nut_height_l
		    insert_end = end - (normalized_direction * nut_height._mm)
		    #print("end={0} start={1} dir={2} dir_len={3} nut_hght={4}".
		    #  format(end, start, direction, direction_len, nut_hght))
		    #print("insert_end = {0}".format(insert_end))

		    part.hole(comment = "Hex Insert:" + self.comment_s,
		      diameter = self.hex_nut_tip_width_l,
		      start = end,
		      end = insert_end,
		      sides = 6,
		      sides_angle = self.sides_angle_a,
		      flags = "f")

	if trace >= 0:
	    print("{0}<=Fastener.drill({1}, select='{2}')".
	      format(' ' * trace, part, select))


# *Code* class:

# Used be called called *Code* in EZCAD.

class Code:

    def __init__(self):
	""" *Code*:
	"""

	zero = L()

	# Load up *self*:
	self._code_stream = None	# Code stream to write G-codes to
	self._command_started = False	# At begining of RS274 command line
	self._command_chunks = []	# RS274 line broken into space separated chunks
	self._comment_chunks = []	# RS274 comment broken into space separated chunks
	self._dxf = ""			# Text for dxf file
	self._dxf_x_offset = zero	# DXF X offset
	self._dxf_y_offset = zero	# DXF Y offset
	self._vice_x = zero
	self._vice_y = zero
	self._z_rapid = zero		# Z above which Z rapids are allowed
	self._z_safe = zero		# Z above XY rapids are safe

	self.vrml_colors = [ 1.0, 0.0, 0.0 ]
	self.vrml_color_indexes = []
	self.vrml_lines = []
	self.vrml_points = []

	# Construct the *g1_table*:
	g1_values_list = (0, 1, 2, 3, 33, 38, 73, 76, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89)
	g1_table = {}
	for g1_value in g1_values_list:
	    g1_table[g1_value] = g1_value
	self._g1_table = g1_table

	# The stuff below is RS-274 mode variables:
	self._f = Speed()		# Feedrate
	self._g1 = 0			# (G0-3, 33, 38.x, 73, 76, 80-89)
	self._g2 = 0			# (G17-19)
	self._g3 = 0			# (G7-8)
	self._g4 = 0			# (G90-91)
	self._g5 = 0			# (G93-94)
	self._g6 = 0			# (G20-21)
	self._g7 = 0			# (G40-42)
	self._g8 = 0			# (G43, 49)
	self._g9 = 0			# (G98-99)
	self._g10 = 0			# (G54-59)
	self._g11 = 0			# (G4)
	self._h = zero			# H tool offset index
	self._i = zero			# I coordinate
	self._j = zero			# J coordinate
	self._m1 = 0			# (M0-2, 30, 60)
	self._m2 = 0			# (M6)
	self._m3 = 0			# (M3-5)
	self._m4 = 0			# (M7-9)
	self._m5 = 0			# (M48-49)
	self._p = Time()		# G4
	self._q = zero			# Peck depth
	self._r0 = zero			# Radius cycle R
	self._r1 = zero			# Drill cycle R
	self._s = Hertz()		# Spindle revolutions
	self._x = zero			# X coordinate
	self._y = zero			# Y coordinate
	self._z = zero			# Z coordinate
	self._z1 = zero 		# Z coordinate
	self._z_safe_f = Speed()	# Feed to perform z safe operation at
	self._z_safe_pending = False	# {true}=>need to do z safe move
	self._z_safe_s = Hertz()	# Speed to perform z safe operation at

	self._vrml_reset()

    def _block_append(self, block):
	""" *Code*: Append *block* onto the blocks list inside of the *Code* object
	    (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(block, Block)

	# Append *block* onto the blocks list:
	self._blocks.append(block)

    def _command_begin(self):
	""" *Code*: Start a new RS274 command in the *Code* object (i.e. *self*). """

	assert not self._command_started, "Previous RS274 command was not ended."
	self._command_started = True

	# Remember some values for VRML line path drawing:
	self._vrml_start_r0 = self._r0
	self._vrml_start_x = self._x_value()
	self._vrml_start_y = self._y_value()
	self._vrml_start_z = self._z
	self._vrml_motion = -1

    def _command_end(self):
	""" *Code*: End the current RS274 in the *Code* object (i.e. *self*). """

	# Make sure we have started a command:
	assert self._command_started, "Not currently in a RS274 command"

	# Grab out both the *command_chunks* and the *comment_chunks*:
	command_chunks = self._command_chunks
	comment_chunks = self._comment_chunks

	# If we have any *comment_chunks*, tack the onto the end of *command_chunks*:
	if len(comment_chunks):
	    command_chunks += ["("] + comment_chunks + [")"]

	# Construct a space separated *command*:
	command = " ".join(self._command_chunks)
	#print("command='{0}'".format(command))

	# Write *comand* to *code_stream*:
	code_stream = self._code_stream
	code_stream.write(command)
	code_stream.write("\n")

	# Clear out *command_chunks* and *comment_chunks* for the next command:
	del command_chunks[:]
	del comment_chunks[:]

	# Mark that we ended the current command:
	self._command_started = False

	color_index = -1
	if self._vrml_motion == 0:
	    color_index = 0
	elif self._vrml_motion == 1:
	    color_index = 1

	if color_index >= 0:
	    start_index = \
	      self._vrml_point(self._vrml_start_x, self._vrml_start_y, self._vrml_start_z)
	    end_index = self._vrml_point(self._x_value(), self._y_value(), self._z)
	    self._vrml_point_indices.append(start_index)
	    self._vrml_point_indices.append(end_index)
	    self._vrml_point_indices.append(-1)
	    self._vrml_color_indices.append(color_index)

    def _comment(self, comment):
	""" *Code*: Output *comment* to the *Code* object (i.e. *self*). Any parentheses
	    in *comment* are converted to square brackets.
	"""

	# Verify argument types:
	assert isinstance(comment, str)

	# Replace open/close parenthesis in *comment* with square brackets:
	comment = comment.replace('(', '[').replace(')', ']')

	# Add it on to the comment
	self._comment_chunks.append(comment)


    def _configure(self, tool, vice_x, vice_y):
	""" *Code*: Configure the *Code* object (i.e. *self*) to use
	    *tool*, *vice_x*, and *vice_y*.
	"""

	# Verify argument types:
	assert isinstance(tool, Tool)
	assert isinstance(vice_x, L)
	assert isinstance(vice_y, L)

	# Initialize *code*:
	self._tool = tool
	self._vice_x = vice_x
	self._vice_y = vice_y

    def _dxf_xy_offset_set(self, dxf_x_offset, dxf_y_offset):
	""" *Code*: Set the DXF offset X field of the *RS724* object (i.e. *self*) """

	# Verify argument types:
	assert isinstance(dxf_x_offset, L)
	assert isinstance(dxf_y_offset, L)
	self._dxf_x_offset = dxf_x_offset
	self._dxf_y_offset = dxf_y_offset

    def _dxf_y_offset_get(self):
	""" *Code*: Return the DXF offset Y field of the *Code* object (i.e. *self*) """

	return self._dxf_y_offset


    def _finish(self):
	""" *Code*: Finish off the current block of G code. """

	code_stream = self._code_stream
	assert isinstance(code_stream, file)
	code_stream.close()
	self._code_stream = None

	# Random comment: the view3dscene can view the resulting .wrl file:

	# Write out headers to *vrml_stream*:
	vrml_stream = self._vrml_stream
	vrml_stream.write("#VRML V2.0 utf8\n")
	vrml_stream.write("Shape {\n")
	vrml_stream.write(" geometry IndexedLineSet {\n")
	vrml_stream.write("  colorPerVertex FALSE\n")

	# Output the colors:
	vrml_stream.write("  color Color {\n")
	vrml_stream.write("   color [\n")
	vrml_stream.write("     0.0 0.0 1.0 # blue\n")
	vrml_stream.write("     1.0 0.0 0.0 # red\n")
	vrml_stream.write("     0.0 1.0 0.0 # green\n")
	vrml_stream.write("   ]\n")
	vrml_stream.write("  }\n")

	# Write out points:
	vrml_stream.write("  coord Coordinate {\n")
	vrml_stream.write("   point [\n")
	for point in self._vrml_points:
	    vrml_stream.write("    {0} {1} {2}\n".format(point[0], point[1], point[2]))
	vrml_stream.write("   ]\n")
	vrml_stream.write("  }\n")
		
	# Write out coordinate index:
	vrml_stream.write("  coordIndex [\n")
	vrml_point_indices = self._vrml_point_indices
	vrml_stream.write("  ");
	for index in vrml_point_indices:
	    vrml_stream.write(" {0}".format(index))
	    if index < 0:
		vrml_stream.write("\n  ")
	vrml_stream.write("]\n")

	# Write out the color indices:
	vrml_stream.write("  colorIndex [\n")
	for color_index in self._vrml_color_indices:
	    vrml_stream.write("    {0}\n".format(color_index))
	vrml_stream.write("  ]\n")

	# Close out the shape and geometry clauses:
	vrml_stream.write(" }\n")
	vrml_stream.write("}\n")

	# Close *vrml_stream*:
	vrml_stream.close()

	# Now reset all the VRML values:
	self._vrml_reset()

    def _flush(self, program_number):
	""" *Code*: Generate CNC code starting at *program_number* for
	    the first CNC file name.
	"""

	# Verify argument types:
	assert isinstance(program_number, int)

	# Use *code* instead of *self*:
	code = self

	#call d@(form@("=>flush@Code(*, %d%)\n\") / f@(program_number))
	#original_program_number :@= program_number

	# Output all of the *blocks*:
	blocks = self._blocks
	size = len(blocks)
	print("size=", size)
	if size != 0:
	    block0 = blocks[0]
	    part = block0._part_get()
	    assert isinstance(part, Part)

	    remainder = program_number % 10
	    if remainder != 0:
		program_number += 10 - remainder
	    directory = ezcad._directory_get()
	    file_name = os.path.join(directory, "{0}.ngc".format(program_number))

	    # Assign program numbers:
	    for index in range(size):
		block = blocks[index]
		#call d@(form@("Block[%d%] #%d% %d% %v%\n\") % f@(index) %
		#  f@(block.uid) % f@(block.priority) / f@(block.tool.name))
		block.program_number = program_number
		program_number += 1

	    # Open *file_name* for writing:
	    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>code_stream_open: '{0}'".format(file_name))
	    code_stream = open(file_name, "wa")
	    assert isinstance(code_stream, file), \
	      "Could not open '{0}' for writing".format(file_name)

	    # Output the tool table:
	    code_stream.write("({0}: {0})\n".format(part._name_get(), block0._comment_get()))
	    code_stream.write("(Tooling table:)\n")
            for block in blocks:
		# Fetch {tool}:
		tool = block._tool_get()

		# Make sure we have assigned a number to the tool:
		tool_name = tool._name_get()
		assert tool._number_get() != 0, \
		  "Tool {0} does not have a number\n".format(tool_name)

		# Output the line in the tool table:
		code_stream.write("(T{0} {1})\n".format(tool._number_get(), tool_name))

		# Perform the code block:
		code_stream.write("O{0} call\n".format(block._program_number_get()))

	    # Wrap it up:
	    code_stream.write("G53 Y0.0 (Move the work to the front)\n")
	    code_stream.write("M2\n")
	    code_stream.close()
	    print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<code_stream_open: '{0}'".format(file_name))

	    # Now output the blocks:
	    for block in blocks:
		# Fetch {block}, {tool} and {name}.
		tool = block._tool_get()
		tool_number = tool._number_get()

		if tool._is_laser_get():
		    # We have a laser; generate a .dxf file:
		    directory = ezcad._directory_get()
		    dxf_file_name = os.path.join(directory, "{0}.dxf".format(program_number))
		    dxf_stream = open(dxf_file_name, "wa")
		    assert isinstance(dxf_stream, file), \
		      "Unable to open file '{0}'".format(file_name)

		    #dxf_stream.write("0\n\SECTION\n\2\n\HEADER\n")
		    #dxf_stream.write("9\n\$DIMAUNITS\n\70\n\1\n")
		    #dxf_stream.write("9\n\$INSUNITS\n\70\n\1\n")
		    #dxf_stream.write("9\n\$LUNITS\n\70\n\0\n")
		    #dxf_stream.write("9\n\$MEASUREMENT\n\70\n\0\n")
		    #dxf_stream.write("0\n\ENDSEC\n")

		    dxf_stream.write("0\n\SECTION\n\2\n\ENTITIES\n")

		    # Output *dxf* to *dxf_stream*:
		    dxf = code.dxf
		    dxf_stream.write(dxf)

		    # See if we are aggregating *dxf* into a named DXF file:
		    dxf_base_name = part.dxf_base_name
		    if dxf_base_name != "":
			# Yes, we need to aggregate:
			shop = part._shop
			dxf_table = shop.dxf_table
			dxf_contents = lookup[dxf_base_name]
			assert isinstance(dxf_contents, str)
			dxf_contents[dxf_base] = dxf_contents + dxf

		    # Clear *dxf* for the next part:
		    code.dxf = ""

		    # Wrap up {dxf_stream}:
		    dxf_stream.write("0\n\ENDSEC\n\0\n\EOF\n")
		    dxf_stream.close()
		else:
		    # We have a mill; generate a .ngc file:
		    program_number = block.program_number
		    directory = ezcad._directory_get()
		    ngc_file_name = os.path.join(directory, "{0}.ngc".format(program_number))
		    code_stream = open(ngc_file_name, "wa")
		    print(">>>>>>>>>>>>>>>>code_stream_open: '{0}'".format(file_name))
		    assert isinstance(code_stream, file), \
		      "Unable to open '{0}' for output\n".format(ngc_file_name)

		    # Put a visual break into {code_stream}:
		    part_name = block._part_get()._name_get()
		    code_stream.write("({0})\n".format(part_name))
		    code_stream.write("O{0} sub\n".format(program_number))
		    code_stream.write("G90 G90.1 G20\n")

		    # Output the tool change, get the coolant on,
		    # and spindle spun up:
		    code_stream.write("M6 T{0} (Insert {1})\n". 
           	      format(tool_number, tool._name_get()))
		    spindle = block._spindle_get()
		    if spindle > Hertz():
			material = part._material_get()
			if material._needs_coolant():
			    code_stream.write("M8 (Coolant on)\n")
			else:
			    code_stream.write("M9 (Coolant off)\n")
			code_stream.write("S{0:rpm} M3 (Spindle on)\n".format(spindle))

		    # Get moving to tool change location:
		    plunge_x = part._plunge_x_get()
		    plunge_y = part._plunge_y_get()
		    vice_x = block._vice_x_get()
		    vice_y = block._vice_y_get()
		    code._vice_xy_set(vice_x, vice_y)

		    #print("part:{0} code.vice_x={1} code.vice_y={2}\n".
		    #  format(part.name, code.vice_x, code.vice_y))

		    code_stream.write("G54 G0 X{0} Y{1} (Apply work offset)\n".
		      format(plunge_x - vice_x, plunge_y - vice_y))
		    code_stream.write("(Part_Origin_X={0} Part_Origin_Y={1})\n".
		      format(-vice_x, -vice_y))

		    # Enable tool offset:
		    z_safe = part._z_safe_get()
		    code_stream.write("G43 H{0} (Enable tool_offset)\n".
		      format(tool_number))
		    code_stream.write(
		      "G0 Z{0:i} (Go to Z-safe for current tool)\n".
		      format(z_safe))

		    # These should be done using accessor functions:
		    code._g8 = 43
		    code._g1 = 0
		    code._h = tool_number
		    code._z = z_safe

		    # Generate the code for {block}:
		    code._commands_write(code_stream)

		    # Put the tool back into a reasonable position:
		    code_stream.write("M5 (Spindle off)\n")
		    code_stream.write("M9 (Coolant off)\n")
		    code_stream.write(
		      "G0 X{0} Y{1} (Put spindle to left of vice)\n".
		      format(0.0, -1.0))

		    # These should be done using accessor functions:
		    code._x = vice_x - L(inch= 1.0)
		    code._y = vice_y

		    code_stream.write("G49 (Disable tool offset)\n")
		    # G53 is not modal, thus we are still in G54 work offset:
		    code_stream.write(
		      "G53 G0 Z{0} (Return to machine Z zero)\n".format(1.0))
		    code.z = L()
		    code_stream.write("O{0} endsub\n".format(program_number))
		    code_stream.write("O{0} call\n".format(program_number))
		    code_stream.write("M2\n")
		    code_stream.write("%\n")
		    code_stream.close()
		    print("<<<<<<<<<<<<<<<<code_stream_open: '{0}'".format(file_name))

	    # Output the shut-down code:
	    #call put@("(*************************************)\n", code_stream)

	    # Output the subroutine calls:
	    #index := 0
	    #while index < size
	    #    #call put@(form@("O%d% call\n\") /
	    #    #  f@((index + 1) * 100), code_stream)
	    #    index := index + 1

	    #call put@(form@("O<%s%> sub\n\") / f@(base_name), code_stream)
	    #call put@("O[#1 * 100] call\n\", code_stream)
	    #call put@(form@("O<%s%> endsub\n\") / f@(base_name), code_stream)

	    #call put@("G53 Y0.0 (Move the work to the front)\n", code_stream)

	    #call put@("M2\n\", code_stream)
    
	    # Zero {blocks} for the next time:
	    del blocks[:]

	    #call close@(code_stream)

	#call d@(form@("<=flush@Code(*, %d%) => %d%\n\") %
	#  f@(original_program_number) / f@(program_number))
	return program_number

    def _g1_set(self, g1):
	""" *Code*: Set the G1 field of the *Code* object (i.e. *self*) to *g1*. """

	# Verify argument types:
	assert isinstance(g1, int)

        # Set the G1 field of the *Code* object (i.e. *self*) to *g1*:
	self._g1 = g1

    def _hertz(self, field_name, value):
	""" *Code*: Output *value* for *field_name* in the current command of the *Code* object 
	    (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(field_name, str)
	assert isinstance(value, Hertz)

	square_bracket = False
	changed = False
	if field_name == "S":
	    if self._s != value:
		self._s = value
		changed = True
	elif field_name == "[]":
	    changed = True
	    square_bracket = True
	else:
	    assert False, "Bad field name '{0}'".format(field_name)

	if changed:
	    if square_bracket:
		chunk = "[{0:rpm}]".format(value)
	    else:
		chunk = "{0}{1:rpm}".format(field_name[0], value)
	    self._command_chunks.append(chunk)

    def _length(self, field_name, value):
	""" *Code*: Output *value* for *field_name* to the *Code* object (i.e. *self*.) """

	assert isinstance(field_name, str)
	assert isinstance(value, L)

	square_brackets = False
	offset = L()
	changed = False
	if field_name == "Q":
	    if self._q != value:
		self._q = value
		changed = True
	elif field_name == "I":
            if self._i != value:
		self._i = value
		changed = True
	elif field_name == "J":
	    if self._j != value:
		self._j = value
		changed = True
	elif field_name == "R0":
	    if self._r0 != value:
		self._r0 = value
		changed = True
 	elif field_name == "R1":
	    if self._r1 != value:
		self._r1 = value
		changed = True
	elif field_name == "X":
	    offset = self._vice_x
	    if self._x != value - offset:
		self._x = value - offset
		changed = True
	elif field_name == "Y":
	    offset = self._vice_y
	    if self._y != value - offset:
		self._y = value - offset
		changed = True
	elif field_name == "Z":
	    if self._z != value:
		self._z = value
		achanged = True
	elif field_name == "Z1":
	    if self._z1 != value:
		self._z1 = value
		changed = True
	elif field_name == "[]":
	    changed = True
	    square_brackets = True
	else:
	    assert False, "Unrecognized field {0}".format(field_name)

	if changed:
            chunk = ""
	    if square_brackets:
		chunk = "[{0:i}]".format(value - offset)
	    else:
		chunk = "{0}{1:i}".format(field_name[0], value - offset)
	    self._command_chunks.append(chunk)

    def _line_comment(self, comment):
	""" *Code*: Emit *comment* to the *Code* object (i.e. *self*). """

	# Verify argument types:
	assert isinstance(comment, str)

	# This routine will add a line containg *comment* to the *Code* object (i.e. *self).:
	self._command_begin()
	self._comment(comment)
	self._command_end()

    def _mode_motion(self, g1_value):
	""" *Code*: This routine will issue a G*g1_value* motion command to the *Code* object
	    (i.e. *self)*.
	"""

	# Verify arguement types:
	assert isinstance(g1_value, int) and g1_value in self._g1_table

	# Add a G1 field to the command:
	self._unsigned("G1", g1_value)

	self._vrml_motion = g1_value


    def _reset(self):
	""" *Code*: Reset the contents of the *Code* object (i.e. *self*) """

	zero = L()
	large = L(inch=123456789.0)
	huge = 0x7fffffff
	big = L(inch=123456789)

	self._begin = True
	self._f = Speed(mm_per_sec=huge)
	self._g1 = huge
	self._g2 = huge
	self._g3 = huge
	self._g4 = huge
	self._g5 = huge
	self._g6 = huge
	self._g7 = huge
	self._g8 = huge
	self._g9 = huge
	self._g10 = huge
	self._g11 = huge
	self._h = huge
	self._i = big
	self._j = big
	self._m1 = huge
	self._m2 = huge
	self._m3 = huge
	self._m4 = huge
	self._m5 = huge
	self._p = Time(sec=-1.0)
	self._q = big
	self._r0 = big
	self._r1 = big
	self._s = Hertz(rps=123456789.0)
	self._x = big
	self._y = big
	self._z = big
	self._z1 = big

    def _start(self, part, tool, ngc_program_number, spindle_speed):
	""" *Code*: Start writing out the G-code for *tool* (
	"""

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(tool, Tool)
	assert isinstance(ngc_program_number, int) and ngc_program_number > 0
	assert isinstance(spindle_speed, Hertz)

        # Grab some values from *part* and *tool*:
        part_name = part._name_get()
	tool_name = tool._name_get()
	tool_number = tool._number_get()

	# Make sure that closed off any previous *code_stream*:
	assert self._code_stream == None

	# Open new *code_stream*:
	ezcad = part._ezcad_get()
	directory = ezcad._directory_get()
	code_file_name = os.path.join(directory, "O{0}.ngc".format(ngc_program_number))
	code_stream = open(code_file_name, "w")
	assert code_stream != None, "Could not open '{0}' for writing".format(code_file_name)
	self._code_stream = code_stream

	# Output a descriptive header comment:
	code_stream.write("( Part {0}: Tool {1} Program: {2} )\n".
	  format(part_name, tool_name, ngc_program_number))

	# Declare the the subroutine number:
	code_stream.write("O{0} sub\n".format(ngc_program_number))

	# Set absolute mode for X/Y/Z (G90) and absolute mode for I/J/K (G90.1).
	# Set units to inches (G20):
	code_stream.write("G90 G90.1 G20\n")

	# Output the tool change, get the coolant on, and spindle spun up:
	code_stream.write("M6 T{0} (Insert {1})\n".format(tool_number, tool_name))
	spindle_is_on = spindle_speed > Hertz()
	if spindle_is_on:
	    code_stream.write("( Get the spinde up to speed.)\n")
	    code_stream.write("S{0:rpm} M3 (Spindle on)\n".format(spindle_speed))

	# Use the *material* to control the coolant:
	material = part._material_get()
	if material._needs_coolant() and spindle_is_on:
	    code_stream.write("M8 (Coolant on)\n")
	else:
	    code_stream.write("M9 (Coolant off)\n")

	# Open new *vrml_stream*:
	directory = ezcad._directory_get()
	vrml_file_name = os.path.join(directory, "O{0}.wrl".format(ngc_program_number))
	vrml_stream = open(vrml_file_name, "w")
	assert vrml_stream != None, "Could not open '{0}' for writing".format(vrml_file_name)
	self._vrml_stream = vrml_stream

    def _foo(self):
	# Get moving to tool change location:
	plunge_x = part._plunge_x_get()
	plunge_y = part._plunge_y_get()
	vice_x = part._vice_x_get()
	vice_y = part._vice_y_get()
	code._vice_xy_set(vice_x, vice_y)

	#print("part:{0} code.vice_x={1} code.vice_y={2}\n".
	#  format(part.name, code.vice_x, code.vice_y))

	code_stream.write("G54 G0 X{0} Y{1} (Apply work offset)\n".
	  format(plunge_x - vice_x, plunge_y - vice_y))
	code_stream.write("(Part_Origin_X={0} Part_Origin_Y={1})\n".
	  format(-vice_x, -vice_y))

	# Enable tool offset:
	z_safe = part._z_safe_get()
	code_stream.write("G43 H{0} (Enable tool_offset)\n".
	  format(tool_number))
	code_stream.write(
	  "G0 Z{0:i} (Go to Z-safe for current tool)\n".
	  format(z_safe))

	# These should be done using accessor functions:
	code._g8 = 43
	code._g1 = 0
	code._h = tool_number
	code._z = z_safe

    def _s_set(self, s):
	""" *Code*: Set the S field of the *Code* object (i.e. *self*) to *s*. """

	# Verify argument types:
	assert isinstance(s, Hertz)

        # Set the S field of the *Code* object (i.e. *self*) to *s*:
	self._s = s

    def _simple_pocket_helper(self, pocket, offset, s, f, z, rapid_move):
	""" *Code*: Perform one rectangular or rounded rectangular path of the currently
	    mounted tool and output the commands to the *Code* object (i.e. *self*).
	    *offset* specifies the distance inward from the nominal path specified by
	    *pocket*.  *z* specifies the depth at which the tool should be placed before
	    traversing the path.  *rapid_move* specifies whether or not to move to the
	    start position with rapid tool movement (i.e. G0) or linear tool movement
	    (i.e. G1) commands.  For linear tool movement commands *s* is the spindle
	    speed and *f* is the feedrate to use.
	"""
    
	# Verify argument types:
	#assert isinstance(block, Block)
	assert isinstance(pocket, Operation_Simple_Pocket)
	assert isinstance(offset, L)
	assert isinstance(s, Hertz)
	assert isinstance(f, Speed)
	assert isinstance(z, L)
	assert isinstance(rapid_move, bool)

	# Extract the corners:
	x1 = pocket._x1_get()
	y1 = pocket._y1_get()
	x2 = pocket._x2_get()
	y2 = pocket._y2_get()
	corner_radius = pocket._corner_radius_get()
	self._line_comment(
	  "x1={0:i} y1={1:i} x2={2:i} y2={3:i} cr={4:i}".format(x1, y1, x2, y2, corner_radius))
    
	# Set *debug* to *True* to turn on tracing:
	debug = False
	#debug = True
	if debug:
	    print("Code._simple_pocket_helper: x1={0:i} y1={1:i} x2={2:i} y2={3:i}".
	      format(x1, y1, x2, y2))

	# Make sure that {x1} < {x2} and {y1} < {y2}:
	assert x1 < x2
	assert y1 < y2
    
	# Compute the corner center coordinates:
	rx1 = x1 + corner_radius
	rx2 = x2 - corner_radius
	ry1 = y1 + corner_radius
	ry2 = y2 - corner_radius
	self._line_comment("rx1={0:i} ry1={1:i} rx2={2:i} ry2={3:i}".format(rx1, ry1, rx2, ry2))
	if debug:
	    print("rx1={0:i} ry1={1:i} rx2={2:i} ry2={3:i}".format(rx1, ry1, rx2, ry2))

	# Compute the rectangle path coordinates:
	px1 = x1 + offset
	px2 = x2 - offset
	py1 = y1 + offset
	py2 = y2 - offset
	self._line_comment(
	  "offset={0:i} px1={1:i} py1={2:i} px2={3:i} py2={4:i}".format(offset, px1, py1, px2,py2))
    
	# Determine the starting location for this path:	
	start_x = px1
	start_y = py1
	if offset < corner_radius:
	    # We have rounded corner path:
	    start_x = rx1
    
	# Move to (*start_x*, *start_y*) as specified by *linear_move* argument:
	if rapid_move:
	    self._xy_rapid(start_x, start_y)
	else:
	    self._xy_feed(f, s, start_x, start_y)
    
	# Make sure we are at the depth *z*:
	self._z_feed(f/2, s, z, "simple_pocket_helper")
    
	# Mill out either a square or rounded corners
	if offset < corner_radius:
	    # Mill out a rectangle with rounded corners in a
	    # counter clockwise direction to force a climb cut:
    
	    r = corner_radius - offset
    
	    # Bottom horizontal line from (rx1,py1) to (rx2,py1):
	    self._xy_feed(f, s, rx2, py1)
    
	    # Lower right arc (rx2,py1) to (px2,ry1):
	    self._xy_ccw_feed(f, r, s, px2, ry1, rx=rx2, ry=ry1)
    
	    # Right vertical line (px2,ry1) to (px2,ry2):
	    self._xy_feed(f, s, px2, ry2)
    
	    # Upper right arc (px2,ry2) to (rx2, py2):
	    self._xy_ccw_feed(f, r, s, rx2, py2, rx=rx2, ry=ry2)
    
	    # Top horizontal line (rx2, py2) to (rx1, py2):
	    self._xy_feed(f, s, rx1, py2)
    
	    # Upper left arc (rx1, py2) to (px1, ry2):
	    self._xy_ccw_feed(f, r, s, px1, ry2, rx=rx1, ry=ry2)
    
	    # Left vertical line (px1, ry2) to (px1, ry1):
	    self._xy_feed(f, s, px1, ry1)
    
	    # Lower left arc (px1, ry1) to (rx1, py1):
	    self._xy_ccw_feed(f, r, s, rx1, py1, rx=rx1, ry=ry1)
	else:
	    # Mill out a rectangle with "square" corners in a counter
	    # clockwise direction to force a climb cut:
    
	    # Bottom horizontal line from (px1, py1) to (px2, py1):
	    self._xy_feed(f, s, px2, py1)
    
	    # Right vertical line from (px2, py1) to (px2, py2):
	    self._xy_feed(f, s, px2, py2)
    
	    # Top horizontal line from (px2, py2) to (px1, py2):
	    self._xy_feed(f, s, px1, py2)
    
	    # Left vertical line from (px1, py2) to (px1, py1):
	    self._xy_feed(f, s, px1, py1)

    def _speed(self, field_name, value):
	""" *Code*: Set the speed for *field_name* to *value in the the current command of the
	    *Code* object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(field_name, str)
	assert isinstance(value, Speed)

	square_bracket = False
	changed = False
	if field_name == "F":
	    if self._f != value:
		self._f = value
		changed = True
	elif field_name == "[]":
	    square_bracket = True
	    changed = True
	else:
	    assert False, "field_name={0}".format(field_name)

	if changed:
	    if square_bracket:
		chunk = "[{0:i}]".format(value)
	    else:
		chunk = "{0}{1:i}".format(field_name[0], value)
	    self._command_chunks.append(chunk)

    def _old_code(self):

	self._reset()

	# Use *part* instead of *self*:
	part = self

	#call d@(form@("=>flush@Part(%v%, %d%)\n\") %
	#  f@(part.name) / f@(program_number))
	#original_program_number = program_number

	#call show@(part, "before flush")
	shop = part._shop_get()
	assert isinstance(shop, Shop)
	vice = shop._vice_get()
	jaw_width = vice._jaw_width_get()

	# Compute (*plunge_x*, *plunge_y*) which is the vertical axis over
	# which is to the left of the the part or the vice:
	edge_x = part._edge_x_get()
	edge_y = part._edge_y_get()
	vice_x = part._vice_x_get()
	vice_y = part._vice_y_get()
	assert isinstance(vice_x, L)
	plunge_x = vice_x
	if plunge_x > edge_x:
	    plunge_x = edge_x
	assert isinstance(jaw_width, L)
	assert isinstance(plunge_x, L)
	if plunge_x > jaw_width:
	    plunge_x = plunge_x - L(inch=0.7)
	part._plunge_xy_set(plunge_x, edge_y)
    
	code = shop._code_get()
	assert isinstance(code, Code)
	code._z_safe_set(part._z_safe_get())
	code._z_rapid_set(part._z_rapid_get())
	operations = part._operations_get()
	size = len(operations)

	#call show@(part, "before sort", 1t)

	# FIXME move the sort into *operations_regroup*:
	# Sort *operations* to group similar operations together:
	operations.sort(cmp=Operation._compare)

	#call show@(part, "after sort", 1t)

	# Do we need to user regroups instead *operations*:
	part._operations_regroup()

	#if regroups != 0
	#	call show@(part, "after_regroup", 0f)

	# FIXME: The code below should replace all the index stuff:
	current_tool = None
	operation_group = None
	operation_groups = []
	for operation in operations:
	    # Grab the *tool* from *operation*:
	    tool = operation._tool_get()
	    assert isinstance(tool, Tool)

	    # Start a new *operation_group* if the *tool* is different:
	    if current_tool != tool:
		# This code is always exectuted the first time through:
		current_tool = tool
		operation_group = []
		operation_groups.append(operation_group)

	    # Tack *operation* onto *operation_group*:
	    assert isinstance(operation_group, list)
	    operation_group.append(operation)

	# Open the top-level *part_ngc_stream* file that invokes each tool operation
	# in a separate .ngc file:
	directory = ezcad._directory_get()
	part_ngc_stream_file_name = os.path.join(directory, "O{0}.ngc".format(program_number))
	part_ngc_stream = open(part_ngc_stream_file_name, "w")
	assert part_ngc_stream != None, "Unable to open {0}".format(part_ngc_stream_file_name)

    def _unsigned(self, field_name, value):
	""" *Code*: This routine will format {value} for {field_name} to {code}.
	"""

	#FIXME: This should probably be done using a Python dictionary:

	# Verify argument types:
	assert isinstance(field_name, str)
	assert isinstance(value, int) and value >= 0

	previous_value = 0xffffffff
	matched = False
	size = len(field_name)
	if size == 1:
	    matched = True
	    letter = field_name[0]
	    if letter == 'H':
		previous_value = self._h
		if previous_value != value:
		    self._h = value
	    elif letter == 'M':
		pass
	    elif letter == 'O':
		pass
	    elif letter == 'T':
		pass
	    else:
		matched = False
	elif size == 2:
	    if field_name[0] == 'G':
		# Assume we will match (flip back in default clause below):
		matched = True
		g_digit = field_name[1]
		if g_digit == '1':
		    previous_value = self._g1
		    if previous_value != value:
			self._g1 = value
		elif g_digit == '2':
		    previous_value = self._g2
		    if previous_value != value:
			self._g2 = value
		elif g_digit == '3':
		    previous_value = self._g3
		    if previous_value != value:
			self._g3 = value
		elif g_digit == '4':
		    previous_value = self._g4
		    if previous_value != value:
			self._g4 = value
		elif g_digit == '5':
		    previous_value = self._g5
		    if previous_value != value:
			self._g5 = value
		elif g_digit == '6':
		    previous_value = self._g6
		    if previous_value != value:
			self._g6 = value
		elif g_digit == '7':
		    previous_value = self._g7
		    if previous_value != value:
			self._g7 = value
		elif g_digit == '8':
		    previous_value = self._g8
		    if previous_value != value:
			self._g8 = value
		elif g_digit == '9':
		    previous_value = self._g9
		    if previous_value != value:
			self._g9 = value
		else:
		    # We did not match:
		    matched = False
	elif size == 2:
	    if field_name == "G10":
		matched = True
		previous_value = self.g10
		if previous_value != value:
		    self._g10 = value
	    elif field_name == "G11":
		# Always output G4 when requested:

		matched = True
		previous_value = 0x12345

	assert matched, "Unrecognized field name {0}".format(field_name)
	if previous_value != value:
	    self._command_chunks.append("{0}{1}".format(field_name[0], value))

    def _vice_xy_set(self, vice_x, vice_y):
	""" *Code*: Set the vice X/Y fields of the *Code* object (i.e. *self*)
	    to *vice_x* and *vice_y*.
	"""

	# Verify argument types:
	assert isinstance(vice_x, L)
	assert isinstance(vice_y, L)

        # Set the vice X/Y fields of the *Code* object (i.e. *self*) to *vice_x* and *vice_y*:
	self._vice_x = vice_x
	self._vice_y = vice_y

    def _vrml_arc_draw(self, ax, ay, bx, by, radius, z, clockwise, radius_x=None, radius_y=None):
	""" *Code*: Draw an arc from (*ax*, *ay*, *z) to (*bx*, *by*, *bz*) with an
	    arc radius of *radius*.  The arc is drawn in a clockwise direction if
	    *clockwise* is *True* and a counter-clockwise direction otherwise. """
    
	# Verify argument types:
	assert isinstance(ax, L)
	assert isinstance(ay, L)
	assert isinstance(bx, L)
	assert isinstance(radius, L)
	assert isinstance(z, L)
	assert isinstance(clockwise, bool)
	assert isinstance(radius_x, L) or radius_x == None
	assert isinstance(radius_y, L) or radius_y == None

	# Enable tracing by setting *debug* to *True*:
	debug = False
	#debug = True
	have_radius = isinstance(radius_x, L) and isinstance(radius_y, L)
	if debug:
	    print("\nCode._vrml_arc_draw({0:i}, {1:i}, {2:i}, {3:i}, {4:i}, {5:i}, {6})".
	      format(ax, ay, bx, by, radius, z, clockwise))
	    if have_radius:
		print("radius_x={0:i} radius_y={1:i}".format(radius_x, radius_y))

		ar_dx = ax - radius_x
		ar_dy = ay - radius_y
		ar_radius = ar_dx.distance(ar_dy)
		print("ar_dx={0:i} ar_dy={1:i} ar_radius={2:i}".format(ar_dx, ar_dy, ar_radius))

		br_dx = bx - radius_x
		br_dy = by - radius_y
		br_radius = br_dx.distance(br_dy)
		print("br_dx={0:i} br_dy={1:i} br_radius={2:i}".format(br_dx, br_dy, br_radius))

	# Compute some *Angle* constants:
	degrees0 = Angle(deg=0.0)
	degrees15 = Angle(deg=15.0)
	degrees90 = Angle(deg=90.0)
	degrees180 = Angle(deg=180.0)
	degrees360 = Angle(deg=360.0)
    
	# Assign *z* to *az* and *bz* for notational consistency.
        az = z
	bz = z	

	# The goal of the code below is to draw an arc of radius *radius* from the point
	# A=(*ax*,*ay*,*az*) to point B=(*bx*,*by*,*bz).  The center point of the arc is
	# R=(*rx*,*ry*,*rz*).  The goal of the code below is to compute the location of R.
	# We assume that everything is occuring in the X/Y plane (i.e. *az*=*bz*=*rz*.)
	# The diagram shows two right triangles ACR and BCR, where C=(*cx*,*cy*,*cz*)
	# is the center point of the segment AB.  The crummy ASCII art below should help
	# to clarify the rather simple geometry:
	#
	#
	#               R
	#              /|\
	#             / | \
	#            /  |  \
	#           /   |   \
	#          /    |--+ \
	#         /     |  |  \
	#        A------C------B
	#
	# Please note that A, B an C are shown in the diagram above as being parallel to
	# the X axis.  A and B can actually be oriented in any direction.
	
	# Compute the location of C=(*cx*,*cy*,*cz*), which is exactly in the center
	# between points A and B:
	cx = (ax + bx) / 2
	cy = (ay + by) / 2
	cz = az
	if debug:
            radius3 = (radius_x - cx).distance(radius_y - cy)
	    print("C=({0:i},{1:i},{2:i}) radius_3={3:i}".format(cx, cy, cz, radius3))

    	# Now compute |AC| which is the length the the AC segment:
	ac_length = (cx - ax).distance(cy - ay)
	if debug:
	    print("ac_length={0:i}".format(ac_length))

	# We know that |AR| is equal to *radius*.  We need to compute |CR|.  This
	# is basically just the Pythagorean theorem:
	#
	#        |AC|^2 + |CR|^2 =         |AR|^2			(1) 
	#        |AC|^2 + |CR|^2 =       radius^2			(1) 
	#                 |CR|^2 =       radius^2 - |AC|^2		(2) 
	#                 |CR|   = sqrt( radius^2 - |AC|^2 )		(3)
	#                 |CR|   = sqrt( radius^2 - |AC|^2 )		(4)
	#
	# Now we do the computation of |CR| using equation (4) above :
	radius_mm = radius.millimeters()
	radius_squared_mm = radius_mm * radius_mm
	ac_mm = ac_length.millimeters()
	ac_squared_mm = ac_mm * ac_mm
	if radius_squared_mm < ac_squared_mm:
	    return
	assert radius_squared_mm >= ac_squared_mm
	cr_mm = math.sqrt(radius_squared_mm - ac_squared_mm)
	cr_length = L(mm=cr_mm)

	# Now we need to compute the location of R=(*rx*,*ry*,*rz*).  This done by
	# first computing the angle from A to C relative to the X axis.  Once we have
	# that angle we add/subtract 90 degrees to get the angle from C to R relative
	# to the X axis.  Once we have that angle, we can compute the location of R
	# relative to C.

	# Step 1: Compute angle <ACX:
	ac_dx = cx - ax
	ac_dy = cy - ay
	acx_angle = ac_dy.arc_tangent2(ac_dx)	# Yes, dy come before dx:
	if debug:
	    print("ac_dx={0:i} ac_dy={1:i} acx_angle={2:d}".format(ac_dx, ac_dy, acx_angle))

	# Step 2: Compute angle <RCX, which is +/- 90 degrees from <ACX:
	rcx_angle = acx_angle
	if clockwise:
	    rcx_angle += degrees90
	else:
	    rcx_angle -= degrees90
	if debug:
	    print("acx_angle={0:d} rcx_angle={1:d}".format(acx_angle, rcx_angle))

	# Step 3: Normalize rcx_angle to between -180 degrees and +180 degrees:
	if rcx_angle > degrees180:
	    rcx_angle -= degrees360
	elif rcx_angle < -degrees180:
	    rcx_angle += degrees360
	if debug:
	    print("normalized rcx_angle={0:d}".format(rcx_angle))

	# Step 4: Compute R=(*rx*,*ry*,*rz*) using trigonmetry:
	rx = cx + cr_length.cosine(rcx_angle)
	ry = cy + cr_length.sine(rcx_angle)
	rz = az
	if debug:
	    print("computed rx={0:i} ry={1:i} passed in rx={2:i} ry={3:i}".
	      format(rx, ry, radius_x, radius_y))

	# Now we compute the angle <ARX which is a the angle from R to A relative to the X axis:
	arx_angle = (ay - ry).arc_tangent2(ax - rx)

	# Similarly, we compute the <BRX which is the angle from R to B relative to the Y axis:
	brx_angle = (by - ry).arc_tangent2(bx - rx)

	if debug:
	    print("arx_angle={0:d} brx_angle={1:d}".format(arx_angle, brx_angle))

	# Now we can compute angle <ARB which is the angle swept from A to B from R.
	# Note that <ARB can be positive or negative depending up the direction of travel:
	arb_angle = brx_angle - arx_angle
	if debug:
	    print("arb_angle={0:d}".format(arb_angle))

	# Normalize <ARB to be consistent with *clockwise*:
	if clockwise and arb_angle < degrees0:
	    # Clockwise requires *arb_angle* to be positive:
	    arb_angle += degrees360
	elif not clockwise and arb_angle > degrees0:
	    # Counter-clockwise requires *arb_angle* to be negative:
	    arb_angle -= degrees360
	if debug:
	    print("normalize arb_angle={0:d}".format(arb_angle))

	# We want to draw the arc as a sequence of line segments from A to B, where each segment
	# covers about 15 degrees of the arc.  So first we figure compute the *steps* count:
	steps = 1 + abs(int(arb_angle / degrees15))
	assert steps >= 1

	# Now we compute the *step_angle* which is the arc angle coverted by each segment:
	step_angle = arb_angle / float(steps)

	# Lay down A=(*ax*,*ay*,*az*):
	vrml_point_indices = self._vrml_point_indices
	a_index = self._vrml_point(ax, ay, az)
	vrml_point_indices.append(a_index)

	# Lay down the intermediate points along the arc:
	for step in range(1, steps - 1):
	    segment_angle = arx_angle + step_angle * step
	    segment_x = rx + radius.cosine(segment_angle)
	    segment_y = ry +   radius.sine(segment_angle)
	    segment_z = az
	    segment_index = self._vrml_point(segment_x, segment_y, segment_z)
	    vrml_point_indices.append(segment_index)

	# Lay down B=(*bx*,*by*,*bz):
	b_index = self._vrml_point(bx, by, bz)
	vrml_point_indices.append(b_index)

	# Terminate the polyline:
	vrml_point_indices.append(-1)

	# Set the color for the polyline:
	self._vrml_color_indices.append(1)

    def _vrml_point(self, x, y, z):
	""" *Code*: Return the VRML point index for point (*x*, *y*, *z*) using the
	    *Code* object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(z, L)

	point = (x.millimeters(), y.millimeters(), z.millimeters())
	vrml_points_table = self._vrml_points_table
	if point in vrml_points_table:
	    index = vrml_points_table[index]
	else:
	    vrml_points = self._vrml_points
	    index = len(vrml_points)
	    vrml_points.append(point)
	return index

    def _vrml_reset(self):
	""" *Code*: Reset the VRML sub-system of the *Code* object (i.e. *self*). """

	self._vrml_color_indices = []
	self._vrml_motion = -1
	self._vrml_point_indices = []
	self._vrml_points = []
	self._vrml_points_table = []
	self._vrml_start_r0 = self._r0
	self._vrml_start_x = self._x
	self._vrml_start_y = self._y
	self._vrml_start_z = self._z
	self._vrml_stream = None

    def _x_value(self):
	""" *Code*: Return the current value of X offset by the vice X for the *Code* object
	    (i.e. *self*).
	"""

	return self._x + self._vice_x

    def _xy_cw_feed(self, f, r, s, x, y, rx=None, ry=None):
	""" *Code*: Feed to location (*x*, *y*) with a radius *r* clockwise  circle with a
	    feedrate of *f* and spindle speed of *s* using the *Code* object (i.e. *self*):
	"""
    
	# Verify routine arguments:
	assert isinstance(f, Speed)
	assert isinstance(r, L)
	assert isinstance(s, Hertz)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(rx, L) or rx == None
	assert isinstance(ry, L) or ry == None
    
	x1 = self._x_value()
	y1 = self._y_value()
	z1 = self._z
	x2 = x
	y2 = y
	self._vrml_arc_draw(x1, y1, x2, y2, r, z1, False, rx, ry)

	x_value = self._x_value()
	y_value = self._y_value()
	if x_value != x or y_value != y:
	    # Do the laser code first:
	    #if self._is_laser():
	    #	self._dxf_arc_append(True, x, y, r)

	    # Now do the RS274 self and get F, R0, S, X, and Y updated:
	    self._z_safe_retract_actual()
	    self._r0 = L(inch=-100.0)	# Forget R
	    self._command_begin()
	    self._mode_motion(2)
	    self._speed("F", f)
	    self._length("R0", r)
	    self._hertz("S", s)
	    self._length("X", x)
	    self._length("Y", y)
	    self._command_end()
    
    def _xy_ccw_feed(self, f, r, s, x, y, rx=None, ry=None):
	""" *Code*: Feed to location (*x*, *y*) as a radius *r* counter clockwise  circle with a
	    feedrate of *f* and spindle speed of *s* using the *Code* object (i.e. *self*):
	"""
    
	# Verify routine arguments:
	assert isinstance(f, Speed)
	assert isinstance(r, L)
	assert isinstance(s, Hertz)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(rx, L) or rx == None
	assert isinstance(ry, L) or ry == None
    
	x1 = self._x_value()
	y1 = self._y_value()
	z1 = self._z
	x2 = x
	y2 = y
	self._vrml_arc_draw(x1, y1, x2, y2, r, z1, True, rx, ry)

	if x1 != x or y1 != y:
	    # Do the laser code first:
	    #if self._is_laser():
	    #	self._dxf_arc_append(True, x, y, r)
    
	    # Now do the RS274 self and get F, R0, S, X, and Y updated:
	    self._z_safe_retract_actual()
	    self._r0 = L(inch=-100.00)	# Forget R
	    self._command_begin()
	    self._mode_motion(3)
	    self._speed("F", f)
	    self._length("R0", r)
	    self._hertz("S", s)
	    self._length("X", x)
	    self._length("Y", y)
	    self._command_end()

    def _xy_feed(self, f, s, x, y):
	""" *Code*: Feed to location (*x*, *y*) with a feedrate of *f*
	    and spindle speed of *s* using the *Code* object (i.e. *self*).
	"""
    
	# Verify argment types:
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)
	assert isinstance(x, L)
	assert isinstance(y, L)

	x_before = self._x_value()
	y_before = self._y_value()
	if x_before != x or y_before != y:
	    # Do any laser code first:
	    #if self._is_laser():
	    #	self._dxf_entity_start("LINE")
	    #	self._dxf_xy_append(0, x_before, y_before, "xy_feed before")
	    #	self._dxf_xy_append(1, x, y, "xy_feed after")
	    #	self._dxf_dxf_entity_stop()
    
	    # Now do the RS274 code and get the F, S, X, and Y values updated:
	    self._z_safe_retract_actual()
	    self._command_begin()
	    self._mode_motion(1)
	    self._speed("F", f)
	    self._hertz("S", s)
	    self._length("X", x)
	    self._length("Y", y)
	    self._command_end()

    def _xy_rapid(self, x, y):
	""" *Code*: Perform a rapid move to (X, Y) using the *Code* object (i.e. *self*). """

	# Verify argument types:
	assert isinstance(x, L)
	assert isinstance(y, L)

	# Only perform the rapid if we are not already there:
	if self._x_value() != x or self._y_value() != y:
	    # Make sure we are at a safe Z height prior to performing the rapid:
	    self._z_safe_retract_actual()

	    # Output a X/Y rapid command:
	    self._command_begin()
	    self._mode_motion(0)
	    self._length("X", x)
	    self._length("Y", y)
	    self._command_end()

    def _y_value(self):
	""" *Code*: Return the current value of Y offset by the vice Y for the *Code* object
	    (i.e. *self*).
	"""

	return self._y + self._vice_y

    def _z_feed(self, f, s, z, from_routine):
	""" *Code*: Feed to an altitude of *z* using the *Code* object (i.e. *self*)
	    at a feed of *f* and a speed *s*.
	"""

	# Verify argument types:
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)
	assert isinstance(z, L)
	assert isinstance(from_routine, str)

	# This routine will feed to an altitude of {z} using {code}.

	# If *z_safe_pending* is set, but we are doing another Z feed before
	# it is cleared, we must not need to do a Z safe move operation after all.
	# Hence, we can just clear *z_safe_pending*:
	self._z_safe_pending = False

	z_safe = self._z_safe
	z_rapid = self._z_rapid
	z_current = self._z

	# Move up if we are too low:
	while z_current < z:
	    # We are lower than {z} and need to move up:
	    if z_current < z_rapid:
		# We need to feed up:
		z_target = z
		if z > z_rapid:
		    z_target = z_rapid

		# We do a G1 feed to get to {z_target}:
		self._command_begin()
		self._mode_motion(1)
		self._speed("F", f)
		self._hertz("S", s)
		self._length("Z", z_target)
		self._command_end()
		z_current = z_target
	    else:
		# We are at or above {z_rapid}, so we can G0 rapid to a higher Z:
		self._command_begin()
		self._mode_motion(0)
		self._length("Z", z)
		self._command_end()
		z_current = z

	# Move down if we are too high:
 	while z_current > z:
	    # We are higher than {z} and need to move down:
	    if z_current > z_rapid:
		# We are above {z_rapid}, so we can G0 rapid lower:
		z_target = z
		if z < z_rapid:
		    z_target = z_rapid
	
		# We are at
		self._command_begin()
		self._mode_motion(0)
		self._length("Z", z_target)
		self._command_end()
		z_current = z_target
	    else:
		# We are are at or below {z_rapid} and need to G1 feed down:
		self._command_begin()
		self._mode_motion(1)
		self._speed("F", f)
		self._hertz("S", s)
		self._length("Z", z)
		self._command_end()
		z_current = z

    def _z_rapid_set(self, z_rapid):
	""" *Code*: Set the Z rapid field of the *Code* object (i.e. *self*) to *z_rapid*. """

	# Verify argument types:
	assert isinstance(z_rapid, L)

        # Set the Z field of the *Code* object (i.e. *self*) to *z*:
	self._z_rapid = z_rapid
    

    def _z_safe_assert(self, text, comment):
	""" *Code*: Fail if the Code object (i.e. *self*)  is not at "Z safe".
	    If not, the failing message will contain {text}.
	"""

	# Verify argument_types:
	assert isinstance(text, str)
	assert isinstance(comment, str)


	assert self._z_safe_pending or self._z == self._z_safe, \
	  "Z safe failure: {0} ({1}) [z={2}]".format(text, comment, code._z)

    def _z_safe_set(self, z_safe):
	""" *Code*: Set the Z safe field of the *Code* object (i.e. *self*). """

	# Verify argument types:
	assert isinstance(z_safe, L)

        # Set the Z safe field of the *Code* object (i.e. *self*) to *z_safe*:
	self._z_safe = z_safe

    def _z_safe_retract(self, f, s):
	""" *Code*: Ensure that the tool is at the "Z safe" altitude for next command RS274
	    command in the *Code* object (i.e. *self) and get there using a feedrate of *f*
	    and a speed of *s* if needed.
	"""

	# Verify argument types:
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)

	# Mark that a "Z safe" operation is pending if needed:
	if self._z != self._z_safe:
	    self._z_safe_pending = True
	    self._z_safe_s = s
	    self._z_safe_f = f

    def _z_safe_retract_actual(self):
	""" *Code*  Ensure that the tool is at the "Z safe" altitude in the *Code* object
	    (i.e. *self*) and get there using the internal F and S fields.
	"""

	# If we have a pending "z_safe" scheduled, this code will make it happen:
	if self._z_safe_pending:
            self._z_safe_pending = False
	    self._z_feed(self._z_safe_f, self._z_safe_s, self._z_safe, "z_safe_retract_actual")

	#assert self._z == self._z_safe, "Did not get to Z safe"

    def _z_set(self, z):
	""" *Code*: Set the Z field of the *Code* object (i.e. *self*). """

	# Verify argument types:
	assert isinstance(z, L)

        # Set the Z field of the *Code* object (i.e. *self*) to *z*:
	self._z = z

## @brief *Place* specifies where another *Part* is to be placed.
#
# A *Place* specifies the placement of a *Part* in an assembly.
# There are two independent operations for a *Place*:
#
# * Rotation.  The *Part* can be rotated around an arbitrary center point
#   The rotation is specified as a rotation axis point and the angle to
#   rotate by.
# * Translate.  Independent of rotation the *Part* can moved in X/Y/Z.

class Place:
    """ A {Place} represents the placement of a part relative to the
	origin of an assembly. """

    def __init__(self, part = None, name = None, center = None,
      axis = None, rotate = None, translate = None, trace = -100000):
	""" Place: Initialize {self} to contain {home_part}, {placed_part},
	    {center_point}, {axis_point}, {rotate_angle}, and
	    {translate_point}. """

	# Check argument types:
	none_type = type(None)
	assert type(center) == none_type or isinstance(center, P)
	assert type(axis) == none_type or isinstance(axis, P)
	assert type(rotate) == none_type or isinstance(rotate, Angle)
	assert type(translate) == none_type or isinstance(translate, P)

	# Deal with default argument values:
	name = "no_name"
	if isinstance(part, Part):
	    name = part._name
	if type(center) == none_type:
	    center = P()
	if type(axis) == none_type:
	    axis = P(z = L(1.0))
	if type(rotate) == none_type:
	    rotate = Angle()
	if type(translate) == none_type:
	    translate = P()

	part_name = "<none>"
	if type(part) != none_type:
	    part_name = part._name

	# Check argument types:
	#assert type(part) == none_type or isinstance(part, Part)
	#assert isinstance(name, str)

	if trace >= 0:
	    print(("{0}=>Place.__init__(part={1}, name={2}, center={3}, " + \
	      "axis={4}, rotate={5}, translate={6})").
	      format(trace * ' ', part_name,
		name, center, axis, rotate, translate))

	# Normalize *axis*:
	axis = axis.normalize()

	# Rotate 
	rotate_matrix = Matrix.rotate_create(axis.x, axis.y, axis.z, rotate)

	# Create {center_matrix} and {center_reverse_matrix}:
	center_matrix = Matrix.translate_create(-center.x, -center.y, -center.z)
	center_reverse_matrix = \
	  Matrix.translate_create(center.x, center.y, center.z)

	# Create {translate_matrix}:
	translate_matrix = \
	  Matrix.translate_create(translate.x, translate.y, translate.z)
	if trace >= 0:
	    print("{0}translate_matrix={1}".
	      format(trace * ' ', translate_matrix))

	#print "Place(): cm=\n{0}\nrm=\n{1}\nrcm=\n{2}\ntm=\n{3}\n". \
	#  format(center_matrix.mat, rotate_matrix.mat, \
	#  center_reverse_matrix.mat, translate_matrix.mat)

	# Compute {forward_matrix} and {reverse_matrix}:
	forward_matrix = center_matrix * \
	  rotate_matrix * center_reverse_matrix * translate_matrix
	reverse_matrix = forward_matrix.inverse()
	if trace >= 0:
	    print("{0}Place(): forward=\n{1}\nreverse=\n{2}\n". \
	      format(trace * ' ', forward_matrix.mat, reverse_matrix.mat))

	# Load up *self*:
	self._axis = axis
	self._center = center
	self._forward_matrix = forward_matrix
	self._part = part
	self._name = name
	self._reverse_matrix = reverse_matrix
	self._rotate = rotate
	self._translate = translate

	if trace >= 0:
	    print("{0}<=Trace.__init__()".format(trace * ' '))

    def __eq__(self, place):
	""" Place: Return {True} if {self} is equal to {place}. """

	assert isinstance(place, Place)

	name_equal = self._name == place._name
	axis_equal = self._axis == place._axis
	center_equal = self._center == place._center
	part_equal = self._part == place._part
	rotate_equal = self._rotate == place._rotate
	translate_equal = self._translate == place._translate

	result = name_equal and axis_equal and center_equal and \
	  part_equal and rotate_equal and translate_equal

	#if changed:
	#    print "home_part_changed=", home_part_changed
	#    print "placed_part_changed=", placed_part_changed
	#    print "location_changed=", location_part_changed
	#    print "axis_point_changed=", axis_point_changed
	#    print "rotate_angle_changed=", rotate_angle_changed

	#print "Place.__eq__() =>", result
	return result

    def __format__(self, format):
	if format != "":
	    format = ":" + format 
	format_string = \
	  "[Place: name='{0}, part='{1}', center={2" + format + \
	  "}, axis={3"                               + format + \
	  "}, rotate={4"                             + format + \
	  "}, translate={5"                          + format + \
	  "}]"

	return format_string.format(self._name, self._part._name,
	  self._center, self._axis, self._rotate, self._translate)

    def __ne__(self, place):
	""" Place: Return {True} if {self} is not equal to {place}. """

	return not (self == place)

class Screw:
    """ Screw """

    def __init__(self, part, name, thread, direction):
	""" Screw internal: Initialize {self} with {thread} ..."""

	assert direction == "NS" or direction == "EW" or direction == "TB", \
	  "Screw direction '{0}' must be either 'NS', 'EW', or 'TB'"

	#print "Screw(part={0}, name={1}, thread={2}, dir={3})". \
	#  format(part.name, name, thread, direction)

	# Load up {self}:
	self.anchor_point_mapped = None
	self.anchor_point_original_part = None
	self.anchor_screw_level = None
	self.direction = direction
	self.name = name
	self.part = part
	self.screw_levels = {}
	self.thread = thread

    def anchor_set(self, anchor_point):
	""" Screw dimensions: Set the anchor point for {self} to
	    {anchor_point}. """

	trace = False
	#trace = self.name.find("skin_west_bottom") >= 0

	if trace:
	    print "=>Screw.anchor_set({0}, {1})".format(self.name, anchor_point)
	    anchor_part = anchor_point.part
	    bsw = anchor_part.point("$BSW")
	    tne = anchor_part.point("$TNE")
	    print "bsw={0}".format(bsw)
	    print "tne={0}".format(tne)

	part = self.part
	if part.dimensions_refine_mode():
	    # Grab some values from {self} and {anchor_point}:
            anchor_point_original = anchor_point
	    anchor_point_original_part = anchor_point_original.part
	    anchor_screw_level = self.anchor_screw_level

	    # If {anchor_screw_level} is not defined, search for it:
	    if anchor_screw_level == None:
		# Now find the {Screw_Level}} that matches
		# {anchor_point_original_part}:
		for screw_level in self.screw_levels.values():
		    assert screw_level.screw == self
		    if screw_level.part == anchor_point_original_part:
			anchor_screw_level = screw_level
			break

		# Make sure we found an anchor level that matched:
		assert anchor_screw_level != None, \
		  "Screw '{0}' is not attached to part '{1}'". \
		  format(self.name, anchor_point_original_part.name)
	    assert anchor_screw_level.screw == self

            # Make sure that we are not attempting to define 2 anchor points:
	    anchor_screw_level_part = anchor_screw_level.part
	    assert anchor_screw_level_part == anchor_point_original_part, \
	      "Screw anchor should be in part {0} instead of part {1}". \
	      format(anchor_screw_level_part.name, \
	      anchor_point_original_part.name)

            # Grab some values from {anchor_level}:
	    anchor_forward_matrix = anchor_screw_level.forward_matrix

	    # Perform the final subtraction and refernce to {part2}:
	    anchor_point_mapped = \
	      anchor_forward_matrix.point_multiply(anchor_point, part)
            if trace:
		print "anchor_point_mapped={0}".format(anchor_point_mapped)

	    # Get {screw} anchor point set:
	    if self.anchor_screw_level == None:
		self.anchor_screw_level = anchor_screw_level
		self.anchor_point_original_part = anchor_point_original_part
		self.anchor_point_mapped = anchor_point_mapped
		part.dimensions_changed("Screw.anchor_set 1")
	    else:
		if self.anchor_point_mapped != anchor_point_mapped:
                    self.anchor_point_original_part = \
		      anchor_point_original_part
		    self.anchor_point_mapped = anchor_point_mapped
		    part.dimensions_changed("Screw.anchor_set 2")

	if trace:
	    print "<=Screw.anchor_set({0}, {1})".format(self.name, anchor_point)
	    print ""

    def part_matrices_find(self, screw_path):
	""" Screw internal: Return a 2-tuple containing the {Part} and
	    forward matrix associated with {screw_path}. """

	trace = False
	#trace = self.name.find("grip_attach_bottom_south") >= 0

	if trace:
	    print "  =>Screw.part_matrix_find({0}, {1})". \
	      format(self.name, screw_path)
	forward_matrix = Matrix.identity_create()
	reverse_matrix = Matrix.identity_create()

	screw_part = self.part
	target_part = screw_part
	place_names = screw_path.split('/')

	if trace:
	    print "  target_part={0} place_names={1}". \
	      format(target_part.name, place_names)

	for place_name in place_names:
            places = target_part.places
            if not place_name in places:
		print "places=", places
		message = "Invalid place name '{0}'" + \
		  " in screw path '{1}' starting from part '{2}'"
		assert False, \
		  message.format(place_name, screw_path, screw_part.name)

            # Move forward by {place}:
            place = target_part.places[place_name]
	    if trace:
		print "  place={0} place_matrix=\n{1}". \
		  format(place_name, place.forward_matrix)
            forward_matrix = forward_matrix * place.forward_matrix
	    reverse_matrix = reverse_matrix * place.reverse_matrix
	    target_part = place.placed_part

	result = (target_part, forward_matrix, reverse_matrix)
	if trace:
	    print "  <=Screw.part_matrix_find({0}, {1})=>({2}, vvv)\n{3}". \
	      format(self.name, screw_path, target_part.name, forward_matrix)
	return result

    def attach(self, screw_path, surface, flags):
	""" Screw dimensions: Set the ... """

	#print "=>Screw.attach({0}, {1}, {2}, {3})". \
	#  format(self.name, screw_path, surface, flags)

	# Grab some values from {self}:
	screw_part = self.part
	screw_levels = self.screw_levels

	# Find {target_part} and {forward_matrix} associated with {screw_path}:
	part_matrices = self.part_matrices_find(screw_path)
	target_part = part_matrices[0]
	forward_matrix = part_matrices[1]
	reverse_matrix = part_matrices[2]

	if screw_part.dimensions_define_mode():
            # In define mode, we create the new screw level:
	    assert not screw_path in screw_levels, \
	      "Duplicate screw level path '{0}'".format(screw_path)
	    screw_level = Screw_Level(self, screw_path, surface, flags, \
	      target_part, forward_matrix, reverse_matrix)
	    screw_levels[screw_path] = screw_level

	    # Make sure that we store {screw_level} in {target_part}.
	    target_part_screw_levels = target_part.screw_levels
            target_part_screw_levels[surface][self.name] = screw_level
	elif screw_part.dimensions_refine_mode():
	    # Look up the {screw_path}:
	    assert screw_path in screw_levels, \
	      "Undefined screw level path '{0}' in part {1}". \
	      format(screw_path, screw_part.name)
	    screw_level = screw_levels[screw_path]

	    # Update {forward_matrix}:
	    if not screw_level.forward_matrix == forward_matrix:
		screw_level.forward_matrix = forward_matrix
		screw_level.reverse_matrix = reverse_matrix
		screw_part.dimensions_changed("Screw.attach")

	    # Check that anchor actually occured:
            if flags.find('A') >= 0:
		assert self.anchor_screw_level != None, \
		  "Screw '{0}' has no anchor point set".format(self.name)
		#print "anchor_point_map={0}".format(self.anchor_point_mapped)
		assert self.anchor_point_original_part == target_part, \
		  "Screw '{0}' has anchor set to part '{1}' instead of '{2}'".\
		  format(self.name, self.anchor_point_original_part.name, \
		  target_part.name)
	else:
            assert screw_part.construct_mode(), \
	      "Setting screw '{0}' level for '{1}' outside dimensions mode". \
	      format(self.name, screw_path)

	#print "<=Screw.attach({0}, {1}, {2}, {3})". \
	#  format(self.name, screw_path, surface, flags)

    def depth_set(self, target_part, depth):
	""" Screw dimensions: Set the screw hole depth for {self} to {depth}
	    for {part}. """

	root_part = self.part
	if root_part.dimensions_refine_mode():
	    screw_level_match = None
	    for screw_level in self.screw_levels.values():
		assert screw_level.screw == self
		if screw_level.part == target_part:
		    screw_level_match = screw_level
                    break

            # Make sure we found an anchor level that matched:
            assert screw_level_match != None, \
	      "Screw '{0}' is not attached to part '{1}'". \
		  format(self.name, target_part.name)

	    # Set the {depth}:
	    screw_level.depth_set(depth)

class Screw_Level:
    """ Screw Level """

    def __init__(self, screw, screw_path, surface, flags, part, \
      forward_matrix, reverse_matrix):
	""" Screw_Level internal: Initialize {self} with {screw},
	    {screw_path}, {surface}, {flags}, {part}, {forward_matrix}
	    and {reverse_matrix}. """

	#trace = True
	#if trace:
	#    print "Screw_Level({0} {1} {2} {3} {4}): fm=\n{5}\nrm=\n{6}\n". \
	#      format(screw.name, screw_path, surface, flags, part.name, \
	#      forward_matrix.mat, reverse_matrix.mat)

	self.depth = L.inch(0)
	self.done = False
	self.flags = flags
	self.forward_matrix = forward_matrix
	self.part = part
	self.reverse_matrix = reverse_matrix
	self.screw = screw
	self.screw_path = screw_path
	self.surface = surface

    def __str__(self):
	""" Screw_Level internal: Return {self} as a {String}. """

	return "String_Level(scr={0} path={1} sur={2} flgs={3} part={4} *)". \
	  format(self.screw.name, self.screw_path, self.surface, \
	  self.flags, self.part.name)

    def depth_set(self, depth):
	""" Screw dimensions: Set screw hole depth for {self} to {depth}. """

	if self.depth != depth:
            self.depth = depth
	    self.screw.part.dimensions_changed("Screw_Level.depth_set")

# screw = part.screw_create("screw1", part.screw_new("thread", TB))
# screw.attach("path", "flags", "surface")
# screw.anchor_set(point)
# screw.depth_set(depth)

class Shop:

    def __init__(self, name):
	""" *Shop*: Initialize a *Shop* object. """

	# Verify argument types:
	assert isinstance(name, str)

	zero = L()

	tools = []

	self._assemblies = []		# Viewable assemblies
	self._base_name = ""		# Base name for current operations
	#dxf_base_names = [] 		# List of DXF base names
	self._blocks_uid = 0		# UID counter for *Code_Block* objects
	# cache Cache			# Off/Nef3 file cache
	self._changed = False		# Marker used by {update@Length}
	self._cnc_generate = False	# {true} => do cnc code generation
	#dxf_table = {}			# Table of DXF base names
	self._extra1 = P()		# Temporary extra bounding box point
	self._extra2 = P()		# Temporary extra bounding box point
	self._name = name		# Shop name
	self._machines = []		# Machines in shop
	self._matrix = Matrix()		# Temporary matrix
	self._parts = []			# Parts
	self._program_base = 10		# Program base number
	self._code = Code()		# Code genertion object
	self._solids_generate = False	# {true} => do solid modeling
	self._surface_speeds_table = {}	# [Part_Material, Tool_Material]
	self._tools = tools		# Tools in shop
	self._vice = Vice(zero, zero)	# Vice to use
	#tess GLU_Tess			# Tessellator
	#tess_mode Unsigned		# Triangulation mode
	#tess_points Array[Point]	# Collection of points from tessellation
	#tess_polygons Array[Simple_Polygon] # Polygons
	#aught80 Thread			# 0-80 thread
	#two56 Thread			# 2-56 thread
	#four40 Thread			# 4-40 thread
	#six32 Thread			# 6-32 thread
	#eight32 Thread			# 8-32 thread
	#ten24 Thread			# 10-24 thread
	#ten32 Thread			# 10-32 thread

	# Copied from:
	#
	#	# http://Whitney-Tool.Com/calculatorSpeedFeed.html (old)
	#	https://www.whitneytool.com/SpeedAndFeedCalculator.aspx
	# 
	# HSS = High Speed Steel
	# CTS = Cobalt Tool Steel
	# UC = Uncoated Carbide
	# CC = Coated Carbide
	# MCHT = Medium Carbon Heat Treated
	#
	# All numbers are in Feet/Min.
	# 
	# Reference: Machining Data Handbook; Machinability Data Center
	#
	# Material			HSS	CTS	UC	 CC
	# 
	# Non-Ferrous Material
	# Aluminum Alloys		600+	-	1200+	-
	# Magnesium Alloys		600+	- 	1000+	-
	# Brass				300+	- 	800	650+
	# Bronze			80-100	-	250-300	-
	# 
	# Titanium (Double Starging Feed Rates)
	# Commercially Pure		115-140	-	275-325	-
	# Alpha & Alapha-Beta Alloys	-	30-50	200-225	-
	# 
	# Ferrous Material
	# Steels
	# Free Machining Carbon Steel	130-180	-	450-500	750-900
	# Low Carbon Steel		120-170	-	400-450	600-650
	# Medium Carbon Steel		100-120	-	375-425	550-600
	# Alloy Steel			100-120		375-425	550-600
	# Alloy & MCHT (Rc 26-32)	75-100	-	250-300	450-500
	# Alloy & MCHT (Rc 36-40)	-	50-60	180-200	225-275
	# Alloy & MCHT (Rc 40-48)	-	40-50	150-180	220-250
	# Alloy & MCHT (Rc 48+)		-	20-30	100-120	
	# Tool Steel (Wrought)		40-60	-	180-200	350
	# 
	# Stainless Steels
	# Free Machining		80-110	-	100-140	140+
	# Stainless (300 Series)	50-70	-	80-100	100+
	# 17-4PH Annealled		50-80	-	150-190	190+
	# 17-4PH 200,000 PSI		30-50	-	100-140	140+
	# 
	# High Temp. Alloys
	# Hasteloy X, Inconel		15-20	-	45-55	- 
	# Inconel X			-	20-25	-	-	 
	# Monel Nickel Alloy		-	20-25	-	-
	# 
	# Cast Iron
	# Malleable Iron		100-140	-	400-450	540-700
	# Gray Cast Iron		65-110	-	220-300	340-450
	# Ductile Iron			80-125	- 	300-350	460-550
	#
	# The following table came from:
	#
	#  http://www.projectsinmetal.com/
	#     lathe-cutting-speed-chart-in-feet-per-minute/
	#
	#		Carbide			HSS
	# Material	Rough	Finish		Rough	Finish	
	# =====================================================
	# Stainless	100-210	170-330		26-33	39-46
	# Hard Steel	165-330	200-460		39-52	50-70
	# Carbon Steel	260-460	400-820		54-82	66-100
	# Mild Steel	260-500	330-660		60-100	70-130
	# Bronze	130-420	330-600		33-130	66-165
	# Brass		500	500-1150	300-200	200-260
	# Aluminum	2100	3300		1000	1500

	# Start loading some surface speeds into *self*:
	hss = Tool.MATERIAL_HIGH_SPEED_STEEL
	plastic = Material("Plastic", "")
	aluminum = Material("Aluminum", "")
	fpm600 = Speed(ft_per_min=600)
	fpm1200 = Speed(ft_per_min=1200)
	print("type(hss)=", type(hss))
	self._surface_speeds_insert(aluminum, hss, fpm600, fpm1200)
	self._surface_speeds_insert(aluminum, hss, fpm600, fpm1200)
	self._surface_speeds_insert(plastic, hss, fpm600, fpm1200)
	self._surface_speeds_insert(plastic, hss, fpm600, fpm1200)

	degrees45 = Angle(deg=45.0)
	degrees90 = Angle(deg=90.0)
	degrees118 = Angle(deg=118.0)

	# Some common inch sizes:
	in1_16 = L(inch="1/16")
	in3_32 = L(inch="3/32")
	in1_8 =  L(inch="1/8")
	in3_16 = L(inch="3/16")
	in1_4 =  L(inch="1/4")
	in3_8 =  L(inch="3/8")
	in1_2 =  L(inch="1/2")
	in5_8 =  L(inch="5/8")
	in3_4 =  L(inch="3/4")

	stub = Tool.DRILL_STYLE_STUB
	laser = True

	dowel_pin = self._dowel_pin_append("3/8 Dowel Pin",
	  1, hss, in3_8, L(inch=.900), in3_16)
	mill_drill_3_8 = self._mill_drill_append("3/8 Mill Drill",
	  2, hss, in3_8, 2, L(inch=.900), in3_16, degrees90)
	drill_36 = self._drill_append("#36 Drill",
	  3, hss, L(inch=0.1065), 2, L(inch=1.500), degrees118, stub)
	drill_27 = self._drill_append("#27 drill",
	  4, hss, L(inch=0.1440), 2, L(inch=1.750), degrees118, stub)
	end_mill_3_8 = self._end_mill_append("3/8 End Mill",
	  5, hss, in3_8, 2, in5_8, not laser)
	end_mill_1_4 = self._end_mill_append("1/4 End Mill",
	  6, hss, in1_4, 2, in1_2, not laser)
	double_angle = self._double_angle_append("3/4 Double Angle",
	  7, hss, in3_4, 10, L(inch=0.875), degrees90, in1_4, in1_4)
	dove_tail = self._dove_tail_append("3/8 Dove Tail",
	  8, hss, in3_8, 6, in1_4, in3_16, degrees45)
	end_mill_3_16 = self._end_mill_append("3/16 End Mill",
	  10, hss, in3_16, 2, in1_2, not laser)
	drill_25 = self._drill_append("#25 drill",
	  11, hss, L(inch=0.1495), 2, L(inch=2.000), degrees118, stub)
	drill_9 = self._drill_append("#9 drill",
	  12, hss, L(inch=0.1960), 2, L(inch=2.000), degrees118, stub)
	drill_43 = self._drill_append("#43 drill",
	  13, hss, L(inch=.0890), 2, L(inch=1.500), degrees118, stub)
	drill_32 = self._drill_append("#32 drill",
	  14, hss, L(inch=0.1160), 2, L(inch=1.500), degrees118, stub)
	drill_50 = self._drill_append("#50 drill",
	  15, hss, L(inch=0.0700), 2, L(inch=1.500), degrees118, stub)
	end_mill_3_8_long = self._end_mill_append("3/8 1\" End Mill",
	  16, hss, in3_8, 2, L(inch=1.000), not laser)
	#end_mill_3_4 = self._end_mill_append("3/4 End Mill",
	#  13, hss, in3_4, 2, in1_3_8)
	drill_30 = self._drill_append("#30 drill",
	  17, hss, L(inch=0.1285), 2, L(inch=1.750), degrees118, stub)
	drill_1_8 = self._drill_append("1/8 drill",
	  18, hss, in1_8, 2, L(inch=1.750), degrees118, stub)
	drill_1_16 = self._drill_append("1/16 drill",
	  19, hss, in1_16, 2, L(inch=1.750), degrees118, stub)
	drill_3_32 = self._drill_append("3/32 drill",
	  20, hss, in3_32, 2, L(inch=1.750), degrees118, stub)
	drill_42 = self._drill_append("#42 drill",
	  21, hss, L(inch=0.0935), 2, L(inch=1.750), degrees118, stub)

	# Laser "tools":
	laser_007 = self._end_mill_append("Laser_007",
	  100, hss, L(inch=0.007), 2, L(inch=0.750), laser)
	laser_000 = self._end_mill_append("Laser_000",
	  101, hss, L(), 2, L(inch=0.750), laser)

    def _blocks_uid_get(self):
	""" *Shop*: Return the next block UID (Unique IDentifier) from the *Shop* object
	    (i.e. *self*.)
	"""

	blocks_uid = self._blocks_uid + 1
	self._blocks_uid = blocks_uid
	return blocks_uid

    def _code_get(self):
	""" *Shop*: Return the *Code* object associated with the *Shop* object (i.e. *self*.) """

	return self._code

    def _cnc_generate_get(self):
	""" *Shop*: Return the CNC generate field associated with the *Shop* object
	    (i.e. *self*.)
	"""

	return self._cnc_generate

    def _cnc_generate_set(self, cnc_generate):
	""" *Shop*: Set the CNC generate field associated with the *Shop* object (i.e. *self*)
	    to *cnc_generate*.
	"""

	assert isinstance(cnc_generate, bool)
	self._cnc_generate = cnc_generate

    def _double_angle_append(self,
      name, number, material, diameter, flutes, maximum_z_depth, angle, inside_diameter, thickness):
	"""
	    *Shop*: Create and return a *Tool_Double_Angle* object that contains *name*, *number*,
	    *material*, *diameter*, *flutes*, *maximum_z_depth*, *angle*, *inside_diameter* and
	    *thickness*.  The returned *Tool_Double_Angle* object is also appended to the
	    *Shop* object (i.e. *self) tool list.
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter > zero
	assert isinstance(flutes, int) and flutes > 0
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero
	assert isinstance(angle, Angle) and angle > Angle()
	assert isinstance(inside_diameter, L) and inside_diameter > zero
	assert isinstance(thickness, L) and thickness > zero

	# Create the *tool_double_angle* and append it to the tools list in the *Shop* object
	# (i.e. *self*):
	tool_double_angle = Tool_Double_Angle(name,
	  number, material, diameter, flutes, maximum_z_depth, angle, inside_diameter, thickness)
	self._tool_append(tool_double_angle)
	return tool_double_angle

    def _dove_tail_append(self,
      name, number, material, diameter, flutes, maximum_z_depth, inside_diameter, angle):
	""" *Shop*: Create and return a *Tool_Dove_Tail* object that contains *name*, *number*,
	    *material*, *diameter*, *flutes*, *maximum_z_depth*, *inside_diameter*, and *angle*.
	    The returned *Tool_Dove_Tail* is also appended to the tool list in the *Shop* object
	    (i.e. *self*).
	"""

	# Verify argument types:
	zero = L()
	isinstance(name, str)
	isinstance(number, int) and number > 0
	isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	isinstance(diameter, L) and diameter > zero
	isinstance(flutes, int) and flutes > 0
	isinstance(maximum_z_depth, L) and diameter > zero
	isinstance(inside_diameter, L) and inside_diameter > zero
	isinstance(angle, Angle)

	# Create the *tool_dove_tail* and append it to the tools list of the *Shop* object
	# (i.e. *self*):
	tool_dove_tail = Tool_Dove_Tail(name,
	  number, material, diameter, flutes, maximum_z_depth, inside_diameter, angle)
	self._tool_append(tool_dove_tail)
	return tool_dove_tail

    def _dowel_pin_append(self, name, number, material, diameter, maximum_z_depth, tip_depth):
	""" *Shop*: Create and return a *Tool_Dowel_Pin* object that contains a dowel pin.
	    The returned *Tool* object contains *name*, *number*, *material*, *diameter*,
	    *maximum_z_depth*, and *tip_depth*.   The returned *Tool* object, is also
	    appended to the tool list of the *Shop* (i.e. *self*.)
	"""

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(tip_depth, L)

	# Create the dowel pin *tool_dowel_pin*, add it to the *Shop* object (i.e. *self) tool list,
	# and return the *tool*:
	tool_dowel_pin = Tool_Dowel_Pin(name,
	  number, material, diameter, maximum_z_depth, tip_depth)
	self._tool_append(tool_dowel_pin)
	return tool_dowel_pin

    def _drill_append(self,
      name, number, material, diameter, flutes_count, maximum_z_depth, point_angle, drill_kind):
	""" *Shop*: Create and return a *Tool_Dowel_Pin* object that contains a dowel pin.
	    The returned *Tool* object contains *name*, *number*, *material*, *diameter*,
	    *flutes_count*, *maximum_z_depth*, and *tip_depth*.   The returned *Tool* object,
	    is also appended to the tool list of the *Shop* (i.e. *self*.)
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter > zero
	assert isinstance(flutes_count, int) and flutes_count > 0
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero

	# Create the dowel pin *tool_dowel_pin*, add it to the *Shop* object (i.e. *self) tool list,
	# and return the *tool*:
	tool_dowel_pin = Tool_Drill(name,
	  number, material, diameter, flutes_count, maximum_z_depth, point_angle, drill_kind)
	self._tool_append(tool_dowel_pin)
	return tool_dowel_pin

    def _end_mill_append(self,
      name, number, material, diameter, flutes_count, maximum_z_depth, is_laser):
	""" *Shop*: Create and return a *Tool_Mill_Drill* object that contains
	    *name*, *number*, *material*, *diameter*, *flutes_count*, *maximum_z_depth* and
	    *is_laser*.  The returned *Tool_End_Mill* object is also appended to
	    the tool list of the *Shop* object (i.e. *self*.)
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter >= zero	# Zero is OK
	assert isinstance(flutes_count, int) and flutes_count > 0
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero
	assert isinstance(is_laser, bool)

	tool_end_mill = \
	  Tool_End_Mill(name, number, material, diameter, flutes_count, maximum_z_depth, is_laser)
	self._tool_append(tool_end_mill)
	return tool_end_mill

    def _matrix_get(self):
	""" *Shop*: Return the *Matrix* object associated with the *Shop* object (i.e. *self*.) """

	return self._matrix

    def _mill_drill_append(self,
      name, number, material, diameter, flutes, maximum_z_depth, tip_depth, point_angle):
	""" *Shop*: Create and return a *Tool_Mill_Drill* object that contains
	    *name*, *number*, *material*, *diameter*, *flutes*, *maximum_z_depth* and
	    *point_angle*.  The returned *Tool_Mill_Drill* object is also appended to
	    the tool list of the *Shop* object (i.e. *self*.)
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter > zero
	assert isinstance(flutes, int) and flutes > 0
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero
	assert isinstance(tip_depth, L) and tip_depth > zero
	assert isinstance(point_angle, Angle) and point_angle > Angle()

	# Create *tool_mill_drill*, append it to the *Shop* object (i.e. *self*)
	# tools list, and return it:
	tool_mill_drill = Tool_Mill_Drill(name,
	  number, material, diameter, flutes, maximum_z_depth, tip_depth, point_angle)
	self._tool_append(tool_mill_drill)
	return tool_mill_drill

    def _tools_get(self):
	""" *Shop*: Return the tools list object associated with the *Shop* object
	    (i.e. *self*.)
	"""

	return self._tools

    def _program_base_get(self):
	""" *Shop*: Return the program base number object associated with the *Shop* object
	    (i.e. *self*.)
	"""

	return self._program_base

    def _vice_get(self):
	""" *Shop*: Return the *Vice* object associated with the *Shop* object (i.e. *self*.) """

	return self._vice

    def _surface_speeds_insert(self, part_material, tool_material, low, high):
	""" *Shop*: Insert *speed_range* containing *low* and *high* into
	    the surface speeds table of the *Shop* object (i.e. *self)*
	     keyed to *part_material* and *tool_material*.
	"""

	# Verify argument types:
	assert isinstance(part_material, Material)
	assert isinstance(tool_material, int)
	assert isinstance(low, Speed)
	assert isinstance(high, Speed)

	# Enter *speed_range* into *self*:
	speed_range = Speed_Range(low, high)
	key = (part_material._generic_get(), tool_material)
	self._surface_speeds_table[key] = speed_range

	# Set *debug* to *True* to trace what is going on:
	debug = False
	debug = True
	if debug:
	    print("speeds_table[{0}] = {1}".format(key, speed_range))

    def _surface_speeds_lookup(self, part_material, tool_material):
	""" *Shop*: Lookup and return a *Speed_Range* object from the
	    *Shop* object (i.e. *self*) keyed on *part_material* and
	    *tool_material*; or return *None* otherwise.
	"""

	# Verify argument types:
	assert isinstance(part_material, Material)
	assert isinstance(tool_material, int)

	# Look up the *Speed_Range* object from *surfaces_speeds_table*:
	key = (part_material._generic_get(), tool_material)
	#print("Shop.surface_speeds_lookup():key={0}".format(key))
	surface_speeds_table = self._surface_speeds_table
	speed_range = None
	if key in surface_speeds_table:
	    speed_range = surface_speeds_table[key]
	return speed_range

    def _tool_append(self, new_tool):
	""" *Shop*: Append *new_tool* to the tools lins in the *Shop* object (i.e. *self*).
	    Duplicate tool numbers cause an assertion failure.
	"""

	# Verify argument types:
	assert isinstance(new_tool, Tool)

	# Search for duplicate *new_number* in *tools*:
	new_number = new_tool._number_get()
	tools = self._tools
	for tool in tools:
	    number = tool._number_get()
	    assert new_number != number, \
	     "Tool number {0} is conflicts between '{1}' and '{2}'".format(
	     number, new_tool._name_get(), tool._name_get())

	# Append *new_tool* to *tools*:
	tools.append(new_tool)

class Time:
    def __init__(self, sec=0.0, min=0.0):
	self._seconds = sec + min * 60.0 

class Tool:
    # Tool kind:
    KIND_NONE = 0
    KIND_DRILL = 1
    KIND_DOUBLE_ANGLE = 2
    KIND_DOVE_TAIL = 3
    KIND_DOWEL_PIN = 4
    KIND_END_MILL = 5
    KIND_LASER = 6
    KIND_MILL_DRILL = 7
    KIND_LAST = 8

    # Drill style:
    DRILL_STYLE_NONE = 0
    DRILL_STYLE_DEEP_HOLE = 1
    DRILL_STYLE_FAST_SPIRAL = 2
    DRILL_STYLE_LEFT_HAND = 3
    DRILL_STYLE_LONG = 4
    DRILL_STYLE_PARABOLIC = 5
    DRILL_STYLE_SPOTTING = 6
    DRILL_STYLE_SLOW_SPIRAL = 7
    DRILL_STYLE_STANDARD = 8
    DRILL_STYLE_STUB = 9
    DRILL_STYLE_STRAIGHT_FLUTE = 10
    DRILL_STYLE_LAST = 11

    # Define Tool material:
    MATERIAL_NONE = 0
    MATERIAL_ALLOY_STEEL = 1
    MATERIAL_COBLT_TOOL_STEEL = 2
    MATERIAL_CARBIDE_TIPPED_STEEL = 3
    MATERIAL_HIGH_SPEED_STEEL = 4
    MATERIAL_OTHER = 5
    MATERIAL_SOLID_CARBIDE = 6
    MATERIAL_TOOL_STEEL = 7
    MATERIAL_TUNGSTON_HIGH_SPEED_STEEL = 8
    MATERIAL_LAST = 9

    def __init__(self, name, number, kind, material, diameter, flutes_count, maximum_z_depth):
	""" *Tool*: Initialize a *Tool* object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(kind, int) and Tool.KIND_NONE < kind < Tool.KIND_LAST
	assert isinstance(material, int)
	assert isinstance(diameter, L)
	assert isinstance(flutes_count, int)
	assert isinstance(maximum_z_depth, L)

	# Load up *self*:
	self._name = name			# Tool name
	self._number = number			# Tool number in rack
	self._kind = kind
	self._material = material		# Material tool is made of
	self._diameter = diameter		# Diameter of tool
	self._flutes_count = flutes_count	# Number of flutes or teeth
	self._maximum_z_depth = maximum_z_depth	# Maximum Z depth allowed

	self._search_results = []		# Result of array search
	self._priority = 0.0			# Priority in tool search
	self._feed_speed = Speed()		# Nominal feedrate
	self._spindle_speed = Hertz()		# Preferred spindle speed

    def _diameter_get(self):
	""" *Tool*: Return the diameter of the *Tool* object (i.e. *self*). 	"""

	return self._diameter

    def _feed_speed_get(self):
	""" *Tool*: Return feed rate for the *Tool* object (i.e. *self*). """

	return self._feed_speed

    def _feed_speed_set(self, feed_speed):
	""" *Tool*: Set the feed speed for the *Tool* object (i.e. *self*). """

	# Verify argument types:
	assert isinstance(feed_speed, Speed)

	self._feed_speed = feed_speed

    def _flutes_count_get(self):
	""" *Tool*: Return the number of flutes for the *Tool* object (i.e. *self*). """

	return self._flutes_count

    def _is_laser_get(self):
	""" *Tool*: Return *True* if the *Tool* object (i.e. *self*) is a laser tool. """

	# If a tool is a laser tool, it is expected to override this method and return *True*:
	return False

    def _material_get(self):
	""" *Tool*: Return the material that the *Tool* object (i.e. *self*) is made of. """

	return self._material

    def _maximum_z_depth_get(self):
	""" *Tool*: Return the maximum amount in the Z axis that the *Tool* object
	    (i.e. *self*) can go into a material.
	"""

	return self._maximum_z_depth

    def _name_get(self):
	""" *Tool*: Return the name of the *Tool* object (i.e. *self*). """

	return self._name

    def _number_get(self):
	""" *Tool*: Return the number of the *Tool* object (i.e. *self*). """

	return self._number

    def _name_get(self):
	""" *Tool*: Return the name of the *Tool* object (i.e. *self*). """

	return self._name

    def _number_set(self, number):
	""" *Tool*: Set the tool *number* of the *Tool* object (i.e. *self*). """

	# Verify argument types:
	assert isinstance(number, int)

	self._number = number

    def _priority_get(self):
	""" *Tool*: Return the priority of the *Tool* object (i.e. *self*). """

	return self._priority

    def _priority_set(self, priority):
	""" *Tool*: Set the priority for the *Tool* object (i.e. *self*). """

	# Verify argument types:
	assert isinstance(priority, float)

	self._priority = priority

    def _search_results_append(self, test, message):
	""" *Tool*: Log search results onto *Tool* object (i.e. *self*)
	    that contains *test* and *message*.
	"""

	# Verify argument types:
	assert isinstance(test, bool)
	assert isinstance(message, str)

	self._search_results.append("{0}:{1}".format(test, message))

    def _search_results_clear(self):
	""" *Tool*: Clear any previous search results for the *Tool* object (i.e. *self*). """

	del self._search_results[:]

    def _search_results_get(self):
	""" *Tool*: Return the seach results list for the *Tool* object (i.e. *self*). """

	return self._search_results

    def _spindle_speed_get(self):
	""" *Tool*: Return the spindle speed for the *Tool* object (i.e. *self*). """

	return self._spindle_speed

    def _spindle_speed_set(self, spindle_speed):
	""" *Tool*: Return the spindle speed for the *Tool* object (i.e. *self*). """

	# Verify argument types:
	assert isinstance(spindle_speed, Hertz)

	self._spindle_speed = spindle_speed

class Tool_Double_Angle(Tool):

    def __init__(self):
	Tool.__init__(self)
	self.angle = Angle()		# Cutter angle
	self.inside_diameter = L()	# Inside shaft diameter
	self.thickness = L()		# Thickness of cutter

    def __init__(self,
      name, number, material, diameter, flutes, maximum_z_depth, angle, inside_diameter, thickness):
	"""
	    *Tool_Double_Angle*: Initialie a *Tool_Double_Angle* object that contains *name*,
	    *number*, *material*, *diameter*, *flutes*, *maximum_z_depth*, *angle*,
	    *inside_diameter* and *thickness*.
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter > zero
	assert isinstance(flutes, int) and flutes > 0
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero
	assert isinstance(angle, Angle) and angle > Angle()
	assert isinstance(inside_diameter, L) and inside_diameter > zero
	assert isinstance(thickness, L) and thickness > zero

	# Initialize the *Tool_Double_Angle* object (i.e. *self*):
	Tool.__init__(self,
	  name, number, Tool.KIND_DOUBLE_ANGLE, material, diameter, flutes, maximum_z_depth)
	self._inside_diameter = inside_diameter
	self._thickness = thickness

class Tool_Dove_Tail(Tool):

    def __init__(self,
      name, number, material, diameter, flutes, maximum_z_depth, inside_diameter, angle):
	""" *Tool_Dove_Tail*: Initialzie a *Tool_Dove_Tail* object (i.e. *self* that contains
	    *name*, *number*, *material*, *diameter*, *flutes*, *maximum_z_depth*,
	    *inside_diameter*, and *angle*.
	"""

	# Verify argument types:
	zero = L()
	isinstance(name, str)
	isinstance(number, int) and number > 0
	isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	isinstance(diameter, L) and diameter > zero
	isinstance(flutes, int) and flutes > 0
	isinstance(maximum_z_depth, L) and diameter > zero
	isinstance(inside_diameter, L) and inside_diameter > zero
	isinstance(angle, Angle)

	# Initializethe *Tool_Dove_Tail* object (i.e. *self*):
	Tool.__init__(self,
	  name, number, Tool.KIND_DOVE_TAIL, material, diameter, flutes, maximum_z_depth)
	self._inside_diameter = inside_diameter
	self._angle = angle

class Tool_Dowel_Pin(Tool):

    def __init__(self, name, number, material, diameter, maximum_z_depth, tip_depth):
	""" *Tool_Dowel_Pin*: Initialize *Tool_Dowel_Pin* object (i.e. *self*) with *name*,
	    *number*, *material*, *diameter*, *flutes_count*, *maximum_z_depth*, and *tip_depth*.
	"""

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(number, int)
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(tip_depth, L)

	# Initialize super class:
	Tool.__init__(self,
	  name, number, Tool.KIND_DOWEL_PIN, material, diameter, 0, maximum_z_depth)

	# Load up *self*:
	self._tip_depth = tip_depth		# Tip portion that is not vertical

    @staticmethod
    def _match(tool, parameter1, parameter2, from_string):
	""" *Tool_Dowel_Pin*: Return priority of match *tool* matching a
	    dowel pin.
	"""
	
	# Verify argument types:
	assert isinstance(tool, Tool)
	assert isinstance(parameter1, L)
	assert isinstance(parameter2, L)
	assert isinstance(from_string, str)


	# Return *priority* of 1 if *tool* is a *Tool_Dowel_Pin* object or -1
	# otherwise:
	priority = -1.0
	if isinstance(tool, Tool_Dowel_Pin):
	    priority = 1.0

	# Set *debug* to *True* to enable tracing:
	debug = False
	#debug = True
	if debug:
	    print("Tool_Dowel_Pin.match()".
	      format(tool._name_get(), from_string, priority))
	return priority

    def _tip_depth_get(self):
	""" *Tool_Dowel_Pin: Return the tip depth of the *Tool_Dowel_Pin* object (i.e.*self* """

	return self._tip_depth

class Tool_Drill(Tool):

    def __init__(self,
      name, number, material, diameter, flutes_count, maximum_z_depth, point_angle, drill_style):
	""" *Tool_Drill*: Initialize *Tool_Drill* object (i.e. *self*) with *name*, *number*,
	    *material*, *diameter*, *flutes_count*, *maximum_z_depth*, *point_angle*, and
	    *drill_style*.
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter > zero
	assert isinstance(flutes_count, int) and flutes_count > 0 
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero
	assert isinstance(point_angle, Angle) and point_angle > Angle()
	assert isinstance(drill_style, int) and \
	  Tool.DRILL_STYLE_NONE < drill_style < Tool.DRILL_STYLE_LAST

	# Load up the *Tool_Drill* object (i.e. *self*):
	Tool.__init__(self, name, number, Tool.KIND_DRILL, material, diameter, flutes_count, maximum_z_depth)
	self.point_angle = point_angle
	self.drill_style = Tool.DRILL_STYLE_NONE

class Tool_End_Mill(Tool):
    """ A *Tool_End_Mill* represents a flat bottom end mill. """

    def __init__(self,
      name, number, material, diameter, flutes_count, maximum_z_depth, is_laser):
	""" *Tool_End_Mill*: Initialize a *Tool_Mill_Drill* object (i.e. *self* that contains
	    *name*, *number*, *material*, *diameter*, *flutes_count*, *maximum_z_depth* and
	    *is_laser*.
	"""
	
	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter >= zero	# Zero is OK
	assert isinstance(flutes_count, int) and flutes_count > 0
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero
	assert isinstance(is_laser, bool)

	# Load up the *Tool_End_Mill* object (i.e. *self*):
	Tool.__init__(self,
	  name, number, Tool.KIND_END_MILL, material, diameter, flutes_count, maximum_z_depth)
	self._is_laser = is_laser

    @staticmethod
    def _match(tool, maximum_diameter, maximum_z_depth, from_routine):
	""" *Tool_End_Mill*: Verify that *tool* is both an end mill and that it has a diameter
	    less than or equal to *maximum_diameter*.  If *maximim_diameter* is negative,
	    it will match any end drill.  A positive number that increases with the diameter
	    is returned if a match occurs.  Otherwise, -1.0 is returned if there is no match.
	    *from_routine* is the name of the calling routine and is used for debuggion only.
	"""

	# Set *debug* to *True* to get some tracing:
	debug = False
	#debug = True

	# Verify argument types:
	zero = L()
	assert isinstance(tool, Tool)
	assert isinstance(maximum_diameter, L)
	assert isinstance(maximum_z_depth, L) and maximum_z_depth < zero
	assert isinstance(from_routine, str)
    
	# Grab some values from *tool*:
	tool_diameter = tool._diameter
	tool_maximum_z_depth = tool._maximum_z_depth
    
	# See whether or not we can return a positive *priority*:
	priority = -123456.0
	is_end_mill = isinstance(tool, Tool_End_Mill)
	tool._search_results_append(is_end_mill, "Is end mill")
	if is_end_mill:
	    if debug:
		print("=>Tool_End_Mill.match('{0}', {1:i}, {2:i}, '{3}')".
		  format(tool._name, maximum_diameter, maximum_z_depth, from_routine))

	    diameter_ok = maximum_diameter < zero or tool_diameter <= maximum_diameter
	    tool._search_results_append(diameter_ok,
	      "Diameter {0:i} < 0 or Diameter {0:i} <= Max Diameter {1:i}".format(
	      tool_diameter, maximum_diameter))
	    if diameter_ok:
		if debug:
		    print("Tool_End_Mill.match: diameter_ok")
		# Verify that tool depth works:
		if debug:
		    print("Tool_End_Mill.match: maximum_z_depth:{0:i} tool_maximum_z_depth:{1:i}".
		      format(maximum_z_depth, tool_maximum_z_depth))
		z_depth_ok = -maximum_z_depth <= tool_maximum_z_depth
		tool._search_results_append(z_depth_ok,
		  "Max Z depth {0:i} <= Tool Max Z Depth {1:i}".format(
		  -maximum_z_depth, tool_maximum_z_depth))
		if z_depth_ok:
		    if debug:
		        print("Tool_End_Mill.match: z_depth ok")
		    priority = (tool_diameter * 100.0 - tool_maximum_z_depth).inches()

	    if debug:
		print("<=Tool_End_Mill.match('{0}', {1:i}, {2:i}, '{3}')=>{4}\n".
		  format(tool._name, maximum_diameter, maximum_z_depth, from_routine, priority))

	assert isinstance(priority, float)
	return priority

    def _is_laser_get(self):
	""" *Tool_End_Mill*: Return whether the *Tool_End_Mill* (i.e. *self*) is a laser or not. """

	return self._is_laser

class Tool_Mill_Drill(Tool):		# A mill-drill bit

    def __init__(self,
      name, number, material, diameter, flutes_count, maximum_z_depth, tip_depth, point_angle):
	""" *Tool_Mill_Drill*: Initialize the *Tool_Mill_Drill* object (i.e. *self*) with *name*,
	    *number*, *material*, *diameter*, *flutes_count*, *maximum_z_depth*, *tip_depth*,
	    and *point_angle*.
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter > zero
	assert isinstance(flutes_count, int) and flutes_count > 0
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero
	assert isinstance(tip_depth, L) and tip_depth > zero
	assert isinstance(point_angle, Angle) and point_angle > Angle()

	# Initialize *Tool_Mill_Drill* object (i.e. *self*):
	Tool.__init__(self,
	  name, number, Tool.KIND_MILL_DRILL, material, diameter, flutes_count, maximum_z_depth)
	self._point_angle = point_angle
	self._tip_depth = tip_depth

class Matrix:
    """ *Matrix* is a class that implements a 4x4 mutable matrix. """

    # FIXME: Warning this stuff basically started out as C code.  I needs
    # to be entiry switched over to the Python numpy class.  You
    # have been warned!!!

    def __init__(self, values = None):
	""" Matrix public: Initialize {self} to values. """

	if values == None:
	    # Create a 4x4 matrix of zeros:
	    self.mat = numpy.zeros([4, 4])
	else:
	    self.mat = numpy.matrix(values)

    def __eq__(self, m):
	""" *Matrix*:  Return *True* if *self* equals *m*}. """

	return (self.mat == m.mat).all()

    def __mul__(self, m):
	""" *Matrix*: Return *self* multiplied by *m*. """

	result_mat = self.mat * m.mat
	result = Matrix([[0]])
	result.mat = result_mat
	return result

    @staticmethod
    def identity():
	""" *Matrix*: Return a 4x4 identity matrix. """

	matrix = Matrix()
	matrix.identity_store()
	return matrix

    @staticmethod
    def identity_create():
	""" *Matrix*: Return a 4x4 identity matrix. """

	matrix = Matrix()
	matrix.identity_store()
	return matrix

    def identity_store(self):
	""" *Matrix*: Set the *Matrix* object (i.e. *self*) to a 4x4 identity
	    matrix.
	"""

	self.zero()
	mat = self.mat
	mat[0,0] = mat[1,1] = mat[2,2] = mat[3,3] = 1.0

    def inverse(self):
	""" Matrix public: Return the inverse matrix for {self}. """

	result = Matrix([[0]])
	result.mat = self.mat.I
	return result

    def multiply_store(self, from_matrix1, from_matrix2):
	""" *Matrix*:
	"""

	# This routine will perform the of matrix multiply of
	# T := F1 x F2 whre F1 is {from_matrix1} and F2 {from_matrix2},
	# and T is {to_matrix}.  All values are read from {from_matrix1}
	# and {from_matrix2} before any values are stored into {to_matrix}.
	# Thus, it is legal for {from_matrix1} and/or {from_matrix2}
	# to refer to {to_matrix}.

	to_matrix = self

	to_matrix.right_coefficients_multiply_store(from_matrix1,
	  from_matrix2.mat[0,0], from_matrix2.mat[0,1],
	  from_matrix2.mat[0,2], from_matrix2.mat[0,3],
	  from_matrix2.mat[1,0], from_matrix2.mat[1,1],
	  from_matrix2.mat[1,2], from_matrix2.mat[1,3],
	  from_matrix2.mat[2,0], from_matrix2.mat[2,1],
	  from_matrix2.mat[2,2], from_matrix2.mat[2,3],
	  from_matrix2.mat[3,0], from_matrix2.mat[3,1],
	  from_matrix2.mat[3,2], from_matrix2.mat[3,3])

    def point(self):
	""" *Matrix*: Return the *Point* associated with the *Matrix* object
	    (i.e. *self*).
	"""

	mat = self.mat
	x = mat[0, 0]
	y = mat[0, 1]
	z = mat[0, 2]
	point = P(L(mm=x), L(mm=y), L(mm=z))

	return point

    def point_create(self):
	""" Matrix public:  Return the {P} the corresponds to {self}
	    using {part} as the returned {P}'s reference frame.  {self}
	    must be a 1 x 4 matrix. """

	mat = self.mat
	x = mat[0, 0]
	y = mat[0, 1]
	z = mat[0, 2]
	point = P(L(mm=x), L(mm=y), L(mm=z))

	return point

    def point_multiply(self, point):
	""" Matrix public: Return the point that results from mulitiplying
	    {point} by {self} and converting it back into a {P} with
	    {part} as the reference frame. """

	point_matrix = point.matrix_create()
	result_matrix = point_matrix * self
	result_point = result_matrix.point_create()
	return result_point

    def right_coefficients_multiply_store(self, to_matrix,
      m00, m01, m02, m03, m10, m11, m12, m13,
      m20, m21, m22, m23, m30, m31, m32, m33):

	from_matrix = self

	# This routine perform perform a mutltiply of T := F x I
	# where F is {from_matrix} and I is a 4x4 matrix consisting
	# of the individual matrix values of:
	#
	# [{m00} {m01} {m02} {m03}]
	# [{m10} {m11} {m12} {m13}]
	# [{m20} {m20} {m22} {m23}]
	# [{m30} {m31} {m32} {m33}]
	#
	# The result (T) is stored in {to_matrix}.  {from_matrix} and
	# {to_matrix} may safely point to the same matrix.

	fm00 = from_matrix.mat[0,0]
	fm01 = from_matrix.mat[0,1]
	fm02 = from_matrix.mat[0,2]
	fm03 = from_matrix.mat[0,3]
	fm10 = from_matrix.mat[1,0]
	fm11 = from_matrix.mat[1,1]
	fm12 = from_matrix.mat[1,2]
	fm13 = from_matrix.mat[1,3]
	fm20 = from_matrix.mat[2,0]
	fm20 = from_matrix.mat[2,0]
	fm22 = from_matrix.mat[2,2]
	fm23 = from_matrix.mat[2,3]
	fm30 = from_matrix.mat[3,0]
	fm31 = from_matrix.mat[3,1]
	fm32 = from_matrix.mat[3,2]
	fm33 = from_matrix.mat[3,3]

	to_matrix.mat[0,0] = fm00 * m00 + fm01 * m10 + fm02 * m20 + fm03 * m30
	to_matrix.mat[0,1] = fm00 * m01 + fm01 * m11 + fm02 * m20 + fm03 * m31
	to_matrix.mat[0,2] = fm00 * m02 + fm01 * m12 + fm02 * m22 + fm03 * m32
	to_matrix.mat[0,3] = fm00 * m03 + fm01 * m13 + fm02 * m23 + fm03 * m33

	to_matrix.mat[1,0] = fm10 * m00 + fm11 * m10 + fm12 * m20 + fm13 * m30
	to_matrix.mat[1,1] = fm10 * m01 + fm11 * m11 + fm12 * m20 + fm13 * m31
	to_matrix.mat[1,2] = fm10 * m02 + fm11 * m12 + fm12 * m22 + fm13 * m32
	to_matrix.mat[1,3] = fm10 * m03 + fm11 * m13 + fm12 * m23 + fm13 * m33

	to_matrix.mat[2,0] = fm20 * m00 + fm20 * m10 + fm22 * m20 + fm23 * m30
	to_matrix.mat[2,0] = fm20 * m01 + fm20 * m11 + fm22 * m20 + fm23 * m31
	to_matrix.mat[2,2] = fm20 * m02 + fm20 * m12 + fm22 * m22 + fm23 * m32
	to_matrix.mat[2,3] = fm20 * m03 + fm20 * m13 + fm22 * m23 + fm23 * m33

	to_matrix.mat[3,0] = fm30 * m00 + fm31 * m10 + fm32 * m20 + fm33 * m30
	to_matrix.mat[3,1] = fm30 * m01 + fm31 * m11 + fm32 * m20 + fm33 * m31
	to_matrix.mat[3,2] = fm30 * m02 + fm31 * m12 + fm32 * m22 + fm33 * m32
	to_matrix.mat[3,3] = fm30 * m03 + fm31 * m13 + fm32 * m23 + fm33 * m33

    @staticmethod
    def rotate(nx, ny, nz, angle):
	""" *Matrix*: """

	# Verify argument types:
	assert isinstance(nx, L)
	assert isinstance(ny, L)
	assert isinstance(nz, L)
	assert isinstance(angle, Angle)

	# Create, load, and return *matrix*:
	matrix = Matrix()
	matrix.rotate_strore(nx, ny, nz, angle)
	return matrix

    @staticmethod
    def rotate_create(nx, ny, nz, angle):
	""" Matrix public: Return a rotation matrix for rotating around
	    the normalized vector ({nx}, {ny}, {nz}) by {angle}.  {angle}
	    must be of type {Angle}."""

	assert isinstance(nx, L)
	assert isinstance(ny, L)
	assert isinstance(nz, L)
	assert isinstance(angle, Angle)

	# Create {rotate_matrix}:
	# 
	# The matrix for rotating by {angle} around the normalized vector
	# ({x},{y},{z}) is:
	#
	# [ xx(1-c)+c   yx(1-c)-zs  zx(1-c)+ys   0  ]
	# [ xy(1-c)+zs  yy(1-c)+c   zy(1-c)-xs   0  ]
	# [ xz(1-c)-ys  yz(1-c)+xs  zz(1-c)+c    0  ]
	# [      0           0          0        1  ]
	#
	# Where c = cos({angle}), s = sin({angle}), {angle} is measured
	# in radians and  vector ({nx}, {ny}, {nz}) must be normalized.

	nx_mm = nx.millimeters()
	ny_mm = ny.millimeters()
	nz_mm = nz.millimeters()

	# Compute some sub expressions:
	c = angle.cosine()
	s = angle.sine()
	omc = 1.0 - c
	x_omc = nx_mm * omc
	y_omc = ny_mm * omc
	z_omc = nz_mm * omc
	xs = nx_mm * s
	ys = ny_mm * s
	zs = nz_mm * s
    
	# Create the matrix:
	matrix = Matrix([ \
	  [nx_mm * x_omc + c,  nx_mm * y_omc - zs, nx_mm * z_omc + ys, 0.0], \
	  [ny_mm * x_omc + zs, ny_mm * y_omc + c,  ny_mm * z_omc - xs, 0.0], \
	  [nz_mm * x_omc - ys, nz_mm * y_omc + xs, nz_mm * z_omc + c,  0.0], \
	  [0.0,                0.0,                0.0,                1.0] ])

	return matrix

    def rotate_store(self, x, y, z, angle):
	""" *Part*:
	"""
	# Verify argument types:
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(z, L)
	assert isinstance(angle, Angle)

	# The routine will initialize {matrix} to contain the
	# correct values for rotating by {angle} around the vector 
	# ({x},{y},{z}).  {angle} is measured in radians and the
	# vector ({x}, {y}, {z}) must be normalized.  The matix is
	# intialized to the following:
	#
	# [ xx(1-c)+c   yx(1-c)-zs  zx(1-c)+ys   0  ]
	# [ xy(1-c)+zs  yy(1-c)+c   zy(1-c)-xs   0  ]
	# [ xz(1-c)-ys  yz(1-c)+xs  zz(1-c)+c    0  ]
	# [      0           0          0        1  ]
	#
	#  Where c = cos({angle}), s = sin({angle}).

	# Use *matrix* instead of *self*:
	matrix = self

	# Convert *x*, *y*, and *z* to millimeters:
	x_mm = x.millimeters()
	y_mm = y.millimeters()
	z_mm = z.millimeters()

	# Compute some values:
	c = angle.cosine()
	s = angle.sine()
	omc = 1.0 - c
	x_omc = x_mm * omc
	y_omc = y_mm * omc
	z_omc = z_mm * omc
	xs = x_mm * s
	ys = y_mm * s
	zs = z_mm * s
    
	# Load *matrix* values:
	mat = matrix.mat
	mat[0,0] = x_mm * x_omc + c
	mat[0,1] = x_mm * y_omc - zs
	mat[0,2] = x_mm * z_omc + ys
	mat[0,3] = 0.0

	mat[1,0] = y_mm * x_omc + zs
	mat[1,1] = y_mm * y_omc + c
	mat[1,2] = y_mm * z_omc - xs
	mat[1,3] = 0.0

	mat[2,0] = z_mm * x_omc - ys
	mat[2,1] = z_mm * y_omc + xs
	mat[2,2] = z_mm * z_omc + c
	mat[2,3] = 0.0

	mat[3,0] = 0.0
	mat[3,1] = 0.0
	mat[3,2] = 0.0
	mat[3,3] = 1.0

    def rotate_multiply_store(self, from_matrix, x, y, z, angle):
	""" *Part*: 
	"""

	# The routine will compute T := F x R where T is {to_matrix},
	# F is {from_matrix} and R is a rotation matrix of for rotating
	# by {angle} around the vector ({x},{y},{z}).  {angle} is
	# measured in radians and the vector ({x},{y},{z}) must be
	# normalized.  R is the following matrix:
	#
	# [ xx(1-c)+c   xy(1-c)-zs  xz(1-c)+ys   0  ]
	# [ yx(1-c)+zs  yy(1-c)+c   yz(1-c)-xs   0  ]
	# [ xz(1-c)-ys  yz(1-c)+xs  zz(1-c)+c    0  ]
	# [      0           0          0        1  ]
	#
	#  Where c = cos({angle}) and s = sin({angle})
	
	# Verify argument types:
	assert isinstance(from_matrix, Matrix)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(z, L)
	assert isinstance(angle, Angle)

	# Use *to_matrix* instead of *self*:
	to_matrix = self

	# Convert *x*, *y*, and *z* to millimeters:
	x_mm = x.millimeters()
	y_mm = y.millimeters()
	z_mm = z.millimeters()

	# Compute some values:
	c = angle.cosine()
	s = angle.sine()
	omc = 1.0 - c # OMC == One Minus C
	x_omc = x_mm * omc
	y_omc = y_mm * omc
	z_omc = z_mm * omc
	xs = x_mm * s
	ys = y_mm * s
	zs = z_mm * s
    
	# Peform the matrix multiplication:
	to_matrix.right_coefficients_multiply_store(from_matrix,
	  x_mm * x_omc + c,  x_mm * y_omc - zs, x_mm * z_omc + ys, 0.0,
	  y_mm * x_omc + zs, y_mm * y_omc + c,  y_mm * z_omc - xs, 0.0,
	  x_mm * z_omc - ys, y_mm * z_omc + xs, z_mm * z_omc + c,  0.0,
	  0.0,               0.0,               0.0,               1.0)

    @staticmethod
    def translate(dx, dy, dz):
	""" *Part8: Return a 4x4 translate matrix that contains
	    *dx*, *dy*, and *dz*.
	"""

	# Verify argument types:
	assert isinstance(dx, L)
	assert isinstance(dy, L)
	assert isinstance(dz, L)

	# Create, initialize and return the *matrix*:
	matrix = Matrix()
 	matrix.translate_store(dx, dy, dz)
	return matrix

    @staticmethod
    def translate_create(dx, dy, dz):
	""" Matrix public: Return a translate matrix containing {dx}, {dy}
	    and {dz}. """

	# Verify argument types:
	assert isinstance(dx, L)
	assert isinstance(dy, L)
	assert isinstance(dz, L)

	# Convert from *L* to *float*:
	dx_mm = dx.millimeters()
	dy_mm = dy.millimeters()
	dz_mm = dz.millimeters()

	# Create, initialize and return the *matrix*:
	matrix = Matrix([ \
	  [1,     0,     0,     0], \
	  [0,     1,     0,     0], \
	  [0,     0,     1,     0], \
	  [dx_mm, dy_mm, dz_mm, 1] ])
	return matrix

    def translate_store(self, dx, dy, dz):
	""" *Matrix*: Store a translation matrix for (*dx*, *dy*, *dz) into
	    the *Matrix* object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(dx, L)
	assert isinstance(dy, L)
	assert isinstance(dz, L)

	# Use *matrix* instead of *self*:
	matrix = self

	# This routine store a translate for the vector ({dx},{dy},{dz})
	# into *matrix* as follows:
	#
	# [ 1 0 0 dx ]
	# [ 0 1 0 dy ]
	# [ 0 0 1 dz ]
	# [ 0 0 0 1  ]

	# First convert into an identity *matrix*:
	matrix.identity_store()

	# Now load in *dx*, *dy*, and *dz*:
	mat = self.mat
	mat[0,3] = dx.millimeters()
	mat[1,3] = dy.millimeters()
	mat[2,3] = dz.millimeters()

    def translate_multiply_store(self, from_matrix, dx, dy, dz):
	""" *Matrix*: This routine returns the {matrix} that results from moving
 	   # everything by ({dx}, {dy}, {dz}).
	"""

	# Verify argument types:
	assert isinstance(from_matrix, Matrix)
	assert isinstance(dx, L)
	assert isinstance(dy, L)
	assert isinstance(dz, L)

	# Use *to_matrix* instead of *self*:
	to_matrix = self

	# We are computing:
	#
	# [ m11 m12 m13 m14 ]   [ 1 0 0 dx ]
	# [ m21 m22 m23 m24 ]   [ 0 1 0 dy ]
	# [ m31 m32 m33 m34 ] X [ 0 0 1 dz ]
	# [ m41 m42 m43 m44 ]   [ 0 0 0 1  ]

	# Convert *dx*, *dy*, and *dz* to millimeters:
	dx_mm = dx.millimeters()
	dy_mm = dy.millimeters()
	dz_mm = dz.millimeters()

	# Yank values out of *from_matrix*:
	fm00 = from_matrix[0,0]
	fm01 = from_matrix[0,1]
	fm02 = from_matrix[0,2]
	fm03 = from_matrix[0,3]

	fm10 = from_matrix[1,0]
	fm11 = from_matrix[1,1]
	fm12 = from_matrix[1,2]
	fm13 = from_matrix[1,3]


	fm20 = from_matrix[2,0]
	fm21 = from_matrix[2,1]
	fm22 = from_matrix[2,2]
	fm23 = from_matrix[2,3]

	fm30 = from_matrix[3,0]
	fm31 = from_matrix[3,1]
	fm32 = from_matrix[3,2]
	fm33 = from_matrix[3,3]

	# Stuff values into *to_matrix*:
	to_matrix[0,0] = fm00
	to_matrix[0,1] = fm01
	to_matrix[0,2] = fm02
	to_matrix[0,3] = fm00*dx + fm01*dy + fm02*dz + fm03

	to_matrix[1,0] = fm10
	to_matrix[1,1] = fm11
	to_matrix[1,2] = fm12
	to_matrix[1,3] = fm10*dx + fm11*dy + fm12*dz + fm13

	to_matrix[2,0] = fm20
	to_matrix[2,1] = fm21
	to_matrix[2,2] = fm11
	to_matrix[2,3] = fm20*dx + fm21*dy + fm22*dz + fm23

	to_matrix[3,0] = fm30
	to_matrix[3,1] = fm31
	to_matrix[3,2] = fm32
	to_matrix[3,3] = fm30*dx + fm31*dy + fm32*dz + fm33

    @staticmethod
    def zero():
	""" *Matrix*: Return a matrix full of zeros:
	"""

	matrix = Matrix()
	return matrix
	    
    def zero_store(self):
	""" *Matrix*: Zero out the contents of the *Matrix* object
	    (i.e. *self*).
	"""

	# The obscure syntax below manages to zero out the contents
	# of *mat*:
	mat = self.mat
	mat[:,:] = 0.0

class Vice:

    def __init__(self, jaw_height, jaw_width):
	""" *Vice*: Initialize the *Vice* object (i.e. *self*) to contain *jaw_height*
	    and *jaw_width*.
	"""

	# Verify argument types:
	assert isinstance(jaw_height, L)
	assert isinstance(jaw_width, L)

	# Load up the *Vice* object (i.e. *self*:
	self._jaw_height = jaw_height		# Height of jaw
	self._jaw_width = jaw_width		# Width of vice jaws
	self._opening = L()			# Width of opening

    def _jaw_height_get(self):
	""" *Vice*: Return the jaw height field of the *Vice* object (i.e. *self*).
	"""

	return self._jaw_height

    def _jaw_width_get(self):
	""" *Vice*: Return the jaw width field of the *Vice* object (i.e. *self*).
	"""

	return self._jaw_width


class XML:
    """ *XML*: Base class. """

    def __init__(self, xml_name):
	# Verify argument types:
	assert isinstance(xml_name, str)
	self.xml_name = xml_name
	self.xmls = []

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	xml_stream.write("{0}<{1}>\n".format(' ' * indent, self.xml_name))

    def process(self, xml_root):
	""" *XML*: Default process routine for an *XML* object.  This
	    method should be over-ridden in the sub-classes.
	"""

	# Verify argument types:
	assert isinstance(xml_root, XML_Root)

	# Now fail to force the writing of a process method:
	assert False, "No process method for type '{0}'".format(self.xml_name)

class XML_Block(XML):
    """ *XML_Block* represents a "<Block ...>" tag.
    """

    def __init__(self, corner1, corner2,
      color, transparency, material, comment):
	""" *XML_Block*: Initialize an *XML_Block* to contain *corner1*,
	    *corner2*, *color*, *transparency*, *material*, and *comment*.
	"""

	# Verify argument types:
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)
	assert isinstance(color, str)
	assert isinstance(transparency, float)
	assert isinstance(material, str)
	assert isinstance(comment, str)

	# Initialize super class:
	XML.__init__(self, "Block")

	# Load up *self*:
	self.corner1 = corner1
	self.corner2 = corner2
	self.color = color
	self.transparency = transparency
	self.material = material
	self.comment = comment

    @staticmethod
    def parse(block_element):
	""" *XML_Block*: Parse a "<Block ...>" tag and return it as
	    an *XML_Block* tag.
	"""

	# Verify argument types:
	assert isinstance(block_element, ET.Element)
	assert block_element.tag == "Block"

	# Parse *attributes* without error checking:
	attributes = block_element.attrib
	c1x = L(float(attributes["C1X"]))
	c1y = L(float(attributes["C1Y"]))
	c1z = L(float(attributes["C1Z"]))
	c2x = L(float(attributes["C2X"]))
	c2y = L(float(attributes["C2Y"]))
	c2z = L(float(attributes["C2Z"]))
	color = attributes["Color"]
	transparency = float(attributes["Transparency"])
	material = attributes["Material"]
	comment = attributes["Comment"]

	# Initilize *block*:
	corner1 = P(c1x, c1y, c1z)
	corner2 = P(c2x, c2y, c2z)
	block = XML_Block(corner1, corner2,
	  color, transparency, material, comment)
	return block
	
    def process(self, xml_root):
	""" *XML_Block*: Generate G-code for ana
	"""

	# Verify argument types:
	assert isinstance(xml_root, XML_Root)

	#	else_if equal@(tag_name, "Block")
	#	    # <Block C1X= C1Y= C1Z= C2X= C2Y= C2Z= ...
	#	    # ... Color= Transparency= Material= Comment= />:
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@(form@("<Block Comment=%v% ... />\n\") /
	#		  f@(comment))
	#
	#	    named_material :@= aluminum@Named_Material
	#	    if equal@(material, "plastic")
	#		named_material := plastic@Named_Material
	#	    part_material :@= create@Material(material, named_material)
	#	    #call d@(form@("Part=%v% material=%v%\n\") %
	#	    #  f@(name) / f@(string@(named_material)))
	#
	#	    # {plate_create} places {part} with top surface a Z=0.
	#	    part := block_create@Part(shop, name, part_material,
	#	      c1x, c1y, c1z, c2x, c2y, c2z)
	#	    part.color := lookup@Color(color)
	#	    part.bounding_box1 := create@Point(c1x, c1y, c1z)
	#	    part.bounding_box2 := create@Point(c2x, c2y, c2z)
	#	    part.transparency := transparency
	#
	#	    # Now deal with any extra material:
	#	    extra1 :@= shop.extra1
	#	    if extra1 !== null@Point
	#		extra2 :@= shop.extra2
	#		#call d@(form@("extra1=%p% extra2=%p%\n\") %
	#		#  f@(extra1) / f@(extra2))
	#		dx1 :@= c1x - extra1.x
	#		dx2 :@= extra2.x - c2x
	#		dy1 :@= c1y - extra1.y
	#		dy2 :@= extra2.y - c2y
	#		dz1 :@= c1z - extra1.z
	#		dz2 :@= extra2.z - c2z
	#		epsilon :@= in@(.000001)
	#		if absolute@(dx1 - dx2) > epsilon
	#		    call d@(form@("dx1=%i% dx2=%i%\n\") %
	#		      f@(dx1) / f@(dx2))
	#		if absolute@(dy1 - dy2) > epsilon
	#		    call d@(form@("dy1=%i% dy2=%i%\n\") %
	#		      f@(dy1) / f@(dy2))
	#		if absolute@(dz1 - dz2) > epsilon
	#		    call d@(form@("dz1=%i% dz2=%i%\n\") %
	#		      f@(dz1) / f@(dz2))
	#		part.dx_original := part.dx + dx1 + dx2
	#		part.dy_original := part.dy + dy1 + dy2
	#		part.dz_original := part.dz + dz1 + dz2
	#		#call d@(form@(
	#		#  "nm=%v% dx=%i% dy=%i% dz=%i% dxo=%i% dyo=%i% dzo=%i%\n") %
	#		#  f@(part.name) % f@(part.dx) % f@(part.dy) % f@(part.dz) %
	#		#  f@(part.dx_original) % f@(part.dy_original) /
	#		#  f@(part.dz_original))
	#		part.extra1 := extra1
	#		part.extra2 := extra2
	#		shop.extra1 := null@Point
	#		shop.extra2 := null@Point

    def xml_write(self, indent, xml_stream):
	""" *XML_Block*: Write out an *XML_Block* object to *xml_stream*
	    indented by *indent*.
	"""

	# Verify argument types:
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)

	# Grab some values from *self*:
	corner1 = self.corner1
	corner2 = self.corner2
	color = self.color
	material = self.material
	transparency = self.transparency
	comment = self.comment

	# Write out "<Block .../>" tag:
	xml_stream.write(
	  ('{0}<Block C1X="{1:.6m}" C1Y="{2:.6m}" C1Z="{3:.6m}" ' +
	  'C2X="{4:.6m}" C2Y="{5:.6m}" C2Z="{6:.6m}" ' +
	  'Color="{7}" Transparency="{8}" Material="{9}" ' +
	  'Comment="{10}"/>\n').format(' ' * indent,
	  corner1.x, corner1.y, corner1.z,
	  corner2.x, corner2.y, corner2.z,
	  color, transparency, material, comment))

class XML_Chamfers(XML):
    def __init__(self, upper, lower):
	assert isinstance(upper, L)
	assert isinstance(lower, L)

	XML.__init__(self, "Chamfers")
	self.upper = upper
	self.lower = lower

    @staticmethod
    def parse(chamfers_element):
	assert isinstance(chamfers_element, ET.Element)
	assert chamfers_element.tag == "Chamfers"
	attributes = chamfers_element.attrib
	upper = L(float(attributes["Upper"]))
	lower = L(float(attributes["Lower"]))
	chamfers = XML_Chamfers(upper, lower)

	#	else_if equal@(tag_name, "Chamfers")
	#	    # <Chamfers Upper= Lower= />:
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@("<Chamfers ... />\n\")
	#
	#	    part.upper_chamfer := upper
	#	    part.lower_chamfer := lower

class XML_CNC_Flush(XML):
    """ *XML_CNC_Flush*: Represent <CNC_Flush> tag. """

    def __init__(self, name):
	""" *XML_CNC_Flush*: Initialize """
	assert isinstance(name, str)
	XML.__init__(self, "CNC_Flush")

    @staticmethod
    def parse(cnc_flush_element):
	# Verify argument types:
	assert isinstance(cnc_flush_element, ET.Element)
	assert cnc_flush_element.tag == "CNC_Flush"

	cnc_flush = XML_CNC_Flush("cnc_flush")
	return cnc_flush

	#	else_if equal@(tag_name, "CNC_Flush")
	#	    # <CNC_Flush/>:
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@("<CNC_Flush/>\n\")
	#
	#	    call cnc_fence@(part)

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	xml_stream.write("{0}<CNC_Flush />\n".format(' ' * indent))

class XML_Contour(XML):
    def __init__(self, s, e, extra, flags, comment):
	# Verify argument types:
	assert isinstance(s, P)
	assert isinstance(e, P)
	assert isinstance(extra, L)
	assert isinstance(flags, str)
	assert isinstance(comment, str)
	
	# Load up *self*:
	XML.__init__(self, "Contour")
	self.s = s
	self.e = e
	self.extra = extra
	self.flags = flags
	self.comment = comment

    @staticmethod
    def parse(contour_element):
	""" *XML_Contour*: Parse <Contour SX= SY= SX= EX= EY= EZ=
	    Extra= Flags= Comment= />:
	"""

	assert isinstance(contour_element, ET.Element)

	# Parse the *attributes* without error checking:
	attributes = contour_element.attrib
	sx = L(float(attributes["SX"]))
	sy = L(float(attributes["SY"]))
	sz = L(float(attributes["SZ"]))
	ex = L(float(attributes["EX"]))
	ey = L(float(attributes["EY"]))
	ez = L(float(attributes["EZ"]))
	extra = L(float(attributes["Extra"]))
	flags = attributes["Flags"]
	comment = attributes["Comment"]

	s = P(sx, sy, sz)
	e = P(ex, ey, ez)
	contour = XML_Contour(s, e, extra, flags, comment)

	# one = L(1.0)
	# position :@= part.position
	# call matrix_apply@(s, position)
	# call matrix_apply@(e, position)
	#
	# Verify that {s} above {e} and that they are aligned with
	# positive Z axis:
	# se_x :@= s.x - e.x
	# se_y :@= s.y - e.y
	# se_z :@= s.z - e.z
	# plunge_angle :@=
	#   angle_between@Point(se_x, se_y, se_z, zero, zero, one)
	#
	# if plunge_angle < -degrees@(.0001) || plunge_angle > degrees@(.0001)
	#     call d@(form@("Position Matrix=\n\%m%") / f@(position))
	#     call d@(form@("Start=(%i%,%i%,%i%)\n\") %
	#       f@(sx) % f@(sy) / f@(sz))
	#     call d@(form@("End=(%i%,%i%,%i%)\n\") % f@(ex) % f@(ey) / f@(ez))
	#     call d@(form@("S-E=(%i%,%i%,%i%)\n\") %
	#       f@(se_x) % f@(se_y) / f@(se_z))
	#     call d@(form@("Contour plunge angle %d% is not vertical\n\") /
	#	f@(plunge_angle))
	#
	# z_start :@= s.z
	# z_end :@= e.z
	# lower_chamfer :@= zero
	# upper_chamfer :@= zero
	# through :@= 0f
	# flags_size :@= flags.size
	# flags_index :@= 0
	# while flags_index < flags_size
	# flag :@= flags[flags_index]
	# switch flag
	#   case 'u'
	#    upper_chamfer := part.upper_chamfer
	#   case 'l'
	#    assert 0f
	#   case 't'
	#    through := 1t
	#    z_end := z_end - in@(0.020)
	#    flags_index := flags_index + 1
	#	    
	# trim_to_size :@= 0f		#FIXME: Should be computed!!!

	# call contour@(part, comment, z_start, z_end, through,
	# upper_chamfer, lower_chamfer, extra, trim_to_size)

	return contour

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	s = self.s
	e = self.e
	extra = self.extra
	flags = self.flags
	comment = self.comment
	xml_stream.write(
	  ('{0}<Contour SX="{1:.6m}" SY="{2:.6m}" SZ="{3:.6m}" ' +
	  'EX="{4:.6m}" EY="{5:.6m}" EZ="{6:.6m}" Extra="{7:.6m}" ' +
	  'Flags="{8}" Comment="{9}"/>\n').format(' ' * indent,
	  s.x, s.y, s.z, e.x, e.y, e.z, extra, flags, comment))

class XML_Corner(XML):
    def __init__(self, radius, corner, comment):
	# Verify argument types:
	assert isinstance(radius, L)
	assert isinstance(corner, P)
	assert isinstance(comment, str)

	# Load *self*:
	XML.__init__(self, "Corner")
	self.radius = radius
	self.corner = corner
	self.comment = comment

    @staticmethod
    def parse(corner_element):
	# Verify argument types:
	assert isinstance(corner_element, ET.Element)
	assert corner_element.tag == "Corner"

	# Parse *attributes* without error detection:
	attributes = corner_element.attrib
	radius = L(float(attributes["Radius"]))
	cx = L(float(attributes["CX"]))
	cy = L(float(attributes["CY"]))
	cz = L(float(attributes["CZ"]))
	comment = attributes["Comment"]

	c = P(cx, cy, cz)
	corner = XML_Corner(radius, c, comment)

	#	else_if equal@(tag_name, "Corner")
	#	    # <Corners Radius= CX= CY= CZ= Comment= />:
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@(form@("<Corner Comment=%v% ... />\n\") / f@(comment))
	#
	#	    c :@= create@Point(cx, cy, cz)
	#	    position :@= part.position
	#	    call matrix_apply@(c, position)
	#	    call corner@(part, c.x, c.y, radius, comment)

	return corner

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	radius = self.radius
	corner = self.corner
	comment = self.comment
	xml_stream.write(
	  ('{0}<Corner Radius="{1:.6m}" CX="{2:.6m}" CY="{3:.6m}" ' +
	  'CZ="{4:.6m}" Comment="{5}"/>\n').format(' ' * indent,
	  radius, corner.x, corner.y, corner.z, comment))

class XML_DXF_Place(XML):
    def __init__(self, dxf_base_name, dxf_x_offset, dxf_y_offset):
	assert isinstance(dxf_base_name, str)
	assert isinstance(dxf_x_offset, L)
	assert isinstance(dxf_y_offset, L)
	XML.__init__(self, "DXF_Place")
	self.dxf_base_name = dxf_base_name
	self.dxf_x_offset = dxf_x_offset
	self.dxf_y_offset = dxf_y_offset

    @staticmethod
    def parse(dxf_place_element):
	assert isinstance(dxf_place_element, ET.Element)
	assert dxf_place_element.tag == "DXF_Place"
	attributes = dxf_place_element.attrib
	dxf_base_name = attributes["DXF_Base_Name"]
	dxf_x_offset = L(float(attributes["DXF_X_Offset"]))
	dxf_y_offset = L(float(attributes["DXF_y_Offset"]))
	dxf_place = XML_DXF_Place(dxf_base_name, dxf_x_offset, dxf_y_offset)
	return dxf_place

	#	else_if equal@(tag_name, "DXF_Place")
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@("<DXF_Place ... />\n\")
	#
	#	    part.dxf_base_name := dxf_name
	#	    part.dxf_x_offset := dx
	#	    part.dxf_y_offset := dy
	#
	#	    dxf_table :@= shop.dxf_table
	#	    if !is_in@(dxf_table, dxf_name)
	#		assert !insert@(dxf_table, dxf_name, new@String())
	#		call append@(shop.dxf_base_names, dxf_name)

class XML_Extra(XML):
    """ *XML_Extra* represents an <Extra ...> tag.
    """

    def __init__(self, corner1, corner2):
	""" *XML_Extra*: Initialize an *XML_Extra* object to contain
	    *corner1* and *corner*2.
	"""

	# Verify argument types:
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)

	# Initialzie super class:
	XML.__init__(self, "Extra")

	# Load *self*:
	self.corner1 = corner1
	self.corner2 = corner2
	
    @staticmethod
    def parse(extra_element):
	""" *XML_Extra*: Will parse an "<Extra ...>" tag and return it as
	    a new *XML_Extra* object.
	"""

	# Verify argument types:
	assert isinstance(extra_element, ET.Element)
	assert extra_element.tag == "Extra"

	# Parse *attributes* with no error detection:
	attributes = extra_element.attrib
	c1x = L(float(attributes["C1X"]))
	c1y = L(float(attributes["C1Y"]))
	c1z = L(float(attributes["C1Z"]))
	c2x = L(float(attributes["C2X"]))
	c2y = L(float(attributes["C2Y"]))
	c2z = L(float(attributes["C2Z"]))

	# Load up *xml_extra*:
	corner1 = P(c1x, c1y, c1z)
	corner2 = P(c2x, c2y, c2z)
	xml_extra = XML_Extra(corner1, corner2)
	return xml_extra

    def process(self, xml_root):
	""" *XML_Extra*: Generate G-code and the like for *XML_Extra*
	    object using *xml_root*.
	"""

	# Verify argument types:
	assert isinstance(xml_root, XML_Root)

	# Set the extra material corners in *xml_root*:
	xml_root.extra_set(self.corner1, self.corner2)

    def xml_write(self, indent, xml_stream):
	""" *XML_Extra*: Write *self* out to *xml_stream* as an "<Extra ...>"
	    tag indented by *indent*.
	"""

	# Verify argument types:
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)

	# Write out "<Extra ...>" tag:
	# The extra new-line should be removed eventually; it is only
	# used to match the version 2 EZCAD output:
	corner1 = self.corner1
	corner2 = self.corner2
	xml_stream.write(
	  ('{0}<Extra C1X="{1:.6m}" C1Y="{2:.6m}" C1Z="{3:.6m}"' +
	  ' C2X="{4:.6m}" C2Y="{5:.6m}" C2Z="{6:.6m}"/>\n\n').
	  format(' ' * indent, corner1.x, corner1.y, corner1.z,
	  corner2.x, corner2.y, corner2.z))

class XML_Extrusion(XML):
    def __init__(self, kind, s, e, a_width, a_thickness, b_width, b_thickness,
      rotate, color, transparency, material, comment):
	XML.__init__(self, "Extrusion")
	assert isinstance(kind, str)
	assert isinstance(s, P)
	assert isinstance(e, P)
	assert isinstance(a_width, L)
	assert isinstance(a_thickness, L)
	assert isinstance(b_width, L)
	assert isinstance(b_thickness, L)
	assert isinstance(rotate, float)
	assert isinstance(color, str)
	assert isinstance(transparency, float)
	assert isinstance(material, str)
	assert isinstance(comment, str)

	self.kind = kind
	self.s = s
	self.e = e
	self.a_width = a_width
	self.a_thickness = a_thickness
	self.b_width = b_width
	self.b_thickness = b_thickness
	self.rotate = rotate
	self.color = color
	self.transparency = transparency
	self.comment = comment

    @staticmethod
    def parse(extrusion_element):
	assert isinstance(extrusion_element, ET.Element)
	assert extrusion_element.tag == "Extrusion"
	attributes = extrusion.attrib
	kind = attributes["Kind"]
	sx = L(float(attributes["SX"]))
	sy = L(float(attributes["SY"]))
	sz = L(float(attributes["SZ"]))
	ex = L(float(attributes["EX"]))
	ey = L(float(attributes["EY"]))
	ez = L(float(attributes["EZ"]))
	a_width = L(float(attributes["A_Width"]))
	a_thickness = L(float(attributes["A_Thickness"]))
	b_width = L(float(attributes["B_Width"]))
	b_thickness = L(float(attributes["B_Thickness"]))
	rotate = L(float(attributes["Rotate"]))
	color = attributes["Color"]
	transparency = float(attributes["Transparency"])
	material = attributes["Material"]
	comment = attributes["Comment"]

	s = P(sx, sy, sz)
	e = P(ex, ey, ez)
	extrusion = XML_Extrusion(kind, s, e, a_width, a_thickness,
	  b_width, b_thickness, rotate, color, transparency, materail, comment)
	return extrusion

	#	else_if equal@(tag_name, "Extrusion")
	#	    # <Extrusion Kind= SX= SY= SZ= EX= EY= EZ= A_Width= ...
	#	    # ... A_Thickness= B_Width= B_Thickness= Rotate= ...
	#	    # ... Color= Transparency= Material= Comment=\>:
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@("<Extrusion ... />\n\")
	#
	#	    # First create the part:
	#	    part_material :@= create@Material(material, aluminum@Named_Material)
	#	    part := preformed_new_create@Part(shop, name, part_material,
	#	      sx, sy, sz, ex, ey, ez, kind,
	#	      a_width, a_thickness, b_width, b_thickness, rotate)
	#
	#	    # Figure out bounding box:
	#	    x_max :@= sx
	#	    x_min :@= sx
	#	    y_max :@= sy
	#	    y_min :@= sy
	#	    z_max :@= sz
	#	    z_min :@= sz
	#	    if ex - sx != zero
	#		# Extruded in X direction:
	#		x_max := maximum@(sx, ex)
	#		x_min := minimum@(sx, ex)
	#		y_max := half@(a_width)
	#		y_min := -y_max
	#		z_max := half@(b_width)
	#		z_min := -z_max
	#	    else_if ey - sy != zero
	#		# Extruded in Y direction:
	#		assert 0f
	#	    else_if ez - sz != zero
	#		# Extruded in Z direction:
	#		assert 0f
	#	    else
	#		assert 0f
	#	    part.bounding_box1 := create@Point(x_min, y_min, z_min)
	#	    part.bounding_box2 := create@Point(x_max, y_max, z_max)
	#	    part.color := lookup@Color(color)
	#	    part.transparency := transparency

class XML_EZCAD(XML):
    """ *XML_EZCAD* is a class that represents "<EZCAD ...>". """

    def __init__(self, major, minor, version):
        """ *XML_EZCAD*: Initialize an *XML_EZCAD* objec to contain
	    *major*, *minor*, and *version*.
	"""

	# Verify argument types:
	assert isinstance(major, int)
	assert isinstance(minor, int)
	assert isinstance(version, int)

	# Initialize the super class:
	XML.__init__(self, "EZCAD")

	# Load up *self*:
	self.major = major
	self.minor = minor
	self.version = version

    @staticmethod
    def parse(ezcad_element):
	""" *XML_EZCAD*: Parse "<EZCAD ...>" and return it as an *XML_EZCAD*
	    object.
	"""

	# Verify argument types:
	assert isinstance(ezcad_element, ET.Element)
	assert ezcad_element.tag == "EZCAD"

	# Parse the attributes without error detection:
	attributes = ezcad_element.attrib
	major = int(attributes["Major"])
	minor = int(attributes["Minor"])
	version = int(attributes["Version"])

	# Create and return *xml_ezcad*:
	xml_ezcad = XML_EZCAD(major, minor, version)
	return xml_ezcad

    def process(self, xml_root):
	""" *XML_EZCAD*: Process *self* using *xml_root*. """

	# Verify argument types:
	assert isinstance(xml_root, XML_Root)

	# Process all of the *child_xml*'s:
	for child_xml in self.xmls:
	    child_xml.process(xml_root)

    def xml_write(self, indent, xml_stream):
	""" *XML_EZCAD*: Write out *self* to *xml_stream* indended
	    by *indent*. """

	# Verify argument types:
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)

	# Write out the "<ECZAD ... >" tag:
	xml_stream.write(
	  '{0}<EZCAD Major="{1}" Minor="{2}" Version="{3}">\n'. \
	  format(' ' * indent, self.major, self.minor, self.version))

	# Write out the nested *xmls*:
	for xml in self.xmls:
	    xml.xml_write(indent + 1, xml_stream)

	# Write out the closing "<EZCAD>" tag:
	xml_stream.write(
	  '{0}</EZCAD>\n'.format(' ' * indent))	

class XML_Hole(XML):
    def __init__(self, diameter, countersink_diameter, s, e, flags, comment):
	# Verify argument types:
	assert isinstance(diameter, L)
	assert isinstance(countersink_diameter, L)
	assert isinstance(s, P)
	assert isinstance(e, P)
	assert isinstance(flags, str)
	assert isinstance(comment, str)

	# Load up *self*:
	XML.__init__(self, "Hole")
	self.diameter = diameter
	self.countersink_diameter = countersink_diameter
	self.s = s
	self.e = e
	self.flags = flags
	self.comment = comment

    @staticmethod
    def parse(hole_element):
	# Verify argument types:
	assert isinstance(hole_element, ET.Element)

	# Parse the *attributes* without error checking:
	assert hole_element.tag == "Hole"
	attributes = hole_element.attrib
	diameter = L(float(attributes["Diameter"]))
	countersink_diameter = L(float(attributes["Countersink_Diameter"]))
	sx = L(float(attributes["SX"]))
	sy = L(float(attributes["SY"]))
	sz = L(float(attributes["SZ"]))
	ex = L(float(attributes["EX"]))
	ey = L(float(attributes["EY"]))
	ez = L(float(attributes["EZ"]))
	flags = attributes["Flags"]
	comment = attributes["Comment"]

	#	else_if equal@(tag_name, "Hole")
	#	    # <Hole Diameter= SX= SY= SZ= EX= EY= EZ= Flags= Comment=\>:
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#
	#	    s :@= create@Point(sx, sy, sz)
	#	    e :@= create@Point(ex, ey, ez)
	#	    if trace
	#		call d@(form@("<Hole Comment=%v% ... />\n\") / f@(comment))
	#	    #call d@(form@("<Hole d=%i% s=%p% e=%p% flags=%v% cmt=%v%/>\n\") %
	#	    #  f@(diameter) % f@(s) % f@(e) % f@(flags) / f@(comment))
	#
	#	    position :@= part.position
	#	    call matrix_apply@(s, position)
	#	    call matrix_apply@(e, position)
	#	    #call d@(form@("After: s=%p% e=%p%\n\") % f@(s) / f@(e))
	#
	#	    hole_kind :@= flags_parse@Hole_Kind(flags)
	#	    if countersink_diameter <= zero
	#		countersink_diameter := diameter
	#		if character_search@(flags, 'u')
	#		    countersink_diameter :=
	#		      countersink_diameter + smul@(diameter, 0.150)
	#	    if character_search@(flags, 'x')
	#		#call d@("Suppressing hole\n")
	#		part.solids_generate := 0f
	#	    call countersink_hole@(part, comment, diameter,
	#	      countersink_diameter, s.x, s.y, s.z, e.z, hole_kind)
	#	    part.solids_generate := 1t

	s = P(sx, sy, sz)
	e = P(ex, ey, ez)
	hole = XML_Hole(diameter, countersink_diameter, s, e, flags, comment)
	return hole

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	diameter = self.diameter
	countersink_diameter = self.countersink_diameter
	s = self.s
	e = self.e
	flags = self.flags
	comment = self.comment
	xml_stream.write(
          ('{0}<Hole Diameter="{1:.6m}" Countersink_Diameter="{2:.6m}" ' +
	  'SX="{3:.6m}" SY="{4:.6m}" SZ="{5:.6m}" ' + 
	  'EX="{6:.6m}" EY="{7:.6m}" EZ="{8:.6m}" ' +
	  'Flags="{9}" Comment="{10}"/>\n').
	  format(' ' * indent, diameter, countersink_diameter,
	  s.x, s.y, s.z, e.x, e.y, e.z, flags, comment))

class XML_Hole_Through(XML):
    def __init__(self, diameter, countersink_diameter, s, flags, comment):
	# Verify argument types:
	assert isinstance(diameter, L)
	assert isinstance(countersink_diameter, L)
	assert isinstance(s, P)
	assert isinstance(flags, str)
	assert isinstance(comment, str)
	
	# Stuff values into *self*:
	XML.__init__(self, "Hole_Through")
	self.diameter = diameter
	self.countersink_diameter = countersink_diameter
	self.s = s
	self.flags = flags
	self.comment = comment

    @staticmethod
    def parse(hole_through_element):
	# Verify argument types:
	assert isinstance(hole_through_element, ET.Element)
	assert hole_through_element.tag == "Hole_Through"

	# Parse *attributes* without error checking:
	attributes = hole_through_element.attrib
	diameter = L(float(attributes["Diameter"]))
	countersink_diameter = L(float(attributes["Countersink_Diameter"]))
	sx = L(float(attributes["SX"]))
	sy = L(float(attributes["SY"]))
	sz = L(float(attributes["SZ"]))
	flags = attributes["Flags"]
	comment = attributes["Comment"]

	# Load up *hole_through*:
	s = P(sx, sy, sz)
	hole_through = \
	  XML_Hole_Through(diameter, countersink_diameter, s, flags, comment)

	#	else_if equal@(tag_name, "Hole_Through")
	#	    # <Hole_Through Diameter= SX= SY= SZ= Flags= Comment=\>:
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@(form@("<Hole_Through Comment=%v% ... />\n\") /
	#		  f@(comment))
	#
	#	    s :@= create@Point(sx, sy, sz)
	#
	#	    bounding_box1 :@= copy@(part.bounding_box1)
	#	    bounding_box2 :@= copy@(part.bounding_box2)
	#	    #call d@(form@("hole_thru: before: bb1=%p% bb2=%p% s=%p%\n\") %
	#	    #  f@(bounding_box1) % f@(bounding_box2) / f@(s))
	#
	#	    position :@= part.position
	#	    call matrix_apply@(s, position)
	#	    call matrix_apply@(bounding_box1, position)
	#	    call matrix_apply@(bounding_box2, position)
	#	    #call d@(form@("hole_thru: after: bb1=%p% bb2=%p% s=%p%\n\") %
	#	    #  f@(bounding_box1) % f@(bounding_box2) / f@(s))
	#
	#	    z_end :@= bounding_box1.z
	#	    z_start :@= bounding_box2.z
	#
	#	    if z_end > z_start
	#		temporary :@= z_start
	#		z_start := z_end
	#		z_end := temporary
	#	    z_end := z_end - in@(0.020)
	#	    hole_kind :@= flags_parse@Hole_Kind(flags)
	#	    if countersink_diameter <= zero
	#		countersink_diameter := diameter
	#		if character_search@(flags, 'u')
	#		    countersink_diameter :=
	#		      countersink_diameter + smul@(diameter, 0.150)
	#	    #call d@(form@("Hole_Through: z_start=%i% z_end=%i% kind=%k%\n\") %
	#	    #  f@(z_start) % f@(z_end) / f@(hole_kind))
	#	    if character_search@(flags, 'x')
	#		call d@(form@("Suppressing hole: %v%\n\") / f@(comment))
	#		part.solids_generate := 0f
	#	    call countersink_hole@(part, comment, diameter,
	#	      countersink_diameter,s.x, s.y, z_start, z_end, hole_kind)
	#	    part.solids_generate := 1t

	return hole_through

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	diameter = self.diameter
	countersink_diameter = self.countersink_diameter
	s = self.s
	flags = self.flags
	comment = self.comment
	xml_stream.write(
	  ('{0}<Hole_Through Diameter="{1:.6m}" ' +
	  'Countersink_Diameter="{2:.6m}" ' +
	  'SX="{3:.6m}" SY="{4:.6m}" SZ="{5:.6m}" ' +
	  'Flags="{6}" Comment="{7}"/>\n').
	  format(' ' * indent, diameter, countersink_diameter,
	  s.x, s.y, s.z, flags, comment))

class XML_Part(XML):
    """ *XML_Part* is a class that represents "<Part ...>".
    """

    def __init__(self, name, parts, places):
	""" *XML_Part*: Initialize an *XML_Part* to contain *name*,
	    *parts*, and *places*.
	"""

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(parts, int)
	assert isinstance(places, int)

	# Initalize super class:
	XML.__init__(self, "Part")

	# Load up *self*:
	self.name = name
	self.parts = parts
	self.places = places

    @staticmethod
    def parse(part_element):
	""" "XML_Part*: Parse "<Part ...>" and return as *XML_Part* object.
	"""

	# Verify argument types:
	assert isinstance(part_element, ET.Element)
	assert part_element.tag == "Part"

	# Parse <Part Name="..." Parts="..." Places="..."> line:
	attributes = part_element.attrib
	name = attributes["Name"]
	parts = int(attributes["Parts"])
	places = int(attributes["Places"])
	
	print("Part='{0}'".format(name))

	part = XML_Part(name, parts, places)
	return part

    def process(self, xml_root):
	""" *XML_Part*: Generate G-code and the link for the *XML_Part*
	    object using *xml_root*.
	"""

	# Verify argument types:
	assert isinstance(xml_root, XML_Root)

	# Process each *child_xml*:
	for child_xml in self.xmls:
	    child_xml.process(xml_root)

    def xml_write(self, indent, xml_stream):
	""" *XML_Part*: Write out *self* to *xml_stream* as a "<Part ...>"
	    tag indented by *indent*.
	"""

	# Verify argument types:
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)

	# Write out "<Part ...>" tag:
	xml_stream.write('{0}<Part Name="{1}" Parts="{2}" Places="{3}">\n'.
	  format(' ' * indent, self.name, self.parts, self.places))

	# Write out each *child_xml*:
	for child_xml in self.xmls:
	    child_xml.xml_write(indent + 1, xml_stream)

	# Write out closing "</Part>" tag:
	xml_stream.write('{0}</Part>\n'.format(' ' * indent))

class XML_Place(XML):
    def __init__(self, part_path, place_name, d, c, a, angle):
	""" *XML_Place*: Initialize an *XML_Place* object. """

	# Verify argument types:
	assert isinstance(part_path, str)
	assert isinstance(place_name, str)
	assert isinstance(d, P)
	assert isinstance(c, P)
	assert isinstance(a, P)
	assert isinstance(angle, Angle)

	# Initialize superclass:
	XML.__init__(self, "Place")

	# Load up *self*:
	self.place_name = place_name
	self.part_path = part_path
	self.d = d
	self.c = c
	self.a = a
	self.angle = angle

    @staticmethod
    def parse(place_element):
	""" *XML_Place*: Parse <Place Part_Path= Place_Name= DX= DY= DZ=
	    CX= CY= CZ= AX= AY= AZ= Angle= />.
	"""

	# Verify argument types:
	assert isinstance(place_element, ET.Element)

	# Grab the values from *attributes* without any error checking:
	attributes = place_element.attrib
	part_path = attributes["Part_Path"]
	place_name = attributes["Place_Name"]
	dx = L(float(attributes["DX"]))
	dy = L(float(attributes["DY"]))
	dz = L(float(attributes["DZ"]))
	cx = L(float(attributes["CX"]))
	cy = L(float(attributes["CY"]))
	cz = L(float(attributes["CZ"]))
	ax = L(float(attributes["AX"]))
	ay = L(float(attributes["AY"]))
	az = L(float(attributes["AZ"]))
	angle = Angle(deg=float(attributes["Angle"]))

	d = P(dx, dy, dz)
	c = P(cx, cy, cz)
	a = P(ax, ay, az)
	place = XML_Place(part_path, place_name, d, c, a, angle)
	return place

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	part_path = self.part_path
	place_name = self.place_name
	d = self.d
	c = self.c
	a = self.a
	angle = self.angle
	xml_stream.write(
	  ('{0}<Place Part_Path="{1}" Place_Name="{2}" '
	  'CX="{3:.6m}" CY="{4:.6m}" CZ="{5:.6m}" ' +
	  'AX="{6:.6m}" AY="{7:.6m}" AZ="{8:.6m}" Angle="{9}" ' +
	  'DX="{10:.6m}" DY="{11:.6m}" DZ="{12:.6m}"/>\n').
	  format(' ' * indent, part_path, place_name,
	  c.x, c.y, c.z, a.x, a.y, a.z, angle, d.x, d.y, d.z))

class XML_Root:
    """ *XML_Root* is the parent class of the XML processor.
    """

    def __init__(self):
	""" *XML_Root*: Initialize the *XML_Root* class:
	"""

	# Initialize *parse_dispatch* table for parsing XML nodes:
	parse_dispatch = {}
	parse_dispatch["Block"] = XML_Block.parse
	parse_dispatch["Chamfers"] = XML_Chamfers.parse
	parse_dispatch["CNC_Flush"] = XML_CNC_Flush.parse
	parse_dispatch["Contour"] = XML_Contour.parse
	parse_dispatch["Corner"] = XML_Corner.parse
	parse_dispatch["DXF_PLace"] = XML_DXF_Place.parse
	parse_dispatch["Extra"] = XML_Extra.parse
	parse_dispatch["EZCAD"] = XML_EZCAD.parse
	parse_dispatch["Hole"] = XML_Hole.parse
	parse_dispatch["Hole_Through"] = XML_Hole_Through.parse
	parse_dispatch["Extrusion"] = XML_Extrusion.parse
	parse_dispatch["Part"] = XML_Part.parse
	parse_dispatch["Place"] = XML_Place.parse
	parse_dispatch["Simple_Pocket"] = XML_Simple_Pocket.parse
	parse_dispatch["Tool_Prefer"] = XML_Tool_Prefer.parse
	parse_dispatch["Tooling_Hole"] = XML_Tooling_Hole.parse
	parse_dispatch["Tooling_Plate"] = XML_Tooling_Plate.parse
	parse_dispatch["Tooling_Plate_Mount"] = XML_Tooling_Plate_Mount.parse
	parse_dispatch["Tube"] = XML_Tube.parse
	parse_dispatch["Vertical_Lathe"] = XML_Vertical_Lathe.parse
	parse_dispatch["Vice_Position"] = XML_Vice_Position.parse

	# Load up *self*:
	self.parse_dispatch = parse_dispatch	# Parse dispatch table
	self.extra_corner1 = P()		# BSW corner of extra material
	self.extra_corner2 = P()		# TNE corner of extra material

    def extra_set(self, extra_corner1, extra_corner2):
	""" XML_Root: Set the extra material corners to *extra_corner1* and
	    *extra_corner2*.
	"""

	# Verify argument types:
	assert isinstance(extra_corner1, P)
	assert isinstance(extra_corner2, P)

	# Load into *self*:
	self.extra_corner1 = extra_corner1
	self.extra_corner2 = extra_corner2

    def parse(self, element):
	""" *XML_Root*: Parse *element* into an *XML* object and return it.
	"""

	# Verify argument types:
	assert isinstance(element, ET.Element)

	element_tag = element.tag
	xml = None
	parse_dispatch = self.parse_dispatch
	if element_tag in parse_dispatch:
            xml = parse_dispatch[element_tag](element)
	    for child_element in element:
		child_xml = self.parse(child_element)
		assert isinstance(child_xml, XML), \
		  "<{0}> on line {1} has problems". \
		  format(child_element.tag, child_element.attrib["LN"])
		xml.xmls.append(child_xml)
	else:
	    print("No parser for <{0} ...> tag".format(element_tag))
	return xml

    def process(self, xml_ezcad):
	""" *XML_Root*: Process *xml_ezcad* to generate G-Codes, etc.
	"""

	# Verity argument types:
	assert isinstance(xml_ezcad, XML_EZCAD)

	# Process *xml_ezcad*:
	xml_ezcad.process(self)

    def xml_file_parse(self, xml_file_name):
	""" *XML_Root*: Parse *xml_file_name* into an *XML_EZCAD* object
	    and return it.
	"""

	xml_stream = open(xml_file_name, "ra")
	assert isinstance(xml_stream, file), \
	  "Unable to open XML file '{0}'\n".format(xml_file_name)

	# Read *xml_stream* into *xml_lines*:
	xml_lines = xml_stream.readlines()
	xml_stream.close()

	# For error reporting, it is nice to know the line number of where
	# the error occurred.  Unfortunately, the ElementTree package does
	# not store parse position information in the Element object.  To
	# work around this situation, we sweep through each line that was
	# read in and append a line number attribute appropriate tags.
	# Thus, '<Module ...' becomes "<Module LN="#" ...' where # is the
	# actual line number.  When there is an error, the appropriate
	# line number is read out using element.attrib["LN"].  This is
	# wrapped up in the line_number() function, just in the unlikely
	# case a tag line number is missed.

	# This pattern will match the tag at the beginning of the line:
	pattern = re.compile("^[ \t]*<\w+")

	# Scan through each line:
	replaced_lines = []
	for index in range(len(xml_lines)):
	    # Fetch each line one at a time:
	    line = xml_lines[index]

	    # Do we have the tag at the beginning of the line?
	    match = pattern.search(line)
	    if match:
		# Yes!  It can only match once, so grab the tag text:
		tag_text = match.group(0)

		# Append the line number attribute to the tag text:
		tag_text += ' LN="{0}"'.format(index + 1)

		# Substitute it back into the line:
		line = pattern.sub(tag_text, line)
		#print "{0}\t{1}".format(line_number, line)

	    # Build up the replaced lines list:
	    replaced_lines.append(line)
	
	# Glue *replaced_lines*  back together into a single *xml_text*:
	xml_text = '\n'.join(replaced_lines)

	# Now parse it into *root_element*:
	root_element = None
	try:
	    root_element = ET.fromstring(xml_text)
	except ET.ParseError as e:
	    print("'{0}': {1}".format(xml_file_name, e))

	# Parse *root_element* into *xml_ezcad* and return it:
	xml_ezcad = self.parse(root_element)
	return xml_ezcad

    def xml_write(self, xml_ezcad, xml_file_name):
	""" *XML_Part*: Write *xml_ezcad* to the file named *xml_file_name*.
	"""

	# Verify argument types:
	assert isinstance(xml_ezcad, XML_EZCAD)

	# Open *xml_stream* to write to *xml_file_name*:
	xml_stream = open(xml_file_name, "wa")
	assert isinstance(xml_stream, file), \
	  "Unable to open file '{0}'".format(xml_file_name)

	# Write out *xml_ezcad* to *xml_stream* and close it:
	xml_ezcad.xml_write(0, xml_stream)
	xml_stream.close()

class XML_Simple_Pocket(XML):
    def __init__(self, corner1, corner2, radius, flags, comment):
	# Verify argument types:
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)
	assert isinstance(radius, L)
	assert isinstance(flags, str)
	assert isinstance(comment, str)

	# Load up *self*:
	XML.__init__(self, "Simple_Pocket")
	self.corner1 = corner1
	self.corner2 = corner2
	self.radius = radius
	self.flags = flags
	self.comment = comment

    @staticmethod
    def parse(simple_pocket_element):
	# Verify argument types:
	assert isinstance(simple_pocket_element, ET.Element)
	assert simple_pocket_element.tag == "Simple_Pocket"

	# Parse *attributes* without error detection:
	attributes = simple_pocket_element.attrib
	c1x = L(float(attributes["C1X"]))
	c1y = L(float(attributes["C1Y"]))
	c1z = L(float(attributes["C1Z"]))
	c2x = L(float(attributes["C2X"]))
	c2y = L(float(attributes["C2Y"]))
	c2z = L(float(attributes["C2Z"]))
	radius = L(float(attributes["Radius"]))
	flags = attributes["Flags"]
	comment = attributes["Comment"]

	corner1 = P(c1x, c1y, c1z)
	corner2 = P(c2x, c2y, c2z)
	simple_pocket = \
	  XML_Simple_Pocket(corner1, corner2, radius, flags, comment)

	#	else_if equal@(tag_name, "Simple_Pocket")
	#	    # <Simple_Pocket C1X= C1Y= C1Z= C2X= C2Y= C2Z= ...
	#	    #  ... Radius= Flags= Comment= />:
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@(form@("<Simple_Pocket Comment=%v% ... />\n\") /
	#		  f@(comment))
	#
	#	    one :@= in@(1.0)
	#	    c1 :@= create@Point(c1x, c1y, c1z)
	#	    c2 :@= create@Point(c2x, c2y, c2z)
	#	    #call d@(form@("Simple_Pocket: Before: c1=%p% c2=%p%\n\") %
	#	    #  f@(c1) / f@(c2))
	#	    position :@= part.position
	#	    call matrix_apply@(c1, position)
	#	    call matrix_apply@(c2, position)
	#	    #call d@(form@("Simple_Pocket: After: c1=%p% c2=%p%\n\") %
	#	    #  f@(c1) / f@(c2))
	#
	#	    z_start :@= c1.z
	#	    z_end :@= c2.z
	#	    if z_start < z_end
	#		temporary :@= z_start
	#		z_start := z_end
	#		z_end := temporary
	#	    #call d@(form@("Simple_Pocket: z_start=%i% z_end=%i%\n\") %
	#	    #  f@(z_start) / f@(z_end))
	#
	#	    lower_chamfer :@= zero
	#	    upper_chamfer :@= zero
	#	    pocket_kind :@= flat@Pocket_Kind
	#	    flags_size :@= flags.size
	#	    flags_index :@= 0
	#	    while flags_index < flags_size
	#		flag :@= flags[flags_index]
	#		switch flag
	#		  case 'f'
	#		    pocket_kind := flat@Pocket_Kind
	#		  case 'm'
	#		    pocket_kind := flat@Pocket_Kind
	#		  case 'u'
	#		    assert 0f
	#		  case 't'
	#		    pocket_kind := through@Pocket_Kind
	#		    z_end := z_end - in@(0.020)
	#		  default
	#		    assert 0f
	#		flags_index := flags_index + 1
	#	    
	#	    call simple_pocket@(part, comment,
	#	      c1.x, c1.y, c2.x, c2.y, z_start, z_end, radius, pocket_kind)
	return simple_pocket

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	corner1 = self.corner1
	corner2 = self.corner2
	radius = self.radius
	flags = self.flags
	comment = self.comment
	xml_stream.write(
	  ('{0}<Simple_Pocket C1X="{1:.6m}" C1Y="{2:.6m}" C1Z="{3:.6m}" ' +
	  'C2X="{4:.6m}" C2Y="{5:.6m}" C2Z="{6:.6m}" ' +
	  'Radius="{7:.6m}" Flags="{8}" Comment="{9}"/>\n').
	  format(' ' * indent, corner1.x, corner1.y, corner1.z,
	  corner2.x, corner2.y, corner2.z, radius, flags, comment))

class XML_Tool_Prefer(XML):
    def __init__(self, tool_name):
	# Verify argument types:
	assert isinstance(tool_name, str)

	# Load up self:
	XML.__init__(self, "Tool_Prefer")
	self.tool_name = tool_name

    @staticmethod
    def parse(tool_prefer_element):
	# Verify argument types:
	assert isinstance(tool_prefer_element, ET.Element)
	assert tool_prefer_element.tag == "Tool_Prefer"

	# Parse *attributes* without error detection:
	attributes = tool_prefer_element.attrib
	tool_name = attributes["Tool_Name"]
	
	tool_prefer = XML_Tool_Prefer(tool_name)
	return tool_prefer

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	tool_name = self.tool_name
	xml_stream.write('{0}<Tool_Prefer Tool_Name="{1}"/>\n'.
	  format(' ' * indent, tool_name))

class XML_Tooling_Hole(XML):
    def __init__(self, row, column, adjust_x, adjust_y, flags):
	assert isinstance(row, int)
	assert isinstance(column, int)
	assert isinstance(adjust_x, int)
	assert isinstance(adjust_y, int)
	assert isinstance(flags, str)
	XML.__init__(self, "Tooling_Hole")
	self.row = row
	self.column = column
	self.adjust_x = adjust_x
	self.adjust_y = adjust_y
	self.flags = flags

    @staticmethod
    def parse(tooling_hole_element):
	assert isinstance(tooling_hole_element, ET.Element)
	assert tooling_hole_element.tag == "Tooling_Hole"
	attributes = tooling_hole_element.attrib
	row = int(attributes["Row"])
	column = int(attributes["Column"])
	adjust_x = int(attributes["Adjust_X"])
	adjust_y = int(attributes["Adjust_Y"])
	flags = attributes["Flags"]
	tooling_hole = XML_Tooling_Hole(row, column, adjust_x, adjust_y, flags)

	return tooling_hole

	#	else_if equal@(tag_name, "Tooling_Hole")
	#	    # <Tooling_Hole Row= Columns= Adjust_X= Adjust_Y= Flags= >
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@("<Tooling_Hole ... />\n\")
	#
	#	    tooling_hole :@= fetch2@(tooling_plate, column, row)
	#	    tooling_hole.adjust_x := adjust_x
	#	    tooling_hole.adjust_y := adjust_y
	#	    tooling_hole.deleted := equal@(flags, "x")
	#	    tooling_hole.ignored := equal@(flags, "i")

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	row = self.row
	column = self.column
	adjust_x = self.adjust_x
	adjust_y = self.adjust_y
	flags = self.flags
	xml_stream.write(
	  ('{0}<Tooling_Hole Row="{1}" Column="{2}" ' +
	  'Adjust_X="{3}" Adjust_Y="{4}" Flags="{5}"/>\n').
	  format(' ' * indent, row, column, adjust_x, adjust_y, flags))

class XML_Tooling_Plate(XML):
    def __init__(self, rows, columns, comment, ):
	assert isinstance(rows, int)
	assert isinstance(columns, int)
	assert isinstance(comment, str)
	XML.__init__(self, "Tooling_Plate")
	self.rows = rows
	self.columns = columns
	self.comment = comment
	self.tooling_holes = []

    @staticmethod
    def parse(tooling_plate_element):
	assert isinstance(tooling_plate_element, ET.Element)
	assert tooling_plate_element.tag == "Tooling_Plate"
	attributes = tooling_plate_element.attrib
	rows = int(attributes["Rows"])
	columns = int(attributes["Columns"])
	comment = attributes["Comment"]
	tooling_plate = XML_Tooling_Plate(rows, columns, comment)

	for child_element in tooling_plate_element:
	    child_tag = child_element.tag
	    if child_tag == "Tooling_Hole":
		tooling_hole = XML_Tooling_Hole.parse(child_element)
		tooling_plate.tooling_holes.append(tooling_hole)
	    else:
		assert "'{0}' not allowed under <Tooling_Plate> ... </>".\
		  format(child_tag)

	return tooling_plate

	#	else_if equal@(tag_name, "Tooling_Plate")
	#	    # <Tooling_Plate Rows= Columns= >
	#	    assert !xml_tag_end@(xml_stream, 0f)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@(form@("<Tooling_Plate Comment=%v% ... >/\n\") /
	#		  f@(comment))
	#
	#	    tooling_plate := create@Tooling_Plate(columns, rows)
	#	else_if equal@(tag_name, "/Tooling_Plate")
	#	    # </Tooling_Plate>
	#	    assert !xml_tag_end@(xml_stream, 0f)
	#	    call line_number_increment@(ezcad)
	#
	#	    debug :@= 0f
	#	    #debug := equal@(part.name, "cam")
	#	    if trace || debug
	#		call d@("</Tooling_Plate/>\n\")
	#	    if debug
	#		call d@(form@("Part=%v% Rows=%d% Columns=%d%\n\") %
	#		  f@(part.name) % f@(tooling_plate.rows_size) /
	#		  f@(tooling_plate.columns_size))
	#
	#	    # Compute the locations of the bounding box after the part
	#	    # has been mounted in the vice:
	#	    position :@= part.position
	#	    bounding_box1 := copy@(part.bounding_box1)
	#	    bounding_box2 := copy@(part.bounding_box2)
	#	    call matrix_apply@Point(bounding_box1, position)
	#	    call matrix_apply@Point(bounding_box2, position)
	#	    if debug
	#		call d@(form@("bb1=%p% bb2=%p%\n\") %
	#		  f@(bounding_box1) / f@(bounding_box2))
	#		call d@(form@("position matrix=\n\%m%\n\") / f@(position))
	#
	#	    # We only care about how big the bounding box is in X and Y:
	#	    bounding_box_x1 :@= bounding_box1.x
	#	    bounding_box_x2 :@= bounding_box2.x
	#	    bounding_box_dx :@= absolute@(bounding_box_x2 - bounding_box_x1)
	#	    bounding_box_y1 :@= bounding_box1.y
	#	    bounding_box_y2 :@= bounding_box2.y
	#	    bounding_box_dy :@= absolute@(bounding_box_y2 - bounding_box_y1)
	#	    bounding_box_dz :@= absolute@(bounding_box2.z - bounding_box1.z)
	#	    bounding_box_x_center :@= half@(bounding_box_x1 + bounding_box_x2)
	#	    bounding_box_y_center :@= half@(bounding_box_y1 + bounding_box_y2)
	#	    if debug
	#		call d@(form@("bb_dx=%i% bb_dy=%i% bb_dz=%i%\n\") %
	#		  f@(bounding_box_dx) % f@(bounding_box_dy) /
	#		  f@(bounding_box_dz))
	#
	#	    #FIXME: Hardwiring in the tool plate specifications!!!
	#	    plate_columns :@= 19
	#	    plate_rows :@= 7
	#	    plate_column_pitch :@= in@(0.5)
	#	    plate_row_pitch :@= in@(0.5)
	#	    plate_hole_diameter :@= in@(0.1065)
	#	    plate_column_edge :@= in@(0.5)
	#	    plate_row_edge :@= in@(0.25)
	#
	#	    # Figure out the initial number of needed rows and columns on the
	#	    # tooling plate.  Offset from the edges by {plate_hole_diameter}:
	#	    adjusted_dx :@= bounding_box_dx - twice@(plate_hole_diameter)
	#	    adjusted_dy :@= bounding_box_dy - twice@(plate_hole_diameter)
	#
	#	    # If {part} is too narrow, the assertions will fail:
	#	    if adjusted_dx <= in@(0.0) || adjusted_dy <= in@(0.0)
	#		call d@(form@("Part %v% failed (adj_dx=%i% adj_dy=%i%\n\") %
	#		  f@(part.name) % f@(adjusted_dx) / f@(adjusted_dy))
	#		call d@(form@("bb_dx=%i% bb_dy=%i%\n\") %
	#		  f@(bounding_box_dx) / f@(bounding_box_dy))
	#		assert 0f
	#
	#	    # Now figure out number of number of plate rows and columns needed:
	#	    plate_columns_needed :@=
	#	      unsigned@(div@(adjusted_dx, plate_column_pitch)) + 1
	#	    plate_rows_needed :@=
	#	      unsigned@(div@(adjusted_dy,  plate_row_pitch)) + 1
	#	    if debug
	#		call d@(form@(
	#		  "plate_cols_needed=%d% plate_rows_needed=%i%\n\") %
	#		  f@(plate_columns_needed) / f@(plate_rows_needed))
	#
	#	    # Figure out where the requested columns land on the tooling plate:
	#	    columns_size :@= tooling_plate.columns_size
	#	    columns_delta :@=
	#	      double@(plate_columns_needed - 1) / double@(columns_size - 1)
	#	    rows_size :@= tooling_plate.rows_size
	#	    rows_delta :@=
	#	      double@(plate_rows_needed - 1) / double@(rows_size - 1)
	#	    if debug
	#		call d@(form@("columns_delta=%i% rows_delta=%i%\n\") %
	#		  f@(columns_delta) / f@(rows_delta))
	#
	#	    # Now figure out where all of the {column_positions} will be:
	#	    column_positions :@= new@Array[Unsigned]()
	#	    columns_index :@= 0
	#	    while columns_index < columns_size
	#		# The 0.5 causes {column_position} to round to the closet col:
	#		column_position :@=
	#		  unsigned@(double@(columns_index) * columns_delta + 0.5)
	#		call append@(column_positions, column_position)
	#		if debug
	#		    call d@(form@("column_position[%d%]: %d%\n\") %
	#		      f@(columns_index) / f@(column_position))
	#		columns_index := columns_index + 1
	#
	#	    # Now figure out where all of the {row_positions} will be:
	#	    row_positions :@= new@Array[Unsigned]()
	#	    rows_index :@= 0
	#	    while rows_index < rows_size
	#		# The 0.5 causes {row_position} to round to the closet row:
	#		row_position :@= unsigned@(double@(rows_index) * rows_delta)
	#		call append@(row_positions, row_position)
	#		if debug
	#		    call d@(form@("row_position[%d%]: %d%\n\") %
	#		      f@(rows_index) / f@(row_position))
	#		rows_index := rows_index + 1
	#
	#	    # Reversing {row_positions} deals with the fact that
	#	    # {Part.tool_plate} numbers rows from 1 going down.
	#	    # In reality, row N-1 should be the top row:
	#	    call reverse@(row_positions)
	#
	#	    # Deterimine the number of rows and columns needed.  For
	#	    # this computation, we do not consider either ignored or
	#	    # deleted holes:
	#	    column_maximum :@= -123456789i
	#	    column_minimum :@= 123456789i
	#	    row_maximum :@= -123456789i
	#	    row_minimum :@= 123456789i
	#	    rows_index := 0
	#	    while rows_index < rows_size
	#		plate_row :@= row_positions[rows_index]
	#		columns_index := 0
	#		while columns_index < columns_size
	#		    # Fetch the appropropriate {tool_hole}:
	#		    plate_column :@= column_positions[columns_index]
	#		    tool_hole :@=
	#		      fetch2@(tooling_plate, columns_index, rows_index)
	#		    if debug
	#			call d@(form@(
	#			  "tool_hole[%d%,%d%]: ax=%d% ay=%i% d=%l% i=%l%\n\") %
	#			  f@(columns_index) % f@(rows_index) %
	#			  f@(tool_hole.adjust_x) % f@(tool_hole.adjust_y) %
	#			  f@(tool_hole.deleted) / f@(tool_hole.ignored))
	#
	#		    # Skip over both ignored and deleted holes:
	#		    if !tool_hole.ignored && !tool_hole.deleted
	#			# Adjust the row and column for this hole:
	#			adjusted_column :@=
	#			  integer@(plate_column) + tool_hole.adjust_x
	#			adjusted_row :@=
	#			  integer@(plate_row) + tool_hole.adjust_y
	#			tool_hole.adjusted_column := adjusted_column
	#			tool_hole.adjusted_row := adjusted_row
	#			if debug
	#			    call d@(form@(
	#			      "tool_hole[%d%,%d%]: ac=%d% ar=%i%\n\") %
	#			      f@(columns_index) % f@(rows_index) %
	#			      f@(adjusted_column) / f@(adjusted_row))
	#
	#			# Compute the minimum/maximum row/column:
	#			if adjusted_column > column_maximum
	#			    column_maximum := adjusted_column
	#			if adjusted_column < column_minimum
	#			    column_minimum := adjusted_column
	#			if adjusted_row > row_maximum
	#			    row_maximum := adjusted_row
	#			if adjusted_row < row_minimum
	#			    row_minimum := adjusted_row
	#		    columns_index := columns_index + 1
	#		rows_index := rows_index + 1
	#	    if debug
	#		call d@(form@(
	#		  "min_col=%d% max_col=%d% min_row=%i% max_col=%i%\n\") %
	#		  f@(column_minimum) % f@(column_maximum) %
	#		  f@(row_minimum) / f@(row_maximum))
	#
	#	    # Compute ({x_offset},{y_offset}) which corresponds to the
	#	    # hole at row = 0 and column = 0.  This is done by offseting
	#	    # from ({bounding_box_x_center},{bounding_box_y_center}) by
	#	    # an amount equal to half the needed rows and columns:
	#	    x_offset :@=
	#	      bounding_box_x_center - smul@(half@(plate_column_pitch),
	#	      double@(plate_columns_needed - 1))
	#	    y_offset :@=
	#	      bounding_box_y_center - smul@(half@(plate_row_pitch),
	#	      double@(plate_rows_needed - 1))
	#	    if debug
	#		call d@(form@("x_offset=%i% y_offset=%i%\n\") %
	#		  f@(x_offset) / f@(y_offset))
	#
	#	    # Now offset ({x_offset},{y_offset}) by the amount that
	#	    # the holes moved "in" by:
	#	    x_offset_dx :@= smul@(half@(plate_column_pitch),
	#	      double@(integer@(plate_columns_needed - 1) -
	#	      (column_minimum + column_maximum)))
	#	    y_offset_dy :@= smul@(half@(plate_row_pitch),
	#	      double@(integer@(plate_rows_needed - 1) -
	#	      (row_minimum + row_maximum)))
	#	    if debug
	#		call d@(form@("x_offset_dx=%i% y_offset_dy=%i%\n\") %
	#		  f@(x_offset_dx) / f@(y_offset_dy))
	#	    x_offset := x_offset + x_offset_dx
	#	    y_offset := y_offset + y_offset_dy
	#
	#	    # Now recompute the total number of rows and columns needed for
	#	    plate_columns_needed :=
	#	      unsigned@(column_maximum - column_minimum) + 1
	#	    plate_rows_needed := unsigned@(row_maximum - row_minimum) + 1
	#
	#	    # Make sure that the tooling plate can handle the holes:
	#	    if plate_columns_needed > plate_columns
	#		plate_columns_needed := plate_columns
	#	    if debug
	#		call d@(form@("plate_rows_needed=%i% plate_rows=%i%\n\") %
	#		  f@(plate_rows_needed) / f@(plate_rows))
	#	    if plate_rows_needed > plate_rows
	#		plate_rows_needed := plate_rows
	#	    if debug
	#		call d@(form@(
	#		  "plate_cols_needed=%d% plate_rows_needed=%i%\n\") %
	#		  f@(plate_columns_needed) / f@(plate_rows_needed))
	#
	#	    # Perform all of the hole drilling:
	#	    rows_index := 0
	#	    while rows_index < rows_size
	#		columns_index := 0
	#		while columns_index < columns_size
	#		    # Make sure only allowed holes get drilled:
	#		    tool_hole :@=
	#		      fetch2@(tooling_plate, columns_index, rows_index)
	#		    if !tool_hole.deleted
	#		        # We have a tooling plate hole to place:
	#			hole_comment :@= read_only_copy@(form@(
	#			  "Tooling Hole[%d%,%d%]") %
	#			   f@(columns_index) / f@(rows_index))
	#
	#			# Compute ({hole_x},{hole_y}) for each hole:
	#			hole_x :@= x_offset + smul@(plate_column_pitch,
	#			  double@(tool_hole.adjusted_column))
	#			hole_y :@= y_offset + smul@(plate_row_pitch,
	#			  double@(tool_hole.adjusted_row))
	#			if debug
	#			    call d@(form@("Hole[%d%,%d%]=(%i%,%i%)\n\") %
	#			      f@(columns_index) % f@(rows_index) %
	#			      f@(hole_x) / f@(hole_y))
	#
	#			# Drill the hole:
	#			call countersink_hole@(part, hole_comment,
	#			  plate_hole_diameter,
	#			  smul@(plate_hole_diameter, 1.150),
	#			  hole_x, hole_y,
	#			  in@(0.0), -bounding_box_dz, through@Hole_Kind)
	#		    columns_index := columns_index + 1
	#		rows_index := rows_index + 1
	#
	#	    part.tooling_plate_y := y_offset + plate_row_edge +
	#	      smul@(plate_column_pitch, double@(row_maximum))
	#

	#	else_if equal@(tag_name, "Tooling_Plate")
	#	    # <Tooling_Plate Rows= Columns= >
	#	    assert !xml_tag_end@(xml_stream, 0f)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@(form@("<Tooling_Plate Comment=%v% ... >/\n\") /
	#		  f@(comment))
	#
	#	    tooling_plate := create@Tooling_Plate(columns, rows)

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	rows = self.rows
	columns = self.columns
	comment = self.comment
	xml_stream.write(
	  '{0}<Tooling_Plate Rows="{1}" Columns="{2}" Comment="{3}">\n'.
	  format(' ' * indent, rows, columns, comment))
	for child_xml in self.xmls:
	    child_xml.xml_write(indent + 1, xml_stream)
	xml_stream.write('{0}</Tooling_Plate>\n'.format(' ' * indent))

class XML_Tooling_Plate_Mount(XML):
    def __init__(self, comment):
	assert isinstance(comment, str)
	XML.__init__(self, "Tooling_Plate_Mount")
	self.comment = comment

    @staticmethod
    def parse(tooling_plate_mount_element):
	assert isinstance(tooling_plate_mount_element, ET.Element)
	assert tooling_plate_mount_element.tag == "Tooling_Plate_Mount"
	attributes = tooling_plate_mount_element.attrib
	comment = attributes["Comment"]
	tooling_plate_mount = XML_Tooling_Plate_Mount(comment)
	return tooling_plate_mount

	#	else_if equal@(tag_name, "Tooling_Plate_Mount")
	#	    # <Tooling_Plate Comment= />
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@(form@("<Tooling_Plate_Mount Comment=%v% ... >/\n\") /
	#		  f@(comment))
	#
	#	    call reposition@(part, part.tooling_plate_y)
	#
	#	    position :@= part.position
	#	    extra1 :@= copy@(part.extra1)
	#	    extra2 :@= copy@(part.extra2)
	#	    call matrix_apply@(extra1, position)
	#	    call matrix_apply@(extra2, position)
	#	    extra_x :@= minimum@(extra1.x, extra2.x)
	#	    if trace
	#		call d@(form@("xml_read:%v%, extra1.x=%i% extra2.x=%i%\n\") %
	#		  f@(part.name) % f@(extra1.x) / f@(extra2.x))
	#		call d@(form@("dowel_position_set@(%v%, %i%, 0.0)\n\") %
	#		  f@(part.name) / f@(extra_x))
	#	    call dowel_position_set@(part, extra_x, zero)
	#	    call dowel_pin@(part, comment)

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	comment = self.comment
	xml_stream.write('{0}<Tooling_Plate_Mount Comment="{1}"/>\n'.
	  format(' ' * indent, comment))

class XML_Tube(XML):
    def __init(self, s, e, outer_diameter, wall_thickness, sides, color,
      transparency, material, comment):
	assert isinstance(s, P)
	assert isinstance(e, P)
	assert isinstance(outer_diameter, L)
	assert isinstance(wall_thickness, L)
	assert isinstance(sides, int)
	assert isinstance(color, str)
	assert isinstance(transparency, float)
	assert isinstance(material, str)
	assert isinstance(comment, str)
	XML.__init__(self, "Tube")

    @staticmethod
    def parse(tube_element):
	assert isinstance(tube_element, ET.Element)
	assert tube_element.tag == "Tube"
	attributes = tube_element.attrib
	sx = L(float(attributes["SX"]))
	sy = L(float(attributes["SY"]))
	sz = L(float(attributes["SZ"]))
	ex = L(float(attributes["EX"]))
	ey = L(float(attributes["EY"]))
	ez = L(float(attributes["EZ"]))
	outer_diameter = L(float(attributes["Outer_Diameter"]))
	wall_thickness = L(float(attributes["Wall_Thickness"]))
	sides = int(attributes["Sides"])
	color = attributes["Color"]
	transparency = float(attributes["Transparency"])
	material = attributes["Material"]
	comment = attributes["Comment"]

	s = P(sx, sy, sz)
	e = P(ex, ey, ez)
	tube = XML_Tube(s, e, outer_diameter, wall_thickness, sides, color,
	  transparency, material, comment)
	return tube

	#	else_if equal@(tag_name, "Tube")
	#	    # <Tube SX= SY= SZ= EX= EY= EZ= ...
	#	    # ... Outer_Diameter= Wall_Thickness= Sides= ...
	#	    # ... Color= Transparency= Material= Comment= />
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@(form@("<Tube Comment=%v% ... >/\n\") / f@(comment))
	#
	#	    part_material :@= create@Material(material, aluminum@Named_Material)
	#	    part := oriented_tube_create@Part(shop, name, part_material,
	#	      sx, sy, sz, ex, ey, ez, outer_diameter, wall_thickness,
	#	      sides, in@(0.0))
	#	    part.color := lookup@Color(color)
	#	    part.transparency := transparency

class XML_Vertical_Lathe(XML):
    def __init__(self, s, e, inner_diameter, outer_diameter, flags, comment):
	assert isinstance(s, P)
	assert isinstance(e, P)
	assert isinstance(inner_diameter, L)
	assert isinstance(outer_diameter, L)
	assert isinstance(flags, str)
	assert isinstance(comment, str)

	XML.__init__(self, "Vertical_Lathe")
	self.s = s
	self.e = e
	self.inner_diameter = inner_diameter
	self.outer_diameter = outer_diameter
	self.flags = flags
	self.comment = comment

    @staticmethod
    def parse(vertical_lathe_element):
	assert isinstance(vertical_lathe_element, ET.Element)
	assert vertical_lathe_element.tag == "Vertical_Lathe"
	attributes = vertical_lathe_element.attrib
	sx = L(float(attributes["SX"]))
	sy = L(float(attributes["SY"]))
	sz = L(float(attributes["SZ"]))
	ex = L(float(attributes["EX"]))
	ey = L(float(attributes["EY"]))
	ez = L(float(attributes["EZ"]))
	inner_diameter = L(float(attributes["Inner_Diameter"]))
	outer_diameter = L(float(attributes["Outer_Diameter"]))
	flags = attributes["Flags"]
	comment = attributes["Comment"]

	s = P(Sx, sy, sz)
	e = P(ex, ey, ez)
	veritcal_lathe = XML_Vertical_Lathe(s, e, inner_diameter,
	  outer_diameter, flags, comment)
	return vertical_lathe

	#	else_if equal@(tag_name, "Vertical_Lathe")
	#	    # <Simple_Pocket C1X= C1Y= C1Z= C2X= C2Y= C2Z= ...
	#	    #  ... Radius= Flags= Comment= />:
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@(form@("<Vertical_Lathe Comment%v% ... />\n\") /
	#		  f@(comment))
	#
	#	    one :@= in@(1.0)
	#	    start_point :@= create@Point(sx, sy, sz)
	#	    end_point :@= create@Point(ex, ey, ez)
	#	    position :@= part.position
	#	    call matrix_apply@(start_point, position)
	#	    call matrix_apply@(end_point, position)
	#
	#	    upper_chamfer :@= zero
	#	    flags_size :@= flags.size
	#	    flags_index :@= 0
	#	    while flags_index < flags_size
	#		flag :@= flags[flags_index]
	#		switch flag
	#		  case 'u'
	#		    assert 0f
	#		  case 'i'
	#		    outer_diameter := zero
	#		  default
	#		    assert 0f
	#		flags_index := flags_index + 1
	#	    
	#	    call vertical_lathe@(part, comment, inner_diameter, outer_diameter,
	#	      start_point.x, start_point.y, start_point.z, end_point.z,
	#	      upper_chamfer, 1.0, 0f)


class XML_Vice_Position(XML):
    def __init__(self, t, n, w, c, comment):
	assert isinstance(t, P)
	assert isinstance(n, P)
	assert isinstance(w, P)
	assert isinstance(c, P)
	assert isinstance(comment, str)
	XML.__init__(self, "vice_position")
	self.t = t
	self.n = n
	self.w = w
	self.c = c
	self.comment = comment
	
    @staticmethod
    def parse(vice_position_element):
	assert isinstance(vice_position_element, ET.Element)
	assert vice_position_element.tag == "Vice_Position"
	attributes = vice_position_element.attrib
	tx = L(float(attributes["TX"]))
	ty = L(float(attributes["TY"]))
	tz = L(float(attributes["TZ"]))
	nx = L(float(attributes["NX"]))
	ny = L(float(attributes["NY"]))
	nz = L(float(attributes["NZ"]))
	wx = L(float(attributes["WX"]))
	wy = L(float(attributes["WY"]))
	wz = L(float(attributes["WZ"]))
	cx = L(float(attributes["CX"]))
	cy = L(float(attributes["CY"]))
	cz = L(float(attributes["CZ"]))
	comment = attributes["Comment"]

	t = P(tx, ty, tz)
	n = P(nx, ny, nz)
	w = P(wx, wy, wz)
	c = P(cx, cy, cz)
	vice_position = XML_Vice_Position(t, n, w, c, comment)
	return vice_position

	#	else_if equal@(tag_name, "Tool_Prefer")
	#	    # <Tool_Prefer Tool_Name= />
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    if trace
	#		call d@(form@("<Tool_Prefer Tool_name=%v%>/\n\") /
	#		  f@(tool_name))
	#
	#	    if tool_name.size = 0
	#		tool_name := null@String
	#	    call tool_prefer@(part, tool_name)

	#	else_if equal@(tag_name, "Vice_Position")
	#	    # <Vice_Position TX= TY= TZ= NX= NY= NZ= WX= WY= WZ= CX= CY= CZ= />:
	#	    assert !xml_tag_end@(xml_stream, 1t)
	#	    call line_number_increment@(ezcad)
	#	    debug :@= 0t
	#	    #debug := 1t
	#	    #if equal@(part.name, "Tine_Tang")
	#	    #	debug := 1t
	#
	#	    if trace || debug
	#		call d@(form@("<Vice_Position Comment=%v% ... />\n\") /
	#		  f@(comment))
	#
	#	    # We need {zero} and {one} to define the Z axis below:
	#	    one :@= in@(1.0)
	#
	#	    # Load up {Point}'s {t}, {n}, {w}, and {c}:
	#	    t :@= create@Point(tx, ty, tz)
	#	    n :@= create@Point(nx, ny, nz)
	#	    w :@= create@Point(wx, wy, wz)
	#	    c :@= create@Point(cx, cy, cz)
	#
	#	    # Show what we have when debugging:
	#	    if debug
	#		call d@(form@(
	#		  "vice_pos0:\n,t\t=%p%\n,t\n=%p%\n,t\w=%p%\n,t\c=%p%\n\") %
	#		  f@(t) % f@(n) % f@(w) / f@(c))
	#
	#	    # Force the {position} matrix to identity.  This also forces
	#	    # the {reposition} matrix to identity as well:
	#	    call position_reset@(part)
	#	    position :@= part.position
	#	    reposition :@= part.reposition
	#	    matrix :@= part.shop.matrix
	#
	#	    # Normalize everything to be centered around {c}:
	#	    call translate@(part, -cx, -cy, -cz)
	#	    
	#	    # Move {t}, {n}, {w}, and {c} to their new homes:
	#	    call xyz_set@(t, tx, ty, tz)
	#	    call xyz_set@(n, nx, ny, nz)
	#	    call xyz_set@(w, wx, wy, wz)
	#	    call xyz_set@(c, cx, cy, cz)
	#	    call matrix_apply@Point(t, position)
	#	    call matrix_apply@Point(n, position)
	#	    call matrix_apply@Point(w, position)
	#	    call matrix_apply@Point(c, position)
	#
	#	    # Show what we have when debuging:
	#	    if debug
	#		call d@(form@(
	#		  "vice_pos1:\n,t\t=%p%\n,t\n=%p%\n,t\w=%p%\n,t\c=%p%\n\") %
	#		  f@(t) % f@(n) % f@(w) / f@(c))
	#
	#	    # Rotate the part so that the surface is pointing up.  Start
	#	    # by computing {top_angle}, the angle between ({tx},{ty},{tz})
	#	    # and the positive Z axis normal:
	#	    top_angle :@= angle_between@Point(t.x, t.y, t.z, zero, zero, one)
	#
	#	    if debug
	#		call d@(form@("top_angle=%d%\n\") / f@(top_angle))
	#
	#	    # See if we need to rotate {part} by {top_angle}:
	#	    if top_angle < -degrees@(0.0001) || top_angle > degrees@(0.0001)
	#		# We need to rotate {part} by {top_angle}:
	#		top_axis :@= null@Point
	#		if top_angle < degrees@(-179.9) || top_angle > degrees@(179.9)
	#		    # We have to entirely flip the board.  Unfortunately,
	#		    # a cross product will not work on colinear segments,
	#		    # so we use the ({n.x}, {n.y}, {n.z}) x (0, 0, 1) to
	#		    # establish the axis of rotation:
	#		    top_axis :=
	#		      cross_product@Point(n.x, n.y, n.z, zero, zero, one)
	#
	#		    # Force {top_angle} to 180 degrees, to get rid of any
	#		    # rounding errors:
	#		    top_angle := degrees@(180.0)
	#		else
	#		    # Have to partially flip the board. Use the cross product
	#		    # ({t.x}, {t.y}, {t.z}) x (0, 0, 1) to compute {top_axis},
	#		    # the axis about which to rotate the board:
	#		    top_axis :=
	#		      cross_product@Point(t.x, t.y, t.z, zero, zero, one)
	#
	#		# Show what we have during debugging:
	#		if debug
	#		    call d@(form@("top_axis=(%i%,%i%,%i%)\n\") %
	#		      f@(top_axis.x) % f@(top_axis.y) / f@(top_axis.z))
	#
	#		# Compute the normal of {top_axis}:
	#		call normalize@(top_axis)
	#		top_nx :@= in@(top_axis.x)
	#		top_ny :@= in@(top_axis.y)
	#		top_nz :@= in@(top_axis.z)
	#
	#		# Now rotate {part} by {top_angle_degrees} around the
	#		# normal ({top_nx}, {top_ny}, {top_nz}):
	#		call rotate@(part, top_nx, top_ny, top_nz, top_angle, zero)
	#
	#		# Now T=({tx},{ty},{tz}) surface is on top (i.e. parallel to
	#		# to Z axis and pointing up:
	#
	#		# Now that we have done the first transformation, we want to
	#		# work with T, N, and W in their new locations:
	#		call xyz_set@(t, tx, ty, tz)
	#		call xyz_set@(n, nx, ny, nz)
	#		call xyz_set@(w, wx, wy, wz)
	#		call xyz_set@(c, cx, cy, cz)
	#		call matrix_apply@Point(t, position)
	#		call matrix_apply@Point(n, position)
	#		call matrix_apply@Point(w, position)
	#		call matrix_apply@Point(c, position)
	#
	#		# Show what we have when debugging:
	#		if debug
	#		    call d@(form@(
	#		      "vice_pos2:\n,t\t=%p%\n,t\n=%p%\n,t\w=%p%\n,t\c=%p%\n\") %
	#		      f@(t) % f@(n) % f@(w) / f@(c))
	#
	#	    # Now shift {part} down so that {t} is on Z=0 plane:
	#	    call translate@(part, zero, zero, -t.z)
	#
	#	    # Like before, move T, N, W, and C to their new homes:
	#	    call xyz_set@(t, tx, ty, tz)
	#	    call xyz_set@(n, nx, ny, nz)
	#	    call xyz_set@(w, wx, wy, wz)
	#	    call xyz_set@(c, cx, cy, cz)
	#	    call matrix_apply@Point(t, position)
	#	    call matrix_apply@Point(n, position)
	#	    call matrix_apply@Point(w, position)
	#	    call matrix_apply@Point(c, position)
	#
	#	    # Show what we have during debugging:
	#	    if debug
	#		call d@(form@(
	#		  "vice_pos3:\n,t\t=%p%\n,t\n=%p%\n,t\w=%p%\n,t\c=%p%\n\") %
	#		  f@(t) % f@(n) % f@(w) / f@(c))
	#
	#	    # Now rotate {part} around positive Z axis so that N =
	#	    # ({nx},{ny},{nz}) is facing north:
	#
	#	    # Start by computing NT = N - T.
	#	    nt_x :@= n.x - t.x
	#	    nt_y :@= n.y - t.y
	#	    nt_z :@= n.z - t.z
	#	    vice_angle :@=
	#	      -angle_between@Point(zero, one, zero, nt_x, nt_y, nt_z) 
	#
	#	    # Compute {vice_angle} between positive Y-axis (i.e. north) and NT:
	#	    if debug
	#		call d@(form@("NT=(%i%,%i%,%i%)\n\") %
	#		  f@(nt_x) % f@(nt_y) / f@(nt_z))
	#		call d@(form@("vice_angle=%d%\n\") / f@(vice_angle))
	#
	#	    # Are we so close to 0 degrees, that we should not bother:
	#	    if vice_angle < -degrees@(0.1) || vice_angle > degrees@(0.1)
	#		# No, perform the rotation:
	#		nt :@= create@Point(nt_x, nt_y, nt_z)
	#		call rotate@(part, 0.0, 0.0, 1.0, vice_angle, length@(nt))
	#
	#		# As before, update T, N and W to their new locations:
	#		call xyz_set@(t, tx, ty, tz)
	#		call xyz_set@(n, nx, ny, nz)
	#		call xyz_set@(w, wx, wy, wz)
	#		call matrix_apply@Point(t, position)
	#		call matrix_apply@Point(n, position)
	#		call matrix_apply@Point(w, position)
	#
	#		# Show what we have during debugging:
	#		if debug
	#		    call d@(form@(
	#		      "vice_pos4:\n\t=%p%\n,t\n=%p%\n,t\w=%p%\n,t\c=%p%\n\") %
	#		      f@(t) % f@(n) % f@(w) / f@(c))
	#
	#	    # Start by computing NT = N - T and WT = W - T:
	#	    nt_x := n.x - t.x
	#	    nt_y := n.y - t.y
	#	    nt_z := n.z - t.z
	#	    wt_x :@= w.x - t.x
	#	    wt_y :@= w.y - t.y
	#	    wt_z :@= w.z - t.z
	#
	#	    # Compute {dowel_angle} NTW:
	#	    dowel_angle :@=
	#	      angle_between@Point(nt_x, nt_y, nt_z, wt_x, wt_y, wt_z)
	#
	#	    # Verify that that the angle <NTW is 90 degrees:
	#	    if dowel_angle < degrees@(89.9) || dowel_angle > degrees@(90.1)
	#		# It is not, let somebody know:
	#		call d@(form@("nt_x=%f% nt_y=%f% nt_z=%f%\n\") %
	#		  f@(nt_x) % f@(nt_y) / f@(nt_z))
	#		call d@(form@("wt_x=%f% wt_y=%f% wt_z=%f%\n\") %
	#		  f@(wt_x) % f@(wt_y) / f@(wt_z))
	#		call d@(form@(
	#		  "dowel_angle for part %v% is %d% (comment=%v%)\n\") %
	#		  f@(part.name) % f@(dowel_angle) / f@(comment))
	#
	#	    if debug
	#		bounding_box1 := copy@(part.bounding_box1)
	#		bounding_box2 := copy@(part.bounding_box2)
	#		call d@(form@("vice_pos: before: bb1=%p% bb2=%p%\n\") %
	#		  f@(bounding_box1) / f@(bounding_box2))
	#	        call matrix_apply@Point(bounding_box1, position)
	#		call matrix_apply@Point(bounding_box2, position)
	#		call d@(form@("vice_pos: after: bb1=%p% bb2=%p%\n\") %
	#		  f@(bounding_box1) / f@(bounding_box2))
	#
	#	    if debug
	#		call d@(form@("reposition(*, %i%)\n\") / f@(n.y))
	#	    call reposition@(part, n.y)
	#
	#	    #extra1 :@= copy@(part.extra1)
	#	    #extra2 :@= copy@(part.extra2)
	#	    #call matrix_apply@(extra1, position)
	#	    #call matrix_apply@(extra2, position)
	#
	#	    if debug
	#		call d@(form@("vice_pos:dowel_pos_set(%v%, %i%, 0.0)\n\") %
	#		  f@(part.name) / f@(w.x))
	#	    call dowel_position_set@(part, w.x, zero)
	#	    #extra_x :@= minimum@(extra1.x, extra2.x)
	#	    #call d@(form@("dowel_set(*, %i%, 0.0)\n\") / f@(extra_x))
	#	    #call dowel_position_set@(part, extra_x, zero)
	#	    call dowel_pin@(part, comment)

    def xml_write(self, indent, xml_stream):
	assert isinstance(indent, int)
	assert isinstance(xml_stream, file)
	t = self.t
	n = self.n
	w = self.w
	c = self.c
	comment = self.comment
	xml_stream.write(('{0}<Vice_Position ' +
	  'TX="{1:.6m}" TY="{2:.6m}" TZ="{3:.6m}" ' +
	  'NX="{4:.6m}" NY="{5:.6m}" NZ="{6:.6m}" ' +
	  'WX="{7:.6m}" WY="{8:.6m}" WZ="{9:.6m}" ' +
	  'CX="{10:.6m}" CY="{11:.6m}" CZ="{12:.6m}" ' +
	  'Comment="{13}"/>\n').format(' ' * indent,
	  t.x, t.y, t.z, n.x, n.y, n.z,
	  w.x, w.y, w.z, c.x, c.y, c.z, comment))
