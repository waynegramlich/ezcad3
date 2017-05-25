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

def float_compare(left, right):
    """ *int*: Return -1, 0, or 1 depending upon whether *left* sorts before,
	at, or after *right*.
    """

    assert isinstance(left, float)
    assert isinstance(right, float)
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
    #	 *mm* + *cm* + *inch* + *ft*.
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
    # of the form "W-N/D" where W is a whole number, N is a numerator, and
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
    #	 angle in the "correct" quadrant.
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

    def compare(self, length):
	""" *L* Return -1, 0, or 1 depending upon whether the *L* object (i.e. *self*)
	    is less than, equal to, or greater than *length*, repectively.
	"""

	# Verify argument types:
	assert isinstance(length, L)

	return float_compare(self._mm, length._mm)


    def cosine(self, angle):
	""" L: Return {self} * cos(angle). """

	assert isinstance(angle, Angle)
	return L(self._mm * angle.cosine())

    def distance(self, dy):
	""" L: Return the length of the line between (*self*,*dy) and the origin.  """
	dx_mm = self._mm
	dy_mm = dy._mm
	distance_mm = math.sqrt(dx_mm * dx_mm + dy_mm * dy_mm)
	return L(mm=distance_mm)

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

    def minimum_maximum(self, length2):
        """ *L*: Return two values with the minimum value coming first. """

	return self.minimum(length2), self.maximum(length2)

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

	if format == "N":
	    return "np.array([[{0}, {1}, {2}, 1.]])".format(
	       self.x.millimeters(), self.y.millimeters(), self.z.millimeters())

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
	  1.0	      ]] )

	return result

    def normalize(self):
	""" *P*: Return *self* normalized to have a length of 1. """

	x = self.x._mm
	y = self.y._mm
	z = self.z._mm
	length = math.sqrt(x * x + y * y + z * z)
	zero = L()
	if length == 0.0:
            normalized = self
	    if EZCAD3.update_count_get() == 0:
		print("P.normailize() was passed all zeros.")
        else:
            normalized = P(L(mm = x / length), L(mm = y / length), L(mm = z / length))

	#print("x={0} y={1} z={2} length={3} nomalized={4:m}".format(x, y, z, length, normalized))
	return normalized
	

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

    @staticmethod
    def polar(angle, distance):
	""" P: """

	assert isinstance(angle, Angle)
	assert isinstance(distance, L)

	x = distance.cosine(angle)
	y = distance.sine(angle)
	return P(x, y, L())

    def triple(self):
        """ *P*: Return the *P* object (i.e. *self*) as an immutable tuple of 3 floats
	    reprecented in millimeters.
	"""
	
	p = self
	return (p.x._mm, p.y._mm, p.z._mm)


    def twice(self):
	""" P dimensions: Return {self} * 2. """

	return self * 2.0

    def x_adjust(self, x):
	""" P dimensions: Return copy of {self} with {x} added to the
	    x field. """

	return P(self.part, self.x + x, self.y, self.z)

    def xy_angle(self):
	""" *P*: """

	return Angle(rad=math.atan2(self.y.millimeters(), self.x.millimeters()))

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
    #	 control formatting.
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

	radians = self.radians
	pi = math.pi
	double_pi = pi + pi
	while (radians > pi):
	    radians -= double_pi
	while (radians < -pi):
	    radians += double_pi
	return Angle(rad=radians)

    def minimum(self, angle2):
	""" *Angle*: Return the minimum of *self* and *angle2*"""

	assert isinstance(angle2, Angle)

	result = self
	if angle2.radians < result.radians:
	    result = angle2
	return result

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

    def xy_direction(self):
	""" *Angle*: Return the unit direction vector in the X/Y plane for the *Angle* object
	    (i.e. *self*).
	"""

	# Compute *x_cosine* and *y_cosine* rouinding values less than *epsilon* to zero:
	epsilon = 1.0e-15
	x_cosine = self.cosine()
	if abs(x_cosine) < epsilon:
	    x_cosine = 0.0
	y_sine = self.sine()
	if abs(y_sine) < epsilon:
	    y_sine = 0.0

	# Create *direction*:
	x = L(mm=x_cosine)
	y = L(mm=y_sine)
	z = L()
	direction = P(x, y, z)

	return direction

class Bend:

    def __init__(self, name, point, radius):
	""" *Bend*: Initialize the *Bend* object (i.e. *self*) to contain *name*,
	    *point*, and *radius*.  *Name* is used for debugging and shows up in
	    generated *G-code*.  *point* is the location of the bend as if it was
	    0 degrees in radius.  *radius* is the actual bend radius.
	"""

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(point, P)
	assert isinstance(radius, L)

	# Use *bend* instead of *self*:
	bend = self

	# Load up the *Bend* object (i.e. *self*) with th:
	bend._name = name
	bend._point = point
	bend._radius = radius

	# We do most of the work in a projected coordinates.  The *projected* point is
	# *point* after it has been rotated and translated so that the points project down
	# down to the X/Y plane.  This value is computed by the *_project*() routine.
	bend._projected_point = None	# Projected point:

	# These fields are computed as a side-effect by *_smallest_inner_radius_compute*():
	bend._is_inside = None		# *True* if *bend* is an inside bend; otherwise *False*
	bend._bearing_change = None	# Set to be the bearing change angle.

	# These fields are computed as a side-effect by *_center_radius_and_tangents_compute*():
	bend._center = None
	bend._incoming_tangent_angle = None
	bend._incoming_tangent_direction = None
	bend._outgoing_tangent_angle = None
	bend._outgoing_tangent_direction = None

    def __format__(self, format):
	""" *Bend*: Return the *Bend* object (i.e. *self*) formated as a string. """

	return "Bend(name='{0}' point={1}, radius={2})".format(
	  self._name, self._point, self._radius)

    def _bearing_change_get(self):
	""" *Bend*: Return the bearing_change field from the *Bend* object (i.e. *self*). """

	return self._bearing_change

    def _center_get(self):
	""" *Bend*: Return the radius center for *Bend* object (i.e. *self*). """

	return self._center

    def _copy(self):
        """ *Bend*: Return a copy of the *Bend* object (i.e. *self*). """

	# Start the copy
	bend = Bend(self._name, self._point, self._radius)

	# Copy the remaining fields:
	bend._projected_point = self._projected_point
        bend._is_inside = self._is_inside
        bend._bearing_change = self._bearing_change
        bend._center = self._center
        bend._incoming_tangent_angle = self._incoming_tangent_angle
        bend._incoming_tangent_direction = self._incoming_tangent_direction
        bend._outgoing_tangent_angle = self._outgoing_tangent_angle
        bend._outgoing_tangent_direction = self._outgoing_tangent_direction

	return bend

    def _incoming_tangent_angle_get(self):
	""" *Bend*: Return the incoming tangent angle for *Bend* object (i.e. *self*)
	    relative to the radius center. """
	
	return self._incoming_tangent_angle

    def _incoming_tangent_compute(self, radius):
	""" *Bend*: Compute and return the incoming tangent point for the *Bend* object
	    (i.e. *self*) with a radius of *radius*.
	"""

	# Verify argument types:
	assert isinstance(radius, L)

	# Compute *incoming_tangent*:
	center = self._center
	assert isinstance(center, P)
	incoming_tangent_direction = self._incoming_tangent_direction
	assert isinstance(incoming_tangent_direction, P)

	incoming_tangent = center + incoming_tangent_direction * radius.millimeters()
	return incoming_tangent

    def _is_inside_get(self):
	""" *Bend*: Set the is_inside field of the *Bend* object (i.e. *self*) to *is_inside. """

	# Store *is_inside* into the *Bend* object (i.e. *self*):
	return self._is_inside

    def _is_inside_set(self, is_inside):
	""" *Bend*: Set the is_inside field of the *Bend* object (i.e. *self*) to *is_inside. """

	# Verify argument types:
	assert isinstance(is_inside, bool)

	# Store *is_inside* into the *Bend* object (i.e. *self*):
	self._is_inside = is_inside

    def _name_get(self):
	""" *Bend*: Return name with the *Bend* object (i.e. *self*). """

	return self._name

    def _outgoing_tangent_angle_get(self):
	""" *Bend*: Return the outgoing tangent angle for *Bend* object (i.e. *self*)
	    relative to the radius center. """
	
	return self._outgoing_tangent_angle

    def _outgoing_tangent_compute(self, radius):
	""" *Bend*: Compute and return the outgoing tangent point for the *Bend* object
	    (i.e. *self*) with a radius of *radius*.
	"""

	# Verify argument types:
	assert isinstance(radius, L)

	# Compute *outgoing_tangent*
	center = self._center
	outgoing_tangent_direction = self._outgoing_tangent_direction
	outgoing_tangent = center + outgoing_tangent_direction * radius.millimeters()
	return outgoing_tangent

    def _radius_center_and_tangents_compute(self, incoming_bend, outgoing_bend, tracing):
	""" *Bend*: Compute the radius center and the radius tangent points for the
	    the *Bend* object (i.e. *self*) where *incoming_bend* is the *Bend* object
	    before and and *outgoing_bend* is the *Bend* object after.  The radius
	    circle center point is the where it is tangent to both the *incoming_bend*
	    and *outgoing_bend*.  The results can be accessed with the *_center_get*(),
	    *_incoming_tangent_get*(), and *_outgoing_tangent_get*().
	"""

	# Verify argument types:
	assert isinstance(incoming_bend, Bend)
	assert isinstance(outgoing_bend, Bend)
	assert isinstance(tracing, int)

	# Use *bend* instead of *self*:
	bend = self

	# Do any requested *tracing*:
	tracing_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    #tracing_detail = 3
	    print("{0}=>Bend._radius_center_and_tangents_compute('{1}', '{2}', '{3}')".format(
	    ' ' * tracing, bend._name, incoming_bend._name, outgoing_bend._name))

	# Below is an ASCII art picture of a *Bend* object (i.e. *self*).  B represents
	# the bend point.  There is circle of radius r centered around the circle center
	# point C.  r is the radius of the bend.  *I* and *O* are the adjecent *Bend* object
	# points.  The circle tangentally touches the the IB line segment at J and the OB
	# line segement at N.  This the length of both JC and NC is r.  The length of BC
	# is given the variable name of d.
	#
	#     
	#     I
	#      \
	#       \	  **|**
	#        \      **  |  **      O
	#	  \   *     |     *   /
	#	   \ *      |      * /
	#     -+--> \*      C      */ <------+-
	#      |     J      |      N	     ^
	#      r      *     |     *	     |
	#      |       \**  |  **/	    
	#     -+------> \ **|** /	     d
	#		 \  |  /	      
	#		  \ | /	             |
	#		   \|/	             V
	#		    B  <-------------+-
	# 
	# Our goal is to compute C using I, B, O and r.  Once we have that we need will also
	# need the direction vectors <<JC>> and <<NC>> for computing the tangent points J and N:
        #
        #        J = C + r * <<JC>>					(1)
	#
        #        N = C + r * <<NC>>					(2)
	#
	# Please note that this geometry is insensitive to whether the bend is an inside
	# bend or an outside bend.  The only time this fails is when the bend is 180 degrees
	# (i.e. no bend) at which point cos(90 - <CBE) evaluates to 0.  This is simply not
	# allowed with a contour.
	#
	# The ASCII art picture above is redrawn with just C, B, N and O shown.
	# Everything has been rotated clockwise by approximate 45-60 degrees.
	# The key improvement is showing that the segement NC is perpendicular
	# to segment NB.  The length of |BC| is d and the length of |CT| is r.
	# Now we can do some geometry and trigonmetry to derive an equation for d.
	#
	#
	#            C
	#           /|
	#          / |
	#  |BC|=d /  | r=|CT|
	#        / +-|
	#       /  | |
	#      B-----N-----O
	#
	# The law of sines lets us write:
	#
	#       d          r
	#    -------- = --------                              (1)
	#    sin <BNC   sin <CBN
	#
        #    <BNC = 90                                        (2)
	#
	#       d        r
	#    ------ = --------                                (3)
	#    sin 90   sin <CBN
	#
	#    d      r
	#    - = --------                                     (4)
	#    1   sin <CBN
	#
	#           r
	#    d = -------- = r / sin( <CBN)                    (5)
	#        sin <CBN
	# 
	#    <CBN = <CBO                                      (6)
	#
	#    d = r / sin(<CBO)                                (7)

	# Compute the IB and OB vectors:
	i = incoming_bend._projected_point
	b = bend._projected_point
	o = outgoing_bend._projected_point
	ib = i - b
	ob = o - b
	if tracing_detail >= 1:
	    indent = ' ' * tracing
	    print("{0}i={1:i} b={2:i} o={3:i}".format(indent, i, b, o))
	    print("{0}ib={1:i} ob={2:i}".format(indent, ib, ob))

	# Compute the angles of IB and OB from the origin:
	ib_angle = ib.xy_angle()
	ob_angle = ob.xy_angle()
	if tracing_detail >= 1:
	    print("{0}ib_angle={1:d} ob_angle={2:d}".format(indent, ib_angle, ob_angle))

	# Compute *ibo_angle* and *cbo_angle*:
	ibo_angle = ib.angle_between(ob)
	cbo_angle = ibo_angle / 2
	if tracing_detail >= 2:
	    print("{0}ibo_angle={1:d} cbo_angle={2:d}".format(indent, ibo_angle, cbo_angle))
	    
	# Compute *d* using *cbo_angle* and *r*:
	r = bend._radius
	d = r / cbo_angle.sine()
	d_mm = d.millimeters()
	if tracing_detail >= 2:
	    print("{0}r={1:i} d={2:i}".format(indent, r, d))

	# We need to find a some point on the line between C and B.
	#
	#         I
	#          \         
	#           \       Q       
	#            \     /|   
	#             \   / C     O
	#              \ /  |    /            
	#        <<IB>> +---M---+ <<OB>>
	#                \  |  /              
	#                 \ | /
	#                  \|/
	#                   B
	#
	# We do this by computing the direction vectors <<IB>> from B to I and <<OB>>
	# form B to O.  By simly adding <<IB>> and <<OB>> together we find another
	# point on the line through B and C.  We show it as Q above.  We normalize the
	# result to get <<CB>>:
	ib_direction = ib.normalize()
	ob_direction = ob.normalize()
	cb_direction = (ib_direction + ob_direction).normalize()

	# Now that we have <<CB>> and d, we can compute C:
	#
	#        C = B + d * <<CB>>
	cb_angle = cb_direction.xy_angle()
	c = b + cb_direction * d_mm
	if tracing_detail >= 3:
	    print("{0}<<ib>>={1:m} <<ob>={2:m} <<cb>>={3:m}". \
	      format(indent, ib_direction, ob_direction, cb_direction))
	if tracing_detail >= 2:
	    print("{0}cb_angle={1:d}".format(indent, cb_angle))
	if tracing_detail >= 1:
	    print("{0}c={1:i}".format(indent, c))

	# The next step is remarkably simple.  The direction vector of <<JC>> from C to J
	# is just 90 degress from the *ib_angle* because the circle is tangent at J.
	# Likewsise, for direction vector <<NC>>.
	degrees90 = Angle(deg=90.0)
	jc_angle = (ib_angle - degrees90).normalize()
	jc_direction = jc_angle.xy_direction()
	nc_angle = (ob_angle + degrees90).normalize()
	nc_direction = nc_angle.xy_direction()
	if tracing_detail >= 2:
	    print("{0}jc_angle={1:d} nc_angle={2:d}".format(indent, jc_angle, nc_angle))
	if tracing_detail >= 3:
	    print("{0}cos(jc_angle)={1} sin(jc_angle){2}".
	      format(indent, jc_angle.sine(), jc_angle.cosine()))
	    print("{0}jc_direction={1:m} nc_direction={2:m}".
	      format(indent, jc_direction, nc_direction))

	# Now we just want to record everything we care about into *bend*:
	if tracing >= 0:
	    print("{0}center={1:i}".format(indent, c))
	bend._center = c
	bend._incoming_tangent_angle = jc_angle
	bend._incoming_tangent_direction = jc_direction
	bend._outgoing_tangent_angle = nc_angle
	bend._outgoing_tangent_direction = nc_direction

	if tracing >= 0:
	    print("{0}<=Bend._radius_center_and_tangents_compute('{1}', '{2}', '{3}')".format(
	    indent, bend._name, incoming_bend._name, outgoing_bend._name))
	    if tracing_detail >= 1:
		print("")

    def _offset_compute(self, before, after, offset, end_mill_radius, tracing=-1000000):
	""" *Bend*: Compute the field values for the *Bend* object (i.e. *self*).
	    *before* is the *Bend* object before and *after* is the *Bend*
	    object after *self* on the contour.  *offset* is the offset
	    from the main contour to mill out.  *offset* is positive to expand
	    the contour and negative to shrink the contour.  *end_mill_radius*
	    is the end mill radius and it can be set to zero deterimine the
	    actual contour locations.
	"""
    
	# Verify argument types:
	assert isinstance(before, Bend)
	assert isinstance(after, Bend)
	assert isinstance(offset, L)
	assert isinstance(end_mill_radius, L)
	assert isinstance(tracing, int)

	# Use *bend* instead of *self*:
	bend = self
	bend_name = bend._name

	# Perform any requested *tracing*:
	trace_detail = 0
	if tracing >= 0:
	    print("{0}=>Bend._offset_compute('{1}', '{2}', '{3}', {4:i})".
	      format(' ' * tracing, bend_name, before._name, after._name, end_mill_radius))
	zero = L()
	smallest_inner_radius = L(mm=-1.0)
    
	# Extract the X/Y coordinates of the {before}, {bend}, and {after}:
	after_projected_point = after._projected_point
	after_x = after_projected_point.x
	after_y = after_projected_point.y
	bend_projected_point = bend._projected_point
	bend_x = bend_projected_point.x
	bend_y = bend_projected_point.y
	before_projected_point = before._projected_point
	before_x = before_projected_point.x
	before_y = before_projected_point.y
	if tracing >= 0 and trace_detail >= 1:
	    print(
	      "{0}c:({1}:{2}), b:({3}:{4}), a:({5}:{6}), o:{7}, emr:{8}".
	      format(' ' * tracing,
	      bend_x, bend_y, before_x, before_y, after_x, after_y, offset, end_mill_radius))
    
	# Compute the distances in X and Y between the bends:
	incoming_dx = before_x - bend_x
	incoming_dy = before_y - bend_y
	outgoing_dx = after_x - bend_x
	outgoing_dy = after_y - bend_y
    
	# Now we compute the absolute angles of the line segments coming into
	# and leaving the bend:
	assert incoming_dx != zero or incoming_dy != zero, \
	  "Bend '{0}' is at the same location as bend '{1}'".format(before._name, bend._name)
	assert outgoing_dx != zero or outgoing_dy != zero, \
	  "Bend '{0}' is at the same location as bend '{1}'".format(bend._name, after._name)
	if True:
	    # Here are some angle constants:
	    degrees0 = Angle(deg=0.0)
	    degrees90 = Angle(deg=90.0)
	    degrees180 = Angle(deg=180.0)

	    # Compute the incoming and outgoing angles as absolute angles:
	    incoming_angle = incoming_dy.arc_tangent2(incoming_dx)
	    incoming_bearing = (incoming_angle + degrees180).normalize()
	    outgoing_bearing = outgoing_dy.arc_tangent2(outgoing_dx)
	    bend._incoming_angle = incoming_angle
	    bend._incoming_bearing = incoming_bearing
	    bend._outgoing_bearing = outgoing_bearing
    
	    # Compute the amount that *incoming_angle* changed to get to
	    # *outgoing_angle*.  If the result is positive, we turned to the
	    # left and we have an outside bend, and if the result is negative,
	    # we turned to the right and we have an inside bend:
	    bearing_change = (outgoing_bearing - incoming_bearing).normalize()
	    bend._bearing_change = bearing_change
    
	    if tracing >= 0 and trace_detail >= 1:
		print("{0}ib=at2(y={1:i}, x={2:i})={3:d}".format(
		  ' ' * tracing, incoming_dy, incoming_dx, incoming_bearing))
		print("{0}ob=at2(y={1:i}, x={2:i})={3:d}".format(
		  ' ' * tracing, outgoing_dy, outgoing_dx, outgoing_bearing))
		print("{0}ib={1:d} + bc({2:d}) = ob({3:d})".format(
		  ' ' * tracing, incoming_bearing, bearing_change, outgoing_bearing))
    
	    half_angle = (outgoing_bearing - incoming_angle).normalize()/2
	    bend._half_angle = half_angle
	    center_angle = (incoming_angle + half_angle).normalize()
	    opposite_center_angle = (center_angle + degrees180).normalize()
	    bend._center_angle = center_angle
	    bend._opposite_center_angle = opposite_center_angle
    
	    # Compute the distance from the bend to the arc center:
	    radius = bend._radius
	    arc_center_offset = (radius / half_angle.sine()).absolute()
	    bend._arc_center_offset = arc_center_offset
    
	    # Compute the location of the arc center:
	    arc_center_x = bend_x + arc_center_offset * center_angle.cosine()
	    arc_center_y = bend_y + arc_center_offset * center_angle.sine()
	    bend._arc_center_x = arc_center_x
	    bend._arc_center_y = arc_center_y
    
	    xxx = degrees90 - half_angle.absolute()
	    arc_after_angle = degrees0
	    arc_before_angle = degrees0
	    inside_bend = bearing_change >= degrees0
	    if inside_bend:
		# Inside bend:
		if radius > zero:
		    smallest_inner_radius = radius
		arc_after_angle = (opposite_center_angle + xxx).normalize()
		arc_before_angle = (opposite_center_angle - xxx).normalize()
		if tracing >= 0 and trace_detail >= 1:
		    print("{0}Inside bend: '{1}'".format(' ' * tracing, bend._name))
	    else:
		# Outside bend:
		arc_after_angle = (opposite_center_angle - xxx).normalize()
		arc_before_angle = (opposite_center_angle + xxx).nomralize()
		if tracing >= 0 and trace_detail >= 1:
		    print("{0}Outside bend: '{1}'", format(' ' * tracing, bend._name))
	    bend._arc_after_angle = arc_after_angle
	    bend._arc_before_angle = arc_before_angle
    
	    if tracing >= 0 and trace_detail >= 1:
		print("{0}ha={1:d} ca={2:d} aaa={3:d} aba={4:d}".format(' ' * tracing,
		  half_angle, center_angle, arc_after_angle, arc_before_angle))
	
	    # Step 1: Figure everything out without *offset* and *end_mill_radius*:
	
	    if tracing and trace_detail >= 1:
		print("{0}Step 1a: r={1:i} o={2:i}".format(' ' * tracing, radius, offset))
		print("{0}Step 1b: aco={1:i} acx={2:i} acy={3:i}".format(
		  ' ' * tracing, arc_center_offset, arc_center_x, arc_center_y))
	
	    # Get these 4 variables defined; we fill them in later:
	    arc_after_x = zero
	    arc_after_y = zero
	    arc_before_x = zero
	    arc_before_y = zero
	
	    if tracing >= 0 and trace_detail >= 1:
		# These would be the values returned if exclude both *offset* and
		# *end_mill_radius* are ignored:
		arc_after_x =  arc_center_x + radius * arc_after_angle.cosine()
		arc_after_y =  arc_center_y + radius * arc_after_angle.sine()
		arc_before_x = arc_center_x + radius * arc_before_angle.cosine()
		arc_before_y = arc_center_y + radius * arc_before_angle.sine()
	
		print("{0}Step 1c: aax={1:i} aay={2:i} abx={3:i} aby={4:i}".format(
		  ' ' * tracing, arc_after_x, arc_after_y, arc_before_x, arc_before_y))
	
	    # Step 2: Now figure everything out with *offset*:
	
	    # Bad things happen when *offset* moves in "behind" the *radius*:
	    adjusted_radius = zero
	    if inside_bend:
		# Inside bend:
		adjusted_radius = radius - offset
	    else:
		# Outside bend:
		adjusted_radius = radius + offset
	    if adjusted_radius < zero:
		# *offset* is too far in and exceeds *radius*.  We need to compute
		# a new center based exclusively on -*offset* and set
		# *adjusted_radius* to zero.  Thus, if the end was rounded
		# (i.e. *radius* > 0.0), we have converted to a sharp bend
		# instead:
		arc_center_offset = -offset / half_angle.sine()
		arc_center_x = bend_x + arc_center_offset * center_angle.cosine()
		arc_center_y = bend_y + arc_center_offset * center_angle.sine()
		adjusted_radius = zero
	
		if tracing >= 0 and trace_detail >= 1:
		    print("{0}Step 2a: ar={1:i} aco={2:i} acx={3:i} acy={4:i}".format(' ' * tracing,
		      adjusted_radius, arc_center_offset, arc_center_x, arc_center_y))
	
	    if tracing >= 0 and trace_detail >= 1:
		# These would be the values return if we include {offset} but
		# exclude *end_mill_radius*:
		arc_after_x =  arc_center_x + adjusted_radius * arc_after_angle.cosine()
		arc_after_y =  arc_center_y + adjusted_radius * arc_after_angle.sine()
		arc_before_x = arc_center_x + adjusted_radius * arc_before_angle.cosine()
		arc_before_y = arc_center_y + adjusted_radius * arc_before_angle.sine()
	
		print("{0}Step 2b: ar={1:i} aax={2:i} aay={3:i} abx={4:i} aby={5:i}".format(
		  ' ' * tracing,
		  adjusted_radius, arc_after_x, arc_after_y, arc_before_x, arc_before_y))
	
	    # Step 3: Figure out everything with both *offset* and *end_mill_radius*:
	    arc_radius = zero
	    if inside_bend:
		# Inside bend:
		arc_radius = adjusted_radius - end_mill_radius
	    else:
		# Outside bend:
		arc_radius = adjusted_radius + end_mill_radius
	    arc_after_x =  arc_center_x + arc_radius * arc_after_angle.cosine()
	    arc_after_y =  arc_center_y + arc_radius * arc_after_angle.sine()
	    arc_before_x = arc_center_x + arc_radius * arc_before_angle.cosine()
	    arc_before_y = arc_center_y + arc_radius * arc_before_angle.sine()
	    
	    if tracing >= 0 and trace_detail >= 1:
		print("{0}Step 3: ar={1:i} aax={2:i} aay={3:i} abx={4:i} aby={5:i}".format(
		  ' ' * tracing,
		  arc_radius, arc_after_x, arc_after_y, arc_before_x, arc_before_y))
	
	    # Stuff the values into the {bend}:
	    bend._arc_after_x = arc_after_x
	    bend._arc_after_y = arc_after_y
	    bend._arc_after_angle = arc_after_angle
	    bend._arc_before_x = arc_before_x
	    bend._arc_before_y = arc_before_y
	    bend._arc_before_angle = arc_before_angle
	    bend._arc_center_x = arc_center_x
	    bend._arc_center_y = arc_center_y
	    bend._arc_radius = arc_radius
    
	if tracing >= 0:
	    print("{0}<=Bend._offset_compute('{1}', '{2}', '{3}', {4:i})=>{5:i}".format(
	       ' ' * tracing, self._name, before._name,
	       after._name, end_mill_radius, smallest_inner_radius))

	return smallest_inner_radius


    def _projected_point_get(self):
	""" *Bend*: Return the projected point of the *Bend* object (i.e. *self*). """

	return self._projected_point

    def _projected_point_set(self, projected_point):
	""" *Bend*: Return the projected point of the *Bend* object (i.e. *self*). """

	# Verify argument types:
	assert isinstance(projected_point, P)

	# Store *projected_point* into *self*:
	self._projected_point = projected_point

    def _point_get(self):
	""" *Bend* Return the corner point associated with the *Bend* object (i.e. *self*). """

	return self._point

    def _radius_get(self):
	""" *Bend* Return the corner radius associated with the *Bend* object (i.e. *self*). """

	return self._radius


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
	""" *Hertz*: Return the frequency of the *Hertz* object (i.e. *self*). """
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

    def _is_steel(self):
        """ *Material*: Return *True* if the *Material* object (i.e. *self*) is steel
	    (or iron) and *False* otherwise.
	"""

	return self._generic.lower() == "steel"

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

	formatting_text = "[bsw={0} tne={1}]"
	if format_text != "":
	    formatting_text = "[bsw={0:" + format_text + "} tne={1:" + format_text + "}]"

	if self._is_empty:
	    result = "[None]"
	else:
	    result = formatting_text.format(self.bsw_get(), self.tne_get())
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
	if not bounding_box._is_empty:
	    self.point_expand(bounding_box.tne_get())
	    self.point_expand(bounding_box.bsw_get())

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


    def volume_get(self):
	""" *Bounding_Box*: Return the volume of the *Bounding_Box* object (i.e. *self*.) """

	return P(self._east - self._west, self._north - self._south, self._top - self._bottom)

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

	return speed1._mm_per_sec == speed2._mm_per_sec

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
	else:
            assert False, "Unrecognized Speed format '{0}'".format(format_text)
	
	# Format the value as a string and return it:
	value = self._mm_per_sec / scale
	speed_format_text = "{0:" + format_text + "f}"
	#print("base={0} format_text='{1}' scaled={2} speed_format_text='{3}'".
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

	if format == "s":
	    result = "[{0:.2f}, {1:.2f}, {2:.2f}, {3:.2f}]".format(
	      self.red, self.green, self.blue, self.alpha)
	else:
	    result = "[r={0:.2f}, g={1:.2f}, b={2:.2f}, a={3:.2f}]".format(
	      self.red, self.green, self.blue, self.alpha)
	return result

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

    def _lookup(self, x, y, label):
	# Verify argument types:
	assert isinstance(x, float)
	assert isinstance(y, float)
	assert isinstance(label, str)

	# Is *point* already in *table*?:
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

    def __init__(self, name):
	""" *Contour*: Initialize the *Contour* object (i.e. *self*) with *name*"""

	# Verify argument types:
	assert isinstance(name, str)

	# Load up the *Contour* object (i.e. *self*):
	self._bends = []		# List of *Bend* objects used by by *Contour* object
	self._inside_bends_identified = False	# *True* if all inside bends have been identified
	self._is_clockwise = False	# *True* if the contour is clockwise; and *False* otherwise
	self._is_projected = False	# *True* if all of the *Bend* objects have been projected
	self._name = name		# The name of the contour
	self._project_transform = Transform()	# The projection *Transform*

    def __format__(self, format):
	""" *Contour*: Format the *Contour* object (i.e. *self*) as a string and return it. """

	# Verify argument types:
	assert isinstance(format, str)

	# Create a list of *chunks* to contentate together into a the final string.  Initilize
	# it with the opening curly brace:
	chunks = ["{"]

	# Visit each *bend* in *bends* to add a new *chunk* to *chunks*:
	bends = self._bends
	for index, bend in enumerate(bends):
	    chunk = "[{0}]: {1}".format(index, bend)
	    chunks.append(chunk)

	# Tack the closing curly brace onto *chunks*:
	chunks.append("}")

	# Convert *chunks* into a space separated *result* string:
	result = '\n'.join(chunks)
	return result

    def _bends_get(self):
	""" *Contour*: Return the list of *Bend* objects associated with the *Contour* object
	    (i.e. *self).
	"""

	return self._bends

    def _bounding_box_expand(self, bounding_box):
	""" *Contour*: Expand *bounding_box* by each point in the *Contour* object (i.e. *self). """

	# Verify argument types:
	assert isinstance(bounding_box, Bounding_Box)

	# Visit each *bend* in *bends* and exp
	bends = self._bends
	for bend in bends:
	    point = bend._point_get()
	    bounding_box.point_expand(point)

    def _copy(self, name):
        """ *Contour*: Return a copy of the *Contour* object (i.e. *self*)
	    with a new name of *name*.
	"""

	# Make a copy of the bends:
	new_bends = []
	for bend in self._bends:
	    new_bends.append(bend._copy())

	# Make a new *contour* object and fill it in:
	contour = Contour(name)
	contour._bends = new_bends
        contour._inside_bends_identify = self._inside_bends_identify
        contour._is_clockwise = self._is_clockwise
        contour._is_isprojected = self._is_projected
        contour._name = self._name
	contour._project_transform = self._project_transform
	
	return contour

    def _inside_bends_identify(self, tracing):
	""" *Contour*: Sweep through each *Bend* object in the *Contour* object (i.e. *self*)
	    and determine if it is an inside bend or an outside bend.  This routine will
	    also set the *is_clockwise* field of the *Contour* object to *True* if the
	    contour is clockwise and *False* otherwise.
	"""

	# How do you figure out whether a *Contour* clockwise vs. counter-clockwise?
	# Basically, the trick is to traverse the contour in a forward direction.
	# Each bend will change the bearing some value between -180 degrees and
	# +180 degrees.  If you compute the sum of all these bearing changes, you
	# will wind up with -360 degrees for a clockwise contour and +360 degrees
	# for a counter-clockwise contour.
	#
	# Once we know the contour direction, we can sweep through all of the bends
	# and determine which bends are inside or outside by inspecting the bearing
	# change.  For clockwise contour, outside bends have a negative bearing change and
	# inside bends have a positive bearing change.  For a counter-clockwise contour,
	# outside bends negative and inside bends are positive.

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform an requested tracing:
	#tracing = 0
	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 3
	    indent = ' ' * tracing
	    print("{0}=>Contour._inside_bends_identify('{1}')".format(indent, self._name))

	# Make srue that we have called *Contour._project*() on the *Contour* object (i.e. *self*):
	assert self._is_projected

	# Some constants:
	degrees0 = Angle()
	degrees180 = Angle(deg=180.0)
	minus_degrees180 = -degrees180
	degrees360 = Angle(deg=360.0)
	large = L(mm=1.e10)
	zero = L()

	# Grab *bends* from the *Contour* object (i.e. *self*) and compute it size:
	bends = self._bends
	bends_size = len(bends)

	# Visit each *bend* in *bends*:
	bearing_changes_sum = degrees0
	bearing_change_angles = []
	for index, at_bend in enumerate(bends):
	    # Get *before_bend* and *after_bend*:
	    before_bend = bends[(index - 1) % bends_size]
	    after_bend = bends[(index + 1) % bends_size]

	    # Grab the X/Y coordinates for *before_bend*, *at_bend*, and *after_bend*:
	    before_projected_point = before_bend._projected_point_get()
            before_x = before_projected_point.x
            before_y = before_projected_point.y
	    at_projected_point = at_bend._projected_point_get()
	    at_x = at_projected_point.x
	    at_y = at_projected_point.y
	    after_projected_point = after_bend._projected_point_get()
	    after_x = after_projected_point.x
	    after_y = after_projected_point.y
	    if tracing_detail >= 3:
		print("{0}bx={1:i} by={2:i} @x={3:i} @y={4:i} ax={5:i} ay={6:i}".format(
		  indent, before_x, before_y, at_x, at_y, after_x, after_y))

	    # Compute the change in dx and dy for entry to and exit from *at_bend*:
	    before_dx = at_x - before_x
	    before_dy = at_y - before_y
	    at_dx = after_x - at_x
	    at_dy = after_y - at_y

	    # Now compute the *bearing_change_angle*:
	    before_to_at_angle = before_dy.arc_tangent2(before_dx)
	    at_to_after_angle = at_dy.arc_tangent2(at_dx)
	    bearing_change_angle = at_to_after_angle - before_to_at_angle

	    # Now adjust *bearing_change_angle* to be between -180 degrees and +180 degrees:
	    if bearing_change_angle < minus_degrees180:
		bearing_change_angle += degrees360
	    elif bearing_change_angle > degrees180:
		bearing_change_angle -= degrees360

	    # Hang onto the *bearing_change_angle* for the inside/outside bend step below:
	    bearing_change_angles.append(bearing_change_angle)	

	    # Update *bearing_changes_sum*:
	    bearing_changes_sum += bearing_change_angle

	# Now we figure out if the contour *is_clockwise* or not:
	is_clockwise = bearing_changes_sum <= degrees0
	self._is_clockwise = is_clockwise
	if tracing >= 0:
	    print("{0}is_clockwise={1}".format(indent, is_clockwise))

	# Now sweep through each of the *bends* and set the determine if it *is_inside* or not:
        outside_bends = 0
	inside_bends = 0
	for index, at_bend in enumerate(bends):
	    bearing_change_angle = bearing_change_angles[index]
	    if is_clockwise:
		is_inside = bearing_change_angle >= degrees0
	    else:
		is_inside = bearing_change_angle <= degrees0
	    at_bend._is_inside_set(is_inside)

	    if is_inside:
		inside_bends += 1
	    else:
		outside_bends += 1

	    if tracing_detail >= 0:
		print("{0}bend['{1}']: is_inside={2} angle_change={3:d}".
		  format(indent, at_bend._name_get(), is_inside, bearing_change_angle))

	# Remember that we have the *inside_bends_identified*:
	self._inside_bends_identified = True

	if tracing >= 0:
	    print(
	      "{0}<=Contour._inside_bends_identify('{1}') clockwise={2} insides={3} outsides={4}".
	      format(indent, self._name, is_clockwise, inside_bends, outside_bends))


    def _is_clockwise_get(self):
	""" *Contour*: Return *True* if the *Contour* object (i.e. *self*) is a clockwise
	    contour or *False* for a counter-clockwise contour.
	"""

	return self._is_clockwise

    def _name_get(self):
	""" *Contour*: Return the name of the *Contour* object (i.e. *self*). """

	return self._name

    def _project(self, transform, tracing = -1000000):
	""" *Contour*: Sweep through the *Contour* object (i.e. *self*) and compute the
	    X/Y coordinates of *Bend* object on an X/Y plane using *transform* each point
	    in the contour prior prior to being projected down to the X/Y plane.
	"""

	# Verify argument types:
	assert isinstance(transform, Transform)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	tracing_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    tracing_detail = 3
	    print("{0}=>Contour._project('{1}', {2:s})".format(indent, self._name, transform))

	# For each *bend* in bends, perform the project from 3D down to 2D:
	zero = L()
	for bend in self._bends:
	    # Grab the X/Y/Z coordinates for *point*:
	    name = bend._name_get()
	    point = bend._point_get()
	    transformed_point = transform * point
	    projected_point = P(transformed_point.x, transformed_point.y, zero)
	    bend._projected_point_set(projected_point)
	    if tracing_detail >= 2:
		print("{0}name='{1}' point={2:m} projected_point={3:m} transform={4:s}".
		  format(indent, name, point, projected_point, transform))

	# Remember that we have performed the projection:
	self._is_projected = True

	# Hang onto the projection *transform* so we can eventually unproject:
	self._project_transform = transform        

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Contour._project('{1}', {2:s})".format(indent, self._name, transform))

    def _radius_center_and_tangents_compute(self, tracing):
	""" *Contour*: """

	# Verify argument types:
	assert isinstance(tracing, int)

	assert self._is_projected
	assert self._inside_bends_identified

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Contour._radius_center_and_tangents_compute('{1}')".
	      format(indent, self._name))

	bends = self._bends
	bends_size = len(bends)
	for index, at_bend in enumerate(self._bends):
	    incoming_bend = bends[(index - 1) % bends_size]
	    outgoing_bend = bends[(index + 1) % bends_size]
	    at_bend._radius_center_and_tangents_compute(incoming_bend, outgoing_bend, tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Contour._radius_center_and_tangents_compute('{1}')".
	      format(indent, self._name))

    def _scad_lines_polygon_append(self, scad_lines, pad, box_enclose, tracing = -1000000):
	""" *Contour*: """

	# Verify argument types:
	assert isinstance(scad_lines, list)
	assert isinstance(pad, str)
	assert isinstance(box_enclose, bool)
	assert isinstance(tracing, int)

	# Use *contour* instead of *self*:
	contour = self

	# Perform any requested *tracing*:
	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 3
	    indent = ' ' * tracing
	    print("{0}=>Contour._scad_lines_polygon_append('{1}', *, '{2}, {3})".
	      format(indent, contour._name, pad, box_enclose))

	# Grab *bends* from *contour*:
	bends = contour._bends_get()
	assert len(bends) >= 1

	# Compute the bounds of a box that will enclose the contour:
	if box_enclose:
	    # Visit each *bend* in *bends8:
	    for index, bend in enumerate(bends):
		# Extract the *x*/*y*/*z* from *bend*:
		point = bend._projected_point_get()
		x = point.x
		y = point.y
		z = point.z

		# The first time through, we just set the miminims/maximums from x/y/z:
		if index == 0:
		    x_maximum = x_minimum = x
		    y_maximum = y_minimum = y
		    z_maximum = z_minimum = z
		else:
		    # All other times, we compute the maximum and minimum:
		    x_minimum = min(x_minimum, x)
		    y_minimum = min(y_minimum, y)
		    z_minimum = min(z_minimum, z)
		    x_maximum = max(x_maximum, x)
		    y_maximum = max(y_maximum, y)
		    z_maximum = max(z_maximum, z)

	    # FIXME: *trim* extra should be computed:
	    # Now tack on extra box area:
	    trim_extra = L(inch=36.0)
	    x1 = x_minimum - trim_extra
	    x2 = x_maximum + trim_extra
	    y1 = y_minimum - trim_extra
	    y2 = y_maximum + trim_extra
	    z1 = z_minimum - trim_extra
	    z2 = z_maximum + trim_extra

	# Start by outputing the polygon directive:
	scad_lines.append("{0}polygon(".format(pad))
            
	# First output each *bend_point* in *bends*::
        contour_pairs = []
	contour_is_clockwise = contour._is_clockwise_get()
	scad_lines.append("{0}  points = [".format(pad))
	bends_size = len(bends)
	for bend_index, bend in enumerate(bends):
	    # Grab some values from *bend*:
	    bend_point = bend._projected_point_get()	# Should this be _projected_bend_point?
	    bend_name = bend._name_get()
	    bend_radius = bend._radius_get()
	    bend_is_inside = bend._is_inside_get()
		
	    # The *arc_radius* is positive or negative depending upon whether
	    # *contour_is_clocwise* or the *bend_is_inside*:
	    if contour_is_clockwise:
		if bend_is_inside:
		    arc_radius = -bend_radius
		else:
		    arc_radius = bend_radius
	    else:
		if bend_is_inside:
		    arc_radius = bend_radius
		else:
		    arc_radius = -bend_radius

	    # We are going to sweep out some line segments from *start_angle* to *end_angle*
	    # (which are angle relative to *center*):
	    center = bend._center_get()
	    start_angle = bend._incoming_tangent_angle_get()
	    start_direction = start_angle.xy_direction()
	    end_angle = bend._outgoing_tangent_angle_get()
	    change_angle = end_angle - start_angle

	    # Make sure that the *bearing_angle_change* is positive for clockwise sweeps
	    # and negative for counter clockwise sweeps:
	    degrees0 = Angle()
	    degrees360 = Angle(deg=360.0)
	    if change_angle >= degrees0:
		if contour_is_clockwise and not bend_is_inside or \
		 not contour_is_clockwise and bend_is_inside:
		    change_angle -= degrees360
		    if tracing_detail >= 1:
			print("{0}decrement".format(indent))
	    elif change_angle < degrees0:
		if contour_is_clockwise and bend_is_inside or \
		  not contour_is_clockwise and not bend_is_inside:
		    if tracing >= 1:
			print("{0}increment".format(indent))
		    change_angle += degrees360
	    if tracing_detail >= 0:
		print("{0}name='{1}' start_angle={2:d} end_angle={3:d} change_angle={4:d}".
		  format(indent, bend_name, start_angle, end_angle, change_angle))

	    # Now we compute *deta_angle*:
	    arc_radius_mm = arc_radius.millimeters()
	    count = max(4, int(arc_radius_mm))
	    delta_angle = change_angle / float(count - 1)
	    if tracing_detail >= 0:
		print("{0}count={1} delta_angle={2:d}".format(indent, count, delta_angle))
	    comma = ","
	    for index in range(count):
		angle = (start_angle + delta_angle * float(index)).normalize()
		point = center + angle.xy_direction() * arc_radius_mm
		if tracing_detail >= 0:
		    print("{0}[{1}/{2}]:angle={3:d} point={4:i}".
		      format(indent, index + 1, count, angle, point))

		if bend_index + 1 == bends_size and index + 1 == count and not box_enclose:
		    comma = ""
		scad_lines.append("{0}   [{1:m}, {2:m}]{3} // {4} {5} of {6}".
		  format(pad, point.x, point.y, comma, bend_name, index + 1, count))
		contour_pairs.append( (point, "{0} {1} of {2}".
		  format(bend_name, index + 1, count) ))

	# Now tack on the 4 points that enclose the contour with a box:
	if box_enclose:
	    scad_lines.append("{0}   [{1:m}, {2:m}], // box SW".format(pad, x1, y1))
	    scad_lines.append("{0}   [{1:m}, {2:m}], // box NW".format(pad, x1, y2))
	    scad_lines.append("{0}   [{1:m}, {2:m}], // box NE".format(pad, x2, y2))
	    scad_lines.append("{0}   [{1:m}, {2:m}]  // box SE".format(pad, x2, y1))

	# Close off all the points:
	scad_lines.append("{0}  ],".format(pad))

	# Now we output the outermost path the encloses the contour:
	contour_pairs_size = len(contour_pairs)
	scad_lines.append("{0}  paths = [".format(pad))
	if box_enclose:
	    scad_lines.append("{0}   [ // box path".format(pad))
	    scad_lines.append("{0}    {1}, // SW box".format(pad, contour_pairs_size))
	    scad_lines.append("{0}    {1}, // NW box".format(pad, contour_pairs_size + 1))
	    scad_lines.append("{0}    {1}, // NE box".format(pad, contour_pairs_size + 2))
	    scad_lines.append("{0}    {1}  // SE box".format(pad, contour_pairs_size + 3))
	    scad_lines.append("{0}   ], // box path".format(pad))

	# Now we output the contour path:
	scad_lines.append("{0}   [ //contour_path".format(pad))
	comma = ","
        for index, contour_pair in enumerate(contour_pairs):
	    contour_point = contour_pair[0]
	    contour_text = contour_pair[1]
	    if index >= contour_pairs_size - 1:
		comma = ""
	    scad_lines.append("{0}    {1}{2} // {3}".format(pad, index, comma, contour_text))
	scad_lines.append("{0}   ] // contour path".format(pad))
	scad_lines.append("{0}  ] // paths".format(pad))
	scad_lines.append("{0} ); // polygon".format(pad))

	# Wrap up any *tracing*:
	if tracing >= 0:
	    print("{0}=>Contour._scad_lines_polygon_append('{1}', *, '{2}, {3})".
	      format(indent, contour._name, pad, box_enclose))

    def _smallest_inner_radius_compute(self, tracing):
	""" *Contour*: Return the smallest bend radius of any the inside bend for the
	    *Contour* object (i.e. *self).  If there are no inside radius, a negative
	    radius is returned.
	"""

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	#tracing = 0
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Contour._smallest_inner_radius_compute('{1}')".format(indent, self._name))

	# Make srue that we have called *Contour._project*() on the *Contour* object (i.e. *self*):
	assert self._is_projected
	assert self._inside_bends_identified

	# Compute the *smallest_inner_radius*:
	zero = L()
	smallest_inner_radius = -L(mm=1.0)	# Start with a big radius:
	for bend in self._bends:
	    # Grab some values from *bend*:
	    point = bend._point_get()
	    radius = bend._radius_get()
	    is_inside = bend._is_inside_get()
	    if tracing >= 0:
		print("{0}bend: point={1:i} radius={2:i} is_inside={3}".
		  format(indent, point, radius, is_inside))

	    if is_inside:
		if smallest_inner_radius < zero:
		    smallest_inner_radius = radius
		else:
		    smallest_inner_radius = smallest_inner_radius.minimum(radius)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Contour._smallest_inner_radius_compute('{1}') => {2}".
	      format(indent, self._name, smallest_inner_radius))

	return smallest_inner_radius

    def _unproject(self, tracing = -1000000):
	""" *Contour*: Unproject the adjusted bend points in the *Contour* to be back
	    in the along the original projection axis. """

	# Verify argument types:
        assert isinstance(tracing, int)

	# Make sure that the *Contour* object (i.e. *self*) was previously projected:
	assert self._is_projected

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = " " * tracing
	    print("{0}=>Contour._unproject()".format(indent))

	# For each *bend* in *bends*, apply the *unproject_transform*:
	unproject_transform = self._project_transform.reverse()
	for bend in self._bends:
	    projected_point = bend._projected_point
	    assert isinstance(projected_point, P)
	    bend._point	= unproject_transform * bend._projected_point

	# Wrap any tracing:
	if tracing >= 0:
	    print("{0}<=Contour._unproject()".format(indent))

    def adjust(self, name, delta, start, end, maximum_radius, minimum_radius, tracing = -1000000):
	""" *Contour*: """

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(delta, L)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(maximum_radius, L)
	assert isinstance(minimum_radius, L)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = " " * tracing
	    print("{0}=>Contour.adjust('{1}', {2}, {3}, {4}, {5}, {6})".
              format(indent, name, delta, start, end, maximum_radius, minimum_radius))

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

	adjusted_contour = self._copy(name)

	degrees0 = Angle(deg=0.0)
	top_surface_transform = Transform.top_surface("top", start, end, degrees0, tracing + 1);
	adjusted_contour._project(top_surface_transform, tracing + 1)

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
	    # (*ax3*, *ay3):

	    # *bend1* occurs *bend2* occurs before *bend3*:
	    bend1 = bends[(bends_index - 1) % bends_size]
	    bend2 = bends[bends_index]
	    bend3 = bends[(bends_index + 1) % bends_size]
	    if tracing >= 0:
		print("{0}bends_index={1}".format(indent, bends_index))
	    
	    # Get the [o]riginal X/Y values:
	    ox1 = bend1._point.x._mm
	    oy1 = bend1._point.y._mm
	    ox2 = bend2._point.x._mm
	    oy2 = bend2._point.y._mm
	    ox3 = bend2._point.x._mm
	    oy3 = bend2._point.y._mm
	    ox4 = bend3._point.x._mm
	    oy4 = bend3._point.y._mm
	    if tracing >= 0:
		print("{0}[{1}]:o1=({2},{3}) o2=({4},{5}) o3=({6},{7}) o4=({8},{9})".
		  format(indent, bends_index, ox1, oy1, ox2, oy2, ox3, oy3, ox4, oy4))

	    # Compute vector (*odx12*, *ody12*) from (*ox1*, *oy1*) to
	    # (*ox2*, *oy2*):
	    odx12 = ox2 - ox1
	    ody12 = oy2 - oy1

	    # Compute the vector (*odx34*, *ody34*) from (*ox3*, *oy3*) to
	    # (*ox4*, *oy4*):
	    odx34 = ox4 - ox3
	    ody34 = oy4 - oy3
	    if tracing >= 0:
		print("{0}[{1}]:od12=({2},{3}) od34=({4},{5})".
		  format(indent, bends_index, odx12, ody12, odx34, ody34))

	    # Compute the lengths of (*dx12*, *dy12*) and (*dx34*, *dy34*):
	    length12 = math.sqrt(odx12 * odx12 + ody12 * ody12)
	    length34 = math.sqrt(odx34 * odx34 + ody34 * ody34)
	    if tracing >= 0:
		print("{0}[{1}]:len12={2} len34={3}".
		  format(indent, bends_index, length12, length34))

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
	    if tracing >= 0:
		print("{0}[{1}]:n12=({2},{3}) n34=({4},{5})".
		  format(indent, bends_index, ndx12, ndy12, ndx34, ndy34))

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
	    if tracing >= 0:
		print("{0}[{1}]:a1=({2},{3}) a2=({4},{5}) a3=({6},{7}) a4=({8},{9})".
		  format(indent, bends_index, ax1, ay1, ax2, ay2, ax3, ay3, ax4, ay4))

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
	    if tracing >= 0:
		print("{0}[{1}]:det12={2} det34={3}".
		  format(indent, bends_index, det12, det34))
		print("{0}[{1}]:x_num={2} y_num={3} denom={4}".
		  format(indent, bends_index, x_numerator, y_numerator, denominator))
		print("{0}[{1}]:ix={2} iy={3}".format(indent, bends_index, ix, iy))

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

	    new_radius = bend2._radius._mm
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
	    adjusted_contour_bend2 = adjusted_contour._bends[bends_index]
	    adjusted_contour_bend2._radius = L(mm=new_radius)
	    adjusted_contour_bend2._projected_point = P(L(mm=ix), L(mm=iy), zero)

	# Now map the bend points back to where they should be:
	adjusted_contour._unproject(tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Contour.adjust('{1}', {2}, {3}, {4}, {5}, {6})".
              format(indent, name, delta, start, end, maximum_radius, minimum_radius))
	return adjusted_contour

    def bend_append(self, name, point, radius):
	""" *Contour*: Create a *Bend* object containing *name*, *point*, and *radius* and
	    append the resulting *Bend* object to the *Contour* object (i.e. *self*).
	    *name* is used for debugging and shows up in generated G-code.  *point*
	    is the location at which a 0 degree radius bend would occur.  *radius* is the
	    desired bend radius.  This routine returns the created *Bend* object.
	    
	"""

	# Check argument types:
	assert isinstance(point, P)
	assert isinstance(radius, L)
	assert isinstance(name, str)

	# Create and append the bend:
	bend = Bend(name, point, radius)
	self._bends.append(bend)
	return bend

    # Used for extrusion:
    def _path_append(self, indexed_points, maximum_angle, tracing = -1000000):
	""" *Contour*: Append a path to *indexed_points*. """

	# Verify argument types:
	assert isinstance(indexed_points, Indexed_Points)
	assert isinstance(maximum_angle, Angle)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    print("{0}=>Contour.path_append(extrude_axis={1}, self={2})".
	          format(trace * ' ', extrude_axis, self))

	bends = self._bends
	for bend in bends:
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
	self._bases_table = {}

	self._bounding_box_dispatch = {
	  "b":      Bounding_Box.b_get,
	  "be":     Bounding_Box.be_get,
	  "bn":     Bounding_Box.bn_get,
	  "bne":    Bounding_Box.bne_get,
	  "bnw":    Bounding_Box.bnw_get,
	  "bs":     Bounding_Box.bs_get,
	  "bse":    Bounding_Box.bse_get,
	  "bsw":    Bounding_Box.bsw_get,
	  "bw":     Bounding_Box.bw_get,
	  "c":      Bounding_Box.c_get,
	  "e":      Bounding_Box.e_get,
	  "n":      Bounding_Box.n_get,
	  "ne":     Bounding_Box.ne_get,
	  "nw":     Bounding_Box.nw_get,
	  "s":      Bounding_Box.s_get,
	  "se":     Bounding_Box.se_get,
	  "sw":     Bounding_Box.sw_get,
	  "t":      Bounding_Box.t_get,
	  "te":     Bounding_Box.te_get,
	  "tn":     Bounding_Box.tn_get,
	  "tne":    Bounding_Box.tne_get,
	  "tnw":    Bounding_Box.tnw_get,
	  "ts":     Bounding_Box.ts_get,
	  "tse":    Bounding_Box.tse_get,
	  "tsw":    Bounding_Box.tsw_get,
	  "tw":     Bounding_Box.tw_get,
	  "volume": Bounding_Box.volume_get,
	  "w":      Bounding_Box.w_get,
	}

	# Make sure various output directories exist:
	self._directory      = directory
	self._dxf_directory  = self._directory_create("dxf", True)
	self._ngc_directory  = self._directory_create("ngc", True)
	self._scad_directory = self._directory_create("scad", False)
	self._stl_directory  = self._directory_create("stl", False)
	self._wrl_directory  = self._directory_create("wrl", True)

	EZCAD3.ezcad = self

    def _directory_create(self, sub_directory, clear):
        """ *EZCAD3*: Force directory named *directory*/*sub_directory* into existence, where
	    *directory* comes for the *EZCAD3* object (i.e. *self.*)  If *clear* is *True*,
	    the directory is cleared.  The resulting directory path is returned. """

	# Create the *full_directory_path*:
	ezcad3 = self
	directory = ezcad3._directory
	full_directory_path = os.path.join(directory, sub_directory)

        # See whether or not the directory already exists:
	if os.path.exists(full_directory_path):
	    # Directory exists:
	    if clear:
		# Clear out  each individual file in the directory using glob style
   		# file matching expression (e.g. `.../wrl/*`):
		for file_name in glob.glob(os.path.join(full_directory_path, "*")):
		    os.remove(file_name)
	else:
	    # Directory does not exist, so create it:
            os.makedirs(full_directory_path)

	return full_directory_path

    def _directory_get(self):
	""" *EZCAD3*: Return the directory to read/write files from/into from the *EZCAD3* object
	    (i.e. *self*).
        """

	return self._directory

    def _dxf_directory_get(self):
	""" *EZCAD3*: Return the directory to read/write DXF files from/into from the *EZCAD3*
	    object (i.e. *self*).
        """

	return self._dxf_directory

    def _ngc_directory_get(self):
	""" *EZCAD3*: Return the directory to read/write NGC files from/into from the *EZCAD3*
	    object (i.e. *self*).
        """

	return self._ngc_directory

    def _scad_directory_get(self):
	""" *EZCAD3*: Return the directory to read/write SCAD files from/into from the *EZCAD3*
	    object (i.e. *self*).
        """

	return self._scad_directory

    def _stl_directory_get(self):
	""" *EZCAD3*: Return the directory to read/write STL files from/into from the *EZCAD3*
	    object (i.e. *self*).
        """

	return self._stl_directory

    def _wrl_directory_get(self):
	""" *EZCAD3*: Return the directory to read/write WRL files from/into from the *EZCAD3*
	    object (i.e. *self*).
        """

	return self._wrl_directory

    def _mode_get(self):
	""" *EZCAD3*: Return the mode of the *EZCAD3* object (i.e. *self*). """

	return self._mode


    def process(self, part):
	""" *EZCAD3*: Perform all of the processing starting at *part*.
	"""

	assert isinstance(part, Part)
	part.process(self)

    @staticmethod
    def update_count_get():
	update_count = 123456789
	try:
	    ezcad = EZCAD3.ezcad
	    if ezcad._mode != EZCAD3.DIMENSIONS_MODE:
		update_count = ezcad._update_count
        except:
            assert False
	return update_count

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
    ORDER_DRILL =                       14
    ORDER_VERTICAL_LATHE =		15
    ORDER_LAST =			16

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

    def __init__(self, name, kind,
      part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed):
	""" *Operation*: Initialize an *Operation* object to contain
	    *name*, *kind*, *part*, *comment*, *sub_priority*, *tool*,
	    *order*, *follows*, *feed_speed*, and *spindle_speed*.
	"""

	# Use *operation* instead of *self*:
	operation = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(kind, int)
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(order, int) and Operation.ORDER_NONE < order < Operation.ORDER_LAST
	assert isinstance(follows, Operation) or follows is None
	assert isinstance(feed_speed, Speed)
	assert isinstance(spindle_speed, Hertz)

	# Load up *operation*:
	operation._name = name
	operation._kind = kind
	operation._part = part
	operation._comment = comment
	operation._follows = follows
	operation._index = -1
	operation._position = 0 # part.position_count
	operation._order = order
	operation._priority = part._priority_get()
	operation._sub_priority = sub_priority
	operation._tool = tool
	operation._feed_speed = feed_speed
	operation._spindle_speed = spindle_speed

    def _compare(self, operation2):
	""" *Operation*: Return -1, 0, or 1 depending upon whether *operation1*
	    (i.e. *self*) should sort before, at, or after *operation2*.
	"""

	# Use *operation1* instead of *self*:
	operation1 = self

	# Verify argument types:
	assert isinstance(operation2, Operation)

	# Sort by *priority* first:
	result = int_compare(operation1._priority, operation2._priority)
	if result != 0:
	    return result

	# Sort by *sub_priority* second:
	result = int_compare(operation1._sub_priority, operation2._sub_priority)
	if result != 0:
	    return result

	# Sort by *tool* third:
	result = operation1._tool._compare(operation2._tool)
	if result != 0:
	    return result

	# Sort by *order* fourth:
	result = int_compare(operation1._order, operation2._order)
	if result != 0:
	    return result

	# Sort by *kind* fifth:
	result = int_compare(operation1._kind, operation2._kind)
	if result != 0:
	    return result

	# Sort by *vice_x* sixth:
	#result = operation1._vice_x.compare(operation2._vice_x)
	#if result != 0:
	#    return result

	# Sort by *vice_y* seventh:
	#result = operation1._vice_y.compare(operation2._vice_y)
	#if result != 0:
	#    return result

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

#    def _vice_x_get(self):
#	""" *Operation*: Return vice Z coordinate for the *Operation* object (i.e. *self*). """
#
#	return self._vice_x
#
#    def _vice_y_get(self):
#	""" *Operation*: Return vice Z coordinate for the *Operation* object (i.e. *self*). """
#
#	return self._vice_y

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
    # both inward and outward.  There is a variable called *contour_offset*,
    # which we will call o, and another variable called *finish_offset*
    # which we will call f.  The contour offset can be positive or negative
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

    def __init__(self,
      part, comment, sub_priority, mill_tool, order, follows, feed_speed, spindle_speed,
      z_start, z_stop, contour, offset, effective_tool_radius, passes, tracing):
	""" *Operation_Contour*: Initialize the *Operation_Contour* object (i.e. *self*)
	    with *part*, *comment*, *sub_priorit*, *mill_tool*, *order*, *follows*,
	    *z_start*, *z_stop*, 
	"""

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(mill_tool, Tool)
	assert isinstance(order, int)
	assert isinstance(follows, Operation) or follows is None
	assert isinstance(feed_speed, Speed)
	assert isinstance(spindle_speed, Hertz)
	assert isinstance(z_start, L)
	assert isinstance(z_stop, L)
	assert isinstance(contour, Contour)
	assert isinstance(offset, L)
	assert isinstance(effective_tool_radius, L)
	assert isinstance(passes, int) and passes > 0
	assert isinstance(tracing, int)

	if tracing >= 1:
	    indent = ' ' * tracing
	    print("{0}=>Operaton_Contour.__init__(*, '{1}', '{2}', {3}, '{4}', {5}, {6},".
	      format(indent, part._name_get(), comment, sub_priority, mill_tool._name_get(),
	      order, isinstance(follows, Operation)))
	    print("{0} s={1:i} f={2:rpm} zs={3:i} ze={4:i} cntr='{5}' off={6:i} tr={7:i} p={8})".
	      format(indent, feed_speed, spindle_speed, z_start, z_stop,
	      contour._name_get(), offset, effective_tool_radius, passes))

	# Initialize super class:
	Operation.__init__(self, "Contour", Operation.KIND_CONTOUR,
	  part, comment, sub_priority, mill_tool, order, follows, feed_speed, spindle_speed)

	# Load up the rest of *self*:
	self._z_start = z_start
	self._z_stop = z_stop
	self._contour = contour
	self._offset = offset
	self._effective_tool_radius = effective_tool_radius
	self._passes = passes

	if tracing >= 1:
	    indent = ' ' * tracing
	    print("{0}<=Operaton_Contour.__init__(*, '{1}', '{2}', {3}, '{4}', {5}, {6},".
	      format(indent, part._name_get(), comment, sub_priority, mill_tool._name_get(),
	      order, isinstance(follows, Operation)))
	    print("{0} s={1:i} f={2:rpm} zs={3:i} ze={4:i} cntr='{5}' off={6:i} tr={7:i} p={8})".
	      format(indent, feed_speed, spindle_speed, z_start, z_stop,
	      contour._name_get(), offset, effective_tool_radius, passes))

    def _cnc_generate(self, tracing=-1000000):
	""" *Operation_Contour*: Generate the CNC code for the *Operation_Contour* object
	    (i.e. *self*).
	"""

	# Use *operaton* instead of *self*:
	operation = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Grab some values from *contour*:
	part = operation._part
	shop = part._shop_get()
	code = shop._code_get()
	tool = operation._tool
	comment = operation._comment

	# Perform an requested *tracing*:
	if tracing >= 0:
	    indent = " " * tracing
	    print("{0}=>Operation_Contour._cnc_generate".format(indent))

	# Record *comment* into *code*:
	code._line_comment(comment)

	# Mark whether *tool* is a laser or not:
	is_laser = tool._is_laser_get()
	code._is_laser_set(is_laser)

	# Grab some values from *tool*:
	s = operation._spindle_speed_get()
	f = operation._feed_speed_get()
	tool_diameter = tool._diameter_get()
	z_feed = f / 2

	# Grap some more values from *operation*:
	z_start = operation._z_start
	z_stop = operation._z_stop
	passes = operation._passes
	contour = operation._contour
	offset = operation._offset
	radius = operation._effective_tool_radius

	# Do the setup step for doing the contour:
	contour._radius_center_and_tangents_compute(tracing + 1)

	z_depth = z_start - z_stop
	depth_per_pass = z_depth / float(passes)

	# Now generate the G-Code.  We visit each corner once, and the
	# first corner twice.  Thus, we need to loop through *size* + 1 times:

	code._z_safe_assert("contour", comment)

	plunge_offset = tool_diameter
	#call d@("plunge_offset temporary set to zero\n")
	zero = L()

	if tracing >= 0:
            print("{0}passes={1}".format(indent, passes))
	for index in range(passes):
	    code._line_comment("Pass {0} of {1}".format(index + 1, passes))

	    # Get cutter down to the correct depth:
	    z = z_start - depth_per_pass * float(index + 1)

	    #call d@(form@("tool=%v% plunge_offset=%i%\n\") %
	    #  f@(tool.name) / f@(plunge_offset))

	    contour = operation._contour
	    code._contour(contour, plunge_offset, offset, radius, True, z, f, s, tracing + 1)

	code._z_safe_retract(z_feed, s)

	if tracing >= 0:
	    print("{0}<=Operation_Contour._cnc_generate".format(indent))

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
      diameter, dowel_point, plunge_point, top_surface_z, xy_rapid_safe_z, tracing=-1000000):
	""" *Operation_Dowel_Pin*: Initialize an *Operation_Dowel_Pin* object (i.e. *self*)
	    to contain *part*, *comment*, *sub_priority*, *tool*, *order*, *follows*,
	    *feed_speed*, *spindle_speed*, *diameter*, *dowel_point*, *plunge_point*,
	    *top_surface_z*, and *xy_rapid_safe_z*.
	"""

	# Use *operation_dowel_pin* instead of *self*:
        operation_dowel_pin = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(follows, Operation) or follows == None
	assert isinstance(feed_speed, Speed) and feed_speed != Speed()
	assert isinstance(spindle_speed, Hertz)
	assert isinstance(diameter, L)
	assert isinstance(dowel_point, P)
	assert isinstance(plunge_point, P)
	assert isinstance(top_surface_z, L)
	assert isinstance(xy_rapid_safe_z, L)

	# Perform any requested tracing:
        if tracing >= 0:
	    pad = ' ' * tracing
	    follows_name = 'None'
	    if isinstance(follows, Operation):
		follows_name = follows._name_get()
	    part_name = part._name_get()

            print(("{0}=>Operation_Dowel_Pin.__init__('{1}', '{2}', {3}, '{4}', '{5}', " +
	      "{6:i}, {7:rpm}, {8:i}, dp={9:i}, pp={10:i}, tz={11:i}, rz={12:i}) ").format(
	      pad, part_name, comment, sub_priority, tool._name_get(), follows_name,
	      feed_speed, spindle_speed, diameter, dowel_point, plunge_point, top_surface_z,
	      xy_rapid_safe_z))

	# Initialize super class:
	Operation.__init__(operation_dowel_pin, "Dowel_Pin", Operation.KIND_DOWEL_PIN,
	  part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed)

	# Load up the rest of *operation_dowel_pin*:
	operation_dowel_pin._diameter = diameter	# Dowel pin diameter
	operation_dowel_pin._dowel_point = dowel_point	# Location to move dowel to
	operation_dowel_pin._plunge_point = plunge_point # Location to plunge dowel down at
	operation_dowel_pin._top_surface_z = top_surface_z
	operation_dowel_pin._xy_rapid_safe_z = xy_rapid_safe_z
	assert isinstance(operation_dowel_pin._top_surface_z, L)

	# Wrap up any *tracing*:
        if tracing >= 0:
            print(("{0}<=Operation_Dowel_Pin.__init__('{1}', '{2}', {3}, '{4}', '{5}', " +
	      "{6:i}, {7:rpm}, {8:i}, dp={9:i}, pp={10:i}, tz={11:i}, rz={12:i}) ").format(
	      pad, part_name, comment, sub_priority, tool._name_get(), follows_name,
	      feed_speed, spindle_speed, diameter, dowel_point, plunge_point, top_surface_z,
	      xy_rapid_safe_z))

    def _cnc_generate(self, tracing):
	""" *Operation_Dowel_Pin*: Generate the CNC G-code for a an
	    *Operation_Dowel_Pin* object (i.e. *self*.)
	"""

	# Use *dowel_pin* instead of *self*.
	operation_dowel_pin = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any *tracing*:
	if tracing >= 0:
	    print("{0}=>Operation_Dowel_Pin._cnc_generate()".format(' ' * tracing))

	# Grab the *part*, *shop*, *code*, *vice*, and *jaw_volume* from *operation_dowel_pin*:
	part = operation_dowel_pin._part_get()
	shop = part._shop_get()
	code = shop._code_get()
	vice = shop._vice_get()
	jaw_volume = vice._jaw_volume_get()
	    
	# Grab some values out of *operation_dowel_pin*:
	comment = operation_dowel_pin._comment
	diameter = operation_dowel_pin._diameter
	cnc_dowel_point = operation_dowel_pin._dowel_point
	cnc_plunge_point = operation_dowel_pin._plunge_point
	cnc_top_surface_z = operation_dowel_pin._top_surface_z
	cnc_xy_rapid_safe_z = operation_dowel_pin._xy_rapid_safe_z

	# Define the speed and feed for these operations:
	ipm10 = Speed(in_per_sec=10.0)
	rpm0 = Hertz()

	# Output the *operation_dowel_pin* comment supplied by the user:
	code._line_comment(comment)

	# Output some dimenions for G-code debugging purposes:
	#code._line_comment(
	#  "{0} Initial Dimensions: {1} x {2} x {3}".format(part._name_get(),
	#  part._dx_original_get(), part._dy_original_get(), part._dz_original_get()))

	#call d@(form@("part:%v% cnc_gen@Op_Dowel_Pin:%i% %i%\n\") %
	#  f@(part.name) % f@(plunge_x) / f@(plunge_y))

	# Rapid over to the plunge point:
	code._xy_rapid_safe_z_force(ipm10, rpm0)
	code._xy_rapid(cnc_plunge_point.x, cnc_plunge_point.y)

	# Now pause to let operator see if Z-safe is at the right height:
	code._command_begin()
	code._unsigned("M", 6)
	code._unsigned("T", 9)
	code._comment("Operator may check that Z-safe is correct")
	code._command_end()

	# Output some information about the *dowel_pin* for G-code debugging:
	#code._line_comment(
	#  "z_stop={0:i} tip_depth={1:i}".format(z_stop, tip_depth))

	# Move slowly down to *z_stop*:
	code._z_feed(ipm10, rpm0, cnc_plunge_point.z, "dowel_pin")

	# Move slowly to (*dowel_x*, *dowel_y*).  This may cause the material in the
	# vice to slide over:
	code._xy_feed(ipm10, rpm0, cnc_dowel_point.x, cnc_dowel_point.y)

	# Now pause again, to let the operator move piece up against the
	# the dowel pin (if it is not already up against there):
	code._command_begin()
	code._unsigned("M", 6)
	code._unsigned("T", 2)
	code._comment("Operator should place part against dowel pin")
	code._command_end()

	# Slowly retract away from the part edge back to *plunge_x* and get
	# back up to Z safe:
	code._xy_feed(ipm10, rpm0, cnc_plunge_point.x, cnc_plunge_point.y)
	code._xy_rapid_safe_z_force(ipm10, rpm0)

	# Wrap-up any *tracing*:
	if tracing >= 0:
	    print("{0}<=Operation_Dowel_Pin._cnc_generate()".format(' ' * tracing))

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

    def _diameter_get(self):
	""" *Operation_Dowel_Pin*: Return the diameter of the *Operation_Dowel_Pin* object
	    (i.e. *self*).
	"""

	return self._diameter

    def _dowel_point_get(self):
	""" *Operation_Dowel_Pin*: Return the dowel push point for the *Operation_Dowel_Pin*
	    object (i.e. *self*).
	"""

	return self._dowel_point


    def _plunge_point_get(self):
	""" *Operation_Dowel_Pin*: Return the plunge point for the *Operation_Dowel_Pin*
	    object (i.e. *self*).
	"""

	return self._plunge_point

    def _top_surface_z_get(self):
	""" *Operation_Dowel_Pin*: Return the top surface Z for the *Operation_Dowel_Pin*
	    object (i.e. *self*).
	"""

	return self._top_surface_z

    def _xy_rapid_safe_z_get(self):
	""" *Operation_Dowel_Pin*: Return the top surface XY rapid safe Z for the
	    *Operation_Dowel_Pin* object (i.e. *self*).
	"""

	return self._xy_rapid_safe_z

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

    def __init__(self, part, comment, sub_priority, tool, order, follows,
      diameter, hole_kind, start, stop, is_countersink):
	""" *Operation_Drill*: Initialize *Operation_Drill* to contain
	    *diameter*, *hole_kind*, *start*, *stop*, and *counter_sink*.
	"""

	# Initialize super class:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(order, int)
	assert isinstance(follows, Operation) or follows == None
	assert isinstance(diameter, L)
	assert isinstance(hole_kind, int)
	assert isinstance(start, P)
	assert isinstance(stop, P)
	assert isinstance(is_countersink, bool)

	# Initialize the super class:
	Operation.__init__(self, "Drill", Operation.KIND_DRILL, part, comment, sub_priority,
	  tool, order, follows, tool._feed_speed_get(), tool._spindle_speed_get())

	# Load up *self*:
	self._diameter = diameter
	self._hole_kind = hole_kind
	self._start = start
	self._stop = stop
	self._is_countersink = is_countersink

    def _cnc_generate(self, tracing):
	""" *Operation_Drill*: Generate the CNC G-code for an
	    *Operation_Drill* object (i.e. *self*).
	"""

	# Use *drill* instead of *self*:
	operation_drill = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Grab some values from *operation_drill*:
	comment = operation_drill._comment
	diameter = operation_drill._diameter
	feed_speed = operation_drill._feed_speed
	hole_kind = operation_drill._hole_kind
	kind = operation_drill._kind
	part = operation_drill._part
	spindle_speed = operation_drill._spindle_speed
	start = operation_drill._start
	stop = operation_drill._stop
	tool = operation_drill._tool

	# Get *shop* and *code*:
	cnc_transform = part._cnc_transform_get()
	shop = part._shop_get()

	code = shop._code_get()

	radius = diameter / 2
	half_radius = radius / 2

	cnc_start = cnc_transform * start
	cnc_stop = cnc_transform * stop


	is_laser = tool._is_laser_get()
	if is_laser:
	    code._dxf_circle(x, y, diameter)
	else:
	    z_start = cnc_start.z
	    z_stop = cnc_stop.z
	    if hole_kind == Part.HOLE_THROUGH:
		assert kind == Operation.KIND_DRILL
		tool_drill = tool
		point_angle = tool_drill._point_angle_get()
		tip_depth = tool_drill._tip_depth_get()
		#call line_comment@(code,
		#  form@("z_stop=%i% diameter=%i% tip_depth=%i%") %
		#  f@(z_stop) % f@(diameter) / f@(tip_depth))
		z_stop -= tip_depth + L(inch=0.040)
	    elif hole_kind == HOLE_KIND_TIP:
		z_stop = cnc_stop.z
	    elif hole_kind == HOLE_KIND_FLAT:
		assert False, "Flat holes can't be done with point drills"
	    else:
		assert False, "Unknown hole kind"

	    # Force the tool to be start from a safe location:
	    code._xy_rapid_safe_z_force(feed_speed, spindle_speed)
	    code._xy_rapid(cnc_start.x, cnc_start.y)
	    code._z_feed(feed_speed, spindle_speed, z_stop, "drill")
	    code._xy_rapid_safe_z_force(feed_speed, spindle_speed)

	    #diameter_divisor = 3.0
	    #material = part._material_get()
	    #switch material.named_material
	    #  case plastic
	    #	# Plastic can really clog up the drill flutes, try to short the
	    #	# the drill spirals:
	    #	diameter_divisor = 0.5    

	    
	
	    #depth = z_start - z_stop
	    #trip_depth = diameter / diameter_divisor
	    #if depth > trip_depth:
	    #	# "Deep" hole:

	    #	# Output G73 G98 F# Q# R# S# X# Y# Z# (comment):
	    #	# Output G83 G98 F# Q# R# S# X# Y# Z# (comment):
	    #	# Output O900 [x] [y] [z_safe] [z_start] [z_stop] [z_step]
	    #	#    [z_back] [feed] [speed] (comment):
	    #	# Output O910 [x] [y] [z_safe] [z_start] [z_stop] [z_step]
	    #	#    [z_back] [feed] [speed] (comment):

	    #	z_rapid = part.z_rapid
	    #	depth = z_rapid - z_stop
	    #	pecks = int(depth / trip_depth) + 1
	    #	# Add just a little to make sure {pecks} * {q} > {depth};
	    #	# The drill will *never* go below the Z value:
	    #	q = (depth / float(pecks)) + L(inch=.005)

	    #	subroutine_code = 900
	    #	if diameter <= L(inch="3/32"):
	    #	    # Use full retract pecking:
	    #	    subroutine_code = 910

	    #	code._command_begin()
	    #	code._subroutine_call(subroutine_code)
	    #	code._length("[]", x - code.vice_x)
	    #	code._length("[]", y - code.vice_y)
	    #	code._length("[]", part.z_safe)
	    #	code._length("[]", z_start)
	    #	code._length("[]", z_stop)
	    #	code._length("[]", radius)
	    #	code._length("[]", half_radius)
	    #	code._speed("[]", f)
	    #	code._hertz("[]", s)
	    #	code._command_comment(comment)
	    #	code._command_end()

	    #	#call begin@(code)
	    #	#call mode_motion@(code, 73)
	    #	#call mode_motion@(code, 83)
	    #	#call mode_canned_cycle_return@(code, 98)
	    #	#call speed@(code, "F", f)
	    #	#call length@(code, "Q", code_length@(q))
	    #	#call length@(code, "R1", z_rapid)
	    #	#call hertz@(code, "S", s)
	    #	#call length@(code, "X", code_length@(x))
	    #	#call length@(code, "Y", code_length@(y))
	    #	#call length@(code, "Z1", code_length@(z_stop))
	    #	#call comment@(code, comment)
	    #	#call end@(code)
	    #else:	
	    #	# Regular hole:

	    #	# Output G81 G98 F# R# S# X# Y# Z# (comment):

	    #	code._command_begin()
	    #	code._mode_motion(81)
	    #	code._mode_canned_cycle_return(98, )
	    #	code._speed("F", f)
	    #	code._length("R1", part.z_rapid)
	    #	code._hertz("S", s)
	    #	code._length("X", x)
	    #	code._length("Y", y)
	    #	code._length("Z1", z_stop)
	    #	code._comment(comment)
	    #	code._command_end()

	    #if is_last:
	    #	# Output a G80 to exit canned cycle mode:
	    #	code._command_begin()
	    #	code._mode_motion(80)
	    #	code._comment("End canned cycle")
	    #	code._commend_end()

	    #	# Forget the G98 or G99 code:
	    #	code.g2 = -1

	    # Do any dwell:
	    #cnc_drill_count = part.cnc_drill_count + 1
	    #part.cnc_drill_count = cnc_drill_count
	    #if not drill.is_countersink and \
	    #      cnc_drill_count % part.cnc_drill_pause == 0:
	    #	# Perform a Spindle Stop, Pause, Spindle_Start, Pause:
	    #	code._command_begin()
	    #	code._unsigned("M", 5)
	    #	code._command_end()

	    #	code._command_begin()
	    #	code._unsigned("G11", 4)
	    #	code._time("P", Code_Time(seconds=3.0))
	    #	code._command_end()

	    #	code._command_begin()
	    #	code._unsigned("M", 3)
	    #	code._hertz("S", s)
	    #	code._command_end()

	    #	code._command_begin()
	    #	code._unsigned("G11", 4)
	    #	code._time("P", Code_Time(seconds=8.0))
	    #	code._command_end()

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

    def __init__(self, part, comment, sub_priority, tool, order, follows, diameter,
      countersink_diameter, hole_kind, start, stop, feed_speed, spindle_speed,
      top_surface_transform, tracing):
	""" *Operation_Round_Pocket* will intitialize an
	    *Operation_Round_Pocket* object (i.e. *self*) to contain
	    *part*, *comment*, *sub_priority*, *tool*, *order*, *follows*,
	    *diameter*, *hole_kind*, *start*, *stop*, *feed_speed*, and *spindle_speed*.
	"""

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(sub_priority, int)
	assert isinstance(tool, Tool)
	assert isinstance(order, int)
	assert isinstance(follows, Operation) or follows == None
	assert isinstance(diameter, L)
	assert isinstance(countersink_diameter, L)
	assert isinstance(hole_kind, int)
	assert isinstance(start, P)
	assert isinstance(stop, P)
	assert isinstance(feed_speed, Speed)
	assert isinstance(spindle_speed, Hertz)
	assert isinstance(tracing, int)
	assert isinstance(top_surface_transform, Transform)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Operation_Round_Pocket.__init__('{1}', '{2}', '{3}', {4}, '{5}', {6}," +
	      " *, {7:i} {8:i}, {9}, {10:i}, {11:i}, {12:i}, {13:rpm}, {14:s}, {15})").
	      format(indent, "Round_Pocket", part._name_get(), comment, sub_priority, tool,
	        order, diameter, countersink_diameter, hole_kind, start, stop, feed_speed,
	        spindle_speed, top_surface_transform, tracing))

	# Initialize superclass:
	Operation.__init__(self, "Round_Pocket", Operation.KIND_ROUND_POCKET,
	  part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed)

	# Load up *self*:
	self._diameter = diameter
	self._countersink_diameter = countersink_diameter
	self.hole_kind = hole_kind
	self._start = start
	self._stop = stop
	self._tracing = -1000000
	#self._tracing = tracing
	self._top_surface_transform = top_surface_transform

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}<=Operation_Round_Pocket.__init__('{1}', '{2}', '{3}', {4}, '{5}', {6}," +
	      " *, {7:i} {8:i}, {9}, {10:i}, {11:i}, {12:i}, {13:rpm}, {14:s})").format(indent,
	      "Round_Pocket", part._name_get(), comment, sub_priority, tool, order, diameter,
	      countersink_diameter, hole_kind, start, stop, feed_speed, spindle_speed,
	      top_surface_transform))

    def _cnc_generate(self, tracing):
	""" *Operation_Round_Pocket*: Generate the CNC G-code for an *Operation_Round_Pocket*
	    object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(tracing, int)

	# Use *round_pocket* instead of *self*:
	round_pocket = self

	# Perform any requested *tracing*:
	if tracing < 0:
	    tracing = round_pocket._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation_Round_Pocket._cnc_generate('{1}')".
	     format(indent, round_pocket._name))

	# Extract some values from *round_pocket*:
	part = round_pocket._part
	diameter = round_pocket._diameter
	start = round_pocket._start
	stop = round_pocket._stop
	tool = round_pocket._tool

	# The *top_surface_transform* has been previously set orient the material correctly for CNC:
	top_surface_transform = round_pocket._top_surface_transform
	mapped_start = top_surface_transform * start
	mapped_stop = top_surface_transform * stop

	# Compute some values based on *diameter*:
	maximum_depth = diameter / 3.0
	radius = diameter / 2

	# Extract some values from *part* and *shop*:
	shop = part._shop
	code = shop._code_get()

	# Extract some values from *tool*:
	tool_diameter = tool._diameter_get()
	comment = round_pocket._comment
	f = tool._feed_speed_get()
	s = tool._spindle_speed_get()
	hole_kind = round_pocket.hole_kind

	# Figure out if *tool* is a laser:
	is_laser = tool._is_laser_get()

	# Compute some values based on {tool_diameter}:
	tool_radius = tool_diameter / 2 
	half_tool_radius = tool_radius / 2

	if is_laser:
	    # We just cut a simple circle:
	    if tracing >= 0:
		print("{0}start={1:i} stop={2:i}".format(indent, mapped_start, mapped_stop))
	    code._dxf_circle(mapped_start.x, mapped_start.y, radius - tool_radius)
	else:
	    # We do all the work to mill out the round_pocket pocket;

	    # Deal with through holes:
	    is_through = False
	    #print("operation.comment='{0}'".format(self._comment))
	    x = mapped_start.x
	    y = mapped_start.y
	    z_start = mapped_start.z
	    z_stop = mapped_stop.z
	    if hole_kind == Part.HOLE_THROUGH:
		is_through = True
		z_stop -= L(inch=0.025)

	    code._comment(
	      "z_start={0} z_stop={1}".format(z_start, z_stop))
	    
	    z_depth = (start - stop).length()
	    passes = int(z_depth / maximum_depth) + 1
	    depth_per_pass = z_depth / float(passes)
	    assert passes < 100

	    #call line_comment@(code,
	    #  form@("z_depth=%i% passes=%i% depth_per_pass=%i%") %
	    #  f@(z_depth) % f@(passes) / f@(depth_per_pass))

	    # Move to position:
	    code._comment(comment)
	    code._z_safe_assert("round_pocket_pocket", comment)
	    code._xy_rapid(x, y)

	    z_feed = f / 4.0
	    shave = L(inch=0.005)
	    radius_remove = radius - shave - tool_radius
	    for depth_pass in range(passes):
		code._comment(
		  "{0} round_pocket pocket [pass {1} of {2}]".
		  format(comment, depth_pass + 1, passes))
		
		# Get to proper depth:
		z = z_start - depth_per_pass * float(depth_pass + 1)

		if is_through:
		    code._xy_feed(f, s, x, y + radius_remove)
		    code._z_feed(z_feed, s, z, "round_pocket_pocket")
		    code._xy_ccw_feed(f, s, radius_remove, x, y - radius_remove)
		    code._xy_ccw_feed(f, s, radius_remove, x, y + radius_remove)
		else:
 		    # We have to mow out all the intervening space:
		    radius_passes = int(radius_remove /  half_tool_radius) + 1
		    pass_remove_delta = radius_remove / float(radius_passes)

		    for radius_index in range(radius_passes):
			pass_remove = pass_remove_delta * float(radius_index + 1)
			code._xy_feed(f, s, x, y + pass_remove)
			code._z_feed(z_feed, s, z, "round_pocket_pocket")
			code._xy_ccw_circle(f, s, pass_remove, x, y - pass_remove)
			code._xy_ccw_circle(f, s, pass_remove, x, y + pass_remove)

	    # Do a "spring pass" to make everybody happy:
	    code._comment("{0} round_pocket pocket 'spring' pass".format(comment))
	    path_radius = radius - tool_radius
	    half_path_radius = path_radius / 2
	    code._xy_feed(f, s, x, y)

	    # Carefully feed the tool to the edge:
	    code._comment("Carefully feed tool to the edge of the hole")
	    code._xy_ccw_feed(f, s, half_path_radius, x + half_path_radius, y + half_path_radius)
	    code._xy_ccw_feed(f, s, half_path_radius, x,                    y + path_radius)

	    # Peform the entire spring cut:
	    code._comment("Peform the entire spring cut")
	    code._xy_ccw_feed(f, s, path_radius, x, y - path_radius)
	    code._xy_ccw_feed(f, s, path_radius, x, y + path_radius)

	    # Carefully remove the tool back to the center:
	    code._comment("Carefully remove the tool back to the center")
	    code._xy_ccw_feed(f, s, half_path_radius, x - half_path_radius, y + half_path_radius)
	    code._xy_ccw_feed(f, s, half_path_radius, x, y)

	    # Safely retract to z safe:
	    code._comment("Safely retract to z safe")
	    code._z_safe_retract(z_feed, s)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Operation_Round_Pocket._cnc_generate('{1}')".
	      format(indent, round_pocket._name))

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

    def _cnc_generate(self, tracing):
	""" *Operation_Simple_Exterior*: Generate the CNC G-code for an
	    *Operation_Simple_Exterior* object (i.e. self).
	"""

	# Verify argument types:
	assert isinstance(tracing, int)

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
      corner1, corner2, corner_radius, tool_radius, pocket_kind):
	""" *Operation_Simple_Pocket*: Initialize an *Operation_Simple_Pocket*
	    object (i.e. *self*) to contain *part*, *comment*, *sub_priority*,
	    *tool*, *order*, *follows*, *corner1*, *corner2*, *corner_radius*,
	    *tool_radius*, and *pocket_kind*.
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
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)
	assert isinstance(corner_radius, L)
	assert isinstance(tool_radius, L)
	assert isinstance(pocket_kind, int)

	# Initialize superclass:
	operation_kind = Operation.KIND_SIMPLE_POCKET
	Operation.__init__(self, "Simple_Pocket", operation_kind,
	  part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed)

	assert corner1.z < corner2.z

	# Load up the rest of *self*:
	self._corner1 = corner1
	self._corner2 = corner2
	self._corner_radius = corner_radius
	self._tool_radius = tool_radius
	self._pocket_kind = pocket_kind

    def _corner_radius_get(self):
	""" *Operation_Simple_Pocket: Return the corner radius for the *Operation_Simple_Pocket*
	    object (i.e. *self*).
	"""

	return self._corner_radius

    def _cnc_generate(self, tracing = -1000000):
	""" *Operation_Simple_Pocket*: Generate the CNC G-code for a
	    *Operation_Simple_Pocket* object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(tracing, int)

	# Use *pocket* instead of *self*:
	pocket = self

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation_Simple_Pocket._cnc_generate('{1}')".
	      format(indent, pocket._comment))

	# Grap some values from *pocket*:
	corner1 = pocket._corner1
	corner2 = pocket._corner2

	tool_radius = pocket._tool_radius
	corner_radius = pocket._corner_radius
	pocket_kind = pocket._pocket_kind
	comment = pocket._comment
	
	# Extract the X/Y/Z coordinates from *corner1* and *corner2*:
	x1 = corner1.x
	y1 = corner1.y
	z1 = corner1.z
	x2 = corner2.x
	y2 = corner2.y
	z2 = corner2.z

	# Determine *z_top* and *zbottom Z*:
        z_top = z1.maximum(z2)
        z_bottom = z1.minimum(z2)

	# Now create *start_point* and *end_point*:
	start_point = P(x1, y1, z_top)
	end_point = P(x1, y1, z_bottom)

	# Now we can create *top_surface_transform* from *start_point* and *end_point*:
	degrees0 = Angle()
	top_surface_transform = \
	  Transform.top_surface("simple_pocket", start_point, end_point, degrees0, tracing + 1)

	transformed_corner1 = top_surface_transform * corner1
	transformed_corner2 = top_surface_transform * corner1

	# Unpack from the corners:
	x1 = corner1.x
	y1 = corner1.y
	z_stop = corner1.z
	x2 = corner2.x
	y2 = corner2.y
	z_start = corner2.z
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
	if tracing >= 0:
            print("{0}is_laser={1}".format(indent, is_laser))

	# Compute the number of depth passes required:
	maximum_operation_depth = tool_radius / 2
	if is_laser:
	    maximum_operation_depth = tool._maximum_z_depth

	total_cut = z_start - z_stop + z_extra
	if tracing >= 0:
	    print("{0}z_start={1:i} z_stop={2:i} z_extra={3:i}".
	      format(indent, z_start, z_stop, z_extra))
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

	is_laser = isinstance(tool, Tool_End_Mill) and tool._is_laser_get()
	if is_laser:
	    code._simple_pocket_helper(pocket, zero, s, f, zero, True)
	else:
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
		    #remaining = remaining - (r + r)
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

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Operation_Simple_Pocket._cnc_generate('{1}')".
	      format(indent, pocket._comment))

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

    def _cnc_generate(self, tracing):
	""" *Operation_Vertical_Lathe*:
	"""

	# Verify argument types:
	assert isinstance(tracing, int)

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

class Parallels:
    """ *Parallels*: A *Parallels* object specifies a set of parallels for use with a *Vice*."""

    def __init__(self, length, thickness, heights):
	""" *Parallels*:  """

	# Use *parallels* instead of *self*:
	parallels = self

        # Verify argument types:
	assert isinstance(length, L)
	assert isinstance(thickness, L)
	assert isinstance(heights, tuple) or isinstance(heights, list)
	zero = L()
	for height in heights:
	    assert isinstance(height, L) and height > zero 

        # Stuff values into *parallels*:
        parallels._length = length
        parallels._thickness = thickness
        parallels._heights = tuple(heights)

    def _length_get(self):
        """ *Parallels: Return the length of each paralle associated with 
	    *Parallels* object (i.e. *self*.) """

	return self._length

    def _heights_get(self):
        """ *Parallels: Return the tuple of heights of the for all the parallels
	    associated with the *Parallels* object (i.e. *self*.) """

	return self._heights

    def _thickness_get(self):
        """ *Parallels: Return the thickness of each parallel associated with the
	    *Parallels* object (i.e. *self*.) """

	return self._heights

class Part:
    """ A *Part* specifies either an assembly of parts or a single physical part. """

    HOLE_THROUGH = 1
    HOLE_TIP = 2
    HOLE_FLAT = 3

    # Flavors of values that can be stored in a {Part}:
    def __init__(self, up, name):
	""" *Part*: Initialize *self* to have a parent of *up*. """

	# Use *part* instead of *self*:
	part = self

        # Check argument types:
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)

	# Some useful abbreviations:
	zero = L()
	one = L(mm=1.0)
	z_axis = P(zero, zero, one)

	ezcad = EZCAD3.ezcad
	part._axis = z_axis
	part._bounding_box = Bounding_Box()
	part._color = None
	part._center = P()
	part._cnc_suppress = False
	part._cnc_transform = None		# Orientation and translate transfrom for CNC vice
	part._dowel_x = zero
	part._dowel_y = zero
	part._dx_original = zero
	part._dxf_x_offset = zero
	part._dxf_y_offset = zero
	part._dy_original = zero
	part._dz_original = zero
	part._dxf_scad_lines = []
	part._ezcad = ezcad
	part._is_part = False	
	part._laser_preferred = False
	part._material = Material("aluminum", "")
	part._name = name
	part._operations = []
	part._plunge_x = zero
	part._plunge_y = zero
	part._position = Matrix.identity()	# Scheduled to go away
	part._position_count = 0
	part._priority = 0
	part._projection_axis = z_axis		# Scheduled to go away
	#part._places = {}
	part._reposition = Matrix()		# Scheduled to go away
	part._rotate = None
	part._shop = ezcad._shop
	part._scad_difference_lines = []
	part._scad_union_lines = []
	part._signature_hash = None
	part._stl_file_name = None
	part._top_surface_transform = Transform()
	part._top_surface_set = False		# Scheduled to go away
	part._tool_preferred = ""
	part._tracing = -1000000
	part._translate = None
	part._visible = True
	part._vice_x = None			# Vice X origin relative to part X origin
	part._vice_y = None			# Vice Y origin relative to part Y origin
	part._z_rapid = None			# Z above which XY rapids are OK
	part._z_safe = None			# Z above which Z rapids are OK
	part.up = up

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

    def cnc_suppress(self):
        """ *Part*: Suppress CNC generation for the *Part* object (i.e. *self*.) """

	# Use *part* insted of *self*
	part = self

	# Mark *part* for CNC suppression:
        part._cnc_suppress = True

    def _cnc_transform_get(self):
        """ *Part*: Return the CNC transform associated with the *Part* object (i.e. *self*.) """

	return self._cnc_transform

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

    def _cnc_flush(self, program_number, tracing):
	""" *Part*: Flush out the CNC code for the *Part* object (i.e. *self*).
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(program_number, int)
	assert isinstance(tracing, int)

	# *part_name* is used in tracing and other places:
	part_name = part._name

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._cnc_flush('{1}', prog_number={2})".
	      format(indent, part_name, program_number))
	    trace_detail = 1

	original_program_number = program_number

	#call show@(part, "before flush")
	shop = part._shop_get()
	assert isinstance(shop, Shop)

	operations = part._operations_get()
	size = len(operations)

	#call show@(part, "before sort", 1t)

	# FIXME move the sort into *operations_regroup*:
	# Sort *operations* to group similar operations together:
	operations.sort(cmp=Operation._compare)

	if trace_detail >= 2:
	    print("{0}len(operations)={1}".format(indent, len(operations)))
	    for operation in operations:
		print("{0}operation name:{1}".format(indent, operation._name_get()))

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
	#print("operation_groups=", operation_groups)

	# Open the top-level *part_ngc_file* that lists both the tool table and
        # calls each tool operation from a single top level .ngc file:
	ezcad = part._ezcad_get()
	ngc_directory = ezcad._ngc_directory_get()
	part_ngc_file_name = os.path.join(ngc_directory, "O{0}.ngc".format(program_number))
	part_ngc_file = open(part_ngc_file_name, "w")
	assert part_ngc_file != None, "Unable to open {0}".format(part_ngc_file_name)
	if trace_detail >= 1:
	    print("{0}part_ngc_file_name='{1}' opened".format(indent, part_ngc_file_name))

	# Output some heading lines for *part_ngc_file*:
	part_ngc_file.write("( Part: {0})".format(part._name_get()))
	part_ngc_file.write("( Tooling table: )\n")

	# Open the top-level *part_wrl_file* file that outputs each visual tool path.
        # Note that these file names ene in .wrl and are written into the .../ngc sub-directory:
	part_wrl_file_name = os.path.join(ngc_directory, "O{0}.wrl".format(program_number))
	part_wrl_file = open(part_wrl_file_name, "w")

	# Output the preceeding VRML code to group group the children shapes:
        part_wrl_file.write("#VRML V2.0 utf8\n")
        #part_wrl_file.write("Viewpoint {description \"Initial view\" position 0 0 150\n")
	#part_wrl_file.write("NavigationInfo { type \"Examine\" }\n")
	#part_wrl_file.write("Def x{0} Group {1}\n".format(part_name, "{"))
	part_wrl_file.write("Group {\n")
	part_wrl_file.write(" children [\n".format(part_name))

	#FIXME: This code belongs in a private routine:
	# Output a coordinate axis to *part_wrl_file*:
	part_wrl_file.write("  #VRML V2.0 utf8\n")
        part_wrl_file.write("  Shape {\n")
	part_wrl_file.write("   geometry IndexedLineSet {\n")
	part_wrl_file.write("    colorPerVertex FALSE\n")
	part_wrl_file.write("    color Color {\n")
	part_wrl_file.write("     color [\n")
	part_wrl_file.write("      1.0 0.0 0.0 # red\n")
	part_wrl_file.write("      0.0 1.0 0.0 # green\n")
	part_wrl_file.write("      0.0 0.0 1.0 # bule\n")
	part_wrl_file.write("     ]\n")
	part_wrl_file.write("    }\n")
	part_wrl_file.write("    coord Coordinate {\n")
	part_wrl_file.write("     point [\n")
	part_wrl_file.write("      0.0 0.0 0.0\n")	# [0]: Origin
	part_wrl_file.write("      25.4 0.0 0.0\n")	# [1]: X axis end-point
	part_wrl_file.write("      0.0 25.4 0.0\n")	# [2]: Y axis end-pont
	part_wrl_file.write("      0.0 0.0 25.4\n")	# [3]: Z axis end-pont
	part_wrl_file.write("      23.0 0.0 2.0\n")	# [4]: X axis arrow head
	part_wrl_file.write("      23.0 0.0 -2.0\n")	# [5]: X axis arrow head
	part_wrl_file.write("      0.0 23.0 2.0\n")	# [6]: Y axis arrow head
	part_wrl_file.write("      0.0 23.0 -2.0\n")	# [7]: Y axis arrow head
	part_wrl_file.write("      2.0 0.0 23.0\n")	# [6]: Z axis arrow head
	part_wrl_file.write("      -2.0 0.0 23.0\n")	# [7]: Z axis arrow head
	part_wrl_file.write("     ]\n")
	part_wrl_file.write("    }\n")
	part_wrl_file.write("    coordIndex [\n")
	part_wrl_file.write("     0 1 -1\n")		# X Axis
	part_wrl_file.write("     1 4 -1\n")		# X Axis Arrow Head
	part_wrl_file.write("     1 5 -1\n")		# X Axis Arrow Head
	part_wrl_file.write("     0 2 -1\n")		# Y Axis
	part_wrl_file.write("     2 6 -1\n")		# Y Axis Arrow Head
	part_wrl_file.write("     2 7 -1\n")		# Y Axis Arrow Head
	part_wrl_file.write("     0 3 -1\n")		# Z Axis
	part_wrl_file.write("     3 8 -1\n")		# Z Axis Arrow Head
	part_wrl_file.write("     3 9 -1\n")		# Z Axis Arrow Head
	part_wrl_file.write("    ]\n")
	part_wrl_file.write("    colorIndex [\n")
	part_wrl_file.write("     0\n")			# X Axis
	part_wrl_file.write("     0\n")			# X Axis Arrow Head
	part_wrl_file.write("     0\n")			# X Axis Arrow Head
	part_wrl_file.write("     1\n")			# Y Axis
	part_wrl_file.write("     1\n")			# Y Axis Arrow Head
	part_wrl_file.write("     1\n")			# Y Axis Arrow Head
	part_wrl_file.write("     2\n")			# Z Axis
	part_wrl_file.write("     2\n")			# Z Axis Arrow Head
	part_wrl_file.write("     2\n")			# Z Axis Arrow Head
	part_wrl_file.write("    ]\n")
	part_wrl_file.write("   }\n")
        part_wrl_file.write("  }\n")

	# For now the first operation must always be an *operation_dowel_pin*:
        assert len(operations) > 0, "No operations present"
	operation_dowel_pin = operations[0]
        assert isinstance(operation_dowel_pin, Operation_Dowel_Pin), \
	  "Got operation '{0}' instead of dowel pin".format(operation.__class__.__name__)

	# Grab some values from *operation_dowel_pin*:
	comment = operation_dowel_pin._comment_get()
	top_surface_z = operation_dowel_pin._top_surface_z_get()
	xy_rapid_safe_z = operation_dowel_pin._xy_rapid_safe_z_get()
	feed_speed = operation_dowel_pin._feed_speed_get()
	spindle_speed = operation_dowel_pin._spindle_speed_get()

	# Now visit each *operation_group* in *operation_groups*:
	for index, operation_group in enumerate(operation_groups):
	    # Compute the *ngc_program_number* for 
	    ngc_program_number = program_number + 1 + index

	    # Grab the *code* object and start the code generation:
	    code = shop._code_get()
	    assert isinstance(code, Code)
	    code._start(part, tool, ngc_program_number, spindle_speed,
	      part_wrl_file, part._stl_file_name, tracing + 1)
	    code._dxf_xy_offset_set(part._dxf_x_offset_get(), part._dxf_y_offset_get())
	    code._top_surface_z_set(top_surface_z)
	    code._xy_rapid_safe_z_set(xy_rapid_safe_z)
	    code._xy_rapid_safe_z_force(feed_speed, spindle_speed)

	    # Output a couple of G-code lines *part_ngc_file*:
	    tool_number = tool._number_get()
	    tool_name = tool._name_get()
	    part_ngc_file.write("( T{0} {1} )\n".format(tool_number, tool_name))
	    part_ngc_file.write("O{0} call\n".format(ngc_program_number))

	    # Sweep through the *operation_group*:
	    part._cnc_flush_helper(operation_group, ngc_program_number, part_wrl_file, tracing + 1)

	    # Close off *code*:
	    code._finish(tracing + 1)

	# Write out the final G-code lines to *part_ngc_file*:
	part_ngc_file.write("G53 Y0.0 ( Move the work to the front )\n")
	part_ngc_file.write("M2\n")
	part_ngc_file.close()

	# Write the part out to the *part_wrl_file*:
	cnc_transform = part._cnc_transform
	part._wrl_write(part_wrl_file, part._cnc_transform, 2,
	  file_name=part._stl_file_name, tracing=tracing + 1)

	# Close out the VRML for *part_wrl_file* and then close it:
	part_wrl_file.write(" ]\n")
	part_wrl_file.write("}\n")
	part_wrl_file.close()


	# Empty out *operations* :
	del operations[:]

	# Compute the next *program number* to be a the next multiple of 10 and return it:
	new_program_number = program_number + len(operation_groups) + 1
	new_program_number = (new_program_number + 9) / 10 * 10

	if tracing >= 0:
	    print("{0}<=Part._cnc__flush('{1}', prog_no={2}) =>{3}".
	      format(indent, part_name, program_number, new_program_number))

	return new_program_number

    def _cnc_flush_helper(self, operations, ngc_program_number, part_wrl_file, tracing):
	""" *Part*: Output the G-code for *operations* to a "On.ngc" file where,
	    N is the *ngc_program_number* using the *Part* object (i.e. *self*).
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(operations, list) and len(operations) > 0
	for operation in operations:
	    assert isinstance(operation, Operation)
	assert isinstance(ngc_program_number, int)
	assert isinstance(part_wrl_file, file)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	part_name = part._name
	indent = ""
	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 1
	    indent = ' ' * tracing
	    debug = True
	    print("{0}=>Part._cnc_flush_helper(part='{1}', len(ops)={2}, ngc_no={3}, *)".
	      format(indent, part_name, len(operations), ngc_program_number))

	zero = L()

	# Grab some values from *part*:
	shop = part._shop_get()

	# Grab the first *operation*:
	operation = operations[0]
	#vice_x = operation._vice_x_get()
	#vice_y = operation._vice_y_get()
	tool = operation._tool_get()

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
	if tracing >= 0:
            print("{0}all_operations_are_drills={1}".format(indent, all_operations_are_drills))

	# FIXME: Move this code to *Operaton_Drill* class!!!
	# Reorder the drill operations to minimize traverses:
	if False and all_operations_are_drills:
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
		    if tracing_detail > 0:
			print("{0}minimum_distance[{1}]={0}\n".
			  format(indent, index, minimum_distance))

			# Swap the two drill operations:
			match_operation = operations[match_index]
			operations[match_index] = operations[index + 1]
			operations[index + 1] = match_operation

	# Now do a cnc generate for each *operaton* in *operations*:
	code = None
	for operation in operations:
	    # Perform any requested *tracing_detail*:
	    if tracing_detail > 0:
		print("{0}operation name:{1} tool={2}".format(indent,
		  operation._name_get(), operation._tool_get()._name_get()))

	    # Grab the *spindle_speed*:
	    spindle_speed = operation._spindle_speed_get()
	    feed_speed = operation._feed_speed_get()

	    # Perform the CNC generation step for *operation*:
	    operation._cnc_generate(tracing + 1)


	if code != None:
	    code._xy_rapid_safe_z_force(feed, speed)
	    code._dxf_xy_offset_set(zero, zero)
	    code._finish()

	# Perform any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part._cnc_flush_helper(part='{1}', len(ops)={2}, ngc_no={3}, *)".
	      format(indent, part_name, len(operations), ngc_program_number))

    def _dimensions_update(self, ezcad, tracing):
	""" *Part*: Update the dimensions of the *Part* object (i.e. *self*)
	    and all of its children *Part*'s.
	"""

	# Verify argument types:
	assert isinstance(ezcad, EZCAD3)
	assert isinstance(tracing, int)

	# Use *part* instead of *self*:
	part = self

	# Do any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * trace
	    print("{0}=>Part._dimensions_update('{1}', *)". \
	      format(indent, part._name))
	    #print("{0}Part._dimensions_update:places={1}". \
	    #  format(' ' * trace, part._places))

	# Start with nothing *changed*:
	changed = 0

	#FIXME: This is weird, why is this necessary???
	part._ezcad = ezcad

	# First record the current values load into *part*.
	sub_parts = []
	before_values = {}
	name = part._name
	for attribute_name in dir(part):
	    attribute = getattr(part, attribute_name)
	    if attribute_name.startswith("__"):
		pass
	    elif attribute_name.endswith("_"):
		assert isinstance(attribute, Part), \
		  "{0}.{1} is not a Part".format(name, attribute_name)
		changed += attribute._dimensions_update(ezcad, tracing + 1)
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

	# Reset the bounding box for *part*:
	before_bounding_box = part._bounding_box
	part._bounding_box = Bounding_Box()

	# Peform dimension updating for *part*:
	
	#part._tracing = tracing + 1
	part.construct()

	# See if anything changed:
	for attribute_name in before_values.keys():
	    before_value = before_values[attribute_name]
	    after_value = getattr(part, attribute_name)
	    assert type(before_value) == type(after_value), \
	      "{0}.{1} before type ({2}) different from after type ({3})". \
	      format(name, attribute_name, type(before_value),
	      type(after_value))
	    if before_value != after_value:
		changed += 1
		if changed > 0:
		    #print("here 2")
		    pass
		if ezcad._update_count > 10 or tracing >= 0:
		    print("{0}Part._dimensions_update:{1}.{2} ({3}=>{4})". \
		      format(' ' * trace, name, attribute_name,
		      before_value, after_value))

	# Now update the *bounding_box* for *part*:
	after_bounding_box = part._bounding_box
	for sub_part in sub_parts:
	    after_bounding_box.bounding_box_expand(sub_part._bounding_box)
	if before_bounding_box != after_bounding_box:
	    changed += 1
	    part._bounding_box = after_bounding_box

	    #part._bounding_box_set(after_bounding_box)

	# Update bounding box with placed *Part* bounding boxes:
	# for place in part._places.values():
	#for sub_part in sub_parts:
	#    # Grab some values from *place* and *part*:
	#    #place_part = place._part
	#    #place_name = place._name
	#    place = Place(center = part._center, axis = part._axis,
	#      rotate = part._rotate, translate = part._translate)
	#    forward_matrix = place._forward_matrix
	#
	#    if trace >= 0:
	#	print("{0}Part._dimensions_update:place={1}". \
	#	  format(' ' * trace, place))
	#	#print("{0}Part._dimensions_update:Merge {1:m} into {2:m}". \
	#	#  format(' ' * trace, place_box, box))
	#
	#    part._box_point_update("[TNE]",
	#      forward_matrix.point_multiply(sub_part.tne))
	#    part._box_point_update("[TNW]",
	#      forward_matrix.point_multiply(sub_part.tnw))
	#    part._box_point_update("[TSE]",
	#      forward_matrix.point_multiply(sub_part.tse))
	#    part._box_point_update("[TSW]",
	#      forward_matrix.point_multiply(sub_part.tsw))
	#    part._box_point_update("BNE]",
	#      forward_matrix.point_multiply(sub_part.bne))
	#    part._box_point_update("[BNW]",
	#      forward_matrix.point_multiply(sub_part.bnw))
	#    part._box_point_update("[BSE]",
	#      forward_matrix.point_multiply(sub_part.bse))
	#    part._box_point_update("[BSW]",
	#      forward_matrix.point_multiply(sub_part.bsw))

	# Determine whether *box* has changed:
	#after_box_changed_count = part._box_changed_count
	#if before_box_changed_count != after_box_changed_count:
	#    if ezcad._update_count > 10:
	#	print("{0} bounding box changed".format(name))
	#    changed += 1

	if tracing >= 0:
	    print("{0}<=Part._dimensions_update('{1}')=>{2}". \
	      format(indent, part._name, changed))

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

    def _dowel_x_get(self):
	""" *Part*: Return the edge X field of the *Part* object (i.e. *self*)
	"""

	return self._dowel_x

    def _dowel_y_get(self):
	""" *Part*: Return the edge Y field of the *Part* object (i.e. *self*)
	"""

	return self._dowel_y

    def _ezcad_get(self):
	""" *Part*: Return the *EZCAD* object from the *Part* object (i.e. *self*)
	"""

	return self._ezcad

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

    def _manufacture(self, ezcad, tracing):
	""" *Part*: Visit the *Part* object (i.e. *self*) and all is the children *Part*'s
	    and perform any manufacturing steps.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(ezcad, EZCAD3)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing < 0 and part._tracing >= 0:
            tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._manufacture:('{1}', *)".format(indent, self._name))

	# Make sure *part* is connected ot *ezcad*:
	part._ezcad = ezcad
	part_name = part._name
	shop = part._shop_get()
	code = shop._code_get()

	# Figure out the *mode_name* and do any requested *tracing*::
	mode = ezcad._mode
	mode_name = "?mode?"
	if mode == EZCAD3.VISUALIZATION_MODE:
	    mode_name = "Visual"
	elif mode == EZCAD3.STL_MODE:
	    mode_name = "STL"
	elif mode == EZCAD3.CNC_MODE:
	    mode_name = "CNC"
	if tracing >= 0:
	    print("{0}mode='{1}'".format(indent, mode_name))

	# First manufacture any child *Part*'s:
	for attribute_name in dir(part):
	    if not attribute_name.startswith("_") and \
	      attribute_name.endswith("_"):
		child_part = getattr(part, attribute_name)
		assert isinstance(child_part, Part), \
		  "{0}.{1} is not a Part".format(part.name, attribute_name)
		child_part._manufacture(ezcad, tracing + 1)

	# Now run construct this *part*
	if tracing >= 0:
	    print("{0}==>Part.construct('{1}')".format(indent, part_name))
	part.construct()
	if tracing >= 0:
	    print("{0}<==Part.construct('{1}')".format(indent, part_name))

	# Do the visualization steps:
	if mode == EZCAD3.VISUALIZATION_MODE:
	    ezcad = part._ezcad_get()
	    wrl_directory = ezcad._wrl_directory_get()
	    wrl_file_name = os.path.join(wrl_directory, "{0}.wrl".format(part._name))
	    wrl_file = open(wrl_file_name, "w")
	    transform = Transform()
	    part._wrl_write(wrl_file, transform, 0, wrl_file_name, tracing = tracing + 1)
	    wrl_file.close()

	# Now generate any CNC files:
	if part._cnc_suppress:
	    print("Part '{0}' is suppressing CNC".format(part._name))
	if mode == EZCAD3.CNC_MODE and not part._cnc_suppress:
	    # Flush out all of the pending CNC operations:
	    if tracing >= 0:
		print("{0}==>Part._manufacture('{1}'):CNC".format(indent, part._name))

	    shop = ezcad._shop
	    program_base = shop._program_base_get()
	    assert program_base % 10 == 0
	    program_number = part._cnc_flush(program_base, tracing + 1)

	    # We want the program base number to start with a mulitple of 10.
	    remainder = program_number % 10
	    if remainder != 0:
		program_number += 10 - remainder
	    shop._program_base_set(program_number)

	    # See if we have a .dxf file to write:
	    if code._dxf_content_avaiable():
		# We do have a .dxf file to write.  Open *dxf_file*:
		dxf_directory = ezcad._dxf_directory_get()
		dxf_file_name = os.path.join(dxf_directory, "{0}.dxf".format(part._name))
		dxf_file = open(dxf_file_name, "w")

		# Output the .dxf file headers:
		#dxf_file.write("0\nSECTION\n2\nHEADER\n")
		#dxf_file.write("9\n$DIMAUNITS\n70\n1\n")
		#dxf_file.write("9\n$INSUNITS\n70\n1\n")
		#dxf_file.write("9\n$LUNITS\n70\n0\n")
		#dxf_file.write("9\n$MEASUREMENT\n70\n0\n")
		#dxf_file.write("0\nENDSEC\n")
		dxf_file.write("0\nSECTION\n2\nENTITIES\n")

		# Output the body of the *dxf_file*:
		code._dxf_write(dxf_file)

		# Close out *dxf_file*:
		dxf_file.write("0\nENDSEC\n0\nEOF\n")
		dxf_file.close()

	    if tracing >= 0:
		print("{0}<==Part._manufacture('{1}'):CNC".format(indent, part._name))

	# Now generate any .stl files:
	if ezcad._mode == EZCAD3.STL_MODE:
	    tracing_detail = -1
	    if tracing >= 0:
		tracing_detail = 3
		print("{0}==>Part._manufacture('{1}'):STL".format(indent, part._name))

	    # Now manufacture this node:
	    #scad_difference_lines = []
	    #scad_union_lines = []
	    #self._scad_difference_lines = scad_difference_lines
	    #self._scad_union_lines = scad_union_lines
	    #self.tracing = tracing + 1
	    #self.construct()
	    #self._scad_difference_lines = None
	    #self._scad_union_lines = None

	    scad_difference_lines = self._scad_difference_lines
	    scad_union_lines = self._scad_union_lines
	    if tracing_detail >= 3:
		print("{0}len(scad_difference_lines)={1} len(scad_union_lines={2}".
		  format(indent, len(scad_difference_lines), len(scad_union_lines)))

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
	    if tracing_detail >= 3:
		print("{0}signature_hash={1}".format(indent, signature_hash))

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

            # Now ensure see whether we have encountered this *signature* before:
	    if signature in base_signatures_table:
		# Yes we have, so reuse the previously assigned name:
		name = base_signatures_table[signature]
		self._name = name
		if tracing_detail >= 3:
		    print("{0}reuse name '{1}'".format(indent, name))
	    else:
		# No we have not so we create a new one *and* write out
		# the associated *name*.scad file:
		name = "{0}{1}".format(base_name, len(base_signatures_table))
		base_signatures_table[signature] = name
		self._name = name
		if tracing_detail >= 3:
		    print("{0}new name '{1}'".format(indent, name))

		# Now we see whether we need to write out the file and
		# run it through openscad:
		stl_directory = ezcad._stl_directory_get()
		stl_file_name = os.path.join(stl_directory,
		  "{0}_{1}.stl".format(name, signature_hash))
		part._stl_file_name = stl_file_name
		if os.path.isfile(stl_file_name):
		    # Since the .stl file already exists, we must have already
		    # written out the .scad file.  Thus, there is nothing
		    # more to do:
		    if tracing_detail >= 3:
			print("{0}{1} already exists".format(indent, stl_file_name))
		    pass
		else:
		    # We need to write out the .scad file and generate the
		    # assocatied .stl file:
		    if tracing_detail >= 3:
			print("{0}{1} does not exist".format(indent, stl_file_name))

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
		    lines.append("    } // union")

		    # Output *scad_difference_lines*:
		    for difference_line in scad_difference_lines:
			lines.append(difference_line)

		    # Close off difference():
		    lines.append("  } // difference")

		    # Perform all the placements:
		    for sub_part in sub_parts:
			#print("Part._manufacture.place={0}".format(place))
			self._scad_transform(lines, 0, center = sub_part._center,
			  axis = sub_part._axis, rotate = sub_part._rotate,
			  translate = sub_part._translate);
			lines.append("{0}();".format(sub_part._name))

		    # Close off the module:
		    lines.append("} // module")
		    lines.append("")

		    # Call the module we just produced:
		    lines.append("{0}();".format(name))
		    lines.append("")

		    # Write out *scad_file*:
		    scad_directory = ezcad._scad_directory_get()
		    scad_file_name = os.path.join(scad_directory, "{0}.scad".format(name))
		    scad_file = open(scad_file_name, "w")
		    scad_file.write("\n".join(lines))
		    scad_file.close()
		    if tracing_detail >= 3:
			print("{0}'{1} written".format(indent, scad_file_name))

		    # Delete any previous *.stl files:
		    stl_directory = ezcad._stl_directory_get()
		    glob_pattern = os.path.join(stl_directory, "{0}_*.stl".format(name))
		    previous_stl_files = glob.glob(glob_pattern)
		    for previous_stl_file in previous_stl_files:
			os.remove(previous_stl_file)
			if tracing_detail >= 3:
			    print("{0}{1} removed".format(indent, previous_stl_file))

		    # Run the command that convert the .scad file into the
		    # associated .stl file:
                    if tracing_detail >= 0:
			print("{0}is_part={1}".format(indent, self._is_part))
		    if self._is_part:
			ignore_file = open("/dev/null", "w")
			scad_directory = ezcad._scad_directory_get()
			scad_file_name = os.path.join(scad_directory, "{0}.scad".format(name))
			command = [ "openscad", "-o", stl_file_name,  scad_file_name ]
			if tracing_detail >= 2:
			    print("{0}command='{1}".format(indent, command))
			subprocess.call(command, stderr=ignore_file) 
			ignore_file.close()

		    # Write out DXF file:
		    dxf_scad_lines = self._dxf_scad_lines
		    if len(dxf_scad_lines) > 0:
			scad_directory = ezcad._scad_directory_get()
			dxf_scad_file_name = os.path.join(scad_directory,
			  "{0}_Dxf.scad".format(name))
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
	    if tracing >= 0:
		print("{0}<==Part._manufacture('{1}'):STL".format(indent, part._name))
    
	if tracing >= 0:
	    print("{0}<=Part._manufacture:('{1}', *)".format(indent, self._name))

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
	""" *Part*: Return the operations list for the *Part* object (i.e. *self*). """
        
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
	      operations[tool1_last_index + 1]._tool_get() == tool:
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
	    if operation._tool_get() == tool:
		# We have another instance of {tool}:
		tool2_first_index = index

		# Find the last operation in the sequence that matches {tool}:
		tool2_last_index = tool2_first_index
		while tool2_last_index + 1 <= last_index and \
		  operations[tool2_last_index + 1]._tool_get() == tool:
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

    def _top_surface_transform_get(self):
	""" *Part*: Return the top surface *Trasform* object from the *Part* object (i.e. *self*).
	"""

	return self._top_surface_transform

    def _tools_dowel_pin_search(self, tracing):
	""" *Part*: Find and return a *Tool_Dowel_Pin* object. """

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._tools_dowel_pin_search('0')".format(indent, self._name))

	# Use *part* instead of *self*:
	part = self

	zero = L()
	dowel_pin = part._tools_search(Tool_Dowel_Pin._match, zero, zero, "dowel_pin", tracing + 1)
	assert isinstance(dowel_pin, Tool_Dowel_Pin)

	if tracing >= 0:
	    tool_name = "NONE"
	    if dowel_pin != None:
		tool_name = dowel_pin._name_get()
	    print("{0}<=Part._tools_dowel_pin_search('0')".format(indent, self._name, tool_name))

	return dowel_pin

    def _tools_drill_search(self, diameter, maximum_z_depth, tracing):
	""" *Part*: Return a drill tool for the *Part* object (i.e. *self*) that has 
	    a diameter of *diameter* and go to a depth of at least *maximum_z_depth*.
	"""

	# Verify argument types:
	assert isinstance(diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(tracing, int)
	
	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._tools_drill_search('{1}', {2:i}, {3:i})".
	      format(indent, self._name, diameter, maximum_z_depth))

	drill_tool = \
	  self._tools_search(Tool_Drill._match, diameter, maximum_z_depth, "drill", tracing + 1)

	if tracing >= 0:
	    tool_name = "NONE"
	    if drill_tool != None:
		toll_name = drill_tool._name_get()
	    print("{0}<=Part._tools_drill_search('{1}', {2:i}, {3:i}) => {4}".
	      format(indent, self._name, diameter, maximum_z_depth, tool_name))

	return drill_tool

    def _tools_end_mill_search(self, maximum_diameter, maximum_z_depth, from_routine, tracing):
	""" *Part*: Search for an end mill with a diameter that is less than or equal to
	    *maximum_diameter* and can mill to a depth of *maximum_z_depth*.  (*from_routine*
	    is used for debugging.
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(maximum_diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(from_routine, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	detail_level = -1
	if tracing >= 0:
	    detail_level = 1
	    indent = ' ' * tracing
	if detail_level >= 0:
	    print("{0}=>Part._tools_end_mill_search('{1}', {2:i}, {3:i}, '{4}')".
	      format(indent, self._name, maximum_diameter, maximum_z_depth, from_routine))

	# Search for a matching *end_mill_tool*:
	search_tracing = -1000000
	if detail_level >= 1:
	    search_tracing = tracing + 1
	end_mill_tool = self._tools_search(Tool_End_Mill._match,
	  maximum_diameter, maximum_z_depth, from_routine, search_tracing)

	# Wrap up any requested *tracing*:
	if detail_level >= 0:
	    tool_name = "NONE"
            if end_mill_tool != None:
		tool_name = end_mill_tool._name_get()
	    print("{0}<=Part._tools_end_mill_search('{1}', {2:i}, {3:i}, '{4}') => '{5}'".format(
	      indent, self._name, maximum_diameter, maximum_z_depth, from_routine, tool_name))

	return end_mill_tool

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

    def _tools_mill_drill_side_search(self, maximum_diameter, maximum_z_depth, tracing):
	""" *Part*: Search a mill drill with a total diameter that is less than or
	    equal to *maximum_diameter* and can mill down to *maximum_z_depth*.
	"""

	# Verify argument types:
	assert isinstance(maximum_diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._tools_mill_drill_side_search('{1}', {2}, {3})".
	      format(indent, self._name, maximum_diameter, maximum_z_depth))

	# Search for a viable *end_mill_tool*:
	end_mill_tool = self._tools_search(Tool_Mill_Drill._mill_drill_side_match,
	  maximum_diameter, maximum_z_depth, "mill drill side", tracing + 1)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    tool_name = "NONE"
	    if end_mill_tool != None:
		tool_name = end_mill_tool._name_get()
	    print("{0}<=Part._tools_mill_drill_side_search('{1}', {2}, {3}) => {4}".
	      format(indent, self._name, maximum_diameter, maximum_z_depth, tool_name))

	return end_mill_tool

    def _tools_mill_drill_tip_search(self, maximum_diameter, maximum_z_depth, tracing):
	""" *Part*: Search a mill drill with a total diameter that is less than or
	    equal to *maximum_diameter* and can mill down to *maximum_z_depth*.
	"""

	# Verify argument types:
	assert isinstance(maximum_diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._tools_mill_drill_tip_search('{1}', {2}, {3})".
	      format(indent, self._name, maximum_diameter, maximum_z_depth))

	# Search for a viable *end_mill_tool*:
	end_mill_tool = self._tools_search(Tool_Mill_Drill._mill_drill_side_match,
	  maximum_diameter, maximum_z_depth, "mill drill side", tracing + 1)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    tool_name = "NONE"
	    if end_mill_tool != None:
		tool_name = end_mill_tool._name_get()
	    print("{0}<=Part._tools_mill_drill_tip_search('{1}', {2}, {3}) => {4}".
	      format(indent, self._name, maximum_diameter, maximum_z_depth, tool_name))

	return end_mill_tool

    def _tools_search(self, match_routine, parameter1, parameter2, from_routine, tracing):
	""" *Part*: Search for the best *Tool* object using *match_routine*,
	    *parameter1*, and *parameter2*.  *from_routine* is used for debugging*:
	"""

	# Verify argument types:
	assert isinstance(parameter1, L)
	assert isinstance(parameter2, L)
	assert isinstance(from_routine, str)
	assert isinstance(tracing, int)

	# Use *part* instead of *self*:
	part = self

	# Perform any requested *tracing*:
	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 1
	    indent = ' ' * tracing
	    print("{0}=>Part._tools_search('{1}', *, {2:i}, {3:i}, '{4}')".
	      format(indent, part._name, parameter1, parameter2, from_routine))

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
	if tracing_detail >= 3:
	    print("{0}available tools={1}".format(indent, len(tools)))

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

	    if tracing_detail >= 3:
		print("{0}Tool:'{1} speed_range_ok={2}'".
		  format(indent, tool._name_get(), speed_range_ok))

	    if speed_range_ok:
		# See whether *tool* has a chance of being a match:
		#if debug:
		#    print("    speed_range={0:.0F}\n".format(speed_range))
		match_tracing = -1000000
		if tracing_detail >= 3:
		    match_tracing = tracing + 1
		priority = match_routine(tool, parameter1, parameter2, from_routine, match_tracing)
		if priority >= 0.0:
		    # Tool is an acceptable match:
		    if tracing_detail >= 2:
			print("{0}Tool: '{1}' priority:{2}".
			  format(indent, tool._name_get(), tool._diameter_get()))

		    tool_preferred = part._tool_preferred
		    if tool_preferred != "":
			if tracing_detail >= 3:
			    print("{0}tool='{1}' preferred='{2}'".
			      format(indent, tool._name, tool_preferred))
			if tool._name == tool_preferred:
			    priority = priority + 100.0
			    if tracing_detail >= 3:
				print("{0}priority={1}".format(indent, priority))

		    # Select the *surface_speed* for *tool_material*:
		    tool._priority = priority
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

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    tool_name = "NONE"
	    if best_tool != None:
		tool_name = best_tool._name_get()
	    print("{0}<=Part._tools_search('{1}', *, {2:i}, {3:i}, '{4}') => '{5}'".
	      format(indent, part._name, parameter1, parameter2, from_routine, tool_name))

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
	
	# Use *part* instead of *self*:
	part = self

	return part._z_safe

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

	# Use *part* instead of *self*:
	part = self

	# Check the *z_rapid* is valid:
	z_rapid = part._z_rapid

	# Return z_rapid regardless:
	return z_rapid

    def _z_rapid_set(self, z_rapid):
	""" *Part*: Set the Z rapid height for horizontal movement above hold-down tooling
	    for *Part* (i.e. *self*) to *z_rapid*
	"""
	
	if not isinstance(z_rapid, L):
            z_rapid = L(inch=0.5)

	# Verify argument types:
	assert isinstance(z_rapid, L)

	# Load *z_rapid* into the *Part* object (i.e. *self*):
	self._z_rapid = z_rapid

    # Public methods come here:

    def block(self, comment, material, color, corner1, corner2, welds, tracing = -100000):
	""" *Part*: Add block of material to the *Part* object (i.e. *self*) that has corners
	    at *corner1* and *corner2*.  The block is made of *material* and rendered in
	    *color*.  Multiple blocks can be "welded" together by specifying that 0 to 6
	    surface layer letters that are weldable (e.g. "tb" would specify that the top
	    and bottom layers are weldable.)
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(material, Material)
	assert isinstance(color, Color)
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)
	assert isinstance(welds, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 1
	    indent = ' ' * tracing
	    print("{0}=>Part.block('{1}', '{2}', '{3}', '{4}', {5:i}, {6:i}, '{7}')".
	      format(indent, part._name, comment, material, color, corner1, corner2, welds))

	# Some constants:
	zero = L()
	one = L(mm=1.0)
	x_axis = P(one, zero, zero)
	y_axis = P(zero, one, zero)
	z_axis = P(zero, zero, one)
	degrees0 = Angle(deg=0.0)
	degrees90 = Angle(deg=90.0)
	degrees180 = Angle(deg=180.0)

	# Record the *color* and *material*:
	part._material = material
	part._color = color

	# Tell the various generators that need to generate a part:
	part._is_part = True

	# Make sure that the corners are diagonal from bottom south west
	# to top north east:
	x1, x2 = corner1.x.minimum_maximum(corner2.x)
	y1, y2 = corner1.y.minimum_maximum(corner2.y)
	z1, z2 = corner1.z.minimum_maximum(corner2.z)

	# Compute the *minimal_corner* and the *maximal_corner*:
        minimal_corner = P(x1, y1, z1)
        maximal_corner = P(x2, y2, z2)

	# Update *bounding_box*:
	bounding_box = part._bounding_box
	bounding_box.point_expand(minimal_corner)
	bounding_box.point_expand(maximal_corner)
	if tracing_detail >= 1:
	    print("{0}corner1={1} corner2={2}".format(indent, minimal_corner, maximal_corner))
	    print("{0}bounding_box={1}".format(indent, bounding_box))

	# Compute the block center points::
	x_center = (x1 + x2) / 2
	y_center = (y1 + y2) / 2
	z_center = (z1 + z2) / 2

	## Do the pre-work needed to compute the *top_surface_transform*.  We will need
	## The *top_surface_point* which is in the center of the top surface.  We need
	## the *top_surface_rotate_axis* to and *top_surface_rotate_angle* to realign
	## the surface with the X/Y plane.
	#
	## Compute *top_surface_point* based on *top* argument:
	#if top.startswith("t"):
	#    # Top surface:
	#    top_surface_point = P(x_center, y_center, z2)
	#    # No rotation needed:
	#    top_surface_rotate_axis = y_axis
	#    top_surface_rotate_angle = degrees0
	#elif top.startswith("b"):
	#    # Bottom surface:
	#    top_surface_point = P(x_center, y_center, z1)
	#    top_surface_rotate_axis = y_axis
	#    top_surface_rotate_angle = degrees180
	#elif top.startswith("n"):
	#    # North surface:
	#    top_surface_point = P(x_center, y2, z_center)
	#    top_surface_rotate_axis = x_axis
	#    top_surface_rotate_angle = degrees90
	#elif top.startswith("s"):
	#    # South surface:
	#    top_surface_point = P(x_center, y1, z_center)
	#    top_surface_rotate_axis = x_axis
	#    top_surface_rotate_angle = -degrees90
	#elif top.startswith("e"):
	#    # East surface:
	#    top_surface_point = P(x2, y_center, z_center)
	#    top_surface_rotate_axis = y_axis
	#    top_surface_rotate_angle = -degrees90
	#elif top.startswith("w"):
	#    # West surface:
	#    top_surface_point = P(x1, y_center, z_center)
	#    top_surface_rotate_axis = y_axis
	#    top_surface_rotate_angle = degrees90
	#else:
        #    assert False, "top argument must start with 't', 'b', 'n', 's', 'e', or 'w'"
	#
	## Compute the *top_surface_transform*:
	#top_surface_transform = Transform()
	#top_surface_transform = top_surface_transform.translate(
	#  "block top to origin", -top_surface_point)
	#top_surface_transform = top_surface_transform.rotate(
	#  "block plane_change to X/Y", top_surface_rotate_axis, top_surface_rotate_angle)
	#
	## Deal with suffix of "90" in *top*:
	#if top.endswith("90"):
	#    assert len(top) == 3, "Too many characters in top argument '{0}'".format(top)
	#    top_surface_transform = top_surface_transform("block rotate 90", z_axis, degrees90)
	#
	## Remember *top_surface_transform* into *part*:
	#part._top_surface_transform = top_surface_transform

	ezcad = part._ezcad
	ezcad_mode = ezcad._mode
	if ezcad_mode == EZCAD3.STL_MODE:
	    # If there are any welds, we make the block a little bigger on each weld surface:
	    adjust = ezcad._adjust * 2
	    weld_extra = adjust.absolute() + L(mm=.01)
	    if 't' in welds:
		z2 += weld_extra
	    if 'b' in welds:
		z1 -= weld_extra
	    if 'n' in welds:
		y2 += weld_extra
	    if 's' in welds:
		y1 -= weld_extra
	    if 'e' in welds:
		x2 += weld_extra	
	    if 'w' in welds:
		x1 -= weld_extra

	    # As with all openscad stuff, the commands go in reverse order onto *union_lines* --
	    # translate first, color second, and cube third:
	    union_lines = part._scad_union_lines
	    pad = ' ' * 6
	    union_lines.append("{0}translate([{1:m}, {2:m}, {3:m}])".
	      format(pad, (x1 + x2)/2, (y1 + y2)/2, (z1 + z2)/2))
	    union_lines.append("{0}color([{1}, {2}, {3}, {4}])".
	      format(pad, color.red, color.green, color.blue, color.alpha))
	    union_lines.append("{0}cube([{1:m}, {2:m}, {3:m}], center=true); // {4}".
	      format(pad, x2 - x1, y2 - y1, z2 - z1, comment))

	if ezcad_mode == EZCAD3.CNC_MODE:
	    # For no do nothing for CNC:
	    pass

	# Wrap up tracing:
	if tracing >= 0:
	    print("{0}<=Part.block('{1}', '{2}', '{3}', '{4}', {5:i}, {6:i}, '{7}')".
	      format(indent, part._name, comment, material, color, corner1, corner2, welds))

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

    def cylinder(self, comment, material, color, start_diameter, end_diameter,
      start, end, sides, sides_angle, welds, flags, tracing = -1000000):
	""" *Part*: Place a *diameter* wide cylinder from *start* to *end*. """

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(material, Material)
	assert isinstance(color, Color)
	assert isinstance(start_diameter, L)
	assert isinstance(end_diameter, L)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(sides, int)
	assert isinstance(sides_angle, Angle)
	assert isinstance(welds, str)
	assert isinstance(flags, str)

        # Use *part* instead of *self*:
	part = self
	part._is_part = True
	part._color = color
	part._material = material

	# Some constants:
        zero = L()

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Part.cylinder('{1}', '{2}', {3}, {4}, {5:i}, {6:i}," +
	      " {7:i}, {8:i}, {9}, {10:d}, '{11}', '{12}')").
	      format(indent, self._name, comment, material, color, start_diameter, end_diameter,
	      start, end, sides, sides_angle, welds, flags))

	# Compute the *maximum_radius*:
	maximum_diameter = start_diameter
        if end_diameter > maximum_diameter:
	    maximum_diameter = end_diameter
	maximum_radius = maximum_diameter / 2

	# Expand the *bounding_box* with 8 points that bracket the maximal cylinder corners:
	thickness = (start - end).length()
	top_transform = Transform.top_surface("top surface", start, end, Angle())
	reverse_transform = top_transform.reverse()
	bounding_box = self._bounding_box
	for x in [-maximum_diameter, maximum_diameter]:
	    for y in [-maximum_diameter, maximum_diameter]:
		for z in [zero, -thickness]:
		    point = reverse_transform * P(x, y, z)
		    bounding_box.point_expand(point)

	if self._ezcad._mode_get() == EZCAD3.STL_MODE:
	    # Output some openscad stuff:
	    pad = ' ' * 4
	    lines = part._scad_union_lines
	    lines.append(("{0}// Part.cylinder('{1}', '{2}', {3}, {4}, {5:i}, {6:i}," +
	      " {7:i}, {8:i}, {9}, {10:d}, '{11}', '{12}')").
	      format(pad, self._name, comment, material, color, start_diameter, end_diameter,
	      start, end, sides, sides_angle, welds, flags))

	    # Perform the cylinder: 
	    part._scad_cylinder(comment, color, start_diameter, end_diameter,
	      start, end, lines, pad, sides, sides_angle, tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print(("{0}<=Part.cylinder('{1}', '{2}', {3}, {4}, {5:i}, {6:i}," +
	      " {7:i}, {8:i}, {9}, {10:d}, '{11}', '{12}')").
	      format(indent, self._name, comment, material, color, start_diameter, end_diameter,
	      start, end, sides, sides_angle, welds, flags))

    def dowel_pin(self, comment, tracing=-1000000):
	""" *Part*: Request that a dowel pin be used to align the *Part* object (i.e. *self*)
	    in the vice with a comment of *comment*.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types
	assert isinstance(comment, str)

	# Perform any requested *tracing*:
	#tracing = part._tracing
	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 1
	    indent = ' ' * tracing
	    print("=>Part.dowel_pin('{0}', '{1}')".format(part._name, comment))

	shop = part._shop
	assert isinstance(shop, Shop)
	if shop._cnc_generate_get():
	    # Find the *tool_dowel_pin* to use :
	    tool_dowel_pin = part._tools_dowel_pin_search(tracing)
	    assert isinstance(tool_dowel_pin, Tool_Dowel_Pin), "No dowel pin tool found"

	    # Immediately grab the *feed_speed* and *spindle_speed* that were computed
	    # for *tool_dowel_pin*:
	    feed_speed = tool_dowel_pin._feed_speed_get()
	    spindle_speed = tool_dowel_pin._spindle_speed_get()

	    if tracing_detail > 0:
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
	    if tracing_detail > 0:
		print("dowel_pin: part={0} bb={1} tip_depth={2}".
		  format(part._name, new_bounding_box, tip_depth))
		print("dowel_pin: before z_depth={0}".format(z_depth))

	    # Figure out how deep the dowel pin actually can go:
	    maximum_z_depth = tool_dowel_pin._maximum_z_depth_get()
	    if z_depth > maximum_z_depth:
		z_depth = maximum_z_depth
	    if tracing_detail > 0:
		print("after z_depth={0}".format(z_depth))

	    if comment == "":
		comment = "Mount {0} x {1} x {2} piece of {3} in vice". \
		  format(part._dx_original, part._dy_original,
		  part._dz_original, _part.material)

	    if tracing_detail > 0:
		print("dowel_pin@({0}): ex={1} ey={2} vx={3} vy={4}".
		  format(part._name, part._dowel_x, part._dowel_y,
		  part._vice_x, part._vice_y))

	    vice = shop._vice_get()
	    jaw_volume = vice._jaw_volume_get()
	    jaw_dx = jaw_volume.x
	    half_jaw_dx = jaw_dx / 2
	    radius = diameter / 2

	    dowel_x = part._dowel_x - radius
	    dowel_y = part._dowel_y
	    plunge_x = dowel_x

	    if tracing_detail > 0:
		print("dowel_pin: {0} jaw_width={1}".
		  format(part._name, jaw_width))

	    # Make sure we always plunge outside of the vise jaws:
	    if plunge_x > -half_jaw_dx:
		plunge_x = -half_jaw_dx

	    # If we are close the jaw edge, move out a little more:
	    if plunge_x > -half_jaw_dx:
		plunge_x = plunge_x - L(inch=0.70)
	    plunge_y = dowel_y
	    plunge_point = P(plunge_x, plunge_y, plunge_z)

	    if tracing_detail > 0:
		print("dowel_pin: {0}: plunge_x={1}".
		  format(part._name, plunge_x))

	    dowel_pin_order = Operation.ORDER_DOWEL_PIN
	    operation = Operation_Dowel_Pin(part,
	      comment, 0, tool_dowel_pin, dowel_pin_order, None, feed_speed, spindle_speed,
	      diameter, dowel_point, plunge_point, top_surface_z, xy_rapid_safe_z, tracing + 1)
	    self._operation_append(operation)

	    if tracing_detail > 0:
		print("dowel_pin('{0}' dowel_point={1:i} plunge_point={2:i}".
		  format(part._name, dowel_point, plunge_point))
	    if tracing_detail > 0:
		print("<=dowel_pin@Part('{0}', '{1}')".
		  format(part._name, comment))

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("<=Part.dowel_pin('{0}', '{1}')".format(part._name, comment))

    def dowel_position_set(self, dowel_point):
	""" *Part*: Set the dowel pin contact point for the *Part* object (i.e. *self*)
	    to *dowel_point*.
	"""

	# This routine will sets the dowel position for {part} to
	# ({dowel_x},{dowel_y}).

	# Verify argument types:
	assert isinstance(dowel_point, P)

	# Grab the *position* matrix and compute the *remapped_dowel_point*:
	position = self._position
	remapped_dowel_point = position.point_multiply(dowel_point)

	# Load up the edge
	self._dowel_x = remapped_dowel_point.x
	self._dowel_y = remapped_dowel_point.y

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

    def extrude(self, comment, material, color,
      contours, start, start_extra, end, end_extra, rotate, tracing = -1000000):
	""" *Part*: Create an extrusion long the axis from *start* to *end* out of *material*
	    (with a render color of *color*.   *contours* is a list of *Contour* objects
	    where the first (required) *Contour* object specifies the outer contour and
	    second and subsequent *Contour* objects (optional) specify internal hole contours.
	    Only the X/Y coordinates of the *Contour* objects are actaully used (i.e. we
	    only use the projection along the Z axis.)  *rotate* specifies how much to
	    rotate the extrusion along the extrusion axis. *comment* will show up in error
	    messages and anycreated G-code.
	"""

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(material, Material)
	assert isinstance(color, Color)
	assert isinstance(contours, list) and len(contours) >= 1
	assert isinstance(start, P)
	assert isinstance(start_extra, L)
	assert isinstance(end, P)
	assert isinstance(end_extra, L)
	assert isinstance(rotate, Angle)
	assert isinstance(tracing, int)

	# Use *part* instead of *self*:
	part = self

	# Perform any requested *tracing*
	#if part._ezcad._mode == EZCAD3.STL_MODE:
	#    tracing = 0
	#if tracing == -1000000:
	#    tracing = part._tracing
	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 1
	    indent = ' ' * tracing
	    print("{0}=>Part.extrude('{1}', {2}, {3}, *, {4:i}, {5:i}, {6:i}, {7:i}, {8:d})".
	      format(indent, comment, material, color, start, start_extra, end, end_extra, rotate))

	# Some constants:
	zero = L()
	one = L(mm=1.0)
	degrees0 = Angle()
	z_axis = P(zero, zero, one)

	# Record the *color* and *material*:
	part._material = material
	part._color = color

	
	extrude_axis = end - start
	extrude_direction = extrude_axis.normalize()
	start_extra_delta = -start_extra.millimeters() * extrude_direction
	start += start_extra_delta
	end_extra_delta = end_extra.millimeters() * extrude_direction
	end += end_extra_delta
	extrude_axis = end - start
	height = extrude_axis.length()

	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 1
	    print("{0}direction={1:m} start_extra_delta={2:i} end_extra_delta={3:i}".
	      format(indent, extrude_direction, start_extra_delta, end_extra_delta))

	# Compute the transform that will place *start* on the X/Y plane with *start*-*end*
	# pointing in the negative Z axis direction:
	top_surface_transform = \
	  Transform.top_surface("extrude", start, end, rotate, tracing + 1)
	if tracing_detail >= 2:
	    print("{0}top_surface_transform={1:s}".format(indent, top_surface_transform))

	# Remember *top_surface_transform* for other operations:
	part._top_surface_transform = top_surface_transform

	# Now compute the transforms that will take a contour in the X/Y plane and map
	# them to appropate planes at *start* and *end*:
	start_contour_transform = top_surface_transform.reverse()
	end_contour_transform = start_contour_transform.translate("move to end", end - start)

	if tracing_detail >= 3:
	    print("{0}before start_contour_transform={1:s}".format(indent, start_contour_transform))

	# Extract the *outer_contour*:
	outer_contour = contours[0]
	outer_contour._project(top_surface_transform, tracing + 1)

	if tracing_detail >= 3:
	    print("{0}middle1 start_contour_transform={1:s}".
	      format(indent, start_contour_transform))

	# Now we can expand the *bounding_box* by visiting each *bend* in *bends*:
	bounding_box = self._bounding_box
	bends = outer_contour._bends_get()
	for bend in bends:
	    # Compute the mapped bend points:
	    projected_point = bend._projected_point_get()
	    bend_start_point = start_contour_transform * projected_point
	    bend_end_point = end_contour_transform * projected_point
	    if tracing_detail >= 2:
		print("{0}point={1:i} bend_start_point={2:i} bend_end_point={3:i}".
		  format(indent, point, bend_start_point, bend_end_point))
	    if tracing_detail >= 3:
		print("{0}start_bend after plane_change={1:i}".
		  format(indent, plane_change.point_multiply(bend_start_point)))

	    # Now expand *bounding_box*:
	    bounding_box.point_expand(bend_start_point)
	    bounding_box.point_expand(bend_end_point)
	if tracing_detail >= 1:
	    print("{0}bounding_box={1:i}".format(indent, bounding_box))
	if tracing_detail >= 3:
	    print("{0}after start_contour_transform={1:s}".format(indent, start_contour_transform))

	part._dowel_x = zero
	part._dowel_y = zero

	self._is_part = True
	ezcad = self._ezcad
	assert isinstance(ezcad, EZCAD3)
	if ezcad._mode == EZCAD3.STL_MODE:
	    if tracing_detail >= 1:
		print("{0}==>Part.extrude:STL".format(indent))
	    # Set up the transform and extrude:
            lines = self._scad_union_lines
	    pad = ' ' * 6

	    lines.append("{0}// extrude('{1}', {2}, {3}, {4:m}, {5:m}, {6:m}, {7:m}, {8:d})".
	      format(pad, comment, material, color, start, start_extra, end, end_extra, rotate))

	    if tracing >= 0:
		print("{0}start_contour_transform={1:s}".format(indent, start_contour_transform))
	    start_contour_transform._scad_lines_append(lines, pad)

	    lines.append(
	      "{0}translate([ 0, 0, {1:m} ]) // Put top of extrusion at origin".
	      format(pad, -height))
            lines.append(
	      "{0}linear_extrude(height = {1:m}, convexity=10)".format(pad, height))
	    
	    # Output the polygon operation:
	    lines.append("{0}polygon(".format(pad))
	    lines.append("{0}  convexity = 10,".format(pad))
	    lines.append("{0}  points = [".format(pad))
	    contours_size = len(contours)
	    for contour_index, contour in enumerate(contours):
		line = "{0}    ".format(pad)
		bends = contour._bends_get()
		bends_size = len(bends)
		for bend_index, bend in enumerate(bends):
		    bend_point = bend._point_get()
		    bend_comma = ","
		    if bend_index + 1 >= bends_size and contour_index + 1 >= contours_size:
			bend_comma = ""
		    line += "[{0:m}, {1:m}]{2} ".format(bend_point.x, bend_point.y, bend_comma)
		line += "// Contour {0}".format(contour_index)
		lines.append(line)
            lines.append("{0}  ],".format(pad))
            lines.append("{0}  paths = [".format(pad))

	    # Output the contour points, staring with the outer contour and followed
	    # by any inner contours:
	    point_index = 0
	    for contour_index, contour in enumerate(contours):
		line = "{0}    [ ".format(pad)
		# Output one set of contour points:
		bends = contour._bends_get()
		bends_size = len(bends)
		for bend_index, bend in enumerate(bends):
 		    bend_point = bend._point_get()
		    bend_comma = ","
		    if bend_index + 1 >= bends_size:
			bend_comma = ""
		    line += "{0}{1} ".format(point_index, bend_comma)
		    point_index += 1

		# Output End of line stuff:
		contour_comma = ","
		if contour_index + 1 >= contours_size:
		    contour_comma = ""
		line += "]{0} // Contour {1}".format(contour_comma, contour_index)
		lines.append(line)

	    # Wrap the polygon up:
            lines.append("{0}  ]".format(pad))
            lines.append("{0}); // polygon".format(pad))

	    # Wrap up any *tracing*:
	    if tracing_detail >= 1:
		print("{0}<==Part.extrude:STL".format(indent))

	# Perform any requested tracing;
	if tracing >= 0:
	    print("{0}<=Part.extrude('{1}', {2}, {3}, *, {4:i}, {5:i}, {6:i}, {7:i}, {8:d})".
	      format(indent, comment, material, color, start, start_extra, end, end_extra, rotate))

    def xhole(self, comment = "NO_COMMENT", diameter = None,
      start = None, end = None, sides = -1, sides_angle = Angle(),
      top = "t", flags = "", start_diameter = None, end_diameter = None,
      trace = -10000):
	""" *Part*: Make a *diameter* hole in the *Part* object (i.e. *self*) with starting 
	    at *start_point* and ending at *end_point*.   *comment will show in any error
	    messages and any generated G-code.  The allowed flag letters in *flags* are:

	      One of 't' (default), 'f', or 'p':
		't'	through hole (i.e. Through)
		'f'	flat hole (i.e. Flat)
		'p'	tip hole (drill tip stops at {end_point} (i.e. tiP)

	      Allowed additional flags:
		'u'	upper hole edge should be chamfered (i.e. Upper)
		'l'	lower hole edge should be chamfered (i.e. Lower)
		'm'	hole is to be milled (i.e. Milled)
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
	#debug = True
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

	# Now visit *part* and all of its children in STL mode:
	ezcad._mode = EZCAD3.STL_MODE
	part._manufacture(ezcad, -1000000)
	ezcad._update_count += 1

	# Now visit *part* and all of its children in visualization mode:
	ezcad._mode = EZCAD3.VISUALIZATION_MODE
	part._manufacture(ezcad, -1000000)
	ezcad._update_count += 1

	# Now visit *part* and all of its children in CNC mode:
	ezcad._mode = EZCAD3.CNC_MODE
	shop = ezcad._shop
	shop._cnc_generate_set(True)
	#part._manufacture(ezcad, 0)
	part._manufacture(ezcad, -1000000)
	shop._cnc_generate_set(False)
	ezcad._update_count += 1

	if debug:
	    print("<=Part.process('{0}')".format(part._name))

    def rectangular_tube_extrude(self, comment, material, color,
      width, height, thickness, start, start_extra, end, end_extra, rotate, tracing = -1000000):
	""" *Part*: Create an extrusion long the axis from *start* to *end* out of *material*
	    (with a render color of *color*.)  The tube will be *height* high, *width* wide,
	    with tube wall thickness of *thickness*.  If *thickness* is zero, the entire
	    tube will be filled in.  The tube starts at *start* and ends at *end*.  *rotate*
	    specifies how much to rotate the extrusion along the extrusion axis. *comment*
	    will show up in error messages and anycreated G-code.
	"""

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(material, Material)
	assert isinstance(color, Color)
	assert isinstance(width, L)
	assert isinstance(height, L)
	assert isinstance(thickness, L) and thickness >= L()
	assert isinstance(start, P)
	assert isinstance(start_extra, L)
	assert isinstance(end, P)
	assert isinstance(end_extra, L)
	assert isinstance(rotate, Angle)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	#tracing = 0
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Part.rectangular_tube_extrude('{1}', '{2}', {3}, {4}, {5:i}, {6:i}," +
	      " {7:i}, {8:i}, {9:i}, {10:i}, {11:i}, {12:d})").
	      format(indent, self._name, comment, material, color,
	      width, height, thickness, start, start_extra, end, end_extra, rotate))

	# Specify the X/Y coordinates of the rectangular tube:
	x1 = -height/2
	x2 = x1 + thickness
	x3 = -x2
	x4 = -x1	
	y1 = -width/2
	y2 = y1 + thickness
	y3 = -y2
	y4 = -y1

	# Construct the *outer_contour* and stuff into *contours* list::
	zero = L()
	outer_contour = Contour("outer tube contour")
	outer_contour.bend_append("outer SW", P(x1, y1, zero), zero)
	outer_contour.bend_append("outer NW", P(x1, y4, zero), zero)
	outer_contour.bend_append("outer NE", P(x4, y4, zero), zero)
	outer_contour.bend_append("outer SE", P(x4, y1, zero), zero)
	contours = [outer_contour]

	# If *thickness* is non-zero, we provide an *inner_contour* to append to *contours*:
	if thickness > zero:
	    # Construct the *inner_contour*:
	    inner_contour = Contour("inner tube contour")
	    inner_contour.bend_append("inner SW", P(x2, y2, zero), zero)
	    inner_contour.bend_append("inner NW", P(x2, y3, zero), zero)
	    inner_contour.bend_append("inner NE", P(x3, y3, zero), zero)
	    inner_contour.bend_append("inner SE", P(x3, y2, zero), zero)
	    contours.append(inner_contour)

	# Perform the extrusion:
	self.extrude(comment, material, color,
	  contours, start, start_extra, end, end_extra, rotate, tracing + 1)

	# Wrap up any requested *tracing*;
	if tracing >= 0:
	    print(("{0}<=Part.rectangular_tube_extrude('{1}', '{2}', {3}, {4}, {5:i}, {6:i}," +
	      " {7:i}, {8:i}, {9:i}, {10:i}, {11:i}, {12:d})\n").
	      format(indent, self._name, comment, material, color,
	      width, height, thickness, start, start_extra, end, end_extra, rotate))

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

    def round_tube_extrude(self, comment, material, color,
      outer_diameter, inner_diameter, start, start_extra, end, end_extra, tracing = -1000000):
	""" *Part*: Create an extrusion long the axis from *start* to *end* out of *material*
	    (with a render color of *color*.)  The tube will be extended in the along the
	    line through *start* and *end* by *start_exta* on the side of *start* and by
	    *end_extra* on the side of *end*.   The tube will have an outside diameter of
	    *outer_diameter* and an inside diameter of *inner_diameter*.  If *inner_diameter*
	    is zero, the tube will be entirely filled with material (i.e. a rod.)  *comment*
	    will show up in error messages and any created G-code.
	"""

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(material, Material)
	assert isinstance(color, Color)
	assert isinstance(outer_diameter, L)
	assert isinstance(inner_diameter, L)
	assert isinstance(start, P)
	assert isinstance(start_extra, L)
	assert isinstance(end, P)
	assert isinstance(end_extra, L)
	assert isinstance(tracing, int)

	# Some constants:
	zero = L()
	degrees0 = Angle(deg=0.0)
	degrees360 = Angle(deg=360.0)

	# Perform any requested *tracing*:
	#tracing = 0
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Part.round_tube_extrude('{1}, '{2}', {3}, {4}, {5:i}, {6:i}," +
	      " {7:i}, {8:i}, {9:i}, {10:i})").
	      format(indent, self._name,  comment, material, color,
	      outer_diameter, inner_diameter, start, start_extra, end, end_extra))

	# Specify the X/Y coordinates of the rectangular tube:
	outer_radius = outer_diameter/2
	outer_contour = Contour("{0} outer tube contour".format(comment))
	count = 32
	delta_angle = degrees360 / float(count)
	for index in range(count):
	    angle = delta_angle * float(index)
	    x = outer_radius.cosine(angle)
	    y = outer_radius.sine(angle)
	    contour_comment = "{0} outer contour {1}".format(comment, index)
	    outer_contour.bend_append(contour_comment, P(x, y, zero), outer_radius)
	contours = [outer_contour]

	# If *inner_diameter* is non-zero, we provide an *inner_contour* to append to *contours*:
	if inner_diameter > zero:
	    # Construct the *inner_contour*:
	    if tracing >= 0:
		print("{0}inner_diameter={1:i}".format(indent, inner_diameter))
	    inner_radius = inner_diameter/2
	    inner_contour = Contour("{0} inner tube contour".format(comment))
	
	    for index in range(count):
		angle = delta_angle * float(index)
		x = inner_radius.cosine(angle)
		y = inner_radius.sine(angle)
		contour_comment = "{0} outer contour {1}".format(comment, index)
		inner_contour.bend_append(contour_comment, P(x, y, zero), inner_radius)
	    contours.append(inner_contour)

	# Perform the extrusion:
	self.extrude(comment, material, color,
	  contours, start, start_extra, end, end_extra, degrees0, tracing + 1)

	# Wrap up any requested *tracing*;
	if tracing >= 0:
	    print(("{0}<=Part.round_tube_extrude('{1}, '{2}', {3}, {4}, {5:i}, {6:i}," +
	      " {7:i}, {8:i}, {9:i}, {10:i})").
	      format(indent, self._name,  comment, material, color,
	      outer_diameter, inner_diameter, start, start_extra, end, end_extra))

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

    def virtual_box(self, comment, corner1, corner2):
	""" *Part*: """

	# Check argument types:
	assert isinstance(comment, str)
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)

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

    def _wrl_write(self,
      wrl_file, transform, wrl_indent, file_name, parts_table = {}, tracing = -1000000):
	""" *Part*: Write *self* to *wrl_file*. """

	# Check argument types:
	assert isinstance(wrl_file, file)
	assert isinstance(transform, Transform)
	assert isinstance(wrl_indent, int)
	assert isinstance(parts_table, dict)
	assert isinstance(file_name, str)
	assert isinstance(tracing, int)

	# Use *part* instead of *self*:
	part = self

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.wrl_write('{1}', file_name='{2}')".
	      format(indent, part._name, file_name))

	# Make sure the top level starts with an empty *parts_table*:
	if wrl_indent == 0:
	    parts_table = {}
	    wrl_file.write("#VRML V2.0 utf8\n")
	    wrl_file.write("# Generated by EZCAD3\n")
	    wrl_file.write("Viewpoint { description \"Initial view\" position 0 0 150 }\n")
	    wrl_file.write("NavigationInfo { type \"EXAMINE\" }\n")

	# Do some preparation work:
	name = self._name
	spaces = ' ' * wrl_indent


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
		    stl_directory = ezcad._stl_directory_get()
		    stl_file_name = os.path.join(
		      stl_directory, "{0}_{1}.stl".format(name, self._signature_hash))
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
			#  Extract *point1*, *point2*, and *point3* from *stl_lines*:
			xlist1 = stl_lines[index + 2].split()
			point1 = P(L(float(xlist1[1])), L(float(xlist1[2])), L(float(xlist1[3])))
			xlist2 = stl_lines[index + 3].split()
			point2 = P(L(float(xlist2[1])), L(float(xlist2[2])), L(float(xlist2[3])))
			xlist3 = stl_lines[index + 4].split()
			point3 = P(L(float(xlist3[1])), L(float(xlist3[2])), L(float(xlist3[3])))

			# Transform *point1*, *point2*, and *point3*:
			point1 = transform * point1
			point2 = transform * point2
			point3 = transform * point3

			# Now convert the transformed *point1*, *point2*, and *point3*, back
			# back into immutable Python tuples that  can be used to index into
                        # a Python dictionary:
			vertex1 = point1.triple()
			vertex2 = point2.triple()
			vertex3 = point3.triple()
    
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
		    spaces = ' ' * wrl_indent
    
		    # Write out "DEF name Shape {":
		    #wrl_file.write(
		    #  "{0}DEF x{1} Shape {{\n".format(spaces, name))
		    wrl_file.write(
		      "{0}#VRML V2.0 utf8\n".format(spaces))
		    wrl_file.write(
		      "{0}Shape {{\n".format(spaces, name))
    
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
		    if color.alpha < 1.0:
			wrl_file.write(
			  "{0}   transparency {1}\n".format(spaces, 1.0 - color.alpha))
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
			    sub_part._wrl_write(wrl_file,
			      transform, wrl_indent + 2, parts_table, file_name)
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
			    sub_part._wrl_write(wrl_file, transform,
			      wrl_indent + 4, parts_table, file_name)
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

	if tracing >= 0:
	    print("{0}<=Part.wrl_write('{1}', file_name='{2}')".
	      format(indent, part._name, file_name))
    
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

    def contour(self, comment, contour, start_point, end_point, extra, flags, tracing = -1000000):
	""" *Part*: Perform an exterior contour operation of the *Part* object (i.e. *self*)
	    using *contour* to describe the path.  The Z start and stop depths are
	    obtained from *start_point.z* and *end_point.z*.  *extra* is the amount
	    of extra material being removed.  (What does *extra* really do?)  *flags*
	    can contain the letter 'u' for an upper chamfer, 'l' for a lower chamfer,
	    and 't' for a contour that cuts entirely through the part entirely.
	    *comment* is used in error messages and any generated G-code.
	"""

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(contour, Contour)
	assert isinstance(start_point, P)
	assert isinstance(end_point, P)
	assert isinstance(extra, L)
	assert isinstance(flags, str)
	for flag in flags:
	    assert flag == 'u' or flag == 'l' or flag == 't', \
	      "Part.contour: Bad flag '{0}' in flags '{1}' for part '{2}'". \
	      format(flag, flags, part._name)

	# Use *part* instead of *self*:
	part = self

	# Perform any requested *tracing*:
	#if tracing == -1000000:
	#    tracing = part._tracing
	tracing_detail = -1
	if tracing >= 0:
	    #tracing_detail = 0
	    indent = ' ' * tracing
	    print("{0}=>Part.contour('{1}, '{2}', {3:i}, {4:i}, {5:i}, '{6}')".
	     format(indent, part._name, comment, start_point, end_point, extra, flags))

	# Before we do anything else, we need to update the bounding box for *part*:
	bounding_box = part._bounding_box
	contour._bounding_box_expand(bounding_box)

	# Figure out if we are in *cnc_mode* or *stl_mode*:
	ezcad = part._ezcad
	cnc_mode = (ezcad._mode == EZCAD3.CNC_MODE)
	stl_mode = (ezcad._mode == EZCAD3.STL_MODE)

	if stl_mode:
	    # Compute the *top_surface_transform*:
	    degrees0 = Angle(deg=0.0)
	    top_surface_transform = \
	      Transform.top_surface("extrude", start_point, end_point, degrees0, tracing + 1)
	    part._top_surface_transform = top_surface_transform

	    # Perform all of the *contour* massaging:
	    top_surface_transform = part._top_surface_transform
	    contour._project(top_surface_transform, tracing + 1)
	    contour._inside_bends_identify(tracing + 1)
	    contour._radius_center_and_tangents_compute(tracing + 1)

	    # Output the *top_surface_transform_reverse* to *difference_lines*:
	    difference_lines = part._scad_difference_lines
	    pad = ' ' * 4
	    top_surface_transform_reverse = top_surface_transform.reverse()
	    top_surface_transform_reverse._scad_lines_append(difference_lines, pad)

	    # Output the linear_extrude to *difference_lines*:
	    extra = L(mm=1.0)
	    difference_lines.append("{0}// Contour {1}".format(pad, contour._name_get()))
	    height = (end_point - start_point).length()
	    difference_lines.append("{0}translate([0, 0, {1:m}])".format(pad, -height - extra))
	    difference_lines.append("{0}linear_extrude(height = {1:m}) {{".
	     format(pad, height + 2.0 * extra))
	    contour._scad_lines_polygon_append(difference_lines, pad, True, tracing + 1)
	    difference_lines.append("{0}}} // linear_extrude".format(pad))

	# Do any requested CNC generation:
	if cnc_mode:
	    # Grab the previously computed *top_surface_transform* that reorients the
            # material properly for the CNC machine:
	    top_surface_transform = part._top_surface_transform

	    # Perform all of the *contour* massaging:
	    contour._project(top_surface_transform, tracing + 1)
	    contour._inside_bends_identify(tracing + 1)
	    contour._radius_center_and_tangents_compute(tracing + 1)

	    # Project all of the points in *contour* onto a plane normal to *projection_axis*
	    # and then rotate that plane to be the X/Y plane:
	    contour._project(top_surface_transform, tracing + 1)

	    # The returned value from *_smallest_inner_diameter_compute*() will be either
	    # negative (for no smallest inner corner radius), or positive for the smallest
	    # inner corner radius.  Either value will find the correct end-mill or mill-drill
	    # tool for searching purposes:
	    zero = L()
	    contour._inside_bends_identify(tracing + 1)
	    smallest_inner_radius = contour._smallest_inner_radius_compute(tracing + 1)
	    smallest_inner_diameter = smallest_inner_radius * 2

	    transformed_start_point = top_surface_transform * start_point
	    transformed_end_point = top_surface_transform * end_point

	    # Find a *Tool* to use for edge milling.  For non-through contours, we should use
	    # a mill drill in preference to an end-mill because that way we can overlap
	    # with any top chamfering hole countersinking (i.e. one fewer tool change).
	    # Otherwise use an end mill:
	    z_stop = transformed_start_point.z - transformed_end_point.z
	    if tracing >= 0:
		print("{0}z_stop={1:i}".format(indent, z_stop))
	    mill_drill_tool = \
	      part._tools_mill_drill_side_search(smallest_inner_diameter, z_stop, tracing + 1)
	    have_mill_drill = isinstance(mill_drill_tool, Tool_Mill_Drill)
	    end_mill_tool = \
	      part._tools_end_mill_search(smallest_inner_diameter, z_stop, "contour", tracing + 1)
	    have_end_mill = isinstance(end_mill_tool, Tool_End_Mill)
	    if tracing >= 0:
		print("{0}mill_drill_tool='{1}'".format(indent, mill_drill_tool))
		print("{0}end_mill_tool='{1}'".format(indent, end_mill_tool))
		print("{0}have_mill_drill={1} have_end_mill={2}".
		  format(indent, have_mill_drill, have_end_mill))

	    # Expand the *flags* into *do_upper_chamfer*, *do_lower_chamfer* and *do_through*:
	    do_upper_chamfer = False
	    do_lower_chamfer = False
	    do_through = False
	    for flag in flags:
		if flag == 't':
		    do_through = True
		elif flag == 'l':
		    do_through = True
		elif flag == 'u':
		    do_through = True

	    # Now select which *mill_tool* to use from either *mill_drill_tool* or *end_mill_tool*:
	    mill_tool = None
	    if do_through:
		mill_drill_tool = None
		have_mill_drill = False
	    if have_end_mill and have_mill_drill:
		# We have both *end_mill_tool* and *mill_drill_tool*.  Select based on priority:
		if mill_drill_tool._priority_get() >= end_mill_tool._priority_get():
		    # The *mill_drill_tool* won on priority:
		    mill_tool = mill_drill_tool
		    have_end_mill = False
		else:
		    # The *end_mill_tool* won on priority:
		    mill_tool = end_mill_tool
		    have_mill_drill = False
	    elif have_mill_drill:
		# The only *mill_tool* available is *mill_drill_tool*:
		mill_tool = mill_drill_tool
	    elif have_end_mill:
		# The only *mill_tool* available is *end_mill_tool*:
		mill_tool = end_mill_tool

	    assert isinstance(mill_tool, Tool), \
	      "{0}:{1} no contour tool: smallest diameter={2:i} z_stop={3:i}".format(
	       part._name, comment, smallest_inner_diameter, z_stop)

	    # Figure out *point_angle* and *tip_depth*:
	    point_angle = Angle(deg=180)	# End-mills are flat
	    tip_depth = zero			# End-mills have no tip
	    if isinstance(mill_tool, Tool_Mill_Drill):
		# Alternatively, mill drills have both a *point_angle* and a *tip_depth*:
		point_angle = mill_drill._point_angle_get()
		tip_depth = mill_tool.tip_depth(point_angle)

	    # Schedule the operation if we found a tool, otherwise error:
	    mill_operation = None
	    if True:
		# Compute *depth_maximum* as a function of *diameter* and *extra*:
		#	{cutter_engagement}	{depth_maximum}
		#	=======================================
		#	1.00 (100%)		{diameter}/4
		#	0.50 (50%)		{diameter}/2
		#	0.25 (25%)		{diameter}/1
		diameter = mill_tool._diameter_get()
		assert diameter > zero

		depth_maximum = zero

		is_laser = mill_tool._is_laser_get()
		if is_laser:
		    # We always cut to maximum depth with the laser:
		    depth_maximum = mill_tool._maximum_z_depth_get()
		else:
		    # When *diameter* equals *extra* the ratio of *diameter*
		    # to *extra* is 1.0.   When *extra* is 1/4th of *diameter*
		    # the ratio of *diameter* over *extra* is 4.0:
		    ratio = diameter / extra

		    # Cap {ratio} between 1.0 and 4.0:
		    if ratio < 1.0:
			ratio = 1.0
		    elif ratio > 4.0:
			ratio = 4.0

		    # *depth_maximum* is *diameter* * *ratio* / 4:
		    depth_maximum = diameter * (ratio / 4.0)
		    if tracing >= 0:
			print("{0}diameter={1:i} extra={2:i} ratio={3} depth_max={4:i}".
			  format(indent, diameter, extra, ratio, depth_maximum))

		# Figure out how many *passes* using *depth_maximum*, *z_start*, *z_stop*,
		# and *z_extra:
		z_start = transformed_start_point.z
		z_end = transformed_end_point.z
		z_extra = zero
		z_depth = (z_start - z_end).absolute() + z_extra

		if tracing >= 0:
		    print("{0}z_depth={1} depth_maximum={2}".format(indent, z_depth, depth_maximum))
		passes = int(z_depth / depth_maximum) + 1

		if tracing >= 0:
		    print("{0}zstt={1:i} zend={2:i} depmax={3:i} zxtr={4:i} zdpth={5:i} pss={6}".
		      format(indent, z_start, z_end, depth_maximum, z_extra, z_depth, passes))
	
		# Set {sub_priority} to 0 to force this operation before
		# top/bottom chamfers:
		sub_priority = 0
		follows = None
		offset = zero
		corners = None
		bends = contour._bends_get()
		mill_feed_speed = mill_tool._feed_speed_get()
		mill_spindle_speed = mill_tool._spindle_speed_get()
		mill_operation_contour = Operation_Contour(part, comment,
		  sub_priority, mill_tool, Operation.ORDER_MILL_DRILL_EXTERIOR,
		  follows, mill_feed_speed, mill_spindle_speed,
		  z_start - tip_depth, z_start - z_depth - tip_depth,
		  contour, offset, diameter/2, passes, tracing + 1)
		part._operation_append(mill_operation_contour)

		if tracing >= 1:
		    print("{0}zstt-td={1:i} zstt-zdep-td={2:i}".format(
		      indent, z_start - tip_depth, z_start - z_depth - tip_depth))

	    # Deal with *do_upper_chamfer*:
	    if do_upper_chamfer:
		# We need to do an upper chamfer that is *top_chamfer* wide:
		# Find a mill drill tool that can do {smallest_inner_radius}:
		top_chamfer_tool = part._tools_mill_drill_chamfer_search(
		  smallest_inner_diameter, z_stop)
		assert not op_chamfer_tool is None, \
		  "No top chamfer tool for '{0}' (diameter={1:i}".format(
		  comment, smallest_inner_diameter)

		if True:
		    # Tool found, schedule the operation:
		    assert isinstance(top_chamfer_tool, Tool_Mill_Drill)
		    tip_depth = mill_drill._tip_depth_get()
		    top_chamfer_comment = "'{0}' top chamfer".format(comment)
		    nominal_tool_radius = top_chamfer_tool._diameter_get() / 4.0
		    z_stop_top_chamfer = z_start - tip_depth/2 - top_chamfer

		    #print("top_ch={0:i} tip_depth={1:i} ntr={2:i} zs={3:i} zstc={4:i}".format(
		    #  top_chamfer, tip_depth, nominal_tool_radius, z_start, z_stop_top_chamfer))

		    # Set {sub_priority} to 1, to be after main contour:
		    sub_priority = 1
		    top_chamfer_operation_contour = Operation_Contour(top_chamfer_comment,
		      sub_priority, top_chamfer_tool,
		      Operation.MILL_DRILL_CHAMFER,
		      mill_operation, z_start, z_stop_top_chamfer, corners,
		      -top_chamfer, nominal_tool_radius, 1)
		    part.operation_append(top_chamfer_operation_contour)

	    # Deal with bottom chamfers:
	    if do_lower_chamfer:
		# We need to do a lower chamfer that is {bottom_chamfer} wide:
		# Find a dove tail tool that can do the {smallest_inner_radius}:
		bottom_chamfer_tool = part.tools_dove_tail_search(smallest_inner_diameter, z_stop)

		if bottom_chamfer_tool == None:
		    # Tool not found, let the user know:
		    assert False, "No bottom chamfer tool for '{0}' (diameter={1:i})".format(
		      comment, smallest_inner_diameter)
		else:
		    # Tool found, schedule the operation:
		    bottom_chamfer_comment = "%f% bottom chamfer".format(comment)
		    nominal_tool_radius = bottom_chamfer_tool.diameter / 4.0
		    # Set {sub_priority} to 1, to be after main contour:
		    sub_priority = 1
		    bottom_chamfer_operation_contour = Operation_Contour(part,
		       bottom_chamfer_comment,
		      sub_priority, bottom_chamfer_tool,
		      Operation.DOVE_TAIL_CHAMFER,
		      z_start, z_stop - nominal_tool_radius + bottom_chamfer,
		      corners, -bottom_chamfer, nominal_tool_radius, 1)

	if tracing >= 0:
	    print("{0}<=Part.contour('{1}, '{2}', '{3:i}', {4:i}', {5:i}, '{6}')".
	     format(' ' * tracing, part._name, comment, start_point, end_point, extra, flags))

    def countersink_hole(self, comment, hole_diameter, countersink_diameter,
      start, stop, flags, sides = -1, sides_angle=Angle(), tracing = -1000000):
	""" *Part*: Put a *hole_diameter* hole into the *Part* object (i.e. *self*) starting
	    at *start* to an end depth of *end*.  If *countersink_diameter* is non-zero,
	    the part will be have a 90 degree countersink of *countersink_diameter* at *start*.
	    *comment* will show up in any  generated RS-274 code.  *hole_kind* specifies the
	    kind of hole.
	"""
    
	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(hole_diameter, L)
	assert isinstance(countersink_diameter, L)
	assert isinstance(start, P)
	assert isinstance(stop, P)
	assert isinstance(flags, str)
	assert isinstance(sides, int)
	assert isinstance(sides_angle, Angle)
	assert isinstance(tracing, int)

	# Use *part* instead of *self*:
	part = self

	# Some constants:
	zero = L()
	degrees0 = Angle(deg=0.0)

	# Perform an requested *tracing*:
	#if self._name[:-1] == "Gear_Box_Shelf" and comment == "Bearing Hole":
	#    tracing = 0
	#if tracing == -1000000:
	#    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.countersink_hole('{1}' '{2}' {3:i} {4:i} {5:i} {6:i} '{7}'".
	      format(indent, part._name, comment, hole_diameter, countersink_diameter,
	      start, stop, flags))
    
	# Compute the radii:
	hole_radius = hole_diameter/2
	countersink_radius = countersink_diameter/2

	# Compute *is_tip_hole*, *is_flat_hole* and *is_through_hole* from *flags*:
	is_through_hole = False
	is_tip_hole = flags.find('p') >= 0
	is_flat_hole = flags.find('f') >= 0
	if not (is_tip_hole or is_flat_hole):
	    is_through_hole = True
	    if part._laser_preferred:
		is_flat_hole = True
		is_through_hole = False
	is_countersink = countersink_diameter > hole_diameter
	if tracing >= 0:
            print("{0}is_through_hole={1} is_tip_hole={2} is_flat_hole={3} is_countersink={4}".
	      format(indent, is_through_hole, is_tip_hole, is_flat_hole, is_countersink))

	# Compute *hole_kind*:
	if is_tip_hole:
	    hole_kind = Part.HOLE_TIP
	elif is_flat_hole:
	    hole_kind = Part.HOLE_TIP
	elif is_through_hole:
	    hole_kind = Part.HOLE_THROUGH
	else:
	    assert False, "No hole kind specified"
	if tracing >= 0:
            print("{0}hole_kind={1}".format(indent, hole_kind))
        
	ezcad = part._ezcad
	mode = ezcad._mode_get()
	if tracing >= 0:
            print("{0}mode == {1}".format(indent, mode))
    	if mode == EZCAD3.CNC_MODE:
	    spot_operation = None
	    try_flat = False
	    if is_through_hole or is_tip_hole:
		# Spot drill and countersink the hole at the same time:
		z_countersink = start.z - countersink_radius
		z_depth = start.z - stop.z

		tool_mill_drill = None
		if is_countersink:
		    tool_mill_drill = \
		      part._tools_mill_drill_tip_search(L(inch=-1.0), z_countersink, tracing + 1)
		tool_drill = part._tools_drill_search(hole_diameter, z_depth, tracing + 1)
    
		if tool_drill != None:
		    # Drill the spot and countersink first:
		    countersink_comment = "{0} [countersink]".format(comment)
    
		    # We want to ensure that countersinks occur first:
		    operation_countersink = None
		    if is_countersink:
			operation_countersink = Operation_Drill(part, counter_sink_comment,
			  sub_priority, tool_mill_drill, Operation.ORDER_MILL_DRILL_COUNTERSINK,
			  None, hole_diameter, Part.HOLE_TIP, start, stop, True)
			part._operation_append(operation_countersink)
    
		    # Now drill the hole:
		    sub_priority = 1
		    operation_drill = Operation_Drill(part, comment, sub_priority, tool_drill,
		      Operation.ORDER_DRILL, operation_countersink, hole_diameter, hole_kind,
		      start, stop, False)
		    part._operation_append(operation_drill)
		else:
		    try_flat = True
	    elif hole_kind == Part.HOLE_FLAT:
		try_flat = True
	    if tracing >= 0:
                print("{0}try_flat={1}".format(indent, try_flat))
    
	    # See if we should try to mill the hole:
	    if try_flat:
		end_mill_tool = part._tools_end_mill_search(hole_diameter,
		  z_depth, "countersink_hole", tracing + 1)
		if end_mill_tool != None:
		    operation_round_pocket = Operation_Round_Pocket(part, comment,
		      0, end_mill_tool, Operation.ORDER_END_MILL_ROUND_POCKET, None,
		      hole_diameter, countersink_diameter, hole_kind, start, end,
		      end_mill_tool._feed_speed_get(), end_mill_tool._spindle_speed_get(),
		      part._top_surface_transform, tracing + 1)
		    part._operation_append(operation_round_pocket)

	if mode == EZCAD3.STL_MODE:
	    drill_start_extra = 0.500
	    if is_flat_hole or is_tip_hole:
		drill_end_extra = 0.0
	    elif is_through_hole:
		drill_end_extra = hole_radius.millimeters() + 0.500
	    else:
		assert False, "Unkown hole type"
	    if tracing >= 0:
		print("{0}drill_start_extra={1} drill_end_extra={2}".
		  format(indent, drill_start_extra, drill_end_extra))

	    lines = part._scad_difference_lines
	    # For debugging:, set *lines* to *_scad_union_lines*:
	    #lines = part._scad_union_lines

	    pad = ' ' * 4
	    lines.append("{0}// countersink_hole('{1}' '{2}' {3:i} {4:i} {5:i} {6:i} '{7}'".
	      format(pad, part._name, comment, hole_diameter, countersink_diameter,
	      start, stop, flags))

	    # Perform the drill with a cylinder:
	    drill_direction = (stop - start).normalize()
	    drill_start = start - drill_direction * drill_start_extra
	    drill_end = stop + drill_direction * drill_end_extra
	    
	    part._scad_cylinder(comment + " drill", Color("black"), hole_diameter, hole_diameter,
	      drill_start, drill_end, lines, pad, sides, sides_angle, tracing + 1)

	    # Peform any requested countersink with a cone:
	    if is_countersink:
		countersink_start_extra = countersink_radius
		countersink_end_extra = countersink_radius
		countersink_start_diameter = countersink_diameter * 2
		countersink_start = start - drill_direction * countersink_start_extra.millimeters()
		countersink_end = start + drill_direction * countersink_end_extra.millimeters()
		zero = L()
		part._scad_cylinder(comment + " countersink", Color("black"),
		  countersink_start_diameter, zero, countersink_start, countersink_end, lines, pad,
		  -1, degrees0 , tracing + 1)

	if tracing >= 0:
	    print("{0}<=Part.countersink_hole('{1}' '{2}' {3:i} {4:i} {5:i} {6:i} '{7}'".
	      format(indent, part._name, comment, hole_diameter, countersink_diameter,
	      start, stop, flags))

    def _scad_cylinder(self, comment, color,
      start_diameter, end_diameter, start, end, lines, pad, sides, sides_angle, tracing = -1000000):
	""" *Part*: Generate an open some openscad code the *Part* object (i.e. *self*)
	    for drawing a cylinder (or cone) that is *start_diameter* on at *start*
	    *end_diameter* at *end*.  *lines* is the list to which the lines are
	    appended with each line prefixed by *pad*.  The cylinder is rendered with *color*.
	"""    

	# Verify argument types:
	assert isinstance(color, Color)
	assert isinstance(comment, str)
	assert isinstance(start_diameter, L)
	assert isinstance(end_diameter, L)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(lines, list)
	assert isinstance(pad, str)
	assert isinstance(sides, int)
	assert isinstance(sides_angle, Angle)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._cylinder('{1}', '{2}', {3:i}, {4:i}, {5:i}, {6:i}, *, {7}, {8})".
	      format(indent, self._name, comment,
	        start_diameter, end_diameter, start, end, pad, color))

	# Only do the append in STL mode:
	if self._ezcad._mode == EZCAD3.STL_MODE:

	    # Append to *lines*:
	    lines.append("{0}// {1} cylinder!".format(pad, comment))

	    # Perform the transform:
	    top_surface_transform = Transform.top_surface(comment, start, end, sides_angle)
	    top_surface_transform_reverse = top_surface_transform.reverse()
	    top_surface_transform_reverse._scad_lines_append(lines, pad)

	    # Compute center axis of the *cylinder*, its *length* and its *center*:
	    zero = L()
	    axis = end - start
	    height = axis.length()
	    assert height > zero, "Cylinder '{0}' has no height".format(comment)

	    lines.append("{0}translate([0, 0, {1:m}])".format(pad, -height))

            r1 = end_diameter / 2
	    r2 = start_diameter / 2
	    assert r1 + r2 > zero, \
	      "sd={0:i} ed={1:i} r1={2:i} r2={3:i}".format(start_diameter, end_diameter, r1, r2)

	    if sides < 3:
		sides = 12
	    if r1 == r2:
		command = "{0}cylinder(r={1:m}, h={2:m}, $fn={3}, $fa=12);".format(
		  pad, r1, height, sides)
	    else:
		command = "{0}cylinder(r1={1:m}, r2={2:m}, h={3:m}, $fn={4}, $fa=12);".format(
		  pad, r1, r2, height, sides)
	    if tracing >= 0:
		print("{0}Part._cylinder: r1={1} r2={2} command='{3}'".
		  format(indent, r1, r2, command))

	    # Output the cylinder in a vertical orientation centered on the
	    # origin.  It is processed before either rotation or translation:
            lines.append(command)
	

	# Wrap up any *tracing*:
	if tracing >= 0:
	    print("{0}<=Part._cylinder('{1}', '{2}', {3:i}, {4:i}, {5:i}, {6:i}, *, {7}, {8})".
	      format(indent, self._name, comment,
	        start_diameter, end_diameter, start, end, pad, color))

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
	""" *Part*: Set the block dimensions increase for the next *block*() operation
	    to increase the block dimentions by *east*, *west*, *north*, *south*, *top*,
	    and *bottom*.  All of these argument must be greater than or equal to zero.
	    These values are immediately reset after the *block*() operation.
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(east, L)   and east   >= zero
	assert isinstance(west, L)   and west   >= zero
	assert isinstance(north, L)  and north  >= zero
	assert isinstance(south, L)  and south  >= zero
	assert isinstance(top, L)    and top    >= zero
	assert isinstance(bottom, L) and bottom >= zero

	# Remember the extra values:
	self._extra_east =   east
	self._extra_west =   west
	self._extra_north =  north
	self._extra_south =  south
	self._extra_top =    top
	self._extra_bottom = bottom

    def extra_xyz(self, dx, dy, dz):
	""" *Part*: Specify the amount of extra materal to be added to the next *block*()
	    operation to be *dx* in the east/west direction, *dy* in the north/south direction,
	    add *dz* in the top/bottom direction.  These value are reset after the next
	    *block*() operation.
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(dx, L) and dx >= zero
	assert isinstance(dy, L) and dy >= zero
	assert isinstance(dz, L) and dz >= zero

	# Pass everything on to *extra_ewnstb*():
	half_dx = dx/2
	half_dy = dy/2
	half_dz = dz/2
	self.extra_ewnstb(half_dx, half_dx, half_dy, half_dy, half_dz, half_dz)

    def extrusion(self, comment, kind, material, color, start, end,
      a_width, a_thickness, b_width, b_thickness, rotate):
	""" Part dimensions: Create a {kind} extrusion manufactured out of
	    {material} that goes from {start_point} to {end_point} rotated
	    by {rotate}.  {a_width}, {a_thickness}, {b_width} and
	    {b_thickness} specify the extra dimentions of the extrusion.
	    The returned part is visualized using {color}. The extrusion
	    must be aligned with one of the X, Y or Z axes. """

	# Check argument types:
	assert isinstance(comment, str)
	assert isinstance(color, Color)
	assert isinstance(material, Material)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(a_width, L)
	assert isinstance(a_thickness, L)
	assert isinstance(b_width, L)
	assert isinstance(b_thickness, L)
	assert isinstance(rotate, Angle)

	# Define soem useful abreviations:
	zero = L()

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

    def hole(self, comment, diameter, start, stop, flags,
      sides=-1, sides_angle=Angle(), tracing = -1000000):
	""" *Part*:  Put a *diameter* hole long the axis from *start* down to *stop*
	    into the *Part* object (i.e. *self*).  *comment* will show up in any
	    generated RS-274 code or error messages:
	"""

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(diameter, L)
	assert isinstance(start, P)
	assert isinstance(stop, P)
	assert isinstance(flags, str)
	assert isinstance(sides, int)
	assert isinstance(sides_angle, Angle)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	#if tracing == -1000000:
	#    tracing = self._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>hole('{1}', '{2}', {3:i}, {4:i}, {5:i}, '{6}')".
	      format(indent, self._name, comment, diameter, start, stop, flags))

	# Perform the hole using the richer *countesink_hole* operation:
	zero = L()
	self.countersink_hole(comment, diameter, zero, start, stop, flags,
          sides=sides, sides_angle=sides_angle, tracing=tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}=>hole('{1}', '{2}', {3:i}, {4:i}, {5:i}, '{6}')".
	      format(indent, self._name, comment, diameter, start, stop, flags))

    def lathe(self, comment, material, color, start, end, contour, faces, tracing = -1000000):
	""" *Part*: Add a a circular symetry object to the *Part* object (i.e. *self*).
	    that starts at *start* and extends towards *end*.  The part is made out of
	    *material* and is rendered in *color*.  *faces* specifies how many faces around
	    the circle to produce.  A *faces* value of 0 will cause a reasonable number
	    of faces to be produced.  The *contour* object specifies the edge profile
	    for the lathe along the Y axis, where all X coordinates specify the radius
	    for the profile.  All X coordinates must be positive, since negative radii
	    are not allowed.  *start* is attached to top-most point in Y direction.
	    Please note, that the object is *NOT* streached to match *end*; *end* only
	    specifies the direction of the lathe contour from *start*.
	"""

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(material, Material)
	assert isinstance(color, Color)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(contour, Contour)
	assert isinstance(faces, int)
	assert isinstance(tracing, int)

	# Use *part* instead of *self*:
	part = self

	# Perform any requested *tracing*:
	#tracing = 0
	tracing_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    tracing_detail = 3
	    print("{0}=>Part.lathe('{1}', '{2}', {3}, {4:i}, {5:i}, *, {6})".
	      format(indent, part._name, comment, material, color, start, end, contour, faces))

	part._is_part = True
	part._color = color
	part._material = material

	ezcad = part._ezcad
	mode = ezcad._mode_get()
	if mode == EZCAD3.STL_MODE:
	    part._is_part = True
	    pad = ' ' * 4

	    identity = Transform()
	    contour._project(identity, tracing + 1)
	    contour._inside_bends_identify(tracing + 1)
	    contour._radius_center_and_tangents_compute(tracing + 1)

	    # Define some constants:
	    zero = L()
	    degrees0 = Angle(deg=0.0)

	    # Make sure *faces* is reasonable:
	    if faces <= 0:
		faces = 16
	    faces = max(faces, 3)

	    # Find *minimum_y* and *maximim_y* for *contour*:
	    bends = contour._bends_get()
	    assert len(bends) > 0, \
	      "Lathe operaton for '{0}:{1}' has no bends in contour".format(part._name, comment)
	    for index, bend in enumerate(bends):
		point = bend._projected_point_get()
		x = point.x
		y = point.y
		assert x >= zero
		if index == 0:
		    y_minimum = y_maximum = y
		else:
		    y_minimum = min(y_minimum, y)
		    y_maximum = max(y_maximum, y)
	    if tracing >= 0:
		print("{0}y_minimum={1:i} y_maximum={2:i}".format(indent, y_minimum, y_maximum))

	    # Do the reverse transform:
	    union_lines = part._scad_union_lines
	    top_surface_transform = Transform.top_surface("lathe", start, end, degrees0)
	    part._top_surface_transform = top_surface_transform
	    top_surface_transform_reverse = top_surface_transform.reverse()
	    top_surface_transform_reverse._scad_lines_append(union_lines, pad)

	    # This is a little stange.  The openscad rotate_extude command takes a profile
	    # in the Y direction and renders it in the Z direction.  Hence, the translate
	    # below uses -*y_maximum* to translate down the Z axis:
	    union_lines.append("{0}translate([0, 0, {1:m}])".format(pad, -y_maximum))

	    # Now we can do a rotate_extrude of the appropriate contour polygon:
	    union_lines.append("{0}color({1:s})".format(pad, color))
	    union_lines.append("{0}rotate_extrude($fn={1})".format(pad, faces))
	    contour._scad_lines_polygon_append(union_lines, pad, False, tracing + 1)

	# Wrap-up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.lathe('{1}', '{2}', {3}, {4:i}, {5:i}, *, {6})".
	      format(indent, part._name, comment, material, color, start, end, contour, faces))

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
	""" *Part*: """
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

    def _scad_transform(self, lines, indent,
       center = None, axis = None, rotate = None, translate = None, tracing = -1000000):
	""" *Part*: """

	none_type = type(None)
	assert isinstance(lines, list)
	assert isinstance(indent, int)
	assert isinstance(center, P) or center == None
	assert isinstance(axis, P) or axis == None
	assert isinstance(rotate, Angle) or rotate == None
	assert isinstance(translate, P) or translate == None
	assert isinstance(tracing, int)

	ezcad = self._ezcad
	if ezcad._mode == EZCAD3.STL_MODE:
            assert isinstance(lines, list)
	    spaces = " " *indent

	    # What we want to do is 4 transforms:
	    #
            #  1) Move the object such that its *center* is at the origin:
	    #  2) Rotate the object around *axis* by *angle*:
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

	    # Steps 1-3 only make sense if exists and is *angle* is non-zero:
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

	print("Part.screw_holes() does not actually do anything" )

	#print "=>Part.screw_holes({0})".format(self.name)

	## Some useful abbreviations:
	#big = L(inch=123456789.0)
	#
	#ezcad = self._ezcad
	#xml_stream = ezcad._xml_stream
	#if xml_stream != None:
	#    # Grab the six surfaces from {self}:
	#    b = self.b
	#    e = self.e
	#    n = self.n
	#    s = self.s
	#    t = self.t
	#    w = self.w
	#
	#    # Extract the bounding box of {self}:
	#    t_z = t.z
	#    b_z = b.z
	#    n_y = n.y
	#    s_y = s.y
	#    e_x = e.x
	#    w_x = w.x
	#
	#    # Figure out {screw_levels_table} based on {top_surface}:
	#    screw_levels = self._screw_levels
	#    top_surface = self._top_surface
	#    screw_levels_table = {}
	#    if self._top_surface_set:
	#	if top_surface == t:
	#	    screw_levels_table = screw_levels["T"]
	#	    top_surface_name = "Top"
	#	    t_z = big
	#	    b_z = -big
	#	elif top_surface == b:
	#	    screw_levels_table = screw_levels["B"]
	#	    top_surface_name = "Bottom"
	#	    t_z = big
	#	    b_z = -big
	#	elif top_surface == n:
	#	    screw_levels_table = screw_levels["N"]
	#	    top_surface_name = "North"
	#	    n_y = big
	#	    s_y = -big
	#	elif top_surface == s:
	#	    screw_levels_table = screw_levels["S"]
	#	    top_surface_name = "South"
	#	    n_y = big
	#	    s_y = -big
	#	elif top_surface == e:
	#	    screw_levels_table = screw_levels["E"]
	#	    top_surface_name = "East"
	#	    e_x = big
	#	    w_x = -big
	#	elif top_surface == w:
	#	    screw_levels_table = screw_levels["W"]
	#	    top_surface_name = "West"
	#	    e_x = big
	#	    w_x = -big
	#	else:
	#	    assert False, \
	#	      "Unexpected top surface for part {0} is {1}". \
	#	      format(self.name, top_surface)
	#
	#    if len(screw_levels_table) == 0:
	#	assert errors_suppress, \
	#	  "Part '{0}' does not have any attached holes on {1} surface".\
	#	  format(self.name, top_surface_name)
	#    else:
	#	for screw_level in screw_levels_table.values():
	#	    #print "screw_level=", screw_level
	#	    screw = screw_level.screw
	#	    #print "screw.name=", screw.name
	#	    trace = False
	#	    #trace = screw.name.find("skin_west_bottom") >= 0
	#
	#	    # Grap {anchor_point_mapped} from {screw}:
	#	    anchor_point_mapped = screw.anchor_point_mapped
	#	    if trace:
	#		print "=>Part.screw_holes({0})".format(self.name)
	#		print "anchor_point_mapped={0}". \
	#		  format(anchor_point_mapped)
	#
	#	    remap_matrix = screw_level.reverse_matrix
	#	    if trace:
	#		print "screw_level_forward_matrix=\n{0}". \
	#		  format(screw_level.forward_matrix.mat	)
	#		print "screw_level_reverse_matrix=\n{0}". \
	#		  format(screw_level.reverse_matrix.mat)
	#		#print "remap_matrix=\n{0}".format(remap_matrix.mat)
	#	    anchor_point_remapped = \
	#	      remap_matrix.point_multiply(anchor_point_mapped, self)
	#	    x = anchor_point_remapped.x
	#	    y = anchor_point_remapped.y
	#	    z = anchor_point_remapped.z
	#	    if trace:
	#		bse = self.point("$BSE")
	#		tnw = self.point("$TNW")
	#		print "bse={0}".format(bse)
	#		print "tnw={0}".format(tnw)
	#		print "anchor_point_remapped={0}". \
	#		  format(anchor_point_remapped)
	#		#print \
	#		#  "w_x={0} e_x={1} s_y={2} n_y={3} b_z={4} t_z={5}". \
	#		#  format(w_x, e_x, s_y, n_y, b_z, t_z)
	#
	#	    # Make sure everything is in on the part:
	#	    assert w_x <= x and x <= e_x, \
	#	      ("X (={0}) for screw {1} not between {2}" + \
	#	      " and {3} (part={4})"). \
	#	      format(x, screw.name, w_x, e_x, self.name)
	#	    assert s_y <= y and y <= n_y, \
	#	      ("Y (={0}) for screw {1} not between {2}" + \
	#	      " and {3} (part={4})"). \
	#	      format(y, screw.name, s_y, n_y, self.name)
	#	    assert b_z <= z and z <= t_z, \
	#	      ("Z (={0}) for screw {1} not between {2}" + \
	#	      " and {3} (part={4})"). \
	#	      format(z, screw.name, b_z, t_z, self.name)
	#
	#	    # Compute {start_point} and {end_point} based on
	#	    # {top_surface} and {depth}:
	#	    depth = screw_level.depth
	#	    if top_surface == t:
	#		start_point = self.point_new(x, y, t.z)
	#		end_point = start_point.z_adjust(-depth)
	#	    elif top_surface == b:
	#		start_point = self.point_new(x, y, b.z)
	#		end_point = start_point.z_adjust(depth)
	#	    elif top_surface == n:
	#		start_point = self.point_new(x, n.y, z)
	#		end_point = start_point.y_adjust(-depth)
	#	    elif top_surface == s:
	#		start_point = self.point_new(x, s.y, z)
	#		end_point = start_point.y_adjust(depth)
	#	    elif top_surface == e:
	#		start_point = self.point_new(e.x, y, z)
	#		end_point = start_point.x_adjust(-depth)
	#	    elif top_surface == w:
	#		start_point = self.point_new(w.x, y, z)
	#		end_point = start_point.x_adjust(depth)
	#	    else:
	#		assert False
	#
	#	    # No do either a through hole or a hole to {depth}:
	#	    if not screw_level.done:
	#		thread = screw.thread
	#		flags = screw_level.flags
	#		#print "Part.screw_holes: part={0} screw={1}" + \
	#		#  " thread={2} flags={3}". \
	#		#  format(self.name, screw_level.screw.name, \
	#		#  thread, flags)
	#		if depth <= L():
	#		    # Drill all the way through:
	#		    self.screw_through(screw.name,thread, \
	#		      start_point, flags)
	#		else:
	#		    # Drill to the specified {end_point}:
	#		    self.screw_hole(screw.name, \
	#		      thread, start_point, end_point, flags)
	#
	#		# Remember that we did this {screw_level}:
	#		screw_level.done = True
	#
	#	if trace:
	#	    print "<=Part.screw_holes({0})".format(self.name)

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

    def simple_pocket(self,
      comment, corner1, corner2, radius, pocket_top, rotate, flags, tracing = -1000000):
	""" *Part*: Create a simple rectangular pocket in the *Part* object (i.e. *self*)
	    bounding corners of *bottom_corner* and *top_corner*, a corner radius if *radius*.
	"""

	# Check argument types:
	assert isinstance(comment, str)
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)
	assert isinstance(radius, L)
	assert isinstance(pocket_top, str)
	assert isinstance(rotate, Angle)
	assert isinstance(flags, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.simple_pocket('{1}', '{2}', {3:i}, {4:i}, {5:i}, '{6}', {7:d}, '{8}')".
	      format(indent, self._name,
	      comment, corner1, corner2, radius, pocket_top, rotate, flags))

	# Use *part* instead of *self:
        part = self

	# Make sure that the corners are diagonal from bottom south west to top north east:
	x1 = min(corner1.x, corner2.x)
	x2 = max(corner1.x, corner2.x)
	y1 = min(corner1.y, corner2.y)
	y2 = max(corner1.y, corner2.y)
	z1 = min(corner1.z, corner2.z)
	z2 = max(corner1.z, corner2.z)

	# Verify that we have properly ordered coordinates:
	update_count = EZCAD3.update_count_get()
	if x1 >= x2 and update_count == 0:
            print("Part.simple_pocket:x1={0} should be less than x2={1}".format(x1, x2))
	if y1 >= y2 and update_count == 0:
	    print("Part.simple_pocket:y1={0} should be less than y2={1}".format(y1, y2))
	if z1 >= z2 and update_count == 0:
            print("Part.simple_pocket:z1={0} should be less than z2={1}".format(z1, z2))

	# Create the *bsw_corner* and *tne_corner* points:
	bsw_corner = P(x1, y1, z1)
	tne_corner = P(x2, y2, z2)

	# Some constants:
	zero = L()
	
	ezcad = self._ezcad
	start_extra = L(mm = 1.0)
	end_extra = zero
	if flags == "t":
	    end_extra = L(mm = 1.0)
	#adjust = ezcad._adjust

	if ezcad._mode == EZCAD3.STL_MODE:
	    if tracing >= 0:
                print("{0}STL_MODE".format(indent))
                #print("bsw_corner={0:m} tne_corner={1:m}".format(bsw_corner, tne_corner))

	    # Compute the center points:
	    x_center = (x1 + x2) / 2
	    y_center = (y1 + y2) / 2
	    z_center = (z1 + z2) / 2

	    # Deal with *pocket_top* argument:
	    if pocket_top == "t":
		start = P(x_center, y_center, z2)
		end =   P(x_center, y_center, z1)
	    elif pocket_top == "b":
		start = P(x_center, y_center, z1)
		end =   P(x_center, y_center, z2)
	    elif pocket_top == "n":
		start = P(x_center, y2, z_center)
		end =   P(x_center, y1, z_center)
	    elif pocket_top == "s":
		start = P(x_center, y1, z_center)
		end =   P(x_center, y2, z_center)
	    elif pocket_top == "e":
		start = P(x2, y_center, z_center)
		end =   P(x1, y_center, z_center)
	    elif pocket_top == "w":
		start = P(x1, y_center, z_center)
		end =   P(x2, y_center, z_center)
	    else:
		assert False, \
		  "pocket_top = '{0}' instead of 't', 'b', 'n', 's', 'e', or 'w'". \
		  format(pocket_top)

	    # Compute *top_transform* and transform the two corners:
	    top_transform = Transform.top_surface(comment, start, end, rotate, tracing + 1)
	    transformed_bsw_corner = top_transform * bsw_corner
	    transformed_tne_corner = top_transform * tne_corner

	    # Now extract the X/Y/Z locatons of the corners:
	    tx1 = transformed_bsw_corner.x
	    ty1 = transformed_bsw_corner.y
	    tz1 = transformed_bsw_corner.z
	    tx2 = transformed_tne_corner.x
	    ty2 = transformed_tne_corner.y
	    tz2 = transformed_tne_corner.z

	    # Now compute the minimum and maximum coordinates:
	    x1 = min(tx1, tx2)
	    y1 = min(ty1, ty2)
	    z1 = min(tz1, tz2) - end_extra
	    x2 = max(tx1, tx2)
	    y2 = max(ty1, ty2)
	    z2 = max(tz1, tz2) + start_extra

	    # Generate the openscad stuff:
	    difference_lines = self._scad_difference_lines
	    pad = ' ' * 6

            # Output the transform that will put everything in the correct location:
	    top_transform_reverse = top_transform.reverse()
	    top_transform_reverse._scad_lines_append(difference_lines, pad)

	    difference_lines.append(
	      "{0}// Part.simple_pocket('{1}', '{2}', {3:i}, {4:i}, {5:i}, '{6}', {7:d}, '{8}')".
	      format(pad, self._name, comment, corner1, corner2, radius, pocket_top, rotate, flags))
	    #difference_lines.append(
	    #  "{0}//start={1:m} end={2:m}".format(pad, start, end))
	    #difference_lines.append(
	    #  "{0}//transformed_bsw_corner={1:m}".format(pad, transformed_bsw_corner))
	    #difference_lines.append(
	    #  "{0}//transformed_tne_corner={1:m}".format(pad, transformed_tne_corner))

	    # Translate the forthcoming polygon:
	    difference_lines.append("{0}translate([0, 0, {1:m}])".
	      format(pad, -(z2 - z1) + start_extra))

	    # Linear extrude the forthcoming polygon:
	    difference_lines.append("{0}linear_extrude(height = {1:m}, center = false)".
              format(pad, z2 - z1))

	    # Output a polygon that represents the pocket.

	    # First put together a list *bend_points* which are the four pocket corners
            # moved inwards by *radius*:
	    bend_points = [
 	      P(x1 + radius, y1 + radius, zero),
 	      P(x1 + radius, y2 - radius, zero),
 	      P(x2 - radius, y2 - radius, zero),
 	      P(x2 - radius, y1 + radius, zero)
            ]

	    # This is where the compute the various angles need to have *corner_sides* sides
            # in each pocket corner:
	    corner_sides = 4
	    degrees90 = Angle(deg=90)
	    angle_delta = degrees90/corner_sides

	    # Assemble *polygon_points* which are the points which are the outline of the pocket:
	    polygon_points = []
	    # The first corner has a *start_angle* that points downward:
	    start_angle = -degrees90

	    # Now iterate through each *bend_point* output the the *corner_sides* + 1 points
            # that make up the pocket bend:
	    for bend_point in bend_points:
		angle = start_angle
		for index in range(corner_sides + 1):
		    adjust = P.polar(angle, radius)
		    corner_point = bend_point + adjust
		    #print("bend_point={0} adjust={1}".format(bend_point, adjust))
		    polygon_points.append(corner_point)
		    angle -= angle_delta
                start_angle -= degrees90

	    # Now create a list of fragments that will be joined together to form the
            # *polygon_command*:
	    polygon_command_parts = []
	    polygon_command_parts.append("{0}polygon(points = [".format(pad))
	    prefix = ""
	    for polygon_point in polygon_points:
		polygon_command_parts.append(
		  "{0}[{1:m}, {2:m}]".format(prefix, polygon_point.x, polygon_point.y))
                prefix = ", "
	    polygon_command_parts.append("], paths = [{0}], convexity = 8);".
              format(range(len(polygon_points))))

	    # Convert *polygon_command_parts* into *polygon_command* and append
            #it to *difference_lines*:
	    polygon_command = "".join(polygon_command_parts)
	    difference_lines.append(polygon_command)

	    # Lower the cube so that it cuts out the pocket.
	    #difference_lines.append(
	    #  "{0}translate([0, 0, {1:m}])".format(pad, -(z2 - z1)/2))

	    # Finally, we can output the cube:
	    #difference_lines.append(
	    #  "{0}cube([{1:m}, {2:m}, {3:m}], center=true);".
	    #  format(pad, x2 - x1, y2 - y1, z2 - z1))

	if ezcad._mode == EZCAD3.CNC_MODE:
	    if tracing >= 0:
                print("{0}Part.simple_pocket: CNC_MODE started".format(indent))

	    # Grab the *top_surface_transform* from *part* that has been previously
            # set to orient the material properly for the CNC machine:
	    top_surface_transform = part._top_surface_transform

	    # Transform the corners
	    transformed_corner1 = top_surface_transform * corner1
	    transformed_corner2 = top_surface_transform * corner2

	    # Extract the X/Y/Z coordinates from *transformed_corner1* and *transformed_corner2*:
	    transformed_x1 = transformed_corner1.x
	    transformed_y1 = transformed_corner1.y
	    transformed_z1 = transformed_corner1.z
	    transformed_x2 = transformed_corner2.x
	    transformed_y2 = transformed_corner2.y
	    transformed_z2 = transformed_corner2.z

	    # Figure out the new corner coordinates:
	    new_x1 = transformed_x1.minimum(transformed_x2)
	    new_x2 = transformed_x1.maximum(transformed_x2)
	    new_y1 = transformed_y1.minimum(transformed_y2)
	    new_y2 = transformed_y1.maximum(transformed_y2)
	    new_z1 = transformed_z1.minimum(transformed_z2)
	    new_z2 = transformed_z1.maximum(transformed_z2)

	    # Now build the 2 new corners:
	    new_corner1 = P(new_x1, new_y1, new_z1)
	    new_corner2 = P(new_x2, new_y2, new_z2)
	    dz = new_z2 - new_z1

	    # Search for a matching *end_mill_tool*:
	    maximum_diameter = 2 * radius
	    end_mill_tool = \
	      self._tools_end_mill_search(maximum_diameter, dz, "simple_pocket", tracing + 1)
	    assert end_mill_tool != None, \
	      "Could not find a end mill to mill {0:i} radius pockets".format(radius)
	    if tracing >= 0:
	        print("{0}end_mill_tool='{1}'".format(indent, end_mill_tool._name_get()))
	    end_mill_diameter = end_mill_tool._diameter_get()
	    end_mill_radius = end_mill_diameter / 2
	    end_mill_feed_speed = end_mill_tool._feed_speed_get()
	    end_mill_spindle_speed = end_mill_tool._spindle_speed_get()

	    # Create the pocket operation:
	    operation_order = Operation.ORDER_END_MILL_SIMPLE_POCKET
	    operation_simple_pocket = Operation_Simple_Pocket(self, comment, 0,
	      end_mill_tool, operation_order, None, end_mill_feed_speed, end_mill_spindle_speed,
	      new_corner1, new_corner2, radius, end_mill_radius, Operation.POCKET_KIND_FLAT)
	    self._operation_append(operation_simple_pocket)

	    if tracing >= 0:
                print("{0}Part.simple_pocket: CNC_MODE done".format(indent))
	# Perform any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.simple_pocket('{1}', '{2}', {3:i}, {4:i}, {5:i}, '{6}', {7:d}, '{8}')".
	      format(indent, self._name,
	      comment, corner1, corner2, radius, pocket_top, rotate, flags))

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

	# Verify argument types:
	assert isinstance(tool_name, str)

	self._tool_preferred = tool_name

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

	    #ezcad = self._ezcad
	    #xml_stream = ezcad._xml_stream
	    #if xml_stream != None:
	    #	xml_stream.write('{0}<Tooling_Plate Rows="{1}" '. \
	    #	  format(" " * ezcad._xml_indent, rows))
	    #	xml_stream.write('Columns="{0}" Comment="{1}">\n'. \
	    #	  format(columns, comment))
	    #    
	    #	ezcad.xml_indent_push()
	    #	for row in range(0, rows):
	    #	    for column in range(0, columns):
	    #		row_column = (row, column)
	    #		adjust = (0, 0)
	    #		if row_column in adjusts:
	    #		    adjust = adjusts[row_column]
	    #		remove = ""
	    #		if row_column in removes:
	    #		    remove = removes[row_column]
	    #
	    #		if trace:
	    #		    print "[{0},{1}] = {2}{3}". \
	    #		      format(row, column, adjust, remove)
	    #
	    #		xml_stream.write('{0}<Tooling_Hole'. \
	    #		  format(" " * ezcad._xml_indent))
	    #		xml_stream.write(' Row="{0}" Column="{1}"'. \
	    #		  format(row, column))
	    #		xml_stream.write(' Adjust_X="{0}" Adjust_Y="{1}"'. \
	    #		  format(adjust[0], adjust[1]))
	    #		xml_stream.write(' Flags="{0}"/>\n'. format(remove))
	    #	    
	    #	ezcad.xml_indent_pop()
	    #	xml_stream.write('{0}</Tooling_Plate>\n'. \
	    #	  format(" " *ezcad._xml_indent))

    def tooling_plate_drill(self, comment, columns, rows, skips, tracing=-1000000):
        """ *Part*: Force the drilling of 
	"""

	# Use *part* instead of *self:
	part = self

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(columns, tuple) or isinstance(columns, list)
	assert isinstance(rows, tuple) or isinstance(rows, list)
	assert isinstance(skips, tuple) or isinstance(skips, list)

	# Do some more argument verification:
	for row in rows:
	    assert isinstance(row, int)
	for column in columns:
            assert isinstance(column, int)
	for skip in skips:
            assert isinstance(skip, tuple) or isinstance(skip, list)
	    assert(len(skip) == 2) and isinstance(skip[0], int) and isinstance(skip[1], int)

	# Only do tracing for the CNC branch:
	tracing_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.tooling_plate_drill('{1}', '{2}', {3}, {4}, {5})".format(
	      indent, part._name, comment, rows, columns, skips))
	    tracing_detail = 1


	shop = part._shop_get()
        #if True or shop._cnc_generate_get():

	# Grab the *tooling_plate* from *shop* and some associated values::
	tooling_plate = shop._tooling_plate_get()
	tooling_plate_dx = tooling_plate._dx_get()
	tooling_plate_dy = tooling_plate._dy_get()
	tooling_plate_dz = tooling_plate._dz_get()
	tooling_plate_rows = tooling_plate._rows_get()
	tooling_plate_columns = tooling_plate._columns_get()
	tooling_plate_drill_diameter =  tooling_plate._drill_diameter_get(part._material)
	tooling_plate_hole_pitch = tooling_plate._hole_pitch_get()

	# Grab the *vice* from *shop* and some associated values:
	vice = shop._vice_get()
	vice_jaw_volume = vice._jaw_volume_get()
	vice_jaw_volume_dz = vice_jaw_volume.z
	parallels = vice._parallels_get()
	    
	# Figure out which *parallel_height* to use from *parallels* to use to mount
        # the *tooling_plate*:
        parallels_heights = sorted(parallels._heights_get())
	selected_parallel_height = parallels_heights[0]
	for parallel_height in parallels_heights:
	    if parallel_height > selected_parallel_height and \
	      parallel_height + tooling_plate_dz <= vice_jaw_volume_dz:
		selected_parallel_height = parallel_height

	# Load up *skips_table* with the (row, column) pairs that are not to be drilled:
	skips_table = {}
	for skip in skips:
	    skip = tuple(skip)
	    skip_column = skip[0]
	    skip_row = skip[1]
	    assert skip_row in rows, "Skip row {0} not in rows {1}".format(skip_row, rows)
	    assert skip_column in columns, \
	      "Skip column {0} not in columns {1}".format(skip_column, columns)
	    skips_table[skip] = skip

	# Sort *rows* and *columns*:
	sorted_rows = sorted(tuple(rows))
	sorted_columns = sorted(tuple(columns))

	# Grab the minimum and maximum from *rows* and *columns*:
	minimum_row = sorted_rows[0]
	maximum_row = sorted_rows[-1]
	minimum_column = sorted_columns[0]
	maximum_column = sorted_columns[-1]

	# Make sure that we do not try to use a hole that does not exist:
	assert minimum_column >= 0, \
	  "Negative column {0} in {1} not allowed".format(minimum_column, columns)
	assert maximum_column < tooling_plate_columns, \
	  "Column {0} in {1} too big".format(maximum_column, columns)
	assert minimum_row >= 0, \
	  "Negative row {0} in {1} not allowed".format(minimum_row, rows)
	assert maximum_row < tooling_plate_rows, \
	  "Row {0} in {1} too big".format(maximum_row, rows)

	# Compute the *rows_spanned* and *columns_spanned*:
        rows_spanned = maximum_row - minimum_row
	columns_spanned = maximum_column - minimum_column

	# What we want is *part_dz* which is the part thickness along the Z axis when
	# it mounted in a vice.  The *top_surface_transform* reorients and translates
	# the part in its original 3 sapce orientation such that the center of the
        # "top" surface is located at the origin (0, 0, 0).  We know that mapping
        # *center_point* of the origin 3 space bounding box with *top_surface_transform*
	# will place the mapped point immedately under origin (i.e. (0, 0, -*part_dz*/2) .)
        # Thus, the following code:
	top_surface_transform = part._top_surface_transform
	center_point = top_surface_transform * part._bounding_box.c_get()
	half_part_dz = -center_point.z
	part_dz = 2 * half_part_dz

	# Figure out where upper left tooling hole is located:
	upper_left_x = -((2 * minimum_column + columns_spanned) * tooling_plate_hole_pitch) / 2
	upper_left_y = -((2 * minimum_row  + rows_spanned) * tooling_plate_hole_pitch) / 2
	zero = L()

	# Drill the *tooling_plate* holes at the intersections of *rows* and *columns*,
        # but ommitting any that are in *skips*:
	reverse_top_surface_transform = top_surface_transform.reverse()
	for row in rows:
	    for column in columns:
		column_row = (column, row)
		if not column_row in skips_table:
		    # Figure out the *start* and *stop* points for the drill as if the
		    # part is centered immediately under the origin (i.e. after
		    # *top_sufrace_transform* is applied to the entire *part*.):
		    x = upper_left_x + tooling_plate_hole_pitch * column
		    y = upper_left_y + tooling_plate_hole_pitch * row
		    start = P(x, y, zero)
		    stop = P(x, y, -part_dz)

		    # Drill here:
		    reversed_start = reverse_top_surface_transform * start
		    reversed_stop = reverse_top_surface_transform * stop
		    if tracing_detail >= 1:
			print(("{0}column={1} row={2} " +
			  "start={3:i} stop={4:i} rstart={5:i} rstop={6:i}").format(
			  indent, column, row, start, stop, reversed_start, reversed_stop))

		    # Drill the hole using *reversted_start* and *reverst_stop*:
		    part.hole("Tooling plate hole ({0}, {1})".format(row, column),
		      tooling_plate_drill_diameter, reversed_start, reversed_stop, "t",
		      tracing=tracing+1)

	# Now remember where to mount the part with the tooling plate:

	# Identify the location of (0, 0) *tooling_plate* hole:
	tooling_plate_columns_dx = tooling_plate_columns * tooling_plate_hole_pitch
	tooling_plate_rows_dy = tooling_plate_rows * tooling_plate_hole_pitch
	tooling_plate_edge_dx = (tooling_plate_dx - tooling_plate_columns_dx) / 2
	tooling_plate_edge_dy = (tooling_plate_dy - tooling_plate_rows_dy) / 2

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Part.tooling_plate_drill('{1}', '{2}', {3}, {4}, {5})".format(
	      indent, part._name, comment, rows, columns, skips))

    def xtooling_plate_mount(self, comment):
	""" *Part*: Cause the mounting plate that holds the *Part* object (i.e. *self*)
	    to be mounted in the vice using a dowel pin operation. """

	# Verify argument types:
	assert isinstance(comment, str)

	# Perform any requested *tracing*;
	tracing = self._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.tooling_plate_mount('{1}'".format(indent, self._name))

	if tracing >= 0:
	    print("{0}!!!!Part.tooling_plate_mount() is not implemented yet!!!!".format(indent))

	#jig_dy :@= part.jig_dy
	#half_jig_dy :@= half@(jig_dy)
	#vice_y :@= half_jig_dy + in@("1/4")
	##call d@(form@("jig_dy=%i% vice_y=%i%\n\") % f@(jig_dy) / f@(vice_y))
	#call reposition@(part, vice_y)
	#call dowel_pin@(part, "Mount on tooling plate")

	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Part.tooling_plate_mount('{1}'".format(indent, self._name))

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

    def vice_mount_with_extra(self,
      comment, top_surface, jaw_surface, extra_dx, extra_dy, tracing=-100000):
	""" *Part*: Cause the *Part* object (i.e. *self*) to be mounted in a vice with
	    *top_surface* facing upwards and *jaw_surface* mounted towards the rear vice jaw.
	    *comment* is the attached to any generated G-code.  *extra_dx* and *extra_dy*
	    is added around the *Part* object to allow for contouring and avoiding drilling
	    holes in to any parallel mounts in the vice under the part.  All transforms are
	    done with respect to the bounding box of the *Part* object.
	"""

	# Use *part* instead of self:
	part = self

	# Verify argument types:
	zero = L()
	assert isinstance(comment, str)
	assert isinstance(top_surface, str) and len(top_surface) == 1 and top_surface in "tbnsew"
	assert isinstance(jaw_surface, str) and len(jaw_surface) == 1 and jaw_surface in "tbnsew"
	assert isinstance(extra_dx, L) and extra_dx >= zero
	assert isinstance(extra_dy, L) and extra_dy >= zero

	part_name = part._name

	# Start any *tracing*:
	trace_detail = -1
	if tracing < 0:
 	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    trace_detail = 1
	    print("{0}=>Part.vice_mount_with_extra('{1}', '{2}', '{3}', '{4}', {5:i}, {6:i})".
	      format(indent, part_name, comment, top_surface, jaw_surface, extra_dx, extra_dy))

	# We only do work if we are in CNC mode:
	shop = part._shop_get()
	if shop._cnc_generate_get():
	    # The *x_axis*, *y_axis*, and *z_axis* are needed for rotations:
	    one = L(mm=1.0)
	    x_axis = P(one, zero, zero)
	    y_axis = P(zero, one, zero)
	    z_axis = P(zero, zero, one)

	    # These are the two rotation constants used to reorient the 
	    degrees90 = Angle(deg=90.0)
	    degrees180 = Angle(deg=180.0)

	    # Grab the 6 possible center surface points (and *center_point*) of the *bounding_box*
	    # for use below:
	    bounding_box = part._bounding_box
	    center_point = bounding_box.c_get()
	    t = bounding_box.t_get()
	    b = bounding_box.b_get()
	    n = bounding_box.n_get()
	    s = bounding_box.s_get()
	    e = bounding_box.e_get()
	    w = bounding_box.w_get()

	    # *top_axis* and *top_rotate* are used to rotate the desired bounding box surface
	    # facing up in the vice.  Leaving as *top_rotate_angle* as *None* means no top
	    # rotation is needed.  *top_point* matches *top_surface*:
	    top_point = None
	    top_axis = None
	    top_rotate_angle = None

	    # *jaw_axis* and *jaw_rotate* are used to rotate the desired bounding box surface
	    # facing towards the rear vice jaw.  Leaving as *jaw_rotate_angle * as *None* means no
	    # jaw rotaton needed.  *jaw_point* matches *jaw_surface*:
	    jaw_point = None
	    jaw_axis = z_axis
	    jaw_rotate_angle = None

	    # *left_dowel_point* specifies the *bounding_box* center surface point to the "left"
	    # of the bounding box center when mounted in the vice.  This is used for the
	    # left dowel pin operation.   If a right dowel pin is needed, it can be computed.
	    left_dowel_point = None

	    # Do a 24 (=6 surfaces x 4 jaw orientations) dispatch on *top_surface* and *jaw_suface*:
	    if 't' in top_surface:
		# Top surface of *bounding_box* facing up from vice:
		top_point = t
		if 'n' in jaw_surface:
		    # No rotation needed
		    jaw_point = n
		    left_dowel_point = w
		elif 's' in jaw_surface:
		    jaw_point = s
		    jaw_rotate_angle = degrees180
		    left_dowel_point = e
		elif 'e' in jaw_surface:
		    jaw_point = e
		    jaw_rotate_angle = degrees90
		    left_dowel_point = n
		elif 'w' in jaw_surface:
		    jaw_point = w
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = s
		else:
		    assert False
	    elif 'b' in top_surface:
		# Bottom surface of *bounding_box* facing up from vice:
		top_point = t
		top_axis = x_axis
		top_rotate_angle = degrees180
		if 'n' in jaw_surface:
		    jaw_point = n
		    jaw_rotate_angle = degrees180
		    left_dowel_point = e
		elif 's' in jaw_surface:
		    jaw_point = s
		    jaw_rotate_angle = None
		    left_dowel_point = w
		elif 'e' in jaw_surface:
		    jaw_point = e
		    jaw_rotate_angle = degrees90
		    left_dowel_point = s
		elif 'w' in jaw_surface:
		    jaw_point = w
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = n
		else:
		    assert False
	    elif 'n' in top_surface:
		# North surface of *bounding_box* facing up from vice:
		top_point = n
		top_axis = x_axis
		top_rotate_angle = degrees90
		if 't' in jaw_surface:
		    jaw_point = t
		    jaw_rotate_angle = -degrees180
		    left_dowel_point = e
		elif 'b' in jaw_surface:
		    jaw_point = b
		    jaw_rotate_angle = None
		    left_dowel_point = w
		elif 'e' in jaw_surface:
		    jaw_point = e
		    jaw_rotate_angle = degrees90
		    left_dowel_point = b
		elif 'w' in jaw_surface:
		    jaw_point = w
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = t
		else:
		    assert False
	    elif 's' in top_surface:
		# South surface of *bounding_box* facing up from vice:
		top_point = s
		top_axis = x_axis
		top_angle = -degrees90
		if 't' in jaw_surface:
		    jaw_point = t
		    jaw_rotate_angle = None
		    left_dowel_point = w
		elif 'b' in jaw_surface:
		    jaw_point = b
		    jaw_rotate_angle = degrees180
		    left_dowel_point = e
		elif 'e' in jaw_surface:
		    jaw_point = e
		    jaw_rotate_angle = degrees90
		    left_dowel_point = t
		elif 'w' in jaw_surface:
		    jaw_point = w
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = b
		else:
		    assert False
	    elif 'e' in top_surface:
		# East surface of *bounding_box* facing up from vice:
		top_point = e
		top_axis = y_axis
		top_angle = -degrees90
		if 't' in jaw_surface:
		    jaw_point = t
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = s
		elif 'b' in jaw_surface:
		    jaw_point = b
		    jaw_rotate_angle = degrees90
		    left_dowel_point = n
		elif 'n' in jaw_surface:
		    jaw_point = n
		    jaw_rotate_angle = None
		    left_dowel_point = t
		elif 's' in jaw_surface:
		    jaw_point = s
		    jaw_rotate_angle = -degrees180
		    left_dowel_point = b
		else:
		    assert False
	    elif 'w' in top_surface:
		# West surface of *bounding_box* facing up from vice:
		top_point = w
		top_axis = y_axis
		top_angle = degrees90
		if 't' in jaw_surface:
		    jaw_point = t
		    jaw_rotate_angle = -degree90
		    left_dowel_point = n
		elif 'b' in jaw_surface:
		    jaw_point = b
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = s
		elif 'n' in jaw_surface:
		    jaw_point = n
		    jaw_rotate_angle = None
		    left_dowel_point = b
		elif 's' in jaw_surface:
		    jaw_point = s
		    jaw_rotate_angle = degrees180
		    left_dowel_point = t
		else:
		    assert False
	    else:
		assert False

	    # Perform any requested *tracing*:
	    if trace_detail >= 2:
		print("{0}center_point={1:i} top_point={2:i} jaw_point={3:i} dowel_point={4:i}".
		  format(indent, center_point, top_point, jaw_point, left_dowel_point))

	    # Grab the *vice* and its associated *vice_jaw volume*:
	    vice = shop._vice_get()
	    vice_jaw_volume = vice._jaw_volume_get()
	    vice_jaw_dx = vice_jaw_volume.x
	    vice_jaw_dy = vice_jaw_volume.y
	    vice_jaw_dz = vice_jaw_volume.z

	    # Grab the *parallels* and associated *parallels_height*:	
	    parallels = vice._parallels_get()
            parallels_heights = parallels._heights_get()

	    # Figure out both the *smallest_parallel_height* and the *largest_parallel_height*:
	    smallest_parallel_height = L(inch=12345.0)
	    largest_parallel_height = L()
            for parallel_height in parallels_heights:
		smallest_parallel_height = smallest_parallel_height.minimum(parallel_height)
		largest_parallel_height  =  largest_parallel_height.maximum(parallel_height)

	    # Now compute the *half_part_dx*, *half_part_dy*, and *half_part_dz* as
	    # augmented by *extra_dx* and *extra_dy*:
	    half_part_dx = center_point.distance(left_dowel_point) + extra_dx / 2
	    half_part_dy = center_point.distance(jaw_point)        + extra_dy / 2
	    half_part_dz = center_point.distance(top_point)

	    # Compute *part_dx*, *part_dy*, and *part_dz* and perform any requested *tracing*:
	    part_dx = 2 * half_part_dx
	    part_dy = 2 * half_part_dy
	    part_dz = 2 * half_part_dz
	    if trace_detail >= 2:
		print("{0}part_dimensions: dx={1:i} dy={2:i} dz={3:i} volume={4:i}".
		  format(indent, part_dx, part_dy, part_dz, bounding_box.volume_get()))

	    # Now figure out *top_surface_z*, which the distance from the vice origin to the
            # center of the top surface of the *part*.  We want to mount such that the top surface
	    # is as near to the top jaw edge as possible.  However, sometimes the parallels will
	    # not let us get that low:
	    selected_parallel_height = smallest_parallel_height
	    for parallel_height in parallels_heights:
		if parallel_height > selected_parallel_height and \
 		  parallel_height + part_dz <= vice_jaw_dz:
                    selected_parallel_height = parallel_height

	    # Now we can compute *cnc_top_surface_z* and *cnc_xy_rapid_safe_z*:
	    cnc_top_surface_z = selected_parallel_height + part_dz
	    cnc_xy_rapid_safe_z = cnc_top_surface_z + L(inch=0.5)

	    # Now we can figure out how we are going to position the part in the vice in Y axis.
            # This allows us to compute *cnc_vice_translate_point*:
	    if part_dx > vice_jaw_dx:
		# The *part* is longer in X axis than vice jaw; center the part in the vice:
		vice_translate_point = P(vice_jaw_dx/2, -half_part_dy, cnc_top_surface_z)
	    else:
		# The *part* is shorter in the X axis than the vice jaw, position part to left
                # edge of vice:
		cnc_vice_translate_point = P(half_part_dx, -half_part_dy, cnc_top_surface_z)

	    # Now we can compute the *cnc_transform*, which rotates the *part* bounding box
	    # so that it is oriented correctly in the vice with the top surface located at
            # *cnc_top_surface_z*:
	    top_surface_transform = Transform()
	    top_surface_transform = top_surface_transform.translate(
	      "move to vice origin", -top_point)
	    if isinstance(top_rotate_angle, Angle):
		top_surface_transform = top_surface_transform.rotate(
		  "vice surface to top", top_axis, top_rotate_angle)
	    if isinstance(jaw_rotate_angle, Angle):
		top_surface_transform = top_surface_transform.rotate(
		  "rotate jaw surface north", jaw_axis, jaw_rotate_angle)
	    part._top_surface_transform = top_surface_transform
	    cnc_transform = top_surface_transform.translate(
	      "position in vice", cnc_vice_translate_point)

	    # Remember *cnc*_transform* in *part* so the rest of the CNC system can use it:
	    part._cnc_transform = cnc_transform

	    # Find the *tool_dowel_pin* and grab some values out of it:
	    tool_dowel_pin = part._tools_dowel_pin_search(tracing + 1)
	    assert isinstance(tool_dowel_pin, Tool_Dowel_Pin), "No dowel pin tool found"
	    tip_depth = tool_dowel_pin._tip_depth_get()
	    maximum_z_depth = tool_dowel_pin._maximum_z_depth_get()
	    feed_speed = tool_dowel_pin._feed_speed_get()
	    spindle_speed = tool_dowel_pin._spindle_speed_get()
	    diameter = tool_dowel_pin._diameter_get()
	    if trace_detail >= 2:
                print("{0}diameter={1:i} maximum_z_depth={2:i} tip_depth={3:i}".format(
		  indent, diameter, maximum_z_depth, tip_depth))

	    # Now we figure out *cnc_dowel_point* which is the location in CNC space where
	    # the tip of the dowel pin is moved to:
	    cnc_dowel_point_x = cnc_vice_translate_point.x - half_part_dx - diameter/2
	    cnc_dowel_point_y = cnc_vice_translate_point.y
	    if part_dz >= maximum_z_depth:
		cnc_dowel_point_z = cnc_top_surface_z - maximum_z_depth
	    else:
		cnc_dowel_point_z = cnc_top_surface_z - part_dz - tip_depth
	    cnc_dowel_point = P(cnc_dowel_point_x, cnc_dowel_point_y, cnc_dowel_point_z)
	    if trace_detail >= 2:
		print("{0}vtp.x={1:i} - hpx={2:i} - rad={3:i} cdpx={4:i}".format(
		  indent, cnc_vice_translate_point.x, half_part_dx, diameter/2, cnc_dowel_point_x))

	    # *cnc_plunge_point* is where the dowel pin initially comes down.:
	    cnc_plunge_point = cnc_dowel_point - P(L(inch=1.000), zero, zero)
	    if trace_detail >= 2:
                print("{0}cnc_dowel_point={1:i} cnc_plunge_point={2:i}".format(
		  indent, cnc_dowel_point, cnc_plunge_point))

	    # Create the *operation_dowel_pin* *Operation* and load in into *part*:
	    order = Operation.ORDER_DOWEL_PIN
	    operation_dowel_pin = Operation_Dowel_Pin(part, "dowel_pin", 0, tool_dowel_pin, order,
	      None, feed_speed, spindle_speed, diameter, cnc_dowel_point, cnc_plunge_point,
	      cnc_top_surface_z, cnc_xy_rapid_safe_z, tracing + 1)
	    part._operation_append(operation_dowel_pin)

	if tracing >= 0:
	    print("{0}<=Part.vice_mount_with_extra('{1}', '{2}', '{3}', '{4}', {5:i}, {6:i})".
	      format(indent, part_name, comment, top_surface, jaw_surface, extra_dx, extra_dy))

	# This is the old code that is no longer used in this routine
	##FIXME: What does this do?!!!
	## Before we get too far, flush any pending screw holes from the previous mounting:
	##if part._top_surface_set:
	##    part.screw_holes(True)
	#
	## Let's create the *z_axis*:
	#one = L(mm=1.0)
	#z_axis = P(zero, zero, one)
	#
	## We need to create a *cnc_transform* matrix that will translate points from the
	## conceptual design orientation into an orientation that is horizontal to the
	## vice plane.  We start by making *top_surface_point* the virtual origin of the
	## vice orientation:
	#cnc_transform = Transform()
	#cnc_transform = cnc_transform.translate("vice_mount top translate", -top_surface_point)
	#
	## Compute the *projection_axis* that is the axis of projection used to project
	## X/Y/Z points in the conceptual design orientation into a planar X/Y orientation
	## for the vice.  This is done be creating a *north_axis* and a *west_axis* and
	## using the cross product operation to create *projection_axis*:
	#north_axis = north_point - surface_point
	#west_axis = west_point - surface_point
	#unnormalized_projection_axis = north_axis.cross_product(west_axis)
	#projection_axis = unnormalized_projection_axis.normalize()
	#if tracing >= 0:
	#    print("{0}north_axis={1:m} west_axis={2:m} unnormalized projection_axis={3:m}".
	#      format(indent, north_axis, west_axis, unnormalized_projection_axis))
	#    print("{0}normailized projection_axis={1:m}".format(indent, projection_axis))
	#
	## Now that we have the *projection_axis* we need to decide whether or not the
	## points in the conceptual design orientation need to be rotated to be in the
	## vice orientation.  This only needs to be done if *projection_axis* is *NOT*
	## aligned with *z_axis*:
	#if projection_axis != z_axis:
	#    # Now we want to adjust the *cnc_transform* to rotate from *projection_axis*
	#    # orientation to *z_axis* orientation.  To do this we need a *rotation_axis*
	#    # a *rotation_angle*.  The *rotation_axis* is computed using a cross product
	#    # between the *projection_axis* and the *z_axis*.  The *rotation_angle* is
	#    # computed using *angle_between*() routine for the *z_axis* and the *projection_axis*.
	#    rotation_axis = z_axis.cross_product(projection_axis).normalize()
	#    rotation_angle = z_axis.angle_between(projection_axis)
	#    cnc_transform = cnc_transfrom.rotate("vice_mount rotate",
	#      rotation_axis, rotation_angle)
	#    if tracing >= 0:
	#	print("{0}rotation_axis={1:m} rotation_angle{2:d}".
	#	  format(' ' * tracing, rotation_axis, rotation_angle))
	#if tracing >= 0 and trace_detail >= 1:
	#    print("{0}position=\n{1}".format(indent, position))
	#
	## Now do the the final translation of *cnc_transform*:
	#half_dx = surface_point.distance(north_point)
	#half_dy = surface_point.distance(west_point)
	#cnc_transform = cnc_transform.translate("vice position translate",
	#  P(half_dx + extra_dx/2, half_dy + extra_dy/2, zero))
	#
	##FIXME: We need to figure out how deep the bounding box is in this orientation!!!
	#bounding_box = part._bounding_box
	#half_dz = surface_point.distance(bounding_box.c_get())
	#
	## Save the *cnc_transform* into *code*:
	#shop = part._shop_get()
	#code = shop._code_get()
	#code._cnc_transform_set(cnc_transform)
	# Wrap up any *tracing*:

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

    def __init__(self, up, name):
	""" *Fastener*: """

	# Verify argument types:
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)

	# Initialize the *Part* super class:
	Part.__init__(self, up, name)
	
	# Disable CNC generation:
	self.cnc_suppress()

	zero = L()
	self.comment_s = name
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
	diameter = self.major_diameter_l
	self.cylinder(self.comment_s, self.material, self.color,
	  diameter, diameter, self.start_p, self.end_p, 16, Angle(deg=0.0), "", "")

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

	    part.hole("'{0} Drill'".format(self.comment_s), diameter, start, end, "t")

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
		if direction.length() <= zero and EZCAD3.update_count_get() == 0:
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
		    part.hole("Flat Head:" + self.comment_s,
		      2 * flat_head_diameter, flat_head_start,
		      flat_head_end, "", tracing = trace + 1)

	    if self.hex_insert_b:
		direction = end - start
		direction_length = direction.length()
		if direction_length <= zero and EZCAD3.update_count_get() == 0:
		    print("Fastener.drill(): part={0}: zero length screw".
		      format(part))
		else:
		    normalized_direction = direction.normalize()
		    nut_height = self.nut_height_l
		    insert_end = end - (normalized_direction * nut_height._mm)
		    #print("end={0} start={1} dir={2} dir_len={3} nut_hght={4}".
		    #  format(end, start, direction, direction_len, nut_hght))
		    #print("insert_end = {0}".format(insert_end))

		    part.hole("Hex Insert:" + self.comment_s,
		      self.hex_nut_tip_width_l, end, insert_end, "f",
                      sides=6, sides_angle=self.sides_angle_a)

	if trace >= 0:
	    print("{0}<=Fastener.drill({1}, select='{2}')".
	      format(' ' * trace, part, select))


# *Code* class:

class Code:
    """ *Code*: The *Code* class outputs G-codes, DXF-files, and VRML based tool paths.
	Most of the focus is on G-codes though.
    """

    def __init__(self):
	""" *Code*: Initialize the new *Code* object (i.e. *self*).
	"""

	# Use *code instead of *self*:
	code = self

	zero = L()

	# Load up *self*:
	code._cnc_transform = None	# Transform for mapping part points into CNC space
	code._code_stream = None	# Code stream to write G-codes to
	code._command_started = False	# At begining of RS274 command line
	code._command_chunks = []	# RS274 line broken into space separated chunks
	code._comment_chunks = []	# RS274 comment broken into space separated chunks
	code._dxf = ""			# Text for dxf file
	code._dxf_lines = []		# List of lines to write out to .dxf file
	code._dxf_x_offset = zero	# DXF X offset
	code._dxf_y_offset = zero	# DXF Y offset
	code._is_laser = False		# True if tool is a laser
	code._part_wrl_file = None	# Place to write unified part path .wrl file
	code._top_surface_z = None	# Top surface of part in CNC coordinates
	code._xy_rapid_safe_z = None	# CNC Z altitude where XY rapids are allowed

	# Construct the *g1_table*:
	g1_values_list = (0, 1, 2, 3, 33, 38, 73, 76, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89)
	g1_table = {}
	for g1_value in g1_values_list:
	    g1_table[g1_value] = g1_value
	code._g1_table = g1_table

	# The stuff below is RS-274 mode variables:
	code._f = Speed()		# Feedrate
	code._g1 = 0			# (G0-3, 33, 38.x, 73, 76, 80-89)
	code._g2 = 0			# (G17-19)
	code._g3 = 0			# (G7-8)
	code._g4 = 0			# (G90-91)
	code._g5 = 0			# (G93-94)
	code._g6 = 0			# (G20-21)
	code._g7 = 0			# (G40-42)
	code._g8 = 0			# (G43, 49)
	code._g9 = 0			# (G98-99)
	code._g10 = 0			# (G54-59)
	code._g11 = 0			# (G4)
	code._h = zero			# H tool offset index
	code._i = zero			# I coordinate
	code._j = zero			# J coordinate
	code._m1 = 0			# (M0-2, 30, 60)
	code._m2 = 0			# (M6)
	code._m3 = 0			# (M3-5)
	code._m4 = 0			# (M7-9)
	code._m5 = 0			# (M48-49)
	code._p = Time()		# G4
	code._q = zero			# Peck depth
	code._r0 = zero			# Radius cycle R
	code._r1 = zero			# Drill cycle R
	code._s = Hertz()		# Spindle revolutions
	code._x = zero			# X coordinate
	code._y = zero			# Y coordinate
	code._z = zero			# Z coordinate
	code._z1 = zero 		# Z coordinate
	code._z_safe_f = Speed()	# Feed to perform z safe operation at
	code._z_safe_pending = False	# {true}=>need to do z safe move
	code._z_safe_s = Hertz()	# Speed to perform z safe operation at

	# Some color index constants:
	code._vrml_motion_color_rapid = 0
	code._vrml_motion_color_cutting = 1
	code._vrml_motion_color_retract = 2

	# Reset VRML fields here.  Note: routine acessses some of the G-code variable,
        # so it must be called last:
	code._vrml_reset()

    def _command_begin(self):
	""" *Code*: Start a new RS274 command in the *Code* object (i.e. *self*). """

	# Use *code* instead of *self*:
	code = self

	assert not code._command_started, "Previous RS274 command was not ended."
	code._command_started = True

	# Remember some values for VRML line path drawing:
	code._vrml_start_r0 = code._r0
	code._vrml_start_x = code._x_value()
	code._vrml_start_y = code._y_value()
	code._vrml_start_z = code._z
	code._vrml_motion_color = -1

    def _cnc_transform_set(self, cnc_transform):
        """ *Code*: Set the CNC transform for the *Code* object (i.e. *self*) to *cnc_transform*.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(cnc_transform, Transform)

	# Stuff *cnc_transform* into *code*:
	code._cnc_transform = cnc_transform

    def _command_end(self):
	""" *Code*: End the current RS274 in the *Code* object (i.e. *self*). """

	# Use *code* instead of *self*:
        code = self

	# Make sure we have started a command:
	assert code._command_started, "Not currently in a RS274 command"

	# Grab out both the *command_chunks* and the *comment_chunks*:
	command_chunks = code._command_chunks
	comment_chunks = code._comment_chunks

	# If we have any *comment_chunks*, tack the onto the end of *command_chunks*:
	if len(comment_chunks):
	    command_chunks += ["("] + comment_chunks + [")"]

	# Construct a space separated *command*:
	command = " ".join(code._command_chunks)
	#print("command='{0}'".format(command))

	# Write *comand* to *code_stream*:
	code_stream = code._code_stream
	code_stream.write(command)
	code_stream.write("\n")

	# Clear out *command_chunks* and *comment_chunks* for the next command:
	del command_chunks[:]
	del comment_chunks[:]

	# Mark that we ended the current command:
	code._command_started = False

	# Figure out what *color_index* to use draw tool path in VRML:
	color_index = code._vrml_motion_color

	# Add some information to draw the tool path in VRML:
	if color_index >= 0:
	    start_index = \
	      code._vrml_point(code._vrml_start_x, code._vrml_start_y, code._vrml_start_z)
	    end_index = code._vrml_point(code._x_value(), code._y_value(), code._z)
	    code._vrml_point_indices.append(start_index)
	    code._vrml_point_indices.append(end_index)
	    code._vrml_point_indices.append(-1)
	    code._vrml_color_indices.append(color_index)

    def _comment(self, comment):
	""" *Code*: Output *comment* to the *Code* object (i.e. *self*). Any parentheses
	    in *comment* are converted to square brackets.
	"""

	# Use *code* instead of *self*:
        code = self

	# Verify argument types:
	assert isinstance(comment, str)

	# Replace open/close parenthesis in *comment* with square brackets:
	comment = comment.replace('(', '[').replace(')', ']')

	# Add it on to the comment
	code._comment_chunks.append(comment)

    def _configure(self, tool, vice_x, vice_y):
	""" *Code*: Configure the *Code* object (i.e. *self*) to use
	    *tool*, *vice_x*, and *vice_y*.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(tool, Tool)
	assert isinstance(vice_x, L)
	assert isinstance(vice_y, L)

	# Initialize *code*:
	code._tool = tool
	code._is_laser = tool._is_laser_get()
	code._vice_x = vice_x
	code._vice_y = vice_y

    def _contour(self, contour,
      plunge_offset, contour_offset, tool_radius, clockwise, z, feed_speed, spindle_speed, tracing):
	""" *Code*: Ggenerate the G-code for *contour* using the *Code* object (i.e. *self*).
	    *plunge_offset* is the offset from contour at which to lower the tool.
	    *contour_offset* is an additional +/- offset from the contour.  *tool_radius*
	    is the tool radius to use.  If *clockwise* is *True*, *bends* are traversed
	    in a clockwise direction  otherwise a counter-clockwise traversal occurs.
	    *z* specifies the depth.  *feed_speed* specifies the feedrate and *spindle_speed*
	    specifies the spindle speed.
	"""

	# Use *code* instead of *self*:
        code = self

	# Verify argument types:
	zero = L()
	assert isinstance(contour, Contour)
	assert isinstance(plunge_offset, L) and plunge_offset >= zero
	assert isinstance(contour_offset, L) and contour_offset >= zero
	assert isinstance(tool_radius, L) and tool_radius >= zero
	assert isinstance(clockwise, bool)
	assert isinstance(z, L)
	assert isinstance(feed_speed, Speed)
	assert isinstance(spindle_speed, Hertz)
	assert isinstance(tracing, int)

	assert feed_speed != Speed()

	# Start performing any *tracing*:
	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 0
	    indent = ' ' * tracing
	    print("{0}=>Code._contour(*, po={1:i} co={2:i} tr={3:i} cl={4} z={5:i} *)".format(
	      indent, plunge_offset, contour_offset, tool_radius, clockwise, z))

	# Grab *bends* from *contour*:
	bends = contour._bends_get()
	bends_size = len(bends)

	# Add in *tool_radius* to the exising *contour_offset*:
	total_offset = contour_offset + tool_radius

	# Search for *west_most_bend* with the X coordinate:
	west_most_index = 0
	west_most_bend = bends[west_most_index]
	west_most_x = west_most_bend._point_get().x
	for index, bend in enumerate(bends):
	    bend_point_x = bend._point_get().x
            if bend_point_x < west_most_x:
		west_most_index = index
		west_most_bend = bend
		west_most_x = bend_point_x
	if tracing_detail >= 1:
	    print("{0}".format(indent, ))

	# Make sure that *west_most_bend* is not an inside bend (i.e. it should be impossible):
	assert not west_most_bend._is_inside_get(), \
	  "west_most_bed of contour '{0}' is inside".format(contour._name_get())

	# Compute the X/Y location to drop the mill bit down:
	contour_is_clockwise = contour._is_clockwise_get()
	west_most_radius = west_most_bend._radius_get()
	start_offset = total_offset + west_most_radius
	plunge_arc_offset = start_offset + plunge_offset
	if contour_is_clockwise == clockwise:
	    start_tangent = west_most_bend._incoming_tangent_compute(start_offset)
	    plunge_tangent = west_most_bend._incoming_tangent_compute(plunge_arc_offset)
	else:
	    start_tangent = west_most_bend._outgoing_tangent_compute(-start_offset)
	    plunge_tangent = west_most_bend._outgoing_tangent_compute(-plunge_arc_offset)
	if tracing_detail >= 1:
	    print("{0}plunge_x={1:i}, plunge_y={2:i}".
	      format(indent, plunge_tangent.x, plunge_tangent.y))
    
	# This routine starts and ends at (*plunge_x*, *plunge_y*).
	# If we are not already at (*plunge_x*, *plunge_y*) we need
	# to get there safely:
	code._xy_rapid(plunge_tangent.x, plunge_tangent.y)

	# Now we get to down to the correct Z level:
	code._z_feed(feed_speed/2, spindle_speed, z, "Contour")

	# Now feed to the start point using a circular motion so that there is no divit
	# at the plunge point:
	if contour_is_clockwise == clockwise:
	    code._xy_ccw_feed(feed_speed, spindle_speed,
	      tool_radius * 1.001, start_tangent.x, start_tangent.y)
	else:
	    code._xy_cw_feed(feed_speed, spindle_speed,
	      tool_radius * 1.001, start_tangent.x, start_tangent.y)

	# We iterate across all of the corners.  We need to visit
	# the first corner one last time at the end; hence we
	# iterate *size* + 1 times through the loop regardless of
	# whether we go clockwise or count-clockwise:
	if contour_is_clockwise == clockwise:
	    # Clockwise (climb) Cut:
	    if tracing_detail >= 0:
		print("{0}clockwise".format(indent))
    
	    # Iterate *bends_size* + 1 times:
	    for bends_index in range(bends_size):

		# Fetch a {bend}:
		index = bends_index % bends_size
		bend = bends[index]
		bend_radius = bend._radius_get()

		if tracing_detail >= 1:
		    bend_name = bend._name_get()
		    bend_point = bend._point_get()
		    print("{0}bend_name='{1}' bend_point='{2}".
		      format(indent, bend_name, bend_point))

		# Compute the *arc_radius*, *arc_start*, and *arc_end* depending upon whether
		# *bend* *is_inside* or not.  Offset everything by *total_offset*:
		is_inside = bend._is_inside_get()
		if is_inside:
		    # For an inside bend, *bend_radius* is decreased by *total_offset*:
		    arc_radius = bend_radius - total_offset

		    # For inside bends, we need to negate the *arc_radius*:
		    arc_start = bend._incoming_tangent_compute(-arc_radius)
		    arc_end = bend._outgoing_tangent_compute(-arc_radius)
		    code._xy_feed(feed_speed, spindle_speed, arc_start.x, arc_start.y)
		    code._xy_ccw_feed(feed_speed, spindle_speed, arc_radius, arc_end.x, arc_end.y)
		else:
		    # For an outside bend, *bend_radius* is increased by *total_offset*:
		    arc_radius = bend_radius + total_offset

		    # For outside bends, no further adjustments are needed:
		    arc_start = bend._incoming_tangent_compute(arc_radius)
		    arc_end = bend._outgoing_tangent_compute(arc_radius)
		    code._xy_feed(feed_speed, spindle_speed, arc_start.x, arc_start.y)
		    code._xy_cw_feed(feed_speed, spindle_speed, arc_radius, arc_end.x, arc_end.y)

		if tracing_detail >= 1:
		    print("{0}arc_radius={1:i} arc_start={1:2} arc_end={3:i}".
		      format(indent, arc_radius, arc_start, arc_end))
		    print("")
	else:
	    # Counter clockwise:
	    if tracing_detail >= 1:
		print("{0}counter-clockwise".format(indent))

	    for bends_index in range(bends_size):
		# Fetch a {bend}:
		index = (bends_size - bends_index) % bends_size
		bend = bends[index]
		bend_radius = bend._radius_get()

		# Compute the *incoming_tangent* and *outgoing_tangent* based on the
		# whether or not the *bend* *is_inside* or not, the *bend_radius* and
		# the *tool_radius*:
		is_inside = bend._is_inside_get()
		if is_inside:
		    # For inside corners, the *bend_radius* is decreased by *total_offset*:
		    arc_radius = bend_radius - total_offset

		    # For counter-clockwise, we swap incoming with outgoing but leave
		    # *arc_radius* unchanged (i.e. positive):
		    arc_end = bend._incoming_tangent_compute(arc_radius)
		    arc_start = bend._outgoing_tangent_compute(arc_radius)
		    code._xy_feed(feed_speed, spindle_speed, arc_start.x, arc_start.y)
		    code._xy_ccw_feed(feed_speed, spindle_speed, arc_radius, arc_end.x, arc_end.y)
		else:
		    # For outside corners, the *bend_radius* is increased by *total_offset*:
		    arc_radius = bend_radius + total_offset

		    # For counter-clockwise, we swap incoming with outgoing and we invert
		    # *arc_radius*:
		    arc_end = bend._incoming_tangent_compute(-arc_radius)
		    arc_start = bend._outgoing_tangent_compute(-arc_radius)
		    code._xy_feed(feed_speed, spindle_speed, arc_start.x, arc_start.y)
		    code._xy_cw_feed(feed_speed, spindle_speed, arc_radius, arc_end.x, arc_end.y)

		if tracing_detail >= 1:
		    print("{0}arc_radius={1:i} arc_start={1:2} arc_end={3:i}".
		      format(indent, arc_radius, arc_start, arc_end))
		    print("")
	    
	# Return back to the (*start_tangent.x*, *start_tangent.y*):
	code._xy_feed(feed_speed, spindle_speed, start_tangent.x, start_tangent.y)

	# Now feed to back to the plunge point in a circular motion to avoid diviting the
	# contour as we retract:
	if contour_is_clockwise == clockwise:
	    code._xy_ccw_feed(feed_speed, spindle_speed,
	      tool_radius * 1.001, plunge_tangent.x, plunge_tangent.y)
	else:
	    code._xy_cw_feed(feed_speed, spindle_speed,
	      tool_radius * 1.001, plunge_tangent.x, plunge_tangent.y)

	if tracing >= 0:
	    print("{0}=>Code._contour(*, po={1:i} co={2:i} tr={3:i} cl={4} z={5:i} *)".format(
	      ' ' * tracing, plunge_offset, contour_offset, tool_radius, clockwise, z))

    def _dxf_angle_append(self, group_code, value):
	""" *Code*: Append {group_code} and {value} to the current DXF entity in the
	    *Code* object (i.e. *self*).
	"""
    
	# Use *code* instead of *self*:
        code = self

	# Verify argument types:
	assert isinstance(group_code, int) and int >= 0
	assert isinstance(value, Angle)

	code._dxf_append("{0}\n{1:d}\n".format(group_code, value), "dxf_angle_append")
    
    def _dxf_append(self, text, from_routine):
	""" *Code*: Append *text* to the DXF component the *Code* object (i.e. *self*).
	    *from_routine* is used for debugging only.
	"""
    
	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(text, str)
	assert isinstance(from_routine, str)

	# Append *text* to *dxf_lines*:
	#assert code._is_laser
	code._dxf_lines.append(text)
    
    def _dxf_arc_append(self, clockwise, end_x, end_y, radius):
	""" *Code*: Generate a DXF arc entity that draws an arc from the current X/Y location
	    the *Code* object (i.e. *self*) to (*end_x*,*end_y*) with radius of *radius*.
	    The arc is drawn clockwise if *clockwise* is *True* and counter-clockwise otherwise.
	"""
    
	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(clockwise, bool)
	assert isinstance(end_x, L)
	assert isinstance(end_y, L)
	assert isinstance(radius, L)

	# Perform any requested *tracing*:
	tracing = -1000000
	#tracing = 0
	detail_level = -1
	if tracing >= 0:
	    #detail_level = 3
	    indent = ' ' * tracing
	    print("{0}=>Code._dxf_arc_append@(cw={1}, end_x={2:i}, end_y={3:i}, radius={4:i})".
	      format(indent, clockwise, end_x, end_y, radius))

	# Compute some constants:
	degrees0 = Angle()
	degrees90 = Angle(deg=90.0)
	degrees180 = Angle(deg=180.0)
	degrees360 = Angle(deg=360.0)
    
	# Grab the start position ({x1}, {y1}):
	start_x = code._x_value()
	start_y = code._y_value()
    
	# What we are going to do here is fit a circle to S=(*start_x*,*start_y) and
	# E=(*end_x*,*end_y*).  We need to find the center of the circle
	# C=(*center_x*, *center_y*).  What we know is that the radius r = |SC| = |EC|.
	# 
	# We draw a line segment from S to E and bisect it at point M.  The line that
	# goes through from MC will be perpendicular to the line through S and E.
	# This gives us a triangle SMC, where the angle <SMC is 90 degrees.
	# Once we have computed |MC| we can compute C.
	#
	# From the Pathagean theorem:
	#
	#        |MS|^2 +|MC|^2 = |SC|^2 = r^2                     (1)
	#
	# Thus,
        #
	#        |M|^2 = r^2C - |MS|^2                             (2)
	#
	#        |MC| = sqrt(r^2 - |MS|^2)                         (3)
	#

	# Let's get *s* and *e* defined:
	zero = L()
	s = P(start_x, start_y, zero)
	e = P(end_x, end_y, zero)
	
	# Find the midpoint *m*:
	m = (s + e) / 2

	# Now compute *ms_length* and *ms_angle*
	ms = m - s
	ms_length = ms.length()
	ms_angle = ms.xy_angle()

	# Now compute *mc_length*
	ms_mm = ms_length.millimeters()
	radius_mm = radius.millimeters()
	mc_mm = math.sqrt(radius_mm * radius_mm - ms_mm * ms_mm)
	mc_length = L(mm=mc_mm)

	# Now compute *cm_angle* being 90 degrees of of *ms_angle*:
	if clockwise:
	    cm_angle = ms_angle - degrees90
	else:
	    cm_angle = ms_angle + degrees90
	cm_angle = cm_angle.normalize()

	# Now compute center *c*:
	center_x = m.x + mc_length.cosine(cm_angle)
	center_y = m.y + mc_length.sine(cm_angle)
	c = P(center_x, center_y, zero)

	# Now we need the *start_angle* and *end_angle* for DXF entity:
	start_angle = (s-c).xy_angle()
	end_angle = (e-c).xy_angle()

	# Generate the required ARC entity:
	code._dxf_entity_start("ARC")
	code._dxf_xy_append(0, center_x, center_y, "_dxf_arc_append")
	code._dxf_length_append(40, radius)
	if clockwise:
	    code._dxf_angle_append(50, end_angle)
	    code._dxf_angle_append(51, start_angle)
	else:
	    code._dxf_angle_append(50, start_angle)
	    code._dxf_angle_append(51, end_angle)
	code._dxf_entity_stop()
    
	if tracing >= 0:
	    print("{0}<=Code._dxf_arc_append@(cw={1}, end_x={2:i}, end_y={3:i}, radius={4:i})".
	      format(indent, clockwise, end_x, end_y, radius))

    def _dxf_circle(self, x, y, radius):
	""" *Code*: Append a circle DXF entity to the *Code* object (i.e. *self*)
	    with a center of (*x*,*y*) and a radius of *radius*.
	"""
    
	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(radius, L)

	code._dxf_entity_start("CIRCLE")
	code._dxf_xy_append(0, x, y, "_dxf_circle")
	code._dxf_length_append(40, radius)
	code._dxf_entity_stop()
    
    def _dxf_entity_start(self, name):
	""" *Code*: Start a DXF entity named *name* with the *Code* object (i.e. *self*).  """

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(name, str)
    
	code._dxf_append("0\n{0}\n".format(name), "_dxf_entity_start")
	code._dxf_integer_append(8, 2)
	code._dxf_integer_append(62, 0)
    
    def _dxf_entity_stop(self):
	""" *Code* Terminate the current DXF entity in the *Code* object (i.e. *self*). """

	# Use *code* instead of *self*:
	code = self

	# Nothing to do right now:
	pass
    
    
    def _dxf_content_avaiable(self):
	""" *Code*: Return *True* if the *Code* object (i.e. *self*) has any DXF file content. """

	# Use *code* instead of *self*:
	code = self

	return len(code._dxf_lines) > 0

    def _dxf_integer_append(self, group_code, value):
	""" *Code*: Append *group_code* and *value* to the current DXF entity the *Code* object
	    (i.e. *self*).
	"""
    
	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(group_code, int) and group_code >= 0
	assert isinstance(value, int)

	code._dxf_append("{0}\n{1}\n".format(group_code, value), "_dxf_integer_append")
    
    def _dxf_length_append(self, group_code, value):
	""" *Code*: Append *group_code* and *value* to the current DXF entity in the
	    *Code* object (i.e. *self*).
	"""
    
	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(group_code, int) and group_code >= 0
	assert isinstance(value, L)

	code._dxf_append("{0}\n{1:m}\n".format(group_code, value), "_dxf_length_append")
    
    def _dxf_write(self, dxf_file):
	""" *Code*: Write the DXF file lines associated with the *Code* object (i.e. *self*)
	    to *dxf_file*.
	"""

	# Use *code* instead of *self*:
	code = self

	for dxf_line in code._dxf_lines:
	    dxf_file.write(dxf_line)
	code._dxf_lines = []
	
    def _dxf_xy_append(self, offset, x, y, from_routine):
	""" *Code*: Append (*x*,*y*) values to the DXF entity of the *Code* object (i.e. *self*),
	    where *offset* specfies the increment added to the 10, 20, 30 record fields.
	"""
    
	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(offset, int) and offset >= 0
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(from_routine, str)

	#call d@(form@("dxf_xy_append@(off=%d%, x+%i%, y=%i%, from=%v%)\n\") %
	#  f@(offset) % f@(x) % f@(y) / f@(from))
    
	code._dxf_length_append(10 + offset, x + code._dxf_x_offset)
	code._dxf_length_append(20 + offset, y + code._dxf_y_offset)
	code._dxf_length_append(30 + offset, L())
    
    def _dxf_xy_offset_set(self, dxf_x_offset, dxf_y_offset):
	""" *Code*: Set the DXF offset X field of the *RS724* object (i.e. *self*) """

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(dxf_x_offset, L)
	assert isinstance(dxf_y_offset, L)
	code._dxf_x_offset = dxf_x_offset
	code._dxf_y_offset = dxf_y_offset

    def _dxf_y_offset_get(self):
	""" *Code*: Return the DXF offset Y field of the *Code* object (i.e. *self*) """

	# Use *code* instead of *self*:
	code = self

	return code._dxf_y_offset

    def _finish(self, tracing=-1000000):
	""" *Code*: Finish off the current block of G code (i.e. *self*.) """

	# Use *code* instead of *self*:
	code = self

	# Close *code_stream*:
	code_stream = code._code_stream
	assert isinstance(code_stream, file)
	code_stream.close()
	code._code_stream = None

	# Perform any requested *tracing*:
	if tracing >= 0:
            indent = ' ' * tracing
            print("{0}=>Code._finish(*)".format(indent))

	# Reset the *top_surface_z* and *xy_rapid_safe_z*:
	code._top_surface_z = None
	code._xy_rapid_safe_z = None

	# Random comment: the view3dscene program can view the resulting .wrl file:

	# Write out headers to *vrml_file*:
	vrml_file = code._vrml_file
	part_wrl_file = code._part_wrl_file

	vrml_file.write(    "  #VRML V2.0 utf8\n")
	part_wrl_file.write("  #VRML V2.0 utf8\n")
	vrml_file.write(    "  Shape {\n")
	part_wrl_file.write("  Shape {\n")
	vrml_file.write(    "   geometry IndexedLineSet {\n")
	part_wrl_file.write("   geometry IndexedLineSet {\n")
	vrml_file.write(    "    colorPerVertex FALSE\n")
	part_wrl_file.write("    colorPerVertex FALSE\n")

	# Output the colors:
	vrml_file.write(    "    color Color {\n")
	part_wrl_file.write("    color Color {\n")
	vrml_file.write(    "     color [\n")
	part_wrl_file.write("     color [\n")
	vrml_file.write(    "       0.0 0.0 1.0 # blue\n")
	part_wrl_file.write("       0.0 0.0 1.0 # blue\n")
	vrml_file.write(    "       1.0 0.0 0.0 # red\n")
	part_wrl_file.write("       1.0 0.0 0.0 # red\n")
	vrml_file.write(    "       0.0 1.0 1.0 # cyan\n")
	part_wrl_file.write("       0.0 1.0 1.0 # cyan\n")
	vrml_file.write(    "     ]\n")
	part_wrl_file.write("     ]\n")
	vrml_file.write(    "    }\n")
	part_wrl_file.write("    }\n")

	# Write out points:
	vrml_file.write(    "    coord Coordinate {\n")
	part_wrl_file.write("    coord Coordinate {\n")
	vrml_file.write(    "     point [\n")
	part_wrl_file.write("     point [\n")
	for point in code._vrml_points:
	    vrml_file.write(    "      {0} {1} {2}\n".format(point[0], point[1], point[2]))
	    part_wrl_file.write("      {0} {1} {2}\n".format(point[0], point[1], point[2]))
	vrml_file.write(    "     ]\n")
	part_wrl_file.write("     ]\n")
	vrml_file.write(    "    }\n")
	part_wrl_file.write("    }\n")
		
	# Write out coordinate index:
	vrml_file.write(    "    coordIndex [\n")
	part_wrl_file.write("    coordIndex [\n")
	vrml_point_indices = code._vrml_point_indices
	vrml_file.write(    "    ");
	part_wrl_file.write("    ");
	for index in vrml_point_indices:
	    vrml_file.write(    "   {0}".format(index))
	    part_wrl_file.write("   {0}".format(index))
	    if index < 0:
		vrml_file.write    ("\n  ")
		part_wrl_file.write("\n    ")
	vrml_file.write(    "  ]\n")
	part_wrl_file.write("  ]\n")

	# Write out the color indices:
	vrml_file.write(    "    colorIndex [\n")
	part_wrl_file.write("    colorIndex [\n")
	for color_index in code._vrml_color_indices:
	    vrml_file.write(    "      {0}\n".format(color_index))
	    part_wrl_file.write("      {0}\n".format(color_index))
	vrml_file.write(    "    ]\n")
	part_wrl_file.write("    ]\n")

	# Close out the shape and geometry clauses:
	vrml_file.write(    "   }\n")
	part_wrl_file.write("   }\n")
	vrml_file.write(    "  }\n")
	part_wrl_file.write("  }\n")

	# Close *vrml_file* but leave *part_wrl_file* open:
	vrml_file.write(" ]\n")
	vrml_file.write("}\n")
	vrml_file.close()

	# Now reset all the VRML values:
	code._vrml_reset()

	# Perform any requested *tracing*:
	if tracing >= 0:
            indent = ' ' * tracing
            print("{0}<=Code._finish(*)".format(indent))

    def _g1_set(self, g1):
	""" *Code*: Set the G1 field of the *Code* object (i.e. *self*) to *g1*. """

	# Use *code* instead *self*:
	code = self

	# Verify argument types:
	assert isinstance(g1, int)

	# Set the G1 field of the *Code* object (i.e. *self*) to *g1*:
	code._g1 = g1

    def _hertz(self, field_name, value):
	""" *Code*: Output *value* for *field_name* in the current command of the *Code* object 
	    (i.e. *self*).
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(field_name, str)
	assert isinstance(value, Hertz)

	square_bracket = False
	changed = False
	if field_name == "S":
	    if code._s != value:
		code._s = value
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
	    code._command_chunks.append(chunk)

    def _is_laser_set(self, is_laser):
	""" *Code*: Set the is_laser field of the *Code* object (i.e. *self*) to *is_laser*. """

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
        assert isinstance(is_laser, bool)

	code._is_laser = is_laser

    def _length(self, field_name, value):
	""" *Code*: Output *value* for *field_name* to the *Code* object (i.e. *self*.) """

	# Use *code* instead of *self*:
        code = self

	# Verify argument tyeps:
	assert isinstance(field_name, str)
	assert isinstance(value, L)

	square_brackets = False
	offset = L()
	changed = False
	if field_name == "Q":
	    if code._q != value:
		code._q = value
		changed = True
	elif field_name == "I":
	    if code._i != value:
		code._i = value
		changed = True
	elif field_name == "J":
	    if code._j != value:
		code._j = value
		changed = True
	elif field_name == "R0":
	    if code._r0 != value:
		code._r0 = value
		changed = True
 	elif field_name == "R1":
	    if code._r1 != value:
		coe._r1 = value
		changed = True
	elif field_name == "X":
	    #offset = code._vice_x
	    offset = L()
	    if code._x != value - offset:
		code._x = value - offset
		changed = True
	elif field_name == "Y":
	    #offset = code._vice_y
	    offset = L()
	    if code._y != value - offset:
		code._y = value - offset
		changed = True
	elif field_name == "Z":
	    if code._z != value:
		code._z = value
		changed = True
	elif field_name == "Z1":
	    if code._z1 != value:
		code._z1 = value
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
	    code._command_chunks.append(chunk)

    def _line_comment(self, comment):
	""" *Code*: Emit *comment* to the *Code* object (i.e. *self*). """

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(comment, str)

	# This routine will add a line containg *comment* to the *Code* object (i.e. *self).:
	code._command_begin()
	code._comment(comment)
	code._command_end()

    def _mode_motion(self, g1_value, motion_color_index):
	""" *Code*: This routine will issue a G*g1_value* motion command to the *Code* object
	    (i.e. *self*).  *motion_color_index* specifies what color index to draw into the tool
	    path .wrl file.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify arguement types:
	assert isinstance(g1_value, int) and g1_value in code._g1_table
	assert isinstance(motion_color_index, int)

	# Add a G1 field to the command:
	code._unsigned("G1", g1_value)
	code._vrml_motion_color = motion_color_index

    def _reset(self):
	""" *Code*: Reset the contents of the *Code* object (i.e. *self*) """

	# Use *code* instead of *self:

	zero = L()
	large = L(inch=123456789.0)
	huge = 0x7fffffff
	big = L(inch=123456789)

	# Reset the *code* object:
	code._begin = True
	code._f = Speed(mm_per_sec=huge)
	code._g1 = huge
	code._g2 = huge
	code._g3 = huge
	code._g4 = huge
	code._g5 = huge
	code._g6 = huge
	code._g7 = huge
	code._g8 = huge
	code._g9 = huge
	code._g10 = huge
	code._g11 = huge
	code._h = huge
	code._i = big
	code._j = big
	code._m1 = huge
	code._m2 = huge
	code._m3 = huge
	code._m4 = huge
	code._m5 = huge
	code._p = Time(sec=-1.0)
	code._q = big
	code._r0 = big
	code._r1 = big
	code._s = Hertz(rps=123456789.0)
	code._x = big
	code._y = big
	code._z = big
	code._z1 = big

    def _start(self,
      part, tool, ngc_program_number, spindle_speed, part_wrl_file, stl_file_name,tracing=-1000000):
	""" *Code*: Start writing out the G-code for *tool* (
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(tool, Tool)
	assert isinstance(ngc_program_number, int) and ngc_program_number > 0
	assert isinstance(spindle_speed, Hertz)
	assert isinstance(part_wrl_file, file)

	# Grab some values from *part* and *tool*:
	part_name = part._name_get()
	tool_name = tool._name_get()
	tool_number = tool._number_get()

	# Do any requested *tracing*:
	tracing_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._start('{1}', '{2}', {3}, {4:rpm}, *)".
	      format(indent, part_name, tool_name, ngc_program_number, spindle_speed))
	    tracing_detail = 1

	# Make sure that closed off any previous *code_stream*:
	assert code._code_stream == None

	# Open new *code_stream*:
	ezcad = part._ezcad_get()
	ngc_directory = ezcad._ngc_directory_get()
	code_file_name = os.path.join(ngc_directory, "O{0}.ngc".format(ngc_program_number))
	code_stream = open(code_file_name, "w")
	assert code_stream != None, "Could not open '{0}' for writing".format(code_file_name)
	code._code_stream = code_stream
	if tracing_detail >= 1:
            print("{0}code_file_name='{1}' opened".format(indent, code_file_name))

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

	# Open new *vrml_file*:
	ngc_directory = ezcad._ngc_directory_get()
	vrml_file_name = os.path.join(ngc_directory, "O{0}.wrl".format(ngc_program_number))
	vrml_file = open(vrml_file_name, "w")
	assert vrml_file != None, "Could not open '{0}' for writing".format(vrml_file_name)

	# Start the top-level group in *vrml_file*:
	vrml_file.write("#VRML V2.0 utf8\n")
        vrml_file.write("Group {\n")
	vrml_file.write(" children [\n")

	# Write out the part visualization:
	part._wrl_write(vrml_file,
	  part._cnc_transform, 2, stl_file_name, parts_table={}, tracing=tracing + 1)

	code._vrml_file = vrml_file
	code._part_wrl_file = part_wrl_file

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Code._start('{1}', '{2}', {3}, {4:rpm}, *)".
	      format(indent, part_name, tool_name, ngc_program_number, spindle_speed))


    def _top_surface_z_set(self, top_surface_z):
        """ *Code*: Set the *top_surface_z* into the *Code* object (i.e. *self*.)
	    Veritical rapid moves are allowed above this distance.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(top_surface_z, L)

	# Save *top_surface_z* into *code*:
	code._top_surface_z = top_surface_z

    #FIXME: This is no longer used!!!
    def _foo(self):
	# Use *code instead of *self*:
        code = self

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

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(s, Hertz)

	# Set the S field of the *Code* object (i.e. *self*) to *s*:
	code._s = s

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
    
	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(pocket, Operation_Simple_Pocket)
	assert isinstance(offset, L)
	assert isinstance(s, Hertz)
	assert isinstance(f, Speed)
	assert isinstance(z, L)
	assert isinstance(rapid_move, bool)

	# Extract the corners:
	corner1 = pocket._corner1
	corner2 = pocket._corner2
	cx1 = corner1.x
	cy1 = corner1.y
	cz1 = corner1.z
	cx2 = corner2.x
	cy2 = corner2.y
	cz2 = corner2.z

	# Now order the arguments:
	x1 = min(cx1, cx2)
	y1 = min(cy1, cy2)
	z1 = min(cz1, cz2)
	x2 = max(cx1, cx2)
	y2 = max(cy1, cy2)
	z2 = max(cz1, cz2)

	corner_radius = pocket._corner_radius_get()
	code._line_comment(
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
	code._line_comment("rx1={0:i} ry1={1:i} rx2={2:i} ry2={3:i}".format(rx1, ry1, rx2, ry2))
	if debug:
	    print("rx1={0:i} ry1={1:i} rx2={2:i} ry2={3:i}".format(rx1, ry1, rx2, ry2))

	# Compute the rectangle path coordinates:
	px1 = x1 + offset
	px2 = x2 - offset
	py1 = y1 + offset
	py2 = y2 - offset
	code._line_comment(
	  "offset={0:i} px1={1:i} py1={2:i} px2={3:i} py2={4:i}".format(offset, px1, py1, px2,py2))
    
	# Determine the starting location for this path:	
	start_x = px1
	start_y = py1
	if offset < corner_radius:
	    # We have rounded corner path:
	    start_x = rx1
    
	# Move to (*start_x*, *start_y*) as specified by *linear_move* argument:
	if rapid_move:
	    code._xy_rapid(start_x, start_y)
	else:
	    code._xy_feed(f, s, start_x, start_y)
    
	# Make sure we are at the depth *z*:
	code._z_feed(f/2, s, z, "simple_pocket_helper")
    
	# Mill out either a square or rounded corners
	if offset < corner_radius:
	    # Mill out a rectangle with rounded corners in a
	    # counter clockwise direction to force a climb cut:
    
	    r = corner_radius - offset
    
	    # Bottom horizontal line from (rx1,py1) to (rx2,py1):
	    code._xy_feed(f, s, rx2, py1)
    
	    # Lower right arc (rx2,py1) to (px2,ry1):
	    code._xy_ccw_feed(f, s, r, px2, ry1, rx=rx2, ry=ry1)
    
	    # Right vertical line (px2,ry1) to (px2,ry2):
	    code._xy_feed(f, s, px2, ry2)
    
	    # Upper right arc (px2,ry2) to (rx2, py2):
	    code._xy_ccw_feed(f, s, r, rx2, py2, rx=rx2, ry=ry2)
    
	    # Top horizontal line (rx2, py2) to (rx1, py2):
	    code._xy_feed(f, s, rx1, py2)
    
	    # Upper left arc (rx1, py2) to (px1, ry2):
	    code._xy_ccw_feed(f, s, r, px1, ry2, rx=rx1, ry=ry2)
    
	    # Left vertical line (px1, ry2) to (px1, ry1):
	    code._xy_feed(f, s, px1, ry1)
    
	    # Lower left arc (px1, ry1) to (rx1, py1):
	    code._xy_ccw_feed(f, s, r, rx1, py1, rx=rx1, ry=ry1)
	else:
	    # Mill out a rectangle with "square" corners in a counter
	    # clockwise direction to force a climb cut:
    
	    # Bottom horizontal line from (px1, py1) to (px2, py1):
	    code._xy_feed(f, s, px2, py1)
    
	    # Right vertical line from (px2, py1) to (px2, py2):
	    code._xy_feed(f, s, px2, py2)
    
	    # Top horizontal line from (px2, py2) to (px1, py2):
	    code._xy_feed(f, s, px1, py2)
    
	    # Left vertical line from (px1, py2) to (px1, py1):
	    code._xy_feed(f, s, px1, py1)

    def _speed(self, field_name, value):
	""" *Code*: Set the speed for *field_name* to *value in the the current command of the
	    *Code* object (i.e. *self*).
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(field_name, str)
	assert isinstance(value, Speed)

	square_bracket = False
	changed = False
	if field_name == "F":
	    assert value != Speed()
	    if code._f != value:
		code._f = value
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
	    code._command_chunks.append(chunk)

    #FIXME: This is old code!!!
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
	dowel_x = part._dowel_x_get()
	dowel_y = part._dowel_y_get()
	vice_x = part._vice_x_get()
	vice_y = part._vice_y_get()
	assert isinstance(vice_x, L)
	plunge_x = vice_x
	if plunge_x > dowel_x:
	    plunge_x = dowel_x
	assert isinstance(jaw_width, L)
	assert isinstance(plunge_x, L)
	if plunge_x > jaw_width:
	    plunge_x = plunge_x - L(inch=0.7)
	part._plunge_xy_set(plunge_x, dowel_y)
    
	code = shop._code_get()
	assert isinstance(code, Code)
	code._z_safe_set(part._z_safe_get())
	code._z_rapid_set(part._z_rapid_get(), part)
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
	ngc_directory = ezcad._ngc_directory_get()
	part_ngc_file_name = os.path.join(ngc_directory, "O{0}.ngc".format(program_number))
	part_ngc_stream = open(part_ngc_file_name, "w")
	assert part_ngc_stream != None, "Unable to open {0}".format(part_ngc_file_name)

    def _unsigned(self, field_name, value):
	""" *Code*: This routine will format {value} for {field_name} to {code}.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(field_name, str)
	assert isinstance(value, int) and value >= 0

	#FIXME: This should probably be done using a Python dictionary:

	previous_value = 0xffffffff
	matched = False
	size = len(field_name)
	if size == 1:
	    matched = True
	    letter = field_name[0]
	    if letter == 'H':
		previous_value = code._h
		if previous_value != value:
		    code._h = value
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
		    previous_value = code._g1
		    if previous_value != value:
			code._g1 = value
		elif g_digit == '2':
		    previous_value = code._g2
		    if previous_value != value:
			code._g2 = value
		elif g_digit == '3':
		    previous_value = code._g3
		    if previous_value != value:
			code._g3 = value
		elif g_digit == '4':
		    previous_value = code._g4
		    if previous_value != value:
			code._g4 = value
		elif g_digit == '5':
		    previous_value = code._g5
		    if previous_value != value:
			code._g5 = value
		elif g_digit == '6':
		    previous_value = code._g6
		    if previous_value != value:
			code._g6 = value
		elif g_digit == '7':
		    previous_value = code._g7
		    if previous_value != value:
			code._g7 = value
		elif g_digit == '8':
		    previous_value = code._g8
		    if previous_value != value:
			code._g8 = value
		elif g_digit == '9':
		    previous_value = code._g9
		    if previous_value != value:
			code._g9 = value
		else:
		    # We did not match:
		    matched = False
	elif size == 2:
	    if field_name == "G10":
		matched = True
		previous_value = code.g10
		if previous_value != value:
		    code._g10 = value
	    elif field_name == "G11":
		# Always output G4 when requested:

		matched = True
		previous_value = 0x12345

	assert matched, "Unrecognized field name {0}".format(field_name)
	if previous_value != value or field_name == "G1":
	    code._command_chunks.append("{0}{1}".format(field_name[0], value))

    def old__vice_xy_set(self, vice_x, vice_y):
	""" *Code*: Set the vice X/Y fields of the *Code* object (i.e. *self*)
	    to *vice_x* and *vice_y*.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(vice_x, L)
	assert isinstance(vice_y, L)

	# Set the vice X/Y fields of the *Code* object (i.e. *self*) to *vice_x* and *vice_y*:
	code._vice_x = vice_x
	code._vice_y = vice_y

    def _vrml_arc_draw(self, ax, ay, bx, by, radius, z, clockwise, radius_x=None, radius_y=None):
	""" *Code*: Draw an arc from (*ax*, *ay*, *z) to (*bx*, *by*, *bz*) with an
	    arc radius of *radius*.  The arc is drawn in a clockwise direction if
	    *clockwise* is *True* and a counter-clockwise direction otherwise. """
    
	# Use *code* instead of *self*:
	code = self

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
	#	       R
	#	      /|\
	#	     / | \
	#	    /  |  \
	#	   /   |   \
	#	  /    |--+ \
	#	 /     |  |  \
	#	A------C------B
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
	#	|AC|^2 + |CR|^2 =	 |AR|^2				(1) 
	#	|AC|^2 + |CR|^2 =       radius^2			(2) 
	#		 |CR|^2 =       radius^2 - |AC|^2		(3) 
	#		 |CR|   = sqrt( radius^2 - |AC|^2 )		(4)
	#		 |CR|   = sqrt( radius^2 - |AC|^2 )		(5)
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
	vrml_point_indices = code._vrml_point_indices
	a_index = code._vrml_point(ax, ay, az)
	vrml_point_indices.append(a_index)

	# Lay down the intermediate points along the arc:
	for step in range(1, steps - 1):
	    segment_angle = arx_angle + step_angle * step
	    segment_x = rx + radius.cosine(segment_angle)
	    segment_y = ry +   radius.sine(segment_angle)
	    segment_z = az
	    segment_index = code._vrml_point(segment_x, segment_y, segment_z)
	    vrml_point_indices.append(segment_index)

	# Lay down B=(*bx*,*by*,*bz):
	b_index = code._vrml_point(bx, by, bz)
	vrml_point_indices.append(b_index)

	# Terminate the polyline:
	vrml_point_indices.append(-1)

	# Set the color for the polyline:
	code._vrml_color_indices.append(code._vrml_motion_color_cutting)

    def _vrml_point(self, x, y, z):
	""" *Code*: Return the VRML point index for point (*x*, *y*, *z*) using the
	    *Code* object (i.e. *self*).
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(z, L)

	point = (x.millimeters(), y.millimeters(), z.millimeters())
	vrml_points_table = code._vrml_points_table
	if point in vrml_points_table:
	    index = vrml_points_table[index]
	else:
	    vrml_points = code._vrml_points
	    index = len(vrml_points)
	    vrml_points.append(point)
	return index

    def _vrml_reset(self):
	""" *Code*: Reset the VRML sub-system of the *Code* object (i.e. *self*). """

	# Use *code* instead of *self*:
	code = self

	# Reset the VRML data structures of *code*:
	code._vrml_color_indices = []
	code._vrml_file = None
	code._vrml_motion_color = -1
	code._vrml_point_indices = []
	code._vrml_points = []
	code._vrml_points_table = []
	code._vrml_start_r0 = code._r0
	code._vrml_start_x = code._x
	code._vrml_start_y = code._y
	code._vrml_start_z = code._z

    def _x_value(self):
	""" *Code*: Return the current value of X offset by the vice X for the *Code* object
	    (i.e. *self*).
	"""

	# Use *code* instead of *self*:
	code = self

	return code._x # + code._vice_x

    def _xy_cw_feed(self, f, s, r, x, y, rx=None, ry=None):
	""" *Code*: Feed to location (*x*, *y*) with a radius *r* clockwise  circle with a
	    feedrate of *f* and spindle speed of *s* using the *Code* object (i.e. *self*):
	"""
    
	# Use *code* instead of *self*:
        code = self

	# Verify routine arguments:
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)
	assert isinstance(r, L)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(rx, L) or rx == None
	assert isinstance(ry, L) or ry == None
    
	x1 = code._x_value()
	y1 = code._y_value()
	z1 = code._z
	x2 = x
	y2 = y
	code._vrml_arc_draw(x1, y1, x2, y2, r, z1, False, rx, ry)

	x_value = code._x_value()
	y_value = code._y_value()
	if x_value != x or y_value != y:
	    # Do the laser code first:
	    if code._is_laser:
	    	code._dxf_arc_append(True, x, y, r)

	    # Now do the RS274 code and get F, R0, S, X, and Y updated:
	    code._z_safe_retract_actual()
	    code._r0 = L(inch=-100.0)	# Forget R
	    code._command_begin()
	    code._mode_motion(2, code._vrml_motion_color_cutting)
	    code._speed("F", f)
	    code._length("R0", r)
	    code._hertz("S", s)
	    code._length("X", x)
	    code._length("Y", y)
	    code._command_end()
    
    def _xy_ccw_feed(self, f, s, r, x, y, rx=None, ry=None):
	""" *Code*: Feed to location (*x*, *y*) as a radius *r* counter clockwise  circle with a
	    feedrate of *f* and spindle speed of *s* using the *Code* object (i.e. *self*):
	"""
    
	# Use *code* instead of *self*:
        code = self

	# Verify routine arguments:
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)
	assert isinstance(r, L)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(rx, L) or rx == None
	assert isinstance(ry, L) or ry == None
    
	x1 = code._x_value()
	y1 = code._y_value()
	z1 = code._z
	x2 = x
	y2 = y
	code._vrml_arc_draw(x1, y1, x2, y2, r, z1, True, rx, ry)

	if x1 != x or y1 != y:
	    # Do the laser code first:
	    if code._is_laser:
	    	code._dxf_arc_append(False, x, y, r)
    
	    # Now do the RS274 code and get F, R0, S, X, and Y updated:
	    code._z_safe_retract_actual()
	    code._r0 = L(inch=-100.00)	# Forget R
	    code._command_begin()
	    code._mode_motion(3, code._vrml_motion_color_cutting)
	    code._speed("F", f)
	    code._length("R0", r)
	    code._hertz("S", s)
	    code._length("X", x)
	    code._length("Y", y)
	    code._command_end()

    def _xy_feed(self, f, s, x, y):
	""" *Code*: Feed to location (*x*, *y*) with a feedrate of *f*
	    and spindle speed of *s* using the *Code* object (i.e. *self*).
	"""
    
	# Use *code* instead of *self*:
        code = self

	# Verify argment types:
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)
	assert isinstance(x, L)
	assert isinstance(y, L)

	x_before = code._x_value()
	y_before = code._y_value()
	if x_before != x or y_before != y:
	    # Do any laser code first:
	    if code._is_laser:
	    	code._dxf_entity_start("LINE")
	    	code._dxf_xy_append(0, x_before, y_before, "xy_feed before")
	    	code._dxf_xy_append(1, x, y, "xy_feed after")
	    	code._dxf_entity_stop()
    
	    # Now do the RS274 code and get the F, S, X, and Y values updated:
	    code._command_begin()
	    code._mode_motion(1, code._vrml_motion_color_cutting)
	    code._speed("F", f)
	    code._hertz("S", s)
	    code._length("X", x)
	    code._length("Y", y)
	    code._command_end()

    def _xy_rapid(self, x, y):
	""" *Code*: Perform a rapid move to (X, Y) using the *Code* object (i.e. *self*). """

	# Use *code* instead of *self:
	code = self

	# Verify argument types:
	assert isinstance(x, L)
	assert isinstance(y, L)

	# Only perform the rapid if we are not already there:
	if code._x_value() != x or code._y_value() != y:
	    # Make sure we are at a safe Z height prior to performing the rapid:
	    assert code._z == code._xy_rapid_safe_z

	    # Output a X/Y rapid command:
	    code._command_begin()
	    code._mode_motion(0, code._vrml_motion_color_rapid)
	    code._length("X", x)
	    code._length("Y", y)
	    code._command_end()

    def _xy_rapid_safe_z_set(self, xy_rapid_safe_z):
        """ *Code*: Save *xy_rapid_safe_z* into the *Code* object (i.e. *self*).
	    *xy_rapid_safe_z* is the Z location above which rapid moves can occur in X and Y.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
        assert isinstance(xy_rapid_safe_z, L)

	# Save *xy_rapid_save_z* into *code*:
	code._xy_rapid_safe_z = xy_rapid_safe_z

    def _xy_rapid_safe_z_force(self, feed, speed):
        """ *Code*: Force the tool for the *Code* object (i.e. *self*) to be at the
	    safe height where rapids in X and Y are allowed.
	"""

	# Use *code* instead of *self*:
        code = self

	xy_rapid_safe_z = code._xy_rapid_safe_z
	code._z_feed(feed, speed, code._xy_rapid_safe_z, "xy_rapid_safe_z_force")

    def _y_value(self):
	""" *Code*: Return the current value of Y offset by the vice Y for the *Code* object
	    (i.e. *self*).
	"""

	# Use *code* instead of *self*:
	code = self

	return code._y # + code._vice_y

    def _z_feed(self, f, s, z, from_routine):
	""" *Code*: Feed to an altitude of *z* using the *Code* object (i.e. *self*)
	    at a feed of *f* and a speed *s*.
	"""

	#print("=>Code._z_feed(*, f={0:i}, s={1:rps}, z={2:i}, '{3}')".
        #  format(f, s, z, from_routine))

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)
	assert isinstance(z, L)
	assert isinstance(from_routine, str)

	# This routine will feed to an altitude of {z} using {code}.

	# If *z_safe_pending* is set, but we are doing another Z feed before
	# it is cleared, we must not need to do a Z safe move operation after all.
	# Hence, we can just clear *z_safe_pending*:
	code._z_safe_pending = False

	top_surface_z = code._top_surface_z
	xy_rapid_safe_z = code._xy_rapid_safe_z
	#print("top_surface_z={0:i} xy_rapid_safe_z={1:i}".format(top_surface_z, xy_rapid_safe_z))

	# Emit G0/G1 commands until we get to *z* (i.e. *z_target*):
	z_current = code._z
	while z_current != z:
	    z_target = z
	    #print("---- z_current={0:i} z_target={1:i}".format(z_current, z_target))

	    if top_surface_z <= z_current <= xy_rapid_safe_z:
		# We are in the retraction portion of the workspace.  G0 rapids
		# are used to move around here:
		
		# Trim *z_target* to be between *top_surface_z* and *xy_rapid_safe_z*:
                if z_target <= top_surface_z:
                    z_target = top_surface_z
		elif z_target >= xy_rapid_safe_z:
                    z_target = xy_rapid_safe_z
		#print("GO trimed z_target={0:i}".format(z_target))

		# Perform a G0 rapid (if necessary):
		if z_current != z_target:
		    #print("G0 performed")
		    code._command_begin()
		    code._mode_motion(0, code._vrml_motion_color_retract)
		    code._length("Z", z_target)
		    code._command_end()
		    z_current = z_target

	    z_target = z
	    if z_current <= top_surface_z:
                # We are in the cutting portion of the workspace.  G1 motion is
                # used to move around here:

		# Make sure *z_target* does not get above *top_surface_z*:
		if z_target >= top_surface_z:
		    z_target = top_surface_z
		#print("G1 trimed z_target={0:i}".format(z_target))

		# Perform G1 move (if necessary):
		if z_current != z_target:
		    #print("G1 performed")
		    code._command_begin()
		    code._mode_motion(1, code._vrml_motion_color_cutting)
		    code._speed("F", f)
		    code._hertz("S", s)
		    code._length("Z", z_target)
		    code._command_end()
		    z_current = z_target

	#print("<=Code._z_feed(*, f={0:i}, s={1:rps}, z={2:i}, '{3}')".
        #  format(f, s, z, from_routine))

    def _xxz_rapid_set(self, z_rapid, part):
	""" *Code*: Set the Z rapid field of the *Code* object (i.e. *self*) to *z_rapid*. """

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(z_rapid, L) or z_rapid == None
	assert isinstance(part, Part)

	if z_rapid == None:
	    print("Part '{0}' does not set z_rapid".format(part._name_get()))
            z_rapid = L(0.5)

	# Set the Z field of the *Code* object (i.e. *self*) to *z*:
	code._z_rapid = z_rapid
    

    def _xxz_safe_assert(self, text, comment):
	""" *Code*: Fail if the Code object (i.e. *self*)  is not at "Z safe".
	    If not, the failing message will contain {text}.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument_types:
	assert isinstance(text, str)
	assert isinstance(comment, str)

	assert code._z_safe_pending or code._z == code._z_safe, \
	  "Z safe failure: {0} ({1}) [z={2}]".format(text, comment, code._z)

    def _xz_safe_set(self, z_safe):
	""" *Code*: Set the Z safe field of the *Code* object (i.e. *self*). """

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	#assert isinstance(z_safe, L)

	# Set the Z safe field of the *Code* object (i.e. *self*) to *z_safe*:
	code._z_safe = z_safe

    def _xz_safe_retract(self, f, s):
	""" *Code*: Ensure that the tool is at the "Z safe" altitude for next command RS274
	    command in the *Code* object (i.e. *self) and get there using a feedrate of *f*
	    and a speed of *s* if needed.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)

	# Mark that a "Z safe" operation is pending if needed:
	if code._z != code._z_safe:
	    code._z_safe_pending = True
	    code._z_safe_s = s
	    code._z_safe_f = f

    def _xz_safe_retract_actual(self):
	""" *Code*  Ensure that the tool is at the "Z safe" altitude in the *Code* object
	    (i.e. *self*) and get there using the internal F and S fields.
	"""

	# Use *code* instead of *self*:
	code = self

	# If we have a pending "z_safe" scheduled, this code will make it happen:
	if code._z_safe_pending:
	    code._z_safe_pending = False
	    code._z_feed(code._z_safe_f, code._z_safe_s, code._z_safe, "z_safe_retract_actual")

	#assert code._z == code._z_safe, "Did not get to Z safe"

    def _z_set(self, z):
	""" *Code*: Set the Z field of the *Code* object (i.e. *self*). """

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(z, L)

	# Set the Z field of the *Code* object to *z*:
	code._z = z

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
	  "}, axis={3"			       + format + \
	  "}, rotate={4"			     + format + \
	  "}, translate={5"			  + format + \
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

	# Use *shop* instead of *self*:
	shop = self

	# Verify argument types:
	assert isinstance(name, str)

	# Define some useful values:
	zero = L()
	tools = []
	vice_volume = P(L(inch=5.0), L(inch=6.0), L(inch=1.0))

	# Tooling plate with #4 holes:
	no4_close_fit_diameter = L(inch=.1160)
	no4_75_percent_thread_diameter = L(inch=.0890)
	no4_50_percent_thread_diameter = L(inch=.0960)
	tooling_plate = Tooling_Plate(L(inch=6.0), L(inch=6.0), L(inch=.250), 11, 11, L(inch=.500),
	  no4_close_fit_diameter, no4_75_percent_thread_diameter, no4_50_percent_thread_diameter)

	# Vice to use:
	parallels = Parallels(L(inch=6.0), L(inch="1/8"), [
	  L(inch="1/2"),
	  L(inch="5/8"),
	  L(inch="3/4"),
	  L(inch="7/8"),
	  L(inch=1.000),
	  L(inch="1-1/8"),
	  L(inch="1-1/4"),
	  L(inch="1/2"),
	  L(inch="1-1/2") ] )
	vice = Vice(vice_volume, parallels)

	# Initialize *shop*:
	shop._assemblies = []		# Viewable assemblies
	shop._base_name = ""		# Base name for current operations
	#dxf_base_names = [] 		# List of DXF base names
	shop._blocks_uid = 0		# UID counter for *Code_Block* objects
	# cache Cache			# Off/Nef3 file cache
	shop._changed = False		# Marker used by {update@Length}
	shop._cnc_generate = False	# {true} => do cnc code generation
	shop._code = Code()		# Code genertion object
	#dxf_table = {}			# Table of DXF base names
	shop._extra1 = P()		# Temporary extra bounding box point
	shop._extra2 = P()		# Temporary extra bounding box point
	shop._name = name		# Shop name
	shop._machines = []		# Machines in shop
	shop._matrix = Matrix()		# Temporary matrix
	shop._parts = []			# Parts
	shop._program_base = 10		# Program base number
	shop._solids_generate = False	# {true} => do solid modeling
	shop._surface_speeds_table = {}	# [Part_Material, Tool_Material]
	shop._tools = tools		# Tools in shop
	shop._tooling_plate = tooling_plate
	shop._vice = vice
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

	# Start loading some surface speeds into *shop*:
	hss = Tool.MATERIAL_HIGH_SPEED_STEEL
	plastic = Material("Plastic", "")
	aluminum = Material("Aluminum", "")
	fpm600 = Speed(ft_per_min=600)
	fpm1200 = Speed(ft_per_min=1200)
	shop._surface_speeds_insert(aluminum, hss, fpm600, fpm1200)
	shop._surface_speeds_insert(aluminum, hss, fpm600, fpm1200)
	shop._surface_speeds_insert(plastic, hss, fpm600, fpm1200)
	shop._surface_speeds_insert(plastic, hss, fpm600, fpm1200)

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

	dowel_pin = shop._dowel_pin_append("3/8 Dowel Pin",
	  1, hss, in3_8, L(inch=.900), in3_16)
	mill_drill_3_8 = shop._mill_drill_append("3/8 Mill Drill",
	  2, hss, in3_8, 2, L(inch=.900), in3_16, degrees90)
	drill_36 = shop._drill_append("#36 Drill",
	  3, hss, L(inch=0.1065), 2, L(inch=1.500), degrees118, stub)
	drill_27 = shop._drill_append("#27 drill",
	  4, hss, L(inch=0.1440), 2, L(inch=1.750), degrees118, stub)
	end_mill_3_8 = shop._end_mill_append("3/8 End Mill",
	  5, hss, in3_8, 2, in5_8, not laser)
	end_mill_1_4 = shop._end_mill_append("1/4 End Mill",
	  6, hss, in1_4, 2, in1_2, not laser)
	double_angle = shop._double_angle_append("3/4 Double Angle",
	  7, hss, in3_4, 10, L(inch=0.875), degrees90, in1_4, in1_4)
	dove_tail = shop._dove_tail_append("3/8 Dove Tail",
	  8, hss, in3_8, 6, in1_4, in3_16, degrees45)
	end_mill_3_16 = shop._end_mill_append("3/16 End Mill",
	  10, hss, in3_16, 2, in1_2, not laser)
	drill_25 = shop._drill_append("#25 drill",
	  11, hss, L(inch=0.1495), 2, L(inch=2.000), degrees118, stub)
	drill_9 = shop._drill_append("#9 drill",
	  12, hss, L(inch=0.1960), 2, L(inch=2.000), degrees118, stub)
	drill_43 = shop._drill_append("#43 drill",
	  13, hss, L(inch=.0890), 2, L(inch=1.500), degrees118, stub)
	drill_32 = shop._drill_append("#32 drill",
	  14, hss, L(inch=0.1160), 2, L(inch=1.500), degrees118, stub)
	drill_50 = shop._drill_append("#50 drill",
	  15, hss, L(inch=0.0700), 2, L(inch=1.500), degrees118, stub)
	end_mill_3_8_long = shop._end_mill_append("3/8 1\" End Mill",
	  16, hss, in3_8, 2, L(inch=1.000), not laser)
	#end_mill_3_4 = shop._end_mill_append("3/4 End Mill",
	#  13, hss, in3_4, 2, in1_3_8)
	drill_30 = shop._drill_append("#30 drill",
	  17, hss, L(inch=0.1285), 2, L(inch=1.750), degrees118, stub)
	drill_1_8 = shop._drill_append("1/8 drill",
	  18, hss, in1_8, 2, L(inch=1.750), degrees118, stub)
	drill_1_16 = shop._drill_append("1/16 drill",
	  19, hss, in1_16, 2, L(inch=1.750), degrees118, stub)
	drill_3_32 = shop._drill_append("3/32 drill",
	  20, hss, in3_32, 2, L(inch=1.750), degrees118, stub)
	drill_42 = shop._drill_append("#42 drill",
	  21, hss, L(inch=0.0935), 2, L(inch=1.750), degrees118, stub)

	# Laser "tools":
	#laser_007 = shop._end_mill_append("Laser_007",
	#  100, hss, L(inch=0.007), 2, L(inch=0.750), laser)
	#laser_000 = shop._end_mill_append("Laser_000",
	#  101, hss, L(), 2, L(inch=0.750), laser)

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
	tool_dowel_pin._feed_speed_set(Speed(in_per_sec = 1.000))
	self._tool_append(tool_dowel_pin)
	#print("tool_dowel_pin=", tool_dowel_pin)
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

    def _program_base_set(self, program_base):
	""" *Shop*: Set program base number associated with the *Shop* object
	    (i.e. *self*) to *program_base*.
	"""

	self._program_base = program_base

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
	#debug = True
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
	new_name = new_tool._name_get()
	tools = self._tools
	for tool in tools:
	    number = tool._number_get()
	    assert new_number != number, \
	      "Tool number {0} is conflicts between '{1}' and '{2}'".format(
	     number, new_tool._name_get(), tool._name_get())
	    name = tool._name_get()
	    assert new_name != name, \
	      "Tool name '{0}' has already been added (tool number {1}".format(name, number)

	# Append *new_tool* to *tools*:
	tools.append(new_tool)

    def _tooling_plate_get(self):
        """ *Shop*: Returns the currently selected tooling plate for the *Shop*object 
	    (i.e. *self*.)
	"""

	return self._tooling_plate

class Time:
    def __init__(self, sec=0.0, min=0.0):
	self._seconds = sec + min * 60.0 

class Tool:
    """ *Tool*: A tool is a tool that can be used to manufacture a part. """

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

    def __format__(self, format):
	""" *Tool*: Return the *Tool* object (i.e. *self*) formatted as a string. """

	return self._name

    def _compare(self, tool2):
	""" *Tool*:  This routine will return -1 or 1 depending upon whether *tool1*
 	    should be used before or after *tool2*.  0 is returned if *tool1*
	    and *tool2* are equal.
	"""

	# Verify argument_types:
	assert isinstance(tool2, Tool)

	# Use *tool1* instead of *self*:
	tool1 = self

	# Do the comparisons:
	result = int_compare(tool1._kind, tool2._kind)
	if result == 0:
	    result = tool1._diameter.compare(tool2._diameter)
	    if result == 0:
		result = int_compare(tool1._number, tool2._number)
	return result

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

	assert feed_speed != Speed()
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
	    *number*, *material*, *diameter*, *maximum_z_depth*, and *tip_depth*.  *material*
	    is the tool material (e.g. HSS, carbide, ...), *diameter* is the dowel pin diameter,
	    *maximum_z_depth* is the maximum distance the entire dowel pin can go down below
	    the part top surface, and *tip_depth* is the amount of the dowel pin at the tip
	    that is not *diameter* round.  (Dowel pins can have a flat, round or conical tip.)
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

	self._feed_speed_set(Speed(in_per_sec=1.0))

	# Load up *self*:
	self._tip_depth = tip_depth		# Tip portion that is not vertical

    @staticmethod
    def _match(tool, parameter1, parameter2, from_string, tracing):
	""" *Tool_Dowel_Pin*: Return priority of match *tool* matching a
	    dowel pin.
	"""
	
	# Verify argument types:
	assert isinstance(tool, Tool)
	assert isinstance(parameter1, L)
	assert isinstance(parameter2, L)
	assert isinstance(from_string, str)
	assert isinstance(tracing, int)

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
	""" *Tool_Dowel_Pin: Return the tip depth of the *Tool_Dowel_Pin* object (i.e. *self*.)
	"""

	return self._tip_depth

class Tool_Drill(Tool):

    def __init__(self,
      name, number, material, diameter, flutes_count, maximum_z_depth, point_angle, drill_style):
	""" *Tool_Drill*: Initialize *Tool_Drill* object (i.e. *self*) with *name*, *number*,
	    *material*, *diameter*, *flutes_count*, *maximum_z_depth*, *point_angle*, and
	    *drill_style*.
	"""

	# 

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
	Tool.__init__(self, name, number, Tool.KIND_DRILL,
	  material, diameter, flutes_count, maximum_z_depth)
	self._point_angle = point_angle
	self._drill_style = drill_style

    def _drill_style_get(self):
	""" *Tool_Drill*: Return the drill style for the *Tool_Drill* object (i.e. *self*). """

	return self._drill_style

    @staticmethod
    def _match(tool, desired_diameter, maximum_z_depth, from_routine, tracing):
	""" *Tool_Drill*:  Return the diameter of *tool*, provided it is a *Drill_Tool*
	    object with a diameter that is very close to *desired_diameter* and
	    can reach at least to *maximum_z_depth*.  If no match occurs, -1.0 is returned.
	    *from_routine* is used for tracing.
	"""
    
	# Verify argument types:
	assert isinstance(tool, Tool)
	assert isinstance(desired_diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(from_routine, str)
	assert isinstance(tracing, int)

	# Perform any requested tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Tool_Drill._match('{1}', {2:i}, {3:i}, '{4}')".
	      format(indent, tool, desired_diameter, maximum_z_depth, from_routine))

	# Priority is negative if *tool* is not a *Tool_Drill* object:
	priority = -1.0
	if isinstance(tool, Tool_Drill):
	    # Make sure *tool* is long enough:
	    if tool._maximum_z_depth >= maximum_z_depth:
		# If the *tool* *diameter* is very close to *desired_diameter* we have match::
		diameter = tool._diameter
		epsilon = L(inch=0.00001)
		if (desired_diameter - diameter).absolute() < epsilon:
		    # We have a match; set *priority* to the diameter:
		    priority = diameter.millimeters()
	assert isinstance(priority, float)

	# Wrap-up any requested tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Tool_Drill._match('{1}', {2:i}, {3:i}, '{4}') => {5}".
	      format(indent, tool, desired_diameter, maximum_z_depth, from_routine, priority))

	return priority

    def _point_angle_get(self):
	""" *Tool_Drill*: Return the point angle for the *Tool_Drill* object (i.e. *self*). """

	return self._point_angle

    def _tip_depth_get(self):
        """ *Tool_Drill*: Return the tip depth for the *Tool_Drill* object (i.e. *self*). """

	# Use *tool_drill* instead of *self*:
	tool_drill = self
        
	point_angle = tool_drill._point_angle
	diameter = tool_drill._diameter
	radius = diameter/2
	tip_depth = radius * (Angle(deg=90.0) - point_angle/2).tangent()
	return tip_depth

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
    def _match(tool, maximum_diameter, maximum_z_depth, from_routine, tracing):
	""" *Tool_End_Mill*: Verify that *tool* is both an end mill and that it has a diameter
	    less than or equal to *maximum_diameter*.  If *maximim_diameter* is negative,
	    it will match any end drill.  A positive number that increases with the diameter
	    is returned if a match occurs.  Otherwise, -1.0 is returned if there is no match.
	    *from_routine* is the name of the calling routine and is used for debuggion only.
	"""

	# Set *debug* to *True* to get some tracing:
	debug = False
	debug = True

	# Verify argument types:
	zero = L()
	assert isinstance(tool, Tool)
	assert isinstance(maximum_diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(from_routine, str)
	assert isinstance(tracing, int)
    
	# Perform any requested *tracing*:
	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 0
	    indent = ' ' * tracing
	if tracing_detail >= 1:
	    print("{0}=>Tool_End_Mill._match('{1}', {2:i}, {3:i}, '{4}')".
	      format(indent, tool._name, maximum_diameter, maximum_z_depth, from_routine))

	# Grab some values from *tool*:
	tool_diameter = tool._diameter
	tool_maximum_z_depth = tool._maximum_z_depth
    
	# See whether or not we can return a positive *priority*:
	priority = -123456.0
	is_end_mill = isinstance(tool, Tool_End_Mill)
	tool._search_results_append(is_end_mill, "Is end mill")
	if is_end_mill:
	    if tracing_detail >= 0:
		print("{0}=>Tool_End_Mill.match('{1}', {2:i}, {3:i}, '{4}')".
		  format(indent, tool._name, maximum_diameter, maximum_z_depth, from_routine))

	    # Somehow, the type *bool_* (from numpy) occasionally gets returned. The cast
            # using bool() works around the problem.  This is just weird:
	    diameter_ok = bool(maximum_diameter < zero) or bool(tool_diameter <= maximum_diameter)
	    assert isinstance(diameter_ok, bool)
	    tool._search_results_append(diameter_ok,
	      "Diameter {0:i} < 0 or Diameter {0:i} <= Max Diameter {1:i}".format(
	      tool_diameter, maximum_diameter))
	    if diameter_ok:
		# Verify that tool depth works:
		if tracing_detail >= 0:
		    print("{0}Tool_End_Mill._match: diameter_ok".format(indent))
		if tracing_detail >= 1:
		    print("{0}Tool_End_Mill._match: max_z_depth:{1:i} tool_maximum_z_depth:{2:i}".
		      format(indent, maximum_z_depth, tool_maximum_z_depth))
		z_depth_ok = -maximum_z_depth <= tool_maximum_z_depth
		#tool._search_results_append(z_depth_ok,
		#  "Max Z depth {0:i} <= Tool Max Z Depth {1:i}".format(
		#  -maximum_z_depth, tool_maximum_z_depth))
		if z_depth_ok:
		    if tracing_detail >= 0:
			print("{0}Tool_End_Mill.match: z_depth ok".format(indent))
		    priority = (tool_diameter * 100.0 - tool_maximum_z_depth).inches()
		    if tool._is_laser:
			priority = tool_diameter.inches()

	    if tracing_detail >= 0:
		print("{0}<=Tool_End_Mill._match('{1}', {2:i}, {3:i}, '{4}') => {5}".format(
		  indent, tool._name, maximum_diameter, maximum_z_depth, from_routine, priority))

	assert isinstance(priority, float)

	# Wrap up any requested *tracing*:
	if tracing_detail >= 1:
	    print("{0}=>Tool_End_Mill._match('{1}', {2:i}, {3:i}, '{4}') => {5}".
	      format(indent, tool._name, maximum_diameter, maximum_z_depth, from_routine, priority))

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

    @staticmethod
    def _mill_drill_side_match(tool, maximum_diameter, maximum_z_depth, from_routine, tracing):
	""" *Tool_Mill_Drill*: Verify that {tool} is both a mill drill and that
	    # it has a diameter less than or equal to *maximum_diameter*.  If
	    # *maximim_diameter* is negative, it will match any mill drill.
	    # A positive number that increases with the diameter is returned
	    # if a match occurs.  Otherwise, -1.0 is returned if there is no match.
	"""

	# Verify argument types:
	assert isinstance(tool, Tool)
	assert isinstance(maximum_diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(from_routine, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("=>{0}Tool_Mill_Drill._mill_drill_side_match('{1}', {2:i}, {3:i}, '{4}')".
	      format(indent, tool._name, maximum_diameter, maximum_z_depth, from_routine))

	priority = -1.0
	if isinstance(tool, Tool_Mill_Drill):
	    tip_depth = tool._tip_depth
	    diameter = tool._diameter
	    zero = L()
	    if maximum_diameter < zero or diameter <= maximum_diameter:
		if maximum_z_depth >= -(tool._maximum_z_depth - tip_depth):
		    priority = diameter.millimeters()

	if tracing >= 0:
	    indent = ' ' * tracing
	    print("<={0}Tool_Mill_Drill._mill_drill_side_match('{1}', {2:i}, {3:i}, '{4}')=>{5}".
	      format(indent, tool._name, maximum_diameter, maximum_z_depth, from_routine, priority))

	return priority

    @staticmethod
    def _mill_drill_tip_match(tool, maximum_diameter, maximum_z_depth, from_routine, tracing):
	""" *Tool_Mill_Drill*: Verify that {tool} is both a mill drill and that
	    # it has a diameter less than or equal to *maximum_diameter*.  If
	    # *maximim_diameter* is negative, it will match any mill drill.
	    # A positive number that increases with the diameter is returned
	    # if a match occurs.  Otherwise, -1.0 is returned if there is no match.
	"""

	# Verify argument types:
	assert isinstance(tool, Tool)
	assert isinstance(maximum_diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(from_routine, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("=>{0}Tool_Mill_Drill._mill_drill_tip_match('{1}', {2:i}, {3:i}, '{4}')".
	      format(indent, maximum_diameter, maximum_z_depth, from_routine))

	priority = -1.0
	if isinstance(tool, Tool_Mill_Drill):
	    tip_depth = tool._tip_depth
	    diameter = tool._diameter
	    zero = L()
	    if maximum_diameter < zero or diameter <= maximum_diameter:
		if maximum_z_depth >= -(tool._maximum_z_depth - tip_depth):
		    priority = diameter.millimeters()

	if tracing >= 0:
	    indent = ' ' * tracing
	    print("<={0}Tool_Mill_Drill._mill_drill_tip_match('{1}', {2:i}, {3:i}, '{4}')=>{5}".
	      format(indent, maximum_diameter, maximum_z_depth, from_routine, priority))

	return priority

class Tooling_Plate:
    """ *Tooling_Plate*: A *Tooling_Plate* object represents a flat plate of holes
	onto which parts can be mounted."""

    def __init__(self, dx, dy, dz, rows, columns, hole_pitch,
      hole_diameter, soft_drill_diameter, steel_drill_diameter):
        """ *Tooling_Plate*:  Initialize the *Tooling_Plate* object to contain
	    *dx*, *dy*, *dz*, *rows*, *columns*, *hole_diameter*, *soft_drill_diameter*,
	    *steel_drill_diameter*.
	"""

	# Use *tooling_plate* instead of *self*:
        tooling_plate = self

        # Verify argument types:
        assert isinstance(dx, L)
        assert isinstance(dy, L)
        assert isinstance(dz, L)
	assert isinstance(rows, int) and rows > 0
	assert isinstance(columns, int) and columns > 0
	assert isinstance(hole_pitch, L)
	assert isinstance(hole_diameter, L)
	assert isinstance(soft_drill_diameter, L)
	assert isinstance(steel_drill_diameter, L)

	# Save everything into *tooling_plate*:
        tooling_plate._columns = columns
        tooling_plate._dx = dx
        tooling_plate._dy = dy
        tooling_plate._dz = dz
        tooling_plate._hole_diameter = hole_diameter
	tooling_plate._hole_pitch = hole_pitch
        tooling_plate._rows = rows
        tooling_plate._soft_drill_diameter = soft_drill_diameter
        tooling_plate._steel_drill_diameter = steel_drill_diameter

    def _columns_get(self):
        """ *Tooling_Plate*: Return the number of columns of mounting holes for the *Tooling_Plate*
	    object (i.e. *self*.)
	"""

	return self._columns

    def _drill_diameter_get(self, material):
        """ *Tooling_Plate*: Return the diameter of the drill holes to be drilled into *material*
	    for using the *Tooling_Plate*  object (i.e. *self*.)  These holes are drilled into the
	    *material* to be held down and supsequently threaded.
	"""

	# Use *tooling_plate* instead of *self*:
        tooling_plate = self

	# Verify argument types:
        assert isinstance(material, Material)

	# Grab the correct drill diameter depending upon *material*.
	drill_diameter = tooling_plate._soft_drill_diameter
	if material._is_steel():
	    drill_diameter = tooling_plate._steel_drill_diameter

	return drill_diameter

    def _dx_get(self):
        """ *Tooling_Plate*: Return the dx dimension of the *Tooling_Plate* object (i.e. *self*.)
	"""

	return self._dx

    def _dy_get(self):
        """ *Tooling_Plate*: Return the dy dimension of the *Tooling_Plate* object (i.e. *self*.)
	"""

	return self._dy

    def _dz_get(self):
        """ *Tooling_Plate*: Return the dz dimension of the *Tooling_Plate* object (i.e. *self*.)
	"""

	return self._dz

    def _hole_diameter_get(self):
        """ *Tooling_Plate*: Return the diameter of the mounting holes for the *Tooling_Plate*
	    object (i.e. *self*.)
	"""

	return self._hole_diameter

    def _hole_pitch_get(self):
        """ *Tooling_Plate*: Return the number distance between holes in X and Y for the
	    *Tooling_Plate* object (i.e. *self*.)
	"""

	return self._hole_pitch

    def _rows_get(self):
        """ *Tooling_Plate*: Return the number of rows of mounting holes for the *Tooling_Plate*
	    object (i.e. *self*.)
	"""

	return self._rows

class Transform:

    # The matrix format is an afine 4x4 matrix in the following format:
    #
    #	[ r00 r01 r02 0 ]
    #   [ r10 r11 r12 0 ]
    #   [ r20 r21 r22 0 ]
    #   [ dx  dy  dz  1 ]
    #
    # The afine point format is a 1x4 matrxi of the following format:
    #
    #   [ x y z 1 ]
    #
    # We multiply with the point on the left (1x4) and the matrix on the right (4x4).
    # This yields a 1x4 point matrix of the same form.
    #
    # Only two transforms are supported:
    # * Transform.translate(Point)
    # * Transform.rotate(Point, Angle)
    #
    # A Transform object is immutable:
    # It should have its inverse matrix computed.
    # It should point to the previous matrix and the transform object.
    # This allows us to use the Matrix object to set up openscad.

    def __init__(self):
	""" *Transform*: Initialize the *Transform* object (i.e. *self*) to be the identity
	    transform.
	"""

	# Create *forward_matrix* and *reverse_matrix* as immutable identity matrices:
	forward_matrix = numpy.eye(4)
	#forward_matrix.flags.writable = False
	reverse_matrix = numpy.eye(4)
	#reverse_matrix.flags.writable = False

	# Load up *self*:
	self._forward_matrix = forward_matrix
	self._reverse_matrix = reverse_matrix
	self._forward_scad_lines = ()
	self._reverse_scad_lines = ()

    def __format__(self, format):
	""" *Transform*: Return a formatted version of the *Transform* object"""

	result = ""
	if format.startswith("s"):
	    result = "{0}".format(self._forward_scad_lines)
	elif format.startswith("m"):
	    result = numpy.array_str(self._forward_matrix, precision=8, suppress_small=True)
	    result = result.replace("\n", " ").replace("  ", " ").replace("[ ", "[")
	    if format.endswith("N"):
		result = "np.array(" + result.replace(" ", ",") + ")"
	else:
	    assert False, "Bad format '{0}'".format(format)
	return result

    def __eq__(self, transform2):
	""" *Transform*: Return *True* if the *Transform* object (i.e. *self*)
	    is equal to *transform2*.
	"""

	#print("transform1._forward_scad_lines=", self._forward_scad_lines)
	#print("transform2._forward_scad_lines=", transform2._forward_scad_lines)
	return self._forward_scad_lines == transform2._forward_scad_lines

    def __ne__(self, transform2):
	""" *Transform*: Return *True* if the *Transform* object (i.e. *self*)
	    is not equal to *transform2*.
	"""

	#print("transform1._forward_scad_lines=", self._forward_scad_lines)
	#print("transform2._forward_scad_lines=", transform2._forward_scad_lines)
	return self._forward_scad_lines != transform2._forward_scad_lines

    def __mul__(self, point):
	""" *Transform*:
	"""

	# Verify argument types:
	assert isinstance(point, P)

	x = point.x.millimeters()
	y = point.y.millimeters()
	z = point.z.millimeters()

	afine_point = numpy.array([ [x, y, z, 1.0] ])

	translated_point = afine_point.dot(self._forward_matrix)
	translated_x = L(mm=translated_point[0, 0])
	translated_y = L(mm=translated_point[0, 1])
	translated_z = L(mm=translated_point[0, 2])

	return P(translated_x, translated_y, translated_z)

    #def _scad_lines_get(self):
    #	""" *Transform*: Return the scad lines for the *Transform* object (i.e. *self*). """
    #
    #	return self._forward_scad_lines

    def _scad_lines_append(self, lines, prefix):
	""" *Transform*: Append the openscad transform lines from the *Transform* object
	    (i.e. *self) to *lines* where each line is preceeded by *prefix*.
	"""

	# Verify argument types:
	assert isinstance(lines, list)
	assert isinstance(prefix, str)
		
	# Append each *line* to *lines* preceeded by *prefix*:
	for line in self._forward_scad_lines:
	    lines.append("{0}{1}".format(prefix, line))

    @staticmethod
    def _zero_fix(value):
	""" *Transform*: Return *value* rounding small values to zero and ensure that there
	    is no -0.0.
	"""

	assert isinstance(value, float)
	if abs(value) < 1.0e-10:
	    value = 0.0
	return value

    def reverse(self):
	""" *Transform*: Return the reverse of the *Transform* object (i.e. *self*). """

	assert len(self._forward_scad_lines) == len(self._reverse_scad_lines)

	result = Transform()
	result._forward_matrix = self._reverse_matrix
	result._reverse_matrix = self._forward_matrix
	result._forward_scad_lines = self._reverse_scad_lines
	result._reverse_scad_lines = self._forward_scad_lines

	assert len(result._forward_scad_lines) == len(result._reverse_scad_lines)

	return result

    def rotate(self, comment, axis, angle, tracing = -1000000):
	""" *Transform*: Return a rotation of the *Transform* object (i.e. *self*) rotated
	    by *angle* around *axis*.  *comment* shows up in the .scad file
	"""

	# Verify argument_types:
	assert isinstance(comment, str)
	assert isinstance(axis, P)
	assert isinstance(angle, Angle)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Transform.rotate('{1}', {2:m}, {3:d})".format(indent, comment, axis, angle))

	# Is this a null rotation?:
	if angle == Angle():
	    # Yes it is a null rotation, so we can return *self*:
	    if tracing >= 0:
		print("{0}null rotation: {1:d}".format(indent, angle))
	    result = self
	else:
	    # This is a non-null rotation, we have to compute both the *forward_matix* and
	    # the *reverse_matrix*:

	    # The matrix for rotating by *angle* around the normalized vector (*x*,*y*,*z*) is:
	    #
	    # [ xx(1-c)+c   yx(1-c)-zs  zx(1-c)+ys   0  ]
	    # [ xy(1-c)+zs  yy(1-c)+c   zy(1-c)-xs   0  ]
	    # [ xz(1-c)-ys  yz(1-c)+xs  zz(1-c)+c    0  ]
	    # [ 0           0           0            1  ]
	    #
	    # Where c = cos(*angle*), s = sin(*angle*), and *angle* is measured in radians.

	    # Convert small values to zero.  Also covnvert -0.0 to 0.0:
	    def zero_clean(value):
		if abs(value) < 1.0e-10:
		    value = 0.0
		return value

	    # First, compute the normalized axis values (*nx*, *ny*, and *nz*):
	    normalized = axis.normalize()
	    zf = Transform._zero_fix
	    nx = zf(normalized.x.millimeters())
	    ny = zf(normalized.y.millimeters())
	    nz = zf(normalized.z.millimeters())

	    # Compute some sub expressions for the *forward_matrix*:
	    c = (-angle).cosine()
	    s = (-angle).sine()
	    omc = 1.0 - c
	    x_omc = nx * omc
	    y_omc = ny * omc
	    z_omc = nz * omc
	    xs = nx * s
	    ys = ny * s
	    zs = nz * s
    
	    # Create the *forward_matrix*:
	    forward_matrix = numpy.array([                                           \
	      [zf(nx * x_omc + c),  zf(nx * y_omc - zs), zf(nx * z_omc + ys), 0.0], \
	      [zf(ny * x_omc + zs), zf(ny * y_omc + c),  zf(ny * z_omc - xs), 0.0], \
	      [zf(nz * x_omc - ys), zf(nz * y_omc + xs), zf(nz * z_omc + c),  0.0], \
	      [0.0,	            0.0,                 0.0,                 1.0] ])

	    if tracing >= 0:
		print("{0}forward_matrix={1}".format(indent, forward_matrix))

	    # Perform the *reverse_matrix* computations using -*angle*:
	    c = angle.cosine()
	    s = angle.sine()
	    omc = 1.0 - c
	    x_omc = nx * omc
	    y_omc = ny * omc
	    z_omc = nz * omc
	    xs = nx * s
	    ys = ny * s
	    zs = nz * s
    
	    # Create the *reverse_matrix*:
	    reverse_matrix = numpy.array([                                          \
	      [zf(nx * x_omc + c),  zf(nx * y_omc - zs), zf(nx * z_omc + ys), 0.0], \
	      [zf(ny * x_omc + zs), zf(ny * y_omc + c),  zf(ny * z_omc - xs), 0.0], \
	      [zf(nz * x_omc - ys), zf(nz * y_omc + xs), zf(nz * z_omc + c),  0.0], \
	      [0.0,	            0.0,                 0.0,                 1.0] ])

	    # Create the new transform *result*:
	    result = Transform()
	    result._forward_matrix = self._forward_matrix.dot(forward_matrix)
	    result._reverse_matrix = numpy.linalg.inv(result._forward_matrix)

	    forward_scad_line =                               \
	      ("rotate(a={0:d}, v=[{1}, {2}, {3}]) //F: {4}". \
	      format( angle, nx, ny, nz, comment), )
	    reverse_scad_line =                               \
	      ("rotate(a={0:d}, v=[{1}, {2}, {3}]) //R: {4}". \
	      format(-angle, nx, ny, nz, comment), )
	    result._forward_scad_lines = forward_scad_line + self._forward_scad_lines
	    result._reverse_scad_lines = self._reverse_scad_lines + reverse_scad_line

	# Wrap-up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Transform.rotate('{1}', {2:m}, {3:d})".format(indent, comment, axis, angle))

	return result

    @staticmethod
    def top_surface(comment, start, end, rotate, tracing = -1000000):
	""" *Transform*: Compute the transform to map *start* - *end* vector to start
	    at the origin and with pointing downwards along the z-axis.  A final
	    rotation of *rotate* is performed at the end.
	"""

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(rotate, Angle)

	# Perform any requested *tracing*:
	#tracing = 0
	tracing_detail = -1
	if tracing >= 0:
	    tracing_detail = 3
	    indent = ' ' * tracing
	    print("{0}=>Transform.top_surface('{1}', {2:i}, {3:i}, {4:d})". \
	      format(indent, comment, start, end, rotate))

	# Some constants:
	epsilon = L(mm=.0000001)
	zero = L()
	one = L(mm=1.0)
	z_axis = P(zero, zero, one)
	negative_z_axis = -z_axis
	degrees0 = Angle()
	
	# Start with the identity *transform*:
	transform = Transform()

	# Move *start* down to origin:
	transform = transform.translate("{0}: start to origin".format(comment), -start)
	if tracing_detail >= 2:
	    print("{0}translate_transform_f={1:s}".format(indent, transform))
	if tracing_detail >= 3:
	    print("{0}translate_transform_r={1:s}\n".format(indent, transform.reverse()))

	direction_axis = end - start
	direction = direction_axis.normalize()
	if tracing_detail >= 1:
	    print("{0}direction={1:m}".format(indent, direction))
	if direction.distance(negative_z_axis).absolute() <= epsilon:
	    # We are already pointed in the negative Z axis direction:
	    if tracing >= 0:
		print("{0}We are already pointed in the negative Z axis direction:".format(indent))
	elif direction.distance(z_axis).absolute() <= epsilon:
	    # We are poinnted in the Z axis direction, just rotate by 180 degrees:
	    y_axis = P(zero, one, zero)
	    degrees180 = Angle(deg=180.0)
	    transform = transform.rotate("{0}: z-axis flip".format(comment),
	      y_axis, degrees180, tracing = tracing + 1)
	    if tracing >= 0:
		print("{0}Z-axis flip={1:m}".format(indent, transform))
	else:
	    # We are pointed in some other direction:

	    # Yes, we need to rotate aound the *rotate_axis* by *rotate_angle*:
	    if tracing_detail >= 3:
		print("{0}negative_z_axis={1:i} direction_axis{2:i}".
		  format(indent, negative_z_axis, direction_axis))
	    #plane_change_axis = negative_z_axis.cross_product(direction_axis).normalize()
	    plane_change_axis = direction_axis.cross_product(negative_z_axis).normalize()
	    rotate_angle = negative_z_axis.angle_between(direction_axis)
	    if tracing_detail >= 1:
		print("{0}plane_change_axis={1:m} rotate_angle={2:d}".
		  format(indent, plane_change_axis, rotate_angle))
	    transform = \
	      transform.rotate("{0}: plane change".format(comment),
	      plane_change_axis, rotate_angle, tracing = tracing + 1)
	    if tracing_detail >= 2:
		print("{0}plane_change transform_f={1:s}".format(indent, transform))
	    if tracing_detail >= 3:
		print("{0}plane_change transform_r={1:s}\n".format(indent, transform.reverse()))

	    # For extrusions, we want to be sure that the contour is aligned with
	    # the plane change.  So we project the *direction_axis* onto the the
	    # X/Y plane and rotate to the projected direction:
	    direction_x = direction_axis.x
	    direction_y = direction_axis.y
	    pre_rotate_angle = direction_y.arc_tangent2(direction_x)
	    assert len(transform._forward_scad_lines) == len(transform._reverse_scad_lines)
	    transform = \
	      transform.rotate("{0}: pre-rotate".format(comment), z_axis, pre_rotate_angle)
	    assert len(transform._forward_scad_lines) == len(transform._reverse_scad_lines)
	    if tracing_detail >= 2:
		print("{0}pre_rotate transform_f={1:s}".format(indent, transform))
	    if tracing_detail >= 3:
		print("{0}pre_rotate transform_r={1:s}\n".format(indent, transform.reverse()))

	# Finally, perform the user requested *rotate*:
	transform = transform.rotate("{0}: user rotate".format(comment),
	  z_axis, rotate, tracing = tracing + 1)
	if tracing_detail >= 2:
	    print("{0}user_rotate transform_f={1:s}".format(indent, transform))
	if tracing_detail >= 2:
	    print("{0}user_rotate transform_r={1:s}\n".format(indent, transform.reverse()))

	# Perform any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Transform.top_surface('{1}', {2:i}, {3:i}, {4:d})=>{5:s}".
	      format(indent, comment, start, end, rotate, transform))

	return transform

    def translate(self, comment, dx_dy_dz):
	""" Transform*: Perform a translate of the *Transform* object (i.e. *self*) by *dx_dy_dz*.
	    *comment* will show up in the .scad file.
	"""

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(dx_dy_dz, P)


	# Extract *dx*, *dy*, and *dz*:
	zf = Transform._zero_fix
	dx = zf(dx_dy_dz.x.millimeters())
	dy = zf(dx_dy_dz.y.millimeters())
	dz = zf(dx_dy_dz.z.millimeters())

	# Compute the negative values:
	ndx = zf(-dx)
	ndy = zf(-dy)
	ndz = zf(-dz)

	# Is the translate null?
	if dx == 0.0 and dy == 0.0 and dz == 0.0:
	    # The translate is null, so we can return *self*:
	    result = self
	else:
	    # The translate is non-null, so we have to compute both *forward_matrix* and
	    # *reverse_matrix*:
	    forward_matrix = numpy.array([
	      [1.0, 0.0, 0.0, 0.0],
	      [0.0, 1.0, 0.0, 0.0],
	      [0.0, 0.0, 1.0, 0.0],
	      [dx,  dy,  dz,  1.0]
	    ])
	    reverse_matrix = numpy.array([
	      [1.0, 0.0, 0.0, 0.0],
	      [0.0, 1.0, 0.0, 0.0],
	      [0.0, 0.0, 1.0, 0.0],
	      [ndx, ndy, ndz, 1.0]
	    ])

	    # Now we can create the new *result*:
	    result = Transform()
	    result._forward_matrix = self._forward_matrix.dot(forward_matrix)
	    result._reverse_matrix = numpy.linalg.inv(result._forward_matrix)
	    forward_scad_line = \
	      ("translate([{0}, {1}, {2}]) //F: {3}".format( dx,  dy,  dz, comment), )
	    reverse_scad_line = \
	      ("translate([{0}, {1}, {2}]) //R: {3}".format(ndx, ndy, ndz, comment), )
	    result._forward_scad_lines = forward_scad_line + self._forward_scad_lines
	    result._reverse_scad_lines = self._reverse_scad_lines + reverse_scad_line
	
	return result

class Matrix:
    """ *Matrix* is a class that implements a 4x4 mutable matrix. """

    # FIXME: Warning this stuff basically started out as C code.  It needs to be entirely
    # switched over to the Python numpy class.  You have been warned!!!


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

    def __format__(self, format):
	""" Matrix*: Return the *Matrix* object (i.e. *self*) as a formated string. """

	matrix = self.mat
	item_format = "   {0:" + format + "}"
	item_format = "   {0:10f}"
	result = ""
	for row in range(4):
	    result += "["
	    for column in range(4):
		result += item_format.format(matrix[row, column])
	    result += "  ]\n"
	return result

    def __mul__(self, m):
	""" *Matrix*: Return *self* multiplied by *m*. """

	result = Matrix()
	result.mat = self.mat.dot(m.mat)
	return result

    def copy(self):
	""" *Matrix*: Return a copy of the *Matrix* object (i.e. *self*). """

	copied_matrix = Matrix()
	numpy.copyto(copied_matrix.mat, self.mat)
	return copied_matrix

    @staticmethod
    def identity():
	""" *Matrix*: Return a 4x4 identity matrix. """

	#FIXME: this should be numby.eye(4):
	matrix = Matrix()
	matrix.mat = numpy.eye(4)
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
	""" *Matrix*: Return the point that results from mulitiplying
	    {point} by {self} and converting it back into a {P} with
	    {part} as the reference frame. """

	point_matrix = numpy.array(
	  [ [point.x.millimeters(), point.y.millimeters(), point.z.millimeters(), 1.0] ])
	#print("point_matix=", point_matrix, "shape=", point_matrix.shape)

	mat = self.mat
	assert isinstance(mat, numpy.ndarray)
	translated_point = point_matrix.dot(mat)
	#translated_matrix = point_matrix.dot(mat)
	#print("point_matrix:")
	#print("{0}".format(point_matrix))
	#print("translated_matrix:")
	#print("{0}:".format(translated_matrix))

	x = L(mm=translated_point[0, 0])
	y = L(mm=translated_point[0, 1])
	z = L(mm=translated_point[0, 2])
	result_point = P(x, y, z)
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
	matrix.rotate_store(nx, ny, nz, angle)
	return matrix

    @staticmethod
    def rotate_create(nx, ny, nz, angle):
	""" *Matrix*: Return a rotation matrix for rotating around the normalized vector
	    (*nx*, *ny*, *nz*) by *angle*."""

	# Verify argument types:
	assert isinstance(nx, L)
	assert isinstance(ny, L)
	assert isinstance(nz, L)
	assert isinstance(angle, Angle)

	# Create *rotate_matrix*:
	# 
	# The matrix for rotating by *angle* around the normalized vector (*x*,*y*,*z*) is:
	#
	# [ xx(1-c)+c   yx(1-c)-zs  zx(1-c)+ys   0  ]
	# [ xy(1-c)+zs  yy(1-c)+c   zy(1-c)-xs   0  ]
	# [ xz(1-c)-ys  yz(1-c)+xs  zz(1-c)+c    0  ]
	# [      0	   0	  0	1  ]
	#
	# Where c = cos(*angle*), s = sin(*angle*), and *angle* is measured in radians.

	# Grab the normal X/Y/Z coordinates:
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
	  [0.0,		0.0,		0.0,		1.0] ])

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
	# [      0	   0	  0	1  ]
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
	# [      0	   0	  0	1  ]
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
	  0.0,	       0.0,	       0.0,	       1.0)

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
	# [ 1  0  0  0 ]
	# [ 0  1  0  0 ]
	# [ 0  0  1  0 ]
	# [ dx dy dz 1 ]

	# First convert into an identity *matrix*:
	matrix.identity_store()

	# Now load in *dx*, *dy*, and *dz*:
	mat = self.mat
	mat[3,0] = dx.millimeters()
	mat[3,1] = dy.millimeters()
	mat[3,2] = dz.millimeters()

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

    def __init__(self, volume, parallels):
	""" *Vice*: Initialize the *Vice* object (i.e. *self*) to contain *volume*,
	    "parallels_thickness", and "parallels.  *volume* is a point (i.e. *P*)
	    that specifies the vice volume dimensions between the vice jaws.
	"""

	# Use *vice* instead of *self*:
	vice = self

	# Verify argument types:
	assert isinstance(volume, P)
	assert isinstance(parallels, Parallels)

	# Load up the *Vice* object (i.e. *self*):
	vice._volume = volume
	vice._parallels = parallels

    def _jaw_volume_get(self):
	""" *Vice*: Return the jaw volume of the *Vice* object (i.e. *self*).
	"""

	return self._volume

    def _parallels_get(self):
	""" *Vice*: Return the *Parallels* object for the *Vice* object (i.e. *self*).
	"""

	return self._parallels



#    #FIXME: This Part routine appears to be no longer used!!!
#    def _flush(self, program_number):
#	""" *Code*: Generate CNC code starting at *program_number* for
#	    the first CNC file name.
#	"""
#
#	# Verify argument types:
#	assert isinstance(program_number, int)
#
#	# Use *code* instead of *self*:
#	code = self
#
#	#call d@(form@("=>flush@Code(*, %d%)\n\") / f@(program_number))
#	#original_program_number :@= program_number
#
#	# Output all of the *blocks*:
#	blocks = code._blocks
#	size = len(blocks)
#	print("size=", size)
#	if size != 0:
#	    block0 = blocks[0]
#	    part = block0._part_get()
#	    assert isinstance(part, Part)
#
#	    remainder = program_number % 10
#	    if remainder != 0:
#		program_number += 10 - remainder
#	    ngc_directory = ezcad._ngc_directory_get()
#	    ngc_file_name = os.path.join(ngc_directory, "{0}.ngc".format(program_number))
#	    print("ngc_file_name='{0}'".format(ngc_file_name))
#	    wrl_directory = ezcad._wrl_directory_get()
#	    wrl_file_name = os.path.join(ngc_directory, "O{0}.wrl".format(program_number))
#	    print("wrl_file_name='{0}'".format(wrl_file_name))
#	    wrl_file = open(wrl_file_name, "wa")
#	    assert isinstance(wrl_file, file)
#
#	    # Assign program numbers:
#	    for index in range(size):
#		block = blocks[index]
#		#call d@(form@("Block[%d%] #%d% %d% %v%\n\") % f@(index) %
#		#  f@(block.uid) % f@(block.priority) / f@(block.tool.name))
#		block.program_number = program_number
#		program_number += 1
#
#	    # Open *ngc_file_name* for writing:
#	    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>code_stream_open: '{0}'".format(ngc_file_name))
#	    code_stream = open(ngc_file_name, "wa")
#	    assert isinstance(code_stream, file), \
#	      "Could not open '{0}' for writing".format(ngc_file_name)
#
#	    # Output the tool table:
#	    code_stream.write("({0}: {0})\n".format(part._name_get(), block0._comment_get()))
#	    code_stream.write("(Tooling table:)\n")
#	    for block in blocks:
#		# Fetch {tool}:
#		tool = block._tool_get()
#
#		# Make sure we have assigned a number to the tool:
#		tool_name = tool._name_get()
#		assert tool._number_get() != 0, \
#		  "Tool {0} does not have a number\n".format(tool_name)
#
#		# Output the line in the tool table:
#		code_stream.write("(T{0} {1})\n".format(tool._number_get(), tool_name))
#
#		# Perform the code block:
#		code_stream.write("O{0} call\n".format(block._program_number_get()))
#
#	    # Wrap it up:
#	    code_stream.write("G53 Y0.0 (Move the work to the front)\n")
#	    code_stream.write("M2\n")
#	    code_stream.close()
#	    print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<code_stream_open: '{0}'".format(ngc_file_name))
#
#	    # Now output the blocks:
#	    for block in blocks:
#		# Fetch {block}, {tool} and {name}.
#		tool = block._tool_get()
#		tool_number = tool._number_get()
#
#		is_laser = tool._is_laser_get()
#		if is_laser:
#		    # We have a laser; generate a .dxf file:
#		    dxf_directory = ezcad._dxf_directory_get()
#		    dxf_file_name = os.path.join(dxf_directory, "{0}.dxf".format(program_number))
#		    dxf_stream = open(dxf_file_name, "wa")
#		    assert isinstance(dxf_stream, file), \
#		      "Unable to open file '{0}'".format(dxf_file_name)
#
#		    #dxf_stream.write("0\n\SECTION\n\2\n\HEADER\n")
#		    #dxf_stream.write("9\n\$DIMAUNITS\n\70\n\1\n")
#		    #dxf_stream.write("9\n\$INSUNITS\n\70\n\1\n")
#		    #dxf_stream.write("9\n\$LUNITS\n\70\n\0\n")
#		    #dxf_stream.write("9\n\$MEASUREMENT\n\70\n\0\n")
#		    #dxf_stream.write("0\n\ENDSEC\n")
#
#		    dxf_stream.write("0\n\SECTION\n\2\n\ENTITIES\n")
#
#		    # Output *dxf* to *dxf_stream*:
#		    dxf = code.dxf
#		    dxf_stream.write(dxf)
#
#		    # See if we are aggregating *dxf* into a named DXF file:
#		    dxf_base_name = part.dxf_base_name
#		    if dxf_base_name != "":
#			# Yes, we need to aggregate:
#			shop = part._shop
#			dxf_table = shop.dxf_table
#			dxf_contents = lookup[dxf_base_name]
#			assert isinstance(dxf_contents, str)
#			dxf_contents[dxf_base] = dxf_contents + dxf
#
#		    # Clear *dxf* for the next part:
#		    code.dxf = ""
#
#		    # Wrap up {dxf_stream}:
#		    dxf_stream.write("0\n\ENDSEC\n\0\n\EOF\n")
#		    dxf_stream.close()
#		else:
#		    # We have a mill; generate a .ngc file:
#		    program_number = block.program_number
#		    ngc_directory = ezcad._ngc_directory_get()
#		    ngc_file_name = os.path.join(ngc_directory, "{0}.ngc".format(program_number))
#		    code_stream = open(ngc_file_name, "wa")
#		    print(">>>>>>>>>>>>>>>>code_stream_open: '{0}'".format(file_name))
#		    assert isinstance(code_stream, file), \
#		      "Unable to open '{0}' for output\n".format(ngc_file_name)
#
#		    # Put a visual break into {code_stream}:
#		    part_name = block._part_get()._name_get()
#		    code_stream.write("({0})\n".format(part_name))
#		    code_stream.write("O{0} sub\n".format(program_number))
#		    code_stream.write("G90 G90.1 G20\n")
#
#		    # Output the tool change, get the coolant on,
#		    # and spindle spun up:
#		    code_stream.write("M6 T{0} (Insert {1})\n". 
#	   	      format(tool_number, tool._name_get()))
#		    spindle = block._spindle_get()
#		    if spindle > Hertz():
#			material = part._material_get()
#			if material._needs_coolant():
#			    code_stream.write("M8 (Coolant on)\n")
#			else:
#			    code_stream.write("M9 (Coolant off)\n")
#			code_stream.write("S{0:rpm} M3 (Spindle on)\n".format(spindle))
#
#		    # Get moving to tool change location:
#		    plunge_x = part._plunge_x_get()
#		    plunge_y = part._plunge_y_get()
#		    vice_x = block._vice_x_get()
#		    vice_y = block._vice_y_get()
#		    code._vice_xy_set(vice_x, vice_y)
#
#		    #print("part:{0} code.vice_x={1} code.vice_y={2}\n".
#		    #  format(part.name, code.vice_x, code.vice_y))
#
#		    code_stream.write("G54 G0 X{0} Y{1} (Apply work offset)\n".
#		      format(plunge_x - vice_x, plunge_y - vice_y))
#		    code_stream.write("(Part_Origin_X={0} Part_Origin_Y={1})\n".
#		      format(-vice_x, -vice_y))
#
#		    # Enable tool offset:
#		    z_safe = part._z_safe_get()
#		    code_stream.write("G43 H{0} (Enable tool_offset)\n".
#		      format(tool_number))
#		    code_stream.write(
#		      "G0 Z{0:i} (Go to Z-safe for current tool)\n".
#		      format(z_safe))
#
#		    # These should be done using accessor functions:
#		    code._g8 = 43
#		    code._g1 = 0
#		    code._h = tool_number
#		    code._z = z_safe
#
#		    # Generate the code for {block}:
#		    code._commands_write(code_stream)
#
#		    # Put the tool back into a reasonable position:
#		    code_stream.write("M5 (Spindle off)\n")
#		    code_stream.write("M9 (Coolant off)\n")
#		    code_stream.write(
#		      "G0 X{0} Y{1} (Put spindle to left of vice)\n".
#		      format(0.0, -1.0))
#
#		    # These should be done using accessor functions:
#		    code._x = vice_x - L(inch= 1.0)
#		    code._y = vice_y
#
#		    code_stream.write("G49 (Disable tool offset)\n")
#		    # G53 is not modal, thus we are still in G54 work offset:
#		    code_stream.write(
#		      "G53 G0 Z{0} (Return to machine Z zero)\n".format(1.0))
#		    code.z = L()
#		    code_stream.write("O{0} endsub\n".format(program_number))
#		    code_stream.write("O{0} call\n".format(program_number))
#		    code_stream.write("M2\n")
#		    code_stream.write("%\n")
#		    code_stream.close()
#		    print("<<<<<<<<<<<<<<<<code_stream_open: '{0}'".format(file_name))
#
#		# Now output the VRML code for *block*:
#                
#	    # Close out the *wrl_file*:
#	    wrl_file.close()
#
#	    # Output the shut-down code:
#	    #call put@("(*************************************)\n", code_stream)
#
#	    # Output the subroutine calls:
#	    #index := 0
#	    #while index < size
#	    #    #call put@(form@("O%d% call\n\") /
#	    #    #  f@((index + 1) * 100), code_stream)
#	    #    index := index + 1
#
#	    #call put@(form@("O<%s%> sub\n\") / f@(base_name), code_stream)
#	    #call put@("O[#1 * 100] call\n\", code_stream)
#	    #call put@(form@("O<%s%> endsub\n\") / f@(base_name), code_stream)
#
#	    #call put@("G53 Y0.0 (Move the work to the front)\n", code_stream)
#
#	    #call put@("M2\n\", code_stream)
#    
#	    # Zero {blocks} for the next time:
#	    del blocks[:]
#
#	    #call close@(code_stream)
#
#	#call d@(form@("<=flush@Code(*, %d%) => %d%\n\") %
#	#  f@(original_program_number) / f@(program_number))
#	return program_number

