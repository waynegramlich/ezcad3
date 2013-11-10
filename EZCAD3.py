#!/usr/bin/python

## @package EZCAD3
#
# More details

import os
import math
import subprocess

## The *Angle* class represents an angle.
#
# An *Angle* object represents an angle.  *Angle*'s can be specified
# in either degrees or radians.
class Angle:
    """ An {Angle} represents an angle. """

    ## *PI* respresents the irrational number &pi;.
    PI = 3.14159265358979323846

    ## *TWO_PI* represents 2&pi;.
    TWO_PI = 2 * PI

    ## Initialize *Angle* to be *scalar_radians*.
    #  @param self *Angle* object to initialize.
    #  @param scalar_radians is the value of the angle in radians.
    #
    # <I>__init__</I>() will initialize *self* (an *Angle* object) to
    # *scalar_radians*.
    def __init__(self, scalar_radians = 0.0):
	""" Angle: Initialize self to contain *scalar_radians*. """

	# Check argument types:
	assert isinstance(scalar_radians, float) or \
	  isinstance(scalar_radians, int)

	# Initliaze the final value:
	self.radians = float(scalar_radians)

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
	return Angle(self.radians + angle.radians)

    ## @brief Divides an *Angle* by *divisor*.
    #  @param self is the *Angle* object to divide.
    #  @param divisor the amount to divide the *Angle* bye.
    #  @returns an *Angle* the corresponds to *self* / *divisor.
    #
    # <I>__div__</I>() will divide *self* by *divisor*.  *divisor* must be
    # either a *float* or an *int*.
    def __div__(self, divisor):
	""" Angle: Return *self* divided by *divisor*. """

	assert isinstance(divisor, float) or isinstance(divisor, int)
	return Angle(self.radians / float(scalar))

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

	assert isinstance(format, str)
	return str(self.radians * 180.0 / Angle.PI)

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
	return Angle(self.radians * multiplier)

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

	assert isinstance(angle, Angle)
	return Angle(-self.radians)

    ## @brief Return *self* &times; *multiplier*.
    #  @param self is the *Angle* object to multiply.
    #  @param multiplier is the number to mulitply *self* by.
    #  @returns *self* &times *multiplier*.
    #
    # <I>__rmul__</I>() returns *self* (an *Angle* object) &times; *multiplier*.
    def __rmul__(self, multiplier):
	""" Angle: Return {self} multiplied by {multiplier}. """

	assert isinstance(multiplier, float) or isinstance(multiplier, int)
	return Angle(self.radians * multiplier)

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
	return Angle(self.radians - angle.radians)

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
	return Angle(scalar_degrees * Angle.PI / 180.0)

    ## @brief Returns *self* converted degrees as a *float*.
    #  @param self is the *Angle* to convert to degrees.
    #  @returns *self* converted to degrees.
    #
    # *degrees*() will return *self* converted to degrees
    def degrees(self):
	""" Angle: Convert {self} back into degrees and return it. """

	return self.radians * 180.0 / Angle.PI

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

    def __init__(self, name = "", red = -1.0, green = -1.0, blue = -1.0, alpha = 1.0):
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

class EZCAD3:
    """ EZCAD3 is the top level engine that executes the design. """

    def __init__(self, major, minor):
	""" {EZCAD}: Initialize the contents of {self} to contain
	    {major} and {minor} version numbers. """

	#print "EZCAD.__init__() called"

	# Check argument types:
	assert major == 3
	assert minor == 0

	# Load up {self}:
	self._major = major
	self._minor = minor
	self._parts_stack = []
	self._xml_indent = 0
	self._xml_stream = None

    def process(self, part):
	assert isinstance(part, Part)
	part.process(self)

class L:
    """ A {L} represents a length. """

    # Currently {L} is represented in inches:

    def __init__(self, mm=0.0, cm=0.0, inch=0.0, ft=0.0):
	""" L: Create {L} with a value of {value}. """

	assert isinstance(mm, float) or isinstance(mm, int)
	assert isinstance(cm, float) or isinstance(cm, int)
	assert isinstance(inch, float) or isinstance(inch, int) or \
	       isinstance(inch, str)
	assert isinstance(ft, float) or isinstance(ft, int)

	if isinstance(inch, str):
	    # We have a string, parse it:
            whole_fraction = inch.split("-")
	    if len(whole_fraction) == 2:
		whole = float(whole_fraction[0])
		faction = whole_fraction[1]
            elif len(whole_fraction) == 1:
		whole = 0.0
		fraction = inches
	    else:
		assert False, "Poorly formed inches string '{0}'".format(inches)
	    numerator_denominator = fraction.split("/")
	    assert len(numerator_denominator) == 2, \
	      "Bad fraction '{0}'".format(fraction)
	    numerator = float(numerator_denominator[0])
	    denominator = float(numerator_denominator[1])
	    inch = whole + numerator / denominator
		
	self._mm = mm + cm * 10.0 + inch * 25.4 + ft * (12.0 * 25.4)

    def __add__(self, length):
	""" L: Add {L} to {self}. """

	assert isinstance(length, L)
	return L(self._mm + length._mm)

    def __div__(self, number):
	""" L: Divide {self} by {scalar}. """
	assert isinstance(number, float) or isinstance(number, int)

	return L(self._mm / number)

    def __eq__(self, length):
	""" L: Return {True} if {self} is equal to {length}. """

	assert isinstance(length, L)
	return self._mm == length._mm

    def __ge__(self, length):
	""" L: Return {True} if {self} is greater than or equal
	    to {length}. """

	assert isinstance(length, L)
	return self._mm >= length._mm

    def __lt__(self, length):
	""" L: Return {True} if {self} is greater than {length}. """

	assert isinstance(length, L)
	return self._mm > length._mm

    def __le__(self, length):
	""" L: Return {True} if {self} is less than or equal
	    to {length}. """

	assert isinstance(length, L)
	return self._mm <= length._mm

    def __lt__(self, length):
	""" L: Return {True} if {self} is less than {length}. """

	assert isinstance(length, L)
	return self._mm < length._mm

    def __format__(self, fmt):
	""" L: Format {self} into a string an return it. """

	assert isinstance(fmt, str)
	mm = self._mm
	if fmt == 'i':
	    mm /=  25.4
	elif fmt == 'f':
            mm /= (25.4 * 12)
	elif fmt == 'c':
	    mm /= 10.0
	return "{0}".format(mm)

    def __mul__(self, number):
	""" L: Multiply {self} by {scalar}. """

	assert isinstance(number, float) or isinstance(number, int)
	return L(self._mm * number)

    def __rmul__(self, number):
	""" L: Multiply {self} by {scalar}. """

	assert isinstance(number, float) or isinstance(number, int)
	return L(self._mm * number)

    def __ne__(self, length):
	""" L: Return {True} if {self} is not equal to {length}. """

	assert isinstance(length, L)
	return self._mm != length._mm

    def __neg__(self):
	""" L: Return the negative of {self}. """

	return L(-self._mm)

    def __str__(self):
	""" L: Return {self} as a formated string. """

	return str(self._mm)

    def __sub__(self, length):
	""" L: Subtract {length} from {self}. """

	assert isinstance(length, L)
	return L(self._mm - length._mm)

    def abs(self):
	""" L: Return absolue value of {self}. """

	assert isinstance(length, L)
	return L(abs(self._mm))

    def arc_tangent2(self, dx):
	""" L: Return the arctangent of {self} over {dx} being careful
	    that the returned angle is in the correct quadrant. """

	assert isinstance(dx, L)
	return Angle.rad(math.atan2(self._mm, dx._mm))

    def centimeters(self):
	""" L: Return {self} as a scalar measured in centimeters. """

	return self._mm / 10.0

    def cosine(self, angle):
	""" L: Return {self} * cos(angle). """

	assert isinstance(angle, Angle)
	return L(self._mm * angle.cosine())

    def half(self):
	""" L: Return half of {self}. """

	return L(self._mm / 2.0)

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

	return self._mm * 25.4

    def sine(self, angle):
	""" L: Return {self} * sin(angle). """

	return L(self._mm * angle.sine())

    def twice(self):
	""" L: Return twice of {self}. """

	return L(self._mm * 2.0)

class Part:
    """ A {Part} specifies either an assembly of parts or a single
	physical part. """

    # Flavors of values that can be stored in a {Part}:
    def __init__(self, up, place = True):

	#print("=>Part.__init__(*, '{0}', *, place={1})".format(name, place))

        # Check argument types:
	none_type = type(None)
	up_type = type(up)
	assert up_type == none_type or isinstance(up, Part)
	assert isinstance(place, bool)

	# Some useful abbreviations:
	zero = L()

	# Find all the child *Part*'s:
	#children = []
	#for attribute_name in dir(self):
	#    attribute = getattr(self, attribute_name)
        #    if isinstance(attribute, Part):
	#	children.append(attribute)
	#	#print("{0}['{1}'] = {2}". \
	#	#  format(name, attribute_name, attribute._name))

	# Load up *self*:
	name = self.__class__.__name__
	self._ezcad = None
	self._name = name
	self._material = None
	self._color = None
	self._places = {}
	self._scad_difference_lines = None
	self._scad_union_lines = None
	self._update_count = 0
	self.up = up

	# Initialize the bounding box information:
	big = L(mm=987654321.0)
	self._box_changed_count = 0
	self._box_points = {}
	self.ex = -big
	self.wx = big
	self.ny = -big
	self.sy = big
	self.tz = -big
	self.bz = big
	self._box_recompute()

	if up_type != none_type and place:
	    #print("Part.__init__: perform place")
	    place = up.place(self, name = name)
	    #print("Part.__init__: place = {0:m}".format(place))

	#print("<=Part.__init__(*, '{0}', *, place={1})".format(name, place))

    ## @brief Formats *self* into a string and returns it.
    #  @param format is the format control string (currently ignored).
    #  @returns a string representation of *self*.
    #
    # *<I>__format__</I>() will return a *self* (a *Box*) object as
    # formatted string.  Currently *format* is ignored, but maybe some
    # time in the future it will contol formaating of the returned string.

    def __format__(self, format):
	""" *Box*: Return formated version of *self*. """
	assert isinstance(format, str)
	if format != "":
	    format = ":" + format
	format_string = "{0}[{1" + format + "}:{2" + format + "},{3" + \
	  format + "}:{4" + format + "}:{5" + format + "}:{6" + format + "}]"
	#print("format_string='{0}'".format(format_string))

	print("[{0}:{1},{2}:{3},{4}:{5}]".format(
	  self.wx._mm, self.ex._mm,
	  self.sy._mm, self.ny._mm,
	  self.sz._mm, self.tz._mm))

	return format_string.format(self._name,
	  self.wx, self.ex,
	  self.sy, self.ny,
	  self.sz, self.tz)

    def __getattr__(self, name):
	if self._update_count == 0:
	    if name.endswith("_"):
		pass
	    elif name.endswith("_a"):
		return Angle()
	    elif name.endswith("_c"):
		return Color()
	    elif name.endswith("_f"):
		return 0.0
	    elif name.endswith("_i"):
		return 0
	    elif name.endswith("_l"):
		return L()
	    elif name.endswith("_m"):
		return Material()
            elif name.endswith("_o"):
		return None
	    elif name.endswith("_s"):
		return ""
	    elif name.endswith("_p"):
		return Point()
	    elif name.endswith("_pl"):
		return Place()
	raise AttributeError(
	  "Part instance has no attribute '{0}'".format(name))

    ## @brief Updates *name*'d *point* is *self* for bounding box calcuation.
    #  @param self is the *Part* to update.
    #  @param name is the name of the *Point* object.
    #  @param point is the *Point* to update.
    #
    # <I>_box_point_update</I>() will update the *Point* named *name* in *self*
    # (a *Part* object) to be *point*.  If the value of *point* has changed
    # from the last time it was updated, the bounding box for *self* is
    # recomputed.
    def _box_point_update(self, name, point):
	""" *Box*: Insert/update the point named *name* to *point*. """

	# Check argument types:
	assert isinstance(name, str)
	assert isinstance(point, Point)

	#print("Box.point_update({0:m}, '{1}', {2:m})".
	#  format(self, name, point))
	# Deterimine if this is an update or the initial insert:
	points = self._box_points
	if name in points.keys():
	    # This is an update; determine if *point* changed:
	    previous = points[name]
	    if previous.x != point.x or \
	      previous.y != point.y or previous.z != point.z:
		# *point* changed, so update everything:
		points[name] = point
		self._box_recompute()
	else:
	    # This is the initial insert; so update everything:
	    points[name] = point
	    self._box_recompute()

    ## @brief Recomputes the bounding box for *self* and any parent *Box*'s.
    #  @param *self* the bounding *Box* to recompute.
    #
    # <I>_recompute</I>() will recompute the bounding box for *self* (a bounding
    # *Box*) and enclosing parent bounding *Box*'s.
    def _box_recompute(self):
	""" *Box*: (Internal use only) Recompute the bounding box corners. """

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
	for point in self._box_points.values():
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
	    # Bounding box changed:
	    self.ex  = ex
	    self.wx = wx
	    self.ny = ny
	    self.sy = sy
	    self.tz = tz
	    self.bz = bz

	    # Compute averate X/Y/Z:
	    cx = (ex + wx).half()
	    cy = (ny + sy).half()
	    cz = (tz + bz).half()

	    # There are 27 points to compute -- 9 on the top slice, 9 on
	    # the middle slice, and 9 on the bottom slice:

	    # Top slice:
	    #	TNW	TN	TNE
	    #	TW	T	TE
	    #	TSW	TS	TSE
	    self.t =   Point(cx, cy, tz)
	    self.te =  Point(ex, cy, tz)
	    self.tn =  Point(cx, ny, tz)
	    self.tne = Point(ex, ny, tz)
	    self.tnw = Point(wx, ny, tz)
	    self.ts =  Point(cx, sy, tz)
	    self.tse = Point(ex, sy, tz)
	    self.tsw = Point(wx, sy, tz)
	    self.tw =  Point(wx, cy, tz)

	    # Middle slice:
	    #	NW	N	NE
	    #	W	C	E
	    #	SW	S	SE
	    self.e =   Point(ex, cy, cz)
	    self.n =   Point(cx, ny, cz)
	    self.ne =  Point(ex, ny, cz)
	    self.nw =  Point(wx, ny, cz)
	    self.c =   Point(cx, cy, cz)
	    self.s =   Point(cx, sy, cz)
	    self.se =  Point(ex, sy, cz)
	    self.sw =  Point(wx, sy, cz)
	    self.w =   Point(wx, cy, cz)

	    # Bottom slice:
	    #	BNW	BN	BNE
	    #	BW	B	BE
	    #	BSW	BS	BSE
	    self.b =   Point(cx, cy, bz)
	    self.be =  Point(ex, cy, bz)
	    self.bn =  Point(cx, ny, bz)
	    self.bne = Point(ex, ny, bz)
	    self.bnw = Point(wx, ny, bz)
	    self.bs =  Point(cx, sy, bz)
	    self.bse = Point(ex, sy, bz)
	    self.bsw = Point(wx, sy, bz)
	    self.bw =  Point(wx, cy, bz)

	    # Recursively update parent *Box*'s until top-most bounding
	    # box is reached:
	    up = self.up
	    if type(up) != type(None):
		up._box_recompute()

            # Keep track if we have changed:
	    self._box_changed_count += 1

    def _dimensions_update(self, trace):

	# Do any requested tracing:
	if trace >= 0:
	    print("{0}=>Part._dimensions_update('{1}')". \
	      format(' ' * trace, self._name))
	    print("{0}Part._dimensions_update:places={1}". \
	      format(' ' * trace, self._places))

	# Start with nothing *changed*:
	changed = 0

	# First record the current values load into *self*.
	before_values = {}
	name = self._name
	for attribute_name in dir(self):
	    attribute = getattr(self, attribute_name)
	    if attribute_name.startswith("__"):
		pass
	    elif attribute_name.endswith("_"):
		assert isinstance(attribute, Part), \
		  "{0}.{1} is not a Part".format(name, attribute_name)
		changed += attribute._dimensions_update(trace + 1)
	    elif attribute_name.endswith("_l"):
		assert isinstance(attribute, L), \
		  "{0}.{1} is not an L (i.e. length)". \
		  format(name, attribute_name)
		before_values[attribute_name] = attribute
	    elif attribute_name.endswith("_p"):
		assert isinstance(attribute, Point), \
		  "{0}.{1} is not a Point (i.e. point)". \
		  format(name, attribute_name)
		before_values[attribute_name] = attribute
	    elif attribute_name.endswith("_a"):
		assert isinstance(attribute, Angle), \
		  "{0}.{1} is not an Angle".format(name, attribute_name)
		before_values[attribute_name] = attribute
	    elif attribute_name.endswith("_f"):
		assert isinstance(attribute, float), \
		  "{0}.{1} is not a float".format(name, attribute_name)
		before_values[attribute_name] = attribute
	    elif attribute_name.endswith("_i"):
		assert isinstance(attribute, int), \
		  "{0}.{1} is not an int".format(name, attribute_name)
		before_values[attribute_name] = attribute
	    #else ignore *attribute_name*:

	# Remember 
	before_box_changed_count = self._box_changed_count

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

		if trace >= 0:
		    print("{0}Part._dimensions_update:{1}.{2} ({3}=>{4})". \
		      format(' ' * trace, name, attribute_name,
		      before_value, after_value))

	# Update bounding box with placed *Part* bounding boxes:
	for place in self._places.values():
	    # Grab some values from *place* and *part*:
	    place_part = place._part
	    place_name = place._name
	    forward_matrix = place._forward_matrix

	    if trace >= 0:
		print("{0}Part._dimensions_update:place={1}". \
		  format(' ' * trace, place))
		#print("{0}Part._dimensions_update:Merge {1:m} into {2:m}". \
		#  format(' ' * trace, place_box, box))

	    self._box_point_update(place_name + "[TNE]",
	      forward_matrix.point_multiply(place_part.tne))
	    self._box_point_update(place_name + "[TNW]",
	      forward_matrix.point_multiply(place_part.tnw))
	    self._box_point_update(place_name + "[TSE]",
	      forward_matrix.point_multiply(place_part.tse))
	    self._box_point_update(place_name + "[TSW]",
	      forward_matrix.point_multiply(place_part.tsw))
	    self._box_point_update(place_name + "[BNE]",
	      forward_matrix.point_multiply(place_part.bne))
	    self._box_point_update(place_name + "[BNW]",
	      forward_matrix.point_multiply(place_part.bnw))
	    self._box_point_update(place_name + "[BSE]",
	      forward_matrix.point_multiply(place_part.bse))
	    self._box_point_update(place_name + "[BSW]",
	      forward_matrix.point_multiply(place_part.bsw))

	# Determine whether *box* has changed:
	after_box_changed_count = self._box_changed_count
	if before_box_changed_count != after_box_changed_count:
	    changed += 1

	if trace >= 0:
	    print("{0}<=Part._dimensions_update('{1}')=>{2}". \
	      format(' ' * trace, self._name, changed))

	return changed

    def _manufacture(self, ezcad):
	""" Part: Output the XML for *self* to *xml_stream*
	    indented by *indent*. """

	#print("=>Part._manufacture:{0}".format(self._name))

	# Check argument types:
	assert isinstance(ezcad, EZCAD3)
	self._ezcad = ezcad

	# First manufacture any child nodes:
	for attribute_name in dir(self):
	    if not attribute_name.startswith("__") and \
	      attribute_name.endswith("_"):
		attribute = getattr(self, attribute_name)
		assert isinstance(attribute, Part), \
		  "{0}.{1} is not a Part".format(self.name, attribute)
		attribute._manufacture(ezcad)

	# Now manufacture this node:
	scad_difference_lines = []
	scad_union_lines = []
	self._scad_difference_lines = scad_difference_lines
	self._scad_union_lines = scad_union_lines
	self.construct()
	self._scad_difference_lines = None
	self._scad_union_lines = None

	# Open *scad_file*:
	name = self._name
	scad_file = open(name + ".scad", "w")

	# Deal with the part *places*:
	lines = []
	place_parts = {}
	places = self._places.values()
	for place in places:
	    place_part = place._part
	    place_parts[place_part._name] = place_part
	for part in place_parts.values():
	    scad_file.write("use <{0}.scad>;\n".format(part._name))

	# Write out the module:
	scad_file.write("module {0}() {{\n".format(name))

	if len(scad_union_lines) != 0:
	    # Get the difference() followed union():
	    scad_file.write("  difference() {\n")
	    scad_file.write("    union() {\n")

	    # Output the *scan_union_lines*:
	    for union_line in scad_union_lines:
		scad_file.write(union_line)
		scad_file.write("\n")

	    # Close off union():
	    scad_file.write("    }\n")

	    # Output *scad_difference_lines*:
	    for difference_line in scad_difference_lines:
		scad_file.write(difference_line)
		scad_file.write("\n")

	    # Close off difference():
	    scad_file.write("  }\n")

        for place in places:
	    self._scad_transform(lines, center = place._center,
	      axis = place._axis, rotate = place._rotate,
	      translate = place._translate);
	    lines.append("{0}();".format(place._part._name))

	for line in lines:
	    scad_file.write(line)
	    scad_file.write("\n")

	# Close off the module:
	scad_file.write("}\n\n")

	# Call the module we just wrote out:
	scad_file.write("{0}();\n".format(name))

	# Close *scad_file*:
	scad_file.close()

	# Run the command:
	ignore_file = open("/dev/null", "w")
	command = [ "openscad",
	  "-o", "{0}.stl".format(name),
	  "{0}.scad".format(name) ]
	print("command=", command)
	subprocess.call(command, stderr=ignore_file) 
	ignore_file.close()

	if False:
	    # For now, write out an offset file:
	    stl_file = open("{0}.stl".format(name), "r")
	    stl_lines = stl_file.readlines()
	    stl_file.close()

	    triangles = []
	    offsets = []
	    vertices = {}
	    size = len(stl_lines)
	    assert stl_lines[0][:5] == "solid"
	    index = 1
	    while index + 4 < size:
		#normal_list = stl_lines[index].split()
		list = stl_lines[index + 2].split()
		vertex1 = (float(list[1]), float(list[2]), float(list[3]))
		list = stl_lines[index + 3].split()
		vertex2 = (float(list[1]), float(list[2]), float(list[3]))
		list = stl_lines[index + 4].split()
		vertex3 = (float(list[1]), float(list[2]), float(list[3]))

		if vertex1 in vertices:
		    offset1 = vertices[vertex1]
		else:
		    offset1 = len(offsets)
                    offsets.append(vertex1)
		    vertices[vertex1] = offset1

		if vertex2 in vertices:
		    offset2 = vertices[vertex2]
		else:
		    offset2 = len(offsets)
                    offsets.append(vertex2)
		    vertices[vertex2] = offset2

		if vertex3 in vertices:
		    offset3 = vertices[vertex3]
		else:
		    offset3 = len(offsets)
                    offsets.append(vertex3)
		    vertices[vertex3] = offset3

		triangles.append( (offset1, offset2, offset3,) )
		index += 7

	    offset_file = open("{0}.off".format(name), "w")
	    offset_file.write("OFF\n\n")
	    offset_file.write(
	      "{0} {1} 0\n".format(len(offsets), len(triangles)))
	    for offset in offsets:
		offset_file.write(
		  "{0} {1} {2}\n".format(offset[0], offset[1], offset[2]))
	    offset_file.write("\n")
	    for triangle in triangles:
		offset_file.write("3 {0} {1} {2}\n". \
		  format(triangle[0], triangle[1], triangle[2]))
	    offset_file.close()

	#print("<=Part._manufacture:{0}".format(self._name))

    def construct(self):
	assert False, \
	  "No construct() method defined for part '{0}'".format(self._name)

    def process(self, ezcad):
	""" {Part}: Generate the XML control file for *self*. """

	assert isinstance(ezcad, EZCAD3)

	# Do the dimensions propogate phase:
	self._update_count = 0
	changed = 1
	while changed != 0:
	    # Find all the child *Part*'s:
            self._update_count += 1
	    print("Dimensions update {0}".format(self._update_count))
	    changed = self._dimensions_update(-1000000)
	    print("Part.process: {0} dimensions_changed\n".format(changed))

	# Open the XML output stream:
	#xml_file_name = self._name + ".xml"
	#xml_stream = open(xml_file_name, "w")
	#assert xml_stream != None, \
	#  "Unable to open XML output file '%s'" % (xml_file_name)

	# Now visit *self* and all of its children:
	self._manufacture(ezcad)

	# Close the XML Stream.
	#ezcad.xml_indent_pop()
	#assert ezcad._xml_indent == 0, \
	#  "XML indentaion bracketing failure"
	#xml_stream.write("</EZCAD>\n")
	#xml_stream.close()
	#ezcad.xml_stream = None

	# Now feed {xml_file_name} into EZCAD_XML:
	#ezcad_xml = os.popen("EZCAD_XML {0}".format(xml_file_name), 'r')
	#ezcad_xml.read()
	#assert ezcad_xml.close() == None, \
	#  "Error running 'EZCAD_XML {0}'".format(xml_file_name)

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

    def block(self, comment = "no comment",  material = None, color = None,
      corner1 = None, corner2 = None,
      center = None, axis = None, rotate = None, translate = None):
	""" {Part} construct: Create a block with corners at {corner1} and
	    {corner2}.  The block is made of {material} and visualized as
	    {color}. """

	#print "block_corners('{0}', {1}, {2}, '{3}', '{4}')".format( \
	#  self.name, corner1, corner2, color, material)

	# Deal with argument defaults:
	none_type = type(None)
	if type(material) == none_type:
	    material = self._material
	    if type(material) == none_type:
		material = Material("plastic", "abs")
	if type(color) == none_type:
	    color = self._color
	    if type(color) == none_type:
		color = Color()
	if type(corner1) == none_type:
	    zero = L(0.0)
	    corner1 = Point(zero, zero, zero)
	if type(corner2) == none_type:
	    one = L.mm(1.0)
	    corner2 = Point(one, one, one)

	# Record the materials:
	self._material = material
	self._color = color

	# Check argument types:
	assert isinstance(comment, str)
	assert isinstance(corner1, Point)
	assert isinstance(corner2, Point)
	assert isinstance(material, Material)
	assert isinstance(color, Color)
	assert type(center) == none_type or isinstance(center, Point)
	assert type(axis) == none_type or isinstance(axis, Point)
	assert type(rotate) == none_type or isinstance(rotate, Angle)
	assert type(translate) == none_type or isinstance(translate, Point)

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

	tne = Point(x2, y2, z2)
	bsw = Point(x1, y1, z1)
	tsw = Point(x1, y1, z2)
	bnw = Point(x1, y2, z1)
	tnw = Point(x1, y2, z2)
	bse = Point(x2, y1, z1)
	tse = Point(x2, y1, z2)
	bne = Point(x2, y2, z1)

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

	union_lines = self._scad_union_lines
	if type(union_lines) != none_type:

	    assert x1 < x2, "x1={0} should be less than x2={1}".format(x1, x2)
	    assert y1 < y2, "y1={0} should be less than y2={1}".format(y1, y2)
	    assert z1 < z2, "xz={0} should be less than z2={1}".format(z1, z2)

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

    def cnc_flush(self):
	""" Part constuct: """

	ezcad = self._ezcad
	xml_stream = ezcad._xml_stream
	if xml_stream != None:
	    xml_stream.write('{0}<CNC_Flush />\n'. \
	      format(" " * ezcad._xml_indent))

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

	# Look up {bsw} and {tne}, the bounding box corner {Point}'s:
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
	ew2 = (e + w).half()
	ns2 = (n + s).half()
	tb2 = (t + b).half()

	# Install 6 bounding box surface {Point}'s into {points}:
	self.xb = self.point_new(ew2, ns2,   b)
	self.xe = self.point_new(  e, ns2, tb2)
	self.xn = self.point_new(ew2,   n, tb2)
	self.xs = self.point_new(ew2,   s, tb2)
	self.xt = self.point_new(ew2, ns2,   t)
	self.xw = self.point_new(  w, ns2, tb2)

	# Install 12 bounding box edge {Point}'s into {points}:
	self.xbe = self.point_new(  e, ns2,   b)
	self.xbn = self.point_new(ew2,   n,   b)
	self.xbs = self.point_new(ew2,   s,   b)
	self.xbw = self.point_new(  w, ns2,   b)
	self.xne = self.point_new(  e,   n, tb2)
	self.xnw = self.point_new(  w,   n, tb2)
	self.xse = self.point_new(  e,   s, tb2)
	self.xsw = self.point_new(  w,   s, tb2)
	self.xte = self.point_new(  e, ns2,   t)
	self.xtn = self.point_new(ew2,   n,   t)
	self.xts = self.point_new(ew2,   s,   t)
	self.xtw = self.point_new(  w, ns2,   t)

	# Install 8 bounding box corner {Point}'s into {points}:
	self.xbne = self.point_new(e, n, b)
	self.xbnw = self.point_new(w, n, b)
	self.xbse = self.point_new(e, s, b)
	self.xbsw = self.point_new(w, s, b)
	self.xtne = self.point_new(e, n, t)
	self.xtnw = self.point_new(w, n, t)
	self.xtse = self.point_new(e, s, t)
	self.xtsw = self.point_new(w, s, t)

	xbsw = self.xbsw
	xtne = self.xtne
	ezcad = self._ezcad
	xml_stream = ezcad._xml_stream
	self._xml_lines.append( ('{0}<Extra C1X="{1}" C1Y="{2}" C1Z="{3}"' + \
	  ' C2X="{4}" C2Y="{5}" C2Z="{6}"/>\n'). \
	  format(" " * ezcad._xml_indent,
	  xbsw.x, xbsw.y, xbsw.z, xtne.x, xtne.y, xtne.z))

    def extra_xyz(self, dx, dy, dz):
	""" {Part}: Add some extra material the block of {self} by {dx},
	    {dy}, and {dz} in the X, Y, and Z dimensions. """

	# Argument type checking:
	assert isinstance(dx, L)
	assert isinstance(dy, L)
	assert isinstance(dz, L)

	# Pass everything on to {extra_ewnstb}:
	half_dx = dx.half()
	half_dy = dy.half()
	half_dz = dz.half()
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
	assert isinstance(start_point, Point)
	assert isinstance(end_point, Point)
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

    def hole(self, comment, diameter, start, end, flags=""):
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

	# Check argument types:
	assert isinstance(comment, str)
	assert isinstance(diameter, L)
	assert isinstance(start, Point)
	assert isinstance(end, Point)

	# Define some useful abbreviations:
	zero = L(0.0)

	# Extract some values from {ezcad}:
	ezcad = self._ezcad

	none_type = type(None)
	difference_lines = self._scad_difference_lines
	if type(difference_lines) != none_type:
	    print("Part.hole:start={0:m} end={1:m}".format(start, end))
	    axis = start - end
	    length = axis.length()
	    assert length > zero, "Cylinder '{0}' has no height".format(comment)

	    print("Part.hole:axis={0:m} length={1:m}".format(axis, length))
	    if axis.x != zero or axis.y != zero:
		# We have to tilt the cylinder.
		z_axis = Point(zero, zero, L(mm=1.0))
		angle = z_axis.angle_between(axis)
		orthogonal_xis = z_axis.dot_product(axis)

	    xxx = start + end
	    print("Part.hole:xxx={0:m}".format(xxx))
	    center = (start + end).half()
	    print("Part.hole:center={0:m}".format(center))
	    difference_lines.append(
	      "    // '{0}' Hole".format(comment))
	    difference_lines.append(
	      "    translate({0:m})".format(center))
            difference_lines.append(
	      "      cylinder(r={0:m}, h={1:m}, center = true, $fn=12);".
	      format(diameter.half(), length + L(mm=1.0)))
	
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

    def place(self, part, name = None,
      center = None, axis = None, rotate = None, translate = None):
	""" Part dimensions: Place {place_part} at {translate_point} relative
	    to {self} with no rotation.  {place_name} is used for point
	    paths. """

	#print("Part.place({0}, part='{1}',name='{2}' ...)". \
	#  format(self._name, part._name, name))

	# Deal with default arguments:
	none_type = type(None)
	if type(name) == none_type:
	    name = part._name
	if type(center) == none_type:
	    center = Point()
	if type(axis) == none_type:
	    axis = Point(z = L(1.0))
	if type(rotate) == none_type:
	    rotate = Angle()
	if type(translate) == none_type:
	    translate = Point()

	# Check argument types:
	assert isinstance(name, str)
	assert isinstance(part, Part)
	assert isinstance(center, Point)
	assert isinstance(axis, Point)
	assert isinstance(rotate, Angle)
	assert isinstance(translate, Point)

	# Create *place* and stuff into *_places*:
	place = Place(part = part, name = name,
	  center = center, axis = axis, rotate = rotate, translate = translate)
	self._places[name] = place

	return place

    #def point(self, point_path):
    #	""" Part dimensions: Return the {Point} associated with {point_path}
    #	    starting from {self}. """
    #
    #	return self.value_lookup(point_path, Part.POINTS)

    def point_map(self, place_path):
	""" Part dimensions: Return the {Point} that corresponds to
	    the point specified by {place_path} in the {self} reference
	    frame. """

	assert not self.part_tree_mode(), \
	  "Part.point_subtract: Not allowed in part tree mode"

	# Wait until after all part names have been defined:
	if self.dimensions_define_mode():
	    zero = L.inch(0)
	    result = Point(self, zero, zero, zero)
	else:
	    # Look up the {Part}/{Point}/{Matrix} triple for {place_path1}:
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

    def _scad_transform(self, lines,
       center = None, axis = None, rotate = None, translate = None):

	none_type = type(None)
	assert type(lines) == none_type or isinstance(lines, list)
	assert type(center) == none_type or isinstance(center, Point)
	assert type(axis) == none_type or isinstance(axis, Point)
	assert type(rotate) == none_type or isinstance(rotate, Angle)
	assert type(translate) == none_type or isinstance(translate, Point)

	if isinstance(lines, list):
	    # Perform any requested translate "last":
	    if type(translate) != none_type:
		lines.append(
		  "      translate([{0:m}, {1:m}, {2:m}])".format(
		  translate.x, translate.y, translate.z))

	    # The rotation only makes sense if there is a *rotate* angle:
	    if type(rotate) != none_type:
		assert isinstance(rotate, Angle)

		# Again, we do the transforms in reverse order.
		# This transform restores the part to its natural location:
		if type(center) == none_type:
		    # Rotate around the *Part* center:
		    lines.append(
		      "    translate([{0:m}, {1:m}, {2:m}])".format(
		      (x1 + x2).half(), (y1 + y2).half(), (z1 + z2).half()))
		else:
		    # Rotate around *center*:
		    assert isinstance(center, Point)
		    lines.append(
		      "      translate([{0:m}, {1:m}, {2:m}])".format(
		      center.x, center.y, center.z))

		# Now do the axis rotation now that the part center is
		# at the origin:
		if type(axis) == none_type:
		    # Rotate around the Z axis:
		    lines.append(
		      "    rotate(a={0}, v=[0, 0, 1])".format(rotate))
		else:
                    # Rotate around *axis*:
		    assert isinstance(axis, Point)
		    lines.append(
		      "      rotate(a={0}, v=[{1:m}, {2:m}, {3:m}])".format(
		      rotate, axis.x, axis.y, axis.z))

		# This is the transform that moves the part to origin
		# such prior to rotation:
		if type(center) == none_type:
		    # Rotate using the center of the cube:
		    lines.append(
		      "    translate([{0:m}, {1:m}, {2:m}])".format(
		      -(x1 + x2).half(), -(y1 + y2).half(), -(z1 + z2).half()))
		else:
		    # Rotate around *center*:
		    lines.append(
		      "      translate([{0:m}, {1:m}, {2:m}])".format(
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
		if self._material.find("steel") == 0:
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
	return L.inch(diameter)

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
	assert isinstance(start_point, Point)
	assert isinstance(end_point, Point)
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
	inch = L.inch
	big = inch(123456789.0)
	
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
			if depth <= L.inch(0):
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

    def simple_pocket(self, comment, \
      corner1_point, corner2_point, radius, flags):
	""" Part construct: Mill a pocket in {self} where {corner1_point}
	    and {corner2_point} specify a diagonal across the pocket.  The
	    radius of the inside corners is {radisu}.  {flags} can have
	    the character 't' for a pocket that goes {self} and 'u' to
	    specify that an upper chamfer is requested.  {comment} will
	    show up in error messages and any generated G-code. """

	# Perform argument type checking:
	assert isinstance(comment, str)
	assert isinstance(corner1_point, Point)
	assert isinstance(corner2_point, Point)
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
	""" Point dimensions: Store Set the value of {string_name} in {self}
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
		value = Point(self, zero, zero, zero)
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
	    flavor_name = "Point"
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

	# Before we get too far, flush any pending screw holes
	# from the previous mounting:
	if self._top_surface_set:
            self.screw_holes(True)

	# Extract some values from {ezcad}:
	ezcad = self._ezcad
	xml_indent = ezcad._xml_indent
	xml_stream = ezcad._xml_stream

	corner1 = self.tne
	corner2 = self.bsw
	cx = (corner1.x + corner2.x).half()
	cy = (corner1.y + corner2.y).half()
	cz = (corner1.z + corner2.z).half()

	# <Vice_Position TX= TY= TZ= NX= NY= NZ= WX= WY= WZ= CX= CY= CZ=/>:
	if xml_stream != None:
	    xml_stream.write( \
	      '{0}<Vice_Position TX="{1}" TY="{2}" TZ="{3}"'. \
	      format(" " * xml_indent, surface_point.x, surface_point.y,
	      surface_point.z))
	    xml_stream.write( \
	      ' NX="{0}" NY="{1}" NZ="{2}"'.format( \
	      north_point.x, north_point.y, north_point.z))
	    xml_stream.write(' WX="{0}" WY="{1}" WZ="{2}"'.format( \
	      west_point.x, west_point.y, west_point.z))
	    xml_stream.write(' CX="{0}" CY="{1}" CZ="{2}"'.format(cx, cy, cz))
	    xml_stream.write(' Comment="{0}"/>\n'.format(comment))
	self._top_surface = surface_point
	self._top_surface_set = True

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

    def __init__(self, part = None, name = None,
      center = None, axis = None, rotate = None, translate = None):
	""" Place: Initialize {self} to contain {home_part}, {placed_part},
	    {center_point}, {axis_point}, {rotate_angle}, and
	    {translate_point}. """

	# Deal with default argument values:
	none_type = type(None)
	if type(name) == none_type:
	    assert type(part) != none_type
	    name = part._name
	if type(center) == none_type:
	    center = Point()
	if type(axis) == none_type:
	    axis = Point(z = L(1.0))
	if type(rotate) == none_type:
	    rotate = Angle()
	if type(translate) == none_type:
	    translate = Point()

	part_name = "<none>"
	if type(part) != none_type:
	    part_name = part._name

	# Check argument types:
	assert type(part) == none_type or isinstance(part, Part)
	assert isinstance(name, str)
	assert isinstance(center, Point)
	assert isinstance(axis, Point)
	assert isinstance(rotate, Angle)
	assert isinstance(translate, Point)

	#print (("Place.__init__(part={0}, name={1}, center={2}, " + \
	#  "axis={3}, rotate={4}, translate={5})"). \
	#  format(part_name, name, center, axis, rotate, translate))

	# Extract some values from {rotate_center}:
	center_x = center.x.inches()
	center_y = center.y.inches()
	center_z = center.z.inches()

	# Extract some values from {axis_point}:
	axis_x = axis.x.inches()
	axis_y = axis.y.inches()
	axis_z = axis.z.inches()
	#print "axis=({0},{1},{2})".format(axis_x, axis_y, axis_z)

	# Normalize rotate axis:
	axis_length = \
	  math.sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z)
	#print "axis_length=", axis_length
	x = axis_x / axis_length
	y = axis_y / axis_length
	z = axis_z / axis_length

	rotate_matrix = Matrix.rotate_create(x, y, z, rotate)

	# Create {center_matrix} and {center_reverse_matrix}:
	center_matrix = Matrix.translate_create(-center_x, -center_y, -center_z)
	center_reverse_matrix = \
	  Matrix.translate_create(center_x, center_y, center_z)

	# Extract some values from {translate_point}:
	translate_x = translate.x.inches()
	translate_y = translate.y.inches()
	translate_z = translate.z.inches()

	# Create {translate_matrix}:
	translate_matrix = \
	  Matrix.translate_create(translate_x, translate_y, translate_z)
	#print "translate_matrix=\n", translate_matrix

	#print "Place(): cm=\n{0}\nrm=\n{1}\nrcm=\n{2}\ntm=\n{3}\n". \
	#  format(center_matrix.mat, rotate_matrix.mat, \
	#  center_reverse_matrix.mat, translate_matrix.mat)

	# Compute {forward_matrix} and {reverse_matrix}:
	forward_matrix = center_matrix * \
	  rotate_matrix * center_reverse_matrix * translate_matrix
	reverse_matrix = forward_matrix.inverse()
	#print "Place(): forward=\n{0}\nreverse=\n{0}\n". \
	#  format(forward_matrix.mat, reverse_matrix.mat)

	# Load up *self*:
	self._axis = axis
	self._center = center
	self._forward_matrix = forward_matrix
	self._part = part
	self._name = name
	self._reverse_matrix = reverse_matrix
	self._rotate = rotate
	self._translate = translate

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

class Point:
    """ {Point} represents a point in 3-space associated with a specific
	{Part} to provide the frame of reference. """

    def __init__(self, x = None, y = None, z = None):
	""" Point: Intialize {self} to contain {part}, {x}, {y}, {z}. """

	# Deal with default arguments:
	if type(x) == type(None):
	    x = L(0.0)
	if type(y) == type(None):
	    y = L(0.0)
	if type(z) == type(None):
	    z = L(0.0)

	# Check argument types:
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(z, L)

	# Load up *self*:
	self.x = x
	self.y = y
	self.z = z

    def __add__(self, point):
	""" Point: Add {point} to {self}. """

	assert isinstance(point, Point)
	return Point(self.x + point.x, self.y + point.y, self.z + point.z)

    def __div__(self, number):
	""" Point: Return the result of dividing {self} by {number}. """

	number = float(number)
	assert isinstance(number, int) or isinstance(number, float)
	return Point(self.x / number, self.y / number, self.z / number)

    def __eq__(self, point):
	""" Point: Return {True} if {self} is equal to {point}. """

	assert isinstance(point, Point)
	return self.x == point.x and self.y == point.y and self.z == point.z

    def __format__(self, format):
	""" *Point*: Return *self* formatted as a string. """

	assert isinstance(format, str)
	if format != "":
            format = ":" + format
	format_string = \
	    "[{0" + format + "}, {1" + format + "}, {2" + format + "}]"

	return format_string.format(self.x, self.y, self.z)

    def __mul__(self, scalar):
	""" Point: Return the result of muliplying {self} by {scalar}. """

	return Point(self.part, \
	  self.x * scalar, self.y * scalar, self.z * scalar)

    def __rmul__(self, scalar):
	""" Point: Return the result of muliplying {self} by {scalar}. """

	return Point(self.part, \
	  self.x * scalar, self.y * scalar, self.z * scalar)

    def __ne__(self, point):
	""" Point: Return {True} if {self} is not equal to {point}. """

	assert isinstance(point, Point)
	return self.x != point.x or self.y != point.y or self.z != point.z

    def __neg__(self):
	""" Point: Return the negative of {self}. """

	return Point(self.part, -self.x, -self.y, -self.z)

    def __str__(self):
	""" Point: Return {self} as a formatted string. """

	return "({0}, {1}, {2})".format(self.x, self.y, self.z)

    def __sub__(self, point):
	""" Point: Subtract {point} from {self}. """

	assert isinstance(point, Point)
	return Point(self.x - point.x, self.y - point.y, self.z - point.z)

    def angle_between(self, point):
	""" Point dimensions: Return the angle between {self} and {point}. """

	# a . b = ||a|| ||b|| cos <AB		(1)
	# (a . b) / (||a|| ||b||) = cos <AB	(2)
	# acos [(a . b) / (||a|| ||b||)] = <AB	(3)

	x1 = self.x.inches()
	y1 = self.y.inches()
	z1 = self.z.inches()
	x2 = point.x.inches()
	y2 = point.y.inches()
	z2 = point.z.inches()
	dot_product = x1 * y1 + x2 * y2 + z1 * z2
	length1 = math.sqrt(x1 * x1 + y1 * y1 + z1 * z1)
	length2 = math.sqrt(x2 * x2 + y2 * y2 + z2 * z2)
	
	return Angle.rad(math.acos(dot_product / (length1 * length2)))

    def cross_product(self, point):
	""" Point dimensions: Return the cross product of {self}
	    with {point}. """

	ux = self.x.inches()
	uy = self.y.inches()
	uz = self.z.inches()
	vx = point.x.inches()
	vy = point.y.inches()
	vz = point.z.inches()

	assert self.part == point.part, \
	  "Point.cross_product: Mis-matched parts"

	inch = L.inch
	return Point(self.part, inch(uy * vz - uz * vy),
	  inch(uz * vx - ux * vz), inch(ux * vy - uy * vx))

    def distance(self, point):
	assert isinstance(point, Point)
	dx = (self.x - point.x).inches()
	dy = (self.y - point.y).inches()
	dz = (self.z - point.z).inches()
	length = Lenght.inch(math.sqrt(dx * dx + dy * dy + dz * dz))
	return length

    def half(self):
	""" Point dimensions: Return {self} / 2. """

	result = self / 2.0
	return result	

    def length(self):
	""" Point dimensions: Return the length of self. """

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

    def points(self, dx, dy, dz):
	""" Part construct: Return a list of points centered around
	    {self} that are separated by {dx} in X, {dy} in Y and {dx}
	    in Z.  If all of {dx}, {dy}, {dz} are non-zero, 8 {Point}'s are
	    returned.  If one of {dx}, {dy}, and {dz} is zero, 4 {Point}'s
	    are returned.  If two of {dx}, {dy}, and {dz} are zero, only
	    2 {Point}'s are returned. """

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
            half_dx = dx.half()
	    x_list = (x - half_dx, x + half_dx)

	if dy == zero:
	    y_list = (y)
	else:
            half_dy = dy.half()
	    y_list = (y - half_dy, y + half_dy)

	if dz == zero:
	    z_list = (z)
	else:
            half_dz = dz.half()
	    z_list = (z - half_dz, z + half_dz)

	# Now iterate over {x_list}, {y_list}, and {z_list}
	# to generate {result}:
	result = []
	for x in x_list:
	    for y in y_list:
		for z in z_list:
		    point = Point(part, x, y, z)
		    result.append(point)

	return result

    def twice(self):
	""" Point dimensions: Return {self} * 2. """

	return self * 2.0

    def x_adjust(self, x):
	""" Point dimensions: Return copy of {self} with {x} added to the
	    x field. """

	return Point(self.part, self.x + x, self.y, self.z)

    def xy_adjust(self, x, y):
	""" Point dimensions: Return copy of {self} with {x} added to the
	    x field and {y} added to the y field. """

	return Point(self.part, self.x + x, self.y + y, self.z)

    def xyz_adjust(self, x, y, z):
	""" Point dimensions: Return copy of {self} with {x} added to the
	    x field, {y} added to the y field, and {z} added to the z field. """

	return Point(self.part, self.x + x, self.y + y, self.z + z)

    def xz_adjust(self, x, z):
	""" Point dimensions: Return copy of {self} with {x} added to the
	    x field and {z} added to the z field. """

	return Point(self.part, self.x + x, self.y, self.z + z)

    def y_adjust(self, y):
	""" Point dimensions: Return copy of {self} with {y} added to the
	    y field. """

	return Point(self.part, self.x, self.y + y, self.z)

    def yz_adjust(self, y, z):
	""" Point dimensions: Return copy of {self} with {y} added to the
	    y field and {z} added to the z field. """

	return Point(self.part, self.x, self.y + y, self.z + z)

    def z_adjust(self, z):
	""" Point dimensions: Return copy of {self} with {z} added to the
	    z field. """

	return Point(self.part, self.x, self.y, self.z + z)

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

from numpy import matrix

class Material:
    def __init__(self, generic, specific):
	""" *Material*: Initialize *self* to contain *generic* and
	    *specific*. """

	self.generic = generic
	self.specific = specific

class Matrix:

    def __init__(self, values):
	""" Matrix public: Initialize {self} to values. """

	self.mat = matrix(values)

    def __eq__(self, m):
	""" Matrix public:  Return {True} if {self} equals {m}. """

	return (self.mat == m.mat).all()

    def __mul__(self, m):
	""" Angle: Return {self} multiplied by {m}. """

	result_mat = self.mat * m.mat
	result = Matrix([[0]])
	result.mat = result_mat
	return result

    @staticmethod
    def identity_create():
	""" Matrix public: Return an identity matrix. """

	result = Matrix([ \
	  [1, 0, 0, 0],   \
	  [0, 1, 0, 0],   \
	  [0, 0, 1, 0],   \
	  [0, 0, 0, 1]  ])
	return result

    def inverse(self):
	""" Matrix public: Return the inverse matrix for {self}. """

	result = Matrix([[0]])
	result.mat = self.mat.I
	return result

    def point_create(self):
	""" Matrix public:  Return the {Point} the corresponds to {self}
	    using {part} as the returned {Point}'s reference frame.  {self}
	    must be a 1 x 4 matrix. """

	mat = self.mat
	x = mat[0, 0]
	y = mat[0, 1]
	z = mat[0, 2]
	point = Point(L(mm=x), L(mm=y), L(mm=z))

	return point

    def point_multiply(self, point):
	""" Matrix public: Return the point that results from mulitiplying
	    {point} by {self} and converting it back into a {Point} with
	    {part} as the reference frame. """

	point_matrix = point.matrix_create()
	result_matrix = point_matrix * self
	result_point = result_matrix.point_create()
	return result_point

    @staticmethod
    def rotate_create(nx, ny, nz, angle):
	""" Matrix public: Return a rotation matrix for rotating around
	    the normalized vector ({nx}, {ny}, {nz}) by {angle}.  {angle}
	    must be of type {Angle}."""

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

	# Compute some sub expressions:
	c = angle.cosine()
	s = angle.sine()
	omc = 1.0 - c
	x_omc = nx * omc
	y_omc = ny * omc
	z_omc = nz * omc
	xs = nx * s
	ys = ny * s
	zs = nz * s
    
	# Create the matrix:
	matrix = Matrix([ \
	  [nx * x_omc + c,  nx * y_omc - zs, nx * z_omc + ys, 0.0], \
	  [ny * x_omc + zs, ny * y_omc + c,  ny * z_omc - xs, 0.0], \
	  [nz * x_omc - ys, nz * y_omc + xs, nz * z_omc + c,  0.0], \
	  [0.0,             0.0,             0.0,             1.0] ])

	return matrix

    @staticmethod
    def translate_create(dx, dy, dz):
	""" Matrix public: Return a translate matrix containing {dx}, {dy}
	    and {dz}. """

	# Create the matrix:
	matrix = Matrix([ \
	  [1,  0,  0,  0], \
	  [0,  1,  0,  0], \
	  [0,  0,  1,  0], \
	  [dx, dy, dz, 1] ])

	return matrix

