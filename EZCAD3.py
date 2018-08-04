####################################################################################################
#<-----------------------------------------100 characters----------------------------------------->#
#
# EZCAD Coding Standards/Style
#
# The overall purpose of the EZCAD Coding standards are to provide legible code that can
# actually be maintained.
#
# The standards are:
#
# * All comments are written in Markdown.  This means that references to variables and class names
#   are marked in *italics* inside all comments.
# * All variable names are in lower case with adverbs, nouns, and verbs, separated by underscore
#   `_` characters.  Fully spelled out words are **strongly** encouraged.  Dropping vowels and
#   and syllables is simply not done.
# * All class names are like variable names, but first letter of each word is capitalized.
#   CamelCase class names are not used.  The only exception is classes that come from imported
#   Python libraries.  The *L* (i.e. length) and *P* classes (i.e. point) classes are an
#   exception to the naming convention of using words because these classes are so heavily used.
# * No lines longer than 100 characters are allowed.  When the code gets too long, is manually
#   formatted to fit within 100 characters.  No exceptions!
# * Indentation levels are 4 spaces.  Continuation lines are typically typically indented to
#   some multiple of 4 with two more spaces.
# * Strings of multiple characters are enclosed in double quotes (`"`) and single characters are
#   enclosed in single quotes (`'`).
# * Classes are listed alpabetically in the code.  The *Angle*, *P* and *L* classes are pulled
#   forward because they sometimes are needed for default argument values in other routines.
# * Within a class, routine definitions are listed alphabetically.  A preceeding underscore `_`
#   is ignored for alphabetical sorting.  An underscore in the middle of the routine definition
#   is treated as a space.
# * Each routine definition is followed by a description of what the routine does in standard
#   Python triple quotes `""" .... """`.  If the routine definition is defined in a class,
#   the class name is listed first in the triple quotes (e.g. `""" *Class_Name*: ..."""`.)
#   This helps the user disambiguate overloaded routine definitions in Python.
# * For longer routines, the *self* variable gets assigined a more descriptive variable name:
#         # Use *class_name* instead of *self*:
#         class_name = self
# * Routine definition variable names are type checked upon routine entry with assertions:
#         assert isinstance(argument_1, Type_1)
#         assert isinstance(argument_2, Type_2)
#         ...
#         assert isinstance(argument_n, Type_N)
# * Routines almost always exit a the end of the routine definition.  `return` statments in
#   in the middle of a routine are strongly discouraged.
# * There is an extensive tracing system used through out the EZCAD code.  There is an optional
#   last argument on most routine definitions of the form `tracing = -1000000`.  Whenever, this
#   argument is positive it causes nested tracing of routine call/returns to occur.  In many
#   routines there is an internal additional variable called *trace_detail* that controls
#   how much additional tracing ocurrs.  *trace_detail* is manually set in the routine.  An
#   example follows:
#
#        # Perform any requested *tracing*:
#        trace_detail = -1
#        if tracing >= 0:
#            indent = ' ' * tracing
#            print("{0}=>Class_Name.Routine_Name({1}, {2}, ...)").format(indent, arg_1, arg2, ...))
#            trace_detail = 2
#
#   There is another matching trace statement at the end of the routine definition.
# * Code within a routine is written in "paragraph style", where there is an
#   English comment that preceeds each chunk of statements.  It should be possible to
#   read the comments before reading the code to figure out generally what is going on.
# * No python global variables are used.  Period.
# * Doxygen was tried for while, but was found to be a lot of typing for not much value.
#   Doxygen is no longer used.
# * Class routine defintions that start with an `_` are "private/protected" and only meant
#   to be used interally by routines in EZCAD.  These routines can be renamed/removed
#   without breaking user code.
# * Class variable definitions are declared in the *__init__()* routine and almost all of
#   start with an underscore (`_`).  The preceeding `_` means "private/protected" and these
#   and their sub-classes may directly access them using "dot notation*".  Accessor functions
#   are provided to allow other classes to these private/protected class variable definitions.
#
####################################################################################################

# Mounting:
#
# This is an overview of the of the various issues associated with part mounting.
# Part mounting is actually a pretty involved topic.
#
# In EZCAD, a part is constructed in two phases:
#
# * Additive phase: The first phase is an additive phase where various solids are
#   "welded" together.  The solids can be blocks, rods (i.e. cylinders), extrusions
#   (i.e. L shaped, I shaped, channel shaped, hollow tubes, etc.)  At the end of
#   the first phase, EZCAD has computed a bounding box (technically it is a prism)
#   that contains all of the solids.  This bounding box/prism is oriented in 3D
#   space such that each surface is parallel to either the X/Y, X/Z, or Y/Z planes.
#
# * Subtractive phase: The second phase is a subtractive phase where material is removed.
#   The mental model is that the user can perform a variety of different opartions, where
#   each operation removes some material.  Example operations are:
#
#   * Exterior Contour.  An exterior contour removes a path of material on the outside
#     surface of a part.
#
#   * Rectangular Pocket:  A rectangular pocket removes a chunk material that is sort
#     rectangular in shape.
#
#   * Round Holes: A round hole is for removing a cylinder of material.
#
#   These operations are *ALWAYS* done in 3D space.  Thus, a hole is removed (i.e. drilled)
#   by specifying a start point an end point and a diameter.
#
# Neither the additive phase nor the subtractive phase actually cares about mounts.
# If the desired output is going to be manufactured by a 3D printer, no mounts are
# actually needed and the generated `.stl` file can be simply fed into the 3D printer
# tool chain.  However, if the part is going to be constructed using a mill or laser
# counter, mounting is very definitely something that needs to be performed by the part
# designer in EZCAD.
#
# A mount specifies a particular orientation of a part.  For a CNC mill, the part is
# is either mounted directly on the CNC table, in a vice that is mounted to the CNC table,
# or on a tooling plate that is mounted in a vice.  A tooling plate is a plate of material
# (usually aluminum) that has rectangular grid of holes cut in it.  Prior to mounting a
# part on a tooling plate, mounting holes have to drilled into the part that match the
# tooling plate rectangular grid.  For a CNC laser cutter, the part simply needs be placed
# on the laser cutter table table.  In short a mount specifies how the part is placed
# or attached to the CNC machine.
#
# In EZCAD, every part has a name.  Each time the user reorients a part for some operations,
# a mount must be specified.  Each mount has a unique name for the part.  Thus, each mount
# can be uniquely specified by a part name and a mount name for the part.
#
# The process of mounting on a tooling plate involves 3 steps:
#
# * Mount the part in a vice ( *vice_mount*(...) )
# * Drill tooling plate mounting holes into the part ( *tooling_plate_drill*(...) )
# * Mount the part on the tooling plate ( *tooling_plate_mount*(...) )
#
# In addition to all of this, there is also a technique call multi-mount that is used
# to combine the mounts of several parts together to make multiple parts at the same
# time using combined CNC paths.  The multi-mount stuff had a greater impact on EZCAD
# than was initially expected.  In order to make the multi-mounts work, they have to
# be queued up to occur *after* all of the regular mounts have occurred.  This is to
# ensure that tool paths for all *Part*'s have been determined.
#
# There are a number of Python classes that are used to support CNC path generation.
# These classes are listed below with a short description of what the class does:
#
# * *Part*: The part represents the part to manufactured.  It has a name, and a table
#   (i.e. Python *dict*) of named *Mount_Operations* objects.  The *Mount_Operations*
#   class is described further below.
#
# * *Operation*: The *Operation* super-class specifies all the stuff that is common to
#   to all operations.  It is sub-classed for specific operations (e.g. *Operation_Drill*,
#   *Operation_Contour*, *Operation_Rectangular_Pocket*, etc.)  Each *Operation* specifies
#   a *Part* and a *Tool* (and a bunch of other less important stuff.)  The *Operation*
#   class does *NOT* specify a *Mount*.  In order to support multi-mounts, the binding
#   between a *Mount* and an *Operation* occurs in a higher level data structure (see
#   *Mount_Operation* below.)
#
# * *Mount*: A mount species a *Part*, a top surface *Transform*, and a mount translate point
#   (i.e. *P* object).  The top surface transform takes the bounding box of the *Part* and
#   translates it so that one surface of the bounding box is centered under the machine
#   origin (0, 0, 0).   The bounding box is rotated so that one of its surfaces is facing
#   "north"; this is the surface that will be facing the north jaw of the vice.  The mount
#   translate point moves the part from machine origin into its final location for CNC tool
#   path generation.
#
# * *Multi_Mount*: A *Multi_Mount* object specifies *Part*, the name of a *Mount* associated
#   with the part, a rotation *Angle*, and an offset (dx, dy).  This is the fundamental
#   unit of the *Multi_Mounts* object.
#
# * *Multi_Mounts*: A *Multi_Mounts* object is basically just a list of *Multi_Mount* objects.
#   A *Multi_Mounts* object is used to generate a coherent CNC path for multiple parts as a
#   single set of combined CNC operations.  Multi-mounts are discussed much more below.
#   With *Multi_Mounts*, a mixture of parts can be machined together.  They can all be the
#   same part, at different locations in the material, or a mixuture of parts.  The user
#   is responsible for keeping the parts from overlapping each other in the material.
#
# * *Mount_Operation*: A *Mount_Operation* specifies a binding between a *Mount* and an
#   *Operation*.
#
# * *Mount_Operations*: A *Mount_Operations* object is a list of *Mount_Operation* objects,
#   where all of the *Mount* objects are compatible with one another.  The *Mount_Operations*
#   object is the fundamental unit of CNC path planning.  Within a *Mount_Operations* list,
#   the *Mount_Operation* objects can be rearranged (with some constraints) to produce more
#   efficient CNC tool paths.
#
# * *Plate*: A *Plate* specifies a bunch of *Multi_Mounts*.  It is used to figure out how
#   to cut the various *Multi_Mounts* chunks from a larger chunk of material.
#
# Now it is time to talk about extra material.  In order to support exterior contouring and
# facing (by the way, facing operations have not been implemented yet) there needs to be some
# extra material on the original part to remove.  There are two solutions to this problem,
# the user can create the original part with extra material to start with *OR* extra material
# can be managed separately.  The later approach has been adopted because the extra material
# does not show up in the bounding box dimensions for each part.  Extra material is specified
# at vice mount time and it is kept track of by EZCAD.  Thus, each mount object contains
# a couple of points that specify where the extra material is.  If fact, there actually
# two pairs of points.
#
# Multi_Mounts:
#
# Frankly, the data structures for supporting multi-mounts evolved over time and are
# probably not entirely correct.  Originally, everything was *Part* centeric.  CNC code
# was generated for a single part.  Eventually, the concept of making multiple parts
# out of a single piece of material was conceieved.  Thus, the concept of a *Multi_Mount*
# was created.  A *Multi_Mount* is basically a *Part* that has a different *Mount* object
# than what is used for doing the CNC for a single *Part*.  The *Multi_Mounts* object can
# be thought of as a super *Part* that makes multiple parts out of one piece of material.
# The code is ugly and hard to understand and should probably be ripped out and redone.

# Imported libraries:
import hashlib				# Used to create unique hash for `.scad` files
import math				# Primarily used for trigonmic functions
import numpy				# Used for 4x4 afine matrix transforms
import os				# Used for operating system functions
import subprocess			# Used to run `openscad` program.
import sys				# System stuff

# Some stand-alone routine definitions:

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
    """ *float*: Return -1, 0, or 1 depending upon whether *left* sorts before,
	at, or after *right*.
    """

    assert isinstance(left, float)
    assert isinstance(right, float)
    if left < right:
	return -1
    elif left == right:
	return 0
    return 1

def string_strip(string, characters):
    """ *str*: Remove all of *characters* from *string*.
    """

    # This code uses the *translate* function that has changed between Python 2 an 3:
    try:
	# Python 2.x:
	result = string.translate(None, characters)
    except TypeError:
	# Python 3.x:
	table = {ord(character): None for character in characters}
	result = string.translate(table)
    return result

def alphanumeric_only(string):
    """ *str*: Remove all characters that are not alphanumeric and underscore.
    """

    return string_strip(string, " !\"#$%&'()*+,-./:;<=>?[\\]^`~")

# The *L*, *P*, and *Angle* classes are listed first because Python does not
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
	assert isinstance(self._mm, float)

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
    # be an *int*, *float*, *L*, or *Speed*.  For a divisor that is an *int* or *float*,
    # the returned result is of type *L*.  For a divisor that is an *L* type,
    # the returned quotient result is of type *float*.  For a divisor that is of type *Speed*:
    # the returned value is of type *Time*:
    def __div__(self, divisor):
	""" *L*: Return {self} / {scalar}. """

	# Check argument types:
	assert isinstance(divisor, float) or \
	  isinstance(divisor, int) or isinstance(divisor, L) or isinstance(divisor, Speed)

	if isinstance(divisor, L):
	    assert divisor._mm != 0.0
	    result = self._mm / divisor._mm
	elif isinstance(divisor, Speed):
	    result = Time(sec = self._mm / divisor._mm_per_sec)
	else:
	    assert divisor != 0.0
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

	# Round small values to zero:
	if -.0000000001 <= value <= .0000000001:
	    value = 0.0

	# Now do the final format:
	if len(format) == 0:
	    result = "{0:.6f}".format(value)
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

	mm = self._mm
	assert isinstance(mm, float)
	return mm

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

	assert isinstance(scalar, float) or isinstance(scalar, int)
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

    def absolute(self):
        """ *P*: Return the the *P* object (i.e. *self*) with all positive coordinate values.
	"""

	return P(self.x.absolute(), self.y.absolute(), self.z.absolute())

    def angle_between(self, point):
	""" *P*: Return the angle between *self*  and *point* with respect to origin. """

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

    def minimum_maximum(self, point2):
	""" *P*: Return two *P* objects where the first returne *P* object has the smaller X/Y/Z
	    coordinates and the second *P* object has the larger X/Y/Z coordinates.
	"""

	# Use *point1* instead of *self*:
	point1 = self

	# Verify argument types:
	assert isinstance(point2, P)

	# Find the minimum and maximum X/Y/Z coordinates for *point1* and *point2*:
	x1, x2 = point1.x.minimum_maximum(point2.x)
	y1, y2 = point1.y.minimum_maximum(point2.y)
	z1, z2 = point1.z.minimum_maximum(point2.z)

	# Return the minimum and maximum points:
	return P(x1, y1, z1,), P(x2, y2, z2)

    def normalize(self):
	""" *P*: Return *self* normalized to have a length of 1. """

	x = self.x._mm
	y = self.y._mm
	z = self.z._mm
	length = math.sqrt(x * x + y * y + z * z)
	zero = L()
	if length == 0.0:
            normalized = self
	    if not EZCAD3.ezcad._dimensions_mode:
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
	""" P: Return a *P* object that corresponds to a point in the X/Y plane that is
	    *distance away* from the origin and is rotated by *angle* from the X axis.
	"""

	assert isinstance(angle, Angle)
	assert isinstance(distance, L)

	x = distance.cosine(angle)
	y = distance.sine(angle)
	zero = L()
	return P(x, y, zero)

    def triple(self):
        """ *P*: Return the *P* object (i.e. *self*) as an immutable tuple of 3 floats
	    represented in millimeters.
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

    def xy_rotate(self, center, angle):
	""" *Point*: Return a new *Point* that has been rotated around the Z access that
	    goes through *center* by *angle*:
	"""

	# Verify argument types:
	assert isinstance(center, P)
	assert isinstance(angle, Angle)

	# Grab *px*, *py*, and *pz*  from the *Point* object (i.e. *self):
	point = self
	px = point.x
	py = point.y
	pz = point.z

	# Grab *cx* and *cy* from the *center*:
	cx = center.x
	cy = center.y

	# Compute *s* (for sine) and *c* (for cosine).
	s = angle.sine()
	c = angle.cosine()

	# Translate (px, py) to "origin* (*cx*, *cy*):
	x = px - cx
	y = py - cy

	# Compute new rotated location (*nx*, *ny*):
	rotated_x = x * c - y * s
	rotated_y = x * s + y * c

	# Return the new *P* object:
	return P(rotated_x + cx, rotated_y + cy, pz)

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
	self._radians = radians
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
	return Angle(rad = self._radians + angle._radians)

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
	    result = Angle(rad = self._radians / float(divisor))
	elif isinstance(divisor, Angle):
	    result = self._radians/divisor._radians
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
	return self._radians == angle._radians

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

	value = self._radians
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
	return self._radians >= angle._radians

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
	return self._radians > angle._radians

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
	return self._radians <= angle._radians

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
	return self._radians < angle._radians

    ## @brief Return *self* &times; *multiplier*.
    #  @param self is the *Angle* object to multiply.
    #  @param multiplier is the number to mulitply *self* by.
    #  @returns *self* &times *multiplier*.
    #
    # <I>__mul__</I>() returns *self* (an *Angle* object) &times; *multiplier*.
    def __mul__(self, multiplier):
	""" Angle: Return *self* divided by *multiplier*. """

	assert isinstance(multiplier, float) or isinstance(multiplier, int)
	return Angle(rad = self._radians * multiplier)

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
	return self._radians != angle._radians

    ## @brief Returns -*self*.
    #  @param self is the *Angle* object to negate.
    #
    # <I>__neg__</I>() will return -*self*.
    def __neg__(self):
	""" Angle: Return negative of {self}. """

	return Angle(rad = -self._radians)

    ## @brief Return *self* &times; *multiplier*.
    #  @param self is the *Angle* object to multiply.
    #  @param multiplier is the number to mulitply *self* by.
    #  @returns *self* &times *multiplier*.
    #
    # <I>__rmul__</I>() returns *self* (an *Angle* object) &times; *multiplier*.
    def __rmul__(self, multiplier):
	""" Angle: Return {self} multiplied by {multiplier}. """

	assert isinstance(multiplier, float) or isinstance(multiplier, int)
	return Angle(rad = self._radians * multiplier)

    ## @brief Returns *self* converted to degrees.
    #  @param self is the *Angle* object to convert to degrees.
    #  @returns *self* converted to degrees.
    #
    # <I>__str__</I>() Returns *self* converted to degrees as a *str* object.
    def __str__(self):
	""" Angle: Return {self} as a string. """

	return str(self._radians * 180.0 / Angle.pi)

    ## @brief Returns *self* - *angle*.
    #  @param self is the *Angle* object to subtract from.
    #  @param angle is the *Angle* to subtract.
    #  @returns *self* - *angle*.
    #
    # <I>__sub__</I>() returns *self* - *angle*.
    def __sub__(self, angle):
	""" Angle: Return {angle} subtracted from {self}. """

	assert isinstance(angle, Angle)
	return Angle(rad = self._radians - angle._radians)

    def absolute(self):
	""" *Angle*: Return absolution value of *Angle* object(i.e. *self*).
	"""

	result = self
	if self._radians < 0.0:
	    result = -self
	return result

    ## @brief Returns the cosine of *self*.
    #  @param self is the *Angle* to compute the cosine of
    #  @returns the cosine of *self*.
    #
    #  *cosine*() returns the cosine of *self* (an *Angle* object.)
    def cosine(self):
	""" Angle: Return the cosine of {self}. """

	return math.cos(self._radians)

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

	return self._radians * 180.0 / math.pi

    ## @brief Returns *self*/2.
    #  @param self is the *Angle* object to divide by 2.
    #  @returns *self/2.
    #
    # *half*() returns the *self*/2.
    def half(self):
	""" Angle: Return half of {self}. """

	return Angle(self._radians / 2.0)

    ## @brief Return a normalized value between -&pi; and &pi;.
    #  @param self is the *Angle* to normalize.
    #  @returns a normalized value for *self*.
    #
    # *normalize*() will return a normalized angle between -&pi; and &pi;.
    def normalize(self):
	""" Angle: Return a normalized angle. """

	radians = self._radians
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
	if angle2._radians < result._radians:
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

	return self._radians

    ## @brief Returns the sine of *self*.
    #  @param self is the *Angle* to compute the sine of
    #  @returns the sine of *self*.
    #
    #  *sine*() returns the sine of *self* (an *Angle* object.)
    def sine(self):
	""" Angle: Return the sine of {self}. """

	return math.sin(self._radians)

    ## @brief Returns the tangent of *self*.
    #  @param self is the *Angle* to compute the tangent of
    #  @returns the tangent of *self*.
    #
    #  *tangent*() returns the tangent of *self* (an *Angle* object.)
    def tangent(self):
	""" Angle: Return the tangent of {self}. """

	return math.tan(self._radians)

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
	# *point* after it has been rotated and translated so that the points transform down
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
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Bend._radius_center_and_tangents_compute('{1}', '{2}', '{3}')".format(
	    ' ' * tracing, bend._name, incoming_bend._name, outgoing_bend._name))
	    trace_detail = 1

	# Below is an ASCII art picture of a *Bend* object (i.e. *self*).  B represents
	# the bend point.  There is circle of radius r centered around the circle center
	# point C.  r is the radius of the bend.  *I* and *O* are the adjecent *Bend* object
	# points.  The circle tangentally touches the the IB line segment at J and the OB
	# line segement at N.  Thus the length of both JC and NC is r.  The length of BC
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
	# Our goal is to compute C using I, B, O and r.  Once we have that, we need will also
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
	if trace_detail >= 1:
	    indent = ' ' * tracing
	    print("{0}i={1:i} b={2:i} o={3:i}".format(indent, i, b, o))
	    print("{0}ib={1:i} ob={2:i}".format(indent, ib, ob))

	# Compute the angles of IB and OB from the origin:
	ib_angle = ib.xy_angle()
	ob_angle = ob.xy_angle()
	if trace_detail >= 1:
	    print("{0}ib_angle={1:d} ob_angle={2:d}".format(indent, ib_angle, ob_angle))

	# Compute *ibo_angle* and *cbo_angle*:
	ibo_angle = ib.angle_between(ob)
	assert ibo_angle != Angle(), \
	  "Bend='{0}' i={1:i} b={2:i} o={3:i} ib={4:i} ob={5:i}".format(bend._name, i, b, o, ib, ob)
	cbo_angle = ibo_angle / 2
	if trace_detail >= 2:
	    print("{0}ibo_angle={1:d} cbo_angle={2:d}".format(indent, ibo_angle, cbo_angle))
	    
	# Compute *d* using *cbo_angle* and *r*:
	r = bend._radius
	d = r / cbo_angle.sine()
	d_mm = d.millimeters()
	if trace_detail >= 2:
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
	if trace_detail >= 3:
	    print("{0}<<ib>>={1:m} <<ob>={2:m} <<cb>>={3:m}". \
	      format(indent, ib_direction, ob_direction, cb_direction))
	if trace_detail >= 2:
	    print("{0}cb_angle={1:d}".format(indent, cb_angle))
	if trace_detail >= 1:
	    print("{0}c={1:i}".format(indent, c))

	# The next step is remarkably simple.  The direction vector of <<JC>> from C to J
	# is just 90 degress from the *ib_angle* because the circle is tangent at J.
	# Likewsise, for direction vector <<NC>>.
	degrees90 = Angle(deg=90.0)
	jc_angle = (ib_angle - degrees90).normalize()
	jc_direction = jc_angle.xy_direction()
	nc_angle = (ob_angle + degrees90).normalize()
	nc_direction = nc_angle.xy_direction()
	if trace_detail >= 2:
	    print("{0}jc_angle={1:d} nc_angle={2:d}".format(indent, jc_angle, nc_angle))
	if trace_detail >= 3:
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
	    if trace_detail >= 1:
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


class BOM:
    """ *BOM*: A *BOM* represents Bill Of Materials.  It is a catch all for summary
	information about a *Part*:
    """

    def __init__(self, part):
	""" *BOM*: Initialize the *BOM* for *part*.
	"""
    
	# Use *bom* instead of *self*:
	bom = self

	# Verify argument types:
	assert isinstance(part, Part)

	# Stuff *part* into *bom* and initialize tables:
	bom._part = part
	bom._blocks = []
	bom._pending_block = None


    def _block_pending(self, part, material, corner1, corner2):
	""" *BOM*: Register a block of *material* with opposing corners of *corner1* and *corner2*
	    with the *BOM* object.
	"""

	# Use *bom* instead of *self*:
	bom = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(material, Material)
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)

	# Hang onto a *pending_block*:
	bom._pending_block = BOM_Block(part, material, corner1, corner2)

    def _block_register(self, transform):
	""" *BOM*: Register any pending *BOM_Block* associated with *BOM*.
	"""

	# Use *bom* instead of *self*:
	bom = self

	# Verify argument types:
	assert isinstance(transform, Transform)

	# Register *pending_bom_block* using to *transform* to reorient it:
	pending_block = bom._pending_block
	if isinstance(pending_block, BOM_Block):
	    pending_block._transform(transform)
	    bom._blocks.append(pending_block)
	    bom._pending_block = None

    def _merge(self, child_bom):
	""" *BOM*: Merge the contents of *child_bom* into to the parent *BOM* object (i.e. *self*.)
	"""

	# Use *parent_bom* instead of *self*:
	parent_bom = self

	# Merge the *BOM_Blocks* from *child_bom* to *parent_bom*:
	parent_bom._blocks.extend(child_bom._blocks)

    def _material_summary(self):
	""" *BOM*: Generate a material summary:
	"""

	# Use *bom* instead of *self*:
	bom = self

	blocks = bom._blocks
	blocks.sort(key=BOM_Block.key)

	thicknesses = {}
	for block in blocks:
	    key = block._key
	    thickness_key = (key[0], key[1], round(key[2], 4))
	    if thickness_key in thicknesses:
		thicknesses[thickness_key].append(block)
	    else:
		thicknesses[thickness_key] = [block]
	
	summary_file = open("/tmp/blocks.txt", "w")
	thicknesses_keys = sorted(thicknesses.keys())
	for thickness_key in thicknesses_keys:
	    grouped_blocks = thicknesses[thickness_key]
	    block0 = grouped_blocks[0]
	    adjusted_dz = block0._adjusted_dz
	    summary_file.write("{0} {1} {2}:\n".
	      format(thickness_key[0], thickness_key[1], adjusted_dz))
	    for block in grouped_blocks:
		dx = block._dx
		dy = block._dy
		dz = block._dz
		part_name = block._part.name_get()
		summary_file.write("    {0:i} x {1:i} x {2:i}: '{3}\n".
		  format(dx, dy, dz, part_name))
	summary_file.close()
	
class BOM_Block:
    """ *BOM_Block*: Represents a block of material.
    """

    # Desired thicknesses in inches:
    HDPE_THICKNESSES = (
      0.0625,
      0.1250,
      0.2500,
      0.3750,
      0.5000,
      0.6250,
      0.7500,
      1.0000,
      1.2500,
      2.0000,
      3.0000,
      4.0000 )
    ALUMINUM_THICKNESSES = (
      0.1250,
      0.2500,
      0.5000,
      1.0000 )

    def __init__(self, part, material, corner1, corner2):
	""" *BOM_Block*: Initialize the *BOM_Block* object to contain *material*, *corner1*, and
	    *corner2*.
	"""

	# Use *bom_block* instead of *self*:
	bom_block = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(material, Material)
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)

	# Load up *bom_block*:
	bom_block._part = part
	bom_block._material = material
	bom_block._corner1 = corner1
	bom_block._corner2 = corner2

    def _transform(self, transform):
	""" *BOM_Block*: Transform the corners of the *BOM_Block* object (i.e. *self*) using
	    *transform*.
	"""

	# Use *bom_block* instead of self
	bom_block = self

	# Verify argument types:
	assert isinstance(transform, Transform)

	# Transform the two corners
	corner1 = transform * bom_block._corner1
	corner2 = transform * bom_block._corner2

	# Compute the volume values:
	dx = (corner2.x - corner1.x).absolute()
	dy = (corner2.y - corner1.y).absolute()
	dz = (corner2.z - corner1.z).absolute()

	# Compute the *key* tuple used for sorting:
	# Determine *adjusted_dz* from the available *material* thicknesses:
	material = bom_block._material
	generic = material.generic_get()
	specific = material.specific_get()
	thicknesses = tuple()
	if specific == "hdpe":
	    # HDPE: High Density Polyethelene:
	    thicknesses = BOM_Block.HDPE_THICKNESSES
	elif generic == "aluminum":
	    # Aluminum:
	    thicknesses = BOM_Block.ALUMINUM_THICKNESSES

	# Grab the first thickness that is thick enough for *adjusted_dz*:
	adjusted_dz_inch = round(dz.inches(), 4)
	for thickness in thicknesses:
	    if thickness >= adjusted_dz_inch:
		adjusted_dz_inch = thickness
		break

	# Compute the *key* tuple:
	dx_inch = dx.inches()
	dy_inch = dy.inches()
	dz_inch = dz.inches()
	area_inch = dx_inch * dy_inch
	part_name = bom_block._part.name_get()
	key = (generic, specific, adjusted_dz_inch, dz_inch, area_inch, dx_inch, dy_inch, part_name)

	# Load up *bom_block*:
	bom_block._adjusted_dz = adjusted_dz_inch
	bom_block._dx = dx
	bom_block._dy = dy
	bom_block._dz = dz
	bom_block._key = key

    def key(self):
	""" *BOM_Block*: Return an immutable tuple for comparing the *BOM_Block* object
	    (i.e. *self*) with other *BOM_Block* objects when sorting.
	"""

	return self._key

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
	elif format == "rpm":
	    result = "{0}".format(frequency * 60.0)
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

    def generic_get(self):
	""" *Material*: Return the generic material name for the *Material* object (i.e. *self*):
	"""

	return self._generic

    def _needs_coolant(self):
	""" *Material*: Return *True* if the *Material* object (i.e. *self*) should
	    be machined with coolant on and *False* otherwise.
	"""

	return self._generic != "plastic"

    def _generic_get(self):
	""" *Material*: Return generic name for *Material* object (i.e. *self*).
	"""

	return self._generic

    def _is_plastic(self):
        """ *Material*: Return *True* if the *Material* object (i.e. *self*) is plastic
	    and *False* otherwise.
	"""

	return self._generic.lower() == "plastic"

    def _is_steel(self):
        """ *Material*: Return *True* if the *Material* object (i.e. *self*) is steel
	    (or iron) and *False* otherwise.
	"""

	return self._generic.lower() == "steel"

    def specific_get(self):
	""" *Material*: Return the specific material name for the *Material* object (i.e. *self*):
	"""

	return self._specific

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
	""" *Bounding_Box*: Return the top/north/east corner of the *Bounding_Box* object
	    (i.e. *self*) as a point (i.e. *P*) object.
	"""

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

    def point_expand(self, point, tracing=-1000000):
	""" *Bounding_Box*: Expand the *Bounding_Box* object (i.e. *self*)
	    to enclose *point*.
	"""

	# Use *bounding_box* instead of *self*:
	bounding_box = self

	# Verify argument types:
	assert isinstance(point, P)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Bounding_Box.point_expand({1:i}, {2:i})".
	      format(indent, bounding_box, point))

	# Extract *x*/*y*/*z* from *point*:
	x = point.x
	y = point.y
	z = point.z

	# Deal with an empty bounding box first:
	if bounding_box._is_empty:
	    # The *bounding_box* is empty, so make it exactly enclose *point*:
	    bounding_box._east = x
	    bounding_box._west = x
	    bounding_box._north = y
	    bounding_box._south = y
	    bounding_box._top = z
	    bounding_box._bottom = z
	else:
	    # Update the minimum and maximum *x*:
	    bounding_box._east =   x.maximum(bounding_box._east)
	    bounding_box._west =   x.minimum(bounding_box._west)
	    bounding_box._north =  y.maximum(bounding_box._north)
	    bounding_box._south =  y.minimum(bounding_box._south)
	    bounding_box._top =    z.maximum(bounding_box._top)
	    bounding_box._bottom = z.minimum(bounding_box._bottom)

	# Make sure the *bounding_box* is not broken:
	assert bounding_box._west <= bounding_box._east
	assert bounding_box._south <= bounding_box._north
	assert bounding_box._bottom <= bounding_box._top

	# Sure that we mark everything as not empty:
	bounding_box._is_empty = False

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Bounding_Box.point_expand({1:i}, {2:i}".
	      format(indent, bounding_box, point))

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

	# Verify argument types:
	assert isinstance(mm_per_sec, float) or isinstance(mm_per_sec, int)
	assert isinstance(in_per_sec, float) or isinstance(in_per_sec, int)
	assert isinstance(cm_per_sec, float) or isinstance(cm_per_sec, int)
	assert isinstance(ft_per_min, float) or isinstance(ft_per_min, int)
	assert isinstance(m_per_sec, float)  or isinstance(m_per_sec, int)

	self._mm_per_sec =			\
	  mm_per_sec +				\
	  cm_per_sec * 10.0 +			\
	  in_per_sec * 25.4 +			\
	  m_per_sec * 1000.0 +			\
	  ft_per_min * 12.0 * 25.4 / 60.0

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
	    # Dividing by length (i.e. *L*) yields a frequency (i.e. *Hertz*):
	    inches_per_second = self.inches_per_second()
	    inches = divisor.inches()
	    rps = inches_per_second / inches
	    hertz =  Hertz(rps=rps)
	    #print("inches_per_second={0}".format(inches_per_second))
	    #print("inches={0}".format(inches))
	    #print("rps={0}".format(rps))
	    #print("hertz={0:rps}rps (={0:rpm}rpm)".format(hertz))
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
	    scale = 60.0 / (25.4 * 12.0)
	    format_text = format_text[:-1]
	elif format_text.endswith('M'):
	    # Meters per second:
	    scale = 1 / 1000.0
	    format_text = format_text[:-1]
	elif format_text.endswith('c'):
	    # Centimeters per second:
	    scale = 1.0 / 10.0
	    format_text = format_text[:-1]
	elif format_text.endswith('i'):
	    # Inches per second:
	    scale =  1.0 / 25.4
	    format_text = format_text[:-1]
	elif format_text.endswith('I'):
	    # Inches per minute:
	    scale =  60.0 / 25.4
	    format_text = format_text[:-1]
	elif format_text.endswith('m'):
	    # Millimeters per second:
	    scale = 1.0
	    format_text = format_text[:-1]
	else:
            assert False, "Unrecognized Speed format '{0}'".format(format_text)
	
	# Format the value as a string and return it:
	value = self._mm_per_sec * scale
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

	return self._mm_per_sec / (25.4 * 12.0) * 60.0

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
	self._red = red
	self._green = green
	self._blue = blue
	self._alpha = alpha
	self._name = name

	#print("Color({0}, {1}, {2}, {3})".
	#  format(self.red, self.green, self.blue, self.alpha))

    def __format__(self, format):
	""" *Color*: Return formatted version of *self*. """

	if format == "s":
	    result = "[{0:.2f}, {1:.2f}, {2:.2f}, {3:.2f}]".format(
	      self._red, self._green, self._blue, self._alpha)
	else:
	    result = "[r={0:.2f}, g={1:.2f}, b={2:.2f}, a={3:.2f}]".format(
	      self._red, self._green, self._blue, self._alpha)
	return result

    def _alpha_get(self):
        """ *Color*: Return the alpha value of the *Color* object (i.e. *self*.)
	"""

	return self._alpha

    def _name_get(self):
	""" *Color*: Return the name associated with the *Color* object (i.e. *self*.)
	"""

	return self._name

    def _red_green_blue_get(self):
        """ *Color*: Return the red, green, and blue components of the *Color* object
	    (i.e. *self*).
	"""

	return self._red, self._green, self._blue

# Are these classes are unused???

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

    def _bends_reverse(self, tracing=-1000000):
	""" *Contour*: Reverse the directions of the bends in the *Contour* object (i.e. *self*.
	"""

	# Use *contour* instead of *self:
	contour = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Contour._bends_reverse(*)".format(indent))

	# Reverse the *bends*:
	contour._bends.reverse()

	# Toggle *is_clockswise_
	contour._is_clockwise = not contour._is_clockwise

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Contour._bends_reverse(*)".format(indent))

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

    def _inside_bends_identify(self, tracing=-1000000):
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
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Contour._inside_bends_identify('{1}')".format(indent, self._name))
	    trace_detail = 1

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
	    if trace_detail >= 3:
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
	if trace_detail >= 1:
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

	    if trace_detail >= 0:
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
	""" *Contour*: Sweep through the *Contour* object (i.e. *self*) and compute the projected
	    X/Y coordinates of *Bend* object on an X/Y plane using *transform* each point
	    in the contour prior to being projected down to the X/Y plane.
	"""

	# Verify argument types:
	assert isinstance(transform, Transform)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Contour._project('{1}', *)".format(indent, self._name))
	    trace_detail = 1

	# For each *bend* in bends, perform the project from 3D down to 2D:
	zero = L()
	if trace_detail >= 2:
	    print("{0}transform={1:v}".format(indent, transform))
	for index, bend in enumerate(self._bends):
	    # Grab the X/Y/Z coordinates for *point*:
	    name = bend._name_get()
	    point = bend._point_get()
	    transformed_point = transform * point
	    projected_point = P(transformed_point.x, transformed_point.y, zero)
	    bend._projected_point_set(projected_point)
	    if trace_detail >= 2:
		print("{0}Bend[{1}]: name='{2}' point={3:i} projected_point={4:i}".
		  format(indent, index, name, point, projected_point, transform))

	# Remember that we have performed the projection:
	self._is_projected = True

	# Hang onto the projection *transform* so we can eventually unproject:
	self._project_transform = transform        

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Contour._project('{1}', *)".format(indent, self._name,))

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
	trace_detail = -1
	if tracing >= 0:
	    trace_detail = 1
	    indent = ' ' * tracing
	    print("{0}=>Contour._scad_lines_polygon_append('{1}', *, '{2}', {3})".
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
	    #trim_extra = L(inch=2.0)
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
		    if trace_detail >= 1:
			print("{0}decrement".format(indent))
	    elif change_angle < degrees0:
		if contour_is_clockwise and bend_is_inside or \
		  not contour_is_clockwise and not bend_is_inside:
		    if tracing >= 1:
			print("{0}increment".format(indent))
		    change_angle += degrees360
	    if trace_detail >= 0:
		print("{0}name='{1}' start_angle={2:d} end_angle={3:d} change_angle={4:d}".
		  format(indent, bend_name, start_angle, end_angle, change_angle))

	    # Now we compute *deta_angle*:
	    arc_radius_mm = arc_radius.millimeters()
	    count = max(4, int(arc_radius_mm))
	    delta_angle = change_angle / float(count - 1)
	    if trace_detail >= 0:
		print("{0}count={1} delta_angle={2:d}".format(indent, count, delta_angle))
	    comma = ","
	    for index in range(count):
		angle = (start_angle + delta_angle * float(index)).normalize()
		point = center + angle.xy_direction() * arc_radius_mm
		if trace_detail >= 0:
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
	    print("{0}<=Contour._scad_lines_polygon_append('{1}', *, '{2}', {3})".
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

class Directory:
    """ *Directory*: A *Directory* object represents a sub-directory that is written to. """

    def __init__(self, parent_directory, sub_directory):
        """ *Directory*: Initialize the *Directory* object (i.e. *self*) to create the *name*
	    *sub_directory* in *parent_directory*.
	"""

	# Use *directory* instead of *self*:
	directory = self

	# Verify argument types:
	assert isinstance(parent_directory, str)
	assert isinstance(sub_directory, str)

	# Create the combined *path*:
	path = os.path.join(parent_directory, sub_directory)

	# Make sure that *sub_directory* for *path* exists:
	try:
	    os.makedirs(path)
	except OSError as exception:
            assert exception.errno == os.errno.EEXIST

	# If the directory previously exists, keep track of all of the *previous_files*
	# so any left overs can be deleted:
	previous_files = dict.fromkeys(os.listdir(path), None)

	# Load up *directory*:
	directory._current_files = {}
	directory._parent_directory = parent_directory
	directory._path = path
	directory._previous_files = previous_files
	directory._sub_directory = sub_directory

    def _clean_up(self, tracing=-1000000):
        """ *Directory*: Clean up the *Directory* object (i.e. *self*) so it has no left
	    over files from previous.
	"""

	# Use *directory* instead of *self*:
	directory = self

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Directory.clean_up('{1}')".format(indent, directory._path))

	# Delete any files that are left over from a previous run:
	path = directory._path
	previous_files = directory._previous_files
	for file_name in previous_files.keys():
	    full_path = os.path.join(path, file_name)
	    os.remove(full_path)
	    del previous_files[file_name]
	    if tracing >= 0:
		print("{0}Removed file '{1}'".format(indent, full_path))

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Directory.clean_up('{1}')".format(indent, directory._path))

    def _exists(self, file_name, tracing=-1000000):
        """ *Directory*: Return *True* if *file_name* is in the *Directory* object (i.e. *self*)
	    and *False* otherwise.
	"""

	# Use *directory* instead of *self*:
	directory = self

	# Verify argument types:
	assert isinstance(file_name, str)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
            print("{0}=>Directory._exists('{1}')".format(indent, file_name))

	# Determine whether *file_name* *exists* or not:
	current_files = directory._current_files
	previous_files = directory._previous_files
	exists = False
	if file_name in current_files:
	    # We have already written out *file_name*:
	    exists = True
	elif file_name in previous_files:
	    # *file_name* exists from a prevoious program run; mark it as current:
            exists = True
            del previous_files[file_name]
	    current_files[file_name] = None

	# Wrap up any requested *tracing* and return *exists* resulst:
	if tracing >= 0:
	    indent = ' ' * tracing
            print("{0}<=Directory._exists('{1}')=>{2}".format(indent, file_name, exists))
	return exists

    def _full_path(self, file_name, tracing=-1000000):
	""" *Directory*: Return the full path for *file_name* using the *Directory* object
	    (i.e. *self*.)
	"""

	# Use *directory* instead of *self*:
	directory = self

	# Verify argument types:
	assert isinstance(file_name, str)
        assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Directory._full_path('{1}', '{2}')".
	      format(indent, directory._path, file_name))

	# Compute *full_path*:
	full_path = os.path.join(directory._path, file_name)

	# Wrap-up any requested *tracing* and return *full_path*:
	if tracing >= 0:
	    print("{0}<=Directory._full_path('{1}', '{2}')".
	      format(indent, directory._path, file_name))
	return full_path

    def _lines_read(self, file_name, tracing=-1000000):
	""" *Directory*: Read *file_name* from the *Directory* object (i.e. *self*) return
	    it as a list of lines.
	"""

	# Use *directory* instead of *self*:
	directory = self

	# Verify argument types:
	assert isinstance(file_name, str)
        assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Directory._lines_read('{1}', '{2}')".
	      format(indent, directory._path, file_name))

	# Open the *file_name* as *input_file*:
	path = directory._path
	full_path = os.path.join(path, file_name)
	input_file = open(full_path, "r")
	assert isinstance(input_file, file), "Unable to open file '{0}'".format(full_path)

	# Read in *input_lines* and close *input_file*:
	input_lines = input_file.readlines()
	input_file.close()

	# Remember that we read *file_name*:
	directory._current_files[file_name] = None

	# Wrap-up any requested *tracing* and return the *input_lines*:
	if tracing >= 0:
	    print("{0}<=Directory._lines_read('{1}', '{2}')".
	      format(indent, directory._path, file_name))
	return input_lines

    def _lines_write(self, file_name, lines, tracing=-1000000):
	""" *Directory*: Write *lines* out to the *Directory* object (i.e. *self*) as *file_name*.
	    The SHA1 hash signature for *lines* is returned.
	"""

	# Use *directory* instead of *self*:
	directory = self

	# Verify argument types:
	assert isinstance(file_name, str)
	assert isinstance(lines, list)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
            print("{0}=>Directory._lines_write('{1}', '{2}', *)".
 	      format(indent, directory._sub_directory, file_name))

	# Create file *content* and associated hash *signature*:
	content = '\n'.join(lines)
	signature = hashlib.sha1(content).hexdigest()

	# Write out the *content* out to *file_name*:
	with directory._write_open(file_name, tracing + 1) as output_file:
	    output_file.write(content)

	# Wrap up any requested *tracing* and return *signature*:
	if tracing >= 0:
	  print("{0}<=Directory._lines_write('{1}', '{2}', *)=>'{3}'".
	    format(indent, directory._sub_directory, file_name, signature))
	return signature

    def _path_get(self):
        """ *Directory*: Return the *Directory* object (i.e. *self*) path. """

	return self._path

    def _write_open(self, file_name, tracing=-1000000):
        """ *Directory*: Open *file* name in the *Directory* object and return the associated
	    Python *file* object.
	"""

	# Use *directory* instead of *self*:
	directory = self

	# Verify argument types:
	assert isinstance(file_name, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Directory.write_open('{1}', '{2}')".
	      format(indent, directory._path, file_name))

	# Open *full_path* for writing:
	path = directory._path
	full_path = os.path.join(path, file_name)
	opened_file = open(full_path, "w")
	assert isinstance(opened_file, file), "Unable to open file '{0}'".format(full_path)

	# Remember that we opened *file_name*
	directory._current_files[file_name] = None

	# Remember to not delete it when when we call *_cleanup*:
	previous_files = directory._previous_files
	if file_name in previous_files:
	    del previous_files[file_name]

	# Wrap up any requested *tracing* and return:
	if tracing >= 0:
	    print("{0}<=Directory.write_open('{1}', '{2}')".
	      format(indent, directory._path, file_name))
	return opened_file

class EZCAD3:
    """ EZCAD3 is the top level engine that executes the design. """

    ezcad = None

    def __init__(self, minor = None, adjust = L(), directory="."):
	""" *EZCAD*: Initialize the contents of {self} to contain
	    *major* and *minor* version numbers. """

	# FIXME: *adjust* should be moved over to the *Part* library!!!
	# It is used for increasing the size of things in STL mode.

	#print "EZCAD.__init__() called"

	# Verify argument types:
	assert minor == 0
	assert isinstance(adjust, L)
	assert isinstance(directory, str)

	# Create the *Directory* objects for each sub-directory:
	dxf_directory  = Directory(directory, "dxf")
	ngc_directory  = Directory(directory, "ngc")
	scad_directory = Directory(directory, "scad")
	stl_directory  = Directory(directory, "stl")
	wrl_directory  = Directory(directory, "wrl")

	# Load up {self}:
	self._adjust = adjust
	self._bases_table = {}
	self._cnc_mode = False
	self._dimensions_mode = False
	self._directory = directory
	self._dxf_directory = dxf_directory
	self._major = 3
	self._minor = minor
	self._mounts = {}
	self._multi_mounts_list = []
	self._multi_mounts_table = {}
	self._ngc_directory = ngc_directory
	self._parts_stack = []
	self._plates_table = {}
	self._scad_directory = scad_directory
	self._shop = Shop("Wayne's Shop")
	self._stl_directory = stl_directory
	self._stl_mode = False
	self._update_count = 0
	self._wrl_directory = wrl_directory
	self._xml_indent = 0

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

	# Write out standard routines:
	self._routines_write()

	# Save a global copy of *ezcad*:
	EZCAD3.ezcad = self

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

    def _mount_register(self, mount):
	""" *EZCAD3*: Register *mount* with the *EZCAD3* object (i.e. *self*):
	"""

	# Use *ezcad* instead of *self*:
	ezcad = self

	# Verify argument types:
	assert isinstance(mount, Mount)

	# Grab the *mounts* table, build a *key* and stuff *mount* into *mounts* using *key*:
	mounts = ezcad._mounts
	mount_name = mount._name_get()
	part = mount._part_get()
	part_name = part._name_get()
	key = part_name + ":" + mount_name
	assert not key in mounts, "Mount {0} has already been created".format(key)
	mounts[key] = mount

    def _multi_mounts_find(self, multi_mounts_name, cell_dx, cell_dy):
	""" *EZCAD3*: Find/create the *Multi_Mounts* object associated with the *multi_mounts_name*
	    and return it.
	"""

	# Use *ezcad* instead of *self*:
	ezcad = self

	# Verify argument types:
	assert isinstance(multi_mounts_name, str) and not ' ' in multi_mounts_name
	assert isinstance(cell_dx, int) and cell_dx >= 0
	assert isinstance(cell_dy, int) and cell_dy >= 0

	# Finde/create the *multi_mounts* object:
	multi_mounts_table = ezcad._multi_mounts_table
	if multi_mounts_name in multi_mounts_table:
	    multi_mounts = multi_mounts_table[multi_mounts_name]
	else:
	    multi_mounts = Multi_Mounts(multi_mounts_name, cell_dx, cell_dy)
	    multi_mounts_table[multi_mounts_name] = multi_mounts
	    #print("Create Multi_Mounts('{0}') id={1}".format(multi_mounts_name, id(multi_mounts)))

	# Return the result:
	return multi_mounts

    def _multi_mounts_list_get(self):
	""" *EZCAD3*: Return the list of *Multi_Mounts* objects from the *EZCAD3* object
	    (i.e. *self*.)
	"""

	return self._multi_mounts_list

    def _ngc_directory_get(self):
	""" *EZCAD3*: Return the directory to read/write NGC files from/into from the *EZCAD3*
	    object (i.e. *self*).
        """

	return self._ngc_directory

    def _parts_stack_pop(self):
	""" *EZCAD3: Pop the top most *Part* from the *Part* stack associated with the *EZCAD3*
	    object (i.e. *self*.)
	"""

	# Pop the top most part from *parts_stack*:
	self._parts_stack.pop()

    def _parts_stack_push(self, part):
	""" *EZCAD3: Push *part* onto the *Part* stack associated with the *EZCAD3* object
	    (i.e. *self*.)
	"""

	# Use *ezcad* instead of *self*:
	ezcad = self

	# Verify argument types:
	assert isinstance(part, Part)

	# Push *part* onto *parts_stack*:
	parts_stack = ezcad._parts_stack
	parts_stack.append(part)

    def _parts_stack_reset(self):
	""" *EZCAD3:* Empty the *Part* stack associated with the *EZCAD3* object (i.e. *self*.)
	"""

	del self._parts_stack[:]

    def _parts_stack_to_text(self):
	""" *EZCAD3*: Return the *Part* stack associated with the *EZCAD3* object (i.e. *self*)
	    as a text string.
	"""

	return " ".join([ "{0}".format(part._name_get()) for part in self._parts_stack ])

    def _plate_lookup(self, name):
	""" *EZCAD3*: Return the *Plate* object associated with *name*.  *None* is returned
	    if there is no such association.
	"""

	result = None
	plates_table = self._plates_table
	if name in plates_table:
	    result = plates_table[name]
	return result

    def _plate_register(self, plate):
	""" *EZCAD3*: Register *plate* with the *EZCAD3* object (i.e. *self*):
	"""
	
	# Use *ezcad* instead of *self*:
	ezcad = self

	# Stuff *plate* into *plates_table* checking for duplicates:
	plates_table = ezcad._plates_table
	plate_name = plate._name_get()
	assert not plate_name in plates_table, "Duplicate create of Plate '{0}'".format(plate_name)
	plates_table[plate_name] = plate

    def _plates_process(self):
	""" *EZCAD3*: Process any queued *Plate* objects associated with the *EZCAD3* object
	    (i.e. *self*.)
	"""

	# Use *ezcad* instead of *self*:
	ezcad = self

	# Process each *plate* in *plates_table*:
	plates_table = ezcad._plates_table
	for plate in plates_table.values():
	    plate._process()

    def _routines_write(self):
        """ *EZCAD3*: Write out the standard routines drill cycles and the like.
	"""

	# Use *ezcad* instead of *self*:
	ezcad = self

	# Grab *ngc_directory*:
	ngc_directory = ezcad._ngc_directory

	# This is a better way to organize the routines:
	lines = []
	lines.append("O19 sub")
	lines.append("  (This is a drill cycle routine that breaks chips by)")
	lines.append("  (doing repeated drill pecks.  The arguments are:)")
	lines.append("  (x y z_safe z_start z_stop z_step z_back feed speed)")
	lines.append("")
	lines.append("  (Put arguments into local variables:)")
	lines.append("  #<x> = #1             (X coordinate to drill at)")
	lines.append("  #<y> = #2             (Y coordinate to drill at)")
	lines.append("  #<z_safe> = #3        (Z location above tooling)")
	lines.append("  #<z_start> = #4       (Z location where drilling starts)")
	lines.append("  #<z_stop> = #5        (Z location where drilling starts)")
	lines.append("  #<z_step> = #6        (Amount to peck down by)")
	lines.append("  #<z_back> = #7        (Amount to retrace on each peck)")
	lines.append("  #<feed> = #8          (Downward drill feed speed)")
	lines.append("  #<speed> = #9         (Spindle speed)")
	lines.append("")
	lines.append("  (Move to just above the hole:)")
	lines.append("  G0 Z#<z_safe>                         (Make sure we are at z_safe)")
	lines.append("  G0 X#<x> Y#<y>                        (Go to [x,y])")
	lines.append("  G0 Z0.10000                           (Rapid down to 0.1)")
	lines.append("  G1 F#<feed> S#<speed> Z#<z_start>     (Feed to z_start)")
	lines.append("")
	lines.append("  (Peck away at hole until done:)")
	lines.append("  #<z> = #<z_start>")
	lines.append("  O901 while [#<z> gt #<z_stop>]")
	lines.append("    G0 Z[#<z> + #<z_back>]              (Rapid retract by z_back)")
	lines.append("    #<z> = [#<z> - #<z_step>]           (Compute new z)")
	lines.append("    O902 if [#<z> lt #<z_stop>]")
	lines.append("      #<z> = #<z_stop>                  (Do not go too deep)")
	lines.append("      O902 endif")
	lines.append("    G1 Z#<z>                            (Drill to new z)")
	lines.append("    O901 endwhile")
	lines.append("  G0 Z#<z_safe>                         (All done, retract to z_safe)")
	lines.append("  O19 endsub")
	lines.append("")
	lines.append("(An example call:)")
	lines.append("O19 call [0.0] [0.0] [0.0] [-0.525] [0.100] [0.010] [10.000] [5000]")
	lines.append("M2")
	ngc_directory._lines_write("19.ngc", lines)

	# Construct subroutine `81.ngc` as a list of *lines* and write out to *ngc_directory*:
	lines = []
	lines.append("( Canned subroutine for drilling holes )")
	lines.append("( o81 call       [#1] [#2] [#3] [#4]    [#5 ]    [#6]      [#7] )")
	lines.append("o81 sub        ( [F]  [X]  [Y]  [Z_TOP] [Z_SAFE] [Z_RETRACT] [Z_BOTTOM] )")
	lines.append("G0 Z[#4]       ( Make sure we are at Z_TOP )")
	lines.append("G0 X[#2] y[#3] ( Make sure we are at [X, Y] )")
	lines.append("G0 Z[#6]       ( Rapid down to Z_RAPID )")
	lines.append("G1 F[#1] Z[#7] ( Drill down to Z_BOTTOM )")
	lines.append("G0 Z[#4]       ( Retract to Z_TOP )")
	lines.append("o81 endsub     ( End of subroutine )")
	lines.append("o81 call [10] [0] [0] [.5] [.1] [-.5] ( Example call )")
	lines.append("m2             ( End of File )")
	ngc_directory._lines_write("81.ngc", lines)

	# Construct subroutine `83.ngc` as a list of *lines* and write out to *ngc_directory*:
	lines = []
	lines.append("( Canned subroutine for drilling holes with pecking )")
	lines.append("( o83 call [#1] [#2] [#3] [#4]    [#5]     [#6]        [#7]       [#8]     )")
	lines.append("o83 sub  ( [F]  [X]  [Y]  [Z_TOP] [Z_SAFE] [Z_RETRACT] [Z_BOTTOM] [Z_PECK] )")
	lines.append("                  ( Local variables [Z]: #9 )")
	lines.append("G0 Z[#4]               ( // Make sure we are at Z_TOP )")
	lines.append("G0 X[#2] y[#3]         ( // Make sure we are at [X, Y] )")
	lines.append("G0 Z[#5]               ( // Rapid to Z_SAFE )")
	lines.append("G1 F[#1] Z[#6]         ( // Feed to Z_RETRACT )")
	lines.append("#9 = #6                ( Z := Z_RETRACT )")
	lines.append("o831 do                ( do { )")
	lines.append("    #9 = [#9 - #8]     (     Z := Z - Z_PECK // Bump Z down by Z_PECK)")
	lines.append("    o832 if [#9 LT #7] (     if Z < Z_BOTTOM { )")
	lines.append("        #9 = #7        (         Z := Z_BOTTOM )")
	lines.append("    o832 endif         (     } )")
	lines.append("    G1 F[#1] Z[#9]     (     // Drill down to Z )")
	lines.append("    G0 Z[#6]           (     // Rapid Retract to Z_RETRACT )")
	lines.append("o831 while [#9 GT #7]  ( } while [#9 > #7] )")
	lines.append("G0 Z[#4]               ( // Rapid retract to Z_TOP )")
	lines.append("o83 endsub         ( // End of subroutine )")
	lines.append("o83 call [10] [0] [0] [.5] [.1] [-.5] [.2] ( Example call )")
	lines.append("m2                 ( // End of File )")
	ngc_directory._lines_write("83.ngc", lines)

    def _scad_directory_get(self):
	""" *EZCAD3*: Return the directory to read/write SCAD files from/into from the *EZCAD3*
	    object (i.e. *self*).
        """

	return self._scad_directory

    def _shop_get(self):
	""" *EZCAD3*: Return the *Shop* object associated with the *EZCAD3* object
	    (i.e. *self*.)
	"""

	shop = self._shop
	assert isinstance(shop, Shop)
	return shop

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

    def process(self, part, tracing = -1000000):
	""" *EZCAD3*: Perform all of the processing starting at *part*.
	"""

	print("=>EZCAD3.process()")

	# Use *ezcad* instead of *self*:
	ezcad = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>EZCAD3.preocess('{1}')".format(indent, part._name_get()))

	# Process *part*:
	part.process(self, tracing = tracing + 1)

	# Process any *plate* objects:
	ezcad._plates_process()

	# Clean up any directories:
	self._dxf_directory._clean_up(tracing = tracing + 1)
	self._ngc_directory._clean_up(tracing = tracing + 1)
	self._scad_directory._clean_up(tracing = tracing + 1)
	self._stl_directory._clean_up(tracing = tracing + 1)
	self._wrl_directory._clean_up(tracing = tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=EZCAD3.preocess('{1}')".format(indent, part._name_get()))

    @staticmethod
    def _update_count_get():
	update_count = 123456789
	try:
	    ezcad = EZCAD3.ezcad
	    update_count = ezcad._update_count
        except:
            assert False
	return update_count

class Mount:
    """ *Mount*: Represents a part mounted on a machine for machining operations.
    """

    def __init__(self, name, part, top_surface_transform, mount_translate_point,
      top_surface_safe_z, xy_rapid_safe_z, is_tooling_plate_mount, selected_parallel_height,
      tracing=-1000000):
	""" *Mount*: Initialize the *Mount* object (i.e. *self*) to contain *name*, *part*,
	    *top_surface_transform*, and *mount_translate_point*.  *name* is provides a
	    mount name for debugging purposes.  *part* is the *Part* object being mounted.
	    *top_surface_translate* is *Transform* object that moves the part from 3D space
	    to located such that entire part is located under the machine origin oriented
	    properly for machining.  The machine origin should just touch the top surface of
	    the part.  *mount_translate_point* is the final part translation to mount in a
	    vice or on a tooling plate.  *top_surface_safe_z* is the machine Z axis altitude
	    above which is safe to prerform rapid moves (i.e. G0) in the Z direction.
	    *xy_rapid_safe_z* is the machine Z axis alitude above which it is safe to perform
	    rapid moves in the X and Y direction.  *is_tooling_plate_mount* is *True* if the
	    part is mounted on a tooling plate and *False* otherwise. *selected_parallel_height*
	    is the height of the parallels to use.
	"""

	# Use *mount* instead of *self*:
	mount = self
        
	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(part, Part)
	assert isinstance(top_surface_transform, Transform)
	assert isinstance(mount_translate_point, P)
	assert isinstance(top_surface_safe_z, L)
	assert isinstance(xy_rapid_safe_z, L)
	assert isinstance(is_tooling_plate_mount, bool)
	assert isinstance(selected_parallel_height, L)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount.__init__('{1}', '{2}', *, {3:i}, {4:i}, {5:i}, {6}, {7:i})".
	      format(indent, name, part._name_get(), mount_translate_point, top_surface_safe_z,
	      xy_rapid_safe_z, is_tooling_plate_mount, selected_parallel_height))
	    trace_detail = 2

	# Create *cnc_transform*:
	if trace_detail >= 2:
	    print("{0}top_surface_transform={1:v} mount_translate_point={2:i}".
	      format(indent, top_surface_transform, mount_translate_point))
	cnc_transform = top_surface_transform.translate("Vice Mount",
	  mount_translate_point, tracing = tracing + 1)
	if trace_detail >= 2:
	    print("{0}top_surface_transform={1:v}".format(indent, top_surface_transform))
	    print("{0}cnc_transform={1:v}".format(indent, cnc_transform))

	# Load up *mount*:
	zero = L()
	one = L(inch=1.000)
	mount._cnc_transform = cnc_transform			# CNC transfrom part to mount
	mount._extra_start_bsw = P(zero, zero, zero)		# BSW extra corner at mount start
	mount._extra_start_tne = P(one, one, one)		# TNE extra corner at mount start
	mount._extra_stop_bsw = P(zero, zero, zero)		# BSW extra corner at mount end
	mount._extra_stop_tne = P(one, one, one)		# TNE extra corner at mount end
	mount._is_tooling_plate_mount = is_tooling_plate_mount	# *True*=>parts are on tooling plate
	mount._mount_translate_point = mount_translate_point	# Translate from CNC origin to mount
        mount._name = name					# Mount name
	mount._part = part					# Part that is mounted
	mount._program_number = -1				# Program num of top level .ngc file
	mount._selected_parallel_height = selected_parallel_height # Height of the parallel to use
	mount._spacers = []					# Tooling plate spacer locations
	mount._stl_vrml = None					# STL *VRML* object for mount
	mount._tooling_plate_holes = []				# Each tooling plate hole (x, y)
	mount._top_surface_safe_z = top_surface_safe_z		# Z altitude above which Z rapids OK
	mount._top_surface_transform = top_surface_transform	# Transfrom from part to CNC orgin
	mount._xy_rapid_safe_z = xy_rapid_safe_z		# Z altitude for XY rapids 
	
	# Make sure that *part*:*name* is only created once:
	ezcad = part._ezcad_get()
	ezcad._mount_register(mount)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Mount.__init__('{1}', '{2}', *, {3:i}, {4:i}, {5:i}, {6}, {7:i})".
	      format(indent, name, part._name_get(), mount_translate_point, top_surface_safe_z,
	      xy_rapid_safe_z, is_tooling_plate_mount, selected_parallel_height))

    def _cnc_transform_get(self):
	""" *Mount*: Return the CNC *Transform* object associated with the *Mount* object
	    (i.e. *self*).
	"""

	return self._cnc_transform

    def _extra_start_get(self, tracing=-1000000):
	""" *Mount*: Return the starting extra material corners for the *Mount* object
	    (i.e. *self*).  These two corners specify of BSW (Bottom South West) and
	    TNE (Top North East) corners of a bounding box that encloses the *Part*
	    associated with the *Mount* object.
	"""
	
	# Use *mount* instead of *self*:
	mount = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount._extra_start_get('{1}')".format(indent, mount._name))
	    trace_detail = 3

	if trace_detail >= 2:
	    print("{0}part='{1}'".format(indent, mount._part._name_get()))

	# Grab the values from *mount*:
	extra_start_bsw = mount._extra_start_bsw
	extra_start_tne = mount._extra_start_tne

	# Verify that they have been set:
	assert isinstance(extra_start_bsw, P)
	assert isinstance(extra_start_tne, P)

	# Wrap up any requested *tracing* and return both *final_bsw* and *final_tne*:
	if tracing >= 0:
	    print("{0}<=Mount._extra_start_get('{1}')=>{2:i},{3:i}".
	      format(indent, mount._name, extra_start_bsw, extra_start_tne))
	return extra_start_bsw, extra_start_tne

    def _extra_start_set(self, extra_start_bsw, extra_start_tne, tracing=-1000000):
	""" *Mount* Store *extra_start_bsw* and *extra_start_tn* into the *Mount* object
	    (i.e. *self*).
	"""

	# Use *mount* instead of *self*:
	mount = self
	
	# Verify argument types:
	assert isinstance(extra_start_bsw, P)
	assert isinstance(extra_start_tne, P)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount._extra_start_store('{1}', {2:i}, {3:i})".
	      format(indent, mount._name, extra_start_bsw, extra_start_tne))
	    trace_detail = 2

	# Compute the *final_bsw* and *final_tne*:
	final_bsw, final_tne = extra_start_bsw.minimum_maximum(extra_start_tne)
	if trace_detail >= 1:
	    print("{0}final_bsw={1:i} final_tne={2:i}".format(indent, final_bsw, final_tne))

	# Fetch the two corners from *mount*:
	mount._extra_start_bsw = final_bsw
	mount._extra_start_tne = final_tne

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}=>Mount._extra_start_store('{1}', {2:i}, {3:i})".
	      format(indent, mount._name, extra_start_bsw, extra_start_tne))

    def _extra_stop_get(self, tracing=-1000000):
	""" *Mount*: Return the stopping extra material corners for the *Mount* object,
	    (i.e. *self*)
	"""
	
	# Use *mount* instead of *self*:
	mount = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount._extra_stop_get('{1}')".format(indent, mount._name))

	# Fetch the two extra stop corners from *mount*:
	extra_stop_bsw = mount._extra_stop_bsw
	extra_stop_tne = mount._extra_stop_tne

	# Wrap up any requested *tracing* and return both *final_bsw* and *final_tne*:
	if tracing >= 0:
	    print("{0}<=Mount._extra_stop_get('{1}', *)=>{2:i},{3:i}".
	      format(indent, mount._name, extra_stop_bsw, extra_stop_tne))
	return extra_stop_bsw, extra_stop_tne

    def _extra_stop_set(self, extra_stop_bsw, extra_stop_tne, tracing=-1000000):
	""" *Mount* Store *extra_stop_bsw* and *extra_stop_tne* into the *Mount* object
	    (i.e. *self*.)  These values must bracke the bounding box of the associated
	    *Part* object.
	"""

	# Use *mount* instead of *self*:
	mount = self
	
	# Verify argument types:
	assert isinstance(extra_stop_bsw, P)
	assert isinstance(extra_stop_tne, P)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Mount._extra_stop_set('{1}', *, {2:i}, {3:i})".
	      format(indent, mount._name, extra_stop_bsw, extra_stop_tne))
	    trace_detail = 2

	# Compute the *final_bsw* and *final_tne*:
	final_bsw, final_tne = extra_stop_bsw.minimum_maximum(extra_stop_tne)
	if trace_detail >= 1:
	    print("{0}final_bsw={1:i} final_tne={2:i}".format(indent, final_bsw, final_tne))

	# Fetch the two corners from *mount*:
	mount._extra_stop_bsw = final_bsw
	mount._extra_stop_tne = final_tne

	# Verify that the two corners were previously set:
	assert isinstance(extra_stop_bsw, P)
	assert isinstance(extra_stop_tne, P)

	# Make sure that *final_bsw*/*final_tne* box subsumes the *part_bsw*/*part_tne* box:
	part = mount._part
	part_bsw = part.bsw
	part_tne = part.tne

	# Floating point errors add up over time, so compensate with *epsilon*.
	ezcad = part._ezcad_get()
	#if ezcad._cnc_mode:
	#FIXME: Need to figure out Z adjustments to *Part* mounting:
	if False:
	    epsilon = L(inch=.000000001)
	    assert final_bsw.x - epsilon <= part_bsw.x, \
	      "final_bsw={0:i} part_bsw={1:i}".format(final_bsw, part_bsw)
	    assert final_bsw.y - epsilon <= part_bsw.y, \
	      "final_bsw={0:i} part_bsw={1:i}".format(final_bsw, part_bsw)
	    assert final_bsw.z - epsilon <= part_bsw.z, \
	      "final_bsw={0:i} part_bsw={1:i}".format(final_bsw, part_bsw)
	    assert part_tne.x <= final_tne.x + epsilon, \
	      "part_tne={0:i} final_tne={1:i}".format(part_tne, final_tne)
	    assert part_tne.y <= final_tne.y + epsilon, \
	      "part_tne={0:i} final_tne={1:i}".format(part_tne, final_tne)
	    assert part_tne.z <= final_tne.z + epsilon, \
	      "part_tne={0:i} final_tne={1:i}".format(part_tne, final_tne)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}=>Mount._extra_stop_set('{1}', *, {2:i}, {3:i})".
	      format(indent, mount._name, extra_stop_bsw, extra_stop_tne))

    def _is_tooling_plate_mount_get(self):
	""" *Mount*: Return *True* if the *Mount* object (i.e. *self*) is being used to mount
	    onto a tooling plate and *False* otherwise.
	"""

	return self._is_tooling_plate_mount

    def _is_vice_mount_get(self):
	""" *Mount*: Return *True* if the *Mount* object (i.e. *self*) is being used to mount
	    directly onto a vice and *False* otherwise.
	"""

	return not self._is_tooling_plate_mount

    def _mount_translate_point_get(self):
        """ *Mount*: Return the mount translate point (i.e. *P*) object associated with the
	    *Mount* object (i.e. *self*.)
	"""

	return self._mount_translate_point

    def _name_get(self):
        """ *Mount*: Return the name string associated with the *Mount* object (i.e. *self*.)
	"""

	return self._name

    def _part_get(self):
        """ *Mount*: Return the *Part* object associated with the *Mount* object (i.e. *self*.)
	"""

	return self._part

    def _selected_parallel_height_get(self):
	""" *Mount*: Return the selected parallel height associated with the *Mount* object
	    (i.e. *self*.)
	"""

	return self._selected_parallel_height

    def _spacers_get(self):
	""" *Mount*: Return the spacer locations associated with the *Mount* object (i.e. *self*.)
	"""

	return self._spacers

    def _spacers_set(self, spacers):
	""" *Mount*: Set the spacer locations for the *Mount* object (i.e. *self*) to *spacers*.
	"""

	# Use *mount* instead of *self*:
	mount = self
        
	# Verify argument types:
        assert isinstance(spacers, list)
	for spacer in spacers:
	    assert len(spacer) == 4
	    assert isinstance(spacer[0], int)
	    assert isinstance(spacer[1], int)
	    assert isinstance(spacer[2], int)
	    assert isinstance(spacer[3], int)

	# Store *spacers* into *mount*:
	mount._spacers = spacers

    def _stl_vrml_get(self, tracing=-1000000):
	""" *Mount*: Return the *VRML* object that represents `.stl` file associated with the
	    *Mount* object (i.e. *self*).
	"""

	# Use *mount* instead of *self*:
	mount = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount._stl_vrml_get('{1}')".format(indent, mount._name))

	# See whether or not we have already computed *stl_vrml*:
	stl_vrml = mount._stl_vrml
	if stl_vrml == None:
	    # It appears we need to compute *stl_vrml*:
	    part = mount._part_get()
	    part_name = part._name_get()
	    part_color = part._color_get()	
	    cnc_transform = mount._cnc_transform
	    stl_vrml_name = "{0}_{1}".format(mount._name, part_name)
	    stl = part._stl_get(tracing = tracing + 1)
	    triangles = stl._triangles_get()
	    stl_vrml = VRML(stl_vrml_name)
	    comment = stl_vrml_name
	    stl_vrml._stl_write(comment, part_color, triangles, tracing = tracing + 1)
	    mount._stl_vrml = stl_vrml

	# Wrap up any requested *tracing* and return *stl_vrml*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Mount._stl_vrml_get('{1}')=>'{2}'".
	      format(indent, mount._name, stl_vrml._name_get()))
	return stl_vrml

    def _tooling_plate_holes_get(self):
	""" *Mount*: Return all of the tooling holes locations with the *Mount* object
	    (i.e. *self*.)
	"""

	return self._tooling_plate_holes

    def _tooling_plate_holes_set(self, tooling_plate_holes, tracing=-1000000):
	""" *Mount*: Set the the tooling holes locations with the *Mount* object
	    (i.e. *self*) to the *tooling_plate_holes* list (i.e. (x, y) tuples).
	"""

	# Use *mount* instead of *self*:
	mount = self

	# Verify argument types:
	assert isinstance(tooling_plate_holes, list)
	for tooling_plate_hole in tooling_plate_holes:
	    assert isinstance(tooling_plate_hole, tuple) and len(tooling_plate_hole) == 2
	    column = tooling_plate_hole[0]
	    row = tooling_plate_hole[1]
	    assert isinstance(column, int) and column >= 0
	    assert isinstance(row, int) and row >= 0

	# Perform any requested *tracing*:
	if tracing > 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount._tooling_plate_holes_set('{1}', {2})".
	      format(indent, mount._name, tooling_plate_holes))

	# Stuff *tooling_plate_holes* into *mount*:
	mount._tooling_plate_holes = tooling_plate_holes

	# Wrap up any requested *tracing*:
	if tracing > 0:
	    print("{0}<=Mount._tooling_plate_holes_set('{1}', {2})".
	      format(indent, mount._name, tooling_plate_holes))

    def _top_surface_transform_get(self):
        """ *Mount*: Return the top surface *Transform* object associated with the *Mount*
	    object (i.e. *self*.)
	"""

	return self._top_surface_transform

    def _top_surface_safe_z_get(self):
        """ *Mount*: Return the top surface safe Z value associated with the *Mount*
	    object (i.e. *self*.)
	"""

	return self._top_surface_safe_z

    def _xy_rapid_safe_z_get(self):
        """ *Mount*: Return the XY rapid safe Z value associated with the *Mount*
	    object (i.e. *self*.)
	"""

	return self._xy_rapid_safe_z

class Mount_Operation:
    """ *Mount_Operation*: A binding between a mount and an operation.
    """

    def __init__(self, name, mount, operation):
	""" *Mount_Operation*: Initialize the *Mount_Operation* object (i.e. *self*) to
	    contain *mount* and *operation*.  *name* is primarily used for debugging.
	"""

	# Use *mount_operation* instead of *self*:
	mount_operation = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(mount, Mount)
	assert isinstance(operation, Operation)

	# Load up *mount_operation*:
	mount_operation._mount = mount
	mount_operation._name = name
	mount_operation._operation = operation

    def _mount_get(self):
	""" *Mount_Operation*: Return the *Mount* object associated with the *Mount_Opeartion*
	    object (i.e. *self):
	"""

	return self._mount

    def _name_get(self):
	""" *Mount_Operation*: Return the name of the *Mount_Opeartion* object (i.e. *self):
	"""

	return self._name

    def _operation_get(self):
	""" *Mount_Operation*: Return the *Operation* object associated with the *Mount_Opeartion*
	    object (i.e. *self):
	"""

	return self._operation

class Mount_Operations:
    """ *Mount_Operations*: Provides a list of *Mount*/*Operation* pairs.
    """

    def __init__(self, name, tracing=-1000000):
	""" *Mount_Operations*: Initialize the *Mount_Operations* (i.e. *self*) to be an empty list.
	"""

	# Use *mount_operations* instead of *self*:
	mount_operations = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(tracing, int)
        
	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount_Operations.__init('{1}')".format(indent, name))

	# Load up *mount_operations*:
	mount_operations._name = name
	mount_operations._pairs = []
	mount_operations._multi_mounts = None

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Mount_Operations.__init('{1}')".format(indent, name))

    def _append(self, mount, operation, prepend=False, tracing=-100000):
	""" *Mount_Operations*: Append the pair of (*mount*, *operation*) to the *Mount_Operations*
	    object (i.e. *self*) list.  If the optional *prepend* argument is *True*, the pair
	    is prepended to the *Mount_Operations* object instead being appended.
	"""

	# Use *mount_operations* instead of *self*:
	mount_operations = self

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(operation, Operation)
	assert isinstance(prepend, bool)
	assert isinstance(tracing, int)

	# Perform any requesting *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount_Operations._append('{1}', '{2}', '{3}', {4})".format(indent,
	      mount_operations._name, mount._name_get(), operation._name_get(), prepend))

	# Grab the fields of *mount_operations*:
	pairs = mount_operations._pairs
	name = mount_operations._name

	# Create the *mount_operation*:
	pairs_size = len(pairs)
	mount_operation_name = "{0}[{1}]".format(name, pairs_size)
	mount_operation = Mount_Operation(mount_operation_name, mount, operation)

	# Append/prepend *mount_operation* to *pairs*:
	if prepend:
	    pairs.insert(0, mount_operation)
	else:
	    pairs.append(mount_operation)
	
	# Wrap up any requesting *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Mount_Operations._append('{1}', '{2}', '{3}', {4})".format(indent,
	      mount_operations._name, mount._name_get(), operation._name_get(), prepend))

    def _cnc_dowel_pin_optimize(self, tracing=-1000000):
	""" *Mount_Operations*: Sweep through the *Mount_Operations* object (i.e. *self*)
	    and replace the *Operation_Dowel_Pin* *Tool* with a better dowel pin.
	"""

	# Use *mount_operations* instead of *self*:
	mount_operations = self
	
	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount_Operations._cnc_dowel_pin_optimize".format(indent))
	    trace_detail = 2
	if trace_detail >= 2:
	    mount_operations._show("cnc_dowel_point_optimize", tracing = tracing + 1)

	# Grab some values from *mount_operations*:
	pairs = mount_operations._pairs

	# Search through *pairs* for a *Operation_Dowel_Pin* followed by an *Operation_Drill*:
	size = len(pairs)
	for index, pair in enumerate(pairs):
	    operation = pair._operation_get()
	    if isinstance(operation, Operation_Dowel_Pin):
		# We have a dowel pin, see if there is a next operation:
		if index + 1 < size:
		    # Grab the *next_operation* after *operation*:
		    next_pair = pairs[index + 1]
		    next_operation = next_pair._operation_get()
		    next_tool = next_operation._tool_get()

		    # Get the *shop* object so we can search it:
		    mount = pair._mount_get()
		    part = mount._part_get()
		    ezcad = part._ezcad_get()
		    shop = ezcad._shop_get()

		    # Search *shop* for *Tool_Dowel_Pin* that matches with *next_tool*:
		    dowel_pin_tool = shop._dowel_pin_search(next_tool)
		    if isinstance(dowel_pin_tool, Tool_Dowel_Pin):
			# Replace the dowel pin in *operation* with *dowel_pin_tool*:
			operation._tool_set(dowel_pin_tool)
			if trace_detail >= 1:
			    print("{0}Replacing dowel pin".format(indent))

			# We did the replacement, so we are done with the search:
			break

		# We only deal with the first dowel pin, so we can stop the search now:
                break

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Mount_Operations._cnc_dowel_pin_optimize".format(indent))

    def _cnc_mount_generate(self, mount_program_number, tracing=-1000000):
        """ *Mount_Operations*: Output the CNC G-code and associated VRML path visualization
	    files for the *Mount_Operations* object (i.e. *self*).  The first G-code file
	    output will start with *mount_program_number*.  The next usable mount program
	    number (which is divisible by 10) is returned.  This routine assumes that each
	    *Mount* object in the *Mount_Operations* object is a compatible with the others.
	"""

	# Use *mount_operations* instead of *self*:
	mount_operations = self

	# Verify argument types:
        assert isinstance(mount_program_number, int) and mount_program_number > 0
        assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount_Operations._cnc_mount_generate('{1}', {2}) id={3}".format(
	      indent, mount_operations._name, mount_program_number, id(mount_operations)))
	    trace_detail = 1

	# Stuff *mount_program_number* into *mount_operations*:
	multi_mounts = mount_operations._multi_mounts_get()
	if isinstance(multi_mounts, Multi_Mounts):
	    multi_mounts._program_number_set(mount_program_number)

	# Grab some values from *mount_operations*:
	pairs = mount_operations._pairs
	assert len(pairs) > 0, "No operations present for '{0}'".format(mount_operations._name)
	pair0 = pairs[0]
	mount0 = pair0._mount_get()
	operation0 = pair0._operation_get()
	part = operation0._part_get()
	part_name = part._name_get()

	# Tack *mount_program_number* ont the program numbers list of *part*:
	part._program_numbers_append(mount_program_number)

	# For now the first *operation0* must always be an *Operation_Mount* or an
	# *Operation_Multi_Mount*: 
	assert isinstance(operation0, Operation_Mount) or \
	  isinstance(operation0, Operation_Multi_Mount), \
	  "Got operation '{0}' instead of Operation_Mount or Operation_Multi_Mount". \
	  format(operation0.__class__.__name__)

	# Create the *mount_vrml* for drawing lines into:
	mount_vrml_lines_name = "{0}_{1}".format(part_name, mount0._name_get())
	mount_vrml_lines = VRML_Lines(mount_vrml_lines_name)
	if trace_detail >= 2:
	    print("{0}mount_vrml_lines={1}".format(indent, mount_vrml_lines))

	# Output a coordinate axis shape to *mount_wrl_file*:
	shop = part._shop_get()
	vice = shop._vice_get()
	vice._coordinates_vrml_append(mount_vrml_lines)

	# Open the top-level *mount_ngc_file* that lists both the tool table and
	# calls each tool operation from a single top level .ngc file:
	ezcad = part._ezcad_get()
	ngc_directory = ezcad._ngc_directory_get()
	mount_ngc_file_name = "{0}.ngc".format(mount_program_number)
	with ngc_directory._write_open(mount_ngc_file_name, tracing + 1) as mount_ngc_file:
	    assert mount_ngc_file != None, "Unable to open {0}".format(mount_ngc_file_name)
	    if trace_detail >= 2:
		print("{0}mount_ngc_file_name='{1}' opened".format(indent, mount_ngc_file_name))

	    # Output some heading lines for *mount_ngc_file*:
	    mount_ngc_file.write("( Part: '{0}'    Program_Number: {1}    Mount: '{2}' )\n".
	      format(part._name_get(), mount_program_number, mount0._name_get()))

	    # Associated with each *Operation* is a priority field that stores the number of times
	    # that *Part.cnc_fence*() has been called for the associated *Part* object.  The rule
	    # is that all operations with a lower priority **MUST** be performed before moving
	    # onto the the next priority.  This completely disallows intermixing of operations
	    # that cross a call to CNC fence.

	    # Create *priority_sorted_pairs* so by the operation priority (i.e. *pair*[1]):
	    priority_sorted_pairs = \
	      sorted(pairs, key=lambda pair: pair._operation_get()._priority_get() )

	    # Now create *priority_mount_operations_list* where each *Mount_Operations* object
	    # in the list has *Operation* objects at the same priority level:
	    current_priority = -1
	    priority_mount_operations_list = []
	    current_priority_mount_operations = None
	    for pair in priority_sorted_pairs:
		# Extract *mount* and *operation* from *pair*.  Get the *priority* from *operation*:
		mount = pair._mount_get()
		operation = pair._operation_get()
		priority = operation._priority_get()

		# If the *priority* does not match *current_priority*, we have to create a new
		# *Mount_Operations* object:
		if priority != current_priority:
		    # Create new *Mount_Operation_Object* and tack it onto th end of
		    # *priority_mount_operations_list*:
		    current_priority = priority
		    current_priority_mount_operations = \
		      Mount_Operations("{0} Priority {1}".format(mount_operations._name, priority))
		    priority_mount_operations_list.append(current_priority_mount_operations)

		# Now we can tack *mount* and *operation* onto the
                # *current_priority_mount_operations*:
		current_priority_mount_operations._append(mount, operation)
	    assert len(priority_mount_operations_list) > 0
	
	    # Now we can generate the *priority_operations* from *priority_groups*:
	    #priority_program_number = mount_program_number
	    #for priority_operations in priority_groups:
	    #    # Output a couple of G-code lines *mount_ngc_file*:
	    #    tool_number = tool._number_get()
	    #    tool_name = tool._name_get()
	    #    mount_ngc_file.write("( T{0} {1} )\n".format(tool_number, tool_name))
	    #    mount_ngc_file.write("O{0} call\n".format(priority_program_number))

	    # Create *mount_vrml_stl* and *mount_vrml_lines*:
	    mount_name = mount0._name_get()
	    mount_vrml_stl_name = "{0}__{1}__STL".format(part_name, mount_name)
	    mount_vrml_stl = VRML_Group(mount_vrml_stl_name)
	    #FIXME: *mount_vrml_lines* is define before this loop!!!
	    #mount_vrml_lines_name = "{0}__{1}__Lines".format(part_name, mount_name)
	    mount_vrml_lines = VRML_Lines(mount_vrml_lines_name)

	    # Now generate all of the tool paths:
	    path_bounding_box = Bounding_Box()
	    priority_program_number = mount_program_number
	    mount_vrml = VRML_Group(part._name_get())
	    total_path_time = Time()
	    for index, priority_mount_operations in enumerate(priority_mount_operations_list):
		if trace_detail >= 2:
		    priority_mount_operations._show("Priority {0}".format(index),
		      tracing = tracing + 1)

		# Now we can generate all the tool operations that are at the same pirority:	
		priority_program_number, path_time = \
		  priority_mount_operations._cnc_priority_generate(priority_program_number + 1,
		  mount_ngc_file, mount_vrml, mount_vrml_lines, mount_vrml_stl, path_bounding_box,
		  tracing = tracing + 1)
		total_path_time += path_time

	    # The *next_mount_program_number* must be divisible by 10:
	    next_mount_program_number = priority_program_number + 9
	    next_mount_program_number -= next_mount_program_number % 10

	    # Write out the final G-code lines to *mount_ngc_file*:
	    # FIXME: Read tool change point from *vice*:
	    mount_ngc_file.write("G49 G0 X-1.5 Y0 Z8 ( Return to tool change point )\n")
	    mount_ngc_file.write("( Estimated time: {0:m} )\n".format(total_path_time))
	    #mount_ngc_file.write("G53 G0 Y0.0 ( Move the work to the front )\n")
	    mount_ngc_file.write("M2\n")
	    mount_ngc_file.write("( Path Bounding Box: {0:i} )\n".format(path_bounding_box))
	    mount_ngc_file.write("( Path Bounding Box Volume: {0:i} )\n".
	      format(path_bounding_box.volume_get()))
	    mount_ngc_file.flush()

	# Now append *mount_vrml_lines* to *mount_vrml*:
	mount_vrml._append(mount_vrml_lines)
	mount_vrml._append(mount_vrml_stl)
	if trace_detail >= 2:
	    print("{0}mount_vrml_lines={1}".format(indent, mount_vrml_lines))

	# Now append *mount_stl_vrml* to *mount_wrl_lines*:
	#mount_stl_vrmls = mount_operations._stl_vrmls_get(tracing = tracing + 1)
	#for mount_stl_vrml in mount_stl_vrmls:
	#    mount_stl_vrml._lines_extend(mount_wrl_lines, 2, tracing = tracing + 1)

	# Now append the tool paths for each tool to *mount_wrl_lines*:
	#for tool_path_vrml in tool_path_vrmls:
	#    tool_path_vrml._lines_extend(mount_wrl_lines, 2, tracing = tracing + 1)

	# Close out the *mount_wrl_lines* with the closing VRML brackets that match the header:
	#mount_wrl_lines.append(" ]\n")
	#mount_wrl_lines.append("}\n")

	# Output *mount_wrl_lines* to the *mount_wrl_file_name* file.  Please note that
	# these file names end in `.wrl` and are written into the *ngc_directory* to
        # coexist with their similarly named `.ngc` counterparts:
	ezcad = part._ezcad_get()
	part_name = part._name_get().replace(' ', '_')
	mount_name = mount0._name_get().replace(' ', '_')
	mount_wrl_file_name = "{0}_{1}_{2}.wrl".format(mount_program_number, part_name, mount_name)
	with ngc_directory._write_open(mount_wrl_file_name, tracing + 1) as mount_wrl_file:
	    mount_wrl_file.write("#VRML V2.0 utf8\n")
	    mount_vrml_text = mount_vrml._text_pad(0)
	    mount_wrl_file.write(mount_vrml_text)

	# Wrap up any requested *tracing* and return the *next_mount_program_number*:
	if tracing >= 0:
	    print("{0}<=Mount_Operations._cnc_mount_generate('{1}', {2})=>{3}".format(
	      indent, mount_operations._name, mount_program_number, next_mount_program_number))
        return next_mount_program_number

    def _cnc_priority_generate(self, priority_program_number, mount_ngc_file,
      cnc_vrml, mount_vrml_lines, mount_vrml_stl, path_bounding_box, tracing=-1000000):
        """ *Mount_Operations*: Generate CNC code for each *Operation* in the *Mount_Operations*
	    object (i.e. *self*), where each *Operation* in the *Mount_Operations* object
	    must have the same priority.  This routine is responsible for regrouping tools to
	    minimize tool changes.  *priority_program_number* is the CNC program number to
	    start using.  When the done, the next program number to use is returned.
	    *mount_ngc_file* is the open file into which mount mount level G codes are written
	    (These are usually just a call to the tool routine.)  The tool path visualization
	    is recorded into *cnc_vrml*.
	"""

	# Use *mount_operations* instead of *self*:
	mount_operations = self
        
	# Verify argument types:
	assert isinstance(priority_program_number, int)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML_Group)
	assert isinstance(mount_vrml_lines, VRML_Lines)
	assert isinstance(mount_vrml_stl, VRML_Group)
	assert isinstance(path_bounding_box, Bounding_Box)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	#tracing = 0
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
            print("{0}=>Part._cnc_priority_generate('{1}', {2}, *, '{3}')".format(indent,
	      mount_operations._name, priority_program_number, cnc_vrml._name_get()))
	    trace_detail = 2

	# This is where we regroup tools in *priority_operations*.  We move tools forward
	# providing we do not violate a *follows* constraint:
	pairs = mount_operations._pairs
	pairs_size = len(pairs)
	assert pairs_size > 0
	pair0 = pairs[0]
	mount0 = pair0._mount_get()
	operation0 = pair0._operation_get()
	pregrouped_pairs = pairs[:]	# Shallow list copy of *pairs*
	regrouped_pairs = []
	while len(pregrouped_pairs) > 0:
	    # Grab the first pair off of *pregrouped_pairs* and stuff onto end
	    # of *regrouped_pairs*.  This loop will terminate since we are shortening:
	    # *pregrouped_pairs* by at least 1 each iteration:
	    lopped_pair = pregrouped_pairs.pop(0)
	    regrouped_pairs.append(lopped_pair)
	    assert len(pregrouped_pairs) + len(regrouped_pairs) == pairs_size

	    # Seach the remains of *pregrouped* operation for an operation that matches
	    # *match_tool*:
	    #lopped_mount = lopped_pair._mount_get()
	    lopped_operation = lopped_pair._operation_get()
	    match_tool = lopped_operation._tool_get()
	    index = 0
	    while index < len(pregrouped_pairs):
		pair = pregrouped_pairs[index]
		#mount = pair._mount_get()
		operation = pair._operation_get()
		tool = operation._tool_get()
		if match_tool == tool:
		    # *tool* matches operation *lopped_tool*; now figure out if we can move
		    # *pair* forward without violating a follows constraint:
		    follows = operation._follows_get()
		    follows_found = False
		    for search_pair in pregrouped_pairs[:index]:
			#search_mount = search_pair._mount_get()
			search_operation = search_pair._operation_get()
			if follows == search_operation:
			    follows_found = True
			    break
		    if follows_found:
			# We can not regroup *pair* because its operation it violates a *follows*
			# contraint somewhere before in *pregrouped_pairs*:
			break

		    # We can move it forward.  Be sure not to incrment *index*:
		    matched_pair = pregrouped_pairs.pop(index)
		    assert isinstance(matched_pair, Mount_Operation)
		    regrouped_pairs.append(matched_pair)
		else:
		    # No match for *match_tool*; try the next *operation*:
		    index += 1
	    assert len(pregrouped_pairs) + len(regrouped_pairs) == pairs_size

	# Do some final sanity checking on the regrouping:
	assert len(regrouped_pairs) == pairs_size
	assert len(pregrouped_pairs) == 0

	# Create *regrouped_mount_operations* by creating and an empty *Mount_Operations* object
	# and replacing its *_pairs* field with *regrouped_pairs*:
	regrouped_mount_operations = \
	  Mount_Operations("Regrouped {0}".format(mount_operations._name))
	regrouped_mount_operations._pairs = regrouped_pairs
	if trace_detail >= 2:
	    regrouped_mount_operations._show("regrouped operations", tracing = tracing + 1)

	regrouped_mount_operations._cnc_dowel_pin_optimize(tracing = tracing + 1)

	# Create *tool_mount_operations_list* where each entry is a *Mount_Operations* object
        # that consists entirely of operations on a single *Tool* object:
	tool_mount_operations_list = []
        current_tool = None
	for pair in regrouped_pairs:
	    # Grab some values from *pair*:
	    mount = pair._mount_get()
            operation = pair._operation_get()
	    tool = operation._tool_get()

	    # Has *tool* has changed?
	    if tool != current_tool:
		# Create a new *Mount_Operations* object for the new *current_tool*:
		current_tool = tool
		current_mount_operations = Mount_Operations("Tool {0}".format(tool._name_get()))
		tool_mount_operations_list.append(current_mount_operations)

	    # No matter what, tack *mount* and *operation* onto the end *current_mount_operarions*:
	    current_mount_operations._append(mount, operation)
	assert len(tool_mount_operations_list) > 0

	# Now sweep through *tool_mount_operations_list* to generate a tool path for each tool:
	total_path_time = Time()
	tool_program_number = priority_program_number
	for index, tool_mount_operations in enumerate(tool_mount_operations_list):
	    # If requested, do a little *trace_detail*:
	    if trace_detail >= 2:
		label = "Tool[{0}] before reorder".format(index)
		tool_mount_operations._show(label, tracing = tracing + 1)

	    # Compute *reordered_tool_mount_operations* from *tool_mount_operations*:
	    pairs = tool_mount_operations._pairs_copy()
	    assert len(pairs) > 0
	    pairs0 = pairs[0]
	    operation0 = pairs0._operation_get()
	    reordered_tool_mount_operations = \
	      operation0._reorder(tool_mount_operations, tracing = tracing + 1)
	    if trace_detail >= 2:
		label = "Tool[{0}] after reorder".format(index)
		tool_mount_operations._show(label, tracing = tracing + 1)

	    # Create *tool_vrml_lines* to record the tool path into:
	    part = operation0._part_get()
	    part_name = part._name_get()
	    mount_name = mount._name_get()
	    tool_name = operation0._tool_get()._name_get()
	    tool_vrml_lines_name = \
	      "{0}_{1}_{2}_Ops{3}]".format(part_name, mount_name, tool_name, index)
	    tool_vrml_lines = VRML_Lines(tool_vrml_lines_name)

	    # Now generate the tool path for *reordered_tool_mount_operations*:
	    tool_program_number, path_time = reordered_tool_mount_operations._cnc_tool_generate(
	      tool_program_number, mount_ngc_file, cnc_vrml, mount_vrml_lines, mount_vrml_stl,
              path_bounding_box, tracing = tracing + 1)
	    total_path_time += path_time

	# Compute *new_priority_program_number* (no change from *tool_program_number*):
	new_priority_program_number = tool_program_number

	# Wrap up any requested *tracing* and return *new_priority_program_number*:
	if tracing >= 0:
            print("{0}<=Part._cnc_priority_generate('{1}', {2}, *, '{3}')".format(indent,
	      mount_operations._name, priority_program_number, cnc_vrml._name_get()))
	return new_priority_program_number, total_path_time

    def _cnc_tool_generate(self, tool_program_number, mount_ngc_file,
      cnc_vrml, mount_vrml_lines, mount_vrml_stl, path_bounding_box, tracing=-1000000):
        """ *Mount_Operations*: Generate CNC code for each *Operation* in *Mount_Operations*
	    object (i.e. *self*.)  Each operation must use the exact same tool.  The generated
	    tool `.ngc` file is numbered with *tool_program_number* and the next usable tool
	    program number is returned.  The generated CNC tool path visualiztion is recorded
	    in *cnc_vrml*.
	"""

	# Use *tool_mount_operations* instead of *self*:
	tool_mount_operations = self

	# Verify argument types:
	assert isinstance(tool_program_number, int)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML_Group)
	assert isinstance(mount_vrml_lines, VRML_Lines)
	assert isinstance(mount_vrml_stl, VRML_Group)
	assert isinstance(path_bounding_box, Bounding_Box)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	deep_tracing = -1000000
	#if tool_program_number == 2303:
	#    print("***********************************************************************")
	#    tracing = 5
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount_Operations._cnc_tool_generate('{1}', {2}, *, '{3}')".format(indent,
	      tool_mount_operations._name, tool_program_number, cnc_vrml._name_get()))
	    trace_detail = 3
	    deep_tracing = tracing + 1

	# Grap some values from *tool_mount_operations*:
	pairs = tool_mount_operations._pairs
	assert len(pairs) > 0
	pair0 = pairs[0]
	mount0 = pair0._mount_get()
	mount0_cnc_transform = mount0._cnc_transform_get()
	operation0 = pair0._operation_get()

	# While *Operation_Mount._cnc_generate*() doe not actually generate any CNC code
	# it *DOES* generate VRML for CNC visualization.  Thus, we must actually call
	# the *_cnc_generate*() method for each *Operation_Mount* object.

	# Do not generate any .ngc file for an *Operation_Mount*:
	new_tool_program_number = tool_program_number
	if True:
	    # Grab some values from *tool_mount_operations*:
	    tool = operation0._tool_get()
	    feed_speed = operation0._feed_speed_get()
	    spindle_speed = operation0._spindle_speed_get()
	    part = operation0._part_get()

	    # Remember that *part* uses *tool*:
	    tool._part_register(part)

	    # Create *cnc_tool_vrml* for recording just this *tool*:
	    part_name = part._name_get()
	    mount_name = mount0._name_get()
	    tool_name = tool._name_get()
	    cnc_tool_vrml_name = "{0}__{1}__{2}__Tool".format(part_name, mount_name, tool_name)
	    cnc_tool_vrml = VRML_Group(cnc_tool_vrml_name)

	    # Create *cnc_vrml_lines* to record the tool path visualation:
	    cnc_vrml_lines_name = "{0}__{1}__{2}__Tool_Path". \
	      format(part_name, mount_name, tool_name)
	    cnc_vrml_lines = VRML_Lines(cnc_vrml_lines_name)
	    cnc_tool_vrml._append(cnc_vrml_lines)

	    # Reset the *code* object and feed it *cnc_vrml_lines* to record tool
	    # path visualization:
	    shop = part._shop_get()
	    code = shop._code_get()
	    code._start(part, tool, tool_program_number,
	      spindle_speed, mount0, cnc_vrml_lines, tracing = tracing + 1)
	    if trace_detail >= 2:
		print("{0}A code._vrml_points={1}".format(indent,
		  ["{0:i}".format(vrml_point) for vrml_point in code._vrml_points]))
	    code._dxf_xy_offset_set(part._dxf_x_offset_get(), part._dxf_y_offset_get())
	    if trace_detail >= 2:
		print("{0}B code._vrml_points={1}".format(indent,
		  ["{0:i}".format(vrml_point) for vrml_point in code._vrml_points]))

	    # We need to output a *Part* VRML for each unique (*Part*, *Mount*) combination.
	    # This is done by keeping track of whenever a new pair comes along.
	    part_mount_table = {}

	    # Perform each *operation* in *tool_operations*:
	    pairs_size = len(pairs)
	    for index, pair in enumerate(pairs):
		# Grab the *tool* and its associated *feed_speed* and *spindle_speed*:
		tool_mount = pair._mount_get()
		tool_mount_name = tool_mount._name_get()
		tool_operation = pair._operation_get()
		tool_operation_name = tool_operation._name_get()
		tool = tool_operation._tool_get()
		tool_name = tool._name_get()
		feed_speed = tool_operation._feed_speed_get()
		spindle_speed = tool_operation._spindle_speed_get()
		assert isinstance(feed_speed, Speed)
		assert isinstance(spindle_speed, Hertz)
		if trace_detail >= 2:
		    print("{0}tool_op[{1}]: name={2} tool={3} feed_speed={4:i} spin_speed={5:rpm}".
		      format(indent, index, tool_operation_name, tool_name,
		      feed_speed, spindle_speed))
		    print("{0}C code._vrml_points={1}".format(indent,
		      ["{0:i}".format(vrml_point) for vrml_point in code._vrml_points]))

		# Make sure we start at a safe height:
		code._xy_rapid_safe_z_force(feed_speed, spindle_speed, tracing = tracing + 1)
		if trace_detail >= 2:
		    print("{0}D code._vrml_points={1}".format(indent,
		      ["{0:i}".format(vrml_point) for vrml_point in code._vrml_points]))

		# Perform the CNC generation step for *operation*:
		is_last = (index + 1 >= pairs_size)
		underscore_part_name = part.name_get().replace(' ', '_')
		underscore_tool_name = tool_name.replace(' ', '_')
		underscore_operation_comment = tool_operation._comment_get().replace(' ', '_')
		code._line_comment("PN={0} TN={1} OC={2}".format(
		  underscore_part_name, underscore_tool_name, underscore_operation_comment))
		if trace_detail >= 2:
		    print("{0}E code._vrml_points={1}".format(indent,
		      ["{0:i}".format(vrml_point) for vrml_point in code._vrml_points]))
		tool_operation._cnc_generate(tool_mount, mount_ngc_file,
		  cnc_tool_vrml, mount_vrml_lines, mount_vrml_stl, is_last, tracing = deep_tracing)

		# Make darn we end at a safe height:
		code._xy_rapid_safe_z_force(feed_speed, spindle_speed)

		# Build a ( *Part*, *Mount* ) pair to use as a key to *part_mount_table*:
		#FIXME: Maybe use ( id(part), id(mount) ) instead???!!!
		part = tool_operation._part_get()
		part_name = part._name_get()
		mount_name = tool_mount._name_get()
		part_mount_key = ( part_name, mount_name )
		if not part_mount_key in part_mount_table:
		    if trace_detail >= 1:
			print("{0}part_mount_key is new".format(indent, part_mount_key))
		    part_mount_table[part_mount_key] = True

	    # Safe the tool and shut down *code*:
	    zero = L()
	    code._xy_rapid_safe_z_force(feed_speed, spindle_speed)
	    code._dxf_xy_offset_set(zero, zero)

	    # Output the call to *mount_ngc_file*:
	    tool_name = tool._name_get()
	    if tool_name != "None":
		mount_ngc_file.write("O{0} call ( T{1}: {2} Priority: {3} Time:{4:m})\n".format(
		  tool_program_number, tool._number_get(), tool_name, operation0._priority_get(),
		  code._time_get()))

	    # Update *path_bounding_box*:
	    path_bounding_box.bounding_box_expand(code._bounding_box_get())
	    path_time = code._finish(tracing + 1)
	    assert isinstance(path_time, Time)

	    cnc_tool_vrml._append(mount_vrml_lines)
	    cnc_tool_vrml._append(mount_vrml_stl)

	    # Now output *cnc_tool_vrml* to the *ngc_directory* if it is not a mount operation:
	    if trace_detail >= 3:
		print("{0}operation0='{1}' is_mount={2}".
		  format(indent, operation0._name_get(), operation0._is_mount_get()))
	    if not operation0._is_mount_get():
		ezcad = part._ezcad_get()
		ngc_directory = ezcad._ngc_directory_get()
		cnc_tool_vrml_padded_text = cnc_tool_vrml._text_pad(0)
		wrl_file_name = "{0}.wrl".format(tool_program_number)
		with ngc_directory._write_open(wrl_file_name) as tool_wrl_file:
		    tool_wrl_file.write("#VRML V2.0 utf8\n")
		    tool_wrl_file.write(cnc_tool_vrml_padded_text)
		if trace_detail >= 2:
		    print("{0}Wrote file '{1}'".format(indent, wrl_file_name))

	    # *cnc_tool_vrml* got closed as a result of calling *_text_pad*.  However, we can still
	    # append it to *cnc_vrml* where all of the tools for given mount are visualized:
	    cnc_vrml._append(cnc_vrml_lines)

	    # Compute *new_tool_program_number* for return; it just increments by one:
	    new_tool_program_number = tool_program_number + 1

	# Wrap up any requested *tracing* and return *new_tool_program_number* and *path_time*:
	if tracing >= 0:
	    print("{0}<=Mount_Operations._cnc_tool_generate('{1}', {2}, *, '{3}')=>{4}, {5:m}".
	      format(indent, tool_mount_operations._name, tool_program_number,
	             cnc_vrml._name_get(), new_tool_program_number, path_time))
	return new_tool_program_number, path_time

    def _extend(self, to_mount, from_mount_operations, flags, tracing=-1000000):
	""" *Mount_Operations*: Append the contents of *from_mount_operations* to the
	    *Mount_Operations* object (i.e. *self*) with substituting in *to_mount* for
	    each append.  If *flags* contains an 'M', mount operations will be stripped
	    out during the copy.  If *flags* contains an 'D', dowel pin operations will
	    be stripped out.  If *flags* contains a 'd', all but the first dowel pin operation
	    will be stripped out.  If *flags* contains an 'F', facing operations will be stripped
	    out of the copy.
	"""

	# Use *to_mount_operations* instead of *self*: 
	to_mount_operations = self

	# Verify argument types:
	assert isinstance(to_mount, Mount)
	assert isinstance(from_mount_operations, Mount_Operations)
	assert isinstance(flags, str)
        assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount_Operations._extend('{1}', '{2}' '{3}')".
	      format(indent, to_mount_operations._name_get(),
	      to_mount._name_get(), from_mount_operations._name_get()))
	    trace_detail = 2

	# Extract values from *to_mount_operations* and *from_mount_operations*:
	to_pairs = to_mount_operations._pairs
	from_pairs = from_mount_operations._pairs

	# Set the stripping flags:
	strip_mounts = 'M' in flags
	strip_dowel_pins = 'D' in flags
	if trace_detail >= 1:
	    print("{0}len(from_pairs)={1} len_to_pairs={2}".
	      format(indent, len(from_pairs), len(to_pairs)))

	# Transfer the from *from_pairs* to *to_pairs*:
	for index, from_pair in enumerate(from_pairs):
	    from_operation = from_pair._operation_get()
	    if strip_mounts and isinstance(from_operation, Operation_Mount) or \
	      strip_dowel_pins and isinstance(from_operation, Operation_Dowel_Pin):
		# Strip out *from_operation*:
		if trace_detail >= 2:
		    print("{0}[{1}]: Skip".format(indent, index))
	    else:
		# Copy *from_operation* over to *to_pairs*:
		to_pair = Mount_Operation(from_pair._name_get(), to_mount, from_operation)
		to_pair_mount = to_pair._mount_get()
		to_pairs.append(to_pair)
		assert to_pairs[-1] == to_pair
		if trace_detail >= 2:
		    print("{0}[{1}]: to_pair_mount='{2}'".
		      format(indent, index, to_pair_mount._name_get()))
		    print("{0}[{1}]: to_mount='{2}', from_operation='{3}'".
		     format(indent, index, to_pair_mount._name_get(), from_operation._name_get()))
		    print("{0}[{1}]: to_mount='{2}', from_operation='{3}'".
		     format(indent, index, to_pairs[-1]._mount_get()._name_get(),
		     to_pairs[-1]._operation_get()._name_get()))
	if trace_detail >= 3:
	    to_mount_operations._show("Mount_Operations._extend")

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Mount_Operations._extend('{1}', '{2}' '{3}')".
	      format(indent, to_mount_operations._name_get(),
	      to_mount._name_get(), from_mount_operations._name_get()))

    def _fetch(self, index):
	""" *Mount_Operations*: Return the *index*'th mount operation pair from the
	    *Mount_Operations* object (i.e. *self*).
	"""
	
	# Use *mount_operations* instead of *self*:
	mount_operations = self

	# Verify argument types:
	assert isinstance(index, int)

	# Return the *index*'th *pair*:
	pairs = mount_operations._pairs
	assert index < len(pairs)
	pair = pairs[index]
	return pair

    def _first_mount_get(self, tracing=-1000000):
	""" *Mount_Operations*: Return the first *Mount* object from the *Mount_Operations* object
	    (i.e. *self*.)
	"""

	# Use mount_operations instead of *self*:
	mount_operations = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' '  * tracing
	    print("{0}=>Mount_Operations._first_mount_get('{1}')".
	      format(indent, mount_operations._name))

	# Search for *first_mount* that is not an *Operation_Multi_Mount*:
	first_mount = None
	pairs = mount_operations._pairs
	for pair in pairs:
	    mount = pair._mount_get()
            operation = pair._operation_get()
	    if not isinstance(operation, Operation_Multi_Mount):
		first_mount = mount
		break

	# Wrap up any requested *tracing* and return *first_mount*:
	if tracing >= 0:
	    print("{0}=>Mount_Operations._first_mount_get('{1}')".
	      format(indent, mount_operations._name))
	return first_mount

    def _multi_mounts_get(self):
	""" *Mount_Operations*: Return any *Multi_Mounts* object associtated with the
	    *Mount_Operations* object (i.e. *self*.)  If there is not associated *Multi_Mounts*
	    object, *None* is returned.
	"""

	return self._multi_mounts


    def _multi_mounts_set(self, multi_mounts):
	""" *Mount_Operations*: Set the *Multi_Mounts* object associated with the
	    *Mount_Operations* object (i.e. *self*) to *multi_mounts*.
	"""

	# Use *mount_operations* instead of *self*:
	mount_operations = self

	# Verify argument types:
	assert isinstance(multi_mounts, Multi_Mounts)

	# Stuff *multi_mounts* into *mount_operations*:
	mount_operations._multi_mounts = multi_mounts


    def _name_get(self):
	""" *Mount_Operations*: Return the name associated with the *Mount_Operarions* object
	    (i.e. *self*.)
	"""

	return self._name

    def _pairs_copy(self):
	""" *Mount_Operations*: Return a copy of the (*Mount*, *Operation*) pairs associated with
	    the *Mount_Operations* object (i.e. *self*).
	"""

	return self._pairs[:]

    def _pairs_extend(self, extend_pairs):
	""" *Mount_Operations*: Extend the pairs list associated with the *Mount_Operations*
	    object (i.e. *self*) with *extend_pairs*:
	"""

	self._pairs.extend(extend_pairs)

    def _prepend(self, mount, operation, tracing=-100000):
	""" *Mount_Operations*: Prepend the pair of (*mount*, *operation*) to the *Mount_Operations*
	    object (i.e. *self*) list.
	"""

	# Use *mount_operations* instead of *self*:
	mount_operations = self

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(operation, Operation)
	assert isinstance(tracing, int)

	# Perform any requesting *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount_Operations._prepend('{1}', '{2}', '{3}')".format(indent,
	      mount_operations._name, mount._name_get(), operation._name_get()))

	# Prepend the *mount*/*operation pair to *mount*:
	mount_operations._append(mount, operation, prepend=True, tracing = tracing+1)

	# Wrap up any requesting *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Mount_Operations._prepend'{1}', '{2}')".
	      format(indent, mount._name_get(), operation._name_get()))

    def _program_number_get(self):
	""" *Mount_Operations*: Return the `.ngc` file program number associated wit the
	    *Mount_Operations* object (i.e. *self*.)
	"""

	return self._program_number

    def _show(self, label, tracing=-1000000):
	""" *Mount_Operations*: Print out the *Mount_Operations* object (i.e. *self*) in
	    a human readable format.
	"""

	# Use *mount_operations* instead of *self*:
        mount_operations = self

	# Verify argument types:
	assert isinstance(label, str)
	assert isinstance(tracing, int)

	# Grab some values from *self*:
	name = mount_operations._name
	pairs = mount_operations._pairs

	# Perform any requested *tracing*:
	indent = ""
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount_Operations._show('{1}', '{2}')".
	      format(indent, mount_operations._name, label))

	# Figure out if all *pairs* have the *same_mount*:
	same_mount = True
	if len(pairs) > 0:
	    pair0 = pairs[0]
	    mount0 = pair0._mount_get()
	    for pair in pairs:
		assert isinstance(pair, Mount_Operation), \
		  "Got {0} instead of Mount_Operation in Mount_Operations '{1}'". \
		  format(pair.__class__.__name__, name)
		if pair._mount_get() != mount0:
		    same_mount = False
		    break

	# All the pairs have the *same_mount*, so do a shorter listing:
	name = mount_operations._name
	if same_mount:
	    # Short listing:
	    print("{0}Mount_Operations name='{1}': label='{2}' Mount0 Name='{3}'".
	      format(indent, name, label, mount0._name_get()))
	else:
	    # Long listing:
	    print("{0}Mount_Operations '{1}': '{2}'".format(indent, name, label))

	# Iterate over all the pairs:
	for index, pair in enumerate(pairs):
	    # Extract *mount* and *operation* form *pair*:
	    assert isinstance(pair, Mount_Operation)
	    mount = pair._mount_get()
	    assert isinstance(mount, Mount)
	    operation = pair._operation_get()
	    assert isinstance(operation, Operation)

	    # Figure out how to print *mount_name*:
	    mount_name = ""
	    if same_mount:
		# Short format lists do not list *mount_name* on every line:
		pass
	    else:
		# Long format lists the *mount_name* name on every line:
		part = mount._part_get()
		mount_name = "p='{0:>12}' m='{1:>8}'".format(part._name_get(), mount._name_get())

	    # Extract values from *operation* and *tool*:
	    #operation_show = operation._show(tracing = tracing + 1)
	    operation_show = operation._show(mount)

	    # Print out the line:
	    print("{0}[{1:>2}]: {2} {3}".
	      format(indent, index, mount_name, operation_show))
	    
	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Mount_Operations._show('{1}', '{2}')".
	      format(indent, mount_operations._name, label))

    def _size_get(self):
	""" *Mount_Operations*: Return the size of the *Mount_Operations* object (i.e. *self*).
	"""

	return len(self._pairs)

    def _stl_vrmls_get(self, tracing=-1000000):
	""" *Mount_Operations*: Return a list of *VRML* objects that corresponds to each
	    different *Mount* object in the *Mount_Operations* object (i.e. *self*.)
	"""

	# Use *mount_operations* instead of *self*:
	mount_operations = self
	
	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Mount_Operations._stl_vrmls_get('{1}')".
	      format(indent, mount_operations._name))

	stl_vrmls_table = {}
	stl_vrmls = []
	pairs = mount_operations._pairs
	for pair in pairs:
	    mount = pair[0]
	    mount_name = mount._name_get()
	    operation = pair[1]
	    if not isinstance(operation, Operation_Multi_Mount):
		part = operation._part_get()
		part_name = part._name_get()
		mount_part_name = "{0}_{1}_STL".format(mount_name, part_name)
		if not mount_part_name in stl_vrmls_table:
		    stl_vrmls_table[mount_part_name] = True
		    stl_vrml = mount._stl_vrml_get(tracing = tracing + 1)
		    stl_vrmls.append(stl_vrml)
	
	# Wrap up any requested *tracing* and return *stl_vrmls*:
	if tracing >= 0:
	    print("{0}<=Mount_Operations._stl_vrmls_get('{1}')".
	      format(indent, mount_operations._name))
	return stl_vrmls

class Multi_Mount:
    """ *Multi_Mount*: Specifies where to mount a part for a multi-part path generation.
    """

    def __init__(self, part, mount_name, rotate, tracing = -1000000):
	""" *Multi_Mount: Initialize the *Multi_Mount* object (i.e. *self*) to
	    start using the *Mount* object specified by *part* and *mount_name* with
	    a rotations of *rotate* which must be a mulitple of 90 degrees.  This
	    initialization is frequently augmented by the *Multi_Mount*.*dx_dy*() method
	    and by the *Multi_Mount*.*tooling_plate*() method.
	"""

	# Use *multi_mount* instead of *self*:
	multi_mount = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(mount_name, str)
	assert isinstance(rotate, Angle)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Multi_Mount.__init__(*, '{1}', '{2}', '{3:d}')")

	# Stuff values into *multi_mount*:
	zero = L()
	big = L(mm=123.4567890)
	multi_mount._column          = -1	
	multi_mount._combined_mount  = None	# The *Mount* used to for overall setup
	multi_mount._dx              = -big
	multi_mount._dy              = -big
	multi_mount._extra_bsw       = None
	multi_mount._extra_tne       = None
	multi_mount._mount_name      = mount_name
	multi_mount._part            = part
	multi_mount._program_numbers = []
	multi_mount._rotate          = rotate
	multi_mount._row             = -1

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Multi_Mount.__init__(*, '{1}', '{2}', '{3:d}')")

    def _column_row_get(self):
	""" *Multi_Mount*: Return the row and column offset for the *Multi_Mount* object
	    (i.e. *self*.)
	"""

	# Use *multi_mount* instead of *self*:
	multi_mount = self
	return multi_mount._column, multi_mount._row

    def _combined_mount_get(self):
	""" *Multi_Mount*: Return the final combined mount for the *Multi_Mount* object
	    (i.e. *self*).
	"""

	# Use *multi_mount* instaed of *self*:
	multi_mount = self
        
	combined_mount = multi_mount._combined_mount
	assert isinstance(combined_mount, Mount)
	return combined_mount

    def _combined_mount_set(self, combined_mount, tracing = -1000000):
	""" *Multi_Mount*: Set the final combined mount for the *Multi_Mount* object (i.e. *self*)
	    to *combined_mount*.  This operation can only be done once.
	"""

	# Use *multi_mount* instead of *self*:
	multi_mount = self

	# Verify argument types:
	assert isinstance(combined_mount, Mount)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Multi_Mount._combined_mount_set('{1}', '{2}')".
	      format(indent, multi_mount._mount_name, combined_mount._name_get()))

	# Stuff *combined_mount* into *multi_mount*:
	assert multi_mount._combined_mount == None, \
	  "Multi_Mount '{0}' was previously set to Mount '{1}'; part_stack={2}". \
	  format(multi_mount._mount_name, combined_mount._name_get(),
	    EZCAD3.ezcad._parts_stack_to_text())
	multi_mount._combined_mount = combined_mount

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Multi_Mount._combined_mount_set('{1}', '{2}')".
	      format(indent, multi_mount._mount_name, combined_mount._name_get()))

    def _copy(self, mount_name, tracing=-100000):
	""" *Multi_Mount*: Return a copy of the *Multi_Mount* object (i.e. *self*) where
	    the mount name his replaced with *mount_name*.
	"""

	# Use *multi_mount* instead of *self*
	multi_mount = self

	# Verify argument types:
	assert isinstance(mount_name, str) and not ' ' in mount_name
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Multi_Mount._copy('{1}', '{2}')".
	      format(indent, multi_mount._mount_name, mount_name))

	new_multi_mount = Multi_Mount(multi_mount._part, mount_name, multi_mount._rotate,
	  tracing = tracing + 1)
	new_multi_mount._column         = multi_mount._column
	new_multi_mount._combined_mount = multi_mount._combined_mount
	new_multi_mount._extra_bsw      = multi_mount._extra_bsw
	new_multi_mount._extra_tne      = multi_mount._extra_tne
	new_multi_mount._row            = multi_mount._row

	# Wrap up any requested *tracing* and return *new_multi_mount*:
	if tracing >= 0:
	    print("{0}=>Multi_Mount._copy('{1}', '{2}')=>'{3}'".
	      format(indent, multi_mount._mount_name, mount_name, new_copy._mount_name))
	return new_multi_mount

    def dx_dy(self, dx, dy):
        """ *Multi_Mount*: Adjust the *Multi_Mount* object (i.e. *self*) by (*dx*, *dy).
	"""

	# Use *multi_mount* instead of *self*
	multi_mount = self

	# Verify argument types:
	assert isinstance(dx, L)
	assert isinstance(dy, L)

	# Adjust *multi_mount* by (*dx, *dy*)
	multi_mount._dx += dx
	multi_mount._dy += dy

	return multi_mount
    def _dx_get(self):
	""" *Multi_Mount*: Return the dx associated with the *Multi_Mount* object (i.e. *self*).
	"""

	assert False
	return self._dx

    def _dy_get(self):
	""" *Multi_Mount*: Return the dy associated with the *Multi_Mount* object (i.e. *self*).
	"""

	return self._dy

    def _extra_get(self):
	""" *Multi_Mount*: Return *extra_bsw* and *extra_tne* from the *Multi_Mount* object
	    (i.e. *self*.)
	"""

	# Use *multi_mount* instead of *self*:
	multi_mount = self

	# Return the results:
	return multi_mount._extra_bsw, multi_mount._extra_tne

    def _extra_set(self, extra_bsw, extra_tne):
	""" *Multi_Mount*: Save *extra_bsw* and *extra_tne* into the *Multi_Mount* object
	    (i.e. *self*.)
	"""

	# Use *multi_mount* instead of *self*:
	multi_mount = self

	# Verify argument types:
	assert isinstance(extra_bsw, P)
	assert isinstance(extra_tne, P)

	# Stuff values into *multi_mount*:
	multi_mount._extra_bsw = extra_bsw
	multi_mount._extra_tne = extra_tne

    def _is_tooling_plate_get(self, tracing=-1000000):
	""" *Multi_Mount*: Return *True* if the *Multi_Mount* object (i.e. *self*) has an
	    offset specificed using tooling plate column and row numbers and *False* otherwise.
	"""
	
	# Use *multi_mount* instead of *self:
	multi_mount = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Multi_Mount._is_tooling_plate_get('{1}')".
	     format(indent, multi_mount._mount_name))

	# Compute the *result*:
	mount = multi_mount._mount_get()
	result = mount._is_tooling_plate_mount_get()

	# Wrap up any requested *tracing* and return:
	if tracing >= 0:
	    print("{0}<=Multi_Mount._is_tooling_plate_get('{1}')=>{2}".
	      format(indent, multi_mount._mount_name, result))
	return result

    def _mount_get(self, tracing=-1000000):
        """ *Multi_Mount*: Return the *Mount* object associate with the *Multi_Mount* object
	    (i.e. *self*.)
	"""

	# Use *multi_mount* instead of *self*:
	multi_mount = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Multi_Mount._mount_get('{1}')".format(indent, multi_mount._mount_name))

	# Grab some values from *multi_mount*:
	part = multi_mount._part
	mount_name = multi_mount._mount_name

	# Now lookup the associatied *mount_operations*:
	mount_operations = part._mount_operations_lookup(mount_name)

	# Now find the associated *mount* object:
	mount_operation = mount_operations._fetch(0)
	mount = mount_operation._mount_get()

	# Wrap up any *tracing* and return mount:
	if tracing >= 0:
	    print("{0}<=Multi_Mount._mount_get('{1}')=>'{2}'".
	      format(indent, multi_mount._mount_name, mount._name_get()))
	return mount

    def _mount_name_get(self):
	""" *Multi_Mount*: Return the mount name of the *Multi_Mount* object (i.e. *self*).
	"""

	return self._mount_name

    def _mount_operations_get(self):
        """ *Multi_Mount*: Return the *Mount* object associate with the *Multi_Mount* object
	    (i.e. *self*.)
	"""

	# Use *multi_mount* instead of *self*:
	multi_mount = self

	# Grab some values from *multi_mount*:
	part = multi_mount._part
	mount_name = multi_mount._mount_name

	# Now lookup the associatied *mount_operations*:
	mount_operations = part._mount_operations_lookup(mount_name)

	return mount_operations

    def _part_get(self):
	""" *Multi_Mount*: Return the *Part* object associated with the *Multi_Mount* object
	    (i.e. *self*.)
	"""

        return self._part

    def _rotate_get(self):
	""" *Multi_Mount*: Return the rotate angle associated with the *Multi_Mount* object
	    (i.e. *self*.)
	"""

	return self._rotate

    def tooling_plate(self, column, row, tracing=-1000000):
	""" *Multi_Mount*: Adjust the *Multi_Mount* object (i.e. *self*) to be mounted at
	    (*column*, *row*) on the tooling plate.  The updated *Multi_Mount* object is returned.
	"""

	# Use *multi_mount* instead of *self*:
	multi_mount = self

	# Verify argument types:
	assert isinstance(column, int)
	assert isinstance(row, int)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Multi_Mount('{1}', {2}, {3})".
	      format(indent, multi_mount._name, column, row))

	part = multi_mount._part
	ezcad = part._ezcad_get()
	shop = ezcad._shop_get()
	tooling_plate = shop._tooling_plate_get()
	hole_pitch = tooling_plate._hole_pitch_get()
	multi_mount._dx += column * hole_pitch
	multi_mount._dy -= row * hole_pitch
	multi_mount._column = column
	multi_mount._row = row

	# Wrap up any requested *tracing* and return the updated *multi_mount*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Multi_Mount('{1}', {2}, {3})=>*".
	      format(indent, multi_mount._name, column, row))
	return multi_mount

class Multi_Mounts:
    """ *Multi_Mounts* is basically just a list of *Multi_Mount* objects.
    """

    def __init__(self, name, cell_dx, cell_dy):
	""" *Multi_Mounts*:  Initialize the *Multi_Mounts* object (i.e. *self*) to be named *name*.
	"""

	# Use *multi_mounts* instead of *self*:
	multi_mounts = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(cell_dx, int) and cell_dx >= 0
	assert isinstance(cell_dy, int) and cell_dy >= 0

	# Load up *multi_mounts*:
	multi_mounts._cell_dx = cell_dx
	multi_mounts._cell_dy = cell_dy
	multi_mounts._combined_multi_mount = None
	multi_mounts._multi_mounts_list = []
	multi_mounts._multi_mounts_table = {}
	multi_mounts._name = name
	multi_mounts._part = None
	multi_mounts._program_number = -1

    def _cell_dx_get(self):
	""" *Multi_Mounts*: Return the cell dx offset associated with the *Multi_Mounts*
	    object (i.e. *self*.)
	"""

	return self._cell_dx

    def _cell_dy_get(self):
	""" *Multi_Mounts*: Return the cell dy offset associated with the *Multi_Mounts*
	    object (i.e. *self*.)
	"""

	return self._cell_dy

    def _combined_multi_mount_get(self):
	""" *Multi_Mounts*: Return the combined multi *Mount* for the *Multi_Mounts* object
	    (i.e. *self*.)
	"""

	# Return *combined_multi_mount* after verifying that it was set properly.
	multi_mounts = self
	combined_multi_mount = multi_mounts._combined_multi_mount
	if combined_multi_mount == None:
	    assert False, "combined_multi_mount not set"
	assert isinstance(combined_multi_mount, Mount)
	return combined_multi_mount

    def _combined_multi_mount_set(self, combined_multi_mount):
	""" *Multi_Mounts*: Set the combined multi *Mount* for the *Multi_Mounts* object
	    (i.e. *self*.)
	"""

	# Use *multi_mounts* instead of *self*:
	multi_mounts = self

	# Verify argument types:
	assert isinstance(combined_multi_mount, Mount)

	# Stuff *combined_multi_mount* into *multi_mounts*:
	multi_mounts._combined_multi_mount = combined_multi_mount

    def _copy(self, mount_name, tracing = -1000000):
	""" *Multi_Mounts*: Return of copy of the *Multi_Mounts* object (i.e. *self*)
	    where *mount_name* has been substituted in for all of the mount names.
	"""

	# Use *multi_mounts* instead of *self*:
	multi_mounts = self

	# Verify argument types:
	assert isinstance(mount_name, str) and not ' ' in mount_name
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Multi_Mounts._copy('{1}', '{2}')".
	      format(indent, multi_mounts._name, mount_name))

	new_multi_mounts = Multi_Mounts(mount_name)
	new_mulit_mounts._program_number = mulit_mounts._program_number
	new_multi_mounts_list = new_multi_mounts._multi_mounts_list
	replace_name = None
	for index, multi_mount in enumerate(multi_mounts._multi_mounts_list):
	    multi_mount_name = multi_mount._mount_name_get()
	    if replace_name == None:
		replace_name = multi_mount_name
	    else:
		assert multi_mount_name == replace_name, \
		  "Multi_Mount[{0}] has of '{1}' which does not match previous name(s) of '{2}'". \
		  format(index, multi_mount_name, replace_name)

	    # Copy *multi_mount* with a new name of *mount_name* and append to
	    # *new_multi_mounts_list*:
	    new_multi_mount = multi_mount._copy(mount_name)
	    new_multi_mounts_list.append(new_multi_mount)

	# Wrap up any requested *tracing* and return *new_multi_mounts*:
	if tracing >= 0:
	    print("{0}<=Multi_Mounts._copy('{1}', '{2}')=>'{3}.".
	      format(indent, multi_mounts._name, mount_name, new_multi_mounts._name))
	return new_multi_mounts

    def dx_dy_append(self, part, mount_name, rotate, dx, dy):
	""" *Multi_Mounts*: Append *part*, *mount_name*, *rotate*, *dx*, and *dy* to the
	    *Multi_Mounts* object (i.e. *self*).
	"""
	
	# Use *multi_mounts* instead of *self*:
	multi_mounts = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(mount_name, str)
	assert isinstance(rotate, Angle)
	assert isinstance(dx, L)
	assert isinstance(dy, L)

	# Create *multi_mount*:
	multi_mount = Multi_Mount(part, mount_name, rotate)
	multi_mount.dx_dy(dx, dy)

	# Append to the *list*:
	multi_mounts._multi_mounts_list.append(multi_mount)

    def _extra_start_dowel_pin_get(self, tooling_plate, tracing=-1000000):
	""" *Multi_Mounts*: Return the extra material corners and dowel pin location for
	    the *Multi_Mounts* object (i.e. *self*.)  The extra material corners are relative
	    to the (0, 0) hole of *tooling_plate*.  The dowel pin result is based relative
	    to the tooling plate (0, 0) mounting hole and will need to be adjusted appropriately.
	"""

	# Use *multi_mounts* instead of *self*:
	multi_mounts = self

	# Verify argument types:
	assert isinstance(tooling_plate, Tooling_Plate)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Multi_Mounts._extra_start_dowel_pin_get('{1}', *)".
	      format(indent, tooling_plate._name_get()))
	    trace_detail = 2

	zero = L()
	big = L(inch=123456789.0)
	bounding_box = Bounding_Box()
	left_most_edge = big
	left_most_index = -1
	multi_mounts_list = multi_mounts._multi_mounts_list
	for index, multi_mount in enumerate(multi_mounts_list):
	    # Grab some values from *multi_mount*:
	    mount       = multi_mount._mount_get()
	    column, row = multi_mount._column_row_get()
	    rotate      = multi_mount._rotate_get()

	    # Do some *trace_detail*:
	    if trace_detail >= 2:
		mount_name = mount._name_get()
		part = mount._part_get()
		part_name = part._name_get()
		print("{0}[{1}]:Part='{2}' Mount='{3}' Col={4} row={5} rotate={6:d}deg".
		  format(indent, index, part_name, mount_name, column, row, rotate))
		part_bsw = part.bsw
		part_tne = part.tne
		print("{0}[{1}]:part_bsw={2:i} part_tne={3:i}".
		  format(indent, index, part_bsw, part_tne))
		cnc_transform = mount._cnc_transform_get()
		cnc_part_bsw = cnc_transform * part_bsw
		cnc_part_tne = cnc_transform * part_tne
		print("{0}[{1}]:cnc_part_bsw={2:i} cnc_part_tne={3:i}".
		  format(indent, index, cnc_part_bsw, cnc_part_tne))

	    # Compute the extra material corners centered on hole (0, 0) of the *tooling_plate*:
	    extra_start_bsw, extra_start_tne = \
	      tooling_plate._extra_material_get(mount, column, row, rotate, tracing = tracing + 1)
	    if trace_detail >=2:
		print("{0}[{1}]: extra_start_bsw={2:i} extra_start_tne={3:i}".
		  format(indent, index, extra_start_bsw, extra_start_tne))
		
	    # Save the results into *multi_mount* (for debugging):
	    multi_mount._extra_set(extra_start_bsw, extra_start_tne)

	    # Expand the *bounding_box* to enclose the extra material assocaiated
	    # with *multi_mount*:
	    bounding_box.point_expand(extra_start_bsw)
	    bounding_box.point_expand(extra_start_tne)

	    # Figure out where the *left_most_edge* is:
	    if extra_start_bsw.x < left_most_edge:
		left_most_edge = extra_start_bsw.x
		left_most_index = index
		left_most_bsw = extra_start_bsw
		left_most_tne = extra_start_tne
	        if trace_detail >=2:
		    print("{0}[{1}]: left_most_bsw={2:i} left_most_tne={3:i}".
		  format(indent, index, left_most_bsw, left_most_tne))
	assert index >= 0, \
	  "No left most edge found for Multi_Mounts '{0}'".format(multi_mounts._name)

	# Grab the final values from *bounding_box*:
	final_bsw = bounding_box.bsw_get()
	final_tne = bounding_box.tne_get()

	# Create the *dowel_pin* point:
	dowel_pin = P(final_bsw.x, (left_most_bsw.y + left_most_tne.y)/2, left_most_bsw.z)

	# Perform any requested *tracing* and return the final corners and the left-most index:
	if tracing >= 0:
	    print("{0}=>Multi_Mounts._extra_start_dowel_pin_get('{1}', *)=>{2:i}, {3:i}, {4:i}".
	      format(indent, tooling_plate._name_get(), final_bsw, final_tne, dowel_pin))
	return final_bsw, final_tne, dowel_pin

    def _is_tooling_plate_get(self, tracing=-1000000):
	""" *Multi_Mounts*: Return *True* if the *Multi_Mounts* object (i.e. *self*) only
	    contains tooling mount *Multi_Mount* objects.
	"""

	# Use *multi_mounts* instead of *self*:
	multi_mounts = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Multi_Mounts._is_tooling_plate_get('{1}')".
	      format(indent, multi_mounts._name))
	    trace_detail = 2

	# Before we do anything we need to compute the *extra_bounding_box* that contains
	# all the *multi_mounts*:
	tooling_plate_mounts_count = 0
	vice_mounts_count = 0
	for index, multi_mount in enumerate(multi_mounts._multi_mounts_get()):
	    # Grab some values from *multi_mount*:
	    if multi_mount._is_tooling_plate_get(tracing = tracing + 1):
		tooling_plate_mounts_count += 1
	    else:
		vice_mounts_count += 1

	# Die if both *tooling_plate_mounts_count* and *vice_mounts_count* are both non-zero:
	if trace_detail >= 1:
	    print("{0}vice_mounts_count={1} tooling_plate_mounts_count={2}".
	      format(indent, vice_mounts_count, tooling_plate_mounts_count))
	assert not (tooling_plate_mounts_count > 0 and vice_mounts_count > 0), \
	  "Multi_Mounts '{0}' has a mixture of tooling_plate and dx_dy entries". \
	  format(multi_mounts_name_get())
		
	# Wrap up any requested *tracing* and return *result*;
	result = tooling_plate_mounts_count > 0
	if tracing >= 0:
	    print("{0}<=Multi_Mounts._is_tooling_plate_get('{1}')=>{2}".
	      format(indent, multi_mounts._name, result))
	return result

    def _multi_mounts_get(self):
	""" *Multi_Mounts*: Return the list of *Multi_Mount* objects from the the *Multi_Mounts*
	    object (i.e. *self*):
	"""

	multi_mounts_list = self._multi_mounts_list
	assert isinstance(multi_mounts_list, list)
	return multi_mounts_list

    def _name_get(self):
	""" *Multi_Mounts*: Return the name of the *Multi_Mounts* object (i.e. *self*).
	"""

	return self._name

    def _name_replace_copy(self, old_mount_name, new_mount_name, tracing = -1000000):
	""" *Multi_Mounts*: Return of copy of the *Multi_Mounts* object (i.e. *self*)
	    where *mount_name* has been substituted in for all of the mount names.
	"""

	# Use *multi_mounts* instead of *self*:
	multi_mounts = self

	# Verify argument types:
	assert isinstance(old_mount_name, str) and not ' ' in old_mount_name
	assert isinstance(new_mount_name, str) and not ' ' in new_mount_name
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Multi_Mounts._copy('{1}', '{2}')".
	      format(indent, multi_mounts._name, mount_name))

	assert old_mount_name in multi_mounts._name, \
	  "'{0}' does not contain '{1}'".format(multi_mounts._name, old_mount_name)
	new_multi_mount_name = multi_mounts._name.replace(old_mount_name, new_mount_name)
	new_multi_mounts = \
	  Multi_Mounts(new_multi_mount_name, multi_mounts._cell_dx, multi_mounts._cell_dy)
	new_multi_mounts_list = new_multi_mounts._multi_mounts_list
	replace_name = None
	for index, multi_mount in enumerate(multi_mounts._multi_mounts_list):
	    multi_mount_name = multi_mount._mount_name_get()
	    assert multi_mount_name.find(old_mount_name), \
	      "Multi_Mount[{0}] is named '{1}' which does contain pattern '{2}'". \
	       format(index, multi_mount_name, new_mount_name)

	    # Copy *multi_mount* with a new name of *mount_name* and append to
	    # *new_multi_mounts_list*:
	    new_multi_mount = \
	      multi_mount._copy(multi_mount_name.replace(old_mount_name, new_mount_name))
	    new_multi_mounts_list.append(new_multi_mount)

	# Wrap up any requested *tracing* and return *new_multi_mounts*:
	if tracing >= 0:
	    print("{0}<=Multi_Mounts._copy('{1}', '{2}')=>'{3}.".
	      format(indent, multi_mounts._name, mount_name, new_multi_mounts._name))
	return new_multi_mounts

    def _part_get(self):
	""" *Multi_Mounts*: Return the *Part* object associated witht the *Multi_Mounts* object
	    (i.e. *self*):
	"""

	return self._part

    def _part_set(self, part):
	""" *Multi_Mounts*: Set the *Part* associated with the *Multi_Mounts* object (i.e. *self*)
	    to *part*.
	"""

	# Use *multi_mounts* instead of *self*:
	multi_mounts = self

	# Verify argument types:
	assert isinstance(part, Part)

	# Make sure we set *part* only once:
	#assert multi_mounts._part == None
	multi_mounts._part = part

    def _program_number_get(self):
	""" *Multi_Mounts*: Return the program number associated with the *Multi_Mounts* object:
	    (i.e. *self*.)	    
	"""

	return self._program_number

    def _program_number_set(self, program_number):
	""" *Multi_Mounts*: Store *program_number* into the associated *Multi_Mounts* object
	    (i.e. *self*.)
	"""
	
	assert isinstance(program_number, int) and program_number > 0
	self._program_number = program_number

    def _size_get(self):
	""" *Multi_Mounts*: Returns the number of *Multi_Mount* objects in the *Multi_Mounts*
	    object (i.e. *self*).
	"""

	return len(self._multi_mounts_list)

    def tooling_plate_append(self, part, mount_name, rotate, column, row):
	""" *Multi_Mounts*: Append *part*, *mount_name*, *rotate*, *column*, and *row* to the
	    *Multi_Mounts* object (i.e. *self*).
	"""
	
	# Use *multi_mounts* instead of *self*:
	multi_mounts = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(mount_name, str)
	assert isinstance(rotate, Angle)
	assert isinstance(column, int) and column >= 0
	assert isinstance(row, int) and row >= 0

	# Create a table key:
	part_name = part._name_get()
	key = (part_name, mount_name, rotate.degrees(), column, row)
	multi_mounts_table = multi_mounts._multi_mounts_table
	multi_mounts_list  = multi_mounts._multi_mounts_list

	# See whether or not this *key* is a repeat:
	if not key in multi_mounts_table:
	    # This is the first time for this *key*:

	    # Create *multi_mount* and register *rotate*, *column* and *row*:
	    multi_mount = Multi_Mount(part, mount_name, rotate)
	    multi_mount.tooling_plate(column, row)

	    # Insert *multi_mount* into *multi_mounts_table* and *multi_mounts_list*:
	    multi_mounts_table[key] = multi_mount
	    multi_mounts_list.append(multi_mount)

	    # Tack the *program_number* for *part* onto the *programs_numbers* list:
	    program_numbers = multi_mount._program_numbers
	    program_number = part._program_number_get()
	    program_numbers.append(program_number)

class Operation:
    """ *Operation* is a class that represents a manufacturing operation.
    """

    # These constants are used to sort the order of operations between
    # operations:
    ORDER_NONE =			0
    ORDER_MOUNT = 			1
    ORDER_DOWEL_PIN =			2
    ORDER_END_MILL_EXTERIOR =		3
    ORDER_MILL_DRILL_EXTERIOR =		4
    ORDER_MILL_DRILL_CHAMFER =		5
    ORDER_MILL_DRILL_COUNTERSINK =	6
    ORDER_MILL_DRILL =			7
    ORDER_END_MILL_DRILL =		8
    ORDER_END_MILL_ROUND_POCKET =	9
    ORDER_END_MILL_SIMPLE_POCKET =	10
    ORDER_MILL_DRILL_POCKET_CHAMFER =	11
    ORDER_MILL_DOVE_TAIL_CHAMFER = 	12
    ORDER_DOUBLE_ANGLE_V_GROOVE =	13
    ORDER_DOUBLE_ANGLE_CHAMFER =	14
    ORDER_DRILL =                       15
    ORDER_VERTICAL_LATHE =		16
    ORDER_LAST =			17

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

    def __init__(self, name, kind, part, comment, sub_priority, tool, order,
     follows, feed_speed, spindle_speed, cnc_start, tracing=-1000000):
	""" *Operation*: Initialize an *Operation* object to contain
	    *name*, *kind*, *part*, *comment*, *sub_priority*, *tool*,
	    *order*, *follows*, *feed_speed*, *spindle_speed*, and *cnc_start*.
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
	assert isinstance(follows, Operation) or follows == None
	assert isinstance(feed_speed, Speed)
	assert isinstance(spindle_speed, Hertz)
	assert isinstance(cnc_start, P)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Operation.__init__('{1}', {2}, '{3}', '{4}', {5}, '{6}', {7}, *," +
              " {8:i} {9:rpm} {10:i} {11})").format(indent, name, kind, part._name_get(), comment,
              sub_priority, tool._name_get(), order, feed_speed, spindle_speed, cnc_start, tracing))

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
	operation._cnc_program_number = -1
	operation._cnc_start = cnc_start		# Point in 3D space where CNC starts.
	operation._tracing = tracing

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print(("{0}<=Operation.__init__('{1}', {2}, '{3}', '{4}', {5}, '{6}', {7}, *," +
              " {8:i} {9:rpm} {10:i} {11})").format(indent, name, kind, part._name_get(), comment,
              sub_priority, tool._name_get(), order, feed_speed, spindle_speed, cnc_start, tracing))

    def _cnc_start_get(self):
        """ *Operation*: Return the CNC start point for the *Operation* object (i.e. *self*).
	"""

	return self._cnc_start

    def _is_mount_get(self):
        """ *Operation*: Return whether or not the *Operation* object (i.e. *self*) is
	    a mount operation or not.  This method returns *False* and is overridden to
	    return *True* in the appropation mount operations.
	"""

	return False

    def _comment_get(self):
	""" *Operation*: Return the comment field of the *Operation* object (i.e. *self*). """

	return self._comment

    def _feed_speed_get(self):
	""" *Operation*: Return the feed_speed field of the *Operation* object (i.e. *self*). """

	return self._feed_speed

    def _follows_get(self):
	""" *Operation*: Return the follows field of the *Operation* object (i.e. *self*). """

	return self._follows

    def _order_get(self):
	""" *Operation*: Return the order field of the *Operation* object (i.e. *self*). """

	return self._order

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

    def _reorder(self, mount_operations, tracing=-1000000):
	""" *Operation*: Reorder *mount_operations* to minimize tool travel time between operations.
	"""

	# Use *operation* instead of *self*:
 	operation = self

	# Verify argument types:
	assert isinstance(mount_operations, Mount_Operations)
	assert isinstance(tracing, int)
	
	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
            indent = ' ' * tracing
	    print("{0}=>Operation._reorder('{1}', *)".format(indent, operation._name))
	    trace_detail = 1

	# Make a copy of *pairs* so we can reorder them without damaging the original:
	pairs = mount_operations._pairs_copy()
	pairs_size = len(pairs)
	if trace_detail >= 1:
	    print("{0}operations size={1}".format(indent, pairs_size))
	if trace_detail >= 2:
	    mount_operations._show("before operations reorder", tracing = tracing + 1)

	# Find the operation that is closest to the machine origin and place it first:
	closest_index = None
	closest_distance = L(inch=123456789.0)
	for index, pair in enumerate(pairs):
	    pair_mount = pair._mount_get()
	    pair_cnc_transform = pair_mount._cnc_transform_get()
	    pair_operation = pair._operation_get()
	    cnc_start = pair_operation._cnc_start_get()
	    transformed_cnc_start = pair_cnc_transform * cnc_start
	    distance = transformed_cnc_start.length()
	    if trace_detail >= 3:
		print("{0}[{1}]: cnc_start={2:i} distance={3:i}".
		  format(indent, index, cnc_start, distance))
	    if distance < closest_distance:
		closest_distance = distance
		pairs[0], pairs[index] = pairs[index], pairs[0]
		if trace_detail >= 3:
		    print("{0}closest_index={1}".format(indent, index))

	if pairs_size >= 3:
	    for current_index in range(pairs_size - 2):
		# Grab the *current_pair* and compute its X/Y start location:
		current_pair = pairs[current_index]
		current_mount = current_pair._mount_get()
		current_operation = current_pair._operation_get()
		current_cnc_transform = current_mount._cnc_transform_get()
		current_cnc_start = current_cnc_transform * current_operation._cnc_start
		current_x = current_cnc_start.x
		current_y = current_cnc_start.y	
		if trace_detail >= 2:
		    print("{0}current[{1}]: x={2:i} y={3:i}".
 		      format(indent, current_index, current_x, current_y))

		minimum_distance = L(mm=123456789.0)
		match_index = -1
		for search_index in range(current_index + 1, pairs_size):
		    search_pair = pairs[search_index]
		    search_mount = search_pair._mount_get()
		    search_operation = search_pair._operation_get()
		    search_cnc_transform = search_mount._cnc_transform_get()
		    search_cnc_start = search_cnc_transform * search_operation._cnc_start
		    search_x = search_cnc_start.x
		    search_y = search_cnc_start.y
		    search_distance = (search_x - current_x).distance(search_y - current_y)
		    if trace_detail >= 2:
			print("{0}search[{1}]: x:{2:i} y:{3:i} dist={4:i}".
			  format(indent, search_index, search_x, search_y, search_distance))
                    if search_distance < minimum_distance:
			minimum_distance = search_distance
			match_index = search_index

		# Move *match_pair* be next after *current_pair* in *pairs*:
		if match_index >= 0:
		    # Swap the two operations:
		    if trace_detail >= 2:
                        print("{0}match_index={1}".format(indent, match_index))
		    match_pair = pairs[match_index]
		    pairs[match_index] = pairs[current_index + 1]
		    pairs[current_index + 1] = match_pair

	# Create *reordered_mount_operations* from *pairs*:
	reordered_mount_operations = \
	  Mount_Operations("reordered {0}".format(operation._name), tracing = tracing + 1)
	reordered_mount_operations._pairs_extend(pairs)
	if trace_detail >= 2:
	    reordered_mount_operations._show("after reorder", tracing = tracing + 1)

	# Wrap up any requested *tracing* and return *reordered_mount_operarions*:
	if tracing >= 0:
            indent = ' ' * tracing
	    print("{0}<=Operation._reorder('{1}', *)=>'{2}'".
	      format(indent, operation._name, reordered_mount_operations._name_get()))
	return reordered_mount_operations

    def _show(self, mount, tracing = -1000000):
        """ *Operation*: Return a one line string represenation of the *Operation* object
	    (i.e. *self*):
	"""

	# Use *operation* instead of *self*:
	operation = self

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation._show('{1}')".format(indent, operation._name))

	# Grab some values from *operation*:
	name = operation._name
        priority = operation._priority
	sub_priority = operation._sub_priority
	tool = operation._tool
	tool_name = tool._name_get()
	tool_number = tool._number_get()
	order = operation._order

	# Format up the *result* string:
	result ="n='{0:<13}' p={1} o={2:>2} sp={3} t#={4:>2} tn='{5:<15}'".format(
	    name, priority, order, sub_priority, tool_number, tool_name)

	# Wrap up any requested *tracing* and return formated *result*:
	if tracing >= 0:
	    print("{0}=>Operation._show('{1}')=>'{2}'".format(indent, operation._name, result))
	return result

    @staticmethod
    def _operations_show(operations, indent, comment, tracing = -1000000):
	""" Part debug: Show {self} indented with {indent}. """

	# Verify argument types:
	assert isinstance(operations, list) or isinstance(operations, tuple)
	assert isinstance(indent, str)
	assert isinstance(comment, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation._operations_show(*, '{1}')".format(indent, comment))

	# Print out each *operation* in *operations*:
	print("{0}{1}".format(indent, comment))
	assert isinstance(operations, list)
	for index, operation in enumerate(operations):
            assert isinstance(operation, Operation)
	    print("{0}[{1:>2}]: {2}".format(indent, index, operation._show()))

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Operation._operations_show(*, '{1}')".format(indent, comment))

    def _spindle_speed_get(self):
	""" *Operation*: Return the spindle_speed field of the *Operation* object (i.e. *self*). """

	return self._spindle_speed

    def _sub_priority_get(self):
	""" *Operation*: Return the sub priority field of the *Operation* object (i.e. *self*). """

	return self._sub_priority

    def _tool_get(self):
	""" *Operation*: Return the tool field of the *Operation* object (i.e. *self*). """

	return self._tool

    def _tool_set(self, tool):
	""" *Operation*: Set the tool field of the *Operation* object (i.e. *self*) to *tool*.
	"""

	assert isinstance(tool, Tool)
	self._tool = tool

class Operation_Contour(Operation):
    """ *Operation_Contour* is a class that represents a contour operation.
    """

    # An *Operation_Contour* corresponds to a shape to be milled out.
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
      z_start, z_stop, contour, offset, effective_tool_radius, passes, tracing = -1000000):
	""" *Operation_Contour*: Initialize the *Operation_Contour* object (i.e. *self*)
	    with *part*, *comment*, *sub_priority*, *mill_tool*, *order*, *follows*,
	    *z_start*, *z_stop*, 
	"""

	# Use *operation_contour* instead of *self*:
	operation_contour = self

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

	# Perform any requesteed *tracing*:
	if tracing >= 1:
	    indent = ' ' * tracing
	    print(("{0}=>Operation_Contour.__init__(*, '{1}', '{2}', {3}, '{4}', {5}, {6}," +
	      " s={7:i} f={8:rpm} zs={9:i} ze={10:i} cntr='{11}' off={12:i} tr={13:i} p={14}," +
              " trc={15})").format(indent, part._name_get(), comment, sub_priority,
	      mill_tool._name_get(), order, isinstance(follows, Operation), feed_speed,
              spindle_speed, z_start, z_stop, contour._name_get(), offset, effective_tool_radius,
              passes, tracing))

	# Initialize super class:
	bends = contour._bends_get()
	assert len(bends) > 0
	bend0 = bends[0]
	cnc_start = bend0._point_get()
	Operation.__init__(operation_contour, "Contour", Operation.KIND_CONTOUR, part, comment,
	  sub_priority, mill_tool, order, follows, feed_speed, spindle_speed, cnc_start,
          tracing + 1)

	# Load up the rest of *self*:
	operation_contour._z_start = z_start
	operation_contour._z_stop = z_stop
	operation_contour._contour = contour
	operation_contour._offset = offset
	operation_contour._effective_tool_radius = effective_tool_radius
	operation_contour._passes = passes

	if tracing >= 1:
	    print(("{0}<=Operation_Contour.__init__(*, '{1}', '{2}', {3}, '{4}', {5}, {6}," +
	      " s={7:i} f={8:rpm} zs={9:i} ze={10:i} cntr='{11}' off={12:i} tr={13:i} p={14}," +
              " trc={15})").format(indent, part._name_get(), comment, sub_priority,
	      mill_tool._name_get(), order, isinstance(follows, Operation), feed_speed,
              spindle_speed, z_start, z_stop, contour._name_get(), offset, effective_tool_radius,
              passes, tracing))
	    assert operation_contour._tracing >= 0

    def _cnc_generate(self,
      mount, mount_ngc_file, cnc_vrml, mount_vrml_lines, mount_vrml_stl, is_last, tracing=-1000000):
	""" *Operation_Contour*: Generate the CNC code for the *Operation_Contour* object
	    (i.e. *self*).
	"""

	# Use *operaton_contour* instead of *self*:
	operation_contour = self

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML_Group)
	assert isinstance(mount_vrml_lines, VRML_Lines)
	assert isinstance(mount_vrml_stl, VRML_Group)
	assert isinstance(is_last, bool)
	assert isinstance(tracing, int)

	# Perform an requested *tracing*:
	trace_detail = -1
	if tracing < 0 and operation_contour._tracing >= 0:
	    tracing = operation_contour._tracing
	if tracing >= 0:
	    indent = " " * tracing
	    print("{0}=>Operation_Contour._cnc_generate('{1}', '{2}', '{3}', *, *, *)".
	      format(indent, operation_contour._name, mount._name_get(), cnc_vrml._name_get()))
	    trace_detail = 3

	# Only trace lower levels when *trace_detail* is high:
	tracing_plus_one = -1000000
	if trace_detail >= 3:
	    tracing_plus_one = tracing + 1

	# Grab some values from *operation_contour*:
	part    = operation_contour._part
	tool    = operation_contour._tool
	comment = operation_contour._comment

	# Grab some values from *tool*:
	is_laser = tool._is_laser_get()
	tool_diameter = tool._diameter_get()

	# Get the *code* from *part* and *shop*:
	shop = part._shop_get()
	code = shop._code_get()

	# Record *comment* into *code*:
	code._line_comment(comment)

	# Mark whether *tool* is a laser or not:
	code._is_laser_set(is_laser)

	# Grab some values from *operation_contour*:
	s = operation_contour._spindle_speed_get()
	f = operation_contour._feed_speed_get()
	z_feed = f / 2
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)

	# Grap some more values from *operation_contour*:
	z_start = operation_contour._z_start
	z_stop = operation_contour._z_stop
	passes = operation_contour._passes
	contour = operation_contour._contour
	offset = operation_contour._offset
	radius = operation_contour._effective_tool_radius

	# Force the *contour* to the *cnc_transform* from *mount* to perform the projection:
	cnc_transform = mount._cnc_transform_get()
	contour._project(cnc_transform, tracing=tracing_plus_one)
	if trace_detail >= 2:
	    print("{0}mount.cnc_transform={1:v}".format(indent, cnc_transform))

	# Do the setup step for doing the contour:
	contour._radius_center_and_tangents_compute(tracing = tracing_plus_one)

	# Start a safe height:
	code._xy_rapid_safe_z_force(f, s)

	# Compute *depth_per_pass*:
	z_depth = z_start - z_stop
	depth_per_pass = z_depth / float(passes)

	plunge_offset = tool_diameter

	# Now generate the G-Code.  We need to perform a number of *passes*:
	for index in range(passes):
	    if trace_detail >= 2:
		print("{0}Pass {1} of {2}".format(indent, index + 1, passes))
	    code._line_comment("Pass {0} of {1}".format(index + 1, passes))

	    # Get cutter down to the correct depth:
	    z = z_start - depth_per_pass * float(index + 1)

	    contour = operation_contour._contour
	    code._contour(contour, plunge_offset, offset, radius,
	      True, z, f, s, tracing_plus_one)
	code._xy_rapid_safe_z_force(f, s)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Operation_Contour._cnc_generate('{1}', '{2}', '{3}', *, *, *)".
	      format(indent, operation_contour._name, mount._name_get(), cnc_vrml._name_get()))

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
      diameter, dowel_point, plunge_point, tracing=-1000000):
	""" *Operation_Dowel_Pin*: Initialize an *Operation_Dowel_Pin* object (i.e. *self*)
	    to contain *part*, *comment*, *sub_priority*, *tool*, *order*, *follows*,
	    *feed_speed*, *spindle_speed*, *diameter*, *dowel_point*, and *plunge_point*,
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

	# Perform any requested tracing:
        if tracing >= 0:
	    pad = ' ' * tracing
	    follows_name = 'None'
	    if isinstance(follows, Operation):
		follows_name = follows._name_get()
	    part_name = part._name_get()
            print(("{0}=>Operation_Dowel_Pin.__init__('{1}', '{2}', {3}, '{4}', '{5}', " +
	      "{6:i}, {7:rpm}, {8:i}, dp={9:i}, pp={10:i})").format(
	      pad, part_name, comment, sub_priority, tool._name_get(), follows_name,
	      feed_speed, spindle_speed, diameter, dowel_point, plunge_point))

	# Initialize super class:
	zero = L()
	cnc_start = plunge_point
	Operation.__init__(operation_dowel_pin, "Dowel_Pin", Operation.KIND_DOWEL_PIN, part,
	  comment, sub_priority, tool, order, follows, feed_speed, spindle_speed, cnc_start)
	# Load up the rest of *operation_dowel_pin*:
	operation_dowel_pin._diameter = diameter	# Dowel pin diameter
	operation_dowel_pin._dowel_point = dowel_point	# Location to move dowel to
	operation_dowel_pin._plunge_point = plunge_point # Location to plunge dowel down at

	# Wrap up any *tracing*:
        if tracing >= 0:
	    operation_dowel_pin._tracing = tracing
            print(("{0}<=Operation_Dowel_Pin.__init__('{1}', '{2}', {3}, '{4}', '{5}', " +
	      "{6:i}, {7:rpm}, {8:i}, dp={9:i}, pp={10:i})").format(
	      pad, part_name, comment, sub_priority, tool._name_get(), follows_name,
	      feed_speed, spindle_speed, diameter, dowel_point, plunge_point))

    def _cnc_generate(self,
      mount, mount_ngc_file, cnc_vrml, mount_vrml_lines, mount_vrml_stl, is_last, tracing=-1000000):
	""" *Operation_Dowel_Pin*: Generate the CNC G-code for a an
	    *Operation_Dowel_Pin* object (i.e. *self*.)
	"""

	# Use *dowel_pin* instead of *self*.
	operation_dowel_pin = self

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML)
	assert isinstance(mount_vrml_lines, VRML_Lines)
	assert isinstance(mount_vrml_stl, VRML_Group)
	assert isinstance(is_last, bool)
	assert isinstance(tracing, int)

	# Perform any *tracing*:
	if tracing < 0 and operation_dowel_pin._tracing >= 0:
	    tracing = operation_dowel_pin._tracing
	    print("***************************")
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation_Dowel_Pin._cnc_generate('{1}', '{2}', '{3}', '{4}', '{5}', {6})".
	      format(indent, operation_dowel_pin._name, mount._name_get(), cnc_vrml._name_get(),
	             mount_vrml_lines._name_get(), mount_vrml_stl._name_get(), is_last))
	    trace_detail = 2

	# Grab some values out of *operation_dowel_pin*:
	cnc_dowel_point  = operation_dowel_pin._dowel_point
	cnc_plunge_point = operation_dowel_pin._plunge_point
	comment          = operation_dowel_pin._comment
	tool             = operation_dowel_pin._tool
	part             = operation_dowel_pin._part_get()

	# Grab some values out of *tool*:
	tool_name       = tool._name_get()
	diameter        = tool._diameter_get()
	tip_depth       = tool._tip_depth_get()
	maximum_z_depth = tool._maximum_z_depth_get()

	# Grab some values out of *part*, *shop*, *vice*, and *mount*:
	shop               = part._shop_get()
	code               = shop._code_get()
	vice               = shop._vice_get()
	jaw_volume         = vice._jaw_volume_get()
	top_surface_safe_z = mount._top_surface_safe_z_get()

	# *cnc_dowel_point* is the desired place to put the bottom circular edge of the dowel pin.
	# *actual_dowel_point* adjusts *cnc_dowel_point* so that the portion that is not
	# *diameter* across is not used (i.e. the bottom *tip_depth* of the dowel pin).
	# Also, we need to ensure that the the *maximum_z_depth* for the dowel pin is not exceeded:

	# Grap desired X/Y/Z coordinates from *cnc_dowel_point* and initialize the actual X/Y/Z
	# coordinates:
	desired_x = cnc_dowel_point.x
	desired_y = cnc_dowel_point.y
	desired_z = cnc_dowel_point.z
	actual_x = desired_x
	actual_y = desired_y
	actual_z = desired_z

	# Adjust the location of *cnc_dowel_point* to *actual_dowel_point* to take into
        # account the dowel point *diameter* and *maximum_z_depth*:
	if cnc_plunge_point.x < desired_x:
	    actual_x -= diameter / 2
	else:
	    actual_x += diameter / 2
	actual_z = (desired_z - tip_depth).maximum(top_surface_safe_z - maximum_z_depth)
	actual_dowel_point = P(actual_x, actual_y, actual_z)

	# The actual plunge point needs to have its Z coordinate adjusted to match the Z
	# coordinate of *actual_dowel_point*:
	actual_plunge_point = P(cnc_plunge_point.x, cnc_plunge_point.y, actual_z)
	if trace_detail >= 1:
	    print("{0}tool_name='{1}'".format(indent, tool_name))
	    print("{0}cnc_dowel_point={1:i}".format(indent, cnc_dowel_point))
	    print("{0}actual_dowel_point={1:i}".format(indent, actual_dowel_point))
	    print("{0}actual_plunge_point={1:i}".format(indent, actual_plunge_point))

	# Define the speed and feed for these operations:
	ipm10 = Speed(in_per_sec=10.0)
	rpm0 = Hertz()

	# Output the *operation_dowel_pin* comment supplied by the user:
	code._line_comment(comment)

	# Rapid over to the plunge point:
	code._xy_rapid_safe_z_force(ipm10, rpm0)
	code._xy_rapid(actual_plunge_point.x, actual_plunge_point.y)

	# Now pause to let operator see if Z-safe is at the right height:
	tool_alternate_number = tool._alternate_number_get()
	code._command_begin()
	code._unsigned("M", 6)
	code._unsigned("T", tool_alternate_number)
	code._comment("Operator may check that Z-safe is correct")
	code._command_end()

	# Output `G43 Htool_number` to enable tool offset:
	code._command_begin()
	# We are already be in G43 mode, so this G-code will not show up; We are just be paranoid:
	code._unsigned("G8", 43)
	code._unsigned("H", tool_alternate_number)
	code._comment("Enable tool offset for tool {0}".format(tool_alternate_number))
	code._command_end()

	# Output some information about the *dowel_point* for G-code debugging:
	#code._line_comment(
	#  "z_stop={0:i} tip_depth={1:i}".format(z_stop, tip_depth))

	# Move slowly down to *z_stop*:
	code._z_feed(ipm10, rpm0, actual_plunge_point.z, "dowel_pin")

	# Move slowly to (*actual_x*, *actual_y*).  This may cause the material in the
	# vice to slide over:
	code._line_comment("actual_dowel_point={0:i}".format(actual_dowel_point))
	code._xy_feed(ipm10, rpm0, actual_x, actual_y)

	# Now pause again, to let the operator move piece up against the
	# the dowel pin (if it is not already up against there):
	code._command_begin()
	code._unsigned("M", 6)
	code._unsigned("T", tool._number_get())
	code._comment("Operator should place part against dowel pin")
	code._command_end()

	# Slowly retract away from the part edge back to *plunge_x* and get
	# back up to Z safe:
	code._xy_feed(ipm10, rpm0, cnc_plunge_point.x, cnc_plunge_point.y)
	code._xy_rapid_safe_z_force(ipm10, rpm0)

	# Visualize *tool* into *mount_vrml*:
	part_name = part._name_get()
	mount_name = mount._name_get()
	tool_name = tool._name_get()
	tip_vrml_lines_name = "{0}__{1}__{2}_Dowel_Tip".format(part_name, mount_name, tool_name)
	tip_vrml_lines = VRML_Lines(tip_vrml_lines_name)
	tool._vrml_append(tip_vrml_lines, 2, actual_dowel_point, tip_depth, tracing + 1)
	cnc_vrml._append(tip_vrml_lines)

	# Wrap-up any *tracing*:
	if tracing >= 0:
	    print("{0}<=Operation_Dowel_Pin._cnc_generate('{1}', '{2}', '{3}', '{4}', '{5}', {6})".
	      format(indent, operation_dowel_pin._name, mount._name_get(), cnc_vrml._name_get(),
	             mount_vrml_lines._name_get(), mount_vrml_stl._name_get(), is_last))

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
      diameter, hole_kind, start, stop, is_countersink, tracing=-1000000):
	""" *Operation_Drill*: Initialize *Operation_Drill* to contain
	    *diameter*, *hole_kind*, *start*, *stop*, and *countersink*.
	"""

	# Use *operation_drill* instead of *self*:
	operation_drill = self

	# Verify argument types:
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
	assert isinstance(tracing, int)

	# Initialize the super class:
	cnc_start = start
	Operation.__init__(operation_drill, "Drill", Operation.KIND_DRILL, part, comment,
	  sub_priority, tool, order, follows, tool._feed_speed_get(), tool._spindle_speed_get(),
          cnc_start)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Operation_Drill.__init__(*, p='{1}', c='{2}', t='{3}', o={4}, *," +
              " d={5:i}, hk={6}, start={7:i}, stop={8:i}, i={9})").
              format(indent, part._name_get(), comment, tool._name_get(),
              order, diameter, hole_kind, start, stop, is_countersink))
	    operation_drill._tracing = tracing

	# Load up *operation_drill*:
	operation_drill._diameter = diameter
	operation_drill._hole_kind = hole_kind
	operation_drill._start = start
	operation_drill._stop = stop
	operation_drill._is_countersink = is_countersink

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}<=Operation_Drill.__init__(*, p='{1}', c='{2}', t='{3}', o={4}, *," +
              " d={5:i}, hk={6}, start={7:i}, stop={8:i}, i={9})").
              format(indent, part._name_get(), comment, tool._name_get(),
              order, diameter, hole_kind, start, stop, is_countersink))

    def _cnc_generate(self,
      mount, mount_ngc_file, cnc_vrml, mount_vrml_lines, mount_vrml_stl, is_last,tracing=-1000000):
	""" *Operation_Drill*: Generate the CNC G-code for an *Operation_Drill* object
	    (i.e. *self*).
	"""

	# Use *drill* instead of *self*:
	operation_drill = self

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML_Group)
	assert isinstance(mount_vrml_lines, VRML_Lines)
	assert isinstance(mount_vrml_stl, VRML_Group)
	assert isinstance(is_last, bool)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing < 0 and operation_drill._tracing >= 0:
	    tracing = operation_drill._tracing
	deep_tracing = -1000000
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation_Drill._cnc_generate('{1}', *, '{2}', '{3}')".
	      format(indent, operation_drill._name, mount._name_get(), cnc_vrml._name_get()))
	    trace_detail = 3
	    deep_tracing = tracing + 1

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
	cnc_transform = mount._cnc_transform_get()
	shop = part._shop_get()
	code = shop._code_get()

	radius = diameter / 2
	half_radius = radius / 2

	cnc_start = cnc_transform * start
	cnc_stop = cnc_transform * stop
	is_laser = tool._is_laser_get()
	if trace_detail >= 1:
	    print("{0}start={1:i} cnc_start={2:i}".format(indent, start, cnc_start))
	    print("{0}stop={1:i} cnc_stop={2:i}".format(indent, stop, cnc_stop))
	    print("{0}is_laser={1}".format(indent, is_laser))
	if trace_detail >= 2:
	    print("{0}E code._vrml_points={1}".
	      format(indent, ["{0:i}".format(vrml_point) for vrml_point in code._vrml_points ]))
	top_surface_safe_z = mount._top_surface_safe_z_get()
	if is_laser:
	    code._dxf_circle(x, y, diameter, tracing=tracing + 1)
	else:
	    cnc_start_z = cnc_start.z
	    cnc_stop_z  = cnc_stop.z
	    z_start = cnc_start_z.maximum(cnc_stop_z)
	    z_stop  = cnc_start_z.minimum(cnc_stop_z)
	    if hole_kind == Part.HOLE_THROUGH:
		assert kind == Operation.KIND_DRILL
		tool_drill = tool
		point_angle = tool_drill._point_angle_get()
		tip_depth = tool_drill._tip_depth_get()
		z_stop -= tip_depth + L(inch=0.040)
	    elif hole_kind == Part.HOLE_TIP:
		z_stop = cnc_stop.z
	    elif hole_kind == Part.HOLE_FLAT:
		assert False, "Flat holes can't be done with point drills"
	    else:
		assert False, "Unknown hole kind"
	    if trace_detail >= 1:
		print("{0}z_start={1:i} z_stop={2:i}".format(indent, z_start, z_stop))

	    # Make darn sure we start a high enough Z:
	    code._xy_rapid_safe_z_force(feed_speed, spindle_speed, tracing=deep_tracing)

	    #drill_depth = top_surface_safe_z - z_stop
	    drill_depth = cnc_start.z - cnc_stop.z
	    trip_depth = 3 * diameter
	    z_safe = top_surface_safe_z
	    if drill_depth > trip_depth:
		# Compute *q* which is the peck distance:
		pecks = int(drill_depth / trip_depth) + 1
		# Add just a little to make sure *pecks* times *q* > *drill_depth*;
		# The drill will *never* go below the Z value:
		q = (drill_depth / float(pecks)) + L(inch=.005)

		# For debugging, show the computation of *pecks* and *q*:
		#code._line_comment("drill_depth={0:i} deep_depth={1:i} pecks={2} q={3:i}".
		#  format(drill_depth, trip_depth, pecks, q))

		# Use a "canned" cycle to peck out the deep hole:
		f = feed_speed		#FIXME: Slow feed down by 20%!!!
		p = None
		#r = top_surface_safe_z + L(inch=0.1000)
		r = z_start
		s = spindle_speed
		x = cnc_start.x
		y = cnc_start.y
		z = z_stop
		code._drill_cycle(98, 83, f, p, q, z_safe, r, s, x, y, z, tracing = deep_tracing)
	    else:
		# Use a "canned" cycle to drill the hole:
		f = feed_speed
		p = None
		q = None
		#r = top_surface_safe_z + L(inch=0.1000)
		r = z_start
		s = spindle_speed      # Leave speed alone, since we are almost always under optimim
		x = cnc_start.x
		y = cnc_start.y
		z = z_stop
		code._drill_cycle(98, 81, f, p, q, z_safe, r, s, x, y, z, tracing = deep_tracing)

	    # Make darn sure we get back to a high enough Z:
	    code._xy_rapid_safe_z_force(feed_speed, spindle_speed, tracing = deep_tracing)

	    if is_last:
		code._command_begin()
		code._unsigned("G9", 80)
		code._command_end()

	# Wrap-up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Operation_Drill._cnc_generate('{1}', *, '{2}', '{3}')".
	      format(indent, operation_drill._name, mount._name_get(), cnc_vrml._name_get()))

	    # OLD CODE:
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

    def _show(self, mount, tracing = -1000000):
        """ *Operation_Drill*: Return text showing the values of the *Operation_Drill* object
	    (i.e. *self*) in a short 1 line format.
	"""

	# Use *operation_drill* instead of *self*:
	operation_drill = self

	# Verify argument types:
	assert isinstance(tracing, int)
	assert isinstance(mount, Mount)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation_Drill._show('{1}')".format(indent, operation_drill._name))

	# Grab some values from *operation_drill*:
        start = operation_drill._start
        stop = operation_drill._stop
        hole_kind = operation_drill._hole_kind
	part = operation_drill._part

	# Grab *cnc_transform* from *mount*:
        cnc_transform = mount._cnc_transform_get()

	# Compute *cnc_start* and *cnc_stop*, the location of *start* and *stop* in CNC coordinates:
        cnc_start = cnc_transform * start
        cnc_stop = cnc_transform * stop

	# Compute *x*, *y*, and *depth*:
	x = cnc_start.x
        y = cnc_start.y
        depth = cnc_start.z - cnc_stop.z
	depth_string = "{0:i}".format(depth)[:6]	# Only show first 6 characters

	# Grab the front half of the text from *Operation* super class:
	prefix = Operation._show(operation_drill, mount)
	suffix = " x={0:i} y={1:i} d={2} hk={3:}".format(x, y, depth_string, hole_kind)
	result = prefix + suffix

	# Wrap up any requested *tracing* and return *result*:
	if tracing >= 0:
	    print("{0}<=Operation_Drill._show('{1}')=>'{2}'".
	      format(indent, operation_drill._name, result))
	return result

class Operation_Mount(Operation):
    """ *Operation_Mount* is a class that indicates that the part mount position has changed.
    """

    def __init__(self, part, comment, vice, jaws_spread,
      parallels_height, tooling_plate, spacers, tooling_plate_present, tracing=-1000000):
	""" *Operation_Mount*: Initialize the *Operation_Mount* object (i.e. *self*) to describe
	    how *part* is mounted for CNC operations.  *comment* shows up in generated G-code.
	    *vice* specifies which vice is being used; `None` is used if there is no vice.
	    *jaws_spread* specifics the distance in Y between the two vice jaws (negative for
	    no vice.)  If a vice is being used, parallels height specifies what height should
	    be selected from the vice parallels (negative for none.)  *tooling_plate* specifies
	    whether or not a tooling plate should be mounted on top of the parallels.  Use `None`
	    if no tooling plate is needed.  *spacers* is a list of tooling plate hole coordinate
	    quadruples, which specify where tooling plate spacers are placed
	    (e.g. [(x1,y1,x2,y2),...,(xN,yN,xN+1,yN+1)] .)  The list can be empty if there is
	    no tooling plate.
	"""

	# Use *operation_mount* instead of *self*:
	operation_mount = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(comment, str)
	assert isinstance(vice, Vice) or vice == None
	assert isinstance(jaws_spread, L)
	assert isinstance(parallels_height, L)
	assert isinstance(tooling_plate, Tooling_Plate) or tooling_plate == None
	assert isinstance(spacers, list)
	assert isinstance(tooling_plate_present, bool)
	assert isinstance(tracing, int)
	for spacer in spacers:
	    assert isinstance(spacer, tuple) and len(spacer) == 4
	    for row_column in spacer:
		assert isinstance(row_column, int) and row_column >= 0

	# Perform any requested *tracing*:
	if tracing >= 1:
	    indent = ' ' * tracing
	    print(("{0}=>Operation_Mount.__init__(*, '{1}', '{2}', '{3}'," +
	      " {4:i} {5:i}, *, {6}, {7}))").format(indent, part._name_get(), comment,
	      vice._name_get(), jaws_spread, parallels_height, spacers, tooling_plate_present))

	# Initialize the *Operation* super class:
	sub_priority = 0
	tool = Tool("None", -1, Tool.KIND_MOUNT, Tool.MATERIAL_NONE, L(), 0, L())
	order = Operation.ORDER_MOUNT
	follows = None
	feed_speed = Speed()
        spindle_speed = Hertz()
	Operation.__init__(operation_mount, "Mount", Tool.KIND_MOUNT, part, comment,
	  sub_priority, tool, order, follows, feed_speed, spindle_speed, part.c)

	# Load up *operation_mount*:
	operation_mount._jaws_spread           = jaws_spread
	operation_mount._parallels_height      = parallels_height
	operation_mount._tooling_plate         = tooling_plate
	operation_mount._tooling_plate_present = tooling_plate_present
	operation_mount._vice                  = vice

	# Wrap up any requested *tracing*:
	if tracing >= 1:
	    print(("{0}<=Operation_Mount.__init__(*, '{1}', '{2}', '{3}'," +
	      " {4:i} {5:i}, *, {6}, {7}))").format(indent, part._name_get(), comment,
	      vice._name_get(), jaws_spread, parallels_height, spacers, tooling_plate_present))

    def _cnc_generate(self,
      mount, mount_ngc_file, cnc_vrml, mount_vrml_lines, mount_vrml_stl, is_last, tracing=-1000000):
        """ *Operation_Mount*: Generate the CNC operations for the *Operation_Mount* object
	    (i.e. *self*).
	"""

	# Use *operation_mount* instead of *self*:
	operation_mount = self

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML_Group)
	assert isinstance(mount_vrml_lines, VRML_Lines)
	assert isinstance(mount_vrml_stl, VRML_Group)
	assert isinstance(is_last, bool)
        assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation_Mount._cnc_generate('{1}', '{2}', *, '{3}')".
	      format(indent, operation_mount._name, mount._name_get(), cnc_vrml._name_get()))
	    trace_detail = 1

	# Write out the part visualization:
	part = operation_mount._part
	top_surface_transform = mount._top_surface_transform_get()
	cnc_transform = mount._cnc_transform_get()
	comment = "Operation_Mount._cnc_generate: Mount='{0}' cnc_transform={1:m}".format(
	  mount._name_get(), cnc_transform)
	if tracing >= 0:
	    print("{0}cnc_transform={1:s}".format(indent, cnc_transform))
				
	# Grab some values out *operation_mount*:
	parallels_height      = operation_mount._parallels_height
	tooling_plate         = operation_mount._tooling_plate
	vice                  = operation_mount._vice

	# Grab the extra material conners from *mount* and render into *mount_vrml_lines*:
	extra_start_bsw, extra_start_tne = mount._extra_start_get()
	cnc_extra_start_bsw = cnc_transform * extra_start_bsw
	cnc_extra_start_tne = cnc_transform * extra_start_tne
	mount_vrml_lines._box_outline("azure", cnc_extra_start_bsw, cnc_extra_start_tne)
	if trace_detail >= 2:
	    print("{0}extra_start_bsw={1:i} extra_start_tne={2:i}".
	      format(indent, extra_start_bsw, extra_start_tne))

	# Write the dimensions of the initial *material* into *mount_ngc_file*:
	cnc_extra_start_dx = cnc_extra_start_tne.x - cnc_extra_start_bsw.x
	cnc_extra_start_dy = cnc_extra_start_tne.y - cnc_extra_start_bsw.y
	cnc_extra_start_dz = cnc_extra_start_tne.z - cnc_extra_start_bsw.z
	material = part._material_get()
	mount_ngc_file.write("( Initial Dimensions X:{0:i}in x Y:{1:i}in x Z:{2:i}in of {3} )\n".
	  format(cnc_extra_start_dx.absolute(), cnc_extra_start_dy.absolute(),
	  cnc_extra_start_dz.absolute(), material))
	selected_parallel_height = mount._selected_parallel_height_get()
	if mount._is_tooling_plate_mount_get():
	    mount_ngc_file.write("( Parallels height {0:i}in for tooling plate )\n".
	      format(selected_parallel_height))
	else:
	    mount_ngc_file.write("( Parallels height {0:i}in )\n".format(selected_parallel_height))

	# If there is a tooling plate, visualize that into *cnc_vrml_lines*:
	jaws_spread = (cnc_extra_start_tne.y - cnc_extra_start_bsw.y).absolute()
	zero = L()
	if tooling_plate != None:
	    corner = P(zero, zero, parallels_height)
	    spacers = mount._spacers_get()
	    tooling_plate._vrml_append(mount_vrml_lines, corner, spacers, mount_ngc_file,
	      tracing = tracing + 1)
	    jaws_spread = tooling_plate._dy_get()
	if trace_detail >= 1:
	    print("{0}jaws_spread={1:i}".format(indent, jaws_spread))

	# Show the *parallels* in *cnc_vrml_linew*:
	parallels = vice._parallels_get()
	parallels._vrml_append(mount_vrml_lines,
	  parallels_height, jaws_spread, tracing = tracing + 1)

	# Show the *vice* jaws:
	vice_jaw_volume = vice._jaw_volume_get()
	corner1 = P(zero,               zero,        zero)
	corner2 = P(vice_jaw_volume.x,  zero,        vice_jaw_volume.z)
	mount_vrml_lines._box_outline("yellow", corner1, corner2)
	corner1 = P(zero,              -jaws_spread, zero)
	corner2 = P(vice_jaw_volume.x, -jaws_spread, vice_jaw_volume.z)
	mount_vrml_lines._box_outline("yellow", corner1, corner2)

	# Show the coordinate axes:
	ezcad = part._ezcad_get()
	shop = ezcad._shop_get()
	vice = shop._vice_get()
	vice._coordinates_vrml_append(mount_vrml_lines)
	
	# See whether or not to visualize the extra material around *part*:
	zero = L()
	if True:
	    # Yes let's visualize the extra material:
	    mount_translate_point = mount._mount_translate_point_get()
	    cnc_bsw = extra_start_bsw + mount_translate_point
	    cnc_tne = extra_start_tne + mount_translate_point

	# Read in the `.stl` file associated with *part* and append to *mount_vrml_stl*:
	stl_vrml = part._stl_vrml_get()
	cnc_transform = mount._cnc_transform_get()
	transformed_stl_vrml = cnc_transform._vrml(stl_vrml, tracing = tracing + 1)
	mount_vrml_stl._append(transformed_stl_vrml)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Operation_Mount._cnc_generate('{1}', '{2}', *, '{3}')".
	      format(indent, operation_mount._name, mount._name_get(), cnc_vrml._name_get()))

    def _is_mount_get(self):
        """ *Operation_Mount*: Return whether or not the *Operation* object (i.e. *self*) is
	    a mount operation or not.  This method returns *True*.
	"""

	return True

    def _parallels_height_get(self):
        """ *Operation_Mount*: Return the vice parallels height from the *Operation_Mount* object
	    (i.e. *self*.)
	"""

	return self._parallels_height

class Operation_Multi_Mount(Operation):
    """ *Operation_Multi_Mount* is a class that manages multiple parts in the same operation.
    """

    def __init__(self, part, multi_mounts, comment, vice, jaws_spread,
      parallels_height, tooling_plate, spacers, tooling_plate_present, tracing=-1000000):
	""" *Operation_Multi_Mount*: Initialize the *Operation_Muili_Mount* object (i.e. *self*)
	    mount each *Multi_Mount* object from *multi_mounts* for a multiple mount CNC
	    path generation.  *part* is the *Part* that is the "parent" for the CNC path
	    generation and the *part* name shows up in the summary `.wrl` file name.
	    *vice* is the *Vice* object to use form mounting and *jaws_spread* specifies
	    the distance between the vice jaws.  *parallels_height* specifies which height
	    from the box of parallels to use.  *tooling_plate* is the tooling plate to use
	    for mounting parts.  *spacers* specifies a list of spacess to mount the parts
	    on the tooling plate.  *spacers* is set to None, if no spacers are to be rendered
	    as part the the CNC path visualization.
	"""

	# Use *operation_multi_mount* instead of *self*:
	operation_multi_mount = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(multi_mounts, Multi_Mounts)
	assert isinstance(comment, str)
	assert isinstance(vice, Vice) or vice == None
	assert isinstance(jaws_spread, L)
	assert isinstance(parallels_height, L)
	assert isinstance(tooling_plate, Tooling_Plate)
	assert isinstance(spacers, list)
	assert isinstance(tooling_plate_present, bool)
	for spacer in spacers:
	    assert isinstance(spacer, tuple) and len(spacer) == 4
	    for row_column in spacer:
		assert isinstance(row_column, int) and row_column >= 0

	# Perform any requested *tracing*:
	if tracing >= 1:
	    indent = ' ' * tracing
	    print(("{0}=>Operation_Mount.__init__(*, p='{1}', fm='{2}', c='{3}'," +
	      " v='{4}' js={5:i} ph={5:i}, *, {6}, {7}))").format(indent,
              part._name_get(),  comment, vice._name_get(), jaws_spread, parallels_height,
	      spacers, tooling_plate_present))

	# Initialize the *Operation* super class:
	sub_priority = 0
	tool = Tool("None", -1, Tool.KIND_MOUNT, Tool.MATERIAL_NONE, L(), 0, L())
	order = Operation.ORDER_MOUNT
	follows = None
	feed_speed = Speed()
        spindle_speed = Hertz()
	cnc_start = part.c
	Operation.__init__(self, "Multi_Mount", Tool.KIND_MOUNT,
	  part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed, cnc_start)

	# Load up *operation_multi_mount*:
	operation_multi_mount._jaws_spread = jaws_spread
	operation_multi_mount._multi_mounts = multi_mounts
	operation_multi_mount._parallels_height = parallels_height
	operation_multi_mount._spacers = spacers
	operation_multi_mount._tooling_plate = tooling_plate
	operation_multi_mount._tooling_plate_present = tooling_plate_present
	operation_multi_mount._vice = vice

	# Wrap up any requested *tracing*:
	if tracing >= 1:
	    print(("{0}<=Operation_Mount.__init__(*, p='{1}', fm='{2}', c='{3}'," +
	      " v='{4}' js={5:i} ph={5:i}, *, {6}))").format(indent,
	      part._name_get(), comment, vice._name_get(), jaws_spread, parallels_height, spacers))

    def _cnc_generate(self, combined_multi_mount,
      mount_ngc_file, cnc_vrml, mount_vrml_lines, mount_vrml_stl, is_last, tracing=-1000000):
        """ *Operation_Multi_Mount*: Perform CNC visualization the *Operation_Multi_Mount* object
	    (i.e. *self*).  This ultimately results in CNC tool paths being rendered into
	    *cnc_vrm", mount tooling being rendered into *mount_vrml_lines* and parts being
	    rendering into *mount_vrml_stl*.  The top level G code is written to *mount_ngc_file*.
	    *is_last* is not used for this operation.
	"""

	# Use *operation_multi_mount* instead of *self*:
	operation_multi_mount = self

	# Verify argument types:
	assert isinstance(combined_multi_mount, Mount)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML_Group)
	assert isinstance(mount_vrml_lines, VRML_Lines)
	assert isinstance(mount_vrml_stl, VRML_Group)
	assert isinstance(is_last, bool)
        assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation_Multi_Mount._cnc_generate('{1}', '{2}', *, '{3}')".format(indent,
	      operation_multi_mount._name, combined_multi_mount._name_get(), cnc_vrml._name_get()))
	    trace_detail = 2

	# Grab some values out *operation_multi_mount*:
	jaws_spread           = operation_multi_mount._jaws_spread
	multi_mounts          = operation_multi_mount._multi_mounts
	parallels_height      = operation_multi_mount._parallels_height
	spacers               = operation_multi_mount._spacers
	tooling_plate         = operation_multi_mount._tooling_plate
	tooling_plate_present = operation_multi_mount._tooling_plate_present
	vice                  = operation_multi_mount._vice

	# Set *debugging* to *True* to add some extra material bounding boxes to the CNC
	# path visulation:
	debugging = False

	# Now visualize each *multi_mount*:
	zero = L()
	first_mount = None
	first_part = None
	for index, multi_mount in enumerate(multi_mounts._multi_mounts_get()):
	    # Unpack *multi_mount* enough to get the associated *mount_operations* as well as the
	    # *dx* and *dy* offsets:
	    #dx = multi_mount._dx_get()
	    #dy = multi_mount._dy_get()
	    combined_mount   = multi_mount._combined_mount_get()
	    mount_name       = multi_mount._mount_name_get()
	    part             = multi_mount._part_get()
	    mount_operations = part._mount_operations_lookup(mount_name)

	    # Create *vrml_triangles* that is associated with *part*:
	    part_color_name = part._color_get()._name_get()
	    part_name = part._name_get()
	    part_stl = part._stl_get()
	    part_triangles = part_stl._triangles_get()
	    vrml_triangles_name = "{0}_Multi_Mount[{1}]".format(part_name, index)
	    vrml_triangles = VRML_Triangles(vrml_triangles_name,
	      part_color_name, part_triangles, tracing = tracing + 1)

	    # Now unpack the first *pair* from *mount_operations*:
	    assert mount_operations._size_get() > 0
	    pair = mount_operations._fetch(0)
	    pair_mount = pair._mount_get()
	    assert isinstance(pair_mount, Mount)

	    # The first time through the loop we remember *first_mount* and *first_part*:
	    if first_mount == None:
		first_mount = pair_mount
		first_part = part

	    # Now translate *mount_cnc_transform* by (*dx*, *dy*, 0):
	    mount_cnc_transform = combined_mount._cnc_transform_get()
	    #mount_dx_dy_cnc_transform = mount_cnc_transform.translate("dx_dy", P(dx, dy, zero))

	    # Now place *vrml_triangles* into the CNC visualization *mount_vrml_stl*
	    # with the correct *dx*/*dy* offset:
	    mount_vrml_triangles = \
	      mount_cnc_transform._vrml(vrml_triangles, tracing = tracing + 1)
	    if not debugging:
		mount_vrml_stl._append(mount_vrml_triangles)

	    if debugging:
		extra_bsw, extra_tne = multi_mount._extra_get()
		mount_vrml_lines._box_outline("purple", extra_bsw, extra_tne)

	assert isinstance(first_mount, Mount) and isinstance(first_part, Part)

	# Create *cnc_vrml_lines* to visualize the multi mount:
	first_part_name = first_part._name_get()
	combined_mount_name = combined_mount._name_get()
	cnc_vrml_lines_name = "{0}__{1}__Multi_Mount".format(first_part_name, combined_mount_name)
	cnc_vrml_lines = VRML_Lines(cnc_vrml_lines_name)
	cnc_vrml._append(cnc_vrml_lines)

	# Show the *parallels* in *mount_vrml_lines*:
	parallels = vice._parallels_get()
	parallels._vrml_append(mount_vrml_lines,
	  parallels_height, jaws_spread, tracing = tracing + 1)

	# Show the vice jaws in *mount_vrml_lines*:
	vice_jaw_volume = vice._jaw_volume_get()
	corner1 = P(zero,               zero,        zero)
	corner2 = P(vice_jaw_volume.x,  zero,        vice_jaw_volume.z)
	mount_vrml_lines._box_outline("yellow", corner1, corner2)
	corner1 = P(zero,              -jaws_spread, zero)
	corner2 = P(vice_jaw_volume.x, -jaws_spread, vice_jaw_volume.z)
	mount_vrml_lines._box_outline("yellow", corner1, corner2)

	# Show the coordinate axes:
	vice._coordinates_vrml_append(mount_vrml_lines)
	
	# Visualize the extra material:
	part = operation_multi_mount._part
	#top_surface_transform = combined_mount._top_surface_transform_get()
	cnc_transform = combined_mount._cnc_transform_get()
	mount_extra_start_bsw, mount_extra_start_tne = combined_mount._extra_start_get()
	mount_extra_start_dx = mount_extra_start_tne.x - mount_extra_start_bsw.x
	mount_extra_start_dy = mount_extra_start_tne.y - mount_extra_start_bsw.y	
	mount_extra_start_dz = mount_extra_start_tne.z - mount_extra_start_bsw.z

	# Draw the tooling plate if *tooling_plate_present*:
	if tooling_plate_present:
	    corner = P(zero, zero, parallels_height)
	    tooling_plate._vrml_append(mount_vrml_lines, corner, spacers, mount_ngc_file,
	      tracing = tracing + 1)

	# Draw the extra material outline:        
	extra_start_bsw, extra_start_tne = combined_multi_mount._extra_start_get()
	mount_vrml_lines._box_outline("azure", extra_start_bsw, extra_start_tne)
	if trace_detail >= 2:
	    print("{0}mount_name='{1}' extra_start_bsw={2:i} extra_start_tne={3:i}".
	      format(indent, combined_mount._name_get(), extra_start_bsw, extra_start_tne))

	# Write the *material* dimensions out to *mount_ngc_file*:
	material = part._material_get()
	extra_start_dx = extra_start_tne.x - extra_start_bsw.x
	extra_start_dy = extra_start_tne.x - extra_start_bsw.y
	extra_start_dz = extra_start_tne.x - extra_start_bsw.z
	mount_ngc_file.write("( Initial dimensions X:{0:i}in x Y:{1:i}in x Z:{2:i}in of {3} )\n".
	  format(extra_start_dx, extra_start_dy, extra_start_dz, material))
	parallels_height = combined_multi_mount._selected_parallel_height_get()
	mount_ngc_file.write("( Parallels Height: {0:i}in)\n".format(parallels_height))

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Operation_Multi_Mount._cnc_generate('{1}', '{2}', *, '{3}')".format(indent,
	      operation_multi_mount._name, combined_multi_mount._name_get(), cnc_vrml._name_get()))

    def _is_mount_get(self):
        """ *Operation_Multi_Mount*: Return whether or not the *Operation* object (i.e. *self*) is
	    a mount operation or not.  This method returns *True*.
	"""

	return True

    def _jaws_spread_set(self, jaws_spread):
	""" *Operation_Multi_Mount*: Set the vice jaws spread for the *Operation_Multi_Mount*
	    object (i.e. *self*) to *jaws_spread*:
	"""

	self._jaws_spread = jaws_spread

class Operation_Round_Pocket(Operation):
    """ *Operation_Round_Pocket* is a class that implements a round pocket
	manufacturing operation.
    """

    def __init__(self, part, comment, sub_priority, tool, order, follows, diameter,
      countersink_diameter, hole_kind, start, stop, feed_speed, spindle_speed, tracing=-100000):
	""" *Operation_Round_Pocket*: will intitialize an
	    *Operation_Round_Pocket* object (i.e. *self*) to contain
	    *part*, *comment*, *sub_priority*, *tool*, *order*, *follows*,
	    *diameter*, *hole_kind*, *start*, *stop*, *feed_speed*, and *spindle_speed*.
	"""

	# Use *operation_round_pocket* instead of *self*:
	operation_round_pocket = self

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

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Operation_Round_Pocket.__init__('{1}', '{2}', '{3}', {4}, '{5}', {6}," +
	           " *, {7:i} {8:i}, {9}, {10:i}, {11:i}, {12:i}, {13:rpm}").
	      format(indent, "Round_Pocket", part._name_get(), comment, sub_priority, tool, order,
	      diameter, countersink_diameter, hole_kind, start, stop, feed_speed, spindle_speed))

	# Initialize superclass:
	cnc_start = start
	Operation.__init__(operation_round_pocket, "Round_Pocket", Operation.KIND_ROUND_POCKET,
	  part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed, cnc_start)
	if tracing >= 0:
	    operation_round_pocket._tracing = tracing

	# Load up *operation_round_pocket*:
	operation_round_pocket._diameter = diameter
	operation_round_pocket._countersink_diameter = countersink_diameter
	operation_round_pocket._hole_kind = hole_kind
	operation_round_pocket._start = start
	operation_round_pocket._stop = stop

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print(("{0}<=Operation_Round_Pocket.__init__('{1}', '{2}', '{3}', {4}, '{5}', {6}," +
	           " *, {7:i} {8:i}, {9}, {10:i}, {11:i}, {12:i}, {13:rpm})").
              format(indent,"Round_Pocket", part._name_get(), comment, sub_priority, tool, order,
	      diameter, countersink_diameter, hole_kind, start, stop, feed_speed, spindle_speed))

    def _cnc_generate(self,
      mount, mount_ngc_file, cnc_vrml, mount_vrml_lines, mount_vrml_stl, is_last, tracing=-1000000):
	""" *Operation_Round_Pocket*: Generate the CNC G-code for an *Operation_Round_Pocket*
	    object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML_Group)
	assert isinstance(mount_vrml_lines, VRML_Lines)
	assert isinstance(mount_vrml_stl, VRML_Group)
	assert isinstance(is_last, bool)
	assert isinstance(tracing, int)

	# Use *operation_round_pocket* instead of *self*:
	operation_round_pocket = self

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing < 0 and operation_round_pocket._tracing >= 0:
	    tracing = operation_round_pocket._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation_Round_Pocket._cnc_generate('{1}', '{2}', *, '{3}')".
	      format(indent, operation_round_pocket._name, mount._name_get(), cnc_vrml._name_get()))
	    trace_detail = 2

	# Extract some values from *operation_round_pocket*:
	comment   = operation_round_pocket._comment
	diameter  = operation_round_pocket._diameter
	hole_kind = operation_round_pocket._hole_kind
	part      = operation_round_pocket._part
	start     = operation_round_pocket._start
	stop      = operation_round_pocket._stop
	tool      = operation_round_pocket._tool

	# The *top_surface_transform* has been previously set orient the material correctly for CNC:
	cnc_transform = mount._cnc_transform_get()
	mapped_start = cnc_transform * start
	mapped_stop = cnc_transform * stop
	if trace_detail >= 2:
	    print("{0}cnc_transform={1:v}".format(indent, cnc_transform))

	# Compute some values based on *diameter*:
	maximum_depth = diameter / 3.0
	radius = diameter / 2

	# Extract some values from *part* and *shop*:
	shop = part._shop
	code = shop._code_get()

	# Extract some values from *tool*:
	tool_diameter = tool._diameter_get()
	f = tool._feed_speed_get()
	s = tool._spindle_speed_get()

	# Figure out if *tool* is a laser:
	is_laser = tool._is_laser_get()

	# Compute some values based on {tool_diameter}:
	tool_radius = tool_diameter / 2 
	half_tool_radius = tool_radius / 2

	if trace_detail >= 1:
	    print("{0}is_laser={1}".format(indent, is_laser))
	if is_laser:
	    # We just cut a simple circle:
	    if trace_detail >= 2:
		print("{0}start={1:i} stop={2:i}".format(indent, mapped_start, mapped_stop))
	    code._dxf_circle(mapped_start.x, mapped_start.y, radius - tool_radius)
	else:
	    # We do all the work to mill out the round_pocket pocket;

	    # Deal with through holes:
	    is_through = False
	    x = mapped_start.x
	    y = mapped_start.y
	    z_start = mapped_start.z
	    z_stop = mapped_stop.z
	    if hole_kind == Part.HOLE_THROUGH:
		is_through = True
		z_stop -= L(inch=0.025)

	    code._line_comment(comment)
	    code._line_comment(
	      "x={0:i} y={1:i} radius={2:i} z_start={3:i} z_stop={4:i} tool_radius={5:i}".
	      format(x, y, radius, z_start, z_stop, tool_radius))
	    
	    z_depth = (start - stop).length()
	    depth_passes = int(z_depth / maximum_depth) + 1
	    depth_per_pass = z_depth / float(depth_passes)
	    assert depth_passes < 100

	    # Move to position:
	    code._xy_rapid_safe_z_force(f, s)
	    code._xy_rapid(x, y)

	    z_feed = f / 4.0
	    shave = L(inch=0.005)
	    radius_remove = radius - shave - tool_radius
	    for depth_pass in range(depth_passes):
		code._line_comment(
		  "{0} round_pocket pocket [depth pass {1} of {2}] radius_remove={3:i}".
		  format(comment, depth_pass + 1, depth_passes, radius_remove))
		
		# Get to proper depth:
		z = z_start - depth_per_pass * float(depth_pass + 1)
		if trace_detail >= 2:
		    print("{0}z={1:i}".format(indent, z))

		if is_through:
		    # We are going all the way through, so we can ignore the material in the middle:
		    code._xy_feed(f, s, x, y + radius_remove)
		    code._z_feed(z_feed, s, z, "round_pocket_pocket", tracing=tracing+1)
		    code._xy_ccw_feed(f, s, radius_remove, x, y - radius_remove)
		    code._xy_ccw_feed(f, s, radius_remove, x, y + radius_remove)
		else:
 		    # We have to mow out all the intervening space:
		    radius_passes = int(radius_remove /  half_tool_radius) + 1
		    pass_remove_delta = radius_remove / float(radius_passes)
		    code._line_comment("radius_passes={1} pass_remove_delta={1:i}".
		      format(radius_passes, pass_remove_delta))

		    #for radius_index in range(radius_passes):
		    for radius_index in range(radius_passes):
			pass_remove = pass_remove_delta * float(radius_index + 1)
			code._line_comment("radius pass {0} of {1}: pass_remove={2:i}".
			  format(radius_index + 1, radius_passes, pass_remove))
			code._xy_feed(f, s, x, y + pass_remove)
			code._z_feed(z_feed, s, z, "round_pocket_pocket[{0}]".format(radius_index),
			  tracing=tracing + 1)
			code._xy_ccw_feed(f, s, pass_remove, x, y - pass_remove)
			code._xy_ccw_feed(f, s, pass_remove, x, y + pass_remove)
		    code._line_comment("radius passes done")
		    code._xy_feed(f, s, x, y)

	    # Do a "spring pass" to make everybody happy:
	    if True:
		code._line_comment("{0} round_pocket pocket 'spring' pass".format(comment))
		path_radius = radius - tool_radius
		half_path_radius = path_radius / 2
		code._xy_feed(f, s, x, y)

		# Carefully feed the tool to the edge:
		code._line_comment("Carefully feed tool to the edge of the hole")
		code._xy_ccw_feed(f, s,
		  half_path_radius, x + half_path_radius, y + half_path_radius)
		code._xy_ccw_feed(f, s, half_path_radius, x, y + path_radius)

		# Peform the entire spring cut:
		code._line_comment("Peform the entire spring cut")
		code._xy_ccw_feed(f, s, path_radius, x, y - path_radius)
		code._xy_ccw_feed(f, s, path_radius, x, y + path_radius)

		# Carefully remove the tool back to the center:
		code._line_comment("Carefully remove the tool back to the center")
		code._xy_ccw_feed(f, s,
		  half_path_radius, x - half_path_radius, y + half_path_radius)
		code._xy_ccw_feed(f, s, half_path_radius, x, y)

	    # Safely retract to z safe:
	    code._line_comment("Safely retract to z safe")
	    code._xy_rapid_safe_z_force(f, s)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Operation_Round_Pocket._cnc_generate('{1}', '{2}', *, '{3}')".
	      format(indent, operation_round_pocket._name, mount._name_get(), cnc_vrml._name_get()))

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

    def _cnc_generate(self,
      mount, mount_ngc_file, cnc_vrml, mount_vrml_lines, mount_vrml_stl, is_last, tracing=-1000000):
	""" *Operation_Simple_Exterior*: Generate the CNC G-code for an
	    *Operation_Simple_Exterior* object (i.e. self).
	"""

	assert False, "This code should be deleted"

	# Verify argument types:
	assert isinstnace(mount, Mount)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML_Group)
	assert isinstance(mount_vrml_lines, VRML_Lines)
	assert isinstance(mount_vrml_stl, VRML_Group)
	assert isinstance(is_last, bool)
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

    ANGLE0 = Angle()

    def __init__(self,
      part, comment, sub_priority, tool, order, follows, feed_speed, spindle_speed,
      corner1, corner2, corner_radius, tool_radius, pocket_kind, rotate=ANGLE0, tracing=-1000000):
	""" *Operation_Simple_Pocket*: Initialize an *Operation_Simple_Pocket*
	    object (i.e. *self*) to contain *part*, *comment*, *sub_priority*,
	    *tool*, *order*, *follows*, *corner1*, *corner2*, *corner_radius*,
	    *tool_radius*, and *pocket_kind*.
	"""

	# Use *operation_simple_pocket* instead of *self*:
	operation_simple_pocket = self

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
	assert isinstance(rotate, Angle)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Operation_Simple_Pocket.__init__(*, p='{1}', c='{2}', sp={3}, t='{4}'," +
	      " o={5}, f={6:i}, s={7:rpm}, c1={8:i}, c2={9:i}, cr={10:i}, tr={11:i}, pk={12}," +
	      " rot={13})").format(indent, part._name_get(), comment, sub_priority,
              tool._name_get(), order, feed_speed, spindle_speed, corner1, corner2,
              corner_radius, tool_radius, pocket_kind, rotate))
	    trace_detail = 2

	# Initialize superclass:
	operation_kind = Operation.KIND_SIMPLE_POCKET
	cnc_start = (corner1 + corner2)/2
	Operation.__init__(operation_simple_pocket, "Simple_Pocket", operation_kind, part, comment,
          sub_priority, tool, order, follows, feed_speed, spindle_speed, cnc_start, tracing + 1)

	# Load up the rest of *operation_simple_pocket*:
	operation_simple_pocket._corner1 = corner1
	operation_simple_pocket._corner2 = corner2
	operation_simple_pocket._corner_radius = corner_radius
	operation_simple_pocket._tool_radius = tool_radius
	operation_simple_pocket._rotate = rotate
	operation_simple_pocket._pocket_kind = pocket_kind

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print(("{0}<=Operation_Simple_Pocket.__init__(*, p='{1}', c='{2}', sp={3}, t='{4}'," +
	      " o={5}, f={6:i}, s={7:rpm}, c1={8:i}, c2={9:i}, cr={10:i}, tr={11:i}, pk={12}," +
              " rot={13})").format(indent, part._name_get(), comment, sub_priority,
	      tool._name_get(), order, feed_speed, spindle_speed, corner1, corner2,
              corner_radius, tool_radius, pocket_kind, rotate))

    def _corner_radius_get(self):
	""" *Operation_Simple_Pocket: Return the corner radius for the *Operation_Simple_Pocket*
	    object (i.e. *self*).
	"""

	return self._corner_radius

    def _corners_get(self):
        """ *Operation_Simple_Pocket*: Return the two corners from the *Operation_Simple_Pocket*
	    object (i.e. *self*.)  These two corners a from the pocket in the original 3D space.
	"""

	return self._corner1, self._corner2

    def _cnc_generate(self,
      mount, mount_ngc_file, cnc_vrml, mount_vrml_lines, mount_vrml_stl, is_last, tracing=-1000000):
	""" *Operation_Simple_Pocket*: Generate the CNC G-code for a
	    *Operation_Simple_Pocket* object (i.e. *self*).
	"""

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML_Group)
	assert isinstance(mount_vrml_lines, VRML_Lines)
	assert isinstance(mount_vrml_stl, VRML_Group)
	assert isinstance(is_last, bool)
	assert isinstance(tracing, int)

	# Use *operation_simple_pocket* instead of *self*:
	operation_simple_pocket = self

	# Perform any requested *tracing*:
	part          = operation_simple_pocket._part
	if tracing < 0:
            if part._tracing >= 0:
		tracing = part._tracing
	    elif operation_simple_pocket._tracing >= 0:
		tracing = operation_simple_pocket._tracing
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Operation_Simple_Pocket._cnc_generate('{1}', '{2}', *, *, *, {3})".
	      format(indent, operation_simple_pocket._comment, mount._name_get(), is_last))
	    trace_detail = 2

	# Grap some values from *operation_simple_pocket*:
	comment       = operation_simple_pocket._comment
	corner_radius = operation_simple_pocket._corner_radius
	corner1       = operation_simple_pocket._corner1
	corner2       = operation_simple_pocket._corner2
	part          = operation_simple_pocket._part
	pocket_kind   = operation_simple_pocket._pocket_kind
	rotate        = operation_simple_pocket._rotate
	tool          = operation_simple_pocket._tool
	tool_radius   = operation_simple_pocket._tool_radius

	# Extract some values from *part* and *shop*:
	shop = part._shop_get()
	code = shop._code_get()

	# Grab some values from *tool*:
	feed_speed           = tool._feed_speed_get()
	spindle_speed        = tool._spindle_speed_get()
	tool_maximum_z_depth = tool._maximum_z_depth_get()
	code._tool_set(tool)

	# Transform *corner1* and *corner2* to mounted CNC coordinates:
	cnc_transform = mount._cnc_transform_get()
	cnc_corner1 = cnc_transform * corner1
	cnc_corner2 = cnc_transform * corner2
	cnc_corner_bsw, cnc_corner_tne = cnc_corner1.minimum_maximum(cnc_corner2)
	if trace_detail >= 3:
	    print("{0}cnc_transform={1:v}".format(indent, cnc_transform))
	if trace_detail >= 2:
	    print("{0}cnc_corner1={1:i} cnc_corner2={2:i}".
	      format(indent, cnc_corner1, cnc_corner2))
	    print("{0}cnc_corner_bsw={1:i} cnc_corner_tne={2:i}".
	      format(indent, cnc_corner_bsw, cnc_corner_tne))

	# Figure out the *cnc_volume*:
	cnc_volume = cnc_corner_tne - cnc_corner_bsw
	cnc_volume_dx = cnc_volume.x
	cnc_volume_dy = cnc_volume.y
	cnc_volume_dz = cnc_volume.z
	if trace_detail >= 2:
	    print("{0}cnc_volume={1:i}".format(indent, cnc_volume))

	# Define some constants:
	zero = L()
	z_extra = zero
	#mount_translate_point = mount._mount_translate_point_get()
	mount_translate_point = P(zero, zero, zero)

	# Start with {comment}:
	code._line_comment(comment)

	# Output the pocket operations:
	code._line_comment(
	  "cnc_corner_bsw={0:i} cnc_corner_tne={1:i} tool_radius={2:i} corner_radius={3:i}".
	  format(cnc_corner_bsw, cnc_corner_tne, tool_radius, corner_radius))

	# Figure out the *minimum_span* across the pocket::
	minimum_span = zero
	if cnc_volume_dx > cnc_volume_dy:
	    minimum_span = cnc_volume_dy
	else:
	    minimum_span = cnc_volume_dx
	half_minimum_span = minimum_span / 2
	if trace_detail >= 2:
	    print("{0}minimum_span={1:i} half_minimum_span={2:i}".
	      format(indent, minimum_span, half_minimum_span))

	# Emit the preparatory G-Code:
	is_laser = tool._is_laser_get()
	if trace_detail >= 1:
            print("{0}is_laser={1}".format(indent, is_laser))

	# Figure out the *maximum_pass_depth*:
	# FIXME: Half tool radius is really conservative.  *tool_radius* should work fine:
	maximum_pass_depth = tool_radius / 2
	material = part._material_get()
	if material._is_plastic():
	    maximum_pass_depth = tool_radius
	if is_laser:
	    maximum_pass_depth = tool_maximum_z_depth
	if trace_detail >= 2:
	    print("{0}tool_radius={1:i} maximum_pass_depth={2:i}".
	      format(indent, tool_radius, maximum_pass_depth))

	# Compute *total_cut* which is the total depth of the pocket including some extra:
	# if the pocket goes all the way through:
	total_cut = cnc_volume_dz + z_extra
	if trace_detail >= 1:
	    print("{0}total_cut={1:i} z_extra={2:i}".format(indent, total_cut, z_extra))

	# Compute the number of required *passes*:
	passes = int(total_cut / maximum_pass_depth)
	while maximum_pass_depth * float(passes) < total_cut:
	    passes = passes + 1
	if trace_detail >= 2:
	    print("{0}passes={1}".format(indent, passes))
	assert passes > 0, \
	  "total_cut={0:i} corner1={1:i} corner2={2:i}".format(total_cut, corner1, corner2)
	assert maximum_pass_depth * float(passes) >= total_cut

	# Compute *z_step* which is the mount remvoed for each pass:
	z_step = total_cut / float(passes)

	# *tool_radius* is the radius of the selected end-mill, which we will henceforth
        # call R. We need to be careful on our inner passes because the end-mill can only
        # cut R/sqrt(2) along the diagnoal.  If each pass is separated by exactly 2R, there
        # will be little triangular islands left behind.
	#
	# We now define T=R/sqrt(2).  To avoid the little island problems, each inner pass
        # must overlap the next pass out by at least (R-T).  Diagramatically:
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
	# where the over lap between two passes is R-T.  In the example above, the offset
        # from the center for pass 1 is R, pass 2 is R+T+R, and for pass 3 is R+T+R+T+R.
        # Of course, the center to edge distance may not be 2R+N*(R+T), for some value of N.
	# The better way to organize this is to offset the end-mill from the edge by R on
        # the pass N, R+T+R on the pass N-1, and R+T+R+T+R on pass N-2, etc.

	r = tool_radius
	# 0.65 < 0.707 = 1/sqrt(2)
	#t :@= smul@(r, 0.65)
	t = r / 2

	code._xy_rapid_safe_z_force(feed_speed, spindle_speed)

	is_laser = isinstance(tool, Tool_End_Mill) and tool._is_laser_get()
	if is_laser:
	    z_end = cnc_corner_tne.z - total_cut
	    code._simple_pocket_helper(cnc_corner_bsw, cnc_corner_tne, corner_radius, z_end,
	    tool_radius, zero, spindle_speed, feed_speed, True, rotate, tracing = tracing + 1)
	else:
	    # Compute the total number of rectangular paths needed:
	    paths = 0
	    remaining = half_minimum_span
	    while remaining > zero:
		if paths == 0:
		    # *r* + *r* does not work if pocket width equals tool width,
		    # so we just use *r* + *t* every time now:
		    remaining -= r + t
		else:
		    remaining -= r + t
		paths += 1

	    # Generate *passes* deep depth passes over the pocket:
	    z_start = cnc_corner_tne.z
	    for depth_pass in range(passes):
		# Output "(Depth Pass # of #)" comment:
		code._line_comment("Depth Pass {0} of {1}".
		  format(depth_pass + 1, passes))

		# Compute the plunge depth for each pass:
		z = z_start - z_step * float(depth_pass + 1)

		if pocket_kind == Operation.POCKET_KIND_THROUGH:
		    # We only need to do the exterior path to a depth of *z_plunge*:
		    code._simple_pocket_helper(cnc_corner_bsw, cnc_corner_tne, corner_radius, z,
		    r, zero, spindle_speed, feed_speed, True, rotate, tracing = tracing + 1)
		elif pocket_kind == Operation.POCKET_KIND_FLAT:
		    # Generate {paths} rectangular passes over the pocket:
		    for path in range(paths):
			# Provide a context comment:
			code._line_comment("Rectangular Path {0} of {1}".format(path + 1, paths))

			# Compute the *offset* for this path:
			offset = r + (r + t) * float((paths - 1) - path)

			# If *offset* exceeds *half_minimum_span* it will cause
			# *px1* > *px2* or *py1* > *py2*, which is bad.
			# We solve the problem by not  letting *offset*
			# exceed *half_minimum_span*:
			if offset > half_minimum_span:
			    offset = half_minimum_span

			# We need to get the tool to the starting point safely:
			rapid_move = False
			if path == 0:
			    # This is the first cut at a new depth:
			    if depth_pass == 0:
				# We are at *z_safe*, so we can move rapidly to the
				# plunge point:
				rapid_move = True
			    else:
				# We are sitting at the outer edge at the old depth,
				# and there is no material between where are now
				# and (*start_x*, *start_y*); a rapid could be too
				# fast, so we do a linear move:
				rapid_move = False
			else:
			    # We are at the previous internal path at this depth
			    # and need to move out and cut material as we go:
			    rapid_move = False

			# Mill out a pocket at level *z*:
		        code._simple_pocket_helper(cnc_corner_bsw, cnc_corner_tne, corner_radius,
			  z, tool_radius, offset, spindle_speed, feed_speed, rapid_move,
			  rotate, tracing = tracing + 1)
		else:
		    assert False, "Unknown pocket kind: {0}".format(pocket_kind)

	# Return the tool to a safe location above the material:
	code._xy_rapid_safe_z_force(feed_speed, spindle_speed)
	code._line_comment("Simple Pocket Done")

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Operation_Simple_Pocket._cnc_generate('{1}', '{2}', *, *, *, {3})".
	      format(indent, operation_simple_pocket._comment, mount._name_get(), is_last))

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

    def _cnc_generate(self,
      mount, mount_ngc_file, cnc_vrml, mount_cnc_lines, mount_cnc_stl, is_last, tracing=-1000000):
	""" *Operation_Vertical_Lathe*:
	"""

	assert False, "This code needs to be fixed"

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(mount_ngc_file, file)
	assert isinstance(cnc_vrml, VRML_Group)
	assert isinstance(mount_cnc_lines, VRML_Lines)
	assert isinstance(mount_cnc_stl, VRML_Group)
	assert isinstance(is_last, bool)
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

    def __init__(self, name, length, thickness, heights):
	""" *Parallels*:  """

	# Use *parallels* instead of *self*:
	parallels = self

	# Verify argument types:
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(length, L)
	assert isinstance(thickness, L)
	assert isinstance(heights, tuple) or isinstance(heights, list)
	zero = L()
	for height in heights:
	    assert isinstance(height, L) and height > zero 

	# Stuff values into *parallels*:
        parallels._heights = tuple(sorted(heights, reverse=True))
        parallels._length = length
	parallels._name = name
        parallels._thickness = thickness

    def _length_get(self):
        """ *Parallels: Return the length of each paralle associated with 
	    *Parallels* object (i.e. *self*.) """

	return self._length

    def _heights_get(self):
        """ *Parallels: Return the tuple of heights of the for all the parallels
	    associated with the *Parallels* object (i.e. *self*.) """

	return self._heights

    def _select(self, dz, vice, tracing=-1000000):
	""" *Parallels*: From the *Parallels* object (i.e. *self*), select a parallel
	    that will ensure that a *Part* of height *dz* will be at or just above the
	    jaw of *vice.
	"""

	# Use *parallels* instead of *vice*:
	parallels = self
	
	# Verify argument types:
	assert isinstance(vice, Vice)
	assert isinstance(dz, L)
	assert isinstance(tracing, int)
	
	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
            print("{0}=>Parallels._select('{1}', {2:i}, '{3}')".
	      format(indent, parallels._name, dz, vice._name_get()))

	# Grab the *jaw_height* from *vice*:
	jaw_volume = vice._jaw_volume_get()
	jaw_height = jaw_volume.z

	# *heights* is sorted from the largest height to the smallest height.  The loop is
	# terminated when the first parallel *height* that is too small is encountered.
	epsilon = L(inch=.000000001)
	heights = parallels._heights
	selected_height = None
	for height in parallels._heights:
	    # Use *epsilon* to deal with rounding errors:
	    if height + dz < jaw_height - epsilon:
		break
	    selected_height = height

	# Make sure we got a *selected_height*:
	assert isinstance(selected_height, L), \
	  "Could not find a parallel that works with part of height {0:i}".format(dz)

	# Wrap up any requested *tracing* and return selected *selected_height*:
	if tracing >= 0:
	    indent = ' ' * tracing
            print("{0}=>Parallels._select('{1}', {2:i}, '{3}')=>{4:i}".
	      format(indent, parallels._name, dz, vice._name_get(), selected_height))
	return selected_height

    def _thickness_get(self):
        """ *Parallels: Return the thickness of each parallel associated with the
	    *Parallels* object (i.e. *self*.) """

	return self._heights

    def _vrml_append(self, vrml, height, dy, tracing=-1000000):
	""" *Parallels*: Write out a VRML visualization of two *height* high parallels separated
	    by *dy* to *wrl_lines* for the *Parallels* object (i.e. *self*).  The VRML lines
	    are indented by *pad* spaces and have *comment* embedded at the front of the VRML.
	"""

	# Use *parallels* instead of *self*:
	parallels = self

	# Verify argument types:
	assert isinstance(vrml, VRML_Lines)
	assert isinstance(height, L)
	assert isinstance(dy, L)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Parallels._wrl_write(*, '{1}', {2:i}, {3:i})".
	      format(indent, vrml._name_get(), height, dy))

	# Grab some values from *parallels*:
	heights = parallels._heights
	length = parallels._length
	thickness = parallels._thickness
	
	# Create the first set of parallels points:
	zero = L()
	corner1 = P(zero, zero, zero)
	corner2 = P(length, -thickness, height)

	# Create the second set of parallels points:
	corner3 = P(zero, thickness - dy, zero)
	corner4 = P(length, -dy, height)
	
	# Draw the two parallels and write them out:
	color = "brown"
	vrml._box_outline(color, corner1, corner2)
	vrml._box_outline(color, corner3, corner4)

	# Wrap-up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Parallels._wrl_write(*, '{1}', {2:i}, {3:i})".
	      format(indent, vrml._name_get(), height, dy))

class Part:
    """ A *Part* specifies either an assembly of parts or a single physical part. """

    HOLE_THROUGH = 1
    HOLE_TIP = 2
    HOLE_FLAT = 3

    ANGLE0 = Angle()

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

	# The *Part* starts off with empty transforms:
	null_transform = Transform()
	root_transform = null_transform
	if isinstance(up, Part):
	    root_transform = up._root_transform
	
	ezcad = EZCAD3.ezcad
	part._axis = z_axis			# Scheduled to go away
	part._bom = BOM(part)
	part._bounding_box = Bounding_Box()
	part._color = None
	part._center = P()			# Is this used anymore
	part._cnc_suppress = False
	part._current_mount = None
	part._dx_original = zero
	part._dxf_x_offset = zero
	part._dxf_y_offset = zero
	part._dy_original = zero
	part._dz_original = zero
	part._dxf_scad_lines = []
	part._extra_start_stop_available = False # *True* when values are in extra start/stop fields
	part._extra_start_bsw = P()		# BSW extra material corner at part construct start
	part._extra_start_tne = P()		# TNE extra material corner at part construct start
	part._extra_stop_bsw = None		# BSW extra material corner at part construct end
	part._extra_stop_tne = None		# TNE extra material corner at part construct end
	part._ezcad = ezcad			# Top level *EZCAD3* object
	part._is_part = False			# *True* => physical part made; *False* => assembly
	part._is_place_only = False		# *True* => Render from *place*() calls only
	part._laser_preferred = False
	part._material = Material("plastic", "abs")
	part._mount_operations_table = {}
	part._mount_operations_list = []	# List of operations that construct part
	part._name = name			# Part name
	part._places = {}
	part._position_count = 0
	part._priority = 0
	part._program_number = -1		# Top level CNC program number for part
	part._program_numbers = []		# All program numers for part
	part._reverse_root_transform = root_transform.reverse()
	part._root_transform = null_transform
	part._scad_difference_lines = []
	part._scad_union_lines = []
	part._shop = ezcad._shop
	part._signature_hash = None
	part._stl_file_name = None
 	part._stl = None			# *STL* object associated with *part* (not assembly)
	part._stl_vrml = None			# *VRML_STL* object for this *part* only
	part._tooling_plate_translate_point = P() # Vice transform for tooling plate.
	part._tool_preferred = ""
	part._tracing = -1000000
	part._visible = True
	part._vrmls = None			# *VRML_Group* objects for this *part* and children
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

    def _bounding_box_get(self, top_surface_transform, tracing=-1000000):
	""" *Part*: Return the for the bounding box of the *Part* object (i.e. *self*)
	    as two corners that have both been normalized to be BSW and TNE after
	    being run through *top_surface_transform*.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(top_surface_transform, Transform)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._bounding_box_get('{1}', *)".format(indent, part._name))

	# Get the *bounding_box* from *part* and extra the two opposing corners:
	bounding_box = part._bounding_box
	bsw = bounding_box.bsw_get()
	tne = bounding_box.tne_get()

	# Now transform the two corners:
	transformed_bsw = top_surface_transform * bsw
	transformed_tne = top_surface_transform * tne

	# Now convert the transformed corners into *final_bsw* and *final_tne*:
	final_bsw, final_tne = transformed_bsw.minimum_maximum(transformed_tne)

	# Wrap up any requested *tracing* and return both *final_bsw* and *final_tne*:
	if tracing >= 0:
	    print("{0}<=Part._bounding_box_get('{1}', *)=>{2:i},{3:i}".
	      format(indent, part._name, final_bsw, final_tne))
	return final_bsw, final_tne

    def _current_mount_get(self):
	""" *Part*: Return the current mount for the *Part* object (i.e. *self*.)
	"""

	# Use *part* instead of *self*:
	part = self

	current_mount = part._current_mount
	assert isinstance(current_mount, Mount)
	return current_mount

    def _color_get(self):
	""" *Part*: Return the color associated with the *Part* object (i.e. *self*.)
	"""

	color = self._color
	assert isinstance(color, Color)
	return color

    def _cnc_manufacture(self, tracing = -1000000):
	""" *Part*: Force the generation of `.ngc` and `.wrl` path files for the *Part* object.
	"""

	# Use *part* instead of *self*
	part = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform an requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._cnc_manufacture('{1}')".format(indent, part._name))
	    trace_detail = 1

	# Process any delayed multi-mounts:
	if part.up == None:
	    ezcad = part._ezcad_get()
	    multi_mounts_list = ezcad._multi_mounts_list_get()
	    for index, multi_mounts in enumerate(multi_mounts_list):
		part == multi_mounts._part_get()
		name = multi_mounts._name_get()
		ezcad._cnc_mode = True
		ezcad._stl_mode = True
		if trace_detail >= 1:
		    print("{0}[{1}]:'{2} (Part:'{3}'".format(indent, index, name, part._name_get()))
		part._multi_mount_process(name, multi_mounts, tracing + 1)
		#part._multi_mount_process(name, multi_mounts, tracing = 5)
		ezcad._cnc_mode = False
		ezcad._stl_mode = False

	# Generate the `.ngc` files for *part*:
	ezcad = part._ezcad
	shop = part._shop_get()
	part_program_number = shop._program_base_get()
	assert part_program_number % 100 == 0
	part_program_number = part._cnc_part_generate(part_program_number, tracing + 1)
	assert part_program_number % 100 == 0
	shop._program_base_set(part_program_number)
	#part._mount = -1

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}=>Part._cnc_manufacture('{1}')".format(indent, part._name))

    def cnc_suppress(self):
        """ *Part*: Suppress CNC generation for the *Part* object (i.e. *self*.) """

	# Use *part* insted of *self*
	part = self

	# Mark *part* for CNC suppression:
        part._cnc_suppress = True

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

    def _cnc_part_generate(self, part_program_number, tracing=-1000000):
	""" *Part*: Flush out the CNC code for the *Part* object (i.e. *self*) using
	    *part_program_number* to for the first .ncg file generated.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(part_program_number, int)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._cnc_part_generate('{1}', {2})".
	      format(indent, part._name, part_program_number))
	    trace_detail = 2
	
	# If requested, show what *operations* looks like before the sort:
	mount_operations_list = part._mount_operations_list
	mount_program_number = part_program_number
	for index, mount_operations in enumerate(mount_operations_list):
	    if trace_detail >= 2:
		print("{0}[{1}]".format(indent, index))
		mount_operations._show("Mount[{0}]".format(index), tracing = tracing + 1)
	    mount_program_number = \
	      mount_operations._cnc_mount_generate(mount_program_number, tracing = tracing + 1)
	
	# Compute the *next_program_number* (which must be divisable by 100):
	next_part_program_number = mount_program_number + 99
        next_part_program_number -= next_part_program_number % 100

	# Wrap up any requested *tracing* and return the *next_part_program_number*:
	if tracing >= 0:
	    print("{0}<=Part._cnc_part_generate('{1}', prog_no={2}) =>{3}".
	      format(indent, part._name, part_program_number, next_part_program_number))
	return next_part_program_number

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
	changed = []

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
	part._priority = 0
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
		changed.append("{0}.{1}".format(part._name, attribute_name))
		if len(changed) > 0:
		    #print("here 2")
		    pass
		#if ezcad._update_count > 20:
		#    assert False, "Dimensions loop:{0}".format(changed)
		if tracing >= 0:
		    print("{0}Part._dimensions_update:{1}.{2} ({3}=>{4})". \
		      format(indent, name, attribute_name,
		      before_value, after_value))

	# Now update the *bounding_box* for *part*:
	after_bounding_box = part._bounding_box
	for sub_part in sub_parts:
	    after_bounding_box.bounding_box_expand(sub_part._bounding_box)
	if before_bounding_box != after_bounding_box:
	    changed.append("{0}.bounding_box".format(part._name))
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

    def _dxf_manufacture(self, tracing=-1000000):
	""" *Part*: Generate a dxf file for the *Part* object (i.e. *self*) if appropriate.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._dxf_manufacture('{1}'):CNC".format(indent, part._name))
	    trace_detail = 2

	# This code may be needed for DXF mode:
	ezcad = part._ezcad
	shop = ezcad._shop_get()
	code = shop._code_get()
	if code._dxf_content_avaiable():
	    # We do have a .dxf file to write.  Open *dxf_file*:
	    dxf_directory = ezcad._dxf_directory_get()
	    assert isinstance(dxf_directory, Directory)

	    # Write out the dxf file content:
	    with dxf_directory._write_open("{0}.dxf".format(part._name)) as dxf_file:
		# Output the .dxf file headers:
		assert isinstance(dxf_file, file)
		dxf_file.write("0\nSECTION\n2\nHEADER\n")
		#original = True
		original = False
		if original:
		    # Old original code (this kind of worked for some unknown reason):
		    dxf_file.write("9\n$DIMAUNITS\n70\n1\n")
		    dxf_file.write("9\n$INSUNITS\n70\n1\n")
		    dxf_file.write("9\n$LUNITS\n70\n0\n")
		    dxf_file.write("9\n$MEASUREMENT\n70\n0\n")
		else:
		    # New code that matches the DXF file spec. for metric:
		    dxf_file.write("9\n$DIMAUNITS\n70\n0\n")
		    dxf_file.write("9\n$INSUNITS\n70\n4\n")
		    dxf_file.write("9\n$LUNITS\n70\n4\n")
		    dxf_file.write("9\n$MEASUREMENT\n70\n1\n")
		dxf_file.write("0\nENDSEC\n")
		dxf_file.write("0\nSECTION\n2\nENTITIES\n")

		# Output the body of the *dxf_file*:
		code._dxf_write(dxf_file)

		# Close out *dxf_file*:
		dxf_file.write("0\nENDSEC\n0\nEOF\n")

	if tracing >= 0:
	    print("{0}<=Part._dxf_manufacture('{1}'):CNC".format(indent, part._name))

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


    def xxx_extra_start_fetch(self, top_surface_transform, tracing=-1000000):
	""" *Part*: Return the starting extra material corners for the *Part* object,
	    (i.e. *self*) where the two corners have been transformed using
	    *top_surface_transform*.  The two returned corners will be respectively
	    at BSW and TNE in the transformed space.
	"""
	
	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(top_surface_transform, Transform)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._extra_start_fetch('{1}', *)".format(indent, part._name))

	# Make sure we have actually got valid values available:
	assert part._extra_start_stop_available, \
	  "Extra start/stop not set for Part '{0}'".format(part._name)

	# Fetch the two extra start corners from *part*:
	extra_start_bsw = part._extra_start_bsw
	extra_start_tne = part._extra_start_tne

	# Verify that both *extra_start_bsw* and *extra_start_tne* have been set:
	assert isinstance(extra_start_bsw, P)
	assert isinstance(extra_start_tne, P)

	# *transform* the the corners:
	transformed_bsw = top_surface_transform * extra_start_bsw
	transformed_tne = top_surface_transform * extra_start_tne

	# Compute the *final_bsw* and *final_tne*:
	final_bsw, final_tne = transformed_bsw.minimum_maximum(transformed_tne)

	# Wrap up any requested *tracing* and return both *final_bsw* and *final_tne*:
	if tracing >= 0:
	    print("{0}<=Part._extra_start_fetch('{1}', *)=>{2:i},{3:i}".
	      format(indent, part._name, final_bsw, final_tne))
	return final_bsw, final_tne

    def xxx_extra_start_store(self,
      reverse_top_surface_transform, extra_start_bsw, extra_start_tne, tracing=-1000000):
	""" *Part* Store *extra_start_bsw* and *extra_start_tn* into the *Part* object
	    (i.e. *self*) after pushing both arguments through *reverse_top_surface_transform*.
	"""

	# Use *part* instead of *self*:
	part = self
	
	# Verify argument types:
	assert isinstance(reverse_top_surface_transform, Transform)
	assert isinstance(extra_start_bsw, P)
	assert isinstance(extra_start_tne, P)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Part._extra_start_store('{1}', *, {2:i}, {3:i})".
	      format(indent, part._name, extra_start_bsw, extra_start_tne))

	# We only get to set the extra start values once per *part*:
	assert part._extra_start_bsw == None
	assert part._extra_start_tne == None

	# *transform* the corners:
	transformed_bsw = reverse_top_surface_transform * extra_start_bsw
	transformed_tne = reverse_top_surface_transform * extra_start_tne

	# Compute the *final_bsw* and *final_tne*:
	final_bsw, final_tne = transformed_bsw.minimum_maximum(transformed_tne)

	# Fetch the two corners from *part*:
	part._extra_start_bsw = final_bsw
	part._extra_start_tne = final_tne

	# Set *extra_start_stop_available* only if *BOTH* fields have been initialzied:
	part._extra_start_stop_available = isinstance(part._extra_stop_bsw, P)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}=>Part._extra_start_store('{1}', *, {2:i}, {3:i})".
	      format(indent, part._name, extra_start_bsw, extra_start_tne))

    def xxx_extra_stop_fetch(self, top_surface_transform, tracing=-1000000):
	""" *Part*: Return the stopping extra material corners for the *Part* object,
	    (i.e. *self*) where the two corners have been transformed using
	    *top_surface_transform*.  The two returned corners will be respectively
	    at BSW and TNE in the transformed space.  The stopping extra material
	    gets smaller as contour and face surfacing operations are performd on *part*.
	"""
	
	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(top_surface_transform, Transform)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._extra_stop_fetch('{1}', *)".format(indent, part._name))

	# Fetch the two extra stop corners from *part*:
	extra_stop_bsw = part._extra_stop_bsw
	extra_stop_tne = part._extra_stop_tne

	# Verify that both *extra_stop_bsw* and *extra_stop_tne* have been set:
	assert isinstance(extra_stop_bsw, P)
	assert isinstance(extra_stop_tne, P)

	# *transform* the to corners:
	transformed_bsw = top_surface_transform * extra_stop_bsw
	transformed_tne = top_surface_transform * extra_stop_tne

	# Compute the *final_bsw* and *final_tne*:
	final_bsw, final_tne = transformed_bsw.minimum_maximum(transformed_tne)

	# Wrap up any requested *tracing* and return both *final_bsw* and *final_tne*:
	if tracing >= 0:
	    print("{0}<=Part._extra_stop_fetch('{1}', *)=>{2:i},{3:i}".
	      format(indent, part._name, final_bsw, final_tne))
	return final_bsw, final_tne

    def xxx_extra_stop_store(self,
      reverse_top_surface_transform, extra_stop_bsw, extra_stop_tne, tracing=-1000000):
	""" *Part* Store *extra_stop_bsw* and *extra_stop_tne* into the *Part* object
	    (i.e. *self*) after pushing both arguments through *reverse_top_surface_transform*.
	    The stopping extra material gets smaller as contour and face surfacing operations
	    are performd on *part*.
	"""

	# Use *part* instead of *self*:
	part = self
	
	# Verify argument types:
	assert isinstance(reverse_top_surface_transform, Transform)
	assert isinstance(extra_stop_bsw, P)
	assert isinstance(extra_stop_tne, P)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Part._extra_stop_store('{1}', *, {2:i}, {3:i})".
	      format(indent, part._name, extra_stop_bsw, extra_stop_tne))

	# It is a really bad idea to set the extra stop values before the extra start values:
	assert isinstance(part._extra_start_bsw, P)
	assert isinstance(part._extra_start_tne, P)

	# *transform* the to corners:
	transformed_bsw = reverse_top_surface_transform * extra_stop_bsw
	transformed_tne = reverse_top_surface_transform * extra_stop_tne

	# Compute the *final_bsw* and *final_tne*:
	final_bsw, final_tne = transformed_bsw.minimum_maximum(transformed_tne)

	# Fetch the two corners from *part*:
	part._extra_stop_bsw = final_bsw
	part._extra_stop_tne = final_tne

	# Set *extra_start_stop_available* only if *BOTH* fields have been initialzied:
	part._extra_start_stop_available = isinstance(part._extra_start_bsw, P)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}=>Part._extra_stop_store('{1}', *, {2:i}, {3:i})".
	      format(indent, part._name, extra_stop_bsw, extra_stop_tne))

    def _ezcad_get(self):
	""" *Part*: Return the *EZCAD* object from the *Part* object (i.e. *self*)
	"""

	return self._ezcad

    def fasten(self, comment, fastener, select, tracing=-1000000):
	""" *Part*: Use *fastener* to drill a hole from the *Part* object (i.e. *self*).
	    *select* is one of "thread", "close", or "free" to specify the desired hole
	    diameter for the hole.  *comment* is used for debugging.
	"""

	# Use *part* instead of *self*:
        part = self

	# Verify argument types:
        assert isinstance(comment, str)
        assert isinstance(fastener, Fastener)
	assert isinstance(select, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.fasten('{1}', '{2}', '{3}')".
	      format(indent, part._name, comment, fastener._name_get()))
	    trace_detail = 2

	# Perform the fasten operation:
	fastener._fasten(comment, part, select, tracing = tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.fasten('{1}', '{2}', '{3}')".
	      format(indent, part._name, comment, fastener._name_get()))

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
		    # Do not return 0, since it trigger division by 0 errors on the first pass:
		    return L(inch=0.000000000000001)
	    elif name.endswith("_m"):
		if first_update:
		    return Material()
	    elif name.endswith("_o"):
		if first_update:
		    return None
	    elif name.endswith("_s"):
		if first_update:
		    return ""
	    elif name.endswith("_t"):
		if first_update:
		    return Transform()
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
	raise AttributeError("Part '{0}' instance has no attribute named '{1}'".
	  format(self._name, name))

    def _manufacture(self, ezcad, tree_depth, tracing=-1000000):
	""" *Part*: Visit the *Part* object (i.e. *self*) and all is the children *Part*'s
	    and perform any manufacturing steps.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(ezcad, EZCAD3)
	assert isinstance(tree_depth, int) and tree_depth >= 0
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing < 0 and part._tracing >= 0:
            tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._manufacture:('{1}', *, {2})".format(indent, self._name, tree_depth))

	# Make sure that *up* is correct:
	up = part.up
	if tree_depth == 0:
	    assert up == None, \
	      "Root Part '{0}' should not have parent (it has one)".format(part._name)
	else:
	    assert up != None, \
	      "Non-Root Part '{0}' does not have a parent (it should have one)".format(part._name)

	# Make sure *part* is connected ot *ezcad*:
	part._ezcad = ezcad
	part_name = part._name
	shop = ezcad._shop_get()
	code = shop._code_get()

	# For debugging, keep the part stack in *ezcad* up to date:
	ezcad._parts_stack_push(part)

	# First manufacture any child *Part*'s:
	for attribute_name in dir(part):
	    if not attribute_name.startswith("_") and attribute_name.endswith("_"):
		child_part = getattr(part, attribute_name)
		assert isinstance(child_part, Part), \
		  "{0}.{1} is not a Part".format(part.name, attribute_name)
		child_part._manufacture(ezcad, tree_depth + 1, tracing + 1)

		# Update the Bill of Materials for *part* with *child_part* BOM:
		part._bom._merge(child_part._bom)

	# Now run *construct* for this *part*
	if tracing >= 0:
	    print("{0}==>Part.construct('{1}') ****************".format(indent, part_name))
	ezcad._stl_mode = True
	ezcad._cnc_mode = True
	part._priority = 0
	if tracing >= 0:
	    part._tracing = tracing + 1
	    part.construct()
	    part._tracing = tracing
	else:
	    part.construct()
	ezcad._stl_mode = False
	ezcad._cnc_mode = False
	if tracing >= 0:
	    print("{0}<==Part.construct('{1}') ****************".format(indent, part_name))

	# Force the `.scad` and `.stl` files in the `scad` and `stl` directories to be generated:
	part._stl_manufacture(tracing = tracing + 1)

	# Force the `.wrl` files in the `wrl` directory to be generated:
	part._wrl_manufacture(tracing = tracing + 1)

	# Force the `.ngc` and `.wrl` files for the `ngc` directory to be generated:
	part._cnc_manufacture(tracing = tracing + 1)

	# Force the `.dxf` files for the `dxf` directory to be generated:
	part._dxf_manufacture(tracing = tracing + 1)

	# For debugging, keep the part stack in *ezcad* up to date:
	ezcad._parts_stack_pop()

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part._manufacture:('{1}', *, {2})".format(indent, self._name, tree_depth))

    def _material_get(self):
	""" *Part*: Return the matrial associated with the *Part* object (i.e. *self*.)"""

	return self._material

    def _mount_operations_lookup(self, mount_name, tracing=-100000):
	""" *Part*: Return the *Mount_Operations* objected named *mount_name* from the
	    *Part* object (i.e. *self*).
	"""
	
	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(mount_name, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._mount_operations_lookup('{1}', '{2})".
	      format(indent, part._name, mount_name))

	# Lookup *
	mount_operations_table = part._mount_operations_table
	assert mount_name in mount_operations_table, \
	  "Could not find '{0}' in operations table of Part '{1}'; should be one of {2}".format(
	  mount_name, part._name, sorted(mount_operations_table.keys()))
	mount_operations = mount_operations_table[mount_name]

	# Wrap up any requested *tracing* and return *mount_operations*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._mount_operations_lookup('{1}', '{2}')=>*".
	      format(indent, part._name, mount_name))
	return mount_operations

    def _mount_register(self, mount, tracing=-1000000):
        """ *Part*: Register *mount* as a new CNC mount for the *Part* object (i.e. *self*).
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
            indent = ' ' * tracing
	    print("{0}=>Part._mount_register('{1}', '{2}')".
	      format(indent, part._name, mount._name_get()))

	# Grab some values from *part* and *mount*:
	mount_operations_table = part._mount_operations_table
	mount_operations_list = part._mount_operations_list
	mount_name = mount._name_get()

	# Make sure we have no duplicates:
	assert not mount_name in mount_operations_table, \
	  "Mount name '{0}' previously used for Part '{1}'!".format(mount_name, part._name)

	# Create a new *mount_operations* object:
	assert mount_name not in mount_operations_table
	mount_operations = Mount_Operations(mount_name)
	mount_operations_table[mount_name] = mount_operations
	mount_operations_list.append(mount_operations)
	part._current_mount = mount

	# Wrap up any requested *tracing* and return *mount_operations*:
	if tracing >= 0:
	    print("{0}<=Part._mount_register('{1}', '{2}')".
	      format(indent, part._name, mount._name_get()))
	return mount_operations

    def _name_get(self):
	""" *Part*: Return the name of the *Part* object (i.e. *self*). """

	return self._name

    def name_get(self):
	""" *Part*: Return the name of the *Part* object (i.e. *self*). """

	return self._name

    def _operation_append(self, operation, tracing=-1000000):
	""" Part*: Append *operation* to the operations list in the *Part* object (i.e. *self*). """

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(operation, Operation)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._operation_append('{1}', '{2}')".
 	      format(indent, part._name, operation._name_get()))

	# Append *operation* to *operations*:
	mount_operations_list = part._mount_operations_list
	assert len(mount_operations_list) > 0, \
 	  "No mount is active for Part '{0}'".format(part._name)
	mount_operations = part._mount_operations_list[-1]
	mount = part._current_mount
	assert isinstance(mount, Mount)
	mount_operations._append(part._current_mount, operation)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part._operation_append('{1}', '{2}')".
 	      format(indent, part._name, operation._name_get()))

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

    def place_fasten(self, comment, transform, fastener, select, tracing=-1000000):
	""" *Part*: Fasten *fastener* to the *Part* object (i.e. *self*) after applying
	    *transform* to move *fastener* to a new location.  "select* must be one of
	    "thread", "close", or "free" to specify the hole diameter.  *comment* is used
	    for debugging.
	"""

	# Use *part* instead of *self*:
        part = self

	# Verify argument types:
        assert isinstance(comment, str)
	assert isinstance(transform, Transform)
        assert isinstance(fastener, Fastener)
	assert isinstance(select, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.place_fasten('{1}', '{2}', *, '{3}')".
	      format(indent, part._name, comment, fastener._name_get()))

	fastener._fasten(comment, part, select, transform=transform, tracing = tracing + 1)


	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.place_fasten('{1}', '{2}', *, '{3}')".
	      format(indent, part._name, comment, fastener._name_get()))

    def plate(self, name, dx, dy, dz, first_letter):
	""" *Part*: Create and return a *Plate* object using the *Part* object (i.e. *self*)
	    register the global plate name.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(dx, L)
	assert isinstance(dy, L)
	assert isinstance(dz, L)
	assert isinstance(first_letter, str) and len(first_letter) == 1

	# Make sure that the returned *plate* is created exacty once:
	ezcad = part._ezcad
	plate = ezcad._plate_lookup(name)
	if plate == None:
	    # Create and register the *plate*:
	    plate = Plate(name, dx, dy, dz, first_letter)
	    ezcad._plate_register(plate)

	assert isinstance(plate, Plate)
	return plate

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
	""" *Part*: Return the priority of the *Part* (i.e. *self*.) """

	return self._priority

    def _program_number_get(self):
	""" *Part*: Return the program number associated with *Part* (i.e. *self*.)
	"""

	return self._program_number

    def _program_numbers_append(self, program_number):
        """ *Part*: Append *program_number* to the program_numbers list associated with
	    the *Part* object (i.e. *self*.)
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(program_number, int) and program_number > 0, \
	  "program_number={0}".format(program_number)

	# Stuff *program_number* into *multi_mount*:
	part._program_numbers.append(program_number)

    def _program_numbers_get(self):
	""" *Multi_Mount*: Return the program numbers list associated with the *Multi_Mount*
	    object (i.e. *self*.)
	"""

	return self._program_numbers

    def _stl_file_name_get(self):
	""" *Part*: Return the STL file name associated with the *Part* object (i.e. *self*.)
	"""

	return self._stl_file_name

    def _stl_get(self, tracing=-1000000):
	""" *Part*: Return the *STL* object associated with the *Part* object (i.e. *self*.)
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._stl_get('{1}')".format(indent, part._name))

	# Read in the *stl* object associated with *part*:
	stl = part._stl
	if stl == None:
	    stl = STL(part, tracing = tracing + 1)
	    part._stl = stl

	# Wrap up any requested *tracing* and return *stl:
	if tracing >= 0:
	    print("{0}<=Part._stl_get('{1}')=>*".format(indent, part._name))
	return stl

    def _stl_vrml_get(self, tracing=-1000000):
	""" *Part*: Return the *VRML_Triangles* object for the *Part* object (i.e. *self*).
	    This only work if the *Part* object is not an assembly.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._stml_vrml_get('{1}')".format(indent, part._name))

	# Verify that *part* is not an assembly:
	assert part._is_part

	# Get the *stl* object:
	stl_vrml = part._stl_vrml
	if stl_vrml == None:
	    # It does not exist yet, so create it:
	    stl = part._stl_get(tracing = tracing + 1)	
	    triangles = stl._triangles_get()
	    color_name = part._color._name_get()
	    stl_vrml = VRML_Triangles(part._name, color_name, colors, tracing = tracing + 1)
	    part._stl_vrml = stl_vrml

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part._stml_vrml_get('{1}')=>'{2}'".
	      format(indent, part._name, stl_vrml._name_get()))

    def _stl_manufacture(self, tracing=-1000000):
	""" *Part*: Force the generation of the `.stl` file for the *Part* object (i.e. *self*).
	"""

	# Use *part* instead of *self*:	
	part = self

	# Verify argument types:
        assert isinstance(tracing, int)

	# Grab some values from *part*:
        part_name = part._name
	part_tracing = part._tracing
	part_scad_union_lines = part._scad_union_lines
	part_scad_difference_lines = part._scad_difference_lines
	ezcad = part._ezcad

	# Perform any requested *tracing*:
	if tracing < 0 and part_tracing >= 0:
	    tracing = part_tracing
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._stl_manufacture('{1}')".format(indent, part_name))
	    trace_detail = 1

	# Write out the module:
	lines = []
	lines.append("module {0}() {{".format(part_name))

	# Get the difference() followed union():
	lines.append("  difference() {")
	lines.append("    union() {")

	# Append the *part_scan_union_lines*:
	lines.extend(part_scad_union_lines)

	# Close off union():
	lines.append("    } // union")

	# Output *scad_difference_lines*:
	lines.extend(part_scad_difference_lines)

	# Close off difference():
	lines.append("  } // difference")

	# Close off the module:
	lines.append("} // module")
	lines.append("")

	# Call the module we just produced:
	lines.append("{0}();".format(part_name))
	lines.append("")

	# Write *lines* out to *scad_file_name*:
	scad_directory = ezcad._scad_directory_get()
	scad_file_name = "{0}.scad".format(part_name)
	signature = scad_directory._lines_write(scad_file_name, lines, tracing + 1)
	part._signature = signature
	if trace_detail >= 1:
	    print("{0}'{1}' written out".format(indent, scad_file_name))

	# Figure out if we can reuse any previously generated `.stl` file:
	stl_file_name = "{0}_{1}.stl".format(part_name, signature)
	part._stl_file_name = stl_file_name
	stl_directory = ezcad._stl_directory_get()
	if stl_directory._exists(stl_file_name, tracing + 1):
	    # Do any requested *trace_detail*:
            if trace_detail >= 1:
                print("{0}stl_file '{1}' exists".format(indent, stl_file_name))
	else:
	    # Run openscad to generate the `.stl` file:
	    stl_path = stl_directory._path_get()
	    stl_full_path = os.path.join(stl_path, stl_file_name)
	    if trace_detail >= 1:
	        print("{0}stl_full_path='{1}'".format(indent, stl_full_path))

	    scad_path = scad_directory._path_get()
	    scad_full_path = os.path.join(scad_path, scad_file_name)
	    if trace_detail >= 1:
	        print("{0}scad_full_path='{1}'".format(indent, scad_full_path))

	    # Execute `openscad` to generate the `.stl` file:
	    ignore_file = open("/dev/null", "w")
	    command = [ "openscad", "-o", stl_full_path, scad_full_path ]
	    if trace_detail >= 1:
	        print("{0}subprocess command='{1}".format(indent, command))
	    subprocess.call(command, stderr=ignore_file) 
	    ignore_file.close()

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part._stl_manufacture('{1}')".format(indent, part._name))

    def _stl_vrml_get(self, tracing=-1000000):
        """ *Part*: Return the *VRML* object that renders the *Part* object (i.e. *self*).
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}Part._stl_vrml_get('{1}')".format(indent, part._name))

	# See whether or not we need to create a *VRML* object:
	stl_vrml = part._stl_vrml
	if stl_vrml == None:
	    # Grap the *STL* object and its associated *stl_triangles*:
	    stl = part._stl_get()
	    stl_triangles = stl._triangles_get()

	    # Create the *stl_vrml* object:
	    part_name = part._name
	    color_name = part._color._name_get()
	    stl_vrml = VRML_Triangles(part_name, color_name, stl_triangles, tracing + 1)
	    part._stl_vrml = stl_vrml

	# Wrap up any requested *tracing* and return *stl_vrml*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Part._stl_vrml_get('{1}')=>'{2}'".
	      format(indent, part._name, stl_vrml._name_get()))
	return stl_vrml

    def _tools_dowel_pin_search(self, preferred_diameter, tracing=-1000000):
	""" *Part*: Find and return a *Tool_Dowel_Pin* object that has as diameter that
	    is closest to *preferred_diameter*. """

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(preferred_diameter, L)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._tools_dowel_pin_search('{1}', {2:i})".
	      format(indent, self._name, preferred_diameter))
	    trace_detail = 1

	# Search for *dowel_pin*:
	deep_tracing = -1000000
	if trace_detail >= 3:
	    deep_tracing = tracing + 1

	zero = L()
	dowel_pin = part._tools_search(Tool_Dowel_Pin._match,
	  preferred_diameter, zero, "dowel_pin", tracing = deep_tracing)
	assert isinstance(dowel_pin, Tool_Dowel_Pin)

	# Wrap-up and requested *tracing* and return *dowel_pin*:
	if tracing >= 0:
	    tool_name = "NONE"
	    if dowel_pin != None:
		tool_name = dowel_pin._name_get()
	    print("{0}<=Part._tools_dowel_pin_search('{1}', {2:i})".
	      format(indent, self._name, preferred_diameter))
	return dowel_pin

    def _tools_drill_search(self, diameter, maximum_z_depth, tracing=-1000000):
	""" *Part*: Return a drill tool for the *Part* object (i.e. *self*) that has 
	    a diameter of *diameter* and go to a depth of at least *maximum_z_depth*.
	"""

	# Verify argument types:
	assert isinstance(diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(tracing, int)
	
	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._tools_drill_search('{1}', {2:i}, {3:i})".
	      format(indent, self._name, diameter, maximum_z_depth))
	    trace_detail = 3

	# Seach for a *drill_tool* that matches our requirements:
	drill_tool = self._tools_search(Tool_Drill._match,
	  diameter, maximum_z_depth, "drill", tracing = tracing + 1)
	if trace_detail >= 2:
	    print("{0}is_drill_tool={1}".format(indent, isinstance(drill_tool, Tool_Drill)))

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    tool_name = "NONE"
	    if isinstance(drill_tool, Tool_Drill):
		tool_name = drill_tool._name_get()
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
	assert isinstance(maximum_z_depth, L) and maximum_z_depth >= zero
	assert isinstance(from_routine, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	detail_level = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._tools_end_mill_search('{1}', {2:i}, {3:i}, '{4}')".
	      format(indent, self._name, maximum_diameter, maximum_z_depth, from_routine))
	    detail_level = 3

	# Search for a matching *end_mill_tool*:
	search_tracing = -1000000
	if detail_level >= 3:
	    search_tracing = tracing + 1
	end_mill_tool = self._tools_search(Tool_End_Mill._match,
	  maximum_diameter, maximum_z_depth, from_routine, search_tracing)

	# Wrap up any requested *tracing* and return *end_mill_tool*:
	if tracing >= 0:
	    tool_name = "NONE"
            if end_mill_tool != None:
		tool_name = end_mill_tool._name_get()
	    print("{0}<=Part._tools_end_mill_search('{1}', {2:i}, {3:i}, '{4}') => '{5}'".format(
	      indent, self._name, maximum_diameter, maximum_z_depth, from_routine, tool_name))
	return end_mill_tool

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
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._tools_mill_drill_side_search('{1}', {2}, {3})".
	      format(indent, self._name, maximum_diameter, maximum_z_depth))
	    trace_detail = 1

	# Search for a viable *end_mill_tool*:
	deep_tracing = -1000000
	if trace_detail >= 3:
	    deep_tracing = tracing + 1
	end_mill_tool = self._tools_search(Tool_Mill_Drill._mill_drill_side_match,
	  maximum_diameter, maximum_z_depth, "mill drill side", deep_tracing)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    tool_name = "NONE"
	    if end_mill_tool != None:
		tool_name = end_mill_tool._name_get()
	    print("{0}<=Part._tools_mill_drill_side_search('{1}', {2}, {3}) => {4}".
	      format(indent, self._name, maximum_diameter, maximum_z_depth, tool_name))

	return end_mill_tool


    # FIXME: tool searching should be part of the *Shop* class:
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

	# Search for a viable *mill_drill_tool*:
	mill_drill_tool = self._tools_search(Tool_Mill_Drill._mill_drill_side_match,
	  maximum_diameter, maximum_z_depth, "mill drill side", tracing + 1)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    tool_name = "NONE"
	    if mill_drill_tool != None:
		tool_name = mill_drill_tool._name_get()
	    print("{0}<=Part._tools_mill_drill_tip_search('{1}', {2}, {3}) => {4}".
	      format(indent, self._name, maximum_diameter, maximum_z_depth, tool_name))

	return mill_drill_tool

    # FIXME: tool searching should be part of the *Shop* class:
    def _tools_search(self, match_routine, parameter1, parameter2, from_routine, tracing=-1000000):
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
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._tools_search('{1}', *, {2:i}, {3:i}, '{4}')".
	      format(indent, part._name, parameter1, parameter2, from_routine))
	    trace_detail = 3

	# For now grab the material of the first part and assume the
	# rest are compatible.  In reality, what we want is to grab
	# all the materials and only accept tools that work with all
	# materials.  That is too much work for now:
	shop = part._shop_get()
	part_material = part._material_get()
	#print("part_material={0}".format(part_material))

	best_tool = None
	best_priority = -1.0

	# For now, hardwire the mill speeds.  We need to actually get the
	# speeds from the *Mill* object:
	maximum_spindle = Hertz(rpm=5000.0)
	minimum_spindle = Hertz(rpm=150.0)
	#print("maximum_spindle={0:rpm}rpm (={0:rps}rps)".format(maximum_spindle))
	#print("minimum_spindle={0:rpm}rpm (={0:rps}rps)".format(minimum_spindle))

	# FIXME: This code probably belongs over in *Shop*:
	# Now search through {tools}:
	zero = L()
	delta = L(inch=.000001)
	tools = shop._tools_get()
	if trace_detail >= 3:
	    print("{0}available tools={1}".format(indent, len(tools)))

	match_tracing = -1000000
	if trace_detail >= 3:
	    match_tracing = tracing + 1
	for tool in tools:
	    # Keep track of search results for this tool:
	    tool._search_results_clear()

	    # Lookup the *speed_range* for *material* and *tool*:
	    tool_material = tool._material_get()
	    speed_range = shop._surface_speeds_lookup(part_material, tool_material, match_tracing)
	    speed_range_ok = isinstance(speed_range, Speed_Range)

	    # Always log *surface_speeds_ok*:
	    tool._search_results_append(speed_range_ok,
	      "Surface speeds are acceptable for part_material '{0}' and tool material '{1}'".
	      format(part_material, tool_material))

	    if trace_detail >= 3:
		print("{0}Tool:'{1} speed_range_ok={2}'".
		  format(indent, tool._name_get(), speed_range_ok))

	    if speed_range_ok:
		# See whether *tool* has a chance of being a match:
		#if debug:
		#    print("    speed_range={0:.0F}\n".format(speed_range))
		priority = match_routine(tool, parameter1, parameter2, from_routine, match_tracing)
		if priority >= 0.0:
		    # Tool is an acceptable match:
		    if trace_detail >= 2:
			print("{0}Tool: '{1}' priority:{2}".
			  format(indent, tool._name_get(), tool._diameter_get()))

		    tool_preferred = part._tool_preferred
		    if tool_preferred != "":
			if trace_detail >= 3:
			    print("{0}tool='{1}' preferred='{2}'".
			      format(indent, tool._name, tool_preferred))
			if tool._name == tool_preferred:
			    priority = priority + 100.0
			    if trace_detail >= 3:
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
			tool_diameter = tool._diameter_get()

			#print("low_speed={0:F}fpm(={0:i}ips)".format(low_speed))
			#print("high_speed={0:F}fpm(={0:i}ips)".format(high_speed))
			#print("Tool diameter={0:i}in".format(tool._diameter))
			#print("maximum_spindle={0:rpm}".format(maximum_spindle))
			#print("minimum_spindle={0:rpm}".format(minimum_spindle))

			pi = 3.14159265358979323846
			desired_spindle = low_speed / (tool_diameter * pi)
			#print("desired_spindle={0:rpm}rpm".format(desired_spindle))

			# Deal with spindle speed limits:
			actual_spindle = desired_spindle
			if desired_spindle > maximum_spindle:
			    # Speed is to fast, clamp it to *maximum_spindle*:
			    actual_spindle = maximum_spindle
			assert actual_spindle >= minimum_spindle,				\
			  "Actual spindle ({0:rpm}) is less than minimum spindle({1:rpm})".	\
			  format(actual_spindle, minimum_spindle)
			#print("actual_spindle={0:rpm}rpm".format(actual_spindle))

			# Record everything we figured out back into *tool*:
			tool._spindle_speed_set(actual_spindle)
			chip_load = L(inch=0.001)
			flutes_count = tool._flutes_count_get()

			# Temporarily store the *feed_speed*, *spindle_speed*  and *priority*
			# in the *tool*.  These values are only good until the end of this
			# routine.  They will possibly get overridden during a subsequent
			# tool search:
			actual_spindle_frequency = actual_spindle.frequency()
			#print("actual_spindle_frequency={0}".format(actual_spindle_frequency))
			feed_speed = Speed(mm_per_sec=
			  ((chip_load * float(flutes_count)).millimeters() *
			  actual_spindle_frequency))
			#print("chip_load={0}".format(chip_load))
			#print("flutes_count={0}".format(flutes_count))
			#print("Set feed_speed for tool {0} to {1:I}".
			#  format(tool._name, feed_speed))
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
		# Check {tool} to see if it is acceptable:
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
	trace_detail = -1
	if tracing >= 0:
	    trace_detail = 1
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
	if trace_detail >= 1:
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
	if ezcad._stl_mode:
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
	    red, green, blue = color._red_green_blue_get()
	    alpha = color._alpha_get()
	    union_lines.append("{0}color([{1}, {2}, {3}, {4}])".
	      format(pad, red, green, blue, alpha))
	    union_lines.append("{0}cube([{1:m}, {2:m}, {3:m}], center=true); // {4}".
	      format(pad, x2 - x1, y2 - y1, z2 - z1, comment))

	if ezcad._cnc_mode:
	    # Add a block record to the *part* BOM (Bill Of Materials):
	    part._bom._block_pending(part, material, corner1, corner2)

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

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Part.Cylinder(c='{1}', p='{2}', m='{3}', clr={4}, sd={5:i}, ed={6:i}," +
	      " s={7:i}, e={8:i}, s={9}, sa={10:d}, w='{11}', f='{12}')").
	      format(indent, part._name, comment, material, color, start_diameter, end_diameter,
              start, end, sides, sides_angle, welds, flags))
	    trace_detail = 2

	comment = comment.replace(' ', '_')

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
	for x in [-maximum_radius, maximum_radius]:
	    for y in [-maximum_radius, maximum_radius]:
		for z in [zero, -thickness]:
		    point = reverse_transform * P(x, y, z)
		    bounding_box.point_expand(point)

	if self._ezcad._stl_mode:
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
	    print(("{0}<=Part.Cylinder(c='{1}', p='{2}', m='{3}', clr={4}, sd={5:i}, ed={6:i}," +
	      " s={7:i}, e={8:i}, s={9}, sa={10:d}, w='{11}', f='{12}')").
	      format(indent, part._name, comment, material, color, start_diameter, end_diameter,
              start, end, sides, sides_angle, welds, flags))

    def dowel_pin(self, comment, dowel_point, plunge_point, preferred_diameter, tracing=-1000000):
	""" *Part*: Request that a dowel pin be used to align the *Part* object (i.e. *self*).
	    The dowel pin will come down at *plunge_point* and move to *dowel_point* and
	    then back.  Only the x and y values of *plunge_point* and *dowel_point* are used.
	    *dowel_point* should be aligned with the actual material point to pushed against.
	    This operation is for experts.  Most people should be using either *vice_mount()*
	    with either the 'l' or 'r' flag set or *tooling_plate_mount()*.  *comment* will
	    show up in any generated G code.  *preferred_diameter* specifies the preferred
	    diameter of the dowel pin.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types
	assert isinstance(comment, str)
	assert isinstance(dowel_point, P)
	assert isinstance(plunge_point, P)
	assert isinstance(preferred_diameter, L)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.dowel_pin('{1}', '{2}', {3:i}, {4:i})".
	      format(indent, part._name, comment, dowel_point, plunge_point))
	    trace_detail = 1

	# Only generate the operation in CNC generate mode:
	ezcad = part._ezcad
	if ezcad._cnc_mode:
	    # Find the *tool_dowel_pin* to use :a
	    tool_dowel_pin = part._tools_dowel_pin_search(preferred_diameter, tracing = tracing + 1)
	    assert isinstance(tool_dowel_pin, Tool_Dowel_Pin), "No dowel pin tool found"

	    # Grab some values from *tool_dowel_pin*:
	    diameter = tool_dowel_pin._diameter_get()
	    feed_speed = tool_dowel_pin._feed_speed_get()
	    maximum_z_depth = tool_dowel_pin._maximum_z_depth_get()
	    spindle_speed = tool_dowel_pin._spindle_speed_get()

	    # This chunk of code figures how deep the part is in the mount position.
	    # It relies on the fact that we know where machine origin lands on the top
	    # surface and where the center of the *part* is from the bounding box:
	    mount = part._current_mount
	    top_surface_transform = mount._top_surface_transform_get()
	    center_point = part.c
	    cnc_center_point = top_surface_transform * center_point
	    half_part_dz = cnc_center_point.length()
	    part_dz = 2 * half_part_dz

	    # Now create the final *cnc_dowel_point* and *cnc_plunge_point*:
	    # *tip_depth* is a total kludge:
	    mount_translate_point = mount._mount_translate_point_get()
	    tip_z = mount_translate_point.z - part_dz
	    cnc_dowel_point = P(dowel_point.x, dowel_point.y, tip_z)
	    cnc_plunge_point = P(plunge_point.x, plunge_point.y, tip_z)
	    if trace_detail >= 1:
		print("{0}dowel_pin('{1}' cnc_dowel_point={2:i} cnc_plunge_point={3:i}".
		  format(indent, part._name, cnc_dowel_point, cnc_plunge_point))

	    # Create the *operation_dowel* pin and append it to the *part* operations list:
	    dowel_pin_order = Operation.ORDER_DOWEL_PIN
	    operation_dowel_pin = Operation_Dowel_Pin(part,
	      comment, 0, tool_dowel_pin, dowel_pin_order, None, feed_speed, spindle_speed,
	      diameter, cnc_dowel_point, cnc_plunge_point, tracing + 1)
	    self._operation_append(operation_dowel_pin)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.dowel_pin('{1}', '{2}', {3:i}, {4:i})".
	      format(indent, part._name, comment, dowel_point, plunge_point))

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
	    (with a render color of *color*).   *contours* is a list of *Contour* objects
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
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.extrude('{1}', {2}, {3}, *, {4:i}, {5:i}, {6:i}, {7:i}, {8:d})".
	      format(indent, comment, material, color, start, start_extra, end, end_extra, rotate))
	    trace_detail = 3

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

	# Do a little tracing:
	if trace_detail >= 1:
	    print("{0}extrude_direction={1:m} start_extra_delta={2:i} end_extra_delta={3:i}".
	      format(indent, extrude_direction, start_extra_delta, end_extra_delta))
	    print("{0}height={1:i}".format(indent, height))

	# Compute the transform that will place *start* on the X/Y plane with *start*-*end*
	# pointing in the negative Z axis direction:
	top_surface_transform = \
	  Transform.top_surface("extrude", start, end, rotate, tracing = tracing + 1)
	if trace_detail >= 2:
	    print("{0}top_surface_transform={1:v}".format(indent, top_surface_transform))

	# Now compute the transforms that will take a contour in the X/Y plane and map
	# them to appropate planes at *start* and *end*:
	reverse_top_surface_transform = top_surface_transform.reverse()
	if trace_detail >= 2:
	    print("{0}reverse_top_surface_transform={1:v}".
	      format(indent, reverse_top_surface_transform))

	# Extract the *outer_contour*:
	outer_contour = contours[0]
	outer_contour._project(top_surface_transform, tracing + 1)

	# Now we can expand the *bounding_box* by visiting each *bend* in *bends*:
	bounding_box = part._bounding_box
	bends = outer_contour._bends_get()
	bottom_offset = P(zero, zero, -height)
	for index, bend in enumerate(bends):
	    # Compute the mapped bend points:
	    top_bend_point = bend._point_get()
	    bottom_bend_point = top_bend_point + bottom_offset
	    if trace_detail >= 3:
		print("{0}[{1}]: top_bend_point={2:i}".
		  format(indent, index, top_bend_point))
		print("{0}[{1}]: bottom_bend_point={2:i}".
		  format(indent, index, bottom_bend_point))

	    transformed_top_bend_point =    reverse_top_surface_transform * top_bend_point
	    transformed_bottom_bend_point = reverse_top_surface_transform * bottom_bend_point
	    if trace_detail >= 3:
		print("{0}[{1}]:  ..transformed_top_bend_point={2:i}".
		  format(indent, index, transformed_top_bend_point))
		print("{0}[{1}]:  ..transformed_bottom_projected_point={2:i}".
		  format(indent, index, transformed_bottom_bend_point))

	    # Now expand *bounding_box*:
	    bounding_box.point_expand(transformed_top_bend_point)
	    bounding_box.point_expand(transformed_bottom_bend_point)
	if trace_detail >= 1:
	    print("{0}bounding_box={1:i}".format(indent, bounding_box))

	# Now render the *contours* in openscad:
	part._is_part = True
	ezcad = part._ezcad
	if ezcad._stl_mode:
	    # Trace STL mode:
	    if trace_detail >= 1:
		print("{0}STL mode entered".format(indent))

	    # Set *lines* to collect OpenSCAD commands.  Due to the natur the of OpenSCAD
	    # language syntax, these lines are appended to *lines in "reverse" order.
	    lines = part._scad_union_lines
	    pad = ' ' * 6

	    # Use *reverse_top_surface_transform* to map the origin to *start*:
	    lines.append("{0}// extrude('{1}', {2}, {3}, {4:m}, {5:m}, {6:m}, {7:m}, {8:d})".
	      format(pad, comment, material, color, start, start_extra, end, end_extra, rotate))
	    reverse_top_surface_transform._scad_lines_append(lines, pad)
	    if trace_detail >= 2:
		print("{0}reverse_top_surface_transform={1:v}".
		  format(indent, reverse_top_surface_transform))

	    # Put top of extrusion at *origin*:
	    lines.append(
	      "{0}translate([ 0, 0, {1:m} ]) // Put top of extrusion at origin".
	      format(pad, -height))
	    lines.append(
	      "{0}linear_extrude(height = {1:m}, convexity=10)".format(pad, height))

	    # Start appending the OpenSCAD `polygon` command to *lines*:
	    lines.append("{0}polygon(".format(pad))
	    lines.append("{0}  convexity = 10,".format(pad))

	    # Output all of the points used by each *contour* in *contours*.
	    # This looks like:
	    # 
	    #        points = [
            #          [x1, y1], ... , [xN, yN], // Contour 0
	    #          [x1, y1], ... , [xN, yN], // Contour 1
	    #          ...
            #          [x1, y1], ... , [xN, yN], // Contour M
            #        ],
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

	    # Output the contour paths with point indices, staring with the outer contour and
	    # followed by any inner contours.  This looks like:
	    #         paths = [
	    #           [0, 1, ..., N], // Contour 0
            #           [N+1, ..., M],  // Contour 1
	    #           ...
            #           [...]           // Contour M
            #         ]
            lines.append("{0}  paths = [".format(pad))
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

	    # Wrap up the `polygon` command:
            lines.append("{0}  ]".format(pad))
            lines.append("{0}); // polygon".format(pad))

	    # Wrap up any *tracing*:
	    if trace_detail >= 1:
		print("{0}STL exited".format(indent))

	# Wrap up any requested *tracing*;
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

	if self._ezcad._cnc_mode:
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
	""" *Part*: Make part visible/invisible.
	"""

	assert isinstance(invisible, bool)
	self._visible = not invisible

    def place(self, name, sub_part, transform):
	""" *Part*: Add a placement of *Part* using *transform* to the *Part* object (i.e. *self*.)
	    The *name* must be unique for *sub_part*.  *name* is also used for debugging.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(sub_part, Part)
	assert isinstance(transform, Transform)

	# Create *place* and insert to *places* indexed by *name*.  The reason why *_places*
        # is a dictionary instead of a list is because because this is called from a
        # construct method which gets called multiple times:
	places = part._places
	if name in places:
	    # We have previously created a *place* for *name*:
	    place = places[name]

	    # Update the *transform* in *place* since it could have changed.
	    place._transform_set(transform)
	else:
	    # This is the first time we have specified this *place* for *name*:
	    place = Place(name, sub_part, transform)
	    places[name] = place

    def place_only_set(self):
	""" *Part*: Set the *Part* object so that it is only rendered from calls to the
	    *Part*.*place*() routine.
	"""

	self._is_place_only = True

    def process(self, ezcad, tracing = -1000000):
	""" *Part*: Perform all the manufacturing processing for a *Part*
	    object (i.e. *self*) and all its sub-*Part*'s.
	"""

	# Verify argument types:
	assert isinstance(ezcad, EZCAD3)
	assert isinstance(tracing, int)

	# Use *part* instead of *self*:
	part = self

	# For debugging, set *debug* to *True*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.process('{1}')".format(indent, part._name))

	# Do the dimensions propogate phase:
	ezcad._dimensions_mode = True
	ezcad._update_count = 0
	changed = ["foo"]
	while len(changed) > 0:
	    if ezcad._update_count > 20:
                print("Dimensions loop cycle:")
		for change_index, change in enumerate(changed):
		    change_indent = ""
		    if change.endswith("bounding_box"):
			change_indent = "        "
		    print("[{0}]:{1}{2}".format(change_index, change_indent, change))
		assert False, "Infinite dimensions loop?"
	    # Find all the child *Part*'s:
	    print("Dimensions update {0}".format(ezcad._update_count + 1))
	    changed = part._dimensions_update(ezcad, -1000000)
	    assert isinstance(changed, list)
	    #changed = part._dimensions_update(ezcad, 0)
	    #part._bounding_box_check(0)
	    print("Part.process: {0} dimension(s) changed\n".format(len(changed)))
	    ezcad._update_count += 1
	ezcad._dimensions_mode = False

	print("*************Part:{0} bounding_box={1:i}".format(self._name, self._bounding_box))

	# Recursively manufacture all of the *Part* objects starting with *part*:
	ezcad._parts_stack_reset()
	part._manufacture(ezcad, 0, tracing + 1)
	ezcad._parts_stack_reset()

	part._bom._material_summary()

	# Clean up any directories:
	ezcad._dxf_directory._clean_up(tracing = tracing + 1)
	ezcad._ngc_directory._clean_up(tracing = tracing + 1)
	ezcad._scad_directory._clean_up(tracing = tracing + 1)
	ezcad._stl_directory._clean_up(tracing = tracing + 1)
	ezcad._wrl_directory._clean_up(tracing = tracing + 1)

	# Process any *plate* objects:
	ezcad._plates_process()

	# Write out the tools summary:
	shop = ezcad._shop_get()
	shop._tools_summary_write()

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.process('{1}')".format(indent, part._name))

    def re_multi_mount(self, multi_mounts, old_name, new_name, tracing=-1000000):
	""" *Part*: Using the *Part* object (i.e. *self*), make a copy of *multi_mounts*
	    substituting *old_name* with *new_name* for all of the mount names.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
        assert isinstance(multi_mounts, Multi_Mounts)
	assert isinstance(old_name, str) and not ' ' in old_name
	assert isinstance(new_name, str) and not ' ' in new_name
	assert isinstance(tracing, int)

	# Perform any requested tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.re_multi_mount('{1}', '{2}', '{3}', '{4}')".
	      format(indent, part._name, multi_mounts._name_get(), old_name, new_name))

	# Only do this operation in *cnc_mode*:
	ezcad = part._ezcad_get()
	if ezcad._cnc_mode:
	    # Make a copy of *multi_mounts* substituting in *mount_name*:
	    new_multi_mounts = \
	      multi_mounts._name_replace_copy(old_name, new_name, tracing = tracing + 1)

	    # Now do the *multi_mount* operation on *new_multi_mounts*:
	    part.multi_mount(new_multi_mounts, None, False, 0)

	# Wrap up any requested tracing*:
	if tracing >= 0:
	    print("{0}<=Part.re_multi_mount('{1}', '{2}', '{3}', '{4}')".
	      format(indent, part._name, multi_mounts._name_get(), old_name, new_name))

    def rectangular_contour(self, comment, radius, start_corner=0, tracing=-1000000):
	""" *Part*: Peform a rectangular exterior contour of the *Part* object (i.e. *self*)
	    where the corners have a radius of *radius*.  *comment* will show up in generated
	    G-code.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(radius, L)
	assert isinstance(start_corner, int) and 0 <= start_corner < 4
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.rectangular_contour('{1}', '{2}', {3:i})".
	      format(indent, part._name, comment, radius))
	    trace_detail = 1

	# We only need to do stuff in STL or CNC mode:
	ezcad = part._ezcad_get()
	if ezcad._stl_mode or ezcad._cnc_mode:
	    # We need to figure out where the four corners are in 3D space.  To do this we first
	    # figure out where everything landed after *top_surface_transform* is applied.  After
	    # that we can reverse the points back into 3D space.  So first we compute the bounding
	    # box coordinates of the *part* after it has been transformed  with the
            # *top_surface_transform* touching the machine origin:

	    # Get the *top_transform* and its reverse from *mount*:
	    mount = part._current_mount
	    assert isinstance(mount, Mount)
	    top_transform = mount._top_surface_transform_get()
	    reverse_top_transform = top_transform.reverse()
	    if trace_detail >= 2:
		print("{0}top_transform={1:v}".format(indent, top_transform))
		print("{0}reverse_top_transform={1:v}".format(indent, reverse_top_transform))

	    # Now grap the bounding box after *top_surface_transform* has been applied:
	    top_bsw, top_tne = \
	      part._bounding_box_get(top_transform, tracing = tracing + 1)
	    top_x1 = top_tne.x
	    top_x0 = top_bsw.x
	    top_y1 = top_tne.y
	    top_y0 = top_bsw.y
	    top_z1 = top_tne.z
	    top_z0 = top_bsw.z
	    if trace_detail >= 2:
		print("{0}top_bsw={1:i} top_tne={2:i}".format(indent, top_bsw, top_tne))
		print("{0}top_x0={1:i} top_x1={2:i}".format(indent, top_x0, top_x1))
		print("{0}top_y0={1:i} top_y1={2:i}".format(indent, top_y0, top_y1))
		print("{0}top_z0={1:i} top_z1={2:i}".format(indent, top_z0, top_z1))

	    # Now we use *reverse_top_surface_transform* to map the points back to 3D space and
	    # append each *corner* to *corners* in a clockswise direction to force climb milling:
	    corners = []
	    corners.append(reverse_top_transform * P(top_x0, top_y0, top_z1))
	    corners.append(reverse_top_transform * P(top_x0, top_y1, top_z1))
	    corners.append(reverse_top_transform * P(top_x1, top_y1, top_z1))
	    corners.append(reverse_top_transform * P(top_x1, top_y0, top_z1))

	    # Rotate *corners* by *start_corner* using some fun Python slice magic:
	    corners = corners[start_corner:] + corners[:start_corner]
	    if trace_detail >= 2:
		print("{0}corner1={1:i}".format(indent, corners[0]))
		print("{0}corner2={1:i}".format(indent, corners[1]))
		print("{0}corner3={1:i}".format(indent, corners[2]))
		print("{0}corner4={1:i}".format(indent, corners[3]))

	    # Create the *exterior_contour*:
	    exterior_contour = Contour("Exterior Contour")
	    for index, corner in enumerate(corners):
		exterior_contour.bend_append("Corner {0}".format(index), corner, radius)
	
	    # Now figure out the *start* and *stop* point back in *part* coordinates:
	    top_start = P(top_x0, (top_y0 + top_y1)/2, top_z1)
	    top_stop  = P(top_x0, (top_y0 + top_y1)/2, top_z0 - L(inch=0.020))
	    start = reverse_top_transform * top_start
	    stop  = reverse_top_transform * top_stop
	    if trace_detail >= 2:
		print("{0}top_start={1:i}".format(indent, top_start))
		print("{0}top_stop={1:i}".format(indent, top_stop))
		print("{0}start={1:i}".format(indent, start))
		print("{0}stop={1:i}".format(indent, stop))

	    # Figure out the bounding box dimensions after *top_transform*:
	    top_dx = top_x1 - top_x0
	    top_dy = top_y1 - top_y0

	    # FIXME: Dealing with extra material in contours is broken!!!
	    part_mount = part._current_mount_get()
	    extra_bsw, extra_tne = part_mount._extra_stop_get(tracing = tracing + 1)
	    top_extra_bsw = top_transform * extra_bsw
	    top_extra_tne = top_transform * extra_tne


	    # *extra* is the material around *part*.  It used to compute how many passes
	    # are needed to machine out the contour:
	    top_extra_dx = top_extra_tne.x - top_extra_bsw.x
	    top_extra_dy = top_extra_tne.y - top_extra_bsw.y
	    extra = (top_extra_dx - top_dx).maximum(top_extra_dy - top_dy)
	    zero = L()
	    if extra <= zero:
		extra = L(mm=1.0)

	    # Trim X/Y coordinates of *extra_stop_bsw* and *extra_stop_tne* down to contour
	    # boundaries.  Leave the Z coordinates alone since no facing operation has occurred.
	    # Stuff the result back into *part*:
	    extra_stop_bsw = reverse_top_transform * P(top_x0, top_y0, top_extra_bsw.z)
	    extra_stop_tne = reverse_top_transform * P(top_x1, top_y1, top_extra_tne.z)
	    part_mount._extra_stop_set(extra_stop_bsw, extra_stop_tne)

	    # Now, actually perform the contour operation:
	    part.contour("Exterior Contour",
	      exterior_contour, start, stop, extra/2, "t", tracing = tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.rectangular_contour('{1}', '{2}', {3:i}".
	      format(indent, part._name, comment, radius))

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

    def round_pocket(self, comment, diameter, start, stop, flags, tracing=-1000000):
	""" *Part*: Manufacture a round pocket in the *Part* object (i.e. *self*) that
	    is *diameter* round, starts at *start* and ends at *end*.  Set *flags* to 't'
	    if the round pocket is supposed to go through *part*.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(diameter, L)
	assert isinstance(start, P)
	assert isinstance(stop, P)
	assert isinstance(flags, str)
	assert isinstance(tracing, int)

	# Perform an requested *tracing*:
	trace_detail = -1
	deep_tracing = -1000000
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
            print("{0}=>Part.round_pocket('{1}', {2:i}, {3:i} {4:i}, '{5}'".
	      format(indent, comment, diameter, start, stop, flags))
	    trace_detail = 2
	    #deep_tracing = tracing + 1
	if trace_detail >= 1:
	    print("{0}start-stop={1:i}".format(indent, start - stop))

	# Cut the round round pocket in the object:
	ezcad = part._ezcad_get()
	if ezcad._stl_mode:
	    lines = part._scad_difference_lines
	    pad = ' ' * 4
	    lines.append("{0}// round_pocket('{1}, {2:i}, {3:i}, {4:i}, '{5}')".
	      format(pad, comment, diameter, start, stop, flags))

	    degrees0 = Angle()
	    hole_axis = (start - stop).normalize()
	    extra = L(inch=0.025).millimeters()
	    #assert isinstance(extra, float)
	    new_start = start + hole_axis * extra
	    new_stop = stop
	    if 't' in flags:
		new_stop = stop - hole_axis * extra
	    if trace_detail >= 1:
		print("{0}start={1:i} new_stop={2:i}".format(indent, start, stop))
		print("{0}hole_axis={1:m} new_start={2:i} new_stop={3:i}".
		  format(indent, hole_axis, new_start, new_stop))
	    part._scad_cylinder(comment, Color("black"), diameter, diameter,
	      new_start, new_stop, lines, pad, 24, degrees0, tracing = deep_tracing)

	# Perform the CNC:
	if ezcad._cnc_mode:
	    hole_z_depth = start.distance(stop)
	    end_mill_tool = part._tools_end_mill_search(diameter,
	      hole_z_depth, "round_pocket", deep_tracing)
	    assert end_mill_tool != None, \
	      "Could not find end mill for '{0}' round pocket".format(comment)
	    if trace_detail >= 1:
		print("{0}hole_z_depth={1:i} end_mill='{2}'".
	          format(indent, hole_z_depth, end_mill_tool._name_get()))

	    # Process the 't' flag to make the pocket go all the way through *part*:
	    new_start = start
	    new_stop = stop
	    hole_kind = Part.HOLE_FLAT
	    extra = L(inch=0.025).millimeters()
	    if 't' in flags:
		# Extend *stop* so that it ends .250inch further "down":
		hole_axis = (start - stop).normalize()
		new_stop = stop - hole_axis * extra
		hole_kind = Part.HOLE_THROUGH
	    if trace_detail >= 1:
		print("{0}start={1:i} new_start={2:i} hole_axis={3:m}".
		  format(indent, start, new_start, hole_axis))
		print("{0}stop={1:i} new_stop={2:i}".
		  format(indent, stop, new_stop))

	    # Create *operation_round_pocket* and stuff it into the *part* operation list:
	    zero = L()
	    operation_round_pocket = Operation_Round_Pocket(part, comment, 0, end_mill_tool,
	      Operation.ORDER_END_MILL_ROUND_POCKET, None, diameter, zero, hole_kind, new_start,
	      new_stop, end_mill_tool._feed_speed_get(), end_mill_tool._spindle_speed_get(),
	      tracing = tracing + 1)
	    part._operation_append(operation_round_pocket)

	# Wrap-up any requested *tracing*:
	if tracing >= 0:
            print("{0}<=Part.round_pocket('{1}', {2:i}, {3:i} {4:i}, '{5}'".
	      format(indent, comment, diameter, start, stop, flags))

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

    def top_face(self, comment, tracing=-1000000):
	""" *Part*: Remove extra material from the top surface of the *Part* object (i.e. *self*.)
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.top_face('{1}')".format(indent, comment))
	    trace_detail = 2

	# Only do stuff if in *ezcad* is in *cnc_mode*:
	ezcad = part._ezcad_get()
	if ezcad._cnc_mode:
	    # Get the current mount:
	    current_mount = part._current_mount
	    assert isinstance(current_mount, Mount), \
	      "Is part '{0}' mounted in vice?".format(part._name)

	    # Convert the extra bounding material bounding box via *top_surface_transform*:
	    extra_stop_bsw, extra_stop_tne = current_mount._extra_stop_get()
	    top_transform = current_mount._top_surface_transform_get()
	    top_extra_stop_bsw, top_extra_stop_tne = \
	      (top_transform * extra_stop_bsw).minimum_maximum(top_transform * extra_stop_tne)
	    if trace_detail >= 2:
		print("{0}Part='{1}'".format(indent, part._name))
		print("{0}part.tne={1:i}".format(indent, part.tne))
		print("{0}part.bsw={1:i}".format(indent, part.bsw))
		print("{0}current_mount='{1}'".format(indent, current_mount._name_get()))
		print("{0}extra_stop_tne={1:i}".format(indent, extra_stop_tne))
		print("{0}extra_stop_bsw={1:i}".format(indent, extra_stop_bsw))
		print("{0}top_transform=\n{1:m}".format(indent, top_transform))
		print("{0}top_extra_stop_tne={1:i}".format(indent, top_extra_stop_tne))
		print("{0}top_extra_stop_bsw={1:i}".format(indent, top_extra_stop_bsw))

	    # Grap the bounding box oriented using *top_surface_transform(*:
	    top_bounding_box_bsw, top_bounding_box_tne = part._bounding_box_get(top_transform)
	    if trace_detail >= 2:
		print("{0}top_bounding_box_tne={1:i}".
		  format(indent, top_bounding_box_tne))
		print("{0}top_bounding_box_bsw={1:i}".
		  format(indent, top_bounding_box_bsw))

	    # Compute *remove_dz* which is the amount mow off the top:
	    remove_dz = top_extra_stop_tne.z - top_bounding_box_tne.z
	    if trace_detail >= 2:
		print("{0}remove_dz={1:i}".format(indent, remove_dz))

	    # Search for an *end_mill_tool* that can do the job:
	    detail_tracing = -1000000
	    if trace_detail >= 3:
		detail_tracing = tracing + 1
	    maximum_diameter = L(inch=1.000)
	    end_mill_tool = self._tools_end_mill_search(maximum_diameter,
	      remove_dz, "simple_pocket", detail_tracing)
	    assert end_mill_tool != None, \
	      "Could not find a end mill to mill off top of part {0}".format(part._name)
	    if trace_detail >= 1:
		 print("{0}end_mill_tool='{1}'".format(indent, end_mill_tool._name_get()))

	    # Grab some values out of *end_mill_tool*:
	    end_mill_diameter = end_mill_tool._diameter_get()
	    end_mill_radius = end_mill_diameter / 2
	    end_mill_feed_speed = end_mill_tool._feed_speed_get()
	    end_mill_spindle_speed = end_mill_tool._spindle_speed_get()

	    # We need the compute the corners of a "simple pocket" that will remove the top face.
	    zero = L()
	    end_mill_offset = P(end_mill_radius, end_mill_radius, zero)
	    top_corner1 = P(top_extra_stop_bsw.x, top_extra_stop_bsw.y, top_bounding_box_tne.z) - \
	      end_mill_offset
	    top_corner2 = top_extra_stop_tne + end_mill_offset

	    # Note that we need to map back from CNC coordinates to actual coordinates
	    # for the *Operation_Simple_Pocket*() below:
	    reverse_top_transform = top_transform.reverse()
	    final_corner1, final_corner2 = (reverse_top_transform * top_corner1). \
	      minimum_maximum(reverse_top_transform * top_corner2)
	    if trace_detail >= 2:
		print("{0}top_corner2={1:i}".format(indent, top_corner2))
		print("{0}top_corner1={1:i}".format(indent, top_corner1))
		print("{0}final_corner2={1:i}".format(indent, final_corner2))
		print("{0}final_corner1={1:i}".format(indent, final_corner1))

	    # Now construct the *operation_simple_packet* and add it to *part*:
	    operation_order = Operation.ORDER_END_MILL_SIMPLE_POCKET
	    pocket_kind = Operation.POCKET_KIND_FLAT
	    operation_simple_pocket = Operation_Simple_Pocket(part, comment, 0,
	      end_mill_tool, operation_order, None, end_mill_feed_speed, end_mill_spindle_speed,
	      final_corner1, final_corner2, end_mill_radius, end_mill_radius, pocket_kind,
	      Part.ANGLE0, tracing = tracing + 1)
	    #operation_simple_pocket._tracing = 5
	    part._operation_append(operation_simple_pocket)

	    # Remember to update the extra material to note that the top face has been removed:
	    new_top_extra_stop_bsw = top_extra_stop_bsw
	    new_top_extra_stop_tne = top_extra_stop_tne - P(zero, zero, remove_dz)
	    new_extra_stop_bsw, new_extra_stop_tne = \
	      (reverse_top_transform * new_top_extra_stop_bsw).minimum_maximum(
	      reverse_top_transform * new_top_extra_stop_tne)
	    if trace_detail >= 2:
		print("{0}new_extra_stop_tne={1:i}".format(indent, new_extra_stop_tne))
		print("{0}new_extra_stop_bsw={1:i}".format(indent, new_extra_stop_bsw))
	    current_mount._extra_stop_set(new_extra_stop_bsw, new_extra_stop_tne)

	    # Make sure that this CNC operation does not get merged with any other operations
            # that follow it:
	    #part.cnc_fence()

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.top_face('{1}')".format(indent, comment))

    def tracing_enable(self):
	""" *Part*: Enable tracing for the *Part* object (i.e. *self*).
	"""

	self._tracing = 0

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

    def _wrl_manufacture(self, tracing = -1000000):
        """ *Part*: Cause the `.wrl` file associated with the *Part* object (i.e. *self*)
	    to be generated.
	"""

	# Use *part* instead of *self*:
        part = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested tracing:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._wrl_manufacture('{1}')".format(indent, part._name))
	    trace_detail = 3

	# OK. This probably the best place to explain some of the gymnastics required to get
	# VRML to work correctly with the *Place* objects associated with a *Part*.  The way
	# this all works is that the Part is usually output as a VRML Group with the following
	# format:
	#
	#        DEF *part_name* Group [
        #         children [
	#          # Children *Part*'s typically only occur in assemblies
	#          *child part1 VRML Node*
	#          ...
        #          *child partN VRML Node*
	#
	#          # STL for *part* (usually *NOT* an assembly):
	#          Shape {
        #           # VRML for Part (usually a bunch of triangles):
	#          }
        #
	#          # Placements typically only occur in assemblies:
	#          Transform { ... children[ USE place1_part_name] }
	#          ...
	#          Transform { ... children[ USE placeN_part_name] }
        #         ]
        #        }
        #
        # The key thing to understand is that the ability to do placements relies on the
	# VRML DEF and USE technology.  While the VRML reference manual does not specify any
	# constraints on the order of DEF and USE, there are many implementations that require
	# the "DEF name" to occur before an "USE name".  This is definitely true with the
        # *view3dscene* program.  For this reason, all children *Part* VRML is placed first in the
	# file to force the child definitions to occur first.  All of the placement transforms
	# occur at the end of the *Part* VRML.
	#
	# In order to work, the Part tree constrains the *Place* object to only be able to
	# reference *Part*'s that are lower in the *Part* tree.  Otherwise, it is possible
	# to generate a USE for a part before it is DEF'ed.
	# 
	# The next issue concerns the *Part.no_automatic_placement*() method which sets the
	# *_is_place_only* flag for a *Part*.  When this flat is set, we do not want to
	# automatically display the *Part* with an empty transform.  The trick for doing
	# this is to use the magic of a VRML Switch statement:
	#
	#         Switch {
	#          choice [
        #           Group DEF *part_name* {
	#            # Rest of VRML for *part_name* here...
	#           }
	#          ]
	#         }
        #
	# Please read the VRML reference manual to find out more about the VRML Switch Node.
	# Basially, by not specifying a `whichChoice` clause, the body of the Switch Node
	# is processed to record the DEF, but it is not automatically rendered.  The node
        # only rendered by an appropriate USE reference.
	#	
	# Lastly, we write out a `.wrl` for each node in the *Part* tree.  For *Part*'s that
	# that have the *_is_place_only* flag set to *True*, this will cause the *Part* to not
	# render.  The kludge to work around this generate an enclosing VRML Group node that
	# REF's the *Part* VRML.  Yes, this is a total kludge but it works.  Enough said.
    
	# Start by creating *part_group_vrmls* which is the workhorse VRML Group node for
	# this *part*:
	part_name = part._name
	part_group_vrmls = VRML_Group(part_name)

	# Put the children *Part*'s into *part_group_vrmls*:
	for attribute_name in dir(self):
	    if not attribute_name.startswith("_") and attribute_name.endswith("_"):
		child_part = getattr(self, attribute_name)
		assert isinstance(child_part, Part)

		# Only call *_wrl_manufacture* for *child_part* if we have not already done so:
		child_part_vrmls = child_part._vrmls
		if child_part_vrmls == None:
		    child_part._wrl_manufacture(tracing = tracing + 1)
		    child_part_vrmls = child_part._vrmls
		    assert isinstance(child_part_vrmls, VRML_Group) or \
		      isinstance(child_part_vrmls, VRML_Switch)

		# Stuff *child_part_vrmls* onto the end *part_group_vrmls*:
		part_group_vrmls._append(child_part_vrmls)

	# Create the *stml_vrml* for *part* if it is an actual part::
	if part._is_part:
	    stl = part._stl_get()
	    stl_triangles = stl._triangles_get()
	    color_name = part._color._name_get()
	    stl_vrml = VRML_Triangles(part._name, color_name, stl_triangles, tracing = tracing + 1)
	    part_group_vrmls._append(stl_vrml)

	# Perform the *places*:
	part_bounding_box = part._bounding_box
	places = part._places
	for place in places.values():
	    # Grab some values out of *place*:
	    assert isinstance(place, Place)
	    place_name = place._name_get()
	    place_part = place._part_get()
	    place_transform = place._transform_get()
	
	    if trace_detail >= 3:
		print("{0}place_name='{1}'".format(indent, place_name))
		print("{0}place_part.bsw={1:i} place_part.tne={2:i}".
		  format(indent, place_part.bsw, place_part.tne))
		print("{0}transformed place_part.bsw={1:i}".
		  format(indent, place_transform * place_part.bsw))
		print("{0}before part_bounding_box={1:i}".format(indent, part_bounding_box))

	    # Expand the *bounding_box*:
	    part_bounding_box.point_expand(place_transform * place_part.bsw)
	    part_bounding_box.point_expand(place_transform * place_part.bnw)
	    part_bounding_box.point_expand(place_transform * place_part.bne)
	    part_bounding_box.point_expand(place_transform * place_part.bse)
	    part_bounding_box.point_expand(place_transform * place_part.tsw)
	    part_bounding_box.point_expand(place_transform * place_part.tnw)
	    part_bounding_box.point_expand(place_transform * place_part.tne)
	    part_bounding_box.point_expand(place_transform * place_part.tse)

	    if trace_detail >= 3:
		print("{0}after part_bounding_box={1:i}".format(indent, part_bounding_box))

	    # Verify that *place_part* occurs "under" *part* in the *Part* tree in order to
	    # ensure that every part is DEF'ed before a subsequent USE:
	    parent_part = place_part
	    while parent_part != None:
		if parent_part == part:
		    break
		parent_part = parent_part.up
	    assert parent_part != None, "Place part '{0}' is not 'under' '{1}' in Part tree. ".\
	      format(place_part._name_get(), part._name)

	    # Now create a VRML USE for *place*:
	    use_comment = "{0}_USE".format(place_part._name_get())
	    use_vrml = VRML_Use(use_comment, place_part, tracing = tracing + 1)
	    transformed_use_vrml = place_transform._vrml(use_vrml, tracing = tracing + 1)
	    part_group_vrmls._append(transformed_use_vrml)

	# If we have a *place_only* part, we hide the definition inside of VRML Switch grouping:
	# The definition is recorded by the VRML viewer, but does not actually get shown until
	# there is a VRML USE reference somewhere else.  Frankly, this is pretty obscure stuff
	# and it causes problems below when we want to write out VRML file:
	is_place_only = part._is_place_only
	if is_place_only:
	    switch_vrmls = VRML_Switch("{0}_Switch".format(part_name), -1)
	    switch_vrmls._append(part_group_vrmls)
	    part_group_vrmls = switch_vrmls
	
	# Stuff *part_group_vrmls* into *part*:
	part._vrmls = part_group_vrmls

	# OK this is ugly.  When the *part* is marked *is_place* only and is the top level
	# VRML file, the VRML Switch grouping causes nothing to be shown.  Bummer!  The solution
	# is to wrap everything in a Group that accesses the part definition:
	final_vrmls = part_group_vrmls
	if is_place_only:
	    use_comment = "{0}_USE".format(part_name)
	    use_vrml = VRML_Use(use_comment, part, tracing = tracing + 1)
	    group_vrml = VRML_Group("TOP_LEVEL")
	    group_vrml._append(part_group_vrmls)
	    group_vrml._append(use_vrml)
	    final_vrmls = group_vrml

	# Write *final_vrmls* out to the file *wrl_file_name*:
	ezcad = part._ezcad
	if not isinstance(part, Fastener):
	    wrl_directory = ezcad._wrl_directory_get()
	    wrl_file_name = "{0}.wrl".format(part._name)
	    wrl_file_text = final_vrmls._text_pad(0)
	    with wrl_directory._write_open(wrl_file_name, tracing = tracing + 1) as wrl_file:
		wrl_file.write("#VRML V2.0 utf8\n")
		wrl_file.write(wrl_file_text)

	# Wrap up any requested *tracing* and return *part_groupo_vrmls*:
	if tracing >= 0:
	    print("{0}<=Part._wrl_manufacture('{1}')".format(indent, part._name))
	return part_group_vrmls

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
	    obtained from *start_point* and *end_point*.  *extra* is the amount
	    of extra material being removed. *flags* can contain the letter 'u' for
	    an upper chamfer, 'l' for a lower chamfer, and 't' for a contour that cuts
	    entirely through the part entirely.  *comment* is used in error messages
	    and any generated G-code.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(contour, Contour)
	assert isinstance(start_point, P)
	assert isinstance(end_point, P)
	assert isinstance(extra, L)
	assert isinstance(flags, str)
	for flag in flags:
	    assert flag in "ult", "Part.contour: Bad flag '{0}' in flags '{1}' for part '{2}'". \
	      format(flag, flags, part._name)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.contour('{1}, '{2}', {3:i}, {4:i}, {5:i}, '{6}')".
	      format(indent, part._name, comment, start_point, end_point, extra, flags))
	    trace_detail = 3

	# *tracing_plus_one* is only set to a positive valued when *trace_detail* is high:
	tracing_plus_one = -1000000
	if trace_detail >= 3:
	    tracing_plus_one = tracing + 1

	# Before we do anything else, we need to update the bounding box for *part*:
	# FIXME: Why???!!!  Removing for now!!!
	#bounding_box = part._bounding_box
	#if trace_detail >= 1:
	#    print("{0}[C0]bb={1:i}".format(indent, part._bounding_box))
	#contour._bounding_box_expand(bounding_box)
	#if trace_detail >= 1:
	#    print("{0}[C1]bb={1:i}".format(indent, part._bounding_box))

	# Figure out if we are in *cnc_mode* or *stl_mode*:
	ezcad = part._ezcad
	cnc_mode = ezcad._cnc_mode and not part._cnc_suppress
	stl_mode = ezcad._stl_mode

	# Only do real work in CNC or STL mode:
	if stl_mode or cnc_mode:
	    # Perform all of the *contour* massaging:
	    mount = part._current_mount
	    assert isinstance(mount, Mount)

	    # Output the linear_extrude to *difference_lines*:
	    if stl_mode:
		# Indicate we have started the STL stuff:
		if trace_detail >= 1:
		    print("{0}STL mode started".format(indent))

		top_surface_transform = mount._top_surface_transform_get()
		transformed_start_point = top_surface_transform * start_point
		transformed_end_point   = top_surface_transform * end_point
		if trace_detail >= 2:
		    print("{0}mount='{1}'".format(indent, mount._name_get() ))
		    print("{0}top_surface_transform={1:v}".format(indent, top_surface_transform))

		contour._project(top_surface_transform, tracing_plus_one)
		contour._inside_bends_identify(tracing_plus_one)
		contour._radius_center_and_tangents_compute(tracing_plus_one)

		# Make sure we cut in a clockwise direction to force climb cutting:
		if not contour._is_clockwise:
		    contour._bends_reverse(tracing_plus_one)
		if trace_detail >= 1:
		    print("{0}is_clockwise={1}".format(indent, contour._is_clockwise))

		# We always contour a little *to_extra* on top.  We contour *bottom_extra*
		# on the bottom if the *flags* has the 't` flag for "through":
		zero = L()
		top_extra = L(mm=1.0)
		bottom_extra = zero
		if 't' in flags:
		    # Through hole flag specified:
		    bottom_extra = L(mm=1.0)
		if trace_detail >= 2:
		    print("{0}top_extra={1:i} bottom_extra={2:i}".
		      format(indent, top_extra, bottom_extra))

		# Now compute the transformed start and end points adjusted by *top_extra*
		# and *bottom_extra*:
		transformed_start_point = top_surface_transform * start_point
		transformed_end_point   = top_surface_transform * end_point
		if trace_detail >= 2:
		    print("{0}transformed_start_point={1:i}".
		      format(indent, transformed_start_point))
		    print("{0}transformed_end_point={1:i}".format(indent, transformed_end_point))
		height = transformed_start_point.z - transformed_end_point.z
		total_height = bottom_extra + height + top_extra
		if trace_detail >= 1:
		    print("{0}bottom_extra={1:i} height={2:i} top_extra={3:i}".
		      format(indent, bottom_extra, height, top_extra))
		    print("{0}total_height={1:i}".format(indent, total_height))

		# Generate the OpenSCAD commands to `linear_extrude` and translate and transform
		# the contour to the correct location.  Remember perform all OpenScad operations
		# are in "reverse" order:
		difference_lines = part._scad_difference_lines
		pad = ' ' * 4
		difference_lines.append("{0}// Contour '{1}' start={2:m} end={3:m}".
		  format(pad, contour._name_get(), start_point, end_point))
		top_surface_transform.reverse()._scad_lines_append(difference_lines, pad)
		difference_lines.append("{0}translate([0, 0, {1:m}])".
		  format(pad, transformed_end_point.z - bottom_extra))
		difference_lines.append("{0}linear_extrude(height = {1:m}) {{".
		  format(pad, total_height.absolute()))
		contour._scad_lines_polygon_append(difference_lines, pad, True, tracing_plus_one)
		difference_lines.append("{0}}} // linear_extrude".format(pad))

		# Indicate we are done with STL stuff:
		if trace_detail >= 1:
		    print("{0}STL mode ended".format(indent))

	    if trace_detail >= 1:
		print("{0}[C0]bb={1:i}".format(indent, part._bounding_box))

	    # Do any requested CNC generation:
	    if cnc_mode:
		# Indicate we have started the CNC stuff:
		if trace_detail >= 1:
		    print("{0}CNC mode started".format(indent))

		# Start with *cnc_transform*:
		cnc_transform = mount._cnc_transform_get()
		if trace_detail >= 2:
		    print("{0}mount='{1}'".format(indent, mount._name_get() ))
		    print("{0}cnc_transform={1:v}".format(indent, cnc_transform))

		# Project all of the points in *contour* onto a plane normal to *projection_axis*
		# and then rotate that plane to be the X/Y plane:
		contour._project(cnc_transform, tracing_plus_one)

		# The returned value from *_smallest_inner_diameter_compute*() will be either
		# negative (for no smallest inner corner radius), or positive for the smallest
		# inner corner radius.  Either value will find the correct end-mill or mill-drill
		# tool for searching purposes:
		zero = L()
		contour._inside_bends_identify(tracing_plus_one)
		smallest_inner_radius = contour._smallest_inner_radius_compute(tracing_plus_one)
		smallest_inner_diameter = smallest_inner_radius * 2
		contour._radius_center_and_tangents_compute(tracing_plus_one)

		cnc_start_point = cnc_transform * start_point
		cnc_end_point   = cnc_transform * end_point
		if trace_detail >= 2:
		    print("{0}cnc_start_point={1:i} cnc_stop_point={2:i}".
		      format(indent, cnc_start_point, cnc_end_point))

		# Find a *Tool* to use for edge milling.  For non-through contours, we should use
		# a mill drill in preference to an end-mill because that way we can overlap
		# with any top chamfering hole countersinking (i.e. one fewer tool change).
		# Otherwise use an end mill:
		z_stop = (cnc_start_point.z - cnc_end_point.z).absolute()
		assert isinstance(z_stop, L)
		assert z_stop >= zero
		#if z_stop > zero:
		#    print("Part '{0}' has problems".format(part._name))
		if trace_detail >= 2:
		    print("{0}z_stop={1:i}".format(indent, z_stop))
		mill_drill_tool = part._tools_mill_drill_side_search(smallest_inner_diameter,
		  z_stop, tracing = tracing_plus_one)
		have_mill_drill = isinstance(mill_drill_tool, Tool_Mill_Drill)
		end_mill_tool = part._tools_end_mill_search(smallest_inner_diameter,
		  z_stop, "contour", tracing = tracing_plus_one)
		have_end_mill = isinstance(end_mill_tool, Tool_End_Mill)
		if trace_detail >= 2:
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
			# Contour all the way through:
			do_through = True
		    elif flag == 'l':
			assert False, "Lower contour chamfers are broken"
			do_through = True
		    elif flag == 'u':
			assert False, "Upper contour chamfers are broken"
			do_through = True
		    else:
			assert False, \
			  "Unrecognized flag '{0}' in contour flags \"{1}\"".format(flag, flags)

		# Now select which *mill_tool* to use from either *mill_drill_tool* or
		# *end_mill_tool*:
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
		point_angle = Angle(deg=180)		# End-mills are flat
		tip_depth = zero			# End-mills have no tip
		if isinstance(mill_tool, Tool_Mill_Drill):
		    # Alternatively, mill drills have both a *point_angle* and a *tip_depth*:
		    tip_depth = mill_tool._tip_depth_get()

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
			if trace_detail >= 0:
			    print("{0}diameter={1:i} extra={2:i} ratio={3} depth_max={4:i}".
			      format(indent, diameter, extra, ratio, depth_maximum))

		    # Figure out how many *passes* using *depth_maximum*, *z_start*, *z_stop*,
		    # and *z_extra:
		    z_start = cnc_start_point.z
		    z_end   = cnc_end_point.z
		    z_extra = zero
		    z_depth = (z_start - z_end).absolute() + z_extra

		    if trace_detail >= 1:
			print("{0}z_depth={1:i} depth_maximum={2:i}".
		          format(indent, z_depth, depth_maximum))
		    passes = int(z_depth / depth_maximum) + 1

		    if trace_detail >= 1:
			print("{0}z_start={1:i} z_end={2:i} dmax={3:i} z_x={4:i} zd={5:i} pss={6}".
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

		    if trace_detail >= 2:
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

			#print("top_ch={0:i} tip_depth={1:i} ntr={2:i} zs={3:i} zstc={4:i}".
                        #  format(top_chamfer, tip_depth, nominal_tool_radius,
                        #  z_start, z_stop_top_chamfer))

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
		    bottom_chamfer_tool = \
		      part.tools_dove_tail_search(smallest_inner_diameter, z_stop)

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

		# Indicate we are done with CNC stuff:
		if trace_detail >= 1:
		    print("{0}CNC mode ended".format(indent))

	if tracing >= 0:
	    print("{0}<=Part.contour('{1}, '{2}', {3:i}, {4:i}, {5:i}, '{6}')".
	     format(' ' * tracing, part._name, comment, start_point, end_point, extra, flags))

    def countersink_hole(self, comment, hole_diameter, countersink_diameter,
      start, stop, flags, sides = -1, sides_angle=Angle(), tracing=-1000000):
	""" *Part*: Put a *hole_diameter* hole into the *Part* object (i.e. *self*) starting
	    at *start* to an end depth of *end*.  If *countersink_diameter* is non-zero,
	    the part will be have a 90 degree countersink of *countersink_diameter* at *start*.
	    *comment* will show up in any  generated RS-274 code.  *hole_kind* specifies the
	    kind of hole.
	"""
    
	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	zero = L()
	assert isinstance(comment, str)
	assert isinstance(hole_diameter, L) and hole_diameter >= zero
	assert isinstance(countersink_diameter, L)
	assert isinstance(start, P)
	assert isinstance(stop, P)
	assert isinstance(flags, str)
	assert isinstance(sides, int)
	assert isinstance(sides_angle, Angle)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	trace_detail = -1
	indent = ""
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.countersink_hole('{1}', '{2}', {3:i}, {4:i}, {5:i}, {6:i}, '{7}')".
	      format(indent, part._name, comment, hole_diameter, countersink_diameter,
	      start, stop, flags))
	    trace_detail = 3

	#FIXME: Why *stl_mode* with *cnc_mode*???!!!
	ezcad = part._ezcad
	cnc_mode = ezcad._cnc_mode
	stl_mode = ezcad._stl_mode
	if cnc_mode or stl_mode:
	    if trace_detail >= 0:
		print("{0}cnc_mode={1} stl_mode={2}".format(indent, cnc_mode, stl_mode))

	    # Specify ome constants:
	    zero = L()
	    degrees0 = Angle(deg=0.0)

	    # Compute *is_tip_hole*, *is_flat_hole* and *is_through_hole* from *flags*:
	    is_through_hole = 't' in flags
	    is_tip_hole = 'p' in flags
	    is_flat_hole = 'f' in flags
	    assert is_through_hole or is_tip_hole or is_flat_hole, \
	      "Specify 't' for through hole, 'p' for tip hole or 'f' for flat hole in flags"
	    assert     is_through_hole and not is_tip_hole and not is_flat_hole or \
	      not is_through_hole and     is_tip_hole and not is_flat_hole or \
	      not is_through_hole and not is_tip_hole and     is_flat_hole,   \
	      "Only specify one of 't', 'p', or 'f' for flags argument"
	
	    # Figure out how deep we want the countsink "cone hole" to go:
	    is_countersink = countersink_diameter > zero
	    if not is_countersink:
		# No countersink is required, so we will countersink to just
		# to just about 3/4 of the drill diameter:
		countersink_diameter = 1.10 * hole_diameter
	    if trace_detail >= 2:
		 print("{0}is_through_hole={1} is_tip_hole={2} is_flat_hole={3} is_countersink={4}".
		   format(indent, is_through_hole, is_tip_hole, is_flat_hole, is_countersink))

	    # Compute the radii:
	    hole_radius = hole_diameter/2
	    countersink_radius = countersink_diameter/2

	    # The hole is drilled in relation to the *top_surface_transform* from *mount*:
	    mount = part._current_mount
	    assert isinstance(mount, Mount)
	    top_surface_transform = mount._top_surface_transform_get()

	    # Compute the transfromed *bsw*, *tne*, *start*, and *stop*:
	    bsw = part.bsw
	    tne = part.tne
	    top_surface_bsw    = top_surface_transform * bsw
	    top_surface_tne    = top_surface_transform * tne
	    top_surface_start  = top_surface_transform * start
	    top_surface_stop   = top_surface_transform * stop

	    # Now we want to be sure that *start* is higher than *stop* in the Z axis.  The
	    # inversion of *start* and *stop* will happen as a side effect of the *fasten*
	    # operation.
	    if top_surface_start.z < top_surface_stop.z:
		start, stop = stop, start
		top_surface_start, top_surface_stop = top_surface_stop, top_surface_start

	    # Make sure that the drill is aligned with the Z axis:
	    slop = L(inch=.0001)
	    if not top_surface_stop.x - slop <= top_surface_start.x <= top_surface_stop.x + slop \
	      or not top_surface_stop.y - slop <= top_surface_start.y <= top_surface_stop.y + slop:
		print("Z axis not aligned for hole '{0}' in part '{1}'".format(comment, part._name))

	    # Make sure that that the drill operation does not start above *top_z*:
	    top_z    = top_surface_tne.z
	    bottom_z = top_surface_bsw.z
	    if top_z < bottom_z:
		top_z, bottom_z = bottom_z, top_z

	    # Compute the *start_z*, *hole_kind*, and *stop_z*:
	    start_z = top_surface_start.z
	    stop_z = top_surface_stop.z
	    if is_tip_hole:
		hole_kind = Part.HOLE_TIP
	    elif is_flat_hole:
		hole_kind = Part.HOLE_FLAT
		try_flat = True
		#hole_z_depth = start_z - stop_z
	    elif is_through_hole:
		hole_kind = Part.HOLE_THROUGH

		# Modifiy *stop* to reflect the new stop location:
		top_surface_stop -= P(zero, zero, L(inch=.100))
		stop_z = top_surface_stop.z
		stop = top_surface_transform.reverse() * top_surface_stop
	    else:
		assert False, "No hole kind specified"
	    hole_z_depth = start_z - stop_z
	    assert hole_z_depth >= zero, \
	      "start_z={0:i} stop_z={1:i} hole_z_depth={2:i}".format(start_z, stop_z, hole_z_depth)

	    # Print out some tracing infromation:
	    if trace_detail >= 1:
		print("{0}hole_kind={1}".format(indent, hole_kind))
	    if trace_detail >= 2 or hole_z_depth <= zero:
		print("{0}start={1:i}".format(indent, start))
		print("{0}stop={1:i}".format(indent, stop))
		print("{0}top_surface_start={1:i}".format(indent, top_surface_start))
		print("{0}top_surface_stop={1:i}".format(indent, top_surface_stop))
		print("{0}start_z={1:i}".format(indent, start_z))
		print("{0}stop_z={1:i}".format(indent, stop_z))
		print("{0}hole_z_depth={1:i}".format(indent, hole_z_depth))
		print("{0}top_surface_tne={1:i} tne={2:i}".format(indent, top_surface_tne, tne))
		print("{0}top_surface_bsw={1:i} bsw={2:i}".format(indent, top_surface_bsw, bsw))

	    success = False
	    spot_operation = None

	    # Search for useful tools:
	    deep_tracing = -1000000
	    if trace_detail >= 3:
		deep_tracing = tracing + 1
	    maximum_mill_drill_diameter = L(inch=1.0)
	    mill_drill_tool = part._tools_mill_drill_tip_search(
	      maximum_mill_drill_diameter, countersink_diameter/2, deep_tracing)
	    drill_tool = part._tools_drill_search(hole_diameter, hole_z_depth, deep_tracing)
	    end_mill_tool = part._tools_end_mill_search(hole_diameter,
	      hole_z_depth, "countersink_hole", deep_tracing)
	    end_mill_tool_is_laser = False
	    if end_mill_tool != None:
		end_mill_tool_is_laser = end_mill_tool._is_laser_get()
	    if trace_detail >= 2:
	        print("{0}mill_drill_tool={1}".format(indent, mill_drill_tool))
	        print("{0}drill_tool={1}".format(indent, drill_tool))
	        print("{0}end_mill_tool={1}".format(indent, end_mill_tool))
		print("{0}end_mill_tool_is_laser={1}".format(indent, end_mill_tool_is_laser))
		print("{0}is_through_hole={1}".format(indent, is_through_hole))
		print("{0}is_tip_hole={1}".format(indent, is_tip_hole))

	    #FIXME: We need to use *Operation_Round_Pocket* when we are generating a .dxf file!!!:
	    if (is_through_hole or is_tip_hole) and mill_drill_tool != None and \
	      drill_tool != None and not end_mill_tool_is_laser:
		# Make sure we have both a *mill_drill_tool* and a *drill_tool*:
		assert isinstance(mill_drill_tool, Tool_Mill_Drill)
		assert isinstance(drill_tool, Tool_Drill)
		if trace_detail >= 1:
		    print("{0}Use a drill".format(indent))

		# Worry about countersinking first:
		operation_countersink = None
		if drill_tool._is_center_cut_get():
		    if trace_detail >= 1:
			print("{0}Drill tool is center cut, so countersink is skipped".
			  format(indent))
		else:
		    # With the *mill_drill_tool* we can actually countersink or at least
		    # spot drill where the drill needs to go:
		    # Compute the depth of the countersink "cone hole" drill using the formula:
		    #
		    #   depth = r * tan(90 - pa/2)
		    #
		    # where r is the hole radius and pa is the point angle of the tool:
		    mill_drill_point_angle = mill_drill_tool._point_angle_get()
		    mill_drill_diameter = mill_drill_tool._diameter_get()
		    mill_drill_radius = mill_drill_diameter / 2 
		    degrees90 = Angle(deg=90.0)
		    countersink_z_depth = \
		      countersink_radius * (degrees90 - mill_drill_point_angle/2).tangent()
				
		    # Figure out where the *countersink_stop* point is:
		    countersink_comment = "{0} [countersink]".format(comment)
		    countersink_stop = start - P(zero, zero, countersink_z_depth)
		    countersink_stop_z = countersink_stop.z

		    # Create *operation_counterink* and append it to the *part* operations:
		    sub_priority = 0
		    operation_countersink = Operation_Drill(part, countersink_comment,
		      sub_priority, mill_drill_tool, Operation.ORDER_MILL_DRILL_COUNTERSINK,
		      None, countersink_diameter, Part.HOLE_TIP,
		      start, countersink_stop, True, tracing = tracing + 1)
		    part._operation_append(operation_countersink, tracing = tracing + 1)

		# Now drill the hole:
		sub_priority = 1
		operation_drill = Operation_Drill(part, comment, sub_priority, drill_tool,
		  Operation.ORDER_DRILL, operation_countersink, hole_diameter, hole_kind,
		  start, stop, False, tracing = tracing + 1)
		part._operation_append(operation_drill, tracing = tracing + 1)
		success = True
	    elif (is_flat_hole or is_through_hole) and end_mill_tool != None:
		mount = part._current_mount
		operation_round_pocket = Operation_Round_Pocket(part, comment,
		  0, end_mill_tool, Operation.ORDER_END_MILL_ROUND_POCKET, None,
		  hole_diameter, countersink_diameter, hole_kind, start, stop,
		  end_mill_tool._feed_speed_get(), end_mill_tool._spindle_speed_get(),
		  tracing = tracing + 1)
		part._operation_append(operation_round_pocket)
		success = True
    	    else:
		print(("Can't drill diameter {0:i} hole {1:i} deep in part '{2}'" +
		  " flags='{3}' md={4} d={5} em={6}").
		  format(hole_diameter, hole_z_depth, part._name, flags, mill_drill_tool != None,
		  drill_tool != None, end_mill_tool != None))

	if ezcad._stl_mode:
	    if trace_detail >= 1:
	         print("{0}STL Mode".format(indent))

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
	    line = "{0}// countersink_hole('{1}' '{2}' {3:i} {4:i} {5:i} {6:i} '{7}'".format(pad,
              part._name, comment, hole_diameter, countersink_diameter, start, stop, flags)
	    lines.append(line)
	    if trace_detail >= 2:
		print("{0}line='{1}'".format(indent, line))

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

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.countersink_hole('{1}', '{2}', {3:i}, {4:i}, {5:i}, {6:i}, '{7}')".
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
	    print("{0}=>Part._scad_cylinder('{1}', '{2}', {3:i}, {4:i}, {5:i}, {6:i}, *, {7}, {8})".
	      format(indent, self._name, comment,
	        start_diameter, end_diameter, start, end, pad, color))

	# Only do the append in STL mode:
	ezcad = self._ezcad
	if ezcad._stl_mode:

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
	    assert height > zero, \
	      "Cylinder '{0}' has no height (start={1:i} end={2:i}".format(comment, start, end)

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
	    print("{0}<=Part._scad_cylinder('{1}', '{2}', {3:i}, {4:i}, {5:i}, {6:i}, *, {7}, {8})".
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

	# Use *part* instead of *self*:
	part = self

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
	#    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>hole('{1}', '{2}', {3:i}, {4:i}, {5:i}, '{6}')".
	      format(indent, part._name, comment, diameter, start, stop, flags))

	# Perform the hole using the richer *countesink_hole* operation:
	zero = L()
	ezcad = part._ezcad_get()
	if ezcad._cnc_mode:
	    assert start.distance(stop) > zero, "start={0:i} stop={1:i}".format(start, stop)
	part.countersink_hole(comment, diameter, zero, start, stop, flags,
	   sides=sides, sides_angle=sides_angle, tracing = tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=hole('{1}', '{2}', {3:i}, {4:i}, {5:i}, '{6}')".
	      format(indent, part._name, comment, diameter, start, stop, flags))

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
	trace_detail = -1
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.lathe('{1}', '{2}', {3}, {4:i}, {5:i}, *, {6})".
	      format(indent, part._name, comment, material, color, start, end, contour, faces))
	    trace_detail = 1

	part._is_part = True
	part._color = color
	part._material = material

	ezcad = part._ezcad
	if ezcad._stl_mode:
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
	    part._top_surface_transform_set(top_surface_transform, "lathe")
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

    def no_automatic_place(self):
	""" *Part*: Disable automatic placement of the *Part* object (i.e. *self*.)
	"""

	# Mark that only *Place* objects are to be used:
	self._is_place_only = True

    #def point(self, point_path):
    #	""" Part dimensions: Return the {P} associated with {point_path}
    #	    starting from {self}. """
    #
    #	return self.value_lookup(point_path, Part.POINTS)

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
	if ezcad._stl_mode:
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
	# Fractional drills in 1/64" increments to 63/64":
	for index in range(63):
	    numerator = index + 1
	    denominator = 64
	    while number % 1 == 0:
		numerator /= 2
		denominator /= 2
	    diameter = float(numerator)/float(denominator)
	    name = "{0}/{1}".format(numerarator, denominator)
	    drill[name] = diameter

	# Create metric screw table *m*:
	m = {}

	# Fill in *m* table:
        # From: https://littlemachineshop.com/reference/TapDrillSizes.pdf
	#               (75% M, 75% I,   50% M, 50% I,   close M, close I, free M, free I)
	m["M1.5x.35"] = (1.15,  "56",    1.25,  "55",    1.60,    "1/16",    1.65,   "52")
	m["M1.6x.35"] = (1.25,  "55",    1.35,  "54",    1.70,    "51",      1.75,   "50")
	m["M1.8x.35"] = (1.45,  "53",    1.55,  "1/16",  1.90,    "49",      2.00,   "5/64")
	m["M2x.45"]   = (1.55,  "1/16",  1.70,  "51",    2.10,    "45",      2.20,   "44")
	m["M2x.40"]   = (1.60,  "52",    1.75,  "50",    2.10,    "45",      2.20,   "44")
	m["M2.2x.45"] = (1.75,  "50",    1.90,  "48",    2.30,    "3/32",    2.40,   "41")
	m["M2.5x.45"] = (2.05,  "46",    2.20,  "44",    2.65,    "37",      2.75,   "7/64")
	m["M3x.60"]   = (2.40,  "41",    2.60,  "37",    3.15,    "1/8",     3.30,   "30")
	m["M3x.50"]   = (2.50,  "39",    2.70,  "36",    3.15,    "1/8",     3.30,   "30")
	m["M3.5x.60"] = (2.90,  "32",    3.10,  "31",    3.70,    "27",      3.85,   "24")
	m["M4x.75"]   = (3.25,  "30",    3.50,  "28",    4.20,    "19",      4.40,   "17")
	m["M4x.70"]   = (3.30,  "30",    3.50,  "28",    4.20,    "19",      4.40,   "17")
	m["M4.5x.75"] = (3.75,  "25",    4.00,  "22",    4.75,    "13",      5.00,   "9")
	m["M5x1.0"]   = (4.00,  "21",    4.40,  "11/64", 5.25,    "5",       5.50,   "7/32")
	m["M5x.90"]   = (4.10,  "20",    4.40,  "17",    5.25,    "5",       5.50,   "7/32")
	m["M5x.80"]   = (4.20,  "19",    4.50,  "16",    5.25,    "5",       5.50,   "7/32")
	m["M5.5x.90"] = (4.60,  "14",    4.90,  "10",    5.80,    "1",       6.10,   "B")
	m["M6x1.0"]   = (5.00,  "8",     5.40,  "4",     6.30,    "E",       6.60,   "G")
	m["M6x.75"]   = (5.25,  "4",     5.50,  "7/32",  6.30,    "E",       6.60,   "G")
	m["M7x1.0"]   = (6.00,  "B",     6.40,  "E",     7.40,    "L",       7.70,   "N")
	m["M7x.75"]   = (6.25,  "D",     6.50,  "F",     7.40,    "L",       7.70,   "N")
	m["M8x1.25"]  = (6.80,  "H",     7.20,  "J",     8.40,    "Q",       8.80,   "S")
	m["M8x1.0"]   = (7.00,  "J",     7.40,  "L",     8.40,    "Q",       8.80,   "S")
	m["M9x1.25"]  = (7.80,  "N",     8.20,  "P",     9.50,    "3/8",     9.90,   "25/64")
	m["M9x1.0"]   = (8.00,  "O",     8.40,  "21/64", 9.50,    "3/8",     9.90,   "25/64")
	m["M10x1.50"] = (8.50,  "R",     9.00,  "T",     10.50,   "Z",       11.00,  "7/16")
	m["M10x1.25"] = (8.80,  "11/32", 9.20,  "23/64", 10.50,   "Z",       11.00,  "7/16")
	m["M10x1.0"]  = (9.00,  "T",     9.40,  "U",     10.50,   "Z",       11.00,  "7/16")
	m["M11x1.50"] = (9.50,  "3/8",   10.00, "X",     11.60,   "29/64",   12.10,  "15/32")
	m["M12x1.75"] = (10.30, "13/32", 10.90, "27/64", 12.60,   "1/2",     13.20,  "33/64")
	m["M12x1.50"] = (10.50, "Z",     11.00, "7/16",  12.60,   "1/2",     13.20,  "33/64")
	m["M12x1.25"] = (10.80, "27/64", 11.20, "7/16",  12.60,   "1/2",     13.20,  "33/64")
	m["M14x2.0"]  = (12.10, "15/32", 12.70, "1/2",   14.75,   "37/64",   15.50,  "39/64")
	m["M14x1.50"] = (12.50, "1/2",   13.00, "33/64", 14.75,   "37/64",   15.50,  "39/64")
	m["M14x1.25"] = (12.80, "1/2",   13.20, "33/64", 14.75,   "37/64",   15.50,  "39/64")
	m["M15x1.50"] = (13.50, "17/32", 14.00, "35/64", 15.75,   "5/8",     16.50,  "21/32")
	m["M16x2.0"]  = (14.00, "35/64", 14.75, "37/64", 16.75,   "21/32",   17.50,  "11/16")
	m["M16x1.50"] = (14.50, "37/64", 15.00, "19/32", 16.75,   "21/32",   17.50,  "11/16")
	m["M17x1.50"] = (15.50, "39/64", 16.00, "5/8",   18.00,   "45/64",   18.50,  "47/64")
	m["M18x2.50"] = (15.50, "39/64", 16.50, "41/64", 19.00,   "3/4",     20.00,  "25/32")
	m["M18x2.0"]  = (16.00, "5/8",   16.75, "21/32", 19.00,   "3/4",     20.00,  "25/32")
	m["M18x1.50"] = (16.50, "21/32", 17.00, "43/64", 19.00,   "3/4",     20.00,  "25/32")
	m["M19x2.50"] = (16.50, "21/32", 17.50, "11/16", 20.00,   "25/32",   21.00,  "53/64")
	m["M20x2.50"] = (17.50, "11/16", 18.50, "23/32", 21.00,   "53/64",   22.00,  "55/64")
	m["M20x2.0"]  = (18.00, "45/64", 18.50, "47/64", 21.00,   "53/64",   22.00,  "55/64")
	m["M20x1.50"] = (18.50, "47/64", 19.00, "3/4",   21.00,   "53/64",   22.00,  "55/64")

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

    def simple_pocket(self,
      comment, corner1, corner2, radius, flags, reverse=False, rotate=ANGLE0, tracing = -1000000):
	""" *Part*: Create a simple rectangular pocket in the *Part* object (i.e. *self*)
	    bounding corners of *bottom_corner* and *top_corner*, a corner radius if *radius*.
	"""

	# Use *part* instead of *self*:
	part = self

	# Check argument types:
	assert isinstance(comment, str)
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)
	assert isinstance(radius, L)
	assert isinstance(flags, str)
	assert isinstance(reverse, bool)
	assert isinstance(rotate, Angle)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	#if tracing < 0 and rotate != Part.ANGLE0:
	#    tracing = 3
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.simple_pocket('{1}', '{2}', {3:i}, {4:i}, {5:i}, '{6}')".
	      format(indent, self._name, comment, corner1, corner2, radius, flags))
	    trace_detail = 2

	# Some constants:
	zero = L()

	ezcad = self._ezcad
	if ezcad._stl_mode or ezcad._cnc_mode:
	    # *top_surface_transform* transform forces the part into the correct orientation
	    # for milling the simple pocket.  So now we need to figure out the *bsw_corner*
	    # and *tne_corner* with the part located immediately under the machine origin:
	    mount = part._current_mount
	    top_surface_transform = mount._top_surface_transform_get()
	    transformed_corner1 = top_surface_transform * corner1
	    transformed_corner2 = top_surface_transform * corner2
	    x1, x2 = transformed_corner1.x.minimum_maximum(transformed_corner2.x)
	    y1, y2 = transformed_corner1.y.minimum_maximum(transformed_corner2.y)
	    z1, z2 = transformed_corner1.z.minimum_maximum(transformed_corner2.z)
	    x_center = (x1 + x2) / 2
	    y_center = (y1 + y2) / 2
	    center = P(x_center, y_center, zero)
	    bsw_corner = P(x1, y1, z1)
	    tne_corner = P(x2, y2, z2)
	    if trace_detail >= 2:
		print("{0}bsw_corner={1:i} tne_corner={2:i}".format(indent, bsw_corner, tne_corner))

	    if ezcad._stl_mode:
		# Compute *top_transform* and transform the two corners:
		if trace_detail >= 2:
		    print("{0}Part.simple_pocket: STL_MODE started".format(indent))

		# As usual with OpenScad, everthing is output in reverse order...

		# Generate the openscad stuff:
		difference_lines = self._scad_difference_lines
		pad = ' ' * 6

		# Determine *start_extra* and *end_extra* which is the amount added to
		# the pocket start and end depth:
		start_extra = L(mm = 1.0)
		end_extra = zero
		if flags == "t":
		    end_extra = L(mm = 1.0)

		# Compute the *start* and *end* points:
		#start = P(x_center, y_center, z2 + start_extra)
		#end =  P(x_center, y_center, z1 - end_extra)

		# Output the transform that will put everything in the correct location:
		top_surface_transform_reverse = top_surface_transform.reverse()
		top_surface_transform_reverse._scad_lines_append(difference_lines, pad)

		difference_lines.append(
		  "{0}// Part.simple_pocket('{1}', '{2}', {3:i}, {4:i}, {5:i}, '{6}', {7:d})".
		  format(pad, self._name, comment, corner1, corner2, radius, flags, rotate))
		#difference_lines.append(
		#  "{0}//start={1:m} end={2:m}".format(pad, start, end))
		#difference_lines.append(
		#  "{0}//transformed_bsw_corner={1:m}".format(pad, transformed_bsw_corner))
		#difference_lines.append(
		#  "{0}//transformed_tne_corner={1:m}".format(pad, transformed_tne_corner))

		# Translate the forthcoming polygon:
		#top_surface_transform.reverse()._scad_lines_append(difference_lines, pad)
		difference_lines.append("{0}translate([0, 0, {1:m}])".
		  format(pad, -(z2 - z1) - end_extra))

		# Linear extrude the forthcoming polygon:
		difference_lines.append("{0}linear_extrude(height = {1:m}, center = false)".
		  format(pad, z2 - z1 + start_extra + end_extra))

		# Output a polygon that represents the pocket.

		# First put together a list *bend_points* which are the four pocket corners
		# moved inwards by *radius*:
		bend_points = [
		  P(x1 + radius, y1 + radius, zero).xy_rotate(center, rotate),
		  P(x1 + radius, y2 - radius, zero).xy_rotate(center, rotate),
		  P(x2 - radius, y2 - radius, zero).xy_rotate(center, rotate),
		  P(x2 - radius, y1 + radius, zero).xy_rotate(center, rotate)
		]
		if reverse:
		    bend_points.reverse()

		# This is where we the compute the various angles need to have *corner_sides* sides
		# in each pocket corner:
		corner_sides = 4
		degrees90 = Angle(deg=90)
		angle_delta = degrees90/corner_sides

		# Assemble *polygon_points* which are the points which are the outline of the
		# pocket:
		polygon_points = []
		# The first corner has a *start_angle* that points downward:
		start_angle = -degrees90

		# Now iterate through each *bend_point* output the the *corner_sides* + 1 points
		# that make up the pocket bend:
		for bend_point in bend_points:
		    angle = start_angle
		    for index in range(corner_sides + 1):
			adjust = P.polar(angle + rotate, radius)
			corner_point = bend_point + adjust
			#print("bend_point={0:i} adjust={1}".format(bend_point, adjust))
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

		# Convert *polygon_command_parts* into a *polygon_command* and append
		# it to *difference_lines*:
		polygon_command = "".join(polygon_command_parts)
		difference_lines.append(polygon_command)

		if trace_detail >= 2:
		    print("{0}Part.simple_pocket: STL_MODE done".format(indent))

	    if ezcad._cnc_mode:
		if trace_detail >= 2:
		    print("{0}Part.simple_pocket: CNC_MODE started".format(indent))

		# Grab the *top_surface_transform* from *part* that has been previously
		# set to orient the material properly for the CNC machine orgin:
		mount = part._current_mount
		cnc_transform = mount._cnc_transform_get()

		# Figure out *dz* which is the depth the pocket:
		cnc_corner1 = cnc_transform * corner1
		cnc_corner2 = cnc_transform * corner2
		dz = (cnc_corner2.z - cnc_corner1.z).absolute()

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
		if trace_detail >= 2:
		    print("{0}corner1={1:i} corner2={2:i}".format(indent, corner1, corner2))
		operation_order = Operation.ORDER_END_MILL_SIMPLE_POCKET
		operation_simple_pocket = Operation_Simple_Pocket(self, comment, 0,
		  end_mill_tool, operation_order, None, end_mill_feed_speed, end_mill_spindle_speed,
		  corner1, corner2, radius, end_mill_radius, Operation.POCKET_KIND_FLAT,
		  rotate, tracing = tracing + 1)
		self._operation_append(operation_simple_pocket)

		if trace_detail >= 2:
		    print("{0}Part.simple_pocket: CNC_MODE done".format(indent))

	# Perform any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.simple_pocket('{1}', '{2}', {3:i}, {4:i}, {5:i}, '{6}')".
	      format(indent, self._name, comment, corner1, corner2, radius, flags))

    def manufacture(self):
	""" Part: Override this method to do any manufacturing steps. """

	print("No manufacture method for '{0}'".format(self))

    def multi_mounts_find(self, multi_mounts_name, cell_dx=0, cell_dy=0):
	""" *Part*: Return the *Multi_Mounts* object that matches *multi_mouns_name* using
	    the *Part* object (i.e. *self*) as an anchor to start the search.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(multi_mounts_name, str) and not ' ' in multi_mounts_name
	assert isinstance(cell_dx, int) and cell_dx >= 0
	assert isinstance(cell_dy, int) and cell_dy >= 0

	# Create/find the *Multi_Mounts* object associated with *multi_mounts_name*:
	ezcad = part._ezcad_get()
	multi_mounts = ezcad._multi_mounts_find(multi_mounts_name, cell_dx, cell_dy)

	return multi_mounts


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

    def tooling_plate_drill(self, plate_mount_name, columns, rows, skips, tracing=-1000000):
	""" *Part*: Force the drilling of tooling plate mounting holes in the *Part* object
	    (i.e. *self*).  *plate_mount_name* specifies the name of the *Mount* object
	    that is created.  *columns* and *rows* specify a grid of holes to be drilled in the
	    *Part* object.  *skips* is a list of grid locations in the grid that will *not* be
	    drilled.  Thus, an empty  *skips* list causes all grid locations to be drilled.
	"""

	# To save on typing, "tooling_plate" is shortened down to "plate*:

	# Use *part* instead of *self:
	part = self

	# Verify argument types:
	assert isinstance(plate_mount_name, str) and not ' ' in plate_mount_name
	assert isinstance(columns, tuple) or isinstance(columns, list)
	assert isinstance(rows, tuple) or isinstance(rows, list)
	assert isinstance(skips, tuple) or isinstance(skips, list)
	for skip in skips:
	    assert isinstance(skip, tuple) and len(skip) == 2
	    for row_column in skip:
		assert isinstance(row_column, int) and row_column >= 0

	# Do some more argument verification:
	for row in rows:
	    assert isinstance(row, int)
	for column in columns:
	     assert isinstance(column, int)
	for skip in skips:
	    assert isinstance(skip, tuple) or isinstance(skip, list)
	    assert(len(skip) == 2) and isinstance(skip[0], int) and isinstance(skip[1], int)

	# Perform any requested *tracing*:
	trace_detail = -1
	deep_tracing = -1000000
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.tooling_plate_drill('{1}', '{2}', {3}, {4}, {5})".format(
	      indent, part._name, plate_mount_name, rows, columns, skips))
	    trace_detail = 2
	    if trace_detail >= 3:
		deep_tracing = tracing + 1

	# Only do stuff if in *ezcad* is in *cnc_mode*:
	ezcad = part._ezcad_get()
	if ezcad._cnc_mode:
	    # Grab the *plate* from *shop* and some associated values.
	    shop = part._shop_get()
	    plate = shop._tooling_plate_get()
	    plate_dx             = plate._dx_get()
	    plate_dy             = plate._dy_get()
	    plate_dz             = plate._dz_get()
	    plate_rows           = plate._rows_get()
	    plate_columns        = plate._columns_get()
	    plate_drill_diameter = plate._drill_diameter_get(part._material)
	    plate_hole_pitch     = plate._hole_pitch_get()
	    plate_spacer_dz      = plate._spacer_dz_get()
	    plate_spacer_width   = plate._spacer_width_get()

	    # Grab the *vice* from *shop* and some associated values:
	    vice = shop._vice_get()
	    vice_jaw_volume    = vice._jaw_volume_get()
	    vice_jaw_volume_dx = vice_jaw_volume.x
	    vice_jaw_volume_dy = vice_jaw_volume.y
	    vice_jaw_volume_dz = vice_jaw_volume.z
	    parallels          = vice._parallels_get()
	    
	    # There are two *Mount* objects in this code.  *vice_mount* refers to the previously
	    # specified mount of the material directly in the *vice*.  Shortly below, *tool_mount*
	    # is created to specify how the material is mounted on the tooling *plate*:
            
	    # We need to know the how the *part* bounding box lays out after the
	    # *top_transform* has been applied:
	    part_bsw = part.bsw
	    part_tne = part.tne
	    vice_mount = part._current_mount_get()
	    top_transform = vice_mount._top_surface_transform_get()
	    top_part_bsw = top_transform * part_bsw
	    top_part_tne = top_transform * part_tne
	    fixed_part_bsw, fixed_part_tne = \
	      top_part_bsw.minimum_maximum(top_part_tne)
	    fixed_part_volume = fixed_part_tne - fixed_part_bsw
	    part_dx = fixed_part_volume.x
	    part_dy = fixed_part_volume.y
	    part_dz = fixed_part_volume.z
	    if trace_detail >= 2:
		print("{0}part_bsw={1:i} part_tne={2:i}".
		  format(indent, part_bsw, part_tne))
		print("{0}top_part_bsw={1:i} top_part_tne.{2:i}".
		  format(indent, top_part_bsw, top_part_tne))
		print("{0}fixed_part_bsw={1:i} fixed_part_tne.{2:i}".
		  format(indent, fixed_part_bsw, fixed_part_tne))
		print("{0}part_dx={1:i} part_dy={2:i} part_dz={3:i}".
		  format(indent, part_dx, part_dy, part_dz))

	    # We will eventually need to know the height of the *part* with extra material:
	    extra_stop_bsw, extra_stop_tne = vice_mount._extra_stop_get()
	    top_extra_stop_bsw, top_extra_stop_tne = \
	      (top_transform * extra_stop_tne).minimum_maximum(top_transform * extra_stop_bsw)
	    extra_dz = top_extra_stop_tne.z - top_extra_stop_bsw.z
	    if trace_detail >= 2:
		print("{0}extra_dz={1:i}".format(indent, extra_dz))
	    zero = L()
	    assert extra_dz >= zero

	    # Figure out which *parallel_height* to use from *parallels* to use to mount
	    # the *tooling_plate*:
	    selected_parallel_height = parallels._select(part_dz, vice, tracing = deep_tracing)

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

	    # Sort *columns* and *rows*:
	    sorted_columns = sorted(tuple(columns))
	    sorted_rows = sorted(tuple(rows))

	    # Grab the minimum and maximum from *rows* and *columns*:
	    minimum_row = sorted_rows[0]
	    maximum_row = sorted_rows[-1]
	    minimum_column = sorted_columns[0]
	    maximum_column = sorted_columns[-1]

	    # Make sure that we do not try to use a hole that does not exist:
	    assert minimum_column >= 0, \
	      "Negative column {0} in {1} not allowed".format(minimum_column, columns)
	    assert maximum_column < plate_columns, \
	      "Column {0} in {1} too big for tooling plate column maximum {2}". \
              format(maximum_column, columns, plate_columns - 1)
	    assert minimum_row >= 0, \
	      "Negative row {0} in {1} not allowed".format(minimum_row, rows)
	    assert maximum_row < plate_rows, \
	      "Row {0} in {1} too big for tooling plate row maximum". \
              format(maximum_row, rows, plate_rows - 1)

	    # Compute the *rows_spanned* and *columns_spanned*:
	    rows_spanned = maximum_row - minimum_row
	    columns_spanned = maximum_column - minimum_column

	    # Figure out where upper left tooling hole is located:
	    #FIXME: Is the `2 * minimum_column` correct?!!!
	    lower_left_x = -((2 * minimum_column + columns_spanned) * plate_hole_pitch) / 2
	    lower_left_y = -((2 * minimum_row  + rows_spanned) * plate_hole_pitch) / 2
	    zero = L()

	    # Drill the *tooling_plate* holes at the intersections of *rows* and *columns*,
	    # but ommitting any that are in *skips*:
	    plate_holes = []
	    reverse_top_transform = top_transform.reverse()
	    if trace_detail >= 1:
		print("{0}top_transform={1:s}".
		  format(indent, top_transform))
		print("{0}reverse_top_surface_transform={1:s}".
		  format(indent, reverse_top_transform))
	    for row in rows:
		#FIXME: We really need to understand maximum tapping depth:
		maximum_drill_depth = L(inch=1.000)
		for column in columns:
		    column_row = (column, row)
		    if not column_row in skips_table:
			# Figure out the *start* and *stop* points for the drill as if the
			# part is centered immediately under the origin (i.e. after
			# *top_sufrace_transform* is applied to the entire *part*.):
			x = lower_left_x + plate_hole_pitch * column
			y = lower_left_y + plate_hole_pitch * row
			start = P(x, y, zero)
			drill_depth = part_dz.minimum(maximum_drill_depth)
			stop =  P(x, y, -drill_depth)

			plate_holes.append(column_row)

			# Drill the hole using *reversed_start* and *reversed_stop*:
			reversed_start = reverse_top_transform * start
			reversed_stop = reverse_top_transform * stop
			part.hole("Tooling plate hole ({0}, {1})".format(row, column),
			  plate_drill_diameter, reversed_start, reversed_stop, "t",
			  tracing = deep_tracing)
			if trace_detail >= 1:
			    print(("{0}column={1} row={2} " +
			      "start={3:i} stop={4:i} rstart={5:i} rstop={6:i}").format(
			      indent, column, row, start, stop, reversed_start, reversed_stop))

	    # We are done with drilling the holes.  Now we move on to figuring out where
	    # everything for the *Part.tooling_plate_mount()* routine.  All the information
	    # computed here is stuffed into the *tooling_plate* object.

	    # Identify the location of (0, 0) *tooling_plate* hole:
	    plate_columns_dx = (plate_columns - 1) * plate_hole_pitch
	    plate_rows_dy    = (plate_rows - 1) * plate_hole_pitch
	    plate_edge_dx    = (plate_dx - plate_columns_dx) / 2
	    plate_edge_dy    = (plate_dy - plate_rows_dy) / 2

	    # Determine which parallel to use to mount the tooling plate such that
	    # the top surface of the tooling plate is near the top surface the vice:
	    selected_parallel_height = parallels._select(plate_dz, vice, tracing = deep_tracing)

	    # The *tooling_plate* is placed so that its upper left corner touchs the upper left
	    # jaw vice corner.  The part is placed so that it is centered in between the
	    # minimum/maximum number of spanned rows/colums        
	    x = plate_edge_dx + minimum_column * plate_hole_pitch + \
	      (columns_spanned * plate_hole_pitch) / 2
	    y = -(plate_edge_dy + minimum_row * plate_hole_pitch +
	      (rows_spanned * plate_hole_pitch) / 2)
	    z = selected_parallel_height + plate_dz + plate_spacer_dz + part_dz
	    plate_translate_point = P(x, y, z)
	    top_surface_safe_z = z + (extra_dz - part_dz).maximum(zero)
	    xy_rapid_safe_z = top_surface_safe_z + L(inch=0.500)

	    # Create the *tooling_plate_mount*:
	    if trace_detail >= 2:
		print("{0}!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!".format(indent))
		print("{0}selected_parallel_height={1:i}".format(indent, selected_parallel_height))
		print("{0}plate_dz={1:i}".format(indent, plate_dz))
		print("{0}plate_spacer_dz={1:i}".format(indent, plate_spacer_dz))
		print("{0}part_dz={1:i}".format(indent, part_dz))
		print("{0}top_surface_safe_z={1:i}".format(indent, top_surface_safe_z))
		print("{0}xy_rapid_safe_z={1:i}".format(indent, xy_rapid_safe_z))
		print("{0}plate_mount.top_transform={1:v}".format(indent, top_transform))
		print("{0}plate_mount.plate_translate_point={1:i}".
		  format(indent, plate_translate_point))
	    plate_mount = Mount(plate_mount_name, part, top_transform,
	      plate_translate_point, top_surface_safe_z, xy_rapid_safe_z, True,
	      selected_parallel_height, tracing = deep_tracing)

	    # Deal with extra material:
	    vice_mount = part._current_mount_get()
	    #vice_mount_cnc_transform = vice_mount._cnc_transform_get()
	    #reverse_vice_cnc_transform = vice_mount_cnc_transform.reverse()
	    #plate_mount_cnc_transform = plate_mount._cnc_transform_get()

	    #vice_extra_start_bsw, vice_extra_start_tne = vice_mount._extra_start_get()
	    #part_extra_start_bsw = reverse_vice_cnc_transform * vice_extra_start_bsw
	    #part_extra_start_tne = reverse_vice_cnc_transform * vice_extra_start_tne
	    #transformed_extra_start_bsw = plate_mount_cnc_transform * part_extra_start_bsw
	    #transformed_extra_start_tne = plate_mount_cnc_transform * part_extra_start_tne
	    #plate_extra_start_bsw, plate_extra_start_tne = \
	    #  transformed_extra_start_bsw.minimum_maximum(transformed_extra_start_tne)
	    #plate_mount._extra_start_set(plate_extra_start_bsw, plate_extra_start_tne)
	    #plate_mount._extra_stop_set( plate_extra_start_bsw, plate_extra_start_tne)

	    extra_start_bsw, extra_start_tne = vice_mount._extra_start_get()
	    plate_mount._extra_start_set(extra_start_bsw, extra_start_tne)
	    plate_mount._extra_stop_set( extra_start_bsw, extra_start_tne)

	    # Remember the *plate_holes* in both *vice_mount* and *plate_mount*.
	    # *vice_mount* is needed for multi-mounting:
	    vice_mount._tooling_plate_holes_set( plate_holes, tracing = tracing + 1)
	    plate_mount._tooling_plate_holes_set(plate_holes, tracing = tracing + 1)

	    # Compute the *dowel_point*:
	    #extra_dx = extra_start_tne.x - extra_start_bsw.x
	    extra_dx = zero

	    plate_cnc_transform = plate_mount._cnc_transform_get()
	    plate_extra_start_bsw = plate_cnc_transform * extra_start_bsw
	    plate_extra_start_tne = plate_cnc_transform * extra_start_tne
	    dowel_point = P(plate_extra_start_bsw.x, \
	      (plate_extra_start_bsw.y + plate_extra_start_tne.y)/2, plate_extra_start_bsw.z)

	    # Specify which *spacers* to visualize:
	    spacers = []
	    if len(rows) >= len(columns):
		# Orient the spacers horizontally:
		for row in rows:
		    spacers.append( (minimum_column, row, maximum_column, row) )
	    else:
		# Orient the spacers vertically:
		for column in columns:
		    spacers.append( (column, minimum_row, column, maximum_row) )
	    plate_mount._spacers_set(spacers)
	    if trace_detail >= 1:
		print("{0}spacers={1}".format(indent, spacers))

	    # Now remember where to mount the part with the tooling plate:
	    plate._mount_reset()
	    plate._dowel_point_set(dowel_point)
	    plate._mount_set(plate_mount)
	    plate._parallels_height_set(selected_parallel_height)
	    if trace_detail >= 1:
		print("{0}parallels_height_set({1:i})".format(indent, selected_parallel_height))
	    plate._spacers_set(spacers)
	    plate._xy_rapid_safe_z_set(xy_rapid_safe_z)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Part.tooling_plate_drill('{1}', '{2}', {3}, {4}, {5})".format(
	      indent, part._name, plate_mount_name, rows, columns, skips))

    def tooling_plate_mount(self, comment, tracing=-100000):
	""" *Part*: Cause the mounting plate that holds the *Part* object (i.e. *self*)
	    to be mounted in the vice using a dowel pin operation. """

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(comment, str)

	# Perform any requested *tracing*;
	trace_detail = -1
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.tooling_plate_mount('{1}', '{2}')".format(indent, self._name, comment))
	    trace_detail = 2

	# Only do this code when CNC generation is on:
	ezcad = part._ezcad_get()
	if ezcad._cnc_mode:
	    shop = part._shop_get()
	    tooling_plate = shop._tooling_plate_get()
	    dowel_point = tooling_plate._dowel_point_get()
	    dy = tooling_plate._dy_get()
	    mount = tooling_plate._mount_get()
	    if trace_detail >= 1:
		print("{0}mount='{1}'".format(indent, mount._name_get()))
	    parallels_height = tooling_plate._parallels_height_get()
	    spacers = tooling_plate._spacers_get()

	    # Register the new *mount* in *part*:
	    part._mount_register(mount, tracing = tracing + 1)

	    # Create *operation_mount* and add it to *part* operations list:
	    vice = shop._vice_get()
	    operation_mount = Operation_Mount(part, comment,
	      vice, dy, parallels_height, tooling_plate, spacers, True, tracing = tracing + 1)
	    part._operation_append(operation_mount)

	    # Perform a dowel pin operation:
	    zero = L()
	    plunge_point = dowel_point - P(L(inch=1.000), zero, zero)
	    preferred_diameter = L(inch="3/8")
	    part.dowel_pin("tooling plate dowel pin",
	      dowel_point, plunge_point, preferred_diameter, tracing = tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Part.tooling_plate_mount('{1}', '{2}')".format(indent, self._name, comment))

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
	    {value} is returned.
	 """

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
	    error messages and any generated G-code.
	"""

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

    def vice_mount(self, name, top_surface, jaw_surface, flags,
      extra_dx=L(), extra_dy=L(), extra_top_dz=L(), extra_bottom_dz=L(), tracing=-100000):
	""" *Part*: Cause the *Part* object (i.e. *self*) to be mounted in a vice with
	    *top_surface* facing upwards and *jaw_surface* mounted towards the rear vice jaw.
	    *top_surface* and *jaw_surface* must be one of 't'" (top), 'b' (bottom), 'n' (north),
	    's' (south), 'e' (east) or 'w' (west).  *comment* is the attached to any
	    generated G-code.  The *Part* dimensions are specified by its bounding box.
	    This routine has a required *flags* string argument can be left empty.
	    The *flags* argument is to support a dowel pin operation.  (See *Part.dowel_pin()*
	    for more information about a dowel pin operation.)  When *flags* contains 'l',
	    "l" it generates a left dowel pin operation, and when it contains 'r' it generates
	    a right dowel pin operation.  The optional *extra_dx*, *extra_dy*, *extra_top_dz*
	    and *extra_bottom_dz* arguments specify how much padding is added to the physical
	    object to prior to being latched into the vice.  (Note that *extra_dx*, *extra_dy*,
	    *extra_top_dz*, extra_bottom_dz* do *NOT* modify the *Part* object bounding_box
	    dimensions in any way.)  The optional *extra_dx* and *extra_dy* specify that the
	    physical part placed into the vice has been padded in the dx and dy vice coordinates.
	    (Half goes on each side.)  Likewise, the optional *extra_top_dz* and *extra_bottom_dz*
	    arguments indicate that the top and bottom of the physical part placed into the vice
	    are padded with material on the top and/or bottom (usually to allow for top facing
	    operations.)
	"""

	# Use *part* instead of self:
	part = self

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(top_surface, str) and len(top_surface) == 1 and top_surface in "tbnsew"
	assert isinstance(jaw_surface, str) and len(jaw_surface) == 1 and jaw_surface in "tbnsew"
	assert isinstance(flags, str)
	for flag in flags:
	    assert flag in "lr", "Flag '{0}' is not allowed for flags argument".format(flag)
	assert isinstance(extra_dx, L) and extra_dx >= zero
	assert isinstance(extra_dy, L) and extra_dy >= zero
	assert isinstance(extra_top_dz, L) # and extra_top_dz >= zero
	assert isinstance(extra_bottom_dz, L) # and extra_bottom_dz >= zero

	# Preform any requested *tracing*:
	detail_tracing = -1000000
	trace_detail = -1
	if tracing < 0 and part._tracing >= 0:
 	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Part.vice_mount('{1}', '{2}', '{3}', '{4}', '{5}, " +
	      "{6:i}, {7:i}, {8:i}, {9:i})").format(indent,
	      part._name, name, top_surface, jaw_surface, flags,
	      extra_dx, extra_dy, extra_top_dz, extra_bottom_dz))
	    trace_detail = 3
	    detail_tracing = tracing + 1

	# We do not need to do anything until we are in CNC mode:
	ezcad = part._ezcad_get()
	if ezcad._cnc_mode:
	    # We need to compute *top_surface_transform* for CNC mode:

	    # The *x_axis*, *y_axis*, and *z_axis* are needed for rotations:
	    one = L(mm=1.0)
	    x_axis = P(one, zero, zero).normalize()
	    y_axis = P(zero, one, zero).normalize()
	    z_axis = P(zero, zero, one).normalize()
	    if trace_detail >= 4:
		print("{0}x_axis={1:m} y_axis={2:m} z_axis={3:m}".
		  format(indent, x_axis, y_axis, z_axis))

	    # These are the two rotation constants used to reorient the bounding box:
	    degrees90 = Angle(deg=90.0)
	    degrees180 = Angle(deg=180.0)

	    # Grab the 6 possible center surface points (and *center_point*) of the *bounding_box*
	    # for use below:
	    c = part.c
	    t = part.t
	    b = part.b
	    n = part.n
	    s = part.s
	    e = part.e
	    w = part.w

	    # *top_axis* and *top_rotate* are used to rotate the desired bounding box surface
	    # facing up in the vice.  Leaving *top_rotate_angle* as *None* means no top
	    # rotation is needed.  *top_point* matches *top_surface*:
	    top_point = None
	    top_axis = None
	    top_rotate_angle = None

	    # *jaw_axis* and *jaw_rotate* are used to rotate the desired bounding box surface
	    # facing towards the rear vice jaw.  Leaving *jaw_rotate_angle * as *None* means no
	    # jaw rotaton needed.  *jaw_point* matches *jaw_surface*:
	    jaw_point = None
	    jaw_axis = z_axis
	    jaw_rotate_angle = None

	    # *left_dowel_point* specifies the *bounding_box* center surface point to the "left"
	    # of the bounding box center when mounted in the vice.  Likewise, for
	    # *right_dowel_point*.  One of these two points is used for a dowel pin operation.
	    left_dowel_point = None
	    right_dowel_point = None

	    # Do a 24 (=6 surfaces x 4 jaw orientations) dispatch on *top_surface* and *jaw_suface*:
	    if 't' in top_surface:
		# Top surface of *bounding_box* facing up from vice:
		top_point = t
		if 'n' in jaw_surface:
		    # No rotation needed
		    jaw_point = n
		    left_dowel_point = w
		    right_dowel_point = e
		elif 's' in jaw_surface:
		    jaw_point = s
		    jaw_rotate_angle = degrees180
		    left_dowel_point = e
		    right_dowel_point = w
		elif 'e' in jaw_surface:
		    jaw_point = e
		    jaw_rotate_angle = degrees90
		    left_dowel_point = n
		    right_dowel_point = s
		elif 'w' in jaw_surface:
		    jaw_point = w
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = s
		    right_dowel_point = n
		else:
		    assert False, "jaw_surface must be one of 'n', 's', 'e', or 'w'"
	    elif 'b' in top_surface:
		# Bottom surface of *bounding_box* facing up from vice:
		top_point = b
		top_axis = x_axis
		top_rotate_angle = degrees180
		if 'n' in jaw_surface:
		    jaw_point = n
		    jaw_rotate_angle = degrees180
		    left_dowel_point = e
		    right_dowel_point = w
		elif 's' in jaw_surface:
		    jaw_point = s
		    jaw_rotate_angle = None
		    left_dowel_point = w
		    right_dowel_point = e
		elif 'e' in jaw_surface:
		    jaw_point = e
		    jaw_rotate_angle = degrees90
		    left_dowel_point = s
		    right_dowel_point = n
		elif 'w' in jaw_surface:
		    jaw_point = w
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = n
		    right_dowel_point = s
		else:
		    assert False, "jaw_surface must be one of 'n', 's', 'e', or 'w'"
	    elif 'n' in top_surface:
		# North surface of *bounding_box* facing up from vice:
		top_point = n
		top_axis = x_axis
		top_rotate_angle = degrees90
		if 't' in jaw_surface:
		    jaw_point = t
		    jaw_rotate_angle = -degrees180
		    left_dowel_point = e
		    right_dowel_point = w
		elif 'b' in jaw_surface:
		    jaw_point = b
		    jaw_rotate_angle = None
		    left_dowel_point = w
		    right_dowel_point = e
		elif 'e' in jaw_surface:
		    jaw_point = e
		    jaw_rotate_angle = degrees90
		    left_dowel_point = b
		    right_dowel_point = t
		elif 'w' in jaw_surface:
		    jaw_point = w
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = t
		    right_dowel_point = b
		else:
		    assert False, "jaw_surface must be one of 't', 'b', 'e', or 'w'"
	    elif 's' in top_surface:
		# South surface of *bounding_box* facing up from vice:
		top_point = s
		top_axis = x_axis
		top_rotate_angle = -degrees90
		if 't' in jaw_surface:
		    jaw_point = t
		    jaw_rotate_angle = None
		    left_dowel_point = w
		    right_dowel_point = e
		elif 'b' in jaw_surface:
		    jaw_point = b
		    jaw_rotate_angle = degrees180
		    left_dowel_point = e
		    right_dowel_point = w
		elif 'e' in jaw_surface:
		    jaw_point = e
		    jaw_rotate_angle = degrees90
		    left_dowel_point = t
		    right_dowel_point = b
		elif 'w' in jaw_surface:
		    jaw_point = w
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = b
		    right_dowel_point = t
		else:
		    assert False, "jaw_surface must be one of 't', 'b', 'e', or 'w'"
	    elif 'e' in top_surface:
		# East surface of *bounding_box* facing up from vice:
		top_point = e
		top_axis = y_axis
		top_rotate_angle = -degrees90
		if 't' in jaw_surface:
		    jaw_point = t
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = s
		    right_dowel_point = n
		elif 'b' in jaw_surface:
		    jaw_point = b
		    jaw_rotate_angle = degrees90
		    left_dowel_point = n
		    right_dowel_point = s
		elif 'n' in jaw_surface:
		    jaw_point = n
		    jaw_rotate_angle = None
		    left_dowel_point = t
		    right_dowel_point = b
		elif 's' in jaw_surface:
		    jaw_point = s
		    jaw_rotate_angle = -degrees180
		    left_dowel_point = b
		    right_dowel_point = t
		else:
		    assert False, "jaw_surface must be one of 't', 'b', 'n', or 's'"
	    elif 'w' in top_surface:
		# West surface of *bounding_box* facing up from vice:
		top_point = w
		top_axis = y_axis
		top_rotate_angle = degrees90
		if 't' in jaw_surface:
		    jaw_point = t
		    jaw_rotate_angle = degrees90
		    left_dowel_point = n
		    right_dowel_point = s
		elif 'b' in jaw_surface:
		    jaw_point = b
		    jaw_rotate_angle = -degrees90
		    left_dowel_point = s
		    right_dowel_point = n
		elif 'n' in jaw_surface:
		    jaw_point = n
		    jaw_rotate_angle = None
		    left_dowel_point = b
		    right_dowel_point = t
		elif 's' in jaw_surface:
		    jaw_point = s
		    jaw_rotate_angle = degrees180
		    left_dowel_point = t
		    right_dowel_point = b
		else:
		    assert False, "jaw_surface must be one of 't', 'b', 'n', or 's'"
	    else:
		assert False, "top_surface must be one of 't', 'b', 'n', 's', 'e', or 'w'"

	    # Perform any requested *trace_detail*:
	    if trace_detail >= 3:
		if isinstance(top_point, P):
		    print("{0}top_point={1:i}".format(indent, top_point))
		else:
		    print("{0}top_point=None".format(indent))
		if isinstance(top_axis, P):
		    print("{0}top_axis={1:m}".format(indent, top_axis))
		else:
		    print("{0}top_axis=None".format(indent))
		print("{0}top_rotate_angle={1}".format(indent, top_rotate_angle))

		if isinstance(jaw_point, P):
		    print("{0}jaw_point={1:i}".format(indent, jaw_point))
		else:
		    print("{0}jaw_point=None".format(indent))
		if isinstance(jaw_axis, P):
		    print("{0}jaw_axis={1:m}".format(indent, jaw_axis))
		else:
		    print("{0}jaw_axis=None".format(indent))
		print("{0}jaw_rotate_angle={1}".format(indent, jaw_rotate_angle))
	    if trace_detail >= 2:
		print("{0}center={1:i} top={2:i} jaw={3:i}".format(indent, c, top_point, jaw_point))
		print("{0}left_dowel={1:i} right_dowel={2:i}".
		  format(indent, left_dowel_point, right_dowel_point))

	    # Now we can compute the *top_surface_transform*, which rotates the *part* bounding box
	    # so that it is oriented correctly for the vice with the top surface located at
	    # *cnc_top_surface_z*:
	    top_surface_transform = Transform()
	    if trace_detail >= 3:
		print("{0}top_surface_transform 0 = {1:v}".format(indent, top_surface_transform))
	    top_surface_transform = top_surface_transform.translate(
	      "move to vice origin", -top_point, tracing = tracing + 1)
	    if trace_detail >= 3:
		print("{0}top_surface_transform 1 = {1:v}".format(indent, top_surface_transform))
	    if isinstance(top_rotate_angle, Angle):
		top_surface_transform = top_surface_transform.rotate(
		"rotate top surface to point up in vice", top_axis, top_rotate_angle,
		tracing = tracing + 1)
	    if trace_detail >= 3:
		print("{0}top_surface_transform 2 = {1:v}".format(indent, top_surface_transform))
	    if isinstance(jaw_rotate_angle, Angle):
		top_surface_transform = top_surface_transform.rotate(
		"rotate jaw surface to point north in vice", jaw_axis, jaw_rotate_angle,
		tracing = tracing + 1)
	    if trace_detail >= 2:
		print("{0}top_surface_transform 3={1:v}".format(indent, top_surface_transform))

	    # Register the block with the *bom* (Bill Of Materials):
	    bom = part._bom
	    bom._block_register(top_surface_transform)

	    # Grab some values from the *vice*:
	    shop = part._shop_get()
	    vice = shop._vice_get()
	    parallels = vice._parallels_get()
	    vice_jaw_volume = vice._jaw_volume_get()
	    vice_jaw_dx = vice_jaw_volume.x
	    vice_jaw_dy = vice_jaw_volume.y
	    vice_jaw_dz = vice_jaw_volume.z
	    if trace_detail >= 1:
		print("{0}vice_jaw_volume={1:i}".format(indent, vice_jaw_volume))

	    # Now compute the *half_part_dx*, *half_part_dy*, and *half_part_dz*:
	    half_part_dx = c.distance(left_dowel_point)
	    half_part_dy = c.distance(jaw_point)
	    half_part_dz = c.distance(top_point)

	    # Compute *part_dx*, *part_dy*, and *part_dz*:
	    part_dx = 2 * half_part_dx
	    part_dy = 2 * half_part_dy
	    part_dz = 2 * half_part_dz
	    if trace_detail >= 2:
	    	print("{0}part dimensions in vice: dx={1:i} dy={2:i} dz={3:i}".
	    	  format(indent, part_dx, part_dy, part_dz))

	    # Grab the associated *parallels_heights* from *parallels* and get
	    # *smallest_parallel_height*:
	    parallels_heights = sorted(parallels._heights_get())
	    smallest_parallel_height = parallels_heights[0]

	    # Now figure out *top_surface_z*, which is the distance from the vice origin to the
	    # center of the top surface of the *part*.  We want to mount such that the top surface
	    # as near to the top jaw edge as possible.  However, sometimes the parallels will
	    # not let us get that low:
	    selected_parallel_height = smallest_parallel_height
	    epsilon = L(inch=.000001)
	    for index, parallel_height in enumerate(parallels_heights):
		# Use *epsilon* to deal with rounding errors:
		trial_top_surface_z = \
		  parallel_height + extra_bottom_dz + part_dz + extra_top_dz - epsilon
		if parallel_height > selected_parallel_height and \
		  trial_top_surface_z <= vice_jaw_dz:
		    selected_parallel_height = parallel_height
		    if trace_detail >= 3:
			print("{0}Parallel[{1}]: parellel_height={2:i} trial_top_surface_z={3:i}".
			  format(indent, index, parallel_height, trial_top_surface_z))
	    if trace_detail >= 1:
		print("{0}parallels_heights={1}".format(indent,
		  [parallel_height.inches() for parallel_height in parallels_heights]))
		print("{0}selected_parallel_height={1:i}".format(indent, selected_parallel_height))

	    # Now we can compute *cnc_top_surface_z* and *cnc_xy_rapid_safe_z*:
	    cnc_top_surface_z = selected_parallel_height + extra_bottom_dz + part_dz
	    cnc_top_surface_safe_z = cnc_top_surface_z + extra_top_dz
	    cnc_xy_rapid_safe_z = cnc_top_surface_safe_z + L(inch=0.5)
	    if trace_detail >= 1:
		print(("{0}cnc_top_surface_z={1:i} " +
		  "cnc_top_surface_safe_z={2:i} cnc_xy_rapid_safe_z={3:i}").
		  format(indent, cnc_top_surface_z, cnc_top_surface_safe_z, cnc_xy_rapid_safe_z))

	    # Process *flags* to figure out if a dowel pin operation is needed:
	    dowel_pin_requested = False
	    is_left_dowel_pin = False
	    is_right_dowel_pin = False
	    if 'l' in flags:
		dowel_pin_requested = True
		is_left_dowel_pin = True
	    if 'r' in flags:
		dowel_pin_requested = True
		is_right_dowel_pin = True
	    assert not (is_left_dowel_pin and is_right_dowel_pin), \
	      "Specifying both 'l' and 'r' in flags argument will not work"

	    # Now we can figure out how we are going to position the *part* in the vice in Y axis.
	    vice_y = -(half_part_dy + extra_dy/2)

	    # Now we figure out where to put the *part* in the vice X axis.  By default, *vice_x*
	    # is set to the dead center of the vice in X:
	    vice_x = vice_jaw_dx/2
	    part_centered_in_vice = True
	    if dowel_pin_requested and part_dx + extra_dx < vice_jaw_dx:
		# We can position the part to the left or the right:
		if is_left_dowel_pin:
		    vice_x = half_part_dx + extra_dx/2
		    part_centered_in_vice = False
		elif is_right_dowel_pin:
		    vice_x = vice_jaw_dx - half_part_dx - extra_dx/2
		    part_centered_in_vice = False
		else:
		    assert False, "Internal error: Neither left or right dowel pin is specified"
	    if trace_detail >= 3:
		print("{0}is_left_dowel_pin={1} is_right_dowel_pin={2}".
		  format(indent, is_left_dowel_pin, is_right_dowel_pin))
		print("{0}part_centered_in_vice={1} vice_x={2:i} vice_y={3:i}".
		  format(indent, part_centered_in_vice, vice_x, vice_y))

	    # Now we can compute *mount_translate_point* which tranlates the *part* from the
	    # machine origin to the correct location in the *vice*:
	    mount_translate_point = P(vice_x, vice_y, cnc_top_surface_z)
	    if trace_detail >= 1:
		print("{0}top_surface_transform={1:s}".format(indent, top_surface_transform))
		print("{0}mount_translate_point={1:i}".format(indent, mount_translate_point))

	    # Now we can issue the mount command for the *part*:
	    cnc_top_surface_safe_z = cnc_top_surface_z + extra_top_dz
	    cnc_xy_rapid_safe_z    = cnc_top_surface_z + extra_top_dz + L(inch=0.500)
	    
	    # Create *mount* and stuff into *part*:
	    mount = Mount(name, part, top_surface_transform, mount_translate_point,
	      cnc_top_surface_safe_z, cnc_xy_rapid_safe_z, False, selected_parallel_height,
	      tracing = tracing + 1)
	    part._mount_register(mount, tracing = tracing + 1)

	    # We can only deal with the extra material if the part is oriented correctly.
	    # Thus, we transform the part bounding box using *top_transform*, add in the
	    # extra material, reverse the transform back, and stuff the results back into *mount*.

	    # Grab the *part* bounding box corners after being mapped using *top_transform*:
	    top_transform = mount._top_surface_transform_get()
	    top_bounding_bsw, top_bounding_tne = \
	      part._bounding_box_get(top_transform, tracing = tracing + 1)

	    # Now add in the extra material:
	    top_extra_start_bsw = top_bounding_bsw - P(extra_dx/2, extra_dy/2, extra_bottom_dz)
	    top_extra_start_tne = top_bounding_tne + P(extra_dx/2, extra_dy/2, extra_top_dz)

	    # Now reverse the extra material back into the original *part* locations:
	    top_reverse_transform = top_transform.reverse()
	    extra_start_bsw = top_reverse_transform * top_extra_start_bsw
	    extra_start_tne = top_reverse_transform * top_extra_start_tne

	    # Now save the values back into *mount* as the initial start and stop values
	    # for extra material:
	    mount._extra_start_set(extra_start_bsw, extra_start_tne)
	    mount._extra_stop_set( extra_start_bsw, extra_start_tne)
	    if trace_detail >= 2:
		print("{0}top_transform={1:s}".format(indent, top_transform))
		print("{0}top_transformed bounding_bsw={1:i} bounding_tne={2:i}".
		  format(indent, top_bounding_bsw, top_bounding_tne))
		print("{0}top_reverse_transformed extra_start_bsw={1:i} extra_start_tne={2:i}".
		  format(indent, extra_start_bsw, extra_start_tne))

	    # Create *operation_mount* and add it to *part* operations list:
	    vice = shop._vice_get()
	    jaws_spread = part_dy + extra_dy
	    tooling_plate = None
	    spacers = []
	    operation_mount = Operation_Mount(part, name,
	      vice, jaws_spread, selected_parallel_height, tooling_plate, spacers, False)
	    part._operation_append(operation_mount)

	    # Generate any needed dowel pin operation:
	    if dowel_pin_requested:
		dowel_x = None
		plunge_x = None
		dowel_pin_comment = "no dowel pin"
		plunge_x_offset = L(inch=0.750)
		if is_left_dowel_pin:
		    dowel_pin_comment = "Left dowel pin"
		    if part_centered_in_vice:
			dowel_x = vice_jaw_dx/2 - (part_dx + extra_dx)/2
		    else:
			dowel_x = zero
		    plunge_x = dowel_x - plunge_x_offset
		elif is_right_dowel_pin:
		    dowel_pin_comment = "Right dowel pin"
		    if part_centered_in_vice:
			dowel_x = vice_jaw_dx/2 + (part_dx + extra_dx)/2
		    else:
			dowel_x = vice_jaw_dx
		    plunge_x = dowel_x + plunge_x_offset
		else:
		    assert False

		# Create *cnc_plunge_point* and *cnc_dowel_point*:
		cnc_plunge_point = P(plunge_x, vice_y, zero)
		cnc_dowel_point = P(dowel_x, vice_y, zero)
		if trace_detail >= 2:
		    print("{0}cnc_dowel_point={1:i} cnc_plunge_point={2:i}".format(
		      indent, cnc_dowel_point, cnc_plunge_point))

		# Perform the dowel pin operation:
		preferred_diameter = L(inch="3/8")

		part.dowel_pin(dowel_pin_comment,
		  cnc_dowel_point, cnc_plunge_point, preferred_diameter, detail_tracing)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print(("{0}<=Part.vice_mount('{1}', '{2}', '{3}', '{4}', '{5}, " +
	      "{6:i}, {7:i}, {8:i}, {9:i})").format(indent,
	      part._name, name, top_surface, jaw_surface, flags,
	      extra_dx, extra_dy, extra_top_dz, extra_bottom_dz))

    def _vice_mount_helper(self, extra_bsw, extra_tne, is_tooling_plate_mount, tracing=-1000000):
	""" *Part*: Return a mount translate point and parallel height for the material bounding
	    box specified by the corners *extra_bsw* and *extra_tne*.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(extra_bsw, P)
	assert isinstance(extra_tne, P)
	assert isinstance(is_tooling_plate_mount, bool)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>_vice_mount_helper({1:i}, {2:i})".format(indent, extra_bsw, extra_tne))
	    trace_detail = 2

	# Now compute volume dimensions of the extra material bounding box:
	extra_dx = extra_tne.x - extra_bsw.x
	extra_dy = extra_tne.y - extra_bsw.y
	extra_dz = extra_tne.z - extra_bsw.z

	# The part is assumed to have its top surface at the machine origin.  Thus,
        # *extra_top_dz* the same as *extra_tne.z*:
	extra_top_dz = extra_tne.z

	# If there is any *extra_top_dz* we will probably be doing a facing operation.
        # For this situataion we want the part top surface to be a bit above vice surface.
	# This is done by setting *extra_facing_dz* to the desired amount above the vice:
	zero = L()
	extra_facing_dz = zero
	if extra_top_dz > zero:
	    # For now, we hard code the extra above the vice jaw to be .025 inches:
	    extra_facing_dz = L(inch=0.025)

	# Grab some values from the *vice*:
	shop = part._shop_get()
	vice = shop._vice_get()
	parallels = vice._parallels_get()
	vice_jaw_volume = vice._jaw_volume_get()
	vice_jaw_dx = vice_jaw_volume.x
	vice_jaw_dy = vice_jaw_volume.y
	vice_jaw_dz = vice_jaw_volume.z
	if trace_detail >= 1:
	    print("{0}vice_jaw_volume={1:i}".format(indent, vice_jaw_volume))

	# Grab the associated *parallels_heights* from *parallels* and get
	# *smallest_parallel_height* and *largest_parallel_height*:
	parallels_heights = sorted(parallels._heights_get())
	smallest_parallel_height = parallels_heights[0]
	selected_parallel_height = smallest_parallel_height
	if is_tooling_plate_mount:
	    # Deal with tooling plate:
	    tooling_plate = shop._tooling_plate_get()
	    tooling_plate_dz = tooling_plate._dz_get()
	    for parallel_height in parallels_heights:
		if parallel_height + tooling_plate_dz <= vice_jaw_dz:
		    selected_parallel_height = parallel_height
	else:
	    # Now we want to select a parallel height that places the top surface of the
	    # *part* at or just above the vice jaws.  Sometimes, the *part* is pretty
	    # thick, so we always start with the smallest parallel height:
	    for parallel_height in parallels_heights:
		if parallel_height > selected_parallel_height and \
		  parallel_height + extra_dz <= vice_jaw_dz + extra_facing_dz:
		    selected_parallel_height = parallel_height
	    if trace_detail >= 1:
		print("{0}selected_parallel_height={1:i}".format(indent, selected_parallel_height))

	# Now we can compute *cnc_top_surface_z* and *cnc_xy_rapid_safe_z*:
	#cnc_top_surface_z = selected_parallel_height + extra_dz
	#cnc_top_surface_safe_z = cnc_top_surface_z + extra_top_dz
	#cnc_xy_rapid_safe_z = cnc_top_surface_safe_z + L(inch=0.5)
	#if trace_detail >= 1:
	#    print(("{0}cnc_top_surface_z={1:i} " +
	#      "cnc_top_surface_safe_z={2:i} cnc_xy_rapid_safe_z={2:i}").
	#      format(indent, cnc_top_surface_z, cnc_top_surface_safe_z, cnc_xy_rapid_safe_z))

	# Now we can figure out how we are going to position the *part* in the vice in Y axis.
	vice_x = extra_dx/2
	vice_y = -extra_dy/2
	vice_z = selected_parallel_height + extra_dz - extra_top_dz + extra_facing_dz

	# Now we can compute *mount_translate_point* which translates the *part* from the
	# machine origin to the correct location in the *vice*:
	mount_translate_point = P(vice_x, vice_y, vice_z)
	if trace_detail >= 1:
	    print("{0}mount_translate_point={1:i}".format(indent, mount_translate_point))

	# Wrap up any requested *tracing* and return both *mount_translate_point* and
        # *selected_parallel_height*:
	if tracing >= 0:
	    print("{0}<=_vice_mount_helper({1:i}, {2:i})=>{3:i}, {4:i}".
	      format(indent, extra_bsw, extra_tne, mount_translate_point, selected_parallel_height))
	return mount_translate_point, selected_parallel_height

    def multi_mount(self, multi_mounts, plate, rotate, sort_order, tracing=-1000000):
	""" *Part*: Set up a multi mounted CNC tool for the *Part* object (i.e. *self*) that
	    consists of *multi_mounts*.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(multi_mounts, Multi_Mounts)
	assert isinstance(plate, Plate) or plate == None
	assert isinstance(tracing, int)
	assert isinstance(rotate, bool)
	assert isinstance(sort_order, int) or isinstance(sort_order, float)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part.multi_mount('{1}', '{2}')".
	      format(indent, part._name, multi_mounts._name_get()))
	    trace_detail = 2

	# Store the *part* into *multi_mounts*:
	multi_mounts._part_set(part)

	# Do not do anthing until *cnc_mode* is set:
	ezcad = part._ezcad_get()
	if ezcad._cnc_mode:
	    ezcad._multi_mounts_list.append(multi_mounts)
	    if trace_detail >= 2:
		print("{0}multi-mounts: part='{1}' name='{2}'".
		  format(indent, part._name, multi_mounts._name_get()))

	    # Deal with *plate*:
	    if isinstance(plate, Plate):
		plate._multi_mounts_add(multi_mounts, rotate, sort_order)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part.multi_mount('{1}', '{2}')".
	      format(indent, part._name, multi_mounts._name_get()))

    def _multi_mount_process(self, name, multi_mounts, tracing=-1000000):
	""" *Part*: Process the *multi_mounts* for the *Part* object (i.e. *self*)
	    using *multi_mounts*.
	"""

	# Use *part* instead of *self*:
	part = self

	# Verify argument types:
	assert isinstance(name, str) and not ' ' in name, "name '{0}' has spaces".format(name)
	assert isinstance(multi_mounts, Multi_Mounts)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing < 0 and part._tracing >= 0:
	    tracing = part._tracing
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Part._multi_mount_process('{1}', '{2}', *)".
	      format(indent, part._name, name))
	    trace_detail = 3

	# Make sure we are in *cnc_mode*:
	ezcad = part._ezcad_get()
	assert ezcad._cnc_mode

	# Grab *tooling_plate* from *ezcad*:
	shop = ezcad._shop_get()
	tooling_plate = shop._tooling_plate_get()

	# Figure out if *multi_mounts* is supposed to be mounted on a *tooling_plate*:
	is_tooling_plate = multi_mounts._is_tooling_plate_get(tracing = tracing + 1)

	# Determine how the extra material is layed out relative to the (0, 0) *tooling_plate* hole:
	hole00_extra_start_bsw, hole00_extra_start_tne, dowel_pin00 = \
	  multi_mounts._extra_start_dowel_pin_get(tooling_plate, tracing = tracing + 1)
	if trace_detail >= 2:
	    print("{0}hole00_extra_start_bsw={1:i}".format(indent, hole00_extra_start_bsw))
	    print("{0}hole00_extra_start_tne={1:i}".format(indent, hole00_extra_start_tne))
	    print("{0}dowel_pin00={1:i}".format(indent, dowel_pin00))

	# Define some constants:
	zero = L()
	origin = P(zero, zero, zero)

	# Grab some values from *tooling_plate*:
	tooling_plate_hole_pitch = tooling_plate._hole_pitch_get()
	tooling_plate_rows       = tooling_plate._rows_get()
	tooling_plate_columns    = tooling_plate._columns_get()
	tooling_plate_dx         = tooling_plate._dx_get()
	tooling_plate_dy         = tooling_plate._dy_get()
	tooling_plate_dz         = tooling_plate._dz_get()
	if trace_detail >= 2:
	    print("{0}tooling_plate: hole_pitch={1:i} columns={2} rows={3}".
	      format(indent, tooling_plate_hole_pitch, tooling_plate_columns, tooling_plate_rows))
	    print("{0}tooling_plate: dx={1:i} dy={2:i} dz={3:i}".
	      format(indent, tooling_plate_dx, tooling_plate_dy, tooling_plate_dz))

	# Now compute some more stuff about the *tooling_plate*:
	tooling_plate_columns_dx = (tooling_plate_columns - 1) * tooling_plate_hole_pitch
	tooling_plate_rows_dy    = (tooling_plate_rows - 1) * tooling_plate_hole_pitch
	tooling_plate_extra_dx   = tooling_plate_dx - tooling_plate_columns_dx
	tooling_plate_extra_dy   = tooling_plate_dy - tooling_plate_rows_dy
	tooling_plate_spacer_dz  = tooling_plate._spacer_dz_get()
	if is_tooling_plate:
	    hole_0_0 = P(tooling_plate_extra_dx / 2, -tooling_plate_extra_dy / 2, zero)
	else:
	    hole_0_0 = P(-hole00_extra_start_bsw.x, -hole00_extra_start_tne.y, zero)
	if trace_detail >= 2:
	    print("{0}tooling_plate: columns_dx={1:i} rows_dy={2:i}".
	      format(indent, tooling_plate_columns_dx, tooling_plate_rows_dy))
	    print("{0}tooling_plate: extra_dx={1:i} extra_dy={1:i}".
	      format(indent, tooling_plate_extra_dx, tooling_plate_extra_dy))
	    print("{0}tooling_plate: hole_0_0={1:i}".format(indent, hole_0_0))
	    print("{0}tooling_plate: spacer_dz={1:i}".format(indent, tooling_plate_spacer_dz))

	# Grab some values from *vice*:
	vice = shop._vice_get()
	parallels = vice._parallels_get()
	jaw_volume = vice._jaw_volume_get()

	# Unpack the *vice* information some more:
	jaw_dx = jaw_volume.x
	jaw_dy = jaw_volume.y
	jaw_dz = jaw_volume.z
	jaws_spread = tooling_plate_dy if is_tooling_plate else zero

	# Initialize some values prior to the loop which fills them in:
	combined_multi_mount = None
	combined_operations = None
	spacers = []
	extra_bounding_box = Bounding_Box()

	# Deal with tooling plate mounting vs. vice mounting:
	tooling_plate_present = multi_mounts._is_tooling_plate_get(tracing = tracing + 1)
        if True:
	    # Scan across each *multi_mount*:
	    for index, multi_mount in enumerate(multi_mounts._multi_mounts_get()):
		# Grab some values from *multi_mount*:
		multi_mount_column, multi_mount_row = multi_mount._column_row_get()
		multi_mount_mount = multi_mount._mount_get()
		multi_mount_mount_operations = multi_mount._mount_operations_get()
		multi_mount_part = multi_mount._part_get()
		multi_mount_rotate = multi_mount._rotate_get()

		# Grab some values from *multi_mount_mount*:
		tooling_plate_holes = multi_mount_mount._tooling_plate_holes_get()
		mount_name = multi_mount_mount._name_get()
		mount_top_surface_transform = multi_mount_mount._top_surface_transform_get()

		# Compute distance between upper left tooling hole and *center* of
		# *multi_mount_part*:
		assert len(tooling_plate_holes) > 0
		sorted_columns = sorted([ column_row[0] for column_row in tooling_plate_holes ])
		sorted_rows =    sorted([ column_row[1] for column_row in tooling_plate_holes ])
		unique_columns = tuple(sorted(tuple(set(sorted_columns))))
		unique_rows =    tuple(sorted(tuple(set(sorted_rows))))
		minimum_column = sorted_columns[0]
		maximum_column = sorted_columns[-1]
		delta_columns  = maximum_column - minimum_column
		minimum_row = sorted_rows[0]
		maximum_row = sorted_rows[-1]
		delta_rows  = maximum_row - minimum_row
		delta_x   = delta_columns      * tooling_plate_hole_pitch
		delta_y   = delta_rows         * tooling_plate_hole_pitch
		column_dx = multi_mount_column * tooling_plate_hole_pitch
		row_dy    = multi_mount_row    * tooling_plate_hole_pitch

		# Create the *rotated_top_surface_transform* which places *multi_mount_part*
		# such that its top surface is centered under the CNC machine origin and
                # rotated in the correct orientation:
		one = L(mm=1.0)
		z_axis = P(zero, zero, one)
		top_surface_transform = multi_mount_mount._top_surface_transform_get()
		rotated_top_surface_transform = \
		  top_surface_transform.rotate("Multi Mount Rotate", z_axis, multi_mount_rotate)

		# Convert *rotate* from an *Angle* to an *int*:
		rotate_degrees = int(multi_mount_rotate.degrees())
		assert rotate_degrees % 90 == 0

		# Now compute both *offset* and *center* that specifies the offset from
		# *hole_0_0* to the upper left mounting hole for *multi_mount_part*.  *center*
		# specifies the distance from the upper left mounting hole to the center of
                # *multi_mount_part*:
		offset = P(column_dx, -row_dy, zero)
		if rotate_degrees % 180 == 0:
		    xy_center = P(delta_x/2, -delta_y/2, zero)
		else:
		    xy_center = P(delta_y/2, -delta_x/2, zero)

		# Figure out the *part* dimensions:
		part_bsw = multi_mount_part.bsw
		part_tne = multi_mount_part.tne
		transformed_part_bsw = rotated_top_surface_transform * part_bsw
		transformed_part_tne = rotated_top_surface_transform * part_tne
		top_surface_bsw, top_surface_tne = \
		  transformed_part_bsw.minimum_maximum(transformed_part_tne)
		part_dx = top_surface_tne.x - top_surface_bsw.x
		part_dy = top_surface_tne.y - top_surface_bsw.y
		part_dz = top_surface_tne.z - top_surface_bsw.z
		if trace_detail >= 2:
		    print("{0}[{1}]: sorted_columns={2} sorted_rows={3}".
		      format(indent, index, sorted_columns, sorted_rows))
		    print("{0}[{1}]: unique_columns={2} unique_rows={3}".
		      format(indent, index, unique_columns, unique_rows))
		    print("{0}[{1}]: mount_top_surface_transform={2:v}".
		      format(indent, index, mount_top_surface_transform))
		    print("{0}[{1}]: part_bsw={2:i} part_tne={3:i}".
		      format(indent, index, part_bsw, part_tne))
		    print("{0}[{1}]: transformed_part_bsw={2:i} transformed_part_tne={3:i}".
		      format(indent, index, transformed_part_bsw, transformed_part_tne))
		    print("{0}[{1}]: top_surface_bsw={2:i} top_surface_tne={3:i}".
		      format(indent, index, top_surface_bsw, top_surface_tne))
		    print("{0}[{1}]: part_dx={2:i} part_dy={3:i} part_dz={2:i}".
		      format(indent, index, part_dx, part_dy, part_dz))

		# Find the smallest parallel that just clears the top of the *vice* for *part_dz*:
		if tooling_plate_present:
		    selected_parallel_height = tooling_plate._parallels_height_get()
		else:
		    selected_parallel_height = parallels._select(part_dz, vice)
		if trace_detail >= 2:
		    print("{0}[{1}]: selected_parallel_height={2:i}".
		      format(indent, index, selected_parallel_height))
		    print("{0}[{1}]: tooling_plate_dz={2:i}".
		      format(indent, index, tooling_plate_dz))
		    print("{0}[{1}]: tooling_plate_spacer_dz={2:i}".
		      format(indent, index, tooling_plate_spacer_dz))
		    print("{0}[{1}]: part_dz={2:i}".
		      format(indent, index, part_dz))

		# Figure out how much extra material there is:
		extra_start_bsw, extra_start_tne = multi_mount_mount._extra_start_get()
		transformed_extra_start_bsw = rotated_top_surface_transform * extra_start_bsw
		transformed_extra_start_tne = rotated_top_surface_transform * extra_start_tne
		fixed_extra_start_bsw, fixed_extra_start_tne = \
		  transformed_extra_start_bsw.minimum_maximum(transformed_extra_start_tne)
		extra_top_z = fixed_extra_start_tne.z - top_surface_tne.z
		if trace_detail >= 2:
		    print("{0}part.name='{1}'".format(indent, part._name_get()))
		    print("{0}part_bsw={1:i}".format(indent, part_bsw))
		    print("{0}part_tne={1:i}".format(indent, part_tne))
		    print("{0}transformed_part_bsw={1:i}".format(indent, transformed_part_bsw))
		    print("{0}transformed_part_tne={1:i}".format(indent, transformed_part_tne))
		    print("{0}extra_start_bsw={1:i}".format(indent, extra_start_bsw))
		    print("{0}extra_start_tne={1:i}".format(indent, extra_start_tne))
		    print("{0}transformed_extra_start_bsw={1:i}".
		      format(indent, transformed_extra_start_bsw))
		    print("{0}transformed_extra_start_tne={1:i}".
		      format(indent, transformed_extra_start_tne))
		    print("{0}fixed_extra_start_bsw={1:i}".format(indent, fixed_extra_start_bsw))
		    print("{0}fixed_extra_start_tne={1:i}".format(indent, fixed_extra_start_tne))
		    print("{0}extra_top_z={1:i}".format(indent, extra_top_z))

		# Figure out where *top_surface_z* is and compute the *top_surface* point:
		if is_tooling_plate:
		    top_surface_z = selected_parallel_height + \
		      tooling_plate_dz + tooling_plate_spacer_dz + part_dz
		else:
		    top_surface_z = selected_parallel_height + part_dz
		top_surface = P(zero, zero, top_surface_z)
		top_surface_z += extra_top_z
		if trace_detail >= 2:
		    print("{0}[{1}]: selected_parallel_height={2:i}".
		      format(indent, index, selected_parallel_height))
		    print("{0}[{1}]: tooling_plate_dz={2:i}".
		      format(indent, index, tooling_plate_dz))
		    print("{0}[{1}]: tooling_plate_spacer_dz={2:i}".
		      format(indent, index, tooling_plate_spacer_dz))
		    print("{0}[{1}]: part_dz={2:i}".
		      format(indent, index, part_dz))
		    print("{0}[{1}]: top_surface_z={2:i}".
		      format(indent, index, top_surface_z))

		# Compute *combined_mount_translate_point* which specifies how to map the
		# part from the machine top surface to the correct location on the multi-mount:
		combined_mount_translate_point = hole_0_0 + offset + xy_center + top_surface
		if trace_detail >= 2:
		    print("{0}[{1}]: multi_mount_part='{2}' mount='{3}'".
		      format(indent, index, multi_mount_part._name_get(), mount_name))
		    print("{0}[{1}]: min_col={2} max_col={3} min_row={4} max_row={5}".format(indent,
		      index, minimum_column, maximum_column, minimum_row, maximum_row))
		    print("{0}[{1}]: delta_rows={2} delta_columns={3} delta_x={4:i} delta_y={5:i}".
		      format(indent, index, delta_rows, delta_columns, delta_x, delta_y))
		    print("{0}[{1}]: column_dx={2:i} row_dy={3:i} xy_center={4:i} offset={5:i}".
		      format(indent, index, column_dx, row_dy, xy_center, offset))
		    print("{0}[{1}]: hole_0_0={2:i}".format(indent, index, hole_0_0))
		    print("{0}[{1}]: offset={2:i}".format(indent, index, offset))
		    print("{0}[{1}]: xy_center={2:i}".format(indent, index, xy_center))
		    print("{0}[{1}]: top_surface={2:i}".format(indent, index, top_surface))
		    print("{0}[{1}]: combined_mount_translate_point={2:i}".
		      format(indent, index, combined_mount_translate_point))

		# Create *combined_mount* which is a *Mount* object that will place
                # *multi_mount_part* properly onto the tool plate:
		combined_mount_name = "{0}[{1}]".format(name, index)
		top_surface_safe_z = top_surface_z
		xy_rapid_safe_z = top_surface_safe_z + L(inch=0.500)
		combined_mount = Mount(combined_mount_name, multi_mount_part,
		  rotated_top_surface_transform, combined_mount_translate_point,
		  top_surface_safe_z, xy_rapid_safe_z, True, selected_parallel_height)
		multi_mount._combined_mount_set(combined_mount)

		# Update the *extra_bounding_box*:
		combined_mount_cnc_transform = combined_mount._cnc_transform_get()
		extra_start_bsw, extra_start_tne = multi_mount_mount._extra_start_get()
		combined_mount._extra_start_set(extra_start_bsw, extra_start_tne)
		cnc_extra_start_bsw = combined_mount_cnc_transform * extra_start_bsw
		cnc_extra_start_tne = combined_mount_cnc_transform * extra_start_tne
		extra_bounding_box.point_expand(cnc_extra_start_bsw)
		extra_bounding_box.point_expand(cnc_extra_start_tne)
		if not is_tooling_plate:
		    extra_dy = extra_bounding_box.tne_get().y - extra_bounding_box.bsw_get().y
		    #extra_dy = (cnc_extra_start_tne.y - cnc_extra_start_bsw.y).absolute()
		    jaws_spread = jaws_spread.maximum(extra_dy)
		if trace_detail >= 2:
		    print("{0}[{1}]: combined_mount.part='{2}' combined_mount.name='{3}'".
		      format(indent, index,
		      combined_mount._part_get()._name_get(), combined_mount._name_get()))
		    print("{0}[{1}]: extra_start_bsw={2:i} extra_start_tne={3:i}".
		      format(indent, index, extra_start_bsw, extra_start_tne))
		    print("{0}[{1}]: combined_mount.cnc_transform={2:v}".
		      format(indent, index, combined_mount_cnc_transform))
		    print("{0}[{1}]: cnc_extra_start_bsw={2:i} cnc_extra_start_tne={3:i}".
		      format(indent, index, cnc_extra_start_bsw, cnc_extra_start_tne))

		# Generate some reasonble looking *spacers* for the CNC visualization:
		if rotate_degrees % 180 == 90:
		    unique_columns, unique_rows = unique_rows, unique_columns
		if len(unique_columns) >= len(unique_rows):
		    for column in unique_columns:
			spacers.append(
			  (column + multi_mount_column, unique_rows[0]  + multi_mount_row,
                           column + multi_mount_column, unique_rows[-1] + multi_mount_row) )
		else:
		    for row in unique_rows:
			spacers.append(
			  (unique_columns[0]  + multi_mount_column, row + multi_mount_row,
                           unique_columns[-1] + multi_mount_column, row + multi_mount_row ) )
		if trace_detail >= 2:
		    print("{0}[{1}]: unique_columns={2} unique_rows={3} spacers={4}".
                          format(indent, index, unique_columns, unique_rows, spacers))
		    
		# The first time through the loop, both *combined_multi_mount* and
                # *combined_operations* are initialized.  *combined_multi_mount* is the
		# *Mount* object that is used for the *Operation_Multi_Mount* and
		# *Operation_Dowel_Pin* operations.  Also, the final extra material is
		# eventually stuffed into *combined_multi_mount*:
		if combined_multi_mount == None:
		    # *null_transform* is not actually used by the either the
		    # *Operation_Mulit_Mount* or *Operation_Dowel_Pin* operations:
		    null_transform = Transform()

		    # Create *combined_multi_mount*:
		    origin = P(zero, zero, zero)
		    combined_multi_mount_name = "{0}_Combined_Multi_Mount".format(name)
		    combined_multi_mount = Mount(combined_multi_mount_name,
		      part, null_transform, origin, top_surface_safe_z, xy_rapid_safe_z, True,
		      selected_parallel_height)
		    multi_mounts._combined_multi_mount_set(combined_multi_mount)

		    # Create *combined_operations*:
		    combined_operations = part._mount_register(combined_multi_mount)

		# Copy the operations from *multi_mount_mount_operations* to *combined_operations*
		# skipping the *Operation_Mount* and *Operation_Dowel_Pin* operations:
		combined_operations._extend(combined_mount,
		  multi_mount_mount_operations, "MD", tracing = tracing + 1)
		combined_operations._multi_mounts_set(multi_mounts)
		#print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
		#combined_operations._extend(combined_mount,
		#  multi_mount_mount_operations, "MD", tracing = 5)

	    # Create *operation_multi_mount* for eventual prepending to *combined_operations*:
	    comment = ""
	    operation_multi_mount = Operation_Multi_Mount(part, multi_mounts,
	      comment, vice, jaws_spread, selected_parallel_height, tooling_plate, spacers,
	      is_tooling_plate, tracing = tracing + 1)

	    # Stuff the bounding box into *combined_multi_mount* for CNC visualization:
	    extra_start_bsw = extra_bounding_box.bsw_get()
	    extra_start_tne = extra_bounding_box.tne_get()
	    combined_multi_mount._extra_start_set(extra_start_bsw, extra_start_tne)
	    if trace_detail >= 3:
		print("{0}mount_name='{1}' extra_start_bsw={2:i} extra_start_tne={3:i}".
		  format(indent, combined_multi_mount._name_get(),extra_start_bsw, extra_start_tne))
		combined_operations._show("Multi Mount operations")
	    if not is_tooling_plate:
		extra_dy = extra_bounding_box.tne_get().y - extra_bounding_box.bsw_get().y
		operation_multi_mount._jaws_spread_set(extra_dy)

	    # Create *operation_dowel_pin* for eventual prepending to *combined_operations*:
	    dowel_pin_point = extra_bounding_box.bw_get()
	    preferred_diameter = L(inch="3/8")
	    dowel_pin_tool = \
	      part._tools_dowel_pin_search(preferred_diameter, tracing = tracing)
	    assert isinstance(dowel_pin_tool, Tool_Dowel_Pin)
	    plunge_point = dowel_pin_point + P(L(inch=-0.700), zero, zero)
	    operation_dowel_pin = Operation_Dowel_Pin(part, "Multi Mount Dowel Pin",
	      0, dowel_pin_tool, Operation.ORDER_DOWEL_PIN, None,
	      dowel_pin_tool._feed_speed_get(), dowel_pin_tool._spindle_speed_get(),
	      dowel_pin_tool._diameter_get(), dowel_pin_point, plunge_point,
	      tracing = tracing + 1)

	    # Prepend *operation_dowel_pin* and *operation_multi_mount* to
	    # *combined_operations* in reverse order so that the mount operation
	    # comes first:
	    combined_operations._prepend(combined_multi_mount,
	      operation_dowel_pin, tracing = tracing + 1)
	    combined_operations._prepend(combined_multi_mount,
	      operation_multi_mount, tracing = tracing + 1)

	    #temporary_bsw, temporary_tne = combined_multi_mount._extra_start_get()
	    #assert temporary_bsw == extra_start_bsw
	    #assert temporary_tne == extra_start_tne
	else:
	    assert False, \
	      "dx_dy mode not supported yet for Multi_Mounts '{0}'".format(multi_mounts._name_get())

	# For debugging:
	#assert False

	# Perform any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Part._multi_mount_process('{1}', '{2}', *)".
	      format(indent, part._name, name))

class Place:
    """ *Place*: Represents a part placement in the visualization. """

    def __init__(self, name, part, transform):
	""" *Place*: ...
	"""

	# Use *place* instead of *self*:
	place = self

	# Verify argument types:
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(part, Part)
	assert isinstance(transform, Transform)

	# Load values into *place*:
	place._name = name
	place._part = part
	place._transform = transform

    def _name_get(self):
	""" *Place*: Return the *name* of the *Place* object (i.e. *self*.)
	"""

        return self._name

    def _part_get(self):
	""" *Place*: Return the *Part* associated with the *Place* object (i.e. *self*.)
	"""

	return self._part

    def _transform_set(self, transform):
	""" *Place*: Replace the current *Transform* associated with the *Place* object
	    (i.e. *self*) with *transform*.
	"""

	# Use *place* instead of *self*:
	place = self

	# Verify argument types:
	assert isinstance(transform, Transform)

	# Update *transform* into *place*:
	place._transform = transform

    def _transform_get(self):
	""" *Place*: Return the *Transform* associated with the *Place* object (i.e. *self*.)
	"""

	return self._transform

class Plate:
    """ *Plate*: Represents a plate of material to be cut up into smaller plates for machining.
    """

    def __init__(self, name, dx, dy, dz, first_letter):
	""" *Plate*: Initialize the *Plate* object to be named *name* and have dimensions of
	    *dx* x *dy* x *dz* where *dx* and *dy* are the length and width and *dz* is the
	    thickness.
	"""

	# Use *plate* instead of *self*:
	plate = self

	# Verify the argument types:
	assert isinstance(name, str)
	assert isinstance(dx, L)
	assert isinstance(dy, L)
	assert isinstance(dz, L)
	assert isinstance(first_letter, str) and len(first_letter) == 1

	# Load *name*, *dx*, *dy* and *dz* into *plate*:
	plate._dx = dx
	plate._dy = dy
	plate._dz = dz
	plate._first_letter = first_letter
	plate._name = name
	plate._triples = []

    def _multi_mounts_add(self, multi_mounts, rotate, sort_order, cell_dx=0, cell_dy=0):
	""" *Plate*: Add a *multi_mount* to the *Plate* (i.e. *self*) multi-mounts list with a
	    a sort order of *sort_order* and a rotate of *True* or *False*.
	"""

	# Use *plate* instead of *self*:
	plate = self

	# Verify argument types:
	assert isinstance(multi_mounts, Multi_Mounts)
	assert isinstance(rotate, bool)
	assert isinstance(sort_order, int) or isinstance(sort_order, float)
	assert isinstance(cell_dx, int) and cell_dx >= 0
	assert isinstance(cell_dy, int) and cell_dy >= 0
	
	# Create a *triple* and append it to *triples*:
	triple = (multi_mounts, rotate, sort_order, cell_dx, cell_dy)
	triples = plate._triples
	triples.append(triple)
	
    def _name_get(self):
        """ *Plate*: Return the *name* associated with the *Plate* object (i.e. *self*):
	"""

	return self._name

    def _process(self, tracing=-1000000):
        """ *Plate*: Write the contents of the *Plate* object (i.e. *self*) out to a file.
	    This involves computing the content as well.
	"""

	# Use *plate* instead of *triples*:
	plate = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Plate.process('{1})".format(indent, plate._name))
	    trace_detail = 1

	# Define some constants:
	zero = L()
	cell_size = L(inch=0.250)
	thousandth = L(inch=0.001)

	# Sort *triples* by the sort order:
	triples = plate._triples
	triples.sort(key=lambda triple: triple[2])
	
	# Grab some values from *plate*:
	dx           = plate._dx
	dy           = plate._dy
	dz           = plate._dz
	name         = plate._name
	first_letter = plate._first_letter
	character_offset = ord(first_letter)

	# Compute the cells grid sizes:
	cells_dx = int((dx + thousandth) / cell_size)
	cells_dy = int((dy + thousandth) / cell_size)
	if trace_detail >= 1:
	    print("{0}cells_dx={1} cells_dy={2}".format(indent, cells_dx, cells_dy))

	
	# Create the *cells_grid* which is *cells_dx* by *cells_dy* in size:
	# Note 1: Start with an initial row of all -1's.
	initial_row = [ -1 for cell_x in range(cells_dx) ]
	# Note 2: The code `list(initial_row)` returns a shallow copy of *initial_row*.
	cells_grid = [ list(initial_row) for cell_y in range(cells_dy) ]
	# Note 3: *cells_grid* is accessed with Y index first (e.g. *cells_grid[y][x]* .)

	# Create the *empty_indices* which contains the X index for the first empty position
	# for each row in *cells_grid*.  Initialize everything to 0:
	empty_indices = [ 0 for cell_y in range(cells_dy) ]
	if trace_detail >= 1:
	    print("{0}len(empty_indices)={1}".format(indent, len(empty_indices)))

	# Open *plate_file* to contain the results:
	plate_name = plate._name
	plate_file_name = "/tmp/{0}_plate.txt".format(plate_name)
	if trace_detail >= 1:
		print("{0}plate_file_name='{1}'".format(indent, plate_file_name))
	with open(plate_file_name, "w") as plate_file:
	    # Write first line of newly opened *plate_file*:
	    plate_file.write("Plate: {0} {1:i} x {2:i} x {3:i}\n".format(name, dx, dy, dz))

	    # Visit each *triple* in *triples*:
	    for index, triple in enumerate(triples):
		# Pull the 5 values out of *triple*:
		multi_mounts = triple[0]
		rotate       = triple[1]
		sort_order   = triple[2]

		# Grab some values from *multi_mount*:
		multi_mounts_name = multi_mounts._name_get()
		program_number    = multi_mounts._program_number_get()
		cell_dx           = multi_mounts._cell_dx_get()
		cell_dy           = multi_mounts._cell_dy_get()
		if trace_detail >= 1:
		    print("{0}[{1}]: multi_mounts_name='{2}'".
		      format(indent, index, multi_mount_name))
		    print("{0}[{1}]: cell_dx={1} cell_dy={2}".
		      format(indent, index, cell_dx, cell_dy))
		    print("{0}[{1}]: program_number={2}".format(indent, index, program_number))
		    print("{0}[{1}]: rotate={2}".format(indent, index, rotate))
		    print("{0}[{1}]: sort_order={2}".format(indent, index, sort_order))
		
		# Figure out the size of the *multi_mounts*:
		combined_multi_mount = multi_mounts._combined_multi_mount_get()
		extra_start_bsw, extra_start_tne = combined_multi_mount._extra_start_get()
		extra_start_dx = extra_start_tne.x - extra_start_bsw.x
		extra_start_dy = extra_start_tne.y - extra_start_bsw.y
		extra_start_dz = extra_start_tne.z - extra_start_bsw.z
		if trace_detail >= 1:
		    print("{0}[{1}]: multi_mount_name='{2}'".
		     format(indent, index, multi_mount._name_get()))
		    print("{0}[{1}]: combined_multi_mount._name={2}".
		      format(indent, index, combined_multi_mount._name_get()))
		    print("{0}[{1}]: extra_start_bsw={2:i}".
		      format(indent, index, extra_start_bsw))
		    print("{0}[{1}]: extra_start_tne={2:i}".
		      format(indent, index, extra_start_tne))

		# Now round the sizes up to be in multiple of cells:
		multi_mount_cells_dx = int(math.ceil((extra_start_dx / cell_size)))
		multi_mount_cells_dy = int(math.ceil((extra_start_dy / cell_size)))
		multi_mount_cells_dz = int(math.ceil((extra_start_dz / cell_size)))

		# Deal with *rotate* by swapping *mutli_mount_cells_dx* and *mulit_mount_cells_dy*:
		if rotate:
		    multi_mount_cells_dx, multi_mount_cells_dy = \
		      multi_mount_cells_dy, multi_mount_cells_dx
		    swapped = " Swapped"
		if trace_detail >= 1:
		    rotated = " Rotated" if rotate else ""
		    print("{0}[{1}]: multi_mount_cells_dx='{2}'{3}".
		      format(indent, index, multi_mount_cells_dx, rotated))
		    print("{0}[{1}]: multi_mount_cells_dy='{2}'{3}".
		      format(indent, index, multi_mount_cells_dy, rotated))
		    print("{0}[{1}]: multi_mount_cells_dz='{2}'".
		      format(indent, index, multi_mount_cells_dz))

		# Now find a place in *cells_grid* to stuff this *multi_mount*:
		if trace_detail >= 3:
		    print("{0}empty_indices={1}".format(indent, empty_indices))
		match_found = False
		for y_search_index in range(cells_dy - multi_mount_cells_dy):
		    #print("y_search_index={0}".format(y_search_index))
		    empty_index = empty_indices[y_search_index]
		    if empty_index + multi_mount_cells_dx <= cells_dx + cell_dx:
			# We have room for at least one row, now check that all rows are available:
			if trace_detail >= 2:
			    print("{0}possible match at {1}".format(indent, y_search_index))
			match_found = True
			for y_index in range(y_search_index, y_search_index + multi_mount_cells_dy):
			    if empty_indices[y_index] > empty_index:
				match_found = False	
				break

			# If *match* found fill in the appropriate portion of *cells_grid*:
			if match_found:
			    if trace_detail >= 1:
				print("{0}We have a match at x={1} y={2}".
				  format(indent, empty_index, y_search_index))

			    # Perform the actual filling of *cells_grid*:
			    match_x_index = empty_index + cell_dx
			    match_y_index = y_search_index + cell_dy
			    new_empty_index = match_x_index + multi_mount_cells_dx
			    if trace_detail >= 2:
				print("{0}empty_index={1} new_empty_index={2} index={3}".
				  format(indent, empty_index, new_empty_index, index))
			    for y_index in range(match_y_index,
			      match_y_index + multi_mount_cells_dy):
				empty_indices[y_index] = new_empty_index
				row = cells_grid[y_index]
				for x_index in range(empty_index + cell_dx,
				  new_empty_index + cell_dx):
				    # *index* corresonds to *triple* index:
				    row[x_index] = index

			    # For extreme tracing, print out the *cells_grid*:
			    if trace_detail >= 3:
				debug_file = sys.stderr
				for yy_index in range(cells_dy):
				    row = cells_grid[yy_index]
				    for xx_index in range(cells_dx):
					grid_value = row[xx_index]
					character = chr(grid_value + character_offset) \
					  if grid_value >= 0 else '.'
					debug_file.write(character)
				    debug_file.write("\n")

			    # Move onto next *multi_mount*:
			    break

		# Generate the initial *summary_line* for this *triple*:
		summary_line = \
		  "[{0}]: {1:>3} {2:>5} {3:5.2f} x {4:5.2f} x {5:3.2f} {6} {7:<30}".format(
		  chr(index + character_offset), sort_order, program_number,
		  (cell_size * multi_mount_cells_dx).inches(),
		  (cell_size * multi_mount_cells_dy).inches(),
		  extra_start_dz.inches(),
		  'R' if rotate else '-', multi_mounts_name)

		# Inform use when no *match_found*:
		if match_found:
		    summary_line += " ({0:5.2f}:{1:5.2f}) ({2:5.2f}:{3:5.2f})".format(
                      (cell_size * match_x_index).inches(),
		      (cell_size * match_y_index).inches(),
                      (cell_size * (match_x_index + multi_mount_cells_dx)).inches(),
		      (cell_size * (match_y_index + multi_mount_cells_dy)).inches())
		else:
		    print("[{0}]: Could not find a place for '{1}' in plate '{2}'".
		      format(chr(index + character_offset), multi_mounts_name, plate_name))

		# Write *summary_line* out to *plate_file*:
		plate_file.write(summary_line)
		plate_file.write('\n')
		if trace_detail >= 1:
		    print("{0}[{1}]: Sumary_line='{2}'".format(indent, index, summary_line))

		# List the *Part*'s included in the *multi_mounts* and their associated
		# extra tasks beyond drilling tooling plate mount holes and doing the initial
		# top level machining:
		multi_mounts_list = multi_mounts._multi_mounts_get()
		assert isinstance(multi_mounts_list, list)
		for multi_mount_index, multi_mount in enumerate(multi_mounts_list):
		    part = multi_mount._part_get()
		    part_name = part._name_get()
		    column, row = multi_mount._column_row_get()
		    part_summary = "  [{0}, {1}] {2}".format(column, row, part_name)

		    # Now sweep the the *program_numbers* of *part* calling out the
		    # additional mounts that need to be performed:
		    program_numbers = part._program_numbers_get()
		    assert len(program_numbers) > 0
		    program_number0 = program_numbers[0]
		    for program_number in program_numbers:
			base_program_number = program_number - program_number0
			if base_program_number % 10 == 0 and \
			  (base_program_number == 0 or base_program_number >= 20):
			    part_summary += " " + str(program_number)
		    if trace_detail >= 2:
			print("{0}{1}".format(indent, part_summary))

		    # Write the *part_sumary* line out to *plate_file*:
		    plate_file.write(part_summary)
		    plate_file.write('\n')

	    # Print out the *cells_grid*:
	    for y_index in range(cells_dy):
		row = cells_grid[y_index]
		for x_index in range(cells_dx):
		    grid_value = row[x_index]
		    character = chr(grid_value + character_offset) if grid_value >= 0 else '.'
		    plate_file.write(character)
		plate_file.write("\n")

	# Wrap-up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Plate.process('{1})".format(indent, plate._name))

class Fastener(Part):
    """ *Fastener*: Represent a fastener (i.e. screw, bolt, etc.) to attach parts.  """
	 
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
	""" *Fastener*: Initialize the *Fastener* object (i.e. *self*) to have a parent of *up*
	    and a name of *name*.
	"""

	# Use *fastener* instead of *self*
	fastener = self

	# Verify argument types:
	assert isinstance(up, Part) or up == None
	assert isinstance(name, str)

	# Initialize the *Part* super class:
	Part.__init__(self, up, name)
	
	# Disable CNC generation:
	fastener.cnc_suppress()

	zero = L()
	fastener.comment_s = name
	fastener.color = Color("red")
	fastener.material = Material("steel", "stainless")
	fastener.start_p = P()
	fastener.end_p = P()
	fastener.flags_s = ""
	fastener.major_diameter_l = zero
	fastener.pitch_l = zero
	fastener.thread75_l = zero
	fastener.thread50_l = zero
	fastener.close_fit_l = zero
	fastener.free_fit_l = zero	
	fastener.flat_head_diameter_l = zero
	fastener.hex_insert_b = False
	fastener.nut_height_l = zero
	fastener.hex_nut_edge_width_l = zero
	fastener.hex_nut_tip_width_l = zero
	fastener.sides_angle_a = Angle()
	fastener.flat_head_point_angle_a = Angle()

    def configure(self, start, end, flags,
      head_washer_diameter = None, tail_washer_diameter = None, sides_angle = None,
      tracing=-1000000):
     	""" *Fastener*: Confitugure the *Fastener* object (i.e. *self*).
	    *comment* shows up in generated CNC.  *material* specifies the fastener material.
	    *color* specifies the rendering color.  *start* and *end* specify the end points
	    of the fastener.  *flags* specifies one or more flag options separated by a
	    colon (':').  The basic fastener class flags are "#0-80", "#1-72", "#2-56",
	    "#4-40", "#6-32", "#10-24", "M3x.05", etc.  {more description needed.}
	"""

	# Use *fastener* instead of *self:
	fastener = self

	# Verify argument types:
	assert isinstance(start, P)
	assert isinstance(end, P)
	assert isinstance(flags, str)
	assert isinstance(head_washer_diameter, L) or head_washer_diameter == None
	assert isinstance(tail_washer_diameter, L) or tail_washer_diameter == None
	assert isinstance(sides_angle, Angle) or sides_angle == None

	# Perform an requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Fastener.configure('{1}', {2:i}, {3:i}, hwd={4}, twd={5}, sa={6})".
	      format(indent, fastener._name_get(), start, end, flags,
	      head_washer_diameter, tail_washer_diameter, sides_angle))
	    trace_detail = 1

	major_diameter = None
	if isinstance(sides_angle, Angle):
	    fastener.sides_angle_a = sides_angle
	fastener.flags_s = flags
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
		fastener.flat_head_point_angle_a = Angle(deg = 81)
	    elif flag == "hi":
		# Hex insert
		fastener.hex_insert_b = True
	    else:
		assert False, "Unrecognized flag ('{0}') in flags('{1}')".format(flag, flags)

	assert isinstance(major_diameter, L), "No fastener size specified (e.g '#4-40' or 'M3x.05')"

	# Load stuff into *fastener*:
	fastener.major_diameter_l = major_diameter
	fastener.pitch_l = pitch
	fastener.thread75_l = thread75
	fastener.thread50_l = thread50
	fastener.close_fit_l = close_fit
	fastener.free_fit_l = free_fit
	fastener.nut_height_l = nut_height
	fastener.hex_nut_edge_width_l = hex_nut_edge_width
	fastener.hex_nut_tip_width_l = hex_nut_tip_width
	fastener.flat_head_diameter_l = flat_head_diameter
	fastener.start_p = start
	fastener.end_p = end

	# Wrap up any *tracing*:
	if tracing >= 0:
	    print("{0}<=Fastener.configure('{1}', {2:i}, {3:i}, hwd={4}, twd={5}, sa={6})".
	      format(indent, fastener._name_get(), start, end, flags,
	      head_washer_diameter, tail_washer_diameter, sides_angle))

    def construct(self):
	""" *Fastener*: """

	# Use *fastener* inested of *self*:
	fastener = self

	assert fastener.comment_s != "NO_COMMENT", \
	  "Fastener name is not set (not configured!)"
	diameter = fastener.thread75_l
	fastener.cylinder(fastener.comment_s, fastener.material, fastener.color,
	  diameter, diameter, fastener.start_p, fastener.end_p, 16, Angle(deg=0.0), "", "")

    def length_get(self):
        """ *Fastener*: Return the length of the *Fastener* object (i.e. *self*.) """

	# Use *fastener* instead of *self*:
	fastener = self

	# Return the length:
	return fastener.start_p.distance(fastener.end_p)


    def nut_ledge(self, part = None, flags = ""):
	""" *Fastener*: Cut out ledge for a screw and a nut. """

	# Check argument types:
	none_type = type(None)
	assert isinstance(part, Part)
	assert isinstance(flags, str)

	if part._ezcad._cnc_mode:
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
	
    def _fasten(self, comment, part, select, transform=None, tracing=-1000000):
	""" *Fastener*: Use  the *Fastener* object (i.e. *self*) to drill a hole in *part*
	    using its current mount for *part*.  *select* should be "thread* for a threaded
	    hole, "close* for a close fit hole, "free" for a looser fit hole, "access"
	    for a free hole with room to snuggle the screw head into, "set_screw" for
	    landing suitable for a set screw.
	"""

	# Use *fastener* instead of *self*:
	fastener = self

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(part, Part)
	assert isinstance(select, str) and \
	  select in ("thread", "close", "loose", "thread_access", "close_access", "set_screw")
	assert isinstance(transform, Transform) or transform == None
	assert isinstance(tracing, int)

	# Perform any requested *tracing*;
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Fastener._fasten('{1}', '{2}', '{3}', '{4}')".
	      format(indent, fastener._name, comment, part._name_get(), select))
	    trace_detail = 1

	# No need to do anything until we are past dimensions mode:
	ezcad = part._ezcad_get()
	if ezcad._stl_mode or ezcad._cnc_mode:
	    # Process the *select* argument:
	    flags = "t"
	    if select == "thread":
		#FIXME: We need to be smart and distinguish between hard and soft materials:
		#material = part._material
		#generic = material._generic
		#if generic == "steel" or generic == "iron":
		#    diameter = fastener.thread50_l
		#else:
		#    diameter = fastener.thread75_l
	        diameter = fastener.thread75_l
	    elif select == "close":
		diameter = fastener.close_fit_l
	    elif select == "free":
		diameter = fastener.free_fit_l
	    elif select == "thread_access":
		diameter = fastener.thread75_l
	    elif select == "close_access":
		diameter = fastener.close_fit_l
	    elif select == "set_screw":
		diameter = 2 * fastener.major_diameter_l
		flags = "f"
	    else:
		assert False, \
		  "select='{0}' not 'thread', 'close', 'loose', 'thread_access' or 'close_access'".\
                  format(select)

	    # Grab *start* and *end* from *fastener*:
	    start = fastener.start_p
	    end = fastener.end_p
	    assert isinstance(start, P), \
	      "No start point specified for fastener '{0}'".fastener._name
	    assert isinstance(end, P), \
	      "No end point specified for fastener '{0}'".fastener._name
	    if trace_detail >= 2:
		print("{0}start={1:i} end={2:i}".format(indent, start, end))

	    # Deal with *Transform* if it is specified:
	    if isinstance(transform, Transform):
		start = transform * start
		end = transform * end
	    if trace_detail >= 1:
		print("{0}start={1:i} end={2:i}".format(indent, start, end))

	    # Clip *start* and *end* to be within the bounding box of *part* to create
	    # *new_start* and *new_end*:
	    part_bsw = part.bsw
	    part_tne = part.tne
	    if trace_detail >= 2:
		print("{0}part_bsw={1:i} part_tne={2:i}".format(indent, part_bsw, part_tne))

	    # Clamp the start/end X coordinate to be between *part_bsw_x* and *part_tne_x*:
	    part_bsw_x = part.bsw.x
	    part_tne_x = part.tne.x
	    new_start_x = start.x.maximum(part_bsw_x).minimum(part_tne_x)
	    new_end_x   =   end.x.maximum(part_bsw_x).minimum(part_tne_x)
	    assert part_bsw_x <= new_start_x <= part_tne_x
	    assert part_bsw_x <= new_end_x   <= part_tne_x
				
	    # Clamp the start/end Y coordinate to be between *part_bsw_y* and *part_tne_y*:
	    part_bsw_y = part_bsw.y
	    part_tne_y = part_tne.y
	    new_start_y = start.y.maximum(part_bsw_y).minimum(part_tne_y)
	    new_end_y   =   end.y.maximum(part_bsw_y).minimum(part_tne_y)
	    assert part_bsw_y <= new_start_y <= part_tne_y
	    assert part_bsw_y <= new_end_y   <= part_tne_y

	    # Clamp the start/end Z coordinate to be between *part_bsw_z* and *part_tne_z*:
	    part_bsw_z = part.bsw.z
	    part_tne_z = part.tne.z
	    new_start_z = start.z.maximum(part_bsw_z).minimum(part_tne_z)
	    new_end_z   =   end.z.maximum(part_bsw_z).minimum(part_tne_z)
	    assert part_bsw_z <= new_start_z <= part_tne_z
	    assert part_bsw_z <= new_end_z   <= part_tne_z

	    # Create *new_start* and *new_end*:
	    new_start   = P(new_start_x, new_start_y, new_start_z)
	    new_end     = P(new_end_x,   new_end_y,   new_end_z)
	    if trace_detail >= 2:
		print("{0}new_start={1:i} new_end={2:i}".format(indent, new_start, new_end))

	    # Now we can drill the hole in *part* with *new_start* and *new_end*:
	    hole_name = "{0} Hole".format(comment)
	    zero = L()
	    #assert new_start.distance(new_end) > zero, \
	    #  "Fastener.fasten: Zero length '{0}': s={1:i} ns={2:i} e={3:i} ne={4:i}".format(
	    #  fastener.comment_s, start, new_start, end, new_end)
	    epsilon = L(inch=.00000001)
	    assert new_start.distance(new_end) > epsilon, \
	      "Fastener.fasten: Part '{0}' (center={1:i}) misses Fastener '{2}' (center={3:i}". \
	      format(fastener.comment_s, fastener.c, part.name_get(), part.c)
	    part.hole(hole_name, diameter, new_start, new_end, flags, tracing = tracing + 1)

	    # An "access" hole is drilled about 1/4 of the way in on to after the first
	    # hole is drilled:
	    if select == "thread_access" or select == "close_access":
		#print("Access hole {0}".format(fastener.comment_s))
		#print("Access hole: ns={0:i} ne={1:i} depth={2:i}".
		#  format(new_start, new_end, new_start.distance(new_end)))
		access_depth = (new_end - new_start) / 4
		access_start = new_start
		access_end   = new_start + access_depth
		#print("Access hole: as={0:i} ae={1:i} ad={2:i}".
		#  format(access_start, access_end, access_start.distance(access_end)))
		major_diameter = fastener.major_diameter_l
		part.hole("Access_{0}".format(hole_name),
		  2 * major_diameter, access_start, access_end, "t", tracing = tracing + 1)

	#FIXME: This code is a little old an needs some work!!!
	if False:
	    zero = L()
	    flat_head_point_angle = fastener.flat_head_point_angle_a
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
		# What we need is the distance CB.  Triangle ABC is a right
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

		flat_head_diameter = fastener.flat_head_diameter_l
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
		    #	 "Part[{0}]:start={1} end={2} normalize={3} fh_end={4}".
		    #    format(part._name,
		    #    start, end, normalized_direction, flat_head_end))
		    part.hole("Flat Head:" + fastener.comment_s,
		      2 * flat_head_diameter, flat_head_start,
		      flat_head_end, "", tracing = trace + 1)

	    if fastener.hex_insert_b:
		direction = end - start
		direction_length = direction.length()
		if direction_length <= zero and EZCAD3.update_count_get() == 0:
		    print("Fastener.drill(): part={0}: zero length screw".
		      format(part))
		else:
		    normalized_direction = direction.normalize()
		    nut_height = fastener.nut_height_l
		    insert_end = end - (normalized_direction * nut_height._mm)
		    #print("end={0} start={1} dir={2} dir_len={3} nut_hght={4}".
		    #  format(end, start, direction, direction_len, nut_hght))
		    #print("insert_end = {0}".format(insert_end))

		    part.hole("Hex Insert:" + fastener.comment_s,
		      fastener.hex_nut_tip_width_l, end, insert_end, "f",
	               sides=6, sides_angle=fastener.sides_angle_a)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Fastener._fasten('{1}', '{2}', '{3}', '{4}')".
	      format(indent, fastener._name, comment, part._name_get(), select))

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
	bounding_box = Bounding_Box()
	code._bounding_box = bounding_box # Bounding box that encloses all path trajectories
	code._code_stream = None	  # Code stream to write G-codes to
	code._command_started = False     # At begining of RS274 command line
	code._command_chunks = []         # RS274 line broken into space separated chunks
	code._comment_chunks = []         # RS274 comment broken into space separated chunks
	code._dxf = ""                    # Text for dxf file
	code._dxf_lines = []              # List of lines to write out to .dxf file
	code._dxf_x_offset = zero         # DXF X offset
	code._dxf_y_offset = zero         # DXF Y offset
	code._is_laser = False            # True if tool is a laser
	code._mount  = None               # *Mount* object to use for current batch of G-code
	code._mount_wrl_lines = None      # The lines to be written into the mount VRML file
	code._rapid_speed = Speed(in_per_sec = 75.0) # Should be read from Mill object
	code._time = Time()		  # Total Tool path time
	code._tool = None		  # Currently selected tool
	code._tool_change_point = None    # Tool change location relative to viceo origin
	code._tool_program_number = -1    # Tool program number
	code._tool_wrl_lines = None       # The lines to be written into the tool VRML file

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
	code._h = -1			# H tool offset index
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
	code._vrml_motion_color_rapid       = "green"
	code._vrml_motion_color_cutting     = "red"
	code._vrml_motion_color_retract     = "cyan"
	code._vrml_motion_color_tool_change = "orange"

	# Reset VRML fields here.  Note: routine acessses some of the G-code variable,
	# so it must be called last:
	code._vrml_reset()

    def _bounding_box_get(self):
	""" *Code*: Return the *Bounding_Box* associated with the *Code* object (i.e. *self*.)
	"""

	return self._bounding_box

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
	code._vrml_motion_color = "yellow"

    def _cnc_transform_set(self, cnc_transform):
	""" *Code*: Set the CNC transform for the *Code* object (i.e. *self*) to *cnc_transform*.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(cnc_transform, Transform)

	# Stuff *cnc_transform* into *code*:
	code._cnc_transform = cnc_transform

    def _command_end(self, vrml_suppress=False, tracing=-1000000):
	""" *Code*: End the current RS274 in the *Code* object (i.e. *self*). """

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(vrml_suppress, bool)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._command_end(*)".format(indent))
	    trace_detail = 2

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

	# Now we can tack the current point onto *vrml_lines_points*:
	if not vrml_suppress:
	    assert -1000000.0 <= code._x.millimeters() <= 1000000.0
	    assert -1000000.0 <= code._y.millimeters() <= 1000000.0
	    assert -1000000.0 <= code._z.millimeters() <= 1000000.0
	    code._vrml_points_append(code._vrml_motion_color, P(code._x, code._y, code._z),
	      tracing = tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Code._command_end(*)".format(indent))

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

    def _contour(self, contour, plunge_offset, contour_offset, tool_radius, clockwise,
      z, feed_speed, spindle_speed, tracing):
	""" *Code*: Generate the G-code for *contour* using the *Code* object (i.e. *self*).
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
	speed0 = Speed()
	hertz0 = Hertz()
	assert isinstance(contour, Contour)
	assert isinstance(plunge_offset, L) and plunge_offset >= zero
	assert isinstance(contour_offset, L) and contour_offset >= zero
	assert isinstance(tool_radius, L) and tool_radius >= zero
	assert isinstance(clockwise, bool)
	assert isinstance(z, L)
	assert isinstance(feed_speed, Speed) and feed_speed > speed0
	assert isinstance(spindle_speed, Hertz) and spindle_speed > hertz0
	assert isinstance(tracing, int)

	# Start performing any *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._contour(*, po={1:i} co={2:i} tr={3:i} cl={4} z={5:i} *)".format(
	      indent, plunge_offset, contour_offset, tool_radius, clockwise, z))
	    trace_detail = 1

	# Grab *bends* from *contour*:
	bends = contour._bends_get()
	bends_size = len(bends)

	# Add in *tool_radius* to the exising *contour_offset*:
	total_offset = contour_offset + tool_radius

	# Find the *first_outside_bend* in *bends*:
	first_outside_bend_index = -1
	first_outside_bend = None
	for index, bend in enumerate(bends):
	    if not bend._is_inside_get():
		first_outside_bend = bend
		first_outside_bend_index = index
		break
	assert first_outside_bend_index >= 0
	if trace_detail >= 1:
	    print("{0}first_outside_bend='{1}' first_outside_bend_index={2}".
	      format(indent, first_outside_bend._name_get(), first_outside_bend_index))

	# Now compute *plunge_tangent* and *start_tangent* where the tool respectively plunges
        # and first starts cutting the contour:
	first_outside_bend_radius = first_outside_bend._radius_get()
	start_offset = total_offset + first_outside_bend_radius
	plunge_arc_offset = start_offset + plunge_offset
	contour_is_clockwise = contour._is_clockwise_get()
	if contour_is_clockwise == clockwise:
	    start_tangent = first_outside_bend._incoming_tangent_compute(start_offset)
	    plunge_tangent = first_outside_bend._incoming_tangent_compute(plunge_arc_offset)
	else:
	    start_tangent = first_outside_bend._outgoing_tangent_compute(-start_offset)
	    plunge_tangent = first_outside_bend._outgoing_tangent_compute(-plunge_arc_offset)
	if trace_detail >= 1:
	    print("{0}plunge_x={1:i}, plunge_y={2:i}".
	      format(indent, plunge_tangent.x, plunge_tangent.y))
    
	# This routine starts and ends at (*plunge_tangent_x*, *plunge_tangent_y*).  If we are
        # not already at (*plunge_tangent_x*, *plunge_tangent_y*) we need to get there safely:
	plunge_tangent_x = plunge_tangent.x
	plunge_tangent_y = plunge_tangent.y
	code._xy_rapid(plunge_tangent_x, plunge_tangent_y)

	# Now we get to down to the correct Z level:
	code._z_feed(feed_speed/2, spindle_speed, z, "Contour")

	# Now feed to the start point using a circular motion so that there is no divit
	# at the plunge point:
	start_tangent_x = start_tangent.x
	start_tangent_y = start_tangent.y
	if contour_is_clockwise == clockwise:
	    code._xy_ccw_feed(feed_speed, spindle_speed,
	      tool_radius * 1.001, start_tangent_x, start_tangent_y)
	else:
	    code._xy_cw_feed(feed_speed, spindle_speed,
 	      tool_radius * 1.001, start_tangent_x, start_tangent_y)

	# Visit each *bend* in the correct order:
	if contour_is_clockwise == clockwise:
	    # Clockwise (climb) Cut:
	    if trace_detail >= 0:
		print("{0}clockwise".format(indent))
    
	    # Iterate *bends_size* + 1 times:
	    for bends_index in range(bends_size):

		# Fetch a {bend}:
		index = (bends_index + first_outside_bend_index) % bends_size
		bend = bends[index]
		bend_radius = bend._radius_get()

		if trace_detail >= 1:
		    bend_name = bend._name_get()
		    bend_point = bend._point_get()
		    print("{0}bend_name='{1}' bend_point='{2:i}".
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
		    code._xy_feed(feed_speed, spindle_speed, arc_start.x, arc_start.y, tracing + 1)
		    code._xy_cw_feed(feed_speed, spindle_speed,
		      arc_radius, arc_end.x, arc_end.y, tracing = tracing + 1)

		if trace_detail >= 1:
		    trace_line = "CW[{0}]: aradius={1:i} astart={2:i} aend={3:i}". \
		      format(bends_index, arc_radius, arc_start, arc_end)
		    print("{0}{1}".format(indent, trace_line))
		    code._comment(trace_line)
	else:
	    # Counter clockwise:
	    if trace_detail >= 1:
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

		if trace_detail >= 1:
		    trace_line = "CCW[{0}]: arc_radius={1:i} arc_start={2:i} arc_end={3:i}". \
		      format(bends_index, arc_radius, arc_start, arc_end)
		    print("{0}{1}".format(indent, trace_line))
		    code._comment(trace_line)
	    
	# Return back to the (*start_tangent.x*, *start_tangent.y*):
	code._xy_feed(feed_speed, spindle_speed, start_tangent_x, start_tangent_y)

	# Now feed to back to the plunge point in a circular motion to avoid diviting the
	# contour as we retract:
	if contour_is_clockwise == clockwise:
	    code._xy_ccw_feed(feed_speed, spindle_speed,
	      tool_radius * 1.001, plunge_tangent_x, plunge_tangent_y)
	else:
	    code._xy_cw_feed(feed_speed, spindle_speed,
	      tool_radius * 1.001, plunge_tangent_x, plunge_tangent_y)

	if tracing >= 0:
	    print("{0}=>Code._contour(*, po={1:i} co={2:i} tr={3:i} cl={4} z={5:i} *)".format(
	      ' ' * tracing, plunge_offset, contour_offset, tool_radius, clockwise, z))

    def _drill_cycle(self, g_complete, g_cycle, f, p, q, z_safe, r, s, x, y, z, tracing=-1000000):
	""" *Code*: Generate one drill cycle for the *Code* object (i.e. *self*.)  *g_complete*
	    must be either 98 or 99 for either a G89 or G99 completion cycle.  *g_cycle* must
	    specify one of the canned drill cycles (i.e. 73, 74, 76, 81, 82, 83, 84, 85, 86, 87,
	    88, or 89.)  *f* is the feed speed, *p* is pause time for some of the cycles.
	    *q* is used for peck drilling.  *r* is the retract level.  *s* is the spindle speed.
	    *x* and *y* specify where to drill in the X/Y plane.  *z* specifies the final
	    drill depth.
	"""

	# Use *code* instead of *self*:
	code = self
			
	# Verify argument types:
	assert isinstance(g_complete, int) and 98 <= g_complete <= 99
	assert isinstance(g_cycle, int) and \
	  g_cycle in (73, 74, 76, 81, 82, 83, 84, 85, 86, 87, 88, 89)
	assert isinstance(f, Speed)
	assert isinstance(p, Time) or p == None
	assert isinstance(q, L) or q == None
	assert isinstance(z_safe, L)
	assert isinstance(r, L) or r == None
	assert isinstance(s, Hertz)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(z, L)

	# Perform any requested *tracing*:
	trace_detail = -1
	deep_tracing = -1000000
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Code._drill_cycle(*, {1}, {2}, " +
	      "f={3:I}, p={4:s}, q={5} r={6:i}, x={7:i}, y={8:i}, z={9:i})").
	      format(indent, g_complete, g_cycle, f, p, q, r, x, y, z))
	    trace_detail = 3
	    deep_tracing = tracing + 1

	# Make sure we are at reasonable Z height:
	code._xy_rapid_safe_z_force(f, s, tracing = tracing + 1)

	# Remember *z_before* so we can set z to the correct height at the end of the cycle:
	z_before = code._z
	point_before = P(x, y, z_before)
	vrml_points = code._vrml_points
	code._vrml_points_append(code._vrml_motion_color_rapid, point_before, tracing = tracing + 1)

	# Grab *xy_rapid_safe_z* from *mount*:
	mount = code._mount
	assert isinstance(mount, Mount)
	xy_rapid_safe_z = mount._xy_rapid_safe_z_get()

	# The routine call can change anything, so force a reset.
	#code._reset()

	# Start the command:
	code._command_begin()

	# Tack on the call:
	command_chunks = code._command_chunks
	command_chunks.append("o{0} call".format(g_cycle))

	# Grab *x_before* and *y_before* needed for time estimation:
	x_before = code._x
	y_before = code._y

	# Tack on the remaining arguments:
	code._speed( "[]", f)			# F
	code._length("[]", x)			# X
	code._length("[]", y)			# Y
	code._length("[]", xy_rapid_safe_z)	# Z_TOP
	code._length("[]", z_safe)		# Z_SAFE
	code._length("[]", r)			# Z_RETRACT
	code._length("[]", z)			# Z_BOTTOM
	if isinstance(q, L):
	    code._length("[]", q)		# Z_PECK
	code._command_end(vrml_suppress=True, tracing=tracing+1)

	# When the drill cycle is done, it will be at (*x*, *y*, *xy_rapid_save_z*):
	code._x = x
	code._y = y
	code._z = xy_rapid_safe_z

	# Compute the X/Y rapid distance:
	dx = (x - x_before).millimeters()
	dy = (y - y_before).millimeters()
	xy_distance_mm = math.sqrt(dx*dx + dy*dy)
	xy_distance = L(mm=xy_distance_mm)

	# Compute an estimate for how long the drill operation will take:
	rapid_speed = code._rapid_speed
	time = code._time
	time += xy_distance / rapid_speed
	time += (xy_rapid_safe_z - r) / rapid_speed
	pecks = int(math.ceil((r - z) / q)) if isinstance(q, L) else 1
	peck_amount = (r - z) / pecks
	for peck_index in range(pecks):
	    peck_dz = (peck_index + 1) * peck_amount
	    time += peck_dz / f
	    time += peck_dz / rapid_speed
	time += (xy_rapid_safe_z - r) / rapid_speed
	code._time = time

	# *g_complete* specifies whether Z ends on *z_before* or *r*:
	if g_complete == 98:
	    z_after = z_before
	else:
            z_after = r
	code._z = z_after

	color_retract = code._vrml_motion_color_retract
	color_cutting = code._vrml_motion_color_cutting
	code._vrml_points_append(color_retract, P(x, y, r),       tracing=deep_tracing)
	code._vrml_points_append(color_cutting, P(x, y, z),       tracing=deep_tracing)
	code._vrml_points_append(color_cutting, P(x, y, r),       tracing=deep_tracing)
	code._vrml_points_append(color_retract, P(x, y, z_after), tracing=deep_tracing)

	if tracing >= 0:
	    print(("{0}<=Code._drill_cycle(*, {1}, {2}, " +
	      "f={3:I}, p={4:s}, q={5} r={6:i}, x={7:i}, y={8:i}, z={9:i})").
	      format(indent, g_complete, g_cycle, f, p, q, r, x, y, z))

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
    
    def _dxf_arc_append(self, clockwise, end_x, end_y, radius, tracing=-1000000):
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
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._dxf_arc_append@(cw={1}, end_x={2:i}, end_y={3:i}, radius={4:i})".
	      format(indent, clockwise, end_x, end_y, radius))
	    trace_detail = 3

	# Compute some constants:
	degrees0 = Angle()
	degrees90 = Angle(deg=90.0)
	degrees180 = Angle(deg=180.0)
	degrees360 = Angle(deg=360.0)
    
	# Grab the start position ({x1}, {y1}):
	start_x = code._x_value()
	start_y = code._y_value()
	if trace_detail >= 1:
	    print("{0}x_before={1:i} {2:i}".format(indent, start_x, start_y))

	# What we are going to do here is fit a circle to S=(*start_x*,*start_y) and
	# E=(*end_x*,*end_y*).  We need to find the center of the circle
	# C=(*center_x*, *center_y*).  Note that only circle arc is shown below
        # when in fact there are two possible solutions.  Other information is ussed
	# to figure out which of the two circle centers to actually use:
        # The crude ASCII art below may be useful:
        #
	#               .....
	#           ..         ..
	#        S-----+--M--------E
	#         \    |  |       /
	#          \   +--+      /
	#           \     |     /
	#            \    |    /
	#            r\   |   /r
	#              \  |  /
	#               \ | /
	#                \|/
	#                 C
	# 
	#
        # What we know is that the radius r = |SC| = |EC|.
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
	#        |MC|^2 = r^2C - |MS|^2                            (2)
	#
	#        |MC| = sqrt(r^2 - |MS|^2)                         (3)
	#

	# Let's get *s* and *e* defined:
	zero = L()
	s = P(start_x, start_y, zero)
	e = P(end_x, end_y, zero)
	
	# Find the midpoint *m* and perform any *tracing*:
	m = (s + e) / 2

	# Now compute *ms_length* and *ms_angle*
	ms = m - s
	ms_length = ms.length()
	ms_angle = ms.xy_angle()
	if trace_detail >= 2:
	    print("{0}m={1:i} m-s={2:i} |ms|={3:i} <ms={4:d}".
	      format(indent, m, ms, ms_length, ms_angle))

	# Now compute *mc_length*:
	ms_mm = ms_length.millimeters()
	radius_mm = radius.millimeters()
	# Unclear why the `float` cast is needed.  Without the cast, the code fails in Python2:
	radius2_minus_ms2 = abs(radius_mm * radius_mm - ms_mm * ms_mm)
	if radius2_minus_ms2 < 0.0:
	    print("Code._dxf_arc_append() has a sign error")
	    radius2_minus_ms2 = -radius2_minus_ms2
	mc_mm = math.sqrt(radius2_minus_ms2)
	mc_length = L(mm=mc_mm)

	# Now compute *cm_angle* being 90 degrees of of *ms_angle*:
	if clockwise:
	    cm_angle = ms_angle - degrees90
	else:
	    cm_angle = ms_angle + degrees90
	mc_angle = cm_angle.normalize()

	# Now compute center *c*:
	center_x = m.x + mc_length.cosine(mc_angle)
	center_y = m.y + mc_length.sine(mc_angle)
	c = P(center_x, center_y, zero)
	if trace_detail >= 3:
	    print("{0}c={1:i} |mc|={2:i} <mc={3:d}".format(indent, c, mc_length, mc_angle))

	# Now we need the *start_angle* and *end_angle* for DXF entity:
	start_angle = (s-c).xy_angle()
	end_angle = (e-c).xy_angle()
	if trace_detail >= 2:
	    print("{0}c={1:i} mc_angle={2:d} start_angle={3:d} end_angle={4:d}".
	      format(indent, c, mc_angle, start_angle, end_angle))

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
    
	# Wrap up any *tracing*:
	if tracing >= 0:
	    if trace_detail >= 1:
		print("{0}x_after={1:i} y_after{2:i}".format(indent, end_x, end_y))
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

    def _finish(self, tracing = -1000000):
	""" *Code*: Finish off the current block of G code (i.e. *self*.)  The estimated time
	    to perform the CNC code is returned.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._finish()".format(indent))

	# Return to *tool_change_point*:
	tool_change_point = code._tool_change_point
	code._command_begin()
	code._unsigned("G8", 49)	# Disable tool offset
	code._mode_motion(0, code._vrml_motion_color_tool_change)
	code._length("X", tool_change_point.x)
	code._length("Y", tool_change_point.y)
	code._length("Z", tool_change_point.z)
	code._comment("Return to tool change point")
	code._command_end()

	# Write out the final commands to wrap up the subroutine and end the 
	code_stream = code._code_stream
	assert isinstance(code_stream, file)
	assert code._tool_program_number % 100 != 0
	code_stream.write("M5 ( Turn spindle off)\n")
	code_stream.write("M9 ( Turn coolant off)\n")
	code_stream.write("O{0} endsub\n".format(code._tool_program_number))
	code_stream.write("O{0} call ( For stand-alone use of NGC file we call {0} here )\n".
	  format(code._tool_program_number))
	code_stream.write("M2 ( End of Program )\n")
	code_stream.write("( Path bounding box: {0:i} )\n".format(code._bounding_box))
	code_stream.write("( Path bounding box volume: {0:i} )\n".
	  format(code._bounding_box.volume_get()))
	code._tool_program_number = -1

	# Close *code_stream*:
	code_stream.close()
	code._code_stream = None

	# Reset the *top_surface_safe_z* and *xy_rapid_safe_z*:
	code._top_surface_safe_z = None
	code._xy_rapid_safe_z = None

	# Write VRML to *tool_wrl_lines* and *mount_wrl_lines*.
	tool_wrl_lines = []
	mount_wrl_lines = code._mount_wrl_lines

	# Flush any remaining points into *vrml_lines*:
	code._vrml_points_flush("")

	# Grab the path time:
	time = code._time

	# Reset *code*:
	code._reset()

	# Now reset all the VRML values:
	code._vrml_reset()

	# Perform any requested *tracing* and return *time*:
	if tracing >= 0:
	     print("{0}<=Code._finish()=>{1:m}".format(indent, time))
	return time

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
	elif field_name == "R":
	    if code._r != value:
		code._r = value
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

    def _mode_motion(self, g1_value, motion_color):
	""" *Code*: This routine will issue a G*g1_value* motion command to the *Code* object
	    (i.e. *self*).  *motion_color* specifies what color to draw into the tool
	    path .wrl file.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify arguement types:
	assert isinstance(g1_value, int) and g1_value in code._g1_table
	assert isinstance(motion_color, str)

	# Add a G1 field to the command:
	code._unsigned("G1", g1_value)
	code._vrml_motion_color = motion_color

    def _reset(self):
	""" *Code*: Reset the contents of the *Code* object (i.e. *self*) """

	# Use *code* instead of *self:
	code = self

	zero = L()
	large = L(inch=123456789.0)
	huge = 0x7fffffff
	big = L(inch=123456789)

	# Reset the *code* object:
	code._begin = True
	code._bounding_box = Bounding_Box()
	code._time = Time()
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
	code._h = -1
	code._i = zero
	code._j = zero
	code._m1 = huge
	code._m2 = huge
	code._m3 = huge
	code._m4 = huge
	code._m5 = huge
	code._p = Time(sec=-1.0)
	code._q = big
	code._r = big
	code._r0 = big
	code._r1 = big
	code._s = Hertz(rps=123456789.0)
	code._x = big
	code._y = big
	code._z = big
	code._z1 = zero

    def _start(self, part, tool, tool_program_number, spindle_speed, mount, vrml, tracing=-1000000):
	""" *Code*: Start writing out the G-code for *tool*...
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(tool, Tool) # and tool._name_get() != "None"
	assert isinstance(tool_program_number, int)
	assert isinstance(spindle_speed, Hertz)
	assert isinstance(mount, Mount)
	assert isinstance(vrml, VRML_Lines)
	assert isinstance(tracing, int)

	# Stuff some values into *code*:
	code._vrml = vrml
	code._tool_program_number = tool_program_number

	# Make darn sure that we are not writing to the top level file:
	assert tool_program_number % 100 != 0

	# Reset *code*:
	code._reset()

	# Do any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._start('{1}', '{2}', {3}, {4:rpm}, '{5}')".
	      format(indent, part._name_get(), tool._name_get(),
	        tool_program_number, spindle_speed, mount._name_get()))
	    trace_detail = 2

	# Grab some values from *part* and *tool*:
	part_name = part._name_get()
	tool_name = tool._name_get()
	tool_number = tool._number_get()

	# Make sure that we closed off any previous *code_stream*:
	assert code._code_stream == None

	# Stuff some values into *code*:
	code._mount = mount

	# Open new *code_stream*:
	ezcad = part._ezcad_get()
	ngc_directory = ezcad._ngc_directory_get()
	# This is pretty ugly.  For a mount "tool", there is no need to generate a `.ngc` file.
	# This is a bit of kludge.  We for a mount "tool", we just write everything to `/dev/null`.
	# Otherwise, we open the *tool_program_number*`.ngc`.
	if tool._name_get() == "None":
	    code_file_name = "/dev/null"
	    code_stream = open(code_file_name, "w")
	else:
	    code_file_name = "{0}.ngc".format(tool_program_number)
	    code_stream = ngc_directory._write_open(code_file_name)
	code._code_stream = code_stream
	if trace_detail >= 2:
	     print("{0}code_file_name='{1}' opened".format(indent, code_file_name))

	# Grap *vice* and *tool_change_point* and save into *code*:
	shop = ezcad._shop_get()
	vice = shop._vice_get()
	tool_change_point = vice._tool_change_point_get()
	code._tool_change_point = tool_change_point

	# Clear out *vrml_points* and intialize the location to be at *tool_change_point*:
	vrml_points = code._vrml_points
	del vrml_points[:]
	vrml_points.append(tool_change_point)
	if trace_detail >= 2:
	    print("{0}SA code._vrml_points={1}".
	      format(indent, [ "{0:i}".format(vrml_point) for vrml_point in vrml_points ]))

	# Output a descriptive header comment:
	code_stream.write("( Part {0}: Tool {1} Program: {2} )\n".
	  format(part_name, tool_name, tool_program_number))

	# Declare the the subroutine number:
	code_stream.write("O{0} sub\n".format(tool_program_number))

	# Set absolute mode for X/Y/Z (G90) and absolute mode for I/J/K (G90.1).
	# Set units to inches (G20):
	code_stream.write("G90 G90.1 G20\n")

	if trace_detail >= 2:
	    print("{0}SB code._vrml_points={1}".
	      format(indent, [ "{0:i}".format(vrml_point) for vrml_point in code._vrml_points ]))

	# Start at *tool_change_point*:
	code._command_begin()
	code._unsigned("G8", 49) # Disable tool offset
	code._mode_motion(0, code._vrml_motion_color_rapid)
	code._length("X", tool_change_point.x)
	code._length("Y", tool_change_point.y)
	code._length("Z", tool_change_point.z)
	code._comment("Start at tool change point with tool offsets disabled")
	code._command_end(vrml_suppress = True, tracing = tracing + 1)
	if trace_detail >= 2:
	    print("{0}SC code._vrml_points={1}".
	      format(indent, [ "{0:i}".format(vrml_point) for vrml_point in code._vrml_points ]))

	# Output the tool change, get the coolant on, and spindle spun up:
	#assert tool_number >= 0, \
	#  "Part:'{0}' Tool:'{1}' num={2}".format(part_name, tool_name, tool_number)

	code_stream.write("M6 T{0} ( Insert tool {0}: {1} )\n".format(tool_number, tool_name))

	# Output `G43 Htool_number` to enable tool offset:
	code._command_begin()
	code._unsigned("G8", 43)
	code._unsigned("H", tool_number)
	code._comment("Enable tool offset for tool {0}".format(tool_number))
	code._command_end()
	if trace_detail >= 2:
	    print("{0}SE code._vrml_points={1}".
	      format(indent, [ "{0:i}".format(vrml_point) for vrml_point in code._vrml_points ]))

	# Turn spindle on or off:
	spindle_is_on = spindle_speed > Hertz()
	if spindle_is_on:
	    code_stream.write("( Get the spindle up to speed )\n")
	    code_stream.write("M3 S{0:rpm} ( Turn spindle on )\n".format(spindle_speed))
	else:
	    code_stream.write("( Spindle is off for this tool )\n")
	    code_stream.write("M5 ( Turn spindle off )\n")

	# Decide whether to turn coolant on or off:
	material = part._material_get()
	if material._needs_coolant() and spindle_is_on:
	    code_stream.write("M8 (Coolant on)\n")
	else:
	    code_stream.write("M9 (Coolant off)\n")

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Code._start('{1}', '{2}', {3}, {4:rpm}, '{5}')".
	      format(indent, part._name_get(), tool._name_get(),
	        tool_program_number, spindle_speed, mount._name_get()))

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

    def _simple_pocket_helper(self, cnc_corner_bsw, cnc_corner_tne, corner_radius, z,
      tool_radius, offset, spindle_speed, feed_speed, rapid_move, rotate, tracing=-1000000):
	""" *Code*: Using the *Code* object (i.e. *self*), generate one rectangular or
	    rounded rectangular path for the currently selected tool:
	    * *cnc_corner_bsw* and *cnc_corner_tne* specify diagonally opposite corners
               of the pocket.
	    * *corner_radius* specifies the inside radius of the pocket corner.
	    * *z* specifies the Z altitude of the tool tip.
            * *tool_radius* specifies the radius of the end mill tool.
	    * *offset* specifies the distance inward from the outermost rectangular path.
	    * *spindle_speed* and and *feed_speed* specify the speed and feed for the tool.
	    * *rapid_move* specifies whether or not to move to the start position with rapid
	      tool movement (i.e. G0) or linear tool movement (i.e. G1) commands.
	    * *rotate* specifies how much to rotate the pocket around its center point.
	"""
    
	# Verify argument types:
	zero = L()
	assert isinstance(cnc_corner_bsw, P)
	assert isinstance(cnc_corner_tne, P)
	assert isinstance(corner_radius, L) and corner_radius >= zero
	assert isinstance(z, L)
	assert isinstance(tool_radius, L)
	assert isinstance(offset, L) and offset >= zero
	assert isinstance(spindle_speed, Hertz)
	assert isinstance(feed_speed, Speed)
	assert isinstance(rapid_move, bool)
	assert isinstance(rotate, Angle)

	# Perform an requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print(("{0}=>Code._simple_pocket_helper(*, bsw={1:i}, tne={2:i}, cr={3:i}" +
	      "z={4:i} tr={5:i}, o={6:i}, ss={7:rpm}, fs={8:i}, rm={9})").format(indent,
	      cnc_corner_bsw, cnc_corner_tne, corner_radius, z, tool_radius, offset,
	      spindle_speed, feed_speed, rapid_move))
	    trace_detail = 2

	# Stuff some debug information into a *line_coment* for *code* (i.e. *self*):
	code = self
	code._line_comment(
	  "bsw={0:i} tne={1:i} cr={2:i} off={3:i} z={4:i}".
	  format(cnc_corner_bsw, cnc_corner_tne, corner_radius, offset, z))
    
	# Yank the X/Y coordinates out of two corners:
	x1 = cnc_corner_bsw.x
	x2 = cnc_corner_tne.x
	y1 = cnc_corner_bsw.y
	y2 = cnc_corner_tne.y
	assert x1 < x2
	assert y1 < y2
	if trace_detail >= 2:
	    print("{0}x1={1:i} y1={2:i} x2={3:i} y2={4:i}".format(indent, x1, y1, x2, y2))
    
	# Compute the arc center locations:
	rx1 = x1 + corner_radius
	rx2 = x2 - corner_radius
	ry1 = y1 + corner_radius
	ry2 = y2 - corner_radius
	assert rx1 < rx2
	assert ry1 < ry2
	if trace_detail >= 2:
	    print("{0}rx1={1:i} ry1={2:i} rx2={3:i} ry2={4:i}".format(indent, rx1, ry1, rx2, ry2))

	# Compute the offset X/Y coordinates:
	ox1 = x1 + offset
	oy1 = y1 + offset
	ox2 = x2 - offset
	oy2 = y2 - offset
	code._line_comment("offset={0:i} ox1={1:i} oy1={2:i} ox2={3:i} oy2={4:i}".
	  format(offset, ox1, oy1, ox2, oy2))
	if trace_detail >= 2:
	    print("{0}offset={1:i} ox1={2:i} oy1={3:i} ox2={4:i} oy2={5:i}".
	      format(indent, offset, ox1, oy1, ox2, oy2))

	# Compute the middle of the pocket:
	middle_x = (x1 + x2) / 2
	middle_y = (y1 + y2) / 2
	middle = P(middle_x, middle_y, zero)

	# We need to compute the following points for the pocket:
        #
        #             NW                NE    
        #              *----------------*    
	#            /                    \
	#          /                        \
        #       WN*    +                +    *EN
	#         |   CNW              CNE   |
	#         |                          |
	#         |   CSW              CSE   |
        #       WS*    +                +    *ES
	#          \                        /
	#            \                    /
        #              *----------------*
        #             SW                SE
        # where:
        #
        # * Cxx stands for the center point for the radius of the 4 pocket corners
        # * xx  stands for one of the end points of a flat surface.
	sw  = P(rx1, oy1, zero)
	se  = P(rx2, oy1,  zero)
	cse = P(rx2, ry1, zero)
	es  = P(ox2, ry1, zero)
	en  = P(ox2, ry2, zero)
	cne = P(rx2, ry2, zero)
	ne  = P(rx2, oy2,  zero)
	nw  = P(rx1, oy2,  zero)
	cnw = P(rx1, ry2, zero)
	wn  = P(ox1, ry2, zero)
	ws  = P(ox1, ry1, zero)
	csw = P(rx1, ry1, zero)

	# Determine the starting location for this path:
	start_x = ox1
	start_y = oy1
	if offset < corner_radius:
	    # We have rounded corner path:
	    start_x = rx1
	start = P(start_x, start_y, zero)

	# Compute rotated values of pocket:
	rsw    = sw.xy_rotate(   middle, rotate)
	rse    = se.xy_rotate(   middle, rotate)
	rcse   = cse.xy_rotate(  middle, rotate)
	res    = es.xy_rotate(   middle, rotate)
	ren    = en.xy_rotate(   middle, rotate)
	rcne   = cne.xy_rotate(  middle, rotate)
	rne    = ne.xy_rotate(   middle, rotate)
	rnw    = nw.xy_rotate(   middle, rotate)
	rcnw   = cnw.xy_rotate(  middle, rotate)
	rwn    = wn.xy_rotate(   middle, rotate)
	rws    = ws.xy_rotate(   middle, rotate)
	rcsw   = csw.xy_rotate(  middle, rotate)
	rstart = start.xy_rotate(middle, rotate)

	# Move to *rstart* as specified by *rapid_move* argument:
	if rapid_move:
	    code._xy_rapid(rstart.x, rstart.y)
	else:
	    code._xy_feed(feed_speed, spindle_speed, rstart.x, rstart.y)
    
	# Make sure we are at the depth *z*:
	code._z_feed(feed_speed/2, spindle_speed, z, "simple_pocket_helper")
    
	# Mill out either a square or rounded corners
	if offset < corner_radius:
	    # Mill out a rectangle with rounded corners in a
	    # counter clockwise direction to force a climb cut:
    
	    # *r* is the radius of the arc portion of the pocket path:
	    r = corner_radius - offset
    
	    # Bottom horizontal line from *rsw* to *rse*:
	    code._xy_feed(feed_speed, spindle_speed, rse.x, rse.y, tracing=tracing+1)
    
	    # Lower right arc *rse* to *res*:
	    code._xy_ccw_feed(feed_speed,
	      spindle_speed, r, res.x, res.y, rx=rcse.x, ry=rcse.y, tracing=tracing+1)
    
	    # Right vertical line *res* to *ren*
	    code._xy_feed(feed_speed, spindle_speed, ren.x, ren.y)
    
	    # Upper right arc *ren* to *rne* with arc center at *rcne*:
	    code._xy_ccw_feed(feed_speed, spindle_speed, r, rne.x, rne.y, rx=rcne.x, ry=rcne.y)

	    # Top horizontal line *rne* to *rnw*:
	    code._xy_feed(feed_speed, spindle_speed, rnw.x, rnw.y)
    
	    # Upper left arc *rnw* to *rwn* with arc center at *rcnw*:
	    code._xy_ccw_feed(feed_speed, spindle_speed, r, rwn.x, rwn.y, rx=rcnw.x, ry=rcnw.y)
    
	    # Left vertical line *rwn* to *rws*:
	    code._xy_feed(feed_speed, spindle_speed, rws.x, rws.y)
    
	    # Lower left arc *rws* to *rsw* with arc center at *rcsw*:
	    code._xy_ccw_feed(feed_speed, spindle_speed, r, rsw.x, rsw.y, rx=rcsw.x, ry=rcsw.y)
	else:
	    # Mill out a rectangle with "square" corners in a counter
	    # clockwise direction to force a climb cut:
    
	    sq_sw = P(ox1, oy1, zero)
	    sq_se = P(ox2, oy1, zero)
	    sq_ne = P(ox2, oy2, zero)
	    sq_nw = P(ox1, oy2, zero)
	    rsq_sw = sq_sw.xy_rotate(middle, rotate)
	    rsq_se = sq_se.xy_rotate(middle, rotate)
	    rsq_ne = sq_ne.xy_rotate(middle, rotate)
	    rsq_nw = sq_nw.xy_rotate(middle, rotate)

	    # Bottom horizontal line from *rsq_sw to* to *rsq_se*
	    code._xy_feed(feed_speed, spindle_speed, rsq_se.x, rsq_se.y, tracing=tracing+1)
    
	    # Right vertical line from *rsq_se* to *rsq_ne*
	    code._xy_feed(feed_speed, spindle_speed, rsq_ne.x, rsq_ne.y)
    
	    # Top horizontal line from *rsq_ne* to *rsq_nw*:
	    code._xy_feed(feed_speed, spindle_speed, rsq_nw.x, rsq_nw.y)
    
	    # Left vertical line from *rsq_nw* to *rsq_sw*:
	    code._xy_feed(feed_speed, spindle_speed, rsq_sw.x, rsq_sw.y)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print(("{0}<=Code._simple_pocket_helper(*, bsw={1:i}, tne={2:i}, cr={3:i}" +
	      "z={4:i} tr={5:i}, o={6:i}, ss={7:rpm}, fs={8:i}, rm={9})").format(indent,
	      cnc_corner_bsw, cnc_corner_tne, corner_radius, z, tool_radius, offset,
	      spindle_speed, feed_speed, rapid_move))

    def _speed(self, field_name, value):
	""" *Code*: Set the speed for *field_name* to *value* in the the current command of the
	    *Code* object (i.e. *self*).
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(field_name, str)

	square_bracket = False
	changed = False
	if field_name == "F":
	    assert isinstance(value, Speed)
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
		chunk = "[{0:I}]".format(value)
	    else:
		chunk = "{0}{1:I}".format(field_name[0], value)
	    code._command_chunks.append(chunk)

    def _time_get(self):
	""" *Code*: Return the time required for the perform the tool path specified by the
	    *Code* object (i.e. *self*:
	"""

	return self._time

    def _tool_set(self, tool):
        """ *Code*: Set the current tool for the *Code* object (i.e. *self*) to *tool*:
	"""

	# Verify argument types:
	assert isinstance(tool, Tool)

	# Load *tool* into *code* (i.e. *self*):
	code = self
	code._tool = tool
	code._is_laser = tool._is_laser_get()

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
	mount_ngc_file_name = os.path.join(ngc_directory, "{0}.ngc".format(program_number))
	part_ngc_stream = open(mount_ngc_file_name, "w")
	assert part_ngc_stream != None, "Unable to open {0}".format(mount_ngc_file_name)

    def _unsigned(self, field_name, value):
	""" *Code*: This routine will format {value} for {field_name} to {code}.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(field_name, str)
	assert isinstance(value, int) # and value >= 0, "value={0}".format(value)

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
	elif size == 3:
	    if field_name == "G10":
		matched = True
		previous_value = code._g10
		if previous_value != value:
		    code._g10 = value
	    elif field_name == "G11":
		# Always output G4 when requested:

		matched = True
		previous_value = 0x12345

	assert matched, "Unrecognized field name {0}".format(field_name)
	if previous_value != value or field_name == "G1" or field_name == "G8":
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

    def _vrml_arc_draw(self,
      ax, ay, bx, by, radius, z, clockwise, radius_x=None, radius_y=None, tracing=-1000000):
	""" *Code*: Draw an arc from (*ax*, *ay*, *z) to (*bx*, *by*, *bz*) with an
	    arc radius of *radius*.  The arc is drawn in a clockwise direction if
	    *clockwise* is *True* and a counter-clockwise direction otherwise.
	"""
    
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

	# Perform an requested_*tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._vrml_arc_draw({1:i}, {2:i}, {3:i}, {4:i}, {5:i}, {6})".
	      format(indent, ax, ay, bx, by, z, clockwise))
	    trace_detail = 3

	# Define some *Angle* constants:
	degrees0 = Angle(deg=0.0)
	degrees15 = Angle(deg=15.0)
	degrees90 = Angle(deg=90.0)
	degrees180 = Angle(deg=180.0)
	degrees360 = Angle(deg=360.0)
    
	# Determine if we *have_radius*:
	have_radius = isinstance(radius_x, L) and isinstance(radius_y, L)
	if have_radius:
	    if trace_detail >= 2:
		print("{0}rx={1:i} ry={2:i}".format(indent, radius_x, radius_y))
		ar_dx = ax - radius_x
		ar_dy = ay - radius_y
		ar_radius = ar_dx.distance(ar_dy)
		print("{0}ar_dx={1:i} ar_dy={2:i} ar_radius={3:i}".
		 format(indent, ar_dx, ar_dy, ar_radius))

		br_dx = bx - radius_x
		br_dy = by - radius_y
		br_radius = br_dx.distance(br_dy)
		print("{0}br_dx={1:i} br_dy={2:i} br_radius={3:i}".
		  format(indent, br_dx, br_dy, br_radius))

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
	if trace_detail >= 2:
	    print("{0}cx={1:i} cy={2:i} cy={3:i}".format(indent, cx, cy, cz))
	#if trace_detail >= 2:
	#	assert isinstance(rx, L)
	#	assert isinstance(cx, L)
	#	assert isinstance(ry, L)
	#	assert isinstance(cy, L)
	#	assert isinstance(cz, L)
	#	radius3 = (radius_x - cx).distance(radius_y - cy)
	#	print("{0}C=({1:i},{2:i},{3:i}) radius_3={4:i}".format(indent, cx, cy, cz, radius3))

	# Now compute |AC| which is the length the the AC segment:
	ac_length = (cx - ax).distance(cy - ay)
	if trace_detail >= 2:
	    print("{0}ac_length={1:i}".format(indent, ac_length))

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
	if trace_detail >= 3:
	    print("{0}radius_mm={1} ac_mm={2}".format(indent, radius_mm, ac_mm))
	#assert radius_squared_mm >= ac_squared_mm, "difference={0}".format(radius_mm - ac_mm)
	cr_mm = math.sqrt(abs(radius_squared_mm - ac_squared_mm))
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
	if trace_detail >= 2:
	    print("{0}ac_dx={1:i} ac_dy={2:i} acx_angle={3:d}".
	      format(indent, ac_dx, ac_dy, acx_angle))

	# Step 2: Compute angle <RCX, which is +/- 90 degrees from <ACX:
	rcx_angle = acx_angle
	if clockwise:
	    rcx_angle += degrees90
	else:
	    rcx_angle -= degrees90
	if trace_detail >= 2:
	    print("{0}acx_angle={1:d} rcx_angle={2:d}".format(indent, acx_angle, rcx_angle))

	# Step 3: Normalize rcx_angle to between -180 degrees and +180 degrees:
	rcx_angle = rcx_angle.normalize()
	if trace_detail >= 2:
	    print("{0}normalized rcx_angle={1:d}".format(indent, rcx_angle))

	# Step 4: Compute R=(*rx*,*ry*,*rz*) using trigonmetry:
	rx = cx + cr_length.cosine(rcx_angle)
	ry = cy + cr_length.sine(rcx_angle)
	rz = az
	if trace_detail >= 2:
	    print("{0}computed rx={1:i} ry={2:i} rz={3:i}".format(indent, rx, ry, rz))
	    if have_radius:
		print("{0}radius_x={1:i} radius_y={2:i}".format(indent, radius_x, radius_y))

	# Now we compute the angle <ARX which is a the angle from R to A relative to the X axis:
	arx_angle = (ay - ry).arc_tangent2(ax - rx)

	# Similarly, we compute the <BRX which is the angle from R to B relative to the Y axis:
	brx_angle = (by - ry).arc_tangent2(bx - rx)

	if trace_detail >= 2:
	    print("{0}arx_angle={1:d} brx_angle={2:d}".format(indent, arx_angle, brx_angle))

	# Now we can compute angle <ARB which is the angle swept from A to B from R.
	# Note that <ARB can be positive or negative depending up the direction of travel:
	arb_angle = brx_angle - arx_angle
	if trace_detail >= 2:
	    print("{0}arb_angle={1:d}".format(indent, arb_angle))

	# Normalize <ARB to be consistent with *clockwise*:
	if clockwise and arb_angle < degrees0:
	    # Clockwise requires *arb_angle* to be positive:
	    arb_angle += degrees360
	elif not clockwise and arb_angle > degrees0:
	    # Counter-clockwise requires *arb_angle* to be negative:
	    arb_angle -= degrees360
	if trace_detail >= 2:
	    print("{0}normalize arb_angle={1:d}".format(indent, arb_angle))

	# We want to draw the arc as a sequence of line segments from A to B, where each segment
	# covers about 15 degrees of the arc.  So first we figure compute the *steps* count:
	steps = 1 + abs(int(arb_angle / degrees15))
	assert steps >= 1
	if trace_detail >= 2:
	    print("{0}steps={1}".format(indent, steps))

	# Now we compute the *step_angle* which is the arc angle coverted by each segment:
	step_angle = arb_angle / float(steps)

	# Lay down A=(*ax*,*ay*,*az*):
	vrml_points = code._vrml_points
	vrml_points.append(P(ax, ay, z))

	# Lay down the intermediate points along the arc:
	for step in range(1, steps - 1):
	    segment_angle = arx_angle + step_angle * step
	    segment_x = rx + radius.cosine(segment_angle)
	    segment_y = ry +   radius.sine(segment_angle)
	    segment_z = az
	    segment = P(segment_x, segment_y, segment_z)
	    vrml_points.append(segment)
	    if trace_detail >= 3:
		print("{0}[{1}]: segment_angle={2:d} segment".
		  format(indent, step, segment_angle, segment))

	# Lay down B=(*bx*,*by*,*bz):
	vrml_points.append(P(bx, by, z))

	# Wrap up any requested_*tracing*:
	if tracing >= 0:
	    print("{0}<=Code._vrml_arc_draw({1:i}, {2:i}, {3:i}, {4:i}, {5:i}, {6})".
	      format(indent, ax, ay, bx, by, z, clockwise))

    def _vrml_point(self, x, y, z, from_routine):
	""" *Code*: Return the VRML point index for point (*x*, *y*, *z*) using the
	    *Code* object (i.e. *self*).
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(z, L)
	assert isinstance(from_routine, str)

	code._vrml_lines_points.append(P(x, y, z))
	point = (x.millimeters(), y.millimeters(), z.millimeters())

    def _vrml_points_append(self, color, point, tracing=-1000000):
	""" *Code*: Append *point* to the VRML points list of the *Code* object (i.e. *self*).
	    The drawn line will be drawn in *color*.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(color, str)
	assert isinstance(point, P)
	assert isinstance(tracing, int)
	assert -1000000.00 <= point.x.millimeters() <= 1000000
	assert -1000000.00 <= point.y.millimeters() <= 1000000
	assert -1000000.00 <= point.z.millimeters() <= 1000000

	# Grab *vrml_points* from *code*:
	vrml_points = code._vrml_points

	# Expand *code* bounding box to include *point*:
	code._bounding_box.point_expand(point)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = tracing * ' '
	    print("{0}=>Code._vrml_points_append(*, {1}, {2:i})".format(indent, color, point))
	    trace_detail = 2

	# A color change forces a flush:
	if trace_detail >= 2:
	    print("{0}before vrm_points={1}".
	      format(indent, [ "{0:i}".format(vrml_point) for vrml_point in vrml_points ]))
	vrml_current_color = code._vrml_current_color
	if color != vrml_current_color:
	    code._vrml_points_flush(color, tracing = tracing + 1)
	if len(vrml_points) == 0 or vrml_points[-1] != point:
	    vrml_points.append(point)
	if trace_detail >= 2:
	    print("{0}after vrm_points={1}".
	      format(indent, [ "{0:i}".format(vrml_point) for vrml_point in vrml_points ]))

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Code._vrml_points_append(*, {1}, {2:i})".format(indent, color, point))

    def _vrml_points_flush(self, new_color, tracing=-1000000):
	""" *Code*: Flush out any vrml lines for the *Code* object (i.e. *self*).
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(new_color, str)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._vrml_points_flush(*)".format(indent))
	    trace_detail = 2

	vrml_current_color = code._vrml_current_color
	vrml_points = code._vrml_points
	if trace_detail >= 2:
	    print("{0}before vrml_current_color='{1}'".format(indent, vrml_current_color))
	    print("{0}before vrml_points={1}".
	      format(indent, [ "{0:i}".format(vrml_point) for vrml_point in vrml_points ]))
	if new_color != vrml_current_color:
	    # The path color has changed, so we need to flush out *vrml_lines_points* using the
	    # current color:
	    if vrml_current_color == "":
	         vrml_current_color = "yellow"
	    vrml = code._vrml
	    assert isinstance(vrml, VRML)
	    vrml._poly_line(vrml_current_color, vrml_points)

	    # We need to clear out *vrml_lines_points* and start it with the *last_point* so
 	    # there will be no gaps in the path:
	    if len(vrml_points) > 0:
		last_point = vrml_points[-1]
		del vrml_points[:]
		vrml_points.append(last_point)
	    else:
		zero = L()
		vrml_points.append(P(zero, zero, zero))
	    code._vrml_current_color = new_color
	if trace_detail >= 2:
	    print("{0}after vrml_current_color='{1}'".format(indent, vrml_current_color))
	    print("{0}after vrml_points={1}".
	      format(indent, [ "{0:i}".format(vrml_point) for vrml_point in vrml_points ]))

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Code._vrml_points_flush(*)".format(indent))

    def _vrml_reset(self):
	""" *Code*: Reset the VRML sub-system of the *Code* object (i.e. *self*). """

	# Use *code* instead of *self*:
	code = self

	zero = L()
	code._mount_wrl_lines = None
	code._wool_wrl_lines = None
	code._vrml = None
	code._vrml_file = None

	code._vrml_points = [ ]
	code._vrml_current_color = ""

    def _x_value(self):
	""" *Code*: Return the current value of X offset by the vice X for the *Code* object
	    (i.e. *self*).
	"""

	# Use *code* instead of *self*:
	code = self

	return code._x # + code._vice_x

    def _xy_cw_feed(self, f, s, r, x, y, rx=None, ry=None, tracing=-1000000):
	""" *Code*: Feed to location (*x*, *y*) with a radius *r* clockwise  circle with a
	    feedrate of *f* and spindle speed of *s* using the *Code* object (i.e. *self*):
	"""
    
	# Verify routine arguments:
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)
	assert isinstance(r, L)
	assert isinstance(x, L)
	assert isinstance(y, L)
	assert isinstance(rx, L) or rx == None
	assert isinstance(ry, L) or ry == None
	assert isinstance(tracing, int)
    
	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._xy_cw_feed(*, f={1:i}, s={2:rpm}, r={3:i} x={4:i} y={5:i}, ...)".
	      format(indent, f, s, r, x, y))
	    trace_detail = 1

	# Grab some values from the *Code* object (i.e. *self*):
	code = self
	x1 = code._x_value()
	y1 = code._y_value()
	z1 = code._z
	x2 = x
	y2 = y
	if trace_detail >= 1:
	    print("{0}x_before={1:i} y_before={2:i}".format(indent, x1, y1))

	code._vrml_arc_draw(x1, y1, x2, y2, r, z1, False, rx, ry)

	x_value = code._x_value()
	y_value = code._y_value()
	if x_value != x or y_value != y:
	    # Do the laser code first:
	    if code._is_laser:
	    	code._dxf_arc_append(True, x, y, r, tracing + 1)

	    # Now do the RS274 code and get F, R0, S, X, and Y updated:
	    #code._xy_rapid_safe_z_force(f, s)
	    code._r0 = L(inch=-100.0)	# Forget R
	    code._command_begin()
	    code._mode_motion(2, code._vrml_motion_color_cutting)
	    code._speed("F", f)
	    code._length("R0", r)
	    code._hertz("S", s)
	    code._length("X", x)
	    code._length("Y", y)
	    code._command_end()

	    # Estimate time it takes to perform the operation using a line segment:
	    zero = L()
	    p1 = P(x1, y1, zero)
	    p2 = P(x,  y,  zero)
	    distance = p1.distance(p2)
	    seconds = distance / f    # feed = distance / sec  ==>  sec = distance / feed
	    code._time += seconds
    
	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    if trace_detail >= 1:
		print("{0}x_after={1:i} y_after={2:i}".format(indent, x2, y2))
	    indent = ' ' * tracing
	    print("{0}<=Code._xy_cw_feed(*, f={1:i}, s={2:rpm}, r={3:i} x={4:i} y={5:i}, ...)".
	      format(indent, f, s, r, x, y))

    def _xy_ccw_feed(self, f, s, r, x, y, rx=None, ry=None, tracing=-1000000):
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
	assert isinstance(tracing, int)
    
	# Peform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._xy_ccw_feed(*, {1:i}, {2:rpm}, {3:i}, {4:i}, {5:i}) ".
	      format(indent, f, s, r, x, y))
	    trace_detail = 2

	# Grab some values from *code*:
	x1 = code._x_value()
	y1 = code._y_value()
	z1 = code._z
	x2 = x
	y2 = y
	code._vrml_arc_draw(x1, y1, x2, y2, r, z1, True, rx, ry, tracing = tracing + 1)

	if x1 != x or y1 != y:
	    if code._is_laser:
	    	code._dxf_arc_append(False, x, y, r, tracing + 1)
    
	    # Now do the RS274 code and get F, R0, S, X, and Y updated:
	    #code._xy_rapid_safe_z_force(f, s)
	    code._r0 = L(inch=-100.00)	# Forget R
	    code._command_begin()
	    code._mode_motion(3, code._vrml_motion_color_cutting)
	    code._speed("F", f)
	    code._length("R0", r)
	    code._hertz("S", s)
	    code._length("X", x)
	    code._length("Y", y)
	    code._command_end()

	    # Estimate time it takes to perform the operation using a line segment:
	    zero = L()
	    p1 = P(x1, y1, zero)
	    p2 = P(x,  y,  zero)
	    distance = p1.distance(p2)
	    seconds = distance / f    # feed = distance / sec  ==>  sec = distance / feed
	    code._time += seconds

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Code._xy_ccw_feed(*, {1:i}, {2:rpm}, {3:i}, {4:i}, {5:i}) ".
	      format(indent, f, s, r, x, y))

    def _xy_feed(self, f, s, x, y, tracing=-1000000):
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
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._xy_feed(f={1:i}, s={2:rpm}, x={3:i}, y={4:i})".
	      format(indent, f, s, x, y))
	    trace_detail = 1

	# Grab the current (*x_before*, *y_before*) values and possibly trace them:
	x_before = code._x_value()
	y_before = code._y_value()
	if trace_detail >= 1:
	    print("{0}x_before={1:i} y_before={2:i}".format(indent, x_before, y_before))

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

	    # Compute the time it takes to perform the operation:
	    zero = L()
	    p1 = P(x_before, y_before, zero)
	    p2 = P(x,        y,        zero)
	    distance = p1.distance(p2)
	    seconds = distance / f    # feed = distance / sec  ==>  sec = distance / feed
	    code._time += seconds

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    if trace_detail >= 1:
		print("{0}x_after={1:i} y_after={2:i}".format(indent, x, y))
	    print("{0}<=Code._xy_feed(f={1:i}, s={2:rpm}, x={3:i}, y={4:i})".
	      format(indent, f, s, x, y))

    def _xy_rapid(self, x, y):
	""" *Code*: Perform a rapid move to (X, Y) using the *Code* object (i.e. *self*). """

	# Use *code* instead of *self:
	code = self

	# Verify argument types:
	assert isinstance(x, L)
	assert isinstance(y, L)

	# Only perform the rapid if we are not already there:
	x1 = code._x_value()
	y1 = code._y_value()
	if x1 != x or y1 != y:
	    # Make sure we are at a safe Z height prior to performing the rapid:
	    mount = code._mount
	    assert code._z == mount._xy_rapid_safe_z_get()

	    # Output a X/Y rapid command:
	    code._command_begin()
	    code._mode_motion(0, code._vrml_motion_color_rapid)
	    code._length("X", x)
	    code._length("Y", y)
	    code._command_end()

	    # Estimate the time it takes to perform the operation:
	    zero = L()
	    p1 = P(x1, y1, zero)
	    p2 = P(x,  y,  zero)
	    distance = p1.distance(p2)
	    rapid_speed = code._rapid_speed
	    seconds = distance / rapid_speed    # feed = distance / sec  ==>  sec = distance / feed
	    code._time += seconds

    def _xy_rapid_safe_z_force(self, feed, speed, tracing=-100000):
	""" *Code*: Force the tool for the *Code* object (i.e. *self*) to be at the
	    safe height where rapids in X and Y are allowed.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(feed, Speed)
	assert isinstance(speed, Hertz)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._xy_rapid_safe_z_force(*, {1:i}, {2:rpm})".format(indent, feed, speed))

	# Make sure we are up at *xy_rapid_z_safe*:
	mount = code._mount
	xy_rapid_safe_z = mount._xy_rapid_safe_z_get()
	code._z_feed(feed, speed, xy_rapid_safe_z, "xy_rapid_safe_z_force", tracing + 1)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Code._xy_rapid_safe_z_force(*, {1:i}, {2:rpm})".format(indent, feed, speed))


    def _y_value(self):
	""" *Code*: Return the current value of Y offset by the vice Y for the *Code* object
	    (i.e. *self*).
	"""

	# Use *code* instead of *self*:
	code = self

	return code._y # + code._vice_y

    def _z_feed(self, f, s, z, from_routine, tracing=-1000000):
	""" *Code*: Feed to an altitude of *z* using the *Code* object (i.e. *self*)
	    at a feed of *f* and a speed *s*.
	"""

	# Use *code* instead of *self*:
	code = self

	# Verify argument types:
	assert isinstance(f, Speed)
	assert isinstance(s, Hertz)
	assert isinstance(z, L)
	assert isinstance(from_routine, str)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Code._z_feed(*, f={1:i}, s={2:rpm}, z={3:i}, '{4}')".
 	      format(indent, f, s, z, from_routine))
	    trace_detail = 1

	# This routine will feed to an altitude of {z} using {code}.

	# If *z_safe_pending* is set, but we are doing another Z feed before
	# it is cleared, we must not need to do a Z safe move operation after all.
	# Hence, we can just clear *z_safe_pending*:
	code._z_safe_pending = False

	mount = code._mount
	top_surface_safe_z = mount._top_surface_safe_z_get()
	xy_rapid_safe_z = mount._xy_rapid_safe_z_get()
	if trace_detail >= 1:
	    print("{0}top_surface_safe_z={1:i} xy_rapid_safe_z={2:i}".
	      format(indent, top_surface_safe_z, xy_rapid_safe_z))

	# Emit G0/G1 commands until we get to *z* (i.e. *z_target*):
	loop_count = 0
	if z > xy_rapid_safe_z:
	    z = xy_rapid_safe_z

	z_current = code._z
	while z_current != z:
	    z_target = z
	    if trace_detail >= 1:
		print(("{0}z_current={1:i} z_target={2:i} " +
		  "top_surface_safe_z={3:i} xy_rapid_safe_z={4:i}").format(indent,
		  z_current, z_target, top_surface_safe_z, xy_rapid_safe_z))

	    loop_count += 1
	    assert loop_count < 10

	    if z_current > top_surface_safe_z or \
	      z_current == top_surface_safe_z and z_target > top_surface_safe_z:
		# We are in the retraction portion of the workspace.  G0 rapids
		# are used to move around here:
		
		# Trim *z_target* to be between *top_surface_safe_z* and *xy_rapid_safe_z*:
		if z_target < top_surface_safe_z:
	             z_target = top_surface_safe_z
		elif z_target > xy_rapid_safe_z:
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
	    if z_current < top_surface_safe_z or \
	      z_current == top_surface_safe_z and z_target < top_surface_safe_z:
	 	# We are in the cutting portion of the workspace.  G1 motion is
	 	# used to move around here:

		# Make sure *z_target* does not get above *top_surface_safe_z*:
		if z_target >= top_surface_safe_z:
		    z_target = top_surface_safe_z
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

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Code._z_feed(*, f={1:i}, s={2:rpm}, z={3:i}, '{4}')".
 	      format(indent, f, s, z, from_routine))

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

class xxx_Place:
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

	# Tooling plate with #4 holes:
	no4_close_fit_diameter         = L(inch=0.1160) # #4-40 close = #32 drill
	no4_75_percent_thread_diameter = L(inch=0.0890) # #4-40 75%   = #43 drill
	no4_50_percent_thread_diameter = L(inch=0.0960) # #4-40 50%   = #41 drill
	no6_close_fit_diameter         = L(inch=0.1440) # #6-32 close = #27 drill
	no6_75_percent_thread_diameter = L(inch=0.1065)	# #6-32 75%   = #36 drill
	no6_50_percent_thread_diameter = L(inch=0.1160)	# #6-32 50%   = #32 drill
	tooling_plate = Tooling_Plate("Wayne's TP", L(inch=9.500), L(inch=4.500), L(inch=.500),
	  19, 9, L(inch=.500), no6_close_fit_diameter, no6_75_percent_thread_diameter,
	  no6_50_percent_thread_diameter, L(inch=0.500), L(inch=0.500))

	# Create *parallels*:
	parallels_length    = L(inch=6.000)
	parallels_thickness = L(inch="1/8")
	parallels = Parallels("Standard", parallels_length, parallels_thickness, [
	  L(inch=  "1/2"),
	  L(inch=  "5/8"),
	  L(inch=  "3/4"),
	  L(inch=  "7/8"),
	  L(inch= 1.000),
	  L(inch="1-1/8"),
	  L(inch="1-1/4"),
	  L(inch="1-3/8"),
	  L(inch="1-1/2") ] )
	tool_change_x = L(inch=-1.500)
	tool_change_y = L()
	tool_change_z = L(inch=8.000)
	tool_change_point = P(tool_change_x, tool_change_y, tool_change_z)

	# Create *vice*:
	vice_volume = P(L(inch=5.1), L(inch=5.0), L(inch=1.5))
	vice = Vice("5in_Vice", vice_volume, parallels, tool_change_point)

	# Initialize *shop*:
	shop._assemblies = []		# Viewable assemblies
	shop._base_name = ""		# Base name for current operations
	#dxf_base_names = [] 		# List of DXF base names
	shop._blocks_uid = 0		# UID counter for *Code_Block* objects
	# cache Cache			# Off/Nef3 file cache
	shop._changed = False		# Marker used by {update@Length}
	shop._cnc_generate = False	# {true} => do cnc code generation
	shop._code = Code()		# Code genertion object
	shop._drills = drills = []	# List of drills
	#dxf_table = {}			# Table of DXF base names
	shop._extra1 = P()		# Temporary extra bounding box point
	shop._extra2 = P()		# Temporary extra bounding box point
	shop._name = name		# Shop name
	shop._machines = []		# Machines in shop
	shop._parts = []		# Parts
	shop._program_base = 100	# Program base number
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
	#print("fpm600={0:F}".format(fpm600))
	fpm1200 = Speed(ft_per_min=1200)
	shop._surface_speeds_insert(aluminum, hss, fpm600, fpm1200)
	shop._surface_speeds_insert(aluminum, hss, fpm600, fpm1200)
	shop._surface_speeds_insert(plastic, hss, fpm600, fpm1200)
	shop._surface_speeds_insert(plastic, hss, fpm600, fpm1200)

	degrees45 = Angle(deg=45.0)
	degrees90 = Angle(deg=90.0)
	degrees118 = Angle(deg=118.0)
	degrees135 = Angle(deg=135.0)

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
	in1   =  L(inch=1.000)
	in4   =  L(inch=4.000)

	center_cut = True
	no_center_cut = not center_cut
	stub = Tool.DRILL_STYLE_STUB
	laser = True

	dowel_pin = shop._dowel_pin_append("3/8 Dowel Pin",
	  1, 9, hss, in3_8, L(inch=.900), in3_16)
	mill_drill_3_8 = shop._mill_drill_append("3/8 Mill Drill [.9in]",
	  2, hss, in3_8, 2, L(inch=.900), degrees90)
	drill_36 = shop._drill_append("#36 [#6-32 75% thread]",
	  3,     hss, L(inch=0.1065), 2, L(inch=1.500), "#36", degrees118, stub, center_cut, drills)
	dowel_36 = shop._dowel_pin_append("#36 Dowel Pin",
	  3, 53, hss, L(inch=0.1065), L(inch=1.500), zero)
	drill_27 = shop._drill_append("#27 [#6-32 close]",
	  4, hss, L(inch=0.1440), 2,  L(inch=1.875), "#27", degrees118, stub, center_cut, drills)
	dowel_27 = shop._dowel_pin_append("#27 Dowel Pin",
	  4, 54, hss, L(inch=0.1440), L(inch=1.875), zero)
	end_mill_3_8 = shop._end_mill_append("3/8 End Mill [5/8in]",
	  5, hss, in3_8, 4, in5_8, not laser)
	dowl_pin_end_mill_3_8 = shop._dowel_pin_append("3/8in Dowel Pin",
	  5, 55, hss, in3_8, in5_8, zero)
	end_mill_1_4 = shop._end_mill_append("1/4 End Mill",
	  6, hss, in1_4, 2, in1_2, not laser)
	double_angle = shop._double_angle_append("3/4 Double Angle",
	  7, hss, in3_4, 10, L(inch=0.875), degrees90, in1_4, in1_4)
	dove_tail = shop._dove_tail_append("3/8 Dove Tail",
	  8, hss, in3_8, 6, in1_4, in3_16, degrees45)
	# Note 9 is the alternate tool dowel pin for tool 1.
	end_mill_3_16 = shop._end_mill_append("3/16 End Mill",
	  10, hss, in3_16, 2, in1_2, not laser)
	drill_25 = shop._drill_append("#25 [#6-32 free]",
          11, hss, L(inch=0.1495), 2, L(inch=2.000), "#25", degrees118, stub, no_center_cut, drills)
	drill_9 = shop._drill_append("#9 [#10 close]",
	  12, hss, L(inch=0.1960), 2, L(inch=2.000), "#9", degrees118, stub, no_center_cut, drills)
	drill_43 = shop._drill_append("#43 [#4-40 75% thread]",
	  13, hss, L(inch=.0890), 2, L(inch=1.250), "#43", degrees118, stub, center_cut, drills)
	dowel_43 = shop._dowel_pin_append("#43 Dowel Pin",
	  13, 63, hss, L(inch=0.0890), L(inch=1.250), zero)
	drill_32 = shop._drill_append("#32 [#4-40 close]",
	  14, hss, L(inch=0.1160), 2,  L(inch=1.625), "#32", degrees118, stub, center_cut, drills)
	dowel_32 = shop._dowel_pin_append("#32 Dowel Pin",
	  14, 64, hss, L(inch=0.1160), L(inch=1.625), zero)
	drill_50 = shop._drill_append("#50 [#0-80 free]",
	  15, hss, L(inch=0.0700), 4, L(inch=0.750), "#50", degrees118, stub, no_center_cut, drills)
	end_mill_3_8_long = shop._end_mill_append("3/8 End Mill [1.6in]",
	  16, hss, in3_8, 4, L(inch=1.600), not laser)
	dowel_3_8_long = shop._dowel_pin_append("3/8 Dowel Pin [1.5in]",
	  16, 66, hss, in3_8, L(inch=1.500), zero)
	#end_mill_3_4 = shop._end_mill_append("3/4 End Mill",
	#  13, hss, in3_4, 2, in1_3_8)
	drill_30 = shop._drill_append("#30 [#4-40 free]",
	  17, hss, L(inch=0.1285), 2, L(inch=1.750), "#30", degrees118, stub, no_center_cut, drills)
	drill_1_8 = shop._drill_append("1/8",
	  18, hss, in1_8, 2, L(inch=1.750), "1/8", degrees118, stub, no_center_cut, drills)
	drill_1_16 = shop._drill_append("1/16",
	  19, hss, in1_16, 2, L(inch=1.750), "1/16", degrees118, stub, no_center_cut, drills)
	drill_3_32 = shop._drill_append("3/32",
	  20, hss, in3_32, 2, L(inch=1.750), "3/32", degrees118, stub, no_center_cut, drills)
	drill_42 = shop._drill_append("#42 [#4-44 75%]",
	  21, hss, L(inch=0.0935), 2, L(inch=1.750), "#42", degrees118, stub, no_center_cut, drills)
	drill_52 = shop._drill_append("#52 [#0-80 close]",
	  22, hss, L(inch=0.0635), 2, L(inch=1.90), "#52", degrees118, stub, center_cut, drills)
	drill_3_64 = shop._drill_append("3/64 [#0-80 75% thread]",
          23, hss, L(inch="3/64"), 2, L(inch=1.05), "3/64", degrees118, stub, center_cut, drills)
	drill_19 = shop._drill_append("#19 [M4x.7 close]",
	  24, hss, L(inch=0.1660), 2, L(inch=2.00), "M4x.7", degrees118, stub, center_cut, drills)
	drill_6mm = shop._drill_append("6mm drill [1.1in]",
	  25, hss, L(mm=6.00),     2, L(mm=28.00), "6mm", degrees135, stub, center_cut, drills)
	drill_8mm = shop._drill_append("8mm drill [3in]",
	  26, hss, L(mm=8.00),    2, L(inch="2-61/64"), "8mm", degrees135, stub, center_cut, drills)

	# Laser "tools":
	laser_007 = shop._end_mill_append("Laser_007",
	  100, hss, L(inch=0.007), 2, L(inch=0.750), laser)
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

    def _dowel_pin_append(self,
      name, number, alternate_number, material, diameter, maximum_z_depth, tip_depth):
	""" *Shop*: Create and return a *Tool_Dowel_Pin* object that contains a dowel pin.
	    The returned *Tool* object contains *name*, *number*, *material*, *diameter*,
	    *maximum_z_depth*, and *tip_depth*.   The returned *Tool* object, is also
	    appended to the tool list of the *Shop* (i.e. *self*.)
	"""

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(alternate_number, int) and alternate_number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(tip_depth, L)

	# Create the dowel pin *tool_dowel_pin*, add it to the *Shop* object (i.e. *self) tool list,
	# and return the *tool*:
	tool_dowel_pin = Tool_Dowel_Pin(name,
	  number, alternate_number, material, diameter, maximum_z_depth, tip_depth)
	tool_dowel_pin._feed_speed_set(Speed(in_per_sec = 1.000))
	self._tool_append(tool_dowel_pin)
	#print("tool_dowel_pin=", tool_dowel_pin)
	return tool_dowel_pin

    def _dowel_pin_search(self, tool):
	""" *Shop*: Return the *Tool_Dowel_Pin* that matches *tool* from the *Shop* object
	    (i.e. *self*) and *None* otherwise.
	"""

	# Use *shop* instead of *self*:
	shop = self

	# Verify argument types:
	assert isinstance(tool, Tool)

	# Search *tools* for a *Tool_Dowel_Pin* that matches *number*:
	dowel_pin_tool = None
	desired_tool_number = tool._number_get()
	for tool in shop._tools:
	    if tool._number_get() == desired_tool_number and isinstance(tool, Tool_Dowel_Pin):
		# We have a match:
		dowel_pin_tool = tool
		break

	return dowel_pin_tool


    def _drill_append(self, name, number, material, diameter,
      flutes_count, maximum_z_depth, drill_name, point_angle, drill_kind, is_center_cut, drills):
	""" *Shop*: Create and return a *Tool_Dowel_Pin* object that contains a dowel pin.
	    The returned *Tool* object contains *name*, *number*, *material*, *diameter*,
	    *flutes_count*, *maximum_z_depth*, *drill_name*, *point_angle*, *drill_kind*, and
	    *is_center_cut*.  *drills* is a list that the created *Tool_Drill* object is
	    appended to.  The returned *Tool_Drill* object, is also appended to the tool
	    list of the *Shop* (i.e. *self*.)
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter > zero
	assert isinstance(flutes_count, int) and flutes_count > 0
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero
	assert isinstance(drill_name, str)
	assert isinstance(point_angle, Angle)
	assert isinstance(drill_kind, int)
	assert isinstance(is_center_cut, bool)
	assert isinstance(drills, list)

	# Create the dowel pin *tool_dowel_pin*, add it to the *Shop* object (i.e. *self) tool list,
	# and return the *tool*:
	tool_drill = Tool_Drill(name, number, material, diameter,
	  drill_name, flutes_count, maximum_z_depth, point_angle, drill_kind, is_center_cut)
	self._tool_append(tool_drill)
	drills.append(tool_drill)
	return tool_drill

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

    def _mill_drill_append(self,
      name, number, material, diameter, flutes, maximum_z_depth, point_angle):
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
	assert isinstance(point_angle, Angle) and point_angle > Angle()

	# Create *tool_mill_drill*, append it to the *Shop* object (i.e. *self*)
	# tools list, and return it:
	tool_mill_drill = Tool_Mill_Drill(name,
	  number, material, diameter, flutes, maximum_z_depth, point_angle)
	self._tool_append(tool_mill_drill)
	return tool_mill_drill

    def _tools_get(self):
	""" *Shop*: Return the tools list object associated with the *Shop* object
	    (i.e. *self*.)
	"""

	return self._tools

    def _tools_summary_write(self):
	""" *Shop*: Write out a summary of CNC tool usage for the *Shop* object (i.e. *self*.)
	"""

	# Use *shop* instead of *self*:
	shop = self

	# Write out the summary file:
	summary_file = open("/tmp/tools.txt", "w")
	tools = shop._tools
	for tool in tools:
	    parts = tool._parts
	    if len(parts) > 0:
		summary_file.write("T{0}: {1}\n".format(tool._number, tool._name))
		parts.sort(key=lambda part: part.name_get())
		for part in parts:
		    summary_file.write("    {0}\n".format(part.name_get()))
	summary_file.close()

	# Write out the tools table:
	color_code = ["BLK", "BRN", "RED", "ORG", "YEL", "GRN", "BLU", "VIO", "GRY", "WHT"]
	with open("/tmp/tools_table.txt", "w") as tools_table_file:
	    tools_table_file.write("Tools_Table:\n")
	    tools_table_file.write("Tool Colors  Diameter Cut Depth Name\n")
	    tools_table_file.write("============================================================\n")
	    for tool in tools:
		number = tool._number_get()
		number_tens = (number / 10) % 10
		number_ones = number % 10
		#print("{0}: {1} {2}".format(number, number_tens, number_ones))
		maximum_z_depth = tool._maximum_z_depth_get()
		diameter = tool._diameter_get()
		name = tool._name_get()
		tools_table_file.write("T{0:02d}: {1}-{2} {3:i} {4:i}  {5}\n".format(
		  number, color_code[number_tens], color_code[number_ones],
		  diameter, maximum_z_depth, name))
		
	# Create a dictionary of collets and fill the with the range of acceptable sizes.
        # These values came off the ER-16 collets listed in the shars.com catalog:
	# We start with *collets* where each triple is ( collet_name, largest, smallest ):
	collets = (
	  ("1/32",  L(inch=0.031), L(inch=0.012)),	# 1/32
	  ("1/16",  L(inch=0.062), L(inch=0.023)),	# 2/32
	  ("3/32",  L(inch=0.093), L(inch=0.054)),	# 3/32
	  ("1/8",   L(inch=0.125), L(inch=0.086)),	# 4/32
	  ("5/32",  L(inch=0.156), L(inch=0.117)),	# 5/32
	  ("3/16",  L(inch=0.187), L(inch=0.148)),	# 6/32
	  ("7/32",  L(inch=0.218), L(inch=0.179)),	# 7/32
	  ("1/4",   L(inch=0.250), L(inch=0.211)),	# 8/32
	  ("9/32",  L(inch=0.281), L(inch=0.242)),	# 9/32
	  ("5/16",  L(inch=0.312), L(inch=0.273)),	# 10/32
	  ("11/32", L(inch=0.343), L(inch=0.304)),	# 11/32
	  ("3/8",   L(inch=0.375), L(inch=0.336)),	# 12/32
	  ("13/32", L(inch=0.406), L(inch=0.367)) )	# 13/32
	
	# Write out the `drills_summary.txt` file:
	drills = shop._drills
	drills.sort(key=lambda drill: drill._diameter.millimeters())
	with open("/tmp/drills_summary.txt", "w") as drill_file:
	    drill_file.write("Tool Drill Diameter Cut Depth Cnt Collets              Name\n")
	    drill_file.write("==================================================================\n")
	    for drill in drills:
		drill._summary_write(drill_file, collets)


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

    def _surface_speeds_lookup(self, part_material, tool_material, tracing=-1000000):
	""" *Shop*: Lookup and return a *Speed_Range* object from the
	    *Shop* object (i.e. *self*) keyed on *part_material* and
	    *tool_material*; or return *None* otherwise.
	"""

	# Verify argument types:
	assert isinstance(part_material, Material)
	assert isinstance(tool_material, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Shop._surface_speeds_lookup(*, '{1}', '{2}')".format(indent,
	      part_material, tool_material))

	# Look up the *Speed_Range* object from *surfaces_speeds_table*:
	key = (part_material._generic_get(), tool_material)
	#print("Shop.surface_speeds_lookup():key={0}".format(key))
	surface_speeds_table = self._surface_speeds_table
	speed_range = None
	if key in surface_speeds_table:
	    speed_range = surface_speeds_table[key]

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Shop._surface_speeds_lookup(*, '{1}', '{2}')=>{3:i}".format(indent,
              part_material, tool_material, speed_range))

	return speed_range

    def _tool_append(self, new_tool):
	""" *Shop*: Append *new_tool* to the tools lins in the *Shop* object (i.e. *self*).
	    Duplicate tool numbers cause an assertion failure.
	"""

	# Use *shop* instead of *self*:
	shop = self

	# Verify argument types:
	assert isinstance(new_tool, Tool)

	# Grab the *tools* list from *shop*:
	tools = shop._tools

	# Only allow dowel pins to be duplicated in *tools*:
	if not isinstance(new_tool, Tool_Dowel_Pin):
	    new_number = new_tool._number_get()
	    new_name = new_tool._name_get()
	    for tool in tools:
		number = tool._number_get()
		assert new_number != number, \
		  "Tool number {0} is conflicts between '{1}' and '{2}'".format(
		 number, new_tool._name_get(), tool._name_get())
		name = tool._name_get()
		#assert new_name != name, \
		#  "Tool name '{0}' has already been added (tool number {1}".format(name, number)

	# Append *new_tool* to *tools*:
	tools.append(new_tool)

    def _tooling_plate_get(self):
	""" *Shop*: Returns the currently selected tooling plate for the *Shop*object 
	    (i.e. *self*.)
	"""

	return self._tooling_plate

class STL:
    """ *STL*: This class represents a `.stl` file.
    """

    def __init__(self, part, tracing=-1000000):
	""" *STL*: Initialize the *STL* object (i.e. *self*) with the contents of the `.stl` file
	    generated for *part*.
	"""

	# Use *stl* instead of *self*:
	stl = self

	# Verify argument types:
	assert isinstance(part, Part)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>STL.__init__(*, '{1}')".format(indent, part._name_get()))
	    trace_detail = 1

	# Read in the .stl file that was generated by OpenSCAD:
	ezcad = part._ezcad_get()
	stl_directory = ezcad._stl_directory_get()
	stl_file_name = "{0}_{1}.stl".format(part._name, part._signature)
	stl_lines = stl_directory._lines_read(stl_file_name, tracing = tracing + 1)
	stl_lines_size = len(stl_lines)
	if trace_detail >= 1:
	    print("{0}stl_file_name='{1}'".format(indent, stl_file_name))
	    print("{0}len(stl_lines)={1}".format(indent, stl_lines_size))
    
	# Extract the *triangles* from *stl_lines*:
	triangles = []
	assert stl_lines[0][:5] == "solid", \
	  "'{0}' does not appear to be a .stl file".format(stl_file_name)
	index = 1
	while index + 4 < stl_lines_size:
		#  Extract *point1*, *point2*, and *point3* from *stl_lines*:
		xlist1 = stl_lines[index + 2].split()
		point1 = ( float(xlist1[1]), float(xlist1[2]), float(xlist1[3]) )
		xlist2 = stl_lines[index + 3].split()
		point2 = ( float(xlist2[1]), float(xlist2[2]), float(xlist2[3]) )
		xlist3 = stl_lines[index + 4].split()
		point3 = ( float(xlist3[1]), float(xlist3[2]), float(xlist3[3]) )

		# Create a *triangle* and append to *triangles*:
		triangle = ( point1, point2, point3)
		triangles.append(triangle)

		# Skip over to the next sequence of 7 lines:
		index += 7

	# Load *triangles* into *stl*:
	stl._stl_file_name = stl_file_name
	stl._triangles = tuple(triangles)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=STL.__init__(*, '{1}')".format(indent, part._name_get()))

    def _triangles_get(self):
	""" *STL*: Return a list of triangles from the *STL* object (i.e. *self*).
	"""

	return self._triangles
            
class Time:
    def __init__(self, sec=0.0, min=0.0):
	self._seconds = sec + min * 60.0 

    def __format__(self, format):
	""" *Time*: Return the *Time* object (i.e. *self) as a string formated by *format*.
	"""

	if format == 'm':
	    result = "{0:.2f}min".format(self._seconds / 60.0)
	else:
	    result = "{0:.2f}".format(self._seconds)
	return result

    def __add__(self, time2):
	""" Return the sum of the *Time* object (i.e. *self*) and *time2*.
	"""

	# Use *time1* instead of *self*:
	time1 = self

	# Verify argument types:
	assert isinstance(time2, Time)

	return Time(sec=time1._seconds + time2._seconds)

class Tool:
    """ *Tool*: A tool is a tool that can be used to manufacture a part. """

    # Tool kind:
    KIND_NONE = 0
    KIND_MOUNT = 1
    KIND_DOWEL_PIN = 2
    KIND_END_MILL = 3
    KIND_MILL_DRILL = 4
    KIND_DRILL = 5
    KIND_DOUBLE_ANGLE = 6
    KIND_DOVE_TAIL = 7
    KIND_LASER = 8
    KIND_LAST = 9

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
	""" *Tool*: Initialize a *Tool* object (i.e. *self*).  *name* is the tool name
	    and *number* is the slot in the tool table to put the tool.  *number* is negative
	    if there is no actual tool (i.e. for a mount operation.)  *kind* specifies
	    what kind of tool it is.  *material* specifies the material tha the tool is made
	    out of (e.g. HSS, Carbide, etc.).  *diameter* is the tool diameter.  Flutes count
	    specifies the number of flutes on the tool.  *maximum_z_depth* is maximum depth
	    that the tool can go into the material being worked on.
	"""

	# Use *tool* instead of *self*:
	tool = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(number, int)
	assert isinstance(kind, int) and Tool.KIND_NONE < kind < Tool.KIND_LAST
	assert isinstance(material, int)
	assert isinstance(diameter, L)
	assert isinstance(flutes_count, int)
	assert isinstance(maximum_z_depth, L)

	# Load up *self*:
	tool._name = name			# Tool name
	tool._number = number			# Tool number in rack
	tool._kind = kind
	tool._material = material		# Material tool is made of
	tool._diameter = diameter		# Diameter of tool
	tool._flutes_count = flutes_count	# Number of flutes or teeth
	tool._maximum_z_depth = maximum_z_depth	# Maximum Z depth allowed
	tool._search_results = []		# Result of array search
	tool._priority = 0.0			# Priority in tool search
	tool._feed_speed = Speed()		# Nominal feedrate
	tool._spindle_speed = Hertz()		# Preferred spindle speed
	tool._parts = []			# List of parts using tool

    def __format__(self, format):
	""" *Tool*: Return the *Tool* object (i.e. *self*) formatted as a string. """

	return self._name

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

    def _part_register(self, part):
	""" *Tool*: Register that *part* uses the *Tool* object (i.e. *self*.) """

	# Use *tool* instead of *self*:
	tool = self

	# Verify argument types:
	assert isinstance(part, Part)

	# Make sure *part* is in *parts* for *tool* without duplicates:
	parts = tool._parts
	if part not in parts:
	    parts.append(part)

    def _parts_get(self):
	""" *Tool*: Return the *Part*'s that use the *Tool* object (i.e. *self*.)
	"""

	return self._parts

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

    def _vrml_append(self, vrml, pad, tip, tip_depth, tracing=-1000000):
	""" *Tool*: Write out a visualization of the *Tool* object in VRML format to
	    *wrl_file* indented by *pad*.  The tip of the tool is located at *tip*.
	"""
	
	# Use *tool* instead of *self*:
	tool = self

	# Verify argument types:
	assert isinstance(vrml, VRML)
	assert isinstance(pad, int)
	assert isinstance(tip, P)
	assert isinstance(tip_depth, L)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Tool._wrl_write('{1}', *, {2}, {3:i}, {4:i})".
   	      format(indent, tool._name, pad, tip, tip_depth))

	# Grab some values from *tool*:
	maximum_z_depth = tool._maximum_z_depth
	diameter = tool._diameter
	radius = diameter/2

	zero = L()
	x0 = tip.x
	y0 = tip.y
	z0 = tip.z
	z1 = z0 + tip_depth
	z2 = z0 + maximum_z_depth

	upper_points = []
	lower_points = []
	count = 16
	sub_angle = 360.0 / count
	for index in range(count):
	    angle = Angle(deg=(sub_angle * float(index)))
	    x = x0 + radius.cosine(angle)
	    y = y0 + radius.sine(angle)
	    upper_point = P(x, y, z2)
	    lower_point = P(x, y, z1)
	    upper_points.append(upper_point)
	    lower_points.append(lower_point)

	points = []
	for index in range(0, count, 2):
	    points.append(tip)
	    points.append(lower_points[index])
	    points.append(upper_points[index])
	    points.append(upper_points[index + 1])
	    points.append(lower_points[index + 1])

	# Create *vrml* lines object to draw the tool outline:
	vrml._poly_line("purple", points)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Tool._wrl_write('{1}', *, {2}, {3:i}, {4:i})".
	      format(indent, tool._name, pad, tip, tip_depth))

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

    def __init__(self,
      name, number, alternate_number, material, diameter, maximum_z_depth, tip_depth):
	""" *Tool_Dowel_Pin*: Initialize *Tool_Dowel_Pin* object (i.e. *self*) with *name*,
	    *number*, *material*, *diameter*, *maximum_z_depth*, and *tip_depth*.  *material*
	    is the tool material (e.g. HSS, carbide, ...), *diameter* is the dowel pin diameter,
	    *maximum_z_depth* is the maximum distance the entire dowel pin can go down below
	    the part top surface, and *tip_depth* is the amount of the dowel pin at the tip
	    that is not *diameter* round.  (Dowel pins can have a flat, round or conical tip.)
	    There are conceptually two dowel pins to work around a LinuxCNC bug.  These are
	    *number* and *alternate_number*.  They must both have the same tool offset in
	    LinuxCNC.
	"""

	# Use *tool_dowel_pin* instead of *self*:
	tool_dowel_pin = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(number, int) and number > 0
	assert isinstance(alternate_number, int) and alternate_number > 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L)
	assert isinstance(maximum_z_depth, L)
	assert isinstance(tip_depth, L)

	# Initialize super class:
	Tool.__init__(self,
	  name, number, Tool.KIND_DOWEL_PIN, material, diameter, 0, maximum_z_depth)

	# Force the speed to be slow:
	tool_dowel_pin._feed_speed_set(Speed(in_per_sec=1.0))

	# Load up *self*:
	tool_dowel_pin._tip_depth = tip_depth		# Tip portion that is not vertical
	tool_dowel_pin._alternate_number = alternate_number

    def _alternate_number_get(self):
	""" *Tool_Dowel_Pin*: Return the alternate tool number for the *Tool_Dowel_Pin* object
	    (i.e. *self*.)
	"""

	return self._alternate_number

    @staticmethod
    def _match(tool, preferred_diameter, parameter2, from_string, tracing=-1000000):
	""" *Tool_Dowel_Pin*: Return priority of match *tool* matching a dowel pin.
	"""
	
	# Verify argument types:
	assert isinstance(tool, Tool)
	assert isinstance(preferred_diameter, L)
	assert isinstance(parameter2, L)
	assert isinstance(from_string, str)
	assert isinstance(tracing, int)

	# Determine the *priority* depending upon whether *tool* is a dowel pin and
	# how close the diameter is to the *preferred_diameter*.  A perfect match returns 100.0
        # and gets smaller as the diameter mismatch grows:
	priority = -1.0
	if isinstance(tool, Tool_Dowel_Pin):
	    priority = 100.0 - (tool._diameter - preferred_diameter).absolute().millimeters()

	    # Set *debug* to *True* to enable tracing:
	    debug = False
	    #debug = True
	    if debug:
		print("Tool_Dowel_Pin.match(): Tool:'{0}' preferred_diameter={1:i} priorit={2}".
	      format(tool._name_get(), preferred_diameter, priority))
	return priority

    def _tip_depth_get(self):
	""" *Tool_Dowel_Pin: Return the tip depth of the *Tool_Dowel_Pin* object (i.e. *self*.)
	"""

	return self._tip_depth

class Tool_Drill(Tool):

    def __init__(self, name, number, material, diameter, drill_name,
      flutes_count, maximum_z_depth, point_angle, drill_style, is_center_cut):
	""" *Tool_Drill*: Initialize *Tool_Drill* object (i.e. *self*) with *name*, *number*,
	    *material*, *diameter*, *flutes_count*, *maximum_z_depth*, *point_angle*,
	    *drill_style*, and *is_center_cut*.
	"""

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter > zero
	assert isinstance(drill_name, str)
	assert isinstance(flutes_count, int) and flutes_count > 0 
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero
	assert isinstance(point_angle, Angle) and point_angle > Angle()
	assert isinstance(drill_style, int) and \
	  Tool.DRILL_STYLE_NONE < drill_style < Tool.DRILL_STYLE_LAST
	assert isinstance(is_center_cut, bool)

	# Load up the *Tool_Drill* object (i.e. *self*):
	Tool.__init__(self, name, number, Tool.KIND_DRILL,
	  material, diameter, flutes_count, maximum_z_depth)
	self._point_angle = point_angle
	self._drill_style = drill_style
	self._drill_name = drill_name
	self._is_center_cut = is_center_cut

    def _drill_style_get(self):
	""" *Tool_Drill*: Return the drill style for the *Tool_Drill* object (i.e. *self*). """

	return self._drill_style

    def _is_center_cut_get(self):
	""" *Tool_Drill*: Return whether the *Tool_Drill* object (i.e. *self*) is a center
	    cut drill.
	"""

	return self._is_center_cut

    @staticmethod
    def _match(tool, desired_diameter, maximum_z_depth, from_routine, tracing=-1000000):
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
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Tool_Drill._match('{1}', {2:i}, {3:i}, '{4}')".
	      format(indent, tool, desired_diameter, maximum_z_depth, from_routine))
	    trace_detail = 3

	# Priority is negative if *tool* is not a *Tool_Drill* object:
	priority = -1.0
	if isinstance(tool, Tool_Drill):
	    # Make sure *tool* is long enough:
	    if tool._maximum_z_depth >= maximum_z_depth:
		# If the *tool* *diameter* is very close to *desired_diameter* we have match::
		diameter = tool._diameter
		epsilon = L(inch=0.00001)
		if trace_detail >= 3:
		    print("{0}desiried_diameter={1:i} diameter={2:i}".
		      format(indent, desired_diameter, diameter))
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

    def _summary_write(self, drill_file, collets):
	""" *Tool_Drill*: Write a summary of the *Tool_Drill* object (i.e. *self*) out to
	    *drill_file*:
	"""

	# Use *tool_drill* instead of *self*:
	tool_drill = self

	# Verify argument types:
	assert isinstance(drill_file, file)
	assert isinstance(collets, tuple)

	# Grab some extra values from *tool_drill*:
	diameter        = tool_drill._diameter
	maximum_z_depth = tool_drill._maximum_z_depth
	drill_name      = tool_drill._drill_name
	name            = tool_drill._name
	number          = tool_drill._number
	parts           = tool_drill._parts

	# Find ER16 collets that can take the *tool_drill*:
	collet_names = []
	extra = L(inch=0.001)
	for collet in collets:
	    collet_name = collet[0]
	    largest = collet[1]  + extra
	    smallest = collet[2] - extra
	    if smallest <= diameter <= largest:
		collet_names.append(collet_name)

	color_code = ["BLK", "BRN", "RED", "ORG", "YEL", "GRN", "BLU", "VIO", "GRY", "WHT"]

	# Write out one line:
	drill_file.write("T{0:02d}  {1:<5} {2:i} {3:i}  {4:>3} {5:<20} {6}\n".
	  format(number, drill_name, diameter, maximum_z_depth, len(parts), collet_names, name))

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
	    it will match any end mill.  A positive number that increases with the diameter
	    is returned if a match occurs.  Otherwise, -1.0 is returned if there is no match.
	    *from_routine* is the name of the calling routine and is used for debugging only.
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
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Tool_End_Mill._match('{1}', {2:i}, {3:i}, '{4}')".
	      format(indent, tool._name, maximum_diameter, maximum_z_depth, from_routine))
	    trace_detail = 1

	# Grab some values from *tool*:
	tool_diameter = tool._diameter
	tool_maximum_z_depth = tool._maximum_z_depth
    
	# See whether or not we can return a positive *priority*:
	priority = -123456.0
	is_end_mill = isinstance(tool, Tool_End_Mill)
	tool._search_results_append(is_end_mill, "Is end mill")
	if is_end_mill:
	    #if trace_detail >= 0:
	    #	print("{0}=>Tool_End_Mill.match('{1}', {2:i}, {3:i}, '{4}')".
	    #	  format(indent, tool._name, maximum_diameter, maximum_z_depth, from_routine))

	    # Somehow, the type *bool_* (from numpy) occasionally gets returned. The cast
	    # using bool() works around the problem.  This is just weird:
	    diameter_ok = bool(maximum_diameter < zero) or bool(tool_diameter <= maximum_diameter)
	    assert isinstance(diameter_ok, bool)
	    tool._search_results_append(diameter_ok,
	      "Diameter {0:i} < 0 or Diameter {0:i} <= Max Diameter {1:i}".format(
	      tool_diameter, maximum_diameter))
	    if diameter_ok:
		# Verify that tool depth works:
		if trace_detail >= 0:
		    print("{0}Tool_End_Mill._match: diameter_ok".format(indent))
		if trace_detail >= 1:
		    print("{0}Tool_End_Mill._match: max_z_depth:{1:i} tool_maximum_z_depth:{2:i}".
		      format(indent, maximum_z_depth, tool_maximum_z_depth))
		z_depth_ok = maximum_z_depth <= tool_maximum_z_depth
		#tool._search_results_append(z_depth_ok,
		#  "Max Z depth {0:i} <= Tool Max Z Depth {1:i}".format(
		#  -maximum_z_depth, tool_maximum_z_depth))
		if z_depth_ok:
		    if trace_detail >= 1:
			print("{0}Tool_End_Mill._match: z_depth ok".format(indent))
		    priority = (tool_diameter * 100.0 - tool_maximum_z_depth).inches()
		    if tool._is_laser:
			priority = tool_diameter.inches()

	    #if trace_detail >= 0:
	    #	print("{0}<=Tool_End_Mill._match('{1}', {2:i}, {3:i}, '{4}') => {5}".format(
	    #	  indent, tool._name, maximum_diameter, maximum_z_depth, from_routine, priority))

	assert isinstance(priority, float)

	# Wrap up any requested *tracing*:
	if trace_detail >= 1:
	    print("{0}<=Tool_End_Mill._match('{1}', {2:i}, {3:i}, '{4}') => {5}".
	      format(indent, tool._name, maximum_diameter, maximum_z_depth, from_routine, priority))

	return priority

    def _is_laser_get(self):
	""" *Tool_End_Mill*: Return whether the *Tool_End_Mill* (i.e. *self*) is a laser or not. """

	return self._is_laser

class Tool_Mill_Drill(Tool):		# A mill-drill bit

    def __init__(self,
      name, number, material, diameter, flutes_count, maximum_z_depth, point_angle):
	""" *Tool_Mill_Drill*: Initialize the *Tool_Mill_Drill* object (i.e. *self*) with *name*,
	    *number*, *material*, *diameter*, *flutes_count*, *maximum_z_depth*, and *point_angle*.
	"""

	# Use *tool_mill_drill* instead of *self*:
	tool_mill_drill = self

	# Verify argument types:
	zero = L()
	assert isinstance(name, str)
	assert isinstance(number, int) and number >= 0
	assert isinstance(material, int) and Tool.MATERIAL_NONE < material < Tool.MATERIAL_LAST
	assert isinstance(diameter, L) and diameter > zero
	assert isinstance(flutes_count, int) and flutes_count > 0
	assert isinstance(maximum_z_depth, L) and maximum_z_depth > zero
	assert isinstance(point_angle, Angle) and point_angle > Angle()

	# Initialize *Tool_Mill_Drill* object (i.e. *self*):
	Tool.__init__(self,
	  name, number, Tool.KIND_MILL_DRILL, material, diameter, flutes_count, maximum_z_depth)
	tool_mill_drill._point_angle = point_angle

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
	    print("{0}=>Tool_Mill_Drill._mill_drill_side_match('{1}', {2:i}, {3:i}, '{4}')".
	      format(indent, tool._name, maximum_diameter, maximum_z_depth, from_routine))

	priority = -1.0
	if isinstance(tool, Tool_Mill_Drill):
	    tip_depth = tool._tip_depth_get()
	    diameter = tool._diameter
	    zero = L()
	    if maximum_diameter < zero or diameter <= maximum_diameter:
		if maximum_z_depth >= -(tool._maximum_z_depth - tip_depth):
		    priority = diameter.millimeters()

	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Tool_Mill_Drill._mill_drill_side_match('{1}', {2:i}, {3:i}, '{4}')=>{5}".
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

    def _point_angle_get(self):
	""" *Tool_Mill_Drill*: Return the point angle for the *Tool_Mill_Drill*.
	"""

	return self._point_angle

    def _tip_depth_get(self):
	""" *Tool_Mill_Drill*: Return the tip depth for the *Tool_Mill_Drill* object
	    (i.e. *self*.)  The tip depth is the vertical distance from the tip to the
	    until the outer diameter of the mill drill is reached.
	"""

	# Use *tool_mill_drill* instead of *self*:
	tool_mill_drill = self

	# tip_depth = r * tan(90 - pa/2) where r is the radius and pa is the point angle:
	diameter = tool_mill_drill._diameter
	point_angle = tool_mill_drill._point_angle
	radius = diameter / 2
	tip_depth = radius * (Angle(deg=90.0) - point_angle/2).tangent()
	return tip_depth

class Tooling_Plate:
    """ *Tooling_Plate*: A *Tooling_Plate* object represents a flat plate of holes
	onto which parts can be mounted."""

    def __init__(self, name, dx, dy, dz, columns, rows, hole_pitch,
      hole_diameter, soft_drill_diameter, steel_drill_diameter, spacer_width, spacer_dz):
	""" *Tooling_Plate*:  Initialize the *Tooling_Plate* object to contain
	    *dx*, *dy*, *dz*, *rows*, *columns*, *hole_diameter*, *soft_drill_diameter*,
	    *steel_drill_diameter*, *spacer_width*, and *spacer_dz*.
	"""

	# Use *tooling_plate* instead of *self*:
	tooling_plate = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(dx, L)
	assert isinstance(dy, L)
	assert isinstance(dz, L)
	assert isinstance(columns, int) and columns > 0
	assert isinstance(rows, int) and rows > 0
	assert isinstance(hole_pitch, L)
	assert isinstance(hole_diameter, L)
	assert isinstance(soft_drill_diameter, L)
	assert isinstance(steel_drill_diameter, L)
	assert isinstance(spacer_width, L)
	assert isinstance(spacer_dz, L)

	# Save everything into *tooling_plate*:
	tooling_plate._columns = columns
	tooling_plate._dx = dx
	tooling_plate._dy = dy
	tooling_plate._dz = dz
	tooling_plate._hole_diameter = hole_diameter
	tooling_plate._hole_pitch = hole_pitch
	tooling_plate._name = name
	tooling_plate._rows = rows
	tooling_plate._soft_drill_diameter = soft_drill_diameter
	tooling_plate._spacer_dz = spacer_dz
	tooling_plate._spacer_width = spacer_width
	tooling_plate._steel_drill_diameter = steel_drill_diameter

	# Reset all of the mount settings:
	tooling_plate._mount_reset()

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

    def _dowel_point_get(self):
	""" *Tooling_Plate*: Return the dowel point stored the *Tooling_Plate* object (i.e. *self*.)
	"""

	dowel_point = self._dowel_point
	assert isinstance(dowel_point, P)
	return dowel_point

    def _dowel_point_set(self, dowel_point):
	""" *Tooling_Plate*: Return the dowel point stored the *Tooling_Plate* object (i.e. *self*.)
	"""

	# Use *tooling_plate* instead of *self*:
	tooling_plate = self

	# Verifiy argument types:
	assert isinstance(dowel_point, P)

	# Stuff *dowel_point* into *tooling_plate*:
	self._dowel_point = dowel_point

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

    def _extra_material_get(self, mount, column, row, rotate, tracing=-1000000):
        """ *Tooling_Plate*: Return the extra material bounding box for the *Part* asociated with
	    *mount* using the *Tooling_Plate* object (i.e. self).  The extra material bounding box
	    is specified a by a BSW (Bottom South West) point and a TNE (Top North East) point.
	    The extra material bounding box is positioned such that the *Part* upper left mounting
	    hole is positioned over the machine origin (i.e. top surface transfrom.)   The *Part*
	    is rotated around a Z axis that goes through the *Part* center and *rotate* must be
	    a multiple of 90 degrees.)
	"""

	# Use *tooling_plate* instead of *self*:
	tooling_plate = self

	# Verify argument types:
	assert isinstance(mount, Mount)
	assert isinstance(rotate, Angle)
	assert isinstance(column, int) # and 0 <= column < tooling_plate._columns
	assert isinstance(row, int)    # and 0 <= row    < tooling_plate._rows
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Tooling_Plate._extra_material_get('{1}', {2}, {3}, {4:d})".
	      format(indent, mount._name_get(), column, row, rotate))
	    trace_detail = 3

	# Convert *rotate* to an *int*:
	rotate_degrees = int(rotate.degrees())
	assert rotate_degrees % 90 == 0, \
	  "Rotate is {0:d} degrees when it needs to be a mulitple of 90 degrees".format(rotate)

	# Figure out the size of the grid of *tooling_plate* mount holes:
	tooling_plate_holes = mount._tooling_plate_holes_get()
	assert len(tooling_plate_holes) > 0, \
	  "Mount '{0}' has no tooling plate holes".format(mount._name_get())
	hole_pitch     = tooling_plate._hole_pitch
	unique_columns = tuple(sorted(set([ column_row[0] for column_row in tooling_plate_holes ])))
	assert len(unique_columns) > 0
	unique_rows    = tuple(sorted(set([ column_row[1] for column_row in tooling_plate_holes ])))
	assert len(unique_rows) > 0
	columns_span   = unique_columns[-1] - unique_columns[0]
	rows_span      = unique_rows[-1]    - unique_rows[0]
	columns_dx     = columns_span * hole_pitch
	rows_dy        = rows_span    * hole_pitch
	if trace_detail >= 3:
	    print("{0}tooling_plate_holes={1}".format(indent, tooling_plate_holes))
	    print("{0}hole_pitch={1:i}".format(indent, hole_pitch))
	    print("{0}unique_columns{1} unique_rows={2}".
	      format(indent, unique_columns, unique_rows))
	    print("{0}columns_span={1} rows_span={2}".format(indent, columns_span, rows_span))
	    print("{0}columns_dx={1:i} rows_dy={2:i}".format(indent, columns_dx, rows_dy))

	# Compute *rotate_transform* that rotates *part* by *rotate* around a Z axis that
        # goes through *mount_translate_point*:
	zero = L()
	one = L(mm=1.000)
	z_axis = P(zero, zero, one)
	rotate_transform = Transform().rotate("rotate@Z_Axis", z_axis, rotate)
	extra_start_bsw, extra_start_tne = mount._extra_start_get()
	if trace_detail >= 3:
	    print("{0}mount='{1}'".format(indent, mount._name_get()))
	    print("{0}rotate_transform={1:v}".format(indent, rotate_transform))
	    print("{0}extra_start_bsw={1:i} extra_start_tne={2:i}".
	      format(indent, extra_start_bsw, extra_start_tne))

	# Figure out where everything is in TS (i.e. top surface) coordinates:
	part = mount._part_get()
	ts_transform       = mount._top_surface_transform_get()
	ts_part_bsw        = ts_transform * part.bsw
	ts_part_center     = ts_transform * part.c
	ts_part_tne        = ts_transform * part.tne
	ts_part_volume     = (ts_part_tne - ts_part_bsw).absolute()
	ts_hole_bsw        = ts_part_center + P(-columns_dx/2, -rows_dy/2, -ts_part_volume.z/2)
	ts_hole_tne        = ts_part_center + P( columns_dx/2,  rows_dy/2,  ts_part_volume.z/2)
	ts_extra_start_bsw = ts_transform * extra_start_bsw
	ts_extra_start_tne = ts_transform * extra_start_tne
	if trace_detail >= 3:
	    print("{0}part='{1}' mount='{2}'".
	      format(indent, part._name_get(), mount._name_get()))
	    print("{0}ts_part_bsw={1:i} ts_part_center={2:i} ts_part_tne={3:i}".
	      format(indent, ts_part_bsw, ts_part_center, ts_part_tne))
	    print("{0}ts_part_volume={1:i}".format(indent, ts_part_volume))
	    print("{0}ts_hole_bsw={1:i} ts_hole_tne={2:i}".
	      format(indent, ts_hole_bsw, ts_hole_tne))
	    print("{0}ts_extra_start_bsw={1:i} ts_extra_start_tne={2:i}".
	      format(indent, ts_extra_start_bsw, ts_extra_start_tne))

	# Push a bunch of interesting points through the *rotate_transform*:
	rotated_ts_part_bsw    = rotate_transform * ts_part_bsw
	rotated_ts_part_center = rotate_transform * ts_part_center
	rotated_ts_part_tne    = rotate_transform * ts_part_tne
	rotated_ts_volume      = (rotated_ts_part_tne - rotated_ts_part_bsw).absolute()
	rotated_ts_hole_bsw    = rotate_transform * ts_hole_bsw
	rotated_ts_hole_tne    = rotate_transform * ts_hole_tne
	rotated_ts_extra_start_bsw = rotate_transform * ts_extra_start_bsw
	rotated_ts_extra_start_tne = rotate_transform * ts_extra_start_tne
	if trace_detail >= 3:
	    print("{0}rotated_ts_part_bsw={1:i} rotated_ts_part_tne={2:i}".
	      format(indent, rotated_ts_part_bsw, rotated_ts_part_tne))
	    print("{0}rotated_ts_volume={1:i}".
	      format(indent, rotated_ts_volume))
	    print("{0}rotated_ts_part_center={1:i} == ts_part_center={2:i}".
	      format(indent, rotated_ts_part_center, ts_part_center))
	    print("{0}rotated_ts_hole_bsw={1:i} rotated_ts_hole_tne={2:i}".
	      format(indent, rotated_ts_hole_bsw, rotated_ts_hole_tne))
	    print("{0}rotated_ts_extra_start_bsw={1:i} rotated_ts_extra_start_tne={2:i}".
	      format(indent, rotated_ts_extra_start_bsw, rotated_ts_extra_start_tne))

	# Fixup the corners after the *rotate_transform*:
	fixed_ts_part_bsw, fixed_ts_part_tne       = \
	  rotated_ts_part_bsw.minimum_maximum(rotated_ts_part_tne)
	fixed_ts_part_volume = (fixed_ts_part_tne - fixed_ts_part_bsw).absolute()
	fixed_ts_hole_bsw, fixed_ts_hole_tne       = \
	  rotated_ts_hole_bsw.minimum_maximum(rotated_ts_hole_tne)
	fixed_extra_start_bsw, fixed_extra_start_tne = \
	  rotated_ts_extra_start_bsw.minimum_maximum(rotated_ts_extra_start_tne)
	if trace_detail >= 3:
	    print("{0}fixed_ts_part_bsw={1:i} fixed_ts_part_tne={2:i}".
	      format(indent, fixed_ts_part_bsw, fixed_ts_part_tne))
	    print("{0}fixed_ts_part_volume={1:i}".
	      format(indent, fixed_ts_part_volume))
	    print("{0}fixed_ts_hole_bsw={1:i} fixed_ts_hole_tne={2:i}".
	      format(indent, fixed_ts_hole_bsw, fixed_ts_hole_tne))
	    print("{0}fixed_extra_start_bsw={1:i} fixed_extra_start_tne={2:i}".
	      format(indent, fixed_extra_start_bsw, fixed_extra_start_tne))
	
	# Compute *final_extra_start_bsw* and *final_extra_start_tne*:
	ts_hole_upper_left = P(fixed_ts_hole_bsw.x, fixed_ts_hole_tne.y, zero)
	corner_offset = P(column * hole_pitch, -row * hole_pitch, zero)
	final_extra_start_bsw = fixed_extra_start_bsw - ts_hole_upper_left + corner_offset
	final_extra_start_tne = fixed_extra_start_tne - ts_hole_upper_left + corner_offset
	if trace_detail >= 3:
	    print("{0}ts_hole_upper_left={1:i} corner_offset={2:i}".
	     format(indent, ts_hole_upper_left, corner_offset))
	    print("{0}final_extra_start_bsw={1:i} final_extra_start_tne={2:i}".
	      format(indent, final_extra_start_bsw, final_extra_start_tne))

	# Wrap up any requested *tracing* and return results:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=Tooling_Plate._extra_material_get('{1}', {2}, {3}, {4:d})=>{5:i}, {6:i}".
	      format(indent, mount._name_get(), column, row, rotate,
	      final_extra_start_bsw, final_extra_start_tne))
	return final_extra_start_bsw, final_extra_start_tne

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

    def _mount_get(self):
	""" *Tooling_Plate*: Return the *Mount* from the *Tooling_Plate* object (i.e. *self*.)
	"""

	mount = self._mount 
	assert isinstance(mount, Mount)
	return mount

    def _mount_reset(self):
	""" *Tooling_Plate*: Reset all the mount settings for the *Tooling_Plate* object
	    (i.e. *self*.)
	"""

	# Use *tooling_plate* instead of *self*:
	tooling_plate = self

	# Reset all of the mount settings to invalid values:
	tooling_plate._dowel_point           = None
	tooling_plate._mount                 = None
	tooling_plate._parallels_height      = None
	tooling_plate._spacers               = None

    def _mount_set(self, mount):
	""" *Tooling_Plate*: Store mount into the *Tooling_Plate* object (i.e. *self*.)
	"""

	# Use *tooling_plate* instead of *self*:
	tooling_plate = self

	# Verify argument types:
	assert isinstance(mount, Mount)

	# Stuff *mount_translate_point* into *tooling_plate*:
	tooling_plate._mount = mount

    def _name_get(self):
	""" *Tooling_Plate*: Return the *Name* from the *Tooling_Plate* object (i.e. *self*.)
	"""

	return self._name

    def _parallels_height_get(self):
	""" *Tooling_Plate*: Return the height of the parallels in the vice for the *Tooling_Plate*
	    object (i.e. *self*.)
	"""

	return self._parallels_height

    def _parallels_height_set(self, parallels_height):
	""" *Tooling_Plate*: Store vice *parallels_heigh*t into the *Tooling_Plate* object
	    (i.e. *self*.)
	"""

	# Use *tooling_plate* instead of *self*:
	tooling_plate = self

	# Verify argument types:
	assert isinstance(parallels_height, L)

	# Stuff *parallels_height* into *tooling_plate*:
	tooling_plate._parallels_height = parallels_height

    def _rows_get(self):
	""" *Tooling_Plate*: Return the number of rows of mounting holes for the *Tooling_Plate*
	    object (i.e. *self*.)
	"""

	return self._rows

    def _spacers_get(self):
	""" *Tooling_Plate*: Return the a list of spacer quadruples used for mounting
	    on the *Tooling_Plate* object (i.e. *self*.)  Each quaduple is a tuple of the
	    form (x1, y1, x2, y2) where x1 and x2 represent a tooling plate column number,
	    and y1 and y2 represent the tooling plate row number.
	"""

	spacers = self._spacers
	assert isinstance(spacers, list)
	return spacers

    def _spacers_set(self, spacers):
	""" *Tooling_Plate*: Set the *spacers* quadruples used for mounting on the
	    *Tooling_Plate* object (i.e. *self*.)  Each quaduple is a tuple of the
	    form (x1, y1, x2, y2) where x1 and x2 represent a tooling plate column number,
	    and y1 and y2 represent the tooling plate row number.
	"""

	# Use *tooling_plate* instead of *self*:
	tooling_plate = self

	# Verify argument types:
	assert isinstance(spacers, list)
	for spacer in spacers:
	    assert isinstance(spacer, tuple) and len(spacer) == 4
	    for row_column in spacer:
		assert isinstance(row_column, int) and row_column >= 0

	# Stuff *spacers* into *tooling_plate*
	tooling_plate._spacers = spacers

    def _spacer_width_get(self):
	""" *Tooling_Plate*: Return the number of spacer width for the *Tooling_Plate*
	    object (i.e. *self*.)
	"""

	return self._spacer_width

    def _spacer_dz_get(self):
	""" *Tooling_Plate*: Return the number of spacer dz for the *Tooling_Plate*
	    object (i.e. *self*.)
	"""

	return self._spacer_dz

    def _vrml_append(self, vrml, corner, spacers, mount_ngc_file, tracing = -1000000):
	""" *Tooling_Plate*: Write out a VRML visulazation of *Tooling_Plate* (i.e. *self*)
	    to *vrml*.  The bottom north west *corner* coordinate is specified.
	"""

	# Use *tooling_plate* instead of *self*:
	tooling_plate = self
	 
	# Verify argument types:
	assert isinstance(vrml, VRML)
	assert isinstance(corner, P)
	assert isinstance(tracing, int)
	assert isinstance(spacers, list)
	assert isinstance(mount_ngc_file, file)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' *tracing
	    print("{0}=>Tooling_Plate._vrml_append('{1}', '{2}', {3:i}, {4})".
	      format(indent, tooling_plate._name, vrml._name_get(), corner, spacers))
	    trace_detail = 1

	# Grab some values from *tooling_plate*
	columns = tooling_plate._columns
	dx = tooling_plate._dx
	dy = tooling_plate._dy
	dz = tooling_plate._dz
	hole_pitch = tooling_plate._hole_pitch
	rows = tooling_plate._rows
	#spacers = tooling_plate._spacers
	spacer_dz = tooling_plate._spacer_dz
	spacer_width = tooling_plate._spacer_width

	extra_dx = dx - hole_pitch * (columns - 1)
	extra_dy = dy - hole_pitch * (rows - 1)
	if trace_detail >=1:
	    print("{0}extra_dx={1:i} extra_dy={2:i}".format(indent, extra_dx, extra_dy))

	zero = L()

	left_column_x = extra_dx/2
	right_column_x = extra_dx/2 + hole_pitch * float(columns - 1)
	top_row_y =    -extra_dy/2
	bottom_row_y = -(extra_dy/2 + hole_pitch * float(rows - 1))
	if trace_detail >= 1:
	    print("{0}left_column_x={1:i} right_column_x={2:i} top_row_y={3:i} bottom_row_y{4:i}".
	      format(indent, left_column_x, right_column_x, top_row_y, bottom_row_y))

	# Draw the tooling plate lines:
	corner1 = corner
	corner2 = corner1 + P(dx, -dy, dz)
	color_name = "orange"
	vrml._box_outline(color_name, corner1, corner2)
	for column in range(columns):
	    x = left_column_x + hole_pitch * float(column)
	    top_point    = corner1 + P(x, top_row_y,    dz)
	    bottom_point = corner1 + P(x, bottom_row_y, dz)
	    vrml._poly_line(color_name, [top_point, bottom_point])
	    if trace_detail >= 2:
		print("{0}column[{1}]: top_point={2:i} bottom_point={3:i}".
		  format(indent, column, top_point, bottom_point))
	for row in range(rows):
	    y = top_row_y - hole_pitch * float(row)
	    left_point =  corner1 + P(left_column_x,  y, dz)
	    right_point = corner1 + P(right_column_x, y, dz)
	    vrml._poly_line(color_name, [left_point, right_point])

	# Now draw the *spacers*:
	color_name = "aqua"
	for index, spacer in enumerate(spacers):
	    # Extract the row/column values from the *spacer* quadruple:
	    column1 = spacer[0]
	    row1 = spacer[1]
	    column2 = spacer[2]
	    row2 = spacer[3]
	    assert column1 <= column2
	    assert row1 <= row2

	    # Provide a hint for the tooling plate spacers:
	    hole_count = (column2 - column1) + (row2 - row1) + 1
	    mount_ngc_file.write("( Spacer[{0}]: {1} holes )\n".format(index, hole_count))

	    # Draw the box outline of the *spacer*:
	    x1 = left_column_x + hole_pitch * float(column1)
	    y1 = top_row_y - hole_pitch * float(row1)
	    x2 = left_column_x + hole_pitch * float(column2)
	    y2 = top_row_y - hole_pitch * float(row2)
	    if row1 == row2:
		# Horizontal spacer:
		corner1 = corner + P(x1 - spacer_width/2, y1 - spacer_width/2, dz)
		corner2 = corner + P(x2 + spacer_width/2, y2 + spacer_width/2, dz + spacer_dz)
	    else:
		# Vertical space:
		corner1 = corner + P(x1 + spacer_width/2, y1 + spacer_width/2, dz)
		corner2 = corner + P(x2 - spacer_width/2, y2 - spacer_width/2, dz + spacer_dz)
	    vrml._box_outline(color_name, corner1, corner2)

	    # Draw some lines to show where the spacer holes belong:
	    p1 = corner + P(x1, y1, dz)
	    p2 = corner + P(x2, y2, dz)
	    p3 = corner + P(x2, y2, dz + spacer_dz)
	    p4 = corner +P(x1, y1, dz + spacer_dz)
	    vrml._poly_line(color_name, [p1, p2, p3, p4, p1])

	    if trace_detail >= 1:
		print("{0}spacers[{1}]: {2}".format(indent, index, spacer))

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Tooling_Plate._vrml_append('{1}', '{2}', {3:i}, {4})".
	      format(indent, tooling_plate._name, vrml._name_get(), corner, spacers))

    def _xy_rapid_safe_z_get(self):
	""" *Tooling_Plate*: Return the top surface transform previoulsy set for the *Tooling_Plate*
	    object (i.e. *self*.)
	"""

	return self._xy_rapid_safe_z

    def _xy_rapid_safe_z_set(self, xy_rapid_safe_z):
	""" *Tooling_Plate*: Store a top surface transform into the *Tooling_Plate* object
	    (i.e. *self*.)
	"""

	# Use *tooling_plate* instead of *self*:
	tooling_plate = self

	# Verify argument types:
	assert isinstance(xy_rapid_safe_z, L)

	# Stuff *xy_rapid_safe_z* into *tooling_plate*:
	tooling_plate._xy_rapid_safe_z = xy_rapid_safe_z

class Transform:

    # The matrix format is an afine 4x4 matrix in the following format:
    #
    #	[ r00 r01 r02 0 ]
    #   [ r10 r11 r12 0 ]
    #   [ r20 r21 r22 0 ]
    #   [ dx  dy  dz  1 ]
    #
    # The afine point format is a 1x4 matrix of the following format:
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
    # This allows us to use the Transform object to set up openscad and VRML.

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
	self._forward_vrmls = ()
	self._reverse_vrmls = ()

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
	elif format.startswith("v"):
	    forward_vrmls = self._forward_vrmls
	    for forward_vrml in forward_vrmls:
		if len(forward_vrml) == 3:
		    # Translate:
		    result += " Translate[{0} {1} {2}]". \
                      format(forward_vrml[0], forward_vrml[1], forward_vrml[2])
		elif len(forward_vrml) == 4:
		    result += " Rotate[{0} {1} {2}  {3}]". \
                      format(forward_vrml[0], forward_vrml[1], forward_vrml[2],
		             forward_vrml[3] * 180.0 / math.pi)
		else:
		    assert False, "Bad forward_vrml: {0}".format(forward_vrml)
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

    def center_rotate(self, comment, center, axis, angle, tracing=-1000000):
	""" *Transform*: Append a rotation transform to the *Transform* object (i.e. *self*).
	    The rotation occurs around *center* along an *axis* by an amount of *angle*.
	"""

	# Use *transform *instead of *self*:
	transform = self

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(center, P)
	assert isinstance(axis, P)
	assert isinstance(angle, Angle)
	assert isinstance(tracing, int)

	# Perform any requested tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Transform.center_rotate('{1}', {2:i}, {3:i}, {4:d}".
	      format(indent, comment, center, axis, angle))
	    trace_detail = 1

	# Move to *center* to origin, rotate around *axis* by *angle*, and move back to *center*:
	transform = transform.translate("{0}: Move center to origin".format(comment), -center)
	transform = transform.rotate("{0}: Rotate along axis".format(comment), axis, angle)
	transform = transform.translate("{0}: Move back to center".format(comment), center)

	# Wrap up any *tracing* and return *transform*:
	if tracing >= 0:
	    print("{0}<=Transform.center_rotate('{1}', {2:i}, {3:i}, {4:d})=>{5:v}".
	      format(indent, comment, center, axis, angle, transform))
	return transform

    def chain(self, next_transform):
	""" *Transform*: Return a new transform that consists of the *Transform* object
	    (i.e. *self*) followed by *next_transform*.
	"""

	# Use *transform* instead of *self*:
	transform = self

	# Verify argument types:
	assert isinstance(next_transform, Transform)

	# Create the *result* *Transform* object:
	result = Transform()

	# Construct the forward and reverse matrices:
	result._forward_matrix = transform._forward_matrix.dot(next_transform._forward_matrix)
	result._reverse_matrix = numpy.linalg.inv(result._forward_matrix)
	
	# Construct the forward and reverse scad_lines:
	result._forward_scad_lines = \
	  next_transform._forward_scad_lines + transform._forward_scad_lines
	result._reverse_scad_lines = \
	  transform._reverse_scad_lines + next_transform._reverse_scad_lines

	# Construct the forward and reverse VRML transforms:
	result._forward_vrmls = transform._forward_vrmls + next_transform._forward_vrmls
	result._reverse_vrmls = next_transform._reverse_vrmls + transform._reverse_vrmls

	# Return the *result*:
	return result

    def reverse(self):
	""" *Transform*: Return the reverse of the *Transform* object (i.e. *self*). """

	assert len(self._forward_scad_lines) == len(self._reverse_scad_lines)

	result = Transform()
	result._forward_matrix = self._reverse_matrix
	result._reverse_matrix = self._forward_matrix
	result._forward_scad_lines = self._reverse_scad_lines
	result._reverse_scad_lines = self._forward_scad_lines
	result._forward_vrmls = self._reverse_vrmls
	result._reverse_vrmls = self._forward_vrmls

	assert len(result._forward_scad_lines) == len(result._reverse_scad_lines)
	assert len(result._forward_vrmls) == len(result._reverse_vrmls)

	return result

    def rotate(self, comment, axis, angle, tracing = -1000000):
	""" *Transform*: Return a rotation of the *Transform* object (i.e. *self*) rotated
	    by *angle* around *axis*.  *comment* shows up in the .scad file
	"""

	# Use *transform* instead of *self*:
	transform = self

	# Verify argument_types:
	assert isinstance(comment, str)
	assert isinstance(axis, P)
	assert isinstance(angle, Angle)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Transform.rotate('{1}', {2:m}, {3:d})".format(indent, comment, axis, angle))
	    trace_detail = 1

	# Is this a null rotation?:
	if angle == Angle():
	    # Yes it is a null rotation, so we can return *transform*:
	    if tracing >= 0:
		print("{0}null rotation: {1:d}".format(indent, angle))
	    result = transform
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
	    result._forward_matrix = transform._forward_matrix.dot(forward_matrix)
	    result._reverse_matrix = numpy.linalg.inv(result._forward_matrix)

	    # Create the *forward_scad_line* and the *reverse_scad_line*
	    forward_scad_line =                               \
	      ("rotate(a={0:d}, v=[{1}, {2}, {3}]) //F: {4}". \
	      format( angle, nx, ny, nz, comment), )
	    reverse_scad_line =                               \
	      ("rotate(a={0:d}, v=[{1}, {2}, {3}]) //R: {4}". \
	      format(-angle, nx, ny, nz, comment), )
	    result._forward_scad_lines = forward_scad_line + transform._forward_scad_lines
	    result._reverse_scad_lines = transform._reverse_scad_lines + reverse_scad_line

	    # Create the *forward_vrml* and the *reverse_vrml*:
	    angle_radians = angle.radians()
	    forward_vrml = (nx, ny, nz,  angle_radians )
	    reverse_vrml = (nx, ny, nz, -angle_radians )
	    if trace_detail >= 3:
		print("{0}before forward_vrmls={1}".format(indent, transform._forward_vrmls))
		print("{0}before reverse_vrmls={1}".format(indent, transform._reverse_vrmls))
	    result._forward_vrmls = transform._forward_vrmls + ( forward_vrml, )
	    result._reverse_vrmls = ( reverse_vrml, ) + transform._reverse_vrmls
	    if trace_detail >= 3:
		print("{0}before forward_vrmls={1}".format(indent, result._forward_vrmls))
		print("{0}before reverse_vrmls={1}".format(indent, result._reverse_vrmls))

	# Wrap-up any requested *tracing* and return *result*:
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
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Transform.top_surface('{1}', {2:i}, {3:i}, {4:d})". \
	      format(indent, comment, start, end, rotate))
	    trace_detail = 1

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
	if trace_detail >= 2:
	    print("{0}translate_transform_f={1:s}".format(indent, transform))
	if trace_detail >= 3:
	    print("{0}translate_transform_r={1:s}\n".format(indent, transform.reverse()))

	direction_axis = end - start
	direction = direction_axis.normalize()
	if trace_detail >= 1:
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
	    if trace_detail >= 3:
		print("{0}negative_z_axis={1:i} direction_axis{2:i}".
		  format(indent, negative_z_axis, direction_axis))
	    #plane_change_axis = negative_z_axis.cross_product(direction_axis).normalize()
	    plane_change_axis = direction_axis.cross_product(negative_z_axis).normalize()
	    rotate_angle = negative_z_axis.angle_between(direction_axis)
	    if trace_detail >= 1:
		print("{0}plane_change_axis={1:m} rotate_angle={2:d}".
		  format(indent, plane_change_axis, rotate_angle))
	    transform = \
	      transform.rotate("{0}: plane change".format(comment),
	      plane_change_axis, rotate_angle, tracing = tracing + 1)
	    if trace_detail >= 2:
		print("{0}plane_change transform_f={1:s}".format(indent, transform))
	    if trace_detail >= 3:
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
	    if trace_detail >= 2:
		print("{0}pre_rotate transform_f={1:s}".format(indent, transform))
	    if trace_detail >= 3:
		print("{0}pre_rotate transform_r={1:s}\n".format(indent, transform.reverse()))

	# Finally, perform the user requested *rotate*:
	transform = transform.rotate("{0}: user rotate".format(comment),
	  z_axis, rotate, tracing = tracing + 1)
	if trace_detail >= 2:
	    print("{0}user_rotate transform_f={1:s}".format(indent, transform))
	if trace_detail >= 2:
	    print("{0}user_rotate transform_r={1:s}\n".format(indent, transform.reverse()))

	# Perform any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Transform.top_surface('{1}', {2:i}, {3:i}, {4:d})=>{5:s}".
	      format(indent, comment, start, end, rotate, transform))

	return transform

    def translate(self, comment, dx_dy_dz, tracing=-1000000):
	""" Transform*: Perform a translate of the *Transform* object (i.e. *self*) by *dx_dy_dz*.
	    *comment* will show up in the .scad file.
	"""

	# Use *transform* instead of *self*:
	transform = self

	# Verify argument types:
	assert isinstance(comment, str)
	assert isinstance(dx_dy_dz, P)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Transform.translate(*, '{1}', {2:i})".format(indent, comment, dx_dy_dz))
	    trace_detail = 1

	if trace_detail >= 1:
	    print("{0}transform={1:v}".format(indent, transform))

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
	    # The translate is null, so we can return *transform*:
	    result = transform
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

	    # Now we can create the new *result*; start with matrix representation:
	    result = Transform()
	    result._forward_matrix = transform._forward_matrix.dot(forward_matrix)
	    result._reverse_matrix = numpy.linalg.inv(result._forward_matrix)

	    # Now do the openscad representation:
	    forward_scad_line = \
	      ("translate([{0}, {1}, {2}]) //F: {3}".format( dx,  dy,  dz, comment), )
	    reverse_scad_line = \
	      ("translate([{0}, {1}, {2}]) //R: {3}".format(ndx, ndy, ndz, comment), )
	    result._forward_scad_lines = forward_scad_line + transform._forward_scad_lines
	    result._reverse_scad_lines = transform._reverse_scad_lines + reverse_scad_line

	    # Now do the VRML representation:
	    forward_vrml = ( dx,  dy,  dz)
	    reverse_vrml = (-dx, -dy, -dz)
	    if trace_detail >= 3:
		print("{0}before forward_vrmls={1}".format(indent, transform._forward_vrmls))
		print("{0}before reverse_vrmls={1}".format(indent, transform._reverse_vrmls))
	    result._forward_vrmls = transform._forward_vrmls + ( forward_vrml, )
	    result._reverse_vrmls = ( reverse_vrml, ) + transform._reverse_vrmls
	    if trace_detail >= 3:
		print("{0}before forward_vrmls={1}".format(indent, result._forward_vrmls))
		print("{0}before reverse_vrmls={1}".format(indent, result._reverse_vrmls))
	
	# Wrap up any requested *tracing* and return *result*:
	if tracing >= 0:
	    print("{0}<=Transform.translate(*, '{1}', {2:i})=>*".format(indent, comment, dx_dy_dz))
	return result

    def _vrml(self, vrml, tracing=-1000000):
	""" *Transform*: Return a new *VRML* object where *vrml* has has the *Transform* object
	    (i.e. *self*) applied to it.
	"""

	# Use *transform* instead of *self*:
	transform = self

	# Verify argument types:
	assert isinstance(vrml, VRML)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Transform._vrml(*, '{1}')".format(indent, vrml._name_get()))

	# Iterate through each of the *forward_transforms*:
	new_vrml = vrml
	for forward_vrml in transform._forward_vrmls:
	    old_vrml = new_vrml
	    old_vrml_name = old_vrml._name_get()
	    if len(forward_vrml) == 3:
		# We have a translate:
		dx = L(mm=forward_vrml[0])
		dy = L(mm=forward_vrml[1])
		dz = L(mm=forward_vrml[2])
		translate = P(dx, dy, dz)
		new_vrml_name = "{0}__Translate".format(old_vrml_name)
		new_vrml = VRML_Translate(new_vrml_name, translate, old_vrml, tracing = tracing + 1)
	    elif len(forward_vrml) == 4:
		# We have a rotation:
		nx = L(mm=forward_vrml[0])
		ny = L(mm=forward_vrml[1])
		nz = L(mm=forward_vrml[2])
		axis = P(nx, ny, nz)
		amount = Angle(rad=forward_vrml[3])
		new_vrml_name = "{0}__Rotate".format(old_vrml_name)
		new_vrml = VRML_Rotate(new_vrml_name, axis, amount, old_vrml, tracing = tracing + 1)
	    else:
		assert False, "Bad VRML transform"
	
	# Wrap up any requested *tracing* and return *new_vrml*:
	if tracing >= 0:
	    print("{0}<=Transform._vrml(*, '{1}')".format(indent, vrml._name_get()))
	return new_vrml

class Vice:

    def __init__(self, name, volume, parallels, tool_change_point):
	""" *Vice*: Initialize the *Vice* object (i.e. *self*) to contain *volume*,
	    "parallels_thickness", and "parallels.  *volume* is a point (i.e. *P*)
	    that specifies the vice volume dimensions between the vice jaws.
	    *tool_change_point* specifies where tool changes occur relative to the vice origin.
	"""

	# Use *vice* instead of *self*:
	vice = self

	# Verify argument types:
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(volume, P)
	assert isinstance(parallels, Parallels)
	assert isinstance(tool_change_point, P)

	# Load up the *Vice* object (i.e. *self*):
	vice._name              = name
	vice._parallels         = parallels
	vice._tool_change_point = tool_change_point
	vice._volume            = volume
	#print("Vice.__init__():Vice._volume={0:i}".format(volume))

    def _coordinates_vrml_append(self, vrml, tracing=-1000000):
	""" *Vice*: Append the vice origin coordinates for the *Vice* object (i.e. *self*)
	    to *vrml*.
	"""

	# Verify argument types:
	assert isinstance(vrml, VRML_Lines)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>Vice._coordinates_vrml_append('{1}', '{2}')".
	      format(indent, vice._name, vrml._name_get()))

	# Define some constants:
	zero = L()
	tip = L(inch=1.000)
	tip_offset = L(inch=0.100)

	# Define the various points that make up the coordinate axes:
	origin = P(zero, zero, zero)
	x_tip = P(tip, zero, zero)
	x_tip1 = P(tip - tip_offset,  tip_offset, zero)
	x_tip2 = P(tip - tip_offset, -tip_offset, zero)
	y_tip = P(zero, tip, zero)
	y_tip1 = P(zero, tip - tip_offset,  tip_offset)
	y_tip2 = P(zero, tip - tip_offset, -tip_offset)
	z_tip = P(zero, zero, tip)
	z_tip1 = P(zero,  tip_offset, tip - tip_offset)
	z_tip2 = P(zero, -tip_offset, tip - tip_offset)

	# Draw the three corordinate axes into *vrml*:
	vrml._poly_line("red",   [origin, x_tip, x_tip1, x_tip2, x_tip])
	vrml._poly_line("green", [origin, y_tip, y_tip1, y_tip2, y_tip])
	vrml._poly_line("blue",  [origin, z_tip, z_tip1, z_tip2, z_tip])

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=Vice._coordinates_vrml_append('{1}', '{2}')".
	      format(indent, vice._name, vrml._name_get()))

    def _jaw_volume_get(self):
	""" *Vice*: Return the jaw volume of the *Vice* object (i.e. *self*).
	"""

	return self._volume

    def _name_get(self):
	""" *Vice*: Return the name of the *Vice* object (i.e. *self*.)
	"""

	return "NO VICE NAME YET"


    def _parallels_get(self):
	""" *Vice*: Return the *Parallels* object for the *Vice* object (i.e. *self*).
	"""

	return self._parallels

    def _tool_change_point_get(self):
	""" *Vice*: Return the tool change point relative to the *Vice* object (i.e. *self*)
	    origin.
	"""

	return self._tool_change_point

class VRML:
    """ *VRML* is used to output VRML (Virtual Reality Markup Language) file.  """

    def __init__(self, name):
	""" *VRML*: Initialize the *VRML* object (i.e. *self*) to prepare for VRML
	    output generation.
	"""

	# Use *vrml* instead *self*:
	vrml = self

	# Verify argument types:
	assert isinstance(name, str)

	# 
	lines = []
	lines.append("\n")
	lines.append("# {0}\n".format(name))

	# Load up *vrml*:
	vrml._is_open = True
	vrml._lines = lines
	vrml._name = alphanumeric_only(name.replace(' ', '_'))
	vrml._pad_table = {}

    def _close(self):
	""" *VRML*: Make sure that the *VRML* object (i.e. *self*) is closed.  This routine
	    is overridden by a sub-class as needed:
	"""

	self._is_open = False

    def _name_get(self):
	""" *VRML*: Return the name of the *VRML* object (i.e. *self*).
	"""

	return self._name

    def _text_pad(self, pad, tracing=-1000000):
	""" *VRML*: Return a string of VRML text for the *VRML* object (i.e. *self*.)
	"""

	assert False, "No _text_pad() routine for {0}".format(self.__class__.__name__)

    def _text_pad_helper(self, pad, tracing=-1000000):
	""" *VRML*: Return a string of VRML text for the *VRML* object (i.e. *self*.)
	"""

	# Use *vrml* instead of *self*:
	vrml = self
	
	# Verify argument types:
	assert isinstance(pad, int)
	assert isinstance(tracing, int)

	# Perform an requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>VRML._text_pad_helper('{1}', {2})".format(indent, vrml._name, pad))
	    trace_detail = 2

	# Make sure all *vrml* is closed:
	if vrml._is_open:
	    vrml._close()
	    assert not vrml._is_open

	# Convert *lines* to *text*:
	pad_table = vrml._pad_table
	if not pad in pad_table:
	    lines = vrml._lines
	    if trace_detail >= 2:
		print("{0}lines[0]='{1}' lines[2]='{2}' lines[3]='{3}".
		  format(indent, lines[0], lines[1], lines[2]))

	    # This is the first time at this padding level; remember the result:
	    padding = ' ' * pad
	    text = padding.join(lines)
	    pad_table[pad] = text
	text = pad_table[pad]

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=VRML._text_pad_helper('{1}', {2})".format(indent, vrml._name, pad))
	return text

class VRML_Group(VRML):
    """ *VRML_Group*: Represents a VRML Group node.
    """

    def __init__(self, name):
	""" *VRML_Group*: Initialize the *VRML_Group* object (i.e. *self*):
	"""

	# Initliaze super-class:
	VRML.__init__(self, name)

	# Create *children* lists:
	self._children = []

    def _append(self, vrml):
	""" *VRML_Group*: Append *vrml* to the *VRML_Group* object (i.e. *self*):
	"""

	# Use *vrml_group* instead of *self*:
	vrml_group = self

	# Verify argument types:
	assert isinstance(vrml, VRML)

	# Make sure that *vrml_group* is still open:
	assert vrml_group._is_open

	# Append *vrml* to *children*:
	vrml_group._children.append(vrml)

    def _text_pad(self, pad, tracing=-1000000):
	""" *VRML_Group*: Return the *VRML_Group* object (i.e. *self*) as a single string 
	"""

	# Use *vrml_group* instead of *self*:
	vrml_group = self

	# Verify argument types:
	assert isinstance(pad, int)	
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>VRML_Group._text_pad('{1}', {2})".format(indent, vrml_group._name, pad))
	    trace_detail = 2

	# Mark this group as *closed*:
	vrml_group._is_open = False

	# Construct the content
	pad_table = vrml_group._pad_table
	if pad in pad_table:
	    # This padding level has been previously computed:
	    text = pad_table[pad]
	else:
	    # Put all the needed text pieces into *lines*:
	    lines = []
	    lines.append("\n")
	    lines.append("### VRML_Group\n")
	    lines.append("DEF {0} Group {{\n".format(vrml_group._name))
	    lines.append(" children [\n")
	    for index, vrml_child in enumerate(vrml_group._children):
		if trace_detail >= 2:
		    print("{0}[{1}] '{2}'".format(indent, index, vrml_child.__class__.__name__))
		vrml_child_text = vrml_child._text_pad(pad + 2, tracing = tracing + 1)
		assert isinstance(vrml_child_text, str), \
		  "Bad return for {0}".format(vrml_child.__class__.__name__)
		lines.append(vrml_child_text)
	    lines.append(" ]\n")
	    lines.append("}\n")

	    # Concatenate all of lines into one big string of *text* padded with *padding*::
	    padding = ' ' * pad
	    text = padding.join(lines)

	    # Stuff *text* into *pad_table*:
	    pad_table[pad] = text

	# Wrap up any requested *tracing* and return *text*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=VRML_Group._text_pad('{1}', {2})".format(indent, vrml_group._name, pad))
	return text

class VRML_Lines(VRML):
    """ *VRML_Lines* represents a bunch of colored VRML lines.
    """

    def __init__(self, name):
	""" *VRML_Lines*: 
	"""

	# Use *vrml_lines* instead of *self*:
	vrml_lines = self

	# Verify argument types:
	assert isinstance(name, str)

	# Initialize *VRML* super-class:
	VRML.__init__(vrml_lines, name)

	# Intialize *vrml_lines*:
	vrml_lines._color_lines = []
	vrml_lines._colors_offset_table = {}
	vrml_lines._point_lines = []
	vrml_lines._points_offset_table = {}
	vrml_lines._poly_line_color_offset_lines = []
	vrml_lines._poly_line_point_offsets_lines = []

    def _close(self):
        """ *VRML_Lines*: Close the *VRML_Lines* to further polyline additions:
	"""

	# Use *vrml_lines* instead of *self*:
        vrml_lines = self

	# Grab *lines* from *VRML* super-class object:
	lines = vrml_lines._lines
	
	# Start the Shape and geometry:
	lines.append("### VRML_Lines\n")
	lines.append("Shape {\n")
	lines.append(" geometry IndexedLineSet {\n")

	# Append the "color" section to *lines*:
	lines.append("  colorPerVertex FALSE\n")
	lines.append("  color Color {\n")
	lines.append("   color [\n")
	color_lines = vrml_lines._color_lines
	for color_line in color_lines:
	    lines.append(color_line)
	lines.append("   ]\n")
	lines.append("  }\n")

	# Append the "coord" section to *lines:
	lines.append("  coord Coordinate {\n")
	lines.append("   point [\n")
	point_lines = vrml_lines._point_lines
	for point_line in point_lines:
	    lines.append(point_line)
	lines.append("   ]\n")
	lines.append("  }\n")
	
	# Append the "coordIndex" section to *lines*:
	poly_line_point_offsets_lines= vrml_lines._poly_line_point_offsets_lines
	lines.append("  coordIndex [\n")
	for poly_line_point_offsets_line in poly_line_point_offsets_lines:
	    lines.append(poly_line_point_offsets_line)
	lines.append("  ]\n")

	# Append "colorIndex" section to *lines*:
	lines.append("  colorIndex [\n")
	poly_line_color_offset_lines = vrml_lines._poly_line_color_offset_lines
	for poly_line_color_offset_line in poly_line_color_offset_lines:
	    lines.append(poly_line_color_offset_line)
	lines.append("  ]\n")

	# Wrap up the Shape and geometry:
	lines.append(" }\n")
	lines.append("}\n")

	# Mark that we are closed for additional polylines:
	vrml_lines._is_open = False

    def _box_outline(self, color_name, corner1, corner2):
	""" *VRML_Lines*: Draw a box outline that touches *corner1* and *corner2* using the
	    *VRML_Lines* object (i.e. *self*).  The lines will be drawn in *color_name*.
	"""

	# Use *vrml_lines* instead of *self*:
	vrml_lines = self

	# Verify argument types:
	assert isinstance(color_name, str)
	assert isinstance(corner1, P)
	assert isinstance(corner2, P)

	# Make sure *vrml_lines* is still open:
	assert vrml_lines._is_open

	# Figure out the X/Y/Z coordinates:
	x1, x2 = corner1.x.minimum_maximum(corner2.x)
	y1, y2 = corner1.y.minimum_maximum(corner2.y)
	z1, z2 = corner1.z.minimum_maximum(corner2.z)

	# Create the 8 box corners:
	p1 = P(x1, y1, z1)	# BSW
	p2 = P(x2, y1, z1)	# BSE
	p3 = P(x2, y2, z1)	# BNE
	p4 = P(x1, y2, z1)	# BNW
	p5 = P(x1, y2, z2)	# TNW
	p6 = P(x2, y2, z2)	# TNE
	p7 = P(x2, y1, z2)	# TSE
	p8 = P(x1, y1, z2)	# TSW

	# Draw the box lines:
	vrml_lines._poly_line(color_name, [p1, p2, p3, p4, p5, p6, p7, p8, p1, p4])
	vrml_lines._poly_line(color_name, [p5, p8])
	vrml_lines._poly_line(color_name, [p2, p7])
	vrml_lines._poly_line(color_name, [p3, p6])

    def _poly_line(self, color_name, points, tracing=-1000000):
	""" *VRML_Lines*: Draw a *color* poly line between *points* in the *VRML_Lines* object
	    (i.e. *self*).
	"""

	# Use *vrml_lines* instead of *self*:
	vrml_lines = self
	 
	# Verify argument types:
	assert isinstance(color_name, str) and color_name != ""
	assert isinstance(points, list)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>VRML_Lines._poly_line('{1}', '{2}', *)".
	      format(indent, vrml_lines._name, color_name))

	# Make sure that *vrml_lines* is still open:
	assert vrml_lines._is_open

	# Figure out the *color_index* for *color*; create one if not there already:
	color_lines = vrml_lines._color_lines
	colors_offset_table = vrml_lines._colors_offset_table
	if color_name in colors_offset_table:
	    # Color was previously created; reuse it:
	    color_offset = colors_offset_table[color_name]
	else:
	    # Color not previously created: create it:
	    color = Color(color_name)
	    red, green, blue = color._red_green_blue_get()
	    color_offset = len(color_lines)
	    colors_offset_table[color_name] = color_offset
	    color_lines.append("    {0} {1} {2} # {3}\n".format(red, green, blue, color_name))
	assert 0 <= color_offset < len(color_lines)

	# This is the padding needed for *poly_line_point_offsets_lines* below:
	line_offsets = [ "  " ]

	# Make sure we have a *coordinate_index* for each *point* in *points*:
	point_lines = vrml_lines._point_lines
	points_offset_table = vrml_lines._points_offset_table
	for index, point in enumerate(points):
	    # Create a *coordinate* tuple key using imutuable *float* types:
	    assert isinstance(point, P), "points[{0}] is not a Point".format(index)
	    point_tuple = point.triple()

	    # Make sure we have a *coordinate_index*:
	    if point_tuple in points_offset_table:
		# The coordinate has already been entered; reuse it: 
		point_offset = points_offset_table[point_tuple]
	    else:
		# The coordinate is new and must be created:
		point_offset = len(point_lines)
		points_offset_table[point_tuple] = point_offset
		point_lines.append("    {0} {1} {2}\n".
		  format(point_tuple[0], point_tuple[1], point_tuple[2]))

	    # Now we can append *coordinate_index* to *line_offsets*:
	    line_offsets.append(" {0}".format(point_offset))

	# Indicate the end of the poly-line with -1:
	line_offsets.append(" -1\n")
	poly_line_point_offsets_line = "".join(line_offsets)
	poly_line_point_offsets_lines = vrml_lines._poly_line_point_offsets_lines
	poly_line_point_offsets_lines.append(poly_line_point_offsets_line)

	# Associate the *color_index* for *color_name* with this poly-line:
	poly_line_color_offset_line = "  {0}\n".format(color_offset)
	poly_line_color_offset_lines = vrml_lines._poly_line_color_offset_lines
	poly_line_color_offset_lines.append(poly_line_color_offset_line)

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=VRML_Lines._poly_line('{1}', '{2}', *)".
	      format(indent, vrml_lines._name, color_name))

    def _text_pad(self, pad, tracing=-1000000):
	""" *VRML_Lines*:
	"""

	# Use *vrml_lines* inestead of *self*:
	vrml_lines = self

	# Verify argument types:
	assert isinstance(pad, int)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>VRML_Lines._text_pad('{1}', {2})".format(indent, vrml_lines._name, pad))

	# Convert *vrml_lines* into one big string of *text*:
	text = vrml_lines._text_pad_helper(pad)

	# Wrap up any requested *tracing* and return *text*:
	if tracing >= 0:
	    print("{0}<=VRML_Lines._text_pad('{1}', {2})".format(indent, vrml_lines._name, pad))
	return text


class VRML_Rotate(VRML):
    """ *VRML_Rotate*: Represents a VRML Rotatiation Transform node.
    """

    def __init__(self, name, axis, amount, child_vrml, tracing=-1000000):
	""" *VRML_Rotate*: Initialize the *VRML_Rotate* node object (i.e. *self*) to
	    represent a rotation of *amount* around *axis*.
	"""
	
	# Use *vrml_rotate* instead of *self*:
	vrml_rotate = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(axis, P)
	assert isinstance(amount, Angle)
	assert isinstance(child_vrml, VRML)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>VRML_Rotate.__init__(*, '{1}', {2}, {3:d}, '{4}')".
	      format(indent, name, axis, amount, child_vrml._name))

	# Intialize the super class:
	VRML.__init__(vrml_rotate, name)

	# Load values into *vrml_rotate*:
	vrml_rotate._axis = axis.normalize()
	vrml_rotate._amount = amount
	vrml_rotate._child_vrml = child_vrml
	vrml_rotate._is_open = False

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=VRML_Rotate.__init__(*, '{1}', {2}, {3:d}, '{4}')".
	      format(indent, name, axis, amount, child_vrml._name))

    def _text_pad(self, pad, tracing=-1000000):
	""" *VRML_Rotate*: Return the *VRML_Rotate* object (i.e. *self*) as a single text
	    string where each line is indented by *pad* spaces.
	"""

	# Use *vrml_rotate* instead of *self*:
	vrml_rotate = self
	
	# Verify argument types:
	assert isinstance(pad, int)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>VRML_Rotate._text_pad('{1}', {2})".format(indent, vrml_rotate._name, pad))

	assert isinstance(vrml_rotate._name, str)

	# Figure out if we have already done this *pad* before:
	pad_table = vrml_rotate._pad_table
	if pad in pad_table:
	    text = pad_table[pad]
	else:
	    # Grab values out of *vrml_rotate*:
	    axis = vrml_rotate._axis
	    axis_x, axis_y, axis_z = axis.triple()
	    amount = vrml_rotate._amount
	    angle = amount.radians()
	    lines = vrml_rotate._lines

	    # Construct all the *lines* that are needed for a VRML rotate:
	    lines.append("### VRML_Rotate\n")
	    lines.append("Transform {\n")
	    lines.append(" rotation {0} {1} {2} {3}\n".format(axis_x, axis_y, axis_z, angle))
	    lines.append(" children [\n")
	    lines.append(   vrml_rotate._child_vrml._text_pad(pad + 2))
	    lines.append(" ]\n")
	    lines.append("}\n")
	    
	    # Convert *lines* into one big *text* string where each lines is prepended by *padding*:
	    padding = ' ' * pad
	    text = padding.join(lines)

	    # Stuff *text* back into *pad_table*:
	    pad_table[pad] = text

	# Wrap up any requested *tracing* and return *text*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=VRML_Rotate._text_pad('{1}', {2})".format(indent, vrml_rotate._name, pad))
	return text

class VRML_Switch(VRML):
    """ *VRML_Switch*: Represents a VRML Switch node.
    """

    def __init__(self, name, which_choice):
	""" *VRML_Switch*: Initialize the *VRML_Switch* object (i.e. *self*):
	"""

	# Use *vrml_switch* instead of *self*:
	vrml_switch = self

	# Verify argument types:
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(which_choice, int)

	# Initialize super-class:
	VRML.__init__(self, name)

	# Create *children* lists:
	vrml_switch._children = []
	vrml_switch._which_choice = which_choice

    def _append(self, vrml):
	""" *VRML_Switch*: Append *vrml* to the *VRML_Switch* object (i.e. *self*):
	"""

	# Use *vrml_switch* instead of *self*:
	vrml_switch = self

	# Verify argument types:
	assert isinstance(vrml, VRML)

	# Make sure that *vrml_switch* is still open:
	assert vrml_switch._is_open

	# Append *vrml* to *children*:
	vrml_switch._children.append(vrml)

    def _text_pad(self, pad, tracing=-1000000):
	""" *VRML_Switch*: Return the *VRML_Switch* object (i.e. *self*) as a single string 
	"""

	# Use *vrml_switch* instead of *self*:
	vrml_switch = self

	# Verify argument types:
	assert isinstance(pad, int)	
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>VRML_Switch._text_pad('{1}', {2})".format(indent, vrml_switch._name, pad))
	    trace_detail = 2

	# Mark this switch as *closed*:
	vrml_switch._is_open = False

	# Construct the content:
	pad_table = vrml_switch._pad_table
	if pad in pad_table:
	    # This padding level has been previously computed:
	    text = pad_table[pad]
	else:
	    # Put all the needed text pieces into *lines*:
	    lines = []
	    lines.append("\n")
	    lines.append("### VRML_Switch\n")
	    lines.append("Switch {\n")
	    which_choice = vrml_switch._which_choice
	    if which_choice >= 0:
		lines.append(" whichChoice {0}\n".format(which_choice))
	    lines.append(" choice [\n")
	    for index, vrml_child in enumerate(vrml_switch._children):
		if trace_detail >= 2:
		    print("{0}[{1}] '{2}'".format(indent, index, vrml_child.__class__.__name__))
		vrml_child_text = vrml_child._text_pad(pad + 2, tracing = tracing + 1)
		assert isinstance(vrml_child_text, str), \
		  "Bad return for {0}".format(vrml_child.__class__.__name__)
		lines.append(vrml_child_text)
	    lines.append(" ]\n")
	    lines.append("}\n")

	    # Concatenate all of lines into one big string of *text* padded with *padding*::
	    padding = ' ' * pad
	    text = padding.join(lines)

	    # Stuff *text* into *pad_table*:
	    pad_table[pad] = text

	# Wrap up any requested *tracing* and return *text*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}<=VRML_Switch._text_pad('{1}', {2})".format(indent, vrml_group._name, pad))
	return text

class VRML_Translate(VRML):
    """ *VRML_Translate*: Represents a VRML Rotatiation Transform node.
    """

    def __init__(self, name, translate, child_vrml, tracing=-1000000):
	""" *VRML_Translate*: Initialize the *VRML_Translate* node object (i.e. *self*) to
	    represent a translation by *translate* for *child_vrml*.
	"""
	
	# Use *vrml_translate* instead of *self*:
	vrml_translate = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(translate, P)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>VRML_Translate.__init__(*, '{1}', {2:i}, '{3}')".
	      format(indent, name, translate, child_vrml._name))

	# Intialize the *VRML* super class:
	VRML.__init__(vrml_translate, name)

	# Load values into *vrml_translate*:
	vrml_translate._child_vrml = child_vrml
	vrml_translate._is_open = False
	vrml_translate._translate = translate

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=VRML_Translate.__init__(*, '{1}', {2:i}, '{3}')".
	      format(indent, name, translate, child_vrml._name))

    def _text_pad(self, pad, tracing=-1000000):
	""" *VRML_Translate* Return the *VRM_Translate* object (i.e. *self*) as a string
	    where each line is preceeded by *pad* spaces.
	"""

	# Use *vrml_translate* instead of *self*:
	vrml_translate = self
	# Verify argument types:
	assert isinstance(pad, int)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}VRML_Translate._text_pad('{1}', {2})".
	      format(indent, vrml_translate._name, pad))
	    trace_detail = 2

	# Figure out if we have already done this *pad* before:
	pad_table = vrml_translate._pad_table
	if pad in pad_table:
	    text = pad_table[pad]
	else:
	    # Grab the values out of *vrml_translate*:
	    translate = vrml_translate._translate
	    dx, dy, dz = translate.triple()
	    if trace_detail >= 2:
		print("{0}dx={1} dy={2} dz={3} translate={4:i}in={4:m}mm".
		  format(indent, dx, dy, dz, translate))

	    # Construct all the *lines* that are needed for a VRML Transform:
	    lines = []
	    lines.append("\n")
	    lines.append("# {0}\n".format(vrml_translate._name))
	    lines.append("### VRML_Translate\n")
	    lines.append("Transform {\n")
	    lines.append(" translation {0} {1} {2}\n".format(dx, dy, dz))
	    lines.append(" children [\n")
	    lines.append(   vrml_translate._child_vrml._text_pad(pad + 2))
	    lines.append(" ]\n")
	    lines.append("}\n")
	    
	    # Convert *lines* into one big *text* string where each lines is prepended by *padding*:
	    padding = ' ' * pad
	    text = padding.join(lines)

	    # Stuff *text* back into *pad_table*:
	    pad_table[pad] = text

	# Wrap up any requested *tracing* and return text:
	if tracing >= 0:
	    print("{0}VRML_Translate._text_pad('{1}', {2})".
	      format(indent, vrml_translate._name, pad))
	return text

class VRML_Triangles(VRML):
    """ *VRML_Triangles*: Represent a VRML shape that consists of colored triangles.
    """

    def __init__(self, name, color_name, triangles, tracing=-1000000):
	""" *VRML_Triangles*: Initialize the *VRML_Triangles* object (i.e. *self*) so that 
	    will render a list of *triangles* in *color*.
	"""

	# Use *vrml_triangles* instead of *self*:
	vrml_triangles = self

	# Verify argument types:
	assert isinstance(name, str)
	assert isinstance(color_name, str)
	assert isinstance(triangles, list) or isinstance(triangles, tuple)
	assert isinstance(tracing, int)

	# Perform any *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>VRML_Triangles.__init__(*, '{1}, '{2}', *)".
	      format(indent, name, color_name))

	# Initialize the super class of *vrml_triangles*:
	VRML.__init__(vrml_triangles, name)

	# In VRML, the triangles are represented as a list of triangle vertices (i.e. a triangle
	# corner) and then triples that index into this array.

	# *corners_table* and *corners* are used to generate a unique indiex for each triangle
	# corner.  *corner_offset_triangle* is collects 3 corner indices together before appended
        # to *corner_offset_triangles* as an immutable tuple.
	corners_table = {}
	corners = []
	corner_offset_triangle = []
	corner_offset_triangles = []
	for triangle in triangles:
	    assert isinstance(triangle, tuple)

	    # Make sure *corner_offset_triangle* is empty:
	    del corner_offset_triangle[:]
	    for corner in triangle:
		assert isinstance(corner, tuple)

		# If *corner* has occurred previously, it will already be in *corners_table*:
		if corner in corners_table:
		    # Reuse the previously allocated *corner_offset*:
		    corner_offset = corners_table[corner]
		else:
		    # Allocate a new *corner_offset* for *corner* an remember it in *corners_table*:
		    corner_offset = len(corners)
		    corners_table[corner] = corner_offset
		    corners.append(corner)

		# Collect three *corner_offset*'s in *corner_offset_triangle*:
		corner_offset_triangle.append(corner_offset)

	    # Record *corner_offset_triangle* into *corner_offset_triangles* as an immutable tuple:
	    corner_offset_triangles.append(tuple(corner_offset_triangle))

	# Mark *vrml_triangles* as closed from the start:
	vrml_triangles._is_open = False

	# Grab the super class *lines* list:
	lines = vrml_triangles._lines

	# Append the Shape and appearance VRML to *lines*:
	lines.append("### VRML_Triangles\n")
	lines.append("Shape {\n")
	lines.append(" appearance Appearance {\n")
	lines.append("  material Material {\n")
	color = Color(color_name)
	red, green, blue = color._red_green_blue_get()
	lines.append("   diffuseColor {0} {1} {2}\n".format(red, green, blue))
	alpha = color._alpha_get()
	if alpha < 1.0:
	    lines.append("   transparency {0}\n".format(1.0 - alpha))
	lines.append("  }\n")
	lines.append(" }\n")

	# Now start the geometry portion of *lines*.  It consists of point coordinates
	lines.append(" geometry IndexedFaceSet {\n")

	# Append the *corners* to *lines*:
	lines.append("  coord Coordinate {\n")
	lines.append("   point [\n")
	for corner in corners:
	    lines.append("    {0} {1} {2}\n".format(corner[0], corner[1], corner[2]))
	lines.append("   ]\n")
	lines.append("  }\n")

	# Append the *corner_offset_triangles* to *lines*:
	lines.append("  coordIndex [\n")
	for corner_offset_triangle in corner_offset_triangles:
	    lines.append("   {0} {1} {2} -1\n".format(corner_offset_triangle[0],
	      corner_offset_triangle[1], corner_offset_triangle[2]))
	lines.append("  ]\n")

	# Wrap up the VRML *shape*:
	lines.append(" }\n")
	lines.append("}\n")

	# Wrap up any *tracing*:
	if tracing >= 0:
	    print("{0}<=VRML_Triangles.__init__(*, '{1}', '{2}', *)".
	      format(indent, name, color_name))

    def _text_pad(self, pad, tracing=-1000000):
	""" *VRML_Triangles*: Return the *VRML_Triangles* object (i.e. *self*) as single
	    string where each line is preceeded by *pad* spaces.
	"""

	# Use *vrml_triangles* instead of *self*:
	vrml_triangles = self
	
	# Verify argument types:
	assert isinstance(pad, int)
	assert isinstance(tracing, int)

	#tracing = 0

	# Perform an requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * pad
	    print("{0}VRML_Triangles._text_pad('{1}', {2})".
	      format(indent, vrml_triangles._name, pad))

	# Convert *vrml_triangles* to padded *text*:
	text = vrml_triangles._text_pad_helper(pad, tracing = tracing + 1)

	# Wrap up any requested *tracing* and return text:
	if tracing >= 0:
	    print("{0}VRML_Triangles._text_pad('{1}', {2})".
	      format(indent, vrml_triangles._name, pad))
	return text

class VRML_Use(VRML):
    """ *VRML_Use*: Represents a VRML Rotatiation Transform node.
    """

    def __init__(self, name, use, tracing=-1000000):
	""" *VRML_Use*: Initialize the *VRML_Use* node object (i.e. *self*) to
	    represent a VRML "use" call.
	"""
	
	# Use *vrml_use* instead of *self*:
	vrml_use = self

	# Verify argument types:
	assert isinstance(name, str) and not ' ' in name
	assert isinstance(use, Part)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}=>VRML_Use.__init__(*, '{1}', '{2}')".
	      format(indent, name, use._name_get()))

	# Intialize the *VRML* super class:
	VRML.__init__(vrml_use, name)

	# Load values into *vrml_use*:
	vrml_use._use = use

	# Wrap up any requested *tracing*:
	if tracing >= 0:
	    print("{0}<=VRML_Use.__init__(*, '{1}', '{2}')".
	      format(indent, name, use._name_get()))

    def _text_pad(self, pad, tracing=-1000000):
	""" *VRML_Use* Return the *VRM_Use* object (i.e. *self*) as a string
	    where each line is preceeded by *pad* spaces.
	"""

	# Use *vrml_use* instead of *self*:
	vrml_use = self

	# Verify argument types:
	assert isinstance(pad, int)
	assert isinstance(tracing, int)

	# Perform any requested *tracing*:
	trace_detail = -1
	if tracing >= 0:
	    indent = ' ' * tracing
	    print("{0}VRML_Use._text_pad('{1}', {2})".
	      format(indent, vrml_use._name, pad))
	    trace_detail = 2

	# Figure out if we have already done this *pad* before:
	pad_table = vrml_use._pad_table
	if pad in pad_table:
	    text = pad_table[pad]
	else:
	    # Grab the values out of *vrml_use*:
	    use = vrml_use._use
	    use_name = use._name_get()

	    # Construct all the *lines* that are needed for a VRML Transform:
	    lines = []
	    lines.append("\n")
	    lines.append("### VRML_Use\n")
	    lines.append("# {0}\n".format(use_name))
	    lines.append("USE {0}\n".format(use_name))
	    
	    # Convert *lines* into one big *text* string where each lines is prepended by *padding*:
	    padding = ' ' * pad
	    text = padding.join(lines)

	    # Stuff *text* back into *pad_table*:
	    pad_table[pad] = text

	# Wrap up any requested *tracing* and return text:
	if tracing >= 0:
	    print("{0}VRML_Use._text_pad('{1}', {2})".
	      format(indent, vrml_use._name, pad))
	return text

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
#	    wrl_file_name = os.path.join(ngc_directory, "{0}.wrl".format(program_number))
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
