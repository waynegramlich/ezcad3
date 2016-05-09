#!/usr/bin/env python

from EZCAD3 import *

def pudb_stringifier(object):
    result = type(object)
    try:
	if isinstance(object, Bounding_Box):
	    result = "BB:{0}".format(object)
	elif isinstance(object, L):
	    result = "L:{0}".format(object)
	elif isinstance(object, P):
	    result = "P:{0}".format(object)
	elif isinstance(object, Angle):
	    result = "Angle:{0:d}"
	elif isinstance(object, Material):
	    result = "Material:{0}".format(object)
	elif isinstance(object, Color):
	    result = "Color:{0}".format(object)
	elif isinstance(object, Part):
	    result = "Color:{0}".format(object)
    except:
	result = type(object)
    return result

