#!/usr/bin/env python

# This is a bunch of unit tests for the *Transform* class.
# It basically verifies that the correct matrices and openscad commands
# are genrated by the *Transform* class.  The matrices are verified with
# assert statements, whereas the openscad lines are verified by viewing
# `translate_rotate.scad` and `rotate_translate.scad`.  The blue cube is
# where the "point" starts (i.e. the location (1, 0, 0)" and the green
# cube is where the point is transformed to>

from EZCAD3 import *

def test_transform():
    print("=>test_transform()")
    # Define some constants
    zero = L()
    one = L(mm=1.0)
    two = L(mm=2.0)
    origin = P(zero, zero, zero)
    x_axis = P(one,  zero, zero)
    y_axis = P(zero, one,  zero)
    z_axis = P(zero, zero, one)
    degrees0 = Angle()
    degrees90 = Angle(deg=90.0)
    
    # Create the identity *Transform*:
    identity = Transform()

    # Check the formatting routines:
    matrix_text = "[[1. 0. 0. 0.] [0. 1. 0. 0.] [0. 0. 1. 0.] [0. 0. 0. 1.]]"
    assert "{0:m}".format(identity) == matrix_text
    assert "{0:s}".format(identity) == "()"

    # Check that equality works:
    assert identity == identity
    assert not identity != identity

    # Do a null translate:
    translate_null = identity.translate("null translate", origin)
    assert translate_null == identity

    # Do a null rotate:
    rotate_null = identity.rotate("null_rotate", z_axis, degrees0)
    assert rotate_null == identity

    # Do a simple *translate* by 1 along the *x_axis*:
    translate = identity.translate("", x_axis)
    matrix_text = "[[1. 0. 0. 0.] [0. 1. 0. 0.] [0. 0. 1. 0.] [1. 0. 0. 1.]]"
    assert "{0:s}".format(translate) == "('translate([1.0, 0.0, 0.0]) //F: ',)"
    assert translate * origin == x_axis
    assert translate * x_axis == P(two, zero, zero)

    # Now make sure we can get back from *translate*:
    translate_reverse = translate.reverse()
    assert translate_reverse * P(two, zero, zero) == x_axis
    assert translate_reverse * x_axis == origin
    assert translate_reverse * P(two, zero, zero) == x_axis
    assert "{0:s}".format(translate_reverse) == "('translate([-1.0, 0.0, 0.0]) //R: ',)"

    # Do a simple rotate around the Z axis:
    rotate = identity.rotate("", z_axis, degrees90)
    matrix_text = "[[0. 1. 0. 0.] [-1. 0. 0. 0.] [0. 0. 1. 0.] [0. 0. 0. 1.]]"
    assert "{0:m}".format(rotate) == matrix_text
    assert "{0:s}".format(rotate) == "('rotate(a=90.0, v=[0.0, 0.0, 1.0]) //F: ',)"
    assert rotate * x_axis == y_axis

    # Now make sure we can get back from *rotate_z_axis*:
    rotate_reverse = rotate.reverse()
    assert rotate_reverse * y_axis == x_axis
    assert "{0:s}".format(rotate_reverse) == "('rotate(a=-90.0, v=[0.0, 0.0, 1.0]) //R: ',)"

    # Do a translate and then a rotate:
    translate_rotate = \
      identity.translate("x_axis translate 2", x_axis).rotate("rotate z_axis 2", z_axis, degrees90)
    #print("translate_rotate={0:m}".format(translate_rotate))
    translate_rotate_p = translate_rotate * x_axis
    #print("translate_rotate_p={0:m}".format(translate_rotate_p))
    assert translate_rotate_p == P(zero, two, zero)

    # Perform the reverse:
    translate_rotate_reverse = translate_rotate.reverse()
    back_p = translate_rotate_reverse * translate_rotate_p
    assert back_p == x_axis
    scad_file_write("translate_rotate.scad", translate_rotate, x_axis)

    # Do a rotate then a translate:
    rotate_translate = \
      identity.rotate("rotate z_axis 3", z_axis, degrees90).translate("x_axis translate 3", x_axis)
    #print("rotate_translate={0:m}".format(translate_rotate))
    rotate_translate_p = rotate_translate * x_axis
    #print("rotate_translate_p={0:m}".format(translate_rotate_p))
    assert rotate_translate_p == P(one, one, zero)

    # Perform the reverse:
    rotate_translate_reverse = rotate_translate.reverse()
    back_p = rotate_translate_reverse * rotate_translate_p
    assert back_p == x_axis
    scad_file_write("rotate_translate.scad", rotate_translate, x_axis)

    top_surface_transform("bottom", -z_axis)
    top_surface_transform("north",   y_axis)
    top_surface_transform("south",  -y_axis)
    top_surface_transform("east",    x_axis)
    top_surface_transform("west",   -x_axis)
    top_surface_transform("top",     z_axis)

    print("<=test_transform()")

def top_surface_transform(comment, end, tracing=-1000000):
    # Verify argument types:
    assert isinstance(comment, str)
    assert isinstance(end, P)

    if tracing >= 0:
	indent = ' ' * tracing
	print("{0}=>top_surface_transform('{1}', {2:m})".format(indent, comment, end))

    # Some constants:
    zero = L()
    one = L(mm=1.00)
    origin = P(zero, zero, zero)
    z_axis = P(zero, zero, one)
    degrees0 = Angle(deg=0.00)
    epsilon = L(mm=.0000001)

    # Perform a z-axis invert:
    transform = Transform.top_surface(comment, origin, end, degrees0, tracing = tracing + 1)
    transform_p = transform * end
    if tracing >= 0:
	print("{0}transform={1:s} transform_p={2:m} -z_axis={3:m}".
	  format(indent, transform, transform_p, -z_axis))
	print("{0}transform={1:m}".format(indent, transform))
    assert transform_p.distance(-z_axis).absolute() < epsilon

    transform_reverse = transform.reverse()
    back_p = transform_reverse * transform_p
    assert back_p == end

    scad_file_write("{0}.scad".format(comment), transform, end)

    if tracing >= 0:
	indent = ' ' * tracing
	print("{0}<=top_surface_transform('{1}', {2:m})".format(indent, comment, end))

def scad_file_write(file_name, transform, point):
    # Verify argument types:
    assert isinstance(file_name, str)
    assert isinstance(transform, Transform)
    assert isinstance(point, P)

    # Extract the X/Y/Z coordinate of *point*:
    x = point.x.millimeters()
    y = point.y.millimeters()
    z = point.z.millimeters()

    # Open *scad_file*:
    scad_file = open(file_name, "w")

    # Generate the forward transform green cube:
    for scad_line in transform._forward_scad_lines:
	scad_file.write(scad_line)
	scad_file.write("\n")

    # The final translate places the *point*:
    scad_file.write("translate([{0}, {1}, {2}])\n".format(x, y, z))
    scad_file.write("color(\"green\", 1)\n")
    scad_file.write("cube([.1, .1, .1], center = true);\n")

    # Openscad does transforms backwards, so we output *reverse_scad_lines* "first":
    for scad_line in transform._reverse_scad_lines:
	scad_file.write(scad_line)
	scad_file.write("\n")

    # Now we output the *forward_scad_lines* "second":
    for scad_line in transform._forward_scad_lines:
	scad_file.write(scad_line)
	scad_file.write("\n")

    # The final translate places the *point* come "third":
    scad_file.write("translate([{0}, {1}, {2}])\n".format(x, y, z))
    scad_file.write("color(\"blue\", 1)\n")
    scad_file.write("cube([.1, .1, .1], center = true);\n")

    scad_file.close()

if __name__ == "__main__":
    test_transform()
