# EZCAD3: Easy (Mechanical) Computer Aided Design -- Version 3

## Introduction

EZCAD stands for Easy Computed Aided Design.  EZCAD is a
mechanical CAD (Computer Aided Design) system with an integrated
mechanical CAM (Computer Aided Manufacture) system.  The CAD
system specifies the geometry of individual parts and how they
are assembled together.  The CAM system specifies the
manufacturing steps required to make a part.

In EZCAD, the designer writes a program using the
(Python)[http://www.python.org/] programming language.
The python program is then executed to generate both the
design visualization and all associated manufacturing files
(e.g. .stl files for 3-printers, G-codes for CNC mills/lathes,
.dxf files for laser printers, etc.)  EZCAD does not include
either a text editor or a visualization tool.  Instead, the designer
uses their preferred text editor (e.g. idle, emacs, vim, etc.) and
they use their own preferred visualization tools  (e.g. meshlab,
etc.)  Other tools are used to visualize generated manufacturing
files.

In short, EZCAD is just a Python library that is executed to generate
the various design files.

## Downloading

This package is distributed in source form from a *git* repository.
To get a copy:

        git clone https://github.com/waynegramlich/ezcad3.git

In addition, this package needs a version of
(OpenSCAD)[http://www.openscad.org/].  For those who use
Linux distributions based on the Debian package system
(e.g. Ubuntu, Mint, etc.), the following will do the trick:

        sudo apt-get install openscad

OpenSCAD is primarily used as a front-end to the wonderful
(CGAL)[http://www.cgal.org/] (Computational Geometry Algorithms
Library).  In particular, OpenSCAD uses the Boolean 3D solids
library embodied in the NEF3 sub-library of CGAL.

Detailed documentation of the Python classes and associated
methods is maintained by (Doxygen)[http://www.doxygen.org/].
For convenience, all of the generated Doxygen documenation
is checked into the repository.  If you want to generated
Doxygen documentation locally, you need to download the program.
I downloaded my copy via:

        sudo apt-get install doxygen

## Licensing

In general, I really like Open Source Licenses.  I have a slight
preference of the GPL open source license, so that is what EZCAD3
is released under.  I emulate the
(Free Software Foundation)[http://www.fsf.org] that any code
contributed to my code branch requires a copyright assignment.
Professor Eben Moglen has written up a short explanation of
(Why the FSF gets copyright assignments from
contributors)[http://www.gnu.org/licenses/why-assign.html].
I strongly feel that having a unified copyright owner provides
maximum protection for the open source software.  While it is
possible to  get
(overly legalistic)[http://ftp.xemacs.org/old-beta/FSF/assign.changes],
I think the following is more that adequate to assign the copyright:

        Hello:

        I am assigning my modifications to EZCAD to
        Wayne C. Gramlich with the requirement that the code
        will continue to be released under the GPL version 3
        license or higher.

        Regards,

        {Your name here}

Enough on this legal stuff!

## Quick Overview

In EZCAD, the most important Python class is the *Part* class.  The
*Part* class is used to both manufacture individual parts *and* to
assemble the individual parts together into sub-assemblies.

In addition to the *Part* class there are bunch of utility classes.
In alphabetical order, the utility classes are:

* *Angle* is used to represent angles.

* *Color* is used to represent a represent the color that an individual
  part is rendered as.

* *L* is used to represent a length  (e.g. mm, cm, inch, ft, etc.)
  (This type is shortened to "L" because is used so much.)

* *Material* specifies what material a part is manufactured out of.

* *P* is used to represent a point in a 3-dimensional Cartesian
   coordinate system.  (This type is shortened to "P" because is used
   so much.)

* *Place* is used to represent a *Part* placement in a sub-assembly.

Each EZCAD class method will validates that its arguments types
are correct.  (If you have ever heard of the "duck typing" design
pattern, EZCAD most definitely does **not** use that particular
programming religion.)

The CAM (Computer Aided Manufacture) portion is based on the
following classes:

* *EZCAD* is the global object that orchestrates the entire process.

* *Shop* is the parent class that specifies all the available
  machines and associated machine tooling.

* *Machine* is a specific machine such as a 3D-printer, laser-cutter,
  mill, lathe, drill press, etc.

* *Tool* is a specific tool that can be loaded into a machine.
  These include drill bits, end-mills, etc.

The shop is specified independently from the part and sub-assembly
file.  The same part will generate different manufacturing files
depending upon what machines and tooling is available in a given shop.
Thus, different G-codes, .dxf files are generated depending upon
what is available in the shop.

EZCAD operates in three distinct phases:

* Part Initialization.  This phase creates a hierarchical name space
  of parts and sub-assemblies.

* Constraint Propagation.  This phase propagates dimensional changes
  around the overall design.  This process also keeps track of the
  bounding box around each individual part and sub-assembly.

* Manufacture and Assembly.  This phase generates the required
  manufacturing files needed for construct the various parts.
  This phase also specifies where *Part*'s are placed in each
  sub-assembly.

Each phase is quickly discussed below:

### Part Initialization Phase ###

Each user individual part is implemented as a Python class object that
sub-classes the *Part* super class.  As expected, the initialization phase
occurs in the *\_\_init\_\_*() method.  Each part has an up-level parent
(abbreviated as "*up*").  The following code fragment shows how it is
done:

        class My_Part(Part):

            def __init__(self, up):
                Part.__init__(self, up)	     # Initialize *Part* super-class
                self.part1_ = Part1(self)    # Initialize *Part1*
                # ...
                self.partN_ = PartN(self)    # Initialize *PartN*

The first line initializes the *Part* super class with *up* as the
up-level parent.  The second and subsequent lines initializes the
various sub-parts.  By convention, each sub-*Part* is stored into
*self* with a trailing "_" suffix.  (There will be more about suffixes
shortly.)

After initialization is complete, each *Part* can access any other
*Part* in the design using standard dot (".") notation.  For example,
*self.up.base_* will access the *base_* *Part* of the parent to *self*.

Due to the nature of Python, the designer will be storing dimensional
values into *Part* member variables.  EZCAD uses member variable suffixes
to prevent accidental name clashes.  Designers are constrained to use
a suffixes listed below:

* "_": *Part*
* "_a": *Angle*
* "_c": *Color*p
* "_i": *int* (i.e. an integer number)
* "_f": *float* (i.e. a floating point number)
* "_l": *L* (i.e. a length)
* "_m": *Material*
* "_o": Some other type
* "_s": *str* (i.e. string)
* "_p": *P* (i.e. a point)
* "_pl": *Place*

All user defined member variables *must* end with one of the suffixes above.

### Constraint Propagation Phase ###

Constraint propagation is what allows the designer to relatively easily
resize a design.  The constraint propagation phase occurs in a
*Part*'s *construct*() method.  Constraint propagation occurs by
repeatably calling the *construct*() method for each *Part* in
the entire design.  The *construct*() methods are repeatably called
until none of the values with approved suffixes (e.g. "\_l", "\_p")
change any more.

To further support constraint propagation, each *Part* maintains a
bounding box for everything it contains.  The designer can easily
access the various bounding box values using short concise member
variable names.  Along the X direction, the two bounding box values
are referred as east ("e") and west ("w").  For the Y direction,
north ("n") and south ("s") are used.  For the Z direction, top ("t")
and bottom ("b") are used.  The bounding box results in 12 length
(i.e. *L*) member variables and 27 point (i.e. *P*) member variables
as listed below:

* *ex*:  East   X length
* *cx*:  Center X length
* *wx*:  West   X length
* *dx*:  Delta  X length (i.e. *ex* - *wx*)
* *ny*:  North  Y length
* *cy*:  Center Y length
* *sy*:  South  Y length
* *dy*:  Delta  Y length (i.e *ny* - *sy*)
* *tz*:  Top    Z length
* *cz*:  Center Z length
* *bz*:  Bottom Z length
* *dz*:  Delta  Z length (i.e. *tz* - *bz*)
* *tne*: Top/North/East    point (i.e. P(ex, ny, tz))
* *tn*:  Top/North         point (i.e. P(cx, ny, tz))
* *tnw*: Top/North/West    point (i.e. P(wx, ny, tz))
* *te*:  Top/East          point (i.e. P(ex, cy, cz))
* *t*:   Top               point (i.e. P(cx, cy, tz))
* *tw*:  Top/West          point (i.e. P(wx, cy, tz))
* *tse*: Top/South/East    point (i.e. P(ex, sy, tz))
* *ts*:  Top/South/East    point (i.e. P(cx, sy, tz))
* *tsw*: Top/South/East    point (i.e. P(wx, sy, tz))
* *ne*:  North/East        point (i.e. P(ex, ny, cz))
* *n*:   North             point (i.e. P(cx, ny, cz))
* *nw*:  North/West        point (i.e. P(wx, ny, cz))
* *e*:   East              point (i.e. P(ex, cy, cz))
* *c*:   Center            point (i.e. P(cx, cy, cz))
* *w*:   West              point (i.e. P(wx, cy, cz))
* *se*:  North/East        point (i.e. P(ex, sy, cz))
* *s*:   North/East        point (i.e. P(cx, sy, cz))
* *sw*:  North/East        point (i.e. P(wx, sy, cz))
* *bne*: Bottom/North/East point (i.e. P(ex, ny, bz))
* *bn*:  Bottom/North      point (i.e. P(cx, ny, bz))
* *bnw*: Bottom/North/West point (i.e. P(wx, ny, bz))
* *be*:  Bottom/East       point (i.e. P(ex, cy, bz))
* *b*:   South             point (i.e. P(cx, cy, bz))
* *bw*:  Bottom/West       point (i.e. P(wx, cy, bz))
* *bse*: Bottom/South/East point (i.e. P(ex, sy, bz))
* *bs*:  Bottom/South/East point (i.e. P(cx, sy, bz))
* *bsw*: Bottom/South/East point (i.e. P(wx, sy, bz))

{Talk about design in assembly space vs individual parts.}

{Short example of constraint propagation goes here.}

### Manufacture and Assembly Phase ###

The manufacturing and assembly phase occurs last.

{More goes here...}

## A Quick Example ##

The example below creates a plastic box with a cover, where the
cover has a lip that keeps the cover centered over the box.
We present the entire code body first then will describe
what is happening on a chunk by chunk basis next:

        #!/usr/bin/env python

        import EZCAD3   # The EZCAD (revision 3) classes:

        class Simple_Box(Part):

            def __init__(self, up, dx=L(mm=100.0), dy=L(mm=50.0),
              dz=L(25.0), wall_thickness=L(mm=10.0),
              material=Material("plastic", "ABS")):
                # Initialize the *Part*:
                Part.__init__(self, up)

                # Remember the initialization values:
                self.dx_l = dx
                self.dy_l = dy
                self.dz_l = dz
                self.wall_thickness_l_ = wall_thickness
                self.material_m = material

                # Instantiate the sub-*Part*'s:
                self.base_ = Simple_Box_Base(self)
                self.cover_ = Simple_Box_Cover(self)

            def construct(self):
                pass

        class Simple_Box_Base(Part):

            def __init__(self, up, place = True):
                Part.__init__(self, up, place)

            def construct(self):
                # Grab some values from *self*:
                box = self.up
                dx = box.dx_l
                dy = box.dy_l
                dz = box.dz_l
                thickness = box.wall_thickness
                material = box.material_m

                # Add another 
                self.height_l = height = dz - thickness
                zero = L()

                # Start with a solid block of the right dimensions:
                height = dz - wall_thickness
                self.block(comment = "Initial block of material",
                  material = material,
                  color = Color("blue"),
                  corner1 = P(-dx/2, -dy/2, zero),
                  corner2 = P( dx/2,  dy/2, height))

                # Pocket out the body of the box:
                self.simple_pocket(comment = "Box Pocket",
                  corner1 = box.bsw + P(thickness, thickness, thickness),
                  corner2 = box.tne - P(thickness, thickness, zero))

        class Simple_Box_Cover(Part):

            def __init__(self, up, place = True):
                Part.__init__(self, up, place)

            def construct(self):
                # Grab some values from *parent* and *base*:
                box = self.up
                dx = box.dx_l
                dy = box.dy_l
                dz = box.dz_l
                material = box.material_m
                thickness = box.wall_thickness_l
                base = box.base_
                base_height = base.height_l

                # Compute local values:
                self.lip_thickness = lip_thickness = thickness/2

                # Do the top part of the cover:
                self.block(comment = "Cover Top",
                  material = material
                  color = Color("green"),
                  corner1 = base.tsw
                  corner2 = base.tne + Point(z = thickness))

                # Do the lip part of the cover:
                self.block(comment = "Cover Lip",
                  corner1 = base.tsw + P(thickness, thickness, -lip_thickness),
                  corner2 = base.tne + Point(-thickness, -thickness, zero))
                 
        ezcad = EZCAD3(0)                # Using EZCAD 3.0
        my_assembly = My_Assembly(None)  # Initialize top-level sub-assembly
        my_assembly.process(ezcad)       # Process the design

### Final Stuff ###

The final step is:

        ezcad = EZCAD3(0)                # Using EZCAD 3.0
        my_assembly = My_Assembly(None)  # Initialize top-level sub-assembly
        my_assembly.process(ezcad)       # Process the design

In this code sequence, an *EZCAD3* object is allocated and initialized.
The major revision number is 3 and the minor revision number is 0.
Next, the top-level sub-assembly is allocated and initialized.
Finally, the *process*() method is invoked to cause the design
to be processed.  For those of you who like to scrunch everything
onto one line:

        My_Assembly(None).process(EZCAD3(0))

will also do the trick.

____________________________________________________________________

#Older stuff#

{This stuff is left over from EZCAD 1 and EZCAD 2.  It will be
largely deleted.}

The designer will spend most of their time manipulating Part's,
Point's, and Length's.

Terminology differs between various programming languages.  For
this document, the term "routine" corresponds to Python methods,
C functions, Pascal procedures, EasyC routines, etc.  While the
examples in this document are written in Python, the more generic
term of "routine" is used instead of the Python specific term
of "method".

When the task of designing an mechanical assembly is
approached, it seems intuitive to use a bottom-up approach
where each physical part is designed, these physical parts
are assembled into sub-assemblies and the sub-assemblies
are assembled into a final assembly.  In practice, this
intuitive bottom-up approach does not work very well.
For example, suppose that the designer wants some holes
for some screws to attach to pieces together.  If the
hole in part is moved, it is necessary to move the
hole in the other part.  The interconnectedness of the
design process makes a pure bottom-up design process
less efficient the long run.

EZCAD deals breaks the design process into the three
distinct phases listed below:

* Phase 1 -- Part Tree Construction:
  All of the physical parts and sub-assemblies are
  arranged into a tree with a root assembly at the
  base of the tree.  The part tree provides the basic
  part naming structure needed for the next phase.

* Phase 2 -- Dimension Discovery:
  In dimension discovery, the dimensions for all of
  physical parts and sub-assemblies are computed.
  This is an iterative process that eventually
  results in stable dimensions that can be used
  for the final phase.

* Phase 3 -- Construction:
  In the final construction phase, the construction
  details for each part are described.  The results
  of this phase are a physical description of the
  geometry of each file, a unified visualization of
  the final assembly, and any appropriate files needed
  to assemble the physical parts.

These three phases are directly visible to the designer
as they write the EZCAD code for each part.

In EZCAD, each part is defined by a single part routine
that has the following rigid code structure:

        # Example boilerplate is in Python syntax:

        def my_part(parent):

            # Phase 1: Part tree construction code goes here:
            part = ...
	
            # Enter dimension discovery:
            if part.dimension_mode():

            # Phase 2: Dimension discovery code occurs here:

            if part.construct_mode():

                # Phase 3: Part construction code goes here:

        # Return the final {Part} object:
        return part.done()

The EZCAD system invokes each part routine at least once
for each phase.  There are two routines called dimension_mode()
and construct_mode() that return different values depending
upon the mode:

<BlockQuote>
  <Table Border="1">
    <TR>
      <TH></TH>
      <TH>dimension\_mode()</TH>
      <TH>construct\_mode()</TH>
    </TR><TR>
      <TD>Phase 1 (Part Tree)</TD>
      <TD>False</TD>
      <TD>False</TD>
    </TR><TR>
      <TD>Phase 2 (Dimensions)</TD>
      <TD>True</TD>
      <TD>False</TD>
    </TR><TR>
      <TD>Phase 3 (Construct)</TD>
      <TD>True</TD>
      <TD>True</TD>
    </TR>
  </Table>
</BlockQuote>

All of the variables defined for phase 1 are available to
phases 2 and 3.  All of the variables defined for phase 2
are also available for phase 3.

The actual call to get EZCAD to process a design is
one line:

    EZCAD.process(root_assembly_routine)

That one call performs all of the work.

That pretty much covers the EZCAD overview.  What comes
next is a basic discussion of the primitive types (i.e.
scalars, Angle, Color, Length, Point, etc.)  After that,
the routines that make up a part routine are introduced.

## History

 Its predecessors
are the experimental EZCAD (version 1) and FunCAD, neither
of which are used any longer.  EZCAD has some similarities
to the OpenSCAD project and RapCAD project.


## Basic Types

This section introduces the basic types are scalars (i.e.
**float** and **int**), Length, Angle, Point, Color,

### Scalars

A scalar is either an Integer or Double.  The EZCAD library
tends to treat Integer and Doubles interchangeably.

An Integer is a signed 32-bit quantity.  For example, 0, 1,
1234, -1, and -123 are all valid Integers.  All of the usual
operations on Integer's are supported by the underlying
programming language (e.g. Python.)

A Double is a 64-bit floating point number.  In general, a
Double constant has a decimal point (e.g. 3.14159265358979323846).
Double constants can be represented in E-format (exponent format)
as well (e.g. -3e-3, 3e3 .31456e1 123.)  Again, the underlying
programming language provides a rich set of operations for
manipulating the Double's.

### Color

Colors are represented as a string name.  The list of color names
are the same as the color names defined for the SVG (Scalable Vector
Graphics) standard.  These color names are:

    alice_blue, antique_white, aqua, aquamarine, azure, beige, bisque,
    black, blanched_almond, blue, blue_violet, brown, burlywood,
    cadet_blue, chartreuse, chocolate, coral, corn_flower_blue,
    corn_silk, crimson, cyan, dark_blue, dark_cyan, dark_goldenrod,
    dark_gray, dark_green, dark_grey, dark_khaki, dark_magenta,
    dark_olive_green, dark_orange, dark_orchid, dark_red, dark_salmon,
    dark_sea_green, dark_slate_blue, dark_slate_gray, dark_slate_grey,
    dark_turquoise, dark_violet, deep_pink, deep_sky_blue, dim_gray,
    dim_grey, dodger_blue, fire_brick, floral_white, forest_green,
    fuchsia, gainsboro, ghost_white, gold, goldenrod, gray, green,
    green_yellow, grey, honey_dew, hot_pink, indian_red, indigo,
    ivory, khaki, lavender, lavender_blush, lawn_green, lemon_chiffon,
    light_blue, light_coral, light_cyan, light_goldenrod_yellow,
    light_gray, light_green, light_grey, light_pink, light_salmon,
    light_sea_green, light_sky_blue, light_slate_gray, light_slate_grey,
    light_steel_blue, light_yellow, lime, lime_green, linen, magenta,
    maroon, medium_aquamarine, medium_blue, medium_orchid, medium_purple,
    medium_sea_green, medium_slate_blue, medium_spring_green,
    medium_turquoise, medium_violet_red, mid_night_blue, mint_cream,
    misty_rose, moccasin, navajo_white, navy, old_lace, olive,
    olive_drab, orange, orange_red, orchid, pale_goldenrod, pale_green,
    pale_turquoise, pale_violet_red, papaya_whip, peach_puff, peru,
    pink, plum, powder_blue, purple, red, rosy_brown, royal_blue,
    saddle_brown, salmon, sandy_brown, sea_green, sea_shell, sienna,
    silver, sky_blue, slate_blue, slate_gray, slate_grey, snow,
    spring_green, steel_blue, tan, teal, thistle, tomato, turquoise,
    violet, wheat, white, white_smoke, yellow, yellow_green

### Material

Materials are represented as colon separated strings of either
the form "general_material" or "general_material:specific_material".
The general materials are:

    aluminum, brass, bronze, copper, iron, plastic, steel, wood

Specific materials are particular to the general material.  For
example aluminum could have "6061", "6063", etc.  Plastic could
be "acrylic", "abs", "pvc", "hdpe", "nylon", "derlin", "teflon",
etc.

### Angle

Angles represent an angle between two intersecting lines.  There
are 2 routines for creating angles:

    deg(scalar_degrees) => Angle
    rad(scalar_radians) => Angle

In Python these routines are static methods and defined as:

    Angle.deg(scalar)
    Angle.rad(scalar)

There are two routines for converting from an {Angle} back to
a Double:

    degrees(Angle) => Double
    radians(Angle) => Double

In Python, these routines are methods in the class {Angle}.

The following operations are available for Angles:

    addition, subtraction, and negation:
	a = Length.deg(90)
	b = Length.rad(3.1415926)
	c = a + b
	d = a - b
	e = -a

    comparison:
	a_lt_b = (a < b)
	a_le_b = (a <= b)
	a_gt_b = (a > b)
	a_ge_b = (a >= b)
	a_eq_b = (a == b)
	a_ne_b = (a != b)

    scalar multiplication and division:
	twice_a = 2.0 * a
	twice_b = b * 2.0
	half_a = a / 2.0

The following trigonometric routines are available:

    sine(angle) => Double
    cosine(angle) => Double
    tangent(angle) => Double

**Length**

A Length is a measurement of distance.  In EZCAD, most values used
are Length's (not scalars).  A scalar is converted to a Length
using some simple routines:

	m(scalar) => Length
	cm(scalar) => Length
	mm(scalar) => Length
	ft(scalar) => Length
	inch(scalar) => Length	# "in" is a reserved word in Python

A Length can be converted back into a scalar using the reverse
routines:

	millimeters(Length) => Double
	centimeters(Length) => Double
	meters(Length) => Double
	inches(Length) => Double
	feet(Length) => Double

In general, the scalar to Length conversion routines use a
shortened abbreviation for the routine name (i.e. m, cm, mm,
ft, etc.)  In Python, "inch" is used instead of "in" because
"in" is a reserved word.  The Length back to scalar routine
uses a fully spelled out name (i.e. millimeters, centimeters,
inches, etc.)

For Python, the scalar to Length conversion routines are
defined as static methods.

	m	Length.m
	cm	Length.cm
	mm	Length.cm
	ft	Length.ft
	inch	Length.inch

Typing out Length.m each time becomes a pain in the rear.  As
a convenience, the designer can put the following statements
at the front of the program to reducing wear and tear on
fingers and keyboards:

	m = Length.m
	cm = Length.cm
	mm = Length.mm
	ft = Length.ft
	inch = Length.inch

The Length to scalar conversion routines are implemented as
regular Python methods.

The following operations are permitted on Length objects:

    addition, subtraction, and negation:
	a = cm(1.0)
	b = cm(2.0)
	c = a + b
	d = a - b
	e = -a

    comparison:
	a_lt_b = (a < b)
	a_le_b = (a <= b)
	a_gt_b = (a > b)
	a_ge_b = (a >= b)
	a_eq_b = (a == b)
	a_ne_b = (a != b)

    scalar multiplication and division:
	twice_a = 2.0 * a
	twice_b = b * 2.0
	half_a = a / 2.0

There are a number of routines for manipulating Length objects:

    sine(length, angle) => Length	# length * sine(angle)
    cosine(length, angle) => Length	# length * cosine(angle)
    arc_tangent2(y_length, x_length) => Angle
    twice(length) => Length		# 2.0 * length
    half(length) => Length		# length / 2.0

**Point**

A Point represents a location relative to the origin of a
specific Part.  While the Part type has not been introduced yet,
suffice to say each Point uses a Part to provide the frame of
reference for interpreting the Point location relative to.
A Point has fields x, y, and z.

A point can be created as in two ways:

	# Python syntax:

	p = Point(part, cm(1), cm(2), cm(3))
	p = part.point_new(cm(1), cm(2), cm(3))

There are a whole slew of operations available Point's:

    addition, subtraction, and negation:
	a = Point(part, cm(1), cm(2), cm(3))
	b = Point(part, cm(2), cm(4), cm(6))
	c = a + b		# Must share common Part
	d = a - b		# Must share common Part
	e = -a

    comparison:

	a_eq_b = (a == b)
	a_ne_b = (a != b)

    scalar multiplication and division:

	ax3 = a * 3.0
	bx4 = 4.0 * b
	a_5 = a / 5.0

An adjustment routine is used to make a copy of a point
by adding to one or more of the x, y, and/or z fields:

	new_point = old_point.x_adjust(cm(1))
	new_point = old_point.y_adjust(cm(2))
	new_point = old_point.z_adjust(cm(3))
	new_point = old_point.xy_adjust(cm(1), cm(2))
	new_point = old_point.xz_adjust(cm(1), cm(3))
	new_point = old_point.yz_adjust(cm(2), cm(3))
	new_point = old_point.xyz_adjust(cm(1), cm(2), cm(3))

Of course, the part, x, y, and z fields can be read out using
standard record notation:

	point_x = point.x
	point_y = point.y
	point_z = point.z
	point_part = point.part

It needs to be mentioned that point addition and subtraction
only work of both Point object have the same Part as the
frame of reference.  It is possible to perform addition
and subtraction of Points from different parts, but it
is more complicated and discussed in a separate section
much further below.


# Part Tree

The first phase of EZCAD is to build the part tree.  The
part tree resembles a hierarchical file system.  Each
part routine does 3 things to construct the Part tree.

     1) The first task is to call the start() routine
	to create the part and give it a name.

     2) The second task is to call any child part routines.

     3) The third task is tor call done()

## Call Tree**

The boilerplate looks as follows:

    def assembly(parent):
	# Create the part:
	a = parent.start("assembly")

	# Invoke the child part routines:
	pr1(a)
	#...
	prN(a)

	# Do the discovery and construct stuff here.

	return a.done()

Basically, the part routines are arranged into a call tree
with a root part routines, that calls other more part
routines until each part routine has been called.  In general,
physical parts are at the leaves of the call tree and assembly
parts are the interior nodes of the call tree.

By using a call tree architecture, a team of designers can
each be responsible for providing a part routine for their
respective part assemblies.  Final integration occurs when
each sub assembly part routine is assembled into a single
root part routine.

For example, let's assume that we have an assembly that
consists one physical part and two other sub-assemblies.
The first sub-assembly has two additional physical parts
and the second sub-assembly has three other physical parts.
These are shown in outline format below:

	assembly
	    part1
	    sub_assembly1
		part1a
		part1b
	    sub_assembly2
		part2a
		part2b
		part2c

The routines for creating these parts are:

    # Python syntax:

    def assembly(parent):
	a = parent.start("assembly")
	p1 = part1(a)
	sa1 = sub_assembly1(a)
	sa2 = sub_assembly2(a)
	#...
	return a.done()

    def part1(parent):
	p1 = parent.start("part1")
	#...
	return p1.done()

    def sub_assembly1(parent):
	sa1 = parent.start("sub_assembly1")
	p1a = part1a(sa1)
	p1b = part1b(sa1)
	#...
	return sa1.done()

    def part1a(parent):
	p1a = parent.start("part1a")
	#...
	return p1a.done()

    def part1b(parent):
	p1b = parent.start("part1b")
	#...
	return p1b.done()

    def sub_assembly2(parent):
	sa2 = parent.start("sub_assembly2")
	p2a = part2a(sa2)
	p2b = part2b(sa2)
	p2c = part2c(sa2)
	#...
	return sa2.done()

    def part2a(parent):
	p2a = parent.start("part2a")
	#...
	return p2a.done()

    def part2b(parent):
	p2b = parent.start("part2b")
	#...
	return p2b.done()

    def part2c(parent):
	p2c = parent.start("part2c")
	#...
	return p2c.done()

After EZCAD has invoked the root assembly part routine the
first time, it has created all of the Part objects and
assembled them into a tree of Part's.

## Part Paths

A part path is a name that allows one part to consistently
reference another part.  A part path looks just like a file
name where the slash character ('/') separates part names.
The special name ".." means go up one level.  Using the part
tree from the example above:

	assembly
	    part1
	    sub_assembly1
		part1a
		part1b
	    sub_assembly2
		part2a
		part2b
		part2c

The part named "part1" can be referenced from part "part2c" as "../part1".
The part named "part1a" can be referenced from part "part2c" as
"../sub_assembly1/part1a".  The part named "part2b" can be referenced
from part "assembly" as "sub_assembly/part2b".  Some code examples
should help to clarify:

    def assembly(parent):
	# Phase 1: Part tree creation:
	a = parent.start("assembly")
	p1 = part1(a)
	sa1 = sub_assembly1(a)
	sa2 = sub_assembly2(a)

	if a.dimensions_mode():
	    # Dimension mode:
	    p1a = a.part("sub_assembly1/part1a")
	    p1b = a.part("sub_assembly1/part1b")
	    p2a = a.part("sub_assembly2/part2a")
	    p2b = a.part("sub_assembly2/part2b")
	    p2c = a.part("sub_assembly2/part2c")

	    # ...

	return a.done()

    def part2c(parent):
	# Phase 1: Part tree creation:
	p2c = parent.start("part2c")

	if p2c.dimensions_mode():
	    # Dimension mode:
	    a = p2c.part("../..")
	    p1 = p2c.part("../../assembly/part1")
	    sa1 = p2c.part("../../assembly/sub_assembly1")
	    p1a = p2c.part("../../assembly/sub_assembly1/part1a")
	    p1b = p2c.part("../../assembly/sub_assembly1/part1b")
	    sa2 = p2c.part("../../assembly/sub_assembly2")
	    p2a = p2c.part("../part2a")
	    p2b = p2c.part("../part2b")

	    # ...

	return p2c.done()

The paths in the part2c() routine above are kind of long since
they all start from the part in variable "p2c".  Part paths
can be computed relative to any part that has been looked up.
The part2c() routine is rewritten below to show how this can
be done:

    def part2c(parent):
	# Phase 1: Part tree creation:
	p2c = parent.start("part2c")

	if p2c.dimensions_mode():
	    # Dimension mode:
	    sa2 = p2c.part("..")
	    p2a = sa2.part("part2a")
	    p2b = sa2.part("part2b")
	    a = sa2.part("..")
	    p1 = a.part("part1")
	    sa1 = a.part("sub_assembly1")
	    p1a = sa1.part("part1a")
	    p1b = sa1.part("part1b")

	    # ...

	return p2c.done()

Ultimately, it is up to the designer to decide which parts
they need to access and which paths to use to reach them
from within a part routine.


## Dimension Discovery

Dimension discovery is a big deal.  It is what allows the
designer to construct a parametric design, where one change
can properly propagate through the design and still wind
up with a correct design.

### Value Storage and Retrieval

The first part of dimension discovery is the ability to store
named values with each part.  The value can be either a scalar,
Angle, Length, or Point.

Each part can be used as a convenient place to store scalar, Length,
Angle, and Point values.  The storage routines are:

	angle_set(part, angle_name, angle_value) => Angle
	length_set(part, length_name, length_value) => Length
	point_set(part, point_name, point_value) => Point
	scalar_set(part, scalar_name, scalar_value) => scalar
	string_set(part, string_name, string_value) => string

Each of these set routines will return the value passed in
as the last argument.  That way these arguments can be further
used in expressions.  These routines should only be called
in the dimension discovery section of a part routine.  A given
angle/length/point/scalar/string name should only be set once;
an error occurs otherwise.

The values can be accessed using the appropriate "get" access
routines listed below:

	angle(part, angle_path) => Angle
	length(part, length_path) => Length
	point(part, point_path) => Point
	scalar(part, scalar_path) => scalar
	string(part, string_path) => string

All four of these routines can take a path.  The path has
the form of:

	{part_path} . {angle/length/point/scalar/string name}

For example:

	"../../sub_assembly1/part1.hole1_diameter"

corresponds to the part path "../../sub_assembly1/part1"
and the length name "hole_diameter".

Dimension discovery mode, EZCAD occurs is two sub-phases.
In sub-phase A, the value names are defined.  In sub-phase B,
the value names are stored.  Sub-phase B is done repeatably
until the stored values stop changing.

An example of all of this is shown below:

        def assembly(parent):
	# Create the assembly {Part}:
	p = parent.start("my_part")

	if p.dimensions_mode():
	    # Compute part dimensions:
	    op = parent.part("other_part")
	    width = op.length("width")
	    corner1 = p.new_point(-width / 2.0, -cm(2), -cm(3))
	    corner2 = p.new_point( width / 2.0,  cm(2),  cm(3))
	    p.point_set("corner1", corner1)
	    p.point_set("corner2", corner2)

	    # Create the part:
	    p.block_create(corner1, corner2)

	    if p.construct_mode():
		# Construct the part here:

In this example, one scalar value called "width" is created
before and the Point values called "corner1" and "corner2"
are created.

## Bounding Box Values

All Parts have some point values that define a bounding box that
enclose all points of the Part.  There are four categories --
1) face values, 2) edge values, 3) corner_values, and
4) miscellaneous values.  They are listed by category and
name below:

    Face Values (6 values):
	$T		# Center of bounding box top (Z axis)
	$B		# Center of bounding box bottom (Z axis)
	$N		# Center of bounding box north face (Y axis)
	$S		# Center of bounding box south face (Y axis)
	$E		# Center of bounding box east face (X axis)
	$W		# Center of bounding box west face (X axis)

    Edge Values (12 values):
	$BN		# Center of bounding box bottom north edge
	$BS		# Center of bounding box bottom south edge
	$BE		# Center of bounding box bottom east edge
	$BW		# Center of bounding box bottom west edge
	$NE		# Center of bounding box north east edge
	$NW		# Center of bounding box north west edge
	$SE		# Center of bounding box south east edge
	$SW		# Center of bounding box south west edge
	$TN		# Center of bounding box top north edge
	$TS		# Center of bounding box top south edge
	$TE		# Center of bounding box top east edge
	$TW		# Center of bounding box top west edge

    Corner values (8 values):
	$BNE		# Bottom north east corner of bounding box
	$BNW		# Bottom north west corner of bounding box
	$BSE		# Bottom south east corner of bounding box
	$BSW		# Bottom south west corner of bounding box
	$TNE		# Top north east corner of bounding box
	$TNW		# Top north west corner of bounding box
	$TSE		# Top south east corner of bounding box
	$TSW		# Top south west corner of bounding box

    Miscellaneous values (3 values):
	$D		# Diagonal from $BSW to $TNE
	$O		# Origin (0,0,0) (Note: may be outside of bounding box)
	$C		# Center of bounding box

It is not legal to set any value that starts with '$'.  Anytime the
bounding box changes, these values are automatically updated.

## Part Assemblies

A Part object is used to either construct a physical part or
an assembly.  This section discusses part assemblies only.
Physical part construction is discussed in a section further
below.

A part assembly consists of one or more Part objects that
are placed in relation to assembly origin.  There are two
placement routines -- place() and place_rotated():

    place(assembly_part, sub_part, place_name, translate_point)
	Place {sub_part} at {translate_point} relative to
	{assembly_part}.  {place_name} is used in place paths
	(place paths are described later) to select amongst
	placement locations.

    place_rotated(assembly_part, sub_part, place_name, 
      center_point, axis_point, angle, translate_point)
	Move {sub_part} from {center_point} to the origin,
	rotate by {angle} about the line from the origin
	through {axis_point}, return the part back to
	{center_point} and finaly translate to {translate_point}
	{place_name} is used in place paths (described later)
	to select amongst placement locations.

Some examples of place() and place_rotated() in action.

    # Python syntax:

    def assembly(parent):
	# Phase 1: Part tree creation:
	a = parent.start("assembly")
	sa1 = sub_assembly1(parent)
	sa2 = sub_assembly2(parent)
	p1 = part1(parent)
	p2 = part2(parent)

	if a.dimensions_mode():
	    # Dimension mode:

	    # Create some locations and rotation axes:
	    sa1_location1 = a.point_new(cm(1), cm(0), cm(-1))
	    sa1_location2 = a.point_new(cm(1), cm(0), cm(1))
	    sa2_location1 = a.point_new(cm(-1), cm(0), cm(-1))
	    sa2_location2 = a.point_new(cm(-1), cm(0), cm(1))
	    p1_location = a.point_new(cm(0), cm(0), cm(-1))
	    p2_location = a.point_new(cm(0), cm(0), cm(-1))
	    p1_axis = a.point_new(cm(1), cm(0), cm(0))
	    p2_axis = a.point_new(cm(-1), cm(0), cm(0))

	    # Perform sub-assembly placement routines:
	    a.place(sa1, "bottom", sa1_location1)
	    a.place(sa1, "top", sa1_location2)
	    a.place(sa2, "bottom", sa2_location1)
	    a.place(sa2, "top", sa2_location2)

	    # Place both {p1} and {p2} {Part}'s:
	    a.place_rotated(p1, p1_location, p1_axis, deg(90))
	    a.place_rotated(p2, p2_location, p2_axis, -deg(90))

	return a.done()

An error occurs if part name and place name in a placement is
not unique for a part assembly.  An assembly part has no need
to call the construct_mode() routine.

That pretty much covers part assemblies.

## Complex Point Arithmetic

Now that place names have been described as part of part
placement for assemblies, it is finally possible to discuss
how to subtract two points that have different Part's as
their reference frame.  This section is what makes it easy
for the designer to align parts, holes, etc. without going
bonkers.  This section is where all of the hassle associated
with dimension discovery pays off.  (Hint: this section is
important!)

The routine listed below is the work horse of complex
point arithmetic:

	point_subtract(part, place_path1, place_path2) => Point
	    Return the {Point} that results from subtracting
	    {place_path2} from {place_path1}, where both
	    {place_path1} and {place_path2} start from a {part}.
	    The returned {Point} will use the {Part} associated
	    with {place_path2} as its reference frame.

A place path has the form:

    place_name1/.../place_nameN.point_name

Starting with part0 (i.e. the first routine argument), place_name1
selects part1, followed by place_name2 selecting part2, and so
forth until place_nameN selected partN.  part0 must be an assembly
part.  Indeed, all parts from part0 through partN-1 must be assembly
parts. The last part can be either an assembly part or a chunk part.
The place_name1 must be a valid place name for part0, and place_name2
must be a valid place name for part1, etc.  Since both place_path1
and place_path2 use place names to specify the final point locations,
EZCAD knows exactly were both parts a located in relation to one
another.  Thus, the subtraction can take place and result in a
consistent value.

A common way of using point_subtract is to find out where a point
on another part is relative to the origin of a given part.  In this
example, a1 is an assembly with two parts called part1 and part2.
The hole in part1 is supposed to align with the hole in part2:

    def part2(parent):
	p2 = parent.start("part2")
	if p2.dimensions_mode():
	    # {a1} is the assembly part that contains part1 and part2:
	    a1 = p2.part("..")

	    # The place path "part1/hole1" specifies the location on part1
	    # where a hole is meant to be drilled on part2.  "part2/$O"
	    # is place path that specifies the origin of part2:
	    hole_location = a1.point_subtract("part1/hole1", "part2/$O")
	    p2.point_set("hole2", hole_location)

	    if p2.construct_mode():
		#...
		mm = Length.mm
		p2.hole_through(mm(3), hole_location, "")
		#...
	return p2.done()	

The point_subtract() routine is what allows the designer to
easily and conveniently propagate locations between parts to
simplify the task of designing the entire assembly.

## Part Construction

Part construction is the process of fabricating the part is
typically done by starting with a chunk of material and removing
material until the desired shape is created.  The alternative
is exercised by 3D printers and is an additive process where
material is assembled into a final shape.  While EZCAD is currently
more focused on the material removal strategy, a part that is
designed around material removal can frequently be done via
material addition as well.  Eventually, some additional operations
that specifically support material addition will be added to
EZCAD for 3D printer support.

The material removal process emulates the steps that can be
done by classical machine tools -- drilling a hole, milling
a pocket, etc.  In this process has three basic steps:

  1. The initial material is provided,

  2. The material is attached to the machine (called a machine
     setup) so that material can be removed.

  3. The various operations are performed to remove the material.

Steps 2 and 3 are repeated until the part arrives at its final
shape.

Sections below go through all of these steps.

## Initial Material

The initial chunk of material is typically a either a
rectilinear prism (i.e. a block) or a material extrusion
like a C-channel, L-channel, T-slot extrusion.  The dimensions
of this initial chunk are defined up in the dimensions discovery
section of a part routine.  The reason why it is done here is
to ensure that the bounding box for the part gets properly
defined during dimension discovery phase.

For defining the rectilinear prism, the two most common routines
are:

    block_corners(part, corner1_point, corner2_point, color, material)
	Create a block of {material} that will be visualized as
	{color} with one corner at {corner1_point} and the other
	at {corner2_point}.  {color} is a string that matches
	an SVG color (e.g. "red", "sky_blue", "blanch_almond", etc.)
	Material is a string that names the material out of which
	the part is manufactured (e.g. "plastic", "aluminum",
	"steel", etc.)

    block_diagonal(part, 
	Create a block of {material} that will be visualized as
	{color} that is centered on the origin with {diagonal}
	equally spanning the two corners equally.  {color} is a
	string that matches an SVG color (e.g. "blue", "blue_green",
	"dodger_blue", etc.)  Material is a string that names the
	material out of which the part is manufactured (e.g.
	"plastic", "aluminum", "steel", etc.)

In addition to blocks, rods and tubes are permitted:

    rod(part, color, material, start_point, end_point, diameter, sides):
	Create a rod for {self} out of {material} that goes from
	{start_point} to {end_point} is {diameter} round.  The rod
	is visualized with {color}.  {sides} specifies how many
	sides to use for visualization; setting {sides} to zero
	gives a reasonable visualization.  The rod must be aligned
	with one of the X, Y or Z axes. """

    tube(part, color, material, 
      start_point, end_point, diameter, wall_thickness, sides):
	""" Create a tube for {part} out of {material} goes from
	    {start_point} to {end_point} is {diameter} round and
	    has a wall thickness of {wall_thickness}.  The tube is
	    visualized with {color}.  {sides} specifies how many sides
	    to use for visualization; setting {sides} to zero gives a
	    reasonable visualization of a tube.  The tube must be
	    aligned with one of the X, Y or Z axes. """

{There are going to be a number of routines for defining
various extrusions.  They have not been written yet, but
when they are, they will be listed here.}

It is extremely common to machine a part down from some
outsize material into its final dimensions.  This is
sufficiently common, that EZCAD has a special routine
for keeping track of it the extra material.  The following
two routines define how much extra material is added
to the 6 basic dimensions of part:

    extra_xyz(part, comment, x, y, z)
	Define extra set of dimensions for extra material
	in the x, y, and z axis directions.  The material
	is split evenly in both sides.

    extra_ewnstb(part, comment, east, west, north, south, top, bottom)
	Define extra set of dimensions for extra material
	for each specific side of {part}.  All values must be
	zero or positive.

These routine will define 26 additional "shadow" points for
the extra material bounding box.  The point names are just
like the bounding box names, but with an extra "X" after
the "$".  The 6 extra surface names are "$XB", "$XE", "$XN",
"$XS", "$XT", and "$XW".  The 12 edge names are "$XBE",
"$XBN", "$XBS", "$XBW", "$XNE", "$XNW", "$XSE", "$XSW",
"$XTE", "$XTN", "$XTS", and "$XTW".  The 8 corner names
are "$XBNE", "$XBNW", "$XBSE", "$XBSW", "$XTNE", "$XTNW",
"$XTSE", and "$XTSW".  One of these two routines should be
called after the routine that defines the initial chunk
material dimensions.

After the initial chunk of material has been created in the
dimension portion the part routine, all remaining operations
occur in the construction section of the part routine.

## Setup Operations

Before a part can operated on by a mill or a lathe it
must be installed on the machine.  For now, only mill
setups are discussed.  For a mill, the part is inserted
into a vice and CNC software is configured to know
exactly where the part is located.

The manufacture of a part may require multiple setups
in order to manufacture the part.  For example, for
a three axis CNC mill, the holes can only be made in
the vertical axis.  Thus, any holes are make by mounting
the part perpendicular to the vertical axis.  If there
are multiple holes is multiple different directions,
the part must be setup before for each different hole
axis direction.  All holes with the same axis direction
are typically done using the same setup.

EZCAD requires that the part be properly mounted prior
to doing a removal operation.  While the visualization
does not need setup operations, all generated G-code does
need the setup operation.

The routine below is used to mount a part in vice:

    vice_position(part, comment, top_point, vice_point, edge_point)
	This routine will cause {part} to be mounted in a vice with
	{top_point} on top, {vice_point} rotated until until it
	is positioned against the vice jaw.  {edge_point} specifies
	the west edge of the part, so it can be used in a "dowel
	pin" operation.

The vice is assumed to be mounted on the mill with its fixed
jaw aligned with the X axis.  Next, the machine origin is set
so that left corner of the fixed jaw is at X=0 and Y=0.  All
G-code that is generated assumes that it knows where the vice
jaw is.  When a part is loaded into a vice, the operator is
responsible setting Z=0 to the top surface of the vice.
Since the part is touching the fixed vice jaw, the machine
"knows" where Y and Z are located on the part.  The final
operation is to find an X location on the part.  There are
many ways of doing this, using an edge finder, a probe tip,
etc.  EZCAD typically uses a technique called a dowel pin.
With a dowel pin, the part is roughly centered in the vice.
The machine brings a round cylinder with a known dimension,
called a dowel pin, to the point that is where the left edge
of the part should be place.  The machine pauses, and the
machine operator pushes the part to the left until it touches
the dowel pin.  The part is clamped into place and the machine
operation resumes.  Once the part is clamped into place, the
CNC software "knows" exactly where the part is and can perform
the correct operations.

Sometimes mounting a part into the vice will not work because
the machine operations will interfere with the vice, thereby
damaging, the vice, the part, the tool, or the machine.
There are many alternative methods of mounting a part to a
machine.  The next method discussed uses something called
a tooling plate.

A tooling plate is typically made out of an aluminum
plate and has a grid of holes drilled in it.  Typically,
the holes are threaded to make it easier to attach parts
to the plate.  Most tooling plates are mounted rigidly on
the mill table.

In contrast to mounting a tooling plate rigidly on the
mill table, a vice mounted tooling plate is one that is
intended to be mounted into a vice.  The advantages of
a vice mounted tooling plate are:

  - The vice can remain mounted on the milling table.
    It is not necessary to remove the vice prior to
    mounting the tooling plate.

  - If the vice is aligned on the milling table, the
    tooling plate becomes aligned as well.  Thus,
    there is no additional alignment step.

  - Multiple tooling plates allow one part to be
    attached to the tooling plate, while the machine
    is operation on another tooling plate.  Thus,
    a vice mounted tooling plate provides a inexpensive
    alternative to machine tool pallet systems.

In order to use a tooling plate, the part is typically
has a number of holes drilled into it.  These holes
are aligned with the tooling plate grid.  The part
is then attached to the tooling plate for further machine
operations.

The steps for using a vice tooling plate are:

    1)	Mount initial chunk of part material into vice.

    2)	Establish the Z coordinate by touching off.

    3)	Establish the X coordinate using a dowel pin.

    4)	Drill the mounting holes in the part.  These
	holes are designed to align with the tooling plate
	hole grid.

    5)	Remove the chunk from the vice.

    6)	Install the chunk onto the tooling plate.

    7)	Install the tooling plate into the vice.

    8)	Touch off the Z coordinate.

    9)	Touch off the X coordinate using a dowel pin.

    10)	Perform the remaining operations on the material chunk.

In step 1, the chunk of material is frequently a bit over sized
so that it can be machined down to the correct size.  The final
machine down to size occurs in step 10.

    tooling_plate(part, comment, flags)
	Cause a grid of tooling plate holes to be drilled
	into {part}.  {flags} are used to adjust the position
	the tooling plate holes.  The grid is 2 row by 2 column
	(i.e. 2 by 2) but {flags} can be modified to specifically
	set the number of rows and/or columns.  For example,
	"3r 4c" establishes a grid that 3 rows by 4 columns.
	The grid holes are initially distributed to just fill
	the {part} bounding box in the X and Y dimensions.
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
	from the center.
	
{Tooling plates need to be discussed as well.}

## Operations

### Holes, Pockets and Screw Holes

The most common operation is to remove a hole of material
from a part.  This is equivalent to drilling a hole with
a twist drill.  For the larger holes, a mill will typically
machine the hole using circular interpolation to provide
the round opening.

The most generic routine for making a round hole is:

    hole(part, comment, diameter, start_point, end_point, flags):
	Part construct: Make a {diameter} hole in {self} with starting at
	{start_point} and ending at {end_point}.  {comment} will show in
	any error messages and any generated G-code.  The allowed flag
	letters in {flags} are:

	  One of 't' (default), 'f', or 'p':
	    't'	through hole (i.e. Through)
	    'f'	flat hole (i.e. Flat)
	    'p'	tip hole (drill tip stops at {end_point} (i.e. tiP)

	  Allowed additional flags:	  
	    'u' upper hole edge should be chamfered (i.e. Upper)
	    'l' lower hole edge should be chamfered (i.e. Lower)
	    'm' hole is to be milled (i.e. Milled)

Since through holes are so common, there is a specific routine
for doing through holes:

    hole_through(part, comment, diameter, start_point, flags):
	Make a {diameter} hole in {part} with starting at {start_point}
	and going all the way through {part}.  {comment} will show in any
	error messages and any generated G-code.  The allowed flag letters
	in {flags} are:

	  Zero, one or more of the following:
	    't'	through hole (i.e. Through)
	    'u' upper hole edge should be chamfered (i.e. Upper)
	    'l' lower hole edge should be chamfered (i.e. Lower)

Related to a hole is a simple pocket which removes a rectangular
pocket of material from the part.  The following routine does a
simple pocket.

    simple_pocket(part, comment, corner1_point, corner2_point, radius, flags)
	Part construction: Mill a pocket in {part} where {corner1_point}
	and {corner2_point} specify a diagonal across the pocket.  The
	radius of the inside corners is {radius}.  {comment} will
	show up in error messages and any generated G-code.  The allowed
	flag letters in {flags} are:

	  Zero, one or more of the following:
	    't'	through hole (i.e. Through)
	    'u' upper hole edge should be chamfered (i.e. Upper)
	    'l' lower hole edge should be chamfered (i.e. Lower)

One of the more common methods for attaching multiple parts of an
assembly together is to use a screw fastener.  When it comes to
screws there are two broad families that are in use world-wide --
imperial (i.e. inches) and metric (i.e. millimeters.)  While there
are religious debates about which screw family is the "the best",
EZCAD is indifferent, since it supports both screw families.

The imperial screw format is of the form "{number}-{tpi}"
or "{fraction}-{tpi}", where tpi stands for threads per inch,
number is a number between 0 and 16, and fraction is a screw
diameter of the form {numerator}/{denominator} or
{whole_number}-{numerator}/{denominator}.  Some examples are
"0-80", "2-56", "4-40", "6-32", "8-32", "10-24", "10-32", "12-24",
"1/4-20", "1-12", "1-1/4-12", etc.  In general, the smaller
screws use the "number" prefix, and the larger screws use the
fractional prefix.  In addition, to "0", there is also "00", and
"000" as in "00-90" and "000-120".  In general, the most common
number drill sizes are 0, 2, 4, 6, 8, 10, and 12.  The metric
screw size have the from "M{diameter}x{pitch}" where both diameter
and pitch are measured in millimeters.  In the United States,
the imperial screws are quite common and the in the rest of the
world the metric screws are more common.  When designing an
assembly, care should be taken to select screws that are locally
available.

When making holes for screws, there are tapping holes and
insertion holes.  Tapping holes are designed to be threaded
with a tapping tool.  Insertion holes are slightly larger than
the screw hole and are designed to have the screw fit through
them with no binding.  When a hole is tapped, the threads are
usually tapped at either 75% for soft materials (e.g. plastic,
aluminum, etc.) or 50% for hard materials (e.g. steel.)
Since EZCAD knows what material is being used, it selects the
correct hole to properly thread the hole.  For insertion,
holes, there is a concept of a "close fit" and a "loose fit".
For a close fit, the screw hole is just a little bit larger
than the outside screw diameter.  For a loose fit, the screw
is holes is a bit larger than a close fit.  This allows for
more "slop" (technical term) in assembly the parts.

When it comes to specifying a screw, either a metric or
imperial screw size can be specified.  Having said that,
there are some imperial and metric that can be interchanged
when it comes to insertion holes.  Here is a short chart
showing some screws that are comparable to one another:

      Imperial		Metric
      ========		======
	0-80		M1.6x.35
	2-56		M2x.4
	4-40		M3x.6
	6-32		M3.5x.7
	8-32		M4x.8
	10-24		M5x.9
	12-24		M6x1.0

For the screw hole routines below, the designer can specify
the screw size in using the preferred screw family.  However,
if the designer adds the 'i' flag to the screw routine, it
allows the EZCAD to adjust the design for either imperial
or metric screw availability.  The preferred screw family
is specified separately from the part design.

The following two routines are used to produce screw holes
in a part:

    screw_hole(part, comment, screw, start_point, end_point, flags):
	Make a hole for {screw} in {part} with starting at
	{start_point} and ending at {end_point}.  {comment} will
	show in any error messages and any generated G-code.
	The allowed flag letters in {flags} are:

	  One of 't' (default), 'f', or 'p':
	    't'	through hole (i.e. Through)
	    'f'	flat hole (i.e. Flat)
	    'p'	tip hole (drill tip stops at {end_point} (i.e. tiP)

	  One of 'd', 'c' (default), or 's':
	    'd'	hole is to be tapped (i.e. threaDed)
	    'c'	hole is a close fit (i.e. Close fit)
	    's'	hole is a loose fit (i.e. LooSe fit)

	  Allowed additional flags:	  
	    'u' upper hole edge should be chamfered (i.e. Upper)
	    'l' lower hole edge should be chamfered (i.e. Lower)
	    'm' hole is to be milled (i.e. Milled)
	    'i'	imperial vs. metric swapping allowed (i.e. Interchangeable)

    screw_through(part, comment, screw, start_point, flags):
	Make a hole for {screw} in {part} with starting at
	{start_point} all the way through {part}.  {comment} will
	show in any error messages and any generated G-code.
	The allowed flag letters in {flags} are:

	  One of 'd', 'c' (default), or 's':
	    'd'	hole is to be tapped (i.e. threaDed)
	    'c'	hole is a close fit (i.e. Close fit)
	    's'	hole is a loose fit (i.e. LooSe fit)

	  Allowed additional flags:	  
	    'u' upper hole edge should be chamfered (i.e. Upper)
	    'l' lower hole edge should be chamfered (i.e. Lower)
	    'm' hole is to be milled (i.e. Milled)
	    'i'	imperial vs. metric swapping allowed (i.e. Interchangeable)

### Other Operations

To mill an exterior path, first the path is specified by repeated
calls to the following routine:

    corner(part, comment, corner_point, radius)
	Part construct: Add a corner with a radius of {radius} to
	{part} using {corner_point} to specify the corner location.
	{comment} will show up any error messages or generated G-code.

The corners have to be specified in either a clockwise or counter
clockwise direction.  For inside corner, {radius} specifies the
radius inside the corner.  For an outside corner, {radius} specifies
the rounding of the corner.

The final routine for exterior milling is:

    contour(part, comment, start_point, end_point, extra, flags):
	Part construct: Cause the current list of corners associated
	with {self} to be removed starting at a depth of {start_point}
	and ending at a depth of {end_point} using an exterior contour.
	{extra} is the amount of extra material being removed. {flags}
	can contain the letter 'u' for an upper chamfer, 'l' for a lower
	chamfer, and 't' for a contour that cuts entirely through {part}.
	{comment} is will show up in error messages and any generated G-code.

The call to this routine also clears the corner list.

Sometimes it is desirable to "lathe" a tube of material off of part
when it is mounted in a vice.  The routine below accomplishes this:

    vertical_lathe(part, comment, axis_start_point, axis_end_point,
      inner_diameter, outer_diameter, flags):
	Part construct: Mill out a tube of material from {part} along
	the axis from {axis_start_point} to {axis_end_point} where
	{inner_diameter} and {outer_diameter} specify the tube
	inside and outside diameter.  {flags} can be 'i' to specify
	that only the inside diameter matters, allowing for an end-mill
	that extends past outside_diameter.  {comment} is used in
	error messages and any generated G-code. """

{More operations will be added here.}

## Fasteners

A fastener is a screw, bolt, rivet, threaded rod, etc.
The basic requirement for fasteners is that two parts
that are to be fastened together require two aligned
holes that are the right size and in the right position
through which the fastener is placed.

Machine screws are a common fastener.  A machine screw
has a diameter, thread pitch, length, and a head.  Pan
head, flat head, round head, socket cap, hex head, square
head, oval head, fillister head, etc.  In addition, each
different head style can have a multitude of different drive
stiles such as Phillips, slotted, star/Torx, Pozidriv(r),
etc. When it comes to fasteners, EZCAD does not really care
about the drive style.  EZCAD fasteners do not visually show
the drive style.

    fastener(part, name, spec, length, start_point, end_point, flags)

