# EZCAD3: Easy (Mechanical) Computer Aided Design -- Version 3

## Introduction

EZCAD stands for Easy Computed Aided Design.  EZCAD is a
mechanical CAD (Computer Aided Design) system with an integrated
mechanical CAM (Computer Aided Manufacture) system.  The CAD
system specifies the geometry of individual parts and how they
are assembled together.  The CAM system specifies the
manufacturing steps required to make a part.

In EZCAD, the designer writes a program using the
[Python](http://www.python.org/) programming language.
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

## Downloading, Documentation, and Licensing

### Downloading

This package is distributed in source form from a *git* repository.
To get a copy:

        git clone https://github.com/waynegramlich/ezcad3.git

In addition, this package needs a version of
[OpenSCAD](http://www.openscad.org/).  For those who use
Linux distributions based on the Debian package system
(e.g. Ubuntu, Mint, etc.), the following will do the trick:

        sudo apt-get install openscad

OpenSCAD is primarily used as a front-end to the wonderful
[CGAL](http://www.cgal.org/)(Computational Geometry Algorithms
Library).  In particular, OpenSCAD uses the Boolean 3D solids
library embodied in the NEF3 sub-library of CGAL.

Detailed documentation of the Python classes and associated
methods is maintained by [Doxygen](http://www.doxygen.org/).
For convenience, all of the generated Doxygen documenation
is checked into the repository.  If you want to generated
Doxygen documentation locally, you need to download the program.
I downloaded my copy via:

        sudo apt-get install doxygen

### Documentation

In addition to this document, there is
[Doxygen generated documentation](html/pages.html).

### Licensing

In general, I really like Open Source Licenses.  I have a slight
preference of the GPL open source license, so that is what EZCAD3
is released under.  I emulate the
[Free Software Foundation](http://www.fsf.org) that any code
contributed to my code branch requires a copyright assignment.
Professor Eben Moglen has written up a short explanation of
[Why the FSF gets copyright assignments from
contributors](http://www.gnu.org/licenses/why-assign.html).
I strongly feel that having a unified copyright owner provides
maximum protection for the open source software.  While it is
possible to get
[overly legalistic](http://ftp.xemacs.org/old-beta/FSF/assign.changes),
I think the following is more that adequate to assign the copyright:

        Hello:

        I am assigning my modifications to EZCAD to
        Wayne C. Gramlich with the requirement that the code
        will continue to be released under the GPL version 3
        license or higher.

        Regards,

        {Your name here}

Enough on this legal stuff!

## EZCAD Overview

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

Each EZCAD class method validates that its arguments types are correct.
(If you have ever heard of the "duck typing" design pattern, EZCAD
most definitely does **not** use that particular programming religion.)

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

### Part Initialization Phase

Each user individual part is implemented as a Python class object that
sub-classes the *Part* super class.  As expected, the initialization phase
occurs in the <I>__init__</I>() method.  Each part has an up-level parent
(abbreviated as "*up*").  The following code fragment shows how it is
done:

        class My_Part(Part):

            def __init__(self, up):
                Part.__init__(self, up)      # Initialize *Part* super-class
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

### Constraint Propagation Phase

Constraint propagation is what allows the designer to relatively easily
resize a design.  The constraint propagation phase occurs in a
*Part*'s *construct*() method.  Constraint propagation occurs by
repeatably calling the *construct*() method for each *Part* in
the entire design.  The *construct*() methods are repeatably called
until none of the values with approved suffixes (e.g. "\_l", "\_p")
change any more.

To further support constraint propagation, each *Part* maintains a
bounding box for everything it contains.  The bounding box uses
an altitude/compass bearing naming convention.  Using this convention,
the X axis is uses East/Center/West names, the Y axis uses
North/Center/South names, and the Z axis uses Top/Center/Bottom name.
The following crude ASCII art should help to convey the concept.

                    T  N (Y axis)
                    | /
                    |/
              W-----*-----E (X axis)
                   /|
                  / |
                 S  B
                   (Z-axis)

The bounding box of a part defines box that consists of 3 by 3 x 3
points.  There are three slices -- the top slice, the center slice,
and bottom slice that are named as follows:

          tnw--tn---tne     nw----n----ne     bnw--bn---bne
          |     |     |     |     |     |     |     |     |
          tw----t----tw     w-----c-----e     bw----b----bw
          |     |     |     |     |     |     |     |     |
          tsw--ts---tse     sw----s----se     bsw--bs---bse
           (Top Slice)      (Center Slice)    (Bottom Slice)

The designer can easily access the various bounding box values
using short concise member variable names.  Thus, the *tne*
member variable corresponds to the Top-North-West corner of
the bounding box, *c* corresponds to the Center of the bounding
box, and *bw* corresponds to the Bottom-West edge of the
bounding box.

Using constraint propagation and bounding boxes, EZCAD encourages
the use of the  *assembly focused* design pattern where each
part is designed to fit precisely into the final assembly from the
beginning.  This is in contrast to the *part focused* design
pattern where each part is individually designed and then
subsequently fitted into the final assembly.

### Manufacture and Assembly Phase

The manufacturing and assembly phase occurs last.  EZCAD can generate
manufacturing files for the following generic classes of machines:

* 3D Printers (.stl file)

* Laser Cutters (.dxf file)

* CNC Mills and Lathes (CNC file)

EZCAD uses a design for manufacture strategy.  For each part, the
designer specifies one or more appropriate machines that could be
used to manufacture the part.  While some parts can only be manufactured
on a single machine class, others can be manufactured on more than
one machine class.  For example, a plastic plate that has some holes,
inter cut outs, and a specific exterior contour can actually be
manufactured on a CNC mill, laser cutter, or 3D printer.  The
decision of which machine to use is deferred until the actual part
needs to be manufactured.

EZCAD uses a two step process for designing a part:

* First, a bunch of material chunks (e.g. blocks, cylinders,
  extrusions, etc.) are "welded" together to provide a rough
  3D chunk of material.

* Second, a bunch of material removal operations (e.g. holes,
  pockets, exterior contour removal, etc.)

The details of these operations are deferred until later.

After each part is designed, there is a final step of placing
each part into the final assembly.  A single part or sub-assembly
can easily be replicated multiple times in the final assembly.

Once the final assembly is present, it can be viewed using the
appropriate 

### A Quick Example

The example below creates a plastic box with a cover, where the
cover has a lip that keeps the cover centered over the box.
We present the entire code body first then will describe
what is happening on a chunk by chunk basis next:

        #!/usr/bin/env python
        
        from EZCAD3 import *   # The EZCAD (revision 3) classes:
        
        class Simple_Box(Part):
        
            def __init__(self, up, dx=L(mm=100.0), dy=L(mm=50.0),
              dz=L(25.0), wall_thickness=L(mm=5.0),
              material=Material("plastic", "ABS")):
                # Initialize the *Part*:
                Part.__init__(self, up)
        
                # Remember the initialization values:
                self.dx_l = dx
                self.dy_l = dy
                self.dz_l = dz
                self.wall_thickness_l = wall_thickness
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
                # Grab some values from *box*:
                box = self.up
                dx = box.dx_l
                dy = box.dy_l
                dz = box.dz_l
                wall_thickness = box.wall_thickness_l
                material = box.material_m
        
                # Add another 
                self.height_l = height = dz - wall_thickness
                zero = L()
        
                # Start with a solid block of the right dimensions:
                height = dz - wall_thickness
                self.block(comment = "Initial block of material",
                  material = material,
                  color = Color("blue", alpha=.5),
                  corner1 = P(-dx/2, -dy/2, zero),
                  corner2 = P( dx/2,  dy/2, height))
        
                # Pocket out the body of the box:
                self.simple_pocket(comment = "Box Pocket",
                 corner1 = self.bsw + P(wall_thickness,
                   wall_thickness, wall_thickness),
                 corner2 = self.tne - P(wall_thickness,
                   wall_thickness, L(mm=-.1)))
        
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
                wall_thickness = box.wall_thickness_l
                base = box.base_
                base_height = base.height_l
                zero = L()
        
                # Compute local values:
                self.lip_thickness = lip_thickness = wall_thickness/2
        
                # Do the top part of the cover:
                self.block(comment = "Cover Top",
                  material = material,
                  color = Color("green", alpha=0.5),
                  corner1 = base.tsw,
                  corner2 = base.tne + P(z = wall_thickness))
        
                # Do the lip part of the cover:
                self.block(comment = "Cover Lip",
                  corner1 = base.tsw + P(wall_thickness,
                    wall_thickness, -lip_thickness),
                  corner2 = base.tne + P(-wall_thickness,
                    -wall_thickness, zero))
                 
        ezcad = EZCAD3(0)                # Using EZCAD 3.0
        simple_box = Simple_Box(None)   # Initialize top-level sub-assembly
        simple_box.process(ezcad)       # Process the design

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

        Simple_Box(None).process(EZCAD3(0))

will also do the trick.

## CAD Class Reference

This section coveres the design classes.  The lower level
classes are covered first, followed by the *Part* class.
The next major section after this covers the manufacturing
classes (e.g. *Shop*, *Machine*, *Tool*, etc.

### The *Angle* Class

The *Angle* class represents an angle.  While most designers
specify their angles in degrees, internally the *Angle* class
represents angles in radians.

The initializer can specify angles in either degrees or radians:

        degrees0 = Angle()                    # No arguments is 0 degrees
        degrees90 = Angle(deg = 90.0)         # Angle specified in degrees
        degrees180 = Angle(rad = 3.14.15926)  # Angle specified in radians

An *Angle* can be convert back into a *float* using a conversion method:

        d = degrees90.deg()                   # Convert to degrees
        r = degress180.rad()                  # Convert to radians

*Angle*'s can be added, subtracted, multiplied, divided, etc.:

        a = degrees90 + degrees180            # Addition
        b = degrees90 - degrees180            # Subtraction
        c = degrees90 * 3                     # Multiplication
        d = degrees90 / 3                     # Division
        e = -degrees90                        # Negation

Comparison of *Angle* objects is supported:

        eq = degrees90 == degrees90     # Equality
        ge = degrees180 >= degrees90    # Greater than or equal
        gt = degrees180 > degrees90     # Greater than
        le = degress0 <= degress180     # Less than or equal
        lt = degress0 < degress180      # Less than
        ne = degrees0 != degress180     # Inequality

An *Angle* object can be formatted using the standard Python format
facility.  The units can be specified by suffix character:

        angle = Angle(deg = 180.0)
        print("degrees={0:d} radians={0:r}".format(angle))

will print:

        degrees=180.0 radians=3.14159265359

There are three trigimetric functions for angles:

        angle.sine()                    # sin(angle)
        angle.cosine()                  # cos(angle)
        angle.tangent()                 # tan(angle)

There is a miscellaneous *Angle* class:

        angle.normalize()               # Angle between -180 and  180 degrees

Read the
[Doxygen generated *Angle* class information](html/classEZCAD3_1_1Angle.html)
to get the detailed method documenation.


### The *Color* Class

*Color* class documentation goes here.

### The *L* (i.e. Length) Class

The *L* class represents a length.  While the length can be specified
in units of centimeters, millimeters, inches, and feet, internally,
the *L* class converts everything to millimeters.  Here are some
examples of creating a length:

        length0 = L()                   # Length of 0.0
        length2mm = L(mm = 2)           # 2.0 millimeters (note *int* is OK)
        length3_5cm = L(cm = 3.5)       # 3.5 centimeters
        length1thou = L(inch = .001)    # One thousandth of an inch
        length6feet = L(ft = 3)         # 3 feet (which happens to be 1 yard)
        length1_3_4 = L(inch = "1-3/4") # 1.75 inches

An *L* object can be converted back into *float* using a conversion method:

        length.mm()                     # Millimeters
        length.cm()                     # Centimeters
        length.inch()                   # Inches
        length.ft()                     # Feet.

As expected, *L* objects can be added, subtracted, multiplied, divided, etc.:

        a = length2mm + length1though   # Addition
        b = length6feet - length3_5cm   # Subtraction
        c = length6feet * 2             # Mulitplication by *int* (or *float*)
        d = length2mm / 2.3             # Division by a *float* (or *int*)
        e = -length3_5cm                # Negation

Comparisons between *L* object are permitted:

        eq = length2mm == length2mm     # Equality
        ge = length2mm >= length0       # Greater than or equal
        gt = length2mm > length1thou    # Greater than
        le = length0 <= length2mm <=    # Less than or equal
        lt = length1thou < length2mm <  # Less than
        ne = length3_5cm != length6feet # Inequality

*L* object can be formatted using the Python formatting system.
For example:

        length = L(inch = 1)
        print("mm={0:m} cm={0:c} inch={0:i} ft={0:f}".format(length))

will print out:

        mm=25.4 cm=2.54 inch=1.0 ft=0.08333333333

If no suffix is provided, millimetes will be used the for output units.
In general, it is best to always specify the units.

In addition, it is possible to use standard float syntax to futher
control the formatting.  Thus:

        print("ft={0:.4f}".format(length)

will print out:

	ft=0.8333

There are trigametric methods for *L* objects that multiply
the a length by a trigimetric function:

        length.sine(angle)        # == length * angle.sine()
        length.cosine(angle)      # == length * angle.cosine()
        length.tangent(angle)     # == length * angle.tangent()

In addtion there is:

        dy.arctangent2(dx)        # == Angle(rad=math.atan2(dy.mm(), dz.mm()))

There are a few miscellaneous *L* methods as well:

        length1.absolute()        # Absolute value of (length1)
        length1.maximim(length2)  # Maximum of length1 and length2
        length2.minimum(length2)  # Minimum of length1 and length2

The detailed method documentation generated by Doxygen is
[available](html/classEZCAD3_1_1L.html).

### The *Material* Class

*Matierial* class documentation goes here.

### The *P* (i.e. point) class

*P* class documenation goes here.

### The *Part* Class

*Part* Class documenation goes here:

Bounding Box:

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
