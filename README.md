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

## Downloading

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

## Documentation

In addition to this document, there is
[Doxygen generated documentation](html/pages.html).

## Licensing

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

