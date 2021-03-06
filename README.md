# PUMI.jl


## What this Library Is

This code provides a way to use PUMI from Julia by wrapping functions in PUMIs
APF API.
The library is composed of several files.  funcs1.cc, funcs1.h form the C++
half of the
interface, PumiInterface.jl forms the Julia half. All the functions in
PumiInterface.jl,
except for `init` and `init2`, mirror apf functions.  I refer to this as the
low level interface.
PumiInterface2.jl provides some higher level function by using the lower level
 functions.
These functions are more convenient but potentially less efficient depending
on the use case.
The files a2.cc and a2.h contain the user supplied reordering function.  
YOU MAY NOT USE THIS CODE FOR ASSIGNMENT 4.  
The files are provided to show how a user can add functions.  If you have
questions about how to interface your reordering function with Julia,
let me know.

funcs1.h has some useful general comments in it.  
funcs1.cc has more descriptive comments for each function.  

## Installation

Now that the library is a Julia package, you can get is using
Pkg.clone(url) followed by Pkg.build().
The build script will check for all needed components and attempt to build
them itself if it does not find
them.  A C++11 compiler is required.

The next two subsections describe the inner workings of the build system.
See the following to sections for the standard installation procedure.

### Install MPI (if needed)
Building Pumi and building the shared library that links to Pumi requires
`mpicxx`.  If `mpicxx` is not found,  the build script will attempt to
install a Debian package of MPICH3.  It will ask for root permission to
install MPICH.  It will also install the MPI.jl package if needed

### Build Pumi (if needed)
The build script attempts to locate Pumi using `cmake --find-package`.
If Pumi is not located in one of CMake's standard search locations,
users can set the `SCOREC_PREFIX` environment variable to Pumi's installation
prefix.
If Pumi is not found, the build system will attempt to build Pumi locally.
The build requires an
MPI implementation (which should be present according to the above process),
and Cmake 3.0.0.  It will install to `deps/core/build/install`.  If you want
to build other software that links to this Pumi build, source the
script script `deps/use_pumi.sh` first.  This will set the necessary
environment variables.

### Build Shared Library
At this point, MPI and Pumi should be be present, so the library itself can
be built, which is done via CMake.  
The build script executes the shell scripts in the `src` directory `config.sh`
and `makeinstall.sh`, which configure and install the shared library,
respectively.  The library is installed to the `repo/install` directory.
If you wish to rebuild the shared library, you should execute these scripts.  
It is also recommended to execute the `cleanup.sh`
script before doing so, to remove any old configuration files.
Previous version of this package required users to source the `use_julialib.sh`
script before each use.  This is no longer necessary.
It is also unnecessary to set `LD_LIBRARY_PATH` or any other environment
variables at runtime.  This package uses `RPATH` to locate the libraries
at runtime (and an `RPATH`-like mechanism for the Julia code).


Thats it, now you are done.

Note Pumi is built without some of its dependencies, including zoltan, so certain features such as load balancing will not work.

### Installation Summary
Here is a summary of the installation process.  
Pre-requisites:
  * CMake version 3.0.0 or higher
  * (Optionally) MPI (the build system will build it for you if needed, but
    you might want to choose a particular implementation yourself)

The build script will:
  * build MPI if it cannot find an existing installation
  * build Pumi if it cannot find an existing installation
  * build the shared library

## Scorec Installation

If you are working on a SCOREC machine (with access to the SCOREC shared file
system), load the Pumi module before building the package.  This will ensure
the Scorec installation of Pumi is found, which is preferable to building
Pumi locally because it links to Zoltan and other optional dependencies that
provide extra features.


## Non Scorec Installation on Scorec machine
If you are working on a Scorec machine but do not want to use the Scorec build
of Pumi, you will have to build it yourself.  Even if the the Pumi module is
not loaded, it seems CMake is still able to find the Pumi installation,
therefore the build system will use it.  

To get around this:
  * Load the compiler and MPI modules you want to use
  * run the build script normally
  * `cd` into `repo/deps` and run `build_local.sh`
  * `source ./use_pumi.sh`
  *  execute `config.sh` and then `makeinstall.sh`

This will build a local copy of Pumi and tell the shared library to
link to it rather than the system one

## Using the Library

The fundamental idea behind the library is that pointers to objects can be passed to and from Julia,
 but not dereferenced by Julia (bad things would happen if your try, so please do not kick the bucket).
  Therefore, a function that in C++ would be called like this:
```
  m_ptr->getDownward(args...)
```

Gets called from Julia like this:
```
  getDownward(m_ptr, args...)
```
Where m_ptr is a pointer to the PUMI object.  The function arguments that are pointers are
currently left untyped in the function declaration.  It is the users responsibility to pass
 it a pointer to the right kind of object.  The name of the arguments will tell you the type
of object the pointer must point to:

```
m_ptr = apf::Mesh2*
mshape_ptr = apf::FieldShape*
eshape_ptr = apf::EntityShape*
entity = apf::MeshEntity*
mel_ptr = apf::MeshElement*
field = apf::FieldShape*
n_ptr = apf::Numbering*
numbering_ptr = apf::Numbering*
numbering = apf::Numbering*
tag_ptr = apf::MeshTag*
```

If you pass an invalid pointer, PUMI will segfault, typically without giving
any useful error data.

The functions init and init2 take in two strings, the names of the mesh and geometry files,
and return pointers to the mesh (and apf::Mesh2 object), and mesh apf::FieldShape object, as
 well as a few other things.  Using these pointers and other functions provided, you can get
pointers to other object, such as Numberings, Tags, EntityShapes, MeshElements etc.



To get pointers to MeshEntitys, the library provides access to 4 iterators, one over each dimension
 entity.  The functions increment* increment these iterators, and the functions get* return a
 MeshEntity pointer.

Most other functions that require a MeshEntity take it as an argument, but there are still one or two
that use the iterators internally.  Any functions that do this are noted as such in funcs1.h.  It is
 the users responsibility to keep track of the location of the iterator.  The reset* functions reset the
iterator to the beginning.


As a final note, the Julia function names mirror (or are extremely similar to) those in PUMI's APF, so
searching PumiInterface.jl for the name of APF functions should get you to the right place quickly.
If you want to find out what function in funcs1.cc a Julia function calls, look at the first argument
 to the ccall function within the Julia funcs, which will have the forms `functionname_name`, and look up
 the value of this global variable at the top of PumiInterface.jl.  The value will be a string, which is
 the name of the C++ function.

## Types of Meshes
`PdePumiInterface` has two kinds of meshes, Continuous Galerkin (CG) and
Discontinuous Galerkin (DG), in two dimensions, 2D and 3D.  Currently,
2D CG, 2D DG and 3D DG are implemented.

The DG outer constructors constructors are

```julia
PumiMeshDG2{T, Tface}(::Type{T}, sbp::AbstractSBP, opts,
                      sbpface::Tface; dofpernode=1, shape_type=2,
                      comm=MPI.COMM_WORLD)

function PumiMeshDG3{T, Tface}(::Type{T}, sbp::AbstractSBP, opts,
                               sbpface::Tface,
                               topo::ElementTopology{3};
                               dofpernode=1, shape_type=2, comm=MPI.COMM_WORLD)

```

where `T` is the element type for all the data arrays (coordinates, metrics
etc), `sbp` is the SBP operator, `opts` is the options dictionary (see below),
`sbpface` is the SBP face operator associated with `sbp`, and `topo` is the
`ElementTopology` object that describes the topology of the reference SBP
element.

The keyword arguments are

 * `dofPerNode`: number of degrees of freedom on each node
 * `shape_type`: this integer deterines the apf::FieldShape that defines the
               node distribution on an element.  It must match the SBP operator
               provided (see below)
 * `comm`: MPI communicator to use.  Currently only MPI.COMM_WORLD is supported.


The `shape_type` values are found in `getFieldShape` in `func1.cc`.  The
commonly used values are:

0: node set associated with Lagrange polynomials (currently no SBP operator equivalent)
1: SBP Gamma CG (2D only)
2: SBP Omega DG
3: SBP Gamma DG
4: SBP Diagonal E (with vertex nodes)
5: SBP Diagonal E (without vertex nodes)



## Options dictionary
The constructors for the various mesh object take an options dictionary
that set the options about how to construct the object, such as boundary
conditions and coloring distance.  The function `set_defaults` takes a
dictonary and supplies default values for any keys not set.  

## Parallelization
Parallel DG meshes are supported.  The use case for parallelization is to
compute the face integrals between faces on different processes.  To do this,
data for elements on a part boundary is copied to the other side of the
interface.

This also requires coloring the mesh in parallel.  Because data is copied
to the other side of the interface, interfacial elements, there are two
independent copies of these elements, and so they can be considered different
elements for the purpose of coloring.  Thus the global coloring problem can be
decomposed into a local coloring problem on each part.  This requires adjacency
information for the remote elements, which is calculated during initialization.

## Curvilinear grids
Curvilinear grids supported in both 2D and 3D.
For curvilinear grids, dxidx and jac are not calculated at the face nodes (ie.
`dxidx_face`, `dxidx_bndry`, `dxidx_sharedface`, `jac_face`, `jac_bndry`, `jac_shareface`
are arrays of size zero).  Instead, users should use the normal vectors
calculated at the face nodes, `nrm_face`, `nrm_bndry`, `nrm_sharedface`.
These arrays are calculated for non-curvilinear meshes as well, and should
be used whenever possible.  See `getAllCoordinatesAndMetrics()` in `metrics.jl`
for more details on which arrays are populated for which modes.

## Periodic Boundary Conditions
DG meshes support periodic boundary conditions.  The .smb file contains
information on which MeshEntities are periodic with which other
MeshEntities.  These entities will automatically be identified as periodic
and treated accordingly.  It is an error to assign a boundary condition to a
geometric entity that has a periodic MeshEntity on it, and an exception will
be thrown.

When creating structured the way MeshCreate does (creating cubes and then
decomposing each cube into 6 tets), there can be a minimum of 2 cubes in
each direction if that direction is periodic.  In this case, for example,
the y direction is periodic if the xz planes are matched.
If this condition is not satisfied, an error is thrown.

### Implementation Details
First some terminology: Matched entities are two entities at different
locations that are periodic with each other.  Remote entities are two
entities that have the same location, and are used to represent entities that
are shared across part boundaries.

If two matched entities are located on the same part, then an Interface is
created between them in `mesh.interfaces`.  If they are on different parts,
then an Interface is created in `mesh.shared_interfaces`.  In order to
determine the relative orientation of the interfaces, vertices that define
the face are required.  In parallel, the processes exchange remote vertices
to calculate the relative orientation.  When there are both remote and
matched vertices present, it is not possible for the sending process to
to uniquely determine which vertices it needs to send.  Consquently each
process must send all matched and remote vertices shared with the given
receiving process and let the receiving process sort out which vertices are
part of the face. The maximum number of vertices that needs to be sent is
bounded by a constant that depends on the topology of the domain.  For example,
a vertex on a cube can have a maximum of 8 matches and 1 remote on a given
process.  The constants `Max_Vert_Matches` and `Max_Vert_Remotes` tell the
code the maximum number of matched and remote vertices that can be shared with
another process, and are used to determine the size of an array that is
sent via MPI during initialization.

## Visualization
All DG meshes interpolate the solution onto the vertices of the mesh for
visualization.  CG meshes higher than second order create a new
mesh that contains all nodes of the orignal mesh as vertices of a new
triangulation and copy the solution values to it.  CG meshes second
order or less can be visualized directly in Paraview.

SBP Omega elements support an option "exact_visualization" (via the options
dictionary) to subtriangulate the mesh such that all nodes plus all vertices
in the origonal mesh are vertices in the new mesh.  This allows the
visualization of the exact nodal values.  The nodal values are interpolated
to the vertices.

## Mesh Warping

Mesh warping is supported for DG meshes.  The function `update_coords()` allows
setting the new coordinates of each element.  The function `commit_coords()`
must be called after all calls to `update_coords()` are complete.


## Derivative calculation

Constructing a mesh object with `Tmsh = Complex128` is supported.  This allows
using the complex step method to calculate derivative with respect to the mesh
coordinate and metrics.  The function `recalcCoordinatesAndMetrics` can be used
to recalculate the volume coordinates and metrics after the `mesh.vert_coords`
field has been updated.

Reverse mode differentiation is supported, both in serial and in parallel,
with keyword arguments to control how non-owned entities are handled.
 See `getAllCoordinatesAndMetrics_rev()`.


## Geometry Interaction

### Geometric Formats

Loading geometric models is supported through Pumi's `gmi` library, which
is wrapped in the Julia `gmi` module.  Currently, three types of models are
supported:

 * null model: if `dmg_name` in the options dictionary is ".null", then
               the null model will be loaded.  This is a fake model used for
               testing only.  Internally, every time a model entity is queried,
               it is created (so by using the model you build the list of
               geometric entities.
 * MDS model: if `dmg_name` specifies a file with extension `.dmg`, then the
              Pumi `mds` model will be used.  This model contains the topology
              of the original CAD model, but no shape information.  See below
              for how to make these.
 * GeomSim model: If Pumi is linked with the (proprietary) Simmetrix libraries,
                  then it is possible to load certain kinds of `.smd` files.
                  These files contain geometry topology, geometry shape, and
                  may contain additional attribute data that was created during
                  the mesh generation process by SimModeler.  See below for
                  how to create these.

It is also possible to load other kinds of geometry files when Pumi is linked
to Simmetrix, although doing so is discouraged.

 * Parasolid (`.x_t` or `.xmt_txt` extension): these files can only be loaded
             on x86 platforms (no BG/Q or Power9).
 * ACIS: (`.sat` extension): limited Linux support


### Conversion

To demonstrate how to create the different files, lets assume a Parasolid CAD
model called `airfoil.x_t` exists (for example, do steps 1-5 [here](http://www.optimaldesignlab.com/optimaldesignlab-wiki/en/#!pages/pdesolver/airfoil_mesh.md))

To make a GeomSim model from a Parasolid or ACIS model (must be done on an
x86 machine):

If the Parasolid model is an assembly model, it must be converted to a
non-manifold model first.  This can be done by:
```
  - Open SimModeler
  - File - Import Geometry - Select the airfoil.x_t
  - If the left pane shows Assembly 1 as the first item with a list of parts
    below it, then this is an assembly model. If this is not an assembly model,
    skip the remaining steps in this box.
  - Modeling - Make NonManifold Model
  - In the window that appears: File - Save As - enter airfoil_nonmanifold as
    the file name
```

The final step will create the files `airfoil_nonmanifold.smd` and
`airfoil_nonmanifold_nat.x_t`.  The second file contains a new, nonmanifold
Parasolid model.  If the original Parasolid model was not an assembly model,
do `cp -v airfoil.x_t airfoil_nonmanifold_nat.x_t` to make the remainder of this
tutorial simipler.


Now make a GeomSim model from the new Parasolid model:
```
  simTranslate airfoil_nonmanifold_nat.x_t airfoil.smd
```

The `airfoil.smd` file can be opened by SimModeler to make a mesh.  After
creating the mesh, SimModeler will offer the user the chance to save both the
model (`.smd`) and mesh `.sms` files.  It is recommended to use this `smd`,
not `airfoil.smd`, when loading the mesh (more on mesh loading below).

Note that `SimModeler` can open the `airfoil.x_t` directly and create a mesh
off of that, however this is not recommended.  Parasolid and GeomSim use
different geometric representations.  The same geometric model should be
use for creating the mesh and loading the mesh with PumiInterface.



To make an MDS model from the Parasolid model:
```
  mdlConvert airfoil_nonmanifold_nat.x_t airfoil.dmg
```

Note that this conversion loses all geometry shape information, as described
above.  `mdlConvert` can also take a GeomSim `.smd` file as the first argument.

As mentioned above, SimModeler saves an `sms` mesh file.  It is recommended
to convert it to the MDS mesh format `smb`.  Assuming SimModeler files were
saved as `airfoil_meshmodel.smd` and `airfoil_mesh.sms`:

```
  convert airfoil_meshmodel.smd airfoil_mesh.sms airfoil_mesh.smb
```

### SimModeler Quirks

There are 2 types of SimModeler model files, and both have the extension `.smd`,
which makes things a bit confusing.  The first type is an attribute file which
does not contain geometry shape information.  The second type is a GeomSim
file, which contains both geometry shape and attribute information.  Only
the second type can be used by PumiInterface.

If SimModeler opens a Parasolid file, generates a mesh, and then saves the
`.smd` and `.sms` files, this creates an `.smd of the first kind.
If `mdlConvert` is used to convert the Parasolid file to a GeomSim `.smd`,
then this `.smd` is opened by SimModeler, the `.smd` that SimModeler saves
will be a `.smd` of the second kind.  I believe this `.smd` has both
geometry shape and attribute information (the instructions for how SimModeler
made the mesh) in it, but this needs to be verified.

### Recommended Workflow

If you want geometry shape information (required Pumi linked to Simmetrix at
runtime):

 1. Make CAD model and save to Parasolid or ACIS
 2. Make non-manifold model if needed
 3. Use `simTranslate to convert to `foo.smd`
 4. Use SimModeler to open `foo.smd`
 5. Make the mesh and save to `bar.smd`. and `bar.sms`
 6. Use `convert` to convert `bar.sms` to `bar.smb`
 7. Load with PumiInterface by setting `smb_name` to `bar.smb` and `dmg_name` to `bar.smd`


If you don't want geometry shape information (or don't have Pumi linked to Simmetrix at runtime)

 1. Make CAD model and save to Parasolid or ACIS
 2. Make non-manifold model if needed
 3. Use `mdlConvert` to convert to `foo.dmg`
 4. Use SimModeler to open the Parasolid or ACIS file
 5. Make the mesh and save to `bar.smd`. and `bar.sms`
 6. Use `convert` to convert `bar.sms` to `bar.smb`
 7. Load with PumiInterface by setting `smb_name` to `bar.smb` and `dmg_name` to `foo.dmg`

### Functionality

The main functionality is the expression of the mesh coordinates in terms of the
CAD parametric coordinates for entities on the boundary.  The conversion of
derivatives wrt xyz to derivatives wrt the parametric coordinates is also
supported (although accuracy varies depending on how the underlying CAD system
computes derivatives).  See the functions in `interface_geo.jl`.


The `PumiMesh` types also keep track of the CAD parametric coordinates for
all coordinate nodes.  In order to keep a consistent data structure, use
`update_coords` or `update_coordsXi` to update the coordinates of the entire
mesh, or `get/setCoords` and `get/setCoordsXi` to node-level granularity.
Do *not* call `apf.setPoint` directly.

Its is generally preferable to work with the CAD parametric coordinates, because
the mapping parametric -> Cartesian is faster and more accurate to compute,
but this only works if the underlying CAD system supports parametric coordinates
(`.dmg` models do not).

# Reference Counting and P-sequencing

It is possible to construct several Julia mesh objects that share an 
underlying `apf::Mesh`.  See the `copy_mesh()` function.  These meshes
may have different types of SBP operators and may be of different degree.

In order to do proper memory managment, the `apf::Field`s and `apf::Numbering`s
associated with each Julia mesh object are stored and reference counting
is used to determine when to free both the fields and numberings, as well
as the underlying mesh object.  The functions in `pde_fields.jl` should be
used to create new fields and numberings (do *not* call
`apf.createPackedField(mesh.m_ptr, ...)` or
`apf.createNumberingJ(mesh.m_ptr, ...)` directly.

By default, when Paraview files are written only fields which are both
associated with the given Julia mesh object and are properly display-able
in VTK format will be written to the Paraview file.  This second restriction
can be removed with the keyword argument `writeall=true` to `writeVisFiles`.
Note that this may substanially increase the file size.

# Utilities

A few C++ executables are built with PumiInterface that are useful as standalone tools.

 * `test_init` this executable mirrors the initialization process used when
               PumiInterface loads a mesh. Useful for experimenting with new
               functionality (before wrapping it in Julia) and for developing
               standalone minimal reproducable examples for filing bug reports.
 * `printGeoClassification`: prints all the geometric entities for a given
                             geometric model, as well as their parametric ranges
                             and whether they are periodic or not.  Useful for
                             quickly inspecting a CAD model.
 * `uniform_refine`: run's Pumi's uniform refinement + snapping + shape
                     correction.  Useful for doing refinement studies on
                     unstructured meshes

If invoked with no arguments, each of these executables prints a message
their arguments.



# Version History

 * v0.1: the old code, before Pumi switched to the CMake build system.  This version is no longer supported

 * v0.2: after the build system switch to CMake

 * v0.3: curvilinear elements (compatable with SBP versions broken_ticon and fixed_ticon)

 * v0.4: 2D linear mesh reverse mode of metrics

 * v0.5: 3D curvilinear metrics and their reverse mode

             Breaking changes:

               Curvilinear metrics are now the default.  interpolateMapping_rev
               and getVertCoords_rev are no longer exported and are hidden
               inside getAllCoordinateAndMetrics_rev.  An options dictionary
               entry can force the use of the linear metrics.
               `commit_coords()` now requires the options dictionary as an
               argument.

               Added default BC.  All faces of all elements now appear somewhere
               in mesh.interfaces, mesh.bndryfaces or mesh.bndries_local.  If
               the options dictionary does not specify a boundary condition for
               any geometric edges, an additional BC is created for them with
               name `defaultBC`.  See `getAllFaceData()`.

 * v0.6: fix bug in in 3D metric calculation, also introduce asserts to make sure
      face orientation can be determined when using periodic boundary conditions.
      Also, removes all global state from libfuncs1, so it is now possible to
      load multiple meshes simultaneously.

 * v0.7: change mesh constructors to be outer constructors, add submesh creation

 * v0.8: upgrade from Julia 0.4 to Julia 0.6
 * v0.9: add parallel reverse mode derivative and gmi interface
 * v0.10: add CAD Xi coordinate interface, derivative, and initial mesh
          adaptation support

[![Build Status](https://travis-ci.org/OptimalDesignLab/PumiInterface.jl.svg?branch=build_system)](https://travis-ci.org/OptimalDesignLab/PumiInterface.jl)
