# PUMI.jl


## What this Library Is

This code provides a way to use PUMI from Julia by wrapping functions in PUMIs APF API.
The library is composed of several files.  funcs1.cc, funcs1.h form the C++ half of the
interface, PumiInterface.jl forms the Julia half. All the functions in PumiInterface.jl,
except for `init` and `init2`, mirror apf functions.  I refer to this as the low level interface.
  PumiInterface2.jl provides some higher level function by using the lower level functions.
These functions are more convenient but potentially less efficient depending on the use case.
The files a2.cc and a2.h contain the user supplied reordering function.  YOU MAY NOT USE THIS
CODE FOR ASSIGNMENT 4.  The files are provided to show how a user can add functions.  If you have
questions about how to interface your reordering function with Julia, let me know.

funcs1.h has some useful general comments in it.  funcs1.cc has more descriptive comments for each
function.  

## Installation

Now that the library is a Julia package, you can get is using Pkg.clone(url) followed by Pkg.build().
The build script will check for all needed components and attempt to build them itself if it does not find
them.  A C++11 compiler is required.

The next two subsections describe the inner workings of the build system. 
See the following to sections for the standard installation procedure.

### Install MPI (if needed)
Building Pumi and building the shared library that links to Pumi requires `mpicxx`.  If `mpicxx` is not found,  the build script will attempt to install a Debian package of MPICH3.  It will ask for root permission to install MPICH.

### Build Pumi (if needed)
The build script
attempts to locate Pumi using `pkg-config`. Pumi must also be locatable in your `LD_LIBRARY_PATH`.
If not found by pkg-config, it will attempt to build Pumi locally.  The build requires an MPI implementation (which
 should be present according to the above process), and Cmake 2.8.6.  It will install to `deps/core/build/install`.  The file `deps/core/build/evars.sh`  indicates the environmental variables needed to build 
the shared library using this installation of Pumi.

### Build Shared Library
At this point, MPI and Pumi should be be present, so the library itself can be built.  The build 
script executes `build_shared.scorec.sh7` in the `src/` directory to build the shared library, and 
creates a file `src/evars.sh`, which exports the environmental variables needed to use the library 
at runtime.  Note that because the library dynamically links to Pumi, Julia needs to be able to find 
both the shared library and Pumi at runtime.  The `evars.sh` adds both the installation of Pumi used 
to build the shared library and the shared library itself to `LD_LIBRARY_PATH`.  This script is 
provided entirely for user convenience.  The package never uses it.  You must  `source` this script 
at runtime to use the library.  

Thats it, now you are done.

Note Pumi is built without some of its dependencies, including zoltan, so certain features such as load balancing will not work.

## Scorec Installation

If you are working on a SCOREC machine (with access to the SCOREC shared file system), load the Pumi module and build the library as follows:

```
Pkg.clone(url)
cd location/of/repo/src
source ./use_julialib.sh
./build_shared.scorec.sh6
```

Now the library is build, using the Scorec installation of Pumi, which is preferable to building Pumi locally.  The file use_julialib.sh sets some environmental variables and loads
some modules needed to make PUMI accessible and to run Julia, and the second script performs the 
actual build.  A successful build produces a file called libfuncs1.so.  The library must be available
 at the runtime of your code.  The use_julialib.sh will set the environmental variables needed to make this true.  You must be in the same
directory as libfuncs1.so when you source the script.

## Non Scorec Installation
If you are not working on a SCOREC machine, or do not want to use the SCOREC 
installation, you can use the following procedure.

1. Unload any existing Pumi installation (if necessary)
2. Load Cmake 2.8.6 or newer
3. Load an MPI implimentation
4. In julia, `Pkg.build("PumiInterface")`

Step 4 will build Pumi and then build the interface library. To *use* the 
library:

```
cd location/of/repo/src  # this step is required
source ./use_julialib2.sh
```

The shell script will set the environmental variables needed for the Julia 
code to find the library.
If you want to rebuild the library manually:

```
source location/of/repo/deps/core/build/evars.sh
cd location/of/repo/src
./build_shared.scorec.sh7
```

This will typically not be necessary because this is done by Step 4 above.

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
 to the ccall function within the Julia funcs, which will have the forms "functionname"_name, and look up
 the value of this global variable at the top of PumiInterface.jl.  The value will be a string, which is
 the name of the C++ function.

## Types of Meshes
`PdePumiInterface` has two kinds of meshes, Continuous Galerkin (CG) and
Discontinuous Galerkin (DG), in two dimensions, 2D and 3D.  Currently, 
2D CG, 2D DG and 3D DG are implemented.

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


## Periodic Boundary Conditions
DG meshes support periodic boundary conditions.  The .smb file contains 
information on which MeshEntities are periodic with which other 
MeshEntities.  These entities will automatically be identified as periodic 
and treated accordingly.  It is an error to assign a boundary condition to a 
geometric entity that has a periodic MeshEntity on it, and an exception will 
be thrown.

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
