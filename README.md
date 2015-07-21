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
them.

Building requires locating Pumi and building a shared library that links to Pumi, using `mpicxx`.  If `mpicxx` is not found,  the build script will attempt to install a Debian package of MPICH3.

The build script
attempts to locate Pumi using `pkg-config`. Pumi must also be locatable in your `LD_LIBRARY_PATH`.
If not found, it will attempt to build Pumi locally.  The build requires an MPI implementation (which
 should be present according to the above process), and Cmake 2.8.6.  It will install to `deps/core/build/install`.  The file `deps/core/build/evars.sh`  indicates the environmental variables needed to build 
the shared library using this installation of Pumi.

At this point, MPI and Pumi should bot be present, so the library itself can be built.  The build 
script executes `build_shared.scorec.sh7` in the `src/` directory to build the shared library, and 
creates a file `src/evars.sh`, which exports the environmental variables needed to use the library 
at runtime.  Note that because the library dynamically links to Pumi, Julia needs to be able to find 
both the shared library and Pumi at runtime.  The `evars.sh` adds both the installation of Pumi used 
to build the shared library and the shared library itself to `LD_LIBRARY_PATH`.  This script is 
provided entirely for user convenience.  The package never uses it.  You must  `source` this script 
at runtime to use the library.  

Thats it, now you are done.

## Scorec Installation

If you are working on a SCOREC machine (with access to the SCOREC shared file system), the library
can be built as follows:

```
Pkg.clone(url)
cd location/of/repo/src
source ./use_julialib.sh
./build_shared.scorec.sh6
```

Now the library is build, using the Scorec installation of Pumi, which is preferable to building Pumi locally.  The file use_julialib.sh sets some environmental variables and loads
some modules needed to make PUMI accessible and to run Julia, and the second script performs the 
actual build.  A successful build produces a file called libfuncs1.so.  The library must be available
 at the runtime of your code.  The use_julialib.sh will ensure this is true.  You must be in the same
directory as libfuncs1.so when you source the script.


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




