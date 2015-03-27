# PUMI.jl
this fork of PUMI.jl is for my Finite Element Programming assignment.  Hopefully
some of these changes are useful and will be merged back into the master branch.

What this Library Is
This code provides a way to use PUMI from Julia by wrapping functions in PUMI's APF API.
The library is composed of several files.  funcs1.cc, funcs1.h form the C++ half of the
interface, PumiInterface.jl forms the Julia half. All the functions in PumiInterface.jl,
except for init and init2, mirror apf functions.  I refer to this as the low level interface.
  PumiInterface2.jl provides some higher level function by using the lower level functions.
these functions are more convienent but potentially less efficient depending on the use case.
The files a2.cc and a2.h contain the user supplied reordering function.  YOU MAY NOT USE THIS
CODE FOR ASSIGNMENT 4.  The files are provided to show how a user can add functions.  If you have
questions about how to interface your reordering function with Julia, let me know.

funcs1.h has some useful general comments in it.  funcs1.cc has more descriptive comments for each
function.  

Building this Library
Currently this code is in a module, not yet a package. You must build the library, and
import the module in your Julia code to access the functions in the module.  All the function
 names are exported.  PUMI must be available to link agains during the building of the library.

If you are working on a SCOREC machine (with access to the SCOREC shared file system), the library
can be built as follows:

cd location/of/your/repo
source ./use_julialib.sh
./build_shared.scorec.sh4

Now the library is build.  The file use_julialib.sh sets some environmental variables and loads
some modules needed to make PUMI accessible and to run Julia, and the second script performs the 
actual build.  A successful build produces a file called libfuncs1.so.  The library must be available
 at the runtime of your code.  The use_julialib.sh will ensure this is true.  You must be in the same
directory as libfuncs1.so when you source the script.

Using Julia
I have built Julia on catan.  Builting it was a hassle, so I recomment just using that build.  It can 
be found here: /users/creanj/build/julia_1_29_15


Using the Library
The fundemental idea behind the library is that pointers to objects can be passed to and from Julia,
 but not dereferenced by Julia (bad things would happen if your try, so please do not kick the bucket).
  Therefore, a function that in C++ would be called like this:

  m_ptr->getDownward(args...)

Gets called from Julia like this:

  getDownward(m_ptr, args...)

Where m_ptr is a pointer to a PUMI object.

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


As a final note, the julia function names mirror (or are extremely similar to) those in PUMI's APF, so
searching PumiInterface.jl for the name of APF functions should get you to the right place quickly. 
If you want to find out what function in funcs1.cc a Julia function calls, look at the first argument
 to the ccall function within the Julia funcs, which will have the forms "functionname"_name, and look up
 the value of this global variable at the top of PumiInterface.jl.  The value will be a string, which is
 the name of the C++ function.




