# this file will build Pumi if it cannot be located by pkg-config

# install PkgFix if not present
if !isdir(joinpath(Pkg.dir(), "PkgFix"))
  Pkg.clone("PkgFix")
end

using PkgFix  # from now on, use PkgFix instead of Pkg for everything

include("./stringmatch.jl")
include("install_pumi.jl")

# package URLs and versions
global const MPI_URL = https://github.com/JuliaParallel/MPI.jl.git
global const MPI_VER = "v0.5.0"

global const SBP_URL = "https://github.com/OptimalDesignLab/SummationByParts.jl.git"
global const SBP_VER = "jcwork"


start_dir = pwd()

run(`./cleanup.sh`)

install_pumi = try run(`cmake --find-package -DNAME=SCOREC -DCOMPILER_ID=GNU -DLANGUAGE=CXX -DMODE=EXIST`) 
          false
          catch 
	  true 
	end
println("install_pumi = ", install_pumi)


install_mpi = try run(`which mpicxx`) 
          false
          catch 
	  true 
	end
println("install_mpi = ", install_mpi)

# install MPI itself
if install_mpi
  cmd_string = "./travis-install-mpi.sh"
  arg_str = "mpich3"
  run(`$cmd_string $arg_str`)
end

# install MPI.jl if needed
if PkgFix.installed("MPI") == nothing
  PkgFix.add(MPI_URL, MPI_VER)
  # because the julia package manager is completely stupid and tries its best to
  # be unusable, the clone command attempts to resolve dependencies, which fails
  # because the version of Compat installed is too old to satisfy the version
  # required by the most recent MPI
  # instead use git directly to install it
  #=
  start_dir2 = pwd()
  cd(Pkg.dir("PumiInterface"))
  run(`git clone https://github.com/JuliaParallel/MPI.jl.git`)
  run(`mv -v ./MPI.jl ./MPI`)
  cd("./MPI")
  run(`git checkout v0.5.0`)
  cd(start_dir2)
  Pkg.build("MPI")
  =#
end

if install_pumi  # did not find pumi
  installPumi()
end

# install non-metadata dependencies
pkg_dict = PkgFix.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "SummationByParts")
  PkgFix.add(SBP_URL, SBP_VER)
  #=
  Pkg.clone("https://github.com/OptimalDesignLab/SummationByParts.jl.git")
  Pkg.checkout("SummationByParts", "jcwork")
  Pkg.build("SummationByParts")
  =#
  # SBP installs ODLCommonTools
end

# now that all dependencies exist, install this package
cd("../src")
run(`./config.sh`)
run(`./makeinstall.sh`)

#=
# build the shared library
cd("../src")
run(`./build_shared.scorec.sh7`)
str5 = string("export LD_LIBRARY_PATH=", pwd(), ":\$LD_LIBRARY_PATH")

# get path to pumi library used to build libfuncs1.so
sonames = readall(`ldd ./libfuncs1.so`)
pumi_path = findWord(sonames, "/libmds.so")
println("pumi_path = ", pumi_path)
str6 = string("export LD_LIBRARY_PATH=", pumi_path, ":\$LD_LIBRARY_PATH")


f = open("evars.sh", "a+")
println(f, str5)
println(f, str6)
close(f)
=#
cd(start_dir)



#ldd ./libfuncs1.so | grep libmds.so | awk -F"=>" '{print $2}' | awk '{print $1}' |  sed 's%/[^/]*$%/%'
