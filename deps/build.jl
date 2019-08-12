# this file will build Pumi if it cannot be located by pkg-config

# install PkgFix if not present
if !isdir(joinpath(Pkg.dir(), "PkgFix"))
  Pkg.clone("https://github.com/OptimalDesignLab/PkgFix.jl.git")
end
start_dir=pwd()
cd(Pkg.dir("PkgFix"))
run(`git checkout upgrade_0.6`)
cd(start_dir)

using PkgFix  # from now on, use PkgFix instead of Pkg for everything

include("./stringmatch.jl")
include("install_pumi.jl")

# package URLs and versions

global const ARRAYVIEWS_URL = "https://github.com/JaredCrean2/ArrayViews.jl.git"
global const ARRAYVIEWS_VER = "work"

global const MPI_URL = "https://github.com/JuliaParallel/MPI.jl.git"
global const MPI_VER = "v0.5.0"

global const SBP_URL = "https://github.com/OptimalDesignLab/SummationByParts.jl.git"
global const SBP_VER = "jcwork"


start_dir = pwd()

run(`./cleanup.sh`)

install_pumi = !haskey(ENV, "SCOREC_PREFIX") || try run(`cmake --find-package -DNAME=SCOREC -DCOMPILER_ID=GNU -DLANGUAGE=CXX -DMODE=EXIST`) 
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

#=
# install MPI.jl if needed
if PkgFix.installed("MPI") == nothing
  PkgFix.add(MPI_URL, branch_ish=MPI_VER)
end
=#

if install_pumi  # did not find pumi
  installPumi()
end

# install non-metadata dependencies
pkg_dict = PkgFix.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "SummationByParts")
  PkgFix.add(SBP_URL, branch_ish=SBP_VER)
  # SBP installs ODLCommonTools
end

if !haskey(pkg_dict, "ArrayViews")
  PkgFix.add(ARRAYVIEWS_URL, branch_ish=ARRAYVIEWS_VER)
end

# now that all dependencies exist, install this package
cd("../src")
run(`./config.sh`)
run(`./makeinstall.sh`)
cd(start_dir)



#ldd ./libfuncs1.so | grep libmds.so | awk -F"=>" '{print $2}' | awk '{print $1}' |  sed 's%/[^/]*$%/%'
