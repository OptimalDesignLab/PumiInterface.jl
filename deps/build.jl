# this file will build Pumi if it cannot be located by pkg-config
include("./stringmatch.jl")
include("install_pumi.jl")
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


if install_mpi
  cmd_string = "./travis-install-mpi.sh"
  arg_str = "mpich3"
  run(`$cmd_string $arg_str`)
  Pkg.clone("MPI")
  Pkg.build("MPI")
end

if install_pumi  # did not find pumi
  installPumi()
end

# install non-metadata dependencies
pkg_dict = Pkg.installed()  # get dictionary of installed package names to version numbers

if !haskey(pkg_dict, "SummationByParts")
  Pkg.clone("https://github.com/OptimalDesignLab/SummationByParts.jl.git")
  Pkg.build("SummationByParts")
  # SBP installs ODLCOmmonTools
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
