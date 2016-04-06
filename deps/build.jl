# this file will build Pumi if it cannot be located by pkg-config
include("./stringmatch.jl")

install_pumi = try run(`pkg-config --exists --libs libmds`) 
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



start_dir = pwd()  # record where we started


if install_mpi
  cmd_string = "./travis-install-mpi.sh"
  arg_str = "mpich3"
  run(`$cmd_string $arg_str`)
end

pumi_version = "e4eabf5"
pumi_version = "HEAD"
if install_pumi  # did not find pumi
  if isdir("./core")
    println("deleting existing Core repo in /deps")
    rm("./core", recursive=true)
  end
  run(`git clone https://github.com/SCOREC/core.git core`) 
  cd("./core")
  run(`git pull`)
  run(`git checkout $pumi_version`)
  mkdir("./build")
  cd("./build")
  mkdir("./install")  # install directory
  run(`../../config.sh`)
  run(`make -j 4`)
  run(`make install`)
  str1 = joinpath( pwd(), "install/lib")
  str2 = string("export LD_LIBRARY_PATH=", str1, ":\$LD_LIBRARY_PATH")
  str3 = joinpath(str1, "pkgconfig")
  str4 = string("export PKG_CONFIG_PATH=", str3, ":\$PKG_CONFIG_PATH")

  # write to file 
  fname = "evars.sh"
  f = open(fname, "a+")

  println(f, str2)
  println(f, str4)

  close(f)

  # update ENV here 

  if haskey(ENV, "LD_LIBRARY_PATH")
    ld_path = ENV["LD_LIBRARY_PATH"]
    ld_path = string(str1, ":", ld_path)
    ENV["LD_LIBRARY_PATH"] = ld_path
  else
    ENV["LD_LIBRARY_PATH"] = str1
  end

  if haskey(ENV, "PKG_CONFIG_PATH")
    pkg_path = ENV["PKG_CONFIG_PATH"]
    pkg_path = string(str3,":", pkg_path)
    ENV["PKG_CONFIG_PATH"] = pkg_path
  else
    ENV["PKG_CONFIG_PATH"] = str3
  end


  cd(start_dir)

end

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

cd(start_dir)



#ldd ./libfuncs1.so | grep libmds.so | awk -F"=>" '{print $2}' | awk '{print $1}' |  sed 's%/[^/]*$%/%'
