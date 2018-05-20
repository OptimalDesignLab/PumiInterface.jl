function installPumi()
  #pumi_version = "e4eabf5"
  #pumi_version = "d1837c5936c28d6ef09abd02c82ea2a4ea9b6f55"
#  pumi_version = "7efdf2319d13ce86994930e867d226ea421c1717"
#  pumi_version = "63fa247358e46352c0d500b1f1735c82447dfddd"
#   pumi_version = "develop"
  pumi_version = "063ad65"
#  pumi_version = "HEAD"

  start_dir = pwd()  # record where we started
  run(`./cleanup.sh`)
#  if isdir("./core")
#    println("deleting existing Core repo in /deps")
#    rm("./core", recursive=true)
#  end
  if !isdir("./core")
    run(`./download.sh`)
  end
  cd("./core")
#  run(`git pull`)
  run(`git checkout $pumi_version`)
  mkdir("./build")
  cd("./build")
  mkdir("./install")  # install directory
  run(`../../config.sh`)
  run(`make -j 4`)
  run(`make install`)
  str1 = joinpath( pwd(), "install/lib")
  str3 = joinpath(str1, "pkgconfig")

  # update ENV in preparatoin for building files in /src
  if haskey(ENV, "LD_LIBRARY_PATH")
    ld_path = ENV["LD_LIBRARY_PATH"]
    ld_path = string(str1, ":", ld_path)
    ENV["LD_LIBRARY_PATH"] = ld_path
  else
    ENV["LD_LIBRARY_PATH"] = str1
  end

  # for backwards compatibility with the pkg-config days
  if haskey(ENV, "PKG_CONFIG_PATH")
    pkg_path = ENV["PKG_CONFIG_PATH"]
    pkg_path = string(str3,":", pkg_path)
    ENV["PKG_CONFIG_PATH"] = pkg_path
  else
    ENV["PKG_CONFIG_PATH"] = str3
  end

  scorec_prefix = joinpath(pwd(), "install")
  println("scorec prefix = ", scorec_prefix)
  ENV["SCOREC_PREFIX"] = scorec_prefix


  cd(start_dir)
  return nothing
end


