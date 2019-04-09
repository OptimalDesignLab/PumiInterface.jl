# reading and writing coordinates to a binary file

const COORDSFILE_MAGIC_BYTE = Int32(26)
const COORDSFILE_VERSION = Int32(1)  # file format version

"""
  Writes the xyz (and xi coordinates if supported) to a binary file.

  The file format is:

   * 4 bytes containing magic number
   * 4 byte integer containing file format version number
   * 4 byte integer specifying if xi coordinates are present
   * 8 byte integer specifying the length of the xyz coordinate vector
   * if xi coordinates are present, and 8 byte integer specifying the length
     of the vector
   * the xyz coordinate vector (of Float64s)
   * if xi coordinates are present, the xi coordinate vector (of Float64s)

  **Inputs**

   * mesh
   * fname: the file name, including extension.  An underscore and MPI rank
            will be appended to the file name (before the extension), and
            one file will be written per process.
"""
function writeCoordsBinary(mesh::PumiMesh, fname::String)

  f = open(get_parallel_fname(fname, mesh.myrank), "w")

  writeMagicByte(mesh, f)
  writeFileVersion(mesh, f)
  writeFileType(mesh, f)
  writeCoords(mesh, f)

  close(f)

  return nothing
end


mutable struct CoordsFileData
  fname::String
  version::Int32
  ftype::Int32
  numCoordDofs::Int64
  numXiDofs::Int64
end

function CoordsFileData(fname::String)
  return CoordsFileData(fname, 0, 0, 0, 0)
end

"""
  Reads the file written by [`writeCoordsBinary`](@ref) and updates the mesh
  object
"""
function readCoordsBinary(mesh::PumiMesh, sbp, opts, fname::String)

  #TODO: write number of MPI ranks to file
  f = open(get_parallel_fname(fname, mesh.myrank), "r")

  obj = CoordsFileData(fname)
  checkMagicByte(mesh, f, obj)
  readFileVersion(mesh, f, obj)
  readFileType(mesh, f, obj)
  readCoords(mesh, sbp, opts, f, obj)

  close(f)

  return nothing
end



"""
  Returns the integer describing if the file contains xyz coordinates and
  the xi coordinates or only the xyz coordinates.
"""
function getCoordsFileType(mesh)

  return Int32(mesh.geoNums.can_eval)
end


#------------------------------------------------------------------------------
# Writer

"""
  Write a magic number to the first 4 bytes of the file
"""
function writeMagicByte(mesh::PumiMesh, f::IO)

  write(f, COORDSFILE_MAGIC_BYTE)
end

"""
  Writes the file format version number
"""
function writeFileVersion(mesh::PumiMesh, f::IO)

  write(f, COORDSFILE_VERSION)
end


"""
  Write some additional information to the file: whether xi coordinates are
  present, and the length of the xyz (and xi, if present) vectors.
"""
function writeFileType(mesh::PumiMesh, f::IO)

  ftype = getCoordsFileType(mesh)
  can_eval = Bool(ftype)

  write(f, ftype)
  # write number of xyz coordinates, and xi coordinate if available
  write(f, Int64(mesh.geoNums.numCoordDofs))
  if can_eval
    write(f, mesh.geoNums.numXiDofs)
  end
end

function writeCoords(mesh::PumiMesh, f::IO)

  xvec = getXCoords(mesh)
  @assert eltype(xvec) == Float64
  write(f, xvec)


  if Bool(getCoordsFileType(mesh))
    xivec = getXiCoords(mesh)
    @assert eltype(xivec) == Float64
    write(f, xivec)
  end
end

#------------------------------------------------------------------------------
# Reader

"""
  Asserts if the file magic byte is not correct
"""
function checkMagicByte(mesh::PumiMesh, f::IO, obj::CoordsFileData)

  val = read(f, Int32)
  if val != COORDSFILE_MAGIC_BYTE
    error("Possible file corruption in coordinate file $(obj.fname)")
  end

  return nothing
end

"""
  Gets the file version number
"""
function readFileVersion(mesh::PumiMesh, f::IO, obj::CoordsFileData)

  val = read(f, Int32)

  if val > COORDSFILE_VERSION
    error("Coordinate file version is from the future: version $val > $COORDSFILE_VERSION")
  end

  obj.version = val

  return nothing
end

"""
  Gets the file type and length of arrays
"""
function readFileType(mesh::PumiMesh, f::IO, obj::CoordsFileData)

  obj.ftype = read(f, Int32)
  can_eval = Bool(obj.ftype)

  obj.numCoordDofs = read(f, Int64)
  if can_eval
    obj.numXiDofs = read(f, Int64)
  end

  return nothing
end

"""
  Read the xyz (and possibly xi) coordinates from the file and update the mesh.
  If the file does not contain xi coordinates and the mesh does, the
  xi coordinates will be calculated from the xyz.

  **Inputs**

   * mesh
   * f: IO object
   * obj: CoordFileData
"""
function readCoords(mesh::PumiMesh, sbp, opts, f::IO, obj::CoordsFileData)

  # read xyz coordinates
  if obj.numCoordDofs != mesh.geoNums.numCoordDofs
    error("Coordinate file $(obj.fname) has $(obj.numCoordDofs) xyz coordinates, which is not equal to $(mesh.geoNums.numCoordDofs), the number of mesh xyz dofs")
  end
  xvec = zeros(Float64, obj.numCoordDofs)
  read!(f, xvec)

  # read xi coordinates
  f_has_xi = Bool(obj.ftype)
  if f_has_xi && mesh.geoNums.can_eval
    if obj.numXiDofs != mesh.geoNums.numXiDofs
      error("Coordinate file $(obj.fname) has $(obj.numXiDofs) xi coordinates, which is not equal to $(mesh.geoNums.numXiDofs), the number of mesh xi dofs")
    end

    xivec = zeros(Float64, obj.numXiDofs)
    read!(f, xivec)
  elseif mesh.geoNums.can_eval && !f_has_xi
    if mesh.myrank == 0
      println(STDERR, "Warning: coordinate file $(obj.fname) does not contain Xi coordinates, yet the mesh supports geometry coordinates")
    end
  elseif f_has_xi && !mesh.geoNums.can_eval
    if mesh.myrank == 0
      println(STDERR, "Warning: coordinate file $(obj.fname) contains Xi coordinates, but the the mesh does not support geometry coordinates")
    end
  end


  # update xi coords, either from the file or by calculating from xyz
  if mesh.geoNums.can_eval && f_has_xi
    update_coordsXi(mesh, sbp, opts, xivec)
  end

  if !f_has_xi && mesh.geoNums.can_eval
    if mesh.myrank == 0
      println("Computing xi coordinates from xyz coordinates")
    end
    xivec = zeros(Float64, mesh.geoNums.numXiDofs)
    coords_xyzToXi(mesh, xvec, xivec)
    update_coordsXi(mesh, sbp, opts, xivec)
  end

  # update xyz coordinates from file
  # The xyz coordinates computed from update_coordsXi are typically only
  # accurate to 1e-8, so set the xyz coordinates here for better accuracy
  update_coords(mesh, sbp, opts, xvec, snap=false)

  return nothing
end

