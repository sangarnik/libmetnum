# - Config file for the Armadillo package
# It defines the following variables
#  ARMADILLO_INCLUDE_DIRS - include directories for Armadillo
#  ARMADILLO_LIBRARY_DIRS - library directories for Armadillo (normally not used!)
#  ARMADILLO_LIBRARIES    - libraries to link against

# Tell the user project where to find our headers and libraries
set(ARMADILLO_INCLUDE_DIRS "/home/nikos/Documents/Tareitas/1er_semestre/analisis_numerico/practica/TP2_EquationSolver/lib/armadillo-12.2.0/build/tmp/include")
set(ARMADILLO_LIBRARY_DIRS "/home/nikos/Documents/Tareitas/1er_semestre/analisis_numerico/practica/TP2_EquationSolver/lib/armadillo-12.2.0/build")

# Our library dependencies (contains definitions for IMPORTED targets)
include("/home/nikos/Documents/Tareitas/1er_semestre/analisis_numerico/practica/TP2_EquationSolver/lib/armadillo-12.2.0/build/ArmadilloLibraryDepends.cmake")

# These are IMPORTED targets created by ArmadilloLibraryDepends.cmake
set(ARMADILLO_LIBRARIES armadillo)

