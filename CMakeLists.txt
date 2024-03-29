cmake_minimum_required(VERSION 3.1)
project(profess LANGUAGES CXX Fortran)

set(BUILD_SHARED_LIBS ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#add_compile_options(-g -Wall -Wextra -pedantic $<$<COMPILE_LANGUAGE:CXX>:-Wfloat-conversion>)
add_compile_options(-g -Ofast -march=native)

# profess library
add_library(profess
    source/ions.cpp
    source/ion_ion.cpp
    source/system.cpp
    source/tools.cpp
    source/minimizers.cpp
    source/functionals/functional.cpp 
    source/functionals/hartree.cpp
    source/functionals/huang_carter.cpp
    source/functionals/ion_electron.cpp
    source/functionals/libxc.cpp
    source/functionals/luo_karasiev_trickey.cpp
    source/functionals/kinetic_class_a.cpp
    source/functionals/kinetic_gga.cpp
    source/functionals/perdew_burke_ernzerhof.cpp
    source/functionals/perdew_zunger.cpp
    source/functionals/perrot.cpp
    source/functionals/smargiassi_madden.cpp
    source/functionals/thomas_fermi.cpp
    source/functionals/thomas_fermi_weizsaecker.cpp
    source/units.cpp
    source/functionals/wang_govind_carter.cpp
    source/functionals/wang_govind_carter_1999_i.cpp
    source/functionals/wang_teter.cpp
    source/functionals/weizsaecker.cpp
)
target_include_directories(profess PUBLIC include include/functionals)

# deft
add_compile_definitions(POCKETFFT_NO_MULTITHREADING)
add_subdirectory(external/deft)
target_link_libraries(profess PUBLIC deft)

# external/lbfgs
add_subdirectory(external/lbfgs)
target_link_libraries(profess PRIVATE lbfgs)

# external/real_space_electrostatic_sum
add_subdirectory(external/real-space-electrostatic-sum)
target_link_libraries(profess PRIVATE real_space_electrostatic_sum)

# external/libkedf
set(JAVA   "OFF" CACHE STRING "")
set(OPENMP "OFF" CACHE STRING "")
set(OPENCL "OFF" CACHE STRING "")
set(MPI    "OFF" CACHE STRING "")
add_subdirectory(external/libKEDF)
target_compile_options(KEDF PRIVATE -Wno-error)
target_compile_options(KEDFClient PRIVATE -Wno-error) # delete eventually
target_include_directories(profess PUBLIC external/libKEDF/include)
target_link_libraries(profess PUBLIC KEDF)

# libxc
target_include_directories(profess PUBLIC ${LIBXC_INCLUDE_DIR})
target_link_directories(profess PUBLIC ${LIBXC_LIBRARY_DIR})
target_link_libraries(profess PUBLIC xc)

# profess-pybind
pybind11_add_module(profess-pybind MODULE pyrofess/pyrofess.cpp)
target_link_libraries(profess-pybind PUBLIC profess)
set_target_properties(profess-pybind PROPERTIES OUTPUT_NAME profess)

# profess-pybind extras
pybind11_add_module(extras MODULE pyrofess/extras.cpp)
target_link_libraries(extras PUBLIC profess)
