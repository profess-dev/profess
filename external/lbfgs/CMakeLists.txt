cmake_minimum_required(VERSION 3.1)
project(lbfgs LANGUAGES CXX Fortran)

if(NOT EXISTS ${CMAKE_CURRENT_LIST_DIR}/lbfgs_um)
    message(STATUS "Downloading Nocedal's LBFGS routine")
    file(DOWNLOAD http://users.iems.northwestern.edu/~nocedal/Software/lbfgs_um.tar.gz ${CMAKE_CURRENT_LIST_DIR}/lbfgs_um.tar.gz)
    execute_process(COMMAND gunzip ${CMAKE_CURRENT_LIST_DIR}/lbfgs_um.tar.gz)
    execute_process(COMMAND tar -xf ${CMAKE_CURRENT_LIST_DIR}/lbfgs_um.tar --directory ${CMAKE_CURRENT_LIST_DIR})
    execute_process(COMMAND rm ${CMAKE_CURRENT_LIST_DIR}/lbfgs_um.tar)
endif()

add_library(lbfgs lbfgs_um/lbfgs.f c_lbfgs.f90)
target_include_directories(lbfgs PUBLIC .)
