
project (fk_mc)  #not sure whether this should be the same name of the executable, but I always see that "convention"
cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 14)

ADD_EXECUTABLE(fk_mc main.cpp some_functions.h observables.h)
link_directories(/usr/include)
target_link_libraries (fk_mc armadillo)
