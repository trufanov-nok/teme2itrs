cmake_minimum_required(VERSION 3.5)

project(teme2itrs LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

######### FIND LIBS ###########

find_package(Eigen3 REQUIRED)
message( " - Eigen3: ${EIGEN3_INCLUDE_DIR}" )


######### ERFA BEGIN ##########
FIND_PATH(
    ERFA_INCLUDE_DIR erfa.h
    PATHS /usr/local/include /usr/include
    DOC "Path to ERFA headers."
)

IF(NOT ERFA_INCLUDE_DIR)
MESSAGE(
    FATAL_ERROR
    "Could not find ERFA library headers.\n"
    "You may need to install a package named liberfa-dev."
)
ENDIF()

FIND_LIBRARY(
    ERFA_LIBRARY NAMES erfa
    PATHS /usr/local/lib /usr/lib
    DOC "Path to ERFA library."
)

IF(NOT ERFA_LIBRARY)
MESSAGE(
    FATAL_ERROR
    "Could not find ERFA library.\n"
    "You may need to install a package named liberfa-dev."
)
ELSE()
    message( " - ERFA: ${ERFA_LIBRARY}" )
ENDIF()
######### ERFA END ###########

######### FIND LIBS END ###########

add_executable(teme2itrs main.cpp
    teme2itrs.h teme2itrs.cpp)
target_link_libraries(teme2itrs Eigen3::Eigen ${ERFA_LIBRARY})
add_dependencies(teme2itrs Eigen3::Eigen)


include(GNUInstallDirs)
install(TARGETS teme2itrs
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
