find_package(ITK REQUIRED)
include(${ITK_USE_FILE})


find_package(CPPUtils REQUIRED)
include_directories(${CPPUtils_INCLUDE})

include_directories(${NRRD_INCLUDE_DIR})

file(GLOB_RECURSE SOURCES "*.cpp")

add_executable(BodyOutlineDelineation ${SOURCES})

message(STATUS "CPPUtils_LIBDIR: ")
message(STATUS "${CPPUtils_LIBDIR}")

message(STATUS "ITK_LIBRARIES: ")
message(STATUS "${ITK_LIBRARIES}")

target_link_libraries (BodyOutlineDelineation "${CPPUtils_LIBDIR}/libCPPUtils.a" ${ITK_LIBRARIES} )

#add_dependencies(BodyOutlineDelineation CPPUtils)
