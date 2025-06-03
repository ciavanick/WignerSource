cmake_minimum_required(VERSION 3.10)
project(WignerUtils)

# Use .so extension on macOS for ROOT compatibility
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT
find_package(ROOT REQUIRED)

# Set paths
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Source files
set(SOURCES
    ${SOURCE_DIR}/CWignerSource.cpp
    ${SOURCE_DIR}/CWignerUtils.cpp
)

# Main library
add_library(WignerUtils SHARED ${SOURCES})
target_include_directories(WignerUtils PRIVATE ${INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(WignerUtils PRIVATE ${ROOT_LIBRARIES})
set_target_properties(WignerUtils PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Install main library and headers
install(TARGETS WignerUtils LIBRARY DESTINATION lib)
install(DIRECTORY ${INCLUDE_DIR}/ DESTINATION include FILES_MATCHING PATTERN "*.h")

# =======================
# ROOT DICTIONARY SECTION
# =======================

# Dictionary headers
set(DICT_HEADERS
    ${INCLUDE_DIR}/CWignerSource.h
    ${INCLUDE_DIR}/CWignerUtils.h
)

# Generate dictionary (note MODULE is needed for .pcm)
ROOT_GENERATE_DICTIONARY(G__WignerUtils
    ${DICT_HEADERS}
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/WignerLinkDef.h
    MODULE WignerUtils
)

# Dictionary library
add_library(WignerDict SHARED G__WignerUtils.cxx)
target_include_directories(WignerDict PRIVATE ${INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(WignerDict PRIVATE WignerUtils ${ROOT_LIBRARIES})
set_target_properties(WignerDict PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Install dictionary lib and PCM
install(TARGETS WignerDict LIBRARY DESTINATION lib)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libWignerUtils_rdict.pcm
    DESTINATION lib
)
