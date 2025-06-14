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

# ========================================
# ROOT Dictionary Generation
# ========================================
set(DICT_HEADERS
    ${INCLUDE_DIR}/CWignerSource.h
    ${INCLUDE_DIR}/CWignerUtils.h
)

ROOT_GENERATE_DICTIONARY(G__WignerUtils
    ${DICT_HEADERS}
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/WignerLinkDef.h
)

# ========================================
# Detect ROOT lib directory for RPATH
# ========================================
execute_process(
  COMMAND root-config --libdir
  OUTPUT_VARIABLE ROOT_LIBRARY_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# ========================================
# Shared library: WignerUtils (with dictionary)
# ========================================
add_library(WignerUtils SHARED ${SOURCES} G__WignerUtils.cxx)
target_include_directories(WignerUtils PRIVATE ${INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(WignerUtils PRIVATE ${ROOT_LIBRARIES})
set_target_properties(WignerUtils PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    INSTALL_RPATH "@loader_path/../lib;${ROOT_LIBRARY_DIR}"
    BUILD_RPATH "@loader_path/../lib;${ROOT_LIBRARY_DIR}"
)

# Install main library and headers
install(TARGETS WignerUtils LIBRARY DESTINATION lib)
install(DIRECTORY ${INCLUDE_DIR}/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libWignerUtils_rdict.pcm
    DESTINATION lib
)

# ========================================
# Executable: wigneroot
# ========================================
add_executable(wigneroot ${SOURCE_DIR}/wigneroot.cpp)
target_include_directories(wigneroot PRIVATE ${INCLUDE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(wigneroot PRIVATE WignerUtils ${ROOT_LIBRARIES})
set_target_properties(wigneroot PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    INSTALL_RPATH "@executable_path/../lib;${ROOT_LIBRARY_DIR}"
    BUILD_RPATH "@executable_path/../lib;${ROOT_LIBRARY_DIR}"
)

# Install the executable
install(TARGETS wigneroot RUNTIME DESTINATION bin)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libWignerUtils_rdict.pcm
    DESTINATION lib
)
