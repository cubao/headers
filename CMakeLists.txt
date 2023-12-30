cmake_minimum_required(VERSION 3.15...3.27)
project(cubao_headers LANGUAGES NONE)

# ---------------------------------------------------------------------------
# Only build tests by default if this is the top-level CMake project
# ---------------------------------------------------------------------------

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CBH_MASTER_PROJECT ON)
else()
  set(CBH_MASTER_PROJECT OFF)
endif()

option(CBH_CREATE_INSTALL_RULES "Create installation rules" ${CBH_MASTER_PROJECT})

# ---------------------------------------------------------------------------
# Do a release build if nothing was specified
# ---------------------------------------------------------------------------

if (CBH_MASTER_PROJECT AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "cubao_headers: setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# ---------------------------------------------------------------------------
# Installation rules
# ---------------------------------------------------------------------------

if(CBH_CREATE_INSTALL_RULES AND NOT CMAKE_SKIP_INSTALL_RULES)
  # Silence warning in GNUInstallDirs due to no enabled languages
  set(CMAKE_INSTALL_LIBDIR "")
  include(GNUInstallDirs)
  set(CBH_INSTALL_DATADIR "${CMAKE_INSTALL_DATADIR}/cubao_headers"
    CACHE PATH "Installation path for read-only architecture-independent cubao_headers data files")

  # Normally these would be configurable by the user, but we can't allow that
  # because the lookup paths are hard-coded in 'cmake/cubao_headers-config.cmake'
  set(CMAKE_INSTALL_INCLUDEDIR "${CBH_INSTALL_DATADIR}/include")
  install(
    DIRECTORY include/
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
endif()