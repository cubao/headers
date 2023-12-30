include_guard(GLOBAL)

# ---------------------------------------------------------------------------
# Helper function to handle undefined CPython API symbols on macOS
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Create shared/static library targets for nanobind's non-templated core
# ---------------------------------------------------------------------------

function (cubao_headers_build_library TARGET_NAME)
  if (TARGET ${TARGET_NAME})
    return()
  endif()

  if (TARGET_NAME MATCHES "-static")
    set (TARGET_TYPE STATIC)
  else()
    set (TARGET_TYPE SHARED)
  endif()

  add_library(${TARGET_NAME} ${TARGET_TYPE}
    EXCLUDE_FROM_ALL
    # TODO
    ${NB_DIR}/include/nanobind/make_iterator.h
  )

  if (TARGET_TYPE STREQUAL "SHARED")
    cubao_headers_link_options(${TARGET_NAME})
    target_compile_definitions(${TARGET_NAME} PRIVATE -DNB_BUILD)
    target_compile_definitions(${TARGET_NAME} PUBLIC -DNB_SHARED)
    cubao_headers_lto(${TARGET_NAME})

    cubao_headers_strip(${TARGET_NAME})
  elseif(NOT WIN32 AND NOT APPLE)
    target_compile_options(${TARGET_NAME} PUBLIC $<${NB_OPT_SIZE}:-ffunction-sections -fdata-sections>)
    target_link_options(${TARGET_NAME} PUBLIC $<${NB_OPT_SIZE}:-Wl,--gc-sections>)
  endif()

  set_target_properties(${TARGET_NAME} PROPERTIES
    POSITION_INDEPENDENT_CODE ON)

  if (MSVC)
    # Do not complain about vsnprintf
    target_compile_definitions(${TARGET_NAME} PRIVATE -D_CRT_SECURE_NO_WARNINGS)
  else()
    # Generally needed to handle type punning in Python code
    target_compile_options(${TARGET_NAME} PRIVATE -fno-strict-aliasing)
  endif()

  if (WIN32)
    if (${TARGET_NAME} MATCHES "abi3")
      target_link_libraries(${TARGET_NAME} PUBLIC Python::SABIModule)
    else()
      target_link_libraries(${TARGET_NAME} PUBLIC Python::Module)
    endif()
  endif()

  # Nanobind performs many assertion checks -- detailed error messages aren't
  # included in Release/MinSizeRel modes
  target_compile_definitions(${TARGET_NAME} PRIVATE
    $<${NB_OPT_SIZE}:NB_COMPACT_ASSERTIONS>)

  target_include_directories(${TARGET_NAME} PUBLIC
    ${Python_INCLUDE_DIRS}
    ${NB_DIR}/include)
  cubao_headers_set_visibility(${TARGET_NAME})
endfunction()

# ---------------------------------------------------------------------------
# Define a convenience function for creating nanobind targets
# ---------------------------------------------------------------------------

  add_library(${name} MODULE ${ARG_UNPARSED_ARGUMENTS})

  cubao_headers_compile_options(${name})
  cubao_headers_link_options(${name})
  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE CXX)

  if (ARG_NB_SHARED AND ARG_NB_STATIC)
    message(FATAL_ERROR "NB_SHARED and NB_STATIC cannot be specified at the same time!")
  elseif (NOT ARG_NB_SHARED)
    set(ARG_NB_STATIC TRUE)
  endif()

  # Stable ABI builds require CPython >= 3.12 and Python::SABIModule
  if ((Python_VERSION VERSION_LESS 3.12) OR
      (NOT Python_INTERPRETER_ID STREQUAL "Python") OR
      (NOT TARGET Python::SABIModule))
    set(ARG_STABLE_ABI FALSE)
  endif()

  set(libname "nanobind")
  if (ARG_NB_STATIC)
    set(libname "${libname}-static")
  endif()

  if (ARG_STABLE_ABI)
    set(libname "${libname}-abi3")
  endif()

  if (ARG_NB_DOMAIN AND ARG_NB_SHARED)
    set(libname ${libname}-${ARG_NB_DOMAIN})
  endif()

  cubao_headers_build_library(${libname})

  target_link_libraries(${name} PRIVATE ${libname})
  cubao_headers_set_visibility(${name})
endfunction()
