set(CUBAO_UTILS_CMAKE_DIRNAME ${CMAKE_CURRENT_LIST_DIR} CACHE INTERNAL "")

macro(activate_common_configuration)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    set(CMAKE_CXX_STANDARD 14)
endmacro()

macro(auto_build_type_and_compile_flags)
    if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
        set(CMAKE_BUILD_TYPE
            "Release"
            CACHE STRING "" FORCE)
        message(STATUS "Set build type to default: ${CMAKE_BUILD_TYPE}")
    else()
        message(STATUS "Your build type: ${CMAKE_BUILD_TYPE}")
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -ggdb")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -ggdb")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    endif()
endmacro()

macro(setup_git_branch)
    execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

macro(setup_git_commit_hash)
    execute_process(
        COMMAND git log -1 --format=%h --abbrev=8
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

macro(setup_git_commit_count)
    execute_process(
        COMMAND git rev-list --count HEAD
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_COUNT
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

macro(setup_git_commit_date)
    execute_process(
        COMMAND bash "-c" "git log -1 --date='format:%Y/%m/%d %H:%M:%S' --format='%cd'"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_DATE
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

macro(setup_git_diff_name_only)
    execute_process(
        COMMAND bash "-c" "git diff --name-only | tr '\n' ','"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_DIFF_NAME_ONLY
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

macro(extract_version_from_changelog)
    execute_process(
        COMMAND bash "-c" "grep -m1 -E '^# v[0-9]+\.[0-9]+\.[0-9]+.*' CHANGELOG.md | cut -c 4-"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE VERSION_FULL
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(VERSION_FULL STREQUAL "")
        set(VERSION_FULL "0.0.1_invalid_version")
    endif()
    string(
        REGEX MATCH
              "^([0-9]+)\\.([0-9]+)\\.([0-9]+)(.*)"
              VERSION_MATCHED
              ${VERSION_FULL})
    set(VERSION_MAJOR ${CMAKE_MATCH_1})
    set(VERSION_MINOR ${CMAKE_MATCH_2})
    set(VERSION_PATCH ${CMAKE_MATCH_3})
    set(VERSION_SUFFIX ${CMAKE_MATCH_4})
endmacro()

macro(setup_username_hostname)
    execute_process(
        COMMAND bash "-c" "echo $(id -u -n)@$(hostname) | tr '\n' ' '"
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE USERNAME_HOSTNAME
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

macro(print_include_directories)
    get_property(
        dirs
        DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        PROPERTY INCLUDE_DIRECTORIES)
    message(STATUS "all include directories:")
    foreach(dir ${dirs})
        message(STATUS "-   ${dir}")
    endforeach()
endmacro()

macro(print_all_linked_libraries target)
    get_target_property(libs ${target} LINK_LIBRARIES)
    message(STATUS "all linked libraries: (against ${target})")
    foreach(lib ${libs})
        message(STATUS "-   ${lib}")
    endforeach()
endmacro()

macro(print_all_variables)
    get_cmake_property(vars VARIABLES)
    list(SORT vars)
    message(STATUS "all variables:")
    foreach(var ${vars})
        message(STATUS "-   ${var}=${${var}}")
    endforeach()
endmacro()

macro(configure_version_h)
    string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)
    set(VERSION_H "${CMAKE_BINARY_DIR}/${PROJECT_NAME}/version.h")
    if(CUBAO_USE_DUMMY_VERSION_H STREQUAL "True")
        configure_file(${CUBAO_UTILS_CMAKE_DIRNAME}/version.h.dummy.in ${VERSION_H} @ONLY)
    else()
        setup_username_hostname()
        setup_git_branch()
        setup_git_commit_hash()
        setup_git_commit_count()
        setup_git_commit_date()
        setup_git_diff_name_only()
        extract_version_from_changelog()
        configure_file(${CUBAO_UTILS_CMAKE_DIRNAME}/version.h.in ${VERSION_H} @ONLY)
    endif()
    install(FILES ${VERSION_H} DESTINATION "include/${PROJECT_NAME}")
endmacro()

macro(configure_output_directories)
    if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
    endif()
    if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
        set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
    endif()
    if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
    endif()
endmacro()
