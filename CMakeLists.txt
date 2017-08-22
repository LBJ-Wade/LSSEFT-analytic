CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(LSSEFT_analytic)

SET(CMAKE_CXX_STANDARD 14)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
INCLUDE(BuildGiNaC)
INCLUDE(CompilerFlags)

# detect toolchain and set compiler flags appropriately
SET_COMPILER_FLAGS()

# RESOLVE DEPENDENCIES

# find required Boost libraries
FIND_PACKAGE(Boost 1.56 REQUIRED COMPONENTS timer)

# find GiNaC libraries
IF(NOT FORCE_BUILD_GINAC)
  FIND_PACKAGE(GiNaC)
ENDIF()

# handle non-standard external dependencies
ADD_CUSTOM_TARGET(DEPS ALL)
SET(EXTERNAL_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/thirdparty/dist")

IF(NOT GINAC_FOUND)
  # GiNaC was not found, so attempt to build it from source
  BUILD_CLN(${CMAKE_BINARY_DIR}/thirdparty/CLN ${EXTERNAL_INSTALL_PREFIX})
  BUILD_GINAC(${CMAKE_BINARY_DIR}/thirdparty/GiNaC ${CLN_INCLUDE_DIRS} ${CLN_LIBRARY_DIRS} ${CLN_LIBRARIES} ${EXTERNAL_INSTALL_PREFIX})
ENDIF()


# add LSSEFT_analytic executable

ADD_EXECUTABLE(LSSEFT_analytic
  main.cpp
  lib/kernel.cpp
  lib/vector.cpp
  services/symbol_factory.cpp
  shared/exceptions.cpp
  )

ADD_DEPENDENCIES(LSSEFT_analytic DEPS)

TARGET_LINK_LIBRARIES(LSSEFT_analytic ${GINAC_LIBRARIES} ${Boost_LIBRARIES})

TARGET_INCLUDE_DIRECTORIES(LSSEFT_analytic PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${GINAC_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  )


# add dummy target for CLion

SET(LIBRARY_FILES
  lib/kernel.cpp
  lib/kernel.h
  lib/vector.cpp
  lib/vector.h
  )

SET(LOCALIZATIONS_FILES
  localizations/en/messages_en.h
  localizations/messages.h
  )

SET(SERVICES_FILES
  services/symbol_factory.cpp
  services/symbol_factory.h
  )

SET(SHARED_FILES
  shared/exceptions.cpp
  shared/exceptions.h
  shared/defaults.h)

SET(TOP_LEVEL_FILES
  main.cpp
  )

ADD_EXECUTABLE(clion_dummy_target
  ${LIBRARY_FILES}
  ${LOCALIZATIONS_FILES}
  ${SERVICES_FILES}
  ${SHARED_FILES}
  ${TOP_LEVEL_FILES}
  )
