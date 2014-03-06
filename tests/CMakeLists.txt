ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
set(BOOST_COMPONENTS unit_test_framework timer system)
SEARCH_FOR_BOOST()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

set(HEADERS TwoDofArm.h)

macro(addUnitTest name)
  if(${UNIT_TESTS})
    add_executable(${name} ${name}.cpp ${HEADERS})
    target_link_libraries(${name} ${Boost_LIBRARIES} SpringEstimator)
    PKG_CONFIG_USE_DEPENDENCY(${name} SpaceVecAlg)
    PKG_CONFIG_USE_DEPENDENCY(${name} RBDyn)
    add_test(${name}Unit ${name})
  endif()
endmacro(addUnitTest)

addUnitTest("SpringEstimatorTest")
