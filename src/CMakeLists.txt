set(SOURCES SpringEstimator.cpp)
set(HEADERS SpringEstimator.h)

add_library(SpringEstimator SHARED ${SOURCES} ${HEADERS})
PKG_CONFIG_USE_DEPENDENCY(SpringEstimator SpaceVecAlg)
PKG_CONFIG_USE_DEPENDENCY(SpringEstimator RBDyn)

set(INSTALL_PATH include/SpringEstimator)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS SpringEstimator DESTINATION "lib")
