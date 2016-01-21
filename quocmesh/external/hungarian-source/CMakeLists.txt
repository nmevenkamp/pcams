# You can download hungarian from:
# http://www.frc.ri.cmu.edu/~lantao/codes/hungarian-v2.0.zip
#
# Lantao Liu, Dylan Shell. "Assessing Optimal Assignment under Uncertainty: An Interval-based Algorithm".
# International Journal of Robotics Research (IJRR). vol. 30, no. 7, pp 936-953. Jun 2011.

# Silence all GCC warnings triggered by wavelet1d.
IF ( CMAKE_COMPILER_IS_GNUCXX )
  ADD_DEFINITIONS ( -w )
  REMOVE_DEFINITIONS ( -pedantic )
ENDIF ( CMAKE_COMPILER_IS_GNUCXX )

SET ( HDRS
  Hungarian.h
  Assignment.h
  BipartiteGraph.h
  HungarianMatrix.h
)

SET ( SRCS
  Hungarian.cpp
  Assignment.cpp
  BipartiteGraph.cpp
)

SET ( EXTERNAL_HUNGARIAN_INCLUDE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL ""
)

ADD_LIBRARY ( hungarian STATIC ${SRCS} ${HDRS} )

SET ( EXTERNAL_HUNGARIAN_LIBRARIES hungarian CACHE INTERNAL "" )