
add_subdirectory(geometryMX)
add_subdirectory(geometryXNeutrals)
add_subdirectory(fluidsolvers)
add_subdirectory(plasmaneutralscoupling)
add_subdirectory(reactionrates)
add_subdirectory(plasmasources)

add_library(spline_definitions INTERFACE)
target_include_directories(spline_definitions
	INTERFACE
	${PROJECT_SOURCE_DIR}/gyselalibxx/simulations/geometryXVx
)
add_library(voicexx::spline_definitions ALIAS spline_definitions)
