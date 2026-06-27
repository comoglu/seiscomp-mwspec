# Unit tests for the framework-independent Brune/Seisan core. These do not
# depend on any SeisComP libraries.

ADD_EXECUTABLE(test_mwspec_brune test_brune.cpp ../brune.cpp)

ADD_TEST(
	NAME test_mwspec_brune
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMAND test_mwspec_brune
)
