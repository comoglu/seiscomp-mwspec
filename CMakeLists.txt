SET(PLUGIN_TARGET mwspec)

SET(
	PLUGIN_SOURCES
		brune.cpp
		config.cpp
		amplitude.cpp
		combiner.cpp
		magnitude.cpp
		plugin.cpp
)

SET(
	PLUGIN_HEADERS
		brune.h
		mwspec.h
		version.h
)

SC_ADD_PLUGIN_LIBRARY(PLUGIN ${PLUGIN_TARGET} "")
SC_LINK_LIBRARIES_INTERNAL(${PLUGIN_TARGET} client)

FILE(GLOB descs "${CMAKE_CURRENT_SOURCE_DIR}/descriptions/*.xml")
INSTALL(FILES ${descs} DESTINATION ${SC3_PACKAGE_APP_DESC_DIR})

IF(SC_GLOBAL_UNITTESTS)
	SUBDIRS(test)
ENDIF(SC_GLOBAL_UNITTESTS)
