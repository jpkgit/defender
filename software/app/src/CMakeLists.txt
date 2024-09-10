set(INSTALL_DEFAULT_BINDIR "bin" CACHE STRING "Appended to CMAKE_INSTALL_PREFIX")

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
get_filename_component(FFTW_LIBRARY_DIRS ${FFTW_LIBRARIES} DIRECTORY)
link_directories(${FFTW_LIBRARY_DIRS})

SET(TOOLS	
	baseline_alerts
	fcc_table
)

if(MSVC)
	add_library(libgetopt_static STATIC
	    ../getopt/getopt.c
	)
	LIST(APPEND TOOLS_LINK_LIBS ${FFTW_LIBRARIES})
else()
    LIST(APPEND TOOLS_LINK_LIBS m fftw3f)
endif()

if(NOT libhackrf_SOURCE_DIR)
	include_directories(${LIBHACKRF_INCLUDE_DIR})
	LIST(APPEND TOOLS_LINK_LIBS ${LIBHACKRF_LIBRARIES})
else()
	LIST(APPEND TOOLS_LINK_LIBS hackrf)
endif()

if(MSVC)
	LIST(APPEND TOOLS_LINK_LIBS libgetopt_static)
endif()

add_executable(baseline_alerts baseline_alerts.c fcc_table.c)
target_link_libraries(baseline_alerts ${TOOLS_LINK_LIBS})
install(TARGETS baseline_alerts RUNTIME DESTINATION ${INSTALL_DEFAULT_BINDIR})

file(COPY fcc_table.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


if( ${WIN32} )
	install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIGURATION>/"
        	DESTINATION ${INSTALL_DEFAULT_BINDIR}
			FILES_MATCHING
			PATTERN "fftw*.dll")
endif( ${WIN32} )
