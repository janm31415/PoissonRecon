
set(HDRS
estimate_normals.h
poisson_recon.h
test_assert.h
)
	
set(SRCS
estimate_normals.cpp
poisson_recon.cpp
test_assert.cpp
test.cpp
)

# general build definitions
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

add_definitions(-DMEMORY_LEAK_TRACKING)

add_executable(libpoisson.test ${HDRS} ${SRCS})
source_group("Header Files" FILES ${hdrs})
source_group("Source Files" FILES ${srcs})

#if (WIN32)
#add_custom_command(TARGET libpoisson.test POST_BUILD 
#   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/data" ${CMAKE_CURRENT_BINARY_DIR}/data)
#endif (WIN32)
#
#if (UNIX)
#add_custom_command(TARGET libpoisson.test POST_BUILD 
#   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/data" "$<TARGET_FILE_DIR:libpoisson.test>/data")
#endif (UNIX)



target_include_directories(libpoisson.test
    PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/../
    )
	
target_link_libraries(libpoisson.test
    PRIVATE
    libpoisson
    )	
