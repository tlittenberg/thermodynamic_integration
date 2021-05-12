include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

add_executable(thermodynamic_integration thermodynamic_integration.c)
target_link_libraries(thermodynamic_integration m)
target_link_libraries(thermodynamic_integration ${GSL_LIBRARIES})

install(TARGETS thermodynamic_integration DESTINATION bin)
