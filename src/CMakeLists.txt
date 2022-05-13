include_directories(commons)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(playground)

add_executable(mmseqs-playground
        ${commons_source_files}
        ${playground_source_files}
        LocalCommandDeclarations.h
        mmseqs-playground.cpp
        )

mmseqs_setup_derived_target(mmseqs-playground)
add_dependencies(mmseqs-playground local-generated)

target_link_libraries(mmseqs-playground version)

install(TARGETS mmseqs-playground DESTINATION bin)
