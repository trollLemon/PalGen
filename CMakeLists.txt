cmake_minimum_required(VERSION 3.31)
project(PalGen C)

set(CMAKE_C_STANDARD 23)




set(LIB lib)
set(INCLUDES include)

include_directories(${INCLUDES})
include_directories(${LIB})

set (
        SOURCE_FILES
        cmd/main.c
        ColorUtils/printColor.c
        Algorithms/K-Means/kmeans.c
)


# Tests
enable_testing()

add_library(ColorUtils STATIC ColorUtils/printColor.c)
add_executable(test_hex tests/test_printColor.c)
target_compile_options(test_hex PRIVATE -fsanitize=address -g)
target_link_libraries(test_hex PRIVATE ColorUtils m -fsanitize=address)

add_library(KMeans STATIC Algorithms/K-Means/kmeans.c)
add_executable(test_kmeans tests/test_kmeans.c)
target_compile_options(test_kmeans PRIVATE -fsanitize=address -g)
target_link_libraries(test_kmeans PRIVATE KMeans m -fsanitize=address)



add_test(NAME test_printColor COMMAND test_hex)
add_test(NAME test_KMeans COMMAND test_kmeans)
# build
add_executable(PalGen ${SOURCE_FILES})
target_compile_options(PalGen PRIVATE -fsanitize=address -g)
target_link_libraries(PalGen PRIVATE m -fsanitize=address )
