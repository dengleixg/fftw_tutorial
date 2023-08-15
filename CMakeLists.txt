cmake_minimum_required(VERSION 3.12)
project(fftw_exp)

include_directories(${PROJECT_SOURCE_DIR}/src_cpp
D:/Tools/fftw-3.3.10/install/include)

link_directories(D:/Tools/fftw-3.3.10/install/lib)

add_executable(test_1d_c2c ${PROJECT_SOURCE_DIR}/src_cpp/test_1d_c2c.cpp)

target_link_libraries(test_1d_c2c fftw3)

#set_property(TARGET test_1d_c2c PROPERTY C_STANDARD 99)

add_executable(test_1d_c2r ${PROJECT_SOURCE_DIR}/src_cpp/test_1d_c2r.cpp)

target_link_libraries(test_1d_c2r fftw3)

add_executable(test_1d_r2c ${PROJECT_SOURCE_DIR}/src_cpp/test_1d_r2c.cpp)

target_link_libraries(test_1d_r2c fftw3)

add_executable(test_2d_c2c ${PROJECT_SOURCE_DIR}/src_cpp/test_2d_c2c.cpp)

target_link_libraries(test_2d_c2c fftw3)