cmake_minimum_required(VERSION 3.12)


project(dustdrift VERSION 0.1
	DESCRIPTION "1D dust drift"
        LANGUAGES CXX)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})


add_library(dustdrift SHARED 
	src/define.cpp include/define.h
	src/funcs.cpp include/funcs.h
	src/initial.cpp include/initial.h
	src/iterate.cpp include/iterate.h
	src/numerics.cpp include/numerics.h
	src/output.cpp include/output.h 
	src/parser.cpp include/parser.h )

target_include_directories(dustdrift PUBLIC include)


add_executable(1D_drift apps/main.cpp)
target_link_libraries(1D_drift PUBLIC dustdrift)
