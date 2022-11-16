cmake_minimum_required(VERSION 3.8)
project(gpue_ramp LANGUAGES CXX CUDA)
find_package(CUDA QUIET REQUIRED)

add_executable(
    ${PROJECT_NAME}

    src/ds.cu
    include/ds.h

    src/evolution.cu
    include/evolution.h

    src/fileIO.cu
    include/fileIO.h

    src/kernels.cu
    include/kernels.h

    src/operators.cu
    include/operators.h

    src/parser.cu
    include/parser.h

    src/split_op.cu
    include/split_op.h

    src/init.cu
    include/init.h
)

target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
	 					-O3 -lcufft 
 					   	-arch=sm_52  
 					   	-gencode=arch=compute_52,code=sm_52 
						-gencode=arch=compute_60,code=sm_60 
						-gencode=arch=compute_61,code=sm_61 
 					   	-gencode=arch=compute_70,code=sm_70 
 					   	-gencode=arch=compute_75,code=sm_75 
						>)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_link_libraries(${PROJECT_NAME} ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} "-Wl,-arch=sm_52")
set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
