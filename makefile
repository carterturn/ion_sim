build:
	nvcc -o ion_sim main.cu input_loader.cpp string_util.cpp -lglfw -lGL -lcublas --gpu-architecture sm_61
