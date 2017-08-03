#define GRAPHICS

#include <iostream>
#include <sys/time.h>
#include <cmath>
#include <cstring>
#include <unistd.h>

#ifdef GRAPHICS
#include <GLFW/glfw3.h>
#include <GL/gl.h>
#endif

#include <cublas_v2.h>
#include <cuda_profiler_api.h>

#include "particle.h"
#include "input_loader.h"

using namespace std;

#define EPSILON 0.000001

__global__ void calculate_force_matrix(particle * particles, double * forces_x, double * forces_y,
				       simulation_config * config){
	int ind = (blockIdx.x * blockDim.x) + threadIdx.x;

	if(ind >= config->number_particles * (config->number_particles + 1) / 2){
		return;
	}

	int i = ind / config->number_particles;
	int j = ind - i * config->number_particles;

	if(j < i){
		ind = config->number_particles * config->number_particles - ind + config->number_particles - 1;
		i = config->number_particles - i;
		j = config->number_particles - j - 1;
	}

	double d_x = (particles[i].x - particles[j].x) * config->meters_per_square;
	double d_y = (particles[i].y - particles[j].y) * config->meters_per_square;
	double F = (K * particles[i].q * particles[j].q) / (d_x*d_x + d_y*d_y);
	double hyp = sqrt(d_x*d_x + d_y*d_y);
	double F_y = (abs(d_y) > EPSILON) * (F * d_y / hyp);
	double F_x = (abs(d_x) > EPSILON) * (F * d_x / hyp);

	forces_x[ind] = F_x;
	forces_y[ind] = F_y;
	forces_x[j * config->number_particles + i] = -1.0 * F_x;
	forces_y[j * config->number_particles + i] = -1.0 * F_y;
}

__global__ void move_particles(particle * particles, double * forces_x, double * forces_y,
			       simulation_config * config){
	int i = (blockIdx.x * blockDim.x) + threadIdx.x;

	if(i >= config->number_particles){
		return;
	}

	particles[i].vx += forces_x[i] * config->dt / particles[i].m;
	particles[i].vy += forces_y[i] * config->dt / particles[i].m;

	particles[i].x += particles[i].vx * config->dt / config->meters_per_square;
	if(particles[i].x > config->width){
		particles[i].x = config->width;
		particles[i].vx = -particles[i].vx * config->wall_elasticity;
	}
	if(particles[i].x < 0){
		particles[i].x = 0;
		particles[i].vx = -particles[i].vx * config->wall_elasticity;
	}
	particles[i].y += particles[i].vy * config->dt / config->meters_per_square;
	if(particles[i].y > config->height){
		particles[i].y = config->height;
		particles[i].vy = -particles[i].vy * config->wall_elasticity;
	}
	if(particles[i].y < 0){
		particles[i].y = 0;
		particles[i].vy = -particles[i].vy * config->wall_elasticity;
	}
}

void check_error(cudaError_t error, int line){
	if(error != cudaSuccess){
		cerr << "Error on line " << line << "\n\t" << cudaGetErrorString(error) << "\n";
	}
}

string cublasGetErrorString(cublasStatus_t status){
	switch(status){
        case CUBLAS_STATUS_SUCCESS:
		return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED:
		return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED:
		return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE:
		return "CUBLAS_STATUS_INVALID_VALUE"; 
        case CUBLAS_STATUS_ARCH_MISMATCH:
		return "CUBLAS_STATUS_ARCH_MISMATCH"; 
        case CUBLAS_STATUS_MAPPING_ERROR:
		return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED:
		return "CUBLAS_STATUS_EXECUTION_FAILED"; 
        case CUBLAS_STATUS_INTERNAL_ERROR:
		return "CUBLAS_STATUS_INTERNAL_ERROR";
	default:
		return "unknown error";
	}
}

void check_error(cublasStatus_t status, int line){
	if(status != CUBLAS_STATUS_SUCCESS){
		cerr << "Error on line " << line << "\n\t" << cublasGetErrorString(status) << "\n";
	}
}

int main(int argc, char* argv[]){

	if(argc < 2){
		cout << "Usage: ./a.out [config file]\n";
		return -1;
	}

	simulation_config config = load_config(argv[1]);

#ifdef GRAPHICS
	glfwInit();
	GLFWwindow * window = glfwCreateWindow(config.width, config.height, "Ion Simulator", NULL, NULL);

	if(!window){
		return -1;
	}

	glfwMakeContextCurrent(window);
	glOrtho(0, config.width, 0, config.height, -1.0, 1.0);
#endif

	cudaDeviceProp props;
	check_error(cudaGetDeviceProperties(&props, 0), __LINE__);

	vector<particle> cpu_particles = load_particles(config.particles_file);
	config.number_particles = cpu_particles.size();

	int move_bsize = props.warpSize;
	int move_blocks;
	while(config.number_particles / move_bsize + 1 > props.maxGridSize[1]){
		move_bsize *= 2;
	}
	if(move_bsize > 1024){
		cerr << "Error: can not simulate " << config.number_particles << " on this platform\n";
	}
	move_blocks = config.number_particles / move_bsize + 1;

	int particle_data_size = config.number_particles * sizeof(particle);

	particle * gpu_particles = NULL;
	check_error(cudaMalloc(&gpu_particles, particle_data_size), __LINE__);
	check_error(cudaMemcpy(gpu_particles, &cpu_particles.front(), particle_data_size, cudaMemcpyHostToDevice),
		    __LINE__);

	simulation_config * gpu_config = NULL;
	check_error(cudaMalloc(&gpu_config, sizeof(simulation_config)), __LINE__);
	check_error(cudaMemcpy(gpu_config, &config, sizeof(simulation_config), cudaMemcpyHostToDevice), __LINE__);

	double * forces_x = NULL;
	double * forces_y = NULL;
	int force_matrix_size = config.number_particles * config.number_particles;
	int force_matrix_memory = force_matrix_size * sizeof(double);
	check_error(cudaMalloc(&forces_x, force_matrix_memory), __LINE__);
	check_error(cudaMalloc(&forces_y, force_matrix_memory), __LINE__);

	int force_bsize;
	int force_blocks;
	{
		force_bsize = props.warpSize;
		int force_threads = config.number_particles * (config.number_particles + 1) / 2;
		while(force_threads / force_bsize + 1 > props.maxGridSize[1]){
			force_bsize *= 2;
		}
		force_blocks = force_threads / force_bsize + 1;
	}

	cout << "Using " << force_blocks << " blocks of size " << force_bsize << " to calculate force matrix\n";
	cout << "Using " << move_blocks << " blocks of size " << move_bsize << " to move particles\n";

        cublasHandle_t blas_handle;
	cublasCreate(&blas_handle);

	double * total_forces_x = NULL;
	double * total_forces_y = NULL;
	int total_force_matrix_memory = config.number_particles * sizeof(double);
	check_error(cudaMalloc(&total_forces_x, total_force_matrix_memory), __LINE__);
	check_error(cudaMalloc(&total_forces_y, total_force_matrix_memory), __LINE__);

	double * forces_ones = NULL;
	check_error(cudaMalloc(&forces_ones, total_force_matrix_memory), __LINE__);
	double * cpu_forces_ones = new double[config.number_particles];
	for(int i = 0; i < config.number_particles; i++){
		cpu_forces_ones[i] = 1.0;
	}
	check_error(cudaMemcpy(forces_ones, cpu_forces_ones, total_force_matrix_memory, cudaMemcpyHostToDevice),
		    __LINE__);
	delete[] cpu_forces_ones;

#ifdef GRAPHICS
	int tick_count = 0;
#endif
	double alpha = 1.0;
	double beta = 0.0;
	for(double t = 0.0; t < config.total_time; t += config.dt){

		// Start stage 1: calculate force matrix
		calculate_force_matrix<<<force_blocks, force_bsize>>>(gpu_particles, forces_x, forces_y, gpu_config);
		check_error(cudaGetLastError(), __LINE__);
		
		check_error(cudaThreadSynchronize(), __LINE__);
		// Synchronize as we end stage 1

		// Start stage 2: sum across rows of force matrix
		check_error(cublasDgemv(blas_handle, CUBLAS_OP_T, config.number_particles, config.number_particles,
					&alpha, forces_x, config.number_particles, forces_ones, 1, &beta,
					total_forces_x, 1), __LINE__);
		check_error(cudaGetLastError(), __LINE__);
		check_error(cublasDgemv(blas_handle, CUBLAS_OP_T, config.number_particles, config.number_particles,
					&alpha, forces_y, config.number_particles, forces_ones, 1, &beta,
					total_forces_y, 1), __LINE__);
		check_error(cudaGetLastError(), __LINE__);

		check_error(cudaThreadSynchronize(), __LINE__);
		// Synchronize as we end stage 2

		// Start stage 3: update particles from forces
		move_particles<<<move_blocks, move_bsize>>>(gpu_particles, total_forces_x, total_forces_y,
							    gpu_config);
		check_error(cudaGetLastError(), __LINE__);

		check_error(cudaThreadSynchronize(), __LINE__);
		// Synchronize as we end stage 3

#ifdef GRAPHICS
		if(tick_count % config.ticks_per_display == 0){
			check_error(cudaMemcpy(&cpu_particles.front(), gpu_particles,
					       particle_data_size, cudaMemcpyDeviceToHost), __LINE__);

			glClear(GL_COLOR_BUFFER_BIT);
			glBegin(GL_POINTS);
			glColor3f(0.0f, 1.0f, 0.0f);
			for(int i = 0; i < config.number_particles; i++){
				glVertex2i((int) cpu_particles[i].x, (int) cpu_particles[i].y);
			}
			glEnd();
			glfwSwapBuffers(window);
			cout << t << "\n";
		}
		tick_count++;
#endif
	}

	cudaFree(forces_x);
	cudaFree(forces_y);
	cudaFree(forces_ones);
	cudaFree(gpu_particles);
	cudaFree(gpu_config);

	cublasDestroy(blas_handle);

	cudaProfilerStop();

	return 0;
}
