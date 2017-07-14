#include <iostream>
#include <sys/time.h>
#include <cmath>
#include <cstring>
#include <unistd.h>

#include <GLFW/glfw3.h>
#include <GL/gl.h>

#include "particle.h"
#include "input_loader.h"

using namespace std;

#define BLOCK_SIZE 32

#define EPSILON 0.000001;

__device__ __host__ bool is_zero(double value){
	return abs(value) < EPSILON;
}

__device__ __host__ double sign(double value){
	return value >= 0 ? 1.0 : -1.0;
}

__device__ __host__ double half_abs(double value){
	return value >= 0 ? value : 0.0;
}

__global__ void calculate_force_matrix(particle * particles, double * forces_x, double * forces_y,
				       simulation_config * config){
	int ind = (blockIdx.x * BLOCK_SIZE) + threadIdx.x;

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
	double F = (K * particles[i].q * particles[j].q) /
		(d_x*d_x + d_y*d_y);
	double F_y;
	double F_x;
	double hyp = sqrt(d_x*d_x + d_y*d_y);
	F_y = is_zero(d_y) ? 0.0 : F * d_y / hyp;
	F_x = is_zero(d_x) ? 0.0 : F * d_x / hyp;

	forces_x[ind] = F_x;
	forces_y[ind] = F_y;
	forces_x[j * config->number_particles + i] = -1.0 * F_x;
	forces_y[j * config->number_particles + i] = -1.0 * F_y;
}

__global__ void move_particles(particle * particles, double * forces_x, double * forces_y,
			       simulation_config * config){
	int i = (blockIdx.x * BLOCK_SIZE) + threadIdx.x;

	if(i >= config->number_particles){
		return;
	}

	for(int j = 0; j < config->number_particles; j++){
		if(j > i){
			int ind = (i + 1) * (2*config->number_particles - i) / 2 + (j - i - 1) -
				config->number_particles;
			particles[i].vx += forces_x[ind] * config->dt / particles[i].m;
			particles[i].vy += forces_y[ind] * config->dt / particles[i].m;
		}
		else if(i > j){
			int ind = (j + 1) * (2*config->number_particles - j) / 2 + (i - j - 1) -
				config->number_particles;
			particles[i].vx += -1.0 * forces_x[ind] * config->dt / particles[i].m;
			particles[i].vy += -1.0 * forces_y[ind] * config->dt / particles[i].m;
		}
	}
	
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

int main(int argc, char* argv[]){

	if(argc < 2){
		cout << "Usage: ./a.out [config file]\n";
		return -1;
	}

	simulation_config config = load_config(argv[1]);

	glfwInit();
	GLFWwindow * window = glfwCreateWindow(config.width, config.height, "Ion Simulator", NULL, NULL);

	if(!window){
		return -1;
	}

	glfwMakeContextCurrent(window);
	glOrtho(0, config.width, 0, config.height, -1.0, 1.0);

	vector<particle> cpu_particles = load_particles(config.particles_file);
	config.number_particles = cpu_particles.size();

	int move_blocks = config.number_particles / BLOCK_SIZE;
	if(config.number_particles % BLOCK_SIZE != 0){
		move_blocks++;
	}

	int particle_data_size = config.number_particles * sizeof(particle);

	particle * gpu_particles = NULL;
	cudaError_t error = cudaMalloc(&gpu_particles, particle_data_size);
	if(error != cudaSuccess){
		cout << cudaGetErrorString(error) << "\n";
	}
	error = cudaMemcpy(gpu_particles, &cpu_particles.front(), particle_data_size, cudaMemcpyHostToDevice);
	if(error != cudaSuccess){
		cout << cudaGetErrorString(error) << "\n";
	}

	simulation_config * gpu_config = NULL;
	error = cudaMalloc(&gpu_config, sizeof(simulation_config));
	if(error != cudaSuccess){
		cout << cudaGetErrorString(error) << "\n";
	}
	error = cudaMemcpy(gpu_config, &config, sizeof(simulation_config), cudaMemcpyHostToDevice);
	if(error != cudaSuccess){
		cout << cudaGetErrorString(error) << "\n";
	}

	double * forces_x = NULL;
	double * forces_y = NULL;
	int force_matrix_size = config.number_particles * config.number_particles;
	int force_matrix_memory = force_matrix_size * sizeof(double);
	error = cudaMalloc(&forces_x, force_matrix_memory);
	if(error != cudaSuccess){
		cout << cudaGetErrorString(error) << "\n";
	}
	error = cudaMalloc(&forces_y, force_matrix_memory);
	if(error != cudaSuccess){
		cout << cudaGetErrorString(error) << "\n";
	}
	int force_blocks = force_matrix_size / BLOCK_SIZE;
	if(force_matrix_size % BLOCK_SIZE != 0){
		force_blocks++;
	}

	pair_index * force_indicies = NULL;
	error = cudaMalloc(&force_indicies, force_matrix_size * sizeof(pair_index));
	if(error != cudaSuccess){
		cout << cudaGetErrorString(error) << "\n";
	}
	pair_index * cpu_force_indicies = new pair_index[force_matrix_size];

	for(int i = 0; i < force_matrix_size; i++){
		cpu_force_indicies[i].i = 0;
		int step = config.number_particles - 1;
		for(int j = 0; j < force_matrix_size; step--){
			j += step;
			if(i < j){
				cpu_force_indicies[i].j = i - j + step + cpu_force_indicies[i].i + 1;
				break;
			}
			cpu_force_indicies[i].i++;
		}
	}
	
	error = cudaMemcpy(force_indicies, cpu_force_indicies, force_matrix_size * sizeof(pair_index),
			   cudaMemcpyHostToDevice);
	if(error != cudaSuccess){
		cout << cudaGetErrorString(error) << "\n";
	}

	int tick_count = 0;
	for(double t = 0.0; t < config.total_time; t += config.dt){
		cudaMemset(forces_x, 0, force_matrix_memory);
		cudaMemset(forces_y, 0, force_matrix_memory);
		
		calculate_force_matrix<<<force_blocks, BLOCK_SIZE>>>(gpu_particles, forces_x, forces_y, gpu_config,
			force_indicies);
		
		move_particles<<<move_blocks, BLOCK_SIZE>>>(gpu_particles, forces_x, forces_y, gpu_config);

		error = cudaMemcpy(&cpu_particles.front(), gpu_particles,
				   particle_data_size, cudaMemcpyDeviceToHost);
		if(error != cudaSuccess){
			cout << cudaGetErrorString(error) << "\n";
		}

		if(tick_count % config.ticks_per_display == 0){
		
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
	}

	cudaFree(forces_x);
	cudaFree(forces_y);
	cudaFree(gpu_particles);
	cudaFree(gpu_config);

	return 0;
}
