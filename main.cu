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

int blocks;

__device__ __host__ bool is_zero(double value){
	return abs(value) < EPSILON;
}

__device__ __host__ double sign(double value){
	return value >= 0 ? 1.0 : -1.0;
}

__global__ void move_particles(particle * particles, simulation_config * config){
	int i = (blockIdx.x * BLOCK_SIZE) + threadIdx.x;

	if(i > config->number_particles) return;
	
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

	for(int j = 0; j < config->number_particles; j++){
		if(i != j){
			double d_x = (particles[i].x - particles[j].x) * config->meters_per_square;
			double d_y = (particles[i].y - particles[j].y) * config->meters_per_square;
			double A = (K * particles[i].q * particles[j].q) /
				(particles[i].m * (d_x*d_x + d_y*d_y));
			double A_y;
			double A_x;
			if(is_zero(d_y)){
				A_x = A * sign(d_x);
				A_y = 0.0;
			}
			else if(is_zero(d_x)){
				A_x = A * sign(d_y);
				A_y = 0.0;
			}
			else{
				double hyp = sqrt(d_x*d_x + d_y*d_y);
				A_y = A * d_y / hyp;
				A_x = A * d_x / hyp;
			}
			particles[i].vx += A_x * config->dt;
			particles[i].vy += A_y * config->dt;
		}
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

	blocks = config.number_particles / BLOCK_SIZE;
	if(config.number_particles % BLOCK_SIZE != 0){
		blocks++;
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

	int tick_count = 0;
	for(double t = 0.0; t < config.total_time; t += config.dt){
		move_particles<<<blocks, BLOCK_SIZE>>>(gpu_particles, gpu_config);

		error = cudaMemcpy(&cpu_particles.front(), gpu_particles, particle_data_size, cudaMemcpyDeviceToHost);
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

	cudaFree(gpu_particles);
	cudaFree(gpu_config);
		
	return 0;
}
