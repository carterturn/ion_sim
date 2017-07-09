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
#define WIDTH 512
#define HEIGHT 512
#define WALL_ELASTICITY 1.0
#define METERS_PER_SQUARE 0.01

#define TOTAL_TIME 50.0
#define DT 0.001
#define TICKS_PER_DISPLAY 50

#define EPSILON 0.000001;

int blocks;

__managed__ int number_particles;

__device__ __host__ bool is_zero(double value){
	return abs(value) < EPSILON;
}

__device__ __host__ double sign(double value){
	return value >= 0 ? 1.0 : -1.0;
}

__global__ void move_particles(particle * particles){
	int i = (blockIdx.x * BLOCK_SIZE) + threadIdx.x;

	if(i > number_particles) return;
	
	particles[i].x += particles[i].vx * DT / METERS_PER_SQUARE;
	if(particles[i].x > WIDTH){
		particles[i].x = WIDTH;
		particles[i].vx = -particles[i].vx * WALL_ELASTICITY;
	}
	if(particles[i].x < 0){
		particles[i].x = 0;
		particles[i].vx = -particles[i].vx * WALL_ELASTICITY;
	}
	particles[i].y += particles[i].vy * DT / METERS_PER_SQUARE;
	if(particles[i].y > HEIGHT){
		particles[i].y = HEIGHT;
		particles[i].vy = -particles[i].vy * WALL_ELASTICITY;
	}
	if(particles[i].y < 0){
		particles[i].y = 0;
		particles[i].vy = -particles[i].vy * WALL_ELASTICITY;
	}

	for(int j = 0; j < number_particles; j++){
		if(i != j){
			double d_x = (particles[i].x - particles[j].x) * METERS_PER_SQUARE;
			double d_y = (particles[i].y - particles[j].y) * METERS_PER_SQUARE;
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
			particles[i].vx += A_x * DT;
			particles[i].vy += A_y * DT;
		}
	}

}

int main(int argc, char* argv[]){

	if(argc < 2){
		cout << "Usage: ./a.out [particle config file]\n";
		return -1;
	}

	glfwInit();
	GLFWwindow * window = glfwCreateWindow((int) WIDTH, (int) HEIGHT, "Ion Simulator", NULL, NULL);

	if(!window){
		return -1;
	}

	glfwMakeContextCurrent(window);
	glOrtho(0, WIDTH, 0, HEIGHT, -1.0, 1.0);

	vector<particle> cpu_particles = load_particles(argv[1]);
	number_particles = cpu_particles.size();

	cout << "A\n";

	blocks = number_particles / BLOCK_SIZE;
	if(number_particles % BLOCK_SIZE != 0){
		blocks++;
	}

	int particle_data_size = number_particles * sizeof(particle);

	cout << "A\n";
	
	particle * gpu_particles = NULL;
	cudaError_t error = cudaMalloc(&gpu_particles, particle_data_size);
	if(error != cudaSuccess){
		cout << cudaGetErrorString(error) << "\n";
	}
	error = cudaMemcpy(gpu_particles, &cpu_particles.front(), particle_data_size, cudaMemcpyHostToDevice);
	if(error != cudaSuccess){
		cout << cudaGetErrorString(error) << "\n";
	}

	cout << "A\n";

	int tick_count = 0;
	for(double t = 0.0; t < TOTAL_TIME; t += DT){
		timespec start;
		clock_gettime(CLOCK_REALTIME, &start);
		move_particles<<<blocks, BLOCK_SIZE>>>(gpu_particles);

		error = cudaMemcpy(&cpu_particles.front(), gpu_particles, particle_data_size, cudaMemcpyDeviceToHost);
		if(error != cudaSuccess){
			cout << cudaGetErrorString(error) << "\n";
		}

		if(tick_count % TICKS_PER_DISPLAY == 0){
			glClear(GL_COLOR_BUFFER_BIT);
			glBegin(GL_POINTS);
			glColor3f(0.0f, 1.0f, 0.0f);
			for(int i = 0; i < number_particles; i++){
				glVertex2i((int) cpu_particles[i].x, (int) cpu_particles[i].y);
			}
			glEnd();
			glfwSwapBuffers(window);
			timespec end;
			clock_gettime(CLOCK_REALTIME, &end);
			cout << t << "\n";
		}
		tick_count++;
	}

	cudaFree(gpu_particles);
		
	return 0;
}
