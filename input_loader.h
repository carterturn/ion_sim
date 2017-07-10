#pragma once

#include <string>
#include <vector>

#include "particle.h"

struct simulation_config {
		int width;
		int height;
		double wall_elasticity;
		double meters_per_square;
		double total_time;
		double dt;
		int ticks_per_display;
		std::string particles_file;
		int number_particles;
};

simulation_config load_config(std::string filename);

std::vector<particle> load_particles(std::string filename);
