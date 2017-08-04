/*
  Copyright 2017 Carter Turnbaugh

  This file is part of Ion Sim.

  Ion Sim is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Ion Sim is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Ion Sim.  If not, see <http://www.gnu.org/licenses/>.
*/

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
