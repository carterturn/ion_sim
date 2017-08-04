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

#include "input_loader.h"
#include "string_util.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

simulation_config load_config(string filename){
		fstream config_file(filename.c_str());

		simulation_config config;

		while(!config_file.eof()){
				string config_line;
				getline(config_file, config_line);

				vector<string> param_value_pair = split_string(config_line, '=');
				if(param_value_pair.size() != 2){
						cerr << "Malformed line\n";
						continue;
				}
				
				string param = trim_whitespace(param_value_pair[0]);
				string value = trim_whitespace(param_value_pair[1]);

				if(param == "width"){
						config.width = atoi(value.c_str());
				}
				else if(param == "height"){
						config.height = atoi(value.c_str());
				}
				else if(param == "wall_elasticity"){
						config.wall_elasticity = strtod(value.c_str(), NULL);
				}
				else if(param == "meters_per_square"){
						config.meters_per_square = strtod(value.c_str(), NULL);
				}
				else if(param == "total_time"){
						config.total_time = strtod(value.c_str(), NULL);
				}
				else if(param == "dt"){
						config.dt = strtod(value.c_str(), NULL);
				}
				else if(param == "ticks_per_display"){
						config.ticks_per_display = atoi(value.c_str());
				}
				else if(param == "particles_file"){
						config.particles_file = value;
				}
		}

		return config;
}

/**
   CSV format: x, y, vx, vy, q, m
 */

vector<particle> load_particles(string filename){
		fstream particle_file(filename.c_str());

		vector<particle> particles;

		int iter = 0;
		while(!particle_file.eof()){
				string particle_config;
				getline(particle_file, particle_config);
				iter++;

				vector<string> columns = split_string(particle_config, ',');
				
				if(columns.size() < 6){
						cerr << "Malformed line " << iter << "\n";
						continue;
				}

				for(int i = 0; i < columns.size(); i++){
						columns[i] = trim_whitespace(columns[i]);
				}

				particle new_particle;
				new_particle.x = strtod(columns[0].c_str(), NULL);
				new_particle.y = strtod(columns[1].c_str(), NULL);
				new_particle.vx = strtod(columns[2].c_str(), NULL);
				new_particle.vy = strtod(columns[3].c_str(), NULL);
				new_particle.q = strtod(columns[4].c_str(), NULL);
				new_particle.m = strtod(columns[5].c_str(), NULL);

				particles.push_back(new_particle);
		}

		return particles;
}
