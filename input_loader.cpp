#include "input_loader.h"
#include "string_util.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

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
