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

#include "string_util.h"

#include <sstream>
#include <vector>

using namespace std;

vector<string> split_string(string data, char splitter){
	vector<string> tokens = vector<string>();
	string tmp = "";
	for(int i = 0; i < data.length(); i++){
		if(data[i] == splitter){
			tokens.push_back(tmp);
			tmp = "";
		}
		else{
			tmp += data[i];
		}
	}
	tokens.push_back(tmp);
	return tokens;
}


string its(int in){
	stringstream a;
	a << in;
	return a.str();
}

string trim_whitespace(string data){
		size_t data_start = data.find_first_not_of(" \t");
		if(data_start == string::npos){
				return "";
		}

		return data.substr(data_start, data.find_last_not_of(" \t") - data_start + 1);
}
