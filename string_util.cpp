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
