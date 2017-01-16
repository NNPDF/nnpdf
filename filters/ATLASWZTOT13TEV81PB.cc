
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <exception>
#include <utility>

#include "ATLASWZTOT13TEV81PB.h"

using namespace std;

//Convert from nb to fb
const auto unit_conversion_factor = 1e6;

const char comment = '#';
const char delimiter = '\t';


void ATLASWZTOT13TEV81PBFilter::ReadData()
{
  const auto data_file = dataPath() + "rawdata/" + fSetName + "/data.txt";
  const auto sys_file = dataPath() + "rawdata/" + fSetName + "/sys_description.txt";

  //Parse Data file
  //We need these two to get an insertion ordered dict...
  auto result = unordered_map<string, vector<double>>{};
  auto keys = vector<string>{};

  auto index_labels = vector<string>{};
  { // Could split this in a function...
  ifstream infile(data_file);
  if(!infile){
    throw runtime_error("Cannot open file: " +  data_file);
  }
  auto line = string{};
  while(getline(infile, line)){
     if (line[line.find_first_not_of(" ")] == comment){
     continue;
     }
     stringstream iss(line);
     auto token = string{};
     while(getline(iss, token, delimiter)){
       if (token.size()){
         keys.push_back(token);
         result[token] = {};
       }
     }
     break;
  }

  while(getline(infile, line)){
     if (line[line.find_first_not_of(" ")] == comment){
     continue;
     }
     stringstream iss(line);
     auto token = string{};
     getline(iss, token, delimiter);
     index_labels.push_back(token);
     for (auto & key : keys){
       auto & values = result.at(key);
       getline(iss, token, delimiter);
       values.push_back(stod(token)*unit_conversion_factor);

     }
  }
  }
  // Parse sysfile
  auto sys_desc = unordered_map<string, pair<string, string>>{};
  {
  ifstream infile(sys_file);
  if(!infile){
    throw runtime_error("Cannot open file: " +  sys_file);
  }

  auto line = string{};
  while(getline(infile, line)){
    stringstream iss(line);
	auto token = string{};
	getline(iss, token, delimiter);
	auto key = token; 
	getline(iss, token, delimiter);
	auto type = token;
	getline(iss, token, delimiter);
	auto name = token;
	sys_desc[key] = {type,name};

  }

  for (auto & key : keys){
    auto & values = result.at(key);
	if(values.size()!=size_t(fNData)){
		throw runtime_error("Size of " + key + " does not coincide with fNdata.");
    }
  }
  }

  auto & k1 = result.at("k1");
  auto & k2 = result.at("k2");
  auto & k3 = result.at("k3");
  copy(k1.begin(), k1.end(), fKin1);
  copy(k2.begin(), k2.end(), fKin2);
  copy(k3.begin(), k3.end(), fKin3);
  result.erase("k1");
  result.erase("k2");
  result.erase("k3");

  auto & stat = result.at("stat");
  copy(stat.begin(), stat.end(), fStat);
  result.erase("stat");

  auto & value = result.at("value");
  copy(value.begin(), value.end(), fData);
  result.erase("value");

  size_t sys_id = 0;
  for(auto & k:keys){
	  if(!result.count(k)){
		  continue;
	  }
	  auto & desc = sys_desc.at(k);
	  auto name = desc.first;
	  auto type = desc.second;
	  for (size_t data_id=0; data_id<size_t(fNData); data_id++){
		  auto & sys = fSys[data_id][sys_id];
		  if (type == "ADD"){
			  sys.type = ADD;
		  }else if (type == "MULT"){
			  sys.type = MULT;
		  }else {
			  throw runtime_error("Unrecognized systematic type: " + type);
		  }
		  auto & sysvalue = result.at(k)[data_id];
		  auto & datavalue = fData[data_id];
		  sys.add = sysvalue;
		  sys.mult = sysvalue/datavalue*100;
      sys.name = name;
	  }
	  result.erase(k);
	  sys_id++;
  }

  if(sys_id != size_t(fNSys)){
	  throw runtime_error("Wrong number of systematics");
  }

 if (result.size()){
	 throw runtime_error("Unused values in data file");
 } 

}


