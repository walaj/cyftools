#ifndef JSON_READER_H
#define JSON_READER_H

#include <iostream>
#include <string>
#include <set>
#include <fstream>
#include <json/json.h>

class JsonReader {
  
public:
 JsonReader(const std::string& filename) : filename_(filename) {}
  
  std::string GetX() const { return x_; }
  std::string GetY() const { return y_; }
  const std::set<std::string>& GetMarkers() const { return markers_; }
  const std::set<std::string>& GetMetaCells() const { return cellmeta_; }  
  
  bool ReadData() {
    std::ifstream json_file(filename_);
    if (!json_file.is_open()) {
      std::cerr << "Error opening json file " << filename_ << std::endl;
      return false;
    }
    
    Json::Value root;
    json_file >> root;
    
    x_ = root["x"].asString();
    y_ = root["y"].asString();
    z_ = root["z"].asString();    
    
    Json::Value markers = root["markers"];
    for (unsigned int i = 0; i < markers.size(); ++i) {
      markers_.insert(markers[i].asString());
    }
    
    Json::Value cellmeta = root["meta_cell"];
    for (unsigned int i = 0; i < cellmeta.size(); ++i) {
      cellmeta_.insert(cellmeta[i].asString());
    }
    
    return true;
  }
  
 private:
  std::string filename_;
  std::string x_;
  std::string y_;
  std::string z_;  
  std::set<std::string> markers_;
  std::set<std::string> cellmeta_;  
};


#endif
