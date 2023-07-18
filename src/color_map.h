#pragma once

#include <vector>

// 255 based colors
struct Color {
  int red;
  int green;  
  int blue;

  float redf()   const { return static_cast<float>(red)   / 255.0f; }
  float greenf() const { return static_cast<float>(green) / 255.0f; }
  float bluef()  const { return static_cast<float>(blue)  / 255.0f; }  
};

typedef std::vector<Color> ColorMap;

extern ColorMap color_map4;
extern ColorMap color_map8;
extern ColorMap color_map12;

// Pick the appropriate colormap 
const ColorMap& getColorMap(int size);
