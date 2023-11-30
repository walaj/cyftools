#pragma once

#include <vector>

// 255 based colors
struct Color {
  int red;
  int green;  
  int blue;
  float alpha = 1;

  float redf()   const { return static_cast<float>(red)   / 255.0f; }
  float greenf() const { return static_cast<float>(green) / 255.0f; }
  float bluef()  const { return static_cast<float>(blue)  / 255.0f; }
  float alphaf()  const { return alpha; }    
};

typedef std::vector<Color> ColorMap;

extern ColorMap color_map4;
extern ColorMap color_map8;
extern ColorMap color_map12;

extern Color color_deep_pink;
extern Color color_yellow;
extern Color color_red;
extern Color color_light_red;
extern Color color_purple;
extern Color color_light_green;
extern Color color_dark_green;
extern Color color_gray;
extern Color color_gray_90;
extern Color color_cyan;
extern Color color_dark_blue;

// Pick the appropriate colormap 
const ColorMap& getColorMap(int size);
