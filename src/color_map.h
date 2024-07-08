#pragma once

#include <vector>
#include <string>
<<<<<<< HEAD
=======
#include <iostream>
#include "cysift.h"
>>>>>>> 723463a3f1330c7129d221ed0a9589a94f517db7

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

  friend std::ostream& operator<<(std::ostream& os, const Color& obj);
  
};

struct CellColor {

  cy_uint pflag;
  cy_uint cflag;

  Color c;

  std::string label;

  friend std::ostream& operator<<(std::ostream& os, const CellColor& obj);
  
};



//using ColorLabelPair = std::pair<Color, std::string>;
typedef std::vector<CellColor> ColorLabelVec;
typedef std::vector<Color> ColorMap;

Color select_color(cy_uint cflag, cy_uint pflag, const ColorLabelVec& palette);

extern ColorMap color_map4;
extern ColorMap color_map8;
extern ColorMap color_map12;

extern Color colorbrewer_3red_dark;
extern Color colorbrewer_3red_medium;
extern Color colorbrewer_3red_light;
extern Color colorbrewer_3green_dark;
extern Color colorbrewer_3green_medium;
extern Color colorbrewer_3green_light;
extern Color colorbrewer_3blue_dark;
extern Color colorbrewer_3blue_medium;
extern Color colorbrewer_3blue_light;

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

ColorLabelVec ColorLabelVecForModule(const std::string& module);
