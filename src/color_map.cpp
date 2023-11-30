#include "color_map.h"
#include <stdexcept>

Color color_yellow = {255,255,0}; 
Color color_deep_pink = {255,20,147};
Color color_red = {255,0,0};
Color color_light_red = {255,102,102};
Color color_cyan = {0,255,255};
Color color_light_green = {144, 238, 144};
Color color_dark_green = {0,100,0};
Color color_gray = {100,100,100};
Color color_gray_90 = {100,100,100,0.1};
Color color_purple = {128,0,128};
Color color_dark_blue = {56,108,176};

  //https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=4
  ColorMap color_map4 = {
    {27,158,119},
    {217,95,2},
    {117,112,179},
    {231,41,138}
  };

  //https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=8
  ColorMap color_map8 = {
    {228,26,28},
    {55,126,184},
    {77,175,74},
    {152,78,163},
    {255,127,0},
    {255,255,51},
    {166,86,40},
    {247,129,191}
  };
  
  // https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
  ColorMap color_map12 = {
    {166,206,227}, //light blue
    {31,120,180},  //dark blue
    {178,223,138}, //light green
    {51,160,44},   //dark green
    {251,154,153}, //pink
    {227,26,28},   //red
    {253,191,111}, //peach
    {255,127,0},   //orange
    {202,178,214}, // light purple
    {106,61,154},  // purple
    {255,255,153}, // yellow
    {177,89,40}    // brown
  };

// get the appropriate color map
const ColorMap& getColorMap(int size) {
    if (size <= 0 || size > 12) {
        throw std::invalid_argument("Size must be between 1 and 12.");
    } else if (size <= 4) {
        return color_map4;
    } else if (size <= 8) {
        return color_map8;
    } else {
        return color_map12;
    }
}
