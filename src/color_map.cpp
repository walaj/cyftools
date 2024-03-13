#include "color_map.h"
#include <stdexcept>
#include <string>
#include <iostream>
#include <cassert>

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

Color colorbrewer_3red_dark     = {222,45 ,38};
Color colorbrewer_3red_medium   = {252,146,114};
Color colorbrewer_3red_light    = {254,224,210};
Color colorbrewer_3green_dark   = {49 ,163,84}; 
Color colorbrewer_3green_medium = {161,217,155};
Color colorbrewer_3green_light  = {229,245,224};
Color colorbrewer_3blue_dark    = {49 ,130,189};
Color colorbrewer_3blue_medium  = {158,202,225};
Color colorbrewer_3blue_light   = {222,235,247};

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

// get color map by module name
ColorLabelMap ColorLabelMapForModule(const std::string& module) {

  ColorLabelMap cm;
  if (module == "tumor") {
    cm = {
      {colorbrewer_3red_light,    "Tumor (Automated)"},
      {colorbrewer_3red_medium,   "Tumor (Manual)"},
      {colorbrewer_3red_dark,     "Tumor (A+M)"},
      {colorbrewer_3blue_light,   "Stroma"},
      {color_cyan,                "Margin"},
      {color_deep_pink,           "Tcell cluster"},
      {colorbrewer_3green_light,  "CD57"},
      {colorbrewer_3green_medium, "CD57"},
      {colorbrewer_3green_dark,   "CD57"}
    };
  }
  else if (module == "tls") {
    cm = {
      {color_deep_pink,         "CD3-only"},
      {color_purple,            "CD20"},
      {color_cyan,              "CD4"},
      {color_red,               "CD8"},
      {color_dark_blue,         "other TLS"}
    };
  }
  else if (module == "margin") {
    cm = {
      {colorbrewer_3red_light,  "Margin (Automated)"},
      {colorbrewer_3red_medium, "Margin (Manual)"},
      {colorbrewer_3red_dark,   "Margin (A+M)"}
    };
  }
  else if (module == "artifact") {
    cm = {
      {color_red, "PanCK CD3"}
    };
  }
  else if (module == "pdl1") {
    cm = {
      {color_light_red,  "PanCK PD-L1 neg"},
      {color_red,        "PanCK PD-L1 pos"},
      {color_light_green,"CD163 PD-L1 neg"},
      {color_dark_green, "CD163 PD-L1 pos"}
    };
  } else if (module == "prostateimmune") {
    cm = {
      {color_light_red,    "CD3+"},
      {color_yellow,       "CD8+"},
      {color_deep_pink,    "CD20+"},
      {color_light_green,  "FOXP3+"},
      {color_purple,       "TLS"}
    };
  } else if (module == "orion") {
      cm = {
	{color_red,        "T-cell PD-1 pos"},
	{color_light_red,  "T-cell PD-1 neg"},
	{color_purple,     "B-cell"},
	{color_dark_green, "PanCK - PD-L1 pos"},
	{color_light_green,"PanCK - PD-L1 neg"},
	{color_cyan,       "Other PD-L1 pos"},
	{color_deep_pink,  "FOXP3 pos"},
	{color_gray,       "Stroma"}
      };
  }
  else if (module == "prostate") {
      cm = {
	{color_light_red, "T-cell PD-1 pos"},
	{color_red,       "T-cell PD-1 neg"},
	{color_purple,    "B-cell"},
	{color_dark_green,"AMACR+"},
	{color_gray,      "Stromal"}
      };
  }
  else if (module == "orionimmune") {
      cm = {
	{color_light_red, "CD3+CD4+"},
	{color_red,       "CD3+CD8+"},
	{color_purple,    "CD3+ only"},
	{color_dark_green,"CD4+ only"},
	{color_dark_blue, "CD8+ only"}
      };
  }
  else if (module == "jhuorion") {
    cm = { 
	{color_red,        "CD3+PD1+"},
	{color_light_red,  "CD3+PD1-"},
	{color_purple,     "CD20+"},
	{color_light_green,"PanCK"},
	{color_cyan,       "CD163+"},
	{color_deep_pink,  "FOXP3+"},
	{color_gray,       "Stroma"}
    };
  }
  else {
    std::cerr << "cyftools png -- color_map.cpp -- error, unknown module " << module << std::endl;
    assert(false);
  }
  return cm;
}

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
