#include "color_map.h"
#include <stdexcept>
#include <string>
#include <iostream>
#include <cassert>

std::ostream& operator<<(std::ostream& os, const Color& obj) {
  os << obj.red << "," << obj.green << "," << obj.blue << "," << obj.alpha;
    return os;
}

std::ostream& operator<<(std::ostream& os, const CellColor& obj) {
    os << "pflag: " << obj.pflag << ", "
       << "cflag: " << obj.cflag << ", "
       << "Color: " << obj.c << ", " 
       << "label: " << obj.label;
    return os;
}


Color select_color(cy_uint cflag, cy_uint pflag, const ColorLabelVec& palette) {
  
  for (const auto& c : palette) {
    //    std::cerr << " cflag " << cflag << " c.cflag " << c.cflag << " ISC " <<
    // IS_FLAG_SET(cflag, c.cflag) << " pflag " << pflag << " c.pflag " <<
    // c.pflag << " ISP " << IS_FLAG_SET(pflag, c.pflag) << " c.c " << c.c << std::endl;
    bool cflag_pass = c.cflag == 0 || IS_FLAG_SET(cflag, c.cflag);
    bool pflag_pass = c.pflag == 0 || IS_FLAG_SET(pflag, c.pflag);
    if (cflag_pass && pflag_pass) {
      return c.c;
    }
  }
  return palette.back().c; // default is the last
}

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
ColorLabelVec ColorLabelVecForModule(const std::string& module) {

  ColorLabelVec cm;
  if (module == "tumor") {
    cm = {
      {TUMOR_FLAG,                     0, colorbrewer_3red_light,  "Tumor (Automated)"},
      {TUMOR_MANUAL_FLAG,              0, colorbrewer_3red_medium, "Tumor (Manual)"},
     {TUMOR_MANUAL_FLAG + TUMOR_FLAG, 0, colorbrewer_3red_dark,   "Tumor (A+M)"},
      {MARGIN_FLAG,                    0, color_cyan,              "Margin"},
      {0,                              0, colorbrewer_3blue_light, "Stroma"}
    };
  }
  else if (module == "tls") {
    cm = {
      {TLS_FLAG, PROSTATE_CD20,color_purple,            "CD20"},
      {TLS_FLAG, PROSTATE_CD4, color_cyan,              "CD4"},
      {TLS_FLAG, PROSTATE_CD8, color_red,               "CD8"},
      {TLS_FLAG, PROSTATE_CD3, color_deep_pink,           "CD3"},
      {TLS_FLAG, 0,            color_dark_blue,         "other TLS"},
      {0,0,color_gray,      "Other"}
    };
  }
  else if (module == "margin") {
    cm = {
      {MARGIN_MANUAL_FLAG + MARGIN_FLAG, 0, colorbrewer_3red_dark,   "Margin (A+M)"},
      {MARGIN_FLAG,                      0, colorbrewer_3red_light,  "Margin (Automated)"},
      {MARGIN_MANUAL_FLAG,               0, colorbrewer_3red_medium, "Margin (Manual)"},
      {0,0,color_gray,      "Other"}
    };
  }
  else if (module == "artifact") {
    cm = {
      {0, ORION_CD3 + ORION_PANCK, color_red, "PanCK CD3"}
    };
  }
  else if (module == "pdl1") {
    cm = {
      {0, ORION_PANCK + ORION_PDL1, color_red,        "PanCK PD-L1 pos"},
      {0, ORION_PANCK, color_light_red,  "PanCK PD-L1 neg"},
      {0, ORION_CD163 + ORION_PDL1, color_dark_green, "CD163 PD-L1 pos"},
      {0, ORION_CD163, color_light_green,"CD163 PD-L1 neg"},
      {0,0,color_gray,      "Stroma"}
    };
  } else if (module == "prostateimmune") {
    cm = {
      {0, PROSTATE_FOXP3, color_light_green,  "FOXP3+"},      
      {0, PROSTATE_CD8, color_red,            "CD8+"},
      {0, PROSTATE_CD3, color_deep_pink,      "CD3+"},      
      {0, PROSTATE_CD20, color_purple,        "CD20+"},
      {0,0,color_gray,      "Other"}
    };
  } else if (module == "orion") {
      cm = {
	{0, ORION_PD1 + ORION_CD3, color_red, "T-cell PD-1 pos"},
	{0, ORION_CD3, color_light_red,       "T-cell PD-1 neg"},
	{0, ORION_CD20, color_purple,     "B-cell"},
	{0, ORION_PANCK + ORION_PDL1, color_dark_green, "PanCK - PD-L1 pos"},
	{0, ORION_PANCK, color_light_green,"PanCK - PD-L1 neg"},
	{0, ORION_PDL1, color_cyan,       "Other PD-L1 pos"},
	{0, ORION_FOXP3, color_deep_pink,  "FOXP3 pos"},
	{0,0,color_gray,      "Stroma"}
      };
  }
  else if (module == "prostate") {
      cm = {
	{0, PROSTATE_PD1 + PROSTATE_CD3, color_light_red, "T-cell PD-1 pos"},
	{0, PROSTATE_CD3, color_red,       "T-cell PD-1 neg"},
	{0, PROSTATE_CD20, color_purple,    "B-cell"},
	{0, PROSTATE_AMCAR, color_dark_green,"AMACR+"},
	{0,0,color_gray,      "Stroma"}
      };
  }
  else if (module == "orionimmune") {
      cm = {
	{0, ORION_CD3 + ORION_CD4, color_light_red, "CD3+CD4+"},
	{0, ORION_CD3 + ORION_CD8, color_red,       "CD3+CD8+"},
	{0, ORION_CD3, color_purple,    "CD3+ only"},
	{0, ORION_CD4, color_dark_green,"CD4+ only"},
	{0, ORION_CD8, color_dark_blue, "CD8+ only"},
	{0,0,color_gray,      "Other"}
      };
  }
  else if (module == "jhuorion") {
    cm = { 
      {0, ORIONJHU_CD3 + ORIONJHU_PD1, color_red,        "CD3+PD1+"},
      {0, ORIONJHU_CD3, color_light_red,  "CD3+PD1-"},
      {0, ORIONJHU_CD20, color_purple,     "CD20+"},
      {0, ORIONJHU_PANCK, color_light_green,"PanCK"},
      {0, ORIONJHU_CD163, color_cyan,       "CD163+"},
      {0, ORIONJHU_FOXP3, color_deep_pink,  "FOXP3+"},
      {0,0,color_gray,       "Stroma"}
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
