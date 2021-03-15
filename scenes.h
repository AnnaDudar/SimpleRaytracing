#ifndef RAYTRACING_D_Programming_Scratchapixel_Raytracing_SCENES_H_INCLUDED
#define RAYTRACING_D_Programming_Scratchapixel_Raytracing_SCENES_H_INCLUDED

#include "lightsClass.h"
#include "objectsClass.h"
#include <vector>




void setSpheres(std::vector<std::unique_ptr<Light>> &lights, std::vector<std::unique_ptr<Object>> &objects, Options &options, const std::string &logFile);




#endif // RAYTRACING_D_Programming_Scratchapixel_Raytracing_SCENES_H_INCLUDED
