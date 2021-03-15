#ifndef RAYTRACING_D_Programming_Scratchapixel_Raytracing_STRUCTS_H_INCLUDED
#define RAYTRACING_D_Programming_Scratchapixel_Raytracing_STRUCTS_H_INCLUDED

#include "vector.h"

const float kInfinity = std::numeric_limits<float>::max();
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

struct Options
{
    int width;
    int height;
    float foV;
    Matrix44_f_t cameraToWorld;
    float imageAspectRatio;
    int maxDepth;
    Vec3_f_t backgroundColor;
    float bias;
};

#endif // RAYTRACING_D_Programming_Scratchapixel_Raytracing_STRUCTS_H_INCLUDED
