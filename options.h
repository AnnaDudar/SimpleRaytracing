#ifndef RAYTRACING_D_Programming_Scratchapixel_Raytracing_OPTIONS_H_INCLUDED
#define RAYTRACING_D_Programming_Scratchapixel_Raytracing_OPTIONS_H_INCLUDED

#include "vector.h"
#include "matrix.h"
#include "constants.h"




struct Options
{
    int width {512};
    int height {512};
    float foV {90.0f};
    float imageAspectRatio {1.f};
    float tgHalfFoV {1.f};
    Matrix44_f_t cameraToWorldMtx {};
    Vec3_f_t backgroundColor {1.f};
    float bias {0.0f};
    int maxDepth {1};

};




struct ObjectOptions
{
    bool culling {true};
    bool smoothShading {true};
};
#endif // RAYTRACING_D_Programming_Scratchapixel_Raytracing_OPTIONS_H_INCLUDED
