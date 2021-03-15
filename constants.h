#ifndef RAYTRACING_D_Programming_Scratchapixel_Raytracing_CONSTANTS_H_INCLUDED
#define RAYTRACING_D_Programming_Scratchapixel_Raytracing_CONSTANTS_H_INCLUDED

#include "vector.h"
#include "matrix.h"
#include <limits>




namespace Consts
{
    constexpr float kEpsilon {1e-8};
    constexpr float kInfinity {std::numeric_limits<float>::max()};
    constexpr float PI {3.141592};
    constexpr float invPI {1 / PI};
    constexpr float PI_2 {2 * PI}; //6.28318530717959;
    const Matrix44_f_t kIdentity {Matrix44_f_t{}};

}




#endif // RAYTRACING_D_Programming_Scratchapixel_Raytracing_CONSTANTS_H_INCLUDED
