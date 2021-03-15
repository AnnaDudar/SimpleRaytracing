#ifndef RAYTRACING_D_Programming_Scratchapixel_Raytracing_FUNCTIONS_H_INCLUDED
#define RAYTRACING_D_Programming_Scratchapixel_Raytracing_FUNCTIONS_H_INCLUDED

#include "constants.h"




inline
float clamp(const float& low, const float &high, const float &v)
{
    return std::max(low, std::min(high, v));
}

inline
float deg2rad(const float &deg)
{
    return deg * Consts::PI / 180;
}

inline
Vec3_f_t mix(const Vec3_f_t& a, const Vec3_f_t& b, const float& mixValue)
{
    return a * mixValue + b * (1 - mixValue);
}

template <typename T>
void ignore (T &/*someVar*/)
{
}




bool solveQuadratic(float a, float b, float c, float& x0, float& x1);

Matrix44_f_t buildRotateXMtx(float rad);

Matrix44_f_t buildRotateYMtx(float rad);

Matrix44_f_t buildRotateZMtx(float rad);

Matrix44_f_t buildTranslateMtx(float x, float y, float z);

Matrix44_f_t buildTranslateMtx(float x);

Matrix44_f_t buildScaleMtx(float scaleX, float scaleY, float scaleZ);

Matrix44_f_t buildScaleMtx(float scaleX);

Vec3_f_t reflect(const Vec3_f_t& I, const Vec3_f_t& N);

Vec3_f_t refract(const Vec3_f_t& I, const Vec3_f_t& N, const float ior);

float fresnel(const Vec3_f_t& I, const Vec3_f_t& N, const float ior);

void setLogName(char *outString);

void writeLog(const std::string &logFile, const std::string &text);

#endif // RAYTRACING_D_Programming_Scratchapixel_Raytracing_FUNCTIONS_H_INCLUDED
