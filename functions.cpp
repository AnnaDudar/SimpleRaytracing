#include "functions.h"
#include "matrix.h"

#include <fstream>
#include <time.h>




/**
solveQuadratic function:
find the roots of the quadratic equation
**/
bool solveQuadratic(float a, float b, float c, float& x0, float& x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0 ) return false;
    else if (discr == 0) x0 = x1 = -0.5 * b / a;
    else
    {
        float q = (b > 0) ?
                -0.5 * (b + sqrt(discr)) :
                -0.5 * (b - sqrt(discr));
        x0 = q / a;
        x1 = c / q;
    }

    return true;
}




/**
buildRotateXMtx, buildRotateYMtx, buildRotateZMtx functions:
return 4x4 matrices used for rotation around X, Y or Z axis correspondingly
**/
Matrix44_f_t buildRotateXMtx(float deg)
{
	float rad = deg2rad(deg);
	Matrix44_f_t rXMtx = Matrix44_f_t(	1.0, 0.0, 0.0, 0.0,
                                        0.0, cos(rad), sin(rad), 0.0,
                                        0.0, -sin(rad), cos(rad), 0.0,
                                        0.0, 0.0, 0.0, 1.0  		);

	return rXMtx;
}

Matrix44_f_t buildRotateYMtx(float deg)
{
	float rad = deg2rad(deg);
	Matrix44_f_t rYMtx = Matrix44_f_t(	cos(rad), 0.0, -sin(rad), 0.0,
                                        0.0, 1.0, 0.0, 0.0,
                                        sin(rad), 0.0, cos(rad), 0.0,
                                        0.0, 0.0, 0.0, 1.0		    );

	return rYMtx;
}

Matrix44_f_t buildRotateZMtx(float deg)
{
	float rad = deg2rad(deg);
	Matrix44_f_t rZMtx = Matrix44_f_t(	cos(rad), sin(rad), 0.0, 0.0,
                                        -sin(rad), cos(rad), 0.0,  0.0,
                                        0.0,  0.0, 1.0, 0.0,
                                        0.0, 0.0, 0.0, 1.0		    );

	return rZMtx;
}

/**
buildTranslateMtx functions:
return 4x4 matrices used for translation
**/
Matrix44_f_t buildTranslateMtx(float Tx, float Ty, float Tz)
{
	Matrix44_f_t tMtx = Matrix44_f_t(	1.0, 0.0, 0.0, 0.0,
                                        0.0, 1.0, 0.0, 0.0,
                                        0.0, 0.0, 1.0, 0.0,
                                        Tx, Ty, Tz, 1.0		);

	return tMtx;
}

Matrix44_f_t buildTranslateMtx(float Tx)
{
	Matrix44_f_t tMtx = Matrix44_f_t(	1.0, 0.0, 0.0, 0.0,
                                        0.0, 1.0, 0.0, 0.0,
                                        0.0, 0.0, 1.0, 0.0,
                                        Tx, Tx, Tx, 1.0		);

	return tMtx;
}

/**
buildScaleMtx functions:
return 4x4 matrices used for scaling
**/
Matrix44_f_t buildScaleMtx(float sX, float sY, float sZ)
{
	Matrix44_f_t sMtx = Matrix44_f_t(	sX, 0.0, 0.0, 0.0,
                                        0.0, sY, 0.0, 0.0,
                                        0.0, 0.0, sZ, 0.0,
                                        0.0, 0.0, 0.0, 1.0		);

	return sMtx;
}

Matrix44_f_t buildScaleMtx(float sX)
{
	Matrix44_f_t sMtx = Matrix44_f_t(	sX, 0.0, 0.0, 0.0,
                                        0.0, sX, 0.0, 0.0,
                                        0.0, 0.0, sX, 0.0,
                                        0.0, 0.0, 0.0, 1.0		);

	return sMtx;
}




/**
reflect function:
computes reflection direction
**/
Vec3_f_t reflect(const Vec3_f_t& I, const Vec3_f_t& N)
{
    return I - 2 * N.dotProduct(I) * N;
}




/**
refract function:
computes refraction direction
**/
Vec3_f_t refract(const Vec3_f_t& I, const Vec3_f_t& N, const float ior)
{
    Vec3_f_t tempN = N;
    float cosI = clamp(-1, 1, I.dotProduct(tempN));
    float etaI = 1;
    float etaT = ior;

    /// incidence ray hits the surface from outside
    if (cosI < 0)
    {
        cosI = - cosI;                                                                              /// I.dotProduct(-N)
    }
    /// incidence ray leaves the surface from inside
    else
    {
        tempN = - tempN;
        std::swap(etaI, etaT);
    }

    float eta = etaI / etaT;
    float ksq = 1 - eta * eta * (1 - cosI * cosI);

    /// total internal reflection
    if (ksq < 0)
        return Vec3_f_t{0};

    float k = sqrtf(std::max(0.f, ksq));
    if (k < 0)
        return Vec3_f_t{0};

    return eta * I + (eta * cosI - k) * tempN;
}




/**
fresnel function:
evaluates Fresnel equation (ratio of reflected light for a given incident view direction, surface normal and surface refractive index)
**/
float fresnel(const Vec3_f_t& I, const Vec3_f_t& N, const float ior)
{
    float cosI = clamp(-1, 1, I.dotProduct(N));
    float etaI = 1;
    float etaT = ior;

    /// incidence ray hits the surface from outside
    if (cosI < 0)
    {
        cosI = - cosI;                                                                              /// I.dotProduct(-N)
    }
    /// incidence ray leaves the surface from inside
    else
    {
        std::swap(etaI, etaT);
    }

    /// compute sinT using Snell's law
    float sinT = etaI / etaT * sqrtf(std::max(0.f, 1 - cosI * cosI));

    /// total internal reflection
    float kr {0.f};
    if (sinT >= 1)
    {
        kr = 1;
    }
    else
    {
        float cosT = sqrtf(std::max(1.0f, 1 - sinT * sinT));
//        cosi = fabsf(cosi);
        float Rpar = ((etaT * cosI) - (etaI * cosT)) / ((etaT * cosI) + (etaI * cosT));
        float Rperp = ((etaI * cosI) - (etaT * cosT)) / ((etaI * cosI) + (etaT * cosT));
        kr = (Rpar * Rpar + Rperp * Rperp) / 2;
        // As a consequence of the conservation of energy, transmittance is given by:
        //kt = 1 - kr;
    }

    return kr;
}




/**
setLogName, writeLog functions:
set the name for a log file and write log info to it correspondingly
**/
void setLogName(char *outString)
{
    time_t t = time(0);                                                                             /// get time now
    struct tm *now = localtime(&t);                                                                 /// localtime() uses the time pointed by t, to fill a tm structure with the values that represent the corresponding local time.

    constexpr int dateBufLength = 31;                                                               /// set buffer for date
    char dateBuffer[dateBufLength];
    strftime (dateBuffer, dateBufLength, "./log/%Y_%m_%d_%H_%M_log.txt", now);                      /// write the values from tm structure to the date buffer

    for (int i = 0; i < dateBufLength; ++i)                                                         /// copy values from date buffer to outString
    {
        outString[i] = dateBuffer[i];
    }
}

void writeLog(const std::string &logFile, const std::string &text)
{
    std::ofstream logOfs;
try
{
    logOfs.open(logFile, std::ios::out | std::ios::binary | std::ios::app);
    if (logOfs.fail())
    {
        throw logFile;
    }

    logOfs << text << '\n';

    logOfs.close();
}
catch (const char *excpt)
{
    std::cerr << "Error: file open failed" << '\n';
    std::cerr << "File: " << excpt << '\n';
    logOfs.close();
}
}
