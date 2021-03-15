#ifndef RAYTRACING_D_Programming_Scratchapixel_Raytracing_LIGHTS_H_INCLUDED
#define RAYTRACING_D_Programming_Scratchapixel_Raytracing_LIGHTS_H_INCLUDED

#include "vector.h"
#include "matrix.h"
#include "constants.h"
#include "functions.h"
#include <tuple>    // for std::ignore




enum class LightType {POINT_LIGHT, DISTANT_LIGHT};




class Light
{
public:
    LightType type {LightType::POINT_LIGHT};
    Vec3_f_t color {1.0f};
    float intensity {1.0f};
    Matrix44_f_t lightToWorldMtx {Consts::kIdentity};

    Light() = default;                                                              /// default constructor: member vars get initialized. If other constructors are provided default constructor is not implicitly created
    virtual ~Light()
    {
    }
/// change position for point light, change direction for distant light
    virtual void transformLightToWorld(const Matrix44_f_t &l2wMtx)                  /// member functions inside the class definition are considered implicitly inline.
    {
        std::cout << "transformLightToWorld function of a Light class, which should be overriden" << '\n';
        std::cout << l2wMtx << '\n';
    }
/// get direction
    virtual Vec3_f_t getDirection(const Vec3_f_t &P)
    {
        ignore(P);
        std::cout << "getDirection function of a Light class, which should be overriden" << '\n';

        return Vec3_f_t{1};
    }
/// get intensity at hit point
    virtual float getIntensityAtPoint(const Vec3_f_t &P)
    {
        ignore(P);
        std::cout << "getIntensityAtPoint function of a Light class, which should be overriden" << '\n';

        return 1.0f;
    }
/// get distance to hit point
    virtual float getDistance(const Vec3_f_t &P)
    {
        ignore(P);
        std::cout << "getDistance function of a Light class, which should be overriden" << '\n';

        return 1.0f;
    }
/// print Light
   virtual void print()
    {
///        std::cout << "type: " << type << '\n';                                   /// no implicit conversion to int
        std::cout << ((type == LightType::POINT_LIGHT) ? "Point light created" : "Distant light created") << '\n';
        std::cout << "color: " << color << '\n';
        std::cout << "intensity: " << intensity << '\n';
    }
};




class PointLight : public Light
{
public:
    Vec3_f_t position {0};

    PointLight() = default;
    PointLight(const Vec3_f_t &c, const float &i, const Vec3_f_t &p): position{p}
    {
        type = LightType::POINT_LIGHT;
        color = c;                                                                  /// can be initialized through parent class constructor, not directly with this->...
        intensity = i;
    }
/// change position
    virtual void transformLightToWorld(const Matrix44_f_t &l2wMtx) override
    {
        lightToWorldMtx = l2wMtx;
        Vec3_f_t tempPosition = lightToWorldMtx.multPointMatrix(position);
        position = tempPosition;
    }
/// get direction
    virtual Vec3_f_t getDirection(const Vec3_f_t &P) override
    {
        return (P - position).normalised();
    }
/// get intensity at hit point
    virtual float getIntensityAtPoint(const Vec3_f_t &P) override
    {
        /// compute the square distance
        float r2 = (P - position).norm();

        return intensity / (4 * Consts::M_PI * r2);
    }
/// get distance to hit point
    virtual float getDistance(const Vec3_f_t &P) override
    {
        return (P - position).length();
    }
/// print Light
    virtual void print() override
    {
        Light::print();
        std::cout << "position: " << position << '\n';
        std::cout << '\n';
    }

};




class DistantLight : public Light
{
public:
    Vec3_f_t direction {0, 0, -1};

    DistantLight() = default;
    DistantLight(const Vec3_f_t &c, const float &i, const Vec3_f_t &d): direction{d}
    {
        type = LightType::DISTANT_LIGHT;
        color = c;                                                            /// can be initialized through parent class constructor, not directly with this->...
        intensity = i;
    }
/// change direction
    virtual void transformLightToWorld(const Matrix44_f_t &l2wMtx) override
    {
        lightToWorldMtx = l2wMtx;
        Vec3_f_t tempDirection = lightToWorldMtx.multVectorMatrix(direction);
        tempDirection.normalise();
        direction = tempDirection;
    }
/// get direction
    virtual Vec3_f_t getDirection(const Vec3_f_t &P) override
    {
        ignore(P);

        return direction;
    }
/// get intensity at hit point, same for distant light
    virtual float getIntensityAtPoint(const Vec3_f_t &P) override
    {
        ignore(P);

        return intensity;
    }
/// get distance to hit point
    virtual float getDistance(const Vec3_f_t &P) override
    {
        ignore(P);

        return Consts::kInfinity;
    }
/// print Light
    virtual void print() override
    {
        Light::print();
        std::cout << "direction: " << direction << '\n';
        std::cout << '\n';
    }

};


#endif // RAYTRACING_D_Programming_Scratchapixel_Raytracing_LIGHTS_H_INCLUDED
