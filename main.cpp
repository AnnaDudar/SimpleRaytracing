#include "vector.h"
#include "matrix.h"
#include "constants.h"
#include "options.h"
#include "lightsClass.h"
#include "objectsClass.h"
#include "raysClass.h"
#include "functions.h"
#include "scenes.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>




/**
trace function:
checks if some of the objects in the scene is intersected by a given ray
**/
Object* trace(const Ray &ray, const std::vector<std::unique_ptr<Object>> &objects, float tNear)
{
    Object *tempPtr {nullptr};

    std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
    for( ; iter != objects.end(); ++iter)
    {
        /// if the ray intersects some object and this object is closer to the viewer than the others then this object is stored in the pointer
        if ((*iter)->intersect(ray) && (*iter)->t < tNear)
        {
            tempPtr = iter->get();
            tNear = (*iter)->t;
        }
    }

    return tempPtr;
}




/**
castRay function:
computes the color at the intersection point of a ray defined by a position and a direction. The function is recursive.
If the material of the intersected object is either reflective or transparent (reflective and refractive), then we compute the reflection/refraction direction and cast one/two new rays into the scene by calling the castRay() function recursively.
When the surface is transparent, we mix the reflection and refraction color using the result of the Fresnel equations.
If the surface is diffuse/glossy we use the Phong illumination model to compute the color at the intersection point.
**/
Vec3_f_t castRay(const Ray &ray,
                 const std::vector<std::unique_ptr<Object>> &objects,
                 const std::vector<std::unique_ptr<Light>> &lights,
                 const Options &options, int depth = 0)
{
    /// 1. find if the ray intersects some of the objects
    /// intersection distance from the ray origin to the hit point
    float tNear = Consts::kInfinity;
    /// find if some object is hit and save in a pointer
    Object *hitObjectPtr = trace(ray, objects, tNear);

    /// 2. if the ray intersects some object find the color of the hit point
    /// color of the hit point
    Vec3_f_t hitColor {0};
    if (hitObjectPtr)
    {
        /// return background color if depth value is reached
        if (depth > options.maxDepth - 1) return options.backgroundColor;                           /// -1 because one more trace check is performed

        /// find properties of the hit point (normal N, texture coordinates tex)
        hitObjectPtr->getHitPointProperties();

        switch (hitObjectPtr->materialProp.matType)
        {
            /// a. diffuse and glossy surface
            case MaterialType::DIFFUSE_AND_GLOSSY:
            {
                Vec3_f_t diffuseHitColor {0};
                Vec3_f_t specularHitColor {0};

                for(uint32_t i = 0; i < lights.size(); ++i)
                {
                    /// direction of light from the light source to P
                    Vec3_f_t L = lights[i]-> getDirection(hitObjectPtr->P);
                    /// shadow ray
                    Ray shadowRay {hitObjectPtr->P + hitObjectPtr->N * options.bias, -L};
                    tNear = lights[i]-> getDistance(hitObjectPtr->P);
                    float isVisible = trace(shadowRay, objects, tNear) ? 0 : 1;               /// pointers can be converted to bool
                    Vec3_f_t light = isVisible * lights[i]-> getIntensityAtPoint(hitObjectPtr->P) * lights[i]-> color;
                    /// facing ratio
                    float facingRatio = std::max(0.f, hitObjectPtr->N.dotProduct(-L));

                    /// shading ratio
                    Vec3_f_t shadingRatio = hitObjectPtr-> materialProp.albedo * Consts::invPI;     /// integral(omega)(albedo * incident light * N.L)* dw <= incident light:
                                                                                                    /// the amount of light that is reflected by the surface can't be greater
                                                                                                    /// than the amount of light that is effectively reaching the surface
                                                                                                    /// minus the amount of light that is absorbed
                    /// pattern
                    //float scale = 5;
                    // sin pattern
                    //float pattern = (0.5 * cos(hitObjectPtr-> tex.m_y * 2 * Consts::M_PI * scale) + 0.5) * (0.5 * sin(hitObjectPtr-> tex.m_x * 2 * Consts::M_PI * scale) + 0.5);
                    // checker board pattern
                    //float angle = deg2rad(0.f);
                    //float s = cos(angle) * hitObjectPtr-> tex.m_x - sin(angle) * hitObjectPtr-> tex.m_y;
                    //float t = sin(angle) * hitObjectPtr-> tex.m_x + cos(angle) * hitObjectPtr-> tex.m_y;
                    float pattern = 1.f; //(fmod(s * scale, 1) > 0.5)^(fmod(t * scale, 1) > 0.5);

                    /// diffuse color
                    diffuseHitColor.m_x += pattern * shadingRatio.m_x * light.m_x * facingRatio;    /// red component
                    diffuseHitColor.m_y += pattern * shadingRatio.m_y * light.m_y * facingRatio;    /// green component
                    diffuseHitColor.m_z += pattern * shadingRatio.m_z * light.m_z * facingRatio;    /// blue component

                    /// the Phong model for the specular reflection
                    /// what would be the reflection direction for this light ray
                    Vec3_f_t R = reflect(L, hitObjectPtr-> N);
                    float specularRatio = std::pow(std::max(0.f, R.dotProduct(-ray.direction)), hitObjectPtr-> materialProp.specularExp);
                    shadingRatio = hitObjectPtr-> materialProp.specularColor;                       /// specular color set to 1 by default

                    /// specular reflection (usually without a color)
                    specularHitColor.m_x += shadingRatio.m_x * light.m_x * specularRatio;
                    specularHitColor.m_y += shadingRatio.m_y * light.m_y * specularRatio;
                    specularHitColor.m_z += shadingRatio.m_z * light.m_z * specularRatio;
                }

                /// hit color is mixed of diffuse and specular colors
                hitColor = hitObjectPtr-> materialProp.Kd * diffuseHitColor + hitObjectPtr-> materialProp.Ks * specularHitColor;

                break;
            }

            /// b. reflective surface
            case MaterialType::REFLECTION:
            {
                /// direction of reflection ray
                Vec3_f_t R = reflect(ray.direction, hitObjectPtr-> N);                              /// reflect() function: for an incident vector, with the end at P!
                R.normalised();
                /// reflection ray
                Ray reflectionRay {hitObjectPtr-> P + hitObjectPtr-> N * options.bias, R};

                /// if reflection ray hits some of the objects then its hit color should be stored to the initial hit color
                if (trace(reflectionRay, objects, tNear))
                {
                    /// shading ratio
                    Vec3_f_t shadingRatio = hitObjectPtr-> materialProp.specularColor;              /// specular color set to 1 by default, shading ratio is optional here

                    /// castRay is called recursively (for the case of several reflective objects in the scene, depth is limited)
                    Vec3_f_t reflColor = castRay(reflectionRay, objects, lights, options, depth + 1);

                    /// hit color is a color of the reflected object
                    hitColor.m_x = shadingRatio.m_x * reflColor.m_x;
                    hitColor.m_y = shadingRatio.m_y * reflColor.m_y;
                    hitColor.m_z = shadingRatio.m_z * reflColor.m_z;
                }

                break;
            }

            /// c. transparent surface
            case MaterialType::REFLECTION_AND_REFRACTION:
            {
                /// ratio of reflected light calculated using Fresnel equation
                float kr = fresnel(ray.direction, hitObjectPtr-> N, hitObjectPtr-> materialProp.ior);

                /// for the case: not internal reflection (incidence ray is completely reflected)
                Vec3_f_t refractColor {};
                if (kr < 1)
                {
                    /// direction of refraction ray
                    Vec3_f_t T =  refract(ray.direction, hitObjectPtr-> N, hitObjectPtr-> materialProp.ior);
                    T.normalised();
                    /// refraction ray
                    Ray refractRay {hitObjectPtr-> P - options.bias * hitObjectPtr-> N, T};

                    /// castRay is called recursively (for the case of several transparent objects in the scene, depth is limited)
                    refractColor = castRay(refractRay, objects, lights, options, depth + 1);
                }
                /// direction of reflection ray
                Vec3_f_t R = reflect(ray.direction, hitObjectPtr-> N);
                R.normalised();
                /// reflection ray
                Ray reflectRay {hitObjectPtr-> P + options.bias * hitObjectPtr-> N, R};

                /// castRay is called recursively (for the case of several transparent objects in the scene, depth is limited)
                Vec3_f_t reflectColor = castRay(reflectRay, objects, lights, options, depth + 1);

                /// hit color is mixed of reflection and refraction colors using the index computed with Fresnel equation
                hitColor = kr * reflectColor + (1 - kr) * refractColor;

                break;
            }
        }
    }
    else
    {
        hitColor = options.backgroundColor;
    }

    return hitColor;
}




/**
render function:
creation of the frame buffer, filling it with the color values, storing frame buffer to the output file
**/
void render(const Options &options, const std::vector<std::unique_ptr<Object>> &objects, const std::vector<std::unique_ptr<Light>> &lights,
            std::string outputFile, const std::string &logFile)
{
    std::cout << "start rendering..." << '\n';
    auto timeStart {std::chrono::high_resolution_clock::now()};
    /// 1. create and fill a frame buffer
    Vec3_f_t *frameBuffer = new Vec3_f_t [options.width * options.height];
    /// first pixel of frame buffer
    Vec3_f_t *pix = frameBuffer;                                                                    /// frame buffer contains Vec3_f_t pixels
    /// set a ray
    Ray ray;
    /// set the ray origin
    ray.origin = options.cameraToWorldMtx.multPointMatrix(ray.origin);                              /// point of origin of a ray = camera position = (0,0)
    std::cout << "Origin world coord: " << ray.origin << '\n';
    std::cout <<  '\n';

    /// fill the frame buffer
    for (int j = 0; j < options.height; ++j)
        for (int i = 0; i < options.width; ++i)
        {
            /// generate the primary ray direction depending on a pixel location
            /// from raster space to NDC space: x E [0, img.width], y E [0, img.height] -> x, y E [0,1]
            Vec2_f_t ndcPnt{(i + 0.5f) / static_cast<float>(options.width), (j + 0.5f) / static_cast<float>(options.height)};
            /// from NDC to screen space: x, y E [0,1] -> x, y E [-1,1]
            Vec2_f_t screenPnt;
            screenPnt.m_x = 2 * ndcPnt.m_x  - 1;
            screenPnt.m_y = 1 - 2 * ndcPnt.m_y;
            /// from screen to camera space: x, y E [-1,1] -> x, y, z E R^3
            Vec3_f_t cameraPnt;
            cameraPnt.m_x = screenPnt.m_x * options.imageAspectRatio * options.tgHalfFoV;
            cameraPnt.m_y = screenPnt.m_y * options.tgHalfFoV;
            cameraPnt.m_z = -1;
            /// from camera to world space
            ray.direction = options.cameraToWorldMtx.multVectorMatrix(cameraPnt);
            ray.direction.normalised();

            /// write the result of the ray cast to a pixel in the frame buffer
            *(pix) = castRay(ray, objects, lights, options);
            ++pix;
        }

    /// 2. save frame buffer to the file
    std::ofstream ofs;
try
{
    ofs.open(outputFile, std::ios::out | std::ios::binary);
    if (ofs.fail())
    {
        throw outputFile;
    }

    ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
    for (int i = 0; i < options.height * options.width; ++i)
    {
        uint8_t r = static_cast<uint32_t>(255 * clamp(0, 1, frameBuffer[i].m_x));
        uint8_t g = static_cast<uint32_t>(255 * clamp(0, 1, frameBuffer[i].m_y));
        uint8_t b = static_cast<uint32_t>(255 * clamp(0, 1, frameBuffer[i].m_z));
        ofs << r << g << b;
    }

    ofs.close();
    delete[] frameBuffer;
}
catch (const char *excpt)
{
    std::cerr << "Error: file open failed" << '\n';
    std::cerr << "File: " << excpt << '\n';
    ofs.close();
    delete[] frameBuffer;
}

    auto timeEnd {std::chrono::high_resolution_clock::now()};
    auto passedTime {std::chrono::duration<double, std::milli> (timeEnd - timeStart).count()};

    std::cout << "end rendering" << '\n';
    std::cerr << "passed time: " << passedTime << '\n';
    std::cout << '\n';

    /// 3. write log file
    std::stringstream ss;
    ss << '\n';
    ss << "rendering time: " << passedTime / 1000 << " s" << '\n';

    int meshTotal {0};
    int sphTotal {0};
    std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
    for( ; iter != objects.end(); ++iter)
    {
        if ((*iter)-> objType == ObjectType::MESH)
            ++ meshTotal;
        else
            ++sphTotal;
    }
    ss << "total number of meshes: " << meshTotal << '\n';
    ss << "total number of spheres: " << sphTotal << '\n';

    writeLog (logFile, ss.str());

}




int main()
{
try
{
    /// 1. creating the scene (adding objects and lights)
    /// dynamic array of pointers to Lights
    std::vector<std::unique_ptr<Light>> lights;
    /// dynamic array of pointers to Objects
    std::vector<std::unique_ptr<Object>> objects;
    /// options
    Options options {};
    /// output file
    std::string outputFile {"outputSpheres.ppm"};
    /// log file (+ 1 for /0 termination)
    constexpr int length = 31;
    char nameBuffer[length];
    setLogName(nameBuffer);
    std::string logFile{nameBuffer};

    /// fill arrays
    setSpheres(lights, objects, options, logFile);

    /// 2. render
    render(options, objects, lights, outputFile, logFile);

}
catch (...)
{
    std::cerr << "Abnormal termination" << '\n';
}

    return 0;
}
