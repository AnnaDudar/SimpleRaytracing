#include "scenes.h"
#include "lightsClass.h"
#include "objectsClass.h"
#include "functions.h"

#include <vector>
#include <fstream>




/**
setMeshFromFile function:
reads a .geo file into the Mesh stucture
**/
std::unique_ptr<Mesh> setMeshFromFile(const char* file, const MaterialProp &matProp, const ObjectOptions &objOpts)
{
    std::ifstream ifs;
    std::cout << "reading mesh from file..." << '\n';
try
{
    ifs.open(file);
    if (ifs.fail())
    {
        throw file;
    }
    std::stringstream ss;
    ss << ifs.rdbuf();

    /// Number of faces
    uint32_t facesNum;
    ss >> facesNum;

    /// reading vertsInFaces array
    std::unique_ptr<uint32_t[]> vertsInFaces {new uint32_t[facesNum]};
    uint32_t vertsNum = 0;
    for (uint32_t i = 0; i < facesNum; ++i)
    {
        ss >> vertsInFaces[i];
        /// compute number of verts
        vertsNum += vertsInFaces[i];
    }

    /// reading vertex index array
    std::unique_ptr<uint32_t[]> vertsIndices {new uint32_t[vertsNum]};
    uint32_t pointsNum = 0;
    for (uint32_t i = 0; i < vertsNum; ++i)
    {
        ss >> vertsIndices[i];
        /// compute number of points
        if (pointsNum < vertsIndices[i])
            pointsNum = vertsIndices[i];
    }
    pointsNum += 1;

    /// reading points array
    std::unique_ptr<Vec3_f_t[]> points {new Vec3_f_t[pointsNum]};
    for (uint32_t i = 0; i < pointsNum; ++i)
    {
        ss >> points[i].m_x >> points[i].m_y >> points[i].m_z;
    }
    /// reading normals array
    std::unique_ptr<Vec3_f_t[]> normals {new Vec3_f_t[vertsNum]};
    for (uint32_t i = 0; i < vertsNum; ++i)
    {
        ss >> normals[i].m_x >> normals[i].m_y >> normals[i].m_z;
    }
    /// reading st coordinates array
    std::unique_ptr<Vec2_f_t[]> texs {new Vec2_f_t[vertsNum]};
    for (uint32_t i = 0; i < vertsNum; ++i)
    {
        ss >> texs[i].m_x >> texs[i].m_y;
    }
    /// reading color array
    std::unique_ptr<Vec3_f_t[]> colors {new Vec3_f_t[vertsNum]};
    for (uint32_t i = 0; i < vertsNum; ++i)
    {
        ss >> colors[i].m_x >> colors[i].m_y >> colors[i].m_z;
    }

    ifs.close();

    return std::unique_ptr<Mesh> {new Mesh(facesNum, std::move(vertsInFaces), std::move(points), std::move(vertsIndices),
                    std::move(normals), std::move(texs), matProp, objOpts)};
}
catch (const char *excpt)
{
    std::cerr << "Error: file open failed" << '\n';
    std::cerr << "File: " << excpt << '\n';
    ifs.close();
}

    return nullptr;
}




/**
setSpheres function:
sets the scene to be rendered
**/
void setSpheres(std::vector<std::unique_ptr<Light>> &lights, std::vector<std::unique_ptr<Object>> &objects, Options &options, const std::string &logFile)
{

    std::stringstream ss;                                                                           /// stringstream for log file
    ss << "Scene: " << '\n';
    std::cout << "Scene: " << '\n';

    /// 1. lights
    /// distant light
    std::unique_ptr<Light> dLightPtr {new DistantLight(Vec3_f_t{0.2f, 0.2f, 0.19f}, 10.0f, Vec3_f_t{0.f, 0.f, -1.f})};
    Matrix44_f_t lightToWorldMtx {Consts::kIdentity};
    dLightPtr-> transformLightToWorld(lightToWorldMtx);
    ss << dLightPtr-> show();
    std::cout << dLightPtr-> show();

    lights.push_back(std::move(dLightPtr));                                                         /// be aware: dLightPtr is moved, cannot be used after this line!

    /// point light
    std::unique_ptr<Light> pLightPtr1 {new PointLight(Vec3_f_t(0.9f, 0.9f, 0.9f), 1000.0f, Vec3_f_t(0.f))};
    lightToWorldMtx = buildTranslateMtx(-5.f, 10.f, 5.f);
    pLightPtr1-> transformLightToWorld(lightToWorldMtx);
    ss << pLightPtr1-> show();
    std::cout << pLightPtr1-> show();

    lights.push_back(std::move(pLightPtr1));                                                        /// be aware: pLightPtr1 is moved, cannot be used after this line!

    /// point light
    std::unique_ptr<Light> pLightPtr2 {new PointLight(Vec3_f_t(0.9f, 0.9f, 0.9f), 1000.0f, Vec3_f_t(0.f))};
    lightToWorldMtx = buildTranslateMtx(8.f, -3.f, 1.f);
    pLightPtr2-> transformLightToWorld(lightToWorldMtx);
    ss << pLightPtr2-> show();
    std::cout << pLightPtr2-> show();

    lights.push_back(std::move(pLightPtr2));                                                        /// be aware: pLightPtr2 is moved, cannot be used after this line!

    /// 2. objects
    /// mesh from sphere
    ss << "Mesh from sphere" << '\n';
    std::cout << "Mesh from sphere" << '\n';
    float position {0.f};
    float radius {1.f};
    MaterialProp sphMeshMatProp {MaterialType::REFLECTION_AND_REFRACTION};                          /// only material type, rest is default
    ObjectOptions sphMeshOpts {true, false};                                                        /// culling and smooth shading

    std::unique_ptr<Sphere> spherePtr {new Sphere{position, radius, sphMeshMatProp, sphMeshOpts}};
    std::unique_ptr<Object> sphereMeshPtr {spherePtr->tessellate(8)};
    Matrix44_f_t objectToWorldMtx = buildScaleMtx(0.6f, 4.75f, 1.f) *  buildTranslateMtx(0.f, 1.75f, -7.f);
    sphereMeshPtr ->transformObjToWorld(objectToWorldMtx);
    ss << sphereMeshPtr-> show();
    std::cout << sphereMeshPtr-> show();

    objects.push_back(std::move(sphereMeshPtr));

    /// spheres
    ss << "Spheres" << '\n';
    std::cout << "Spheres" << '\n';
    int objTotal = 7;                                                                               /// total number of objects
    int objNum = 5;                                                                                 /// number of objects in outer loop
    int sphNum = 3;                                                                                 /// number of spheres in inner loop
    for (int i = 0; i < objTotal; ++i)
    {
        for (int j = 0; j < sphNum; ++j)
        {
            /// sphere 1
            int index = i * objNum + j;
            float r {0.75f};                                                                        /// radius of a sphere
            float rReduceFactor {0.02f};
            float radius {std::max(r * (1.f - rReduceFactor * index), 0.01f)};
            float outr {4.f};                                                                       /// radius of a circle
            float outrReduceFactor {0.025f};
            float outerRadius {std::max(outr * (1.f - outrReduceFactor * index), 0.01f)};
            float density {0.85f};
            float yReduceFactor{0.0005f};
            Vec3_f_t position { Vec3_f_t(static_cast<float>(outerRadius * cos(density * index)),
                                        index * (0.25f - yReduceFactor * index),
                                        static_cast<float>(outerRadius * sin(density * index)))};

            MaterialProp sphere1MatProp {MaterialType::DIFFUSE_AND_GLOSSY};
            sphere1MatProp.Kd = std::max(1.f - 0.35f * j, 0.f);
            sphere1MatProp.Ks = std::max(1.f - sphere1MatProp.Kd, 0.f);
            sphere1MatProp.albedo = Vec3_f_t{0.999f, 0.806f, 0.401f};
            sphere1MatProp.specularColor = Vec3_f_t{1.f};
            sphere1MatProp.specularExp = std::pow(15, j);
            ObjectOptions sphere1Opts {true, true};                                                 /// culling and smooth shading

            std::unique_ptr<Object> spherePtr1 {new Sphere(position, radius, sphere1MatProp, sphere1Opts)};
            Matrix44_f_t objectToWorldMtx = buildTranslateMtx(0.f, -1.85f, -7.f);                   /// sphere: only translation, no scaling or rotation!
            spherePtr1->transformObjToWorld(objectToWorldMtx);
            ss << spherePtr1-> show();
            std::cout << spherePtr1-> show();

            objects.push_back(std::move(spherePtr1));                                               /// unique_ptr when passed to a function as an argument uses move constructor which accepts only &&rref as an argument, so only temporary objects can be passed to a function

        }

        /// sphere 2
        int index = i * objNum + sphNum + 0;
        float r {0.75f};                                                                            /// radius of a sphere
        float rReduceFactor {0.02f};
        float radius {std::max(r * (1.f - rReduceFactor * index), 0.01f)};
        float outr {4.f};                                                                           /// radius of a circle
        float outrReduceFactor {0.025f};
        float outerRadius {std::max(outr * (1.f - outrReduceFactor * index), 0.01f)};
        float density {0.85f};
        float yReduceFactor{0.0005f};
        Vec3_f_t position { Vec3_f_t(static_cast<float>(outerRadius * cos(density * index)),
                                    index * (0.25f - yReduceFactor * index),
                                    static_cast<float>(outerRadius * sin(density * index)))};

        MaterialProp sphere2MatProp {MaterialType::DIFFUSE_AND_GLOSSY};
        sphere2MatProp.Kd = 0.5f;
        sphere2MatProp.Ks = 0.5f;
        sphere2MatProp.albedo = Vec3_f_t{0.75164f + 0.24725f, 0.60648f + 0.1995f, 0.22648f + 0.0745f}; /// gold: diffuse color
        sphere2MatProp.specularColor = Vec3_f_t{0.62828f, 0.5558f, 0.36607f};                       /// gold: specular color
        sphere2MatProp.specularExp = 2.15f;
        ObjectOptions sphere2Opts {true, true};                                                     /// culling and smooth shading

        std::unique_ptr<Object> spherePtr2 {new Sphere(position, radius, sphere2MatProp, sphere2Opts)};
        objectToWorldMtx = buildTranslateMtx(0.f, -1.85f, -7.f);                                    /// sphere: only translation, no scaling or rotation!
        spherePtr2->transformObjToWorld(objectToWorldMtx);
        ss << spherePtr2-> show();
        std::cout << spherePtr2-> show();

        objects.push_back(std::move(spherePtr2));                                                   /// unique_ptr when passed to a function as an argument uses move constructor which accepts only &&rref as an argument, so only temporary objects can be passed to a function

        /// sphere 3
        index = i * objNum + sphNum + 1;
        r = 0.75f;                                                                                  /// radius of a sphere
        rReduceFactor = 0.02f;
        radius = std::max(r * (1.f - rReduceFactor * index), 0.01f);
        outr = 4.f;                                                                                 /// radius of a circle
        outrReduceFactor = 0.025f;
        outerRadius = std::max(outr * (1.f - outrReduceFactor * index), 0.01f);
        density = 0.85f;
        yReduceFactor = 0.0005f;
        position = Vec3_f_t(static_cast<float>(outerRadius * cos(density * index)),
                                    index * (0.25f - yReduceFactor * index),
                                    static_cast<float>(outerRadius * sin(density * index)));

        MaterialProp sphere3MatProp {MaterialType::DIFFUSE_AND_GLOSSY};
        sphere3MatProp.Kd = 0.5f;
        sphere3MatProp.Ks = 0.5f;
        sphere3MatProp.albedo = Vec3_f_t{0.714f, 0.4284f, 0.1814f};                                 /// grey: diffuse color
        sphere3MatProp.specularColor = Vec3_f_t{0.3936f, 0.2719f, 0.1667f};                         /// white: specular color
        sphere3MatProp.specularExp = 6.52f;
        ObjectOptions sphere3Opts {true, true};                                                     /// culling and smooth shading

        std::unique_ptr<Object> spherePtr3 {new Sphere(position, radius, sphere3MatProp, sphere3Opts)};
        objectToWorldMtx = buildTranslateMtx(0.f, -1.85f, -7.f);                                    /// sphere: only translation, no scaling or rotation!
        spherePtr3->transformObjToWorld(objectToWorldMtx);
        ss << spherePtr3-> show();
        std::cout << spherePtr3-> show();

        objects.push_back(std::move(spherePtr3));                                                   /// unique_ptr when passed to a function as an argument uses move constructor which accepts only &&rref as an argument, so only temporary objects can be passed to a function
    }

    /// plain/mesh
    ss << "Plain/ mesh" << '\n';
    std::cout << "Plain/ mesh" << '\n';

    uint32_t plFacesNum = 2;
    uint32_t plPointsNum = 4;
    uint32_t plVertsNum = 8;
    /// vertsInFaces array
    std::unique_ptr<uint32_t[]> plVertsInFaces {new uint32_t[plFacesNum]    {4, 4}    };
    /// points array
    std::unique_ptr<Vec3_f_t[]> plPoints {new Vec3_f_t[plPointsNum]         {{-100, -3, 0}, {100, -3, 0}, {100, -3, -40}, {-100, -3, -40}}    };        /// center of the sphere plus radius
    /// vertex index array
    std::unique_ptr<uint32_t[]> plVertsIndices {new uint32_t[plVertsNum]    {0, 1, 2, 3, 1, 0, 3, 2}    };
    /// normals array
    std::unique_ptr<Vec3_f_t[]> plNormals {new Vec3_f_t[plVertsNum]         {{0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}}    };
    /// st coordinates array
    std::unique_ptr<Vec2_f_t[]> plTexs {new Vec2_f_t[plVertsNum]            {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}}    };
    MaterialProp plMeshMatProp {MaterialType::REFLECTION};
    plMeshMatProp.specularColor = Vec3_f_t{0.75f};
    ObjectOptions plMeshOpts {true, false};                                                         /// culling and smooth shading

    std::unique_ptr<Object> plainMeshPtr = std::unique_ptr<Mesh> {new Mesh(plFacesNum, std::move(plVertsInFaces), std::move(plPoints), std::move(plVertsIndices),
                                   std::move(plNormals), std::move(plTexs), plMeshMatProp, plMeshOpts)};
    objectToWorldMtx = Consts::kIdentity;
    plainMeshPtr ->transformObjToWorld(objectToWorldMtx);
    ss << plainMeshPtr-> show();
    std::cout << plainMeshPtr-> show();

    objects.push_back(std::move(plainMeshPtr));

    ss << "Wall/ mesh" << '\n';
    std::cout << "Wall/ mesh" << '\n';
    MaterialProp wallMatProp {MaterialType::DIFFUSE_AND_GLOSSY};
    wallMatProp.Kd = 1.f;
    wallMatProp.Ks = 1.f;
    ObjectOptions wallOpts {true, false};                                                            /// culling and smooth shading

    std::unique_ptr<Object> wallMeshPtr = setMeshFromFile("./scenes/inputWall.geo", wallMatProp, wallOpts);
    if (wallMeshPtr)
    {
        objectToWorldMtx = Consts::kIdentity;
        wallMeshPtr-> transformObjToWorld(objectToWorldMtx);
        ss << wallMeshPtr-> show();
        std::cout << wallMeshPtr-> show();

        objects.push_back(std::move(wallMeshPtr));
    }

    /// 3. options
    options.backgroundColor = Vec3_f_t{1.f};
    options.width = 1920;                                                                           // 1920 * 1080 = 2073600 pixel; 2073600 pixel * 3 bytes = 6220800 ~ 6075 Kb
    options.height = 1080;
    options.foV = 60.f;
    options.imageAspectRatio = options.width /  static_cast<float>(options.height);
    options.tgHalfFoV = tan(deg2rad(options.foV * 0.5f));
    options.cameraToWorldMtx = buildTranslateMtx(0.f, 0.f, 5.f);
    options.bias = 0.0001f;
    options.maxDepth = 5;
    ss << "Options:" << '\n';
    ss << "image width: " << options.width << '\n';
    ss << "image height: " << options.height << '\n';
    ss << "image aspect ratio: " << options.imageAspectRatio << '\n';
    ss << "field of view: " << options.foV << '\n';

    /// 4. write log file
    writeLog (logFile, ss.str());
}
