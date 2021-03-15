#ifndef RAYTRACING_D_Programming_Scratchapixel_Raytracing_OBJECTS_H_INCLUDED
#define RAYTRACING_D_Programming_Scratchapixel_Raytracing_OBJECTS_H_INCLUDED

#include "options.h"
#include "constants.h"
#include "raysClass.h"




enum class ObjectType {SPHERE, TRIANGLE, MESH};
enum class MaterialType {DIFFUSE_AND_GLOSSY, REFLECTION_AND_REFRACTION, REFLECTION};




struct MaterialProp
{
    MaterialType matType {MaterialType::DIFFUSE_AND_GLOSSY};
    /// diffuse and glossy
    float Kd {0.8f};
    float Ks {0.2f};
    Vec3_f_t albedo {0.18f};
    Vec3_f_t specularColor {1.f};
    float specularExp {25.f};
    /// reflection and refraction
    float ior {1.3f};

};




class Object
{
public:
    ObjectType objType {ObjectType::MESH};
    /// object-to-world transformation matrix
    Matrix44_f_t objectToWorldMtx {Consts::kIdentity};
    /// material properties
    MaterialProp materialProp {MaterialType::DIFFUSE_AND_GLOSSY, 0.8f, 0.2f, 0.18f, 1.f, 25.f, 1.3f};
    /// object options: if culling is applied, if smooth shading is applied
    ObjectOptions objectOpts {};
    /// hit point and its properties: uv, n, N, tex
    float t {0};
    Vec3_f_t P {0};
    /// barycentric coordinates
    Vec2_f_t uv {0};
    /// normal not normalised
    Vec3_f_t n {0};
    /// normal normalised
    Vec3_f_t N {0};
    /// texture coordinates
    Vec2_f_t tex {0};

public:
    Object() {}
    virtual ~Object()
    {
    }
    /// transform object with the object-to-world matrix
    virtual void transformObjToWorld (const Matrix44_f_t &o2wMtx);
    /// find a hit point
    virtual bool intersect(const Ray& ray);
    /// get hit point properties: N, tex
    virtual void getHitPointProperties();
    /// print object
    virtual std::string show();
};




class Triangle : public Object
{
public:
    Vec3_f_t v0 {0};
    Vec3_f_t v1 {0};
    Vec3_f_t v2 {0};
    Vec3_f_t n0 {0};
    Vec3_f_t n1 {0};
    Vec3_f_t n2 {0};
    Vec2_f_t tx0 {0};
    Vec2_f_t tx1 {0};
    Vec2_f_t tx2 {0};
//    Vec3_f_t clr0 {0};                                                                            /// only for investigation reasons: different colours for different vertices
//    Vec3_f_t clr1 {0};
//    Vec3_f_t clr2 {0};

public:
    Triangle()
    {
        objType = ObjectType::TRIANGLE;
    }
    Triangle (const Vec3_f_t &vtx0, const Vec3_f_t &vtx1, const Vec3_f_t &vtx2,
            const Vec3_f_t &norm0, const Vec3_f_t &norm1, const Vec3_f_t &norm2,
            const Vec2_f_t &tex0, const Vec2_f_t &tex1, const Vec2_f_t &tex2,
            const MaterialProp &matProp, const ObjectOptions &objOpts);
    virtual ~Triangle()
    {
    }
    /// transform object with the object-to-world matrix
    virtual void transformObjToWorld (const Matrix44_f_t &o2wMtx) override;
    /// find a hit point
    virtual bool intersect(const Ray &ray) override;
    /// get hit point properties: N, tex
    virtual void getHitPointProperties() override;
    /// print triangle
    void print()
    {
        std::cout << v0 << " " << v1 << " " << v2 << '\n';
    }
};




class Mesh : public Object
{
public:
    static uint32_t counter;
    uint32_t id {0};
    uint32_t trisNum {0};
    uint32_t pointsNum {0};
    uint32_t vertsNum {0};

//    std::unique_ptr<uint32_t[]> vertsInFaces {nullptr};
    std::unique_ptr<Vec3_f_t[]> points {nullptr};
    std::unique_ptr<uint32_t[]> vertsIndices {nullptr};
    std::unique_ptr<Vec3_f_t[]> normals {nullptr};
    std::unique_ptr<Vec2_f_t[]> texs {nullptr};

    Triangle hitTriangle {};                                                                        /// used by several member-functions

public:
    Mesh()
    {
        id = ++counter;
        std::cout << "Mesh " << id << " created" << '\n';
        objType = ObjectType::MESH;

    }
    Mesh(uint32_t fNum, std::unique_ptr<uint32_t[]> vInFArr,                                        /// using move semantics
           std::unique_ptr<Vec3_f_t[]> pArr, std::unique_ptr<uint32_t[]> vIndArr,
           std::unique_ptr<Vec3_f_t[]> nArr, std::unique_ptr<Vec2_f_t[]> stArr,
           const MaterialProp &matProp, const ObjectOptions &objOpts);

    ~Mesh()
    {
        std::cout << "Mesh " << id << " deleted" << '\n';
        std::cout << '\n';
    }
    Mesh (const Mesh &mesh) = delete;
    Mesh& operator = (const Mesh &mesh) = delete;
    /// transform object with the object-to-world matrix
    virtual void transformObjToWorld (const Matrix44_f_t &o2wMtx) override;
    /// find a hit point
    virtual bool intersect(const Ray &ray) override;
    /// get hit point properties: N, tex
    virtual void getHitPointProperties() override;
    /// print mesh
    virtual std::string show() override;
};




class Sphere : public Object
{
public:
    static uint32_t counter;
    uint32_t id {0};
    Vec3_f_t center {Vec3_f_t{0}};
    float radius {1.f};
    float radius2 {1.f};

public:
    Sphere()
    {
        id = ++counter;
        std::cout << "Sphere " << id << " created" << '\n';
        objType = ObjectType::SPHERE;
    }
    Sphere(const Vec3_f_t &c, const float r, const MaterialProp &matProp, const ObjectOptions &objOpts);
    virtual ~Sphere()
    {
        std::cout << "Sphere " << id << " deleted" << '\n';
        std::cout << '\n';
    }
    /// transform object with the object-to-world matrix
    virtual void transformObjToWorld (const Matrix44_f_t &o2wMtx) override;
    /// find a hit point
    virtual bool intersect (const Ray& ray) override;
    /// get hit point properties: N, tex
    virtual void getHitPointProperties() override;
    /// convert sphere to mesh
    std::unique_ptr<Mesh> tessellate(uint32_t divs);
    /// print sphere
    virtual std::string show() override;
};




#endif // RAYTRACING_D_Programming_Scratchapixel_Raytracing_OBJECTS_H_INCLUDED
