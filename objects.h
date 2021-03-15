#ifndef RAYTRACING_D_Programming_Scratchapixel_Raytracing_OBJECTS_H_INCLUDED
#define RAYTRACING_D_Programming_Scratchapixel_Raytracing_OBJECTS_H_INCLUDED

#include "options.h"
#include "constants.h"
#include "rays.h"




enum class ObjectType {SPHERE, TRIANGLE, MESH};
enum class MaterialType {DIFFUSE_AND_GLOSSY, REFLECTION_AND_REFRACTION, REFLECTION};




struct MaterialInfo
{
    MaterialType matType {MaterialType::DIFFUSE_AND_GLOSSY};
    /// diffuse and glossy
    float albedo {0.18f};
    float Kd {0.8};
    float Ks {0.2};
    float specularExponent {25};
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
    MaterialType matType {MaterialType::DIFFUSE_AND_GLOSSY};
    float albedo {0.18f};
    float Kd {0.8};
    float Ks {0.2};
    float specularExponent {25};
    float ior {1.3f};
    /// object options: if culling is applied, if smooth shading is applied
    ObjectOptions objectOpts {};
    /// hit point and its properties: uv, n, N, tex, color
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
    /// color of a hit point
    Vec3_f_t color {0};

public:
// mojno sdelat protected, chtoby snaruji dostupa ne bylo, no byl dostup klassam-naslednikam, no togda net dostupa iz main
//   Object(const ObjectType &tp, const Vec3_f_t &alb, bool cull): type{tp}, albedo{alb}, culling{cull}
    Object() {}
    virtual ~Object()                                                       // not necessary
    {
    }
    /// transform object with the object-to-world matrix
    virtual void transformObjToWorld (const Matrix44_f_t &o2wMtx);
    /// find a hit point
    virtual bool intersect(const Ray& ray); //, int& index
    /// get hit point properties: N, tex
    virtual void getHitPointProperties();  //const Vec3_f_t& I, const int& index
//    virtual Vec3_f_t evaluateDiffuseColor(const Vec2_f_t& ) const {return diffuseColor; }
    void printObjType(ObjectType objType)
    {
        switch (static_cast<int>(objType))
        {
        case 0:
            std::cout << "SPHERE" << '\n';
            break;
        case 1:
            std::cout << "TRIANGLE" << '\n';
            break;
        case 2:
            std::cout << "MESH" << '\n';
            break;
        }
    }
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
    Vec3_f_t clr0 {0};
    Vec3_f_t clr1 {0};
    Vec3_f_t clr2 {0};

public:
    Triangle()
    {
        objType = ObjectType::TRIANGLE;
    }
    Triangle (const Vec3_f_t &vtx0, const Vec3_f_t &vtx1, const Vec3_f_t &vtx2,
            const Vec3_f_t &norm0, const Vec3_f_t &norm1, const Vec3_f_t &norm2,
            const Vec2_f_t &tex0, const Vec2_f_t &tex1, const Vec2_f_t &tex2,
            const Vec3_f_t &col0, const Vec3_f_t &col1, const Vec3_f_t &col2,
            const MaterialInfo &matInfo, const ObjectOptions &objOpts);
    virtual ~Triangle()                                                       // not necessary
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
//    bool triangulated {true};
//    uint32_t facesNum {0};
    uint32_t trisNum {0};
    uint32_t pointsNum {0};
    uint32_t vertsNum {0};

//    std::unique_ptr<uint32_t[]> vertsInFaces {nullptr};
    std::unique_ptr<Vec3_f_t[]> points {nullptr};
    std::unique_ptr<uint32_t[]> vertsIndices {nullptr};
    //uint32_t *vertsIndices {nullptr};
    std::unique_ptr<Vec3_f_t[]> normals {nullptr};
    std::unique_ptr<Vec2_f_t[]> texs {nullptr};
    std::unique_ptr<Vec3_f_t[]> colors {nullptr};

    Triangle hitTriangle {};                                                /// used by several member-functions

public:
    Mesh()
    {
        id = ++counter;
        std::cout << "Mesh " << id << " created" << '\n';
        objType = ObjectType::MESH;

    }
    Mesh(uint32_t fNum, std::unique_ptr<uint32_t[]> vInFArr,                /// using move semantics
         std::unique_ptr<Vec3_f_t[]> pArr, std::unique_ptr<uint32_t[]> vIndArr,
         std::unique_ptr<Vec3_f_t[]> nArr, std::unique_ptr<Vec2_f_t[]> stArr, std::unique_ptr<Vec3_f_t[]> cArr,
         const MaterialInfo &matInfo, const ObjectOptions &objOpts);
    ~Mesh()                                                                 // not necessary?!
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
    void print();
/*
    Vec3_f_t evaluateDiffuseColor(Vec2_f_t& st) const
    {
        float scale = 5;
        float pattern = (fmodf (st.x * scale, 1) > 0.5) ^ (fmodf (st.y * scale, 1) > 0.5);

        return mix(Vec3_f_t(0.815, 0.235, 0.031), Vec3_f_t(0.937, 0.937, 0.231), pattern);
    }
*/
};




class Sphere : public Object
{
public:
    static uint32_t counter;
    uint32_t id {0};
    Vec3_f_t center {Vec3_f_t{0}};
    float radius {1.f};
    float radius2 {1.f};
    Vec3_f_t col {Vec3_f_t{1}};                                                       // this is a color of an object!

public:
    Sphere()
    {
        id = ++counter;
        std::cout << "Sphere " << id << " created" << '\n';
        objType = ObjectType::SPHERE;
    }
    Sphere(const Vec3_f_t &c, const float r, const Vec3_f_t &cls,
           const MaterialInfo &matInfo, const ObjectOptions &objOpts);
    virtual ~Sphere()
    {
        std::cout << "Sphere " << id << " deleted" << '\n';
        std::cout << '\n';
    }
    /// transform object with the object-to-world matrix
    virtual void transformObjToWorld (const Matrix44_f_t &o2wMtx) override;
    /// find a hit point
    virtual bool intersect (const Ray& ray) override; //, int& index, Vec2_f_t& uv
    /// get hit point properties: N, tex
    virtual void getHitPointProperties() override; //, const Vec3_f_t& I, const int& index, const Vec2_f_t& uv, , Vec2_f_t& st
    /// convert sphere to mesh
    std::unique_ptr<Mesh> tessellate(uint32_t divs);                                 // generatePolyShphere(float rad, uint32_t divs)
};




#endif // RAYTRACING_D_Programming_Scratchapixel_Raytracing_OBJECTS_H_INCLUDED
