#include "objects.h"
#include "constants.h"
#include "functions.h"
#include "options.h"




/// class Object
/// transform object with the object-to-world matrix
void Object::transformObjToWorld (const Matrix44_f_t &o2wMtx)
{
    std::cout << "Object class should not be instantiated" << '\n';
    std::cout << o2wMtx << '\n';
}
/// find a hit point
bool Object::intersect(const Ray& ray) //, int& index
{
    std::cout << "Object class should not be instantiated" << '\n';
    std::cout << ray.origin << " " << ray.direction << '\n';

    return true;
}
/// get hit point properties: N, tex
void Object::getHitPointProperties()  //const Vec3_f_t& I, const int& index
{
}




/// class Triangle
/// constructor with parameters
Triangle::Triangle (const Vec3_f_t &vtx0, const Vec3_f_t &vtx1, const Vec3_f_t &vtx2,
                    const Vec3_f_t &norm0, const Vec3_f_t &norm1, const Vec3_f_t &norm2,
                    const Vec2_f_t &tex0, const Vec2_f_t &tex1, const Vec2_f_t &tex2,
                    const Vec3_f_t &col0, const Vec3_f_t &col1, const Vec3_f_t &col2,
                    const MaterialInfo &matInfo, const ObjectOptions &objOpts):
                    v0{vtx0}, v1{vtx1}, v2{vtx2},
                    n0{norm0}, n1{norm1}, n2{norm2},
                    tx0{tex0}, tx1{tex1}, tx2{tex2},
                    clr0{col0}, clr1{col1}, clr2{col2}
{
//    static int triangleConstrCount {0};
//    std::cout << "triangle constructor: " << triangleConstrCount++ << '\n';
//    std::cout <<  '\n';

    objType = ObjectType::TRIANGLE;
    matType = matInfo.matType;
    albedo = matInfo.albedo;
    Kd = matInfo.Kd;
    Ks = matInfo.Ks;
    specularExponent = matInfo.specularExponent;
    ior = matInfo.ior;
    objectOpts = objOpts;

//        std::cout << "Triangle is created" << '\n';
//        std::cout << '\n';
}

/// transform object with the object-to-world matrix
void Triangle::transformObjToWorld (const Matrix44_f_t &o2wMtx)
{
//    static int triangleTransform {0};
//    std::cout << "triangle transform: " << triangleTransform++ << '\n';
//    std::cout <<  '\n';

    objectToWorldMtx = o2wMtx;
    /// transform vertices
    Vec3_f_t tempVtx0 = objectToWorldMtx.multPointMatrix(v0);
    Vec3_f_t tempVtx1 = objectToWorldMtx.multPointMatrix(v1);
    Vec3_f_t tempVtx2 = objectToWorldMtx.multPointMatrix(v2);
    v0 = tempVtx0;
    v1 = tempVtx1;
    v2 = tempVtx2;
    /// transform normals
    Matrix44_f_t transNormMtx = objectToWorldMtx.invert().transpose();
    Vec3_f_t tempNorm0 = transNormMtx.multVectorMatrix(n0);
    Vec3_f_t tempNorm1 = transNormMtx.multVectorMatrix(n1);
    Vec3_f_t tempNorm2 = transNormMtx.multVectorMatrix(n2);
    n0 = tempNorm0;
    n1 = tempNorm1;
    n2 = tempNorm2;
}

/// find a hit point
bool Triangle::intersect(const Ray &ray)
{
//    static int triangleIntersect {0};
//    std::cout << "triangle intersect: " << triangleIntersect++ << '\n';
//    std::cout <<  '\n';

    /// compute plane's normal
    Vec3_f_t v0v1 = v1 - v0;
    Vec3_f_t v0v2 = v2 - v0;
    Vec3_f_t v1v2 = v2 - v1;

    n = v0v1.crossProduct(v0v2);
    float area2 = n.length();    /// / 2;        /// can be optimized

    /// Step 1: finding P
    /// check if ray and plane are parallel
    /// if n * d > 0 the triangle is backfacing
    /// if n * d ~ 0, the ray misses the triangle
    float nd = n.dotProduct(ray.direction);

    if ((objectOpts.culling) & (nd > Consts::kEpsilon))
        return false;

    if (fabs(nd) < Consts::kEpsilon)
        return false;

    /// compute t
    Vec3_f_t ov0 = v0 - ray.origin;
    t = ov0.dotProduct(n) / nd;
    /// check if the triangle is behind the ray
    if (t < 0)
        return false;
    /// compute the intersection point
    P = ray.origin + t * ray.direction;

    /// Step 2: inside-outside test
    Vec3_f_t insOut;
    /// edge 0
    Vec3_f_t v0P = P - v0;
    insOut = v0v1.crossProduct(v0P);
    if (n.dotProduct(insOut) < 0)
        return false;
    /// edge 1
    Vec3_f_t v1P = P - v1;
    insOut = v1v2.crossProduct(v1P);
    uv[0] = insOut.length()  / area2;    // / 2
    if (n.dotProduct(insOut) < 0)
        return false;
    /// edge 2
    Vec3_f_t v2P = P - v2;
    insOut = (-v0v2).crossProduct(v2P);
    uv[1] = insOut.length() / area2;  // / 2
    if (n.dotProduct(insOut) < 0)
        return false;

    return true;

}

/// get hit point properties: N, tex
void Triangle::getHitPointProperties()
{
//    static int triangleHitPointProp {0};
//    std::cout << "triangle get hit point properties: " << triangleHitPointProp++ << '\n';
//    std::cout <<  '\n';

    N = n.normalised();
//        N = uv.m_x * n0 + uv.m_y * n1 + (1 - uv.m_x - uv.m_y) * n2;
    tex = uv[0] * tx0 + uv[1] * tx1 + (1 - uv[0] - uv[1]) * tx2;
    color = uv[0] * clr0 + uv[1] * clr1 + (1 - uv[0] - uv[1]) * clr2;

//        std::cerr << "get hit point props Triangle: " << color << '\n';
}




/// class Mesh
uint32_t Mesh::counter{0};

/// constructor with parameters
Mesh::Mesh(uint32_t fNum, std::unique_ptr<uint32_t[]> vInFArr,                                  /// move constructor of a smart pointer is called, copy constructor is deleted, ownership transferred
           std::unique_ptr<Vec3_f_t[]> pArr, std::unique_ptr<uint32_t[]> vIndArr,
           std::unique_ptr<Vec3_f_t[]> nArr, std::unique_ptr<Vec2_f_t[]> stArr, std::unique_ptr<Vec3_f_t[]> cArr,
           const MaterialInfo &matInfo, const ObjectOptions &objOpts)
{
    id = ++counter;
    std::cout << "Mesh " << id << " created" << '\n';

    objType = ObjectType::MESH;
    matType = matInfo.matType;
    albedo = matInfo.albedo;
    Kd = matInfo.Kd;
    Ks = matInfo.Ks;
    specularExponent = matInfo.specularExponent;
    ior = matInfo.ior;
    objectOpts = objOpts;

    // vertsNum > pointsNum, for cube: 24 verts (6 faces with 4 verts), but 8 points
    // vInFArr.length = fNum; vIndArr.length = vertsNum; pArr.length = pointsNum; nArr.length = sArr.length = cArr.length = vertsNum)
//        facesNum = fNum;
    //memcpy(vertsInFaces.get(), vInFArr, sizeof(int) * facesNum);
    // find out how many triangles we need to create for this mesh
    uint32_t k = 0;
    if (vInFArr)
    {
        for (uint_fast32_t i = 0; i < fNum; ++i)
        {
    //            if ((vInFArr[i] > 3) && (triangulated = true))
    //                triangulated = false;
    //            vertsNum += vInFArr[i];
            trisNum += vInFArr[i] - 2;
            for (uint32_t j = 0; j < vInFArr[i]; ++j)                                           /// vnutri kajdogo face smotrim na index of a vertex i ishem max index
            {
                if (pointsNum < vIndArr[k + j])
                    pointsNum = vIndArr[k + j];
            }
            k += vInFArr[i];
        }
        pointsNum += 1;
        vertsNum = trisNum * 3;
    }

    /// transfer ownership to points, no changes needed
    if (pArr)
    {
         points = std::move(pArr);
    /// allocate memory for points array
//    points = std::unique_ptr<Vec3_f_t[]>{new Vec3_f_t[pointsNum]};
//    for (uint_fast32_t i = 0; i < pointsNum; ++i)
//        points[i] = pArr[i];

    }

    /// allocate memory for other arrays
    /// mesh is triangulated
    if (vIndArr && nArr && stArr && cArr)
    {
        vertsIndices = std::unique_ptr<uint32_t[]> {new uint32_t[vertsNum]};
        normals = std::unique_ptr<Vec3_f_t[]>{new Vec3_f_t[vertsNum]};
        texs = std::unique_ptr<Vec2_f_t[]>{new Vec2_f_t[vertsNum]};
        colors = std::unique_ptr<Vec3_f_t[]>{new Vec3_f_t[vertsNum]};

        uint32_t l = 0;
        for (uint32_t i = 0, k = 0; i < fNum; ++i)                                              /// for each face
        {
            for (uint32_t j = 0; j < vInFArr[i] - 2; ++j)                                       /// for each triangle in the face
            {
                vertsIndices[l] = vIndArr[k];
                vertsIndices[l + 1] = vIndArr[k + j + 1];
                vertsIndices[l + 2] = vIndArr[k + j + 2];
                if (nArr)
                {
                    normals[l] = nArr[k];
                    normals[l + 1] = nArr[k + j + 1];
                    normals[l + 2] = nArr[k + j + 2];
                }
                if (stArr)
                {
                    texs[l] = stArr[k];
                    texs[l + 1] = stArr[k + j + 1];
                    texs[l + 2] = stArr[k + j + 2];
                }
                if (cArr)
                {
                    colors[l] = cArr[k];
                    colors[l + 1] = cArr[k + j + 1];
                    colors[l + 2] = cArr[k + j + 2];
                }
                l += 3;
            }
            k += vInFArr[i];
        }
        // you can use move if the input geometry is already triangulated
        // N = std::move(normals); // transfer ownership
        // sts = std::move(st); // transfer ownership
    }
}                                                                                               /// parameters pointing to resource are destroyed here, resource does not exist any more

/// transform object with the object-to-world matrix
void Mesh::transformObjToWorld (const Matrix44_f_t &o2wMtx)
{
//    static int meshTransform {0};
//    std::cout << "mesh Transform: " << meshTransform++ << '\n';
//    std::cout <<  '\n';


    objectToWorldMtx = o2wMtx;
    // transform vertices
    for (uint_fast32_t i = 0; i < pointsNum; ++i)
    {
        Vec3_f_t tempVtx = objectToWorldMtx.multPointMatrix(points[i]);
        points[i] = tempVtx;
    }
    // transform normals
    Matrix44_f_t transNormMtx = objectToWorldMtx.invert().transpose();
    for (uint_fast32_t i = 0; i < vertsNum; ++i)
    {
        Vec3_f_t tempNorm = transNormMtx.multVectorMatrix(normals[i]);
        normals[i] = tempNorm;
        normals[i].normalise();
    }
}

/// find a hit point
bool Mesh::intersect(const Ray &ray)
{
//    static int meshIntersect {0};
//    std::cout << "mesh Intersect: " << meshIntersect++ << '\n';
//    std::cout <<  '\n';

    bool intersect {false};
    float tNear {Consts::kInfinity};                                                                    /// optimize? tNear can be passed with parameter, so no need to start with kInfinity!!!
//        uint32_t hitTriIndex {0};
    for (uint32_t k = 0; k < trisNum; ++k)
    {
        const Vec3_f_t &v0 = points[vertsIndices[k * 3]];
        const Vec3_f_t &v1 = points[vertsIndices[k * 3 + 1]];
        const Vec3_f_t &v2 = points[vertsIndices[k * 3 + 2]];
        const Vec3_f_t &n0 = normals[k * 3];
        const Vec3_f_t &n1 = normals[k * 3 + 1];
        const Vec3_f_t &n2 = normals[k * 3 + 2];
        const Vec2_f_t &tx0 = texs[k * 3];
        const Vec2_f_t &tx1 = texs[k * 3 + 1];
        const Vec2_f_t &tx2 = texs[k * 3 + 2];
        const Vec3_f_t &clr0 = colors[k * 3];
        const Vec3_f_t &clr1 = colors[k * 3 + 1];
        const Vec3_f_t &clr2 = colors[k * 3 + 2];

        Triangle temp {v0, v1, v2, n0, n1, n2, tx0, tx1, tx2, clr0, clr1, clr2, MaterialInfo {matType, albedo, ior}, objectOpts};    /// vertices will be used further, so no need to transfer ownership, so passed by reference
        if (temp.intersect(ray) && temp.t < tNear)
        {
            tNear = temp.t;
            //hitTriIndex = k;
            //tempTriPtr = &tri;                                                                        /// bad idea, tri will be out of scope at the end of iteration, so tempTriPtr will be dangling
            hitTriangle = temp;

            intersect |= true;                                                                          /// bit XOR 0^1=1, 1^0=1, 0^0=0, 1^1=0
        }
    }
    if (intersect)
    {
        t = hitTriangle.t;
        P = hitTriangle.P;
        uv = hitTriangle.uv;
        n = hitTriangle.n;

/*      //slower with index and object is lost, should be created once again and intersect() called once more
        const Vec3_f_t &v0 = points[vertsIndices[hitTriIndex * 3]];
        const Vec3_f_t &v1 = points[vertsIndices[hitTriIndex * 3 + 1]];
        const Vec3_f_t &v2 = points[vertsIndices[hitTriIndex * 3 + 2]];
        const Vec3_f_t &n0 = normals[hitTriIndex * 3];
        const Vec3_f_t &n1 = normals[hitTriIndex * 3 + 1];
        const Vec3_f_t &n2 = normals[hitTriIndex * 3 + 2];
        const Vec2_f_t &tx0 = texs[hitTriIndex * 3];
        const Vec2_f_t &tx1 = texs[hitTriIndex * 3 + 1];
        const Vec2_f_t &tx2 = texs[hitTriIndex * 3 + 2];
        const Vec3_f_t &clr0 = colors[hitTriIndex * 3];
        const Vec3_f_t &clr1 = colors[hitTriIndex * 3 + 1];
        const Vec3_f_t &clr2 = colors[hitTriIndex * 3 + 2];

        hitTriangle = Triangle {v0, v1, v2, n0, n1, n2, tx0, tx1, tx2, clr0, clr1, clr2};
        hitTriangle.intersect(ray);                 // write hit t, P and uv coords
        hitTriangle.getHitPointProperties();        // write hit normal, tex coords and color
*/
    }

    return intersect;
}

/// get hit point properties: N, tex
void Mesh::getHitPointProperties()
{
//    static int meshHitPointProp {0};
//    std::cout << "mesh get hit point properties: " << meshHitPointProp++ << '\n';
//    std::cout <<  '\n';

    hitTriangle.getHitPointProperties();
    N = hitTriangle.N;
    tex = hitTriangle.tex;
    color = hitTriangle.color;

//        std::cerr << "get hit point props Mesh: " << color << '\n';
}

/// print mesh
void Mesh::print()
{
    std::cerr << "printing mesh..." << '\n';
    std::cerr << "id: " << id << '\n';
    std::cerr << "trisNum: " << trisNum << '\n';
    std::cerr << "vertsNum: " << vertsNum << '\n';
    std::cerr << "pointsNum: " << pointsNum << '\n';
/*
    std::cerr << "vertsInFaces: " << '\n';
    for (uint32_t i = 0; i < facesNum; ++i)
        std::cerr << vertsInFaces[i] << '\t';
    std::cerr << '\n';
*/
    std::cerr << "points: " << '\n';
    for (uint32_t i = 0; i < pointsNum; ++i)
        std::cerr << points[i] << '\t';
    std::cerr << '\n';

    std::cerr << "vertsIndices: " << '\n';
    for (uint32_t i = 0; i < vertsNum; ++i)
        std::cerr << vertsIndices[i] << '\t';
    std::cerr << '\n';

    std::cerr << "normals: " << '\n';
    for (uint32_t i = 0; i < vertsNum; ++i)
        std::cerr << normals[i] << '\t';
    std::cerr << '\n';

    std::cerr << "texs: " << '\n';
    for (uint32_t i = 0; i < vertsNum; ++i)
        std::cerr << texs[i] << '\t';
    std::cerr << '\n';

    std::cerr << "colors: " << '\n';
    for (uint32_t i = 0; i < vertsNum; ++i)
        std::cerr << colors[i] << '\t';
    std::cerr << '\n';

    std::cout << '\n';
}




/// class Sphere
uint32_t Sphere::counter{0};

/// constructor with parameters
Sphere::Sphere(const Vec3_f_t& c, const float r, const Vec3_f_t &cls,
       const MaterialInfo &matInfo, const ObjectOptions &objOpts):
        center{c}, radius{r}, radius2{r * r}, col{cls}
{
    id = ++counter;
    std::cout << "Sphere " << id << " created" << '\n';

    objType = ObjectType::SPHERE;
    matType = matInfo.matType;
    albedo = matInfo.albedo;
    Kd = matInfo.Kd;
    Ks = matInfo.Ks;
    specularExponent = matInfo.specularExponent;
    ior = matInfo.ior;
    objectOpts = objOpts;
//        std::cout << "Sphere is created" << '\n';
//        std::cout << '\n';
}

/// transform object with the object-to-world matrix
void Sphere::transformObjToWorld (const Matrix44_f_t &o2wMtx)
{
//    static int sphereTransform {0};
//    std::cout << "sphere Transform: " << sphereTransform++ << '\n';
//    std::cout <<  '\n';

    objectToWorldMtx = o2wMtx;
    // transform center
    Vec3_f_t tempCenter = objectToWorldMtx.multPointMatrix(center);
    center = tempCenter;
}

/// find a hit point
bool Sphere::intersect (const Ray& ray) //, int& index, Vec2_f_t& uv
{
//    static int sphereIntersect {0};
//    std::cout << "sphere intersect: " << sphereIntersect++ << '\n';
//    std::cout <<  '\n';

    //analytic solution: |O + tD|^2 - R^2 = 0 - vektorno-parametricheskoe uravnenie, iz P^2 - R^2 = 0
    float t0, t1;

    Vec3_f_t O = ray.origin - center;
    float a = (ray.direction).dotProduct(ray.direction);
    float b = 2 * (ray.direction).dotProduct(O);
    float c = O.dotProduct(O) - radius2;
//        std::cout << '\n';
//        std::cout << "a: " << a << '\n';
//        std::cout << "b: " << b << '\n';
//        std::cout << "c: " << c << '\n';

    if (!solveQuadratic(a, b, c, t0, t1))
        return false;
    if (t0 > t1)
        std::swap(t0, t1);
    // if t0 is negative, let's use t1 instead
    if (t0 < 0)
        t0 = t1;
    // both t0 and t1 are negative
    if (t0 < 0)
        return false;   ///!!

    t = t0;
    // orig is a point, dir is a vector, tnear is a scalar, result is a vector (dir * tnear.x - orig.x, dir * tnear.y - orig.y, dir * tnear.z - orig.z)
    P = ray.origin + ray.direction * t;

    return true;
}

/// get hit point properties: N, tex
void Sphere::getHitPointProperties() //, const Vec3_f_t& I, const int& index, const Vec2_f_t& uv, , Vec2_f_t& st
{
//    static int sphereHitPointProp {0};
//    std::cout << "sphere get hit point properties: " << sphereHitPointProp++ << '\n';
//    std::cout <<  '\n';

    n = (P - center);
    N = n.normalised();
    // In this particular case, the normal is simular to a point on a unit sphere
    // centred around the origin. We can thus use the normal coordinates to compute
    // the spherical coordinates of Phit.
    // atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
    // acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
    tex.m_x = atan2(N.m_z, N.m_x);
    tex.m_y = acosf(N.m_y);
    tex.m_x = (1 + tex.m_x / Consts::M_PI) * 0.5;
    tex.m_y = tex.m_y / Consts::M_PI;
    color = col;
}

/// convert sphere to mesh                                                                      /// ochen dolgo!
std::unique_ptr<Mesh> Sphere::tessellate(uint32_t divs)     // generatePolyShphere(float rad, uint32_t divs)
{
//    static int sphereTesselate {0};
//    std::cout << "sphere tesselate: " << sphereTesselate++ << '\n';
//    std::cout <<  '\n';

    ///generate points, norms, texs and colors
    uint32_t pointsNum = (divs - 1) * divs + 2;                                                 /// 2 sloya (polnye secheniya) * 3 dolki (polsecheniya) + 2 vershiny

    std::unique_ptr<Vec3_f_t[]> points{new Vec3_f_t[pointsNum]};
    std::unique_ptr<Vec3_f_t[]> norms {new Vec3_f_t[pointsNum]};
    std::unique_ptr<Vec2_f_t[]> txs {new Vec2_f_t[pointsNum]};

    float th = -Consts::M_PI;
    float ph = -Consts::M_PI_2;
    float dTh = Consts::M_PI / divs;
    float dPh = Consts::M_PI_2 / divs;

    points[0] = norms[0] = Vec3_f_t(0, -radius, 0) + center;                                    /// lower vertex
    uint32_t p = 1;
    for (uint32_t i = 0; i < divs - 1; i++)                                                     /// po gorizontali = shirota
    {
        th += dTh;
        ph = -Consts::M_PI;                          ///mojno ubrat?!
        for (uint32_t j = 0; j < divs; j++)                                                     /// po vertikali = dolgota
        {
            float x = radius * sin(th) * sin(ph) + center.m_x;
            float y = radius * cos(th) + center.m_y;
            float z = radius * sin(th) * cos(ph) + center.m_z;
            //float x = radius * cos(th) * cos(ph);
            //float y = radius * sin(th);
            //float z = radius * cos(th) * sin(ph);

            points[p] = norms[p] = Vec3_f_t(x, y, z);
            txs[p].m_x = th / Consts::M_PI + 0.5;                                               /// perevod [0; 2Pi] v interval [0; 1]
            txs[p].m_y = ph * 0.5 / Consts::M_PI + 0.5;                                         /// perevod [-Pi; Pi] v interval [0; 1]

            ph += dPh;
            p++;
        }
    }
    points[p] = norms[p] = Vec3_f_t(0, radius, 0) + center;

    //create connectivity list, verts in faces and verts indices arrays
    uint32_t facesNum = divs * divs;
    std::unique_ptr<uint32_t[]> vertsInFaces{new uint32_t[facesNum]};

    uint32_t vertsNum = (6 + (divs - 2) * 4) * divs;                                            /// 2 fans sverhu i snizu + quads, vse eto umnojit na kolichestvo vertikalnh polos
    std::unique_ptr<uint32_t[]> vertsIndices{new uint32_t[vertsNum]};                           /// obshee kolichestvo vershin
    std::unique_ptr<Vec3_f_t[]> colors {new Vec3_f_t[vertsNum]};

    uint32_t vid = 1, vNum = 0, l = 0;
    p = 0;
    for (uint32_t i = 0; i < divs; i++)                                                         /// kolichestvo gorizontalnyh polosok
    {
        for (uint32_t j = 0; j < divs; j++)                                                     /// kolichestvo polinomov v poloske
        {
        if (i == 0)                                                                             /// esli verhnyaa poloska
        {
            vertsInFaces[p++] = 3;
            vertsIndices[l] = 0;
            vertsIndices[l + 1] = (j == (divs - 1)) ? vid : j + vid + 1;
            vertsIndices[l + 2] = j + vid;
            colors[l] = col;
            colors[l + 1] = col;    // + Vec3_f_t(0.1);
            colors[l + 2] = col;    // - Vec3_f_t(0.1);
            l += 3;
        }
        else if (i == divs - 1)                                                                 /// esli nijnyaya poloska
        {
            vertsInFaces[p++] = 3;
            vertsIndices[l] = j + vid + 1 - divs;
            vertsIndices[l + 1] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs;
            vertsIndices[l + 2] = vid + 1;
            colors[l] = col;
            colors[l + 1] = col;    // + Vec3_f_t(0.1);
            colors[l + 2] = col;    // - Vec3_f_t(0.1);

            l += 3;
        }
        else                                                                                    /// esli centralnye poloski
        {
            vertsInFaces[p++] = 4;
            vertsIndices[l] = j + vid + 1 - divs;
            vertsIndices[l + 1] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs;
            vertsIndices[l + 2] = (j == (divs - 1)) ? vid + 1 : j + vid + 2;
            vertsIndices[l + 3] = j + vid + 1;
            colors[l] = col;
            colors[l + 1] = col;    // + Vec3_f_t(0.1);
            colors[l + 2] = col;    // - Vec3_f_t(0.1);
            colors[l + 3] = col;    // - Vec3_f_t(0.2);

            l += 4;
        }
        vNum++;
        }
    vid = vNum;
    }
    // rewrite arrays, now norm, tex and color info for each vert, not for each point
    std::unique_ptr<Vec3_f_t[]> normals {new Vec3_f_t[vertsNum]};
    std::unique_ptr<Vec2_f_t[]> texs {new Vec2_f_t[vertsNum]};

    for (uint32_t i = 0; i < vertsNum; ++i)
    {
        normals[i] = norms[vertsIndices[i]];
        texs[i] = txs[vertsIndices[i]];
    }

    return std::unique_ptr<Mesh> {new Mesh(facesNum, std::move(vertsInFaces), std::move(points), std::move(vertsIndices),
                    std::move(normals), std::move(texs), std::move(colors), MaterialInfo {matType, albedo, ior}, objectOpts)};

}
