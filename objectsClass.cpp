#include "objectsClass.h"
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
bool Object::intersect(const Ray& ray)
{
    std::cout << "Object class should not be instantiated" << '\n';
    std::cout << ray.origin << " " << ray.direction << '\n';

    return true;
}
/// get hit point properties: N, tex
void Object::getHitPointProperties()
{
}
/// print object
std::string Object::show()
{
    return std::string{};
}




/// class Triangle
/// constructor with parameters
Triangle::Triangle (const Vec3_f_t &vtx0, const Vec3_f_t &vtx1, const Vec3_f_t &vtx2,
                    const Vec3_f_t &norm0, const Vec3_f_t &norm1, const Vec3_f_t &norm2,
                    const Vec2_f_t &tex0, const Vec2_f_t &tex1, const Vec2_f_t &tex2,
                    const MaterialProp &matProp, const ObjectOptions &objOpts):
                    v0{vtx0}, v1{vtx1}, v2{vtx2},
                    n0{norm0}, n1{norm1}, n2{norm2},
                    tx0{tex0}, tx1{tex1}, tx2{tex2}
{
    objType = ObjectType::TRIANGLE;
    materialProp = matProp;
    objectOpts = objOpts;
}

/// transform object with the object-to-world matrix
void Triangle::transformObjToWorld (const Matrix44_f_t &o2wMtx)
{
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
    /// compute plane's normal
    Vec3_f_t v0v1 = v1 - v0;
    Vec3_f_t v0v2 = v2 - v0;
    Vec3_f_t v1v2 = v2 - v1;

    n = v0v1.crossProduct(v0v2);
    float area2 = n.length();                // / 2;        /// can be optimized

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
    uv[0] = insOut.length()  / area2;   // / 2
    if (n.dotProduct(insOut) < 0)
        return false;
    /// edge 2
    Vec3_f_t v2P = P - v2;
    insOut = (-v0v2).crossProduct(v2P);
    uv[1] = insOut.length() / area2;    // / 2
    if (n.dotProduct(insOut) < 0)
        return false;

    return true;

}

/// get hit point properties: N, tex
void Triangle::getHitPointProperties()
{
    /// smooth shading: vertex normals are used and then interpolated
    if (objectOpts.smoothShading)
    {
        N = uv.m_x * n0 + uv.m_y * n1 + (1 - uv.m_x - uv.m_y) * n2;
        N.normalise();
    }
    /// faceted shading: plane normal is used
    else
    {
        N = n.normalised();
    }
    tex = uv[0] * tx0 + uv[1] * tx1 + (1 - uv[0] - uv[1]) * tx2;
}




/// class Mesh
uint32_t Mesh::counter{0};

/// constructor with parameters
Mesh::Mesh(uint32_t fNum, std::unique_ptr<uint32_t[]> vInFArr,                                      /// move constructor of a smart pointer is called, copy constructor is deleted, ownership transferred
           std::unique_ptr<Vec3_f_t[]> pArr, std::unique_ptr<uint32_t[]> vIndArr,
           std::unique_ptr<Vec3_f_t[]> nArr, std::unique_ptr<Vec2_f_t[]> stArr,
           const MaterialProp &matProp, const ObjectOptions &objOpts)
{
    id = ++counter;
    objType = ObjectType::MESH;
    materialProp = matProp;
    objectOpts = objOpts;

    /// find out how many triangles are needed to create for this mesh
    uint32_t k = 0;
    if (vInFArr)
    {
        for (uint_fast32_t i = 0; i < fNum; ++i)
        {
            trisNum += vInFArr[i] - 2;
            /// inside each face check the index of a vertex and find the max index
            for (uint32_t j = 0; j < vInFArr[i]; ++j)
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
    if (vIndArr && nArr && stArr)
    {
        vertsIndices = std::unique_ptr<uint32_t[]> {new uint32_t[vertsNum]};
        normals = std::unique_ptr<Vec3_f_t[]>{new Vec3_f_t[vertsNum]};
        texs = std::unique_ptr<Vec2_f_t[]>{new Vec2_f_t[vertsNum]};

        uint32_t l = 0;
        for (uint32_t i = 0, k = 0; i < fNum; ++i)                                                  /// for each face
        {
            for (uint32_t j = 0; j < vInFArr[i] - 2; ++j)                                           /// for each triangle in the face
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
                l += 3;
            }
            k += vInFArr[i];
        }
    }
}                                                                                                   /// parameters pointing to resource are destroyed here, resource does not exist any more

/// transform object with the object-to-world matrix
void Mesh::transformObjToWorld (const Matrix44_f_t &o2wMtx)
{
    objectToWorldMtx = o2wMtx;
    /// transform vertices
    for (uint_fast32_t i = 0; i < pointsNum; ++i)
    {
        Vec3_f_t tempVtx = objectToWorldMtx.multPointMatrix(points[i]);
        points[i] = tempVtx;
    }
    /// transform normals
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
    bool intersect {false};
    float tNear {Consts::kInfinity};                                                                /// optimize? tNear can be passed with parameter, so no need to start with kInfinity!!!

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

        Triangle temp {v0, v1, v2, n0, n1, n2, tx0, tx1, tx2, materialProp, objectOpts};            /// vertices will be used further, so no need to transfer ownership, so passed by reference
        if (temp.intersect(ray) && temp.t < tNear)
        {
            tNear = temp.t;
            hitTriangle = temp;

            intersect |= true;                                                                      /// bit XOR 0^1=1, 1^0=1, 0^0=0, 1^1=0
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
    hitTriangle.getHitPointProperties();
    N = hitTriangle.N;
    tex = hitTriangle.tex;
}

/// print mesh
std::string Mesh::show()
{
    std::stringstream ss;
    ss << "Mesh " << id << '\n';
    ss << "number of triangles: " << trisNum << '\n';
    ss << "number of vertices: " << vertsNum << '\n';
    ss << "number of points: " << pointsNum << '\n';
    ss << "culling: " << (objectOpts.culling ? "on" : "off") << '\n';
    ss << "smooth shading: " << (objectOpts.smoothShading ? "on" : "off") << '\n';
    ss << '\n';

    return std::string(ss.str());
}




/// class Sphere
uint32_t Sphere::counter{0};

/// constructor with parameters
Sphere::Sphere(const Vec3_f_t& c, const float r, const MaterialProp &matProp, const ObjectOptions &objOpts):
        center{c}, radius{r}, radius2{r * r}
{
    id = ++counter;
    objType = ObjectType::SPHERE;
    materialProp = matProp;
    objectOpts = objOpts;

}

/// transform object with the object-to-world matrix
void Sphere::transformObjToWorld (const Matrix44_f_t &o2wMtx)                                       /// only translation!!! rotation does not make sense, scaling - via radius!
{
    objectToWorldMtx = o2wMtx;
    /// transform center
    Vec3_f_t tempCenter = objectToWorldMtx.multPointMatrix(center);
    center = tempCenter;
}

/// find a hit point
bool Sphere::intersect (const Ray& ray)
{
    /// analytic solution: |O + tD|^2 - R^2 = 0 - vektorno-parametricheskoe equation, from P^2 - R^2 = 0
    float t0, t1;

    Vec3_f_t O = ray.origin - center;
    float a = (ray.direction).dotProduct(ray.direction);
    float b = 2 * (ray.direction).dotProduct(O);
    float c = O.dotProduct(O) - radius2;

    if (!solveQuadratic(a, b, c, t0, t1))
        return false;
    if (t0 > t1)
        std::swap(t0, t1);
    /// if t0 is negative, then use t1 instead
    if (t0 < 0)
        t0 = t1;
    /// both t0 and t1 are negative
    if (t0 < 0)
        return false;

    t = t0;

    P = ray.origin + ray.direction * t;

    return true;
}

/// get hit point properties: N, tex
void Sphere::getHitPointProperties()
{
    n = (P - center);
    N = n.normalised();

    tex.m_x = atan2(N.m_z, N.m_x);                                                                  /// atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
    tex.m_y = acosf(N.m_y);                                                                         /// acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
    tex.m_x = (1 + tex.m_x / Consts::PI) * 0.5;
    tex.m_y = tex.m_y / Consts::PI;
}

/// convert sphere to mesh                                                                          /// optimize? ochen dolgo!
std::unique_ptr<Mesh> Sphere::tessellate(uint32_t divs)
{
    ///generate points, norms and texs
    uint32_t pointsNum = (divs - 1) * divs + 2;                                                     /// 2 sloya (polnye secheniya) * 3 dolki (polsecheniya) + 2 vershiny

    std::unique_ptr<Vec3_f_t[]> points{new Vec3_f_t[pointsNum]};
    std::unique_ptr<Vec3_f_t[]> norms {new Vec3_f_t[pointsNum]};
    std::unique_ptr<Vec2_f_t[]> txs {new Vec2_f_t[pointsNum]};

    float th = -Consts::PI;
    float ph = -Consts::PI_2;
    float dTh = Consts::PI / divs;
    float dPh = Consts::PI_2 / divs;

    points[0] = norms[0] = Vec3_f_t(0, -radius, 0) + center;                                        /// lower vertex
    uint32_t p = 1;
    for (uint32_t i = 0; i < divs - 1; i++)                                                         /// po gorizontali = shirota
    {
        th += dTh;
        ph = -Consts::PI;
        for (uint32_t j = 0; j < divs; j++)                                                         /// po vertikali = dolgota
        {
            float x = radius * sin(th) * sin(ph) + center.m_x;
            float y = radius * cos(th) + center.m_y;
            float z = radius * sin(th) * cos(ph) + center.m_z;

            points[p] = norms[p] = Vec3_f_t(x, y, z);                                               /// vertex normals!
            txs[p].m_x = th / Consts::PI + 0.5;                                                     /// perevod [0; 2Pi] v interval [0; 1]
            txs[p].m_y = ph * 0.5 / Consts::PI + 0.5;                                               /// perevod [-Pi; Pi] v interval [0; 1]

            ph += dPh;
            p++;
        }
    }
    points[p] = norms[p] = Vec3_f_t(0, radius, 0) + center;

    /// create verts in faces and verts indices arrays
    uint32_t facesNum = divs * divs;
    std::unique_ptr<uint32_t[]> vertsInFaces{new uint32_t[facesNum]};

    uint32_t vertsNum = (6 + (divs - 2) * 4) * divs;                                                /// 2 fans sverhu i snizu + quads, vse eto umnojit na kolichestvo vertikalnh polos
    std::unique_ptr<uint32_t[]> vertsIndices{new uint32_t[vertsNum]};                               /// obshee kolichestvo vershin

    uint32_t vid = 1, vNum = 0, l = 0;
    p = 0;
    for (uint32_t i = 0; i < divs; i++)                                                             /// kolichestvo gorizontalnyh polosok
    {
        for (uint32_t j = 0; j < divs; j++)                                                         /// kolichestvo polinomov v poloske
        {
        if (i == 0)                                                                                 /// esli verhnyaa poloska
        {
            vertsInFaces[p++] = 3;
            vertsIndices[l] = 0;
            vertsIndices[l + 1] = (j == (divs - 1)) ? vid : j + vid + 1;
            vertsIndices[l + 2] = j + vid;
            l += 3;
        }
        else if (i == divs - 1)                                                                     /// esli nijnyaya poloska
        {
            vertsInFaces[p++] = 3;
            vertsIndices[l] = j + vid + 1 - divs;
            vertsIndices[l + 1] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs;
            vertsIndices[l + 2] = vid + 1;
            l += 3;
        }
        else                                                                                        /// esli centralnye poloski
        {
            vertsInFaces[p++] = 4;
            vertsIndices[l] = j + vid + 1 - divs;
            vertsIndices[l + 1] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs;
            vertsIndices[l + 2] = (j == (divs - 1)) ? vid + 1 : j + vid + 2;
            vertsIndices[l + 3] = j + vid + 1;
            l += 4;
        }
        vNum++;
        }
    vid = vNum;
    }
    /// rewrite arrays, now normals and tex info for each vert, not for each point
    std::unique_ptr<Vec3_f_t[]> normals {new Vec3_f_t[vertsNum]};
    std::unique_ptr<Vec2_f_t[]> texs {new Vec2_f_t[vertsNum]};

    for (uint32_t i = 0; i < vertsNum; ++i)
    {
        normals[i] = norms[vertsIndices[i]];
        texs[i] = txs[vertsIndices[i]];
    }

    return std::unique_ptr<Mesh> {new Mesh(facesNum, std::move(vertsInFaces), std::move(points), std::move(vertsIndices),
                    std::move(normals), std::move(texs), materialProp, objectOpts)};
}

/// print sphere
std::string Sphere::show()
{
    std::stringstream ss;
    ss << "Sphere " << id << '\n';
    ss << "center: " << center << '\n';
    ss << "radius: " << radius << '\n';
    ss << "culling: " << (objectOpts.culling ? "on" : "off") << '\n';
    ss << "smooth shading: " << (objectOpts.smoothShading ? "on" : "off") << '\n';    ss << '\n';

    return std::string(ss.str());
}
