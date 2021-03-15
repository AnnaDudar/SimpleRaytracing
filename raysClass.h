#ifndef RAYTRACING_D_Programming_Scratchapixel_Raytracing_RAYS_H_INCLUDED
#define RAYTRACING_D_Programming_Scratchapixel_Raytracing_RAYS_H_INCLUDED




class Ray
{
public:
    Vec3_f_t origin {} ;
    Vec3_f_t direction {Vec3_f_t(1).normalised()};
public:
    Ray() = default;
    Ray(Vec3_f_t o, Vec3_f_t d): origin{o}, direction{d}
    {}
    void print() const
    {
        std::cout << origin << '\n';
        std::cout << direction << '\n';
        std::cout << '\n';
    }
};




#endif // RAYTRACING_D_Programming_Scratchapixel_Raytracing_RAYS_H_INCLUDED
