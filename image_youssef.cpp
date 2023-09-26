struct Sphere{
    vec3 c;// Center
    float r;// Radius
    int i;// Texture Id
};

struct Ellipsoide{
    vec3 c; // Center
    vec3 r; // Radius
    int i; // Texture Id
};

struct Cylinder{
    vec3 a;// Bottom Center
    vec3 b;// Middle Center
    float r;// Radius
    int i;// Texture Id
};

struct Plane{
    vec3 n;// Normal
    vec3 p;// Point
    int i;// Texture Id
};

struct Hit{
    float t;// Intersection depth
    vec3 n;// Normal
    int i;// Texture Id
};

struct Ray{
    vec3 o;// Origin
    vec3 d;// Direction
};

struct Material
{
    vec3 d;// Diffuse
};

float Checkers(in vec2 p)
{
    // Filter kernel
    vec2 w=fwidth(p)+.001;
    // Box box filter
    vec2 i=2.*(abs(fract((p-.5*w)*.5)-.5)-abs(fract((p+.5*w)*.5)-.5))/w;
    // xor pattern
    return.5-.5*i.x*i.y;
}

// Compute point on ray
vec3 Point(Ray ray,float t)
{
    return ray.o+t*ray.d;
}

// Compute color
// i : Texture index
// p : Point
Material Texture(vec3 p,int i)
{
    if(i==1)
    {
        return Material(vec3(.8,.5,.4));
    }
    else if(i==0)
    {
        // compute checkboard
        float f=Checkers(.5*p.xy);
        vec3 col=vec3(.4,.5,.7)+f*vec3(.1);
        return Material(col);
    }
    return Material(vec3(0));
}

// Sphere intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectSphere(Ray ray,Sphere sph,out Hit x)
{
    vec3 oc=ray.o-sph.c;
    float b=dot(oc,ray.d);
    float c=dot(oc,oc)-sph.r*sph.r;
    float d=b*b-c;
    if(d>0.)
    {
        float t=-b-sqrt(d);
        if(t>0.)
        {
            vec3 p=Point(ray,t);
            x=Hit(t,normalize(p-sph.c),sph.i);
            
            return true;
        }
    }
    return false;
    
}

// Ellipsoide intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectEllipsoide(Ray ray,Ellipsoide ellip,out Hit x)
{
    vec3 oc=ray.o-ellip.c;
    float a = dot((ray.d/ellip.r),(ray.d/ellip.r));
    float b = 2.0*dot((oc/ellip.r), (ray.d/ellip.r));
    float c = dot((oc/ellip.r), (oc/ellip.r)) - dot(ellip.r, ellip.r);

    float d = b*b - 4. *a *c;
    
    if(d>0.)
    {
        float t1 = (-b - sqrt(d)) / (2.0*a);
        float t2 = (-b + sqrt(d)) / (2.0*a);
        float t = min(t1, t2);
        if(t>0.)
        {
            vec3 p=Point(ray,t);
            //vec3 normal = normalize( (p-ellip.c) / dot(ellip.r, ellip.r));
            
            //x=Hit(t,normal,ellip.i);
            x=Hit(t,normalize(p-ellip.c),ellip.i);
            return true;
        }
    }
    return false;
}

// Cylinder intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectCylinder(Ray ray,Cylinder cyl,out Hit x)
{
    vec3 oc=ray.o-cyl.c;
    
    float b=dot(oc,ray.d);
    float c=dot(oc,oc)-cyl.r*cyl.r;
    
    float d=b*b-c;
    
    if(d>0.)
    {
        float t=-b-sqrt(d);
        if(t>0.)
        {
            vec3 p=Point(ray,t);
            x=Hit(t,normalize(p-cyl.c),cyl.i);
            
            return true;
        }
    }
    return false;
    
}

// Plane intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectPlane(Ray ray,Plane pl,out Hit x)
{
    float t=-dot(ray.o-pl.p,pl.n)/dot(ray.d,pl.n);
    if(t>0.)
    {
        
        x=Hit(t,vec3(0,0,1),0);
        return true;
    }
    return false;
}

// Scene intersection
// ray : The ray
//   x : Returned intersection information
bool Intersect(Ray ray,out Hit x)
{
    // Spheres
    const Sphere sph1=Sphere(vec3(0.,0.,1.),1.,1);
    const Sphere sph2=Sphere(vec3(2.,0.,2.),1.,1);
    const Plane pl=Plane(vec3(0.,0.,1.),vec3(0.,0.,0.),0);
    
    const Ellipsoide ellip1 = Ellipsoide(vec3(-4., 3., 2.), vec3(1.6,1.,0.5), 1);
    
    x=Hit(1000.,vec3(0),-1);
    Hit current;
    bool ret=false;
    if(IntersectSphere(ray,sph1,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
    
    if(IntersectSphere(ray,sph2,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
    if(IntersectPlane(ray,pl,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
    if(IntersectEllipsoide(ray,ellip1,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
    
    return ret;
}

vec3 Background(vec3 rd)
{
    return mix(vec3(.8,.8,.9),vec3(.7,.7,.8),rd.z);
}

// Camera rotation matrix
// ro : Camera origin
// ta : Target point
mat3 setCamera(in vec3 ro,in vec3 ta)
{
    vec3 cw=normalize(ta-ro);
    vec3 cp=vec3(0,0,1);
    vec3 cu=-normalize(cross(cw,cp));
    vec3 cv=-normalize(cross(cu,cw));
    return mat3(cu,cv,cw);
}

// Apply color model
// m : Material
// n : normal
vec3 Color(Material m,vec3 n)
{
    vec3 light=normalize(vec3(1,1,1));
    
    float diff=clamp(dot(n,light),0.,1.);
    vec3 col=m.d*diff+vec3(.2,.2,.2);
    return col;
}

// Rendering
vec3 Shade(Ray ray)
{
    // Intersect contains all the geo detection
    Hit x;
    bool idx=Intersect(ray,x);
    
    if(idx)
    {
        vec3 p=Point(ray,x.t);
        Material mat=Texture(p,x.i);
        
        return Color(mat,x.n);
    }
    else
    {
        return Background(ray.d);
    }
    
    return vec3(0);
}

void mainImage(out vec4 fragColor,in vec2 fragCoord)
{
    // From uv which are the pixel coordinates in [0,1], change to [-1,1] and apply aspect ratio
    vec2 uv=(-iResolution.xy+2.*fragCoord.xy)/iResolution.y;
    
    // Mouse control
    vec2 mouse=iMouse.xy/iResolution.xy;
    
    // Ray origin
    vec3 ro=12.*normalize(vec3(sin(2.*3.14*mouse.x),cos(2.*3.14*mouse.x),1.4*(mouse.y-.1)));
    vec3 ta=vec3(0.,0.,1.5);
    mat3 ca=setCamera(ro,ta);
    
    // Ray
    vec3 rd=ca*normalize(vec3(uv.xy*tan(radians(22.5)),1.));
    
    // Render
    vec3 col=Shade(Ray(ro,rd));
    
    fragColor=vec4(col,1.);
}
