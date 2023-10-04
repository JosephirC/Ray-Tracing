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

struct Cube{
    vec3 a;
    vec3 b;
    int i;
};

struct Plane{
    vec3 n;// Normal
    vec3 p;// Point
    int i;// Texture Id
};

struct disc{
    vec3 n;
    vec3 p;
    float r;
    int i ;
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

//fonction résolution equation du second degré:
float solvRoots(float a, float b, float c){
    float t;
    //on calcul le déterminant d pour trouver les racine de l'équation
    float d = b*b - 4.*a*c;
    if(d>0.)
    {
        float t1 = (-b - sqrt(d)) / (2.0*a);
        float t2 = (-b + sqrt(d)) / (2.0*a);
        float t = min(t1, t2);
        return t;
    }
}

// Plane intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectPlane(Ray ray,Plane pl,out Hit x){
    float t=-dot(ray.o-pl.p,pl.n)/dot(ray.d,pl.n);
    if(t>0.)
    {
        
        x=Hit(t,vec3(0,0,1),0);
        return true;
    }
    return false;
}

bool IntersectCube(Ray ray, Cube cb, out Hit x){
    vec3 t0, t1;
    //défini tout les cas de t0 et t1 pour x puis y puis z dans tous les angles de la cam
    if (ray.d.x >= 0.) { 
        t0.x = (cb.a.x - ray.o.x) / ray.d.x;
        t1.x = (cb.b.x - ray.o.x) / ray.d.x;
    } 
    else { 
        t0.x = (cb.b.x - ray.o.x) / ray.d.x;
        t1.x = (cb.a.x - ray.o.x) / ray.d.x;
    } 
    
    if (ray.d.y >= 0.) { 
        t0.y = (cb.a.y - ray.o.y) / ray.d.y;
        t1.y = (cb.b.y - ray.o.y) / ray.d.y;
    } 
    else { 
        t0.y = (cb.b.y - ray.o.y) / ray.d.y;
        t1.y = (cb.a.y - ray.o.y) / ray.d.y;
    } 
    
    if (ray.d.z >= 0.) { 
        t0.z = (cb.a.z - ray.o.z) / ray.d.z;
        t1.z = (cb.b.z - ray.o.z) / ray.d.z;
    } 
    else { 
        t0.z = (cb.b.z - ray.o.z) / ray.d.z;
        t1.z = (cb.a.z - ray.o.z) / ray.d.z;
    }
    //filtre des cas ne touchant pas la boite
    if (t0.x > t1.y || t0.y > t1.x) return false;
    //filtre pour trouver tmin et tmax
    float tmin = max(t0.x, t0.y);
    float tmax = min(t1.x, t1.y);
    //même opération étendu à z
    if (tmin > t1.z || t0.z > tmax) return false;
    tmin = max(tmin, t0.z);
    tmax = min(tmax, t1.z);
    //sélection du t:
    float t = (tmin < 0.) ? tmax : tmin;
    if(t>0.){
        vec3 p=Point(ray,t);
        //le centre du cube
        vec3 c;
        if (length(cb.a) < length(cb.b))
            c = cb.a + cb.b/2.;
        else
            c = cb.b + cb.a/2.;
        x=Hit(t, normalize(p-c),cb.i);
        return true;
    }
    return false;
    
}

// disc intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectDisc(Ray ray,disc disc,out Hit x){
    bool pl = IntersectPlane(ray, Plane(disc.n, disc.p, disc.i), x);
    if(pl){
        vec3 p=Point(ray, x.t);
        if(length(p-disc.p) < disc.r){
            x = Hit(x.t, disc.n, disc.i);
            return true;
        }
    }
    return false;
}

// Sphere intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectSphere(Ray ray,Sphere sph,out Hit x){
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

bool IntersectEllipsoide(Ray ray,Ellipsoide ellip,out Hit x){
    vec3 oc=ray.o-ellip.c;
    float a = dot((ray.d/ellip.r),(ray.d/ellip.r));
    float b = 2.0*dot((oc/ellip.r), (ray.d/ellip.r));
    float c = dot((oc/ellip.r), (oc/ellip.r)) - dot(ellip.r, ellip.r);

    float t = solvRoots(a, b, c);
    if(t>0.)
    {
        vec3 p=Point(ray,t);
        //vec3 normal = normalize( (p-ellip.c) / dot(ellip.r, ellip.r));
          
        //x=Hit(t,normal,ellip.i);
        x=Hit(t,normalize(p-ellip.c),ellip.i);
        return true;
    }
    return false;
}

// Cylinder intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectCylinderBase(Ray ray,Cylinder cyl,out Hit x)
{
    vec3 oa = ray.o-cyl.a;
    vec3 u = normalize(cyl.b - cyl.a);
    
    float a = length(ray.d) - dot(ray.d, u)*dot(ray.d, u); 
    float b = 2. * ( dot(oa, ray.d) - dot(oa, u) * dot(ray.d, u));
    float c = dot(oa, oa) - dot(oa, u)*dot(oa, u) - cyl.r*cyl.r;
    
    float t = solvRoots(a, b, c);
    if(t>0.){
        vec3 p=Point(ray,t);
        float v = dot(p-cyl.a, u);
        vec3 h = cyl.a + v*u;
        if (v>=0. && v<= length(cyl.b-cyl.a)){
            x=Hit(t,normalize(p-h),cyl.i);
            return true;
        }
    }
    return false;
}

bool IntersectCapsule(Ray ray, Cylinder cyl, out Hit x){
    x=Hit(1000.,vec3(0),-1);
    Hit x_;
    Sphere sph1 = Sphere(cyl.a, cyl.r, cyl.i);
    Sphere sph2 = Sphere(cyl.b, cyl.r, cyl.i);
    bool b1 = IntersectSphere(ray, sph1, x_);
    if(b1 && x_.t < x.t){
        x=x_;
    }
    bool b2 = IntersectSphere(ray, sph2, x_);
    if(b2 && x_.t < x.t){
        x=x_;
    }
    bool corp = IntersectCylinderBase(ray, cyl, x_);
    if(corp && x_.t < x.t){
        x=x_;
    }
    return b1 || corp || b2;  
}

bool IntersectCylinder(Ray ray, Cylinder cyl, out Hit x){
    x=Hit(1000.,vec3(0),-1);
    Hit x_;
    
    disc ds1 = disc(normalize(cyl.b-cyl.a), cyl.a, cyl.r, cyl.i);
    disc ds2 = disc(normalize(cyl.b-cyl.a), cyl.b, cyl.r, cyl.i);
    Sphere sph2 = Sphere(cyl.b, cyl.r, cyl.i);
    bool b1 = IntersectDisc(ray, ds1, x_);
    if(b1 && x_.t < x.t){
        x.t=x_.t;
    }
    bool b2 = IntersectDisc(ray, ds2, x_);
    if(b2 && x_.t < x.t){
        x.t=x_.t;
    }
    bool corp = IntersectCylinderBase(ray, cyl, x_);
    if(corp && x_.t < x.t){
        x=x_;
    }
    return b1 || corp || b2;  
}

// OPERATORS

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
    
    const Cylinder cyll1 = Cylinder(vec3(3.,2.,2.), vec3(2., 4., 5.), 1. ,1);
    
    const disc ds = disc(normalize(vec3(0.,0.,1.)), vec3(3.,3.,1.), 4.,1);
    
    const Cube cb = Cube(vec3(-5., -5., 0.), vec3(-2., 0., 4.), 1);
    
    x=Hit(1000.,vec3(0),-1);
    Hit current;
    bool ret=false;
    if(IntersectSphere(ray,sph1,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
    
    if(IntersectCube(ray,cb,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
    
    if(IntersectDisc(ray,ds,current)&&current.t<x.t){
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
    if(IntersectCylinder(ray,cyll1,current)&&current.t<x.t){
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
