// IMAGE - M1 INFO
/* VERY IMPORTANT !!

    - TO SEE THE AMBIENT OCCLUSION PLEASE COMMENT THE LINE `x = Hit(1000., vec3(0), -1);` in Shade()
    AND DECOMMENT IT IN INTERSECT().

    - IN COLOR() YOU CAN PLAY AROUND WITH THE NUMBER OF SCENE LIGHTS.

    - TO MAKE GOURSAT STOP GROWING REPLACE `+ 3.5 * cos(iTime) + 14.5` WITH  ` + 11.8`

    - TO MAKE OCTAEDER AND PIPECONNECTOR STOP GROWING REPLACE `+ cos(iTime) * 67.5 - 82.5` WITH `-25.`
*/ 

#define MAX_REFLECTION 5

// Definition of primitives
struct Sphere {
    vec3 center;// Center
    float radius;// Radius
    int id;// Texture Id
};

struct Ellipsoid {
    vec3 center;// Center
    vec3 radius;// Radius
    int id;// Texture Id
};

struct Cylinder {
    vec3 bottomCenter;// Bottom Center
    vec3 topCenter;// Top Center
    float radius;// Radius
    int id;// Texture Id
};

struct Disc {
    vec3 normal;// Normal
    vec3 point;// Center
    float radius;// Radius
    int id;// Texture Id
};

struct Capsule {
    vec3 bottomCenter;// Bottom Center
    vec3 topCenter;// Top Center
    float radius;// Radius
    int id;// Texture Id
};

struct Box {
    vec3 bottomCorner;// First Corner
    vec3 topCorner;// Second Corner
    int id;// Texture Id
};

struct Torus {
    vec3 center;// Center
    float bigRadius;// Radius between center and the donut cirlce
    float smallRadius;// Radius between the donut circle and the outer circle
    int id;// Texture Id
};

struct Goursat {
    vec3 center;// Center
    int id;// Texture Id
};

struct Octaeder {
    vec3 center;// Center
    int id;// Texture Id
};

struct PipeConnector {
    vec3 center;// Center
    int id;// Texture Id
};

struct Plane {
    vec3 normal;// Normal
    vec3 point;// Point
    int id;// Texture Id
};

struct Hit {
    float t;// Intersection depth
    vec3 normal;// Normal
    int id;// Texture Id
};

struct Ray {
    vec3 origin;// Origin
    vec3 direction;// Direction
    bool isHomo;// Homothetie boolean
    bool isRot;// Rotation boolean
};

struct Light {
    vec3 lightColor;//color
    vec3 lightPos;//position
};

struct Scene {
    //for transformation ray( translation, homo, rotation )
    Ray tabRay[10];
    vec3 tabScale[10];
    vec3 tabAngle[10];

    int nbSphere;// Number of displayed Spheres in the Scene
    Sphere tabSphere[10];// Array of Spheres
    Plane plane;// Plane of the Scene

    int nbEllipsoid;// Number of displayed Ellipsoide in the Scene
    Ellipsoid tabEllipsoid[10];// Array of Ellipsoides

    int nbCylinder;// Number of displayed Cylinder in the Scene
    Cylinder tabCylinder[10];// Array of Cylinders

    int nbCapsule;// Number of displayed Capsule in the Scene
    Capsule tabCapsule[10];// Array of Capsules

    int nbBox;// Number of displayed Box in the Scene
    Box tabBox[10];// Array of Boxes

    int nbTorus;// Number of displayed Torus in the Scene
    Torus tabTorus[10];// Array of Torus

    int nbGoursat;// Number of displayed EL Goursat in the Scene
    Goursat tabGoursat[10];// Array of EL Goursat

    int nbOctaeder;// Number of displayed Octaeder in the Scene
    Octaeder tabOctaeder[10];// Array of Octaeders

    int nbPipeConnector;// Number of displayed PipeConnector in the Scene 
    PipeConnector tabPipeConnector[10];// Array of PipeConnector

    int nbLight;// Number of displayed Light sources in the Scene
    Light tabLight[25];// Array of Light sources
};

struct Material {
    vec3 diffuse;// Diffuse
    vec3 ambient;//ambiant
    vec3 specular;//specular
    float coefShininess; // reflection coef
    float reflexivity; // reflexivity du materiel
    vec3 mirrorColor; // couleur du miroir
};

//////////////////////////////////////////////////////////////////////////
//Texturing

/**
 * @brief fonctions using for checkerboard Texture
 * 
 * @param vec2 p 
 * @return float 
 */
float Checkers(in vec2 p) {
    // Filter kernel
    vec2 w = fwidth(p)+.001;
    // Box box filter
    vec2 i = 2. * (abs(fract((p-.5*w)*.5)-.5)-abs(fract((p+.5*w)*.5)-.5))/w;
    // xor pattern
    return 0.5 - 0.5 * i.x * i.y;
}

/**
 * @brief Compute point on ray
 * 
 * @param Ray ray 
 * @param float t 
 * @return vec3 
 */
vec3 Point(Ray ray,float t) {
    return ray.origin + t * ray.direction;
}

/**
 * @brief Hashing function
 * 
 * @param p : Vector in space
 * @return float : a random number in [-1,1]
 */
float Hash(in vec3 p)  
{
    p  = fract( p*0.3199+0.152 );
	p *= 17.0;
    return fract( p.x*p.y*p.z*(p.x+p.y+p.z) );
}

/**
 * @brief Procedural value noise with cubic interpolation
 * 
 * @param p : Point 
 * @return float : a random coefiscient
 */
float Noise(in vec3 p) {
    vec3 i = floor(p);
    vec3 f = fract(p);
  
    f = f*f*(3.0-2.0*f);
    // Could use quintic interpolation instead of cubic
    // f = f*f*f*(f*(f*6.0-15.0)+10.0);

    return mix(mix(mix( Hash(i+vec3(0,0,0)), 
                        Hash(i+vec3(1,0,0)),f.x),
                   mix( Hash(i+vec3(0,1,0)), 
                        Hash(i+vec3(1,1,0)),f.x),f.y),
               mix(mix( Hash(i+vec3(0,0,1)), 
                        Hash(i+vec3(1,0,1)),f.x),
                   mix( Hash(i+vec3(0,1,1)), 
                        Hash(i+vec3(1,1,1)),f.x),f.y),f.z);
}

/**
 * @brief Sum of Noize
 * 
 * @param p : Point
 * @param waveLength : Controls the size of noise patterns
 * @param coef : Coefficient for noise calculation
 * @param iterations : Number of iterations for noise summation
 * @return float : Resulting turbulence value
 */
float Turbulence(in vec3 p, float waveLength, float coef, int iterations) {// somme de bruits 
    float turbulence = coef * Noise(p / waveLength);
    for (int i = 0; i < iterations; i++) { //boucle pour calculer la somme de bruit
        coef = coef * 0.5;
        waveLength = waveLength * 0.5;
        turbulence = turbulence + coef * Noise(p / waveLength);
    }
    return turbulence;
}

/**
 * @brief convert a float into the type int
 * 
 * @param a : Float
 * @return int 
 */
int convert(float a){
    if (a < 0.0) {
        a= a - 1.0;
    }
    return int(a);
}

/**
 * @brief Texture of a checkerboard
 * 
 * @param point : Point in the space
 * @param color1 : First color of the checkerboard
 * @param color2 : Second color of the checkerboard
 * @return vec3 : Resulting texture
 */
vec3 Checkboard(vec3 point, vec3 color1, vec3 color2 ){
    vec3 w = fwidth(point) +.001;
    int x = convert(point.x);    //convertion en entier
    int y = convert(point.y);    //correction partie entière des valeurs négatifs pour evité le cas de int(0.4)=int(-0.4)
    int z = convert(point.z);

        if ((x + y + z) % 2 == 0) {//modulo 2 pour faire 2 cas: coordonées paire, coordonées impaire
        return color1 / w;
    } else {
        return color2 / w;
    }   
}

/**
 * @brief Texture of a marble
 * 
 * @param point : Point in the space
 * @param color1 : First material of the marble
 * @param color2 : Second material of the marble
 * @return Material : Resulting Marble
 */
Material MarbleTexture(vec3 point, Material color1, Material color2) {
    point = point + Turbulence(point, 1., 20., 10);
    float t = cos(point.x);
    Material finalColor = Material(color1.diffuse * t + color2.diffuse * (0.6 - t), vec3(0.2), color1.specular * t + color2.specular * (1.0 - t), 50., 0., vec3(0.));
    return finalColor;
}

/**
 * @brief Texture of a Veins marble
 * 
 * @param point : Point in the space
 * @param color1 : Vein material of the marble
 * @param color2 : back material of the marble
 * @return Material : Resulting Veins Marble
 */
Material VeinMarbleTexture(vec3 point, Material color1, Material color2) {
    point = point + Turbulence(point, 5., 5., 10);
    float t = tan(point.x*3.);
    if(t < 4.) {
        return color1;
    } 
    else {
        return color2;
    }
}

// Source : https://www.shadertoy.com/view/wl3czM
float n21(vec2 p) {
	const vec3 s = vec3(7, 157, 0);
	vec2 h,
	     ip = floor(p);
	p = fract(p);
	p = p * p * (3. - 2. * p);
	h = s.zy + dot(ip, s.xy);
	h = mix(fract(sin(h) * 43.5453), fract(sin(h + s.x) * 43.5453), p.x);
	return mix(h.x, h.y, p.y);
}

// Source : https://www.shadertoy.com/view/wl3czM
float n11(float p) {
	float ip = floor(p);
	p = fract(p);
	vec2 h = fract(sin(vec2(ip, ip + 1.) * 12.3456) * 43.5453);
	return mix(h.x, h.y, p * p * (3. - 2. * p));
}

/**
 * @brief Texture of a Wood 
 * 
 * @param point : Point in the space
 * @return Material : Resulting a Wood Texture
 */
float WoodTexture(vec2 p) {
	p.x *= 71.;
	p.y *= 1.9;
	return n11(n21(p) * 30.);
}

/**
 * @brief Texture of a Wood 
 * 
 * @param point : Point in the space
 * @return Material : Resulting a Wood Texture
 */
vec3 WoodTexture2(vec3 p, vec3 color1, vec3 color2) {
    p = p + Turbulence(p, 5., 0.5, 5);
    float r = sqrt (p.x * p.x + p.y * p.y);
    float pattern = 0.5 + 0.5 * cos (10. * 3.1415927 * r) ;   
    return mix(color1, color2, pattern); 
}

/**
 * @brief Definition of all texture in the scene
 * 
 * @param p : Point in the space
 * @param i : Texture index
 * @return Material : Specific material per index
 */
Material Texture(vec3 p, int i) {
    if (i == 1) { // uniform
         // return Material(vec3(.8,.0,.1),vec3(0.2,0.2,0.2), vec3(0.2, 0.2, 0.2), 50., 0., vec3(0.));
        return Material(vec3(.8,.5,.4), vec3(0.1), vec3(0.7, 0.7, 0.7), 50., 0., vec3(0.));
    }
    else if (i == 2) { // variation
        vec3 colorA = vec3(0.1,0.1,0.9);
        vec3 colorB = vec3(1.,0.8,0.2);
        return Material(mix(colorA, colorB, sin(iTime*0.5)),vec3(0.4,0.4,0.4), vec3(0.2, 0.2, 0.2), 50., 0., vec3(0.) );
    }
    else if (i == 3) { // Checkboard
        vec3 texture = Checkboard(p, vec3(1., 1., 1.), vec3(0., 0., 0.));
        return Material(texture, vec3(0.4,0.4,0.4), vec3(0.2, 0.2, 0.2), 50., 0., vec3(0.) );
    }
    else if (i == 4) { // Marble
        Material color1 = Material(vec3(1.2, 1.2, 1.2), vec3(0.4,0.4,0.4), vec3(0.9, 0.9, 0.9), 50., 0., vec3(0.));
        Material color2 = Material(vec3(0.83, .83, .83), vec3(0.4,0.4,0.4), vec3(0.2, 0.2, 0.2), 50., 0., vec3(0.));
        Material marble = MarbleTexture(p, color1, color2);
        return marble;
    }
    else if (i == 5) { // Marble
        Material color1 = Material(vec3(0.3, 0.3, 0.3), vec3(0.2,0.2,0.2), vec3(0, 0, 0), 1., 0., vec3(0.));
        Material color2 = Material(vec3(0.9, 0.9, 0.2), vec3(0.7,0.7,0.7), vec3(3, 3, 3), 100., 0., vec3(0.));
        Material marble = VeinMarbleTexture(p, color1, color2);
        return marble;
    }
    else if (i == 6) { // Mirror
        return Material(vec3(.8,.5,.4), vec3(.7), vec3(0.2), 50., 1., vec3(0,0, 0.));
    }
    else if (i == 7) { // Wood
        float woodValue = WoodTexture(p.xy);
        vec3 color1 = vec3(0.4, 0.2, 0.1);
        vec3 color2 = vec3(0.2, 0.1, 0.07);
        vec3 texture = mix(color1, color2, woodValue);
        return Material(texture, vec3(.2), vec3(0.1), 50., 0., vec3(0, 0, 0));
    }
    else if (i == 8) { // Wood 2
        vec3 color1 = vec3(.17, .1, .05);
        vec3 color2 = vec3(.08, .05, .03);
        vec3 texture = WoodTexture2(p, color1, color2);
        return Material(texture, vec3(.7), vec3(0.2), 50., 0., vec3(0, 0, 0));
    }
    else if (i == 0) { // classic checkboard
        float f = Checkers(.5*p.xy);
        vec3 col = vec3(.4,.5,.7) + f * vec3(.1);
        return Material(col, vec3(0.2, 0.2, 0.2), vec3(0.9, 0.9, 0.9), 50., 0., vec3(0.));
    }
    return Material(vec3(0),vec3(0.,0.,0.), vec3(0.2, 0.2, 0.2), 50., 0., vec3(0.));
}

/**
 * @brief resolv second order equation : ax² + bx + c = 0
 * 
 * @param a : float
 * @param b : float
 * @param c : float
 * @return float : min of the 2 solution
 */
float solvRoots(float a, float b, float c) {
    float t;
    //on calcul le déterminant d pour trouver les racine de l'équation
    float d = b * b - 4. * a * c;
    if (d > 0.){
        float t1 = (-b - sqrt(d)) / (2.0 * a);
        float t2 = (-b + sqrt(d)) / (2.0 * a);
        float t = min(t1, t2);
        return t;
    }
}

//////////////////////////////////////////////////////////////////////
//Intersection functions

/**
 * @brief Classic sphere intersection
 * 
 * @param ray : The ray
 * @param sph : Structure information
 * @param x : Returned intersection information
 * @return true
 * @return false
 */
bool IntersectSphere(Ray ray, Sphere sph, out Hit x) {
    vec3 oc = ray.origin - sph.center;
    float b = dot(oc, ray.direction);
    float c = dot(oc, oc) - sph.radius * sph.radius;
    float d = b * b - c;
    if (d>0.) {
        float t = -b - sqrt(d);
        if (t>0.) {
            vec3 p = Point(ray,t);
            x = Hit(t, normalize(p-sph.center), sph.id); 
            return true;
        }
    }
    return false;
}

/**
 * @brief Classic Plane intersection
 * 
 * @param ray : The ray
 * @param pl : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectPlane(Ray ray, Plane pl, out Hit x) {
    pl.normal = normalize(pl.normal);
    float t = -dot(ray.origin - pl.point, pl.normal) / dot(ray.direction, pl.normal);
    if (t > 0.) {
        vec3 p = Point(ray, t);
        x = Hit(t, pl.normal, pl.id);
        return true;
    }
    return false;
}

/**
 * @brief intersect between a Ray and an Ellipsoid
 * 
 * @param ray : The ray
 * @param ellip : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectEllipsoid(Ray ray, Ellipsoid ellip, out Hit x) {
    vec3 oc = ray.origin - ellip.center;
    float a = dot((ray.direction / ellip.radius), (ray.direction / ellip.radius));
    float b = 2. * dot((oc / ellip.radius), (ray.direction / ellip.radius));
    float c = dot((oc/ellip.radius), (oc/ellip.radius)) - dot(ellip.radius, ellip.radius);

    float t = solvRoots(a, b, c);
    if (t>0.) {
        vec3 p = Point(ray,t);
        x = Hit(t,normalize(p-ellip.center),ellip.id);
        return true;
    }
    return false;
}

/**
 * @brief intersect between a Ray and a Cylinder body
 * 
 * @param ray : The ray
 * @param cyl : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectCylinderBody(Ray ray, Cylinder cyl, out Hit x) {
    vec3 oa = ray.origin - cyl.bottomCenter;
    vec3 u = normalize(cyl.topCenter - cyl.bottomCenter);
    
    float a = dot(ray.direction, ray.direction) - dot(ray.direction, u) * dot(ray.direction, u); 
    float b = 2. * ( dot(oa, ray.direction) - dot(oa, u) * dot(ray.direction, u));
    float c = dot(oa, oa) - dot(oa, u) * dot(oa, u) - cyl.radius * cyl.radius;
    float t = solvRoots(a, b, c);
    if (t>0.) {
        vec3 p = Point(ray, t);
        float v = dot(p - cyl.bottomCenter, u);
        vec3 h = cyl.bottomCenter + v * u;
        if (v >= 0. && v <= length(cyl.topCenter - cyl.bottomCenter)) {
            x = Hit(t, normalize(p - h), cyl.id);
            return true;
        }
    }
    return false;
}

/**
 * @brief intersect between a Ray and a Disc
 * 
 * @param ray : The ray
 * @param disc : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectDisc(Ray ray,Disc disc,out Hit x) {
    disc.normal = normalize(disc.normal);
    bool pl = IntersectPlane(ray, Plane(disc.normal, disc.point, disc.id), x);
    if (pl) {
        vec3 p = Point(ray, x.t);
        if (length(p - disc.point) < disc.radius) {
            // x = Hit(x.t,disc.n, disc.i);
            return true;
        }
    }
    return false;
}

/**
 * @brief intersect between a Ray and a Cylinder
 * 
 * @param ray : The ray
 * @param cyl : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectCylinder(Ray ray, Cylinder cyl, out Hit x) {
    x = Hit(1000., vec3(0), -1);
    Hit x_;
    
    Disc ds1 = Disc(-normalize(cyl.topCenter - cyl.bottomCenter), cyl.bottomCenter, cyl.radius, cyl.id);
    Disc ds2 = Disc(normalize(cyl.topCenter - cyl.bottomCenter), cyl.topCenter, cyl.radius, cyl.id);
    bool b1 = IntersectDisc(ray, ds1, x_);
    if (b1 && x_.t < x.t) {
        x = x_;
    }
    bool b2 = IntersectDisc(ray, ds2, x_);
    if (b2 && x_.t < x.t) {
        x = x_;
    }
    bool body = IntersectCylinderBody(ray, cyl, x_);
    if (body && x_.t < x.t) {
        x = x_;
    }
    return b1 || b2 || body;  
}

/**
 * @brief intersect between a Ray and a Capsule
 * 
 * @param ray : The ray
 * @param cap : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectCapsule(Ray ray, Capsule cap, out Hit x) {
    x = Hit(1000., vec3(0), -1);
    Hit x_;
    Sphere sph1 = Sphere(cap.bottomCenter, cap.radius, cap.id);
    Sphere sph2 = Sphere(cap.topCenter, cap.radius, cap.id);
    bool b1 = IntersectSphere(ray, sph1, x_);
    if (b1 && x_.t < x.t) {
        x = x_;
    }
    bool b2 = IntersectSphere(ray, sph2, x_);
    if (b2 && x_.t < x.t) {
        x = x_;
    }
    Cylinder cyl;
    cyl.bottomCenter = cap.bottomCenter;
    cyl.topCenter = cap.topCenter;
    cyl.radius = cap.radius;
    cyl.id = cap.id;
    bool body = IntersectCylinderBody(ray, cyl, x_);
    if (body && x_.t < x.t) {
        x = x_;
    }
    return b1 || body || b2;  
}

// Source : https://www.shadertoy.com/view/fsB3Wt
float cbrt(in float x) { return sign(x) * pow(abs(x), 1.0 / 3.0); }

int solveQuartic(in float a, in float b, in float c, in float d, in float e, inout vec4 roots) {
    b /= a; c /= a; d /= a; e /= a; // Divide by leading coefficient to make it 1

    // Depress the quartic to x^4 + px^2 + qx + r by substituting x-b/4a
    // This can be found by substituting x+u and the solving for the value
    // of u that makes the t^3 term go away
    float bb = b * b;
    float p = (8.0 * c - 3.0 * bb) / 8.0;
    float q = (8.0 * d - 4.0 * c * b + bb * b) / 8.0;
    float r = (256.0 * e - 64.0 * d * b + 16.0 * c * bb - 3.0 * bb * bb) / 256.0;
    int n = 0; // Root counter

    // Solve for a root to (t^2)^3 + 2p(t^2)^2 + (p^2 - 4r)(t^2) - q^2 which resolves the
    // system of equations relating the product of two quadratics to the depressed quartic
    float ra =  2.0 * p;
    float rb =  p * p - 4.0 * r;
    float rc = -q * q;

    // Depress using the method above
    float ru = ra / 3.0;
    float rp = rb - ra * ru;
    float rq = rc - (rb - 2.0 * ra * ra / 9.0) * ru;

    float lambda;
    float rh = 0.25 * rq * rq + rp * rp * rp / 27.0;
    if (rh > 0.0) { // Use Cardano's formula in the case of one real root
        rh = sqrt(rh);
        float ro = -0.5 * rq;
        lambda = cbrt(ro - rh) + cbrt(ro + rh) - ru;
    }

    else { // Use complex arithmetic in the case of three real roots
        float rm = sqrt(-rp / 3.0);
        lambda = -2.0 * rm * sin(asin(1.5 * rq / (rp * rm)) / 3.0) - ru;
    }

    // Newton iteration to fix numerical problems (using Horners method)
    // Suggested by @NinjaKoala
    for (int i=0; i < 2; i++) {
        float a_2 = ra + lambda;
        float a_1 = rb + lambda * a_2;
        float b_2 = a_2 + lambda;

        float f = rc + lambda * a_1; // Evaluation of λ^3 + ra * λ^2 + rb * λ + rc
        float f1 = a_1 + lambda * b_2; // Derivative

        lambda -= f / f1; // Newton iteration step
    }

    // Solve two quadratics factored from the quartic using the cubic root
    if (lambda < 0.0) return n;
    float t = sqrt(lambda); // Because we solved for t^2 but want t
    float alpha = 2.0 * q / t, beta = lambda + ra;

    float u = 0.25 * b;
    t *= 0.5;

    float z = -alpha - beta;
    if (z > 0.0) {
        z = sqrt(z) * 0.5;
        float h = +t - u;
        roots.xy = vec2(h + z, h - z);
        n += 2;
    }

    float w = +alpha - beta;
    if (w > 0.0) {
        w = sqrt(w) * 0.5;
        float h = -t - u;
        roots.zw = vec2(h + w, h - w);
        if (n == 0) roots.xy = roots.zw;
        n += 2;
    }

    return n;
}

/**
 * @brief intersect between a Ray and a Box
 * 
 * @param ray : The ray
 * @param bx : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectBox(Ray ray, Box bx, out Hit x) {
    vec3 t0, t1;
    //défini tout les cas de t0 et t1 pour x puis y puis z dans tous les angles de la cam
    if (ray.direction.x >= 0.) { 
        t0.x = (bx.bottomCorner.x - ray.origin.x) / ray.direction.x;
        t1.x = (bx.topCorner.x - ray.origin.x) / ray.direction.x;
    } 
    else { 
        t0.x = (bx.topCorner.x - ray.origin.x) / ray.direction.x;
        t1.x = (bx.bottomCorner.x - ray.origin.x) / ray.direction.x;
    } 
    
    if (ray.direction.y >= 0.) { 
        t0.y = (bx.bottomCorner.y - ray.origin.y) / ray.direction.y;
        t1.y = (bx.topCorner.y - ray.origin.y) / ray.direction.y;
    } 
    else { 
        t0.y = (bx.topCorner.y - ray.origin.y) / ray.direction.y;
        t1.y = (bx.bottomCorner.y - ray.origin.y) / ray.direction.y;
    } 
    
    if (ray.direction.z >= 0.) { 
        t0.z = (bx.bottomCorner.z - ray.origin.z) / ray.direction.z;
        t1.z = (bx.topCorner.z - ray.origin.z) / ray.direction.z;
    } 
    else { 
        t0.z = (bx.topCorner.z - ray.origin.z) / ray.direction.z;
        t1.z = (bx.bottomCorner.z - ray.origin.z) / ray.direction.z;
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
    if (t > 0.) {
        vec3 p = Point(ray, t);
        //le centre du box
        vec3 c = (bx.bottomCorner + bx.topCorner) / 2.;
        x = Hit(t, normalize(p - c), bx.id);
        return true;
    }
    return false;
}

/**
 * @brief intersect between a Ray and a Torus
 * 
 * @param ray : The ray
 * @param tor : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectTorus(Ray ray, Torus tor, out Hit x) { 

    // http://heigeas.free.fr/laure/ray_tracing/tore.html

    if (tor.bigRadius > tor.smallRadius) {
        vec3 oc = ray.origin - tor.center;

        float f = dot(ray.direction, ray.direction);
        float g = 2. * dot(oc, ray.direction); 
        float h = dot(oc, oc) + (tor.bigRadius * tor.bigRadius) - (tor.smallRadius * tor.smallRadius);

        float i = - 4. * (tor.bigRadius * tor.bigRadius) * ((ray.direction.x * ray.direction.x) + (ray.direction.y * ray.direction.y));
        float j = - 8. * (tor.bigRadius * tor.bigRadius) * ((oc.x * ray.direction.x) + (oc.y * ray.direction.y)); 
        float k = - 4. * (tor.bigRadius * tor.bigRadius) * ((oc.x * oc.x) + (oc.y * oc.y)); 

        float a = f * f;
        float b = 2. * f * g;
        float c = (2. * f * h) + (g * g) + i;
        float d = (2. * g * h) + j;
        float e = (h * h) + k;

        vec4 roots;
        int nroots = solveQuartic(a, b, c, d, e, roots); 

        if (nroots > 0) {
            float t1, t2, t;
            if (nroots == 2) {
                t =  min(roots.x, roots.y);
            }

            if (nroots == 4) {
                t1 = min(roots.x, roots.y);
                t2 = min(roots.z, roots.w);
                t = min (t1, t2);
            }

            if (t > 0.) {
                vec3 p = Point(ray, t);
                vec3 normale;

                vec3 odt = oc + ray.direction * t;
                float oxdt = oc.x + ray.direction.x * t;
                float oydt = oc.y + ray.direction.y * t;
                float ozdt = oc.z + ray.direction.z * t;

                normale.x = 4. * oxdt * (dot(oc, oc) + 2. * dot(oc, ray.direction * t) + dot(ray.direction * t, ray.direction * t) + dot(tor.bigRadius, tor.bigRadius) - dot(tor.smallRadius, tor.smallRadius)) - 8. * dot(tor.bigRadius, tor.bigRadius) * oxdt;
                normale.y = 4. * oydt * (dot(oc, oc) + 2. * dot(oc, ray.direction * t) + dot(ray.direction * t, ray.direction * t) + dot(tor.bigRadius, tor.bigRadius) - dot(tor.smallRadius, tor.smallRadius)) - 8. * dot(tor.bigRadius, tor.bigRadius) * oydt;
                normale.z = 4. * ozdt * (dot(oc, oc) + 2. * dot(oc, ray.direction * t) + dot(ray.direction * t, ray.direction * t) + dot(tor.bigRadius, tor.bigRadius) - dot(tor.smallRadius, tor.smallRadius));
                
                x = Hit(t, normalize(normale),tor.id);
                return true;
            }
        }
        return false;
    }
}


/**
 * @brief intersect between a Ray and a Goursat Box
 * 
 * @param ray : The ray
 * @param goursat : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectGoursat(Ray ray, Goursat goursat, out Hit x) {

    // http://heigeas.free.fr/laure/ray_tracing/cubetroue.html

    vec3 oc = ray.origin - goursat.center;
    vec3 rayD2 = ray.direction * ray.direction;
    vec3 rayD3 = rayD2 * ray.direction;
    vec3 rayO2 = oc * oc;
    vec3 rayO3 = rayO2 * oc;

    float f = dot(rayD2, rayD2);
    float g = 4.* dot(oc, rayD3);
    float h = 6. * dot(rayO2, rayD2);
    float i = 4. * dot(rayO3, ray.direction);
    float j = dot(rayO2, rayO2);

    float k = -5. * dot(ray.direction, ray.direction);
    float l = -10. * dot(oc, ray.direction);
    float m = -5. * dot(oc, oc) + 3.5 * cos(iTime) + 14.5;

    float a = f;
    float b = g;
    float c = h + k;
    float d = i + l;
    float e = j + m;

    vec4 roots;
    int nroots = solveQuartic(a, b, c, d, e, roots);

    if (nroots > 0) {
        float t1, t2, t;
        if (nroots == 2) {
            t =  min(roots.x, roots.y);
        }

        if (nroots == 4) {
            t1 = min(roots.x, roots.y);
            t2 = min(roots.z, roots.w);
            t = min (t1, t2);
        }

        if (t > 0.) {
            vec3 p = Point(ray, t);
            vec3 normale;

            float ocX2 = oc.x * oc.x;
            float ocX3 = ocX2 * oc.x;
            
            float ocY2 = oc.y * oc.y;
            float ocY3 = ocY2 * oc.y;
            
            float ocZ2 = oc.z * oc.z;
            float ocZ3 = ocZ2 * oc.z;


            float dXt = ray.direction.x * t;
            float dYt = ray.direction.y * t;
            float dZt = ray.direction.z * t;

            float dXt2 = dXt * dXt;
            float dYt2 = dYt * dYt;
            float dZt2 = dZt * dZt;

            float dXt3 = dXt2 * dXt;
            float dYt3 = dYt2 * dYt;
            float dZt3 = dZt2 * dZt;

            normale.x = 4. * ( ocX3 + 3.*dot(ocX2, dXt) + 3.*dot(oc.x, dXt2) + dXt3) - 10. * (oc.x + ray.direction.x * t );
            normale.y = 4. * ( ocY3 + 3.*dot(ocY2, dYt) + 3.*dot(oc.y, dYt2) + dYt3) - 10. * (oc.y + ray.direction.y * t );
            normale.z = 4. * ( ocZ3 + 3.*dot(ocZ2, dZt) + 3.*dot(oc.z, dZt2) + dZt3) - 10. * (oc.z + ray.direction.z * t );

            x = Hit(t, normalize(normale - p), goursat.id);
            return true;
        }
    }
    return false;
}

/**
 * @brief intersect between a Ray and an Octaeder
 * 
 * @param ray : The ray
 * @param octa : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectOctaeder(Ray ray, Octaeder octa, out Hit x) {

    // https://mathcurve.com/surfaces/goursat/goursat.shtml

    vec3 oc = ray.origin - octa.center;
    vec4 roots;

    float rdx = ray.direction.x;
    float rdy = ray.direction.y;
    float rdz = ray.direction.z;

    float rdx2 = rdx * rdx;
    float rdy2 = rdy * rdy;
    float rdz2 = rdz * rdz;

    float odx = oc.x * rdx;
    float ody = oc.y * rdy;
    float odz = oc.z * rdz;

    float ocx2 = oc.x * oc.x;
    float ocy2 = oc.y * oc.y;
    float ocz2 = oc.z * oc.z;

    float l = 8. * (dot(rdx2, rdy2) + dot(rdy2, rdz2) + dot(rdz2, rdx2)); 
    float m = 8. * (2. * (dot(odx, rdy2) + dot(ody, rdx2) 
                        + dot(ody, rdz2) + dot(odz, rdy2) 
                        + dot(odz, rdx2) + dot(odx, rdz2)));

    float n =  8. * (dot(ocx2, rdy2) + dot(ocy2, rdx2)
                    + dot(ocy2, rdz2) + dot(ocz2, rdy2)
                    + dot(ocz2, rdx2) + dot(ocx2, rdz2)
                    + 4. * (dot(odx, ody) + dot(ody, odz) + dot(odz, odx)));
    float o = 8. * (2. * (dot(ocx2, ody) + dot(ocy2, odx) 
                        + dot(ocy2, odz) + dot(ocz2, ody) 
                        + dot(ocz2, odx) + dot(ocx2, odz) ));

    float p = 8. * (dot(ocx2, ocy2) + dot(ocy2, ocz2) + dot(ocz2, ocx2));

    float q = 5. * (dot(ray.direction, ray.direction));
    float r = 5. * (2. * (dot(oc, ray.direction)));
    float s = 5. * (dot(oc, oc)) + cos(iTime) * 67.5 - 82.5;

    float f = l;
    float g = m;
    float h = n + q;
    float i = o + r;
    float j = p + s;

    int nroots = solveQuartic(f, g, h, i, j, roots);

    if (nroots > 0) {
        float t1, t2, t;
        if (nroots == 2) {
            t =  min(roots.x, roots.y);
        }

        if (nroots == 4) {
            t1 = min(roots.x, roots.y);
            t2 = min(roots.z, roots.w);
            t = min (t1, t2);
        }

        if (t > 0.) {
            vec3 p = Point(ray, t);
            vec3 normale;

            normale.x = 16. * (oc.x + ray.direction.x * t) * ( dot(oc.y, oc.y) + 2. * dot(oc.y, ray.direction.y * t) + dot(ray.direction.y * t, ray.direction.y * t) 
                            + dot(oc.z, oc.z) + 2. * dot(oc.z, ray.direction.z * t) + dot(ray.direction.z * t, ray.direction.z * t)) 
                            + 10. * (oc.x + ray.direction.x * t);
            normale.y = 16. * (oc.y + ray.direction.y * t) * ( dot(oc.x, oc.x) + 2. * dot(oc.x, ray.direction.x * t) + dot(ray.direction.x * t, ray.direction.x * t) 
                            + dot(oc.z, oc.z) + 2. * dot(oc.z, ray.direction.z * t) + dot(ray.direction.z * t, ray.direction.z * t)) 
                            + 10. * (oc.y + ray.direction.y * t);
            normale.z = 16. * (oc.z + ray.direction.z * t) * ( dot(oc.x, oc.x) + 2. * dot(oc.x, ray.direction.x * t) + dot(ray.direction.x * t, ray.direction.x * t) 
                            + dot(oc.y, oc.y) + 2. * dot(oc.y, ray.direction.y * t) + dot(ray.direction.y * t, ray.direction.y * t)) 
                            + 10. * (oc.z + ray.direction.z * t);

            x = Hit(t, normalize(normale), octa.id);
            return true;
        }
    }
    return false;
}

/**
 * @brief intersect between a Ray and a PipeConnector 
 * 
 * @param ray : The ray
 * @param pip : Structure information
 * @param x : Returned intersection information
 * @return true 
 * @return false 
 */
bool IntersectPipeConnector(Ray ray, PipeConnector pip, out Hit x) {

    // https://mathcurve.com/surfaces/goursat/goursat.shtml

    vec3 oc = ray.origin - pip.center;
    vec4 roots;

    float rdx = ray.direction.x;
    float rdy = ray.direction.y;
    float rdz = ray.direction.z;

    float rdx2 = rdx * rdx;
    float rdy2 = rdy * rdy;
    float rdz2 = rdz * rdz;

    float odx = oc.x * rdx;
    float ody = oc.y * rdy;
    float odz = oc.z * rdz;

    float ocx2 = oc.x * oc.x;
    float ocy2 = oc.y * oc.y;
    float ocz2 = oc.z * oc.z;

    float l = 2. * (dot(rdx2, rdy2) + dot(rdy2, rdz2) + dot(rdz2, rdx2)); 
    float m = 2. * (2. * (dot(odx, rdy2) + dot(ody, rdx2) 
                        + dot(ody, rdz2) + dot(odz, rdy2) 
                        + dot(odz, rdx2) + dot(odx, rdz2)));

    float n =  2. * (dot(ocx2, rdy2) + dot(ocy2, rdx2)
                    + dot(ocy2, rdz2) + dot(ocz2, rdy2)
                    + dot(ocz2, rdx2) + dot(ocx2, rdz2)
                    + 4. * (dot(odx, ody) + dot(ody, odz) + dot(odz, odx)));
    float o = 2. * (2. * (dot(ocx2, ody) + dot(ocy2, odx) 
                        + dot(ocy2, odz) + dot(ocz2, ody) 
                        + dot(ocz2, odx) + dot(ocx2, odz) ));

    float p = 2. * (dot(ocx2, ocy2) + dot(ocy2, ocz2) + dot(ocz2, ocx2));

    float q = -3. * (dot(ray.direction, ray.direction));
    float r = -3. * (2. * (dot(oc, ray.direction)));
    float s = -3. * (dot(oc, oc)) + cos(iTime) * 67.5 - 82.5;

    float f = l;
    float g = m;
    float h = n + q;
    float i = o + r;
    float j = p + s;

    int nroots = solveQuartic(f, g, h, i, j, roots);

    if (nroots > 0) {
        float t1, t2, t;
        if (nroots == 2) {
            t =  min(roots.x, roots.y);
        }

        if (nroots == 4) {
            t1 = min(roots.x, roots.y);
            t2 = min(roots.z, roots.w);
            t = min (t1, t2);
        }

        if (t > 0.) {
            vec3 p = Point(ray, t);
            vec3 normale;

            normale.x = 4. * (oc.x + ray.direction.x * t) * ( dot(oc.y, oc.y) + 2. * dot(oc.y, ray.direction.y) 
                            + dot(ray.direction.y * t, ray.direction.y * t) + dot(oc.z, oc.z) 
                            + 2. * dot(oc.z, ray.direction.z) + dot(ray.direction.z * t, ray.direction.z * t)) 
                            + 6. * (oc.x + ray.direction.x * t);
            normale.y = 4. * (oc.y + ray.direction.y * t) * ( dot(oc.x, oc.x) + 2. * dot(oc.x, ray.direction.x) 
                            + dot(ray.direction.y * t, ray.direction.y * t) + dot(oc.z, oc.z) 
                            + 2. * dot(oc.z, ray.direction.z) + dot(ray.direction.z * t, ray.direction.z * t)) 
                            + 6. * (oc.y + ray.direction.y * t);
            normale.z = 4. * (oc.z + ray.direction.z * t) * ( dot(oc.x, oc.x) + 2. * dot(oc.x, ray.direction.x) 
                            + dot(ray.direction.y * t, ray.direction.y * t) + dot(oc.z, oc.z) 
                            + 2. * dot(oc.z, ray.direction.z) + dot(ray.direction.z * t, ray.direction.z * t))
                            + 6. * (oc.z + ray.direction.z * t);

            x = Hit(t, normalize(normale), pip.id);
            return true;
        }
    }
    return false;
}

/**
 * @brief Translates the origin of a ray by a specified point.
 * 
 * @param ray : The ray
 * @param p : point
 * @return Ray the translated ray 
 */
Ray Translation(Ray ray, vec3 p) {
    return Ray(ray.origin - p, ray.direction, false, false);
}

/**
 * @brief Translates the origin of a ray by a specified point.
 * 
 * @param ray : The ray vector in space
 * @param p : point
 * @return Ray the translated ray in space
 */
vec3 Translation(vec3 ray, vec3 p) {
    return ray - p;
}

/**
 * @brief Applies rotation to a given normal vector around a specified point.
 * 
 * This function applies rotations to a normal vector around the specified point 'tr'.
 * 
 * @param normal The normal vector to be rotated.
 * @param rot A vector representing the rotation angles around X, Y, and Z axes.
 * @param tr The center point around which the rotation is performed.
 * @return vec3 The rotated normal vector.
 */
vec3 RotationNormal(vec3 normal, vec3 rot, vec3 tr) {
    //construire les matrices de rotations
    mat3 rotationX = mat3(
        1., 0.      , 0.       ,
        0., cos(rot.x), -sin(rot.x),
        0., sin(rot.x), cos(rot.x)
    );
    mat3 rotationY = mat3(
        cos(rot.y), 0., -sin(rot.y),
        0.      , 1., 0.       ,
        sin(rot.y), 0., cos(rot.y)
    );
    mat3 rotationZ = mat3(
        cos(rot.z), -sin(rot.z), 0.,
        sin(rot.z), cos(rot.z) , 0.,
        0.      , 0.       , 1.
    );
    normal = rotationX * rotationY * rotationZ * normal;
    return normal;
}

/**
 * @brief Applies rotation to a ray around a specified point.
 * 
 * This function applies rotations to a given ray using the rotation angle around the specified point 'tr'.
 * 
 * @param ray The original ray to be rotated.
 * @param rot A vector representing the rotation angles around X, Y, and Z axes.
 * @param tr The center point around which the rotation is performed.
 * @return Ray The rotated ray.
 */
Ray Rotation(Ray ray, vec3 rot, vec3 tr) {
    //construire les matrices de rotations
    mat3 rotationX = mat3(
        1., 0.      , 0.       ,
        0., cos(rot.x), -sin(rot.x),
        0., sin(rot.x), cos(rot.x)
    );
    mat3 rotationY = mat3(
        cos(rot.y), 0., -sin(rot.y),
        0.      , 1., 0.       ,
        sin(rot.y), 0., cos(rot.y)
    );
    mat3 rotationZ = mat3(
        cos(rot.z), -sin(rot.z), 0.,
        sin(rot.z), cos(rot.z) , 0.,
        0.      , 0.       , 1.
    );
    //ramener à 0
    ray.direction = rotationZ * rotationY * rotationX * ray.direction;
    ray.origin = Translation(ray.origin, tr);
    //effectuer la rotation
    ray.origin = rotationZ * rotationY * rotationX * ray.origin;
    //ramener à où c'était
    ray.origin = Translation(ray.origin, -tr);
    return Ray(ray.origin, ray.direction, ray.isHomo, true);
}

/**
 * @brief Applies a homothety transformation to a ray around a specified point.
 * 
 * This function applies a homothety transformation to a given ray by changing its direction and origin,
 * relative to the specified point 'tr'.
 * 
 * @param ray The original ray to be transformed.
 * @param homo A vector representing scaling factors for the ray's direction.
 * @param tr The center point around which the homothety transformation is performed.
 * @return Ray The transformed ray.
 */
Ray Homothetie(Ray ray, vec3 homo, vec3 tr) {
    // change direction
    ray.direction = ray.direction / homo;

    // move to 0 origin
    ray.origin = Translation(ray.origin, tr);
    // change origin
    ray.origin = ray.origin / homo;
    // get back the origin
    ray.origin = Translation(ray.origin, -tr);
    // normalize the direction
    ray.direction = normalize(ray.direction);
    return Ray(ray.origin, ray.direction, true, ray.isRot);
}

/**
 * @brief Transforms a hit point obtained in homothetic space back to the original space.
 * 
 * This function takes a hit point 'homoHit' obtained in homothetic space and transforms it back
 * to the original space using the given baseRay and homoRay for reference. The 'scale' vector
 * is used to scale the transformed point.
 * 
 * @param homoHit The hit point in homothetic space.
 * @param baseRay The original ray used as a reference for transformation.
 * @param homoRay The transformed ray used as a reference for transformation.
 * @param scale A vector representing the scaling factors.
 * @return Hit The hit point in the original space.
 */
Hit Homothetie(Hit homoHit, Ray baseRay, Ray homoRay, vec3 scale) {
    // point d'intersection dans le repère non transformé
    vec3 homoPoint = Point(homoRay, homoHit.t);
    // point d'intersection sur le vrai objet transformé
    homoPoint *= scale;
    
    Hit hit = homoHit;
    // nouveau t à partir du vrai point transformé (on résout juste homoPoint = baseRay.o + t * baseRay.d, 
    // en choisissant n'importe quelle composante, ici x)
    hit.t = (homoPoint.x - baseRay.origin.x) / baseRay.direction.x;
    
    return hit;
}

/**
 * @brief Computes the center point between two given points.
 * 
 * @param cornerA The first point.
 * @param cornerB The second point.
 * @return vec3 The center point between 'cornerA' and 'cornerB'.
 */
vec3 FindTheCenter(vec3 cornerA, vec3 cornerB) {
    return (cornerA + cornerB)/2.;
}

/**
 * @brief Creates a scene of all the objects.
 * 
 * @param ray The ray
 * @return Dynamic scene.
 */
Scene sceneAllObjects(Ray ray){
    Scene scene;
    scene.plane = Plane(vec3(0.,0.,1.), vec3(0.,0.,0.),0);

    for (int i = 0; i < 10; i++) {
        scene.tabRay[i] = ray;
    }

    scene.nbSphere = 1;
    scene.nbEllipsoid = 1;
    scene.nbCylinder = 1;
    scene.nbCapsule = 1;
    scene.nbBox = 1;
    scene.nbTorus = 1;
    scene.nbGoursat = 1;
    scene.nbOctaeder = 1;
    scene.nbPipeConnector = 1;
    
    scene.tabSphere[0] = Sphere(vec3(-4.,5.,4.),1.,1);
    scene.tabEllipsoid[0] = Ellipsoid(vec3(-4., 0., 4.), vec3(1., 1.2, 0.5), 1);
    scene.tabCylinder[0] = Cylinder(vec3(0.,5.,0.), vec3(0., 5., 2.1), .75, 1);
    scene.tabCapsule[0] = Capsule(vec3(0.,2.,0.), vec3(0., 2., 3.1), 0.5 ,1);
    scene.tabBox[0] = Box(vec3(-1.1, -1.1, 0.1), vec3(2.5, 0.1, 2.1), 1);
    
    scene.tabRay[5] = Translation(scene.tabRay[0], vec3(0, -3, 2)); 
    scene.tabTorus[0] = Torus(vec3(0., 0., 0.), 1., .5, 1);
    
    scene.tabGoursat[0] = Goursat(vec3(0., -8., 2.1), 1);
    
    scene.tabOctaeder[0] = Octaeder(vec3(-10., -8., 4), 1);
    scene.tabPipeConnector[0] = PipeConnector(vec3(-23., -3, 10), 1);
    
    return scene;
}

/**
 * @brief Creates a scene where we test Translation, Homotethy and Rotation functions.
 * 
 * @param ray The ray
 * @return Dynamic scene.
 */
Scene sceneTransformation(Ray ray) {
    Scene scene;
    scene.plane = Plane(vec3(0.,0.,1.), vec3(0.,0.,0.),4);
    
    for (int i = 0; i < 9; i++) {
        scene.tabRay[i] = ray;
    }

    // //Box homothethy
    scene.nbBox = 1;
    scene.tabBox[0] = Box(vec3(-4., -1., 2.), vec3(-2, 1., 4.), 2);
    scene.tabScale[0] = vec3(1., sin(iTime+1.5*3.1415)+2., cos(iTime)+2.);
    scene.tabRay[0] = Homothetie(scene.tabRay[0], scene.tabScale[0], FindTheCenter(scene.tabBox[0].bottomCorner, scene.tabBox[0].topCorner));

    scene.nbTorus = 1;
    scene.tabTorus[0] = Torus(vec3(0., 0., 0.), 1., .2, 7);
    scene.tabRay[1] = Translation(scene.tabRay[1], vec3(0., 0., 4.));
    scene.tabScale[1] = vec3(1., 1., 4.);
    scene.tabAngle[1] = vec3(0.5*-iTime, 0., 0.);
    scene.tabRay[1] = Rotation(scene.tabRay[1], scene.tabAngle[1], scene.tabTorus[0].center);
    scene.tabRay[1] = Homothetie(scene.tabRay[1], scene.tabScale[1], scene.tabTorus[0].center);

    scene.nbGoursat = 1;
    scene.tabGoursat[0] = Goursat(vec3(0, 6, 1), 5); 
    scene.tabScale[2] = vec3(1.,2.,1.);
    scene.tabRay[2] = Homothetie(scene.tabRay[2], scene.tabScale[2], scene.tabGoursat[0].center);

    scene.nbOctaeder = 1;
    scene.tabOctaeder[0] = Octaeder(vec3(0, -6, 2.2), 1);
    scene.tabAngle[3] = vec3(0., iTime, 10.*iTime);
    scene.tabRay[3] = Rotation(scene.tabRay[3], scene.tabAngle[3], scene.tabOctaeder[0].center);

    return scene;
}

/**
 * @brief Creates a scene with spheres.
 * 
 * @param ray The ray
 * @return Dynamic scene.
 */
Scene sceneSphere(Ray ray) {
    Scene scene;
    scene.plane = Plane(vec3(0.,0.,1.), vec3(0.,0.,0.),0);
    
    for (int i = 0; i < 10; i++) {
        scene.tabRay[i] = ray;
    }

     scene.nbSphere = 4;
     scene.tabSphere[0] = Sphere(vec3(-3, 0, 1), 1.5, 6);
     scene.tabSphere[1] = Sphere(vec3(0.,-3.,1.),1.5,8);
     scene.tabSphere[2] = Sphere(vec3(3.,0.,1.),1.5,6);
     scene.tabSphere[3] = Sphere(vec3(0.,3.,1.),1.5,7);

    return scene;
}

/**
 * @brief Creates a scene where we display our Textures.
 * 
 * @param ray The ray
 * @return Dynamic scene.
 */
Scene sceneTexture(Ray ray) {
    Scene scene;
    scene.plane = Plane(vec3(0.,0.,1.), vec3(0.,0.,0.),0);

    for (int i = 0; i < 10; i++) {
        scene.tabRay[i] = ray;
    }

    scene.nbSphere = 1;
    scene.nbEllipsoid = 1;
    scene.nbCylinder = 1;
    scene.nbCapsule = 1;
    scene.nbBox = 1;
    scene.nbTorus = 1;
    scene.nbGoursat = 0;
    scene.nbOctaeder = 0;
    scene.nbPipeConnector = 0;
    
    scene.tabSphere[0] = Sphere(vec3(0.,-5.,1.), 2.,7);
    scene.tabEllipsoid[0] = Ellipsoid(vec3(-4., -5., 4.5), vec3(1.,1.,0.5), 4);
    scene.tabCylinder[0] = Cylinder(vec3(0.,3.,0.), vec3(0., 3., 2.1), .75, 3);
    scene.tabCapsule[0] = Capsule(vec3(-2.,3.,2.), vec3(-4., 3., 3.), 0.5 ,2);
    scene.tabBox[0] = Box(vec3(-4., -1.5, 0.), vec3(-2.1, 1.5, 5.), 5);
    scene.tabTorus[0] = Torus(vec3(0., 0., 0.), 1., .5, 6);
    scene.tabOctaeder[0] = Octaeder(vec3(0., 5., 3.), 1);
    
    return scene;
}

Scene sceneMultiLight(Ray ray) {
     /* RAPPEL : 
        ACTIVER MULTILIGHT DANS COLOR()
        LES POSITIONS DE LUMIERES ET LEURS COULEURS SONT GERES DANS LA FONCTION COLOR()
    */
    Scene scene;
    scene.plane = Plane(vec3(0.,0.,1.), vec3(0.,0.,0.),0);

    for (int i = 0; i < 10; i++) {
        scene.tabRay[i] = ray;
    }

    scene.nbSphere = 2;
    scene.tabSphere[0] = Sphere(vec3(0.,-2.,1.), 1.,1);
    scene.tabSphere[1] = Sphere(vec3(0.,2.,1.), 1.,1);

    return scene;
}

Scene sceneOA(Ray ray) {
    /* RAPPEL : 
        POUR TESTER L'OA IL FAUT COMMENTER LA LIGNE `x = Hit(1000.,vec3(0),-1);`
        DE LA FONCTION SHADE ET IL FAUT LA DECOMMENTER DE LA FONCTION INTERSECT
    */
    Scene scene;
    scene.plane = Plane(vec3(0.,0.,1.), vec3(0.,0.,0.),0);

    for (int i = 0; i < 10; i++) {
        scene.tabRay[i] = ray;
    }

    scene.nbGoursat = 1;
    scene.tabGoursat[0] = Goursat(vec3(0.,0.,2.5),1);

    return scene;
}

/**
 * @brief Performs ray intersection with the scene.
 *
 * @param ray The ray to intersect with the scene object.
 * @param x Reference to the Hit structure, stores intersection information.
 * @return True if there is an intersection, false otherwise.
 */
bool Intersect(Ray ray, inout Hit x) {

    /* RAPPEL : 
    ON DECOMENTE ICI ET ON COMMENTE DANS LA FONCTION SHADE() ET ON PEUT VOIR L'OA 
    NOUBLIEZ DE MODIFER LA FONCTON COLOR POUR L'OA AUSSI
    */
    // x = Hit(1000., vec3(0.), -1);

    Scene scene = sceneAllObjects(ray);
    Hit current;
    bool ret=false;
    int idR = 0;

    if (IntersectPlane(ray,scene.plane,current) && current.t<x.t) {
        x = current;
        ret=true;
    }

    for (int i = 0; i < scene.nbSphere; i++) {
        if (IntersectSphere(scene.tabRay[i + idR], scene.tabSphere[i], current)) {
            if (scene.tabRay[i+idR].isHomo)
                current = Homothetie(current, ray, scene.tabRay[i+idR], scene.tabScale[i+idR]);
            if (current.t < x.t){
                x = current;
                if (scene.tabRay[i+idR].isRot)
                    x.normal = RotationNormal(x.normal, -scene.tabAngle[i+idR], scene.tabSphere[i].center);
                ret=true;
            }
        }
    }
    idR += scene.nbSphere;
    
    for (int i = 0; i < scene.nbEllipsoid; i++) {
        if (IntersectEllipsoid(scene.tabRay[i + idR],scene.tabEllipsoid[i],current) && current.t<x.t) {
            if (scene.tabRay[i+idR].isHomo)
                current = Homothetie(current, ray, scene.tabRay[i+idR], scene.tabScale[i+idR]);
            if (current.t < x.t){
                x = current;
                if (scene.tabRay[i+idR].isRot)
                    x.normal = RotationNormal(x.normal, -scene.tabAngle[i+idR], scene.tabEllipsoid[i].center);
                ret=true;
            }
        }
    }
    idR += scene.nbEllipsoid;

    for (int i = 0; i < scene.nbCylinder; i++) {
        if (IntersectCylinder(scene.tabRay[i + idR], scene.tabCylinder[i],current) && current.t<x.t) {
            if (scene.tabRay[i+idR].isHomo)
                current = Homothetie(current, ray, scene.tabRay[i+idR], scene.tabScale[i+idR]);
            if (current.t < x.t){
                x = current;
                if (scene.tabRay[i+idR].isRot)
                    x.normal = RotationNormal(x.normal, -scene.tabAngle[i+idR], FindTheCenter(scene.tabCylinder[i].bottomCenter, scene.tabCylinder[i].topCenter));
                ret=true;
            }
        }
    }
    idR += scene.nbCylinder;

    for (int i = 0; i < scene.nbCapsule; i++) {
        if (IntersectCapsule(scene.tabRay[i + idR], scene.tabCapsule[i],current)&&current.t<x.t) {
            if (scene.tabRay[i+idR].isHomo)
                current = Homothetie(current, ray, scene.tabRay[i+idR], scene.tabScale[i+idR]);
            if (current.t < x.t){
                x = current;
                if (scene.tabRay[i+idR].isRot)
                    x.normal = RotationNormal(x.normal, -scene.tabAngle[i+idR], FindTheCenter(scene.tabCapsule[i].bottomCenter, scene.tabCapsule[i].topCenter));
                ret=true;
            }
        }
    }
    idR += scene.nbCapsule;

    for(int i = 0; i < scene.nbBox; i++) {
        if (IntersectBox(scene.tabRay[i+idR], scene.tabBox[i], current)) {
            if (scene.tabRay[i+idR].isHomo)
                current = Homothetie(current, ray, scene.tabRay[i+idR], scene.tabScale[i+idR]);
            if (current.t < x.t) {
                x = current;
                if(scene.tabRay[i+idR].isRot)
                    x.normal = RotationNormal(x.normal, -scene.tabAngle[i+idR], FindTheCenter(scene.tabBox[i].bottomCorner, scene.tabBox[i].topCorner));
                ret=true;
            }
        }
    }
    idR += scene.nbBox;

    for (int i = 0; i < scene.nbTorus; i++) {
        if (IntersectTorus(scene.tabRay[i + idR], scene.tabTorus[i], current)) {
            if(scene.tabRay[i + idR].isHomo)
                current = Homothetie(current, ray, scene.tabRay[i + idR], scene.tabScale[i + idR]);
            if(current.t < x.t){
                x = current;
                if(scene.tabRay[i + idR].isRot)
                    x.normal = RotationNormal(x.normal, -scene.tabAngle[i + idR], scene.tabTorus[i].center);
                ret=true;
            }
        }
    }
    idR += scene.nbTorus;

    for (int i = 0; i < scene.nbGoursat; i++) {
        if (IntersectGoursat(scene.tabRay[i + idR], scene.tabGoursat[i], current)) {
            if (scene.tabRay[i+idR].isHomo)
                current = Homothetie(current, ray, scene.tabRay[i + idR], scene.tabScale[i+idR]);
            if (current.t < x.t){
                x = current;
                if (scene.tabRay[i+idR].isRot)
                    x.normal = RotationNormal(x.normal, -scene.tabAngle[i + idR], scene.tabGoursat[i].center);
                ret=true;
            }
        }
    }
    idR += scene.nbGoursat;

    for (int i = 0; i < scene.nbOctaeder; i++) {
        if (IntersectOctaeder(scene.tabRay[i + idR], scene.tabOctaeder[i], current)) {
            if (scene.tabRay[i + idR].isHomo)
                current = Homothetie(current, ray, scene.tabRay[i + idR], scene.tabScale[i+idR]);
            if (current.t < x.t){
                x = current;
                if (scene.tabRay[i + idR].isRot)
                    x.normal = RotationNormal(x.normal, -scene.tabAngle[i + idR], scene.tabOctaeder[i].center);
                ret=true;
            }
        }
    }
    idR += scene.nbOctaeder;

    for (int i = 0; i < scene.nbPipeConnector; i++) {
        if (IntersectPipeConnector(ray, scene.tabPipeConnector[i], current)) {
            if (scene.tabRay[i + idR].isHomo)
                current = Homothetie(current, ray, scene.tabRay[i + idR], scene.tabScale[i + idR]);
            if (current.t < x.t){
                x = current;
                if (scene.tabRay[i + idR].isRot)
                    x.normal = RotationNormal(x.normal, -scene.tabAngle[i + idR], scene.tabPipeConnector[i].center);
                ret=true;
            }
        }
    }
    idR += scene.nbPipeConnector;
    return ret;
}

/**
 * @brief Generates a random direction within a hemisphere oriented by a given normal.
 * 
 * @param seed An integer seed value to determine the randomness of the direction.
 * @param normal The normal vector representing the orientation of the hemisphere.
 * @return vec3 A random direction vector within the hemisphere.
 */
vec3 Hemisphere(int seed,vec3 normal) {
    float a, b;
    a = fract(sin(176.19 * float(seed)));// Uniform randoms
    b = fract(sin(164.19 * float(seed)));
    float u = 2. * 3.1415 * a;// Random angle
    float v = acos(2. * b - 1.);// Arcosine distribution to compensate for poles
    vec3 d = vec3(cos(u) * cos(v), sin(u) * cos(v), sin(v));// Direction
    if (dot(d, normal) < 0.) { // Hemishper
        d = -d; 
    }
    return d;
}

/**
 * @brief Computes ambient occlusion at a given point on a surface.
 * 
 * @param point The point on the surface for which ambient occlusion is calculated.
 * @param normal The surface normal at the specified point.
 * @param numSamples The number of samples to use for occlusion calculation.
 * @return float The ambient occlusion value in the range [0, 1].
 */
float AmbientOcclusion(vec3 point, vec3 normal, int numSamples) {
    if (numSamples == 0) {
        return 1.;
    }

    float occ = 0.0;
    
    for (int i = 0; i < numSamples; i++) {
        // Generate a random direction in a hemisphere around the normal
        vec3 hemisphereDir = Hemisphere(i, normal);
        
        // Create a ray from the point in the direction of the hemisphereDir
        Ray occRay = Ray(point + normal * 0.001, hemisphereDir, false, false);
        
        // Check for intersections with scene objects
        Hit occHit;
        occHit = Hit(.1, vec3(0.), -1);
        bool hit = Intersect(occRay, occHit);
        
        // If no intersection, increase occlusion
        if (hit && occHit.t < 3. ) {
            occ += 1.0;
        }
    }

    // Normalize occlusion value
    occ /= float(numSamples) * 4.5;
    
    return occ;
}

/**
 * @brief Calculates the background color based on the viewing direction.
 * 
 * @param rd The viewing direction vector.
 * @return vec3 The background color for the given viewing direction.
 */
vec3 Background(vec3 rd) {
    return mix(vec3(.8,.8,.9),vec3(.7,.7,.8),rd.z);
}

/**
 * @brief Creates a camera rotation matrix based on the camera origin and target point.
 * 
 * @param ro The camera origin.
 * @param ta The target point that the camera is looking at.
 * @return mat3 The camera rotation matrix.
 */
mat3 setCamera(in vec3 ro, in vec3 ta) {
    vec3 cw = normalize(ta - ro);
    vec3 cp = vec3(0, 0, 1);
    vec3 cu = -normalize(cross(cw, cp));
    vec3 cv = -normalize(cross(cu, cw));
    return mat3(cu, cv, cw);
}

/**
 * @brief Sets up an array of lights with evenly distributed positions and colors.
 *
 * @param l an array of Light structures representing the lights. It is modified in place.
 * @param n The number of lights in one dimension (total number of lights = n * n).
 */
void MultiLight(inout Light l[25], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int index = i * n + j;
            l[index].lightColor = vec3(1.0 / float(n * n), 1.0 / float(n * n), 1.0 / float(n * n));
            l[index].lightPos = vec3(1.0 + float(i) / 10.0, 1.0 + float(j) / 10.0, 4.0);
        }
    }
}

/**
 * @brief Apply Phong Model to illuminate the scene and add object Shadows. 
 * Can also permit you to see the OA with `if(true)` and `m.ambient` returned
 *
 * @param m Object material.
 * @param n The surface normal at the specified point.
 * @param p The point on the surface.
 * @param camera The camera ray.
 * @return The computed color for the specified point on the surface.
 */
vec3 Color(Material m, vec3 n, vec3 p, Ray camera) {
    /*RAPPEL : 
        Eclairage :     
            UTILISES MULTILIGHT POUR LES OMBRES DEGRADES OU BIEN UTILISES 2 LUMIERES 
            POUR ECLAIRCIR LA SCENE

        Voir OA : N'OUBLIEZ PAS D'AJOUTER `x = Hit(1000., vec3(0), -1);` dans Intersect()
            et l'enlever de Shade(). 
            Decommenter le `if (true)` et commenter l'autre 
            `if (!Intersect(r, randomHit) ...)` et retourner m.ambient au lieu 
            de finalColor
    */
    Scene scene;

    // OMBRE DEGRADE
    // scene.nbLight = 25;
    // MultiLight(scene.tabLight, 5);

    // NOMBRE DE LUMIERES INITIALISES AVEC LES COORDOONES DE LUMIERES
    scene.nbLight = 2;
    scene.tabLight[0].lightPos = vec3(0,0,4);
    scene.tabLight[0].lightColor = vec3(1,1,1);
    scene.tabLight[1].lightPos = vec3(4,4,2);
    scene.tabLight[1].lightColor = vec3(1,1,1);

    // Rotation de lumiere sur l'axe z
    Ray rotLight = Rotation(Ray(scene.tabLight[0].lightPos, vec3(0), false, false), vec3(0, 0, iTime), vec3(1, 0, 0));
    // Pour faire une rotation sur la lumiere decommente cette ligne de code
    // scene.tabLight[0].lightPos = rotLight.origin; 
    // scene.tabLight[1].lightPos = rotLight.origin; 

    Hit randomHit;
    vec3 finalColor;

    // Je dois calculer la lumiere spectrale et la lumiere diffus           
    vec3 camDirection = normalize(camera.origin - p); // Direction de la caméra  
    for (int i = 0; i < scene.nbLight; i++) {
         vec3 lightDirection = normalize(scene.tabLight[i].lightPos - p);
         Ray r = Ray(p + n * 0.001, lightDirection, false, false); // le rayon que j'envoie de point d'intersect de mon objet

        randomHit = Hit( length(scene.tabLight[i].lightPos - p), vec3(0.), -1);
    
        if (!Intersect(r, randomHit) && randomHit.t >= length(scene.tabLight[i].lightPos - p)) {
        //if (true) {
            vec3 reflectDirection = reflect(-lightDirection, n);// Direction de réflexion de lumiere depuis mon objet   
            
            float spec = pow(max(dot(camDirection, reflectDirection), 0.0),  m.coefShininess);// coefficient de shininess contrôle la netteté du reflet             
            vec3 specularColor = m.specular * spec * scene.tabLight[i].lightColor;// Éclairage speculaire             
            
            float diff = max(dot(n, lightDirection), 0.0);     
            vec3 diffuseColor = m.diffuse * diff * scene.tabLight[i].lightColor;// Éclairage diffus

            // Couleur a retourner
            finalColor += (specularColor * (1. / float(scene.nbLight)) + diffuseColor);// / float(scene.nbLight);
            // finalColor = m.ambient;
        } else {
            finalColor += vec3(0,0,0); 
        }   
    }

    float ao = AmbientOcclusion(p, n, 64); 
    return finalColor * (1. - ao );
}

// Rendering
vec3 Shade(Ray ray) {
    vec3 resultColor = vec3(0.0);

    for (int reflection = 0; reflection < MAX_REFLECTION; reflection++) {
        Hit x;
        x = Hit(1000., vec3(0), -1);
        bool idx = Intersect(ray, x);

        if (idx) {
            vec3 p = Point(ray, x.t);
            Material mat = Texture(p, x.id);

            // return x.normal;//débug normale

            // Si le matériau est un miroir, calculez la couleur réfléchie
            if (mat.reflexivity > 0.) {
                vec3 n = x.normal;
                vec3 reflectDir = reflect(ray.direction, n);
                ray = Ray(p + n * 0.001, reflectDir, false, false);
                resultColor += (1.0 - mat.reflexivity) * Color(mat, n, p, ray) + mat.mirrorColor;
            }
            else {
                // Si ce n'est pas un miroir, calculez la couleur normale et terminez la boucle
                resultColor += Color(mat, x.normal, p, ray);
                break;
            }
        }
        else {
            // Si aucune intersection n'est trouvée, retournez la couleur d'arrière-plan
            resultColor += Background(ray.direction);
            break;
        }
    }

    return resultColor;
}

vec3 ColorOriginal(Material m,vec3 n) {
    vec3 light = normalize(vec3(10, 2, 5));
    float diff = clamp(dot(n, light), 0., 1.);
    vec3 col = m.diffuse * diff + vec3(.2,.2,.2);
    return col;
}

vec3 ShadeOriginal(Ray ray) {
    // Intersect contains all the geo detection
    Hit x;

    //DECOMMENTER SI VOUS ACTIVER L'OA SEULEMENT, A REMETTRE PAR LA SUITE
    x = Hit(1000., vec3(0), -1);
    bool idx=Intersect(ray,x);
    
    if(idx) {
        vec3 p=Point(ray,x.t);
        Material mat = Texture(p, x.id);
        
        return ColorOriginal(mat, x.normal);
    }
    else {
        return Background(ray.direction);
    }
    return vec3(0);
}

void mainImage(out vec4 fragColor,in vec2 fragCoord) {
    // From uv which are the pixel coordinates in [0,1], change to [-1,1] and apply aspect ratio
    vec2 uv=(-iResolution.xy + 2.* fragCoord.xy) / iResolution.y;
    
    // Mouse control
    vec2 mouse = iMouse.xy / iResolution.xy;
    
    // Ray origin
    //défini la position de la cam
    vec3 ro=13. * normalize(vec3(sin(2. * 3.14 * mouse.x), cos(2. * 3.14 * mouse.x), 1.4 * (mouse.y - .1)));
    vec3 ta=vec3(0., 0., 1.5);
    mat3 ca=setCamera(ro, ta);
    
    // Ray
    vec3 rd = ca*normalize(vec3(uv.xy * tan(radians(22.5)), 1.));
    
    // Render
    vec3 col = Shade(Ray(ro, rd, false, false));
    
    fragColor = vec4(col,1.);
}
