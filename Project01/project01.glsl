// -----------------------------------------------------------------------------------

/*

MIT License

Copyright (c) 2021 Gerard Marcos Freixas

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

// -----------------------------------------------------------------------------------

/*
 *    For clarifying which variables are meant to be output parameters from a function.
 */
#define OUT

/*
 *   Undefine for no anti-aliasing.
 */
#define ANTIALIASING

/*
 *   Octaves number for the fractional brownian motion.
 */
#define NUM_OCTAVES 5

/*
 *   Maximum number of objects, for each type.
 */
#define MAX_PLANES          1
#define MAX_SPHERES         2
#define MAX_POINTLIGHTS     1

/*
 *   Object type.
 */
#define TYPE_PLANE          0
#define TYPE_SPHERE         1

/*
 *   Plane type.
 */
#define PLANE_CHECKERS      1
#define PLANE_CLOUDS        2

/*
 *   Material type.
 */
#define MATERIAL_REFLECTIVE 0
#define MATERIAL_REFRACTIVE 1

// -----------------------------------------------------------------------------------

/*
 *   The camera allows to see the scene from any point of view.
 */
struct Camera
{
    vec3 position;
    
    vec3 look_at;
    vec3 forward;
    vec3 right;
    vec3 up;

    float zoom;
};

/*
 *   Either a plane or a sphere, have its own material.
 */
struct Material 
{
    float diffuse;
    float specular;
    float shininess;
    float ambience;
    float reflection;
    
    int type; // 0 for reflective, 1 for refractive
};

/*
 *   The plane object.
 */
struct Plane 
{
    vec3 position;
    vec3 normal;
    vec3 color;

    int type;
    
    Material material;
};

/*
 *   The sphere object.
 */
struct Sphere
{
    vec3 position;
    float radius;
    vec3 color;
    
    Material material;
};

/*
 *   The light object.
 */
struct PointLight
{
    vec3 position;
    vec3 color; // TODO
    float intensity;
};

// -----------------------------------------------------------------------------------

Camera camera = Camera(vec3(0.0, 0.0, -2.0), vec3(0.0), vec3(0.0), vec3(0.0), vec3(0.0), 1.0);

const Material material01 = Material(1.0, 1.0, 100.0, 1.0, 1.0, MATERIAL_REFLECTIVE);
const Material material02 = Material(0.4, 1.0, 120.0, 0.3, 1.0, MATERIAL_REFLECTIVE);
const Material material03 = Material(0.2, 1.0, 100.0, 0.4, 1.0, MATERIAL_REFLECTIVE);
const Material material04 = Material(0.2, 1.0, 100.0, 0.4, 1.0, MATERIAL_REFLECTIVE);
const Material material05 = Material(0.3, 1.0, 90.0, 0.5, 1.0, MATERIAL_REFLECTIVE);
const Material material06 = Material(0.5, 0.5, 77.0, 0.7, 1.0, MATERIAL_REFLECTIVE);
const Material material07 = Material(0.4, 0.25, 124.0, 0.2, 1.0, MATERIAL_REFLECTIVE);

Plane plane01 = Plane(vec3(0.0, -0.2, 0.0), vec3(0.0, 1.0, 0.0), vec3(0.0), PLANE_CHECKERS, material01);
//Plane plane02 = Plane(vec3(0.0, 0.5, 0.0), vec3(0.0, -1.0, 0.0), vec3(0.0, 0.2, 0.8), PLANE_CLOUDS, material02);

Sphere sphere01 = Sphere(vec3(0.1, 0.0, 0.0), 0.07, vec3(1.0, 0.5, 0.3), material06);
Sphere sphere02 = Sphere(vec3(-0.1, 0.0, 0.0), 0.09, vec3(0.5, 0.3, 0.5), material07);

PointLight pointlight01 = PointLight(vec3(0.0, 0.2, -0.1), vec3(1.0, 1.0, 1.0), 10.0);

Plane planes[MAX_PLANES];
Sphere spheres[MAX_SPHERES];
PointLight pointlights[MAX_POINTLIGHTS];

const int bounces = 2;

// -----------------------------------------------------------------------------------

void setup()
{
    planes[0] = plane01;
    //planes[1] = plane02;

    spheres[0] = sphere01;
    spheres[1] = sphere02;
    
    pointlights[0] = pointlight01;
}

/*
 *   Solves a quadratic equation.
 *
 *   @param[IN] 'a' parameter of the quadratic equation.
 *   @param[IN] 'b' parameter of the quadratic equation.
 *   @param[IN] 'c' parameter of the quadratic equation.
 *   @param[OUT] First solution of the equation.
 *   @param[OUT] Second solution of the equation.
 *
 *   @return True if the equation is solvable, false if the equation is unsolvable.
 */
bool solve_quadratic(in float a, in float b, in float c, out float t0, out float t1)
{
    float disc = b * b - 4.0 * a * c;
    
    if (disc < 0.0)
    {
        return false;
    }
    
    if (disc == 0.0)
    {
        t0 = t1 = -b / (2.0 * a);
        return true;
    }
    
    t0 = (-b + sqrt(disc)) / (2.0 * a);
    t1 = (-b - sqrt(disc)) / (2.0 * a);
    
    return true;
}

float rand(in vec2 n)
{
	return fract(sin(dot(n, vec2(12.9898, 4.1414))) * 43758.5453);
}

float noise(in vec2 n)
{
    const vec2 d = vec2(0.0, 1.0);
    vec2 b = floor(n), f = smoothstep(vec2(0.0), vec2(1.0), fract(n));
	
    return mix(mix(rand(b), rand(b + d.yx), f.x), mix(rand(b + d.xy), rand(b + d.yy), f.x), f.y);
}

float fbm(in vec2 x)
{
	float v = 0.0;
	float a = 0.5;
	vec2 shift = vec2(100);
	// Rotate to reduce axial bias
    mat2 rot = mat2(cos(0.5), sin(0.5), -sin(0.5), cos(0.5));
	for (int i = 0; i < NUM_OCTAVES; ++i)
    {
		v += a * noise(x);
		x = rot * x * 2.0 + shift;
		a *= 0.5;
	}
    
	return v;
}

/*
 *   Get the material of an object.
 *
 *   @param[IN] Index of the object.
 *   @param[IN] Type of object.
 *
 *   @return The material of the desired object.
 */
Material get_material(in int index, in int type)
{
    switch(type)
    {
        case TYPE_PLANE:
        {
            return planes[index].material;
        }
        break;

        case TYPE_SPHERE:
        {
            return spheres[index].material;
        }
        break;
        
        default:
        {} break;
    }
}

vec3 calculate_fresnel(in vec3 surface_normal, in vec3 surface_point_position, in float intensity)
{
    return vec3(intensity + (1.0 - intensity) * pow(1.0 - dot(-surface_normal, normalize(surface_point_position - camera.position)), 5.0));
}

vec3 get_color(in vec3 viewDir, in vec3 surfacePointPosition, in vec3 objectColor, in PointLight pointLight, in vec3 surfaceNormal, in Material material)
{
    vec3 lightVector = surfacePointPosition - pointLight.position;
    vec3 lightDir = normalize(lightVector);   
    
   	float lightIntensity = (pow(0.1, 2.0) / pow(sqrt(dot(lightVector, lightVector)), 2.0)) * pointLight.intensity;
    
    float coeff = -dot(lightDir, surfaceNormal);
    
    vec3 ambient = material.ambience * objectColor;

    vec3 diffuse = material.diffuse * max(coeff, 0.0) * objectColor * lightIntensity;
       
    vec3 halfwayDir = normalize(lightDir + viewDir);
    vec3 specular = pow(max(-dot(surfaceNormal, halfwayDir), 0.0), material.shininess) * material.specular * objectColor * lightIntensity;
    
    vec3 fresnel = calculate_fresnel(surfaceNormal, surfacePointPosition, 0.01);

    vec3 color = diffuse + specular + ambient + fresnel;
    //color = pow(color, vec3(0.4545454545));
    
    return color;
}

bool intersect_plane(in Plane plane, in vec3 origin, in vec3 rayDirection, out float hitDistance, out vec3 Phit) 
{ 
    float denom = dot(plane.normal, rayDirection);
    if (denom < 1e-6)
    {
        vec3 p0l0 = plane.position - origin;

        hitDistance = dot(p0l0, plane.normal) / denom;

        if(hitDistance >= 0.0)
        {
            Phit = origin + rayDirection * hitDistance;

            return true;
        }
        
        return (hitDistance >= 0.0); 
    }

    return false;
}

bool intersect_sphere(in vec3 origin, in vec3 direction, in Sphere sphere, out float dist, out vec3 surfaceNormal, out vec3 Phit)
{
    vec3 L = origin - sphere.position;
    
    OUT float a = dot(direction, direction);
    OUT float b = 2.0 * dot(direction, L);
    OUT float c = dot(L, L) - pow(sphere.radius, 2.0);
    
    OUT float t0;
    OUT float t1;
    
    if (solve_quadratic(a, b, c, t0, t1))
    {        
        if (t0 > t1) 
        {
        	float temp = t0;
            t0 = t1;
            t1 = temp;
        } 
 
        if (t0 < 0.0)
        { 
            t0 = t1;
            if (t0 < 0.0) return false;
        }  
             
        dist = t0;
       
        Phit = origin + dist * direction;
        surfaceNormal = normalize(Phit - sphere.position);               
        
        return true;
    }  
     
    return false;
}

void calculateShadow(in vec3 pHit, out vec3 finalColor, in float ambient, in int type, in int index)
{
    vec3 shadowSurfaceNormal;
    vec3 shadowRay = pointlights[0].position - pHit;
    vec3 shadowRayDirection = normalize(shadowRay);
    float distanceToLight = sqrt(dot(shadowRay, shadowRay));
    vec3 shadowPhit;
    
    float dist;
    
    for(int i = 0; i < planes.length(); ++i)
	{
 		if (type == TYPE_PLANE && index == i)
        {
            continue;
        }
        
        if (intersect_plane(planes[i], pHit, shadowRay, dist, shadowPhit))
        {    
            if (dist < distanceToLight)
            {                
 				finalColor *= 0.25 * ambient;
            }
        }
    }
    
    for(int i = 0; i < spheres.length(); ++i)
	{
        if (type == TYPE_SPHERE && index == i)
        {
            continue;  
        }
    
        if (intersect_sphere(pHit, shadowRay, spheres[i], dist, shadowSurfaceNormal, shadowPhit))
        {
            if (dist > 0.0 && dist < distanceToLight)
            {
            	finalColor *= 0.25 * ambient;
            }
        }
    }
}

vec3 get_reflection(in vec3 direction, in vec3 surface_normal)
{
    return reflect(direction, surface_normal);
}

vec3 get_refraction(in vec3 direction, in vec3 surface_normal)
{
    float eta1 = 1.0;
    float eta2 = 1.9;
    float eta = eta1 / eta2;
    
    float c1 = dot(surface_normal, direction);
    if (c1 < 0.0)
    {
        c1 = -c1;
    }
    else
    {
        surface_normal = -surface_normal;
        eta = 1.0 / eta;
    }    
    
   	float theta = acos(c1);
    
    float k = 1.0 - eta * eta * sin(theta) * sin(theta);
    if (k < 0.0) 
        return vec3(0.0);
    
    float c2 = sqrt(k);
    
    vec3 ret = eta * direction + surface_normal * (eta * c1 - c2);
    
    return ret; 
}

vec3 create_ray(in vec3 origin, in vec3 direction)
{
    vec3 result = vec3(0.0);

    int previous_type = -1;
    int previous_index = -1;

    OUT vec3 hit = origin;
    vec3 bounce_pass_hit;
    
    bool has_intersected = false;
    
    for(int bounce = 0; bounce < bounces; ++bounce)
    {
        OUT vec3 surface_normal;

        float dist = 1e+10;
        OUT float object_hit_distance = dist;

        int current_type = -1;
        int current_index = -1;

        vec3 bounce_pass_color;
        
        for(int i = 0; i < planes.length(); ++i)
        {
            if(previous_type == TYPE_PLANE && previous_index == i)
            {
                continue;
            }

            if(intersect_plane(planes[i], origin, direction, object_hit_distance, hit))
            {
                if(object_hit_distance < dist)
                {
                    //result = planes[i].color;
                    dist = object_hit_distance;

                    switch(planes[i].type)
                    {
                        case PLANE_CHECKERS:
                        {
                            if(mod(floor(hit.x) - floor(hit.z), 2.0) == 0.0)
                            {
                                bounce_pass_color = get_color(direction, hit, planes[i].color, pointlights[0], planes[i].normal, planes[i].material);
                            }
                            else
                            {
                                bounce_pass_color = get_color(direction, hit, vec3(1.0), pointlights[0], planes[i].normal, planes[i].material);
                            }
                        }
                        break;

                        case PLANE_CLOUDS:
                        {

                        }
                        break;

                        default:
                        {} break;
                    }

                    surface_normal = planes[i].normal;

                    calculateShadow(hit, bounce_pass_color, planes[i].material.ambience, TYPE_PLANE, i);

                    current_type = TYPE_PLANE;
                    current_index = i;

                    bounce_pass_hit = hit;
                }
                
                has_intersected = true;
            }
        }
        
        for(int i = 0; i < spheres.length(); ++i)
        {
            if(previous_type == TYPE_SPHERE && previous_index == i)
            {
                continue;
            }

            if(intersect_sphere(origin, direction, spheres[i], object_hit_distance, surface_normal, hit))
            {
                if(object_hit_distance < dist)
                {
                    //result = spheres[i].color;
                    dist = object_hit_distance;
                    bounce_pass_color = get_color(direction, hit, spheres[i].color, pointlights[0], surface_normal, spheres[i].material);
                
                    calculateShadow(hit, bounce_pass_color, spheres[i].material.ambience, TYPE_SPHERE, i);

                    current_type = TYPE_SPHERE;
                    current_index = i;

                    bounce_pass_hit = hit;
                }
                
                has_intersected = true;
            }
        }
        
        Material mat = get_material(previous_type, previous_index);
        
        if(bounce == 0)
        {
            result += bounce_pass_color;
        }
        else
        {
            result += mat.specular * bounce_pass_color;
        }

        if(current_type < 0)
        {
            break;
        }

        origin = bounce_pass_hit/* + surface_normal * 1e-3*/;

        /*if(mat.reflective)
        {*/
            direction = get_reflection(direction, surface_normal);
        /*}*/

        /*if(mat.refractive)
        {*/
            //direction = get_refraction(direction, surface_normal);
        /*}*/

        previous_type = current_type;
        previous_index = current_index;
    }
    
    if(has_intersected)
    {
        return result;
    }
    else
    {
        /*vec3 top_color = vec3(0.8, 0.4, 1.0);
        vec3 bot_color = vec3(0.85, 0.9, 1.0);
        return mix(bot_color, top_color, direction.y);*/
        vec3 top_color = vec3(0.8, 0.4, 1.0);
        vec3 bot_color = vec3(0.85, 0.9, 1.0);
        vec3 m = mix(bot_color, top_color, direction.y);
        
        float f = fbm(vec2(direction.x - origin.x, direction.y * 0.5 - origin.y * 0.5));
        vec3 c = m + f * 0.25;
        return c;
    }
}

/* ---------------------------- */

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    setup();

    float aspect_ratio = iResolution.x / iResolution.y;
    vec2 uv = fragCoord/iResolution.xy - 0.5;
    uv.x *= aspect_ratio;

    camera.position = vec3(sin(iTime) * 1.0, 0.0, -cos(iTime) * 1.0);
    
    vec3 ray_origin = camera.position;

    camera.forward = normalize(camera.look_at - ray_origin);
    camera.right = cross(vec3(0.0, 1.0, 0.0), camera.forward);
    camera.up = cross(camera.forward, camera.right);

    vec3 center = ray_origin + camera.forward * camera.zoom;
    
#ifdef ANTIALIASING
    vec3 intersection1 = center + (uv.x - 0.001) * camera.right + uv.y * camera.up;
    vec3 intersection2 = center + (uv.x + 0.001) * camera.right + uv.y * camera.up;
    vec3 intersection3 = center + uv.x * camera.right + (uv.y - 0.001) * camera.up;
    vec3 intersection4 = center + uv.x * camera.right + (uv.y + 0.001) * camera.up;

    vec3 ray_direction1 = normalize(intersection1 - ray_origin);
    vec3 ray_direction2 = normalize(intersection2 - ray_origin);
    vec3 ray_direction3 = normalize(intersection3 - ray_origin);
    vec3 ray_direction4 = normalize(intersection4 - ray_origin);
    
    vec3 color1 = create_ray(ray_origin, ray_direction1);
    vec3 color2 = create_ray(ray_origin, ray_direction2);
    vec3 color3 = create_ray(ray_origin, ray_direction3);
    vec3 color4 = create_ray(ray_origin, ray_direction4);
    
    vec3 color = (color1 + color2 + color3 + color4) / 4.0;
#else
    vec3 intersection = center + (uv.x) * camera.right + uv.y * camera.up;
    
    vec3 ray_direction = normalize(intersection - ray_origin);
    
    vec3 color = create_ray(ray_origin, ray_direction);
#endif

    fragColor = vec4(color, 1.0);
}

// -----------------------------------------------------------------------------------