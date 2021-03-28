#define OUT

// Undefine for no anti-aliasing
#define ANTIALIASING

#define MAX_PLANES          1
#define MAX_SPHERES         3
#define MAX_POINTLIGHTS     1

#define TYPE_PLANE          0
#define TYPE_SPHERE         1

#define PLANE_CHECKERS      1
#define PLANE_CLOUDS        2

/* ---------------------------- */

struct Camera
{
    vec3 position;
    
    vec3 look_at;
    vec3 forward;
    vec3 right;
    vec3 up;

    float zoom;
};

struct Material 
{
    float diffuse;
    float specular;
    float shininess;
    float ambience;
    float reflection;
};

struct Plane 
{
    vec3 position;
    vec3 normal;
    vec3 color;

    int type;
    
    Material material;
};

struct Sphere
{
    vec3 position;
    float radius;
    vec3 color;
    
    Material material;
};

struct PointLight
{
    vec3 position;
    vec3 color; // TODO
    float intensity;
};

/* ---------------------------- */

Camera camera = Camera(vec3(0.0, 0.0, -2.0), vec3(0.0), vec3(0.0), vec3(0.0), vec3(0.0), 2.0);

const Material material01 = Material(1.0, 1.0, 100.0, 1.0, 1.0);
const Material material02 = Material(0.4, 1.0, 120.0, 0.3, 1.0);
const Material material03 = Material(0.2, 1.0, 100.0, 0.4, 1.0);
const Material material04 = Material(0.2, 1.0, 100.0, 0.4, 0.0);

Plane plane01 = Plane(vec3(0.0, -0.2, 0.0), vec3(0.0, 1.0, 0.0), vec3(0.0), PLANE_CHECKERS, material01);

Sphere sphere01 = Sphere(vec3(0.1, 0.0, 0.0), 0.07, vec3(1.0, 0.5, 0.3), material02);
Sphere sphere02 = Sphere(vec3(-0.1, 0.0, 0.0), 0.09, vec3(0.5, 0.3, 0.5), material03);
Sphere sphere03 = Sphere(vec3(0.5, 0.0, 0.0), 0.2, vec3(0.8, 0.6, 0.0), material04);

PointLight pointlight01 = PointLight(vec3(0.0, 0.2, -0.1), vec3(1.0, 1.0, 1.0), 10.0);

Plane planes[MAX_PLANES];
Sphere spheres[MAX_SPHERES];
PointLight pointlights[MAX_POINTLIGHTS];

const int bounces = 2;

/* ---------------------------- */

void setup()
{
    planes[0] = plane01;

    spheres[0] = sphere01;
    spheres[1] = sphere02;
    spheres[2] = sphere03;
    
    pointlights[0] = pointlight01;
}

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

float fresnel(in vec3 I, in vec3 N, in float ior) 
{
    float kr;

    float cosi = clamp(-1.0, 1.0, dot(I, N)); 
    
    float etai = 1.0;
    float etat = ior; 
    
    if (cosi > 0.0) 
    {
        float temp = etai;
        etai = etat;
        etat = temp;        
    } 
    
    float sint = etai / etat * sqrt(max(0.0, 1.0 - cosi * cosi)); 

    if (sint >= 1.)
    { 
        kr = 1.0; 
    } 
    else 
    { 
        float cost = sqrt(max(0.0, 1.0 - sint * sint)); 
        cosi = abs(cosi); 
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)); 
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)); 
        kr = (Rs * Rs + Rp * Rp) / 2.0; 
    }

    return kr;
}

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
    }
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
    
    //float f = fresnel(viewDir, surfaceNormal, coeff);
    //float f = 1. - dot(normalize(camera.position), normalize(surfaceNormal));
    float f = pow(1.0 - max(dot(surfaceNormal, -lightDir), 0.0), 10.0); 
    
    //vec3 color = ambient + diffuse + specular;
    vec3 color = ambient + mix(diffuse, vec3(1.0), specular) + f;
    
    color = pow(color, vec3(0.4545454545));
    
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
 				finalColor *= 2.0 * ambient;        
                //finalColor = vec3(0.0);
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
            	finalColor *= 2.0 * ambient;
                //finalColor = vec3(0.0);
            }
        }
    }
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
        
        if(bounce == 0)
        {
            result += bounce_pass_color;
        }
        else
        {
            Material mat = get_material(previous_type, previous_index);
            result += mat.specular * bounce_pass_color;
        }

        if(current_type < 0)
        {
            break;
        }

        origin = bounce_pass_hit/* + surface_normal * 1e-3*/;
        direction = reflect(direction, surface_normal);

        previous_type = current_type;
        previous_index = current_index;
    }
    
    if(has_intersected)
    {
        return result / float(bounces);
    }
    else
    {
        vec3 top_color = vec3(0.8, 0.4, 1.0);
        vec3 bot_color = vec3(0.85, 0.9, 1.0);
        return mix(bot_color, top_color, direction.y);
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