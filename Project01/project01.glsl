#define OUT

/* ---------------------------- */

struct Camera
{
    vec3 position;
    float focal_distance;
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

Camera camera = Camera(vec3(0.0, 0.0, -0.2), 0.6);

const Material material01 = Material(0.5, 0.4, 70.0, 0.6, 1.0);
const Material material02 = Material(0.4, 0.2, 120.0, 0.3, 1.0);

Plane plane = Plane(vec3(0.0, -0.2, 0.0), vec3(0.0, 1.0, 0.0), vec3(0.0, 0.2, 0.8), material01);

Sphere sphere01 = Sphere(vec3(0.1, 0.0, 0.0), 0.07, vec3(1.0, 0.5, 0.3), material01);
Sphere sphere02 = Sphere(vec3(-0.1, 0.0, 0.0), 0.09, vec3(0.5, 0.3, 0.5), material02);

PointLight pointlight01 = PointLight(vec3(0.0, 0.2, -0.1), vec3(1.0, 1.0, 1.0), 10.0);

Sphere spheres[2];
PointLight pointlights[1];

/* ---------------------------- */

void setup()
{
    spheres[0] = sphere01;
    spheres[1] = sphere02;
    
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
    
    vec3 color = ambient + diffuse + specular;
    //color = pow(color, vec3(0.4545454545));
    
    return color;
}

bool intersect_plane(in Plane plane, in vec3 origin, in vec3 rayDirection,
                     out float hitDistance, out vec3 Phit) 
{ 
    float denom = dot(plane.normal, rayDirection);
    if (denom < 1e-6)
    {
        vec3 p0l0 = plane.position - origin;
        hitDistance = dot(p0l0, plane.normal) / denom;
        Phit = origin + rayDirection * hitDistance;
        
        return (hitDistance >= 0.0); 
    }

    return false;
}

bool intersect_sphere(in vec3 origin, in vec3 direction, in Sphere sphere,
                      out float dist, out vec3 surfaceNormal, out vec3 Phit)
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

vec3 create_ray(in vec3 origin, in vec3 direction)
{
    vec3 result = vec3(0.2);

    float dist = 10000000.0;; 
    
    OUT vec3 hit;
    OUT float object_hit_distance = 0.0;
    
    OUT vec3 sphere_normal;
    
    if(intersect_plane(plane, origin, direction, object_hit_distance, hit))
    {
        if(object_hit_distance < dist)
        {
            //result = plane.color;
            dist = object_hit_distance;
            result = get_color(direction, hit, plane.color, pointlights[0], plane.normal, plane.material);
        }
    }
    
    for(int i = 0; i < spheres.length(); ++i)
    {
        if(intersect_sphere(origin, direction, spheres[i], object_hit_distance, sphere_normal, hit))
        {
            if(object_hit_distance < dist)
            {
                //result = spheres[i].color;
                dist = object_hit_distance;
                result = get_color(direction, hit, spheres[i].color, pointlights[0], sphere_normal, spheres[i].material);
            }
        }
    }
    
    return result;
}

/* ---------------------------- */

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    setup();

    float aspect_ratio = iResolution.x / iResolution.y;
    vec2 uv = fragCoord/iResolution.xy - 0.5;
    uv.x *= aspect_ratio;
    
    vec3 near_clip_plane_position = vec3(uv.x, uv.y, camera.position.z + camera.focal_distance);

    vec3 ray_origin = camera.position;
    vec3 ray_direction = normalize(near_clip_plane_position - ray_origin);
    
    vec3 color = create_ray(ray_origin, ray_direction);

    fragColor = vec4(color, 1.0);
}