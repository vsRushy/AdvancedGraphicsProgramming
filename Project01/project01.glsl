#define OUT

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

Camera camera = Camera(vec3(0.0, 0.0, -1.0), vec3(0.0),
                       vec3(0.0), vec3(0.0), vec3(0.0),
                       2.0);

const Material material01 = Material(0.5, 0.4, 70.0, 0.6, 1.0);
const Material material02 = Material(0.4, 0.2, 120.0, 0.3, 1.0);

Plane plane01 = Plane(vec3(0.0, -0.2, 0.0), vec3(0.0, 1.0, 0.0), vec3(0.0, 0.2, 0.8), material01);

Sphere sphere01 = Sphere(vec3(0.1, 0.0, 0.0), 0.07, vec3(1.0, 0.5, 0.3), material01);
Sphere sphere02 = Sphere(vec3(-0.1, 0.0, 0.0), 0.09, vec3(0.5, 0.3, 0.5), material02);

PointLight pointlight01 = PointLight(vec3(0.0, 0.2, -0.1), vec3(1.0, 1.0, 1.0), 10.0);

Plane planes[1];
Sphere spheres[2];
PointLight pointlights[1];

const int bounces = 2;

/* ---------------------------- */

void setup()
{
    planes[0] = plane01;

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

        if(hitDistance >= 0.0)
        {
            Phit = origin + rayDirection * hitDistance;

            return true;
        }
        
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
 		if (type == 0 && index == i)
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
        if (type == 1 && index == i)
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
            if(previous_type == 0 && previous_index == i)
            {
                continue;
            }

            if(intersect_plane(planes[i], origin, direction, object_hit_distance, hit))
            {
                if(object_hit_distance < dist)
                {
                    //result = planes[i].color;
                    dist = object_hit_distance;
                    bounce_pass_color = get_color(direction, hit, planes[i].color, pointlights[0], planes[i].normal, planes[i].material);
                    
                    bounce_pass_color += vec3(mod(floor(hit.x) - floor(hit.z), 2.0) == 0.0) + vec3(0.15);
                    
                    surface_normal = planes[i].normal;

                    calculateShadow(hit, bounce_pass_color, planes[i].material.ambience, 0, i);

                    current_type = 0;
                    current_index = i;

                    bounce_pass_hit = hit;
                }
                
                has_intersected = true;
            }
        }
        
        for(int i = 0; i < spheres.length(); ++i)
        {
            if(previous_type == 1 && previous_index == i)
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
                
                    calculateShadow(hit, bounce_pass_color, spheres[i].material.ambience, 1, i);

                    current_type = 1;
                    current_index = i;

                    bounce_pass_hit = hit;
                }
                
                has_intersected = true;
            }
        }
        
        if(bounce == 0)
        {
            result = bounce_pass_color;
        }
        else
        {
            result += bounce_pass_color;
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
    vec3 intersection = center + uv.x * camera.right + uv.y * camera.up;

    vec3 ray_direction = normalize(intersection - ray_origin);
    
    vec3 color = create_ray(ray_origin, ray_direction);

    fragColor = vec4(color, 1.0);
}