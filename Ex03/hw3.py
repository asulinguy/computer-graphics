from helper_classes import *
import matplotlib.pyplot as plt

def render_scene(camera, ambient, lights, objects, screen_size, max_depth):
    I_A = ambient
    width, height = screen_size
    ratio = float(width) / height
    screen = (-1, 1 / ratio, 1, -1 / ratio)  # left, top, right, bottom

    image = np.zeros((height, width, 3))

    for i, y in enumerate(np.linspace(screen[1], screen[3], height)):
        for j, x in enumerate(np.linspace(screen[0], screen[2], width)):
            # screen is on origin
            pixel = np.array([x, y, 0])
            origin = camera
            direction = normalize(pixel - origin)
            ray = Ray(origin, direction)
            nearest_object, _ = ray.nearest_intersected_object(objects)
            if nearest_object is None:
                color = np.zeros(3)
            color = calc_color(ray, nearest_object, objects, lights, I_A, max_depth)
            # We clip the values between 0 and 1 so all pixel values will make sense.
            image[i, j] = np.clip(color,0,1)
                
            

    return image


# Write your own objects and lights
# TODO
def your_own_scene():
    camera = np.array([0,0,1])
    lights = []
    objects = []
    return camera, lights, objects



def calc_color(ray, nearest_object, objects, lights, I_A, max_depth, level=1, epsilon=0.01):
    color = np.zeros(3)
    color +=  I_A * nearest_object.ambient
    intersection_point = ray.origin + ray.direction * ray.nearest_intersected_object(objects)[1] 
    normal = find_normal(nearest_object, intersection_point)
    intersection_point += normal * epsilon
    for light in lights:
        # Check if the light is blocked and if not, calculate the color
        if not light.is_light_blocked(intersection_point, objects, nearest_object):
            diffuse_color = calc_diffuse(nearest_object, light, intersection_point, normal)
            specular_color = calc_specular(ray, light, intersection_point, nearest_object, normal)
            color = color + (diffuse_color + specular_color) * light.get_intensity(intersection_point)
    
    # Reflection and refraction
    if level < max_depth:
        # Reflection
        reflected_ray_direction = reflected(ray.direction, normal)
        reflected_ray = Ray(intersection_point, reflected_ray_direction)
        objects_without_nearest = [obj for obj in objects if obj != nearest_object]
        hit_obj,_ = reflected_ray.nearest_intersected_object(objects_without_nearest)
        if hit_obj is not None:
            color += nearest_object.reflection * calc_color(reflected_ray, hit_obj, objects, lights, I_A, max_depth, level + 1)

        # Refraction
        
        
    return color





def calc_diffuse(obj, light, intersection_point, normal):
    light_ray = light.get_light_ray(intersection_point)
    return np.multiply(obj.diffuse, np.dot(normal, light_ray.direction))

def calc_specular(camera, light, point, object, normal):
    N = normal
    L = -1 * light.get_light_ray(point).direction
    V_hat = normalize(camera.origin - point)
    R_hat = reflected(L, N)
    return np.multiply(object.specular, (np.dot(V_hat, R_hat)) ** (object.shininess))