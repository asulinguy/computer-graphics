from helper_classes import *
import matplotlib.pyplot as plt

def render_scene(camera, ambient, lights, objects, screen_size, max_depth):
    width, height = screen_size
    ratio = float(width) / height
    screen = (-1, 1 / ratio, 1, -1 / ratio)  # left, top, right, bottom
    image = np.zeros((height, width, 3))

    for i, y in enumerate(np.linspace(screen[1], screen[3], height)):
        for j, x in enumerate(np.linspace(screen[0], screen[2], width)):
            pixel = np.array([x, y, 0])
            origin = camera
            direction = normalize(pixel - origin)
            ray = Ray(origin, direction)
            color = np.zeros(3)
            intersection = find_intersection(ray, objects)

            if intersection is not None:
                intersection_point, nearest_object = intersection
                intersection_point += get_point_bias(nearest_object, intersection_point, ray)
                color = get_color(camera, ambient, lights, nearest_object, ray, intersection_point, 0, max_depth, objects)

                image[i, j] = np.clip(color, 0, 1)

    return image

def get_color(camera, ambient, lights, nearest_object, ray, intersection_point, depth, max_depth, objects):
    color = calc_ambient(nearest_object, ambient)
    for light in lights:
        shadow_light = shadow_coefficient(light, intersection_point, objects)
        diffuse_color = calc_diffuse(nearest_object, light, intersection_point)
        specular_color = calc_specular(camera, light, intersection_point, nearest_object)

        color = color + (diffuse_color + specular_color) * shadow_light * light.get_intensity(intersection_point)

    depth += 1
    if depth >= max_depth:
        return color

    r_ray = Ray(intersection_point, reflected(ray.direction, object_norm(nearest_object, intersection_point)))

    intersection = find_intersection(r_ray, objects)
    if intersection is None:
        return color

    r_hit, hit_object = intersection
    r_hit += get_point_bias(hit_object, r_hit, r_ray)
    color += nearest_object.reflection * get_color(camera, ambient, lights, hit_object, r_ray, r_hit, depth + 1, max_depth, objects)

    return color

def find_intersection(ray: Ray, objects):
    intersection = ray.nearest_intersected_object(objects)
    if intersection is None:
        return None

    min_distance, nearest_object = intersection
    intersection_point = ray.origin + (min_distance * ray.direction)
    return intersection_point, nearest_object

def shadow_coefficient(light: LightSource, point, objects):
    light_ray = light.get_light_ray(point)
    intersection = light_ray.nearest_intersected_object(objects)

    if intersection is None:
        return 1

    distance, _ = intersection
    if distance < light.get_distance_from_light(point):
        return 0

    return 1

def calc_ambient(obj, ambient):
    return obj.ambient * ambient

def calc_diffuse(obj, light, intersection_point):
    normal = object_norm(obj, intersection_point, None)  
    light_ray = light.get_light_ray(intersection_point)
    return np.multiply(obj.diffuse, np.dot(normal, light_ray.direction))

def calc_specular(camera, light, point, object):
    N = object_norm(object, point)
    L = -1 * light.get_light_ray(point).direction
    V_hat = normalize(camera - point)
    R_hat = reflected(L, N)

    return np.multiply(object.specular, (np.dot(V_hat, R_hat)) ** (object.shininess))

def object_norm(obj, intersection_point, ray=None):
    """
    Returns the normal at the intersection point based on the object type.

    Args:
        obj (Object3D): The intersected object.
        intersection_point: The intersection point on the object.
        ray (Ray): The ray that intersects the object (optional, used for certain types).

    Returns:
        np.array: The normal vector at the intersection point.
    """
    if isinstance(obj, Plane):
        return obj.normal
    elif isinstance(obj, Sphere):
        return normalize(intersection_point - obj.center)
    elif isinstance(obj, Triangle):
        return obj.compute_normal()
    elif isinstance(obj, Pyramid):
        if hasattr(obj, 'last_intersected_triangle') and obj.last_intersected_triangle is not None:
            return obj.last_intersected_triangle.compute_normal()
    else:
        print(type(obj))
        print(obj)
        raise NotImplementedError("Normal retrieval for this object type is not implemented.")

def get_point_bias(intersected_object, intersection_point, ray=None, epsilon=0.01):
    """
    Shifts an intersection point slightly away from the surface.

    Args:
        intersected_object (Object3D): The object the intersection occurs with.
        intersection_point (np.array): The point of intersection.
        ray (Ray): The ray that intersects the object.
        epsilon (float): The amount to shift the intersection point along the normal vector.

    Returns:
        np.array: The adjusted intersection point.
    """
    return epsilon * object_norm(intersected_object, intersection_point, ray)


# Write your own objects and lights
def your_own_scene():
    camera = np.array([0,0,1])
    lights = []
    objects = []
    return camera, lights, objects
