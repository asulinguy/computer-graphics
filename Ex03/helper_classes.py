import numpy as np


# This function gets a vector and returns its normalized form.
def normalize(vector):
    return vector / np.linalg.norm(vector)


# This function gets a vector and the normal of the surface it hit
# This function returns the vector that reflects from the surface
def reflected(vector, axis):
    return vector - 2 * np.dot(vector, axis) * axis

## Lights


class LightSource:
    def __init__(self, intensity):
        self.intensity = intensity

    # This function determines if the light is blocked by any object between the light and the intersection
    def is_light_blocked(self, intersection, objects, nearest_object):
        objects_without_nearest = objects
        ray_to_light = self.get_light_ray(intersection)
        distance_to_light = self.get_distance_from_light(intersection)
        _, min_distance = ray_to_light.nearest_intersected_object(objects_without_nearest)
        return min_distance < distance_to_light

class DirectionalLight(LightSource):

    def __init__(self, intensity, direction):
        super().__init__(intensity)
        self.direction = normalize(np.array(direction)) 

    # This function returns the ray that goes from the light source to a point
    def get_light_ray(self,intersection_point):
        # For a directional light, the direction is constant and does not depend on the intersection point
        return Ray(intersection_point, -self.direction)

    # This function returns the distance from a point to the light source
    def get_distance_from_light(self, intersection):
        return np.inf

    # This function returns the light intensity at a point
    def get_intensity(self, intersection):
        return self.intensity


class PointLight(LightSource):
    def __init__(self, intensity, position, kc, kl, kq):
        super().__init__(intensity)
        self.position = np.array(position)
        self.kc = kc
        self.kl = kl
        self.kq = kq

    # This function returns the ray that goes from a point to the light source
    def get_light_ray(self,intersection):
        return Ray(intersection, normalize(self.position - intersection))

    # This function returns the distance from a point to the light source
    def get_distance_from_light(self,intersection):
        return np.linalg.norm(intersection - self.position)
    
    

    # This function returns the light intensity at a point
    def get_intensity(self, intersection):
        d = self.get_distance_from_light(intersection)
        return self.intensity / (self.kc + self.kl*d + self.kq * (d**2))


class SpotLight(LightSource):
    def __init__(self, intensity, position, direction, kc, kl, kq):
        super().__init__(intensity)
        # TODO

    # This function returns the ray that goes from the light source to a point
    def get_light_ray(self, intersection):
        #TODO
        pass

    def get_distance_from_light(self, intersection):
        #TODO
        pass

    def get_intensity(self, intersection):
        #TODO
        pass


class Ray:
    def __init__(self, origin, direction):
        self.origin = origin
        self.direction = direction

    # The function is getting the collection of objects in the scene and looks for the one with minimum distance.
    # The function should return the nearest object and its distance (in two different arguments)
    def nearest_intersected_object(self, objects):
        nearest_object = None
        min_distance = np.inf

        for obj in objects:
            intersection_data = obj.intersect(self)

            if intersection_data is not None:
                distance, _ = intersection_data

                if distance > 0 and distance < min_distance:
                    min_distance = distance
                    nearest_object = obj

        return nearest_object,min_distance


class Object3D:
    def set_material(self, ambient, diffuse, specular, shininess, reflection):
        self.ambient = ambient
        self.diffuse = diffuse
        self.specular = specular
        self.shininess = shininess
        self.reflection = reflection


class Plane(Object3D):
    def __init__(self, normal, point):
        self.normal = np.array(normal)
        self.point = np.array(point)

    def compute_normal(point,normal):
        return np.array(normal)
        
    def intersect(self, ray: Ray):
        v = self.point - ray.origin
        t = np.dot(v, self.normal) / (np.dot(self.normal, ray.direction) + 1e-6)
        if t > 0:
            return t, self
        else:
            return None


class Triangle(Object3D):
    """
        C
        /\
       /  \
    A /____\ B

    The fornt face of the triangle is A -> B -> C.
    
    """
    def __init__(self, a, b, c):
        self.a = np.array(a)
        self.b = np.array(b)
        self.c = np.array(c)
        self.normal = self.compute_normal()

    # computes normal to the trainagle surface. Pay attention to its direction!
    def compute_normal(self):
        edge1 = self.b - self.a
        edge2 = self.c - self.a
        normal = np.cross(edge1, edge2)
        return normalize(normal)

    def intersect(self, ray: Ray):
        # Ray-triangle intersection algorithm (Möller-Trumbore)
        vertex0 = self.a
        vertex1 = self.b
        vertex2 = self.c
        edge1 = vertex1 - vertex0
        edge2 = vertex2 - vertex0
        h = np.cross(ray.direction, edge2)
        a = np.dot(edge1, h)
        if -1e-6 < a < 1e-6:
            return None  # Ray is parallel to the triangle

        f = 1.0 / a
        s = ray.origin - vertex0
        u = f * np.dot(s, h)
        if u < 0.0 or u > 1.0:
            return None
        q = np.cross(s, edge1)
        v = f * np.dot(ray.direction, q)
        if v < 0.0 or u + v > 1.0:
            return None

        t = f * np.dot(edge2, q)
        if t > 1e-6:  # Check for a real intersection
            return t, self
        else:
            return None
        

class Pyramid(Object3D):
    """     
            D
            /\*\
           /==\**\
         /======\***\
       /==========\***\
     /==============\****\
   /==================\*****\
A /&&&&&&&&&&&&&&&&&&&&\ B &&&/ C
   \==================/****/
     \==============/****/
       \==========/****/
         \======/***/
           \==/**/
            \/*/
             E 
    
    Similar to Traingle, every from face of the diamond's faces are:
        A -> B -> D
        B -> C -> D
        A -> C -> B
        E -> B -> A
        E -> C -> B
        C -> E -> A
    """
    def __init__(self, v_list):
        self.v_list = v_list
        self.triangle_list = self.create_triangle_list()

    def create_triangle_list(self):
        l = []
        t_idx = [
                [0,1,3],
                [1,2,3],
                [0,3,2],
                 [4,1,0],
                 [4,2,1],
                 [2,4,0]]
        # TODO
        return l

    def apply_materials_to_triangles(self):
        # TODO
        pass

    def intersect(self, ray: Ray):
        # TODO
        pass

class Sphere(Object3D):
    def __init__(self, center, radius: float):
        self.center = center
        self.radius = radius

    def intersect(self, ray: Ray):
        #TODO
        pass




def find_normal(obj, intersection_point = None):
    if obj.__class__.__name__ == "Plane":
        return obj.normal
    if obj.__class__.__name__ == "Triangle":
        return obj.normal
    else:
        raise KeyError("Normal not implemented for this object type")



