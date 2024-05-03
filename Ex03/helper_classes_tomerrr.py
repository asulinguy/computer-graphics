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


class DirectionalLight(LightSource):

    def __init__(self, intensity, direction):
        super().__init__(intensity)
        self.direction = direction
        self.intensity = intensity

    # This function returns the ray that goes from a point to the light source
    def get_light_ray(self,intersection_point):
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

        return min_distance, nearest_object



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

    def intersect(self, ray: Ray):
        v = self.point - ray.origin
        t = np.dot(v, self.normal) / (np.dot(self.normal, ray.direction) + 1e-6)
        if t > 0:
            return t, self
        else:
            return None
        
        # denominator = np.dot(self.normal, ray.direction)
        # if abs(denominator) < 1e-6:  # Check if ray is parallel to the plane
        #     return None

        # v = self.point - ray.origin
        # t = np.dot(v, self.normal) / denominator
        # if t > 0:  # Ensure the intersection is in front of the ray origin
        #     return t, self
        # else:
        #     return None


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

    # computes normal to the trianagle surface. Pay attention to its direction!
    def compute_normal(self):
        return normalize(np.cross(self.b - self.a, self.c - self.a))

    def intersect(self, ray):
        origin, direction = ray.origin, ray.direction
        edge1, edge2 = self.b - self.a, self.c - self.a
        p = np.cross(direction, edge2)
        det = np.dot(edge1, p)
        
        if abs(det) < 1e-8:
            return None  # The ray is parallel to the plane of the triangle

        inv_det = 1.0 / det
        t = origin - self.a
        u = np.dot(t, p) * inv_det
        if not (0 <= u <= 1):
            return None

        q = np.cross(t, edge1)
        v = np.dot(direction, q) * inv_det
        if v < 0 or u + v > 1:
            return None

        t = np.dot(edge2, q) * inv_det
        if t < 0:
            return None  # The triangle is behind the ray

        intersection_point = origin + t * direction
        return t, intersection_point
        

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
        for i in t_idx:
            l.append(Triangle(self.v_list[i[0]], self.v_list[i[1]], self.v_list[i[2]]))
        return l

    def apply_materials_to_triangles(self):
        for t in self.triangle_list:
            t.set_material(self.ambient, self.diffuse, self.specular, self.shininess, self.reflection)

    def intersect(self, ray: Ray):
        closest_intersection = None
        closest_triangle = None
        for triangle in self.triangle_list:
            result = triangle.intersect(ray)
            if result is not None:
                if closest_intersection is None or result[0] < closest_intersection[0]:
                    closest_intersection = result
                    closest_triangle = triangle
        if closest_intersection:
            self.last_intersected_triangle = closest_triangle
            return closest_intersection
        return None

class Sphere(Object3D):
    def __init__(self, center, radius: float):
        self.center = center
        self.radius = radius

    def intersect(self, ray: Ray):
        #TODO
        pass
