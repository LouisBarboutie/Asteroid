"""
Author: Konrad Barboutie
Date: 2025/10/31
Description: vector class definition for solar system 3d visualisation
Url: https://thepythoncodingbook.com/2021/12/11/simulating-3d-solar-system-python-matplotlib/
"""
import math

class Vector:
    def __init__(self, x=0.0,y=0.0,z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"Vector({self.x}, {self.y}, {self.z})"

    def __str__(self):
        return f"{self.x}i + {self.y}j + {self.z}k"

    def __getitem__(self, item):
        if item == 0:
            return self.x
        elif item == 1:
            return self.y
        elif item == 2:
            return self.z
        else:
            raise IndexError("There are only three elements in the vector ")

    def __add__(self, other):
        return Vector(
            self.x + other.x,
            self.y + other.y,
            self.z + other.z
        )
    def __sub__(self, other):
        return Vector(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z
        )
    def __mul__(self, other):
        if isinstance(other, Vector):
            return (
                self.x * other.x
                + self.y * other.y
                + self.z * other.z
            )
        elif isinstance(other, int) or isinstance(other, float):
            return Vector(
                self.x * other,
                self.y * other,
                self.z * other
            )
        else:
            raise TypeError("Operand must be vector, int or float")
    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Vector(
                self.x / other,
                self.y / other,
                self.z / other,
            )
        else:
            raise TypeError("Operand must be int or float")
    def get_magnitude(self):
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
    def normalize(self):
        magnitude = self.get_magnitude()
        return Vector(
            self.x / magnitude,
            self.y / magnitude,
            self.z / magnitude,
        )
test =  Vector(3,6,9)#Vector(1,1,1)
print(test.get_magnitude())
print(test.normalize())
print(test.normalize().get_magnitude())
