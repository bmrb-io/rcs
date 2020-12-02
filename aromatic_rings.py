from atoms import Atom
import numpy as np


class Aromatic_ring:

    def __init__(self, atoms_list):

        self.atoms_list = atoms_list
        self.center_position = self.find_center() 

    def find_distance(self, other_atom):

        dist_sum = 0
        for atom in self.atoms_list:
            dist_vec = [
                a - b for a, b in zip(other_atom.position, atom.position)
            ]
        return dist_sum / 6

    def find_center(self):

        center_x, center_y, center_z = 0, 0, 0
        for atom in self.atoms_list:
            center_x += atom.position[0]
            center_y += atom.position[1]
            center_z += atom.position[2]
        center_position = [center_x / 6, center_y / 6, center_z / 6]

        return center_position

    def find_angle(self, other_atom):

        pos0 = self.atoms_list[0].position
        pos1 = self.atoms_list[1].position

        vec0 = [pos0[i] - self.center_position[i] for i in range(3)]
        vec1 = [pos1[i] - self.center_position[i] for i in range(3)]
        vec_h = [
            other_atom.position[i] - self.center_position[i] for i in range(3)
        ]

        n = np.cross(vec0, vec1)

        phi = np.arccos(
            np.dot(n, vec_h) / (np.linalg.norm(n) * np.linalg.norm(vec_h))
        )
        phi = phi * 180 / np.pi

        return phi
