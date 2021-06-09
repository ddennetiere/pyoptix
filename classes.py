# coding: utf-8
import numpy as np


class Beamline(object):
    def __init__(self):
        super().__init__()
        self.elements = []
        self.chains = []

    def align(self):
        pass

    def chain(self):
        for element in self.elements:
            print(element)
            if element.next is None:
                new_chain = True
                for chain in self.chains:
                    if chain[0] == element:
                        new_chain = False
                if new_chain:
                    self.chains.append([element])
                    print("new_chain")

        def get_previous(an_element):
            for element in self.elements:
                if element.next == an_element.name:
                    return element

        for chain in self.chains:
            previous = 1
            while previous is not None:
                    previous = get_previous(chain[0])
                    if previous:
                        chain.insert(0, previous)

    def add_element(self, new_element):
        if new_element not in self.elements:
            self.elements.append(new_element)
        else:
            print("Warning : element already in beamline")
        for element in self.elements:
            if element.name == new_element.next:
                element.previous = new_element.name


class OpticalElement(object):
    def __init__(self, name="", phi=0, psi=0, theta=0, d_phi=0, d_psi=0, d_theta=0, x=0, y=0, z=0,
                 d_x=0, d_y=0, d_z=0, next=None, previous=None, distance_from_previous=0):
        super().__init__()
        self.name = name
        self.phi = phi
        self.psi = psi
        self.theta = theta
        self.d_phi = d_phi
        self.d_psi = d_psi
        self.d_theta = d_theta
        self.x = x
        self.y = y
        self.z = z
        self.d_x = d_x
        self.d_y = d_y
        self.d_z = d_z
        self.next = next
        self.previous = previous
        self.distance_from_previous = distance_from_previous

    def __repr__(self):
        description = f"Element {self.name} de classe {self.__class__}"
        description += f"\n\t at {self.distance_from_previous} m from {self.previous}"
        description += f"\n\t pointing to {self.next}"
        return description
