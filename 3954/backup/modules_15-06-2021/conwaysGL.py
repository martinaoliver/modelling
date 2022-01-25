import pickle
from class_circuit_eq import *

import os
import sys
modulepath = os.path.expanduser('~/Documents/modelling/6eq/modules')  # os.path.expanduser(path) : return the argument with an initial component of ~ or ~user replaced by that user’s home directory.
sys.path.append(modulepath)
path = os.path.expanduser('~/Documents/modelling/6eq/parameter_space_search') #path of project folder: folder where code, results and parameterfiles are found.
from sympy import *

import pickle
from class_circuit_eq import *
import itertools

import numpy as np
from cellular_automaton import CellularAutomaton, MooreNeighborhood, CAWindow, EdgeRule
import random
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from cellular_automaton import CellularAutomaton, MooreNeighborhood, CAWindow, EdgeRule

ALIVE = [1.0]
DEAD = [0]


class ConwaysCA(CellularAutomaton):
    """ Cellular automaton with the evolution rules of conways game of life """

    def __init__(self):
        super().__init__(dimension=[100, 100],
                         neighborhood=MooreNeighborhood(EdgeRule.FIRST_AND_LAST_CELL_OF_DIMENSION_ARE_NEIGHBORS))

    def init_cell_state(self, __):  # pylint: disable=no-self-use
        rand = random.randrange(0, 16, 1)
        init = max(.0, float(rand - 14))
        return [init]

    def evolve_rule(self, last_cell_state, neighbors_last_states):
        new_cell_state = last_cell_state
        alive_neighbours = self.__count_alive_neighbours(neighbors_last_states)
        if last_cell_state == DEAD and alive_neighbours == 3:
            new_cell_state = ALIVE
        if last_cell_state == ALIVE and alive_neighbours < 2:
            new_cell_state = DEAD
        if last_cell_state == ALIVE and 1 < alive_neighbours < 4:
            new_cell_state = ALIVE
        if last_cell_state == ALIVE and alive_neighbours > 3:
            new_cell_state = DEAD
        return new_cell_state

    @staticmethod
    def __count_alive_neighbours(neighbours):
        alive_neighbors = []
        for n in neighbours:
            if n == ALIVE:
                alive_neighbors.append(1)
        return len(alive_neighbors)


if __name__ == "__main__":
    CAWindow(cellular_automaton=ConwaysCA(),
             window_size=(1000, 830)).run(evolutions_per_second=40)
