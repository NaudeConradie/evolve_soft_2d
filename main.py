##  Main program

#   Imports
import importlib

from evolve_soft_2d import classes, file_paths, log, plotting, utility
from evolve_soft_2d.evolve import cppns, gen_alg, lsystems
from evolve_soft_2d.result import analyse, obtain
from evolve_soft_2d.unit import create, inspect, modify, rep_grid

from py_mentat import py_connect, py_disconnect

################################################################################

#   Main function

def main():

    #   Reload the modules
    importlib.reload(classes)
    importlib.reload(file_paths)
    importlib.reload(plotting)
    importlib.reload(utility)
    importlib.reload(cppns)
    importlib.reload(gen_alg)
    importlib.reload(lsystems)
    importlib.reload(analyse)
    importlib.reload(obtain)
    importlib.reload(create)
    importlib.reload(inspect)
    importlib.reload(modify)
    importlib.reload(rep_grid)

    #   Initialisations
    #   The template case identifier
    case = 1
    #   The number of elements in each axis direction
    x_e = 15
    y_e = 15
    #   The length of each side in mm
    e_s = 10
    #   The thickness of the unit boundary
    b = 3
    #   The number of increments per second to analyse
    n_steps = 5
    #   The text name of the table used for the applied displacement and load
    table_name = "ramp_input"
    #   The magnitude of the applied displacement in mm
    d_mag = y_e*e_s/2
    #   The magnitude of the applied internal pressure in MPa
    p_mag = 0.025

    #   The unit generation method
    g_meth = "r"
    #   The analysis method
    a_meth = "m"

    #   Genetic algorithm parameters
    gen = 15
    prob = [0.5, 0.1, 0.5]
    point = [1, 2, 2]

    #   Prepare the unit parameters
    temp = classes.template(case, x_e, y_e, e_s, b, classes.mold_star_15, n_steps, table_name, d_mag, p_mag)

    #   Create the template
    create.temp_create(temp)

    #   Check if the analysis method is the Monte Carlo method
    if a_meth == "m":
        analyse.monte_carlo(temp, g_meth)

    #   Check if the analysis method is the Genetic Algorithm method
    elif a_meth == "g":
        gen_alg.g_a(temp, gen, prob, point, g_meth)

    #   View the boundary conditions of the template
    inspect.view_bc()

    return

if __name__ == "__main__":

    py_connect("", 40007)

    main()
    
    py_disconnect