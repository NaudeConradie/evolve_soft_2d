##  Functions used with the representative grids

#   Imports
import numpy
import itertools

from scipy.ndimage import measurements

from evolve_soft_2d.log import m_log
from evolve_soft_2d.unit import inspect

################################################################################

def create_grid(
    x: int,
    y: int,
    d: int,
    ) -> list:
    """Create a grid of digits

    Parameters
    ----------
    x : int
        The number of elements in the x-direction
    y : int
        The number of elements in the y-direction
    d : int
        The digit to create the grid of

    Returns
    -------
    list
        The representative grid
    """
    grid = [[d]*(x) for i in range(y)]

    if y == 1:

        grid = grid[0]

    return grid
 
################################################################################

def find_cluster(grid: list) -> [bool, list]:
    """Find all clusters of elements using the representative grid

    Parameters
    ----------
    grid : list
        The representative grid of ones

    Returns
    -------
    [bool, list]
        True if free clusters were found, false otherwise
        The grid with clusters incrementally labelled
    """

    #   Initialisations
    s = [[0, 1, 0], [1, 1, 1], [0, 1, 0]]

    #   Identify the clusters in the grid
    grid_label, cluster = measurements.label(grid, structure = s)

    #   Check if more than one cluster is found
    if cluster > 1:

        #   Set flag
        found = True

        # if cluster == 2:
        #     m_log.warning("{} free element cluster found!".format(cluster - 1))
        # else:
        #     m_log.warning("{} free element clusters found!".format(cluster - 1))

    else:

        #   Set flag
        found = False

    return found, grid_label

################################################################################

def find_dia(
    template,
    grid: list,
    ) -> [list, list]:
    """Find diagonally connected elements

    Parameters
    ----------
    template : template
        The unit template parameters
    grid : list
        The representative grid

    Returns
    -------
    [list, list]
        The lists of diagonally connected element pairs
    """    

    #   Initialisations
    dia_pairs_pos = []
    dia_pairs_neg = []

    #   Prepare the grid for analysis
    grid_flat = list(reversed(grid))
    grid_flat = list(itertools.chain.from_iterable(grid_flat))

    #   Set the list of elements to be checked as all internal elements except the last row
    e_check = [i - 1 for i in template.e_internal]
    e_check = e_check[:-(template.x_e - 2*template.b)]

    #   Loop through the list of elements to be checked
    for i in e_check:

        #   Check if the current element exists
        if grid_flat[i] == 1:

            #   Check if the element only has a diagonal neighbour in the positive direction
            if grid_flat[i + 1] == 0 and grid_flat[i + template.x_e] == 0 and grid_flat[i + template.x_e + 1] == 1:

                #   Add the diagonally connected element pair coordinates to the list
                dia_pairs_pos.append([i + 1, i + template.x_e + 1 + 1])

        #   Check if the current element does not exist
        elif grid_flat[i] == 0:

            #   Check if the next element only has a diagonal neighbour in the negative direction
            if grid_flat[i + 1] == 1 and grid_flat[i + template.x_e] == 1 and grid_flat[i + template.x_e + 1] == 0:

                #   Add the diagonally connected element pair coordinates to the list
                dia_pairs_neg.append([i + 1 + 1, i + template.x_e + 1])

    return dia_pairs_pos, dia_pairs_neg

################################################################################

def rem_el_grid(
    template,
    rem: list,
    ) -> numpy.array:
    """Removes elements from the representative grid

    Parameters
    ----------
    template : template
        The unit template parameters
    rem : list
        The element IDs of the elements to be removed

    Returns
    -------
    numpy.array
        The grid with zeros in the places of the removed elements
    """

    #   Initialisations
    grid_rem = numpy.array(template.grid)

    #   Loop through the number of elements to be removed
    for i in range(0, len(rem)):

        #   Remove the element from the grid
        grid_rem[template.x_e - (rem[i] - 1)//template.x_e - 1][rem[i]%template.x_e - 1] = 0

    return grid_rem

################################################################################

def rem_el_free_grid(
    template,
    grid_label: list,
    ) -> [numpy.array, list]:
    """Remove free element clusters from the representative grid

    Parameters
    ----------
    template : template
        The unit template parameters
    grid_label : list
        Representative grid with clusters incrementally labelled

    Returns
    -------
    [numpy.array, list]
        The grid with zeros in place of the removed elements
        The list of removed elements
    """
    
    #   Initialisations
    rem_i = 1
    rem = []

    grid_rem = numpy.array(template.grid)

    #   Loop through the elements in the x-direction
    for i in range(0, template.x_e):

        #   Loop through the elements in the y-direction
        for j in range(0, template.y_e):

            #   Check if the labelled grid has an element numbered greater than 1
            if grid_label[template.x_e - i - 1][j] > 1:

                #   Remove the element from the grid
                grid_rem[template.x_e - (rem_i - 1)//template.x_e - 1, rem_i%template.x_e - 1] = 0

                #   Add the index of the element to the list of removed elements
                rem.append(rem_i)

            elif grid_label[template.x_e - i - 1][j] == 0:

                grid_rem[template.x_e - (rem_i - 1)//template.x_e - 1, rem_i%template.x_e - 1] = 0

            #   Increment the removed element counter
            rem_i = rem_i + 1

    return grid_rem, rem