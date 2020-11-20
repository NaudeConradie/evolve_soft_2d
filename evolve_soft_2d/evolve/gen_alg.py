##  Genetic Algorithms

#   Imports

import numpy

from evolve_soft_2d import utility
from evolve_soft_2d.evolve import lsystems
from evolve_soft_2d.result import analyse
from evolve_soft_2d.unit import create

################################################################################

def g_a(
    template,
    gen: int,
    prob: list,
    point: list,
    meth: str,
    ) -> None:

    #   Initialisations
    pop_all = []
    pop_best = []

    #   Check if the unit generation method is set to L-Systems
    if meth == "l":

        all_max = ls_all_max
        all_min = ls_all_min
        
    #   Check if the unit generation method is set to CPPNs
    elif meth == "c":

        all_max = cppn_all_max
        all_min = cppn_all_min

    #   Check if the unit generation method is set to random
    else:

        all_max = [n_u + 1, len(template.e_internal) + 1]
        all_min = [1, 0]

    print(all_max)
    print(all_min)

    #   Obtain the number of chromosomes
    chrom = len(all_max)

    print(chrom)
        
    #   Generate a list of unit parameters
    pop_i, param_i = create.gen_init_units(template, n_u, meth, [all_max, all_min])

    print(pop_i)
    print(param_i)

    pop_all.append(pop_i)

    #   Loop through the number of generations
    for _ in range(0, gen):

        #   Initialisations
        par_i = []

        #   Run the current population of units
        fp_lu, empty_id, full_id = create.run_units(template, pop_i, meth)
        analyse.rank_u(template, meth, fp_lu)
        
        #   Rank the current population of units
        param_rank_i = analyse.rank_pop(fp_lu, empty_id, full_id, param_i)

        print(param_rank_i)

        #   Loop through the population size in increments of 2 
        for _ in range(0, n_u, 2):

            #   Select two parents
            par_1 = sel_par(n_u, param_rank_i)
            par_2 = sel_par(n_u, param_rank_i)

            print(par_1)
            print(par_2)

            #   Potentially apply crossover between the parents to produce two children
            chi_1, chi_2 = crossover(chrom, prob[0], point[0], par_1, par_2)

            #   Potentially randomly mutate the children
            chi_1 = mut(chrom, all_max, all_min, prob[1], point[1], chi_1)
            chi_2 = mut(chrom, all_max, all_min, prob[1], point[1], chi_2)

            #   Potentially apply biased mutation to the children
            chi_1 = mut_bias(chrom, all_max, all_min, prob[2], point[2], chi_1)
            chi_2 = mut_bias(chrom, all_max, all_min, prob[2], point[2], chi_2)

            print(chi_1)
            print(chi_2)

            #   Add the children to the list of parameters
            par_i.append(chi_1)
            par_i.append(chi_2)

        #   Store the best member of the current generation
        pop_best.append(param_rank_i[0])

        #   Check if the unit generation method is CPPNs
        if meth == "c":

            #   Create the list of unit class objects from the list of parameters
            pop_i = create.gen_bred_units(template, meth, par_i, c_mod_max = all_max[0])

        else:

            #   Create the list of unit class objects from the list of parameters
            pop_i = create.gen_bred_units(template, meth, par_i)

        #   Store the new list of parameters
        param_i = par_i

        #   Add the current population to the list of all populations
        pop_all.append(pop_i)

    return

################################################################################

def gen_par(
    pop_n: int,
    chrom: int,
    all_max: list,
    all_min: list,
    ) -> list:
    """Generate a random initial population

    Parameters
    ----------
    pop_n : int
        The population size
    chrom : int
        The number of chromosomes a member has
    all_max : list
        The list of the maximum allele values
    all_min : list
        The list of the minimum allele values

    Returns
    -------
    list
        The population
    """    

    #   Initialisations
    pop_i = []

    #   Loop through the number of members to generate
    for i in range(0, pop_n):

        pop_i.append([])

        #   Loop through the number of chromosomes
        for j in range(0, chrom):

            #   Add a new random allele value
            pop_i[i].append(numpy.random.randint(all_min[j], all_max[j]))

    return pop_i

################################################################################

def sel_par(
    pop_n: int,
    param_i: list,
    ) -> list:
    """Select a parent to produce the new generation

    Parameters
    ----------
    pop_n : int
        The population size
    fit_i : list
        The fitness indices
    pop : list
        The population
    c : float
        The fitness function constant

    Returns
    -------
    list
        The parent
    """    

    #   Initialisations
    i = 0
    par_found = False
    par = param_i[0]

    #   Generate a random threshold between 0 and 1
    x = numpy.random.uniform(size = 1)

    #   Loop until a parent is found
    while i < pop_n and not par_found:

        #   Calculate the upper and lower bounds for comparison with the threshold
        x_bel = sum([fit_weight(pop_n, j) for j in range(1, i + 1)])
        x_abo = sum([fit_weight(pop_n, j) for j in range(1, i + 2)])

        #   Check if the threshold is within the bounds
        if x_bel < x and x <= x_abo:

            #   Assign the parent
            par = param_i[i]

            #   Set the exit flag
            par_found = True

        else:

            #   Increment the counter
            i += 1

    return par

################################################################################

def crossover(
    chrom: int,
    prob_co: float,
    point_co: int,
    par_1: list,
    par_2: list
    ) -> [list, list]:
    """Apply genetic crossover to two children from two parents

    Parameters
    ----------
    chrom : int
        The number of chromosomes a member has
    prob_co : float
        The probability of crossover occurring
    point_co : int
        The potential number of crossover points
    par_1 : list
        The first parent
    par_2 : list
        The second parent

    Returns
    -------
    [list, list]
        The two children
    """    

    #   Initialisations
    chi_1 = par_1[:]
    chi_2 = par_2[:]

    #   Loop through the number of potential crossover points
    for _ in range(0, point_co):

        #   Generate a random value between 0 and 1
        x = numpy.random.uniform(size = 1)

        #   Check if the random value is below the probability
        if x < prob_co:

            #   Select a random index for crossover to occur
            co_point = numpy.random.randint(1, chrom - 2)

            #   Create temporary copies of the children
            chi_1_c = chi_1[:]
            chi_2_c = chi_2[:]

            #   Apply the crossover
            chi_1 = chi_1_c[:co_point] + chi_2_c[co_point:]
            chi_2 = chi_2_c[:co_point] + chi_1_c[co_point:]

    return chi_1, chi_2

################################################################################

def mut(
    chrom: int,
    all_max: list,
    all_min: list,
    prob_m: float,
    point_m: int,
    chi: list,
    ) -> list:
    """Apply random mutation to a member

    Parameters
    ----------
    chrom : int
        The number of chromosomes a member has
    all_max : list
        The list of the maximum allele values
    all_min : list
        The list of the minimum allele values
    prob_m : float
        The probability of random mutation occurring
    point_m: int
        The number of potential random mutations
    chi : list
        The member

    Returns
    -------
    list
        The mutated member
    """

    #   Loop for the number of potential random mutations
    for _ in range(0, point_m):

        #   Generate a random value between 0 and 1
        x = numpy.random.uniform(size = 1)

        #   Check if the value is less than the probability for random mutations
        if x < prob_m:

            #   Select the point of random mutation
            m_point = numpy.random.randint(0, chrom - 1)

            #   Apply the mutation
            chi[m_point] = numpy.random.randint(all_min[m_point], all_max[m_point])

    return chi

################################################################################

def mut_bias(
    chrom: int,
    all_max: list,
    all_min: list,
    prob_bm: float,
    point_bm: float,
    chi: list,
    ) -> list:
    """Apply biased mutation to a member

    Parameters
    ----------
    chrom : int
        The number of chromosomes a member has
    all_max : list
        The list of the maximum allele values
    all_min : list
        The list of the minimum allele values
    prob_bm : float
        The probability of biased mutation occurring
    point_bm : float
        The number of potential biased mutations
    chi : list
        The member

    Returns
    -------
    list
        The mutated member
    """    

    #   Loop for the number of potential biased mutations
    for _ in range(0, point_bm):

        #   Generate a random value between 0 and 1
        x = numpy.random.uniform(size = 1)

        #   Check if the value is less than the probability for biased mutations
        if x < prob_bm:

            #   Select the point of biased mutation
            bm_point = numpy.random.randint(0, chrom - 1)

            #   Generate a random value between 0 and 1
            y = numpy.random.uniform(size = 1)

            #   Check if the value is less than 0.5
            if y < 0.5:

                #   Decrease the allele value by 1
                chi[bm_point] -= 1

            else:

                #   Increase the allele value by 1
                chi[bm_point] += 1

            #   Ensure that the allele value is not out of bounds
            chi[bm_point] = utility.clean_int(chi[bm_point], all_max[bm_point], lb = all_min[bm_point])

    return chi

################################################################################

def fit_weight(
    p: int,
    r: int,
    c: int = 1,
    ) -> float:
    """The fitness weighting function

    Parameters
    ----------
    p : int
        The population size
    c : int
        The fitness constant
    r : int
        The rank of the selected member

    Returns
    -------
    float
        The fitness weight
    """

    y = (2*((p + 1 - r)**c))/(p**2 + p)

    return y

################################################################################

#   The number of units per generation
n_u = 1000

#   The L-System allele ranges
ls_all_max = [n_u + 1, len(lsystems.a_all), len(lsystems.e_var) + 1, 6, 6]
ls_all_min = [1, 0, 1, 2, 1]

#   The CPPN allele ranges
cppn_all_max = [n_u + 1, 2, 2, 11, 32, 101]
cppn_all_min = [1, 1, 1, 2, 2, 0]