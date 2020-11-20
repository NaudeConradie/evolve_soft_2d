##  Functions used for obtaining and inspecting results

#   Imports
import csv
import fileinput
import linecache
import numpy
import pandas

from scipy.optimize import curve_fit
from scipy.spatial.distance import directed_hausdorff

from evolve_soft_2d import plotting, utility
from evolve_soft_2d.evolve import gen_alg, lsystems
from evolve_soft_2d.file_paths import create_fp_file
from evolve_soft_2d.result import obtain
from evolve_soft_2d.unit import create, modify

################################################################################

def monte_carlo(
    template,

    meth: str,
    ) -> None:
    """Perform a Monte Carlo analysis on a population of units

    Parameters
    ----------
    template : template
        The unit template parameters
    meth : str
        The unit generation method
        l : L-Systems
        c : CPPNs
        r : Random generation
    """

    #   Check if the unit generation method is set to L-Systems
    if meth == "l":
        
        #   Generate a list of unit parameters
        l_u = create.gen_init_units(template, gen_alg.n_u, meth, [gen_alg.ls_all_max, gen_alg.ls_all_min])[0]

    #   Check if the unit generation method is set to CPPNs
    elif meth == "c":

        #   Generate a list of unit parameters
        l_u = create.gen_init_units(template, gen_alg.n_u, meth, [gen_alg.cppn_all_max, gen_alg.cppn_all_min])[0]

    #   Check if the unit generation method is set to random
    else:

        #   Generate a list of unit parameters
        l_u = create.gen_init_units(template, gen_alg.n_u, meth, [[gen_alg.n_u + 1, len(template.e_internal) + 1], [1, 0]])[0]

    #   Run the population of units
    fp_lu, empty_id, full_id = create.run_units(template, l_u, meth)

    #   Rank the population of units
    rank_u(template, meth, empty_id, full_id, fp_lu)

    return

################################################################################

def rank_u(
    template,
    g_meth: str,
    empty_id: str,
    full_id: str,
    fp_lu: list,
    ) -> None:
    """Rank units according to their performance

    Parameters
    ----------
    template
        The unit template parameters
    fp_lu : str
        The file path of the log file of units created during the last simulation
    """

    #   Initialisations
    v = []
    data = pandas.DataFrame()

    label = []
    label.append("Constraint Energy X")
    label.append("Constraint Energy Y")
    label.append("Constraint Energy")
    label.append("Internal Energy X")
    label.append("Internal Energy Y")
    label.append("Internal Energy")

    #   Read the list of units created during the last simulation
    lu = ranking(fp_lu, empty_id, full_id)

    #   Loop through all labels
    for i in label:

        #   Append the data to the list
        v.append(obtain.read_all(template, lu, i))

    #   Add the Hausdorff distance to the list of labels
    label.append("Hausdorff Distance")

    #   Read the Hausdorff distance variables
    v.append(obtain.read_all_hd(template, lu))

    label.append("Number of Elements Removed")

    #   Create a list of the number of elements removed from every element
    v.append([utility.find_int_in_str(i) for i in lu])

    #   Check if the template case is 1
    if template.case == 1:

        #   Add the necessary labels for comparison
        label.append("Height to Width Ratio")
        label.append("Absolute Change In Width")

        #   Read the fitness data
        v_fit = obtain.read_all_fit(template, lu)

        #   Append the fitness data to the list
        for i in v_fit:
            v.append(i)

    if g_meth == "l":

        label.append("Axiom ID")
        label.append("Number of Rules")
        label.append("Rule Length")
        label.append("Number of Iterations")

        l_axiom = []
        l_rule_n = []
        l_rule_l = []
        l_i = []

        for i in lu:

            curr_mod = utility.open_v(template, i)

            l_axiom.append(lsystems.a_all.index(curr_mod.ls.axiom))
            l_rule_n.append(len(curr_mod.ls.gramm))
            l_rule_l.append(ls_avg_rule_l(curr_mod.ls))
            l_i.append(curr_mod.ls.n)

        v.append(l_axiom)
        v.append(l_rule_n)
        v.append(l_rule_l)
        v.append(l_i)

    elif g_meth == "c":

        label.append("Model ID")
        label.append("Scale")
        label.append("Number of Hidden Layers")
        label.append("Size of the Initial Hidden Layer")
        label.append("Element Removal Threshold")

        c_m_id = []
        c_scale = []
        c_hl_n = []
        c_hl_s = []
        c_thresh = []

        for i in lu:

            curr_mod = utility.open_v(template, i)

            c_m_id.append(curr_mod.cp.mod_id)
            c_scale.append(curr_mod.cp.cppn.scale)
            c_hl_n.append(curr_mod.cp.cppn.hl_n)
            c_hl_s.append(curr_mod.cp.cppn.hl_s)
            c_thresh.append(curr_mod.cp.cppn.thresh)

        v.append(c_m_id)
        v.append(c_scale)
        v.append(c_hl_n)
        v.append(c_hl_s)
        v.append(c_thresh)

    #   Store the list of units
    data["Unit ID"] = lu

    #   Loop through the values
    for i in range(0, len(v)):

        #   Add the values to the dataframe
        data[label[i]] = v[i]

    #   Check if the template case is 1
    if template.case == 1:

        #   Studentize the fitness values
        data["Height to Width Ratio"] = (data["Height to Width Ratio"] - data["Height to Width Ratio"].mean())/data["Height to Width Ratio"].std()
        data["Absolute Change In Width"] = (data["Absolute Change In Width"] - data["Absolute Change In Width"].mean())/data["Absolute Change In Width"].std()

        #   Calculate a single fitness value
        data["Fitness"] = data["Height to Width Ratio"] + data["Absolute Change In Width"]

    #   Read the timestamp of the simulation
    tm = utility.read_str(fp_lu[0], -25, -4)

    # plot_data = 

    # #   Plot the desired graphs from the results
    # plotting.plot_all(template, v, n_e, label, tm)

    data = data.replace(0, numpy.nan)

    # plotting.hist_all(template, tm, data)
    # plotting.scat_all(template, tm, data)

    #   Sort the dataframe according to the fitness values
    data.sort_values(by = ["Fitness"], ascending = False, inplace = True, ignore_index = True)
    
    #   Save the list of best performing units
    # data["Unit ID"].to_csv(fp_lu[5], header = False, index = False)

    return data

################################################################################

def ranking(
    fp_lu: list,
    empty_id: str,
    full_id: str,
    ) -> list:

    # #   Read the list of all units created
    # lu = obtain.read_lu(fp_lu[0])

    #   Read the ranked list of all units created
    with open(fp_lu[5], 'r') as f:
        lu_rank = f.read()

    try:
        
        #   Read the list of empty units created
        with open(fp_lu[3], 'r') as f:
            lu_empty = f.read()

        #   Replace the placeholder empty unit ID in the ranked list with all generated empty units
        lu_rank = lu_rank.replace(empty_id, lu_empty)

    except:
        
        #   Replace the placeholder empty unit ID in the ranked list with a blank space
        lu_rank = lu_rank.replace(empty_id, "")

    try:

        #   Read the list of full units created
        with open(fp_lu[4], 'r') as f:
            lu_full = f.read()

        #   Replace the placeholder full unit ID in the ranked list with all generated empty units
        lu_rank = lu_rank.replace(full_id, lu_full)

    except:

        #   Replace the placeholder full unit ID in the ranked list with a blank space
        lu_rank = lu_rank.replace(full_id, "")

    # try:

    #     #   Read the list of failed units created
    #     with open(fp_lu[2], 'r') as f:
    #         lu_fail = f.read()

    #     #   Append the list of failed units to the ranked list
    #     lu_rank += lu_fail

    # except:
    #     pass

    #   Format the list of ranked units
    lu_rank = list(lu_rank.split("\n"))
    while "" in lu_rank:
        lu_rank.remove("")

    return lu_rank

################################################################################

def rank_pop(
    fp_lu: list,
    empty_id: str,
    full_id: str,
    par: list,
    ) -> list:

    #   Initialisations
    data = pandas.DataFrame()

    #   Read the list of all units created
    lu = obtain.read_lu(fp_lu[0])

    #   Read the ranked list of all units created
    with open(fp_lu[5], 'r') as f:
        lu_rank = f.read()

    try:
        
        #   Read the list of empty units created
        with open(fp_lu[3], 'r') as f:
            lu_empty = f.read()

        #   Replace the placeholder empty unit ID in the ranked list with all generated empty units
        lu_rank = lu_rank.replace(empty_id, lu_empty)

    except:
        
        #   Replace the placeholder empty unit ID in the ranked list with a blank space
        lu_rank = lu_rank.replace(empty_id, "")

    try:

        #   Read the list of full units created
        with open(fp_lu[4], 'r') as f:
            lu_full = f.read()

        #   Replace the placeholder full unit ID in the ranked list with all generated empty units
        lu_rank = lu_rank.replace(full_id, lu_full)

    except:

        #   Replace the placeholder full unit ID in the ranked list with a blank space
        lu_rank = lu_rank.replace(full_id, "")

    try:

        #   Read the list of failed units created
        with open(fp_lu[2], 'r') as f:
            lu_fail = f.read()

        #   Append the list of failed units to the ranked list
        lu_rank += lu_fail

    except:
        pass

    #   Format the list of ranked units
    lu_rank = list(lu_rank.split("\n"))
    while "" in lu_rank:
        lu_rank.remove("")

    #   Add the relevant data to the dataframe
    data["Unit ID"] = lu
    data["Parameters"] = par

    #   Create the sorting index according to the ranked unit list
    lu_rank_index = dict(zip(lu_rank, range(0, len(lu_rank))))

    #   Add the sorting index to the dataframe
    data["Rank"] = data["Unit ID"].map(lu_rank_index)

    #   Sort the dataframe according to the sorting index
    data.sort_values(["Rank"], ascending = [True], inplace = True)

    #   Remove the sorting index from the dataframe
    data.drop("Rank", 1, inplace = True)

    par_sort = data["Parameters"].tolist()

    return par_sort

################################################################################

def ls_avg_rule_l(ls) -> float:

    ls_rules = [i[1] for i in ls.gramm]

    avg_rule_l = utility.avg_str_l(ls_rules)

    return avg_rule_l

################################################################################

def constraint_energy(
    template,
    l: str,
    ) -> None:
    """Calculate the constraint energy for a unit

    Parameters
    ----------
    template : template
        The unit template parameters
    l : str
        The label for the results file
        Either a template or unit identifier
    """    
    
    #   Initialisations
    label = []
    label.append("Displacement X")
    label.append("Displacement Y")
    label.append("Displacement")
    label.append("Reaction Force X")
    label.append("Reaction Force Y")
    label.append("Reaction Force")

    label_c_e = []
    label_c_e.append("Constraint Energy X")
    label_c_e.append("Constraint Energy Y")
    label_c_e.append("Constraint Energy")

    v = numpy.zeros((len(label), template.n_steps + 1, template.n_n))

    #   Loop through all the variable labels
    for i in range(0, len(label)):

        #   Create the file path of the results file
        fp_r_f = create_fp_file(template, label[i] + "_" + l, "r")

        #   Store the results from the file
        v[i] = numpy.genfromtxt(fp_r_f, delimiter = ",")

    #   Decrement the node IDs of the external nodes by 1 to be used as array indices
    n_external_i = [i - 1 for i in template.n_external]

    #   Store only the external node values
    v_ex = v[:, :, n_external_i]

    #   Loop through every constraint energy label
    for i in range(0, len(label_c_e)):

        #   Initialise the constraint energy array
        c_e = numpy.zeros(len(v_ex[i]))

        #   Loop through every step in the unit
        for j in range(0, len(v_ex[i])):

            #   Loop through every external node
            for k in range(0, len(v_ex[i, i])):

                #   Calculate the constraint energy for the current step
                c_e[j] += v_ex[i, j, k]*v_ex[i + len(label_c_e), j, k]/1000

        #   Save the constraint energy to a .csv file
        save_numpy_array_to_csv(template, label_c_e[i] + "_" + l, c_e)

    return

################################################################################

def internal_energy(
    template,
    l: str,
    ) -> None:
    """Calculate the internal energy for a unit

    Parameters
    ----------
    template 
        The unit template parameters
    l : str
        The label for the results file
        Either a template or unit identifier
    """    
    
    #   Initialisations
    label = []
    label.append("Comp 11 of Total Strain")
    label.append("Comp 22 of Total Strain")
    label.append("Total Strain Energy Density")

    label_i_e = []
    label_i_e.append("Internal Energy X")
    label_i_e.append("Internal Energy Y")
    label_i_e.append("Internal Energy")

    s = numpy.zeros((len(label), template.n_steps + 1, template.n_n))

    #   Loop through all the variable labels
    for i in range(0, len(label)):

        #   Create the file path of the results file
        fp_r_f = create_fp_file(template, label[i] + "_" + l, "r")

        #   Store the results from the file
        s[i] = numpy.genfromtxt(fp_r_f, delimiter = ",")

    #   Loop through all the internal energy labels
    for i in range(0, len(label_i_e)):

        #   Initialise the internal energy array
        i_e = numpy.zeros(len(s[i]))

        #   Loop through every step in the unit
        for j in range(0, len(s[i])):

            #   Loop through every external node
            for k in range(0, len(s[i, i])):

                #   Calculate the internal energy for the current step
                i_e[j] += s[i, j, k]

        #   Save the internal energy to a .csv file
        save_numpy_array_to_csv(template, label_i_e[i] + "_" + l, i_e)

    return

################################################################################

def disp(
    template,
    l: str,
    ) -> list:
    """Read the displacement values for a unit

    Parameters
    ----------
    template : template
        The unit template parameters
    l : str
        The label for the results file
        Either a template or unit identifier

    Returns
    -------
    list
        The list of displacement values
    """    

    #   Initialisations
    label = []
    label.append("Displacement X")
    label.append("Displacement Y")

    d = numpy.zeros((len(label), template.n_steps + 1, template.n_n))

    #   Loop through the list of labels
    for i in range(0, len(label)):

        #   Create the file path of the results file
        fp_r_f = create_fp_file(template, label[i] + "_" + l, "r")

        #   Store the results from the file
        d[i] = numpy.genfromtxt(fp_r_f, delimiter = ",")

    #   Decrement the node IDs of the external nodes by 1 to be used as array indices
    n_external_i = [i - 1 for i in template.n_external]

    #   Store only the external node values
    d_ex = d[:, :, n_external_i]

    #   Store only the final values in the simulation
    d_ex = d_ex[:, template.n_steps]

    return d_ex

################################################################################

def disp_fit(
    template,
    l: str,
    ) -> None:
    """Calculate and save the displacement fitness values

    Parameters
    ----------
    template : template
        The unit template parameters
    l : str
        The label for the results file
        Either a template or unit identifier
    """

    #   Obtain the deformed boundary displacement values
    d_def = disp(template, l)

    #   Check if the template case is 1
    if template.case == 1:

        #   Obtain the fitness measures for case 1
        fit_m = disp_fit_1(template, d_def)

        #   Define the label
        label = "Case 1 Fitness Measures_"

    #   Save the internal energy to a .csv file
    save_numpy_array_to_csv(template, label + l, fit_m)

    return

################################################################################

def disp_fit_1(
    template,
    d: list,
    ) -> list:
    """Obtain the displacement fitness values for case 1

    Parameters
    ----------
    template : template
        The unit template parameters
    d : list
        The displacement values

    Returns
    -------
    list
        The fitness values
    """

    #   Initialisations
    d_a = []

    #   Split the displacement values according to the four sides of the unit
    d_b, d_t, d_l, d_r = split_d(template, d)

    #   Swap the x- and y-coordinates of the left and right sides
    d_l = utility.list_swap_i(d_l, 0, 1)
    d_r = utility.list_swap_i(d_r, 0, 1)

    #   Add all side values to one list
    d_split = [d_b, d_t, d_l, d_r]

    #   Loop through the side values
    for i in d_split:

        #   Fit a horizontal curve through the side coordinates
        a = curve_fit(utility.f_const, i[0], i[1])[0]

        #   Add the curve height to the list of side displacements
        d_a.append(a[0])

    #   Calculate the height to width ratio
    h_w = (d_a[1] - d_a[0])/(d_a[3] - d_a[2])

    #   Check if the deformed width is larger than the original width
    if abs(d_a[3] - d_a[2]) >= template.x_s:

        #   Calculate the deformed width to original width ratio
        w = abs((d_a[3] - d_a[2])/template.x_s)

    else:

        #   Calculate the deformed width to original width ratio
        w = abs(template.x_s/(d_a[3] - d_a[2]))

    fit_m = [h_w, w]
    fit_m = numpy.array(fit_m)

    return fit_m

################################################################################

def hausdorff_d(
    template,
    l: str,
    ) -> None:
    """Calculate the Hausdorff distance values for a unit

    Parameters
    ----------
    template : template
        The unit template parameters
    l : str
        The label for the results file
        Either a template or unit identifier
    """    

    #   Read the displacement values of the unit
    d_def = disp(template, l)
    d_def = numpy.transpose(d_def)

    d_des = numpy.transpose(template.d)

    #   Calculate the Hausdorff distance between the unit and the template displacements
    h_d = directed_hausdorff(d_def, d_des)

    #   Save the internal energy to a .csv file
    save_numpy_array_to_csv(template, "Hausdorff Distance_" + l, h_d)

    return

################################################################################

def split_d(
    template,
    d: numpy.array,
    ) -> [numpy.array, numpy.array, numpy.array, numpy.array]:
    """Split the displacement values according to the four sides of the unit, including corner nodes

    Parameters
    ----------
    template : template
        The unit template parameters
    d : numpy.array
        The displacement values

    Returns
    -------
    [numpy.array, numpy.array, numpy.array, numpy.array]
        The split displacement values
    """    

    d_b = d[:, 0:template.x_n]
    d_t = d[:, -(template.x_n):len(template.n_external)]
    d_l = [[d[0, 0]] + [d[0, i] for i in range(template.x_n, len(template.n_external) - (template.x_n - 1), 2)], [d[1, 0]] + [d[1, i] for i in range(template.x_n, len(template.n_external) - (template.x_n - 1), 2)]]
    d_r = [[d[0, i] for i in range(template.x_n - 1, len(template.n_external) - (template.x_n - 1), 2)] + [d[0, len(template.n_external) - 1]], [d[1, i] for i in range(template.x_n - 1, len(template.n_external) - (template.x_n - 1), 2)] + [d[1, len(template.n_external) - 1]]]

    return d_b, d_t, d_l, d_r

################################################################################

def save_numpy_array_to_csv(
    template,
    l: str,
    data: numpy.array,
    ) -> None:
    """Write the results to .csv files

    Parameters
    ----------
    template 
        The unit template parameters
    l : str
        The label of the data
    data : numpy.array
        The results to be stored
    """    
    
    #   Create the file path of the results file
    fp_r_f = create_fp_file(template, l, "r")

    #   Write the data to the results file
    numpy.savetxt(fp_r_f, data, delimiter = ",")

    print("{}.csv saved".format(l))

    return