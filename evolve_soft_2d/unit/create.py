##  Functions used with the Marc Mentat units

#   Imports
import numpy
import time

from py_mentat import py_send

from evolve_soft_2d import classes, utility
from evolve_soft_2d.evolve import cppns, gen_alg, lsystems
from evolve_soft_2d.unit import inspect, modify, rep_grid
from evolve_soft_2d.result import analyse, obtain
from evolve_soft_2d.file_paths import create_fp_file

################################################################################

def gen_init_units(
    template,
    n: int,
    meth: str,
    a: list,
    ) -> [list, list]:
    """Generate an initial population of units

    Parameters
    ----------
    template
        The unit template parameters
    n : int
        The number of units to generate
    meth : str
        The unit generation method
        l:  L-Systems
        c:  CPPNs
        r:  Random
    a : list
        The maximum and minimum parameters

    Returns
    -------
    list
        The population
    """    

    #   Initialisations
    rem = []

    #   Check if the generation method is specified as L-Systems
    if meth == "l":

        #   Initialisations
        ls = []

        #   Generate the random parameters
        par = gen_alg.gen_par(n, 5, a[0], a[1])

        #   Loop through the list of parameters
        for i in par:

            #   Create the current L-System instance
            ls_i = lsystems.gen_lsystem(i[0], lsystems.e_vocabulary, i[1], i[2], i[3], i[4])

            #   Add the L-System instance to the list
            ls.append(ls_i)

            #   Obtain the list of elements to be removed
            rem.append(lsystems.interpret_word(template, ls_i.word))

        #   Save the list of units
        lu = numpy.array([rem, ls]).T.tolist()

    #   Check if the generation method is specified as CPPNs
    elif meth == "c":

        #   Initialisations
        cp = []

        #   Generate the random parameters
        par = gen_alg.gen_par(n, 6, a[0], a[1])

        #   Loop through the list of parameters
        for i in par:

            #   Generate the CPPN
            cp_c = cppns.cppn(i[0], a[0][1], i[2], i[3], i[4], i[5], template.x_e - 2*template.b, template.y_e - 2*template.b)

            #   Save the specific CPPN model
            cp_i = cppns.cppn_i(cp_c, i[1])

            #   Add the CPPN model to the list
            cp.append(cp_i)

            #   Obtain the list of elements to be removed
            rem.append(cppns.cppn_rem(template, cp_i.grid))

        #   Save the list of units
        lu = numpy.array([rem, cp]).T.tolist()

    else:

        #   Generate the random parameters
        par = gen_alg.gen_par(n, 2, a[0], a[1])

        #   Loop through the list of parameters
        for i in par:

            #   Add the list of elements to be removed
            rem.append(utility.sel_random(template.e_internal, f = i[1], seed = i[0]))

        #   Save the list of units
        lu = numpy.array([rem, par]).T.tolist()

    return lu, par

################################################################################

def gen_bred_units(
    template,
    meth: str,
    par: list,
    c_mod_max: int = None,
    ) -> [list, list]:
    """Generate a bred population of units

    Parameters
    ----------
    template
        The unit template parameters
    meth : str
        The unit generation method
        l:  L-Systems
        c:  CPPNs
        r:  Random
    par : list
        The list of unit parameters
    c_mod_max : int, optional
        The maximum number of CPPN models per network, by default None

    Returns
    -------
    list
        The population
    """    

    #   Initialisations
    rem = []

    #   Check if the generation method is specified as L-Systems
    if meth == "l":

        #   Initialisations
        ls = []

        #   Loop through the list of parameters
        for i in par:

            #   Create the current L-System instance
            ls_i = lsystems.gen_lsystem(i[0], lsystems.e_vocabulary, i[1], i[2], i[3], i[4])

            #   Add the L-System instance to the list
            ls.append(ls_i)

            #   Obtain the list of elements to be removed
            rem.append(lsystems.interpret_word(template, ls_i.word))

        #   Save the list of units
        lu = numpy.array([rem, ls]).T.tolist()

    #   Check if the generation method is specified as CPPNs
    elif meth == "c":

        #   Initialisations
        cp = []

        #   Loop through the list of parameters
        for i in par:

            #   Generate the CPPN
            cp_c = cppns.cppn(i[0], c_mod_max, i[2], i[3], i[4], i[5], template.x_e - 2*template.b, template.y_e - 2*template.b)

            #   Save the specific CPPN model
            cp_i = cppns.cppn_i(cp_c, i[0])

            #   Add the CPPN model to the list
            cp.append(cp_i)

            #   Obtain the list of elements to be removed
            rem.append(cppns.cppn_rem(template, cp_i.grid))

        #   Save the list of units
        lu = numpy.array([rem, cp]).T.tolist()

    else:

        #   Loop through the list of parameters
        for i in par:

            #   Add the list of elements to be removed
            rem.append(utility.sel_random(template.e_internal, f = i[1], seed = i[0]))

        #   Save the list of units
        lu = numpy.array([rem, par]).T.tolist()

    return lu

################################################################################

def run_units(
    template,
    l_u: list,
    meth: str,
    ) -> [list, str, str]:
    """Run a list of units

    Parameters
    ----------
    template : template
        The unit template parameters
    l_u : list
        The list of unit parameters
    meth : str
        The unit generation method
        l:  L-Systems
        c:  CPPNs
        r:  Random

    Returns
    -------
    [str, str]
        The file paths of the list of units generated and ranked
    """    

    #   The time the unit generation starts as a string label
    t = time.strftime("_%Y-%m-%d--%H-%M-%S", time.gmtime())

    #   Create the file paths of the log files of units created during the simulation
    #   0:  all
    #   1:  successful
    #   2:  failed
    #   3:  empty
    #   4:  full
    #   5:  ranked
    fp_lu = [create_fp_file(template, t, "l"), create_fp_file(template, t + "_success", "l"), create_fp_file(template, t + "_failed", "l"), create_fp_file(template, t + "_empty", "l"), create_fp_file(template, t + "_full", "l"), create_fp_file(template, t + "_ranked", "l")]

    #   Run the empty and full units
    empty_id, full_id = empty_full(template, fp_lu)    

    #   Loop through the list of units generated
    for i in l_u:

        #   Open the template file
        modify.open_model(template.fp_t_mud)

        #   Generate the grid with elements removed, search for any free element clusters, and update the grid and list
        grid_rem, rem = gen_grid_rem_free(template, i[0])

        #   Check if all internal elements were removed from the current unit
        if len(rem) == len(template.e_internal):

            #   Add an exception for an empty unit
            add_exc(template, fp_lu, rem, grid_rem, empty_id, meth, i[1])

        elif len(rem) == 0:

            #   Add an exception for an empty unit
            add_exc(template, fp_lu, rem, grid_rem, full_id, meth, i[1])

        else:

            #   Check if the unit generation method is set to L-Systems
            if meth == "l":

                #   Remove the elements, save and run the model and obtain the desired results
                rem_el_run_results(template, rem, grid_rem, fp_lu, ls = i[1])

            #   Check if the unit generation method is set to CPPNs
            elif meth == "c":

                #   Remove the elements, save and run the model and obtain the desired results
                rem_el_run_results(template, rem, grid_rem, fp_lu, cp = i[1])

            #   Check if the unit generation method is set to random
            else:

                #   Remove the elements, save and run the model and obtain the desired results
                rem_el_run_results(template, rem, grid_rem, fp_lu)

    return fp_lu, empty_id, full_id

################################################################################

def gen_grid_rem_free(
    template,
    rem: list,
    ) -> [list, list]:
    """Generate the representative grid with the desired elements removed, search for any free-floating element clusters, add them to the list of elements to be removed and update the representative grid

    Parameters
    ----------
    template : template
        The unit template parameters
    rem : list
        The list of elements to be removed

    Returns
    -------
    [list, list]
        The representative grid and the updated list of elements to be removed
    """    

    #   Remove the elements from the representative grid
    grid_rem = rep_grid.rem_el_grid(template, rem)

    #   Search for clusters
    (found_free, grid_label) = rep_grid.find_cluster(grid_rem)

    #   Check if free clusters were found
    if found_free:

        #   Remove the free clusters
        (grid_rem, rem_free) = rep_grid.rem_el_free_grid(template, grid_label)

        #   Update the list of removed elements
        rem = utility.add_sort_list(rem, rem_free)

    return grid_rem, rem

################################################################################

def rem_el_run_results(
    template,
    rem: list,
    grid_rem: list,
    fp_lu: list,
    log_flag: bool = True,
    ls = None,
    cp = None,
    ) -> None:
    """Remove the elements from the unit model, save and run the model, obtain the desired results from the model and reopen the template model

    Parameters
    ----------
    template : template
        The unit template parameters
    rem : list
        The list of elements to be removed
    grid_rem : list
        The representative grid with the elements removed
    fp_lu : str
        The file path of the log file of the list of all units generated
    ls : int, optional
        The L-system rule length, by default None
    cp : list, optional
        The CPPN parameters, by default None
    """    

    #   Save the current unit parameters
    curr_mod = classes.unit_p(template, rem, grid_rem, ls = ls, cp = cp)

    #   Remove the elements from the unit
    modify.rem_el(rem)

    #   Remove connections between diagonally connected elements
    modify.rem_connect(template, grid_rem)

    #   Add the internal pressure boundary conditions
    modify.add_bc_p_internal(curr_mod)

    #   Copy the unit to create neighbours
    modify.copy_neighbours(template)

    #   Add the boundary conditions for the neighbouring units
    corner_bc(curr_mod)

    #   Add the loadcase and job
    modify.add_lcase(template, 2, bc_ip + bc_fd)
    modify.add_job(2, bc_ip + bc_fd)

    #   Save the altered unit
    modify.save_model(curr_mod.fp_u_mud)

    #   Run the model
    curr_mod.run_success = modify.run_model(template, 2, curr_mod.u_id, curr_mod.fp_u_mud, curr_mod.fp_u_log[1], curr_mod.fp_u_t16[1], template.d)

    #   Check if the model run was successful
    if curr_mod.run_success:

        #   Obtain the constraint and internal energies of the current model
        curr_mod.c_e = obtain.read_xym(template, curr_mod.u_id, "Constraint Energy")
        curr_mod.i_e = obtain.read_xym(template, curr_mod.u_id, "Internal Energy")
    
        #   Log the current unit ID
        with open(fp_lu[1], "a") as f:
            print(curr_mod.u_id, file = f)

    else:

        #   Log the current unit ID
        with open(fp_lu[2], "a") as f:
            print(curr_mod.u_id, file = f)

    if log_flag:
        with open(fp_lu[0], "a") as f:
            print(curr_mod.u_id, file = f)

    #   Log the current unit parameters
    with open(curr_mod.fp_u_l, "w") as f:
        print(curr_mod, file = f)

    #   Save the current unit class object
    utility.save_v(template, curr_mod, curr_mod.u_id)

    #   Reopen the template
    modify.open_model(template.fp_t_mud)

    return

################################################################################

def empty_full(
    template,
    fp_lu: list,
    ) -> [str, str]:
    """Run the empty and full unit cases

    Parameters
    ----------
    template : template
        The unit template parameters
    fp_lu : str
        The file path of the log file of the list of all units generated

    Returns
    -------
    [str, str]
        The unit IDs of the empty and full units
    """    

    #   Generate the grid with all elements removed
    grid_rem_e, rem_e = gen_grid_rem_free(template, template.e_internal)

    #   Remove the elements, save and run the model and obtain the desired results
    rem_el_run_results(template, rem_e, grid_rem_e, fp_lu, log_flag = False)

    #   Generate the unit ID of the empty unit
    empty_id = str(len(rem_e)) + "_" + utility.gen_hash(utility.list_to_str(rem_e, "_"))

    #   Generate the grid with all elements removed
    grid_rem_f, rem_f = gen_grid_rem_free(template, [])

    #   Remove the elements, save and run the model and obtain the desired results
    rem_el_run_results(template, rem_f, grid_rem_f, fp_lu, log_flag = False)

    #   Generate the unit ID of the full unit
    full_id = str(len(rem_f)) + "_" + utility.gen_hash(utility.list_to_str(rem_f, "_"))

    return empty_id, full_id

################################################################################

def add_exc(
    template,
    fp_lu: str,
    rem: list,
    grid_rem: list,
    u_id: str,
    meth: str,
    meth_v: list,
    ) -> None:
    """Add an exception case to the list of units

    Parameters
    ----------
    template : template
        The unit template parameters
    fp_lu : str
        The file path of the list of exception units
    rem : list
        The list of elements to be removed
    grid_rem : list
        The representative grid with elements removed
    u_id : str
        The relevant unit ID
    meth : str
        The unit generation method
    meth_v : list
        The unit generation method class object
    """    

    #   Check if the unit generation method is set to L-Systems
    if meth == "l":

        #   Save the unit parameters
        curr_mod = classes.unit_p(template, rem, grid_rem, ls = meth_v)

    #   Check if the unit generation method is set to CPPNs
    elif meth == "c":

        #   Save the unit parameters
        curr_mod = classes.unit_p(template, rem, grid_rem, cp = meth_v)

    #   Check if the unit generation method is set to random
    else:

        #   Save the unit parameters
        curr_mod = classes.unit_p(template, rem, grid_rem)

    #   Read the constraint and internal energy values from the prerun case
    curr_mod.c_e = obtain.read_xym(template, u_id, "Constraint Energy")
    curr_mod.i_e = obtain.read_xym(template, u_id, "Internal Energy")

    #   Log the unit parameters
    with open(curr_mod.fp_u_l, "w") as f:
        print(curr_mod, file = f)

    #   Save the unit class object
    utility.save_v(template, curr_mod, curr_mod.u_id)

    if len(rem) != 0:

        #   Log the current unit ID
        with open(fp_lu[3], "a") as f:
            print(curr_mod.u_id, file = f)

    else:

        #   Log the current unit ID
        with open(fp_lu[4], "a") as f:
            print(curr_mod.u_id, file = f)

    with open(fp_lu[0], "a") as f:
        print(curr_mod.u_id, file = f)

    return

################################################################################

def temp_create(template) -> list:
    """Create a template from which elements may be removed

    Case:   1
    Unidirectional extension
    Case:   2
    Bidirectional extension
    Case:   2
    Pure shear strain
    Case:   3
    Elongation of one side

    Parameters
    ----------
    template : template
        The unit template parameters
    """

    #   Clear the workspace
    py_send("*new_model yes")

    #   Construct the grid
    modify.create_nodes(template)
    modify.create_elements(template)

    #   Add the application graph
    modify.add_ramp(template)

    #   Determine which case is selected and apply the appropriate boundary conditions
    if template.case == 1:
        template_1_bc(template)
    elif template.case == 2:
        template_2_bc(template)
    elif template.case == 3:
        template_3_bc(template)
    else:
        template_4_bc(template)

    #   Add template properties
    modify.add_geom_prop()
    modify.add_mat_ogden(template.ogd_mat)
    modify.add_contact_body()

    #   Add the loadcase
    modify.add_lcase(template, 1, template_bc_fd[template.case - 1])

    #   Add the job
    modify.add_job(1, template_bc_fd[template.case - 1])
    modify.save_model(template.fp_t_mud)

    #   Run the model
    template.run_success = modify.run_model(template, 1, template.t_id, template.fp_t_mud, template.fp_t_log, template.fp_t_t16)

    #   Check if the model run was successful
    if template.run_success:

        template.d = analyse.disp(template, template.t_id + "_1")

        #   Obtain the constraint energy of the current model
        template.c_e = obtain.read_xym(template, template.t_id, "Constraint Energy")
        template.i_e = obtain.read_xym(template, template.t_id, "Internal Energy")

    #   Save the template
    modify.save_model(template.fp_t_mud)

    #   Log the template parameters
    with open(template.fp_t_l, "w") as f:
        print(template, file = f)
        
    utility.save_v(template, template, template.t_id)

    return

################################################################################

def template_1_bc(template) -> None:
    """Apply the displacement boundary conditions for case 1

    Parameters
    ----------
    template : template
        The unit template parameters
    """    

    #   Add the boundary conditions
    modify.add_bc_fd_edge("yy1", "y", "y", 0, 0)
    modify.add_bc_fd_edge("yy2", "y", "y", template.y_s, template.d_mag, tab_nam = template.tab_nam)
    modify.add_bc_fd_edge("xx1", "x", "x", 0, 0)
    modify.add_bc_fd_edge("xx2", "x", "x", template.x_s, 0)

    return


################################################################################

def template_2_bc(template) -> None:
    """Apply the displacement boundary conditions for case 2

    Parameters
    ----------
    template : template
        The unit template parameters
    """

    #   Add the boundary conditions
    modify.add_bc_fd_edge("yy1", "y", "y", 0, 0)
    modify.add_bc_fd_edge("yy2", "y", "y", template.y_s, template.d_mag, tab_nam = template.tab_nam)
    modify.add_bc_fd_edge("xx1", "x", "x", 0, 0)
    modify.add_bc_fd_edge("xx2", "x", "x", template.x_s, template.d_mag, tab_nam = template.tab_nam)

    return

################################################################################

def template_3_bc(template) -> None:
    """Apply the displacement boundary conditions for case 3

    Parameters
    ----------
    template : template
        The unit template parameters
    """

    #   Add the boundary conditions
    modify.add_bc_fd_edge("xy1", "x", "y", 0, 0)
    modify.add_bc_fd_edge("xy2", "x", "y", template.y_s, template.d_mag, tab_nam = template.tab_nam)
    modify.add_bc_fd_edge("yx1", "y", "x", 0, 0)
    modify.add_bc_fd_edge("yx2", "y", "x", template.x_s, 0)

    return

################################################################################

def template_4_bc(template) -> None:
    """Apply the displacement boundary conditions for case 3

    Parameters
    ----------
    template : template
        The unit template parameters
    """

    #   Add the boundary conditions
    modify.add_bc_fd_edge("yy1", "y", "y", 0,            0)
    modify.add_bc_fd_edge("yy2", "y", "y", template.y_s, 0)
    
    modify.add_bc_fd_node("xf1", "x", 1,            0)
    modify.add_bc_fd_node("xf2", "x", template.x_n, 0)
    modify.add_bc_fd_node("xn",  "x", template.n_n - template.x_n + 1, -template.d_mag, tab_nam = template.tab_nam)
    modify.add_bc_fd_node("xp",  "x", template.n_n, template.d_mag, tab_nam = template.tab_nam)

    return

################################################################################

def corner_bc(unit_p) -> None:
    """Apply the fixed displacement boundary conditions to the corner nodes

    Parameters
    ----------
    unit_p : unit_p
        The unit parameters
    """

    #   Find the corner node IDs
    n_corn = inspect.find_n_corn(unit_p)

    #   Apply the fixed displacement boundary conditions
    modify.add_bc_fds_nodes("xy_corn", n_corn, 0)

    return

################################################################################

#   Lists of template boundary condition labels

template_bc_fd = [
    ["bc_fd_yy1", "bc_fd_yy2", "bc_fd_xx1", "bc_fd_xx2"],
    ["bc_fd_yy1", "bc_fd_yy2", "bc_fd_xx1", "bc_fd_xx2"],
    ["bc_fd_xy1", "bc_fd_xy2", "bc_fd_yx1", "bc_fd_yx2"],
    ["bc_fd_yy1", "bc_fd_yy2", "bc_fd_xf1", "bc_fd_xf2", "bc_fd_xn", "bc_fd_xp"],
    ]

bc_ip = ["bc_load_yp", "bc_load_xn", "bc_load_yn", "bc_load_xp"]

bc_fd = ["bc_fd_xy_corn"]