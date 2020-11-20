##  Functions used for obtaining and inspecting results

#   Imports
import csv
import itertools
import linecache
import numpy
import os.path
import re
import time

from evolve_soft_2d import plotting, utility
from evolve_soft_2d.file_paths import create_fp_file
from evolve_soft_2d.log import m_log, en_log
from evolve_soft_2d.unit import modify

from py_mentat import py_send, py_get_float

################################################################################

def check_out(
    fp_mud: str,
    fp_log: str,
    fp_t16: str,
    j_id: int
    ) -> bool:
    """Check if the updated output files exist

    Parameters
    ----------
    fp_mud : str
        The file path of the model file
    fp_log : str
        The file path of the model log file
    fp_t16 : str
        The file path of the model t16 file
    j_id : int
        The job ID

    Returns
    -------
    bool
        True if the updated output files exist, false if otherwise
    """

    #   Initialisations
    t0 = time.time()
    success = False
    decided = False

    #   Time to wait for files
    t = 1

    #   Time allowed for files to run
    t_max = 300
    t_out = 0

    #   Text to look for when searching the log files
    exit_number_str = re.compile("exit number", re.IGNORECASE)
    access_viol_str = re.compile("access violation", re.IGNORECASE)

    #   Obtain the timestamp of the last time the unit file was modified
    t_mud = os.path.getmtime(fp_mud)
    
    #   Wait until the log file exists and has been updated
    utility.wait_file_exist(fp_log, "log", t)
    utility.wait_file_update(fp_log, t_mud, "log", t)

    #   Loop until an exit number is detected
    while 1:

        t_out += 1

        #   Search the log file for an exit number
        found_e_n, e_n = utility.search_text_file(fp_log, exit_number_str)
        found_a_v, _ = utility.search_text_file(fp_log, access_viol_str)

        #   Check if an exit number was found
        if found_e_n:

            #   Output the exit number
            exit_number = utility.find_int_in_str(e_n)

            #   Exit the loop
            break

        elif found_a_v:

            #   Set the exit number to the loss of connection to the license server exit number
            exit_number = 67

            #   Exit the loop
            break

        elif t_out >= t_max:

            exit_number = 0
            py_send("*kill_job *monitor_job")


        #   Wait and check again
        else:
            utility.wait(t, "exit number to be found")

    #   Check if the exit number indicates a successful run
    if exit_number == 3004:

        #   Set the success flag
        success = True

        m_log.info("Unit run successfully")

    #   Check if the exit number indicates a loss of connection to the license server
    elif exit_number == 67:

        m_log.error("License server connection timed out or failed")

        # #   Loop until a valid decision is made
        # while not decided:

        #     #   Request a user response
        #     dec = input("Would you like to run again? (y/n) ")

        #     #   Check if the response was yes
        #     if dec == "y":

        #         #   Set flag that a decision was made
        #         decided = True

        #   Rerun the job
        modify.run_job(j_id)

        #   Check if the updated output files exist
        success = check_out(fp_mud, fp_log, fp_t16, j_id)

            # #   Check if the response was no
            # elif dec == "n":

            #     #   Set flag that a decision was made
            #     decided = True

            #     m_log.warning("Results cannot be analysed")

            # #   Check if response was invalid
            # else:
            #     m_log.error("Invalid input received")
            #     print("Please either type a single \"y\" for yes or \"n\" for no")

    elif exit_number == 0:

        en_log.error("Model took too long to run. Job terminated")

    #   Output a warning
    else:
        en_log.error("Model run unsuccessfully with exit number {}".format(exit_number))
        m_log.warning("Results cannot be analysed. Check Mentat log file and exit number for details")

    #   Check if the unit was run successfully without a loss of connection to the license server
    if success and not decided:

        #   Wait until the t16 file exists and has been updated
        utility.wait_file_exist(fp_t16, "t16", t)
        utility.wait_file_update(fp_t16, t_mud, "t16", t)

        t1 = time.time()
        m_log.info("Results generated in approximately {:.3f}s".format(t1 - t0))
    
    #   Check if the unit was run successfully with a loss of connection to the license server
    elif success and decided:

        m_log.info("Results generated after initial connection failure")

    #   Output an error message
    else:

        t1 = time.time()
        m_log.warning("Results failed to generate after approximately {:.3f}s".format(t1 - t0))
        
    return success

################################################################################

def all_val(
    template,
    l: str,
    fp_t16: str,
    ) -> None:
    """Obtain values for all nodes from results

    Parameters
    ----------
    template 
        The unit template parameters
    l : str
        The label for the results file
        Either a template or unit identifier
    fp_t16 : str
        The file path of the model t16 file
    """

    #   The labels of the desired results
    label = []
    label.append("Displacement X")
    label.append("Displacement Y")
    label.append("Reaction Force X")
    label.append("Reaction Force Y")
    label.append("Total Strain Energy Density")
    label.append("Comp 11 of Total Strain")
    label.append("Comp 22 of Total Strain")
    label.append("Displacement")
    label.append("Reaction Force")

    #   Open the results file
    py_send("@main(results) @popup(modelplot_pm) *post_open \"{}\"".format(fp_t16))
    py_send("*post_numerics")

    #   Loop through all given labels
    for i in range(0, len(label)):

        #   Initialise list for the current label
        v = []

        #   Rewind the post file to the initial step
        py_send("*post_rewind")

        #   Set the post file to the current label
        py_send("*post_value {}".format(label[i]))

        #   Loop through all steps of the post file
        for j in range(0, template.n_steps + 1):

            #   Append an empty list to create a new index for every step
            v.append([])

            #   Loop through all 
            for k in range(1, template.n_n + 1):

                #   Append the current node's value to the list at the current step
                v[j].append(py_get_float("scalar_1({})".format(k)))

            #   Increment the post file
            py_send("*post_next")

        #   Save the values to a .csv file
        save_2d_list_to_csv(template, label[i] + "_" + l, v)

    #   Rewind and close the post file
    py_send("*post_rewind")
    py_send("*post_close")

    return

################################################################################

def internal_edges(
    unit_p,
    ) -> list:
    """Find all internal edges

    Parameters
    ----------
    unit_p : unit_p
        The unit parameters class object

    Returns
    -------
    list
        The internal edges sorted in order of boundary conditions to be applied
    """    

    #   Initialisations
    ip_yp = []
    ip_yn = []
    ip_xp = []
    ip_xn = []

    #   Format the grid for the search
    grid_flat = list(reversed(unit_p.grid))
    grid_flat = list(itertools.chain.from_iterable(grid_flat))

    #   Format the element indices of the elements removed
    rem_i = [i - 1 for i in unit_p.rem]

    #   Loop through all removed elements
    for i in rem_i:

        #   Check if the element above has not been removed
        if grid_flat[i + unit_p.template.x_e] == 1:

            #   Add the element to the appropriate list
            ip_yp.append(i + unit_p.template.x_e + 1)

        #   Check if the element below has not been removed
        if grid_flat[i - unit_p.template.x_e] == 1:

            #   Add the element to the appropriate list
            ip_yn.append(i - unit_p.template.x_e + 1)

        #   Check if the element to the right has not been removed
        if grid_flat[i + 1] == 1:

            #   Add the element to the appropriate list
            ip_xp.append(i + 1 + 1)

        #   Check if the element to the left has not been removed
        if grid_flat[i - 1] == 1:

            #   Add the element to the appropriate list
            ip_xn.append(i - 1 + 1)

    #   Combine all lists in the appropriate order
    ip_e = [ip_yp, ip_xn, ip_yn, ip_xp]

    return ip_e

################################################################################

def read_val(
    template,
    l: str,
    ) -> float:
    """Read a specific value from a results file

    Parameters
    ----------
    template 
        The unit template parameters
    l : str
        The label for the results file
        Either a template or unit identifier

    Returns
    -------
    float
        The found value
    """

    #   Create the file path of the results file
    fp_r_f = create_fp_file(template, l, "r")

    #   Read the value as a float if possible
    try:
        val = float(linecache.getline(fp_r_f, template.n_steps + 1))
    except:
        val = None

    return val

################################################################################

def read_xym(
    template,
    l: str,
    t: str,
    ) -> list:
    """Read the X, Y, and magnitude values for a type of result

    Parameters
    ----------
    template
        The unit template parameters
    l : str
        The label for the results file
        Either a template or unit identifier
    t : str
        The type of result

    Returns
    -------
    list
        The list of result values
    """    

    #   Initialisations
    label = []
    label.append(t + " X_")
    label.append(t + " Y_")
    label.append(t + "_")

    val = []

    #   Loop through the list of labels
    for i in label:

        #   Add the result value to the list
        val.append(read_val(template, i + l + "_1"))

    return val
    
################################################################################

def read_all(
    template,
    lu: list,
    l: str,
    ) -> list:
    """Read all result values from a list of units and plot them

    Parameters
    ----------
    l : str
        The result label
    lu : list
        The list of units
    template : template
        The unit template parameters

    Returns
    -------
    list
        The list of values
    """    

    #   Initialisations
    v = []

    #   Loop through all units
    for i in range(0, len(lu)):

        #   Initialisations
        v.append([])

        #   Create the file path of the results file
        fp_r_f = create_fp_file(template, l + "_" + lu[i] + "_2", "r")

        #   Store the results from the file
        v[i] = linecache.getline(fp_r_f, template.n_steps + 1)

    #   Clean the list
    v = [i.rstrip() for i in v]

    #   Cast the list items to floats and check for any failed results
    (v, v_f) = utility.list_to_float(v)

    #   Output a warning if any models failed to deliver results
    if v_f != 0:
        print("Warning: {} models failed to deliver results!".format(v_f))

    v = list(map(abs, v))

    return v

################################################################################

def read_all_hd(
    template,
    lu: list,
    ) -> list:

    #   Initialisations
    hd = []

    #   Loop through all units
    for i in range(0, len(lu)):

        #   Create the file path of the results file
        fp_r_f = create_fp_file(template, "Hausdorff Distance_" + lu[i] + "_2", "r")

        #   Store the results from the file
        hd.append(linecache.getline(fp_r_f, 1))

    #   Clean the list
    hd = [i.rstrip() for i in hd]

    #   Cast the list items to floats and check for any failed results
    (hd, hd_f) = utility.list_to_float(hd)

    #   Output a warning if any models failed to deliver results
    if hd_f != 0:
        print("Warning: {} models failed to deliver results!".format(hd_f))

    hd = list(map(abs, hd))

    return hd

################################################################################

def read_all_fit(
    template,
    lu: list,
    ) -> list:

    if template.case == 1:

        #   Initialisations
        r = []
        w = []

        #   Loop through all units
        for i in range(0, len(lu)):

            #   Create the file path of the results file
            fp_r_f = create_fp_file(template, "Case 1 Fitness Measures_" + lu[i] + "_2", "r")

            #   Store the results from the file
            r.append(linecache.getline(fp_r_f, 1))
            w.append(linecache.getline(fp_r_f, 2))

        #   Clean the list
        r = [i.rstrip() for i in r]
        w = [i.rstrip() for i in w]

        #   Cast the list items to floats and check for any failed results
        (r, r_f) = utility.list_to_float(r)
        (w, w_f) = utility.list_to_float(w)  

        #   Output a warning if any models failed to deliver results
        if r_f != 0:
            print("Warning: {} models failed to deliver results!".format(r_f))
        if w_f != 0:
            print("Warning: {} models failed to deliver results!".format(w_f))

        r = list(map(abs, r))
        w = list(map(abs, w))

        fit_m = [r, w]

    return fit_m

################################################################################

def read_lu(fp_lu: str) -> list:
    """Read the list of units

    Parameters
    ----------
    fp_lu : str
        The file path of the log file of units created during the last simulation

    Returns
    -------
    list
        The list of units
    """

    lu = []

    #   Read the list from the file
    for i in open(fp_lu, "r"):

        i = i.replace('"', '').strip()

        lu.append(i)

    #   Clean the list
    lu = [i.rstrip() for i in lu]

    return lu

################################################################################

def save_2d_list_to_csv(
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
    with open(fp_r_f, 'w') as f:
        wr = csv.writer(f)
        wr.writerows(data)

    print("{}.csv saved".format(l))

    return