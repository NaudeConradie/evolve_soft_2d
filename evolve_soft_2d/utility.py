##  Utility functions

#   Imports
import hashlib
import math
import numpy
import os.path
import pandas
import pickle
import re
import time

from pathlib import Path
from scipy.special import comb

from evolve_soft_2d import file_paths

################################################################################

def wait(
    t: float,
    f: str,
    ) -> None:
    """Wait for a specified time

    Parameters
    ----------
    t : float
        The time in seconds to wait
    f : str
        The name of the object being waited for
    """

    print("Waiting for %s..." % f)

    time.sleep(t)

    return

################################################################################
  
def wait_file_exist(
    file_name: str,
    label: str,
    t: float,
    ) -> None:
    """Wait until a specified file exists

    Parameters
    ----------
    file_name : str
        The name of the file to be waited for
    label : str
        The label for the output message
    t : float
        The time in seconds to wait per loop
    """
    
    #   Loop until the file exists
    while 1:

        #   Check if the file exists
        exists = if_file_exist(file_name)

        #   Break out of the loop if the file exists
        if exists:

            break

        #   Wait and check again
        else:
            wait(t, "%s file to be created..." % label)

    return

################################################################################
  
def wait_file_update(
    file_name: str,
    t0: float,
    label: str,
    t: float,
    ) -> None:
    """Wait until a specified file is updated

    Parameters
    ----------
    file_name : str
        The name of the file to be waited for
    t0 : float
        The time since which the file should have been updated
    label : str
        The label for the output message
    t : float
        The time in seconds to wait per loop
    """

    #   Loop until the file has been updated
    while 1:

        #   See how recently the file has been updated
        t_mod = os.path.getmtime(file_name)

        #   Break out of the loop if the file has been updated since the given time
        if t_mod > t0:

            break

        #   Wait and check again
        else:
            wait(t, "%s file to be updated..." % label)

    return

################################################################################

def if_file_exist(file_name: str) -> bool:
    """Check if a specified file exists

    Parameters
    ----------
    file_name : str
        The name of the file to be checked

    Returns
    -------
    bool
        True if the file exists, False otherwise
    """

    exists = os.path.exists(file_name)

    return exists

################################################################################

def make_folder(l: str) -> None:
    """Make a folder if it does not exist

    Parameters
    ----------
    l : str
        The folder name and path
    """

    Path(l).mkdir(parents=True, exist_ok=True)

    return

################################################################################

def f_const(
    x,
    const,
    ):
    """A constant function

    Parameters
    ----------
    x
        x data
    const
        A constant

    Returns
    -------
    
        A constant
    """    

    return const

################################################################################

def sigmoid(x: numpy.array) -> numpy.array:
    """Sigmoid function

    Parameters
    ----------
    x : numpy.array
        The input data

    Returns
    -------
    numpy.array
        The output data
    """

    s = 1/(numpy.exp(-x) + 1)

    return s

################################################################################

def smooth_relu(x: numpy.array) -> numpy.array:
    """Smooth ReLu function

    Parameters
    ----------
    x : numpy.array
        The input data

    Returns
    -------
    numpy.array
        The output data
    """

    s_r = numpy.log(numpy.exp(x) + 1)

    return s_r

################################################################################

def nCr(
    n: int,
    r: int = 0,
    l: list = [],
    ) -> int:
    """Calculate the n Choose r value for a single, multiple or all r

    Parameters
    ----------
    n : int
        The number to choose from
    r : int, optional
        The numbers to choose, by default 0
    l : list, optional
        The list of numbers to choose, by default []

    Returns
    -------
    int
        The total possible choices
    """

    #   Initialisations
    c = 0

    #   Check if the list is not empty
    if l != []:

        #   Loop through the list
        for i in l:

            #   Cumulatively calculate the n Choose r value for each list item
            c += comb(n, i, exact = True)

    #   Check if r is not zero
    elif r != 0:

        #   Calculate the n Choose r value
        c = comb(n, r, exact = True)

    else:

        #   Loop through all possible r
        for i in range(0, n):

            #   Cumulatively calculate the n Choose r value
            c += comb(n, i, exact = True)

    return c

################################################################################

def sel_random(
    l: list,
    f: int = None,
    seed: int = None,
    ) -> list:
    """Randomly select a random number of numbers from a given list of numbers

    Parameters
    ----------
    l : list
        The given list of numbers
    f : int, optional
        How many numbers to select, by default None
    seed : int, optional
        The seed for random generation, by default None

    Returns
    -------
    list
        The randomly selected numbers
    """

    #   Initialisations
    sel = []
    l_temp = l[:]

    #   Check if how many numbers to select is not zero
    if f is not None:

        #   Determine how many numbers to select from the given number
        n_sel = f

    else:

        if seed is not None:

            numpy.random.seed(seed = seed)

        #   Determine how many numbers to select from the given list
        n_sel = numpy.random.randint(low = 1, high = len(l_temp))

    #   Loop through the amount of numbers to be selected
    for i in range(0, n_sel):

        if seed is not None:

            numpy.random.seed(seed = seed + i)
        
        #   Select a random number from the list of numbers
        sel.append(numpy.random.choice(numpy.asarray(l_temp)))

        #   Remove the selected number from the list of numbers to prevent the same number from being selected more than once
        l_temp.remove(sel[i])

    #   Sort the list of selected numbers
    sel.sort()

    return sel

################################################################################

def gen_random(
    l: list,
    n: int,
    seed: int = None
    ) -> str:
    """Generate a random string from a list of characters

    Parameters
    ----------
    l : list
        The list of characters
    n : int
        The length of the string

    Returns
    -------
    str
        The randomly generated string
    """

    #   Initialisations
    s = ""

    #   Loop for the desired length of the string
    for i in range(0, n):

        if seed is not None:

            numpy.random.seed(seed + i)

        #   Append the next random character
        s += numpy.random.choice(l)

    return s

################################################################################

def gen_hash(s: str) -> str:
    """Generate a random hash code from a given string

    Parameters
    ----------
    s : str
        The string to be hashed

    Returns
    -------
    str
        The hash code
    """

    m = hashlib.md5()
    m.update(bytes(s, encoding = 'utf8'))
    hash_code = str(m.hexdigest())

    return hash_code

################################################################################

def search_text_file(
    file_name: str,
    find_text: str,
    ) -> (bool, str):
    """Search a text file for the first occurrence of a given text string

    Parameters
    ----------
    file_name : str
        The name of the file to be searched through
    find_text : str
        The text to be searched for

    Returns
    -------
    (bool, str)
        True if the text was found, false otherwise
        The entire line containing the text
    """

    #   Initialisations
    found_text = ""
    found = False

    #   Open the file to be read
    with open(file_name, "rt") as f:

        #   Loop through every line in the file until the text string
        for line in f:

            #   Check if the text is in the current line of the file
            if find_text.search(line) != None:

                #   Save the entire line of text containing the desired text
                found_text = line.rstrip("\n")

                #   Set the found flag to true
                found = True

                #   Exit the loop
                break

    return (found, found_text)

################################################################################

def unique_list(
    l1: list,
    l2: list,
    ) -> list:
    """Find unique members among two lists

    Parameters
    ----------
    l1 : list
        The first list
    l2 : list
        The second list

    Returns
    -------
    list
        The unique list members
    """

    l = list((set(l1) | set(l2)) - (set(l1) & set(l2)))

    return l

################################################################################

def add_sort_list(
    l1: list,
    l2: list,
    ) -> list:
    """Add two lists and sort them

    Parameters
    ----------
    l1 : list
        The first list to be added
    l2 : list
        The second list to be added

    Returns
    -------
    list
        The added and sorted list
    """
    
    #   Add two lists
    l = l1 + l2

    #   Sort the added lists
    l.sort()

    return l

################################################################################

def normalise_list(
    l: list,
    d: int,
    f: int = 1,
    ) -> list:
    """Normalise a list according to a given displacement

    Parameters
    ----------
    l : list
        The list of values
    n : int, optional
        The displacement, by default 0

    Returns
    -------
    list
        The normalised list
    """    

    l = [i/f + d for i in l]

    return l

################################################################################

def list_swap_i(
    l: list,
    i1: int,
    i2: int,
    ) -> list:
    """Swap list values at specified indices

    Parameters
    ----------
    l : list
        The list
    i1 : int
        The first index
    i2 : int
        The second index

    Returns
    -------
    list
        The list with swapped indices
    """    

    l[i1], l[i2] = l[i2], l[i1]

    return l

################################################################################

def clean_list(
    l: list,
    ub: int,
    lb: int = 0,
    ) -> list:
    """Clean a list according to a given upper and lower bound

    Parameters
    ----------
    l : list
        The list to be cleaned
    ub : int
        The inclusive upper bound
    lb : int, optional
        The exclusive lower bound, by default 0

    Returns
    -------
    list
        The cleaned list
    """

    #   Initialisations
    l_temp = l[:]

    #   Replace all out-of-bounds values with the boundary values
    l_temp = [clean_int(i, ub, lb = lb) for i in l_temp]

    #   Sort the list and remove duplicates
    l_temp = list(set(l_temp))

    return l_temp

################################################################################

def clean_int(
    i: int,
    ub: int,
    lb: int = 0,
    ) -> int:
    """Clean an integer according to a given upper and lower bound

    Parameters
    ----------
    i : int
        The integer to be cleaned
    ub : int
        The inclusive upper bound
    lb : int, optional
        The exclusive lower bound, by default 0

    Returns
    -------
    int
        The cleaned integer
    """    

    #   Initialisations
    i_temp = i

    #   Check if the integer is above the upper bound
    if i_temp > ub:

        #   Set it to the upper bound
        i_temp = ub

    #   Check if the integer is below or equal to the lower bound
    elif i_temp <= lb:

        #   Set it to one above the lower bound
        i_temp = lb + 1

    return i_temp

################################################################################

def clean_str(
    s: str,
    l: list,
    r: list,
    ) -> str:
    """Replace substrings within a given string

    Parameters
    ----------
    s : str
        The string
    l : list
        The list of substrings
    r : list
        The list of replacement substrings

    Returns
    -------
    str
        The string with the substrings removed
    """    

    #   Loop through every substring in the list
    for i in range(0, len(l)):

        #   Remove all occurrences of the substring
        s = s.replace(l[i], r[i])

    return s

################################################################################

def list_to_float(l: list) -> [list, int]:
    """Extract a list of floats from a given list

    Parameters
    ----------
    l : list
        The list to be extracted from

    Returns
    -------
    [list, int]
        The list of floats
        The number of list of items that failed to convert to string
    """

    #   Initialisations
    l_o = []
    l_f = 0

    #   Loop through the list
    for i in l:

        #   Add the current list item as a float
        try:
            l_o.append(float(i))

        #   Increment the failure counter
        except:
            l_o.append(0)
            l_f += 1

    return l_o, l_f

################################################################################

def list_to_str(
    l: list,
    c: str,
    ) -> str:
    """Convert a list into a string connected by a given symbol

    Parameters
    ----------
    l : list
        The list to be converted
    c : str
        The symbol to be inserted between list items

    Returns
    -------
    str
        The output string
    """

    s = c.join(map(str, l))

    return s

################################################################################

def find_int_in_str(s: str) -> int:
    """Finds the first integer in a string

    Parameters
    ----------
    s : str
        The string to be searched

    Returns
    -------
    int
        The integer found
    """
    
    i = int(re.search(r'\d+', s).group())

    return i

################################################################################

def read_str(
    s: str,
    i1: int,
    i2: int,
    ) -> str:
    """Read a string segment

    Parameters
    ----------
    s : str
        The entire string
    i1 : int
        The initial index
    i2 : int
        The last index

    Returns
    -------
    str
        The string segment
    """    

    s_last = s[i1:i2]

    return s_last

################################################################################

def save_v(
    template,
    v,
    l: str,
    ) -> None:
    """Save a variable as a file

    Parameters
    ----------
    v : 
        The variable to be saved
    l : str
        The file path
    """

    fp = file_paths.create_fp_file(template, l, "v")

    #   Open the file to be written to
    f = open(fp, 'wb')

    #   Save the variable
    pickle.dump(v, f)

    #   Close the file
    f.close()

    return

################################################################################

def open_v(
    template,
    l: str,
    ):
    """Open a variable saved to a file

    Parameters
    ----------
    fp : str
        The file path

    Returns
    -------
    
        The saved variable
    """

    fp = file_paths.create_fp_file(template, l, "v")

    #   Open the file to be read from
    f = open(fp, "rb")

    #   Store the variable
    v = pickle.load(f)

    #   Close the file
    f.close()

    return v

################################################################################

def avg_str_l(l: list) -> float:

    l_len = []

    for i in l:

        l_len.append(len(i))

    avg = sum(l_len)/len(l_len)

    return avg

################################################################################

def coord_to_grid(c: pandas.DataFrame) -> numpy.array:

    #   Normalise the coordinates
    c.x = normalise_list(c.x, max(c.x))
    c.y = normalise_list(c.y, max(c.y))

    c = c.reset_index(drop = True)
    c = c.astype(int)

    grid = numpy.zeros((max(c.x) + 1, max(c.y) + 1))

    for _, r in c.iterrows():

        grid[r.x, r.y] = 1

    grid = array_flip_h(grid)

    grid = grid.T

    return grid

################################################################################

def array_flip_h(a):

    a_temp = a.copy()

    a_temp = numpy.flip(a_temp, axis = 1)

    return a_temp
