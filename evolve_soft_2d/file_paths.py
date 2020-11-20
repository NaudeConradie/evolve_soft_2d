##  File paths

# astroid==2.3.3
# bleach==3.1.4
# certifi==2020.4.5.1
# chardet==3.0.4
# colorama==0.4.3
# cycler==0.10.0
# docutils==0.16
# idna==2.9
# importlib-metadata==1.6.0
# isort==4.3.21
# keyring==21.2.0
# kiwisolver==1.2.0
# lazy-object-proxy==1.4.3
# matplotlib==3.2.1
# mccabe==0.6.1
# mpmath==1.1.0
# numpy==1.15.1
# pandas==1.0.3
# pkginfo==1.5.0.1
# pyDOE==0.3.8
# Pygments==2.6.1
# pylint==2.4.4
# pyparsing==2.4.7
# python-dateutil==2.8.1
# pytz==2019.3
# pywin32-ctypes==0.2.0
# readme-renderer==25.0
# requests==2.23.0
# requests-toolbelt==0.9.1
# scipy==1.4.1
# seaborn==0.10.1
# six==1.14.0
# sympy==1.5.1
# tqdm==4.45.0
# twine==3.1.1
# typed-ast==1.4.1
# urllib3==1.25.9
# webencodings==0.5.1
# wrapt==1.11.2
# zipp==3.1.0

#   Imports
from evolve_soft_2d import utility

################################################################################

#   Fixed file paths
# fp = r'C:\Users\Naude Conradie\Desktop\Repository\Masters-Project\Models\MarcMentat'
# fp = r'C:\Users\19673418\Documents\Masters-Project\Models\MarcMentat'
fp = r'D:\Models'
# fp = r'G:\Models'
fp_g = fp + r'\Graphs'
fp_r = fp + r'\Results'
fp_t = fp + r'\Templates'
fp_u = fp + r'\Units'
fp_v = fp + r'\Variables'

################################################################################

def create_fp_folder(
    template,
    p: str,
    ) -> str:
    """Create the file path of the results folder up to the template level

    Parameters
    ----------
    template : class
        The unit template parameters
    p : str
        The specific folder
        "g" - The graph files
        "r" - The results folder
        "t" - The unit template folder
        "u" - The unit folder
        
    Returns
    -------
    str
        The folder path
    """

    #   Generate the folder name according to the current case and template
    fp_folder = r'\grid_' + template.t_id

    #   Determine which folder path to generate the folder name along
    if p == "g":
        fp_folder = fp_g + fp_folder
    elif p == "r":
        fp_folder = fp_r + fp_folder
    elif p == "t":
        fp_folder = fp_t + fp_folder
    elif p == "u":
        fp_folder = fp_u + fp_folder
    elif p == "v":
        fp_folder = fp_v + fp_folder

    #   Create the folder if it does not exist
    utility.make_folder(fp_folder)

    return fp_folder

################################################################################

def create_fp_file(
    template,
    l: str,
    p: str,
    unit = None,
    ) -> str:
    """Create the file path of the desired file

    Parameters
    ----------
    template : class
        The unit template parameters
    l : str
        The unique label of the file
    p : str
        The specific file path
        "l" - The log file of units generated during the current run of the programme
        "g" - The graph files
        "r" - The results files
        "t" - The unit template files
        "u" - The unit files
    unit : class, optional
        The unit template parameters, by default None

    Returns
    -------
    str
        The file path
    """    
    
    #   Determine which folder path to generate the file name along
    if p == "g":
        fp_folder = create_fp_folder(template, p)
        fp_file = fp_folder + "\\" + l + ".png"
    if p == "l":
        fp_file = fp_u + r'\grid_' + template.t_id + l + ".log"
    elif p == "r":
        fp_folder = create_fp_folder(template, p)
        fp_file = fp_folder + "\\" + l + ".csv"
    elif p == "t":
        fp_folder = create_fp_folder(template, p)
        fp_file = fp_folder + r'\grid_' + template.t_id + l
    elif p == "u":
        fp_folder = create_fp_folder(template, p)
        fp_file = fp_folder + r'\grid_' + unit.u_id

        #   Create the folder if it does not exist
        utility.make_folder(fp_file)
        
        fp_file = fp_file + r'\grid_' + unit.u_id + l
    elif p == "v":
        fp_folder = create_fp_folder(template, p)
        fp_file = fp_folder + "\\" + l

    return fp_file