##  Logging settings

#   Imports
from evolve_soft_2d.file_paths import fp

import logging

################################################################################

#   File paths for logs
m_log_path = fp + r'\units.log'
en_log_path = fp + r'\exit_number.log'

#   Basic logging configuration settings
logging.basicConfig(level = logging.DEBUG)
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

#   Logging handler for terminal
t_han = logging.StreamHandler()
t_han.setLevel(logging.WARNING)
t_for = logging.Formatter("%(name)s - %(levelname)s: %(message)s")
t_han.setFormatter(t_for)

#   Format for messages logged to file
f_for = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s: %(message)s", datefmt = "%d %b %Y - %H:%M:%S")

#   Logging handler for log file for units
m_f_han = logging.FileHandler(m_log_path)
m_f_han.setLevel(logging.DEBUG)
m_f_han.setFormatter(f_for)

#   Logging handler for log file for exit numbers
en_f_han = logging.FileHandler(en_log_path)
en_f_han.setLevel(logging.DEBUG)
en_f_han.setFormatter(f_for)

#   Logger for units
m_log = logging.getLogger("Unit Logger")
m_log.addHandler(t_han)
m_log.addHandler(m_f_han)

#   Logger for exit numbers
en_log = logging.getLogger("Exit Number Logger")
en_log.addHandler(t_han)
en_log.addHandler(en_f_han)