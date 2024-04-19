import ctypes as ct
import numpy as np
from enum import Enum
import sys
import os

class OS(Enum):
    WIN = 1
    LINUX = 2
    MAC = 3
    UNKNOWN = 4

def determine_os():
    """
    Get the operating system of the current machine.
    """
    import platform
    if platform.system() == "Windows":
        return OS.WIN
    elif platform.system() == "Linux":
        return OS.LINUX
    elif platform.system() == "Darwin":
        return OS.MAC
    return OS.UNKNOWN

def get_shared_lib_name():
    if determine_os() == OS.WIN:
        return "ausaxs.dll"
    elif determine_os() == OS.LINUX:
        return "libausaxs.so"
    elif determine_os() == OS.MAC:
        return "libausaxs.dylib"
    return ""

path = os.path.join(sys.argv[1], get_shared_lib_name())
ausaxs = ct.CDLL(str(path))
ausaxs.evaluate_sans_debye.argtypes = [
    ct.POINTER(ct.c_double), # q vector
    ct.POINTER(ct.c_double), # x vector
    ct.POINTER(ct.c_double), # y vector
    ct.POINTER(ct.c_double), # z vector
    ct.POINTER(ct.c_double), # w vector
    ct.c_int,                # nq (number of points in q)
    ct.c_int,                # nc (number of points in x, y, z, w)
    ct.POINTER(ct.c_int),    # status (0 = success, 1 = q range error, 2 = other error)
    ct.POINTER(ct.c_double)  # Iq vector for return value
]
ausaxs.evaluate_sans_debye.restype = None # don't expect a return value

q = np.arange(1e-4, 1, 0.01)
x = np.array([0, 1, 2, 3, 4])
y = np.array([1, 2, 3, 4, 5])
z = np.array([2, 3, 4, 5, 6])
w = np.array([1, 1, 1, 1, 1])
nq = len(q)
nc = len(x)
status = ct.c_int(0)
Iq = np.zeros(nq)
ausaxs.evaluate_sans_debye(
    q.ctypes.data_as(ct.POINTER(ct.c_double)), 
    x.ctypes.data_as(ct.POINTER(ct.c_double)), 
    y.ctypes.data_as(ct.POINTER(ct.c_double)), 
    z.ctypes.data_as(ct.POINTER(ct.c_double)), 
    w.ctypes.data_as(ct.POINTER(ct.c_double)), 
    nq, 
    nc, 
    ct.byref(status), 
    Iq.ctypes.data_as(ct.POINTER(ct.c_double))
)
print("OK" if status.value == 0 else "ERROR")
sys.exit(status.value)