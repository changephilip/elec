import numpy as np
import os

if __name__ == "__main__":
    lib = np.ctypeslib.load_library("dielec.so",".")
    func1 = getattr(lib,'test')
    func1.restype = None
    func1.argtypes = (np.ctypeslib.ndpointer(np.float32 ,ndim=2,flags="aligned"),
                        np.ctypeslib.ctypes.c_ulong,
                        np.ctypeslib.ndpointer(np.float32 ,ndim=2, flags ="aligned,contiguous,writeable"))
    
    a = np.arange(9).reshape(3,3)
    a=a.astype('float32')
    b = np.zeros_like(a)

    func1(a,3,b)
    print(b)

