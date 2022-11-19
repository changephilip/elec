import numpy as np
import os

lib = np.ctypeslib.load_library("external.so",".")
#cdistance = getattr(lib,'calc_wrap')
lib.calc_wrap.restype = None
lib.calc_wrap.argtypes = (np.ctypeslib.ndpointer(np.float32, ndim=2, flags="aligned"),
                    np.ctypeslib.ndpointer(np.float32,ndim=2,flags="aligned"),
                    np.ctypeslib.ndpointer(np.float32,ndim=1,flags="aligned"),
                    np.ctypeslib.ndpointer(np.float32, ndim=1, flags="aligned"),
                    np.ctypeslib.ctypes.c_ulong,
                    np.ctypeslib.ctypes.c_ulong,
                    np.ctypeslib.ndpointer(np.float32,ndim=2,flags="aligned,contiguous,writeable"))

lib.euclidean_distance.restype = None
lib.euclidean_distance.argtypes = (np.ctypeslib.ndpointer(np.float32, ndim=2, flags="aligned"),
                    np.ctypeslib.ndpointer(np.float32,ndim=2,flags="aligned"),
                    np.ctypeslib.ndpointer(np.float32,ndim=1,flags="aligned"),
                    np.ctypeslib.ndpointer(np.float32, ndim=1, flags="aligned"),
                    np.ctypeslib.ctypes.c_ulong,
                    np.ctypeslib.ctypes.c_ulong,
                    np.ctypeslib.ndpointer(np.float32,ndim=2,flags="aligned,contiguous,writeable"))



if __name__ == "__main__":
    lib = np.ctypeslib.load_library("dielec.so",".")
    func1 = getattr(lib,'test')
    func1.restype = None
    func1.argtypes = (np.ctypeslib.ndpointer(np.float32 ,ndim=2,flags="aligned"),
                        np.ctypeslib.ctypes.c_ulong,
                        np.ctypeslib.ndpointer(np.float32 ,ndim=2, flags ="aligned,contiguous,writeable"))
    
    #a = np.arange(9).reshape(3,3)
    #a=a.astype('float32')
    #b = np.zeros_like(a)

    #func1(a,3,b)
    #print(b)

    cdistance = getattr(lib,'calc_wrap')
    cdistance.restype = None
    cdistance.argtypes = (np.ctypeslib.ndpointer(np.float32, ndim=2, flags="aligned"),
                        np.ctypeslib.ndpointer(np.float32,ndim=2,flags="aligned"),
                        np.ctypeslib.ndpointer(np.float32,ndim=2,flags="aligned"),
                        np.ctypeslib.ndpointer(np.float32, ndim=2, flags="aligned"),
                        np.ctypeslib.ctypes.c_ulong,
                        np.ctypeslib.ctypes.c_ulong,
                        np.ctypeslib.ndpointer(np.float32,ndim=2,flags="aligned,contiguous,writeable"))
    
    a = np.arange(12).reshape(4,3)
    a = a.astype('float32')
    b = np.arange(9).reshape(3,3)
    b = b.astype('float32')
    c = np.zeros_like(np.arange(12).reshape(4,3)).astype('float32')
    am = np.array([0.0,0.0,0.0]).astype('float32').reshape(3,1)
    bm = np.array([0.0,0.0,0.0]).astype('float32').reshape(3,1)

    cdistance(a,b,am,bm,4,4,c)
    print(c)

