"""
xqmaps,py

Define transforms that return an array of x and an array of Q²(GeV²) points
given some input kinematics. The arrays have to be of the same length, but not
the same as the input arrays.
"""

def DIS(k1, k2, k3, **extra):
    return k1, k2