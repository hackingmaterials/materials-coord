"""
This module contains supplemental methods used by CNBase derived classes in materialscord.cn_methods
"""

def Brunner(structure, n, mode="reciprocal", tol=1.0e-4, radius=8.0):
    """
    Helper function to compute Brunner's reciprocal gap and realspace gap CN.
    """
    nl = structure.get_neighbors(structure.sites[n], radius)
    ds = [i[-1] for i in nl]
    ds.sort()

    if mode == "reciprocal":
        ns = [1.0/ds[i] - 1.0/ds[i+1] for i in range(len(ds) - 1)]
    elif mode == "relative":
        ns = [ds[i]/ds[i+1] for i in range(len(ds) - 1)]
    elif mode == "real":
        ns = [ds[i] - ds[i+1] for i in range(len(ds) - 1)]
    else:
        raise ValueError("Unknown Brunner CN mode.")

    d_max = ds[ ns.index(max(ns)) ]
    cn = {}
    for i in nl:
        if i[-1] < d_max + tol:
            el = i[0].species_string
            if el in cn:
                cn[el] += 1.0
            else:
                cn[el] = 1.0
    return cn


def getDict(stringIn):
    """
    Helper function for Delaunay CNs by David Mrdjenovich et al. (LBL)
    """
    if stringIn[0] == '{' :
        stringIn = stringIn[1:len(stringIn)]
    if stringIn[len(stringIn) - 1] == '}' :
        stringIn = stringIn[0:len(stringIn) - 1]
    entries = stringIn.split(",")
    toReturn = {}
    for s in entries :
        keyVal = s.split(":")
        toReturn[keyVal[0]] = int(keyVal[1])
    return toReturn