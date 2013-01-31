#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os, csg
import struct
import numpy as np


def loadABA(filename):
    fp = open(filename,"rb")

    # Read the vertex count
    vCount = struct.unpack("I", fp.read(4))[0]

    print "vertex Count:", vCount
    vertexDataSize = vCount * 4 * (3+3)

    rawArray = np.fromstring(fp.read(vertexDataSize), dtype=np.float32)
    xn = rawArray[0::6]
    yn = rawArray[1::6]
    zn = rawArray[2::6]
    
    xs = rawArray[3::6]
    ys = rawArray[4::6]
    zs = rawArray[5::6]

    # from pylab import *
    # plot(xs, zs, "bo")
    # show()
    # 1/0

    polygons = []
    polyCount = struct.unpack("I", fp.read(4))[0]
    for polygonIndex in xrange(polyCount):
        pointCount = struct.unpack("H", fp.read(2))[0]
        if polygonIndex % 100 == 0:
            print polygonIndex, polyCount

    
        vertices = []
        for pointIndex in xrange(pointCount):
            index = struct.unpack("I", fp.read(4))[0]
            vertices.append(csg.Vertex.fromXYZ(xs[index], ys[index], zs[index]))
            
        for indx in xrange(pointCount-2):
            even = (indx % 2 == 0)
            v1 = vertices[indx + 0]
            v2 = vertices[indx + 1]
            v3 = vertices[indx + 2]

            # Triange strip winding angle hoot-nanny
            if even:
                polygons.append(csg.Polygon([v1, v2, v3]))
            else:
                polygons.append(csg.Polygon([v2, v1, v3]))

    return csg.PolygonMesh(polygons)

