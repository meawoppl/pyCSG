import itertools
import numpy as np

#CSG Library for Python
# Modeled after: https://github.com/evanw/csg.js/blob/master/csg.js

def unit(array):
    return array / np.linalg.norm(array)

def lerp(b, a, t):
    return b + ((a - b) * t)

class Vertex(object):
    def __init__(self, position, normal):
        self.pos = np.array(position, dtype=np.float64, copy=True).squeeze()
        self.normal = unit( np.array(normal, dtype=np.float64, copy=True).squeeze() )

    def flip(self):
        self.normal *= -1

    def interpolate(self, other, t):
        return Vertex(lerp(self.pos, other.pos, t), lerp(self.normal, other.normal, t))

    def copy(self):
        return Vertex(self.pos.copy(), self.normal.copy())

    def __str__(self):
        return "pyCSG Vertex:" + str(self.pos)

    
class Plane(object):
    PLANE_eps = 1e-5
    def __init__(self, normal, w):
        self.normal = np.array(normal, dtype=np.float64, copy=True).squeeze()
        self.w = np.array(w, dtype=np.float64, copy=True).squeeze()

    @staticmethod
    def fromPoints(a, b, c):
        normal = unit( np.cross(b-a, c-a) ) 
        return Plane(normal, np.dot( normal, a ) )

    def flip(self):
        self.normal *= -1

    def copy(self):
        return Plane(self.normal, self.w)
     
    def splitPolygon(self, polygon):
        COPLANAR = 0
        FRONT = 1
        BACK = 2
        SPANNING = 3

        polyType = 0
        vertTypes = []
         
        for vert in polygon.vertices:
            t = np.dot( self.normal, vert.pos) - self.w
            typ = COPLANAR
            if t < -self.PLANE_eps:
                typ = BACK
            elif t > self.PLANE_eps:
                typ = FRONT
            else:
                typ = COPLANAR
                
            polyType |= typ
            vertTypes.append(typ)

        # Return lists . . .
        cf, cb, f, b = [], [], [], []

        # If the polygon is in the back or front of a plane, append it as such
        if   polyType == FRONT:
            f.append(polygon)
        elif polyType == BACK:
            b.append(polygon)
        elif polyType == COPLANAR:
            if np.dot(self.normal, polygon.plane.normal) > 0:
                cf.append(polygon)
            else:
                cb.append(polygon)
                
        elif polyType == SPANNING:
            newFront = []
            newBack  = []

            vtxCount = len(polygon.vertices)
            for vertIndx1 in range(vtxCount):
                vertIndx2 = (vertIndx1 + 1) % vtxCount
        
                # Get the front/back/coplanar type of the vertex
                ti = vertTypes[vertIndx1]
                tj = vertTypes[vertIndx2]

                # Get the assocaited vertices
                vi = polygon.vertices[vertIndx1]
                vj = polygon.vertices[vertIndx2]

                # If the vertex is in front or back, stick it there.
                if ti != BACK: newFront.append(vi.copy())
                if ti != FRONT: newBack.append(vi.copy())
                
                # If the current edge spans the cut plane
                if (ti | tj) == SPANNING:
                    # Calculate where it hits it
                    t = (self.w - self.normal.dot(vi.pos)) / (self.normal.dot(vj.pos - vi.pos))
                    v = vi.interpolate(vj, t)
                    
                    print "The edge between", vi
                    print "\tand", vj
                    print "\tis cut.", v
                    
                    # Add that new vertex to the front and back vertex lists
                    newFront.append( v.copy() )
                    newBack.append(  v.copy() )

            if len(newFront) >= 3: 
                f.append(Polygon(newFront, polygon.shared))
            if len(newBack) >= 3: 
                b.append(Polygon(newBack, polygon.shared))

        return cf, cb, f, b

class Polygon(object):
    def __init__(self, vertices, shared=None):
        self.vertices = vertices
        self.shared = shared
        self.plane = Plane.fromPoints(vertices[0].pos, vertices[1].pos, vertices[2].pos)

    def fromArray(ndarray):
        vertices = [Vertex(tri) for tri in ndarray]
        return Polygon(vertices)

    def copy(self):
        return Polygon(vertices.copy(), shared.copy())

    def flip(self):
        # Flip the order (and hence computed normal)
        self.vertices = self.vertices[::-1]
        
        # Flip the per/vertex normal as well
        for  v in self.vertices: v.flip()

    def __str__(self):
        return "pyCSG Polygon with %i vertices" % len(self.vertices)


class UniqueVertexCollector(object):
    def __init__(self):
        self.vertexIndex = {}
        self.indexVertex = {}
        self.indexCount = 0

    def addVertex(self, vertex):
        stringRep = vertex.tostring()
        if stringRep in self.vertexIndex:
            return
        
        self.vertexIndex[stringRep] = self.indexCount
        self.indexVertex[self.indexCount] = vertex
        self.indexCount += 1

    def getVertex(self, index):
        return self.indexVertex[index]

    def getIndex(self, vertex):
        return self.vertexIndex[vertex.tostring()]

    def __len__(self):
        return len(self.vertexIndex)


class PolygonMesh(object):
    def __init__(self, polygonList):
        self.polygons = polygonList

    # TODO!
    def writeSTL(self, filename):
        pass
    def writePLY(self, filename):
        # Go through, and determine the unique vertices
        vc = UniqueVertexCollector()
        for polygon in self.polygons:
            for vertex in polygon.vertices:
                vc.addVertex(vertex.pos)


        vertexCount = len(vc)
        faceCount   = len(self.polygons)
        # Open the file and stary writing the header
        f = open(filename, "w")
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("element vertex %i\n" % vertexCount )

        for axName in ["x", "y", "z"]:
            f.write("property float32 %s\n" % axName)
    
        f.write("element face %i\n" % faceCount )
        f.write("property list uint8 int32 vertex_index\n")
        f.write("end_header\n")

        # Write the vertices
        for indx in xrange(len(vc)):
            f.write( "%f %f %f\n" % tuple(vc.getVertex(indx)) )
        
        # Write the indices
        for poly in self.polygons:
            nVert = len(poly.vertices)
            indxs = [vc.getIndex(vert.pos) for vert in poly.vertices]
            f.write(str(nVert) + " " + " ".join([str(i) for i in indxs]) + "\n")


        print "Wrote {0} vertices, and {1} faces to file '{2}'".format(vertexCount, faceCount, filename)
        f.close()


    def __add__(self, otherMesh):
        t1 = BSPNode(   self.polygons + []    )
        t2 = BSPNode( otherMesh.polygons + [] )
        print "len0", len(t1.allPolygons()), len(t2.allPolygons())

        t1.clipTo(t2)
        print "len1", len(t1.allPolygons()), len(t2.allPolygons())

        # t2.clipTo(t1)
        # print "len3", len(t1.allPolygons()), len(t2.allPolygons())

        # t2.invert()
        # print "len4", len(t1.allPolygons()), len(t2.allPolygons())

        # t2.clipTo(t1)
        # print "len5", len(t1.allPolygons()), len(t2.allPolygons())

        # t2.invert()
        # print "len6", len(t1.allPolygons()), len(t2.allPolygons())
        
        # t1.build(t2.allPolygons())
        # print "len7", len(t1.allPolygons()), len(t2.allPolygons())

        return PolygonMesh(t1.allPolygons())


class BSPNode(object):
    def __init__(self, polygons=[]):
        self.plane = None
        self.front = None
        self.back  = None
        self.polygons = []

        if len(polygons) != 0: self.build(polygons)

    def invert(self):
        # Flip all the polygons on this node
        for poly in self.polygons: poly.flip()

        # Flip my own plane
        self.plane.flip()

        # Invert children
        if self.front: self.front.invert()
        if self.back : self.back.invert()
        
        # Invert front and back branches
        self.front, self.back = self.back, self.front


    def allPolygons(self):
        p = []
        if self.front: p += self.front.allPolygons()
        if self.back:  p += self.back.allPolygons()

        return  self.polygons + p

    
    def clipTo(self, otherBSP):
        self.polygons = otherBSP.clipPolygons(self.polygons)

        # Trim the front and back as well
        if self.front: self.front.clipTo(otherBSP)
        if self.back:  self.back.clipTo(otherBSP)


    def clipPolygons(self, polygons):
        # If a node lacks a splitting plane, then it is a leaf, so no clipping
        if not self.plane: return polygons[:]
        
        fr, bk = [], []
        for polygon in polygons:
            cFront, cBack, front, back =  self.plane.splitPolygon(polygon)
            fr += cFront + front
            bk += cBack  + back

        if self.front: 
            fr = self.front.clipPolygons(fr)
        if self.back : 
            bk = self.back.clipPolygons(bk)
        
        return fr + bk

        
    def build(self, polygons):
        if len(polygons) == 0: return

        if not self.plane: self.plane = polygons[0].plane.copy()

        newFront = []
        newBack  = []
        for polygon in polygons:
            cFront, cBack, front, back = self.plane.splitPolygon(polygon)
            self.polygons += cFront + cBack
            newFront += front
            newBack  += back
            
        if len(front) != 0 and not self.front:
            self.front = BSPNode()
            self.front.build(newFront)

        if len(back) != 0 and not self.back:
            self.back = BSPNode()
            self.back.build(newBack)


def makeCircle(center, radius, nTheta):
    pass


def makeTetrahedron(center):
    oort = 1.0 / np.sqrt(2)
    pts = np.array([[+1.0,  0.0, -oort],
                    [-1.0,  0.0, -oort],
                    [ 0.0, -1.0,  oort],
                    [ 0.0, +1.0,  oort]])
    # Make it axis centered in x, y, z
    mean_xyz = pts.mean(axis=0)
    pts -= mean_xyz
    pts += np.array(center)

    print "Tet points:", pts

    tetVerts = []
    for pt in pts:
        ver = Vertex(pt, unit(pt))
        tetVerts.append(ver)

    polys = []
    polys.append( Polygon( [tetVerts[0], tetVerts[2], tetVerts[1]] ) )
    polys.append( Polygon( [tetVerts[1], tetVerts[2], tetVerts[3]] ) ) 
    polys.append( Polygon( [tetVerts[0], tetVerts[1], tetVerts[3]] ) ) 
    polys.append( Polygon( [tetVerts[0], tetVerts[3], tetVerts[2]] ) ) 

    return PolygonMesh(polys)




if __name__ == "__main__":
    t1 = makeTetrahedron(np.array([0.0,0.0,0.0]))
    t2 = makeTetrahedron(np.array([0.0,0.0,-0.5]))

    t3 = t1 + t2

    t3.writePLY("unioned.ply")

    t4 = PolygonMesh( t1.polygons + t2.polygons )
    t4.writePLY("superimposed.ply")
