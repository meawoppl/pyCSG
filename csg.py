import itertools, md5
import numpy as np

# CSG Library for Python
# Modeled after: https://github.com/evanw/csg.js/blob/master/csg.js


def unit(array):
    return array / np.linalg.norm(array)

def lerp(b, a, t):
    return b + ((a - b) * t)

class Vertex(object):
    def __init__(self, position, normal):
        self.pos = np.array(position, dtype=np.float64, copy=True).squeeze()
        self.normal = unit( np.array(normal, dtype=np.float64, copy=True).squeeze() )

    @staticmethod
    def fromXYZ(x, y, z):
        pos = np.array([x,y,z])
        return Vertex(pos, np.ones(3))

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
        self.w *= -1

    def copy(self):
        return Plane(self.normal, self.w)
     
    def splitPolygon(self, polygon):
        '''Take a polygon in, and split it by my plane into coplanar and front/back categories.'''
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
            f.append( polygon )
        elif polyType == BACK:
            b.append( polygon )
        elif polyType == COPLANAR:
            if np.dot(self.normal, polygon.plane.normal) > 0:
                cf.append( polygon )
            else:
                cb.append( polygon )
                
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
                if ti != BACK: newFront.append( vi )
                if ti != FRONT: newBack.append( vi )
                
                # If the current edge spans the cut plane
                if (ti | tj) == SPANNING:
                    # Calculate where it hits it
                    t = (self.w - self.normal.dot(vi.pos)) / (self.normal.dot(vj.pos - vi.pos))
                    v = vi.interpolate(vj, t)
                                
                    # Add that new vertex to both the front and back vertex lists
                    newFront.append( v.copy() )
                    newBack.append(  v.copy() )

            if len(newFront) >= 3: 
                f.append(Polygon(newFront, polygon.shared))
            if len(newBack) >= 3: 
                b.append(Polygon(newBack, polygon.shared))

        return cf, cb, f, b


def planeCoalesce(polygons):
    edgeCount = {}

    vc = UniqueVertexCollector()

    for polygon in polygons:
        nVertices = len(polygon.vertices)
        for i in xrange(nVertices):
            j = (i + 1) % nVertices
            v1 = polygon.vertices[i]
            v2 = polygon.vertices[j]
            
            # Put the vertices them in a arbitrary to index 
            vi1 = vc.addVertex(v1.pos)
            vi2 = vc.addVertex(v2.pos)
            
            edge_tuple = tuple(sorted([vi1, vi2]))

            edgeCount[edge_tuple] = edgeCount.get(edge_tuple, 0) + 1




class Polygon(object):
    def __init__(self, vertices, shared={}):
        self.vertices = vertices + []
        self.shared = shared.copy()
        if "color" not in shared: 
            self.shared["color"] = (0.5, 0.5, 0.5) # Grey if otherwise unlabeled
        

        self.plane = Plane.fromPoints(vertices[0].pos, vertices[1].pos, vertices[2].pos)

        m = md5.new()
        m.update("".join([v.pos.tostring() for v in self.vertices]))
        self._hash = abs(long(m.hexdigest(), 16))

    def fromArray(ndarray):
        vertices = [Vertex(tri) for tri in ndarray]
        return Polygon(vertices)

    def copy(self):
        return Polygon(self.vertices + [], self.shared.copy())

    def flip(self):
        # Flip the order (and hence computed normal)
        self.vertices = self.vertices[::-1]
        
        # Flip the per/vertex normal as well
        for  v in self.vertices: v.flip()

    def __str__(self):
        printable = ["\t" + str(v.pos) for v in self.vertices]
        return "\n".join(["pyCSG Polygon with %i vertices" % len(self.vertices)] + printable )

    def __eq__(self, other):
        # Make sure a pair of polygons have the same point count
        if len(self.vertices) != len(other.vertices):
            return False
        
        # Make sure each of the points is the same
        for vtx1, vtx2 in zip(self.vertices, other.vertices):
            if not all(vtx1.pos == vtx2.pos):
                return False
        return True

    def __hash__(self):        
        return self._hash


class UniqueVertexCollector(object):
    def __init__(self):
        # Bidirectional dictionary.  
        # ndarray -> vertexID and vice versa
        self.vertexIndex = {}
        self.indexVertex = {}

        # Internal counter for vertex insertion
        self.indexCount = 0

    def addVertex(self, vertex):
        # ndarrays are not hashable, but their string rep is
        stringRep = vertex.tostring()
        if stringRep in self.vertexIndex:
            return
        
        # The bidirectional dict (2 dicts)
        self.vertexIndex[stringRep] = self.indexCount
        self.indexVertex[self.indexCount] = vertex
        self.indexCount += 1

        return self.indexCount - 1

    def getVertex(self, index):
        return self.indexVertex[index]

    def getIndex(self, vertex):
        return self.vertexIndex[vertex.tostring()]

    def __len__(self):
        return len(self.vertexIndex)

    @staticmethod
    def fromPolygonList(lst):
        vc = UniqueVertexCollector()
        for polygon in lst:
            for vertex in polygon.vertices:
                vc.addVertex(vertex.pos)
        return vc
        

class PolygonMesh(object):
    def __init__(self, polygonList = []):
        self.polygons = polygonList

    def writeOBJ(self, filename):
        m = open("colors.mtl", "w")

        # m.write(""
        # Figure out what all colors we are using
        def tupleToColorString(colorTuple):
            return "color" + "_".join(["%0.2f" % c for c in colorTuple])

        colorsWritten = set()
        for polygon in self.polygons:
            polygonColor = polygon.shared["color"]
            if polygonColor in colorsWritten:
                continue
            m.write("newmtl " + tupleToColorString(polygonColor) + "\n")
            m.write("Ka %f %f %f\n" % polygonColor)
            m.write("Kd %f %f %f\n" % polygonColor)
            m.write("illum 2\n")
            m.write("\n")
            colorsWritten.add(polygonColor)
        
        m.close()

        f = open(filename, "w")
        f.write("# Created by pyCSG\n")
        f.write("mtllib colors.mtl\n")
        
        # Get the unique vertices
        vc = UniqueVertexCollector.fromPolygonList(self.polygons)

        # Write the vertices (in order!)
        for indx in xrange(len(vc)):
            f.write("v %f %f %f\n" % tuple(vc.getVertex(indx)))

        # Write the face indices
        for polygon in self.polygons:
            # Tell it to use the right color
            f.write("usemtl " + tupleToColorString(polygon.shared["color"]) + "\n")
            # .obj is 1 based hence the +1
            vIndices = [str(vc.getIndex(v.pos) + 1) for v in polygon.vertices]
            f.write( "f " + " ".join(vIndices) + "\n" )

        f.close()
            
    # TODO!
    def writeSTL(self, filename):
        pass

    def writePLY(self, filename):
        # Go through, and determine the unique vertices
        vc = UniqueVertexCollector.fromPolygonList(self.polygons)

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
        for clr in ["red", "blue", "green"]:
            f.write("property uchar %s\n" % clr)

        f.write("end_header\n")

        # Write the vertices
        for indx in xrange(len(vc)):
            vtx = vc.getVertex(indx)
            f.write( "%f %f %f\n" % tuple( vtx ) )
        
        # Write the indices
        for poly in self.polygons:
            nVert = len(poly.vertices)
            indxs = [vc.getIndex(vert.pos) for vert in poly.vertices]
            count_nIndexes = str(nVert) + " " + " ".join([str(i) for i in indxs])
            # Cast the 0-1 color into uint8 land, and stringify
            colors = "%i %i %i" % tuple([int(clr*255) for clr in poly.shared.get("color", (0.5, 0.5, 0.5) )])
            f.write(count_nIndexes + " " + colors + "\n")


        print "Wrote {0} vertices, and {1} faces to file '{2}'".format(vertexCount, faceCount, filename)
        f.close()


    def __or__(self, otherMesh):
        # Build trees with copies of the two polygon sets
        t1 = BSPNode(   self.polygons  )
        t2 = BSPNode( otherMesh.polygons )
        t1.clipTo(t2)
        t2.clipTo(t1)
        t2.invert()
        t2.clipTo(t1)
        t2.invert()
        t1.build(t2.allPolygons())

        return PolygonMesh(t1.allPolygons())

    def __and__(self, otherMesh):
        # MRG TODO: not workign 100% correctly.
        t1 = BSPNode(   self.polygons    )
        t2 = BSPNode( otherMesh.polygons )
        t1.invert()
        t2.clipTo(t1)
        t2.invert()
        t1.clipTo(t2)
        t2.clipTo(t1)
        t1.build(t2.allPolygons())
        t1.invert()
        return PolygonMesh(t1.allPolygons())


    def __sub__(self, otherMesh):
        t1 = BSPNode(   self.polygons    )
        t2 = BSPNode( otherMesh.polygons )
        
        t1.invert()
        t1.clipTo(t2)
        t2.clipTo(t1)
        t2.invert()
        t2.clipTo(t1)
        t2.invert()
        t1.build(t2.allPolygons())
        t1.invert()

        return PolygonMesh(t1.allPolygons())


def union( *meshes ):
    print meshes

    bsps = [BSPNode(m.polygons) for m in meshes]


    for bsp1, bsp2 in itertools.combinations(bsps, 2):
        bsp1.clipTo(bsp2)
        bsp2.clipTo(bsp1)

        print bsp1, bsp2
    

    mergeBSP = bsps[0]
    for bsp in bsps[1:]:
        bsp.invert()
        bsp.clipTo(mergeBSP)
        bsp.invert()
        mergeBSP.build(bsp.allPolygons())

    return PolygonMesh(mergeBSP.allPolygons())



class BSPNode(object):
    def __init__(self, polygons=[]):
        self.plane = None
        self.front = None
        self.back  = None
        self.polygons =  [] 

        if len(polygons) != 0: self.build(polygons + [])

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
        '''Return a shallow copy to all the polygons upstream of this node.'''
        p = []
        if self.front: p += self.front.allPolygons()
        if self.back:  p += self.back.allPolygons()
        return  self.polygons + p

    
    def clipTo(self, otherBSP):
        '''Clip all the polygons in this tree with respect to another tree.'''
        # My new polygons are those clipped against the tree passed
        self.polygons = otherBSP.clipPolygons(self.polygons)

        # Trim the front and back as well
        if self.front: self.front.clipTo(otherBSP)
        if self.back:  self.back.clipTo(otherBSP)

    def clipPolygons(self, polygons):
        # If a node lacks a splitting plane, then it is a leaf, so no clipping
        if not self.plane: return polygons

        # A list that keys a single polygon to its child (split) polygons
        parentChildList = {}

        fr, bk = [], []
        # For each polygon
        for polygon in polygons:
            # Split the polygon by my plane
            cFront, cBack, front, back =  self.plane.splitPolygon(polygon)
            fr += cFront + front
            bk += cBack  + back

            # If a polygon gets split (becomes more than one polygon)
            if len(fr + bk) > 1: 
                parentChildList[polygon] = fr + bk
                print "Cut detected", hash(polygon)

        # The generated polygons that are in front/back of  me might
        # need to further be clipped by other obects in front of me
        if self.front: 
            fr = self.front.clipPolygons(fr)
        if self.back: 
            bk = self.back.clipPolygons(bk)
        else:
            # This is the actual clip.  
            # The polygons were partitioned above, 
            # but now there is nothing behind me in the BPS they can be discarded
            bk = []

        survivingPolygons = fr + bk
        # Lastly we check if all the child polygons are still in this bsp after slicing 
        # for parent, children in parentChildList.iteritems():
        #     figure()
        #     print hash(parent), [hash(child) for child in children]
        #     # Make sure that all of the child polygons are still around
        #     allChildrenAlive = all([child in survivingPolygons for child in children])
        #     if not allChildrenAlive: 
        #         parent.shared["color"] = (1,0,0)
        #         continue

        #     print "All Children Alive", hash(parent)
        #     for child in children: 
        #         survivingPolygons.remove(child)
        #         # for vtx in child.vertices:

        #     parent.shared["color"] = (0,1,0)
        #     survivingPolygons.append(parent)

        return survivingPolygons

        
    def build(self, polygons):
        '''Build up this bsp tree based on the passed polygons.'''
        if len(polygons) == 0: return

        # Assign this node a cutting plane if it does not already exist.
        if not self.plane: self.plane = polygons[0].plane.copy()

        # Collect polygons that go in back and front of this plane
        newFront = []
        newBack  = []
        for polygon in polygons:
            # Split (or just allocate) the polygon into coplanar/front/back parts 
            cFront, cBack, front, back = self.plane.splitPolygon(polygon)
            
            # Collect the coplanar ones on this node
            self.polygons += cFront + cBack

            # Pile up the front and back ones to push down the front/back BSP's
            newFront += front
            newBack  += back
            
        if len(newFront) != 0:
            # Init a new node if necessary
            if not self.front: self.front = BSPNode()
            # Push Polys up the tree
            self.front.build(newFront)

        if len(newBack) != 0:
            # Init a new node if necessary
            if not self.back: self.back = BSPNode()
            # Push Polys up the tree
            self.back.build(newBack)


def makeSphere(center, radius, nTheta, nPhi):
    phis   = np.linspace(-np.pi/2, np.pi/2, nPhi)
    thetas = np.linspace(0, 2*np.pi, nTheta)

    def computePos(theta, phi, radius):
        x = np.cos(theta) * np.cos(phi) * radius
        y = np.sin(theta) * np.cos(phi) * radius
        z = np.sin(phi) * radius
        return np.array([x, y, z], dtype=np.float64)

    clr = randomColor()
    polyCollection = []
    for phiIndx in range(nPhi)[:-1]:
        for thetaIndx in range(nTheta)[:-1]:
            theta1 = thetas[thetaIndx]
            phi1   = phis  [phiIndx]

            theta2 = thetas[thetaIndx+1]
            phi2   = phis  [phiIndx+1]

            p1 = computePos(theta1, phi1, radius) + center
            p2 = computePos(theta2, phi1, radius) + center
            p3 = computePos(theta1, phi2, radius) + center
            p4 = computePos(theta2, phi2, radius) + center

            v1 = Vertex( p1, unit(p1) )
            v2 = Vertex( p2, unit(p2)  )
            v3 = Vertex( p3, unit(p3)  )
            v4 = Vertex( p4, unit(p4)  )

            dup123 = (np.allclose(v1.pos, v2.pos) or np.allclose(v2.pos, v3.pos) or np.allclose(v1.pos, v2.pos)  )
            if not dup123:
                polyCollection.append( Polygon([v1, v2, v3], {"color":clr} ))
                
            dupe324 = (np.allclose(v3.pos, v2.pos) or np.allclose(v2.pos, v4.pos) or np.allclose(v3.pos, v4.pos)  )
            if not dupe324:
                polyCollection.append( Polygon([v3, v2, v4], {"color":clr} ))

    return PolygonMesh(polyCollection)


def randomColor():
    import random
    # MRG FIXME
    return (0.5, 0.5, 0.5) # tuple([random.uniform(0,1) for r in range(3)])


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

    clr = randomColor()
    polys = []
    polys.append( Polygon( [tetVerts[0], tetVerts[2], tetVerts[1]] , {"color":clr} ) )
    polys.append( Polygon( [tetVerts[1], tetVerts[2], tetVerts[3]] , {"color":clr} ) ) 
    polys.append( Polygon( [tetVerts[0], tetVerts[1], tetVerts[3]] , {"color":clr} ) ) 
    polys.append( Polygon( [tetVerts[0], tetVerts[3], tetVerts[2]] , {"color":clr} ) ) 

    return PolygonMesh(polys)


def testOperators(m1, m2, name):
    # sup = PolygonMesh(m1.polygons + m2.polygons)
    # sup.writePLY(name + "-super.ply")
    # sup.writeOBJ(name + "-super.obj")

    union = m1 | m2
    # union.writePLY(name + "-union.ply")
    union.writeOBJ(name + "-union.obj")

    # sub = m1 - m2
    # sub.writePLY(name + "-sub.ply")
    # sub.writeOBJ(name + "-sub.obj")

    # intersect = m1 & m2
    # intersect.writePLY(name + "-int.ply")
    # intersect.writeOBJ(name + "-int.obj")



def doVasculature(filename):
    xs, ys, zs, rs = np.load(filename).T

    sortr = rs.argsort()[::-1]

    xyz = np.c_[xs, ys, zs]

    xyz = xyz[sortr,:]
    rs  = rs[sortr]

    msh = PolygonMesh()
    for i in range(len(rs)):
        msh |= makeSphere(xyz[i], rs[i], 4, 4)
        print i, rs[i]
        
        if (i != 0) and (i % 25 == 0):
            msh.writeOBJ("vasculature-test-%i.obj" % i)




if __name__ == "__main__":
    # doVasculature("spheres.npy")
    # t5 = makeSphere(np.zeros(3), 1.0, 7, 7)
    # t6 = makeSphere(np.ones(3), 1, 7, 7)

    import meshIO
    m1 = meshIO.loadABA("grey.msh")
    m2 = makeSphere(np.array([2941, 3325, 2549]), 1000, 10, 10)

    over = PolygonMesh(m1.polygons + m2.polygons)
    over.writeOBJ("brain-super.obj")
    
    m3 = m1 - m2
    m3.writeOBJ("brain-minus-sphere.obj")

    t1 = makeTetrahedron(np.array([0,0,0]))
    t2 = makeTetrahedron(np.array([0,0,np.sqrt(2)+0.1]))
    testOperators(t1, t2, "tet")

    # meshes = []
    # meshes = [makeSphere(np.ones(3)*0.5*x, 1, 5, 5) for x in range(10)]

    # result = union(*tuple(meshes))

    # result.writeOBJ("sphere-merge.obj")
