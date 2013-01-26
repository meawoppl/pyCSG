import csg, unittest
import numpy as np

originPt = np.zeros(3)

unitXPoint = np.array([1.0, 0.0, 0.0], dtype=np.float64)
unitYPoint = np.array([0.0, 1.0, 0.0], dtype=np.float64)
unitZPoint = np.array([0.0, 0.0, 1.0], dtype=np.float64)

smallX = np.array([1.0e-7, 0.0, 0.0], dtype=np.float64)


xyPlane = csg.Plane.fromPoints( originPt, unitXPoint, unitYPoint )
xzPlane = csg.Plane.fromPoints( originPt, unitXPoint, unitZPoint )
yzPlane = csg.Plane.fromPoints( originPt, unitYPoint, unitZPoint )

bigVert = [csg.Vertex((xyz*2) - 1 , xyz) for xyz in [unitXPoint, unitYPoint, unitZPoint] ]
bigPoly = csg.Polygon(bigVert)


binNames = ["Coplanar Front", "Coplanar Back", "Front", "Back"]

class TestPlaneSplitting(unittest.TestCase):
    def test_obvious_split(self):
        for plane in [xyPlane, xzPlane, yzPlane]:
            cf, cb, f, b = plane.splitPolygon(bigPoly)
            
            self.assertNotEqual(len(f), 0, msg="The front bin should not be empty")
            self.assertNotEqual(len(b), 0, msg="The back bin should not be empty" )

            self.assertEqual(len(cb), 0)

if __name__ == "__main__":
    unittest.main()
        


