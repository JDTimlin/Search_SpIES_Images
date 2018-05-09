import numpy as np

class point(object):
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def dot(self, b):
        """Compute the dot product of two vectors"""
        return self.x*b.x + self.y*b.y + self.z*b.z

    def cross(self,  b):
        """Compute the cross product of two vectors"""
        return point((self.y*b.z - b.y*self.z),
                     (self.z*b.x - b.z*self.x),
                     (self.x*b.y - b.x*self.y)
                    )

    def normed(self):
        """Return a normalized copy of the vector"""
        d = np.sqrt(self.dot(self))
        return point(self.x/d, self.y/d, self.z/d)

    def dist_to(self, b):
        diff = self - b
        return np.sqrt(diff.dot(diff))

    def __add__(self, b):
        return point(self.x + b.x,
                     self.y + b.y,
                     self.z + b.z
                    )

    def __sub__(self, b):
        return point(self.x - b.x,
                     self.y - b.y,
                     self.z - b.z
                    )

    def __radd__(self, b):
        return self + b

    def __repr__(self):
        return str((self.x, self.y, self.z))

def midpoint(a, b):
    return (a+b).normed()

def opposites(pts):
    maxdist = 0.0
    for a in pts:
        for b in pts:
            if a is b: continue

            if a.dist_to(b) > maxdist:
                pair = (a,b)
                maxdist = a.dist_to(b)

    return pair

def center(pts):
    s = sum(pts, point(0.0,0.0,0.0))
    return s.normed()

class greatcap:
    # a spherical cap covering exactly half of the sphere, defined by the unit
    # vector of the north polar axis
    def __init__(self, a, b):
        self.north = a.cross(b).normed()
        self.defining_points = (a,b)
        
    def contains(self, pt):
    # return whether the cap contains another point
        return self.north.dist_to(pt) < np.sqrt(2)

    def invert(self):
        self.north = point(-self.north.x, -self.north.y, -self.north.z)

    def __repr__(self):
        return str(self.north)

def ccw(a, b, c):
    n = (b - a).cross(c - b)
    return n.dot(a + b + c) > 0.0

# alternate

#def ccw(a, b, c):
#    n = (a.cross(b)).cross(b.cross(c))
#    return n.dot( hemisphere.center ) > 0.0

def is_poly_convex(pts):
    triplets = zip(pts, offset(pts, 1), offset(pts, 2))

    for a, b, c in triplets:
        if not ccw(a,b,c):
            return False
    return True

def offset(lst, i):
    """Returns a cyclically shifted list where the ith element becomes the
    ith element
    """
    return lst[i:] + lst[:i]

class sphpoly:
    def __init__(self, pts):
        self.pts = pts

        # TODO 

        # make bounding circle
        self.center = center(pts)
        self.radius = 2 * np.arcsin(max([self.center.dist_to(x) for x in pts]) / 2.)

        # if points are not initially CCW, reverse and check again
        if not is_poly_convex(self.pts):
            #print "reversing"
            self.pts = list(reversed(self.pts))

        if not is_poly_convex(self.pts):
            raise ValueError("points are not ordered to form a convex polygon")

        # make pairs of points going around the polygon
        pairs = zip(self.pts, offset(self.pts, 1))
 
        # each pair makes a great circle edge, which is defined by a spherical
        # cap
        self.caps = [greatcap(a,b) for (a,b) in pairs]

    def contains(self, pt):
        for cap in self.caps:
            if not cap.contains(pt):
                return False
        return True

    def overlaps(self, poly):
        """Return whether input polygon overlaps this one, i.e. whether any
        point which defines one polygon is contained by the other
        """
        for pt in self.pts:
            if poly.contains(pt):
                return True

        for pt in poly.pts:
            if self.contains(pt):
                return True

        return False

    def is_convex(self):
        triplets = zip(self.pts, offset(self.pts, 1), offset(self.pts, 2))

        for a, b, c in triplets:
            if not ccw(a,b,c):
                return False
        return True

    # Alternate method: check that all points except those which define
    # each cap lie inside the cap
    #def is_convex(self):
    #    for cap in self.caps:
    #        for pt in list(set(self.pts) - set(cap.defining_points)):
    #            if not cap.contains(pt):
    #                return False
    #    return True

def main():
    delta = 1e-1
    pts = [
            point(1.0,   0.0,   0.0).normed(),
            point(1.0,   0.0, delta).normed(),
            point(1.0, delta, delta).normed(),
            point(1.0, delta,   0.0).normed(),
            ]

    poly = sphpoly(pts)

    print "Convex:", poly.is_convex()

    testpts = [
            point( 1.0,   delta/2.,  -delta/2.).normed(),
            point( 1.0,  -delta/2.,   delta/2.).normed(),
            point( 1.0, 3*delta/2.,   delta/2.).normed(),
            point( 1.0,   delta/2., 3*delta/2.).normed(),
            point(-1.0,  -delta/2.,  -delta/2.).normed(),
            point( 1.0,   delta/2.,   delta/2.).normed(),
            ]

    for p in testpts:
        print poly.contains(p)

    print poly.radius


if __name__ == "__main__":
    main()
