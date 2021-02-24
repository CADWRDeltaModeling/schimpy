# -*- coding: utf-8 -*-

""" invdisttree.py: inverse-distance-weighted interpolation using KDTree
    fast, solid, local
"""
from __future__ import division
import numpy as np
from scipy.spatial import cKDTree as KDTree
from shapely.geometry import LineString
# http://docs.scipy.org/doc/scipy/reference/spatial.html

__date__ = "2010-11-09 Nov"  # weights, doc

class Invdisttree:
    """ Inverse-distance-weighted interpolation using KDTree
    
    Examples
    --------   
    
    >>> invdisttree = Invdisttree( X, z )  -- data points, values
    >>> interpol = invdisttree( q, nnear=3, eps=0, p=1, weights=None, stat=0 )
    # interpolates z from the 3 points nearest each query point q;
    >>> interpol(q)
    
    Finds the 3 data points nearest point q, at distances d1 d2 d3
    and returns the IDW average of the values z1 z2 z3    
    .. math:: (z1/d1 + z2/d2 + z3/d3) / (1/d1 + 1/d2 + 1/d3) = .55 z1 + .27 z2 + .18 z3

    How many nearest neighbors should one take?
    1. start with 8 11 14 .. 28 in 2d 3d 4d .. 10d; see Wendel's formula
    2. make 3 runs with nnear= e.g. 6 8 10, and look at the results
    
    There is also a parameter p that weights nearer points more, farther points less.
    In 2d, the circles around query points have :math:`areas ~ distance^2`
    So p=2 is essentially inverse-area weighting:
    
    .. math::
    (z1/area1 + z2/area2 + z3/area3)/ (1/area1 + 1/area2 + 1/area3) = .74 z1 + .18 z2 + .08 z3
    
    Notes
    -----
        If the components of the X coordinates measure different things, Euclidean distance
        can be way off.  For example, if X0 is in the range 0 to 1
        but X1 0 to 1000, the X1 distances will swamp X0;
        rescale the data, i.e. make X0.std() ~= X1.std() .
    """
    def __init__( self, X, z, leafsize=10, stat=0 ):
        """Constructor using coordinates and data
        
        Parameters
        ----------
        X : np.ndarray 
            Coordinates of data points
            
        z : np.ndarray
            Data values at data points
        
        leafsize: int
            tree parameter
            
        stat : bool
            accumulate wsum, wn for average weights
        """
            
        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))
        self.tree = KDTree( X, leafsize=leafsize )  # build the tree
        self.x = X
        self.z = z
        self.stat = stat
        self.wn = 0
        self.wsum = None;

    def __call__( self, q, nnear=6, eps=0, p=1, weights=None, gridboundary=None ):
        """ Apply the interpolator to find nearest neighbors of each query point 

            Parameters
            ----------
            q : Nx2 
                Destination points. may be one point, or a batch of points.
   
            eps: float 
                approximate nearest, dist <= (1 + eps) * true nearest
    
            p: float
                power for decay with distance `
    
            weights: float
                optional multipliers for :math:`1 / distance^p`, of the same shape as q
            
            gridboundary: True 
                avoid crossing land
            
            Returns:
            Interpolated values
        
        """
        if gridboundary is not None: raise NotImplementedError("not implemented with gridboundary option")
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
        if self.wsum is None:
            self.wsum = np.zeros(nnear)


        # The reason we use KDTree here is for faster computation and also
        # we want to find the nearest N points and their indices.
        self.distances, self.ix = self.tree.query( q, k=nnear, eps=eps )
        interpol = np.zeros( (len(self.distances),) + np.shape(self.z[0]) )
        jinterpol = 0

        for i, (dist, ix) in enumerate(zip( self.distances, self.ix )):
            if nnear == 1:
                wz = self.z[ix]
            elif dist[0] < 1e-10:
                wz = self.z[ix[0]]
            else:  # weight z s by 1/dist --
                w = 1 / dist**p

                if weights is not None:
                    w *= weights[ix]  # >= 0

                # when the line joining the two points intersects with the boundary
                # give that point a very low weight.
                # The following code is too slow, so it needs more development.
#                if gridboundary is not None:
#                    if nnear ==1:
#                        line = LineString([q[i],self.x[ix]])
#                        if gridboundary.intersects(line):
#                            weightg = 1e-6
#                    else:
#                        weightg = []
#                        for ixi in ix:
#                            line = LineString([q[i],self.x[ixi]])
#                            if gridboundary.intersects(line):
#                                weightg.append(1e-6)
#                            else:
#                                weightg.append(1.0)
#                    w *= weightg

                w /= np.sum(w)
                wz = np.dot( w, self.z[ix] )
                if self.stat:
                    self.wn += 1
                    self.wsum += w
            interpol[jinterpol] = wz
            jinterpol += 1
        return interpol if qdim > 1  else interpol[0]

#
if __name__ == "__main__":
    import sys

    N = 10000
    Ndim = 2
    Nask = N  # N Nask 1e5: 24 sec 2d, 27 sec 3d on mac g4 ppc
    Nnear = 8  # 8 2d, 11 3d => 5 % chance one-sided -- Wendel, mathoverflow.com
    leafsize = 10
    eps = .1  # approximate nearest, dist <= (1 + eps) * true nearest
    p = 1  # weights ~ 1 / distance**p
    cycle = .25
    seed = 1

    exec( "\n".join( sys.argv[1:] ) ) # python this.py N= ...
    np.random.seed(seed )
    np.set_printoptions( 3, threshold=100, suppress=True )  # .3f

    print( "\nInvdisttree:  N %d  Ndim %d  Nask %d  Nnear %d  leafsize %d  eps %.2g  p %.2g" % (
        N, Ndim, Nask, Nnear, leafsize, eps, p) )

    def terrain(x):
        """ ~ rolling hills """
        return np.sin( (2*np.pi / cycle) * np.mean( x, axis=-1 ))

    known = np.random.uniform( size=(N,Ndim) ) ** .5  # 1/(p+1): density x^p
    z = terrain( known )
    ask = np.random.uniform( size=(Nask,Ndim) )

#...............................................................................
    invdisttree = Invdisttree( known, z, leafsize=leafsize, stat=1 )
    interpol = invdisttree( ask, nnear=Nnear, eps=eps, p=p )

    print( "average distances to nearest points: %s" % \
        np.mean( invdisttree.distances, axis=0 ) )
    print( "average weights: %s" % (invdisttree.wsum / invdisttree.wn) )
        # see Wikipedia Zipf's law
    err = np.abs( terrain(ask) - interpol )
    print( "average |terrain() - interpolated|: %.2g" % np.mean(err) )

    # print "interpolate a single point: %.2g" % \
    #     invdisttree( known[0], nnear=Nnear, eps=eps )