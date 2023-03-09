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
    def __init__(self, x, nnear=10, leafsize=10):
        """Constructor using coordinates and data
        
        Parameters
        ----------
        x : np.ndarray, shape (n,2)
            Coordinates of data points 
            n is the number of data points
        
        leafsize: positive int, optional
            The number of points at which the algorithm switches over to brute-force. 
            Default: 10

        nnear : positive int, optional
            The number of nearest neighbors to be included, the default is 10.
            
        stat : bool
            accumulate wsum, wn for average weights
        """
            
        self.tree = KDTree( x, leafsize=leafsize )  # build a tree
        self.x = x
        self.nnear = nnear
        
    def weights(self, node_xy, eps=0, p=2, weights=None):
        """
        Find the nearest neighbors and assign weights to each neibors 

        Parameters
        ----------
        node_xy : np.ndarray, shape (n,2)
            Coordinates of the query points 
            n is the number of query points

        eps : positive float, optional
            approximate nearest, dist <= (1 + eps) * true nearest. The default is 0.
        p : positive int, optional
            power. The default is 2.
        weights : float, optional
            optional multiplier for weights, same dimension as node_xy
            The default is None.

        Returns
        -------
        None.

        """

        node_xy = np.asarray(node_xy)
        qdim = node_xy.ndim
        if qdim == 1:
            node_xy = np.array([node_xy])

        # calculate distances and indices of the nearest points
        self.distances, self.ix = self.tree.query(node_xy, k=self.nnear, eps=eps)
        w = 1/self.distances**p
        if weights is not None:
            w *= weights
        self.w = w      
    
    def interp(self, obs, ninterp=6):
        """
        Perform interpolation

        Parameters
        ----------
        obs : np.ndarray
            observational values; nan values are allowed and will be ignored.
            
        ninterp : positive int, ninterp<=nnear
            The number of nearest points to perform interpolation
            The default is 6. 

        Returns
        -------
        val_interp : np.ndarray
            interpolated values
        """
        assert ninterp+np.count_nonzero(np.isnan(obs))<=self.nnear, "not enough data points for nnear"
        assert len(self.x) == len(obs), "len(x) %d != len(z) %d" % (len(self.x), len(obs))
        obs[obs==-9999.0] = np.nan
        idm = np.argwhere(np.isnan(obs))  # check if any value is missing
        if len(idm)==0: # all data points are available
            val = obs[self.ix][:,:ninterp].copy()  # obs data points to perform interpolation
            wn = self.w[:,:ninterp].copy()                     # weights to perform interpolation
        else: # some data is missing
            wc = self.w.copy()                     
            wc[np.isin(self.ix,idm)] = np.nan
            wn = wc[:,:ninterp].copy()  
            idr = np.argwhere(np.isnan(wn).any(axis=1))[:,0] # at least one of the first ninterp obs data is missing in these rows
            if len(idr)>0:
                idc = np.array([np.argwhere(~np.isnan(c))[:ninterp,0].T for c in wc[idr]])
                ix_new = self.ix[:,:ninterp].copy()
                ix_new[idr] = np.array([self.ix[r,c] for r, c in zip(idr,idc)])
                wn[idr] = np.array([wc[r,c] for r,c in zip(idr,idc)])
                val = obs[ix_new]
            else:
                val = obs[self.ix][:,:ninterp].copy()  # obs data points to perform interpolation
                wn = self.w[:,:ninterp].copy()                     # weights to perform interpolation
        nw = np.divide(wn.T,wn.sum(axis=1)).T  # normalized weights
        #assert (abs(nw.sum(axis=1)-1)<1e-4).all()  # make sure that it sums up to very close to one        
        val_interp = nw*val
        return val_interp.sum(axis=1)        

#if __name__ == "__main__":
    ## example usage: 
    # tree = interp_2d.Invdisttree(obs_xy)
    # tree.weights(node_xy)
    # values_v = tree.interp(obs.temperature.values)   