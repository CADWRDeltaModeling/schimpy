# -*- coding: utf-8 -*-
"""invdisttree.py: inverse-distance-weighted interpolation using KDTree"""

import numpy as np
from scipy.spatial import cKDTree as KDTree


class Invdisttree:
    """Inverse-distance-weighted interpolation using KDTree

    Examples
    --------
    tree = interp_2d.Invdisttree(obs_xy) # initialize KDTree with observational points
    tree.weights(node_xy)  # calculate weights for each node.
    values_v = tree.interp(obs.temperature.values)   # perform spatial interpolation with the calculated weights from the previous step
    """

    def __init__(self, x, nnear=10, leafsize=10):
        """Constructor using coordinates and data

        Parameters
        ----------
        x : np.ndarray, shape (n,2)
            Coordinates of data points
            n is the number of data points

        nnear : positive int, optional
            The number of nearest neighbors to be included, the default is 10.

        leafsize: positive int, optional
            The number of points at which the algorithm switches over to brute-force.
            Default: 10

        """

        self.tree = KDTree(x, leafsize=leafsize)  # build a tree
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
        w = 1 / self.distances**p
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
        assert (
            ninterp + np.count_nonzero(np.isnan(obs)) <= self.nnear
        ), "not enough data points for nnear"
        assert len(self.x) == len(obs), "len(x) %d != len(z) %d" % (
            len(self.x),
            len(obs),
        )
        obs = obs.copy()
        obs[obs == -9999.0] = np.nan
        idm = np.argwhere(np.isnan(obs))  # check if any value is missing
        if len(idm) == 0:  # all data points are available
            # obs data points to perform interpolation
            val = obs[self.ix][:, :ninterp].copy()
            # weights to perform interpolation
            wn = self.w[:, :ninterp].copy()
        else:  # some data is missing
            wc = self.w.copy()
            wc[np.isin(self.ix, idm)] = np.nan
            wn = wc[:, :ninterp].copy()
            # at least one of the first ninterp obs data is missing in these rows
            idr = np.argwhere(np.isnan(wn).any(axis=1))[:, 0]
            if len(idr) > 0:
                idc = np.array(
                    [np.argwhere(~np.isnan(c))[:ninterp, 0].T for c in wc[idr]]
                )
                ix_new = self.ix[:, :ninterp].copy()
                ix_new[idr] = np.array([self.ix[r, c] for r, c in zip(idr, idc)])
                wn[idr] = np.array([wc[r, c] for r, c in zip(idr, idc)])
                val = obs[ix_new]
            else:
                # obs data points to perform interpolation
                val = obs[self.ix][:, :ninterp].copy()
                # weights to perform interpolation
                wn = self.w[:, :ninterp].copy()
        nw = np.divide(wn.T, wn.sum(axis=1)).T  # normalized weights
        # assert (abs(nw.sum(axis=1)-1)<1e-4).all()  # make sure that it sums up to very close to one
        val_interp = nw * val
        return val_interp.sum(axis=1)
