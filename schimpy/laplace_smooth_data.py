import numpy as np


def laplace_smooth_data2(mesh,data,kappa=0.05,dt=1.0,iter_total=150):
    """
    Apply Laplacian smoothing to a scalar field on a mesh using explicit diffusion.

    This function performs iterative smoothing by updating each node's value
    toward the average of its neighbors. It corresponds to forward Euler
    integration of the diffusion equation.

    Parameters
    ----------
    mesh : object
        A mesh object with a `nodes` array and a method `get_neighbor_nodes(ndx)`
        returning neighbor node indices for a given node.
    data : ndarray of shape (N,)
        Initial scalar field values at each mesh node.
    kappa : float, optional
        Diffusivity coefficient. Controls the rate of smoothing. Default is 0.05.
    dt : float, optional
        Timestep used for each smoothing iteration. Default is 1.0.
    iter_total : int, optional
        Number of smoothing iterations to perform. Default is 150.

    Returns
    -------
    smooth_data : ndarray of shape (N,)
        The smoothed scalar field after `iter_total` iterations.
    """    
    nodes = mesh.nodes
    edges = mesh.edges
    rate = dt * kappa
    x = nodes[:, 0:2]
    iter = 0
    smooth_data = data.copy()
    while iter < iter_total:
        zz = smooth_data.copy()
        for ndx in range(nodes.shape[0]):
            nds = mesh.get_neighbor_nodes(ndx)
            vnode = smooth_data[ndx]
            vneighbor = [smooth_data[n] for n in nds]
            vave = np.mean(vneighbor)
            zz[ndx] = vnode + rate * (vave - vnode)

        iter += 1
        smooth_data = zz
    return smooth_data


 

def laplace_smooth_with_vel2(mesh,data,vel,kappa=0.05,dt=1.,iter_total=1):
    """
    Apply Laplacian smoothing with additive forcing to a scalar field on a mesh.

    This function performs iterative diffusion of `data` while incorporating
    a source or forcing term `vel` at each node. It evolves the field
    using explicit Euler integration of the PDE:
        ∂u/∂t = κ ∇²u + vel

    Parameters
    ----------
    mesh : object
        A mesh object with a `nodes` array and a method `get_neighbor_nodes(ndx)`
        returning neighbor node indices for a given node.
    data : ndarray of shape (N,)
        Initial scalar field values at each mesh node.
    vel : ndarray of shape (N,)
        Additive source or forcing term to apply at each node.
    kappa : float, optional
        Diffusivity coefficient. Default is 0.05.
    dt : float, optional
        Timestep used for each smoothing iteration. Default is 1.0.
    iter_total : int, optional
        Number of smoothing iterations to perform. Default is 1.

    Returns
    -------
    smooth_data : ndarray of shape (N,)
        The smoothed and forced scalar field after `iter_total` iterations.
    """    
    nodes = mesh.nodes
    # edges = mesh.edges
    # x =nodes[:,0:2]
    iter = 0
    smooth_data = data.copy()
    while iter < iter_total:
        zz = smooth_data.copy()
        # print "iter = {} smooth_data = {} ".format(iter,smooth_data[100000:100005])
        for ndx in range(nodes.shape[0]):
            nds = mesh.get_neighbor_nodes(ndx)
            vnode = smooth_data[ndx]
            vneighbor = [smooth_data[n] for n in nds]
            vave = np.mean(vneighbor)
            zz[ndx] = vnode + dt * kappa * (vave - vnode) + dt * vel[ndx]

        iter += 1
        smooth_data = zz
    return smooth_data



def laplace_smooth_with_vel3(mesh,nlayer,data,vel,kappa=0.05,dt=1.,iter_total=1):
    """
    Apply layer-aware Laplacian smoothing with forcing to a scalar field on a mesh.

    This function diffuses values across mesh nodes but restricts neighbor contributions
    to those from the same or deeper layers. Shallower neighbors are excluded from
    smoothing influence, preserving stratification or flow directionality.

    Parameters
    ----------
    mesh : object
        A mesh object with a `nodes` array and a method `get_neighbor_nodes(ndx)`
        returning neighbor node indices for a given node.
    nlayer : ndarray of shape (N,)
        Integer layer index for each node. Smoothing ignores neighbors with shallower layer values.
    data : ndarray of shape (N,)
        Initial scalar field values at each mesh node (used as fallback for excluded neighbors).
    vel : ndarray of shape (N,)
        Additive source or forcing term applied at each node.
    kappa : float, optional
        Diffusivity coefficient. Default is 0.05.
    dt : float, optional
        Timestep used for each smoothing iteration. Default is 1.0.
    iter_total : int, optional
        Number of smoothing iterations to perform. Default is 1.

    Returns
    -------
    smooth_data : ndarray of shape (N,)
        The smoothed scalar field after `iter_total` iterations with layer-aware filtering and forcing.
    """    
    nodes = mesh.nodes
    # edges = mesh.edges
    # x =nodes[:,0:2]
    iter = 0
    smooth_data = data.copy()
    while iter < iter_total:
        zz = smooth_data.copy()
        for ndx in range(nodes.shape[0]):
            nds = mesh.get_neighbor_nodes(ndx)
            nlayer_node = nlayer[ndx]
            vnode = smooth_data[ndx]
            vneighbor = smooth_data[nds]
            vbase = data[nds]
            neighbor = np.where(nlayer[nds] >= nlayer_node, vneighbor, vbase)
            vave = np.mean(neighbor)
            zz[ndx] = vnode + dt * kappa * (vave - vnode) + dt * vel[ndx]

        iter += 1
        smooth_data = zz
    return smooth_data
