import numpy as np
from scipy.linalg import block_diag

def get_l_transformation(l):
    """
    Returns the transformation matrix T (N_sph, N_cart) for angular momentum l.
    Matches the Cartesian order: _CART[l]
    Spherical order: m = -l, ..., 0, ..., +l
    """
    if l == 0: return np.array([[1.0]])
    if l == 1: return np.eye(3)

    if l == 2: # d-orbitals (5x6)
        rt3 = np.sqrt(3.0)
        T = np.zeros((5, 6))
        T[0, 3] = 1.0                                # dxy (1,1,0)
        T[1, 5] = 1.0                                # dyz (0,1,1)
        T[2, 0], T[2, 1], T[2, 2] = -0.5, -0.5, 1.0  # dz2 (2,0,0), (0,2,0), (0,0,2)
        T[3, 4] = 1.0                                # dxz (1,0,1)
        T[4, 0], T[4, 1] = 0.5*rt3, -0.5*rt3         # dx2-y2
        return T

    if l == 3: # f-orbitals (7x10)
        T = np.zeros((7, 10))
        rt5, rt3 = np.sqrt(5.0), np.sqrt(3.0)
        # m=0: fz3
        T[3, 2], T[3, 4], T[3, 6] = 1.0, -1.5, -1.5 
        # m=+1: fxz2
        T[4, 7], T[4, 0], T[4, 5] = 1.0, -0.25, -0.25
        # m=-1: fyz2
        T[2, 8], T[2, 1], T[2, 3] = 1.0, -0.25, -0.25
        # m=+2: fz(x2-y2)
        T[5, 4], T[5, 6] = 0.5*rt3, -0.5*rt3
        # m=-2: fxyz
        T[1, 9] = 1.0
        # m=+3: fx(x2-3y2)
        T[6, 0], T[6, 5] = 0.25*rt5/rt3, -0.75*rt5/rt3
        # m=-3: fy(3x2-y2)
        T[0, 3], T[0, 1] = 0.75*rt5/rt3, -0.25*rt5/rt3
        return T

    raise NotImplementedError(f"l={l} not implemented.")


def spherical_to_cartesian_matrix(C_sph, shell_l_list):
    """
    Transforms MO coefficients from Spherical back to Cartesian.
    
    Args:
        C_sph (np.ndarray): MO coefficients in Spherical basis (N_sph, N_mo)
        shell_l_list (list): List of l-values for each shell.
        
    Returns:
        np.ndarray: Matrix in Cartesian basis (N_cart, N_mo)
    """
    # 1. Build the same T matrix used for the forward transform
    blocks = [get_l_transformation(l) for l in shell_l_list]
    T = block_diag(*blocks)
    
    # 2. Use the transpose to map Spherical -> Cartesian
    # C_cart = T.T @ C_sph
    return np.dot(T.T, C_sph)
