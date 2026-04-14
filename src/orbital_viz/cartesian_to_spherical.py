import numpy as np
from scipy.linalg import block_diag

def get_l_transformation(l):
    """
    Returns the transformation matrix (2l+1, N_cart) for angular momentum l.
    Matches the Cartesian order: _CART[l]
    Spherical order: m = -l, ..., 0, ..., +l
    """
    if l == 0: # s: (0,0,0) -> s
        return np.array([[1.0]])
    
    if l == 1: # p: (1,0,0), (0,1,0), (0,0,1) -> px, py, pz
        # Standard ordering for m=(-1, 0, 1) is usually (y, z, x) or (x, y, z)
        # Here we map directly: px=x, py=y, pz=z
        return np.eye(3)

    if l == 2: # d: (2,0,0), (0,2,0), (0,0,2), (1,1,0), (1,0,1), (0,1,1)
        # Target Spherical: dxy, dyz, dz2, dxz, dx2-y2 (m = -2, -1, 0, 1, 2)
        # Note: coefficients include Cartesian normalization factors
        rt3 = np.sqrt(3.0)
        # Rows: m=-2, -1, 0, 1, 2
        # Cols: x2, y2, z2, xy, xz, yz
        T = np.zeros((5, 6))
        T[0, 3] = 1.0                                # dxy
        T[1, 5] = 1.0                                # dyz
        T[2, 0], T[2, 1], T[2, 2] = -0.5, -0.5, 1.0  # dz2 (2z2 - x2 - y2) / 2
        T[3, 4] = 1.0                                # dxz
        T[4, 0], T[4, 1] = 0.5*rt3, -0.5*rt3         # dx2-y2 (sqrt(3)/2 * (x2 - y2))
        return T

    if l == 3: # f: 10 Cartesian components -> 7 Spherical
        # Cartesian: (300), (030), (003), (210), (201), (120), (021), (102), (012), (111)
        # Target Spherical: m = -3, -2, -1, 0, 1, 2, 3
        T = np.zeros((7, 10))
        rt5, rt3, rt6 = np.sqrt(5.0), np.sqrt(3.0), np.sqrt(6.0)
        
        # m=0: fz3 = z(5z2 - 3r2)/2 = (5*z3 - 3*x2z - 3*y2z)/2
        T[3, 2], T[3, 4], T[3, 6] = 1.0, -1.5, -1.5 
        # m=+1: fxz2 = x(5z2 - r2)*c = (4*xz2 - x3 - xy2)*c
        T[4, 7], T[4, 0], T[4, 5] = 1.0, -0.25, -0.25
        # m=-1: fyz2 = y(5z2 - r2)*c
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

    raise NotImplementedError(f"Hand-coded l={l} transformation not provided.")


def cartesian_to_spherical_matrix(C_cart, shell_l_list):
    """
    Transforms MO coefficients from Cartesian to Spherical.
    
    Args:
        C_cart: array (N_cart, N_mo)
        shell_l_list: list of l for each shell in the basis
    """
    blocks = [get_l_transformation(l) for l in shell_l_list]
    T = block_diag(*blocks)
    
    # C_sph = T @ C_cart
    return np.dot(T, C_cart)
