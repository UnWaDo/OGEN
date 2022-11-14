import numpy as np


"""
A sript to generate van der Waals surface of molecules.
"""

# Van der Waals radii (in angstrom) from 10.1021/jp8111556
vdw_r = {'H': 1.10, 'HE': 1.40,
         'LI': 1.81, 'BE': 1.53, 'B': 1.92, 'C': 1.70,
         'N': 1.55, 'O': 1.42, 'F': 1.47, 'NE': 1.54,
         'NA': 2.27, 'MG': 1.73, 'AL': 1.84, 'SI': 2.10,
         'P': 1.80, 'S': 1.80, 'CL': 1.75, 'AR': 1.88, 'BR': 1.83, 'I': 1.98}

bohr_to_angstrom = 0.52917721092

def surface(n):
    """Computes approximately n points on unit sphere. Code adapted from GAMESS.

    Parameters
    ----------
    n : int
        approximate number of requested surface points

    Returns
    -------
    ndarray
        numpy array of xyz coordinates of surface points
    """

    u = []
    eps = 1e-10
    nequat = int(np.sqrt(np.pi*n))
    nvert = int(nequat/2)
    nu = 0
    for i in range(nvert+1):
        fi = np.pi*i/nvert
        z = np.cos(fi)
        xy = np.sin(fi)
        nhor = int(nequat*xy+eps)
        if nhor < 1:
            nhor = 1
        for j in range(nhor):
            fj = 2*np.pi*j/nhor
            x = np.cos(fj)*xy
            y = np.sin(fj)*xy
            if nu >= n:
                return np.array(u)
            nu += 1
            u.append([x, y, z])
    return np.array(u)

def vdw_surface(coordinates, elements, scale_factor, density, input_radii):
    """Computes points outside the van der Waals surface of molecules.

    Parameters
    ----------
    coordinates : ndarray
        cartesian coordinates of the nuclei, in units of angstrom
    elements : list
        The symbols (e.g. C, H) for the atoms
    scale_factor : float
        The points on the molecular surface are set at a distance of
        scale_factor * vdw_radius away from each of the atoms.
    density : float
        The (approximate) number of points to generate per square angstrom
        of surface area. 1.0 is the default recommended by Kollman & Singh.
    input_radii : dict
        dictionary of user's defined VDW radii

    Returns
    -------
    radii : dict
        A dictionary of scaled VDW radii
    surface_points : ndarray
        array of the coordinates of the points on the surface

    """
    radii = {}
    surface_points = []
    # scale radii
    for i in elements:
        if i in radii.keys():
            continue
        if i in input_radii.keys():
            radii[i] = input_radii[i] * scale_factor
        elif i in vdw_r.keys():
            radii[i] = vdw_r[i] * scale_factor
        else:
            raise KeyError('%s is not a supported element; ' %i
                         + 'use the "VDW_RADII" option to add '
                         + 'its van der Waals radius.')
    # loop over atomic coordinates
    for i in range(len(coordinates)):
        # calculate approximate number of ESP grid points
        n_points = int(density * 4.0 * np.pi* np.power(radii[elements[i]], 2))
        # generate an array of n_points in a unit sphere around the atom
        dots = surface(n_points)
        # scale the unit sphere by the VDW radius and translate
        dots = coordinates[i] + radii[elements[i]] * dots
        for j in range(len(dots)):
            save = True
            for k in range(len(coordinates)):
                if i == k:
                    continue
                # exclude points within the scaled VDW radius of other atoms
                d = np.linalg.norm(dots[j] - coordinates[k])
                if d < radii[elements[k]]:
                    save = False
                    break
            if save:
                surface_points.append(dots[j])

    return np.array(surface_points), radii

def gen_points(
    coordinates,
    elements,
    scale_factors = [1.2, 1.4, 1.6, 1.8, 2.0],
    density = 1.0,
    input_radii = vdw_r
):
    points = []
    for scale_factor in scale_factors:
        shell, radii = vdw_surface(coordinates, elements, scale_factor,
                            density, input_radii)
        points.append(shell / bohr_to_angstrom)
    points = np.concatenate(points)
    return points
