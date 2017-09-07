"""
    materials.py is a subroutine of NumBAT that defines Material objects,
    these represent dispersive lossy refractive indices and possess
    methods to interpolate n from tabulated data.

    Copyright (C) 2017  Bjorn Sturmberg, Kokou Dossou
"""


data_location = '../backend/material_data/'


class Material(object):
    """ Represents a material with:

            Refractive index []
            Density [kg/m3]
            Stiffness tensor component [Pa]
            Photoelastic tensor component []
            Acoustic loss tensor component [Pa s]

    """
 

# Array that converts between 4th rank tensors in terms of x,y,z and Voigt notation
#               [[xx,xy,xz], [yx,yy,yz], [zx,zy,zz]]
to_Voigt = np.array([[0,5,4], [5,1,3], [4,3,2]]) 


def rotation_matrix_sum(i, j, k, l, tensor_orig, mat_R):
    """
    Inner loop of rotation matrix summation.
    """
    tensor_prime_comp = 0
    for q in range(3):
        for r in range(3):
            V1 = to_Voigt[q,r]
            for s in range(3):
                for t in range(3):
                    V2 = to_Voigt[s,t]
                    tensor_prime_comp += mat_R[i,q] * mat_R[j,r] * mat_R[k,s] * mat_R[l,t] * tensor_orig[V1,V2]

    return tensor_prime_comp


def rotate_tensor(tensor_orig, theta, rotation_axis):
    """
    Rotate all acoustic material tensor by theta radians around chosen
    rotation_axis.

    Args:
        tensor_orig  (array): Tensor to be rotated.

        theta  (float): Angle to rotate by in radians.

        rotation_axis  (str): Axis around which to rotate.
    """
    if rotation_axis == 'x-axis':
        mat_R = np.array([[1,0,0], [0,np.cos(theta),-np.sin(theta)], [0,np.sin(theta),np.cos(theta)]])
    if rotation_axis == 'y-axis':
        mat_R = np.array([[np.cos(theta),0,np.sin(theta)], [0,1,0], [-np.sin(theta),0,np.cos(theta)]])
    if rotation_axis == 'z-axis':
        mat_R = np.array([[np.cos(theta),-np.sin(theta),0], [np.sin(theta),np.cos(theta),0], [0,0,1]])

    tensor_prime = np.zeros((6,6))
    for i in range(3):
        for j in range(3):
            V1 = to_Voigt[i,j]
            for k in range(3):
                for l in range(3):
                    V2 = to_Voigt[k,l]
                    tensor_prime[V1,V2] = rotation_matrix_sum(i,j,k,l,tensor_orig,mat_R)

    return tensor_prime


def isotropic_stiffness(E, v):
    """
    Calculate the stiffness matrix components of isotropic 
    materials, given the two free parameters.

    Ref: www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_isotropic.cfm

    Args:
        E  (float): Youngs_modulus

        v  (float): Poisson_ratio
    """
    c_11 = E*(1-v)/((1+v)*(1-2*v))
    c_12 = E*(v)/((1+v)*(1-2*v))
    c_44 = (E*(1-2*v)/((1+v)*(1-2*v)))/2

    return c_11, c_12, c_44


Vacuum = Material("Vacuum")

Si_2016_Smith = Material("Si_2016_Smith")
Si_2015_Van_Laer = Material("Si_2015_Van_Laer")
Si_2013_Laude = Material("Si_2013_Laude")
Si_test_anisotropic = Material("Si_test_anisotropic")

SiO2_2016_Smith = Material("SiO2_2016_Smith")
SiO2_2015_Van_Laer = Material("SiO2_2015_Van_Laer")
SiO2_2013_Laude = Material("SiO2_2013_Laude")

As2S3_2017_Morrison = Material("As2S3_2017_Morrison")
As2S3_2016_Smith = Material("As2S3_2016_Smith")

GaAs_2016_Smith = Material("GaAs_2016_Smith")
