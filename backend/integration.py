"""
    mode_calcs.py is a subroutine of NumBAT that contains methods to
    calculate the EM and Acoustic modes of a structure.

    Copyright (C) 2016  Bjorn Sturmberg, Kokou Dossou 
"""

import numpy as np

def quad_triangle(nquad):
    """ Implementation of quad_triangle

                   Integration numerique
                   ---------------------
    Evalue les integrales elementaires des composantes convectives
    sur chaque triangle. on utilise ici la methode de hammer a
    seize points de gauss qui integre exactement des polynomes du
    huitieme degre.

    Google Translate:
    Evaluates integrals elementary convective components
    on each triangle. here we use the method of a hammer
    sixteen points of gauss that integrates exactly the polynomials
    eighth degree.

    Reference
    J. N. Lyness and D. Jespersen
    "Moderate Degree Symmetric Quadrature Rules for the Triangle"
    J. Inst. Math. Appl., 1975, 15(1), pp. 19-32
    "J. Inst. Math. Appl." is now Continued as "IMA J. Appl. Math."
    J. Inst. Math. Appl. = Journal of the Institute of Mathematics and its Applications
    IMA J. Appl. Math.   = IMA Journal of Applied Mathematics
    """

    if nquad is not 16:
        raise ValueError, 'integration.quad_triangle nquad must == 16.'
    wq = np.zeros(nquad)
    xq = np.zeros(nquad)
    yq = np.zeros(nquad)

    poidbar = 1.443156076777862e-1 / 2.0
    poid1   = 2.852749028018549e-1 / 6.0
    poid2   = 9.737549286959440e-2 / 6.0
    poid3   = 3.096521116041552e-1 / 6.0
    poid4   = 1.633818850466092e-1 / 12.0

    coorbar = 1.0 / 3.0

    coor1grp1 = 4.592925882927229e-1
    coor2grp1 = 8.141482341455413e-2
    coor1grp2 = 5.054722831703103e-2
    coor2grp2 = 8.989055433659379e-1
    coor1grp3 = 1.705693077517601e-1
    coor2grp3 = 6.588613844964797e-1
    coor1grp4 = 7.284923929554041e-1
    coor2grp4 = 2.63112829634638689e-1
    coor3grp4 = 8.394777409957211e-3

    i = 0
    xq[i] = coorbar
    yq[i] = coorbar
    wq[i] = poidbar
    i = 1
    xq[i] = coor1grp1
    yq[i] = coor1grp1
    wq[i] = poid1
    i = 2
    xq[i] = coor1grp1
    yq[i] = coor2grp1
    wq[i] = poid1
    i = 3
    xq[i] = coor2grp1
    yq[i] = coor1grp1
    wq[i] = poid1
    i = 4
    xq[i] = coor1grp2
    yq[i] = coor1grp2
    wq[i] = poid2
    i = 5
    xq[i] = coor1grp2
    yq[i] = coor2grp2
    wq[i] = poid2
    i = 6
    xq[i] = coor2grp2
    yq[i] = coor1grp2
    wq[i] = poid2
    i = 7
    xq[i] = coor1grp3
    yq[i] = coor1grp3
    wq[i] = poid3
    i = 8
    xq[i] = coor1grp3
    yq[i] = coor2grp3
    wq[i] = poid3
    i = 9
    xq[i] = coor2grp3
    yq[i] = coor1grp3
    wq[i] = poid3
    i = 10
    xq[i] = coor1grp4
    yq[i] = coor2grp4
    wq[i] = poid4
    i = 11
    xq[i] = coor2grp4
    yq[i] = coor1grp4
    wq[i] = poid4
    i = 12
    xq[i] = coor2grp4
    yq[i] = coor3grp4
    wq[i] = poid4
    i = 13
    xq[i] = coor3grp4
    yq[i] = coor2grp4
    wq[i] = poid4
    i = 14
    xq[i] = coor3grp4
    yq[i] = coor1grp4
    wq[i] = poid4
    i = 15
    xq[i] = coor1grp4
    yq[i] = coor3grp4
    wq[i] = poid4

    return [wq, xq, yq]