"""
OXASL_OPTPCASL - Cost functions for optimizing ASL protocols

Copyright 2019 University of Oxford
"""
import numpy as np

class CostMeasure(object):

    def cost(self, cov):
        """
        Calculate cost

        :param cov: Sensitivity covariance matrix [..., 2, 2]
        :return: cost [...]
        """
        raise NotImplementedError()

class LOptimalCost(CostMeasure):
    """
    Optimize CBF or ATT
    """
    def __init__(self, A):
        self.A = A
        self.name = 'L-optimal'

    def cost(self, cov):
        """
        Cost is taken from a subset of the covariance
        matrix (e.g. just the CBF or ATT parts)
        """
        cost = np.abs(np.matmul(self.A, cov))
        # Force trace function to batch across leading dimensions
        return np.trace(cost, axis1=-1, axis2=-2)

class CBFCost(LOptimalCost):
    """
    Optimize CBF
    """
    def __init__(self):
        LOptimalCost.__init__(self, [[1, 0],  [0, 0]])
        self.name = 'L-optimal (CBF)'

class ATTCost(LOptimalCost):
    """
    Optimize ATT
    """
    def __init__(self):
        LOptimalCost.__init__(self, [[0, 0],  [0, 1]])
        self.name = 'L-optimal (ATT)'

class DOptimalCost(CostMeasure):
    """
    Optimize for both CBF and ATT variance
    """
    def __init__(self):
        self.name = 'D-optimal'

    def cost(self, cov):
        """
        Cost is determinant of covariance matrix
        """
        return np.abs(np.linalg.det(cov))
