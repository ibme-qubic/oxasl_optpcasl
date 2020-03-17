"""
OXASL_OPTPCASL - Cost functions for optimizing ASL protocols

Copyright 2019 University of Oxford
"""
import numpy as np

class CostMeasure(object):

    def cost(self, hessian):
        raise NotImplementedError()

class LOptimalCost(CostMeasure):
    """
    Optimize CBF or ATT
    """
    def __init__(self, A):
        self.A = A
        self.name = 'L-optimal'

    def cost(self, hessian):
        """
        :param hessian: Sensitivity Hessian matrix [..., 2, 2]
        :return: cost [...]
        """
        det = np.linalg.det(hessian)
        #print("det", det)
        cov = np.zeros(hessian.shape)
        cov[det != 0] = np.linalg.inv(hessian[det != 0])
        cov[det == 0] = np.inf
        #print("cov\n", cov)
        # Correct for inf*0 errors in A*inverse
        #cov[np.isnan(cov)] = np.inf

        # Change into (ml/100g/min)
        cov[..., 0, 0] = cov[..., 0, 0] * 6000 * 6000
        cov[..., 0, 1] = cov[..., 0, 1] * 6000
        cov[..., 1, 0] = cov[..., 1, 0] * 6000

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

    def cost(self, hessian):
        det_h = np.linalg.det(hessian)
        return 1.0/np.abs(det_h)
