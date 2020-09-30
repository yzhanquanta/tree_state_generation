import numpy as np

class tree:
	"""Define a tree state class to compute the matrix representation."""
	def __init__(self, photonNum: int, photonDim: int):
		"""
		Create a tree state object.

		Args:
			photonNum: The total number of photons in the tree;
			photonDim: The dimension of each photonic Hilbert space, {|0>, |1>} or {|vac>, |0>, |1>};
			_stateVec: The state vector of the tree.
		"""
		self._photonNum = int(photonNum)
		self._photonDim = int(photonDim)

		if (self._photonDim == 2):
			initEachPhoton = np.array([[1], [1]])/np.sqrt(2)
		elif (self._photonDim == 3):	# include the vaccuum
			initEachPhoton = np.array([[0], [1], [1]])/np.sqrt(2)
		else:
			print("What is the Hilbert space dimension for each photon?")
		self._stateVec = initEachPhoton
		for ii in range(self._photonNum - 1):
			self._stateVec = np.kron(self._stateVec, initEachPhoton)


	def entangle(self, photonInd1: int, photonInd2: int):
		"""
		Apply a CZ gate to entangle photons "photonInd1" and "photonInd2".

		Args:
			photonInd1: The index of the first photon;
			photonInd2: The index of the second photon, self._photonNum >= photonInd2 > photonInd1 > 0.
		"""
		dimPhotonBetween = np.power(self._photonDim, photonInd2 - photonInd1 + 1)
		dimPhotonBefore = np.power(self._photonDim, photonInd1 - 1)
		dimPhotonAfter = np.power(self._photonDim, self._photonNum - photonInd2)

		idenBefore = np.identity(dimPhotonBefore, dtype=complex)
		idenAfter = np.identity(dimPhotonAfter, dtype=complex)
		CZBetween = np.identity(dimPhotonBetween, dtype=complex)

		for ii in range(np.power(self._photonDim, photonInd2 - photonInd1 - 1)):
			indexTemp = dimPhotonBetween - self._photonDim*ii - 1
			CZBetween[indexTemp][indexTemp] = -1

		CZMat = np.kron(idenBefore, np.kron(CZBetween, idenAfter))
		self._stateVec = CZMat @ self._stateVec


	def getStateVec(self) -> np.array:
		"""Return the state vector."""
		return self._stateVec


