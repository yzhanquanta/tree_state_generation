import numpy as np

class densityMatrix:
	"""Define relevant parameters for the density matrix."""
	def __init__(self, photonNum: int, errorProb: float, tcohOvertph: float, errorLabel: str):
		"""
		Create a densityMatrix object.

		Args:
			photonNum: The number of photons in the tree;
			errorProb: The probability a pulse gets wrong;
			errorLabel: The error type ('g'-green pulse, 'r'-red pulse, 'b'-blue pulse, 'c'-coherence time, 's'-rescattering error, 'n'-no error);
			tcohOvertph: t_coh/t_ph;
			_atomDim: The dimension of the emitter, {|g0>, |g1>, |g2>, |eL>};
			_photonDim: The dimension of each photon, {|vac>, |0>, |1>};
			_dim: The dimension of the density matrix;
			_rho: The density matrix of the emitter-photon state;
			_sigmax: The Pauli sigma X matrix;
			_sigmay: The Pauli sigma Y matrix;
			_sigmaz: The Pauli sigma Z matrix.
		"""
		self._photonNum = int(photonNum)
		self._errorProb = errorProb
		self._errorLabel = errorLabel
		self._atomDim = 4
		self._photonDim = 2
		self._dim = int(np.power(self._photonDim, self._photonNum)*self._atomDim)	# dimension of emitter = 4 (|g0>, |g1>, |g2>, |eL>)
		self._dimTree = int(self._dim/self._atomDim)	# dimension of tree photons
		self._rho = np.array([[1/2, 1/2, 0, 0], [1/2, 1/2, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])	# initialization of density matrix

		self._tcohOvertph = tcohOvertph

		# for throwing away the vaccum part of newly generated photon
		self._transBlock = np.array([[0, 1, 0], [0, 0, 1]])


	@property
	def extendAPhoton(self):
		"""Extend the Hilbert space by a photon."""
		photonState = np.array([[1], [0], [0]])	# start in vaccuum
		self._rho = np.kron(self._rho, photonState @ photonState.conj().T)

	def pi_01(self, currentPhotonNum: int):
		"""Apply a green pi-pulse to the current density matrix."""
		rotAtom = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
		# rotAtom = self._sigmax	# rotation matrix when we only consider erroneous atom space
		# rotAtom = np.array([[0, -1.0j, 0, 0], [-1.0j, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])	# rotation matrix in the atom space
		identityPhoton = np.identity(3*self._photonDim**(currentPhotonNum - 1), dtype=complex)

		rotTot = np.kron(rotAtom, identityPhoton)	# direct product to get rotation matrix
		if (self._errorLabel == 'g'):	# possibly imperfect green pulse
			rotX = np.kron(np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), identityPhoton)
			rotY = np.kron(np.array([[0, -1.0j, 0, 0], [1.0j, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), identityPhoton)
			rotZ = np.kron(np.array([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), identityPhoton)
			# rotX = np.kron(np.array([[0, 1], [1, 0]]), identityPhoton)	# rotation matrix when we only consider erroneous atom space
			# rotY = np.kron(np.array([[0, -1.0j], [1.0j, 0]]), identityPhoton)
			# rotZ = np.kron(np.array([[1, 0], [0, -1]]), identityPhoton)
			mixedRho = 1/4*(self._rho + rotX @ self._rho @ rotX.conj().T + rotY @ self._rho @ rotY.conj().T + rotZ @ self._rho @ rotZ.conj().T)
			self._rho = (1 - self._errorProb)*(rotTot @ self._rho @ rotTot.conj().T) + self._errorProb*mixedRho
		else:
			self._rho = rotTot @ self._rho @ rotTot.conj().T

	def hadamard_01(self, currentPhotonNum: int):
		"""Apply a pi/4 Y-rotation to reset the emitter."""
		rotAtom = np.array([[1/np.sqrt(2), 1/np.sqrt(2), 0, 0], [-1/np.sqrt(2), 1/np.sqrt(2), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
		identityPhoton = np.identity(self._photonDim**currentPhotonNum, dtype=complex)

		rotTot = np.kron(rotAtom, identityPhoton)	# direct product to get rotation matrix
		if (self._errorLabel == 'g'):	# possibly imperfect green pulse
			rotX = np.kron(np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), identityPhoton)
			rotY = np.kron(np.array([[0, -1.0j, 0, 0], [1.0j, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), identityPhoton)
			rotZ = np.kron(np.array([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), identityPhoton)
			mixedRho = 1/4*(self._rho + rotX @ self._rho @ rotX.conj().T + rotY @ self._rho @ rotY.conj().T + rotZ @ self._rho @ rotZ.conj().T)
			self._rho = (1 - self._errorProb)*(rotTot @ self._rho @ rotTot.conj().T) + self._errorProb*mixedRho
		else:
			self._rho = rotTot @ self._rho @ rotTot.conj().T

	def pi_12(self, currentPhotonNum: int, aim: str):
		"""Apply a red pi-pulse to the current density matrix."""
		rotAtom = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
		# rotAtom = np.array([[1, 0, 0, 0], [0, 0, -1.0j, 0], [0, -1.0j, 0, 0], [0, 0, 0, 1]])	# rotation matrix in the atom space
		if (aim == 'e'):
			identityPhoton = np.identity(3*self._photonDim**(currentPhotonNum - 1), dtype=complex)
		elif (aim == 'c'):
			identityPhoton = np.identity(self._photonDim**currentPhotonNum, dtype=complex)
		else:
			print("Aim input wrong!")

		rotTot = np.kron(rotAtom, identityPhoton)	# direct product to get rotation matrix
		if (self._errorLabel == 'r'):	# possibly imperfect red pulse
			rotX = np.kron(np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]), identityPhoton)
			rotY = np.kron(np.array([[1, 0, 0, 0], [0, 0, -1.0j, 0], [0, 1.0j, 0, 0], [0, 0, 0, 1]]), identityPhoton)
			rotZ = np.kron(np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]), identityPhoton)
			mixedRho = 1/4*(self._rho + rotX @ self._rho @ rotX.conj().T + rotY @ self._rho @ rotY.conj().T + rotZ @ self._rho @ rotZ.conj().T)
			self._rho = (1 - self._errorProb)*(rotTot @ self._rho @ rotTot.conj().T) + self._errorProb*mixedRho
		else:
			self._rho = rotTot @ self._rho @ rotTot.conj().T

	def pi3_12(self, currentPhotonNum: int, aim: str):
		"""Apply a red 3pi-pulse to the current density matrix."""
		rotAtom = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
		# rotAtom = np.array([[1, 0, 0, 0], [0, 0, 1.0j, 0], [0, 1.0j, 0, 0], [0, 0, 0, 1]])	# rotation matrix in the atom space
		if (aim == 'e'):
			identityPhoton = np.identity(3*self._photonDim**(currentPhotonNum - 1), dtype=complex)
		elif (aim == 'c'):
			identityPhoton = np.identity(self._photonDim**currentPhotonNum, dtype=complex)
		else:
			print("Aim input wrong!")

		rotTot = np.kron(rotAtom, identityPhoton)	# direct product to get rotation matrix
		if (self._errorLabel == 'r'):	# possibly imperfect red pulse
			rotX = np.kron(np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]), identityPhoton)
			rotY = np.kron(np.array([[1, 0, 0, 0], [0, 0, -1.0j, 0], [0, 1.0j, 0, 0], [0, 0, 0, 1]]), identityPhoton)
			rotZ = np.kron(np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]), identityPhoton)
			mixedRho = 1/4*(self._rho + rotX @ self._rho @ rotX.conj().T + rotY @ self._rho @ rotY.conj().T + rotZ @ self._rho @ rotZ.conj().T)
			self._rho = (1 - self._errorProb)*(rotTot @ self._rho @ rotTot.conj().T) + self._errorProb*mixedRho
		else:
			self._rho = rotTot @ self._rho @ rotTot.conj().T

	def pi_2e(self, currentPhotonNum: int):
		"""Apply a blue pi-pulse to the current density matrix."""
		rotAtom = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
		# rotAtom = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, -1.0j], [0, 0, -1.0j, 0]])	# rotation matrix in the atom space
		identityPhoton = np.identity(3*self._photonDim**(currentPhotonNum - 1), dtype=complex)

		rotTot = np.kron(rotAtom, identityPhoton)	# direct product to get rotation matrix
		if (self._errorLabel == 'b'):	# possibly imperfect blue pulse
			rotX = np.kron(np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]), identityPhoton)
			rotY = np.kron(np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, -1.0j], [0, 0, 1.0j, 0]]), identityPhoton)
			rotZ = np.kron(np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]), identityPhoton)
			mixedRho = 1/4*(self._rho + rotX @ self._rho @ rotX.conj().T + rotY @ self._rho @ rotY.conj().T + rotZ @ self._rho @ rotZ.conj().T)
			self._rho = (1 - self._errorProb)*(rotTot @ self._rho @ rotTot.conj().T) + self._errorProb*mixedRho
		else:
			self._rho = rotTot @ self._rho @ rotTot.conj().T
		
	def emission(self, photonIndex: int, timeBin: int):
		"""
		Spontaneously emit the 'timeBin' time-bin of the 'photonIndex'-th photon.

		Arg:
			photonIndex: The index of the emitted photon;
			timeBin: The time-bin of the emitted photon, 0-earlier, 1-later.
		"""
		swapBlock1 = np.zeros((3, 3))
		swapBlock1[timeBin + 1][0] = 1
		swapBlock2 = np.identity(3, dtype=complex)
		swapBlock2[0][0] = 0
		for ii in range(photonIndex - 1):
			sizeTemp = 3*np.power(2, ii)
			zeroBlock = np.zeros((sizeTemp, sizeTemp))
			swapBlock1 = np.block([[swapBlock1, zeroBlock], [zeroBlock, swapBlock1]])
			swapBlock2 = np.block([[swapBlock2, zeroBlock], [zeroBlock, swapBlock2]])

		zeroMat = np.zeros((3*2**(photonIndex - 1), 3*2**(photonIndex - 1)))
		idenMat = np.identity(3*2**(photonIndex - 1), dtype=complex)
		emiMat = np.block([[idenMat, zeroMat, zeroMat, zeroMat],
			[zeroMat, idenMat, zeroMat, swapBlock1],
			[zeroMat, zeroMat, idenMat, zeroMat],
			[zeroMat, zeroMat, zeroMat, swapBlock2]])

		self._rho = emiMat @ self._rho @ emiMat.conj().T

		# possible coherence time error, updated after the emission of each time-bin
		if (self._errorLabel == 'c'):
			identityPhoton = np.identity(3*self._photonDim**(photonIndex - 1), dtype=complex)
			rotZ = np.kron(np.array([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), identityPhoton)
			decohProb = (1 - np.exp(-1/self._tcohOvertph))/2
			self._rho = (1 - decohProb)*self._rho + decohProb*(rotZ @ self._rho @ rotZ.conj().T)


		# throw away the vaccuum part
		if (timeBin == 1):
			idenFormer = np.identity(self._atomDim*2**(photonIndex - 1), dtype=complex)
			transMat = np.kron(idenFormer, self._transBlock)
			self._rho = transMat @ self._rho @ transMat.conj().T

		
	def decohWhenWaiting(self, numSlots: int):
		"""Possible decoherence between two CZ gates. But not a problem for depth-2 trees."""
		identityPhoton = np.identity(self._dimTree, dtype=complex)
		if (self._errorLabel == 'c'):
			rotZ = np.kron(np.array([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), identityPhoton)
			decohProb = (1 - np.exp(-numSlots/self._tcohOvertph))/2
			self._rho = (1 - decohProb)*self._rho + decohProb*(rotZ @ self._rho @ rotZ.conj().T)

	
	def rescattering(self, currentPhotonNum: int, photonIndex: int, timeBin: int):
		"""
		The emitter rescatters with the 'timeBin' time-bin of the 'photonIndex'-th photon.

		Arg:
			photonIndex: The index of the rescattering photon;
			timeBin: The time-bin of the rescattering photon, 1-earlier, 0-later.
		"""
		# construct transformation matrix for CZ gate
		if (timeBin == 0):
			czBlock = np.array([[-1, 0], [0, 1]])
		elif (timeBin == 1):
			czBlock = np.array([[1, 0], [0, -1]])
		else:
			czBlock = np.identity(2, dtype=complex)
			print("Time-bin index error!")

		czMat = czBlock
		for ii in range(photonIndex - 1):
			sizeTemp = np.power(2, ii + 1)
			zeroBlock = np.zeros((sizeTemp, sizeTemp))
			czMat = np.block([[czMat, zeroBlock], [zeroBlock, czMat]])
		czMat = np.kron(czMat, np.identity(np.power(2, currentPhotonNum - photonIndex), dtype=complex))
		idenMat = np.identity(2**currentPhotonNum, dtype=complex)
		zeroMat = np.zeros((2**currentPhotonNum, 2**currentPhotonNum))

		czMat = np.block([[idenMat, zeroMat, zeroMat, zeroMat],
			[zeroMat, czMat, zeroMat, zeroMat],
			[zeroMat, zeroMat, idenMat, zeroMat],
			[zeroMat, zeroMat, zeroMat, idenMat]])

		# transform the density matrix
		# if (self._errorLabel == 's'):
		# 	self._rho = (1 - self._errorProb)*(czMat @ self._rho @ czMat.conj().T) + self._errorProb/2*(self._rho + (czMat @ self._rho @ czMat.conj().T))
		# else:
		# 	self._rho = czMat @ self._rho @ czMat.conj().T
		self._rho = czMat @ self._rho @ czMat.conj().T 	# we consider the rescattering error in a different way
		
	
	def getRho(self) -> np.array:
		"""Return the density matrix."""
		return self._rho
		
	
	def EGate(self, photonIndex: int):
		"""
		The E gate sequence to generate the 'photonIndex'-th photon.
		"""
		self.extendAPhoton
		self.pi_01(photonIndex)
		self.pi_12(photonIndex, 'e')
		self.pi_2e(photonIndex)
		self.emission(photonIndex, 0)
		self.pi_01(photonIndex)
		self.pi_12(photonIndex, 'e')
		self.pi_01(photonIndex)
		self.pi_2e(photonIndex)
		self.emission(photonIndex, 1)
		self.hadamard_01(photonIndex)

	
	def CZGate(self, currentPhotonNum: int, photonIndex: int):
		"""
		The CZ gate sequence to entangle the 'photonIndex'-th photon.
		"""
		self.pi_12(currentPhotonNum, 'c')
		self.rescattering(currentPhotonNum, photonIndex, 0)
		self.pi3_12(currentPhotonNum, 'c')
		self.rescattering(currentPhotonNum, photonIndex, 1)


	def partialTracedRho(self) -> np.array:
		"""Return the photonic density matrix by doing partial trace over atomic Hilbert space."""
		g0Part = self._rho[0*self._dimTree:1*self._dimTree, 0*self._dimTree:1*self._dimTree]
		g1Part = self._rho[1*self._dimTree:2*self._dimTree, 1*self._dimTree:2*self._dimTree]
		g2Part = self._rho[2*self._dimTree:3*self._dimTree, 2*self._dimTree:3*self._dimTree]
		eLPart = self._rho[3*self._dimTree:4*self._dimTree, 3*self._dimTree:4*self._dimTree]

		return g0Part + g1Part + g2Part + eLPart


