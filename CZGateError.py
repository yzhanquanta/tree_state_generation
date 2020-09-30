import numpy as np


def overlap(gammaBOverGammaR: float) -> float:
	"""Compute the overlap between ideal scattered wavepacket and distorted wavepacket."""
	return -1 + 4*gammaBOverGammaR**2


def fidelity(photonNum: int, parentRelation: np.array, gammaBOverGammaR: float) -> float:
	"""
	Compute the fidelity of tree state due to CZ gate error.

	Args:
		photonNum: The total number of photons in the tree;
		stateVec: The state vector of the tree;
		parentRelation: A 1D array restoring the parent node information, 0 means no parent node.
	"""
	dimTree = np.power(2, photonNum)	# we only consider {|0>, |1>} here without |vac>
	normEachPhoton = 1/photonNum	# normalization for each photon's wavepacket

	count = 0	# count the number of distorted wavepackets
	for ii in range(dimTree):
		binaryStr = "{0:20b}".format(ii)
		for jj in range(photonNum - 1):
			digitTemp = binaryStr[20 - photonNum + jj]
			parentTemp = binaryStr[20 - photonNum + parentRelation[jj] - 1]
			# if (digitTemp == ' '):
			# 	digitTemp == '0'
			if (digitTemp == '1') and (parentTemp == '1'):
				count = count + 1

	overlapOfWP = overlap(gammaBOverGammaR)
	F = (1/dimTree*(dimTree - count/photonNum - count*overlapOfWP/photonNum))**2

	return F
