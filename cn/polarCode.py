import numpy as np 
import functions as pcfun
import math
import copy
import time
from time import perf_counter

class Node: #designed for progressive bit-flipping algorithm
	def __init__(self):
		'''self.curU = -1 
		self.preU = []'''
		self.sequenceU = []

class Path:

	def __init__(self, N=128, m=6):

		self.N = N 
		self.n = int(math.log2(N))
		self.LLRs = np.zeros(2 * self.N - 1, dtype=float)
		self.BITs = np.zeros((2, self.N - 1), dtype=int)
		self.vHat = np.zeros(self.N, dtype=int)
		self.uHat = np.zeros(self.N, dtype=int)
		#self.dHat = np.zeros(self.K)
		self.pathMetric = 0
		self.curState = [0 for i in range(m)]
		self.omiga = 0
		#self.unitCal = 0

	def updateLLRs(self, position):
		
		n = self.n 
		self.unitCal = 0
		if position == 0:
			nextlevel = n 
		else:
			lastlevel = (bin(position)[2:].zfill(n)).find('1') + 1
			start = int(np.power(2, lastlevel - 1)) - 1
			end = int(np.power(2, lastlevel) - 1) - 1
			
			for i in range(start, end + 1):
				self.unitCal += 1
				self.LLRs[i] = pcfun.lowerconv(self.BITs[0][i],
									self.LLRs[end + (i - start) + 1],
									self.LLRs[end + (i - start) + 1 + 2 ** (lastlevel - 1)])
				
			nextlevel = lastlevel - 1
			
		for lev in range(nextlevel, 0, -1):

			start = int(np.power(2, lev - 1)) - 1
			end = int(np.power(2, lev) - 1) - 1
			for indx in range(start, end + 1):
				exp1 = end + (indx - start)
				llr1 = self.LLRs[exp1 + 1]
				llr2 = self.LLRs[exp1 + 1 + 2 ** (lev - 1)]
				self.unitCal += 1
				self.LLRs[indx] = np.sign(llr1) * np.sign(llr2) * min(abs(llr1), abs(llr2))
	

	def updateBits(self, position):

		N = self.N
		latestbit = self.uHat[pcfun.bitreversed(position, self.n)]

		if position == N - 1:
			return
		elif position < N // 2:
			self.BITs[0][0] = latestbit
		else:
			lastlevel = (bin(position)[2:].zfill(self.n)).find('0') + 1
			self.BITs[1][0] = latestbit
			for lev in range(1, lastlevel - 1):
				st = int(np.power(2, lev - 1)) - 1
				ed = int(np.power(2, lev) - 1) - 1
				for i in range(st, ed + 1):
					self.BITs[1][ed + (i - st) + 1] = (self.BITs[0][i] + self.BITs[1][i]) % 2
					self.BITs[1][ed + (i - st) + 1 + 2 ** (lev - 1)] = self.BITs[1][i]

			lev = lastlevel - 1
			st = int(np.power(2, lev - 1)) - 1
			ed = int(np.power(2, lev) - 1) - 1
			for i in range(st, ed + 1):
				self.BITs[0][ed + (i - st) + 1] = (self.BITs[0][i] + self.BITs[1][i]) % 2
				self.BITs[0][ed + (i - st) + 1 + 2 ** (lev - 1)] = self.BITs[1][i]


	def updateBits1(self, position):

		N = self.N
		latestbit = self.uHat[pcfun.bitreversed(position, self.n)]

		if position == N - 1:
			return
		elif position < N // 2:
			self.BITs[0][0] = latestbit
		else:
			lastlevel = (bin(position)[2:].zfill(self.n)).find('0')
			self.BITs[1][0] = latestbit
			for lev in range(1, lastlevel):
				st = int(np.power(2, lev - 1)) - 1
				ed = int(np.power(2, lev) - 1) - 1
				for i in range(st, ed + 1):
					self.BITs[1][ed + (i - st) + 1] = (self.BITs[0][i] + self.BITs[1][i]) % 2
					self.BITs[1][ed + (i - st) + 1 + 2 ** (lev - 1)] = self.BITs[1][i]

			lev = lastlevel
			st = int(np.power(2, lev - 1)) - 1
			ed = int(np.power(2, lev) - 1) - 1
			for i in range(st, ed + 1):
				self.BITs[0][ed + (i - st) + 1] = (self.BITs[0][i] + self.BITs[1][i]) % 2
				self.BITs[0][ed + (i - st) + 1 + 2 ** (lev - 1)] = self.BITs[1][i]

class PolarCode:

	def __init__(self, N, K, construct, dSNR, stackSize=0, listSize=0, gen=[]):

		self.codewordLength = N
		self.infoLen = K
		self.n = int(math.log2(N))
		self.designSNR = dSNR
		self.indicesRev = [pcfun.bitreversed(i, self.n) for i in range(self.codewordLength)]
		self.polarMask = pcfun.degaBuild(N, K, dSNR) if construct.lower() == 'dega' \
					else pcfun.rmPolarBuild(N, K, dSNR) if construct.lower() == 'rmpolar' \
					else pcfun.WSConstruction(N, K, dSNR, gen) if construct.lower() == 'ws' \
					else pcfun.pacConstruction(N, K, dSNR, gen) if construct.lower() == 'my' \
					else None
		self.frozenMask = (self.polarMask + 1) % 2
		self.A = np.array(range(N), dtype=int)[self.polarMask == 1]
		self.Fn = np.array([[1, 0], [1, 1]])
		self.GN = pcfun.GN(N, self.Fn)
		self.criticalSet = pcfun.generateCriticalSet(self.frozenMask)
		self.LLRs = np.zeros(2 * self.codewordLength - 1, dtype=float)
		self.BITs = np.zeros((2, self.codewordLength - 1), dtype=int)
		self.gen = gen
		self.m = len(gen) - 1
		#self.curState = []
		#self.ANV
		self.cutoff_rate = np.zeros(N, dtype=float)
		self.pe = np.zeros(N, dtype=float)

		self.snrType = ''
		self.modu = ''

		self.D = stackSize
		self.L = listSize

		self.crcPoly = 0x0 # in normal form
		self.crcPolyArrayForm = []#np.array([int(bin(self.crcPoly)[2:][i]) for i in range(len(bin(self.crcPoly)[2:]))], dtype=int)
		self.crcWidth = 0
		self.crc8Table = []
		self.isCRC = False

		self.threshold = 0
		self.delta = 1

		self.bitMetric = 0;

# ---------------- Functions ----------------

	def PCRs(self):
		KI = self.infoLen - self.crcWidth
		GC = pcfun.getGC(KI, self.crcPolyArrayForm)
		#KC = len(g) - 1
		Rs = [[] for i in range(self.crcWidth)]
		#Ru = [[] for i in range(KC)]
		Ru = []

		for l in range(self.crcWidth):
			for i in range(KI):
				if GC[i, l + KI] == 1:
					Rs[l].append(i)
			Rs[l].append(KI + l)
	
		for l in range(self.crcWidth):
			Ru.append([self.A[i] for i in Rs[l]])
		return Ru

	def TransformPCRs(self, Ru):
		D = np.zeros([self.crcWidth, self.codewordLength], dtype=int)
		Q = [[] for i in range(self.crcWidth)]
		for l in range(self.crcWidth):
			for j in range(self.codewordLength):
				if j in Ru[l]:
					D[l, j] = 1
		tmp = [[0, 0] for _ in range(len(D))]
		for i, row in enumerate(D):
			lock = True
			for column in range(self.codewordLength):
				if row[column] == 1 and lock:
					tmp[i][0] = column
					lock = False
				elif row[column] == 1:
					tmp[i][1] += 1
				elif lock:
					continue
				else:
					break
		#print(D[[2, 5, 1, 8, 3, 6, 9, 4, 7, 10, 0]])
		#D = D[[2, 5, 1, 8, 3, 6, 9, 4, 7, 10, 0]]

		D = pcfun.rowEchelonForm(D)

		for l in range(self.crcWidth):
			for j in range(self.codewordLength):
				if D[l, j] == 1:
					Q[l].append(j)
		return Q

	def punctureUnit(self, x=[], M=0):
		Np = self.codewordLength - M
		self.p = np.ones(self.codewordLength, dtype=int)
		self.p[:Np] = 0
		self.p = self.p[self.indicesRev]
		print(self.A)
		for i, value in enumerate(self.p):
			if value == 0 and i in self.A:
				print(i)
		input('')
		x = x[self.p == 1]
		return x

# ---------------- Encode ----------------

	def rateProfile(self, info):

		u = np.zeros(self.codewordLength, dtype=int)
		u[self.polarMask == 1] = info

		return u 

	def mul_matrix(self, u):

	    n = self.n 
	    N = self.codewordLength
	    polarcoded = u

	    for i in range(n):
	        if i == 0:
	            polarcoded[0:N:2] = (polarcoded[0:N:2] + polarcoded[1:N:2]) % 2
	        elif i == (n - 1):
	            polarcoded[0:int(N/2)] = (polarcoded[0:int(N/2)] + polarcoded[int(N/2):N]) % 2
	        else:
	            enc_step = int(np.power(2, i))
	            for j in range(enc_step):
	                polarcoded[j:N:(2 * enc_step)] = (polarcoded[j:N:(2 * enc_step)]
	                                                    + polarcoded[j + np.power(2, i):N:(2 * enc_step)]) % 2
	    return polarcoded

	def polarEncode(self, info):

		u = self.rateProfile(info)
		x = self.mul_matrix(u)

		return x 

	'''def CRCPolarEncode3(self, info):

		precoded = pcfun.crcTableEncode(info, self.crcPoly, self.crcWidth)
		u = self.rateProfile(precoded)
		x = self.mul_matrix(u)

		return x'''

	def CRCPolarEncode2(self, info):

		precoded = pcfun.crcEncode(info, self.crcPoly)
		u = self.rateProfile(precoded)
		x = self.mul_matrix(u)

		return x

	def CRCPolarEncode(self, info): # CRC polynomial is in non-standard form

		#precoded = pcfun.crcEncode(info, self.crcPoly)
		precoded = np.dot(info, pcfun.getGC(self.infoLen - self.crcWidth, self.crcPolyArrayForm)) % 2
		u = self.rateProfile(precoded)
		x = self.mul_matrix(u)

		return x

	def RCPPEncode(self, info, M=0):

		precoded = np.dot(info, pcfun.getGC(self.infoLen - self.crcWidth, self.crcPolyArrayForm)) % 2
		u = self.rateProfile(precoded)
		x = self.mul_matrix(u)
		x1 = self.punctureUnit(x, M)
		return x1

	'''def RCPPEncode(self, info, M=0):
		self.punctureUnit(M=M)
		precoded = np.dot(info, pcfun.getGC(self.infoLen - self.crcWidth, self.crcPolyArrayForm)) % 2
		self.GR = self.GN[:, self.p == 1]
		self.GR = self.GR[self.A, :]
		#print(np.dot(precoded, self.GR) % 2)
		return np.dot(precoded, self.GR) % 2'''

	def pacEncode(self, info):

		v = self.rateProfile(info)
		u = pcfun.convEncode(v, self.gen, self.m)
		x = self.mul_matrix(u)
		return x

	def extract(self, uHat):

		decodedInfo = np.array(list(), dtype=int)
		mask = self.polarMask
		for i in range(len(mask)):
			if mask[i] == 1:
				decodedInfo = np.append(decodedInfo, uHat[i])
		return decodedInfo

# ---------------- Decode ----------------

	def decode(self, softMess, y=[], LMax=0):
		if self.decoderType.lower() == 'scfano':
			return self.SCFano(softMess)
		elif self.decoderType.lower() == 'sc':
			return self.scDecoder(softMess)
		elif self.decoderType.lower() == 'scl':
			return self.sclDecoder(softMess)
		elif self.decoderType.lower() == 'scs':
			return self.scsDecoder(softMess)
		elif self.decoderType.lower() == 'lva':
			return self.pacListViterbiDecoder(softMess)
		elif self.decoderType.lower() == 'sva':
			return self.pacStackViterbiDecoder(softMess)
		elif self.decoderType.lower() == 'paclist':
			return self.pacSCLDecoder(softMess)
		elif self.decoderType.lower() == 'pacstack':
			return self.pacStackDecoder2(softMess)
		elif self.decoderType.lower() == 'pacfano':
			return self.PACFano(softMess)
		elif self.decoderType.lower() == 'cs-aided':
			return self.pacStackDecoder1(softMess)
		elif self.decoderType.lower() == 'rowshanlistviterbi':
			return self.pac_viterbi_decoder(softMess)
		elif self.decoderType.lower() == 'polargenie':
			return self.oracleAssistPolarDecoder(softMess)
		elif self.decoderType.lower() == 'pbf':
			return self.progressiveBitFlipping(softMess)
		elif self.decoderType.lower() == 'polarsd':
			return self.sphereDecoderDynamic(y)
		elif self.decoderType.lower() == 'ca_sd':
			return self.CA_SD(softMess)
		elif self.decoderType.lower() == 'ca_hd':
			return self.CA_HD(y, softMess, LMax)
		else:
			input('unexpected decoder type')

	def IandD(self): #dynamic bound
		GHat = self.GN
		for i in range(self.codewordLength):
			if self.polarMask[i] == 0:
				GHat[i, :] = 0
		I = [pcfun.transform4Sphere(i, GHat) for i in range(self.codewordLength)]
		d = np.zeros(self.codewordLength, dtype=int)
		for i in range(self.codewordLength):
			if self.polarMask[i] == 1:
				d[i] = len(I[i])
		return (I, d)

	def moveBack(self, beta, j, T, gama):
		while True:

			# muPre = -np.inf if j == -1 else pm[ii] if j == 0 else beta[j - 1]
			muPre = -np.inf if j == -1 else 0 if j == 0 else beta[j - 1]
			if muPre >= T:
				if gama[j] == 0:
					j = j - 1
					B = 1
					return (T, j, B)
				else:
					j = j - 1
			else:
				B = 0
				T = T - self.delta
				return (T, j, B)

	'''def moveBack1(self, beta, j, T, gama):
		while True:

			# muPre = -np.inf if j == -1 else pm[ii] if j == 0 else beta[j - 1]
			muPre = -np.inf if j == -1 else 0 if j == 0 else beta[j - 1]
			if muPre >= T:
				if gama[j] == 0 and (self.A[j] in self.criticalSet):
					j = j - 1
					B = 1
					return (T, j, B)
				else:
					j = j - 1
			else:
				B = 0
				T = T - self.delta
				return (T, j, B)'''

	def pacStackPathFork(self, currPath):

		idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
		currPath.updateLLRs(idxRev)
		edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
		edgeValue1 = pcfun.conv1Bit(1, currPath.curState, self.gen)
		pathState0 = pcfun.getNextState(0, currPath.curState, self.m)
		pathState1 = pcfun.getNextState(1, currPath.curState, self.m)
		'''
		branchMetric0 = np.log(1 + np.exp(-(1 - 2 * edgeValue0) * currPath.LLRs[0]))
		branchMetric1 = np.log(1 + np.exp(-(1 - 2 * edgeValue1) * currPath.LLRs[0]))
		'''
		'''
		penalty = np.abs(currPath.LLRs[0])
		if currPath.LLRs[0] > 0:
			branchMetric0 = 0 if (edgeValue0 == 0) else penalty
			branchMetric1 = 0 if (edgeValue1 == 0) else penalty
		elif currPath.LLRs[0] < 0:
			branchMetric0 = 0 if (edgeValue0 == 1) else penalty
			branchMetric1 = 0 if (edgeValue1 == 1) else penalty
		else:
			input('warning')
		'''
		penalty = np.abs(currPath.LLRs[0])
		penalty /= np.log(2)
		if currPath.LLRs[0] > 0:
			branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.cutoff_rate[currPath.currPosition]
			branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.cutoff_rate[currPath.currPosition]
		elif currPath.LLRs[0] < 0:
			branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.cutoff_rate[currPath.currPosition]
			branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.cutoff_rate[currPath.currPosition]
		else:
			input('warning')
		copyPath = copy.deepcopy(currPath)
		currPath.pathMetric += branchMetric0
		currPath.uHat[currPath.currPosition] = edgeValue0
		currPath.vHat[currPath.currPosition] = 0
		currPath.curState = pathState0
		currPath.updateBits(idxRev)
		copyPath.pathMetric += branchMetric1
		copyPath.uHat[copyPath.currPosition] = edgeValue1
		copyPath.vHat[copyPath.currPosition] = 1 
		copyPath.curState = pathState1
		copyPath.updateBits(idxRev)

		if self.T <= (self.D - 2):

			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

		else:

			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
			del self.trellisPathStack[0]
			self.T -= 1
		
	def stackViterbiFork(self, currPath):

		if self.T <= (self.D - 2):

			idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
			currPath.updateLLRs(idxRev)
			edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
			edgeValue1 = pcfun.conv1Bit(1, currPath.curState, self.gen)
			pathState0 = pcfun.getNextState(0, currPath.curState, self.m)
			pathState1 = pcfun.getNextState(1, currPath.curState, self.m)
			'''
			branchMetric0 = np.log(1 + np.exp(-(1 - 2 * edgeValue0) * currPath.LLRs[0]))
			branchMetric1 = np.log(1 + np.exp(-(1 - 2 * edgeValue1) * currPath.LLRs[0]))
			'''
			'''
			penalty = np.abs(currPath.LLRs[0])
			if currPath.LLRs[0] > 0:
				branchMetric0 = 0 if (edgeValue0 == 0) else penalty
				branchMetric1 = 0 if (edgeValue1 == 0) else penalty
			elif currPath.LLRs[0] < 0:
				branchMetric0 = 0 if (edgeValue0 == 1) else penalty
				branchMetric1 = 0 if (edgeValue1 == 1) else penalty
			else:
				input('warning')
			'''
			'''
			likelihoodRatio = np.exp(currPath.LLRs[0])
			branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
			branchMetric1 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue1))) - self.cutoff_rate[currPath.currPosition]
			'''
			'''
			metric0 = (np.exp((1 - edgeValue0) * currPath.LLRs[0])) / (np.exp(currPath.LLRs[0]) + 1)
			metric1 = (np.exp((1 - edgeValue1) * currPath.LLRs[0])) / (np.exp(currPath.LLRs[0]) + 1)
			branchMetric0 = np.log(metric0 / (1 - self.pe[currPath.currPosition]))
			branchMetric1 = np.log(metric1 / (1 - self.pe[currPath.currPosition]))
			'''
			'''
			likelihoodRatio = np.exp(currPath.LLRs[0])
			branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
			branchMetric1 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue1))) - self.cutoff_rate[currPath.currPosition]
			'''
			
			penalty = np.abs(currPath.LLRs[0])
			penalty /= np.log(2)
			if currPath.LLRs[0] > 0:
				branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
				branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[currPath.currPosition]
			elif currPath.LLRs[0] < 0:
				branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
				branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[currPath.currPosition]
			else:
				input('warning')
			

			copyPath = copy.deepcopy(currPath)
			currPath.pathMetric += branchMetric0
			currPath.uHat[currPath.currPosition] = edgeValue0
			currPath.vHat[currPath.currPosition] = 0
			currPath.curState = pathState0
			currPath.updateBits(idxRev)
			copyPath.pathMetric += branchMetric1
			copyPath.uHat[copyPath.currPosition] = edgeValue1
			copyPath.vHat[copyPath.currPosition] = 1 
			copyPath.curState = pathState1
			copyPath.updateBits(idxRev)
			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

		else:
			#print('hi')
			currPath.currPosition -= 1
			self.trellisPathStack.append(currPath)
			self.T += 1
			numStates = int(2 ** self.m)
			tmpTrellisPathStack = []
			pathStateMap = [[] for i in range(numStates)]

			for i in range(self.D):
				pathStateMap[pcfun.bin2dec(self.trellisPathStack[i].curState)].append(self.trellisPathStack[i])

			for i in range(numStates):

				numBranches = len(pathStateMap[i])
				if numBranches > 0:

					if pathStateMap[i][-1].currPosition < (self.codewordLength - 1):

						if self.polarMask[pathStateMap[i][-1].currPosition + 1] == 1:

							pathStateMap[i][-1].currPosition += 1
							idxRev = pcfun.bitreversed(pathStateMap[i][-1].currPosition, self.n)
							pathStateMap[i][-1].updateLLRs(idxRev)
							edgeValue0 = pcfun.conv1Bit(0, pathStateMap[i][-1].curState, self.gen)
							edgeValue1 = pcfun.conv1Bit(1, pathStateMap[i][-1].curState, self.gen)
							pathState0 = pcfun.getNextState(0, pathStateMap[i][-1].curState, self.m)
							pathState1 = pcfun.getNextState(1, pathStateMap[i][-1].curState, self.m)
							'''
							branchMetric0 = np.log(1 + np.exp(-(1 - 2 * edgeValue0) * currPath.LLRs[0]))
							branchMetric1 = np.log(1 + np.exp(-(1 - 2 * edgeValue1) * currPath.LLRs[0]))
							'''
							'''
							penalty = np.abs(pathStateMap[i][-1].LLRs[0])
							if pathStateMap[i][-1].LLRs[0] > 0:
								branchMetric0 = 0 if (edgeValue0 == 0) else penalty
								branchMetric1 = 0 if (edgeValue1 == 0) else penalty
							elif pathStateMap[i][-1].LLRs[0] < 0:
								branchMetric0 = 0 if (edgeValue0 == 1) else penalty
								branchMetric1 = 0 if (edgeValue1 == 1) else penalty
							else:
								input('warning')
							'''
							'''
							likelihoodRatio = np.exp(pathStateMap[i][-1].LLRs[0])
							branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[pathStateMap[i][-1].currPosition]
							branchMetric1 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue1))) - self.cutoff_rate[pathStateMap[i][-1].currPosition]
							'''
							'''
							metric0 = (np.exp((1 - edgeValue0) * pathStateMap[i][-1].LLRs[0])) / (np.exp(pathStateMap[i][-1].LLRs[0]) + 1)
							metric1 = (np.exp((1 - edgeValue1) * pathStateMap[i][-1].LLRs[0])) / (np.exp(pathStateMap[i][-1].LLRs[0]) + 1)
							branchMetric0 = np.log(metric0 / (1 - self.pe[pathStateMap[i][-1].currPosition]))
							branchMetric1 = np.log(metric1 / (1 - self.pe[pathStateMap[i][-1].currPosition]))
							'''
							'''
							likelihoodRatio = np.exp(pathStateMap[i][-1].LLRs[0])
							branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[pathStateMap[i][-1].currPosition]
							branchMetric1 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue1))) - self.cutoff_rate[pathStateMap[i][-1].currPosition]
							'''
							
							penalty = np.abs(pathStateMap[i][-1].LLRs[0])
							penalty /= np.log(2)
							if pathStateMap[i][-1].LLRs[0] > 0:
								branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
								branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
							elif pathStateMap[i][-1].LLRs[0] < 0:
								branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
								branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
							else:
								input('warning')
							
							copyPath = copy.deepcopy(pathStateMap[i][-1])
							pathStateMap[i][-1].pathMetric += branchMetric0
							pathStateMap[i][-1].uHat[pathStateMap[i][-1].currPosition] = edgeValue0
							pathStateMap[i][-1].vHat[pathStateMap[i][-1].currPosition] = 0
							pathStateMap[i][-1].curState = pathState0
							pathStateMap[i][-1].updateBits(idxRev)
							copyPath.pathMetric += branchMetric1
							copyPath.uHat[copyPath.currPosition] = edgeValue1
							copyPath.vHat[copyPath.currPosition] = 1 
							copyPath.curState = pathState1
							copyPath.updateBits(idxRev)
							pathStateMap[i].append(copyPath)
							pathStateMap[i].sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
							del pathStateMap[i][0]

						else:

							pathStateMap[i][-1].currPosition += 1
							idxRev = pcfun.bitreversed(pathStateMap[i][-1].currPosition, self.n)
							pathStateMap[i][-1].updateLLRs(idxRev)
							edgeValue0 = pcfun.conv1Bit(0, pathStateMap[i][-1].curState, self.gen)
							pathStateMap[i][-1].curState = pcfun.getNextState(0, pathStateMap[i][-1].curState, self.m)
							pathStateMap[i][-1].uHat[pathStateMap[i][-1].currPosition] = edgeValue0
							pathStateMap[i][-1].vHat[pathStateMap[i][-1].currPosition] = self.polarMask[pathStateMap[i][-1].currPosition]
							'''
							branchMetric0 = np.log(1 + np.exp(-(1 - 2 * edgeValue0) * currPath.LLRs[0]))
							'''
							'''
							penalty = np.abs(pathStateMap[i][-1].LLRs[0])
							if pathStateMap[i][-1].LLRs[0] > 0:
								branchMetric0 = 0 if (edgeValue0 == 0) else penalty
								
							elif pathStateMap[i][-1].LLRs[0] < 0:
								branchMetric0 = 0 if (edgeValue0 == 1) else penalty
								
							else:
								input('warning')
							'''
							'''
							likelihoodRatio = np.exp(pathStateMap[i][-1].LLRs[0])
							branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0)))
							'''
							'''
							metric0 = (np.exp((1 - edgeValue0) * pathStateMap[i][-1].LLRs[0])) / (np.exp(pathStateMap[i][-1].LLRs[0]) + 1)
							branchMetric0 = np.log(metric0 / (1 - self.pe[pathStateMap[i][-1].currPosition]))
							'''
							'''
							likelihoodRatio = np.exp(pathStateMap[i][-1].LLRs[0])
							branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[pathStateMap[i][-1].currPosition]
							'''
							
							penalty = np.abs(pathStateMap[i][-1].LLRs[0])
							penalty /= np.log(2)
							if pathStateMap[i][-1].LLRs[0] > 0:
								branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
							elif pathStateMap[i][-1].LLRs[0] < 0:
								branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
							else:
								input('warning')
							
							pathStateMap[i][-1].pathMetric += branchMetric0
							pathStateMap[i][-1].updateBits(idxRev)
							pathStateMap[i].sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

					for k in range(numBranches):
						tmpTrellisPathStack.append(pathStateMap[i][k])

			self.trellisPathStack = tmpTrellisPathStack[:]
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

	def stackViterbiFork1(self, currPath):

		if self.T <= (self.D - 2):

			idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
			currPath.updateLLRs(idxRev)
			edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
			edgeValue1 = pcfun.conv1Bit(1, currPath.curState, self.gen)
			pathState0 = pcfun.getNextState(0, currPath.curState, self.m)
			pathState1 = pcfun.getNextState(1, currPath.curState, self.m)
			
			penalty = np.abs(currPath.LLRs[0])
			penalty /= np.log(2)
			if currPath.LLRs[0] > 0:
				branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
				branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[currPath.currPosition]
			elif currPath.LLRs[0] < 0:
				branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
				branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[currPath.currPosition]
			else:
				input('warning')
			'''likelihoodRatio = np.exp(currPath.LLRs[0])
			branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
			branchMetric1 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue1))) - self.bias[currPath.currPosition]'''

			'''if np.abs(currPath.LLRs[0]) >= 50:

				penalty = np.abs(currPath.LLRs[0])
				penalty /= np.log(2)
				if currPath.LLRs[0] > 0:
					branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
					branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[currPath.currPosition]
				elif currPath.LLRs[0] < 0:
					branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
					branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[currPath.currPosition]
			elif currPath.LLRs[0] == 0:
				input('warning')
			else:
				branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
				branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[currPath.currPosition]'''

			copyPath = copy.deepcopy(currPath)
			currPath.pathMetric += branchMetric0
			currPath.uHat[currPath.currPosition] = edgeValue0
			currPath.vHat[currPath.currPosition] = 0
			currPath.curState = pathState0
			currPath.updateBits(idxRev)
			copyPath.pathMetric += branchMetric1
			copyPath.uHat[copyPath.currPosition] = edgeValue1
			copyPath.vHat[copyPath.currPosition] = 1 
			copyPath.curState = pathState1
			copyPath.updateBits(idxRev)
			self.ANV += 1
			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

		else:
			#print('hi')
			currPath.currPosition -= 1
			self.trellisPathStack.append(currPath)
			self.T += 1
			numStates = int(2 ** self.m)
			tmpTrellisPathStack = []
			pathStateMap = [[] for i in range(numStates)]
			pruneNum = 0

			for i in range(self.D):
				pathStateMap[pcfun.bin2dec(self.trellisPathStack[i].curState)].append(self.trellisPathStack[i])

			for i in range(numStates):

				numBranches = len(pathStateMap[i])
				if numBranches > 0:

					if pathStateMap[i][-1].currPosition < (self.codewordLength - 1):

						if self.polarMask[pathStateMap[i][-1].currPosition + 1] == 1:

							pathStateMap[i][-1].currPosition += 1
							idxRev = pcfun.bitreversed(pathStateMap[i][-1].currPosition, self.n)
							pathStateMap[i][-1].updateLLRs(idxRev)
							edgeValue0 = pcfun.conv1Bit(0, pathStateMap[i][-1].curState, self.gen)
							edgeValue1 = pcfun.conv1Bit(1, pathStateMap[i][-1].curState, self.gen)
							pathState0 = pcfun.getNextState(0, pathStateMap[i][-1].curState, self.m)
							pathState1 = pcfun.getNextState(1, pathStateMap[i][-1].curState, self.m)
							
							penalty = np.abs(pathStateMap[i][-1].LLRs[0])
							penalty /= np.log(2)
							if pathStateMap[i][-1].LLRs[0] > 0:
								branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
								branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
							elif pathStateMap[i][-1].LLRs[0] < 0:
								branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
								branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
							else:
								input('warning')
							'''likelihoodRatio = np.exp(pathStateMap[i][-1].LLRs[0])
							branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.bias[pathStateMap[i][-1].currPosition]
							branchMetric1 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue1))) - self.bias[pathStateMap[i][-1].currPosition]'''

							'''if np.abs(pathStateMap[i][-1].LLRs[0]) >= 50:

								penalty = np.abs(pathStateMap[i][-1].LLRs[0])
								penalty /= np.log(2)
								if pathStateMap[i][-1].LLRs[0] > 0:
									branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
									branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
								elif pathStateMap[i][-1].LLRs[0] < 0:
									branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
									branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
							elif pathStateMap[i][-1].LLRs[0] == 0:
								input('warning')
							else:
								branchMetric0 = 1 - np.log2(1 + np.exp(pathStateMap[i][-1].LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[pathStateMap[i][-1].currPosition]
								branchMetric1 = 1 - np.log2(1 + np.exp(pathStateMap[i][-1].LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[pathStateMap[i][-1].currPosition]'''
							
							copyPath = copy.deepcopy(pathStateMap[i][-1])
							pathStateMap[i][-1].pathMetric += branchMetric0
							pathStateMap[i][-1].uHat[pathStateMap[i][-1].currPosition] = edgeValue0
							pathStateMap[i][-1].vHat[pathStateMap[i][-1].currPosition] = 0
							pathStateMap[i][-1].curState = pathState0
							pathStateMap[i][-1].updateBits(idxRev)
							copyPath.pathMetric += branchMetric1
							copyPath.uHat[copyPath.currPosition] = edgeValue1
							copyPath.vHat[copyPath.currPosition] = 1 
							copyPath.curState = pathState1
							copyPath.updateBits(idxRev)
							pathStateMap[i].append(copyPath)
							pruneNum += 1
							self.ANV += 1
							#pathStateMap[i].sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
							#del pathStateMap[i][0]

						else:

							pathStateMap[i][-1].currPosition += 1
							idxRev = pcfun.bitreversed(pathStateMap[i][-1].currPosition, self.n)
							pathStateMap[i][-1].updateLLRs(idxRev)
							edgeValue0 = pcfun.conv1Bit(0, pathStateMap[i][-1].curState, self.gen)
							pathStateMap[i][-1].curState = pcfun.getNextState(0, pathStateMap[i][-1].curState, self.m)
							pathStateMap[i][-1].uHat[pathStateMap[i][-1].currPosition] = edgeValue0
							pathStateMap[i][-1].vHat[pathStateMap[i][-1].currPosition] = self.polarMask[pathStateMap[i][-1].currPosition]
							
							penalty = np.abs(pathStateMap[i][-1].LLRs[0])
							penalty /= np.log(2)
							if pathStateMap[i][-1].LLRs[0] > 0:
								branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
							elif pathStateMap[i][-1].LLRs[0] < 0:
								branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
							else:
								input('warning')
							'''likelihoodRatio = np.exp(pathStateMap[i][-1].LLRs[0])
							branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.bias[pathStateMap[i][-1].currPosition]'''

							'''if np.abs(pathStateMap[i][-1].LLRs[0]) >= 50:
								penalty = np.abs(pathStateMap[i][-1].LLRs[0])
								penalty /= np.log(2)
								if pathStateMap[i][-1].LLRs[0] > 0:
									branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
								elif pathStateMap[i][-1].LLRs[0] < 0:
									branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[pathStateMap[i][-1].currPosition]
							elif pathStateMap[i][-1].LLRs[0] == 0:
								input('warning')
							else:
								branchMetric0 = 1 - np.log2(1 + np.exp(pathStateMap[i][-1].LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[pathStateMap[i][-1].currPosition]'''
							
							self.ANV += 1
							pathStateMap[i][-1].pathMetric += branchMetric0
							pathStateMap[i][-1].updateBits(idxRev)
							#pathStateMap[i].sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
					
					for paths in pathStateMap[i]:
						tmpTrellisPathStack.append(paths)

			self.trellisPathStack = tmpTrellisPathStack[:]
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
			self.trellisPathStack = self.trellisPathStack[pruneNum:]
			if len(self.trellisPathStack) != self.D:
				
				print(self.D)
				print(len(self.trellisPathStack))
				print(pruneNum)
				print(len(tmpTrellisPathStack))
				
				input('pause')

	def pacSCLPathFork(self, pos):

		tmpTrellisPathList = []

		for i in range(len(self.trellisPathList)):

			edgeValue0 = pcfun.conv1Bit(0, self.trellisPathList[i].curState, self.gen)
			edgeValue1 = pcfun.conv1Bit(1, self.trellisPathList[i].curState, self.gen)
			pathState0 = pcfun.getNextState(0, self.trellisPathList[i].curState, self.m)
			pathState1 = pcfun.getNextState(1, self.trellisPathList[i].curState, self.m)
			'''penalty = np.abs(self.trellisPathList[i].LLRs[0])
			if self.trellisPathList[i].LLRs[0] > 0:

				branchMetric0 = 0 if (edgeValue0 == 0) else penalty
				branchMetric1 = 0 if (edgeValue1 == 0) else penalty

			elif self.trellisPathList[i].LLRs[0] < 0:

				branchMetric0 = 0 if (edgeValue0 == 1) else penalty
				branchMetric1 = 0 if (edgeValue1 == 1) else penalty

			else:
				input('warning')'''
			Li = self.trellisPathList[i].LLRs[0] / np.log(2)
			branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
			branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue1)))
			copyPath = copy.deepcopy(self.trellisPathList[i])
			self.trellisPathList[i].uHat[pos] = edgeValue0
			self.trellisPathList[i].vHat[pos] = 0
			self.trellisPathList[i].curState = pathState0
			self.trellisPathList[i].pathMetric += branchMetric0
			copyPath.uHat[pos] = edgeValue1
			copyPath.vHat[pos] = 1 
			copyPath.curState = pathState1
			copyPath.pathMetric += branchMetric1

			tmpTrellisPathList.append(copyPath)

		if (len(self.trellisPathList) + len(tmpTrellisPathList)) <= self.L:
			self.trellisPathList += tmpTrellisPathList
		else:
			
			tmpTrellisPathList += self.trellisPathList
			#tmpTrellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=False)
			tmpTrellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=True)
			
			self.trellisPathList = tmpTrellisPathList[:self.L]
	
	def listViterbiFork(self, pos):

		tmpTrellisPathList = []

		for i in range(len(self.trellisPathList)):

			edgeValue0 = pcfun.conv1Bit(0, self.trellisPathList[i].curState, self.gen)
			edgeValue1 = pcfun.conv1Bit(1, self.trellisPathList[i].curState, self.gen)
			pathState0 = pcfun.getNextState(0, self.trellisPathList[i].curState, self.m)
			pathState1 = pcfun.getNextState(1, self.trellisPathList[i].curState, self.m)
			penalty = np.abs(self.trellisPathList[i].LLRs[0])
			if self.trellisPathList[i].LLRs[0] > 0:

				branchMetric0 = 0 if (edgeValue0 == 0) else penalty
				branchMetric1 = 0 if (edgeValue1 == 0) else penalty

			elif self.trellisPathList[i].LLRs[0] < 0:

				branchMetric0 = 0 if (edgeValue0 == 1) else penalty
				branchMetric1 = 0 if (edgeValue1 == 1) else penalty

			else:
				input('warning')

			copyPath = copy.deepcopy(self.trellisPathList[i])
			self.trellisPathList[i].uHat[pos] = edgeValue0
			self.trellisPathList[i].vHat[pos] = 0
			self.trellisPathList[i].curState = pathState0
			self.trellisPathList[i].pathMetric += branchMetric0
			copyPath.uHat[pos] = edgeValue1
			copyPath.vHat[pos] = 1 
			copyPath.curState = pathState1
			copyPath.pathMetric += branchMetric1

			tmpTrellisPathList.append(copyPath)

		if (len(self.trellisPathList) + len(tmpTrellisPathList)) <= self.L:
			self.trellisPathList += tmpTrellisPathList
		else:
			
			numStates = 2 ** self.m
			#localListSize = int(self.L / numStates)
			tmpTrellisPathList += self.trellisPathList
			pathStateMap = [[] for i in range(numStates)]
			tmp1 = []

			for path in tmpTrellisPathList:
				pathStateMap[pcfun.bin2dec(path.curState)].append(path)

			for i in range(numStates):

				if len(pathStateMap[i]) > 0:

					pathStateMap[i].sort(key=lambda path_List: path_List.pathMetric, reverse=False)
					pathStateMap[i] = pathStateMap[i][:int(len(pathStateMap[i]) / 2)]

					for path in pathStateMap[i]:
						tmp1.append(path)
			
			self.trellisPathList = tmp1
			if len(self.trellisPathList) != self.L:
				input('warning')
			
	def polarSCLFork(self, pos):

		tmpPathList = []
		for path in self.pathList:

			#penalty = np.abs(path.LLRs[0])
			if path.LLRs[0] == 0:
				input('warning')
				#print('warning')
			copyPath = copy.deepcopy(path)
			path.uHat[pos] = 0
			#branchMetric0 = 0 if path.LLRs[0] > 0 else penalty
			#branchMetric0 = (0 != 1 / 2 * (1 - np.sign(path.LLRs[0]))) * penalty
			Li = path.LLRs[0] / np.log(2)
			branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
			path.pathMetric += branchMetric0
			copyPath.uHat[pos] = 1 
			#branchMetric1 = 0 if copyPath.LLRs[0] < 0 else penalty
			#branchMetric1 = (1 != 1 / 2 * (1 - np.sign(path.LLRs[0]))) * penalty
			branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
			copyPath.pathMetric += branchMetric1
			#self.ANV += 1

			tmpPathList.append(copyPath)

		if (len(tmpPathList) + len(self.pathList)) <= self.L:
			self.pathList += tmpPathList
		else:

			tmpPathList += self.pathList
			#tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=False)
			tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
			self.sortNum += 1

			self.pathList = tmpPathList[:self.L]
			if len(self.pathList) != self.L:
				input('warning')

	def scsPathFork(self, currPath):

		idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
		currPath.updateLLRs(idxRev)
		penalty = np.abs(currPath.LLRs[0])
		if currPath.LLRs[0] == 0:
			input('warning')
		branchMetric0 = 0 if (currPath.LLRs[0] > 0) else penalty
		branchMetric1 = 0 if (currPath.LLRs[0] < 0) else penalty
		copyPath = copy.deepcopy(currPath)
		currPath.pathMetric += branchMetric0
		currPath.uHat[currPath.currPosition] = 0
		currPath.updateBits1(idxRev)
		copyPath.pathMetric += branchMetric1
		copyPath.uHat[copyPath.currPosition] = 1
		copyPath.updateBits1(idxRev)
		
		
		if self.T <= (self.D - 2):
		
			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=True)
		
		else:

			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=True)
			del self.trellisPathStack[0]
			self.T -= 1

	def polarSCLFork2(self, pos):

		numPaths = len(self.pathList)
		if 2 * numPaths <= self.L:
			for ii in range(numPaths):
				copyPath = copy.deepcopy(self.pathList[ii])
				Li = copyPath.LLRs[0] / np.log(2)
				branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
				self.pathList[ii].pathMetric += branchMetric0
				copyPath.uHat[pos] = 1
				branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
				copyPath.pathMetric += branchMetric1
				self.pathList.append(copyPath)
		else:
			PM = np.zeros(2 * numPaths, dtype=float)
			discardTag = np.ones(2 * numPaths, dtype=int)
			pathStates = np.zeros(numPaths, dtype=int)
			for ii in range(numPaths):
				Li = self.pathList[ii].LLRs[0] / np.log(2)
				branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
				branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
				PM[ii] = self.pathList[ii].pathMetric + branchMetric0
				PM[ii + numPaths] = self.pathList[ii].pathMetric + branchMetric1
			self.sortNum += 1
			indices = sorted(range(2 * numPaths), key=lambda x: PM[x], reverse=True)
			indices = np.array(indices)
			discardTag[indices[:self.L]] = 0
			disIdxRecord = []
			for ii in range(numPaths):
				if discardTag[ii] == 1 and discardTag[ii + numPaths] == 1:
					pathStates[ii] = 1
					disIdxRecord.append(ii)
				elif discardTag[ii] == 1:
					pathStates[ii] = 2
				elif discardTag[ii + numPaths] == 1:
					pathStates[ii] = 3
				elif discardTag[ii] == 0 and discardTag[ii + numPaths] == 0:
					pathStates[ii] = 4
				else:
					input('state warning!')
			for ii in range(len(pathStates)):
				if pathStates[ii] == 2:
					self.pathList[ii].pathMetric = PM[ii + numPaths]
					self.pathList[ii].uHat[pos] = 1
				elif pathStates[ii] == 3:
					self.pathList[ii].pathMetric = PM[ii]
				elif pathStates[ii] == 4:
					'''copyPath = copy.deepcopy(self.pathList[ii])
					self.pathList[ii].pathMetric = PM[ii]
					self.pathList[disIdxRecord[0]].pathMetric = PM[ii + numPaths]
					copyPath.pathMetric = PM[ii + numPaths]
					copyPath.uHat[pos] = 1'''
					self.pathList[disIdxRecord[0]] = copy.deepcopy(self.pathList[ii])
					self.pathList[ii].pathMetric = PM[ii]
					self.pathList[disIdxRecord[0]].pathMetric = PM[ii + numPaths]
					self.pathList[disIdxRecord[0]].uHat[pos] = 1
					del disIdxRecord[0]

	def sclDecoderN(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.sortNum = 0
		self.unitCal = 0
		self.mT = -12
		idxRecord = []
		peSC = 1
		for i in range(self.codewordLength):
			'''if i >= self.criticalSet[-1]:
				break'''
			if self.polarMask[i] == 1 and self.GA[i] < 18:
				idxRecord.append(i)
		'''for i in range(self.codewordLength):
			if self.polarMask[i] == 1 and i not in idxRecord:
				peSC *= (1 - self.pe[i])
		print(peSC, 1 - peSC, f'{1 - peSC:.4e}')
		print(len(idxRecord), idxRecord)
		print(len(self.criticalSet))
		input('21q4e')'''

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1
			#self.unitCal += (self.pathList[0].unitCal * len(self.pathList))
			if self.polarMask[i] == 1:
				if i in idxRecord:
					self.PSCLFork2(i)
				else:
					for path in self.pathList:
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if branchMetric0 > branchMetric1:
							path.uHat[i] = 0
							path.pathMetric += branchMetric0
						else:
							path.uHat[i] = 1
							path.pathMetric += branchMetric1
			else:

				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					#branchMetric0 = np.log(1 + np.exp(-(1 - 2 * path.uHat[i]) * path.LLRs[0]))
					#penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] == 0:
						print(path.LLRs)
						input('warning')
						#print('warning')
					#branchMetric0 = 0 if (path.LLRs[0] > 0) else penalty
					#branchMetric0 = (0 != 1 / 2 * (1 - np.sign(path.LLRs[0]))) * penalty
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0
					#self.ANV += 1

			for path in self.pathList:
				path.updateBits(idxRev)


		#self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=False)
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			#self.crc8Table = pcfun.buildCRC8Table(self.crcPoly)
			best = self.extract(self.pathList[0].uHat)
			candidate = pcfun.crcEncode(best, self.crcPoly)
			#candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			#if (np.sum(pcfun.crc8TableMethod(best, self.crc8Table))) == 0:
			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = pcfun.crcEncode(nextBest, self.crcPoly)
					#nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def sclDecoderOracleGCA(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.sortNum = 0
		self.unitCal = 0
		idxRecord = []
		for i in range(self.codewordLength):
			if self.polarMask[i] == 1 and self.GA[i] < 40:
				idxRecord.append(i)
		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1
			#self.unitCal += (self.pathList[0].unitCal * len(self.pathList))
			if self.polarMask[i] == 1:
				if i in idxRecord:
					self.polarSCLFork2(i)
				else:
					for path in self.pathList:
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if branchMetric0 > branchMetric1:
							path.uHat[i] = 0
							path.pathMetric += branchMetric0
						elif branchMetric1 > branchMetric0:
							path.uHat[i] = 1
							path.pathMetric += branchMetric1
						else:
							input('b0 = b1')
			else:

				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					#branchMetric0 = np.log(1 + np.exp(-(1 - 2 * path.uHat[i]) * path.LLRs[0]))
					#penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] == 0:
						print(path.LLRs)
						input('warning')
						#print('warning')
					#branchMetric0 = 0 if (path.LLRs[0] > 0) else penalty
					#branchMetric0 = (0 != 1 / 2 * (1 - np.sign(path.LLRs[0]))) * penalty
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0
					#self.ANV += 1

			for path in self.pathList:
				path.updateBits(idxRev)


		#self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=False)
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		for path in self.pathList:
			if (path.uHat == self.message).all():
				return self.extract(path.uHat)
		return self.extract(self.pathList[0].uHat)

	def sclDecoderOracle(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.sortNum = 0
		self.unitCal = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1
			#self.unitCal += (self.pathList[0].unitCal * len(self.pathList))
			if self.polarMask[i] == 1:
				self.polarSCLFork2(i)
			else:

				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					#branchMetric0 = np.log(1 + np.exp(-(1 - 2 * path.uHat[i]) * path.LLRs[0]))
					#penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] == 0:
						print(path.LLRs)
						input('warning')
						#print('warning')
					#branchMetric0 = 0 if (path.LLRs[0] > 0) else penalty
					#branchMetric0 = (0 != 1 / 2 * (1 - np.sign(path.LLRs[0]))) * penalty
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0
					#self.ANV += 1

			for path in self.pathList:
				path.updateBits(idxRev)


		#self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=False)
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		for path in self.pathList:
			if (path.uHat == self.message).all():
				return self.extract(path.uHat)
		return self.extract(self.pathList[0].uHat)

	def sclDecoder2(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.sortNum = 0
		self.unitCal = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1
			#self.unitCal += (self.pathList[0].unitCal * len(self.pathList))
			if self.polarMask[i] == 1:
				self.polarSCLFork2(i)
			else:

				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					#branchMetric0 = np.log(1 + np.exp(-(1 - 2 * path.uHat[i]) * path.LLRs[0]))
					#penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] == 0:
						print(path.LLRs)
						input('warning')
						#print('warning')
					#branchMetric0 = 0 if (path.LLRs[0] > 0) else penalty
					#branchMetric0 = (0 != 1 / 2 * (1 - np.sign(path.LLRs[0]))) * penalty
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0
					#self.ANV += 1

			for path in self.pathList:
				path.updateBits(idxRev)


		#self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=False)
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			#self.crc8Table = pcfun.buildCRC8Table(self.crcPoly)
			best = self.extract(self.pathList[0].uHat)
			candidate = pcfun.crcEncode(best, self.crcPoly)
			#candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			#if (np.sum(pcfun.crc8TableMethod(best, self.crc8Table))) == 0:
			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = pcfun.crcEncode(nextBest, self.crcPoly)
					#nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def sclDecoder3(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.sortNum = 0
		self.unitCal = 0
		idxRecord = []
		for i in range(self.codewordLength):
			if self.polarMask[i] == 1 and self.GA[i] < 51:
				idxRecord.append(i)

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1
			
			if self.polarMask[i] == 1:
				if i in idxRecord:
					self.polarSCLFork2(i)
				else:
					for path in self.pathList:
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if branchMetric0 > branchMetric1:
							path.uHat[i] = 0
							path.pathMetric += branchMetric0
						elif branchMetric1 > branchMetric0:
							path.uHat[i] = 1
							path.pathMetric += branchMetric1
						else:
							input('b0 = b1')
			else:
				#counter = 0
				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					if path.LLRs[0] == 0:
						print(path.LLRs)
						input('warning')
						#print('warning')
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)


		#self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=False)
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			#self.crc8Table = pcfun.buildCRC8Table(self.crcPoly)
			best = self.extract(self.pathList[0].uHat)
			candidate = pcfun.crcEncode(best, self.crcPoly)
			#candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			#if (np.sum(pcfun.crc8TableMethod(best, self.crc8Table))) == 0:
			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = pcfun.crcEncode(nextBest, self.crcPoly)
					#nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def scsDecoder(self, softMess):

		self.trellisPathStack = [Path(self.codewordLength, self.m)]
		self.T = 1 
		self.trellisPathStack[0].LLRs[self.codewordLength - 1:] = softMess
		self.trellisPathStack[0].currPosition = -1
		self.ANV = 0
		flag = True

		while flag:

			currPath = self.trellisPathStack.pop()
			self.T -= 1

			if self.polarMask[currPath.currPosition + 1] == 1:

				self.ANV += 1
				currPath.currPosition += 1
				self.scsPathFork(currPath)

			else:

				self.ANV += 1
				currPath.currPosition += 1
				idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
				currPath.updateLLRs(idxRev)
				currPath.uHat[currPath.currPosition] = self.polarMask[currPath.currPosition]
				
				penalty = np.abs(currPath.LLRs[0])
				#penalty /= np.log(2)
				
				currPath.pathMetric += 0 if (currPath.LLRs[0] > 0) else penalty
				if currPath.LLRs[0] == 0:
					input('warning')
				
				currPath.updateBits1(idxRev)
				self.trellisPathStack.append(currPath)
				self.T += 1
				self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=True)

			if self.trellisPathStack[-1].currPosition == (self.codewordLength - 1):

				best = self.trellisPathStack[-1].uHat
				flag = False
				

		return self.extract(best)


	def scDecoder(self, softMess, sequenceU=[]):

		trellisPath = Path(self.codewordLength)
		trellisPath.LLRs[self.codewordLength - 1:] = softMess
		self.llrs = np.zeros(self.codewordLength, dtype=float)
		self.uHat = np.zeros(self.codewordLength, dtype=int)
		self.ANV = 0
		self.mT = -6
		self.predictIdx = []

		for i in range(self.codewordLength):

			ii = pcfun.bitreversed(i, self.n)
			trellisPath.updateLLRs(ii)

			self.ANV += 1
			self.llrs[i] = trellisPath.LLRs[0]

			if self.polarMask[i] == 1:

				if trellisPath.LLRs[0] > 0:
					trellisPath.uHat[i] = 0
				elif trellisPath.LLRs[0] < 0:
					trellisPath.uHat[i] = 1
				else:
					input('Warning')

			else:
				trellisPath.uHat[i] = 0

			if i in sequenceU:
				'''if len(sequenceU) == self.maxLevel:
					print('hi')'''
				trellisPath.uHat[i] = 1 - trellisPath.uHat[i]

			'''Li = trellisPath.LLRs[0] / (np.log(2))
			branchMetric0 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 0)))
			branchMetric1 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 1)))
			if (i not in sequenceU) and (branchMetric0 > self.mT and branchMetric1 > self.mT):
				if i in self.criticalSet[:-1]:
					self.predictIdx.append(i)
				#if self.polarMask[i] == 1:
				#	self.predictIdx2.append(i)'''
			trellisPath.updateBits(ii)

		self.uHat = trellisPath.uHat
		return self.extract(trellisPath.uHat)

	def progressiveBitFlipping(self, softMess): #caution

		maxLevel = self.maxLevel
		l = 0
		S = [[] for i in range(maxLevel)]
		nodes = [[] for i in range(maxLevel)]
		mu = pcfun.GA(self.codewordLength, self.infoLen, self.designSNR)
		self.scDecoder(softMess)
		best = self.extract(self.uHat)

		if np.sum(pcfun.crcEncode(best, self.crcPoly)[-self.crcWidth:]) != 0: 
			
			S[l] = pcfun.generateCriticalSet(self.frozenMask)
			M = np.abs(self.llrs[S[l]] / np.sqrt(mu[S[l]]))
			#M = np.abs(self.llrs[S[l]] * np.sqrt(mu[S[l]]))
			#M = np.abs(self.llrs[S[l]] / self.pe[S[l]])
			sortedIndex = M.argsort()
			S[l] = S[l][sortedIndex]
			for u in S[l]:
				node = Node()
				node.sequenceU.append(u)
				nodes[l].append(node)

			while (l <= maxLevel - 1) and np.sum(pcfun.crcEncode(best, self.crcPoly)[-self.crcWidth:]) != 0:
				
				curNode = 0
				while curNode < len(S[l]):

					self.scDecoder(softMess, nodes[l][curNode].sequenceU)
					best = self.extract(self.uHat)
					if np.sum(pcfun.crcEncode(best, self.crcPoly)[-self.crcWidth:]) != 0 and (l + 1) <= (maxLevel - 1):

						SPrime = pcfun.modifyCriticalSet(self.polarMask, S[l][curNode])
						M = np.abs(self.llrs[SPrime] / np.sqrt(mu[SPrime]))
						#M = np.abs(self.llrs[SPrime] * np.sqrt(mu[SPrime]))
						#M = np.abs(self.llrs[SPrime] / self.pe[SPrime])
						sortedIndex = M.argsort()
						SPrime = SPrime[sortedIndex]
						for u in SPrime:
							node = Node()
							node.sequenceU += nodes[l][curNode].sequenceU
							node.sequenceU.append(u)
							nodes[l + 1].append(node)
							S[l + 1].append(u)
						curNode += 1
					elif np.sum(pcfun.crcEncode(best, self.crcPoly)[-self.crcWidth:]) != 0 and (l + 1) == maxLevel:
						curNode += 1
						continue
					elif np.sum(pcfun.crcEncode(best, self.crcPoly)[-self.crcWidth:]) == 0:
						break
					else:
						input('warning')
					#curNode += 1 
				l += 1
		
		return best[:len(best) - self.crcWidth]

	def scFlip(self, softMess, sequenceU=[]):

		trellisPath = Path(self.codewordLength)
		trellisPath.LLRs[self.codewordLength - 1:] = softMess
		self.llrs = np.zeros(self.codewordLength, dtype=float)
		self.uHat = np.zeros(self.codewordLength, dtype=int)
		self.ANV = 0
		#self.mT = -10
		self.predictIdx = []

		for i in range(self.codewordLength):

			ii = pcfun.bitreversed(i, self.n)
			trellisPath.updateLLRs(ii)

			self.ANV += 1
			self.llrs[i] = trellisPath.LLRs[0]

			if self.polarMask[i] == 1:

				if trellisPath.LLRs[0] > 0:
					trellisPath.uHat[i] = 0
				elif trellisPath.LLRs[0] < 0:
					trellisPath.uHat[i] = 1
				else:
					input('Warning')

			else:
				trellisPath.uHat[i] = 0

			if i in sequenceU:
				'''if len(sequenceU) == self.maxLevel:
					print('hi')'''
				trellisPath.uHat[i] = 1 - trellisPath.uHat[i]

			Li = trellisPath.LLRs[0] / (np.log(2))
			branchMetric0 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 0)))
			branchMetric1 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 1)))
			'''if sequenceU and (i > max(sequenceU)) and (branchMetric0 > self.mT and branchMetric1 > self.mT):
				if i in self.criticalSet[:-1]:
					self.predictIdx.append(i)
				#if self.polarMask[i] == 1:
				#	self.predictIdx2.append(i)'''
			if not sequenceU:
				if i in self.criticalSet[:-1] and (branchMetric0 > self.mT[self.currL] and branchMetric1 > self.mT[self.currL]):
					self.predictIdx.append(i)
			elif (i > max(sequenceU)) and (branchMetric0 > self.mT[self.currL] and branchMetric1 > self.mT[self.currL]):
				
				if i in self.criticalSet[:-1]:
					self.predictIdx.append(i)
			trellisPath.updateBits(ii)

		self.uHat = trellisPath.uHat
		return self.extract(trellisPath.uHat)

	def scFlip2(self, softMess, sequenceU=[]):

		trellisPath = Path(self.codewordLength)
		trellisPath.LLRs[self.codewordLength - 1:] = softMess
		self.llrs = np.zeros(self.codewordLength, dtype=float)
		self.uHat = np.zeros(self.codewordLength, dtype=int)
		self.ANV = 0
		#self.mT = -10
		self.predictIdx = []

		for i in range(self.codewordLength):

			ii = pcfun.bitreversed(i, self.n)
			trellisPath.updateLLRs(ii)

			self.ANV += 1
			self.llrs[i] = trellisPath.LLRs[0]

			if self.polarMask[i] == 1:

				if trellisPath.LLRs[0] > 0:
					trellisPath.uHat[i] = 0
				elif trellisPath.LLRs[0] < 0:
					trellisPath.uHat[i] = 1
				else:
					input('Warning')

			else:
				trellisPath.uHat[i] = 0

			if i in sequenceU:
				trellisPath.uHat[i] = 1 - trellisPath.uHat[i]

			#Li = trellisPath.LLRs[0] / (np.log(2))
			#branchMetric0 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 0)))
			#branchMetric1 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 1)))

			if self.polarMask[i] == 1 and self.currL < self.maxLevel:
				Li = trellisPath.LLRs[0] / (np.log(2))
				branchMetric0 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 0)))
				branchMetric1 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 1)))
				if not sequenceU:
					if i in self.criticalSet[:-1] and min(branchMetric0, branchMetric1) > self.mT[self.currL]:
						self.predictIdx.append(i)
				else:
					if (i > max(sequenceU)) and min(branchMetric0, branchMetric1) > self.mT[self.currL]:
						if i in self.criticalSet[:-1]:
							self.predictIdx.append(i)

			'''if not sequenceU:
				if i in self.criticalSet[:-1] and (branchMetric0 > self.mT[self.currL] and branchMetric1 > self.mT[self.currL]):
					self.predictIdx.append(i)
			elif (i > max(sequenceU)) and (branchMetric0 > self.mT[self.currL] and branchMetric1 > self.mT[self.currL]):
				
				if i in self.criticalSet[:-1]:
					self.predictIdx.append(i)'''
			trellisPath.updateBits(ii)

		self.uHat = trellisPath.uHat
		return self.extract(trellisPath.uHat)

	def scFlip3(self, softMess, sequenceU=[]):

		if sequenceU:
			trellisPath = Path(self.codewordLength)
			#trellisPath.LLRs[self.codewordLength - 1:] = softMess
			trellisPath.LLRs = self.intermediateLLRs[self.currL - 1][self.curNode]#copy.deepcopy()
			trellisPath.BITs = self.intermediateBITs[self.currL - 1][self.curNode]#copy.deepcopy()
			trellisPath.uHat = self.intermediateU[self.currL - 1][self.curNode]
			self.llrs = np.zeros(self.codewordLength, dtype=float)
			self.uHat = np.zeros(self.codewordLength, dtype=int)
			counter = 0

			self.ANV = 0
			self.predictIdx = []

			for i in range(max(sequenceU) + 1, self.codewordLength):

				ii = pcfun.bitreversed(i, self.n)
				trellisPath.updateLLRs(ii)

				self.ANV += 1
				self.llrs[i] = trellisPath.LLRs[0]

				if self.polarMask[i] == 1:

					if trellisPath.LLRs[0] > 0:
						trellisPath.uHat[i] = 0
					elif trellisPath.LLRs[0] < 0:
						trellisPath.uHat[i] = 1
					else:
						input('Warning')

				else:
					trellisPath.uHat[i] = 0

				if self.polarMask[i] == 1 and self.currL < self.maxLevel:
					if counter < self.attempts[self.currL]:
						Li = trellisPath.LLRs[0] / (np.log(2))
						branchMetric0 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 1)))
						if i in self.criticalSet[:-1] and min(branchMetric0, branchMetric1) > self.mT[self.currL]:
							self.intermediateLLRs[self.currL].append(copy.deepcopy(trellisPath.LLRs))
							trellisPath.uHat[i] = 1 - trellisPath.uHat[i]
							self.intermediateU[self.currL].append(copy.deepcopy(trellisPath.uHat))
							trellisPath.updateBits(ii)
							trellisPath.uHat[i] = 1 - trellisPath.uHat[i] # recover
							self.intermediateBITs[self.currL].append(copy.deepcopy(trellisPath.BITs))
							self.predictIdx.append(i)
							counter += 1

				trellisPath.updateBits(ii)

			self.uHat = trellisPath.uHat
			return self.extract(trellisPath.uHat)

		else:

			trellisPath = Path(self.codewordLength)
			trellisPath.LLRs[self.codewordLength - 1:] = softMess
			self.llrs = np.zeros(self.codewordLength, dtype=float)
			self.uHat = np.zeros(self.codewordLength, dtype=int)
			self.ANV = 0
			self.predictIdx = []
			counter = 0

			for i in range(self.codewordLength):

				ii = pcfun.bitreversed(i, self.n)
				trellisPath.updateLLRs(ii)

				self.ANV += 1
				self.llrs[i] = trellisPath.LLRs[0]

				if self.polarMask[i] == 1:

					if trellisPath.LLRs[0] > 0:
						trellisPath.uHat[i] = 0
					elif trellisPath.LLRs[0] < 0:
						trellisPath.uHat[i] = 1
					else:
						input('Warning')

				else:
					trellisPath.uHat[i] = 0

				if self.polarMask[i] == 1 and self.currL < self.maxLevel:
					if counter < self.attempts[self.currL]:
						Li = trellisPath.LLRs[0] / (np.log(2))
						branchMetric0 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 1)))
						if i in self.criticalSet[:-1] and min(branchMetric0, branchMetric1) > self.mT[self.currL]:
							self.intermediateLLRs[self.currL].append(copy.deepcopy(trellisPath.LLRs))
							trellisPath.uHat[i] = 1 - trellisPath.uHat[i]
							self.intermediateU[self.currL].append(copy.deepcopy(trellisPath.uHat))
							trellisPath.updateBits(ii)
							trellisPath.uHat[i] = 1 - trellisPath.uHat[i] # recover
							self.intermediateBITs[self.currL].append(copy.deepcopy(trellisPath.BITs))
							self.predictIdx.append(i)
							counter += 1

				trellisPath.updateBits(ii)

			self.uHat = trellisPath.uHat
			return self.extract(trellisPath.uHat)

	def lowComplexityPBF3(self, softMess): #caution

		maxLevel = self.maxLevel
		l = 0
		self.currL = l
		S = [[] for _ in range(maxLevel + 1)]
		nodes = [[] for _ in range(maxLevel + 1)]
		self.intermediateLLRs = [[] for _ in range(maxLevel)]
		self.intermediateBITs = [[] for _ in range(maxLevel)]
		self.intermediateU = [[] for _ in range(maxLevel)]
		self.mT = [-12, -10, -6]
		attempts = [40, 35, 8]
		self.attempts = attempts
		self.scNum = 0
		while l <= maxLevel:
			if l == 0:
				self.scFlip3(softMess)
				self.scNum += 1
				best = self.extract(self.uHat)
				check = pcfun.crcEncode(best, self.crcPoly)
				if np.sum(check[-self.crcWidth:]) == 0:
					print(f'succeed in level {l}')
					return best[:len(best) - self.crcWidth]
				else:
					if l < maxLevel:
						for u in self.predictIdx[:attempts[l]]:
							node = Node()
							node.sequenceU.append(u)
							nodes[l + 1].append(node)
							S[l + 1].append(u)
			else:
				curNode = 0
				self.curNode = curNode
				while curNode < len(S[l]):
					self.scFlip3(softMess, nodes[l][curNode].sequenceU)
					self.scNum += 1
					best = self.extract(self.uHat)
					check = pcfun.crcEncode(best, self.crcPoly)
					if np.sum(check[-self.crcWidth:]) == 0:
						print(f'succeed in level {l}')
						return best[:len(best) - self.crcWidth]
					else:
						if l < maxLevel:
							for u in self.predictIdx[:attempts[l]]:
								node = Node()
								node.sequenceU += nodes[l][curNode].sequenceU
								node.sequenceU.append(u)
								nodes[l + 1].append(node)
								S[l + 1].append(u)
					curNode += 1
					self.curNode = curNode
			l += 1
			self.currL = l
				
		return best[:len(best) - self.crcWidth]

	def lowComplexityPBF2(self, softMess): #caution

		maxLevel = self.maxLevel
		l = 0
		self.currL = l
		S = [[] for i in range(maxLevel + 1)]
		nodes = [[] for i in range(maxLevel + 1)]
		self.mT = [-12, -10, -6]
		attempts = [40, 35, 8]
		self.scNum = 0
		while l <= maxLevel:
			if l == 0:
				self.scFlip2(softMess)
				self.scNum += 1
				best = self.extract(self.uHat)
				check = pcfun.crcEncode(best, self.crcPoly)
				if np.sum(check[-self.crcWidth:]) == 0:
					print(f'succeed in level {l}')
					return best[:len(best) - self.crcWidth]
				else:
					if l < maxLevel:
						for u in self.predictIdx[:attempts[l]]:
							node = Node()
							node.sequenceU.append(u)
							nodes[l + 1].append(node)
							S[l + 1].append(u)
			else:
				curNode = 0
				self.curNode = curNode
				while curNode < len(S[l]):
					self.scFlip2(softMess, nodes[l][curNode].sequenceU)
					self.scNum += 1
					best = self.extract(self.uHat)
					check = pcfun.crcEncode(best, self.crcPoly)
					if np.sum(check[-self.crcWidth:]) == 0:
						print(f'succeed in level {l}')
						return best[:len(best) - self.crcWidth]
					else:
						if l < maxLevel:
							for u in self.predictIdx[:attempts[l]]:
								node = Node()
								node.sequenceU += nodes[l][curNode].sequenceU
								node.sequenceU.append(u)
								nodes[l + 1].append(node)
								S[l + 1].append(u)
					curNode += 1
					self.curNode = curNode
			l += 1
			self.currL = l
				
		return best[:len(best) - self.crcWidth]

	def lowComplexityPBF(self, softMess): #caution

		maxLevel = self.maxLevel
		l = 0
		self.currL = l
		S = [[] for i in range(maxLevel)]
		nodes = [[] for i in range(maxLevel)]
		mu = pcfun.GA(self.codewordLength, self.infoLen, self.designSNR)
		self.mT = [-10, -7]
		attempts = [30, 16]
		self.scFlip(softMess)
		best = self.extract(self.uHat)

		if np.sum(pcfun.crcEncode(best, self.crcPoly)[-self.crcWidth:]) != 0: 

			'''S[l] = pcfun.generateCriticalSet(self.frozenMask)
			M = np.abs(self.llrs[S[l]] / np.sqrt(mu[S[l]]))
			sortedIndex = M.argsort()
			S[l] = S[l][sortedIndex]
			for u in S[l]:
				node = Node()
				node.sequenceU.append(u)
				nodes[l].append(node)'''
			for u in self.predictIdx[:attempts[l]]:
			#for u in self.predictIdx:
				node = Node()
				node.sequenceU.append(u)
				nodes[l].append(node)
				S[l].append(u)

			while (l <= maxLevel - 1) and np.sum(pcfun.crcEncode(best, self.crcPoly)[-self.crcWidth:]) != 0:
				
				curNode = 0
				#self.mT = mT[l]
				self.currL = l
				while curNode < len(S[l]):

					self.scFlip(softMess, nodes[l][curNode].sequenceU)
					best = self.extract(self.uHat)
					if np.sum(pcfun.crcEncode(best, self.crcPoly)[-self.crcWidth:]) != 0 and (l + 1) <= (maxLevel - 1):

						#SPrime = pcfun.modifyCriticalSet(self.polarMask, S[l][curNode])
						#M = np.abs(self.llrs[SPrime] / np.sqrt(mu[SPrime]))
						#M = np.abs(self.llrs[SPrime] * np.sqrt(mu[SPrime]))
						#M = np.abs(self.llrs[SPrime] / self.pe[SPrime])
						#sortedIndex = M.argsort()
						#SPrime = SPrime[sortedIndex]
						'''for u in SPrime:
							node = Node()
							node.sequenceU += nodes[l][curNode].sequenceU
							node.sequenceU.append(u)
							nodes[l + 1].append(node)
							S[l + 1].append(u)'''
						for u in self.predictIdx[:attempts[l + 1]]:
							node = Node()
							node.sequenceU += nodes[l][curNode].sequenceU
							node.sequenceU.append(u)
							nodes[l + 1].append(node)
							S[l + 1].append(u)
						curNode += 1
					elif np.sum(pcfun.crcEncode(best, self.crcPoly)[-self.crcWidth:]) != 0 and (l + 1) == maxLevel:
						curNode += 1
						continue
					elif np.sum(pcfun.crcEncode(best, self.crcPoly)[-self.crcWidth:]) == 0:
						print(f'succeed in level {l + 1}')
						break
					else:
						input('warning')
					#curNode += 1 
				l += 1
		
		return best[:len(best) - self.crcWidth]

	def sclDecoder(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.sortNum = 0
		self.unitCal = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1
			#self.unitCal += (self.pathList[0].unitCal * len(self.pathList))
			if self.polarMask[i] == 1:
				self.polarSCLFork(i)
			else:

				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					#branchMetric0 = np.log(1 + np.exp(-(1 - 2 * path.uHat[i]) * path.LLRs[0]))
					#penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] == 0:
						print(path.LLRs)
						input('warning')
						#print('warning')
					#branchMetric0 = 0 if (path.LLRs[0] > 0) else penalty
					#branchMetric0 = (0 != 1 / 2 * (1 - np.sign(path.LLRs[0]))) * penalty
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0
					#self.ANV += 1

			for path in self.pathList:
				path.updateBits(idxRev)


		#self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=False)
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			#self.crc8Table = pcfun.buildCRC8Table(self.crcPoly)
			best = self.extract(self.pathList[0].uHat)
			candidate = pcfun.crcEncode(best, self.crcPoly)
			#candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			#if (np.sum(pcfun.crc8TableMethod(best, self.crc8Table))) == 0:
			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = pcfun.crcEncode(nextBest, self.crcPoly)
					#nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def PSCLFork2(self, pos):

		numPaths = len(self.pathList)
		numAfter = numPaths
		PM = np.zeros(2 * numPaths, dtype=float)
		pathStates = np.zeros(numPaths, dtype=int)
		for ii in range(numPaths):
			Li = self.pathList[ii].LLRs[0] / np.log(2)
			branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
			branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
			if branchMetric0 >= self.mT and branchMetric1 >= self.mT:
				pathStates[ii] = 1
				PM[ii] = self.pathList[ii].pathMetric + branchMetric0
				PM[ii + numPaths] = self.pathList[ii].pathMetric + branchMetric1
				numAfter += 1
			elif branchMetric0 < self.mT:
				pathStates[ii] = 2
				PM[ii] = -np.inf
				PM[ii + numPaths] = self.pathList[ii].pathMetric + branchMetric1
			elif branchMetric1 < self.mT:
				pathStates[ii] = 3
				PM[ii] = self.pathList[ii].pathMetric + branchMetric0
				PM[ii + numPaths] = -np.inf

		if numAfter <= self.L:
			for ii in range(len(pathStates)):
				if pathStates[ii] == 1:
					copyPath = copy.deepcopy(self.pathList[ii])
					copyPath.pathMetric = PM[ii + numPaths]
					copyPath.uHat[pos] = 1
					self.pathList[ii].pathMetric = PM[ii]
					self.pathList.append(copyPath)
				elif pathStates[ii] == 2:
					self.pathList[ii].pathMetric = PM[ii + numPaths]
					self.pathList[ii].uHat[pos] = 1
				elif pathStates[ii] == 3:
					self.pathList[ii].pathMetric = PM[ii]
		else:
			discardTag = np.ones(2 * numPaths, dtype=int)
			pathStates = np.zeros(numPaths, dtype=int)
			self.sortNum += 1
			indices = sorted(range(2 * numPaths), key=lambda x: PM[x], reverse=True)
			indices = np.array(indices)
			discardTag[indices[:self.L]] = 0
			disIdxRecord = []
			for ii in range(numPaths):
				if discardTag[ii] == 1 and discardTag[ii + numPaths] == 1:
					#print('hi')
					pathStates[ii] = 1
					disIdxRecord.append(ii)
				elif discardTag[ii] == 1:
					pathStates[ii] = 2
				elif discardTag[ii + numPaths] == 1:
					pathStates[ii] = 3
				elif discardTag[ii] == 0 and discardTag[ii + numPaths] == 0:
					pathStates[ii] = 4
				else:
					input('state warning!')
			for ii in range(len(pathStates)):
				if pathStates[ii] == 2:
					self.pathList[ii].pathMetric = PM[ii + numPaths]
					self.pathList[ii].uHat[pos] = 1
				elif pathStates[ii] == 3:
					self.pathList[ii].pathMetric = PM[ii]
				elif pathStates[ii] == 4:
					if disIdxRecord:
						self.pathList[disIdxRecord[0]] = copy.deepcopy(self.pathList[ii])
						self.pathList[ii].pathMetric = PM[ii]
						self.pathList[disIdxRecord[0]].pathMetric = PM[ii + numPaths]
						self.pathList[disIdxRecord[0]].uHat[pos] = 1
						del disIdxRecord[0]
					else:
						copyPath = copy.deepcopy(self.pathList[ii])
						self.pathList[ii].pathMetric = PM[ii]
						copyPath.pathMetric = PM[ii + numPaths]
						copyPath.uHat[pos] = 1
						self.pathList.append(copyPath)
			if len(self.pathList) != self.L:
				input('list size warning')

	def PSCLN(self, softMess, isCRC=False, isRCPP=False): #若u1可靠，则整个subblock被SC；反之前N类SCL

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.mT = -19
		self.sortNum = 0
		self.unitCal = 0
		mySum = 0
		'''for i in self.criticalSet:
			count = 1
			print(i, end=' ')
			while self.polarMask[i + 1] == 1 and (i + 1) not in self.criticalSet:
				count += 1
				i += 1
				if i == self.codewordLength - 1:
					break
			print(count)'''
			#mySum += count
		#print(mySum)

		step4CS = []
		for i in self.criticalSet[:-1]:
			count = 0
			while self.polarMask[i + 1] == 1 and (i + 1) not in self.criticalSet:
				count += 1
				i += 1
				step4CS.append(i)
				if count == 2: #step size
					break
		print(len(step4CS), len(self.criticalSet))
		peSC = 1
		for i in range(self.codewordLength):
			if self.polarMask[i] == 1 and i not in step4CS and i not in self.criticalSet[:-1]:
				peSC *= (1 - self.pe[i])
		print(peSC, f'{1 - peSC:.4e}')
		input('234d')


		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1

			if self.polarMask[i] == 1:
				if i in self.criticalSet[:-1]:
					tmp = []
					for path in self.pathList:
						path.isSC = True 
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if min(branchMetric0, branchMetric1) >= self.mT:
							path.isSC = False
							copyPath = copy.deepcopy(path)
							path.pathMetric += branchMetric0
							copyPath.uHat[i] = 1
							copyPath.pathMetric += branchMetric1
							tmp.append(copyPath)
						else:
							if branchMetric0 > branchMetric1:
								path.pathMetric += branchMetric0
							else:
								path.pathMetric += branchMetric1
								path.uHat[i] = 1
					if (len(tmp) + len(self.pathList)) <= self.L:
						self.pathList += tmp
					else:
						tmp += self.pathList
						tmp.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
						self.sortNum += 1
						self.pathList = tmp[:self.L]

						if len(self.pathList) != self.L:
							input('warning')

				elif i in step4CS:
					tmp = []
					for path in self.pathList:
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if path.isSC:
							if branchMetric0 > branchMetric1:
								path.pathMetric += branchMetric0
							else:
								path.pathMetric += branchMetric1
								path.uHat[i] = 1
						else:
							if min(branchMetric0, branchMetric1) >= self.mT:
								copyPath = copy.deepcopy(path)
								path.pathMetric += branchMetric0
								copyPath.uHat[i] = 1
								copyPath.pathMetric += branchMetric1
								tmp.append(copyPath)
							else:
								if branchMetric0 > branchMetric1:
									path.pathMetric += branchMetric0
								else:
									path.pathMetric += branchMetric1
									path.uHat[i] = 1
					if (len(tmp) + len(self.pathList)) <= self.L:
						self.pathList += tmp
					else:
						tmp += self.pathList
						tmp.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
						self.sortNum += 1
						self.pathList = tmp[:self.L]

						if len(self.pathList) != self.L:
							input('warning')
				else:
					for path in self.pathList:
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if branchMetric0 > branchMetric1:
							path.pathMetric += branchMetric0
						else:
							path.pathMetric += branchMetric1
							path.uHat[i] = 1

			else:
				for path in self.pathList:
					path.uHat[i] = self.polarMask[i]
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)

		if len(self.pathList) == 0:
			input('decoding failure')
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			best = self.extract(self.pathList[0].uHat)
			candidate = pcfun.crcEncode(best, self.crcPoly)
			#candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = pcfun.crcEncode(nextBest, self.crcPoly)
					#nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def PSCL3(self, softMess, isCRC=False, isRCPP=False): #若u1可靠，则整个subblock被SC；反之类SCL

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.mT = -10
		self.sortNum = 0
		self.unitCal = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1

			if self.polarMask[i] == 1:
				if i in self.criticalSet[:-1]:
					tmp = []
					for path in self.pathList:
						path.isSC = True 
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if min(branchMetric0, branchMetric1) >= self.mT:
							path.isSC = False
							copyPath = copy.deepcopy(path)
							path.pathMetric += branchMetric0
							copyPath.uHat[i] = 1
							copyPath.pathMetric += branchMetric1
							tmp.append(copyPath)
						else:
							if branchMetric0 > branchMetric1:
								path.pathMetric += branchMetric0
							else:
								path.pathMetric += branchMetric1
								path.uHat[i] = 1
					if (len(tmp) + len(self.pathList)) <= self.L:
						self.pathList += tmp
					else:
						tmp += self.pathList
						tmp.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
						self.sortNum += 1
						self.pathList = tmp[:self.L]

						if len(self.pathList) != self.L:
							input('warning')

				else:
					if i >= self.criticalSet[-1]:
						for path in self.pathList:
							Li = path.LLRs[0] / np.log(2)
							branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
							branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
							if branchMetric0 > branchMetric1:
								path.pathMetric += branchMetric0
							else:
								path.pathMetric += branchMetric1
								path.uHat[i] = 1
					else:
						tmp = []
						for path in self.pathList:
							Li = path.LLRs[0] / np.log(2)
							branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
							branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
							if path.isSC:
								if branchMetric0 > branchMetric1:
									path.pathMetric += branchMetric0
								else:
									path.pathMetric += branchMetric1
									path.uHat[i] = 1
							else:
								if min(branchMetric0, branchMetric1) >= self.mT:
									copyPath = copy.deepcopy(path)
									path.pathMetric += branchMetric0
									copyPath.uHat[i] = 1
									copyPath.pathMetric += branchMetric1
									tmp.append(copyPath)
								else:
									if branchMetric0 > branchMetric1:
										path.pathMetric += branchMetric0
									else:
										path.pathMetric += branchMetric1
										path.uHat[i] = 1
						if (len(tmp) + len(self.pathList)) <= self.L:
							self.pathList += tmp
						else:
							tmp += self.pathList
							tmp.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
							self.sortNum += 1
							self.pathList = tmp[:self.L]

							if len(self.pathList) != self.L:
								input('warning')

			else:
				for path in self.pathList:
					path.uHat[i] = self.polarMask[i]
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)

		if len(self.pathList) == 0:
			input('decoding failure')
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			best = self.extract(self.pathList[0].uHat)
			candidate = pcfun.crcEncode(best, self.crcPoly)
			#candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = pcfun.crcEncode(nextBest, self.crcPoly)
					#nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def PSCL2(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.mT = -20
		self.sortNum = 0
		self.unitCal = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1

			if self.polarMask[i] == 1:
				self.PSCLFork2(i)
			else:

				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)

		if len(self.pathList) == 0:
			input('decoding failure')
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			best = self.extract(self.pathList[0].uHat)
			candidate = pcfun.crcEncode(best, self.crcPoly)
			#candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = pcfun.crcEncode(nextBest, self.crcPoly)
					#nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def PSCLFork(self, pos):

		tmpPathList = []
		discardIdx = []
		#path0Survive = False
		#path1Survive = False
		for i, path in enumerate(self.pathList):

			copyPath = copy.deepcopy(path)
			Li = path.LLRs[0] / np.log(2)
			branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
			if branchMetric0 >= self.mT:
				path.uHat[pos] = 0
				path.pathMetric += branchMetric0
				#path0Survive = True
			else:
				discardIdx.append(i)
				self.pathList[i] = 1

			branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
			if branchMetric1 >= self.mT:
				#copyPath = copy.deepcopy(path)
				copyPath.uHat[pos] = 1 
				copyPath.pathMetric += branchMetric1
				tmpPathList.append(copyPath)
				#path1Survive = True
			#if branchMetric0 >= self.mT and branchMetric1 >= self.mT:
			#	self.count += 1

		while 1 in self.pathList:
			self.pathList.remove(1)

		if (len(tmpPathList) + len(self.pathList)) <= self.L:
			self.pathList += tmpPathList
		else:

			tmpPathList += self.pathList
			tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
			self.sortNum += 1

			self.pathList = tmpPathList[:self.L]
			if len(self.pathList) != self.L:
				input('warning')

	def PSCL(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.mT = -12
		self.sortNum = 0
		self.unitCal = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1

			if self.polarMask[i] == 1:
				self.PSCLFork(i)
			else:

				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)

		if len(self.pathList) == 0:
			input('decoding failure')
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			#self.crc8Table = pcfun.buildCRC8Table(self.crcPoly)
			best = self.extract(self.pathList[0].uHat)
			candidate = pcfun.crcEncode(best, self.crcPoly)
			#candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			#if (np.sum(pcfun.crc8TableMethod(best, self.crc8Table))) == 0:
			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = pcfun.crcEncode(nextBest, self.crcPoly)
					#nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def CAPSCL2(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.mT = -20
		self.sortNum = 0
		self.unitCal = 0
		self.count = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1

			if self.polarMask[i] == 1:
				if i in self.criticalSet[:-1]:
					self.PSCLFork2(i)
				else:
					for path in self.pathList:
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if branchMetric0 > branchMetric1:
							path.uHat[i] = 0
							path.pathMetric += branchMetric0
						elif branchMetric1 > branchMetric0:
							path.uHat[i] = 1
							path.pathMetric += branchMetric1
						else:
							input('b0 = b1')
			else:

				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)

		if len(self.pathList) == 0:
			input('decoding failure')
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			#self.crc8Table = pcfun.buildCRC8Table(self.crcPoly)
			best = self.extract(self.pathList[0].uHat)
			candidate = pcfun.crcEncode(best, self.crcPoly)
			#candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			#if (np.sum(pcfun.crc8TableMethod(best, self.crc8Table))) == 0:
			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = pcfun.crcEncode(nextBest, self.crcPoly)
					#nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def CAPSCL(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.mT = -14
		self.sortNum = 0
		self.unitCal = 0
		self.count = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1

			if self.polarMask[i] == 1:
				if i in self.criticalSet[:-1]:
					self.PSCLFork(i)
				else:
					for path in self.pathList:
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if branchMetric0 > branchMetric1:
							path.uHat[i] = 0
							path.pathMetric += branchMetric0
						elif branchMetric1 > branchMetric0:
							path.uHat[i] = 1
							path.pathMetric += branchMetric1
						else:
							input('b0 = b1')
			else:

				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)

		if len(self.pathList) == 0:
			input('decoding failure')
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			#self.crc8Table = pcfun.buildCRC8Table(self.crcPoly)
			best = self.extract(self.pathList[0].uHat)
			candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			#if (np.sum(pcfun.crc8TableMethod(best, self.crc8Table))) == 0:
			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def DPSCLN(self, softMess, isCRC=False, isRCPP=False): #若u1可靠，则整个subblock被SC；反之前N类SCL

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.mT = -15
		self.sortNum = 0
		self.unitCal = 0
		alpha = 1.0

		step4CS = []
		for i in self.criticalSet[:-1]:
			count = 0
			while self.polarMask[i + 1] == 1 and (i + 1) not in self.criticalSet:
				count += 1
				i += 1
				step4CS.append(i)
				if count == 2: #step size
					break

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)
				self.unitCal += path.unitCal
				self.ANV += 1

			if self.polarMask[i] == 1:
				if i in self.criticalSet[:-1]:
					tmp = []
					for path in self.pathList:
						path.isSC = True 
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						#if min(branchMetric0, branchMetric1) > -self.GA[i]:
						if min(branchMetric0, branchMetric1) > 1 - np.log2(1 + 2 ** (-self.GA[i] * ((-1) ** 1))):
							path.isSC = False
							copyPath = copy.deepcopy(path)
							path.pathMetric += branchMetric0
							copyPath.uHat[i] = 1
							copyPath.pathMetric += branchMetric1
							tmp.append(copyPath)
						else:
							if branchMetric0 > branchMetric1:
								path.pathMetric += branchMetric0
							else:
								path.pathMetric += branchMetric1
								path.uHat[i] = 1
					if (len(tmp) + len(self.pathList)) <= self.L:
						self.pathList += tmp
					else:
						tmp += self.pathList
						tmp.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
						self.sortNum += 1
						self.pathList = tmp[:self.L]

						if len(self.pathList) != self.L:
							input('warning')

				elif i in step4CS:
					tmp = []
					for path in self.pathList:
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if path.isSC:
							if branchMetric0 > branchMetric1:
								path.pathMetric += branchMetric0
							else:
								path.pathMetric += branchMetric1
								path.uHat[i] = 1
						else:
							if min(branchMetric0, branchMetric1) > self.mT:
								copyPath = copy.deepcopy(path)
								path.pathMetric += branchMetric0
								copyPath.uHat[i] = 1
								copyPath.pathMetric += branchMetric1
								tmp.append(copyPath)
							else:
								if branchMetric0 > branchMetric1:
									path.pathMetric += branchMetric0
								else:
									path.pathMetric += branchMetric1
									path.uHat[i] = 1
					if (len(tmp) + len(self.pathList)) <= self.L:
						self.pathList += tmp
					else:
						tmp += self.pathList
						tmp.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
						self.sortNum += 1
						self.pathList = tmp[:self.L]

						if len(self.pathList) != self.L:
							input('warning')
				else:
					for path in self.pathList:
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if branchMetric0 > branchMetric1:
							path.pathMetric += branchMetric0
						else:
							path.pathMetric += branchMetric1
							path.uHat[i] = 1

			else:
				for path in self.pathList:
					path.uHat[i] = self.polarMask[i]
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)

		if len(self.pathList) == 0:
			input('decoding failure')
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			best = self.extract(self.pathList[0].uHat)
			candidate = pcfun.crcEncode(best, self.crcPoly)
			#candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = pcfun.crcEncode(nextBest, self.crcPoly)
					#nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def DPSCL(self, softMess, isCRC=False, isRCPP=False):

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.mT = -12
		self.dT = 1 * 1e-5
		self.sortNum = 0
		alpha = 1.5
		count = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)

			if self.polarMask[i] == 1:
				
				tmpPathList = []
				discardIdx = []
				if i in self.criticalSet[:-1]:
					for ii, path in enumerate(self.pathList):
						mllr = abs(path.LLRs[0])
						if mllr < alpha * abs(np.log((1- self.pe[i]) / self.pe[i])):
							
							copyPath = copy.deepcopy(path)
							Li = path.LLRs[0] / np.log(2)
							branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
							if branchMetric0 >= self.mT:
								path.uHat[i] = 0
								path.pathMetric += branchMetric0
							else:
								discardIdx.append(ii)
								self.pathList[ii] = 1

							branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
							if branchMetric1 >= self.mT:
								#copyPath = copy.deepcopy(path)
								copyPath.uHat[i] = 1 
								copyPath.pathMetric += branchMetric1
								tmpPathList.append(copyPath)

						else:
							Li = path.LLRs[0] / np.log(2)
							branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
							branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
							if branchMetric0 > branchMetric1:
								path.uHat[i] = 0
								path.pathMetric += branchMetric0
							elif branchMetric1 > branchMetric0:
								path.uHat[i] = 1
								path.pathMetric += branchMetric1
							else:
								input('b0 = b1')
				else:
					for path in self.pathList:
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						if branchMetric0 > branchMetric1:
							path.uHat[i] = 0
							path.pathMetric += branchMetric0
						elif branchMetric1 > branchMetric0:
							path.uHat[i] = 1
							path.pathMetric += branchMetric1
						else:
							input('b0 = b1')

				while 1 in self.pathList:
					self.pathList.remove(1)

				if (len(tmpPathList) + len(self.pathList)) <= self.L:
					self.pathList += tmpPathList
				else:

					tmpPathList += self.pathList
					tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
					self.sortNum += 1

					self.pathList = tmpPathList[:self.L]
					if len(self.pathList) != self.L:
						input('warning')
			else:

				for path in self.pathList:

					path.uHat[i] = self.polarMask[i]
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)

		if len(self.pathList) == 0:
			input('decoding failure')
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			#self.crc8Table = pcfun.buildCRC8Table(self.crcPoly)
			best = self.extract(self.pathList[0].uHat)
			candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			#if (np.sum(pcfun.crc8TableMethod(best, self.crc8Table))) == 0:
			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def SRSCL2(self, softMess, isCRC=False, isRCPP=False): #dot caution

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.sortNum = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)

			if self.polarMask[i] == 1:
				
				tmpPathList = []
				for path in self.pathList:

					if path.LLRs[0] > np.log((1 - self.pe[i]) / self.pe[i]):
						path.uHat[i] = 0
						Li = path.LLRs[0] / np.log(2)
						branchMetric = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						path.pathMetric += branchMetric
					elif path.LLRs[0] < -np.log((1 - self.pe[i]) / self.pe[i]):
						path.uHat[i] = 1
						Li = path.LLRs[0] / np.log(2)
						branchMetric = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						path.pathMetric += branchMetric
					else:
						
						copyPath = copy.deepcopy(path)
						path.uHat[i] = 0
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						path.pathMetric += branchMetric0
						copyPath.uHat[i] = 1 
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						copyPath.pathMetric += branchMetric1

						tmpPathList.append(copyPath)

				if (len(tmpPathList) + len(self.pathList)) <= self.L:
					self.pathList += tmpPathList
				else:

					tmpPathList += self.pathList
					#tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=False)
					tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
					self.sortNum += 1
					self.pathList = tmpPathList[:self.L]

			else:

				for path in self.pathList:

					#path.omiga += 1
					path.uHat[i] = self.polarMask[i]
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)

		if len(self.pathList) == 0:
			input('decoding failure')
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			best = self.extract(self.pathList[0].uHat)
			candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def SRSCL(self, softMess, isCRC=False, isRCPP=False): #dot caution

		self.pathList = [Path(self.codewordLength)]
		if isRCPP:
			(self.pathList[0].LLRs[self.codewordLength - 1:])[self.p == 1] = softMess
		else:
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.sortNum = 0
		omiga = np.inf
		alpha = 1.5

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.pathList:
				path.updateLLRs(idxRev)

			if self.polarMask[i] == 1:
				
				tmpPathList = []
				for path in self.pathList:
					absLLR = abs(path.LLRs[0])
					if absLLR > abs(alpha * np.log((1 - self.pe[i]) / self.pe[i])):
						path.omiga += 1
						path.uHat[i] = 0 if path.LLRs[0] > alpha * np.log((1 - self.pe[i]) / self.pe[i]) else 1
						Li = path.LLRs[0] / np.log(2)
						branchMetric = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0))) if path.LLRs[0] > alpha * np.log((1 - self.pe[i]) / self.pe[i]) else 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						path.pathMetric += branchMetric
					else:
						
						copyPath = copy.deepcopy(path)
						path.uHat[i] = 0
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						path.pathMetric += branchMetric0
						path.omiga = 0
						copyPath.uHat[i] = 1 
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 1)))
						copyPath.pathMetric += branchMetric1
						copyPath.omiga = 0

						tmpPathList.append(copyPath)

				if (len(tmpPathList) + len(self.pathList)) <= self.L:
					self.pathList += tmpPathList
				else:

					tmpPathList += self.pathList
					surviveNum = 0
					for i, path in enumerate(tmpPathList):
						if path.omiga > omiga:
							surviveNum += 1
					if surviveNum == 0:
						#tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=False)
						tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
						self.sortNum += 1
						self.pathList = tmpPathList[:self.L]
					else:
						for i, path in enumerate(tmpPathList):
							if path.omiga < omiga:
								del tmpPathList[i]
						if len(tmpPathList) > self.L:
							#tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=False)
							tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
							self.sortNum += 1
							self.pathList = tmpPathList[:self.L]
						else:
							self.pathList = tmpPathList[:]
			else:

				for path in self.pathList:

					#path.omiga += 1
					path.uHat[i] = self.polarMask[i]
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
					path.pathMetric += branchMetric0

			for path in self.pathList:
				path.updateBits(idxRev)

		if len(self.pathList) == 0:
			input('decoding failure')
		self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
		self.sortNum += 1

		if isCRC:
			
			best = self.extract(self.pathList[0].uHat)
			candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]
				return best[:len(best) - self.crcWidth]

		else:
			best = self.pathList[0].uHat

		return self.extract(best)

	def pacSCDecoder(self, softMess):

		'''
		self.LLRs = np.zeros(2 * self.codewordLength - 1, dtype=float)
		self.BITs = np.zeros((2, self.codewordLength - 1), dtype=int)
		self.LLRs[self.codewordLength - 1:] = softMess
		self.vHat = np.zeros(self.codewordLength, dtype=int)
		self.curState = [0 for i in range(self.m)]
		'''
		path = Path(self.codewordLength, self.m)
		path.LLRs[self.codewordLength - 1:] = softMess
		self.errCount = 0

		for i in range(self.codewordLength):

			ii = pcfun.bitreversed(i, self.n)
			path.updateLLRs(ii)

			if self.polarMask[i] == 1:

				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				edgeValue1 = pcfun.conv1Bit(1, path.curState, self.gen)
				'''
				branchMetric0 = np.log(1 + np.exp(-(1 - 2 * edgeValue0) * path.LLRs[0]))
				branchMetric1 = np.log(1 + np.exp(-(1 - 2 * edgeValue1) * path.LLRs[0]))
				'''
				
				penalty = np.abs(path.LLRs[0])
				if path.LLRs[0] > 0:

					branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					branchMetric1 = 0 if (edgeValue1 == 0) else penalty
				elif path.LLRs[0] < 0:

					branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					branchMetric1 = 0 if (edgeValue1 == 1) else penalty
				else:
					input('warning')
				
				if branchMetric0 < branchMetric1:

					path.vHat[i] = 0
					edgeValue = edgeValue0
					path.curState = pcfun.getNextState(0, path.curState, self.m)

				elif branchMetric1 < branchMetric0:

					path.vHat[i] = 1
					edgeValue = edgeValue1
					path.curState = pcfun.getNextState(1, path.curState, self.m)

				else:
					input('warning')

			else:

				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				edgeValue = edgeValue0
				path.vHat[i] = self.polarMask[i]
				path.curState = pcfun.getNextState(0, path.curState, self.m)

			path.uHat[i] = edgeValue
			path.updateBits(ii)
		
		return self.extract(path.vHat)

	def oracleAssistPolarDecoder(self, softMess):

		trellisPath = Path(self.codewordLength)
		trellisPath.LLRs[self.codewordLength - 1:] = softMess
		mu = pcfun.GA(self.codewordLength, self.infoLen, self.designSNR)
		idxRecord = []
		for i in range(self.codewordLength):
			if self.polarMask[i] == 1 and self.GA[i] < 32:
				idxRecord.append(i)
		print(idxRecord, len(idxRecord))
		input('2')
		self.errCount = 0
		self.errorIndex = []
		self.predictIdx = []
		self.predictIdx2 = []
		self.recordLLRs = []
		self.mT = [-11] + [-8 for _ in range(self.maxLevel - 1)]
		self.pauseFlag = False
		
		print('new line ---------------')
		for i in range(self.codewordLength):

			ii = pcfun.bitreversed(i, self.n)
			trellisPath.updateLLRs(ii)
			self.recordLLRs.append(trellisPath.LLRs[0])

			if self.polarMask[i] == 1:

				if trellisPath.LLRs[0] > 0:
					trellisPath.uHat[i] = 0
				elif trellisPath.LLRs[0] < 0:
					trellisPath.uHat[i] = 1
				else:
					input('Warning')

			else:
				trellisPath.uHat[i] = 0

			if trellisPath.uHat[i] != self.message[i]:
				trellisPath.uHat[i] = self.message[i]
				Li = trellisPath.LLRs[0] / (np.log(2))
				branchMetric0 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 0)))
				branchMetric1 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 1)))
				print(i, min(branchMetric0, branchMetric1), i in self.criticalSet, i in idxRecord)
				self.errCount += 1
				if i not in idxRecord:
					self.pauseFlag = True

			'''Li = trellisPath.LLRs[0] / (np.log(2))
			branchMetric0 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 0)))
			branchMetric1 = 1 - np.log2(1 + pow(2, -Li * ((-1) ** 1)))
			#if branchMetric0 > self.mT and branchMetric1 > self.mT:
			if branchMetric0 > self.mT[self.errCount] and branchMetric1 > self.mT[self.errCount]:
				if i in self.criticalSet[:-1]:
					self.predictIdx.append(i)'''
				#if self.polarMask[i] == 1:
				#	self.predictIdx2.append(i)

			'''if self.errCount < self.maxLevel and trellisPath.uHat[i] != self.message[i]:
				self.errCount += 1
				self.errorIndex.append(i)
				trellisPath.uHat[i] = self.message[i]
				print(i, min(branchMetric0, branchMetric1), i in self.criticalSet)'''

			'''if self.errCount < self.maxLevel and trellisPath.uHat[i] != self.message[i]:
				self.errCount += 1
				self.errorIndex.append(i)
				trellisPath.uHat[i] = self.message[i]
				print(i, min(branchMetric0, branchMetric1))'''

			trellisPath.updateBits(ii)
		#input('4231')
		if self.pauseFlag:
			for i in self.criticalSet:
				count = 1
				print(i, end=' ')
				while self.polarMask[i + 1] == 1 and (i + 1) not in self.criticalSet:
					count += 1
					i += 1
					if i == self.codewordLength - 1:
						break
				print(count)
			input('pause')
		'''self.recordLLRs = np.array(self.recordLLRs)
		print(self.errorIndex)
		print(self.predictIdx)
		print(self.criticalSet)
		input('pp')
		sortedIndex = sorted(self.criticalSet, key=lambda x: np.abs(self.recordLLRs[x] / np.sqrt(mu[x])))
		sortedIndex = np.array(sortedIndex, dtype=int)
		print(sortedIndex)
		#input('pp')
		for err in self.errorIndex:
			if err in self.predictIdx:
				for i, idx in enumerate(self.predictIdx):
					if idx == err:
						print(i + 1)
						break
			else:
				print(-1)
		#print(len(self.predictIdx2))
		print(len(self.predictIdx))
		#print(self.recordLLRs[np.array(self.predictIdx)])
		input('p')'''
		#return self.extract(trellisPath.uHat)
		return self.extract(trellisPath.uHat)[:-self.crcWidth]

	def oracleAssistPACDecoder(self, softMess):

		path = Path(self.codewordLength, self.m)
		path.LLRs[self.codewordLength - 1:] = softMess
		self.errCount = 0

		for i in range(self.codewordLength):

			ii = pcfun.bitreversed(i, self.n)
			path.updateLLRs(ii)

			if self.polarMask[i] == 1:

				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				edgeValue1 = pcfun.conv1Bit(1, path.curState, self.gen)
				
				branchMetric0 = np.log(1 + np.exp(-(1 - 2 * edgeValue0) * path.LLRs[0]))
				branchMetric1 = np.log(1 + np.exp(-(1 - 2 * edgeValue1) * path.LLRs[0]))
				
				'''
				penalty = np.abs(path.LLRs[0])
				if path.LLRs[0] > 0:

					branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					branchMetric1 = 0 if (edgeValue1 == 0) else penalty
				elif path.LLRs[0] < 0:

					branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					branchMetric1 = 0 if (edgeValue1 == 1) else penalty
				else:
					input('warning')
				'''
				if branchMetric0 < branchMetric1:

					path.vHat[i] = 0
					edgeValue = edgeValue0
					#path.curState = pcfun.getNextState(0, path.curState, self.m)

				elif branchMetric1 < branchMetric0:

					path.vHat[i] = 1
					edgeValue = edgeValue1
					#path.curState = pcfun.getNextState(1, path.curState, self.m)

				else:
					input('warning')

			else:

				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				edgeValue = edgeValue0
				path.vHat[i] = self.polarMask[i]
				

			
			if path.vHat[i] != self.message[i]:

				path.vHat[i] = self.message[i]
				path.uHat[i] = 1 - edgeValue
				path.curState = pcfun.getNextState(path.vHat[i], path.curState, self.m)
				path.updateBits(ii)
				self.errCount += 1

			else:
				path.uHat[i] = edgeValue
				path.curState = pcfun.getNextState(path.vHat[i], path.curState, self.m)
				path.updateBits(ii)
			
		return self.extract(path.vHat)

	def pacSCLPathFork2(self, pos):

		numPaths = len(self.trellisPathList)
		if 2 * numPaths <= self.L:
			for ii in range(numPaths):
				edgeValue0 = pcfun.conv1Bit(0, self.trellisPathList[ii].curState, self.gen)
				edgeValue1 = pcfun.conv1Bit(1, self.trellisPathList[ii].curState, self.gen)
				pathState0 = pcfun.getNextState(0, self.trellisPathList[ii].curState, self.m)
				pathState1 = pcfun.getNextState(1, self.trellisPathList[ii].curState, self.m)
				copyPath = copy.deepcopy(self.trellisPathList[ii])
				Li = copyPath.LLRs[0] / np.log(2)
				branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
				self.trellisPathList[ii].pathMetric += branchMetric0
				self.trellisPathList[ii].uHat[pos] = edgeValue0
				self.trellisPathList[ii].curState = pathState0
				
				branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue1)))
				copyPath.pathMetric += branchMetric1
				copyPath.vHat[pos] = 1
				copyPath.uHat[pos] = edgeValue1
				copyPath.curState = pathState1
				self.trellisPathList.append(copyPath)
		else:
			edgeValues = np.zeros(2 * numPaths, dtype=int)
			nextStates = [[] for _ in range(2 * numPaths)]
			PM = np.zeros(2 * numPaths, dtype=float)
			discardTag = np.ones(2 * numPaths, dtype=int)
			pathStates = np.zeros(numPaths, dtype=int)
			for ii in range(numPaths):
				edgeValue0 = pcfun.conv1Bit(0, self.trellisPathList[ii].curState, self.gen)
				edgeValue1 = pcfun.conv1Bit(1, self.trellisPathList[ii].curState, self.gen)
				#pathState0 = pcfun.getNextState(0, self.trellisPathList[ii].curState, self.m)
				#pathState1 = pcfun.getNextState(1, self.trellisPathList[ii].curState, self.m)
				Li = self.trellisPathList[ii].LLRs[0] / np.log(2)
				branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
				branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue1)))
				PM[ii] = self.trellisPathList[ii].pathMetric + branchMetric0
				PM[ii + numPaths] = self.trellisPathList[ii].pathMetric + branchMetric1
				edgeValues[ii] = edgeValue0
				edgeValues[ii + numPaths] = edgeValue1
				nextStates[ii] = pcfun.getNextState(0, self.trellisPathList[ii].curState, self.m)
				nextStates[ii + numPaths] = pcfun.getNextState(1, self.trellisPathList[ii].curState, self.m)
			self.sortNum += 1
			indices = sorted(range(2 * numPaths), key=lambda x: PM[x], reverse=True)
			indices = np.array(indices)
			discardTag[indices[:self.L]] = 0
			disIdxRecord = []
			for ii in range(numPaths):
				if discardTag[ii] == 1 and discardTag[ii + numPaths] == 1:
					pathStates[ii] = 1
					disIdxRecord.append(ii)
				elif discardTag[ii] == 1:
					pathStates[ii] = 2
				elif discardTag[ii + numPaths] == 1:
					pathStates[ii] = 3
				elif discardTag[ii] == 0 and discardTag[ii + numPaths] == 0:
					pathStates[ii] = 4
				else:
					input('state warning!')
			for ii in range(len(pathStates)):
				if pathStates[ii] == 2:
					self.trellisPathList[ii].pathMetric = PM[ii + numPaths]
					self.trellisPathList[ii].vHat[pos] = 1
					self.trellisPathList[ii].uHat[pos] = edgeValues[ii + numPaths]
					self.trellisPathList[ii].curState = nextStates[ii + numPaths]
				elif pathStates[ii] == 3:
					self.trellisPathList[ii].pathMetric = PM[ii]
					self.trellisPathList[ii].uHat[pos] = edgeValues[ii]
					self.trellisPathList[ii].curState = nextStates[ii]
				elif pathStates[ii] == 4:
					self.trellisPathList[disIdxRecord[0]] = copy.deepcopy(self.trellisPathList[ii])
					self.trellisPathList[ii].pathMetric = PM[ii]
					self.trellisPathList[ii].uHat[pos] = edgeValues[ii]
					self.trellisPathList[ii].curState = nextStates[ii]
					self.trellisPathList[disIdxRecord[0]].pathMetric = PM[ii + numPaths]
					self.trellisPathList[disIdxRecord[0]].vHat[pos] = 1
					self.trellisPathList[disIdxRecord[0]].uHat[pos] = edgeValues[ii + numPaths]
					self.trellisPathList[disIdxRecord[0]].curState = nextStates[ii + numPaths]
					del disIdxRecord[0]

	def pacSCLDecoder2(self, softMess):

		self.trellisPathList = [Path(self.codewordLength, self.m)]
		self.trellisPathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.sortNum = 0
		self.runningTime = 0

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.trellisPathList:

				path.updateLLRs(idxRev)
				#if self.polarMask[i] == 1:
				#	self.ANV += 1

			if self.polarMask[i] == 1:
				self.pacSCLPathFork2(i)
			else:

				for path in self.trellisPathList:

					edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
					path.uHat[i] = edgeValue0
					path.vHat[i] = self.polarMask[i]
					path.curState = pcfun.getNextState(0, path.curState, self.m)
					'''penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] > 0:
						branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					elif path.LLRs[0] < 0:
						branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					else:
						input('warning')'''
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
					path.pathMetric += branchMetric0
					#self.ANV += 1

			for path in self.trellisPathList:
				path.updateBits(idxRev)

		#self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=False)
		self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=True)
		self.sortNum += 1
		best = self.trellisPathList[0].vHat

		return self.extract(best)

	def PACPSCLFork(self, pos):

		tmpPathList = []
		for i, path in enumerate(self.trellisPathList):

			copyPath = copy.deepcopy(path)
			edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
			Li = path.LLRs[0] / np.log(2)
			branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
			if branchMetric0 >= self.mT:
				path.uHat[pos] = edgeValue0
				path.curState = pcfun.getNextState(0, path.curState, self.m)
				path.pathMetric += branchMetric0
			else:
				self.trellisPathList[i] = 1

			edgeValue1 = pcfun.conv1Bit(1, copyPath.curState, self.gen)
			branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue1)))
			if branchMetric1 >= self.mT:
				#copyPath = copy.deepcopy(path)
				copyPath.vHat[pos] = 1 
				copyPath.uHat[pos] = edgeValue1
				copyPath.curState = pcfun.getNextState(1, copyPath.curState, self.m)
				copyPath.pathMetric += branchMetric1
				tmpPathList.append(copyPath)

		while 1 in self.trellisPathList:
			self.trellisPathList.remove(1)

		if (len(tmpPathList) + len(self.trellisPathList)) <= self.L:
			self.trellisPathList += tmpPathList
		else:

			tmpPathList += self.trellisPathList
			tmpPathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
			self.sortNum += 1

			self.trellisPathList = tmpPathList[:self.L]
			if len(self.trellisPathList) != self.L:
				input('warning')

	def PACCAPSCL(self, softMess):

		self.trellisPathList = [Path(self.codewordLength, self.m)]
		self.trellisPathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.sortNum = 0
		self.runningTime = 0
		self.mT = -12

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.trellisPathList:

				path.updateLLRs(idxRev)

			if self.polarMask[i] == 1:
				if i in self.criticalSet[:-1]:
					self.PACPSCLFork(i)
				else:
					for path in self.trellisPathList:
						edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
						edgeValue1 = pcfun.conv1Bit(1, path.curState, self.gen)
						curState0 = pcfun.getNextState(0, path.curState, self.m)
						curState1 = pcfun.getNextState(1, path.curState, self.m)
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue1)))
						if branchMetric0 > branchMetric1:
							path.uHat[i] = edgeValue0
							path.curState = curState0
							path.pathMetric += branchMetric0
						else:
							path.vHat[i] = 1
							path.uHat[i] = edgeValue1
							path.curState = curState1
							path.pathMetric += branchMetric1
			else:

				for path in self.trellisPathList:

					edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
					path.uHat[i] = edgeValue0
					path.vHat[i] = self.polarMask[i]
					path.curState = pcfun.getNextState(0, path.curState, self.m)
					'''penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] > 0:
						branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					elif path.LLRs[0] < 0:
						branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					else:
						input('warning')'''
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
					path.pathMetric += branchMetric0

			for path in self.trellisPathList:
				path.updateBits(idxRev)

		#self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=False)
		self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=True)
		self.sortNum += 1
		best = self.trellisPathList[0].vHat

		return self.extract(best)

	def PACGCAPSCL(self, softMess):

		self.trellisPathList = [Path(self.codewordLength, self.m)]
		self.trellisPathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.sortNum = 0
		self.runningTime = 0
		self.mT = -11
		idxRecord = []
		for i in range(self.codewordLength):
			if self.polarMask[i] == 1 and self.GA[i] < 11:
				idxRecord.append(i)

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.trellisPathList:

				path.updateLLRs(idxRev)

			if self.polarMask[i] == 1:
				if i in idxRecord:
					self.PACPSCLFork(i)
				else:
					for path in self.trellisPathList:
						edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
						edgeValue1 = pcfun.conv1Bit(1, path.curState, self.gen)
						curState0 = pcfun.getNextState(0, path.curState, self.m)
						curState1 = pcfun.getNextState(1, path.curState, self.m)
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
						branchMetric1 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue1)))
						if branchMetric0 > branchMetric1:
							path.uHat[i] = edgeValue0
							path.curState = curState0
							path.pathMetric += branchMetric0
						else:
							path.vHat[i] = 1
							path.uHat[i] = edgeValue1
							path.curState = curState1
							path.pathMetric += branchMetric1
			else:

				for path in self.trellisPathList:

					edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
					path.uHat[i] = edgeValue0
					path.vHat[i] = self.polarMask[i]
					path.curState = pcfun.getNextState(0, path.curState, self.m)
					'''penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] > 0:
						branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					elif path.LLRs[0] < 0:
						branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					else:
						input('warning')'''
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
					path.pathMetric += branchMetric0

			for path in self.trellisPathList:
				path.updateBits(idxRev)

		#self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=False)
		self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=True)
		self.sortNum += 1
		best = self.trellisPathList[0].vHat

		return self.extract(best)

	def PACPSCL(self, softMess):

		self.trellisPathList = [Path(self.codewordLength, self.m)]
		self.trellisPathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.sortNum = 0
		self.runningTime = 0
		self.mT = -11

		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.trellisPathList:

				path.updateLLRs(idxRev)

			if self.polarMask[i] == 1:
				self.PACPSCLFork(i)
			else:

				for path in self.trellisPathList:

					edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
					path.uHat[i] = edgeValue0
					path.vHat[i] = self.polarMask[i]
					path.curState = pcfun.getNextState(0, path.curState, self.m)
					'''penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] > 0:
						branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					elif path.LLRs[0] < 0:
						branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					else:
						input('warning')'''
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
					path.pathMetric += branchMetric0

			for path in self.trellisPathList:
				path.updateBits(idxRev)

		#self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=False)
		self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=True)
		self.sortNum += 1
		best = self.trellisPathList[0].vHat

		return self.extract(best)

	def pacSCLDecoder(self, softMess):

		self.trellisPathList = [Path(self.codewordLength, self.m)]
		self.trellisPathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.runningTime = 0

		#t1 = perf_counter()
		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.trellisPathList:

				path.updateLLRs(idxRev)
				#if self.polarMask[i] == 1:
				#	self.ANV += 1

			if self.polarMask[i] == 1:
				self.pacSCLPathFork(i)
			else:

				for path in self.trellisPathList:

					edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
					path.uHat[i] = edgeValue0
					path.vHat[i] = self.polarMask[i]
					path.curState = pcfun.getNextState(0, path.curState, self.m)
					'''penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] > 0:
						branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					elif path.LLRs[0] < 0:
						branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					else:
						input('warning')'''
					Li = path.LLRs[0] / np.log(2)
					branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** edgeValue0)))
					path.pathMetric += branchMetric0
					#self.ANV += 1

			for path in self.trellisPathList:
				path.updateBits(idxRev)

		#self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=False)
		self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=True)
		best = self.trellisPathList[0].vHat
		#t2 = perf_counter()
		#self.runningTime = (t2 - t1) * 1000

		return self.extract(best)

	def pacStackDecoder(self, softMess):

		self.trellisPathStack = [Path(self.codewordLength, self.m)]
		self.T = 1 
		self.trellisPathStack[0].LLRs[self.codewordLength - 1:] = softMess
		self.trellisPathStack[0].currPosition = -1
		self.visitsCount = np.zeros(self.codewordLength, dtype=int)
		flag = True

		while flag:

			currPath = self.trellisPathStack.pop()
			self.T -= 1

			if self.polarMask[currPath.currPosition + 1] == 1:

				currPath.currPosition += 1
				self.visitsCount[currPath.currPosition] += 2
				self.pacStackPathFork(currPath)

			else:

				currPath.currPosition += 1
				self.visitsCount[currPath.currPosition] += 1
				idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
				currPath.updateLLRs(idxRev)
				edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
				currPath.curState = pcfun.getNextState(0, currPath.curState, self.m)
				currPath.uHat[currPath.currPosition] = edgeValue0
				currPath.vHat[currPath.currPosition] = self.polarMask[currPath.currPosition]
				'''
				branchMetric0 = np.log(1 + np.exp(-(1 - 2 * edgeValue0) * currPath.LLRs[0]))
				'''
				'''
				penalty = np.abs(currPath.LLRs[0])
				if currPath.LLRs[0] > 0:
					branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					
				elif currPath.LLRs[0] < 0:
					branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					
				else:
					input('warning')
				'''
				penalty = np.abs(currPath.LLRs[0])
				penalty /= np.log(2)
				if currPath.LLRs[0] > 0:
					branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.cutoff_rate[currPath.currPosition]
				elif currPath.LLRs[0] < 0:
					branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.cutoff_rate[currPath.currPosition]
				else:
					input('warning')
				currPath.pathMetric += branchMetric0
				currPath.updateBits(idxRev)
				self.trellisPathStack.append(currPath)
				self.T += 1
				self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

			if self.trellisPathStack[-1].currPosition == (self.codewordLength - 1):

				best = self.trellisPathStack[-1].vHat
				flag = False
				

		return self.extract(best)

	def pacListViterbiDecoder(self, softMess):

		self.trellisPathList = [Path(self.codewordLength, self.m)]
		self.trellisPathList[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.runningTime = 0

		t1 = perf_counter()
		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			for path in self.trellisPathList:

				path.updateLLRs(idxRev)
				if self.polarMask[i] == 1:
					self.ANV += 1

			if self.polarMask[i] == 1:
				self.listViterbiFork(i)
			else:

				for path in self.trellisPathList:

					edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
					path.uHat[i] = edgeValue0
					path.vHat[i] = self.polarMask[i]
					path.curState = pcfun.getNextState(0, path.curState, self.m)
					penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] > 0:
						branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					elif path.LLRs[0] < 0:
						branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					else:
						input('warning')
					path.pathMetric += branchMetric0
					self.ANV += 1

			for path in self.trellisPathList:
				path.updateBits(idxRev)

		self.trellisPathList.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=False)
		#print(len(self.trellisPathList))
		best = self.trellisPathList[0].vHat
		t2 = perf_counter()
		self.runningTime = (t2 - t1) * 1000
		return self.extract(best)

	def pacStackViterbiDecoder(self, softMess):

		self.trellisPathStack = [Path(self.codewordLength, self.m)]
		self.T = 1 
		self.trellisPathStack[0].LLRs[self.codewordLength - 1:] = softMess
		self.trellisPathStack[0].currPosition = -1
		self.criticalMask = np.zeros(self.codewordLength, dtype=int)
		self.criticalMask[self.criticalSet] = 1
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35
		self.ANV = 0
		self.runningTime = 0
		flag = True

		t1 = perf_counter()
		while flag:

			currPath = self.trellisPathStack.pop()
			self.T -= 1

			if self.polarMask[currPath.currPosition + 1] == 1:

				currPath.currPosition += 1
				self.stackViterbiFork1(currPath)
				'''
				else:

					currPath.currPosition += 1
					idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
					currPath.updateLLRs(idxRev)
					edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
					edgeValue1 = pcfun.conv1Bit(1, currPath.curState, self.gen)
					penalty = np.abs(currPath.LLRs[0])
					penalty /= np.log(2)
					if currPath.LLRs[0] > 0:
						branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.cutoff_rate[currPath.currPosition]
						branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.cutoff_rate[currPath.currPosition]
					elif currPath.LLRs[0] < 0:
						branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.cutoff_rate[currPath.currPosition]
						branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.cutoff_rate[currPath.currPosition]
					else:
						input('warning')
					currPath.pathMetric += branchMetric0 if (branchMetric0 > branchMetric1) else branchMetric1
					currPath.vHat[currPath.currPosition] = 0 if (branchMetric0 > branchMetric1) else 1 
					currPath.uHat[currPath.currPosition] = edgeValue0 if (branchMetric0 > branchMetric1) else edgeValue1
					currPath.curState = pcfun.getNextState(currPath.vHat[currPath.currPosition], currPath.curState, self.m)
					currPath.updateBits(idxRev)
					self.trellisPathStack.append(currPath)
					self.T += 1
					self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
				'''
			else:

				currPath.currPosition += 1
				idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
				currPath.updateLLRs(idxRev)
				edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
				currPath.curState = pcfun.getNextState(0, currPath.curState, self.m)
				currPath.uHat[currPath.currPosition] = edgeValue0
				currPath.vHat[currPath.currPosition] = self.polarMask[currPath.currPosition]
				'''
				branchMetric0 = np.log(1 + np.exp(-(1 - 2 * edgeValue0) * currPath.LLRs[0]))
				'''
				'''
				penalty = np.abs(currPath.LLRs[0])
				if currPath.LLRs[0] > 0:
					branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					
				elif currPath.LLRs[0] < 0:
					branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					
				else:
					input('warning')
				'''
				'''
				likelihoodRatio = np.exp(currPath.LLRs[0])
				branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0)))
				'''
				
				'''likelihoodRatio = np.exp(currPath.LLRs[0])
				branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]'''
				
				
				penalty = np.abs(currPath.LLRs[0])
				penalty /= np.log(2)
				if currPath.LLRs[0] > 0:
					branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
				elif currPath.LLRs[0] < 0:
					branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
				else:
					input('warning')

				'''if np.abs(currPath.LLRs[0]) >= 50:
					penalty = np.abs(currPath.LLRs[0])
					penalty /= np.log(2)
					if currPath.LLRs[0] > 0:
						branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
					elif currPath.LLRs[0] < 0:
						branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
				elif currPath.LLRs[0] == 0:
					input('warning')
				else:
					branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]'''
				
				currPath.pathMetric += branchMetric0
				currPath.updateBits(idxRev)
				self.ANV += 1
				self.trellisPathStack.append(currPath)
				self.T += 1
				self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

			if self.trellisPathStack[-1].currPosition == (self.codewordLength - 1):


				best = self.trellisPathStack[-1].vHat
				flag = False
				
		t2 = perf_counter()
		self.runningTime = (t2 - t1) * 1000
		return self.extract(best)

	def pacStackDecoder1(self, softMess):

		self.trellisPathStack = [Path(self.codewordLength, self.m)]
		self.T = 1 
		self.trellisPathStack[0].LLRs[self.codewordLength - 1:] = softMess
		self.trellisPathStack[0].currPosition = -1
		self.criticalMask = np.zeros(self.codewordLength, dtype=int)
		self.criticalMask[self.criticalSet] = 1
		self.visitsCount = np.zeros(self.codewordLength, dtype=int)
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35
		self.ANV = 0
		self.runningTime = 0
		flag = True
		
		t1 = perf_counter()
		while flag:

			currPath = self.trellisPathStack.pop()
			self.T -= 1

			if self.polarMask[currPath.currPosition + 1] == 1:

				if self.criticalMask[currPath.currPosition + 1] == 1:

					self.ANV += 1
					currPath.currPosition += 1
					self.visitsCount[currPath.currPosition] += 2
					self.pathFork1(currPath)

				else:

					self.ANV += 1
					currPath.currPosition += 1
					self.visitsCount[currPath.currPosition] += 1
					idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
					currPath.updateLLRs(idxRev)
					edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
					edgeValue1 = pcfun.conv1Bit(1, currPath.curState, self.gen)
					penalty = np.abs(currPath.LLRs[0])
					penalty /= np.log(2)
					if currPath.LLRs[0] > 0:
						branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
						branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[currPath.currPosition]
					elif currPath.LLRs[0] < 0:
						branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
						branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[currPath.currPosition]
					else:
						input('warning')
					currPath.pathMetric += branchMetric0 if (branchMetric0 > branchMetric1) else branchMetric1
					currPath.vHat[currPath.currPosition] = 0 if (branchMetric0 > branchMetric1) else 1 
					currPath.uHat[currPath.currPosition] = edgeValue0 if (branchMetric0 > branchMetric1) else edgeValue1
					currPath.curState = pcfun.getNextState(currPath.vHat[currPath.currPosition], currPath.curState, self.m)
					currPath.updateBits(idxRev)
					self.trellisPathStack.append(currPath)
					self.T += 1
					self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

			else:

				self.ANV += 1
				currPath.currPosition += 1
				self.visitsCount[currPath.currPosition] += 1
				idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
				currPath.updateLLRs(idxRev)
				edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
				currPath.curState = pcfun.getNextState(0, currPath.curState, self.m)
				currPath.uHat[currPath.currPosition] = edgeValue0
				currPath.vHat[currPath.currPosition] = self.polarMask[currPath.currPosition]
				'''
				branchMetric0 = np.log(1 + np.exp(-(1 - 2 * edgeValue0) * currPath.LLRs[0]))
				'''
				'''
				penalty = np.abs(currPath.LLRs[0])
				if currPath.LLRs[0] > 0:
					branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					
				elif currPath.LLRs[0] < 0:
					branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					
				else:
					input('warning')
				'''
				'''
				likelihoodRatio = np.exp(currPath.LLRs[0])
				branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
				'''
				
				penalty = np.abs(currPath.LLRs[0])
				penalty /= np.log(2)
				if currPath.LLRs[0] > 0:
					branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
				elif currPath.LLRs[0] < 0:
					branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
				else:
					input('warning')
				
				currPath.pathMetric += branchMetric0
				currPath.updateBits(idxRev)
				self.trellisPathStack.append(currPath)
				self.T += 1
				self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

			if self.trellisPathStack[-1].currPosition == (self.codewordLength - 1):

				best = self.trellisPathStack[-1].vHat
				flag = False
				
		t2 = perf_counter()
		self.runningTime = (t2 - t1) * 1000
		return self.extract(best)

	def pathFork1(self, currPath):

		idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
		currPath.updateLLRs(idxRev)
		edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
		edgeValue1 = pcfun.conv1Bit(1, currPath.curState, self.gen)
		pathState0 = pcfun.getNextState(0, currPath.curState, self.m)
		pathState1 = pcfun.getNextState(1, currPath.curState, self.m)
		'''
		branchMetric0 = np.log(1 + np.exp(-(1 - 2 * edgeValue0) * currPath.LLRs[0]))
		branchMetric1 = np.log(1 + np.exp(-(1 - 2 * edgeValue1) * currPath.LLRs[0]))
		'''
		'''
		penalty = np.abs(currPath.LLRs[0])
		if currPath.LLRs[0] > 0:
			branchMetric0 = 0 if (edgeValue0 == 0) else penalty
			branchMetric1 = 0 if (edgeValue1 == 0) else penalty
		elif currPath.LLRs[0] < 0:
			branchMetric0 = 0 if (edgeValue0 == 1) else penalty
			branchMetric1 = 0 if (edgeValue1 == 1) else penalty
		else:
			input('warning')
		'''
		'''
		likelihoodRatio = np.exp(currPath.LLRs[0])
		branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.bias
		branchMetric1 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue1))) - self.bias
		'''
		'''
		likelihoodRatio = np.exp(currPath.LLRs[0])
		branchMetric0 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
		branchMetric1 = 1 - np.log2(1 + likelihoodRatio ** (-(1 - 2 * edgeValue1))) - self.cutoff_rate[currPath.currPosition]
		'''
		
		penalty = np.abs(currPath.LLRs[0])
		penalty /= np.log(2)
		if currPath.LLRs[0] > 0:
			branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
			branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[currPath.currPosition]
		elif currPath.LLRs[0] < 0:
			branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
			branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[currPath.currPosition]
		else:
			input('warning')
		
		copyPath = copy.deepcopy(currPath)
		currPath.pathMetric += branchMetric0
		currPath.uHat[currPath.currPosition] = edgeValue0
		currPath.vHat[currPath.currPosition] = 0
		currPath.curState = pathState0
		currPath.updateBits(idxRev)
		copyPath.pathMetric += branchMetric1
		copyPath.uHat[copyPath.currPosition] = edgeValue1
		copyPath.vHat[copyPath.currPosition] = 1 
		copyPath.curState = pathState1
		copyPath.updateBits(idxRev)
		
		if self.T <= (self.D - 2):
		
			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
		
		else:

			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
			del self.trellisPathStack[0]
			self.T -= 1
		
	'''if np.abs(currPath.LLRs[0]) >= 50:

			penalty = np.abs(currPath.LLRs[0])
			penalty /= np.log(2)
			if currPath.LLRs[0] > 0:
				branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
				branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[currPath.currPosition]
			elif currPath.LLRs[0] < 0:
				branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
				branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[currPath.currPosition]
		elif currPath.LLRs[0] == 0:
			input('warning')
		else:

			branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
			branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[currPath.currPosition]'''

	def pacStackDecoder2(self, softMess):

		self.trellisPathStack = [Path(self.codewordLength, self.m)]
		self.T = 1 
		self.trellisPathStack[0].LLRs[self.codewordLength - 1:] = softMess
		self.trellisPathStack[0].currPosition = -1
		self.visitsCount = np.zeros(self.codewordLength, dtype=int)
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35
		self.ANV = 0
		self.runningTime = 0
		self.numPaths = 0
		flag = True

		#t1 = perf_counter()
		while flag:

			currPath = self.trellisPathStack.pop()
			self.T -= 1

			if self.polarMask[currPath.currPosition + 1] == 1:

				currPath.currPosition += 1
				#self.visitsCount[currPath.currPosition] += 2
				self.pathFork2(currPath)
				self.ANV += 1

			else:

				currPath.currPosition += 1
				#self.visitsCount[currPath.currPosition] += 1
				idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
				currPath.updateLLRs(idxRev)
				edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
				currPath.curState = pcfun.getNextState(0, currPath.curState, self.m)
				currPath.uHat[currPath.currPosition] = edgeValue0
				currPath.vHat[currPath.currPosition] = self.polarMask[currPath.currPosition]

				'''penalty = np.abs(currPath.LLRs[0])
				penalty /= np.log(2)
				if currPath.LLRs[0] > 0:
					#branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
				elif currPath.LLRs[0] < 0:
					#branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
				else:
					input('warning')'''
				branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
				#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[currPath.currPosition]
				#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]

				currPath.pathMetric += branchMetric0
				currPath.updateBits(idxRev)
				self.ANV += 1
				self.trellisPathStack.append(currPath)
				self.T += 1
				self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

			if self.trellisPathStack[-1].currPosition == (self.codewordLength - 1):

				self.numPaths = len(self.trellisPathStack)
				best = self.trellisPathStack[-1].vHat
				flag = False
				
		#t2 = perf_counter()
		#self.runningTime = (t2 - t1) * 1000
		return self.extract(best)

	def pathFork2(self, currPath):

		idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
		currPath.updateLLRs(idxRev)
		edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
		edgeValue1 = pcfun.conv1Bit(1, currPath.curState, self.gen)
		pathState0 = pcfun.getNextState(0, currPath.curState, self.m)
		pathState1 = pcfun.getNextState(1, currPath.curState, self.m)
		
		'''penalty = np.abs(currPath.LLRs[0])
		penalty /= np.log(2)
		if currPath.LLRs[0] > 0:
			branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[currPath.currPosition]
			branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[currPath.currPosition]
		elif currPath.LLRs[0] < 0:
			branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[currPath.currPosition]
			branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[currPath.currPosition]
		else:
			input('warning')'''
		branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
		branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[currPath.currPosition]
		#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[currPath.currPosition]
		#branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.I[currPath.currPosition]
		#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
		#branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.cutoff_rate[currPath.currPosition]

		copyPath = copy.deepcopy(currPath)
		currPath.pathMetric += branchMetric0
		currPath.uHat[currPath.currPosition] = edgeValue0
		currPath.vHat[currPath.currPosition] = 0
		currPath.curState = pathState0
		currPath.updateBits(idxRev)
		copyPath.pathMetric += branchMetric1
		copyPath.uHat[copyPath.currPosition] = edgeValue1
		copyPath.vHat[copyPath.currPosition] = 1 
		copyPath.curState = pathState1
		copyPath.updateBits(idxRev)
		
		
		'''if self.T <= (self.D - 2):
		
			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
		
		else:

			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
			del self.trellisPathStack[0]
			self.T -= 1
			if len(self.trellisPathStack) != self.D:
				input('warning')'''
		self.trellisPathStack.append(currPath)
		self.trellisPathStack.append(copyPath)
		self.T += 2
		self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
	
	def PSCS(self, softMess):

		self.trellisPathStack = [Path(self.codewordLength, self.m)]
		self.T = 1 
		self.trellisPathStack[0].LLRs[self.codewordLength - 1:] = softMess
		self.trellisPathStack[0].currPosition = -1
		self.visitsCount = np.zeros(self.codewordLength, dtype=int)
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35
		self.ANV = 0
		self.runningTime = 0
		self.numPaths = 0
		self.mT = -11
		flag = True

		#t1 = perf_counter()
		while flag:

			currPath = self.trellisPathStack.pop()
			self.T -= 1

			if self.polarMask[currPath.currPosition + 1] == 1:

				currPath.currPosition += 1
				#self.visitsCount[currPath.currPosition] += 2
				self.PSCSFork(currPath)
				self.ANV += 1

			else:

				currPath.currPosition += 1
				#self.visitsCount[currPath.currPosition] += 1
				idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
				currPath.updateLLRs(idxRev)
				edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
				currPath.curState = pcfun.getNextState(0, currPath.curState, self.m)
				currPath.uHat[currPath.currPosition] = edgeValue0
				currPath.vHat[currPath.currPosition] = self.polarMask[currPath.currPosition]

				#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
				branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[currPath.currPosition]
				#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
				currPath.pathMetric += branchMetric0
				currPath.updateBits(idxRev)
				self.ANV += 1
				self.trellisPathStack.append(currPath)
				self.T += 1
				self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

			if self.trellisPathStack[-1].currPosition == (self.codewordLength - 1):

				self.numPaths = len(self.trellisPathStack)
				best = self.trellisPathStack[-1].vHat
				flag = False
				
		#t2 = perf_counter()
		#self.runningTime = (t2 - t1) * 1000
		return self.extract(best)

	def PSCSFork(self, currPath):

		idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
		currPath.updateLLRs(idxRev)
		edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
		edgeValue1 = pcfun.conv1Bit(1, currPath.curState, self.gen)
		pathState0 = pcfun.getNextState(0, currPath.curState, self.m)
		pathState1 = pcfun.getNextState(1, currPath.curState, self.m)
		retainBoth = False
		
		#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
		#branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[currPath.currPosition]
		branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[currPath.currPosition]
		branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.I[currPath.currPosition]
		#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
		#branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.cutoff_rate[currPath.currPosition]

		#copyPath = copy.deepcopy(currPath)
		if branchMetric0 >= self.mT and branchMetric1 >= self.mT:
			copyPath = copy.deepcopy(currPath)
			currPath.pathMetric += branchMetric0
			currPath.uHat[currPath.currPosition] = edgeValue0
			currPath.vHat[currPath.currPosition] = 0
			currPath.curState = pathState0
			currPath.updateBits(idxRev)
			copyPath.pathMetric += branchMetric1
			copyPath.uHat[copyPath.currPosition] = edgeValue1
			copyPath.vHat[copyPath.currPosition] = 1 
			copyPath.curState = pathState1
			copyPath.updateBits(idxRev)
			retainBoth = True
		elif branchMetric0 >= self.mT:
			currPath.pathMetric += branchMetric0
			currPath.uHat[currPath.currPosition] = edgeValue0
			currPath.vHat[currPath.currPosition] = 0
			currPath.curState = pathState0
			currPath.updateBits(idxRev)
		elif branchMetric1 >= self.mT:
			currPath.pathMetric += branchMetric1
			currPath.uHat[currPath.currPosition] = edgeValue1
			currPath.vHat[currPath.currPosition] = 1
			currPath.curState = pathState1
			currPath.updateBits(idxRev)
		
		if retainBoth:
			self.trellisPathStack.append(currPath)
			self.trellisPathStack.append(copyPath)
			self.T += 2
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
		else:
			self.trellisPathStack.append(currPath)
			self.T += 1
			self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

	def GCAPSCS(self, softMess):

		self.trellisPathStack = [Path(self.codewordLength, self.m)]
		self.T = 1 
		self.trellisPathStack[0].LLRs[self.codewordLength - 1:] = softMess
		self.trellisPathStack[0].currPosition = -1
		self.visitsCount = np.zeros(self.codewordLength, dtype=int)
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35
		self.ANV = 0
		self.runningTime = 0
		self.numPaths = 0
		self.mT = -11
		flag = True
		idxRecord = []
		for i in range(self.codewordLength):
			if self.polarMask[i] == 1 and self.GA[i] < 18:
				idxRecord.append(i)
		#t1 = perf_counter()

		while flag:

			currPath = self.trellisPathStack.pop()
			self.T -= 1

			if self.polarMask[currPath.currPosition + 1] == 1:
				if (currPath.currPosition + 1) in idxRecord:

					currPath.currPosition += 1
					#self.visitsCount[currPath.currPosition] += 2
					self.PSCSFork(currPath)
					self.ANV += 1
				else:
					currPath.currPosition += 1
					#self.visitsCount[currPath.currPosition] += 1
					idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
					currPath.updateLLRs(idxRev)
					edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
					edgeValue1 = pcfun.conv1Bit(1, currPath.curState, self.gen)
					#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
					#branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[currPath.currPosition]
					branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[currPath.currPosition]
					branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.I[currPath.currPosition]
					#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
					#branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.cutoff_rate[currPath.currPosition]
					if branchMetric0 > branchMetric1:
						currPath.curState = pcfun.getNextState(0, currPath.curState, self.m)
						currPath.uHat[currPath.currPosition] = edgeValue0
						currPath.vHat[currPath.currPosition] = 0
						currPath.pathMetric += branchMetric0
					else:
						currPath.curState = pcfun.getNextState(1, currPath.curState, self.m)
						currPath.uHat[currPath.currPosition] = edgeValue1
						currPath.vHat[currPath.currPosition] = 1
						currPath.pathMetric += branchMetric1

					currPath.updateBits(idxRev)
					self.ANV += 1
					self.trellisPathStack.append(currPath)
					self.T += 1
					self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
			else:

				currPath.currPosition += 1
				#self.visitsCount[currPath.currPosition] += 1
				idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
				currPath.updateLLRs(idxRev)
				edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
				currPath.curState = pcfun.getNextState(0, currPath.curState, self.m)
				currPath.uHat[currPath.currPosition] = edgeValue0
				currPath.vHat[currPath.currPosition] = self.polarMask[currPath.currPosition]

				#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
				branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[currPath.currPosition]
				#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
				currPath.pathMetric += branchMetric0
				currPath.updateBits(idxRev)
				self.ANV += 1
				self.trellisPathStack.append(currPath)
				self.T += 1
				self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

			if self.trellisPathStack[-1].currPosition == (self.codewordLength - 1):

				self.numPaths = len(self.trellisPathStack)
				best = self.trellisPathStack[-1].vHat
				flag = False
				
		#t2 = perf_counter()
		#self.runningTime = (t2 - t1) * 1000
		return self.extract(best)

	def CAPSCS(self, softMess):

		self.trellisPathStack = [Path(self.codewordLength, self.m)]
		self.T = 1 
		self.trellisPathStack[0].LLRs[self.codewordLength - 1:] = softMess
		self.trellisPathStack[0].currPosition = -1
		self.visitsCount = np.zeros(self.codewordLength, dtype=int)
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35
		self.ANV = 0
		self.runningTime = 0
		self.numPaths = 0
		self.mT = -11
		flag = True

		#t1 = perf_counter()
		while flag:

			currPath = self.trellisPathStack.pop()
			self.T -= 1

			if self.polarMask[currPath.currPosition + 1] == 1:
				if (currPath.currPosition + 1) in self.criticalSet[:-1]:

					currPath.currPosition += 1
					#self.visitsCount[currPath.currPosition] += 2
					self.PSCSFork(currPath)
					self.ANV += 1
				else:
					currPath.currPosition += 1
					#self.visitsCount[currPath.currPosition] += 1
					idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
					currPath.updateLLRs(idxRev)
					edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
					edgeValue1 = pcfun.conv1Bit(1, currPath.curState, self.gen)
					branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
					branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[currPath.currPosition]
					#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[currPath.currPosition]
					#branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.I[currPath.currPosition]
					#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
					#branchMetric1 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.cutoff_rate[currPath.currPosition]
					if branchMetric0 > branchMetric1:
						currPath.curState = pcfun.getNextState(0, currPath.curState, self.m)
						currPath.uHat[currPath.currPosition] = edgeValue0
						currPath.vHat[currPath.currPosition] = 0
						currPath.pathMetric += branchMetric0
					else:
						currPath.curState = pcfun.getNextState(1, currPath.curState, self.m)
						currPath.uHat[currPath.currPosition] = edgeValue1
						currPath.vHat[currPath.currPosition] = 1
						currPath.pathMetric += branchMetric1

					currPath.updateBits(idxRev)
					self.ANV += 1
					self.trellisPathStack.append(currPath)
					self.T += 1
					self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)
			else:

				currPath.currPosition += 1
				#self.visitsCount[currPath.currPosition] += 1
				idxRev = pcfun.bitreversed(currPath.currPosition, self.n)
				currPath.updateLLRs(idxRev)
				edgeValue0 = pcfun.conv1Bit(0, currPath.curState, self.gen)
				currPath.curState = pcfun.getNextState(0, currPath.curState, self.m)
				currPath.uHat[currPath.currPosition] = edgeValue0
				currPath.vHat[currPath.currPosition] = self.polarMask[currPath.currPosition]

				branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[currPath.currPosition]
				#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[currPath.currPosition]
				#branchMetric0 = 1 - np.log2(1 + np.exp(currPath.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[currPath.currPosition]
				currPath.pathMetric += branchMetric0
				currPath.updateBits(idxRev)
				self.ANV += 1
				self.trellisPathStack.append(currPath)
				self.T += 1
				self.trellisPathStack.sort(key=lambda trellis_path_stack: trellis_path_stack.pathMetric, reverse=False)

			if self.trellisPathStack[-1].currPosition == (self.codewordLength - 1):

				self.numPaths = len(self.trellisPathStack)
				best = self.trellisPathStack[-1].vHat
				flag = False
				
		#t2 = perf_counter()
		#self.runningTime = (t2 - t1) * 1000
		return self.extract(best)

	def SCFano(self, softMess):

		flag = True
		path = Path(self.codewordLength)
		path.LLRs[self.codewordLength - 1:] = softMess
		path.pm = np.zeros(self.codewordLength, dtype=float)
		beta = np.zeros(self.infoLen, dtype=float)
	 	#betaCut = np.zeros(self.infoLen, dtype=float)
		gama = np.zeros(self.infoLen, dtype=int)
		B = 0
		C = 0
		j = -1
		i = -1
		iPrevious = -2
		T = self.threshold
		self.iterations = 0
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35
	    # self.threshold = 0
	    # self.ANV = 0
	    #self.pe = pcfun.PEDega(self.codewordLength, self.infoLen, self.designSNR)
		'''for k in range(self.codewordLength):
			C += np.log(1 - self.pe[k])'''

		while flag:
		
			'''if i != iPrevious:
				idxRev = pcfun.bitreversed(i + 1, self.n)
				path.updateLLRs(idxRev)
				iPrevious = i'''

			if self.polarMask[i + 1] == 1:

				idxRev = pcfun.bitreversed(i + 1, self.n)
				path.updateLLRs(idxRev)
				branchMetric0 = np.log(((np.exp(path.LLRs[0]) / (np.exp(path.LLRs[0]) + 1))) / (1 - self.pe[i + 1]))
				branchMetric1 = np.log(((1 / (np.exp(path.LLRs[0]) + 1))) / (1 - self.pe[i + 1]))
				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * 0))) - self.bias[i + 1]
				#branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * 1))) - self.bias[i + 1]
				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * 0))) - self.cutoff_rate[i + 1]
				#branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * 1))) - self.cutoff_rate[i + 1]
				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * 0))) - self.I[i + 1]
				#branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * 1))) - self.I[i + 1]

				if branchMetric0 == branchMetric1:
					input('warning')
				'''pathMetric0 = path.pm[i + 1 - 1] + branchMetric0
				pathMetric1 = path.pm[i + 1 - 1] + branchMetric1'''
				pathMetric0 = (0 + branchMetric0) if ((i + 1) == 0) else (path.pm[i + 1 - 1] + branchMetric0)
				pathMetric1 = (0 + branchMetric1) if ((i + 1) == 0) else (path.pm[i + 1 - 1] + branchMetric1)

				if max(pathMetric0, pathMetric1) >= T:
					if B == 0:

						i += 1
						j += 1
						path.pm[i] = pathMetric0 if (pathMetric0 > pathMetric1) else pathMetric1
						path.uHat[i] = 0 if (pathMetric0 > pathMetric1) else 1
						beta[j] = path.pm[i]
						# betaCut[j] = min(pathMetric0, pathMetric1)
						gama[j] = 0
						# muPre = 0 if ((i - 1) == -1) else path.pm[i - 1]
						muPre = 0 if j == 0 else beta[j - 1]
						if muPre < T + self.delta:
							while T + self.delta <= path.pm[i]:
								T += self.delta

						path.updateBits(idxRev)

					else:
						if min(pathMetric0, pathMetric1) > T:

							i += 1
							j += 1
							path.pm[i] = pathMetric0 if (pathMetric0 < pathMetric1) else pathMetric1
							path.uHat[i] = 0 if (pathMetric0 < pathMetric1) else 1
							beta[j] = path.pm[i]
							gama[j] = 1
							'''muPre = 0 if j == 0 else beta[j - 1]
							if muPre < T + self.delta:
								while T + self.delta <= path.pm[i]:
									T += self.delta'''
							path.updateBits(idxRev)
							B = 0

						else:
							(T, j, B) = self.moveBack(beta, j, T, gama)
							i = self.A[0] - 1 if j == -1 else self.A[j]

							for ii in range(0, i + 1):
								idxRev = pcfun.bitreversed(ii, self.n)
								path.updateLLRs(idxRev)
								path.updateBits(idxRev)

				else:
					(T, j, B) = self.moveBack(beta, j, T, gama)
					i = self.A[0] - 1 if j == -1 else self.A[j]

					for ii in range(0, i + 1):
						idxRev = pcfun.bitreversed(ii, self.n)
						path.updateLLRs(idxRev)
						path.updateBits(idxRev)
	
			else:

				i += 1
				idxRev = pcfun.bitreversed(i, self.n)
				path.updateLLRs(idxRev)
				branchMetric0 = np.log(((np.exp(path.LLRs[0]) / (np.exp(path.LLRs[0]) + 1))) / (1 - self.pe[i]))
				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * 0))) - self.bias[i]
				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * 0))) - self.cutoff_rate[i]
				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * 0))) - self.I[i]
				path.pm[i] = 0 + branchMetric0 if i == 0 else path.pm[i - 1] + branchMetric0
				#path.pm[i] = C + branchMetric0 if i == 0 else path.pm[i - 1] + branchMetric0
				#path.pm[i] = 0 if i == 0 else path.pm[i - 1]

				path.uHat[i] = self.polarMask[i]
				path.updateBits(idxRev)

			if i == (self.codewordLength - 1):
				self.gama = gama
				self.uHat = path.uHat
				flag = False
		return self.extract(path.uHat)

	def PACFano(self, softMess):

		flag = True
		path = Path(self.codewordLength, self.m)
		path.LLRs[self.codewordLength - 1:] = softMess
		path.pm = np.zeros(self.codewordLength, dtype=float)
		beta = np.zeros(self.infoLen, dtype=float)
		# betaCut = np.zeros(self.infoLen, dtype=float)
		gama = np.zeros(self.infoLen, dtype=int)
		stateStorage = [[] for i in range(self.infoLen)]
		B = 0
		#C = 0
		j = -1
		i = -1
		#iPre = -2
		self.iterations = 0
		T = self.threshold
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35
		'''self.pe = pcfun.PEDega(self.codewordLength, self.infoLen, self.designSNR)
		for k in range(self.codewordLength):
			C += np.log(1 - self.pe[k])'''

		while flag:

			'''print(i)
			time.sleep(0.1)'''
			if self.polarMask[i + 1] == 1:

				idxRev = pcfun.bitreversed(i + 1, self.n)
				path.updateLLRs(idxRev)
				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				edgeValue1 = pcfun.conv1Bit(1, path.curState, self.gen)
				nextState0 = pcfun.getNextState(0, path.curState, self.m)
				nextState1 = pcfun.getNextState(1, path.curState, self.m)
				'''penalty = np.abs(path.LLRs[0])
				penalty /= np.log(2)

				if path.LLRs[0] > 0:
					branchMetric0 = 1 - (0 if (edgeValue0 == 0) else penalty) - self.bias[i + 1]
					branchMetric1 = 1 - (0 if (edgeValue1 == 0) else penalty) - self.bias[i + 1]
				elif path.LLRs[0] < 0:
					branchMetric0 = 1 - (0 if (edgeValue0 == 1) else penalty) - self.bias[i + 1]
					branchMetric1 = 1 - (0 if (edgeValue1 == 1) else penalty) - self.bias[i + 1]
				else:
					input('Warning')'''

				'''branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[
					i + 1]
				branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.cutoff_rate[
					i + 1]'''
				branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[
					i + 1]
				branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[
					i + 1]
				'''branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[
					i + 1]
				branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.I[
					i + 1]'''
				'''branchMetric0 = np.log(np.exp(path.LLRs[0]) / (np.exp(path.LLRs[0]) + 1)) - np.log(1.0 - self.pe[i + 1])
				branchMetric1 = np.log(1 / (np.exp(path.LLRs[0]) + 1)) - np.log(1.0 - self.pe[i + 1])'''
				# pathMetric0 = path.pm[path.currPosition - 1] + branchMetric0
				# pathMetric1 = path.pm[path.currPosition - 1] + branchMetric1
				pathMetric0 = (0 + branchMetric0) if ((i + 1) == 0) else (path.pm[i + 1 - 1] + branchMetric0)
				pathMetric1 = (0 + branchMetric1) if ((i + 1) == 0) else (path.pm[i + 1 - 1] + branchMetric1)

				if max(pathMetric0, pathMetric1) >= T:

					if B == 0:

						self.iterations += 1
						i += 1
						j += 1
						path.pm[i] = pathMetric0 if (pathMetric0 > pathMetric1) else pathMetric1
						path.vHat[i] = 0 if (pathMetric0 > pathMetric1) else 1
						path.uHat[i] = edgeValue0 if (pathMetric0 > pathMetric1) else edgeValue1
						beta[j] = path.pm[i]
						# betaCut[j] = min(pathMetric0, pathMetric1)
						gama[j] = 0

						path.updateBits(idxRev)

						# stateStorage[j] = path.curState
						path.curState = nextState0 if (pathMetric0 > pathMetric1) else nextState1
						stateStorage[j] = path.curState #deep copy?

						muPre = path.pm[i - 1] if (j == 0) else beta[j - 1]
						#muPre = 0 if (j == 0) else beta[j - 1]
						if muPre < (T + self.delta):
							while (T + self.delta) <= path.pm[i]:
								T += self.delta
					else:
						if min(pathMetric0, pathMetric1) > T:

							i += 1
							j += 1
							path.pm[i] = pathMetric0 if (pathMetric0 < pathMetric1) else pathMetric1
							path.vHat[i] = 0 if (pathMetric0 < pathMetric1) else 1
							path.uHat[i] = edgeValue0 if (pathMetric0 < pathMetric1) else edgeValue1
							beta[j] = path.pm[i]
							# betaCut[j] = max(pathMetric0, pathMetric1)
							gama[j] = 1
							# stateStorage[j] = path.curState
							path.curState = nextState0 if (pathMetric0 < pathMetric1) else nextState1
							stateStorage[j] = path.curState #deep copy?
							path.updateBits(idxRev)

							#miu = -math.inf if path.currPosition == 0 else path.pm[path.currPosition - 1]
							muPre = path.pm[i - 1] if (j == 0) else beta[j - 1]
							#muPre = 0 if (j == 0) else beta[j - 1]
							if muPre < (T + self.delta):
								while (T + self.delta) <= path.pm[i]:
									T = T + self.delta
							B = 0
						else:
							jj = j
							# stateStorage[j] = path.curState
							(T, j, B) = self.moveBack(beta, j, T, gama)
							i = self.A[0] - 1 if j == -1 else self.A[j]
							path.curState = [0 for mm in range(self.m)] if j == -1 else stateStorage[
								j] if jj != j else path.curState
							for ii in range(0, i + 1):
								idxRev = pcfun.bitreversed(ii, self.n)
								path.updateLLRs(idxRev)
								path.updateBits(idxRev)
				else:
					jj = j
					#stateStorage[j] = path.curState
					(T, j, B) = self.moveBack(beta, j, T, gama)
					i = self.A[0] - 1 if j == -1 else self.A[j]
					path.curState = [0 for mm in range(self.m)] if j == -1 else stateStorage[j] if jj != j else path.curState
					for ii in range(0, i + 1):
						idxRev = pcfun.bitreversed(ii, self.n)
						path.updateLLRs(idxRev)
						path.updateBits(idxRev)
			else:

				i += 1
				idxRev = pcfun.bitreversed(i, self.n)
				path.updateLLRs(idxRev)

				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				nextState0 = pcfun.getNextState(0, path.curState, self.m)
				'''penalty = np.abs(path.LLRs[0])
				penalty /= np.log(2)

				if path.LLRs[0] > 0:
					branchMetric0 = 1 - (0 if edgeValue0 == 0 else penalty) - self.bias[i]
				elif path.LLRs[0] < 0:
					branchMetric0 = 1 - (0 if edgeValue0 == 1 else penalty) - self.bias[i]
				else:
					input('warning')'''

				'''branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.cutoff_rate[
					i] #????'''
				branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[
					i] #?????? why i + 1?
				'''branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[
					i]'''
				#branchMetric0 = np.log(np.exp(path.LLRs[0]) / (np.exp(path.LLRs[0]) + 1)) - np.log(1.0 - self.pe[i])
				#path.pm[i] = C + branchMetric0 if i == 0 else path.pm[i - 1] + branchMetric0
				path.pm[i] = 0 + branchMetric0 if i == 0 else path.pm[i - 1] + branchMetric0
				# path.pm[path.currPosition] = path.pm[path.currPosition - 1] + branchMetric0

				path.vHat[i] = self.polarMask[i]
				path.uHat[i] = edgeValue0
				path.curState = nextState0 #deep copy?
				path.updateBits(idxRev)


			if i == (self.codewordLength - 1):
				flag = False

		return self.extract(path.vHat)

	def moveBack222(self, i, T, pm, indicator):
		while True:
			muPre = -np.inf if i == -1 else 0 if i == 0 else pm[i - 1]
			if muPre < T:
				T -= self.delta
				indicator[i + 1] = 0
				return (i, T, indicator)
			else:
				#indicator[i] += 1
				if (indicator[i] + 1) == 2 or (i not in self.criticalSet):
					i = i - 1
					continue
				else:
					indicator[i] += 1
					i = i - 1
					return (i, T, indicator)

	def fano222(self, softMess):
		flag = True
		N = self.codewordLength
		n = int(np.log2(N))
		#m = self.m
		path = Path(N, self.m)
		path.LLRs[N - 1:] = softMess
		pm = np.zeros(N, dtype=float)
		state = [[] for i in range(N)]
		indicator = np.zeros(N, dtype=int)
		i = -1
		pointer = -1
		self.iterations = 0
		T = self.threshold
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35

		while flag:
			'''print(i)
			time.sleep(0.1)'''
			pointer = i + 1 
			if self.polarMask[pointer] == 1:
				if pointer in self.criticalSet[:-1]:
					path.updateLLRs(pcfun.bitreversed(pointer, n))
					edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
					nextState0 = pcfun.getNextState(0, path.curState, self.m)
					edgeValue1 = pcfun.conv1Bit(1, path.curState, self.gen)
					nextState1 = pcfun.getNextState(1, path.curState, self.m)

					#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[pointer]
					#branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[pointer]
					branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[pointer]
					branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.I[pointer]
					pathMetric0 = 0 + branchMetric0 if pointer == 0 else pm[pointer - 1] + branchMetric0
					pathMetric1 = 0 + branchMetric1 if pointer == 0 else pm[pointer - 1] + branchMetric1
					metricMax = pathMetric0 if (pathMetric0 > pathMetric1) else pathMetric1
					metricMin = pathMetric0 if (pathMetric0 < pathMetric1) else pathMetric1
					if indicator[pointer] == 0:
						pm[pointer] = metricMax
					else:
						pm[pointer] = metricMin
				else:
					path.updateLLRs(pcfun.bitreversed(pointer, n))
					edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
					nextState0 = pcfun.getNextState(0, path.curState, self.m)
					edgeValue1 = pcfun.conv1Bit(1, path.curState, self.gen)
					nextState1 = pcfun.getNextState(1, path.curState, self.m)

					#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[pointer]
					#branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[pointer]
					branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[pointer]
					branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.I[pointer]
					pathMetric0 = 0 + branchMetric0 if pointer == 0 else pm[pointer - 1] + branchMetric0
					pathMetric1 = 0 + branchMetric1 if pointer == 0 else pm[pointer - 1] + branchMetric1

					pm[pointer] = pathMetric0 if (pathMetric0 > pathMetric1) else pathMetric1
			else:

				path.updateLLRs(pcfun.bitreversed(pointer, n))
				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				nextState0 = pcfun.getNextState(0, path.curState, self.m)

				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[pointer]
				branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[pointer]
				pathMetric0 = 0 + branchMetric0 if pointer == 0 else pm[pointer - 1] + branchMetric0
				#pathMetric1 = inf
				pm[pointer] = pathMetric0

			if pointer >= self.criticalSet[-1] or pm[pointer] >= T:
				self.iterations += 1
				muPre = -np.inf if pointer == -1 else 0 if pointer == 0 else pm[pointer - 1]
				if muPre < (T + self.delta):
					while (T + self.delta) <= pm[pointer]:
						T = T + self.delta
				i += 1
				path.vHat[i] = 0 if pm[i] == pathMetric0 else 1 
				path.uHat[i] = edgeValue0 if pm[i] == pathMetric0 else edgeValue1
				path.curState = copy.deepcopy(nextState0) if pm[i] == pathMetric0 else copy.deepcopy(nextState1)
				state[i] = copy.deepcopy(path.curState)
				path.updateBits(pcfun.bitreversed(i, n))

				if (i + 1) == N:
					return self.extract(path.vHat)
				indicator[i + 1] = 0
			else:
				j = i
				(i, T, indicator) = self.moveBack222(i, T, pm, indicator)
				if i == -1:
					path.curState = copy.deepcopy([0 for mm in range(self.m)])
				elif i != j:
					path.curState = copy.deepcopy(state[i])
				if j != i:
					for ii in range(0, i + 1):
						idxRev = pcfun.bitreversed(ii, n)
						path.updateLLRs(idxRev)
						path.updateBits(idxRev)

	def moveBack2(self, i, T, pm, indicator):
		while True:
			muPre = -np.inf if i == -1 else 0 if i == 0 else pm[i - 1]
			if muPre < T:
				T -= self.delta
				indicator[i + 1] = 0
				return (i, T, indicator)
			else:
				#indicator[i] += 1
				if (indicator[i] + 1) == 2 or (self.polarMask[i] == 0):
					i = i - 1
					continue
				else:
					indicator[i] += 1
					i = i - 1
					return (i, T, indicator)

	def fano(self, softMess):
		flag = True
		N = self.codewordLength
		n = int(np.log2(N))
		#m = self.m
		path = Path(N, self.m)
		path.LLRs[N - 1:] = softMess
		pm = np.zeros(N, dtype=float)
		state = [[] for i in range(N)]
		indicator = np.zeros(N, dtype=int)
		i = -1
		pointer = -1
		self.iterations = 0
		T = self.threshold
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35

		while flag:
			'''print(i)
			time.sleep(0.1)'''
			pointer = i + 1 
			if self.polarMask[pointer] == 1:

				path.updateLLRs(pcfun.bitreversed(pointer, n))
				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				nextState0 = pcfun.getNextState(0, path.curState, self.m)
				edgeValue1 = pcfun.conv1Bit(1, path.curState, self.gen)
				nextState1 = pcfun.getNextState(1, path.curState, self.m)

				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[pointer]
				#branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[pointer]
				branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[pointer]
				branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.I[pointer]
				pathMetric0 = 0 + branchMetric0 if pointer == 0 else pm[pointer - 1] + branchMetric0
				pathMetric1 = 0 + branchMetric1 if pointer == 0 else pm[pointer - 1] + branchMetric1
				metricMax = pathMetric0 if (pathMetric0 > pathMetric1) else pathMetric1
				metricMin = pathMetric0 if (pathMetric0 < pathMetric1) else pathMetric1
				if indicator[pointer] == 0:
					pm[pointer] = metricMax
				else:
					pm[pointer] = metricMin
			else:

				path.updateLLRs(pcfun.bitreversed(pointer, n))
				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				nextState0 = pcfun.getNextState(0, path.curState, self.m)

				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[pointer]
				branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[pointer]
				pathMetric0 = 0 + branchMetric0 if pointer == 0 else pm[pointer - 1] + branchMetric0
				#pathMetric1 = inf
				pm[pointer] = pathMetric0

			if pm[pointer] >= T:
				self.iterations += 1
				muPre = -np.inf if pointer == -1 else 0 if pointer == 0 else pm[pointer - 1]
				if muPre < (T + self.delta):
					while (T + self.delta) <= pm[pointer]:
						T = T + self.delta
				i += 1
				path.vHat[i] = 0 if pm[i] == pathMetric0 else 1 
				path.uHat[i] = edgeValue0 if pm[i] == pathMetric0 else edgeValue1
				path.curState = copy.deepcopy(nextState0) if pm[i] == pathMetric0 else copy.deepcopy(nextState1)
				state[i] = copy.deepcopy(path.curState)
				path.updateBits(pcfun.bitreversed(i, n))

				if (i + 1) == N:
					return self.extract(path.vHat)
				indicator[i + 1] = 0
			else:
				j = i
				(i, T, indicator) = self.moveBack2(i, T, pm, indicator)
				if i == -1:
					path.curState = copy.deepcopy([0 for mm in range(self.m)])
				elif i != j:
					path.curState = copy.deepcopy(state[i])
				if j != i:
					for ii in range(0, i + 1):
						idxRev = pcfun.bitreversed(ii, n)
						path.updateLLRs(idxRev)
						path.updateBits(idxRev)

	def moveBack3(self, i, T, pm, indicator):
		while True:
			muPre = -np.inf if i == -1 else 0 if i == 0 else pm[i - 1]
			if muPre < T:
				T -= self.delta
				indicator[i + 1] = 0
				return (i, T, indicator)
			else:
				#indicator[i] += 1
				if (indicator[i] + 1) == 2 or (self.polarMask[i] == 0): #or (self.branchMetrics[i] < self.mT):
					i = i - 1
					continue
				else:
					indicator[i] += 1
					i = i - 1
					return (i, T, indicator)

	def fano1(self, softMess):
		flag = True
		N = self.codewordLength
		n = int(np.log2(N))
		#m = self.m
		path = Path(N, self.m)
		path.LLRs[N - 1:] = softMess
		pm = np.zeros(N, dtype=float)
		state = [[] for i in range(N)]
		indicator = np.zeros(N, dtype=int)
		i = -1
		pointer = -1
		self.iterations = 0
		T = self.threshold
		self.branchMetrics = np.zeros(self.codewordLength, dtype=float)
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35
		self.mT = -14

		while flag:

			pointer = i + 1 
			if self.polarMask[pointer] == 1:

				path.updateLLRs(pcfun.bitreversed(pointer, n))
				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				nextState0 = pcfun.getNextState(0, path.curState, self.m)
				edgeValue1 = pcfun.conv1Bit(1, path.curState, self.gen)
				nextState1 = pcfun.getNextState(1, path.curState, self.m)

				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[pointer]
				#branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[pointer]
				branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[pointer]
				branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.I[pointer]
				pathMetric0 = 0 + branchMetric0 if pointer == 0 else pm[pointer - 1] + branchMetric0
				pathMetric1 = 0 + branchMetric1 if pointer == 0 else pm[pointer - 1] + branchMetric1
				metricMax = pathMetric0 if (pathMetric0 > pathMetric1) else pathMetric1
				metricMin = pathMetric0 if (pathMetric0 < pathMetric1) else pathMetric1
				if indicator[pointer] == 0:
					pm[pointer] = metricMax
				else:
					pm[pointer] = metricMin
			else:

				path.updateLLRs(pcfun.bitreversed(pointer, n))
				edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
				nextState0 = pcfun.getNextState(0, path.curState, self.m)

				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[pointer]
				branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[pointer]
				branchMetric1 = -np.inf
				pathMetric0 = 0 + branchMetric0 if pointer == 0 else pm[pointer - 1] + branchMetric0
				pm[pointer] = pathMetric0

			if pm[pointer] >= T:
				self.iterations += 1
				muPre = -np.inf if pointer == -1 else 0 if pointer == 0 else pm[pointer - 1]
				if muPre < (T + self.delta):
					while (T + self.delta) <= pm[pointer]:
						T = T + self.delta
				i += 1
				self.branchMetrics[i] = branchMetric1 if pm[i] == pathMetric0 else branchMetric0
				path.vHat[i] = 0 if pm[i] == pathMetric0 else 1 
				path.uHat[i] = edgeValue0 if pm[i] == pathMetric0 else edgeValue1
				path.curState = copy.deepcopy(nextState0) if pm[i] == pathMetric0 else copy.deepcopy(nextState1)
				state[i] = copy.deepcopy(path.curState)
				path.updateBits(pcfun.bitreversed(i, n))

				if (i + 1) == N:
					return self.extract(path.vHat)
				indicator[i + 1] = 0
			else:
				j = i
				(i, T, indicator) = self.moveBack3(i, T, pm, indicator)
				if i == -1:
					path.curState = copy.deepcopy([0 for mm in range(self.m)])
				elif i != j:
					path.curState = copy.deepcopy(state[i])
				if j != i:
					for ii in range(0, i + 1):
						idxRev = pcfun.bitreversed(ii, n)
						path.updateLLRs(idxRev)
						path.updateBits(idxRev)

	def fanoPolar(self, softMess):
		flag = True
		N = self.codewordLength
		K = self.infoLen
		n = int(np.log2(N))
		path = Path(N, self.m)
		path.LLRs[N - 1:] = softMess
		pm = np.zeros(N, dtype=float)
		state = [[] for i in range(N)]
		indicator = np.zeros(N, dtype=int)
		i = -1
		pointer = -1
		self.iterations = 0
		T = self.threshold
		self.delta = 0.1
		self.bias = np.zeros(self.codewordLength, dtype=float)
		self.bias[self.polarMask == 1] = 1.35

		while flag:
			'''print(i)
			time.sleep(0.1)'''
			pointer = i + 1 
			if self.polarMask[pointer] == 1:

				path.updateLLRs(pcfun.bitreversed(pointer, n))
				edgeValue0 = 0
				edgeValue1 = 1

				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.bias[pointer]
				#branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.bias[pointer]
				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[pointer]
				#branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - self.I[pointer]
				branchMetric0 = np.log(((np.exp(path.LLRs[0]) / (np.exp(path.LLRs[0]) + 1))) / (1 - self.pe[pointer]))
				branchMetric1 = np.log(((1 / (np.exp(path.LLRs[0]) + 1))) / (1 - self.pe[pointer]))
				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - (K / N)
				#branchMetric1 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue1))) - (K / N)
				pathMetric0 = 0 + branchMetric0 if pointer == 0 else pm[pointer - 1] + branchMetric0
				pathMetric1 = 0 + branchMetric1 if pointer == 0 else pm[pointer - 1] + branchMetric1
				metricMax = pathMetric0 if (pathMetric0 > pathMetric1) else pathMetric1
				metricMin = pathMetric0 if (pathMetric0 < pathMetric1) else pathMetric1
				if indicator[pointer] == 0:
					pm[pointer] = metricMax
				else:
					pm[pointer] = metricMin
			else:

				path.updateLLRs(pcfun.bitreversed(pointer, n))
				edgeValue0 = 0

				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - self.I[pointer]
				branchMetric0 = np.log(((np.exp(path.LLRs[0]) / (np.exp(path.LLRs[0]) + 1))) / (1 - self.pe[pointer]))
				#branchMetric0 = 1 - np.log2(1 + np.exp(path.LLRs[0]) ** (-(1 - 2 * edgeValue0))) - (K / N)
				pathMetric0 = 0 + branchMetric0 if pointer == 0 else pm[pointer - 1] + branchMetric0
				pathMetric1 = -np.inf
				pm[pointer] = pathMetric0

			if pm[pointer] >= T:

				if self.polarMask[pointer] == 1:
					self.iterations += 1
				muPre = -np.inf if pointer == -1 else 0 if pointer == 0 else pm[pointer - 1]
				if muPre < (T + self.delta):
					while (T + self.delta) <= pm[pointer]:
						T = T + self.delta
				i += 1
				path.uHat[i] = 0 if pm[i] == pathMetric0 else 1 
				path.updateBits(pcfun.bitreversed(i, n))

				if (i + 1) == N:
					return self.extract(path.uHat)
				indicator[i + 1] = 0
			else:
				j = i
				(i, T, indicator) = self.moveBack2(i, T, pm, indicator)
				if j != i:
					for ii in range(0, i + 1):
						idxRev = pcfun.bitreversed(ii, n)
						path.updateLLRs(idxRev)
						path.updateBits(idxRev)

	def FanoANDSphere(self, y, softMess):
		self.iterations = 0
		N = self.codewordLength
		uHat = np.zeros(N, dtype=int)
		dHat = self.fanoPolar(softMess)
		uHat[self.polarMask == 1] = dHat
		xHat = self.polarEncode(dHat)
		xCorrect = self.polarEncode(self.correct)
		eucDistance = 0
		eucDistance2 = 0
		for i in range(N):
			eucDistance += ((1 - y[i]) / 2 -  xHat[i]) ** 2
			eucDistance2 += ((1 - y[i]) / 2 - xCorrect[i]) ** 2
		squareRadius = eucDistance
		#print(squareRadius, eucDistance2)
		return self.sphereDecoderDynamic(y, squareRadius, uHat)

	'''def computeP(self, K, y):
		if K == 0:
			i = 0  
		else:
			i = self.A[K]
		#tmp1 = 0
		tmp2 = 0
		for k in range(i, self.codewordLength):
			tmp1 = 0
			for j in range(k, self.codewordLength):
				tmp1 += self.GN[j][k] * self.uHat[j]
			tmp1 %= 2
			tmp2 += (y[k] - (1 - 2 * tmp1)) ** 2
		return tmp2

	def spherePolar(self, K, y, counter=0): #y = x + n
		if counter == 0:
			self.uHat = np.zeros(self.codewordLength, dtype=int)
			self.uOptimal = np.zeros(self.codewordLength, dtype=int)
			P = np.zeros(self.infoLen, dtype=float)
			self.squareRadius = np.inf
			self.LB = np.array([(y[i] - 1) ** 2 if ((y[i] - 1) ** 2 < (y[i] + 1) ** 2) else (y[i] + 1) ** 2 for i in range(self.codewordLength)])
			K = K - 1
			#print('hi')
		#print(K)
		#time.sleep(0.1)
		#print(self.extract(self.uOptimal))
		for i in range(2):
			if i == 0:
				self.uHat[self.A[K]] = 0
				#self.P[K] = self.computeP(K, y) if (K + 1) == self.infoLen else self.P[K + 1] + self.computeP(K, y)
				P = self.computeP(K, y)
			else:
				self.uHat[self.A[K]] = 1
				#self.P[K] = self.computeP(K, y) if (K + 1) == self.infoLen else self.P[K + 1] + self.computeP(K, y)
				P = self.computeP(K, y)

			if P > self.squareRadius:
				continue
			else:
				if K == 0:
					self.uOptimal = copy.deepcopy(self.uHat) #caution
					self.squareRadius = P
					#print(self.squareRadius)
					#print(self.extract(self.uOptimal))
				else:
					self.spherePolar(K - 1, y, counter + 1)

		if K == self.infoLen - 1:
			return self.extract(self.uOptimal)'''

	def computeMetric2(self, uHat, i, y):
		tmp = 0  
		for j in range(i, self.codewordLength):
			tmp += self.GN[j][i] * uHat[j]
		tmp %= 2
		return (((1 - y[i]) / 2) - tmp) ** 2 

	def computeMetric(self, uHat, i, y):
		tmp = 0  
		for j in range(i, self.codewordLength):
			tmp += self.GN[j][i] * uHat[j]
		tmp %= 2
		return (y[i] - (1 - 2 * tmp)) ** 2

	def lowerBound(self, LB, l):
		tmp = 0  
		for i in range(l):
			tmp += LB[i]
		return tmp

	def sphereDecoderFixed(self, y, squareRadius=np.inf, uOptimal=[]):

		N = self.codewordLength
		K = self.infoLen
		#LB = np.array([(y[i] - 1) ** 2 if ((y[i] - 1) ** 2 < (y[i] + 1) ** 2) else (y[i] + 1) ** 2 for i in range(self.codewordLength)])
		LB = np.array([(((1 - y[i]) / 2) - 0) ** 2 if (((1 - y[i]) / 2) - 0) ** 2 < (((1 - y[i]) / 2) - 1) ** 2 else (((1 - y[i]) / 2) - 1) ** 2 for i in range(N)])
		counter = np.zeros(N, dtype=int)
		uHat = np.zeros(N, dtype=int)
		d = np.zeros(N, dtype=float)
		#squareRadius = np.inf
		omiga = [0, 1]
		flag = True
		i = N - 1
		self.iterations = 0

		while flag:

			self.iterations += 1
			if self.polarMask[i] == 1:
				uHat[i] = omiga[counter[i]]
				d[i] = 0 + self.computeMetric2(uHat, i, y) if i == N - 1 else d[i + 1] + self.computeMetric2(uHat, i, y)
			else:
				uHat[i] = 0
				d[i] = 0 + self.computeMetric2(uHat, i, y) if i == N - 1 else d[i + 1] + self.computeMetric2(uHat, i, y)

			if d[i] + self.lowerBound(LB, i) > squareRadius:

				if self.polarMask[i] == 1:
					counter[i] += 1
					if counter[i] == len(omiga):
						counter[i] %= 2

						i += 1
						if i == N:
							'''print(d[i - 1], squareRadius)
							input('0')'''
							return self.extract(uOptimal)
						while self.polarMask[i] == 0 or (counter[i] + 1 == len(omiga)):
							i += 1
							if i == N:
								return self.extract(uOptimal)
						counter[:i] = 0
						counter[i] += 1
				else:

					i += 1
					if i == N:
						input('1')
					while self.polarMask[i] == 0 or (counter[i] + 1 == len(omiga)):
						i += 1
						if i == N:
							return self.extract(uOptimal)
					counter[:i] = 0
					counter[i] += 1
			else:

				if i == 0:
					uOptimal = copy.deepcopy(uHat)
					squareRadius = d[0]
					#print(squareRadius)

					i = self.A[0]
					counter[i] += 1
					if counter[i] == len(omiga):
						counter[i] %= 2

						i += 1
						if i == N:
							input('2')
						while self.polarMask[i] == 0 or (counter[i] + 1 == len(omiga)):
							i += 1
							if i == N:
								return self.extract(uOptimal)
						counter[:i] = 0
						counter[i] += 1
				else:
					i -= 1

	def dynamicBound(self, l, y, I, uHat):
		IL = I[l]
		m = 0
		for I in IL:
			m0 = 0  
			m1 = 0
			for i in I:
				t = 0
				for j in range(l, self.codewordLength):
					t += self.GN[j, i] * uHat[j]
				m0 += ((1 - y[i]) / 2 - ((0 + t) % 2)) ** 2
				m1 += ((1 - y[i]) / 2 - ((1 + t) % 2)) ** 2
			m += min(m0, m1)
		return m

	def sphereDecoderDynamic(self, y, squareRadius=np.inf, uOptimal=[]):

		N = self.codewordLength
		K = self.infoLen
		(I, d1) = self.IandD() 
		LB = np.array([(((1 - y[i]) / 2) - 0) ** 2 if (((1 - y[i]) / 2) - 0) ** 2 < (((1 - y[i]) / 2) - 1) ** 2 else (((1 - y[i]) / 2) - 1) ** 2 for i in range(N)])
		counter = np.zeros(N, dtype=int)
		uHat = np.zeros(N, dtype=int)
		d = np.zeros(N, dtype=float)
		omiga = [0, 1]
		flag = True
		i = N - 1
		self.iterations = 0

		while flag:

			self.iterations += 1
			if self.polarMask[i] == 1:
				uHat[i] = omiga[counter[i]]
				d[i] = 0 + self.computeMetric2(uHat, i, y) if i == N - 1 else d[i + 1] + self.computeMetric2(uHat, i, y)
			else:
				uHat[i] = 0
				d[i] = 0 + self.computeMetric2(uHat, i, y) if i == N - 1 else d[i + 1] + self.computeMetric2(uHat, i, y)

			#pathMetric = d[i] + (self.lowerBound(LB, i) if d1[i] == 0 else self.dynamicBound(i, y, I, uHat))
			pathMetric = d[i] + (self.dynamicBound(i, y, I, uHat) if d1[i] > 0 and self.polarMask[i] == 1 else self.lowerBound(LB, i))
			if pathMetric > squareRadius:

				if self.polarMask[i] == 1:
					counter[i] += 1
					if counter[i] == len(omiga):
						counter[i] %= 2

						i += 1
						if i == N:
							return self.extract(uOptimal)#[:-self.crcWidth]
						while self.polarMask[i] == 0 or (counter[i] + 1 == len(omiga)):
							i += 1
							if i == N:
								return self.extract(uOptimal)#[:-self.crcWidth]
						counter[:i] = 0
						counter[i] += 1
				else:

					i += 1
					if i == N:
						input('1')
					while self.polarMask[i] == 0 or (counter[i] + 1 == len(omiga)):
						i += 1
						if i == N:
							return self.extract(uOptimal)#[:-self.crcWidth]
					counter[:i] = 0
					counter[i] += 1
			else:

				if i == 0:
					uOptimal = copy.deepcopy(uHat)
					squareRadius = d[0]
					#print(squareRadius)

					i = self.A[0]
					counter[i] += 1
					if counter[i] == len(omiga):
						counter[i] %= 2

						i += 1
						if i == N:
							input('2')
						while self.polarMask[i] == 0 or (counter[i] + 1 == len(omiga)):
							i += 1
							if i == N:
								return self.extract(uOptimal)#[:-self.crcWidth]
						counter[:i] = 0
						counter[i] += 1
				else:
					i -= 1

	def CA_SD(self, y, squareRadius=np.inf, uOptimal=[], softMess=[]):

		N = self.codewordLength
		K = self.infoLen
		(I, d1) = self.IandD() 
		LB = np.array([(((1 - y[i]) / 2) - 0) ** 2 if (((1 - y[i]) / 2) - 0) ** 2 < (((1 - y[i]) / 2) - 1) ** 2 else (((1 - y[i]) / 2) - 1) ** 2 for i in range(N)])
		counter = np.zeros(N, dtype=int)
		uHat = np.zeros(N, dtype=int)
		d = np.zeros(N, dtype=float)
		omiga = [0, 1]
		flag = True
		i = N - 1
		Q = self.TransformPCRs(self.PCRs())
		P = [min(Q[l]) for l in range(len(Q))]

		'''if softMess.size > 0:
			vPrime = self.scDecoder(softMess)
			check = np.dot(vPrime, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			if np.sum(check) == 0:
				return vPrime[:-self.crcWidth]
			else:

				vPrime[-self.crcWidth:] = (np.dot(vPrime[:-self.crcWidth], pcfun.getGC(self.infoLen - self.crcWidth, self.crcPolyArrayForm)) % 2)[-self.crcWidth:]
				uOptimal = np.zeros(self.codewordLength, dtype=int)
				uOptimal[self.polarMask == 1] = vPrime
				eucDistance = 0
				for i in range(N):
					eucDistance += ((1 - y[i]) / 2 -  self.polarEncode(vPrime)[i]) ** 2
				squareRadius = eucDistance
				#print(squareRadius)'''

		#self.iterations = 0

		while flag:

			#self.iterations += 1
			if self.polarMask[i] == 1 and i not in P:
				uHat[i] = omiga[counter[i]]
				d[i] = 0 + self.computeMetric2(uHat, i, y) if i == N - 1 else d[i + 1] + self.computeMetric2(uHat, i, y)
			elif self.polarMask[i] == 0:
				uHat[i] = 0
				d[i] = 0 + self.computeMetric2(uHat, i, y) if i == N - 1 else d[i + 1] + self.computeMetric2(uHat, i, y)
			elif i in P:
				for l, kl in enumerate(P):
					if kl == i:
						break
				tmp = 0
				for t in Q[l]:
					tmp += uHat[t] if t != kl else 0
				tmp %= 2
				uHat[i] = tmp
				d[i] = 0 + self.computeMetric2(uHat, i, y) if i == N - 1 else d[i + 1] + self.computeMetric2(uHat, i, y)
			else:
				input('warning')

			#pathMetric = d[i] + (self.lowerBound(LB, i) if d1[i] == 0 else self.dynamicBound(i, y, I, uHat))
			pathMetric = d[i] + (self.dynamicBound(i, y, I, uHat) if d1[i] > 0 and self.polarMask[i] == 1 else self.lowerBound(LB, i))
			#pathMetric = d[i] +  self.dynamicBound(i, y, I, uHat)
			#pathMetric = d[i] + (self.dynamicBound(i, y, I, uHat) if d1[i] > 0 and self.polarMask[i] == 1 else 0)

			if self.polarMask[i] == 0 and d1[i] > 0:
				input('error warning, CA-SD')
			
			if pathMetric > squareRadius:

				if self.polarMask[i] == 1 and i not in P:
					counter[i] += 1
					if counter[i] == len(omiga):
						counter[i] %= 2

						i += 1
						if i == N:
							return self.extract(uOptimal)[:-self.crcWidth]
						while self.polarMask[i] == 0 or (counter[i] + 1 == len(omiga)) or i in P:
							i += 1
							if i == N:
								#print(uOptimal)
								return self.extract(uOptimal)[:-self.crcWidth]
						counter[:i] = 0
						counter[i] += 1
				else:

					i += 1
					if i == N:
						input('1')
					while self.polarMask[i] == 0 or (counter[i] + 1 == len(omiga)) or i in P:
						i += 1
						if i == N:
							#print(uOptimal)
							return self.extract(uOptimal)[:-self.crcWidth]
					counter[:i] = 0
					counter[i] += 1
			else:

				if i == 0:
					uOptimal = copy.deepcopy(uHat)
					squareRadius = d[0]
					#print(squareRadius)

					i = self.A[0]
					counter[i] += 1
					if counter[i] == len(omiga):
						counter[i] %= 2

						i += 1
						if i == N:
							input('2')
						while self.polarMask[i] == 0 or (counter[i] + 1 == len(omiga)) or i in P:
							i += 1
							if i == N:
								return self.extract(uOptimal)[:-self.crcWidth]
						counter[:i] = 0
						counter[i] += 1
				else:
					i -= 1

	def CA_HD(self, y, softMess, LMax=0, squareRadius=np.inf): #to be determined
		self.L = 1
		self.sortNum = 0
		while self.L <= LMax:
			print(f'L = {self.L}, LMax = {LMax}, AD-SCL')
			self.pathList = [Path(self.codewordLength)]
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
			#self.ANV = 0

			for i in range(self.codewordLength):

				idxRev = pcfun.bitreversed(i, self.n)
				for path in self.pathList:
					path.updateLLRs(idxRev)	

				if self.polarMask[i] == 1:
					self.polarSCLFork(i)
				else:

					for path in self.pathList:

						path.uHat[i] = self.polarMask[i]
						'''penalty = np.abs(path.LLRs[0])
						if path.LLRs[0] == 0:
							input('warning')
						branchMetric0 = 0 if (path.LLRs[0] > 0) else penalty'''
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						path.pathMetric += branchMetric0

				for path in self.pathList:
					path.updateBits(idxRev)


			self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
				
			'''best = self.extract(self.pathList[0].uHat)
			candidate = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

			if (np.sum(candidate[-self.crcWidth:])) == 0:
				return best[:len(best) - self.crcWidth]
			else:
				for path in self.pathList[1:]:
					nextBest = self.extract(path.uHat)
					nextCandidate = np.dot(nextBest, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2

					if (np.sum(nextCandidate[-self.crcWidth:])) == 0:
						return nextBest[:len(nextBest) - self.crcWidth]'''
			for l in range(len(self.pathList)):
				best = self.extract(self.pathList[l].uHat)
				check = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2
				if (np.sum(check[-self.crcWidth:])) == 0:
					print(f'AD-SCL succeeded in L = {self.L}, CA-HD')
					return best[:-self.crcWidth]
			self.L *= 2

		for l in range(LMax):
			u = self.extract(self.pathList[l].uHat)
			b = u[:-self.crcWidth]
			s = np.dot(b, pcfun.getGC(self.infoLen - self.crcWidth, self.crcPolyArrayForm)) % 2
			u[-self.crcWidth:] = s[-self.crcWidth:]
			xTilde = self.polarEncode(u)

			eucDistance = 0
			for i in range(self.codewordLength):
				eucDistance += (((1 - y[i]) / 2 ) - xTilde[i]) ** 2
			if eucDistance < squareRadius:
				squareRadius = eucDistance
		print(f'ready to start CA-SD, squared radius = {squareRadius}, CA-HD')
		return
		return self.CA_SD(y, squareRadius=squareRadius, softMess=np.array([]))

	def CA_HD1(self, y, softMess, LMax=0, squareRadius=np.inf):
		self.L = 1
		self.sortNum = 0
		self.mT = -14
		while self.L <= LMax:
			print(f'L = {self.L}, LMax = {LMax}, AD-SCL')
			self.pathList = [Path(self.codewordLength)]
			self.pathList[0].LLRs[self.codewordLength - 1:] = softMess
			#self.ANV = 0

			for i in range(self.codewordLength):

				idxRev = pcfun.bitreversed(i, self.n)
				for path in self.pathList:
					path.updateLLRs(idxRev)	

				if self.polarMask[i] == 1:
					self.PSCLFork(i)
				else:

					for path in self.pathList:

						path.uHat[i] = self.polarMask[i]
						Li = path.LLRs[0] / np.log(2)
						branchMetric0 = 1 - np.log2(1 + 2 ** (-Li * ((-1) ** 0)))
						path.pathMetric += branchMetric0

				for path in self.pathList:
					path.updateBits(idxRev)


			self.pathList.sort(key=lambda path_List: path_List.pathMetric, reverse=True)
				
			for l in range(len(self.pathList)):
				best = self.extract(self.pathList[l].uHat)
				check = np.dot(best, pcfun.getGC(self.infoLen, self.crcPolyArrayForm)) % 2
				if (np.sum(check[-self.crcWidth:])) == 0:
					print(f'AD-SCL succeeded in L = {self.L}, CA-HD1')
					return best[:-self.crcWidth]
			self.L *= 2

		for l in range(LMax):
			u = self.extract(self.pathList[l].uHat)
			b = u[:-self.crcWidth]
			s = np.dot(b, pcfun.getGC(self.infoLen - self.crcWidth, self.crcPolyArrayForm)) % 2
			u[-self.crcWidth:] = s[-self.crcWidth:]
			xTilde = self.polarEncode(u)

			eucDistance = 0
			for i in range(self.codewordLength):
				eucDistance += (((1 - y[i]) / 2 ) - xTilde[i]) ** 2
			if eucDistance < squareRadius:
				squareRadius = eucDistance
		print(f'ready to start CA-SD, squared radius = {squareRadius}, CA-HD1')
		return
		return self.CA_SD(y, squareRadius=squareRadius, softMess=np.array([]))

# ---------------- Rowshan's realization ----------------
	
	def trellisFork(self, pos):

		edgeValue = [0 for i in range(2*self.currNumPaths)]
		msgValue = [0 for i in range(2*self.currNumPaths)]
		pathMetric = [0.0 for i in range(2*self.currNumPaths)]
		pathState = [[] for i in range(2*self.currNumPaths)]
		pathStateMap =  [[] for i in range(int(2 ** self.m))]

		for i in range(self.currNumPaths):
			i2 = i+self.currNumPaths
			if self.trellisPath[i].LLRs[0] > 0:
				edgeValue[i] = pcfun.conv1Bit(0, self.trellisPath[i].curState, self.gen)
				edgeValue[i2] = 1 - edgeValue[i]
				pathMetric[i] = self.trellisPath[i].pathMetric + (0 if edgeValue[i]==0 else 1) * np.abs(self.trellisPath[i].LLRs[0])
				pathMetric[i2] = self.trellisPath[i].pathMetric + (0 if edgeValue[i2]==0 else 1) * np.abs(self.trellisPath[i].LLRs[0])
				if pathMetric[i2] > pathMetric[i]:
					msgValue[i] = 0
					msgValue[i2] = 1
					pathState[i] = pcfun.getNextState(0, self.trellisPath[i].curState, self.m)
					pathState[i2] = pcfun.getNextState(1, self.trellisPath[i].curState, self.m)
				else:
					edgeValue[i] = 1 - edgeValue[i]
					edgeValue[i2] = 1 - edgeValue[i2]
					tempPM = pathMetric[i]
					pathMetric[i] = pathMetric[i2]
					pathMetric[i2] = tempPM
					msgValue[i] = 1
					msgValue[i2] = 0
					pathState[i] = pcfun.getNextState(1, self.trellisPath[i].curState, self.m)
					pathState[i2] = pcfun.getNextState(0, self.trellisPath[i].curState, self.m)
	
			else:
				edgeValue[i] = pcfun.conv1Bit(1, self.trellisPath[i].curState, self.gen)
				edgeValue[i2] = 1 - edgeValue[i]
				pathMetric[i] = self.trellisPath[i].pathMetric + (0 if edgeValue[i]==1 else 1) * np.abs(self.trellisPath[i].LLRs[0])
				pathMetric[i2] = self.trellisPath[i].pathMetric + (0 if edgeValue[i2]==1 else 1) * np.abs(self.trellisPath[i].LLRs[0])
				if pathMetric[i2] > pathMetric[i]:
					msgValue[i] = 1
					msgValue[i2] = 0
					pathState[i] = pcfun.getNextState(1, self.trellisPath[i].curState, self.m)
					pathState[i2] = pcfun.getNextState(0, self.trellisPath[i].curState, self.m)
				else:
					edgeValue[i] = 1 - edgeValue[i]
					edgeValue[i2] = 1 - edgeValue[i2]
					tempPM = pathMetric[i]
					pathMetric[i] = pathMetric[i2]
					pathMetric[i2] = tempPM
					msgValue[i] = 0
					msgValue[i2] = 1
					pathState[i] = pcfun.getNextState(0, self.trellisPath[i].curState, self.m)
					pathState[i2] = pcfun.getNextState(1, self.trellisPath[i].curState, self.m)

			pathStateMap[pcfun.bin2dec(pathState[i])].append(i)
			pathStateMap[pcfun.bin2dec(pathState[i2])].append(i2)

		if 2*self.currNumPaths <= self.L:
			for i in range(self.currNumPaths):
				i2 = i+self.currNumPaths
				copy_path = copy.deepcopy(self.trellisPath[i])
				self.trellisPath[i].pathMetric = pathMetric[i]
				self.trellisPath[i].vHat[pos] = msgValue[i]
				self.trellisPath[i].uHat[pos] = edgeValue[i]
				self.trellisPath[i].curState = pathState[i]
				copy_path.pathMetric = pathMetric[i2]
				copy_path.vHat[pos] = msgValue[i2]
				copy_path.uHat[pos] = edgeValue[i2]
				copy_path.curState = pathState[i2]
				self.trellisPath.append(copy_path)
		else:
			discarded_paths = np.zeros(2*self.currNumPaths, dtype=np.int8)
			survived_paths = np.zeros(2*self.currNumPaths, dtype=np.int8)
			duplicated_paths = np.zeros(self.currNumPaths, dtype=np.int8)
			swapping_paths = np.zeros(self.currNumPaths, dtype=np.int8)
			retaining_paths = np.zeros(self.currNumPaths, dtype=np.int8)
			survivors = []
			deleted_paths = []
			num_states = int(2 ** self.m)

			for i in range(num_states):
				num_branches = len(pathStateMap[i])
				prune_start_idx = int(num_branches/2)
				if num_branches > 0:
					branchMetrics =  np.zeros(num_branches, dtype=float)
					for k in range(num_branches):
						branchMetrics[k] = pathMetric[pathStateMap[i][k]]
					bM_sorted_idx = np.argsort(branchMetrics)

					for k in range(prune_start_idx,num_branches):
						if pathStateMap[i][bM_sorted_idx[k]] < self.currNumPaths: #if ? > self.curr_num_paths-1, it will be discarded anyway and we don't need to copy another path to its place.
							discarded_paths[pathStateMap[i][bM_sorted_idx[k]]] = 1
					for k in range(prune_start_idx):
						if pathStateMap[i][bM_sorted_idx[k]] > self.currNumPaths-1:
							survived_paths[pathStateMap[i][bM_sorted_idx[k]]] = 1
							survivors.append(pathStateMap[i][bM_sorted_idx[k]])

			for i in range(self.currNumPaths):
				if discarded_paths[i] == 1 and survived_paths[i+self.currNumPaths] == 1:
					swapping_paths[i] = 1
					discarded_paths[i] = 0
				elif discarded_paths[i] == 1 and survived_paths[i+self.currNumPaths] == 0:
					deleted_paths.append(i) 
				elif discarded_paths[i] == 0 and survived_paths[i+self.currNumPaths] == 1:
					duplicated_paths[i] = 1 #1; duplicated, 0: deleted 
				elif discarded_paths[i] == 0 and survived_paths[i+self.currNumPaths] == 0:
					retaining_paths[i] = 1

			for i in range(self.currNumPaths): # In order of the paths stored in the memory
				if  swapping_paths[i] == 1: # Swapping the i-th path with i2-th path
					self.trellisPath[i].vHat[pos] = msgValue[i+self.currNumPaths]
					self.trellisPath[i].uHat[pos] = edgeValue[i+self.currNumPaths]
					self.trellisPath[i].curState = pathState[i+self.currNumPaths]
					self.trellisPath[i].pathMetric = pathMetric[i+self.currNumPaths]
					survivors.remove(i+self.currNumPaths)
				elif  retaining_paths[i] == 1: #the i-th path retained, not the i2-th path
					self.trellisPath[i].vHat[pos] = msgValue[i]
					self.trellisPath[i].uHat[pos] = edgeValue[i]
					self.trellisPath[i].curState= pathState[i]
					self.trellisPath[i].pathMetric = pathMetric[i]
				elif  duplicated_paths[i] == 1: #Issue: when duplicating, if there is no deleted path? Can't be
					self.trellisPath[i].vHat[pos] = msgValue[i]
					self.trellisPath[i].uHat[pos] = edgeValue[i]
					self.trellisPath[i].curState = pathState[i]
					self.trellisPath[i].pathMetric = pathMetric[i]
					self.trellisPath[deleted_paths[0]] = copy.deepcopy(self.trellisPath[i]) #in v2, the index of deleted_paths[0] might be > self.curr_num_paths
					self.trellisPath[deleted_paths[0]].vHat[pos] = msgValue[i+self.currNumPaths]
					self.trellisPath[deleted_paths[0]].uHat[pos] = edgeValue[i+self.currNumPaths]
					self.trellisPath[deleted_paths[0]].curState = pathState[i+self.currNumPaths]
					self.trellisPath[deleted_paths[0]].pathMetric = pathMetric[i+self.currNumPaths]
					deleted_paths.remove(deleted_paths[0])
					survivors.remove(i+self.currNumPaths)

	def pac_viterbi_decoder(self, softMess, issystematic=False):

		self.trellisPath = [Path(self.codewordLength, self.m)]
		self.trellisPath[0].LLRs[self.codewordLength - 1:] = softMess
		self.ANV = 0
		self.runningTime = 0

		t1 = perf_counter()
		for i in range(self.codewordLength):

			idxRev = pcfun.bitreversed(i, self.n)
			self.currNumPaths = len(self.trellisPath)
			for path in self.trellisPath:

				path.updateLLRs(idxRev)
				if self.polarMask[i] == 1:
					self.ANV += 1

			if self.polarMask[i] == 1:
				self.trellisFork(i)
			else:

				for path in self.trellisPath:

					edgeValue0 = pcfun.conv1Bit(0, path.curState, self.gen)
					path.uHat[i] = edgeValue0
					path.vHat[i] = self.polarMask[i]
					path.curState = pcfun.getNextState(0, path.curState, self.m)
					penalty = np.abs(path.LLRs[0])
					if path.LLRs[0] > 0:
						branchMetric0 = 0 if (edgeValue0 == 0) else penalty
					elif path.LLRs[0] < 0:
						branchMetric0 = 0 if (edgeValue0 == 1) else penalty
					else:
						input('warning')
					path.pathMetric += branchMetric0
					self.ANV += 1

			for path in self.trellisPath:
				path.updateBits(idxRev)

		self.trellisPath.sort(key=lambda trellisPath: trellisPath.pathMetric, reverse=False)
		#print(len(self.trellisPathList))
		best = self.trellisPath[0].vHat
		t2 = perf_counter()
		self.runningTime = (t2 - t1) * 1000
		return self.extract(best)