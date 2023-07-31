from operator import itemgetter
import matplotlib.pyplot as plt
import numpy as np
import math
import copy

class Node: #designed for progressive bit-flipping algorithm
	def __init__(self):
		'''self.curU = -1 
		self.preU = []'''
		self.sequenceU = []

# ---------------- Auxiliary Functions ----------------

def isAllZero(x):
	if x.size > 0:
		for entry in x:
			if entry != 0:
				return False
		return True
	else:
		return True

def rowEchelonForm(A): #caution: only suitable for int type
	columnNum = len(A[0])
	rowNum = len(A)
	indexes = []
	for column in range(columnNum):
		lock = False
		for row in range(rowNum):
			if isAllZero(A[row, :column]) and A[row, column] != 0 and lock == False:
				tmp = copy.deepcopy(A[row, :])
				indexes.append(row)
				lock = True
			elif isAllZero(A[row, :column]) and A[row, column] != 0 and lock == True:
				common = A[row, column] * tmp[column] / math.gcd(A[row, column], tmp[column])
				A[row, :] *= int(common / A[row, column])
				A[row, :] += tmp * int(-common / tmp[column])
				#A[row, :] += tmp
				A[row, :] %= 2
				#head = A[row, column]
				#A[row, :] = -tmp[column] * A[row, :]
				#A[row, :] += tmp * head

				#A[row, :] += -(A[row, column] / tmp[column]) * tmp
			else:
				continue
	if len(indexes) != rowNum:
		for row in range(rowNum):
			if row not in indexes:
				indexes.append(row)
	return A[indexes]

def polynomialMod(a, b):
	lenA = len(a)
	lenB = len(b)
	lenC = lenA - lenB + 1
	#c = [0 for i in range(lenC)]
	if lenC < 0:
		input('Warning')
	else:
		c = np.zeros(lenC, dtype=int)
		tmp = a[:lenB]
		for i in range(lenC):
			c[i] = 0 if tmp[0] == 0 else 1
			tmp = (tmp + b * c[i]) % 2
			if i != lenC - 1:
				tmp = np.append(tmp, a[lenB + i])
				tmp = tmp[-lenB:]
	return tmp[1:]
	#return tmp

def JFunction(t):
	return (1 - 2 ** (-0.3073 * (t ** (2 * 0.8935)))) ** 1.1064

def JFunctionInv(t):
	return (-1 / 0.3073 * np.log2(1 - (t ** (1 / 1.1064)))) ** (1 / (2 * 0.8935))

def int_to_binlist(num: int, size: int):
    """"""
    return [int(bit) for bit in bin(num)[2:].zfill(size)]

def standardForm(crcPoly):
	
	size = len(bin(crcPoly)[3:])
	arrayForm = np.array([int(bit) for bit in bin(crcPoly)[3:].zfill(size)])
	num = 0

	for i in range(len(arrayForm)):
		num += arrayForm[i] * 2 ** (len(arrayForm) - 1 - i)
		
	return arrayForm, '0x' + hex(num)[2:].zfill(len(hex(crcPoly)) - len(hex(num)))

def Factorial(n: int):
	if n == 0:
		return 1 
	else:
		mul = 1
		for i in range(1, n + 1):
			mul *= i
		return mul

def combination(n: int, r: int):

	return int(Factorial(n) / (Factorial(r) * Factorial(n - r)))

def fails(list1, list2):
    """returns number of bit errors"""
    return np.sum(np.absolute(list1 - list2))

def bitreversed(num: int, n) -> int:
    """"""
    return int(''.join(reversed(bin(num)[2:].zfill(n))), 2)

def lowerconv(upperdecision, upperllr, lowerllr):
	if upperdecision == 0:
		return lowerllr + upperllr
	else:
		return lowerllr - upperllr

def countOnes(num: int):

	ones = 0
	binary = bin(num)[2:]
	lenBin = len(binary)

	for i in range(lenBin):
		if binary[i] == '1':
			ones += 1
	return ones

def GRowsWeight(N: int, K: int):

	w = np.zeros(N, dtype=int)

	for i in range(N):
		w[i] = countOnes(i)

	return w

def phi_inv(x: float):
    if (x>12):
        return 0.9861 * x - 2.3152
    elif (x<=12 and x>3.5):
        return x*(0.009005 * x + 0.7694) - 0.9507
    elif (x<=3.5 and x>1):
        return x*(0.062883*x + 0.3678) - 0.1627
    else:
        return x*(0.2202*x + 0.06448)

def bin2dec(binary): 
    decimal = 0
    for i in range(len(binary)): 
        decimal = decimal + binary[i] * pow(2, i) 
    return decimal

def PEDega(N: int, K: int, dsnr_db: float):

	mllr = GA(N, K, dsnr_db)
	pe = np.zeros(N, dtype=float)

	for i in range(N):
		pe[i] = (1 / 2) * math.erfc(1 / math.sqrt(2) * math.sqrt(mllr[i] / 2))

	return pe


def GA1(N: int, K: int, dsnr_db:float):

	n = int(math.log2(N))
	sigma_sq = 1/(2*K/N*np.power(10,dsnr_db/10))
	mllr = np.zeros(2 * N - 1, dtype=float)
	mllr[N - 1:] = 2 / sigma_sq
	mllr1 = []

	for i in range(N):
		ii = bitreversed(i, n)
		if ii == 0:
			nextlevel = n 
		else:
			lastlevel = (bin(ii)[2:].zfill(n)).find('1') + 1
			start = int(np.power(2, lastlevel - 1)) - 1
			end = int(np.power(2, lastlevel) - 1) - 1

			for i in range(start, end + 1):
				mllr[i] = 2 * mllr[end + (i - start) + 1]

			nextlevel = lastlevel - 1

		for lev in range(nextlevel, 0, -1):

			start = int(np.power(2, lev - 1)) - 1
			end = int(np.power(2, lev) - 1) - 1
			for indx in range(start, end + 1):
				exp1 = end + (indx - start)
				mllr[indx] = phi_inv(mllr[exp1 + 1])
		mllr1.append(mllr[0])
	return mllr1


def GA(N: int, K: int, dsnr_db: float): #Gaussian Approximation
	sigma_sq = 1/(2*K/N*np.power(10,dsnr_db/10))
	mllr = [2 / sigma_sq]
	while len(mllr) != N:
		tmp = []
		for i in range(len(mllr)):
			tmp.append(phi_inv(mllr[i]))
			tmp.append(2 * mllr[i])
		mllr = tmp
	mllr = np.array(mllr)
	return mllr

def symmetricCapacity(N: int, K: int, dsnr_db: float):
	sigma_sq = 1/(2*K/N*np.power(10,dsnr_db/10))
	t = [JFunction(2 / sigma_sq)]
	while len(t) != N:
		tmp = []
		for i in range(len(t)):
			tmp.append(1 - JFunction(np.sqrt(2) * JFunctionInv(1 - t[i])))
			tmp.append(JFunction(np.sqrt(2) * JFunctionInv(t[i])))
		t = tmp
	t = np.array(t)
	return t

def bhatta(N: int, K: int, dsnr_db: float):

	mllr = GA(N, K, dsnr_db)
	sigmaSq = 2 / mllr
	bhattaPara = np.zeros(N, dtype=float)

	for i in range(N):
		bhattaPara[i] = np.exp(-1 / (2 * sigmaSq[i]))

	return bhattaPara

def cutoffRate(N:int, K:int, dsnr_db:float):

	bhattaPara = bhatta(N, K, dsnr_db)
	cutoff_rate = np.zeros(N, dtype=float)

	for i in range(N):
		cutoff_rate[i] = np.log2(2 / (1 + bhattaPara[i]))

	return cutoff_rate

def GA_Es(N: int, dsnr_db: float):
	sigma_sq = 1/(2*np.power(10,dsnr_db/10))
	mllr = [2 / sigma_sq]
	while len(mllr) != N:
		tmp = []
		for i in range(len(mllr)):
			tmp.append(phi_inv(mllr[i]))
			tmp.append(2 * mllr[i])
		mllr = tmp
	mllr = np.array(mllr)
	return mllr

def bhatta_Es(N: int, dsnr_db: float):

	mllr = GA_Es(N, dsnr_db)
	sigmaSq = 2 / mllr
	bhattaPara = np.zeros(N, dtype=float)

	for i in range(N):
		bhattaPara[i] = np.exp(-1 / (2 * sigmaSq[i]))

	return bhattaPara

def cutoffRate_Es(N:int, dsnr_db:float):

	bhattaPara = bhatta_Es(N, dsnr_db)
	cutoff_rate = np.zeros(N, dtype=float)

	for i in range(N):
		cutoff_rate[i] = np.log2(2 / (1 + bhattaPara[i]))

	return cutoff_rate

def generateCriticalSet(frozenBits):

	N = frozenBits.size
	n = int(np.log2(N))
	hw = []
	criticalSet = []
	A = -1 * np.ones((n + 1, N), dtype=int)

	for i in range(N):
		A[-1, i] = frozenBits[i]

	for i in range(n - 1, -1, -1):
		for j in range(0, np.power(2, i)):
			A[i, j] = A[i + 1, 2 * j] + A[i + 1, 2 * j + 1]

	for i in range(0, n + 1):
		for j in range(0, np.power(2, i)):
			if A[i, j] == 0:

				index_1 = j 
				index_2 = j 

				for k in range(1, n - i + 1):

					index_1 = 2 * index_1
					index_2 = 2 * index_2 + 1

					for p in range(index_1, index_2 + 1):
						A[i + k, p] = -1
				criticalSet.append(index_1)
	criticalSet.sort()
	return np.array(criticalSet)

def modifyCriticalSet(polarMask, iMax):

	infoSet = np.zeros(len(polarMask), dtype=int)
	for i in range(len(polarMask)):
		if (i > iMax) and (polarMask[i] == 1):
			infoSet[i] = 1 
	frozenSet = (infoSet + 1) % 2
	return generateCriticalSet(frozenSet)

def constructTree(polarMask, maxLevel):

	S = [[] for i in range(maxLevel)]
	nodes = [[] for i in range(maxLevel)]
	frozenMask = (1 + polarMask) % 2
	l = 0 

	while l <= (maxLevel - 1):
		if l == 0:
			SPrime = generateCriticalSet(frozenMask)
			for u in SPrime:
				node = Node()
				node.sequenceU.append(u)
				nodes[l].append(node)
				S[l].append(u)

		else:
			curNode = 0
			while curNode < len(S[l - 1]):
				SPrime = modifyCriticalSet(polarMask, S[l - 1][curNode])
				for u in SPrime:
					node = Node()
					node.sequenceU += nodes[l - 1][curNode].sequenceU
					node.sequenceU.append(u)
					nodes[l].append(node)
					S[l].append(u)
				curNode += 1
		l += 1
	return nodes

def getFn(n :int, F):
	
	if n == 1:
		return np.kron(F, 1)
	else:
		return np.kron(getFn(n - 1, F), F)

def GN(N, F): #GN = BN * Fn or Fn
	#N = int(2 ** n)
	#indicesRev = [bitreversed(i, n) for i in range(N)]
	n = int(np.log2(N))
	return getFn(n, F)

def transform4Sphere(l, GHat):
	if l == 0:
		return [] #None
	else:
		columns = [[GHat[:l, i], i] for i in range(l)]
		I = []
		while columns:

			test1 = np.zeros(len(columns[0][0]), dtype=int)
			#print()
			if (columns[0][0] == test1).all():
				del columns[0]
				continue
			else:
				column = columns[0]
				del columns[0]
				tmp = []

			for i in range(len(columns)):
				if (column[0] == columns[i][0]).all():
					tmp.append(columns[i][1])
			if tmp:
				tmp.append(column[1])
				for i in tmp:
					for index, column in enumerate(columns):
						if i == column[1]:
							del columns[index]
				tmp.sort(reverse=True)
				I.append(tmp)
				#I.reverse()
		I.reverse()
		return I

# ---------------- Reed-Muller Polar Construction ----------------

def rmPolarBuild(N: int, K: int, design_snr=0):

	mask = [[i, 0, 0.0, 1] for i in range(N)]
	values = GRowsWeight(N, K)
	values2 = GA(N, K, design_snr)

	for i in range(N):
		mask[i][1] = values[i]
		mask[i][2] = values2[i]

	weightCount = np.zeros(int(math.log2(N)) + 1, dtype=int)

	for i in range(N):
		weightCount[values[i]] += 1

	bitCount = 0
	k = 0

	while bitCount + weightCount[k] <= N - K:

		for i in range(N):

			if values[i] == k:
				mask[i][3] = 0 
				bitCount += 1
		k += 1

	mask2 = []

	for i in range(N):
		if mask[i][1] == k:
			mask2.append(mask[i])

	mask2 = sorted(mask2, key=itemgetter(2))
	remainder = (N - K) - bitCount
	#print(remainder)
	for i in range(remainder):
		mask[mask2[i][0]][3] = 0

	rateProfile = np.array([i[3] for i in mask])
	return rateProfile

def rmPolarBuild_Es(N: int, K: int, design_snr=0):

	mask = [[i, 0, 0.0, 1] for i in range(N)]
	values = GRowsWeight(N, K)
	values2 = GA_Es(N, design_snr)

	for i in range(N):
		mask[i][1] = values[i]
		mask[i][2] = values2[i]

	weightCount = np.zeros(int(math.log2(N)) + 1, dtype=int)

	for i in range(N):
		weightCount[values[i]] += 1

	bitCount = 0
	k = 0

	while bitCount + weightCount[k] <= N - K:

		for i in range(N):

			if values[i] == k:
				mask[i][3] = 0 
				bitCount += 1
		k += 1

	mask2 = []

	for i in range(N):
		if mask[i][1] == k:
			mask2.append(mask[i])

	mask2 = sorted(mask2, key=itemgetter(2))
	remainder = (N - K) - bitCount

	for i in range(remainder):
		mask[mask2[i][0]][3] = 0

	rateProfile = np.array([i[3] for i in mask])
	return rateProfile
# ---------------- Density Evolution With Gaussian Approximation Construction ----------------

def degaBuild(N: int, K: int, design_snr=0):
	mask = [[i, 0.0, 1] for i in range(N)]
	values = GA(N, K, design_snr)
	for i in range(N):
		mask[i][1] = values[i]
	mask = sorted(mask, key=itemgetter(1), reverse=False)
	for i in range(N - K):
		mask[i][2] = 0
	mask = sorted(mask, key=itemgetter(0))
	
	return np.array([i[2] for i in mask])

# ---------------- Convolution Encoder ----------------

def conv1Bit(inBit, curState, gen):

	gLen = len(gen)
	gBit = inBit * gen[0]

	for i in range(1, gLen):

		if gen[i] == 1:
			gBit ^= curState[i - 1]

	return gBit

def getNextState(inBit, curState, m):

	if inBit == 0:
		nextState = [0] + curState[0:m - 1]
	else:
		nextState = [1] + curState[0:m - 1]

	return nextState

def convEncode(inCode, gen, m):

	curState = [0 for i in range(m)]
	inCodeLen = len(inCode)
	convCode = np.zeros(inCodeLen, dtype=int)
	
	for i in range(inCodeLen):

		inBit = inCode[i]
		output = conv1Bit(inBit, curState, gen)
		curState = getNextState(inBit, curState, m)
		convCode[i] = output

	return convCode

# ---------------- CRC Encoder ----------------

def getGC(KI, g): #get generator matrix for CRC
	KC = len(g) - 1
	K = KC + KI
	identityMatrix = np.eye(KI)
	matrixRight = np.zeros([KI, KC])
	for i in range(KI):
		x = np.zeros(K - i, dtype=int)
		x[0] = 1
		matrixRight[i, :] = polynomialMod(x, g)
	GC = np.concatenate((identityMatrix, matrixRight), axis=1)
	GC = np.array(GC, dtype=int)
	return GC

def PCRs(KI, g):
	A = np.array(range(N), dtype=int)[degaBuild(64, 12, 2.5) == 1]
	GC = getGC(KI, g)
	print(GC)
	KC = len(g) - 1
	Rs = [[] for i in range(KC)]
	#Ru = [[] for i in range(KC)]
	Ru = []
	for l in range(KC):
		for i in range(KI):
			if GC[i, l + KI] == 1:
				Rs[l].append(i)
		Rs[l].append(KI + l)
	for l in range(KC):
		Ru.append([A[i] for i in Rs[l]])
	#print(Rs)
	return Ru

def crcEncode(message, crcPoly): # bitwise implementation caution!!! wrong!!!!
	
	crcPoly, hexForm = standardForm(crcPoly)
	registerLength = len(crcPoly)
	augmentedMessage = np.append(message, [0 for i in range(registerLength)])
	register = np.zeros(registerLength, dtype=int)

	for i in range(len(augmentedMessage)):

		topBit = register[0]
		register = np.append(register[1:registerLength], augmentedMessage[i])

		if topBit == 1:
			register ^= crcPoly 
	
	augmentedMessage[len(augmentedMessage) - registerLength:] = register 
	return augmentedMessage

def buildCRC8Table(crcPoly: int):

	generator = np.uint16(crcPoly)
	crc8Table = list()

	for div in range(256):
		curByte = div
		for bit in range(8):

			if np.bitwise_and(curByte, np.uint8(0x80)) != np.uint8(0x00):

				curByte <<= 1
				curByte ^= generator
			else:
				curByte <<= 1
		crc8Table.append(np.uint8(curByte))
	return crc8Table

def crc8TableMethod(info, crcTable):

	crc = 0 
	if (info.size % 8) != 0:
		pad0 = np.zeros((info.size // 8 * 8 + 8) - info.size, dtype=np.int8)
		info = np.append(pad0, info)
	coef = np.array([128, 64, 32, 16, 8, 4, 2, 1])
	for b in range(0, len(info), 8):
		pos = np.uint8((crc) ^ np.sum(info[b:b + 8] * coef))
		crc = crcTable[pos]
	return int_to_binlist(crc, 8)

def crcTableEncode(message, crcPoly, width):

	if width == 8:
		return np.append(message, crc8TableMethod(message, buildCRC8Table(crcPoly)))
	else:
		pass

# ---------------- Weighted Sum Based Construction ----------------

def getAandS(N: int, K: int, g):

	t = 0
	n = int(math.log2(N))
	for t in range(0, n):
		lowerBound = 0 
		upperBound = 0
		for p in range(t + 1, n + 1):

			lowerBound += combination(n, p)
		for p in range(t, n + 1):

			upperBound += combination(n, p)
		if lowerBound <= K < upperBound:
			t = t
			break

	GWeight = GRowsWeight(N, K)
	A = []
	S = []
	for i in range(len(GWeight)):
		if GWeight[i] > t:
			A.append(i)
		if GWeight[i] == t:
			S.append(i)

	return A, S 


def updateBandTau(A, N, g):
	b = np.zeros(N, dtype=int)
	b[A] = 1 #????
	'''
	for i in range(len(b)):
		if i in A:
			b[i] = 1 
	'''

	tau = np.zeros(N, dtype=int)
	for i in range(N):
		for j in range(len(g)):
			if (i - j) >= 0:
				tau[i] += b[i - j] * g[j]

	return tau 

def updateTheta(w, tau, g, N):

	theta = np.zeros(N, dtype=float)
	for i in range(len(theta)):
		for j in range(len(g)):
			if (i + j) < N:
				theta[i] += (g[j] * w[i + j]) / (tau[i + j] + 1)

	return theta

def updateAandS(A, S, theta):

	iStar = [[0.0, 0] for i in range(len(S))]
	for i in range(len(S)):
		iStar[i][0], iStar[i][1] = theta[S[i]], S[i]
	
	iStar = sorted(iStar, key=itemgetter(0), reverse=True)
	
	iStar = iStar[0][1]

	S.remove(iStar)
	A.append(iStar)
	A.sort()
	

	return A, S 

def WSConstruction(N, K, dsnr_db, g):

	A, S = getAandS(N, K, g)
	if len(A) < K:

		cutoff_rate = cutoffRate(N, K, dsnr_db)
		w = np.ceil(cutoff_rate / 0.1)

		while len(A) < K:

			tau = updateBandTau(A, N, g)
			theta = updateTheta(w, tau, g, N)
			A, S = updateAandS(A, S, theta)

	rateProfile = np.zeros(N, dtype=int)
	A = np.array(A)
	rateProfile[A] = 1

	return rateProfile

# ---------------- Draft ----------------

def findMin(theta, S):
	lock = True
	for i in range(len(theta)):
		if i not in S and lock:
			minNum = theta[i]
			lock = False
		elif i in S:
			continue
		if theta[i] < minNum and i not in S:
			minNum = theta[i]
	return minNum

def pacConstruction(N, K, dsnr_db, g):
	I = symmetricCapacity(N, K, dsnr_db)
	rmPolar = rmPolarBuild(N, K, dsnr_db)
	dega = degaBuild(N, K, dsnr_db)
	tau = np.zeros(N, dtype=int)
	theta = np.zeros(N, dtype=float)
	
	for i in range(len(tau)):
		for j in range(len(g)):
			if i >= j:
				tau[i] += dega[i - j] * g[j]
		
	w = I / (tau + 1)

	for i in range(N):
		for j in range(len(g)):
			if i + j < N:
				theta[i] += g[j] * w[i + j]
	print(tau)
	#print(theta)
	print(dega)
	c = 0
	for i, value in enumerate(tau):
		if dega[i] == 1:
			c += I[i] / value
	print(c)

def IandD(n, F, mask): #dynamic bound
	N = int(pow(2, n))
	GHat = getFn(n, F)
	for i in range(N):
		if mask[i] == 0:
			GHat[i, :] = 0
	#print(GHat)
	I = [transform4Sphere(i, GHat) for i in range(N)]
	d = np.zeros(N, dtype=int)
	for i in range(N):
		if mask[i] == 1:
			d[i] = len(I[i])
	#print(I)
	print(d)
	return (I, d)

def GA2(N: int, K: int, dsnr_db: float): #Gaussian Approximation
	sigma_sq = 1/(2*K/N*np.power(10,dsnr_db/10))
	mllr = [2 / sigma_sq]
	while len(mllr) != N:
		tmp = []
		for i in range(len(mllr)):
			tmp.append(phiInv(1 - (1 - phi(mllr[i])) ** 2))
			tmp.append(2 * mllr[i])
		mllr = tmp
	mllr = np.array(mllr)
	return mllr

def degaBuild2(N: int, K: int, design_snr=0):
	mask = [[i, 0.0, 1] for i in range(N)]
	values = GA2(N, K, design_snr)
	for i in range(N):
		mask[i][1] = values[i]
	mask = sorted(mask, key=itemgetter(1), reverse=False)
	for i in range(N - K):
		mask[i][2] = 0
	mask = sorted(mask, key=itemgetter(0))
	
	return np.array([i[2] for i in mask])

def phi(x):
	if (x >= 0) and (x <= 10):
		return np.exp(-0.4527 * x ** 0.859 + 0.0218)
	else:
		return np.sqrt(np.pi / x) * np.exp(-x / 4) * (1 - 10 / 7 / x)

def phiInv(y):
	if (y <= 1.0221) and (y >= 0.0388):
		x = ((0.0218 - np.log(y)) / 0.4527) ** (1 / 0.86)
		return x
	else:
		x0 = 0.0388
		x1 = x0 - (phi(x0) - y) / derivativePhi(x0)
		delta = np.abs(x1 - x0)
		epsilon = 1e-3
		while delta >= epsilon:
			x0 = x1
			x1 = x1 - (phi(x1) - y) / derivativePhi(x1)
			if x1 > 1e2:
				epsilon = 10
			delta = np.abs(x1 - x0)
		return x1

def derivativePhi(x):
	if (x >= 0) and (x <= 10):
		return -0.4527 * 0.86 * x ** (-0.14) * phi(x)
	else:
		return np.exp(-x / 4) * np.sqrt(np.pi / x) * (-1 / 2 / x * (1 - 10/ 7 / x) - 1 / 4 * (1 - 10/ 7 / x) + 10 / 7 / x / x)

def PEDega2(N: int, K: int, dsnr_db: float):

	mllr = GA2(N, K, dsnr_db)
	pe = np.zeros(N, dtype=float)

	for i in range(N):
		pe[i] = (1 / 2) * math.erfc(1 / math.sqrt(2) * math.sqrt(mllr[i] / 2))

	return pe
import scipy.special as ss
F = np.array([[1, 0], [1, 1]])
N = 2 ** 10
R = 512 / 1024
K = int(N * R)
n = np.log2(N)
dsnr = 2.25
g = [1, 0, 1, 1, 0, 1, 1]
#pacConstruction(N ,K, dsnr, g)
#mask = degaBuild(N, K, dsnr)
mask2 = degaBuild(N, K, dsnr)
frozenMask2 = (1 - mask2) % 2
pe = PEDega(N, K, dsnr)
mu = GA(N, K, dsnr)
idxRecord = []
for i in range(N):
	if mask2[i] == 1 and mu[i] < 59:
		idxRecord.append(i)

peSC = 1
CS = generateCriticalSet(frozenMask2)
for i in range(N):
	if mask2[i] == 1 and i not in idxRecord:
		peSC *= (1 - pe[i])
print(f'{1 - peSC:.4e}', len(idxRecord), len(CS))
'''def QFunction(x):
	return 1 / 2 * math.erfc(math.sqrt(2) / 2 * x)
def QFunctionInv(x):
	return math.sqrt(2) * ss.erfcinv(2 * x)
#print(QFunction((np.log(1 - pe[127]) - np.log(pe[127]) - mu[127]) / np.sqrt(2 * mu[127])))
corrPro = 1
count = 0
for i in range(N):
	if mask2[i] == 1:
		if not np.isnan(QFunction(QFunctionInv(pe[i]) + (np.log(1 / pe[i] - 1)) / (2 * QFunctionInv(pe[i])))):
			corrPro *= (1 - 1 * QFunction(QFunctionInv(pe[i]) + (np.log(1 / pe[i] - 1)) / (2 * QFunctionInv(pe[i]))))
print(f'{1 - corrPro:.4e}', count)'''
# ---------------- Novel Construction ----------------
