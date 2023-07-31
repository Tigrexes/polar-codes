from channel import Channel
from polarCode import PolarCode
from time import time, perf_counter
import numpy as np
import functions as pcf
import matplotlib.pyplot as plt
import platform
import time

#crcPoly = 0x107 # 0b100000111 0x107, 0x07 is standard form
#crcPoly = 0x10001E5
crcPoly = 0x1800063 # for progressive bit flipping
#crcPoly = 0x1621
#crcPoly = 0xA22E5
#1010 0010 0010 1110 0101
#crcWidth = len(hex(crcPoly)[3:]) * 4
crcWidth = 24#len([int(bin(crcPoly)[2:][i]) for i in range(len(bin(crcPoly)[2:]))]) - 1

N = 2 ** 10 # block length
R = 512 / 1024 # code rate
K = int(N * R) # number of information bits
K2 = K + crcWidth
R2 = K2 / N
construct = "dega" # 'dega': density evolution with gaussian approximation, 'rmPolar': Reed-Muller Polar, 'WS': Weighted Sum
dsnr = 2.25 # design SNR for dega
convGen = [1, 0, 1, 1, 0, 1, 1]#[1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1], [1, 0, 1, 1, 0, 1, 1]
m = len(convGen) - 1

#localStackSize = 20
stackSize = 100
localListSize = 2 ** 0
listSize = localListSize * 2 ** m
maxLevel = 3


iterations = 10 ** 8
errCnt = 100 # total number of error to count at each SNR 

snrType = 'SNRb' # 'SNRb':Eb/N0, 'SNR':Es/N0
modu = 'BPSK'
'''
polarCoders = []

sva1 = PolarCode(N, K, construct, dsnr, gen=[1, 0, 1, 1, 0, 1, 1])
sva1.convGenName = str([1, 0, 1, 1, 0, 1, 1])
sva1.decoderType = 'sva'
sva1.D = 100

sva2 = PolarCode(N, K, construct, dsnr, gen=[1, 0, 1, 1, 0, 1, 1])
sva2.convGenName = str([1, 0, 1, 1, 0, 1, 1])
sva2.decoderType = 'pacstack'
sva2.D = (100 + 2 ** m)
print(sva2.D)

sva3 = PolarCode(N, K, construct, dsnr, gen=[1, 0, 1, 1, 0, 1, 1])
sva3.convGenName = str([1, 0, 1, 1, 0, 1, 1])
sva3.decoderType = 'cs-aided'
sva3.D = 100

sva4 = PolarCode(N, K, construct, dsnr, gen=[1, 1, 1])
sva4.convGenName = str([1, 1, 1])
sva4.decoderType = 'rowshanlistviterbi'
sva4.L = 32

sva5 = PolarCode(N, K, construct, dsnr, gen=[1, 0, 1, 1, 0, 1, 1])
sva5.convGenName = str([1, 0, 1, 1, 0, 1, 1])
sva5.decoderType = 'paclist'
sva5.L = 32

polarCoders.append(sva1)
polarCoders.append(sva2)
polarCoders.append(sva3)
polarCoders.append(sva4)
polarCoders.append(sva5)
'''

pcode = PolarCode(N, K, construct, dsnr, gen=convGen)
pcode.snrType = snrType
pcode.modu = modu
pcode.gen = convGen
pcode.m = m
pcode.D = stackSize
pcode.L = 256
pcode.crcPoly = crcPoly
pcode.crcWidth = crcWidth
#pcode.crcPolyArrayForm = np.array([1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1], dtype=int)#np.array([int(bin(crcPoly)[2:][i]) for i in range(len(bin(crcPoly)[2:]))], dtype=int)
#pcode.crcPolyArrayForm = np.array([1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0,
#									1, 1, 1, 0, 0, 1, 0, 1], dtype=int)
#1010 0010 0010 1110 0101
pcode.crcPolyArrayForm = np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1], dtype=int) #CRC-24
pcode.cutoff_rate = pcf.cutoffRate(N, K, dsnr)
pcode.pe = pcf.PEDega(N, K, dsnr)
pcode.I = pcf.symmetricCapacity(N, K, dsnr)
pcode.GA = pcf.GA(N, K, dsnr)
pcode.maxLevel = maxLevel

pcode2 = PolarCode(N, K, 'dega', dsnr, gen=convGen)
pcode2.snrType = snrType
pcode2.modu = modu
pcode2.gen = convGen
pcode2.m = m
pcode2.D = stackSize
#pcode2.cutoff_rate = pcf.cutoffRate(N, K2, dsnr)
pcode2.I = pcf.symmetricCapacity(N, K, dsnr)


snrRange = np.arange(2.25, 7.0, 0.5)

class BERFER():
	def __init__(self):
		self.snrRange = list()
		self.ber = list()
		self.fer = list()

result = BERFER()

for snr in snrRange:

	print(f'\nSNR = {snr} dB')
	ber = 0
	fer = 0 

	ber2 = 0 
	fer2 = 0

	errCount = 0
	failureCount = 0
	ANV = 0
	ANV2 = 0
	errTotal = []

	#ch = Channel(pcode.modu, snr, pcode.snrType, R)
	ch = Channel(modu, snr, snrType, R)
	#ch = Channel(modu, snr, snrType, 1 / 2)

	for t in range(iterations):

		message = np.random.randint(0, 2, size=K, dtype=int)

		pcode.message = np.zeros(N, dtype=int)
		pcode.message[pcode.polarMask == 1] = message
		#pcode.message = np.zeros(N, dtype=int)
		#pcode.message[pcode.polarMask == 1] = pcf.crcEncode(message, crcPoly)
		#pcode.correct = np.zeros(N, dtype=int)
		#pcode.correct = message
		

		x = pcode.polarEncode(message)
		#x2 = pcode.polarEncode(message)
		modulatedX = ch.modulate(x)
		#modulatedX2 = ch.modulate(x2)
		y = ch.addNoise(modulatedX)
		#y, y2 = ch.addNoise2(modulatedX, modulatedX2)
		llr_ch = ch.calcLLR(y)
		#llr_ch2 = ch.calcLLR(y2)
		decoded = pcode.sclDecoderOracle(llr_ch)
		ANV += pcode.sortNum
		print(pcode.unitCal, pcode.sortNum)
		'''if pcode.errCount > maxLevel:
			print(f'Errors in channel = {pcode.errCount}')
			#errBits = 1
			fer += 1
			fer2 += 1
			print(f'Error # {fer} t = {t}, FER = {fer / (t + 1):0.2e}, sort num = {ANV / (t + 1):.2f}, genie')
			print(f'Error # {fer2} t = {t}, FER = {fer2 / (t + 1):0.2e}, scNum = {scNum}, lowPBF')
			continue
		else:
			print(f'Errors in channel = {pcode.errCount}')
			errBits = 0'''

		decoded2 = pcode.sclDecoderOracleGCA(llr_ch)
		ANV2 += pcode.sortNum
		print(pcode.unitCal, pcode.sortNum)
		input('dd')
		'''scNum = pcode.scNum
		if (message == decoded2).all():
			print(f'Correct # {fer2} t = {t}, FER = {fer2 / (t + 1):0.2e}, scNum = {scNum}, lowPBF')'''
		
		
		#ANV += pcode.iterations
		'''if (decoded == message).all():
			pass#print('correct\n')
		else:
			correct = pcode.polarEncode(message)
			error = pcode.polarEncode(decoded)
			eucDistance1 = 0
			eucDistance2 = 0
			for i in range(N):
				eucDistance1 += (((1 - y[i]) / 2 ) - correct[i]) ** 2
				eucDistance2 += (((1 - y[i]) / 2 ) - error[i]) ** 2
				#eucDistance1 += (y[i] - (1 - 2 * correct[i])) ** 2
				#eucDistance2 += (y[i] - (1 - 2 * error[i])) ** 2
			print(eucDistance1, eucDistance2)
			#print(message)
			#input('pause')
		
		if pcode.errCount > 0:

			errTotal.append(pcode.errCount)
			failureCount += 1

			if failureCount == 4000:

				errCount = [0 for i in range(max(errTotal))]

				for i in range(len(errTotal)):
					errCount[errTotal[i] - 1] += 1

				errCount = np.array(errCount)
				y = (errCount / len(errTotal))
				x = np.array(range(len(errCount))) + 1
				plt.bar(x, y)
				plt.show()'''
		
		
		errBits = pcf.fails(message, decoded)
		ber += errBits

		
		errBits2 = pcf.fails(message, decoded2)
		ber2 += errBits2
		

		if errBits > 0:
			fer += 1
			print(f'Error # {fer} t = {t}, FER = {fer / (t + 1):0.2e}, sort num = {ANV / (t + 1):.2f}, SCL')
			#print(f'Error # {fer} t = {t}, FER = {fer / (t + 1):0.2e}, ANV = {ANV / (t + 1)}')
			'''correct = pcode.polarEncode(message)
			error = pcode.polarEncode(decoded)
			eucDistance1 = 0
			eucDistance2 = 0
			for i in range(N):
				eucDistance1 += (1 / 2 * (y[i] + 1) - correct[i]) ** 2
				eucDistance2 += (1 / 2 * (y[i] + 1) - error[i]) ** 2
			if eucDistance2 < eucDistance1:
				print('warning')
				input('pause')
			print(eucDistance1, eucDistance2)
			print('\n')'''
		
		'''if errBits > 0:
			correct = pcode.pacEncode(message)
			error = pcode.pacEncode(decoded)
			eucDistance1 = 0
			eucDistance2 = 0
			for i in range(N):
				eucDistance1 += (1 / 2 * (y[i] + 1) - correct[i]) ** 2
				eucDistance2 += (1 / 2 * (y[i] + 1) - error[i]) ** 2
			print(eucDistance1, eucDistance2)
			input('pause')'''

		if errBits2 > 0:
			fer2 += 1
			print(f'Error # {fer2} t = {t}, FER = {fer2 / (t + 1):0.2e}, sort num = {ANV2 / (t + 1):.2f}, GCASCL')
		
		
		if (fer >= errCnt) and (t >= 99999):
			print(f'\nAt {snr} dB FER is {fer / (t + 1):0.2e}')
			break
		
		
		
		#if fer2 == errCnt:
			#print(f'\nAt {snr} dB FER is {fer2 / (t + 1):0.2e}')
			#break
		

		

		#if min(fer, fer2) == errCnt:
			#print(f'\nAt {snr} dB FER is {fer / (t + 1):0.2e}, n')
			#print(f'\nAt {snr} dB FER is {fer2 / (t + 1):0.2e}, v')
			#print(f'\nt = {t} ANV Normal = {ANVNormal / (t + 1)}, normal')
			#print(f'\nt = {t} ANV CriticalSet = {ANVCriticalSet / (t + 1)}, criticalSet')
			#break
		

		if (t % 2000 == 0) and (t != 0):
			print(f'\nt = {t} FER = {fer / (t + 1):0.2e}')
			

	#result.snrRange.append(snr)
	#result.ber.append(ber / ((t + 1) * K))
	#result.fer.append(fer / (t + 1))


'''
snrRange = np.arange(1, 3.5, 0.5)

class BERFER():
	def __init__(self):
		self.snrRange = list()
		self.ber = []
		self.fer = []

#result = BERFER()
bySNRFER = []
results = []
for i in range(len(polarCoders)):
	results.append(BERFER())

for snr in snrRange:

	print(f'\nSNR = {snr} dB')
	ber = np.array([0 for i in range(len(polarCoders))])
	fer = np.array([0 for i in range(len(polarCoders))])

	failureCount = 0
	errTotal = []

	ANVCriticalSet = 0
	ANVNormal = 0

	ch = Channel(modu, snr, snrType, R)

	for t in range(iterations):

		xs = []
		message = np.random.randint(0, 2, size=K, dtype=int)

		for pcode in polarCoders:
			xs.append(pcode.pacEncode(message))

		modulatedXs = ch.modulateMultiple(xs)
		ys = ch.addNoiseMultiple(modulatedXs)
		llr_chs = ch.calcLLRMultiple(ys)

		decodeds = []
		for i in range(len(polarCoders)):
			if i == len(polarCoders) - 1:
				decodeds.append(polarCoders[i].pacStackDecoder2(llr_chs[i]))
			else:
				decodeds.append(polarCoders[i].pacStackViterbiDecoder(llr_chs[i]))	
		

		errBits = [0 for i in range(len(decodeds))]

		for i in range(len(decodeds)):
			errBits[i] = pcf.fails(message, decodeds[i])
			ber[i] += errBits[i]

		for i in range(len(errBits)):

			if errBits[i] > 0:
				fer[i] += 1
				print(f'Error # {fer[i]} t = {t}, FER = {fer[i] / (t + 1):0.2e}, {polarCoders[i].convGenName}')
		
		
		if min(fer) == errCnt:
			for i in range(len(fer)):
				print(f'\nAt {snr} dB FER is {fer[i] / (t + 1):0.2e}, {polarCoders[i].convGenName}')
			break

	
	#result.ber.append(ber / ((t + 1) * K))
	#result.fer.append(fer / (t + 1))
	
	for i in range(len(polarCoders)):
		results[i].fer.append(fer[i] / (t + 1))
		results[i].ber.append(ber[i] / ((t + 1) * K))

	bySNRFER.append(fer / (t + 1))

'''

'''
snrRange = np.arange(0.0, 3.5, 0.5)

class BERFER():
	def __init__(self):
		self.snrRange = list()
		self.ber = []
		self.fer = []
		self.complexity = np.array([])
		self.runningTime = np.array([])

#result = BERFER()
bySNRFER = []
results = []
for i in range(len(polarCoders)):
	results.append(BERFER())

for snr in snrRange:

	print(f'\nSNR = {snr} dB')
	ber = np.array([0 for i in range(len(polarCoders))])
	fer = np.array([0 for i in range(len(polarCoders))])

	failureCount = 0
	errTotal = []

	complexity = np.array([0 for i in range(len(polarCoders))])
	runningTime = np.array([0 for i in range(len(polarCoders))])

	ch = Channel(modu, snr, snrType, R)

	for t in range(iterations):

		xs = []
		message = np.random.randint(0, 2, size=K, dtype=int)

		for pcode in polarCoders:
			xs.append(pcode.pacEncode(message))

		modulatedXs = ch.modulateMultiple(xs)
		ys = ch.addNoiseMultiple(modulatedXs)
		llr_chs = ch.calcLLRMultiple(ys)

		decodeds = []
		
		for i in range(len(polarCoders)):
			decodeds.append(polarCoders[i].decode(llr_chs[i]))

		for i in range(len(polarCoders)):
			complexity[i] += polarCoders[i].ANV
			runningTime[i] += polarCoders[i].runningTime
			#print(polarCoders[i].runningTime, polarCoders[i].ANV, polarCoders[i].decoderType)
		#input('pause')
		
		errBits = [0 for i in range(len(decodeds))]

		for i in range(len(decodeds)):
			errBits[i] = pcf.fails(message, decodeds[i])
			ber[i] += errBits[i]

		for i in range(len(errBits)):

			if errBits[i] > 0:
				fer[i] += 1
				#print(f'Error # {fer[i]} t = {t}, FER = {fer[i] / (t + 1):0.2e}, {polarCoders[i].convGenName}')
				print(f'Error # {fer[i]} t = {t}, FER = {fer[i] / (t + 1):0.2e}, decoder index = {i}, decoder type = {polarCoders[i].decoderType}')
		
		
		if min(fer) >= errCnt and t >= 9999: #careful
			for i in range(len(fer)):
				#print(f'\nAt {snr} dB FER is {fer[i] / (t + 1):0.2e}, {polarCoders[i].convGenName}')
				print(f'\nAt {snr} dB FER is {fer[i] / (t + 1):0.2e}, decoder index = {i}, decoder type = {polarCoders[i].decoderType}')
			break

	
	#result.ber.append(ber / ((t + 1) * K))
	#result.fer.append(fer / (t + 1))
	
	for i in range(len(polarCoders)):
		results[i].fer.append(fer[i] / (t + 1))
		results[i].ber.append(ber[i] / ((t + 1) * K))
		#results[i].complexity.append(complexity[i] / (t + 1))
		results[i].complexity = np.append(results[i].complexity, complexity[i] / (t + 1))
		results[i].runningTime = np.append(results[i].runningTime, runningTime[i] / (t + 1))

	bySNRFER.append(fer / (t + 1))
'''
#-----------------------complexity plot-----------------------#
'''
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.tick_params(top='on', right='on', which='both')
plt.xticks(snrRange)
markers = ['o', 'v', 's', 'd', '*', 'x']
colors = ['red', 'violet', 'lawngreen', 'blue', 'black', 'yellow']

for i in range(len(results)):
	if i == 3:
		plt.plot(snrRange, results[i].complexity / results[3].complexity, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
			color=colors[i], lw='1.8', clip_on=False, alpha=1, 
			label=f'SC')
	elif i == 4:
		plt.plot(snrRange, results[i].complexity / results[3].complexity, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
			color=colors[i], lw='1.8', linestyle='dotted', clip_on=False, alpha=1, 
			label=f'SCL (L = {polarCoders[i].L})')
	elif i == 5:
		plt.plot(snrRange, results[i].complexity / results[3].complexity, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
			color=colors[i], lw='1.8', clip_on=False, alpha=1, 
			label=f'SCS (D = {polarCoders[i].D})')
	else:
		plt.plot(snrRange, results[i].complexity / results[3].complexity, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
			color=colors[i], lw='1.8', clip_on=False, alpha=1, 
			label=f'SC-Fano (Δ = {polarCoders[i].delta})')

plt.legend(framealpha=1, edgecolor='dimgray')
#plt.title('PAC(128,64), LVA vs SVA')
plt.xlabel('${E_b}$/${N_0}$[dB]', fontsize=12)
plt.ylabel('Normalized Complexity', fontsize=12)

plt.grid(True, which='major', linestyle='-')
plt.grid(True, which='minor', linestyle='-')
plt.ylim(0, max(results[0].complexity / results[3].complexity))
plt.xlim(min(snrRange), max(snrRange))
plt.show()
'''
#-----------------------FER plot-----------------------#
'''
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.tick_params(top='on', right='on', which='both')
plt.xticks(snrRange)
markers = ['o', 'v', 's', 'd', '*', 'x']
colors = ['steelblue', 'darkorange', 'forestgreen', 'crimson', 'mediumorchid', 'sienna']
#colors = ['red', 'violet', 'lawngreen', 'blue', 'black', 'yellow']

for i in range(len(results)):
	if i == 0:
		plt.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
			color=colors[i], lw='1.8', linestyle='dotted', clip_on=False, alpha=1, 
			label=f'SVD, D={polarCoders[i].D}, g={polarCoders[i].convGenName}')
	elif i == 1:
		plt.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
			color=colors[i], lw='1.8', linestyle='dotted', clip_on=False, alpha=1, 
			label=f'LVA, L={polarCoders[i].L}, g={polarCoders[i].convGenName}')
	elif i == 2:
		plt.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
			color=colors[i], lw='1.8', linestyle='dotted', clip_on=False, alpha=1, 
			label=f'SVD, D={polarCoders[i].D}, g={polarCoders[i].convGenName}')
	elif i == 3:
		plt.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
			color=colors[i], lw='1.8', linestyle='dotted', clip_on=False, alpha=1, 
			label=f'LVA, L={polarCoders[i].L}, g={polarCoders[i].convGenName}')
	elif i == 4:
		plt.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
			color=colors[i], lw='1.8', linestyle='dotted', clip_on=False, alpha=1, 
			label='LVA, ' + '${L_G}$=' + f'{polarCoders[i].L}, g={polarCoders[i].convGenName}')

plt.legend(loc='lower left', framealpha=1, edgecolor='dimgray')
#plt.title('PAC(128,64), LVA vs SVA')
plt.xlabel('${E_b}$/${N_0}$[dB]', fontsize=12)
plt.ylabel('FER', fontsize=12)


count = 0
if min(bySNRFER[-1]) >= (10 ** -4):
	for i in str(min(bySNRFER[-1])):
		#print('hi')
		if i == '0':
			count += 1
		elif i != '.':
			break
else:
	count = str(min(bySNRFER[-1]))[-1]
count = int(count)

plt.yscale('log')
plt.grid(True, which='major', linestyle='-')
plt.grid(True, which='minor', linestyle='dotted')
plt.ylim(10 ** (-count), 1)
plt.xlim(min(snrRange), max(snrRange))
plt.show()
'''
#-----------------------Double YY Axes plot-----------------------#

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
fig = plt.figure()
ax1 = fig.add_subplot()

markers = ['o', 'v', 's', 'd', '*', 'x']
colors = ['steelblue', 'darkorange', 'forestgreen', 'crimson', 'mediumorchid', 'sienna']

for i in range(len(results)):
	if i == 0:
		ax1.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='-', clip_on=False, alpha=1, 
				label=f'SVD, D={polarCoders[i].D}, FER')
	elif i == 1:
		ax1.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='-', clip_on=False, alpha=1, 
				label=f'CSD, D={polarCoders[i].D}, FER')
	elif i == 2:
		ax1.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='-', clip_on=False, alpha=1, 
				label=f'CS-Aided, D={polarCoders[i].D}, FER')
	elif i == 3:
		ax1.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='-', clip_on=False, alpha=1, 
				label=f'LVA, L={polarCoders[i].L}, FER, in [3]')
	elif i == 4:
		ax1.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='-', clip_on=False, alpha=1, 
				label=f'LD, ' + '${L_G}$=' + f'{polarCoders[i].L}, FER, in [2]')

ax2 = ax1.twinx()
for i in range(len(results)):
	if i == 0:
		ax2.plot(snrRange, results[i].runningTime, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='--', clip_on=False, alpha=1, 
				label=f'SVD, D={polarCoders[i].D}, Running Time')
	elif i == 1:
		ax2.plot(snrRange, results[i].runningTime, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='--', clip_on=False, alpha=1, 
				label=f'CSD, D={polarCoders[i].D}, Running Time')
	elif i == 2:
		ax2.plot(snrRange, results[i].runningTime, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='--', clip_on=False, alpha=1, 
				label=f'CS-Aided, D={polarCoders[i].D}, Running Time')
	elif i == 3:
		ax2.plot(snrRange, results[i].runningTime, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='--', clip_on=False, alpha=1, 
				label=f'LVA, L={polarCoders[i].L}, Running Time, in [3]')
	elif i == 4:
		ax2.plot(snrRange, results[i].runningTime, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='--', clip_on=False, alpha=1, 
				label=f'LD, ' + '${L_G}$=' + f'{polarCoders[i].L}, Running Time, in [2]')

plt.xlim(min(snrRange), max(snrRange))

ax1.set_yscale('log')
ax1.grid()
#ax1.tick_params(direction='in')
ax1.set_xlabel('${E_b}$/${N_0}$[dB]')
ax1.set_ylabel('FER')

#ax2.tick_params(direction='in')
ax2.set_ylabel('Running Time (ms)')
fig.legend(loc=1, bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)

plt.show()

#-----------------------Double YY Axes plot2-----------------------#

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
fig = plt.figure()
ax1 = fig.add_subplot()

markers = ['o', 'v', 's', 'd', '*', 'x']
colors = ['steelblue', 'darkorange', 'forestgreen', 'crimson', 'mediumorchid', 'sienna']

for i in range(len(results)):
	if i == 0:
		ax1.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='-', clip_on=False, alpha=1, 
				label=f'SVD, D={polarCoders[i].D}, FER')
	elif i == 1:
		ax1.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='-', clip_on=False, alpha=1, 
				label=f'CSD, D={polarCoders[i].D}, FER')
	elif i == 2:
		ax1.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='-', clip_on=False, alpha=1, 
				label=f'CS-Aided, D={polarCoders[i].D}, FER')
	elif i == 3:
		ax1.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='-', clip_on=False, alpha=1, 
				label=f'LVA, L={polarCoders[i].L}, FER, in [3]')
	elif i == 4:
		ax1.plot(snrRange, results[i].fer, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='-', clip_on=False, alpha=1, 
				label=f'LD, ' + '${L_G}$=' + f'{polarCoders[i].L}, FER, in [2]')

ax2 = ax1.twinx()
for i in range(len(results)):
	if i == 0:
		ax2.plot(snrRange, results[i].complexity, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='--', clip_on=False, alpha=1, 
				label=f'SVD, D={polarCoders[i].D}, ANV')
	elif i == 1:
		ax2.plot(snrRange, results[i].complexity, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='--', clip_on=False, alpha=1, 
				label=f'CSD, D={polarCoders[i].D}, ANV')
	elif i == 2:
		ax2.plot(snrRange, results[i].complexity, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='--', clip_on=False, alpha=1, 
				label=f'CS-Aided, D={polarCoders[i].D}, ANV')
	elif i == 3:
		ax2.plot(snrRange, results[i].complexity, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='--', clip_on=False, alpha=1, 
				label=f'LVA, L={polarCoders[i].L}, ANV, in [3]')
	elif i == 4:
		ax2.plot(snrRange, results[i].complexity, marker=markers[i], markersize='6', markerfacecolor='None', mew=1.8,
				color=colors[i], lw='1.8', linestyle='--', clip_on=False, alpha=1, 
				label=f'LD, ' + '${L_G}$=' + f'{polarCoders[i].L}, ANV, in [2]')

plt.xlim(min(snrRange), max(snrRange))

ax1.set_yscale('log')
ax1.grid()
#ax1.tick_params(direction='in')
ax1.set_xlabel('${E_b}$/${N_0}$[dB]')
ax1.set_ylabel('FER')

#ax2.tick_params(direction='in')
ax2.set_ylabel('ANV')
fig.legend(loc=1, bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)

plt.show()

for i in range(len(results)):
	print(results[i].runningTime, results[i].complexity)