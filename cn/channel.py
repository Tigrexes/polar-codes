import numpy as np 

class Channel:
	def __init__(self, modulation, snrdB, snrType, Rc):
		self.modulation = modulation
		self.M = 2 if modulation.lower() == 'bpsk' else None
		self.noisePower = self.calc_N0(snrdB, snrType, Rc)

	def calc_N0(self, snrdB, snrType, Rc):

		if snrType.lower() == 'snr':
			return 1 / pow(10, snrdB / 10)
		else:
			return 1 / (np.log2(self.M) * Rc * pow(10, snrdB / 10))

	def modulate(self, m):
		modulated = []
		if self.modulation.lower() == 'bpsk':
			modulated = [(1 - 2 * x) for x in m]
			return modulated
		else:
			pass

	def modulateMultiple(self, ms):
		modulated = []
		if self.modulation.lower() == 'bpsk':
			for m in ms:
				modulated.append([(1 - 2 * x) for x in m])
			return modulated
		else:
			pass

	def addNoise(self, signal):
		if self.modulation.lower() == 'bpsk':
			return signal + np.sqrt(self.noisePower / 2) * np.random.standard_normal(len(signal))

	def addNoise2(self, signal1, signal2):
		if self.modulation.lower() == 'bpsk':

			noise = np.sqrt(self.noisePower / 2) * np.random.standard_normal(len(signal1))
			return signal1 + noise, signal2 + noise

	def addNoiseMultiple(self, signals):
		if self.modulation.lower() == 'bpsk':

			signalLength = len(signals[0])
			noise = np.sqrt(self.noisePower / 2) * np.random.standard_normal(signalLength)

			for i in range(len(signals)):
				signals[i] += noise 
			return signals

	def calcLLR(self, c):
		llr = []
		if self.modulation.lower() == 'bpsk':
			llr = [(4 / self.noisePower * y) for y in c]
			llr = np.array(llr)
			return llr 
		else:
			pass

	def calcLLRMultiple(self, cs):
		llr = []
		if self.modulation.lower() == 'bpsk':
			for c in cs:
				llr.append(np.array([4 / self.noisePower * y for y in c]))
			return llr 
		else:
			pass