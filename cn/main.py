import numpy as np
from operator import itemgetter

def my_build_mask(N: int, K: int, design_snr=0):
	mask = [[i, 0.0, 1] for i in range(N)]
	values = ga(N, K, design_snr)
	for i in range(N):
		mask[i][1] = values[i]
	mask = sorted(mask, key=itemgetter(1), reverse=False)
	for i in range(N - K):
		mask[i][2] = 0
	mask = sorted(mask, key=itemgetter(0))
	
	return np.array([i[2] for i in mask])

def ga(N: int, K: int, dsnr_db: float):
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

def phi_inv(x: float):
    if (x>12):
        return 0.9861 * x - 2.3152
    elif (x<=12 and x>3.5):
        return x*(0.009005 * x + 0.7694) - 0.9507
    elif (x<=3.5 and x>1):
        return x*(0.062883*x + 0.3678)- 0.1627
    else:
        return x*(0.2202*x + 0.06448)

def dega_construct(N: int, K: int, dsnr_db: float):
    # bhattacharya_param = [0.0 for i in range(N)]
    mllr = np.zeros(N, dtype=float)
    # snr = pow(10, design_snr / 10)
    #dsnr = np.power(10, dsnr_db / 10)
    sigma_sq = 1/(2*K/N*np.power(10,dsnr_db/10))
    mllr[0] = 2/sigma_sq
    #mllr[0] = 4 * K/N * dsnr
    for level in range(1, int(np.log2(N)) + 1):
        B = np.power(2, level)
        for j in range(int(B / 2)):
            T = mllr[j]
            mllr[j] = phi_inv(T)
            mllr[int(B / 2 + j)] = 2 * T
    return mllr

def build_mask(N: int, K: int, design_snr=0):
    """Generates mask of polar code
    in mask 0 means frozen bit, 1 means information bit"""
    # each bit has 3 attributes
    # [order, bhattacharyya value, frozen / imformation position]
    # 0 - frozen, 1 - information
    mask = [[i, 0.0, 1] for i in range(N)]
    # Build mask using Bhattacharya values
    #values = G_rows_wt(N, K)
    values = dega_construct(N, K, design_snr)
    #values = bhattacharyya_count(N, design_snr)
    # set bhattacharyya values
    for i in range(N):
        mask[i][1] = values[i]
    # sort channels due to bhattacharyya values
    mask = sorted(mask, key=itemgetter(1), reverse=False)   #DEGA, RM
    #mask = sorted(mask, key=itemgetter(1), reverse=True)    #bhattacharyya
    # set mask[i][2] in 1 for channels with K lowest bhattacharyya values
    for i in range(N-K):
        mask[i][2] = 0
    # sort channels due to order
    mask = sorted(mask, key=itemgetter(0))
    
    # return positions bits
    return np.array([i[2] for i in mask])

def mul_matrix(precoded, N):

    n = int(np.log2(N))
    polarcoded = precoded
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

def bitreversed(num: int, n) -> int:
    """"""
    return int(''.join(reversed(bin(num)[2:].zfill(n))), 2)

def lowerconv(upperdecision, upperllr, lowerllr):
	if upperdecision == 0:
		return lowerllr + upperllr
	else:
		return lowerllr - upperllr

def update_llrs(position, llrs, bits):
	

	if position == 0:
		nextlevel = n 
	else:
		lastlevel = (bin(position)[2:].zfill(n)).find('1') + 1
		start = int(np.power(2, lastlevel - 1)) - 1
		end = int(np.power(2, lastlevel) - 1) - 1
		for i in range(start, end + 1):
			llrs[i] = lowerconv(bits[0][i],
								llrs[end + 2 * (i - start) + 1],
								llrs[end + 2 * (i - start) + 2])

		nextlevel = lastlevel - 1
	for lev in range(nextlevel, 0, -1):

		start = int(np.power(2, lev - 1)) - 1
		end = int(np.power(2, lev) - 1) - 1
		for indx in range(start, end + 1):
			exp1 = end + 2 * (indx - start) # exp1 = end + (indx - start)
			llr1 = llrs[exp1 + 1] #llr1 = llrs[exp1 + 1]
			llr2 = llrs[exp1 + 2] #llr2 = llrs[exp1 + 1 + 2 ** (lev - 1)]
			llrs[indx] = np.sign(llr1) * np.sign(llr2) * min(abs(llr1), abs(llr2))
	if position == 0:
		print(llrs)
		input('pause')
	return llrs

def update_bits(position, latestbit, bits):
	N = codewordLength
	if position == N - 1:
		return
	elif position < N // 2:
		bits[0][0] = latestbit
	else:
		lastlevel = (bin(position)[2:].zfill(n)).find('0') + 1
		bits[1][0] = latestbit
		for lev in range(1, lastlevel - 1):
			st = int(np.power(2, lev - 1)) - 1
			ed = int(np.power(2, lev) - 1) - 1
			for i in range(st, ed + 1):
				bits[1][ed + 2 * (i - st) + 1] = (bits[0][i] + bits[1][i]) % 2
				bits[1][ed + 2 * (i - st) + 2] = bits[1][i]

		lev = lastlevel - 1
		st = int(np.power(2, lev - 1)) - 1
		ed = int(np.power(2, lev) - 1) - 1
		for i in range(st, ed + 1):
			bits[0][ed + 2 * (i - st) + 1] = (bits[0][i] + bits[1][i]) % 2
			bits[0][ed + 2 * (i - st) + 2] = bits[1][i]
	return bits

def my_update_llrs(position, llrs, bits):
	
	

	if position == 0:
		nextlevel = n 
	else:
		lastlevel = (bin(position)[2:].zfill(n)).find('1') + 1
		start = int(np.power(2, lastlevel - 1)) - 1
		end = int(np.power(2, lastlevel) - 1) - 1
		
		for i in range(start, end + 1):
			
			llrs[i] = lowerconv(bits[0][i],
								llrs[end + (i - start) + 1],
								llrs[end + (i - start) + 1 + 2 ** (lastlevel - 1)])
			
		nextlevel = lastlevel - 1
	for lev in range(nextlevel, 0, -1):

		start = int(np.power(2, lev - 1)) - 1
		end = int(np.power(2, lev) - 1) - 1
		for indx in range(start, end + 1):
			exp1 = end + (indx - start)
			llr1 = llrs[exp1 + 1]
			llr2 = llrs[exp1 + 1 + 2 ** (lev - 1)]
			llrs[indx] = np.sign(llr1) * np.sign(llr2) * min(abs(llr1), abs(llr2))
	
	return llrs

def my_update_bits(position, latestbit, bits):
	N = codewordLength
	if position == N - 1:
		return
	elif position < N // 2:
		bits[0][0] = latestbit
	else:
		lastlevel = (bin(position)[2:].zfill(n)).find('0') + 1
		bits[1][0] = latestbit
		for lev in range(1, lastlevel - 1):
			st = int(np.power(2, lev - 1)) - 1
			ed = int(np.power(2, lev) - 1) - 1
			for i in range(st, ed + 1):
				bits[1][ed + (i - st) + 1] = (bits[0][i] + bits[1][i]) % 2
				bits[1][ed + (i - st) + 1 + 2 ** (lev - 1)] = bits[1][i]

		lev = lastlevel - 1
		st = int(np.power(2, lev - 1)) - 1
		ed = int(np.power(2, lev) - 1) - 1
		for i in range(st, ed + 1):
			bits[0][ed + (i - st) + 1] = (bits[0][i] + bits[1][i]) % 2
			bits[0][ed + (i - st) + 1 + 2 ** (lev - 1)] = bits[1][i]
	return bits

codewordLength = 1024#len(codeword)
log2N = int(np.log2(codewordLength))
indicesRev = [bitreversed(i, log2N) for i in range(codewordLength)]
R = 0.5
K = int(codewordLength * R) 


#polarMask = build_mask(codewordLength, K, 2.5)[indicesRev] #careful
polarMask = my_build_mask(codewordLength, K, 2.5)


message = np.random.randint(0, 2, size=K, dtype=int)
#message = np.array([0, 0, 0, 0, 0, 0, 1, 1])

u = np.zeros(codewordLength, dtype=int)
u[polarMask == 1] = message

'''print(message)
print(u)'''
x = mul_matrix(u[indicesRev], codewordLength)

'''print(x)'''
modulated_x = [(1 - 2 * c) for c in x]
'''print(modulated_x)'''
Rc = K / codewordLength
snrdB = 1
noisePower = 1 / (Rc * pow(10, snrdB / 10))
y = modulated_x + np.sqrt(noisePower / 2) * np.random.standard_normal(len(modulated_x))
'''print(y)'''
llrCh = np.array([(4 / noisePower * y) for y in y])


llrs = np.zeros(2 * codewordLength - 1, dtype=float)
bits = np.zeros((2, codewordLength - 1), dtype=int)
llrs[codewordLength - 1:] = llrCh

llrOb = []
n = int(np.log2(codewordLength))

t = 1
uHat = np.zeros(codewordLength, dtype=int)
for i in range(codewordLength):
	ii = bitreversed(i, n)
	llrs = my_update_llrs(ii, llrs, bits)
	llrOb.append(llrs[0])
	

	if polarMask[ii] == 1:
		if llrs[0] > 0:
			uHat[ii] = 0
		elif llrs[0] < 0:
			uHat[ii] = 1
		else:
			input('Warning')
		
		
		
	else:
		uHat[ii] = 0
	bits = my_update_bits(ii, uHat[ii], bits)

#uHat = uHat[indicesRev]
decoded = []
for i in range(len(polarMask)):
	if polarMask[i] == 1:
		decoded.append(uHat[i])
decoded = np.array(decoded, dtype=int)
print(f"decoded =   {decoded}")
print(f'message =   {message}')
print(f'polarMask = {polarMask}')
print(f'u =         {u}')
print(f'uHat =      {uHat}')
print((message == decoded).all())


'''
def my_update_llrs(position, llrs, bits):
	
	

	if position == 0:
		nextlevel = n 
	else:
		lastlevel = (bin(position)[2:].zfill(n)).find('1') + 1
		start = int(np.power(2, lastlevel - 1)) - 1
		end = int(np.power(2, lastlevel) - 1) - 1
		
		for i in range(start, end + 1):
			
			llrs[i] = lowerconv(bits[0][i],
								llrs[end + (i - start) + 1],
								llrs[end + (i - start) + 1 + 2 ** (lastlevel - 1)])
			
		nextlevel = lastlevel - 1
	for lev in range(nextlevel, 0, -1):

		start = int(np.power(2, lev - 1)) - 1
		end = int(np.power(2, lev) - 1) - 1
		for indx in range(start, end + 1):
			exp1 = end + (indx - start)
			llr1 = llrs[exp1 + 1]
			llr2 = llrs[exp1 + 1 + 2 ** (lev - 1)]
			llrs[indx] = np.sign(llr1) * np.sign(llr2) * min(abs(llr1), abs(llr2))
	
	return llrs

def my_update_bits(position, latestbit, bits):
	N = codewordLength
	if position == N - 1:
		return
	elif position < N // 2:
		bits[0][0] = latestbit
	else:
		lastlevel = (bin(position)[2:].zfill(n)).find('0') + 1
		bits[1][0] = latestbit
		for lev in range(1, lastlevel - 1):
			st = int(np.power(2, lev - 1)) - 1
			ed = int(np.power(2, lev) - 1) - 1
			for i in range(st, ed + 1):
				bits[1][ed + (i - st) + 1] = (bits[0][i] + bits[1][i]) % 2
				bits[1][ed + (i - st) + 1 + 2 ** (lev - 1)] = bits[1][i]

		lev = lastlevel - 1
		st = int(np.power(2, lev - 1)) - 1
		ed = int(np.power(2, lev) - 1) - 1
		for i in range(st, ed + 1):
			bits[0][ed + (i - st) + 1] = (bits[0][i] + bits[1][i]) % 2
			bits[0][ed + (i - st) + 1 + 2 ** (lev - 1)] = bits[1][i]
	return bits

llrs[codewordLength - 1:] = (np.array([-2.0, -2.5, -4.0, 1.0, -6.5, 6.0, 16.6, 3.5]))

for i in range(codewordLength):
	ii = bitreversed(i, n)
	llrs = my_update_llrs(ii, llrs, bits)
	llrOb.append(llrs[0])
	if polarMask[ii] == 1:
		bits = my_update_bits(ii, 1, bits)
	else:
		bits = my_update_bits(ii, 0, bits)
print(llrOb)
'''
