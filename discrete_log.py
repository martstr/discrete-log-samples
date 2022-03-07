#!/usr/bin/env python -i

from math import sqrt, ceil
from time import perf_counter

def bruteforce(generator, group_element, modulus):
	# Solves a discrete log problem using a bruteforce approach, i.e.
	# finds a number exponent such that
	#
	#	generator^exponent = group_element (mod modulus)

	tic = perf_counter()

	value = 1
	for exponent in range(1, modulus):
		value = value * generator % modulus
		if value == group_element:
			toc = perf_counter()
			print(f"Duration: {toc - tic:0.4f} seconds")
			return exponent
	toc = perf_counter()
	print(f"Duration: {toc - tic:0.4f} seconds")

	return 0

def bsgs(generator, group_element, modulus, order = None, standalone = True):
	# Solves a discrete log problem using Shank's algorithm.
	# (also known as baby-step, giant-step)

	if standalone:
		tic = perf_counter()
	if order is None:
		order = modulus-1
	m = ceil(sqrt(order))

	L1 = []
	for j in range(0, m):
		L1.append((j, pow(generator, (m*j), modulus)))
	L1 = sorted(L1, key=lambda x: abs(x[1]))

	L2 = []
	for i in range(0, m):
		L2.append((i, (group_element * pow(generator, -i, modulus)) % modulus))
	L2 = sorted(L2, key=lambda x: abs(x[1]))

	found = False
	i = 0
	j = 0

	while (not found) and (i < m) and (j < m):
		if L1[j][1] == L2[i][1]:
			found = True
		elif abs(L1[j][1]) > abs(L2[i][1]):
			i = i + 1
		else:
			j = j + 1

	if found:
		if standalone:
			toc = perf_counter()
			print(f"Duration: {toc - tic:0.4f} seconds")
		return m*L1[j][0]+L2[i][0] % order

	if standalone:
		toc = perf_counter()
		print(f"Duration: {toc - tic:0.4f} seconds")

	return 0

def pi1(x, order, modulus):
	# The projection to the other subgroup than the one with the given
	# order
	return pow(abs(x), order, modulus)

def pohlig_hellman(generator, group_element, modulus, Q, C, verbose=True):
	# For each factor, push the element into a subgroup by multiplying
	# by all the other factors.

	tic = perf_counter()
	D = []
	for q, c in zip(Q, C):
		if verbose:
			print("Pohlig-Hellman in the subgroup "+str(q)+"^"+str(c))
		m = 1
		for a, b in zip(Q, C):
			if q != a:
				m = m * pow(a,b)
		h = pi1(generator, m, modulus)
		y = pi1(group_element, m, modulus)
		d = pohlig_hellman_prime(q, c, h, y, modulus, verbose)
		D.append(d)

	n = [pow(q, c) for q, c in zip(Q, C)]

	res = crt(D, n)
	toc = perf_counter()
	print(f"Duration: {toc - tic:0.4f} seconds")

	return res

def pohlig_hellman_prime(q, e, g, h, modulus, verbose = True):
	# Here, we push the elements further down into smaller subgroups
	# of prime order p instead of p to some power.
	# Finally, we use Shank's algorithm to solve the problem there.

	x = 0
	inv_g = modinv(g, modulus)
	gamma = pow(g, pow(q, e-1), modulus)

	for k in range(e):
		if verbose:
			print("    BSGS in the subgroup "+str(q)+"^"+str(k))
		y = (pow(inv_g, x)*h) % modulus
		hk = pow(y, pow(q, e-1-k), modulus)
		d = bsgs(gamma, hk, modulus, q, False)
		if verbose:
			print("    d such that "+str(gamma)+"^d = "+str(hk)+":", d)
		x = x + pow(q, k)*d

	return x

def crt(a, n):
	# Chinese remainder theorem
	# Both a and n are lists.
	l = len(a)

	m = 1
	for p in n:
		m = m*p

	sum = 0
	for ai, p in zip(a, n):
		si = modinv(m//p, p)
		sum = (sum + ai*si*m//p) % m

	return sum

def egcd(a, b):
	# Extended Euclidean algorithm. Returns g, x, y such that
	#
	#	g = gcd(a, b)
	#	ax + by = gcd(a, b)

    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
	# Returns the inverse of a modulo m, i.e. the number b such that
	#
	#	ab = 1 (mod m)

    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse for '+str(a)+' modulo '+str(m)+
        	' does not exist')
    else:
        return x % m

print("""
Usage: python -i discrete_log.py
   or: ./discrete_log.py (if you have set the executable flag on a *nix box)

The file implements some example algorithms to break the discrete log
	o bruteforce
	o Shank's algorithm (Baby-step, giant-step, abbr. bsgs)
	o Pohlig-Hellman

	bruteforce(generator, group_element, modulus):
		generator is the base element
		group_element is the element of which we want to find the dlog
		modulus is the (main) prime we're working with

	bsgs(generator, group_element, modulus, order):
		generator: the base element
		group_element: the element of which we want to find the dlog
		modulus: the (main) prime we're working with
		order: the parameter that defines the upper bound for the lists

	pohlig_hellman(generator, group_element, modulus, Q, C):
		generator: the base element
		group_element: the element of which we want to find the dlog
		modulus: the (main) prime we're working with
		Q is a list of prime factors of modulus-1
		C is a list of prime multiplicities for Q

Example

	modulus = 20201227
	generator = 7
	group_element = 75735546
	Q = [2, 3, 29, 116099]
	C = [1, 1, 1, 1]

	# modulus-1 = 20201226 = 2 * 3 * 29 * 116099
""")
modulus = 20201227
generator = 7
group_element = 7573546
Q = [2, 3, 29, 116099]
C = [1, 1, 1, 1]

large_modulus = 3467882422082372663
large_group_element = 1400024179417283732
large_generator = 5
large_Q = [2, 1196813, 1448798777287]
large_C = [1, 1, 1]
