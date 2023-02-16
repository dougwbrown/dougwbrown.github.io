from sympy import *
import numpy as np
import functools
import random
import traceback
import ast

"""
The constants I X Y Z H S Sdg T Tdg SX SXdg SWAP are gates without args.
The functions P Rx Ry Rz R U are gates with args.
"""

def U(th,ph,la): return Matrix(
	[
		[cos(th/2), -exp(J*la)*sin(th/2)],
		[exp(J*ph)*sin(th/2), exp(J*(ph+la))*cos(th/2)]
	])

J = sqrt(-1) # use J for sqrt(-1) instead of I ...
I = Matrix([[1,0],[0,1]]) # ... because QC wants I = the identity gate
X = Matrix([[0,1],[1,0]])
Y = Matrix([[0,-J],[J,0]])
Z = Matrix([[1,0],[0,-1]])
H = U(pi/2,0,pi)
S = Matrix([[1,0],[0,J]])
Sdg = Matrix([[1,0],[0,-J]])
T = Matrix([[1,0],[0,exp(J*pi/4)]])
Tdg = Matrix([[1,0],[0,exp(-J*pi/4)]])
SX = Matrix([[1/2+J/2,1/2-J/2],[1/2-J/2,1/2+J/2]])
SXdg = Matrix([[1/2-J/2,1/2+J/2],[1/2+J/2,1/2-J/2]])
SWAP = Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])

def P(th): return U(0,0,th)
def Rx(th): return U(th,-pi/2,pi/2)
def Ry(th): return U(th,0,0)
def Rz(th): return exp(-J*th/2)*U(0,th,0)

def R(th,ph): return U(th,ph-pi/2,-ph+pi/2)

def kron(*args):
	""" Return the kronecker product of one or more matrices. """
	assert len(args) >= 1
	return Matrix(functools.reduce(np.kron,args))

_ketmap = {
	'0': [[1],[0]],
	'1': [[0],[1]],
	'+': [[1/sqrt(2)],[1/sqrt(2)]],
	'-': [[1/sqrt(2)],[-1/sqrt(2)]]
}

def ket(st):
	"""
	Return ket of arg.

	Arg is single string of 0s, 1s, +s, or -s, no whitespace.
	e.g.: ket('0+') returns |01>
	"""
	assert isinstance(st, str)
	assert set(st) <= set('01+-')
	return Matrix(kron(*list(map(lambda s: _ketmap[s],st))))

def bra(st):
	"""
	Return bra of arg.

	Arg is single string of 0s, 1s, +s, or -s, no whitespace.
	e.g.: bra('0+') returns <01|
	"""
	assert isinstance(st, str)
	assert set(st) <= set('01+-')
	return Matrix(transpose(ket(st)))

def ckt(*args):
	"""
	Build a gate from one or more gates.

	The first arg is the number of qubits.
	The rest of the args are each (gate, controls, targets), where
		gate is a gate matrix, like X or Rz(pi/2).
		controls is a list of zero or more control qubits.
		targets is a list of one or more target qubits.
	There must be a target for each qubit in the gate (e.g. 1 for X, 2 for SWAP) and they must be sequential and ascending.
	Use perm() if you want to swizzle the qubits.
	See examples of ckt in the self-test code.
	"""
	assert len(args) > 1
	assert isinstance(args[0], int)
	qubits = args[0]
	ms = [_gate(qubits, g, c, t) for g,c,t in args[1:]]
	r = Matrix(functools.reduce(lambda a,b:a*b, reversed(ms)))
	return r

def _gate(qubits, gate, controls, targets):
	assert qubits > 0
	assert isinstance(gate, Matrix)
	assert isinstance(controls, list) or isinstance(controls, tuple)
	assert isinstance(targets, list) or isinstance(targets, tuple)
	tc = set(controls)
	tt = set(targets)
	tq = set(range(qubits))
	assert tc <= tq
	assert tt <= tq
	assert tc.isdisjoint(tt)
	d,dd = gate.shape
	assert d == dd
	assert targets == list(range(targets[0], targets[0]+len(targets)))
	assert 2**len(targets) == d
	nc = 2**len(controls)
	R = []
	kb = [ket('0')*bra('0'), ket('1')*bra('1')]
	for inputState in range(nc):
		k = []
		for qubit in range(qubits):
			if qubit in targets:
				if inputState == nc-1:
					if qubit == targets[0]:
						k.append(gate)
				else:
					k.append(I)
			elif qubit in controls:
				wcb = controls.index(qubit)
				t = (inputState >> wcb) & 1
				k.append(kb[t])
			else:
				k.append(I)
		R.append(functools.reduce(kron, reversed(k)))
	return functools.reduce(lambda a,b:a+b, R) # sum() does not work here

def perm(*args):
	"""
	Rearrange the qubits.

	Args (ints) are in the new order.
	See examples in the self-test code.
	"""
	assert isinstance(args, list) or isinstance(args, tuple)
	assert all(isinstance(t,int) and t >= 0 for t in args)
	assert len(args) == len(set(args))
	n = len(args)
	n2 = 2**n
	r = [[0 for _ in range(n2)] for _ in range(n2)]
	for fr in range(n2):
		to = 0
		for i in range(n):
			b = (fr >> n-args[i]-1) & 1
			to += b << (n-i-1)
		r[to][fr] = 1
	return Matrix(r)

_varz_cache = {}

def varz(str):
	"""
	Define a symbol.

	Arg is one string, a whitespace-separated list of symbols.
	Each symbol is name:type, where
		name is the name.
		type is a,s,or g, for angle, state, or gate.
	You specify the type so symbolic representation includes appropriate constraints,
	and numeric evaluation in eqn does an appropriate substitution.
	"""
	r = []
	for thing in str.split():
		varName,varType = thing.split(':')
		_varz_cache[varName] = varType
		if varType == 'a':
			r.append(var(varName))
		elif varType == 's':
			a = var(varName + "_0")
			b = var(varName + "_1")
			k = sqrt(abs(a)**2 + abs(b)**2)
			r.append(Matrix([[a/k],[b/k]]))
		elif varType == 'g':
			a = var(varName + "_0")
			b = var(varName + "_1")
			c = var(varName + "_2")
			r.append(U(a,b,c))
		else:
			assert False
	return r

def _isMatrix(a):
	return isinstance(a, MutableDenseMatrix) or isinstance(a, ImmutableDenseMatrix)

def eqs(a,b):
	"""
	Return whether a exactly equals b.

	May fail to simplify enough to decide equality.
	All the examples in the self-test code correctly return true with just a simplify of the two args,
	except the one noted below.  rewrite(exp) worked for that case.
	You're on your own for any other simplification problems.
	"""
	assert _isMatrix(a) and _isMatrix(b)
	if simplify(a) == simplify(b):
		return True
	# this is needed for eqs(R(theta,phi), exp(-J * (theta/2) * (cos(phi)*X + sin(phi)*Y)))
	return (a-b).as_immutable().rewrite(exp).simplify().is_zero_matrix

def eqn(a,b):
	""" Return whether a==b, or close enough (1e-14), after random numeric substition for symbols."""
	assert _isMatrix(a) and _isMatrix(b)
	def rn(): return random.random()*2-1
	def ra(): return rn()*2*pi
	diff = a-b
	sublist = []
	for varName in _varz_cache:
		varType = _varz_cache[varName]
		if varType == 'a':
			sublist.append((varName, ra()))
		elif varType == 's':
			sublist.append((varName + "_0", rn()+J*rn()))
			sublist.append((varName + "_1", rn()+J*rn()))
		elif varType == 'g':
			sublist.append((varName + "_0", ra()))
			sublist.append((varName + "_1", ra()))
			sublist.append((varName + "_2", ra()))
		else:
			assert False
	diff = diff.subs(sublist)
	diff = diff.evalf()
	t = diff.shape
	assert len(t) == 2
	for a in range(t[0]):
		for b in range(t[1]):
			if abs(diff[a,b]) > 1e-14:
				return False
	return True

if __name__ == '__main__':

	# the self-test code.

	theta,phi,g1,g2,g3,psi = varz('theta:a phi:a g1:g g2:g g3:g psi:s')

	def runtest(f):
		assert f(I, ket('0')*bra('0')+ket('1')*bra('1'))
		assert f(X, ket('1')*bra('0')+ket('0')*bra('1'))
		assert f(X*X, I)
		assert f(Y, J*ket('1')*bra('0')-J*ket('0')*bra('1'))
		assert f(Z, ket('0')*bra('0')-ket('1')*bra('1'))
		assert f(H*H, I)
		assert f(S*S, Z)
		assert f(Sdg, adjoint(S))
		assert f(S*Sdg, I)
		assert f(T*T, S)
		assert f(Tdg, adjoint(T))
		assert f(T*Tdg, I)
		assert f(SX*SX, X)
		assert f(SXdg, adjoint(SX))
		assert f(SX*SXdg,I)
		assert f(SWAP, ket('00')*bra('00') + ket('01')*bra('10') + ket('10')*bra('01') + ket('11')*bra('11'))
		assert f(SWAP*SWAP, kron(I,I))
		assert f(kron(g1,g2), SWAP * kron(g2,g1) * SWAP)
		assert f(P(pi), Z)
		assert f(P(pi/2), S)
		assert f(P(pi/4), T)
		assert f(P(theta), exp(J*(theta/2))*Rz(theta))
		assert f(Rx(theta), exp(-J * theta/2 * X))
		assert f(Ry(theta), exp(-J * theta/2 * Y))
		assert f(Rz(theta), exp(-J * theta/2 * Z))
		assert f(R(theta,phi), exp(-J * (theta/2) * (cos(phi)*X + sin(phi)*Y)))
		assert f(ket('0+'), Matrix([[1/sqrt(2)],[1/sqrt(2)],[0],[0]]))
		assert f(bra('-0'), Matrix([[1/sqrt(2),0,-1/sqrt(2),0]]))
		assert f(bra('01')*ket('01'), Matrix([1]))
		assert f(ket('01')*bra('01'), Matrix([[0,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]]))
		assert f(g1, ckt(1,(g1,[],[0])))
		assert f(g1*g2, ckt(1,(g2,[],[0]),(g1,[],[0])))
		assert f(kron(g1,g2), ckt(2, (g2,[],[0]), (g1,[],[1])))
		assert f(ckt(2, (X,[0],[1])), Matrix([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]]))
		assert f(ckt(2, (X,[1],[0])), Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]))
		assert f(kron(I,I), perm(0,1))
		assert f(SWAP, perm(1,0))
		assert f(Rx(phi)*Ry(theta)*Ry(-theta)*Rx(-phi), I)
		# 1.4 4
		assert f(X, H*Z*H)
		# 2.2 2
		assert f(X*ket('0'), ket('1'))
		assert f(X*ket('1'), ket('0'))
		# "now compare them to mine"
		assert f(X*ket('-'), -ket('-'))
		CNOT = ckt(2,(X,[0],[1]))
		assert f(CNOT*ket('-0'),kron(ket('-'),ket('0')))
		assert f(CNOT*ket('-0'),ket('-0'))
		assert f(CNOT*ket('-1'),kron(X*ket('-'),ket('1')))
		assert f(CNOT*ket('-1'),kron(-ket('-'),ket('1')))
		assert f(CNOT*ket('-1'),-ket('-1'))
		assert f(CNOT*ket('-+'),1/sqrt(2)*(CNOT*ket('-0')+CNOT*ket('-1')))
		assert f(CNOT*ket('-+'),1/sqrt(2)*(ket('-0')+kron(X,I)*ket('-1')))
		assert f(CNOT*ket('-+'),1/sqrt(2)*(ket('-0')-ket('-1')))
		assert f(CNOT*ket('-+'),kron(ket('-'),(1/sqrt(2)*(ket('0')-ket('1')))))
		assert f(CNOT*ket('-+'),ket('--'))
		# section 1.Making a Controlled-Z from a CNOT
		assert f(H*X*H, Z)
		assert f(H*Z*H, X)
		assert f(ckt(2,(Z,[0],[1])),ckt(2,(H,[],[1]),(X,[0],[1]),(H,[],[1])))
		assert f(ckt(2,(Y,[0],[1])),ckt(2,(Sdg,[],[1]),(X,[0],[1]),(S,[],[1])))
		assert f(ckt(2,(H,[0],[1])),ckt(2,(Ry(pi/4),[],[1]),(X,[0],[1]),(Ry(-pi/4),[],[1])))
		# section 2.Swapping Qubits
		assert f(ckt(2,(X,[0],[1]),(X,[1],[0]))*ket('01'),ket('10'))
		assert f(kron(I,I),ckt(2,(X,[0],[1]),(X,[1],[0]),(X,[1],[0]),(X,[0],[1])))
		assert f(SWAP,ckt(2,(X,[1],[0]),(X,[0],[1]),(X,[1],[0])))
		assert f(SWAP,ckt(2,(X,[0],[1]),(X,[1],[0]),(X,[0],[1])))
		# section 3.Controlled Rotations
		assert f(ckt(2,(Ry(theta),[0],[1])),ckt(2,(Ry(theta/2),[],[1]),(X,[0],[1]),(Ry(-theta/2),[],[1]),(X,[0],[1])))
		assert f(ckt(2,(Rz(theta),[0],[1])),ckt(2,(Rz(theta/2),[],[1]),(X,[0],[1]),(Rz(-theta/2),[],[1]),(X,[0],[1])))
		# section 4.The Toffoli
		assert f(ckt(3,(X,[0,1],[2])),ckt(3,(SX,[1],[2]),(X,[0],[1]),(SXdg,[1],[2]),(X,[0],[1]),(SX,[0],[2])))
		assert f(ckt(3,(X,[0,1],[2])),ckt(3,(H,[],[2]),(X,[1],[2]),(Tdg,[],[2]),(X,[0],[2]),(T,[],[2]),(X,[1],[2]),(Tdg,[],[2]),(X,[0],[2]),(T,[],[1]),(T,[],[2]),(H,[],[2]),(SWAP,[],[1,2]),(X,[0],[2]),(T,[],[0]),(Tdg,[],[2]),(X,[0],[2]),(SWAP,[],[1,2])))
		assert f(ckt(3,(X,[0,1],[2])),ckt(3,(H,[],[2]),(X,[1],[2]),(Tdg,[],[2]),(X,[0],[2]),(T,[],[2]),(X,[1],[2]),(Tdg,[],[2]),(X,[0],[2]),(T,[],[1]),(T,[],[2]),(H,[],[2]),(X,[0],[1]),(T,[],[0]),(Tdg,[],[1]),(X,[0],[1])))
		assert f(ckt(3,(H,[0],[2]),(Z,[1],[2]),(H,[0],[2]))*kron(psi,ket('00')),kron(psi,ket('00')))
		assert f(ckt(3,(H,[0],[2]),(Z,[1],[2]),(H,[0],[2]))*kron(psi,ket('11')),kron(X*psi,ket('11')))
		assert f(ckt(3,(H,[0],[2]),(Z,[1],[2]),(H,[0],[2]))*kron(psi,ket('01')),kron(psi,ket('01')))
		assert f(ckt(3,(H,[0],[2]),(Z,[1],[2]),(H,[0],[2]))*kron(psi,ket('10')),kron(Z*psi,ket('10')))
		# Dirac Notation
		assert f(ket('0+'), Matrix([[1/sqrt(2)],[1/sqrt(2)],[0],[0]]))
		assert f(bra('-0'), Matrix([[1/sqrt(2),0,-1/sqrt(2),0]]))
		assert f(bra('01')*ket('01'), Matrix([1]))
		assert f(ket('01')*bra('01'), Matrix([[0,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]]))
		# Gate Generation: ckt (is qcc in online calc.js)
		CNOT1 = ckt(2,(X,[0],[1]))
		CNOT2 = ckt(2,(X,[1],[0]))
		assert f(ckt(2,(H,[],[0]),(H,[],[1]),(CNOT1,[],[0,1]),(H,[],[0]),(H,[],[1])),CNOT2)
		assert f(ckt(2,(H,[],[0]),(H,[],[1]),(CNOT2,[],[0,1]),(H,[],[0]),(H,[],[1])),CNOT1)
		assert f(ckt(3,(g1,[],[0]),(g2,[],[1]),(g3,[],[2])),kron(g3,g2,g1))
		assert f(ckt(2,(Rx(pi/2),[],[0])),kron(I,Rx(pi/2)))
		assert f(ckt(2,(SWAP,[],[0,1]),(g1,[],[0]),(g2,[],[1]),(SWAP,[],[0,1])),kron(g1,g2))
		# Gate Generation: perm
		assert f(perm(0,1),kron(I,I))
		assert f(perm(1,0),SWAP)
		assert f(perm(0,1,2),kron(I,I,I))
		assert f(perm(0,2,1),kron(I,SWAP))
		assert f(perm(1,0,2),kron(SWAP,I))
		assert f(perm(1,2,0),kron(I,SWAP)*kron(SWAP,I))
		assert f(perm(2,0,1),kron(SWAP,I)*kron(I,SWAP))
		assert f(perm(2,1,0),kron(I,SWAP)*kron(SWAP,I)*kron(I,SWAP))
		assert f(kron(I,I,I),perm(0,1,2))
		assert f(kron(I,I,I),perm(2,1,0)*perm(2,1,0))
		assert f(kron(I,I,I),perm(1,2,0)*perm(1,2,0)*perm(1,2,0))
		# Random things
		assert f(I, Rx(theta)*Ry(phi)*Ry(-phi)*Rx(-theta))
		assert f(abs(bra('0')*psi)**2+abs(bra('1')*psi)**2,Matrix([1]))

	runtest(eqn)
	print('numeric tests ok')
	runtest(eqs)
	print('symbolic tests ok')
