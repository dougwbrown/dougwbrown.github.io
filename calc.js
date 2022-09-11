/*
	(C) Copyright Doug Brown 2022

	This code is licensed under the Apache License, Version 2.0. You may
	obtain a copy of this license at http://www.apache.org/licenses/LICENSE-2.0.

	Any modifications or derivative works of this code must retain this
	copyright notice, and modified files need to carry a notice indicating
	that they have been altered from the originals.
*/

// Setting CatchButton true adds a radio button that selects here (the results) or the console gets error messages.
// It's not useful unless you're familiar with the source code, so the choice resides here in the source code.
// If CatchButton is false, then there is no button to redirect error messages, and they always go to results.
var CatchButton = false

// Scope contains constants and functions that go into context.scope, and from there to evaluate(context.scope)
var Scope = {}

// build gates
{
	let i = math.complex(0,1)
	let mi = math.complex(0,-1)
	let mm = (a,b) =>  math.multiply(a,b)
	let mm3 = (a,b,c) => math.multiply(a, math.multiply(b,c))
	let sqrt = math.sqrt
	let exp = math.exp
	let mexp = (x) => math.subtract(0, exp(x))
	let sin = math.sin
	let cos = math.cos
	let pi = math.pi
	let ohpio2 = math.add(1/2,math.divide(math.complex(0,1),2))
	let ohmio2 = math.add(1/2,math.divide(math.complex(0,-1),2))

	Scope.I = [[1,0],[0,1]]
	Scope.X = [[0,1],[1,0]]
	Scope.Y = [[0,mi],[i,0]]
	Scope.Z = [[1,0],[0,-1]]
	Scope.H = [[1/sqrt(2),1/sqrt(2)],[1/sqrt(2),-1/sqrt(2)]]
	Scope.S = [[1,0],[0,i]]
	Scope.Sdg = [[1,0],[0,mi]]
	Scope.T = [[1,0],[0,exp(mm(i,pi/4))]]
	Scope.Tdg = [[1,0],[0,exp(mm(mi,pi/4))]]
	Scope.SX = [[ohpio2,ohmio2],[ohmio2,ohpio2]]
	Scope.SXdg = [[ohmio2,ohpio2],[ohpio2,ohmio2]]
	Scope.SWAP = [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]

	Scope.P = (th) => [[1,0],[0,exp(mm(i,th))]]
	Scope.Rx = (th) => [[cos(th/2),mm(mi,sin(th/2))],[mm(mi,sin(th/2)),cos(th/2)]]
	Scope.Ry = (th) => [[cos(th/2),-sin(th/2)],[sin(th/2),cos(th/2)]]
	Scope.Rz = (th) => [[exp(mm(mi,th/2)),0],[0,exp(mm(i,th/2))]]

	Scope.R = (th,ph) => [[cos(th/2), mm3(mi,exp(mm(mi,ph)),sin(th/2))],[mm3(mi,exp(mm(i,ph)),sin(th/2)), cos(th/2)]]

	Scope.U = (th,ph,la) => [[cos(th/2), mm(mexp(mm(i,la)),sin(th/2))],[mm(exp(mm(i,ph)),sin(th/2)),mm(exp(mm(i,(math.add(ph,la)))),cos(th/2))]]

}

// external use
function createCalc(div, text) {
	// div can be div id or div
	if (typeof(div) == 'string') div = document.getElementById(div)
	$(document).ready(function() {
		createCalc1(div, text)
	})
}

// internal use
function createCalc1(div, text) {
	let trio = document.createElement("div")
	trio.style = "border-style:hidden hidden hidden solid; margin-top: 10px"
	div.appendChild(trio)
	let expr = document.createElement("textarea")
	expr.style = "margin-left: 5px"
	let buttons = document.createElement("div")
	let result = document.createElement("div")
	result.style = "margin-left: 5px; margin-top: 5px"
	trio.appendChild(expr)
	trio.appendChild(buttons)
	trio.appendChild(result)
	expr.type = 'text'
	expr.cols = '80'
	expr.rows = '10'
	expr.value = text
	expr.placeholder = "enter expression here"
	let context = {div:div, expr:expr, buttons:buttons, result:result, options:{}}
	context.scope = Object.assign({}, Scope) // shallow copy sufficient
	expr.oninput = () => updateScreen(context)
	expr.style.display = 'none'
	result.addEventListener('click', function (event) {
		toggleDisplay(context)
		if (expr.style.display == 'block') expr.focus()
	})
	buttons.style.display = 'none'
	if (CatchButton) addRadio(context, ['Catch Errors:', 'here', 'console'], context, 'catching')
	addRadio(context, ['Parenthesis:', 'keep', 'auto', 'all'], context.options, 'parenthesis')
	addRadio(context, ['Implicit:', 'hide', 'show'], context.options, 'implicit')
	addButton(context, 'Re-evaluate', (event) => updateScreen(context))
	addButton(context, 'User Manual', (event) => window.open("calc.doc.html", '_blank') )
	// addButton(context, 'User Manual', (event) => window.open("calc.doc.html", '_blank') )
	expr.focus()
	updateScreen(context)
}

function addButton(context, label, action) {
	let t = document.createElement('button')
	t.textContent = label
	t.addEventListener('click', action)
	context.buttons.appendChild(t)
}

function addRadio(context, labels, flag, flagele) {
	context.buttons.append(document.createTextNode(labels[0]))
	let ts = []
	for (let i = 1; i < labels.length; i++) {
		ts.push(document.createElement('input'))
		let t = ts[ts.length-1]
		t.type = 'radio'
		if (i == 1) t.checked = true
		context.buttons.append(t)
		context.buttons.append(document.createTextNode(labels[i]))
	}
	flag[flagele] = labels[1]
	let spacer = document.createElement('span')
	spacer.style = "display:inline-block;width:15px"
	context.buttons.append(spacer)
	for (let listener = 0; listener < ts.length; listener++) {
		let t = ts[listener]
		t.addEventListener('click', (event) => {
			for (let other = 0; other < ts.length; other++) {
				if (other != listener) ts[other].checked = false
			}
			flag[flagele] = labels[listener+1]
			updateScreen(context)
		})
	}
}

function toggleDisplay(context) {
	let t = (context.expr.style.display == 'block') ? 'none' : 'block'
	context.expr.style.display = t
	context.buttons.style.display = t
}

function perm(arg) {
	if (arguments.length != 1) throw "wrong number of arguments for function perm"
	if (typeof(arg) != 'string') throw "function perm takes a string argument: " + arg
	if (!arg.match(/^\d+$/)) throw "argument for function perm must be digits: " + arg
	let theperm = arg.split('').map(bit => parseInt(bit))
	let n = theperm.length
	let n2 = 2**n
	let r = math.zeros(n2,n2)._data
	for (let fr = 0; fr < n2; fr++) {
		let to = 0
		for (let i = 0; i < n; i++) {
			let b = (fr >> n-theperm[i]-1) & 1
			to += b << (n-i-1)
		}
		r[to][fr] = 1
	}
	return r
}

function qc(args0, math, scope) {
	if (args0.length < 1 || args0.length % 2 == 0) throw "wrong number of arguments for function qc"
	args = args0.map(arg => arg.evaluate(scope))
	let qubits = args[0]
	args0[0].qcHack= qubits
	let r = math.identity(2**qubits)
	for (let i = 1; i < args.length; i+=2) {
		let U = args[i]
		if (scalar(U)) throw "bad gate in qc: " + U
		if (!args[i+1].match(/^(\d+>)?\d+$/)) throw "bad syntax for control/targets: " + args[i+1]
		let ct = args[i+1].split('>')
		let controls
		let targets
		if (ct.length == 1) {
			controls = []
			targets = ct[0].split('').map(d=>parseInt(d))
		} else if (ct.length == 2) {
			controls = ct[0].split('').map(d=>parseInt(d))
			targets = ct[1].split('').map(d=>parseInt(d))
		} else {
			throw "bad syntax for controls/targets: " + args[i+1]
		}
		while (2**targets.length != ('_data' in U ? U._data.length : U.length)) {
			throw "wrong number of target qubits: " + args[i+1]
		}
		for (let j = 1; j < targets.length; j++) {
			if (targets[j] != 1+targets[j-1]) throw "target qubits must be consequtive and ascending: " + targets
		}
		args0[i+1].qcHack = {controls:controls, targets:targets}
		qcCheck(targets,controls,qubits)
		r = math.multiply(r, ncGate(qubits, U, controls, targets))
	}
	return r
}

// I don't know how to evaluate the qubits and controls/targets args in qcLatex(),
// because I don't have access to context.scope,
// so I save the values during qc() in node.args[whatever].qcHack

function qcLatex(node, options) {
	let r = []
	let qubits = node.args[0].qcHack
	for (let i = 1; i < node.args.length; i+=2) {
		let Ud = node.args[i].toTex(options)
		let targets = node.args[i+1].qcHack.targets
		let controls = node.args[i+1].qcHack.controls
		let col = [] // even though it's going in as a row for now
		for (let row = 0; row < qubits; row++) {
			if (targets.includes(row)) {
				col.push(Ud)
			} else if (controls.includes(row)) {
				col.push('c')
			} else {
				col.push('-')
			}
		}
		r.push(col)
	}
	return pretty(math.transpose(r), mt='Bmatrix')
}

function qcCheck(targets, controls, qubits) {
	targets.forEach(function(t) {if (controls.includes(t)) throw "targets and controls cannot overlap"})
	controls.forEach(function(c) {if (targets.includes(c)) throw "targets and controls cannot overlap"})
	targets.forEach(function(t) {if (t >= qubits) throw "targets and controls cannot exceed the qubits"})
	controls.forEach(function(c) {if (c >= qubits) throw "targets and controls cannot exceed the qubits"})
}

function eq() {
	if (arguments.length < 2) throw "function eq must have at least 2 arguments"
	let r = true
	for (let i = 1; i < arguments.length; i++) {
		let b = eq1(arguments[0], arguments[i])
		r &&= b
	}
	return r
}

function eq1(a,b) {
	if (scalar(a) || scalar(b)) {
		// compareNatural takes care of many checks
		return math.compareNatural(a,b) == 0
	}
	// compareNatural gets confused about array vs object
	let r = true
	let t = math.matrix(math.subtract(a,b))._data
	t.forEach(function (value, index, matrix) {
		value.forEach(function(v) {
			if (!isclose(v,0)) { r = false }
		})
	})
	return r
}

let ketmap = {
	'0': [[1],[0]],
	'1': [[0],[1]],
	'+': [[1/math.sqrt(2)],[1/math.sqrt(2)]],
	'-': [[1/math.sqrt(2)],[-1/math.sqrt(2)]]
}

function checkDiracArgs(args, fun, n) {
	if (args.length != n) throw `wrong number of args to function ${fun}`
	for (let i = 0; i < n; i++) {
		if (typeof(args[i]) != 'string' || !args[i].match(/^[01\+\-]+$/)) throw `${n>1?'each':''} argument of function ${fun} must be a string with only 0's, 1's, +'s, and -'s in it: ` + args[i]
	}
}

function ket(s) {
	checkDiracArgs(arguments, 'ket', 1)
	let r = ketmap[s[0]]
	for (let i = 1; i < s.length; i++) {
		r = math.kron(r, ketmap[s[i]])
	}
	return r
}

function bra(s) {
	checkDiracArgs(arguments, 'bra', 1)
	return math.transpose(ket(s))
}

function braket(s0,s1) {
	checkDiracArgs(arguments, 'braket', 2)
	if (typeof(s0) != "string" || typeof(s1) != "string") throw "function braket wants string arguments"
	return math.multiply(bra(s0), ket(s1))
}

function ketbra(s0,s1) {
	checkDiracArgs(arguments, 'ketbra', 2)
	if (typeof(s0) != "string" || typeof(s1) != "string") throw "function ketbra wants string arguments"
	return math.multiply(ket(s0), bra(s1))
}

function ncGate(qubits, gate, controls, targets) {
	let nc = 2**controls.length
	let I = Scope['I']
	let R = 0
	for (let i = 0; i < nc; i++) {
		let k = null
		for (let j = qubits-1; j >= 0; j--) {
			let kk = null
			if (targets.includes(j)) {
				if (i == nc-1) {
					kk = (j == targets[0]) ? gate : null
				} else {
					kk = I
				}
			} else if (controls.includes(j)) {
				let wcb = controls.indexOf(j)
				let t = (i >> wcb) & 1 
				kk = t ? ketbra('1','1') : ketbra('0','0')
			} else {
				kk = I
			}
			if (kk != null) {
				k = k == null ? kk : math.kron(k, kk)
			}
		}
		R = math.add(R, k)
	}
	return R
}

function scalar(x) {
	let t = typeof(x)
	return t == 'number' || t == 'boolean' || t == 'string' || x.isComplex
}

function isclose(a,b) {
	const eps = 1e-14
	return math.abs(a-b) < eps
}

function mj(tex) {
	return MathJax.tex2svg(tex, {em: 16, ex: 6, display: false});
}

function updateScreen(context) {
	if (!CatchButton || context.catching == 'here') {
		try {
			gogo(context)
		} catch (err) {
			context.result.innerHTML = `<span style="color: red;"> ${err.toString()} </span>`
		}
	} else {
		gogo(context)
	}
}

function randomState() {
	if (arguments.length > 0) throw "function randomState does not take any arguments"
	// return random number between -1 and 1
	function rn() {
		return (Math.random() * 2) - 1 
	}
	// return random state
	function rq() {
		let j1 = math.complex(0,1)
		let t1 = math.complex(rn(), rn())
		let t2 = math.complex(rn(), rn())
		let k = math.sqrt(1/(math.abs(t1)**2 + math.abs(t2)**2))
		return [[math.multiply(t1,k)],[math.multiply(t2,k)]]
	}
	return rq()
}

function randomAngle() {
	if (arguments.length > 0) throw "function randomAngle does not take any arguments"
	let hi = 2*math.pi
	let low = -2*math.pi
	return Math.random() * (hi-low) + low
}

function atToDotStarToKron() { 
	let r = arguments[0]
	for (let i = 1; i < arguments.length; i++) {
		r = math.kron(r, arguments[i])
	}
	return r
}

const customFunctions = {
	ket: ket,
	bra: bra,
	braket: braket,
	ketbra: ketbra,
	eq: eq,
	qc: qc,
	randomAngle: randomAngle,
	randomState: randomState,
	atToDotStarToKron: atToDotStarToKron,
	perm: perm,
	assert: function () {
		if (eq.apply({}, arguments)) {
			return true
		} else {
			throw "FAIL"
		}
	}
}

qc.rawArgs = true

customFunctions.ket.toTex = '\\left| ${args[0]} \\right\\rangle'
customFunctions.bra.toTex = '\\left\\langle ${args[0]} \\right|'
customFunctions.ketbra.toTex = '\\left| ${args[0]} \\right\\rangle\\left\\langle ${args[1]} \\right|'
customFunctions.braket.toTex = '\\left\\langle ${args[0]} | ${args[1]} \\right\\rangle'

customFunctions.qc.toTex = function(node,options) {
	return qcLatex(node, options)
}

customFunctions.atToDotStarToKron.toTex = function(node,options) {
	return node.args.map(t=>t.toTex(options)).join('\\otimes ')
}

customFunctions.eq.toTex = function(node,options) {
	let t = node.args.map(t=>t.toTex(options))
	let r = ['', t[0], '='+t[1], '']
	t = t.splice(2)
	while (t.length > 0) {
		r = r.concat(['', '', '='+t[0],''])
		t = t.splice(1)
	}
	return r
}

customFunctions.assert.toTex = customFunctions.eq.toTex

math.import(customFunctions)

function lhsLaTeX(node, options) {
	if (node.type == 'ConstantNode') { return node.value }  // get rid of quotes around a string
	if (node.type == 'SymbolNode') {
		switch (node.name) {
			case 'Sdg': return ' S^{\\dagger} '
			case 'Tdg': return ' T^{\\dagger} '
			case 'SXdg': return ' SX^{\\dagger} '
		}
	}
	if (node.type == 'AssignmentNode') {
		let object = node.object.toTex(options)
		let r = node.value.toTex(options)
		if (typeof(value) == 'object') {
			r[0] = object
		} else {
			r = [object, '', r, '']
		}
		r[0] += ':='
		return r
	}
}

function rhsLaTeX(node, options) {
	if (node.type != 'ConstantNode') throw "internal error"
	return pretty(node.value)
}

function gogo(context) {
	context.result.innerHTML = ''
	let expr = context.expr.value
	if (expr.trim() == '') throw "empty expression"
	if (expr.indexOf('.*') >= 0) throw ".* is being used for something else"
	expr = expr.replace(/@/g, '.*')
	expr = '\n'+expr // force BlockNode
	let parsed = math.parse(expr)
	parsed = fixDotStar(parsed)
	let compiled = parsed.compile()
	let results = compiled.evaluate(context.scope).entries
	parsed.blocks = parsed.blocks.filter(block => block.visible)
	let mjs = []
	// there are four columns
	// col 0  = var :=
	// col 3  = --> result
	// col 1&2 are for eq()
	// normal expressions would like to go into a merged col 1&2,
	// but I couldn't figure out how to merge columns, so I put them in col 2
	for (let i = 0; i < results.length; i++) {
		let node = new math.ConstantNode(results[i])
		let lhs = parsed.blocks[i].node.toTex(Object.assign({handler: lhsLaTeX}, context.options))
		let rhs = node.toTex(Object.assign({handler: rhsLaTeX}, context.options))
		if (typeof(lhs) != 'object') {
			lhs = ['','',lhs,'']
		}
		lhs[lhs.length-1] = '\\longrightarrow ' + rhs
		while (lhs.length > 0) {
			mjs.push(lhs.slice(0,4).join('&'))
			lhs = lhs.slice(4)
		}
	}
	let mjin = ' \\begin{array}{rrll} '
	mjin += mjs.join(' \\\\ ')
	mjin += ' \\end{array} '
	context.result.appendChild(mj(mjin))
}

function fixDotStar(node) {
	return node.transform(function (node, path, parent) {
		if (node.type == 'OperatorNode' && node.op == '.*') {
			return new math.FunctionNode(new math.SymbolNode('atToDotStarToKron'), node.args.map(fixDotStar))
		} else {
			return node
		}
	})
}

function pretty(x,mt='bmatrix') {
	if (typeof(x) == 'undefined') {
		return 'undefined '
	}
	if (typeof(x) == 'boolean') {
		return x ? 'true ' : 'false '
	}
	if (typeof(x) == 'string') {
		return x
	}
	if (scalar(x)) {
		return num_to_latex(x)
	}
	if (math.size(x).length == 1) { // 1D
		x = x.map(t => typeof(t) == 'string' ? t : num_to_latex(t))
		return x.join('&')
	}
	let r = '\\begin{' + mt + '}'
	x = math.matrix(x)._data
	for (let i = 0; i < x.length; i++) {
		r += pretty(x[i])
		r += '\\\\'
	}
	r += '\\end{' + mt + '}'
	return r
}
