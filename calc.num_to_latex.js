/*
	This code came from Qiskit.

	(C) Copyright IBM 2017, 2020.

	This code is licensed under the Apache License, Version 2.0. You may
	obtain a copy of this license at http://www.apache.org/licenses/LICENSE-2.0.

	Any modifications or derivative works of this code must retain this
	copyright notice, and modified files need to carry a notice indicating
	that they have been altered from the originals.
	
	num_to_latex and proc_value were converted from python to javascript
	and further changes were made.
*/


function num_to_latex(num) {
	/*
		Takes a complex number as input and returns a latex representation

			Args:
				num (numerical): The number to be converted to latex.

			Returns:
				str: Latex representation of num

		Result is combination of maximum 4 strings in the form:
			{common_facstring} ( {realstring} {operation} {imagstring}i )
		common_facstring: A common factor between the real and imaginary part
		realstring: The real part (inc. a negative sign if applicable)
		operation: The operation between the real and imaginary parts ('+' or '-')
		imagstring: Absolute value of the imaginary parts (i.e. not inc. any negative sign).
		This function computes each of these strings and combines appropriately.
	*/

	let r = num.isComplex ? num.re : num
	let i = num.isComplex ? num.im : 0
    let common_factor = null

    // try to factor out common terms in imaginary numbers
    if (isclose(math.abs(r), math.abs(i)) && ! isclose(r, 0) && ! isclose(i, 0)) {
        common_factor = math.abs(r)
        r = r/common_factor
        i = i/common_factor
	}

    common_terms = [
        [math.sqrt(2), ' \\sqrt{2} '],
        [1/math.sqrt(2), ' \\tfrac{1}{\\sqrt{2}} '],
        [1/math.sqrt(3), ' \\tfrac{1}{\\sqrt{3}} '],
        [math.sqrt(2/3), ' \\sqrt{\\tfrac{2}{3}} '],
        [math.sqrt(3/4), ' \\sqrt{\\tfrac{3}{4}} '],
        [1/math.sqrt(8), ' \\tfrac{1}{\\sqrt{8}} '],
		[math.pi, ' \\pi '],
		[math.e, ' e '],
		[math.pi/2, ' \\frac{\\pi}{2} '],
		[math.pi/3, ' \\frac{\\pi}{3} '],
		[2*math.pi/3, ' \\frac{2\\pi}{3} '],
		[math.pi/4, ' \\frac{\\pi}{4} '],
		[3*math.pi/4, ' \\frac{3\\pi}{4} '],
    ]

    // Get string (or None) for common factor between real and imag
	let common_facstring = common_factor == null ? null : proc_value(common_factor)

    // Get string for real part
    realstring = proc_value(r)

    // Get string for both imaginary part and operation between real and imaginary parts
	let operation = i > 0 ? "+" : "-"
	let imagstring = proc_value(i > 0 ? i : -i)
    if (imagstring == "1") {
        imagstring = ""  // Don't want to return '1i', just 'i'
	}

    // Now combine the strings appropriately:
    if (imagstring == "0") {
        return realstring  // realstring already contains the negative sign (if needed)
	}
    if (realstring == "0") {
        // imagstring needs the negative sign adding
		return (operation == '-' ? '-' : '') + imagstring + 'i'
	}
    if (common_facstring != null) {
		return common_facstring + '(' + realstring + operation + imagstring + 'i)'
	} else {
		return realstring + operation + imagstring + 'i'
	}
}

function proc_value(val) {
	// This function converts a real value to a latex string
	// First, see if val is close to an integer:
	val_mod = val % 1
	if (isclose(val_mod, 0) || isclose(val_mod, 1)) {
		// If so, return that integer
		return math.round(val).toString() // return str(int(np.round(val)))
	}
	// Otherwise, see if it matches one of the common terms
	for (let i = 0; i < common_terms.length; i++) {
		if (isclose(math.abs(val), common_terms[i][0])) {
			let latex_str = common_terms[i][1]
			return val > 0 ? latex_str : '-' + latex_str
		}
	}
	// try to factorise val nicely
	let t = math.fraction(val)
	let num = t.n
	let denom = t.d
	if (num + denom < 20) {
		if (denom == 1) {
			return num
		}
		// If fraction is 'nice' 
		return (val > 0 ? '' : '-') + ` \\tfrac{ ${math.abs(num)} }{ ${math.abs(denom)} }`
	} else {
		// Failing everything else, return val as a decimal
		return val.toFixed(5)
	}
}
