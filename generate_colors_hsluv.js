var colors = (function () {
	function cca(s, index) {
		var x = s.charCodeAt(index);
		return x == x ? x : undefined;
	}

	function HxSubStr(s, pos, len) {
		if (pos != null && pos != 0 && len != null && len < 0) return "";

		if (len == null) len = s.length;

		if (pos < 0) {
			pos = s.length + pos;
			if (pos < 0) pos = 0;
		} else if (len < 0) {
			len = s.length + len - pos;
		}

		return s.substr(pos, len);
	}

	function parseInt(x) {
		var v = parseInt(x, 10);
		if (v == 0 && (cca(x, 1) == 120 || cca(x, 1) == 88)) v = parseInt(x);
		if (isNaN(v)) return null;
		return v;
	}

	function StringHex(n, digits) {
		var s = "";
		var hexChars = "0123456789ABCDEF";
		do {
			s = hexChars.charAt(n & 15) + s;
			n >>>= 4;
		} while (n > 0);
		if (digits != null) while (s.length < digits) s = "0" + s;
		return s;
	}

	var Geometry = {
		intersectLineLine: function (a, b) {
			var x = (a.intercept - b.intercept) / (b.slope - a.slope);
			var y = a.slope * x + a.intercept;
			return { x: x, y: y };
		},

		distanceFromOrigin: function (point) {
			return Math.sqrt(Math.pow(point.x, 2) + Math.pow(point.y, 2));
		},

		distanceLineFromOrigin: function (line) {
			return Math.abs(line.intercept) / Math.sqrt(Math.pow(line.slope, 2) + 1);
		},

		perpendicularThroughPoint: function (line, point) {
			var slope = -1 / line.slope;
			var intercept = point.y - slope * point.x;
			return { slope: slope, intercept: intercept };
		},

		angleFromOrigin: function (point) {
			return Math.atan2(point.y, point.x);
		},

		normalizeAngle: function (angle) {
			var m = 2 * Math.PI;
			return (angle % m + m) % m;
		},

		lengthOfRayUntilIntersect: function (theta, line) {
			return line.intercept / (Math.sin(theta) - line.slope * Math.cos(theta));
		},

	};

	var Hsluv = {
		getBounds: function (L) {
			var result = [];
			var sub1 = Math.pow(L + 16, 3) / 1560896;
			var sub2;
			if (sub1 > Hsluv.epsilon) {
				sub2 = sub1;
			} else {
				sub2 = L / Hsluv.kappa;
			}
			var _g = 0;
			while (_g < 3) {
				var c = _g++;
				var m1 = Hsluv.m[c][0];
				var m2 = Hsluv.m[c][1];
				var m3 = Hsluv.m[c][2];
				var _g1 = 0;
				while (_g1 < 2) {
					var t = _g1++;
					var top1 = (284517 * m1 - 94839 * m3) * sub2;
					var top2 = (838422 * m3 + 769860 * m2 + 731718 * m1) * L * sub2 - 769860 * t * L;
					var bottom = (632260 * m3 - 126452 * m2) * sub2 + 126452 * t;
					result.push({ slope: top1 / bottom, intercept: top2 / bottom });
				}
			}
			return result;
		},

		maxSafeChromaForL: function (L) {
			var bounds = Hsluv.getBounds(L);
			var min = 1.7976931348623157e+308;
			var _g = 0;
			while (_g < 2) {
				var i = _g++;
				var length = Geometry.distanceLineFromOrigin(bounds[i]);
				min = Math.min(min, length);
			}
			return min;
		},

		maxChromaForLH: function (L, H) {
			var hrad = H / 360 * Math.PI * 2;
			var bounds = Hsluv.getBounds(L);
			var min = 1.7976931348623157e+308;
			var _g = 0;
			while (_g < bounds.length) {
				var bound = bounds[_g];
				++_g;
				var length = Geometry.lengthOfRayUntilIntersect(hrad, bound);
				if (length >= 0) min = Math.min(min, length);
			}
			return min;
		},

		dotProduct: function (a, b) {
			var sum = 0;
			var _g1 = 0;
			var _g = a.length;
			while (_g1 < _g) {
				var i = _g1++;
				sum += a[i] * b[i];
			}
			return sum;
		},
		fromLinear: function (c) {
			if (c <= 0.0031308) return 12.92 * c; else return 1.055 * Math.pow(c, 0.416666666666666685) - 0.055;
		},

		toLinear: function (c) {
			if (c > 0.04045) return Math.pow((c + 0.055) / 1.055, 2.4); else return c / 12.92;
		},

		xyzToRgb: function (tuple) {
			return [Hsluv.fromLinear(Hsluv.dotProduct(Hsluv.m[0], tuple)), Hsluv.fromLinear(Hsluv.dotProduct(Hsluv.m[1], tuple)), Hsluv.fromLinear(Hsluv.dotProduct(Hsluv.m[2], tuple))];
		},

		rgbToXyz: function (tuple) {
			var rgbl = [Hsluv.toLinear(tuple[0]), Hsluv.toLinear(tuple[1]), Hsluv.toLinear(tuple[2])];
			return [Hsluv.dotProduct(Hsluv.minv[0], rgbl), Hsluv.dotProduct(Hsluv.minv[1], rgbl), Hsluv.dotProduct(Hsluv.minv[2], rgbl)];
		},

		yToL: function (Y) {
			if (Y <= Hsluv.epsilon) return Y / Hsluv.refY * Hsluv.kappa; else return 116 * Math.pow(Y / Hsluv.refY, 0.333333333333333315) - 16;
		},

		lToY: function (L) {
			if (L <= 8) return Hsluv.refY * L / Hsluv.kappa; else return Hsluv.refY * Math.pow((L + 16) / 116, 3);
		},
		xyzToLuv: function (tuple) {
			var X = tuple[0];
			var Y = tuple[1];
			var Z = tuple[2];
			var divider = X + 15 * Y + 3 * Z;
			var varU = 4 * X;
			var varV = 9 * Y;
			if (divider != 0) {
				varU /= divider;
				varV /= divider;
			} else {
				varU = NaN;
				varV = NaN;
			}
			var L = Hsluv.yToL(Y);
			if (L == 0) return [0, 0, 0];
			var U = 13 * L * (varU - Hsluv.refU);
			var V = 13 * L * (varV - Hsluv.refV);
			return [L, U, V];
		},

		luvToXyz: function (tuple) {
			var L = tuple[0];
			var U = tuple[1];
			var V = tuple[2];
			if (L == 0) return [0, 0, 0];
			var varU = U / (13 * L) + Hsluv.refU;
			var varV = V / (13 * L) + Hsluv.refV;
			var Y = Hsluv.lToY(L);
			var X = 0 - 9 * Y * varU / ((varU - 4) * varV - varU * varV);
			var Z = (9 * Y - 15 * varV * Y - varV * X) / (3 * varV);
			return [X, Y, Z];
		},

		luvToLch: function (tuple) {
			var L = tuple[0];
			var U = tuple[1];
			var V = tuple[2];
			var C = Math.sqrt(U * U + V * V);
			var H;
			if (C < 0.00000001) H = 0; else {
				var Hrad = Math.atan2(V, U);
				H = Hrad * 180.0 / 3.1415926535897932;
				if (H < 0) H = 360 + H;
			}
			return [L, C, H];
		},

		lchToLuv: function (tuple) {
			var L = tuple[0];
			var C = tuple[1];
			var H = tuple[2];
			var Hrad = H / 360.0 * 2 * Math.PI;
			var U = Math.cos(Hrad) * C;
			var V = Math.sin(Hrad) * C;
			return [L, U, V];
		},

		hsluvToLch: function (tuple) {
			var H = tuple[0];
			var S = tuple[1];
			var L = tuple[2];
			if (L > 99.9999999) return [100, 0, H];
			if (L < 0.00000001) return [0, 0, H];
			var max = Hsluv.maxChromaForLH(L, H);
			var C = max / 100 * S;
			return [L, C, H];
		},

		lchToHsluv: function (tuple) {
			var L = tuple[0];
			var C = tuple[1];
			var H = tuple[2];
			if (L > 99.9999999) return [H, 0, 100];
			if (L < 0.00000001) return [H, 0, 0];
			var max = Hsluv.maxChromaForLH(L, H);
			var S = C / max * 100;
			return [H, S, L];
		},

		hpluvToLch: function (tuple) {
			var H = tuple[0];
			var S = tuple[1];
			var L = tuple[2];
			if (L > 99.9999999) return [100, 0, H];
			if (L < 0.00000001) return [0, 0, H];
			var max = Hsluv.maxSafeChromaForL(L);
			var C = max / 100 * S;
			return [L, C, H];
		},

		lchToHpluv: function (tuple) {
			var L = tuple[0];
			var C = tuple[1];
			var H = tuple[2];
			if (L > 99.9999999) return [H, 0, 100];
			if (L < 0.00000001) return [H, 0, 0];
			var max = Hsluv.maxSafeChromaForL(L);
			var S = C / max * 100;
			return [H, S, L];
		},

		rgbToHex: function (tuple) {
			var h = "#";
			var _g1 = 0;
			var _g = tuple.length;
			while (_g1 < _g) {
				var i = _g1++;
				var chan = tuple[i];
				h += StringHex(Math.round(chan * 255), 2).toLowerCase();
			}
			return h;
		},

		hexToRgb: function (hex) {
			hex = hex.toUpperCase();
			return [parseInt("0x" + HxSubStr(hex, 1, 2)) / 255.0, parseInt("0x" + HxSubStr(hex, 3, 2)) / 255.0, parseInt("0x" + HxSubStr(hex, 5, 2)) / 255.0];
		},

		lchToRgb: function (tuple) {
			return Hsluv.xyzToRgb(Hsluv.luvToXyz(Hsluv.lchToLuv(tuple)));
		},

		rgbToLch: function (tuple) {
			return Hsluv.luvToLch(Hsluv.xyzToLuv(Hsluv.rgbToXyz(tuple)));
		},

		hsluvToRgb: function (tuple) {
			return Hsluv.lchToRgb(Hsluv.hsluvToLch(tuple));
		},

		rgbToHsluv: function (tuple) {
			return Hsluv.lchToHsluv(Hsluv.rgbToLch(tuple));
		},

		hpluvToRgb: function (tuple) {
			return Hsluv.lchToRgb(Hsluv.hpluvToLch(tuple));
		},

		rgbToHpluv: function (tuple) {
			return Hsluv.lchToHpluv(Hsluv.rgbToLch(tuple));
		},

		hsluvToHex: function (tuple) {
			return Hsluv.rgbToHex(Hsluv.hsluvToRgb(tuple));
		},

		hpluvToHex: function (tuple) {
			return Hsluv.rgbToHex(Hsluv.hpluvToRgb(tuple));
		},

		hexToHsluv: function (s) {
			return Hsluv.rgbToHsluv(Hsluv.hexToRgb(s));
		},

		hexToHpluv: function (s) {
			return Hsluv.rgbToHpluv(Hsluv.hexToRgb(s));
		},

		m: [[3.240969941904521, -1.537383177570093, -0.498610760293], [-0.96924363628087, 1.87596750150772, 0.041555057407175], [0.055630079696993, -0.20397695888897, 1.056971514242878]],

		minv: [[0.41239079926595, 0.35758433938387, 0.18048078840183], [0.21263900587151, 0.71516867876775, 0.072192315360733], [0.019330818715591, 0.11919477979462, 0.95053215224966]],

		refY: 1.0,

		refU: 0.19783000664283,

		refV: 0.46831999493879,

		kappa: 903.2962962,

		epsilon: 0.0088564516,

	};

	/**
	 *
	 * Generate colors. Eyeballing, distinct ranges seem to be:
	 * 0-360  hues
	 * 40-100 saturation
	 * 20-80  lightness (below and above this, colors fade too much)
	 *
	 * However, eyeballing with my own protanomally as a
	 * "worst case baseline", assuming other color-vision
	 * deficiencies are similar except with other hues,
	 * the following rought rules of thumb seem to hold:
	 *
	 *   - changing hues with steps of 60 is a minimum
	 *     at all ranges, so six variations for hue
	 *
	 *   - lightness brings relatively high shifts in
	 *     perception for all saturations and hues,
	 *     steps of 6 seem ok, giving eleven options in total:
	 *     [20, 26, 32, 38, 44, 50, 56, 62, 68, 74, 80]
	 *
	 *   - color distinction peaks with 62 and 68 lightness,
	 *     so these have a preference when dealing with
	 *     fewer categories
	 *
	 *   - for 44, 50 and 56 and 62 lightness, saturation steps
	 *     40, 70 and 100 are required for distinction (3 options)
	 *
	 *   - 32 and 74 lightness, saturation should be either
	 *     60 or 100 to avoid overly similar colors (2 options)

	 *   - 20, 26, and 80 lightness should only have 100 saturation
	 *
	 *   - improving distinction between close colors is most
	 *     important. Since hues are cyclic, and we only care
	 *     about the distance between hues for a given L and S,
	 *     we can make a "checkerboard pattern" when varying
	 *     these values by alternating a +30 hue offset.
	 *     That means that identical hues are separated
	 *     by at least two steps of L+S, which hopefully
	 *     improves things a bit.
	 *
	 * That leads to a max total of 114 HSLuv colors, many of which
	 * will be very close of course. The trick is to pick less colors
	 * by default
	*/

	// include actual HSLuv values, and "normalised" values
	// used later for picking colors that are the furthest apart
	let hslTuples = [], normTuples = [];

	// later on we need a starting color, a "seed"
	// Let's take [60, 100, 80], it's bright yellow-ish.
	// The furthest away color is [240, 100, 20],
	// dark-blue/purplish
	let seedColorIndex = 0,
		seedH = 60,
		seedS = 100,
		seedL = 80;

	// Eleven lightness steps
	for (let i = 0; i < 11; i++) {
		// in order of "most perceptually important"
		let l, s, h;

		// perceived color distance is not weighed equally,
		// even with the step sizes defined above.
		// For all-else-equal values to me the following holds
		// when I compare single-step changes:
		//
		// - lightness has a bigger impact on perceptual changes,
		//   and does not depend on colorvision
		// - saturation second
		// - hue the least (especially at lower saturation/lightness)
		//
		// So we weigh l > s > h
		// (I'll probably tweak this to see the effect on the output)
		const lWeight = 2, sWeight = 1.25, hWeight = 1;

		l = 80 - i * 6;

		switch (l) {
			case 44:
			case 50:
			case 56:
			case 62:
				// three saturation steps
				for (let j = 0; j < 3; j++) {
					// go from most to least saturated
					s = 100 - j * 30;

					// six hue steps
					for (let k = 0; k < 6; k++) {
						// bit-twiddling hack to make
						// hue-offset "checkerboard"
						// dependent on lightness and
						// saturation steps
						h = k * 60 + ((i + j) & 1) * 30;

						// save tuple
						hslTuples.push([h, s, l]);

						// save position in "normalised" color space
						normTuples.push([i * lWeight, j * sWeight, k * hWeight]);

					}
				}
				break;
			case 32:
			case 74:
				// two saturation steps
				for (let j = 0; j < 2; j++) {
					// go from most to least saturated
					s = 100 - j * 40;

					// six hue steps
					for (let k = 0; k < 6; k++) {
						// bit-twiddling hack to make
						// hue-offset "checkerboard"
						// dependent on lightness and
						// saturation steps
						h = k * 60 + ((i + j) & 1) * 30;

						// save tuple
						hslTuples.push([h, s, l]);

						// save position in "normalised" color space
						normTuples.push([i * lWeight, j * sWeight, k * hWeight]);

					}
				}
				break;
			case 20:
			case 26:
			case 80:
				// no saturation loop, always 100
				s = 100;
				for (let k = 0; k < 6; k++) {
					// bit-twiddling hack to make
					// hue-offset "checkerboard"
					// dependent on lightness and
					// saturation steps
					h = k * 60 + (i & 1) * 30;

					// If we match the chosen seed-color,
					// store its the hslTuples index
					if (h == seedH && s == seedS && l == seedL) {
						seedColorIndex = hslTuples.length;
					}

					// save tuple
					hslTuples.push([h, s, l]);

					// save position in "normalised" color space
					normTuples.push([i * lWeight, 0, k * hWeight]);
				}

		}
	}

	// convert to RGB
	const hexColors = hslTuples.map(hsl => {
		return Hsluv.hsluvToHex(hsl);
	});
	console.log(hexColors);

	// Make "distance" matrix according to my own human metric (again)
	const totalColors = normTuples.length;
	let distance = new Float64Array(totalColors * totalColors);

	const { sqrt } = Math;
	for (let i = 0; i < totalColors; i++) {
		for (let j = i + 1; j < totalColors; j++) {
			const [hi, si, li] = normTuples[i];
			const [hj, sj, lj] = normTuples[j];
			// hue is cyclical
			const dh = Math.min(Math.abs(hi - hj), 6 - Math.abs(hi - hj));
			// will  be squared, no need to take absolute value
			const ds = si - sj;
			const dl = li - lj;

			distance[i + j * totalColors] =
				distance[j + i * totalColors] =
				sqrt(dh * dh + ds * ds + dl * dl);
		}
	}

	// actual color selection.
	// Ideally: maximise the minimal color distance
	// Practically: unfeasible. Unless I can somehow implement
	// TSP in JavaScript.

	// First approach tried is Tim Holy's approximation:
	// start with one color, find most distant color from it,
	// find next most distant color from those two, etc.
	// Color distance will quickly decrease as pallette does
	// let there be N total colors (in our case, N = 198)

	let,
		colors = [hexColors[seedColorIndex]],
		colorIndices = [seedColorIndex],
		nColors = [colors],
		allIndices = hexColors.map((_, index) => {
		// skip the seedColorIndex, note that
		// we'll end up with an extra index
		return index + (index < seedColorIndex ? 0 : 1);
	});
	// remove extra index
	allIndices.pop();


	for (let i = 1; i < totalColors; i++) {
		colors = colors.slice(0);
		colorIndices = colorIndices.slice(0);
		let maxMinDist = 0;
		let maxMindDistIdx = -1;
		// for each new color, go through the distance
		// matrix and find the color with the biggest
		// minimum distance from the current color
		for (let j = 0; j < i; j++) {
			for (let k = 0; k < allIndices.length; k++) {
				const kIdx = allIndices[k]*totalColors;

				// start with distance from first color
				let dist = distance[colorIndices[0] + kIdx];

				for (let l = 1; l < colorIndices.length; l++) {
					let nDist = distance[colorIndices[l] + kIdx];
					// keep track of smallest distance
					// compared to previous colors
					dist = nDist < dist ? nDist : dist;
				}

			}
		}
	}
	// add white in front of each color
	colors = colors.map(arr => arr.unshift('#FFFFFF'));
})();