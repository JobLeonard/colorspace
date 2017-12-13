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

},

var Hsluv = {
	getBounds: function (L) {
		var result = [];
		var sub1 = Math.pow(L + 16, 3) / 1560896;
		var sub2;
		if (sub1 > Hsluv.epsilon) sub2 = sub1; else sub2 = L / Hsluv.kappa;
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
			var length = distanceLineFromOrigin(bounds[i]);
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
			var length = lengthOfRayUntilIntersect(hrad, bound);
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
 * 20-80  lightness
 *
 * However, eyeballing with my own protanomally as a
 * "worst case baseline", assuming other color-vision
 * deficiencies are similar but with other hues,
 * the following rought rules of thumb seem to hold:
 *
 *   - changing hues with steps of 60 is a minimum
 *     at all ranges, so six variations for hue
 *
 *   - lightness brings relatively high shifts in
 *     perception for all saturations and hues,
 *     steps of 10 seem ok, giving seven options in total:
 *     [20, 30, 40, 50, 60, 70, 80]
 *
 *   - Color distinction peaks at 60 and 70 lightness,
 *     so these have a preference when dealing with
 *     fewer categories
 *
 *   - for 40, 50 and 60 lightness, saturation steps
 *     40, 70 and 100 are required for distinction (3 options)
 *
 *   - 30 and 70 lightness, saturation should be either
 *     60 or 100 to avoid overly similar colors (2 options)

 *   - 20 and 80 lightness should only have 100 saturation
 *
 *
 * That leads to the 90 following HSLuv colors, written in [L, S, H]
 * because I'm an idiot and can't be bothered to type this out again,
 * so I'll just use `.map(lsh => [lsh[2], lsh[1], lsh[0]])` at the end:


	// 60 lightness, 100 saturation
	[60, 100, 0], [60, 100, 60], [60, 100, 120],
	[60, 100, 180], [60, 100, 240], [60, 100, 300],
	// 60 lightness, 70 saturation
	[60, 70, 0], [60, 70, 60], [60, 70, 120],
	[60, 70, 180], [60, 70, 240], [60, 70, 300],
	// 60 lightness, 40 saturation
	[60, 40, 0], [60, 40, 60], [60, 40, 120],
	[60, 40, 180], [60, 40, 240], [60, 40, 300],

	// 40 lightness, 100 saturation
	[40, 100, 0], [40, 100, 60], [40, 100, 120],
	[40, 100, 180], [40, 100, 240], [40, 100, 300],
	// 40 lightness, 70 saturation
	[40, 70, 0], [40, 70, 60], [40, 70, 120],
	[40, 70, 180], [40, 70, 240], [40, 70, 300],
	// 40 lightness, 40 saturation
	[40, 40, 0], [40, 40, 60], [40, 40, 120],
	[40, 40, 180], [40, 40, 240], [40, 40, 300],

	// 50 lightness, 100 saturation
	[50, 100, 0], [50, 100, 60], [50, 100, 120],
	[50, 100, 180], [50, 100, 240], [50, 100, 300],
	// 50 lightness, 70 saturation
	[50, 70, 0], [50, 70, 60], [50, 70, 120],
	[50, 70, 180], [50, 70, 240], [50, 70, 300],
	// 50 lightness, 40 saturation
	[50, 40, 0], [50, 40, 60], [50, 40, 120],
	[50, 40, 180], [50, 40, 240], [50, 40, 300],

	// 70 lightness, 100 saturation
	[70, 100, 0], [70, 100, 60], [70, 100, 120],
	[70, 100, 180], [70, 100, 240], [70, 100, 300],
	// 70 lightness, 60 saturation
	[70, 60, 0], [70, 60, 60], [70, 60, 120],
	[70, 60, 180], [70, 60, 240], [70, 60, 300],

	// 30 lightness, 100 saturation
	[30, 100, 0], [30, 100, 60], [30, 100, 120],
	[30, 100, 180], [30, 100, 240], [30, 100, 300],
	// 30 lightness, 60 saturation
	[30, 60, 0], [30, 60, 60], [30, 60, 120],
	[30, 60, 180], [30, 60, 240], [30, 60, 300],

	// 80 lightness, 100 saturation
	[80, 100, 0], [80, 100, 60], [80, 100, 120],
	[80, 100, 180], [80, 100, 240], [80, 100, 300],

	// 20 lightness, 100 saturation
	[20, 100, 0], [20, 100, 60], [20, 100, 120],
	[20, 100, 180], [20, 100, 240], [20, 100, 300],

 */

let tuples = [
	// 60 lightness, 100 saturation
	[60, 100, 0], [60, 100, 60], [60, 100, 120],
	[60, 100, 180], [60, 100, 240], [60, 100, 300],
	// 60 lightness, 70 saturation
	[60, 70, 0], [60, 70, 60], [60, 70, 120],
	[60, 70, 180], [60, 70, 240], [60, 70, 300],
	// 60 lightness, 40 saturation
	[60, 40, 0], [60, 40, 60], [60, 40, 120],
	[60, 40, 180], [60, 40, 240], [60, 40, 300],

	// 40 lightness, 100 saturation
	[40, 100, 0], [40, 100, 60], [40, 100, 120],
	[40, 100, 180], [40, 100, 240], [40, 100, 300],
	// 40 lightness, 70 saturation
	[40, 70, 0], [40, 70, 60], [40, 70, 120],
	[40, 70, 180], [40, 70, 240], [40, 70, 300],
	// 40 lightness, 40 saturation
	[40, 40, 0], [40, 40, 60], [40, 40, 120],
	[40, 40, 180], [40, 40, 240], [40, 40, 300],

	// 50 lightness, 100 saturation
	[50, 100, 0], [50, 100, 60], [50, 100, 120],
	[50, 100, 180], [50, 100, 240], [50, 100, 300],
	// 50 lightness, 70 saturation
	[50, 70, 0], [50, 70, 60], [50, 70, 120],
	[50, 70, 180], [50, 70, 240], [50, 70, 300],
	// 50 lightness, 40 saturation
	[50, 40, 0], [50, 40, 60], [50, 40, 120],
	[50, 40, 180], [50, 40, 240], [50, 40, 300],

	// 70 lightness, 100 saturation
	[70, 100, 0], [70, 100, 60], [70, 100, 120],
	[70, 100, 180], [70, 100, 240], [70, 100, 300],
	// 70 lightness, 60 saturation
	[70, 60, 0], [70, 60, 60], [70, 60, 120],
	[70, 60, 180], [70, 60, 240], [70, 60, 300],

	// 30 lightness, 100 saturation
	[30, 100, 0], [30, 100, 60], [30, 100, 120],
	[30, 100, 180], [30, 100, 240], [30, 100, 300],
	// 30 lightness, 60 saturation
	[30, 60, 0], [30, 60, 60], [30, 60, 120],
	[30, 60, 180], [30, 60, 240], [30, 60, 300],

	// 80 lightness, 100 saturation
	[80, 100, 0], [80, 100, 60], [80, 100, 120],
	[80, 100, 180], [80, 100, 240], [80, 100, 300],

	// 20 lightness, 100 saturation
	[20, 100, 0], [20, 100, 60], [20, 100, 120],
	[20, 100, 180], [20, 100, 240], [20, 100, 300],
].map(lsh => [lsh[2], lsh[1], lsh[0]]);

// needed for distance matrix
let luvTuples = tuples.map(hsl => {
	Hsluv.lchToLuv(Hsluv.hsluvToLch(hsl));
});

const { sqrt, abs } = Math;

let distance = new Array(luvTuples.length);
for (let i = 0; i < distance.length; i++){
	distance[i] = new Float64Array(distance.length);
}
for (let i = 0; i < luvTuples.length; i++) {
	const [l1, u1, v1] = luvTuples[i];
	for (let j = 0; j < luvTuples.length; j++) {
		const [l2, u2, v1] = luvTuples[j];
		distance[i][j] = sqrt(
			abs(l1 * l1 - l2 * l2) +
			abs(u1 * u1 - u2 * u2) +
			abs(v1 * v1 - v2 * v2)
		);
	}
}