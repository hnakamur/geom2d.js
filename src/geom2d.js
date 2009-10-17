var geom2d = function() {

function Bezier(points) {
  this.points = points;
}
Bezier.epsilon = 1e-6;
Bezier.intersections = function(bezier1, bezier2) {
  var c1 = bezier1.coefficients(), c2 = bezier2.coefficients(),
      n = Math.max(bezier1.degree(), bezier2.degree()), cx = [], cy = [],
      i, matrix, dcmp;
  for (i = 0; i <= n; ++i) {
    cx[i] = (c1[i].x() || 0) - (c2[i].x() || 0);
console.log('i=' + i + ', c1x=' + c1[i].x() + ', c2x=' + c2[i].x() + ', cx=' + cx[i]);
    cy[i] = (c1[i].y() || 0) - (c2[i].y() || 0);
console.log('i=' + i + ', c1y=' + c1[i].y() + ', c2y=' + c2[i].y() + ', cy=' + cy[i]);
  }
  matrix = new numeric.Matrix([
    [cx[3], cx[2], cx[1], cx[0], 0, 0],
    [0, cx[3], cx[2], cx[1], cx[0], 0],
    [0, 0, cx[3], cx[2], cx[1], cx[0]],
    [cy[3], cy[2], cy[1], cy[0], 0, 0],
    [0, cy[3], cy[2], cy[1], cy[0], 0],
    [0, 0, cy[3], cy[2], cy[1], cy[0]]
  ]);
  console.log('matrix=' + matrix);
  dcmp = new numeric.LUDecomposition(matrix);
  var matrixLU = dcmp.matrixL().multiply(dcmp.matrixU());
  var matrixPA = dcmp.matrixP().multiply(matrix);
  console.log('PA=' + matrixPA);
  console.log('LU=' + matrixLU);
  console.log('det=' + dcmp.determinant());
  console.log('L=' + dcmp.matrixL());
  console.log('U=' + dcmp.matrixU());
  console.log('P=\n' + dcmp.matrixP());
  var vectorB = numeric.Vector.zero(matrix.rowSize);
  var vectorBp = dcmp.solve(vectorB);
  console.log('bp=' + vectorBp);
  var ge = new numeric.GaussElimination(matrix, vectorB);
  console.log('singular=' + ge.isSingular());
  console.log('solution=' + ge.solution());
};
var proto = Bezier.prototype;
proto.degree = function() {
  return this.points.length - 1;
};
proto.pointAtT = function(t) {
  return numeric.Vector.polynomial(t, this.coefficients());
};
proto.pointAtT2 = function(t) {
  var n = this.degree(), t1 = 1 - t, p = this.points, q = [], i, j;
  for (i = n; i >= 1; --i) {
    for (j = 0; j < i; ++j)
      q[j] = p[j].scalarMult(t1).add(p[j + 1].scalarMult(t));
    p = q;
    q = [];
  }
  return p[0];
};
proto.coefficients = function() {
  if (!this._coefficients) {
    var p = this.points;
    switch (this.degree()) {
    case 1:
      var p0 = p[0], p1 = p[1];
      this._coefficients = [p0, p1.subtract(p0)];
      break;
    case 2:
      var p0 = p[0], p1 = p[1], p2 = p[2],
          p10 = p1.subtract(p0), p21 = p2.subtract(p1);
      this._coefficients = [p0, p10.scalarMult(2), p21.subtract(p10)];
      break;
    case 3:
      var p0 = p[0], p1 = p[1], p2 = p[2], p3 = p[3],
          p10 = p1.subtract(p0), p21 = p2.subtract(p1), p32 = p3.subtract(p2),
          p21_10 = p21.subtract(p10), p32_21 = p32.subtract(p21);
      this._coefficients =
        [p0, p10.scalarMult(3), p21_10.scalarMult(3), p32_21.subtract(p21_10)];
      break;
    }
  }
  return this._coefficients;
};
proto.derivativeAtT = function(t) {
  return numeric.Vector.polynomial(t, this.derivativeCoefficients());
};
proto.derivativeCoefficients = function() {
  if (!this._derivativeCoefficients) {
    var c = this.coefficients();
    switch (this.degree()) {
    case 1:
      this._derivativeCoefficients = [c[1]];
      break;
    case 2:
      this._derivativeCoefficients = [c[1], c[2].scalarMult(2)];
      break;
    case 3:
      this._derivativeCoefficients = [c[1], c[2].scalarMult(2), c[3].scalarMult(3)];
      break;
    }
  }
  return this._derivativeCoefficients;
};
proto.secondDerivativeAtT = function(t) {
  return numeric.Vector.polynomial(t, this.secondDerivativeCoefficients());
};
proto.secondDerivativeCoefficients = function() {
  if (!this._secondDerivativeCoefficients) {
    var c = this.coefficients();
    switch (this.degree()) {
    case 1:
      this._secondDerivativeCoefficients = [0];
      break;
    case 2:
      this._secondDerivativeCoefficients = [c[2].scalarMult(2)];
      break;
    case 3:
      this._secondDerivativeCoefficients = [c[2].scalarMult(2), c[3].scalarMult(6)];
      break;
    }
  }
  return this._secondDerivativeCoefficients;
};
proto.inflectionPointTs = function(includeEnds) {
  // http://www.caffeineowl.com/graphics/2d/vectorial/cubic-inflexion.html
  var coef = this.coefficients();
  switch (this.degree()) {
  case 1:
  case 2:
    return [];
  case 3:
    var oneThird = 1 / 3,
      a = coef[1].scalarMult(oneThird), b = coef[2].scalarMult(oneThird), c = coef[3];
      roots = new numeric.QuadraticEquation(
        b.crossProduct(c).z(), a.crossProduct(c).z(), a.crossProduct(b).z()).
        realRoots();
    var predicate = includeEnds ?
      function(t) { return 0 <= t && t <= 1; } :
      function(t) { return Bezier.epsilon <= t && t <= 1 - Bezier.epsilon; };
    return filterValues(roots, predicate);
  }
};
proto.subdivideAtInfectionPoints = function() {
  var ts = this.inflectionPointTs(), n = ts.length;
  switch (n) {
  case 0:
    return [this];
  case 1:
    return this.subdivideAtT(ts[0]);
  case 2:
    var t1 = ts[1], segments = this.subdivideAtT(t1);
    return segments[0].subdivideAtT(ts[0] / t1).concat(segments[1]);
  }
};
proto.intermediatePointsAtT = function(t) {
  var newPoints = [], points = this.points, n = points.length;
  for (var i = 0; i < n - 1; ++i)
    newPoints.push(new Bezier([points[i], points[i + 1]]).pointAtT(t));
  return newPoints;
};
proto.subdivideAtT = function(t) {
  var p = this.points;
  switch (this.degree()) {
  case 1:
    var q = this.pointAtT(t);
    return [new Bezier([p[0], q]), new Bezier(q, p[1])];
  case 2:
    var q = this.intermediatePointsAtT(t),
        r = new Bezier(q).pointAtT(t);
    return [new Bezier([p[0], pp[0], r]), new Bezier([r, pp[1], p[2]])];
  case 3:
    var q = this.intermediatePointsAtT(t),
        r = new Bezier(q).intermediatePointsAtT(t),
        s = new Bezier(r).pointAtT(t);
    return [new Bezier([p[0], q[0], r[0], s]),
      new Bezier([s, r[1], q[2], p[3]])];
  }
};
proto.calcCurveLength = function() {
//    if (this.inflectionPointTs().length > 0) {
//      var segments = this.subdivideAtInfectionPoints();
//      var length = 0;
//      for (var i = 0, n = segments.length; i < n; ++i)
//        length += segments[i].calcCurveLength();
//      return length;
//    }

  var segments = this.segments();
  return segments[segments.length - 1].accLen;
};
proto.calcCurveLength2 = function() {
  var self = this;
  function f(t) {
    return self.derivativeAtT(t).length();
  }
  return new numeric.integral(f, 0, 1, Bezier.epsilon);
};
proto.getTAtLength = function(length) {
  var segments = this.segments();
  var index = binarySearch(segments, 'accLen', length);
  if (index == -1)
    return undefined;
  var segment = segments[index];
  var accLen = segment.accLen;
  if (length == accLen || index == segments.length - 1)
    return segment.t;
  var nextSegment = segments[index + 1];
  return calcLinearInterpolation(length, accLen, nextSegment.accLen,
      segment.t, nextSegment.t);
};
proto.segments = function() {
  if (!this._segments) {
    var segments = [];
    var self = this;
    var length = calcIntegrationBySimpson(function(t) {
      var deriv = self.derivativeAtT(t);
      var derivLen = deriv.length();
      segments.push({t: t, deriv: deriv, derivLen: derivLen});
      return derivLen;
    }, 0, 1, Bezier.epsilon);

    segments.sort(function(a, b) { return a.t - b.t; });

    var n = segments.length, h = 1 / (n - 1);
    var derivLen0 = segments[0].derivLen;
    var sumOdd = 0, sumEven = 0;
    var i = 0;
    segments[i].len = 0;
    while (++i < n) {
      var derivLen = segments[i].derivLen;
      segments[i].accLen =
        h / 3 * (derivLen0 + 4 * sumOdd + 2 * sumEven + derivLen);
      if (i % 2)
        sumOdd += derivLen;
      else
        sumEven += derivLen;
    }
    this._segments = segments;
  }
  return this._segments;
};

function filterValues(values, predicate) {
  if (values) {
    var ret = [], value;
    for (var i = 0, n = values.length; i < n; ++i) {
      value = values[i];
      if (predicate(value))
        ret.push(value);
    }
    return ret;
  }
  else
    return values;
}

function binarySearch(mappings, searchIndex, searchValue) {
  // http://en.wikipedia.org/wiki/Binary_search_algorithm
  var iMin = 0, iMax = mappings.length - 1, iMid, midValue;
  if (searchValue < mappings[iMin][searchIndex] ||
      searchValue > mappings[iMax][searchIndex])
    return -1;
  do {
    iMid = Math.floor((iMin + iMax) / 2);
    midValue = mappings[iMid][searchIndex];
    if (searchValue > midValue)
      iMin = iMid + 1;
    else
      iMax = iMid - 1;
  } while (midValue !== searchValue && iMin <= iMax);
  return iMid;
}

function calcLinearInterpolation(x, x0, x1, y0, y1) {
  // http://en.wikipedia.org/wiki/Linear_interpolation
  return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

return {
  Bezier: Bezier
};
  
}();
