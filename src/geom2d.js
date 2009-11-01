var geom2d = function() {

var sum = numeric.sum, numberEquals = numeric.numberEquals;

// immutable two dimensional vector.
function Vector2d(x, y) {
  this.x = x;
  this.y = y;
}
mix(Vector2d, {
  vectorsFromNestedXYArray: function vectorsFromNestedXYArray(xyArray) {
    var vectors = [], n = xyArray.length, i, xy;
    for (i = 0; i < n; ++i) {
      xy = xyArray[i];
      vectors.push(new Vector2d(xy[0], xy[1]));
    }
    return vectors;
  },
  sum: function sum(vectors) {
    var xs = [], ys = [], n = vectors.length, i, v;
    for (i = 0; i < n; ++i) {
      v = vectors[i];
      xs.push(v.x);
      ys.push(v.y);
    }
    return new Vector2d(sum(xs), sum(ys));
  },
  polynomial: function polynomial(t, coefficients) {
    var xs = [], ys = [], n = coefficients.length, ti = 1, i, v;
    for (i = 0; i < n; ++i, ti *= t) {
      v = coefficients[i];
      xs.push(v.x * ti);
      ys.push(v.y * ti);
    }
    return new Vector2d(sum(xs), sum(ys));
  },
  ZERO: new Vector2d(0, 0)
});
mix(Vector2d.prototype, {
  angleTo: function angleTo(vector) {
    // p. 1120, Numerical Recipe 3rd ed.
    return Math.atan2(this.zComponentOfCrossProduct(vector),
        this.dotProduct(vector));
  },
  toString: function toString() {
    return '[' + this.x + ' ' + this.y + ']';
  },
  unit: function unit() {
    return this.scalarDiv(this.length());
  },
  length: function length() {
    if (isUndefined(this._length))
      this._length = Math.sqrt(this.dotProduct(this));
    return this._length;
  },
  dotProduct: function dotProduct(vector) {
    return sum([this.x * vector.x, this.y * vector.y]);
  },
  zComponentOfCrossProduct: function zComponentOfCrossProduct(vector) {
    return this.x * vector.y - this.y * vector.x;
  },
  equals: function equals(vector, epsilon) {
    return numberEquals(this.x, vector.x, epsilon) &&
      numberEquals(this.y, vector.y, epsilon);
  },
  equalsToOneOf: function equalsToOneOf(vectors, epsilon) {
    for (var i = 0, n = vectors.length; i < n; ++i) {
      if (this.equals(vectors[i], epsilon))
        return true;
    }
    return false;
  },
  add: function add(vector) {
    return new Vector2d(this.x + vector.x, this.y + vector.y);
  },
  subtract: function subtract(vector) {
    return new Vector2d(this.x - vector.x, this.y - vector.y);
  },
  negate: function negate() {
    return new Vector2d(-this.x, -this.y);
  },
  scalarMult: function scalarMult(factor) {
    return new Vector2d(this.x * factor, this.y * factor);
  },
  scalarDiv: function scalarDiv(factor) {
    return new Vector2d(this.x / factor, this.y / factor);
  }
});

function Line(p0, p1) {
  this.p0 = p0;
  this.p1 = p1;
}
mix(Line.prototype, {
  bbox: function bbox() {
    if (!this._bbox)
      this._bbox = Rectangle.fromTwoCornerPoints(this.p0, this.p1);
    return this._bbox;
  },
  coefficients: function coefficients() {
    if (!this._coefficients)
      this._coefficients = [this.p0, this.dv()];
    return this._coefficients;
  },
  distance: function distance(point) {
    return Math.abs(
      point.subtract(this.p0).zComponentOfCrossProduct(this.dv().unit()));
  },
  dv: function dv() {
    if (!this._dv)
      this._dv = this.p1.subtract(this.p0);
    return this._dv;
  },
  intersection: function intersection(line) {
    var ts = this.tsAtIntersection(line);
    return ts ? this.pointAtT(ts[0]) : ts;
  },
  length: function length() {
    return this.dv().length();
  },
  tAtLength: function tAtLength(length) {
    return length / this.length();
  },
  tAtX: function tAtX(x) {
    return (x - this.p0.x) / this.dv().x;
  },
  tAtY: function tAtY(y) {
    return (y - this.p0.y) / this.dv().y;
  },
  tsAtIntersection: function tsAtIntersection(line) {
    var dp0 = line.p0.subtract(this.p0), dv1 = this.dv(), dv2 = line.dv(),
        a = dv1.x, b = -dv2.x, c = dv1.y, d = -dv2.y, e = dp0.x, f = dp0.y,
        det = a * d - b * c, t, u;
    if (det !== 0) {
      t = (d * e - b * f) / det, u = (a * f - c * e) / det;
      return [t, u];
    }
    else
      return undefined;
  },
  pointAtT: function pointAtT(t) {
    return Vector2d.polynomial(t, this.coefficients());
  }
});

function Rectangle(x, y, width, height) {
  this.x = x;
  this.y = y;
  this.width = width;
  this.height = height;
}
mix(Rectangle, {
  fromTwoCornerPoints: function fromTwoCornerPoints(p0, p1) {
    var x0 = p0.x, y0 = p0.y, x1 = p1.x, y1 = p1.y;
    return new Rectangle(Math.min(x0, x1), Math.min(y0, y1),
      Math.abs(x1 - x0), Math.abs(y1 - y0));
  }
});
mix(Rectangle.prototype, {
  intersects: function intersects(rect) {
    return this.x <= rect.x + rect.width && this.x + this.width >= rect.x &&
      this.y <= rect.y + rect.height && this.y + this.height >= rect.y;
  }
});

function QuadraticBezier(p0, p1, p2) {
  this.p0 = p0;
  this.p1 = p1;
  this.p2 = p2;
}
mix(QuadraticBezier.prototype, {
  coefficients: function coefficients() {
    if (!this._coefficients) {
      var p0 = this.p0, p1 = this.p1, p2 = this.p2,
          p10 = p1.subtract(p0), p21 = p2.subtract(p1);
      this._coefficients = [p0, p10.scalarMult(2), p21.subtract(p10)];
    }
    return this._coefficients;
  },
  pointAtT: function pointAtT(t) {
    return Vector2d.polynomial(t, this.coefficients());
  },
  tAtXY: function tAtPoint(x, y) {
    var c = this.coefficients(), c0 = c[0], c1 = c[1], c2 = c[2];
    var a0 = c0.x, a1 = c1.x, a2 = c2.x;
    var b0 = c0.y, b1 = c1.y, b2 = c2.y;
    var d21 = a2 * b1 - a1 * b2;
    var t = -sum([b2 * x, -a2 * y, a2 * b0 - a0 * b2]) / d21;
    if (isNaN(t)) {
      t = -sum([b1 * x, -a1 * y, a1 * b0 - a0 * b1]) /
        sum([b2 * x, -a2 * y, a2 * b0 - a0 * b2]);
    }
    return t;
  },
  tPairsAtIntersections: function tPairsAtIntersections(bezier) {
    var tc = this.coefficients();
    var a0 = tc[0].x, a1 = tc[1].x, a2 = tc[2].x;
    var b0 = tc[0].y, b1 = tc[1].y, b2 = tc[2].y;
    var a2b1 = a2 * b1 - a1 * b2;
    var a2b0 = a2 * b0 - a0 * b2;
    var a1b0 = a1 * b0 - a0 * b1;
    // implicitized curve 1: a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
    var a = -b2 * b2;
    var b = 2 * a2 * b2;
    var c = -a2 * a2;
    var d = b1 * a2b1 - 2 * b2 * a2b0;
    var e = -a1 * a2b1 + 2 * a2 * a2b0;
    var f = a2b1 * a1b0 - a2b0 * a2b0;

    var bc = bezier.coefficients();
    var c0 = bc[0].x, c1 = bc[1].x, c2 = bc[2].x;
    var d0 = bc[0].y, d1 = bc[1].y, d2 = bc[2].y;

    var poly = new numeric.Polynomial([
      a*c0*c0+b*c0*d0+c*d0*d0+d*c0+e*d0+f,
      2*a*c0*c1+b*(c1*d0+c0*d1)+2*c*d0*d1+d*c1+e*d1,
      a*(2*c0*c2+c1*c1)+b*(c2*d0+c0*d2+c1*d1)+c*(2*d0*d2+d1*d1)+d*c2+e*d2,
      2*a*c1*c2+b*(c2*d1+c1*d2)+2*c*d1*d2,
      a*c2*c2+b*c2*d2+c*d2*d2
    ]);
    console.log('poly.coefficients=' + JSON.stringify(poly.coefficients));
    var us = poly.realRootsBetween(0, 1, 1e-3);
    console.log('us=' + JSON.stringify(us));

    var tPairs = [];
    for (var i = 0, n = us.length; i < n; ++i) {
      var u = us[i];
      var p = bezier.pointAtT(u);
      var t = this.tAtXY(p.x, p.y);
      tPairs.push([t, u]);
    }
    return tPairs;
  },
  tsAtX: function tsAtX(x) {
    var c = this.coefficients;
    return realRootsOfQuadraticEquation(c[2].x, c[1].x, c[0].x - x,
      function(t) { return 0 <= t && t <= 1; }
    );
  },
  tsAtY: function tsAtY(y) {
    var c = this.coefficients;
    return realRootsOfQuadraticEquation(c[2].y, c[1].y, c[0].y - y,
      function(t) { return 0 <= t && t <= 1; }
    );
  }
});

function Bezier() {
  if (arguments.length === 1)
    this.points = arguments[0];
  else
    this.points = arguments;
}
Bezier.epsilon = 1e-6;
Bezier.intersections2 = function(bezier1, bezier2) {
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
mix(Bezier.prototype, {
  degree: function degree() {
    return this.points.length - 1;
  },
  pointAtT: function pointAtT(t) {
    return Vector2d.polynomial(t, this.coefficients());
  },
  coefficients: function coefficients() {
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
          [p0, p10.scalarMult(3), p21_10.scalarMult(3),
          p32_21.subtract(p21_10)];
        break;
      }
    }
    return this._coefficients;
  },
  derivativeAtT: function derivativeAtT(t) {
    return Vector2d.polynomial(t, this.derivativeCoefficients());
  },
  derivativeCoefficients: function derivativeCoefficients() {
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
        this._derivativeCoefficients =
          [c[1], c[2].scalarMult(2), c[3].scalarMult(3)];
        break;
      }
    }
    return this._derivativeCoefficients;
  },
  inflectionPointTs: function inflectionPointTs() {
    // http://www.caffeineowl.com/graphics/2d/vectorial/cubic-inflexion.html
    var coef = this.coefficients();
    switch (this.degree()) {
    case 1:
    case 2:
      return [];
    case 3:
      var a = coef[1].scalarDiv(3), b = coef[2].scalarDiv(3), c = coef[3];
      return realRootsOfQuadraticEquation(
        b.zComponentOfCrossProduct(c),
        a.zComponentOfCrossProduct(c),
        a.zComponentOfCrossProduct(b),
        function(t) { return 0 <= t && t <= 1; }
      );
    }
  },
  secondDerivativeAtT: function secondDerivativeAtT(t) {
    return Vector2d.polynomial(t, this.secondDerivativeCoefficients());
  },
  secondDerivativeCoefficients: function secondDerivativeCoefficients() {
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
        this._secondDerivativeCoefficients =
          [c[2].scalarMult(2), c[3].scalarMult(6)];
        break;
      }
    }
    return this._secondDerivativeCoefficients;
  },
  intermediatePointsAtT: function intermediatePointsAtT(t) {
    var newPoints = [], points = this.points, n = points.length;
    for (var i = 0; i < n - 1; ++i)
      newPoints.push(new Line(points[i], points[i + 1]).pointAtT(t));
    return newPoints;
  },
  tsOfZeroXYDerivPoints: function tsOfZeroXYDerivPoints() {
    var cv = this.derivativeCoefficients();
    var c0 = cv[0], c1 = cv[1], c2 = cv[2] || Vector2d.ZERO;
    var filterFunc = function (t) { return 0 < t && t < 1; };
    var ts = numeric.uniqAndSort(
      [0, 1].concat(
        realRootsOfQuadraticEquation(c2.x, c1.x, c0.x, filterFunc),
        realRootsOfQuadraticEquation(c2.y, c1.y, c0.y, filterFunc)
      )
    );
    return ts;
  },
  segmentsDividedAtTs: function segmentsDividedAtTs(ts) {
    var segments = [];
    for (var i = 0, n = ts.length - 1; i < n; ++i) {
      var t1 = ts[i], t2 = ts[i + 1];
      var p1 = i == 0 ? this.pointAtT(t1) : segments[i - 1].p2;
      var p2 = this.pointAtT(t2);
      segments.push(new BezierSegment(this, t1, t2, p1, p2));
    }
    return segments;
  },
  beziersDividedAtT: function beziersDividedAtT(t) {
    var p = this.points;
    switch (this.degree()) {
    case 1:
      var q = this.pointAtT(t);
      return [new Bezier(p[0], q), new Bezier(q, p[1])];
    case 2:
      var q = this.intermediatePointsAtT(t),
          r = new Bezier(q).pointAtT(t);
      return [new Bezier(p[0], pp[0], r), new Bezier(r, pp[1], p[2])];
    case 3:
      var q = this.intermediatePointsAtT(t),
          r = new Bezier(q).intermediatePointsAtT(t),
          s = new Bezier(r).pointAtT(t);
      return [new Bezier(p[0], q[0], r[0], s), new Bezier(s, r[1], q[2], p[3])];
    }
  },
  calcCurveLength: function calcCurveLength() {
  //    if (this.inflectionPointTs().length > 0) {
  //      var segments = this.subdivideAtInfectionPoints();
  //      var length = 0;
  //      for (var i = 0, n = segments.length; i < n; ++i)
  //        length += segments[i].calcCurveLength();
  //      return length;
  //    }

    var segments = this.segments();
    return segments[segments.length - 1].accLen;
  },
  calcCurveLength2: function calcCurveLength2() {
    var self = this;
    function f(t) {
      return self.derivativeAtT(t).length();
    }
    return new numeric.integral(f, 0, 1, Bezier.epsilon);
  },
  getTAtLength: function getTAtLength(length) {
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
  },
  segments: function segments() {
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

      var n = segments.length;
      var h = 1 / (n - 1);
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
  }
});

// 1. 2本の曲線をX軸またはY軸に平行な点でそれぞれ分割する。
// 2. 分割された曲線のバウンディングボックスが重なる組み合わせを見つける。
// 3. 分割された区間で曲線が直線と見なせるか判定する。
//   3-1. 区間内で曲線の傾きがバウンディングボックスの対角線の傾きと同じに
//        なる点を探す。そのような点は1つまたは2つ存在する。
//   3-2. その点とバウンディングボックスの対角線との距離を調べる。
//   3-3. その距離が指定の誤差範囲内なら直線と見なす。
// 4. バウンディングボックスが重なる組み合わせのうち、直線と見なせないほう
//    （片方または両方）を2分割して2.に戻る。
// 5. バウンディングボックスが重なる組み合わせの両方が直線と見なせるようなら
//    対角線の2線分の交点を求めればそれが元の2曲線の交点。
function BezierIntersectionBBoxAlgorithm(
  bezier1, bezier2, pointTolerance, lineTolerance) {
  this.bezier1 = bezier1;
  this.bezier2 = bezier2;
  this.pointTolerance = pointTolerance;
  this.lineTolerance = lineTolerance;
  this.points = [];
}
mix(BezierIntersectionBBoxAlgorithm, {
  execute: function execute(bezier1, bezier2, pointTolerance, lineTolerance) {
    return new BezierIntersectionBBoxAlgorithm(
      bezier1, bezier2, pointTolerance, lineTolerance).execute();
  }
});
mix(BezierIntersectionBBoxAlgorithm.prototype, {
  execute: function execute() {
    var segmentPairs = this.calcInitialSegmentPairs();
    this.iterationCount = 0;
    while (segmentPairs.length > 0) {
      segmentPairs = this.processSegmentPairs(segmentPairs);
      ++this.iterationCount;
    }
console.log('iterationCount=' + this.iterationCount);
    return this.points;
  },
  calcInitialSegmentPairs: function calcInitialSegmentPairs() {
    var bezier1 = this.bezier1, bezier2 = this.bezier2;
    var segments1 =
      bezier1.segmentsDividedAtTs(bezier1.tsOfZeroXYDerivPoints());
    var segments2 =
      bezier2.segmentsDividedAtTs(bezier2.tsOfZeroXYDerivPoints());
    return this.pairsOfSegmentsWithOverlappedBbox(segments1, segments2);
  },
  findIntersectionsOrDivideSegments:
  function findIntersectionsOrDivideSegments(segmentPair) {
    var segment1 = segmentPair[0], segment2 = segmentPair[1];
    var isLine1 = segment1.isAssumableAsLine(this.lineTolerance);
    var isLine2 = segment2.isAssumableAsLine(this.lineTolerance);
    if (isLine1 && isLine2) {
      var point = segment1.toLine().intersection(segment2.toLine());
      if (point && !point.equalsToOneOf(this.points, this.pointTolerance))
        this.points.push(point);
      return [];
    }
    else {
      var segments1 = isLine1 ? [segment1] : segment1.divideAtMidT();
      var segments2 = isLine2 ? [segment2] : segment2.divideAtMidT();
      return this.pairsOfSegmentsWithOverlappedBbox(segments1, segments2);
    }
  },
  processSegmentPairs: function processSegmentPairs(segmentPairs) {
    var newSegmentPairs = [];
    for (var i = 0, n = segmentPairs.length; i < n; ++i) {
      newSegmentPairs = newSegmentPairs.concat(
        this.findIntersectionsOrDivideSegments(segmentPairs[i])
      );
    }
    return newSegmentPairs;
  },
  pairsOfSegmentsWithOverlappedBbox:
  function pairsOfSegmentsWithOverlappedBbox(segments1, segments2) {
    var n1, n2, i1, i2, segmentPairs, segment1, segment2;
    n1 = segments1.length;
    n2 = segments2.length;
    segmentPairs = [];
    for (i1 = 0; i1 < n1; ++i1) {
      segment1 = segments1[i1];
      for (i2 = 0; i2 < n2; ++i2) {
        segment2 = segments2[i2];
        if (segment1.toLine().bbox().intersects(segment2.toLine().bbox())) {
          segmentPairs.push([segment1, segment2]);
        }
      }
    }
    return segmentPairs;
  }
});

function BezierSegment(bezier, t1, t2, p1, p2) {
  if (!(0 <= t1 && t1 <= 1))
    throw new Error('t1 must be between 0 and 1.');
  if (!(0 <= t2 && t2 <= 1))
    throw new Error('t2 must be between 0 and 1.');
  if (!(t1 < t2))
    throw new Error('t1 must be less than t2.');

  this.bezier = bezier;
  this.t1 = t1;
  this.t2 = t2;
  this.p1 = p1 || bezier.pointAtT(t1);
  this.p2 = p2 || bezier.pointAtT(t2);
}
mix(BezierSegment.prototype, {
  divideAtMidT: function divideAtMidT() {
    var bezier = this.bezier, t1 = this.t1, t2 = this.t2;
    var tMid = (t1 + t2) / 2;
    var pMid = bezier.pointAtT(tMid);
    return [
      new BezierSegment(bezier, t1, tMid, this.p1, pMid),
      new BezierSegment(bezier, tMid, t2, pMid, this.p2)
    ];
  },
  getTsAtSameGradientAsLine: function getTsAtSameGradientAsLine() {
    var cv = this.bezier.derivativeCoefficients();
    var c0 = cv[0];
    var c1 = cv[1];
    var c2 = cv[2] || Vector2d.ZERO;
    var dv = this.toLine().dv();
    var dx = dv.x, dy = dv.y;
    var t1 = this.t1, t2 = this.t2;
    return realRootsOfQuadraticEquation(
      c2.x * dy - c2.y * dx,
      c1.x * dy - c1.y * dx,
      c0.x * dy - c0.y * dx,
      function(t) { return t1 <= t && t <= t2; }
    );
  },
  isAssumableAsLine: function isAssumableAsLine(maxDistance) {
    var bezier = this.bezier;
    switch (bezier.degree()) {
    case 1:
      return true;
    case 2:
    case 3:
      var ts = this.getTsAtSameGradientAsLine();
      var n = ts.length;
      if (n === 0)
        throw new Error('Unexpected error');
      var dist = 0;
      var line = this.toLine();
      for (var i = 0; i < n; ++i) {
        var t = ts[i], d = line.distance(bezier.pointAtT(t));
        dist = Math.max(dist, d);
      }
      return dist <= maxDistance;
    }
  },
  toBezier: function toLine() {
    if (this.t1 === 0) {
      if (this.t2 === 1)
        return this.bezier;
      else
        return this.bezier.beziersDividedAtT(t2)[0];
    }
    else {
      if (this.t2 === 1)
        return this.bezier.beziersDividedAtT(t1)[1];
      else {
        var beziers = this.bezier.beziersDividedAtT(this.t1);
        return beziers[0].beziersDividedAtT(this.t2 / this.t1)[1];
      }
    }
  },
  toLine: function toLine() {
    return new Line(this.p1, this.p2);
  }
});

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

function realRootsOfQuadraticEquation(a, b, c, filterFunc) {
  return new numeric.QuadraticEquation(a, b, c).realRoots(filterFunc);
}

function isUndefined(obj) {
  return obj === undefined;
}

function mix(dest, src) {
  for (var k in src)
    dest[k] = src[k];
  return dest;
}

return {
  Vector2d: Vector2d,
  Line: Line,
  Rectangle: Rectangle,
  QuadraticBezier: QuadraticBezier,
  Bezier: Bezier,
  BezierSegment: BezierSegment,
  BezierIntersectionBBoxAlgorithm: BezierIntersectionBBoxAlgorithm
};
  
}();
