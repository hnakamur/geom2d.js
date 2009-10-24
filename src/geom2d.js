var geom2d = function() {

var sum = numeric.sum, numberEquals = numeric.numberEquals,
    QuadraticEquation = numeric.QuadraticEquation;

// immutable two dimensional vector.
function Vector2d(x, y) {
  this.x = x;
  this.y = y;
}
Vector2d.ZERO = new Vector2d(0, 0);
Vector2d.vectorsFromNestedXYArray = function(xyArray) {
  var vectors = [], n = xyArray.length, i, xy;
  for (i = 0; i < n; ++i) {
    xy = xyArray[i];
    vectors.push(new Vector2d(xy[0], xy[1]));
  }
  return vectors;
};
Vector2d.sum = function(vectors) {
  var xs = [], ys = [], n = vectors.length, i, v;
  for (i = 0; i < n; ++i) {
    v = vectors[i];
    xs.push(v.x);
    ys.push(v.y);
  }
  return new Vector2d(sum(xs), sum(ys));
};
Vector2d.polynomial = function(t, coefficients) {
  var xs = [], ys = [], n = coefficients.length, ti = 1, i, v;
  for (i = 0; i < n; ++i, ti *= t) {
    v = coefficients[i];
    xs.push(v.x * ti);
    ys.push(v.y * ti);
  }
  return new Vector2d(sum(xs), sum(ys));
};
proto = Vector2d.prototype;
proto.toString = function() {
  return '[' + this.x + ' ' + this.y + ']';
};
proto.unit = function() {
  return this.scalarDiv(this.length());
};
proto.length = function() {
  return Math.sqrt(this.dotProduct(this));
};
proto.dotProduct = function(vector) {
  return sum([this.x * vector.x, this.y * vector.y]);
};
proto.zComponentOfCrossProduct = function(vector) {
  return this.x * vector.y - this.y * vector.x;
};
proto.equals = function(vector, epsilon) {
  return numberEquals(this.x, vector.x, epsilon) &&
    numberEquals(this.y, vector.y, epsilon);
};
proto.equalsToOneOf = function(vectors, epsilon) {
  for (var i = 0, n = vectors.length; i < n; ++i) {
    if (this.equals(vectors[i], epsilon))
      return true;
  }
  return false;
};
proto.add = function(vector) {
  return new Vector2d(this.x + vector.x, this.y + vector.y);
};
proto.subtract = function(vector) {
  return new Vector2d(this.x - vector.x, this.y - vector.y);
};
proto.negate = function() {
  return new Vector2d(-this.x, -this.y);
};
proto.scalarMult = function(factor) {
  return new Vector2d(this.x * factor, this.y * factor);
};
proto.scalarDiv = function(factor) {
  return new Vector2d(this.x / factor, this.y / factor);
};

function Line(p0, p1) {
  this.p0 = p0;
  this.p1 = p1;
}
Line.intersection = function(line1, line2) {
  var params = this.paramsAtIntersection(line1, line2);
  return params ? line1.pointAtT(params[0]) : params;
};
Line.paramsAtIntersection = function(line1, line2) {
  var dp0 = line2.p0.subtract(line1.p0), dv1 = line1.dv(), dv2 = line2.dv(),
      a = dv1.x, b = -dv2.x, c = dv1.y, d = -dv2.y, e = dp0.x, f = dp0.y,
      det = a * d - b * c, t, u;
  if (det) {
    t = (d * e - b * f) / det, u = (a * f - c * e) / det;
//console.log('direct t=' + t + ', u=' + u);
//var p = line1.pointAtT(t), q = line2.pointAtT(u);
//console.log('direct: p[t]=' + p + ', q[u]=' + q + ', difLen=' + p.subtract(q).length());
    return [t, u];
  }
  else
    return undefined;
}
Line.paramsAtIntersection2 = function(line1, line2) {
  var dp0 = line2.p0.subtract(line1.p0), dv1 = line1.dv(), dv2 = line2.dv(),
      matrix = new numeric.Matrix([
        [dv1.x, -dv2.x],
        [dv1.y, -dv2.y]]),
      vector = new numeric.Vector(dp0.x, dp0.y),
      ge = new numeric.GaussElimination(matrix, vector);
console.log('dv1=' + dv1 + ', dv2=' + dv2);
console.log('matrix=' + matrix + ', vector=' + vector);
var a = matrix.elements;
console.log('det=' + numeric.Matrix.det2by2(a[0][0], a[0][1], a[1][0], a[1][1]));
  if (ge.isSingular())
    return undefined;
  var sol = ge.solution(), t = sol.component(0), u = sol.component(1);
console.log('t=' + t + ', u=' + u);
console.log('point=' + matrix.multiply(sol) + ', difLen=' + matrix.multiply(sol).subtract(vector).length());
var p = line1.pointAtT(t), q = line2.pointAtT(u);
console.log('p[t]=' + p + ', q[u]=' + q + ', difLen=' + p.subtract(q).length());

  var dcmp = new numeric.LUDecomposition(matrix);
  sol = dcmp.solve(vector);
  t = sol.component(0), u = sol.component(1);
console.log('LUDecomp t=' + t + ', u=' + u);
console.log('LU: point=' + matrix.multiply(sol) + ', difLen=' + matrix.multiply(sol).subtract(vector).length());
var p = line1.pointAtT(t), q = line2.pointAtT(u);
console.log('LU: p[t]=' + p + ', q[u]=' + q + ', difLen=' + p.subtract(q).length());

  //if (0 <= t && t <= 1 & 0 <= u && u <= 1) {
    return [t, u];
  //}
  //else
  //  return undefined;
};
var proto = Line.prototype;
proto.dv = function() {
  if (!this._dv)
    this._dv = this.p1.subtract(this.p0);
  return this._dv;
}
proto.length = function() {
  return this.dv().length();
};
proto.pointAtT = function(t) {
  return Vector2d.polynomial(t, [this.p0, this.dv()]);
};
proto.toImplicitLine = function() {
  var p0 = this.p0, p1 = this.p1, x0 = p0.x, y0 = p0.y, x1 = p1.x, y1 = p1.y;
  return new ImplicitLine(y1 - y0, x0 - x1, x1 * y0 - x0 * y1);
};
proto.distance = function(point) {
  return Math.abs(this.signedDistance(point));
};
proto.signedDistance = function(point) {
  return point.subtract(this.p0).zComponentOfCrossProduct(this.dv().unit());
};

// Implicit line ax + by + c = 0
function ImplicitLine(a, b, c) {
  this.a = a;
  this.b = b;
  this.c = c;
}
proto = ImplicitLine.prototype;
proto.distance = function(point) {
  return sum([this.a * point.x, this.b * point.y, this.c]);
}

function Rectangle(x, y, width, height) {
  this.x = x;
  this.y = y;
  this.width = width;
  this.height = height;
}
Rectangle.fromTwoCornerPoints = function(points) {
  var p0 = points[0], p1 = points[1],
      x0 = p0.x, y0 = p0.y, x1 = p1.x, y1 = p1.y;
  return new Rectangle(Math.min(x0, x1), Math.min(y0, y1),
    Math.abs(x1 - x0), Math.abs(y1 - y0));
}
var proto = Rectangle.prototype;
proto.intersects = function(rect) {
  return this.x <= rect.x + rect.width && this.x + this.width >= rect.x &&
    this.y <= rect.y + rect.height && this.y + this.height >= rect.y;
};

function Bezier(points) {
  this.points = points;
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

//
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

Bezier.intersections = function(bezier1, bezier2, tolerance) {
  var segments1 =
        Bezier.intersections.segmentsDividedByZeroXYDerivPoints(bezier1),
      segments2 =
        Bezier.intersections.segmentsDividedByZeroXYDerivPoints(bezier2),
      segmentPairs =
        Bezier.intersections.findSegmentPairs(segments1, segments2),
      points = [];
console.log('initial segmentPairs:');
console.log(segmentPairs);
  for (var i = 0, n = segmentPairs.length; i < n; ++i)
    findIntersectionInSegmentPair(segmentPairs[i]);

  function findIntersectionInSegmentPair(segmentPair) {
console.log('findIntersectionInSegmentPair start. segmentPair:');
console.log(segmentPair);
    var segment1 = segmentPair[0], segment2 = segmentPair[1],
        isLine1 = segment1.isAssumableAsLine(tolerance),
        isLine2 = segment2.isAssumableAsLine(tolerance)
    if (isLine1 && isLine2) {
      var point = Line.intersection(segment1.toLine(), segment2.toLine());
console.log('point');
console.log(point);
      if (point && !point.equalsToOneOf(points))
        points.push(point);
    }
    else {
      var segments1 = isLine1 ? [segment1] : segment1.divideInHalves(),
          segments2 = isLine2 ? [segment2] : segment2.divideInHalves(),
          segmentPairs =
            Bezier.intersections.findSegmentPairs(segments1, segments2);
      for (var i = 0, n = segmentPairs.length; i < n; ++i)
        findIntersectionInSegmentPair(segmentPairs[i]);
    }
  }
console.log('points');
console.log(points);
  return points;
};
Bezier.intersectionSegmentPairs = function(bezier1, bezier2) {
  var segments1, segments2, segmentPairs;
  segments1 = Bezier.intersections.segmentsDividedByZeroXYDerivPoints(bezier1);
  segments2 = Bezier.intersections.segmentsDividedByZeroXYDerivPoints(bezier2);
  segmentPairs = Bezier.intersections.findSegmentPairs(segments1, segments2);
//segmentPairs = segmentPairs.slice(0, 1);
  for (var i = 0; i < 1; ++i)
    segmentPairs = Bezier.intersections.narrowSegmentPairs(bezier1, bezier2, segmentPairs);
  return segmentPairs;
};
Bezier.intersections.segmentsDividedByZeroXYDerivPoints = function(bezier) {
  var cv = bezier.derivativeCoefficients(), c0 = cv[0], c1 = cv[1],
      c2 = cv[2] || Vector2d.ZERO, ts = [0, 1],
      filterFunc = function (t) { return 0 < t && t < 1; },
      ts, segments;
  ts = numeric.uniqAndSort(
    [0, 1].concat(
      new QuadraticEquation(c2.x, c1.x, c0.x).realRoots(filterFunc),
      new QuadraticEquation(c2.y, c1.y, c0.y).realRoots(filterFunc))
  );
//console.log('ts');
//console.log(ts);

  segments = [];
  for (i = 0, n = ts.length - 1; i < n; ++i) {
    var tStart = ts[i], tEnd = ts[i + 1],
        pStart = i == 0 ? bezier.pointAtT(tStart) : segments[i - 1].pEnd,
        pEnd = bezier.pointAtT(tEnd)
    segments.push(
      new BezierMonotonousSegment(bezier, tStart, tEnd, pStart, pEnd));
  }
//console.log('segments');
//console.log(segments);
  return segments;
};
Bezier.intersections.findSegmentPairs = function(segments1, segments2) {
  var n1, n2, i1, i2, segmentPairs, segment1, segment2;
  n1 = segments1.length;
  n2 = segments2.length;
  segmentPairs = [];
  for (i1 = 0; i1 < n1; ++i1) {
    segment1 = segments1[i1];
    for (i2 = 0; i2 < n2; ++i2) {
      segment2 = segments2[i2];
      if (segment1.box().intersects(segment2.box())) {
        segmentPairs.push([segment1, segment2]);
      }
    }
  }
  return segmentPairs;
};
Bezier.intersections.segmentLength = function(segment) {
  return new Line(segment.pStart.x, segment.pStart.y, segment.pEnd.x, segment.pEnd.y).length();
};
Bezier.intersections.narrowSegmentPairs = function(bezier1, bezier2, segmentPairs) {
  var newSegmentPairs = [];
  for (var i = 0, n = segmentPairs.length; i < n; ++i) {
    var segmentPair = segmentPairs[i];
    var segment1 = segmentPair[0];
    var segment2 = segmentPair[1];
    var length1 = Bezier.intersections.segmentLength(segment1);
    var length2 = Bezier.intersections.segmentLength(segment2);
    var segments1, segments2;
    if (length1 > length2) {
      segments1 = Bezier.intersections.divideSegment(bezier1, segment1);
      segments2 = [segment2];
    }
    else {
      segments1 = [segment1];
      segments2 = Bezier.intersections.divideSegment(bezier2, segment2);
    }
    newSegmentPairs = newSegmentPairs.concat(
      Bezier.intersections.findSegmentPairs(segments1, segments2));
  }
  return newSegmentPairs;
};
Bezier.intersections.findIntersectionInSegmentPair = function(bezier1, bezier2, segmentPair) {
  var point = Bezier.intersections.intersectionOfSegmentPairLines(segmentPair);
  var i = 1;
  while (true) {
    var segments1 = Bezier.intersections.divideSegment(bezier1, segmentPair[0]);
    var segments2 = Bezier.intersections.divideSegment(bezier2, segmentPair[1]);
    var segmentPairs = Bezier.intersections.findSegmentPairs(segments1, segments2);

    var newPoints = [];
    var newSegmentPairs = [];
    for (var i = 0, n = segmentPairs.length; i < n; ++i) {
      segmentPair = segmentPairs[i];
      var newPoint = Bezier.intersections.intersectionOfSegmentPairLines(segmentPair);
      if (newPoint) {
        newSegmentPairs.push(segmentPair);
        newPoints.push(newPoint);
      }
    }
    if (newSegmentPairs.length > 1) {
console.log('segmentPairs.length=' + newSegmentPairs.length);
console.log('segmentPairs');
console.log(newSegmentPairs);
      throw new Error('Unexpected segment pair count = ' + newSegmentPairs.length);
    }

    if (newPoints.length === 0)
      return undefined;
    if (newPoints[0].equals(point, 1e-10)) {
console.log('iteration count=' + i);
      return point;
    }

    point = newPoints[0];
    segmentPair = newSegmentPairs[0];
    ++i;
  }

};
Bezier.intersections.intersectionOfSegmentPairLines = function(segmentPair) {
  var line1 = new Line(segmentPair[0].pStart.x, segmentPair[0].pStart.y,
    segmentPair[0].pEnd.x, segmentPair[0].pEnd.y);
  var line2 = new Line(segmentPair[1].pStart.x, segmentPair[1].pStart.y,
    segmentPair[1].pEnd.x, segmentPair[1].pEnd.y);
  return Line.intersection(line1, line2);
};
Bezier.intersections.divideSegment = function(bezier, segment) {
  var tStart = segment.tStart,
      tEnd = segment.tEnd,
      pStart = segment.pStart,
      pEnd = segment.pEnd,
      tMid = (tStart + tEnd) / 2,
      pMid = bezier.pointAtT(tMid),
      box0 = Rectangle.fromTwoCornerPoints([pStart, pMid]),
      box1 = Rectangle.fromTwoCornerPoints([pMid, pEnd]);
  return [
    {tStart: tStart, tEnd: tMid, pStart: pStart, pEnd: pMid, box: box0},
    {tStart: tMid, tEnd: tEnd, pStart: pMid, pEnd: pEnd, box: box1}
  ];
};
//Bezier.intersections.calcBoxes = function(bezier, segments) {
//  var n = segments.length, i, segment, t, box;
//  for (i = 0; i < n; ++i) {
//    segment = segments[i];
//    t = segment.t;
//    if (!segment.p)
//      segment.p = bezier.pointAtT(t);
//  }
//
//  for (i = 0; i < n - 1; ++i) {
//    segment = segments[i];
//    if (!segment.box)
//      segment.box = Rectangle.fromTwoCornerPoints(segment.p, segments[i + 1].p);
//  }
//};
var proto = Bezier.prototype;
proto.degree = function() {
  return this.points.length - 1;
};
proto.pointAtT = function(t) {
  return Vector2d.polynomial(t, this.coefficients());
};
proto.pointAtT2 = function(t) {
  var n = this.degree(), t1 = 1 - t, p = this.points, q = [], i, j;
  for (i = n; i >= 1; --i) {
    for (j = 0; j < i; ++j)
      q[j] = p[j].scalarMult(t1).add(p[j + 1].scalarMult(t));
    p = q;
    q = [];
  }
  return Point.fromVector(p[0]);
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
  return Vector2d.polynomial(t, this.derivativeCoefficients());
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
      this._derivativeCoefficients =
        [c[1], c[2].scalarMult(2), c[3].scalarMult(3)];
      break;
    }
  }
  return this._derivativeCoefficients;
};
proto.secondDerivativeAtT = function(t) {
  return Vector2d.polynomial(t, this.secondDerivativeCoefficients());
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
      this._secondDerivativeCoefficients =
        [c[2].scalarMult(2), c[3].scalarMult(6)];
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
      roots = new QuadraticEquation(
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
proto.segmentsDividedByZeroXYDerivPoints = function() {
  var cp = this.derivativeCoefficients(), ts = [0, 1], segments = [], i;
  for (i = 0; i < 2; ++i) {
    var eqn = new QuadraticEquation(
          cp[2].component(i), cp[1].component(i), cp[0].component(i)),
        roots = filterValues(eqn.realRoots(),
          function(t) { return 0 < t && t < 1});
    ts = ts.concat(roots);
  }
  var ts = numeric.uniqAndSort(ts);
  for (var i = 0, n = ts.length - 1; i < n; ++i) {
//console.log('i=' + i + ', ts[i]=' + ts[i] + ', ts[i+1]=' + ts[i + 1]);
    segments.push(new BezierMonotonousSegment(this, ts[i], ts[i + 1]));
  }
  return segments;
};

function BezierMonotonousSegment(bezier, t1, t2, p1, p2) {
  this.bezier = bezier;
  this.t1 = t1;
  this.t2 = t2;
  this.p1 = p1 || bezier.pointAtT(t1);
  this.p2 = p2 || bezier.pointAtT(t2);
}
proto = BezierMonotonousSegment.prototype;
proto.box = function() {
  return Rectangle.fromTwoCornerPoints([this.p1, this.p2]);
};
proto.toLine = function() {
  return new Line(this.p1, this.p2);
};
proto.isAssumableAsLine = function(maxDistance) {
  var bezier = this.bezier;
  switch (bezier.degree()) {
  case 1:
    return true;
  case 2:
  case 3:
    var ts = this.getTsAtSameGradientAsLine(), n = ts.length,
        line = this.toLine(), i;
    if (n === 0)
      throw new Error('Unexpected error');
    dist = 0;
    for (i = 0; i < n; ++i) {
      var t = ts[i];
      var d = line.distance(bezier.pointAtT(ts[i]));
//console.log('t=' + t + ', d=' + d);
      dist = Math.max(dist, d);
    }
    return dist <= maxDistance;
  }
};
proto.getTsAtSameGradientAsLine = function() {
  var cv = this.bezier.derivativeCoefficients(), c0 = cv[0], c1 = cv[1],
      c2 = cv[2] || Vector2d.ZERO,
      dv = this.toLine().dv(), dx = dv.x, dy = dv.y,
      t1 = this.t1, t2 = this.t2;
  return new QuadraticEquation(
    c2.x * dy - c2.y * dx,
    c1.x * dy - c1.y * dx,
    c0.x * dy - c0.y * dx
  ).realRoots(function(t) { return t1 <= t && t <= t2; });
};
proto.divideInHalves = function() {
  var bezier = this.bezier, t1 = this.t1, t2 = this.t2,
      p1 = this.p1, p2 = this.p2,
      tMid = (t1 + t2) / 2, pMid = bezier.pointAtT(tMid);
  return [
    new BezierMonotonousSegment(bezier, t1, tMid, p1, pMid),
    new BezierMonotonousSegment(bezier, tMid, t2, pMid, p2)
  ];
};

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
  Vector2d: Vector2d,
  Line: Line,
  Rectangle: Rectangle,
  Bezier: Bezier,
  BezierMonotonousSegment: BezierMonotonousSegment
};
  
}();
