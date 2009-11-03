var numeric = function() {

// References:
// [1] What Every Computer Scientist Should Know About Floating-Point Arithmetic
//     http://docs.sun.com/source/806-3568/ncg_goldberg.html
// [2] Numerical Recipe 3rd ed. (ISBN: 978-0-521-88068-8)
//     NOTE: To avoid copyright violation, I never looked at source codes in
//           this book. I read only explanations of algorithms.
//           All implementations in this source code are all my own works.

var MACHINE_EPSILON = 2.220446049250313e-16;

// a * x^2 + b * x + c = 0
function QuadraticEquation(a, b, c) {
  this.a = a;
  this.b = b;
  this.c = c;
}
mix(QuadraticEquation.prototype, {
  valueAt: function valueAt(x) {
    return sum([this.a * x * x, this.b * x, this.c]);
    //return (this.a * x + this.b) * x + this.c; <- this is unstable.
  },
  realRoots: function realRoots(predicate) {
    // http://en.wikipedia.org/wiki/Quadratic_equation
    var a = this.a, b = this.b, c = this.c;
    var roots = [];
    if (a !== 0) {
      var dis = b * b - 4 * a * c;
      if (dis > 0) {
        var t = -0.5 * (b + sgn(b) * Math.sqrt(dis)), x1 = t / a, x2 = c / t;
        roots = x1 < x2 ? [x1, x2] : [x2, x1];
      }
      else if (dis === 0)
        roots = [-0.5 * b / a];
      else
        roots = [];
    }
    else if (b !== 0)
      roots = [-c / b];
    else if (c === 0)
      throw new Error('Infinite number of roots exist (Any value is ok).');

    return predicate ? filterValues(roots, predicate) : roots;
  }
});

// x^3 + a x^2 + b x + c = 0
function CubicEquation(a, b, c) {
  this.a = a;
  this.b = b;
  this.c = c;
}
mix(CubicEquation.prototype, {
  derivativeValueAt: function valueAt(x) {
    return sum([3 * x * x, 2 * this.a * x, this.b]);
    //return (3 * x + 2 * this.a) * x + this.b;
  },
  realRoots: function realRoots(predicate) {
    // see p.228 of Numerical Recipe 3rd ed (ISBN: 978-0-521-88068-8)
    var a = this.a, b = this.b, c = this.c;
    var q = (a * a - 3 * b) / 9;
    var r = (2 * a * a * a - 9 * a * b + 27 * c) / 54;
    var r2 = r * r;
    var q3 = q * q * q;
    var adiv3 = a / 3;
    var roots;
    if (r2 < q3) {
      var theta = Math.acos(r / Math.sqrt(q3));
      var m2sqrtq = -2 * Math.sqrt(q);
      roots = [
        m2sqrtq * Math.cos(theta / 3) - adiv3,
        m2sqrtq * Math.cos((theta + 2 * Math.PI) / 3) - adiv3,
        m2sqrtq * Math.cos((theta - 2 * Math.PI) / 3) - adiv3
      ];
    }
    else {
      var sgnr = r > 0 ? 1 : -1;
      var la = -sgnr * Math.pow(Math.abs(r) + Math.sqrt(r2 - q3), 1 / 3);
      var lb = la !== 0 ? q / la : 0;
      roots = [(la + lb) - adiv3];
      if (la === lb)
        roots.push(-(la + lb) / 2 - adiv3);
    }

    roots = uniqAndSort(roots);

    var self = this;
    var newtonMethodEpsilon = 1e-12;
    roots = new NewtonRaphsonMethod(
      function(x) { return self.valueAt(x); },
      function(x) { return self.derivativeValueAt(x); }
    ).polishRoots(roots, newtonMethodEpsilon);

    return predicate ? filterValues(roots, predicate) : roots;
  },
  valueAt: function valueAt(x) {
    var x2 = x * x, x3 = x2 * x;
    return sum([x3, this.a * x2, this.b * x, this.c]);
    //return ((x + this.a) * x + this.b) * x + this.c;
  }
});

// sum(c[i] * t^i) (i=0..n)
function Polynomial(coefficients) {
  var n = coefficients.length - 1;
  while (n >= 0 && coefficients[n] === 0)
    --n;
  this.coefficients = coefficients.slice(0, n + 1);
}
Polynomial.resultant = function(f, g) {
  var degree = Math.max(f.degree(), g.degree()),
      a = f.coefficients, b = g.coefficients;
  function ab(i, j) {
    return (a[i] || 0) * (b[j] || 0) - (a[j] || 0) * (b[i] || 0);
  }

  switch (degree) {
  case 1:
    return ab(1, 0);
  case 2:
    var a2b1 = ab(2, 1), a2b0 = ab(2, 0),
        a1b0 = ab(1, 0);
    return Matrix.det2by2(a2b1, a2b0, a2b0, a1b0);
  case 3:
    var a3b2 = ab(3, 2), a3b1 = ab(3, 1), a3b0 = ab(3, 0),
        a2b1 = ab(2, 1), a2b0 = ab(2, 0),
        a1b0 = ab(1, 0);
    return Matrix.det3by3(
      a3b2, a3b1, a3b0,
      a3b1, a3b0 + a2b1, a2b0,
      a3b0, a2b0, a1b0
    );
  default:
    throw new Error('Only polynomials of 3 or lower degree are supported.');
  }
}
mix(Polynomial.prototype, {
  degree: function degree() {
    return this.coefficients.length - 1;
  },
  derivative: function derivative() {
    var c = this.coefficients, dc = [];
    for (i = 1, degree = this.degree(); i <= degree; ++i)
      dc.push(c[i] * i);
//console.log('derivative: dc=' + JSON.stringify(dc));
    return new Polynomial(dc);
  },
  realRootsBetween: function realRootsBetween(xMin, xMax, epsilon) {
    if (this.degree() <= 2) {
      var c = this.coefficients;
      return new QuadraticEquation(c[2], c[1], c[0]).realRoots(
        function(x) { return xMin <= x && x <= xMax; });
    }
    else {
      var deriv = this.derivative();
      var derivRoots = deriv.realRootsBetween(xMin, xMax, epsilon);
//console.log('derivRoots=' + JSON.stringify(derivRoots));
      var xs = uniqAndSort([xMin, xMax].concat(derivRoots));
      var n = xs.length;
      var sgnYs = [];
      var roots = [];
      for (var i = 0; i < n; ++i) {
        var x = xs[i];
        var sgnY = sgnYs[i] = sgn(this.valueAt(x));
        if (sgnY === 0)
          roots.push(x);
      }

      var f = bind(this.valueAt, this);
      var df = bind(deriv.valueAt, deriv);
      for (var i = 0; i < n - 1; ++i) {
        if (sgnYs[i] * sgnYs[i + 1] === -1) {
          var finder = new OneRootFinder(f, df);
          var root = finder.findOneRootBetween(xs[i], xs[i + 1], epsilon);
//console.log('findOneRoot iter=' + finder.iter);
          roots.push(root);
        }
      }

      return uniqAndSort(roots);
    }
  },
  valueAt: function valueAt(t) {
    var terms = [], c = this.coefficients, ti = 1;
    for (var i = 0, degree = this.degree(); i <= degree; ++i) {
      terms.push(c[i] * ti);
      ti *= t;
    }
    return sum(terms);

    //var s = 0, c = this.coefficients, i = c.length;
    //while (--i >= 0) {
    //  s = s * t + c[i];
    //}
    //return s;
  }
});

function OneRootFinder(func, derivFunc) {
  this.func = func;
  this.derivFunc = derivFunc;
}
mix(OneRootFinder.prototype, {
  calcEpsilonFromMinAndMax: function calcEpsilonFromMinAndMax(xMin, xMax) {
    return MACHINE_EPSILON * (Math.abs(xMin) + Math.abs(xMax)) / 2;
    //return MACHINE_EPSILON * (Math.abs(this.func(xMin)) + Math.abs(this.func(xMax))) / 2;
  },
  findOneRootBetween: function findOneRootBetween(xMin, xMax, epsilon) {
//console.log('findOneRootBetween(): epsilon=' + epsilon);
    var isEpsProvided = epsilon !== undefined;
    if (!isEpsProvided) {
      epsilon = this.calcEpsilonFromMinAndMax(xMin, xMax);
//console.log('findOneRootBetween(): calculated epsilon=' + epsilon);
    }
    var yAtXMin = this.func(xMin);
    if (numberEquals(yAtXMin, 0, epsilon))
      return xMin;
    var yAtXMax = this.func(xMax);
    if (numberEquals(yAtXMax, 0, epsilon))
      return xMax;
//console.log('findOneRootBetween(): yAtXMin=' + yAtXMin + ', yAtXMax=' + yAtXMax);
    var sgnYAtXMin = sgn(yAtXMin), sgnYAtXMax = sgn(yAtXMax);
    if (sgnYAtXMin * sgnYAtXMax >= 0)
      throw new Error('Values at xMin and xMax must have different signs.');

    var xNewton, yAtXNewton;
    for (this.iter = 0; this.iter < this.maxIteration; ++this.iter) {
//console.log('findOneRootBetween(): xMin=' + xMin + ', xMax=' + xMax);
      // First, we proceed one step with the bisection method.
      var xMid = (xMin + xMax) / 2, yAtXMid = this.func(xMid);
//console.log('findOneRootBetween(): yAtXMid=' + yAtXMid + ', xMid=' + xMid);
      if (numberEquals(yAtXMid, 0, epsilon))
        return xMid;
      if (xMid === xMin || xMid === xMax) {
        if (isEpsProvided)
          throw new Error('Specified epsilon is too small for this root.');
        else
          return Math.abs(yAtXMin) < Math.abs(yAtXMax) ? xMin : xMax;
      }

      if (sgn(yAtXMid) === sgnYAtXMin)
        xMin = xMid;
      else
        xMax = xMid;
//console.log('findOneRootBetween(): after bisection xMin=' + xMin + ', xMax=' + xMax);

      // We use the middle point as the starting point of the Newton method
      // in cases:
      // 1) in the first iteration, 2) the previous result with the Newton
      // method was out of the range determined by the bisection method.
      //
      // Otherwise, the previous result with the Newton method was in the
      // range, we use that point.
      if (xNewton === undefined) {
        xNewton = xMid;
        yAtXNewton = yAtXMid;
      }

      // Now, we proceed one step with the Newton method. If the result
      // is in the range determined with the bisection, we use it to narrow
      // the range. Otherwise, we discard it.
      //
      // Note we do not need to check the derivative value is zero.
      // In that case, the value of the division below will be one of NaN,
      // Infinity, and -Infinity. The result with the Newton method will be
      // out of the range determined with the bisection, and hence will
      // be discarded.
      xNewton -= yAtXNewton / this.derivFunc(xNewton);
//console.log('findOneRootBetween(): xNewton=' + xNewton);
      if (xMin <= xNewton && xNewton <= xMax) {
        yAtXNewton = this.func(xNewton);
//console.log('findOneRootBetween(): yAtXNewton=' + yAtXNewton);
        if (numberEquals(yAtXNewton, 0, epsilon))
          return xNewton;

        if (sgn(yAtXNewton) == sgnYAtXMin)
          xMin = xNewton;
        else
          xMax = xNewton;
//console.log('findOneRootBetween(): after newton xMin=' + xMin + ', xMax=' + xMax);
      }
      else
        xNewton = undefined;
    }
    throw new Error('Maximum number of iterations exceeded.');
  },
  maxIteration: 100
});

function NewtonRaphsonMethod(func, derivFunc) {
  this.func = func;
  this.derivFunc = derivFunc;
}
mix(NewtonRaphsonMethod.prototype, {
  findOneRoot: function findOneRoot(guess, epsilon) {
    log('findOneRoot initial guess=' + guess);
    var f = this.func, df = this.derivFunc, i = -1;
    while (++i < this.maxIteration) {
      var value = f(guess);
      if (value === 0)
        return guess;

      var derivValue = df(guess);
      if (derivValue === 0)
        throw new Error('A derivative value must not be zero.');

      var nextGuess = guess - value / derivValue;
      log('guess=' + guess + ' next=' + nextGuess);
      if (numberEquals(nextGuess, guess, epsilon))
        return nextGuess;
      guess = nextGuess;
    }
    throw new Error('Iteration count reached the limit.');
  },
  maxIteration: 100,
  polishRoots: function polishRoots(roots, epsilon) {
    var polishedRoots = [];
    for (var i = 0, n = roots.length; i < n; ++i) {
      var root = roots[i];
      var value = this.func(root);
      if (value !== 0) {
        var newRoot = this.findOneRoot(root, epsilon);
        var newValue = this.func(newRoot);
        if (Math.abs(newValue) < Math.abs(value))
          polishedRoots.push(newRoot);
        else
          polishedRoots.push(root);
      }
      else
        polishedRoots.push(root);
    }
    return polishedRoots;
  }
});

function numberEquals(x, y, epsilon) {
  return epsilon ? Math.abs(x - y) <= epsilon : x === y;
}

function sum(values, addAlgorithm) {
  return (addAlgorithm || sum.defaultAddAlgorithm)(values);
}
sum.addWithKahanAlgorithm = function(values) {
  // OK: Firefox 3.5.3
  // NG: Safari 4.0.3 (6531.9), Chrome 4.0.221.8
  var n = values.length;
  if (n === 0)
    return 0;
  var s = values[0],
      c = 0,
      y, t;
  for (var j = 1; j < n; ++j) {
    y = values[j] - c;
    t = s + y;
    c = (t - s) - y;
    s = t;
  }
  return s;
}
sum.addWithSortAndAddAlgorithm = function(values) {
  // When we add floating point numbers, we should avoid undesireable effects
  // which cause loss of significant bits:
  // - Round-off error
  //   - Caused by additions of two values whose exponents are very different.
  //     The bits of the value with smaller exponent are lost.
  // - Cancellation
  //   - Caused by subtractions of nearly equal values or additions of 
  //     values with different signs and near absolutes.
  //
  // The algorithm implemented here:
  // - Calculate the sum of negative values and positive values separately
  //   in order to avoid cancellations.
  // - When adding the values of same sign, first sort them by the absolute
  //   values in ascending order.
  //
  // http://en.wikipedia.org/wiki/Floating_point
  // http://en.wikipedia.org/wiki/Numerical_stability
  // http://en.wikipedia.org/wiki/Loss_of_significance

  var n = values.length;
  switch (n) {
  case 0:
    return 0;
  case 1:
    return values[0];
  case 2:
    return values[0] + values[1];
  default:
    values = values.concat().sort(sum.compareFunc);

    var positiveStart = 0;
    while (positiveStart < n && values[positiveStart] < 0)
      ++positiveStart;

    var negativeSum = 0;
    for (var i = positiveStart - 1; i >= 0; --i)
      negativeSum += values[i];

    var positiveSum = 0;
    for (i = positiveStart; i < n; ++i)
      positiveSum += values[i];

    return positiveSum + negativeSum;
  }
};
sum.compareFunc = function(a, b) { return a - b; };
sum.addWithNaiveAlgorithm = function(values) {
  var total = 0;
  for (var i = 0, n = values.length; i < n; ++i)
    total += values[i];
  return total;
};
sum.defaultAddAlgorithm = sum.addWithSortAndAddAlgorithm;
//sum.defaultAddAlgorithm = sum.addWithKahanAlgorithm;

//function filter(values, func) {
//  var ret = [];
//  for (var i = 0, n = values.length; i < n; ++i) {
//    var value = values[i];
//    if (func(value))
//      ret.push(value);
//  }
//}

function integral(f, a, b, relativeEps, algorithm) {
  (algorithm || integral.defaultAlgorithm)(f, a, b, relativeEps);
}
integral.simpsonAlgorithm = function(f, a, b, relativeEps) {
  // http://en.wikipedia.org/wiki/Simpson%27s_rule
  var n = 2;
  var h = (b - a) / n;
  var endValues = [f(a), f(b)];
  var oddValuesTimes4 = [4 * f(a + h)];
  var oldResult = h / 3 * sum(endValues.concat(oddValuesTimes4));
  var evenValuesTimes2 = [];
  while (true) {
    n *= 2;
    h /= 2;
    for (var i = 0, m = oddValuesTimes4.length; i < m; ++i)
      evenValuesTimes2.push(oddValuesTimes4[i] / 2);
    oddValuesTimes4 = [];
    for (var i = 1; i < n; i += 2)
      oddValuesTimes4.push(4 * f(a + h * i));
    newResult = h / 3 * sum(
      endValues.concat(oddValuesTimes4, evenValuesTimes2));
    if (Math.abs((newResult - oldResult) / oldResult) < relativeEps)
      break;
    oldResult = newResult;
  }
//log('simpson n=' + n + ', Result=' + newResult);
  return newResult;
};
integral.defaultAlgorithm = integral.simpsonAlgorithm;

function sortBy(values, func) {
  var n = values.length;
  if (n <= 1)
    return values;

  var buf = [];
  for (var i = 0; i < n; ++i) {
    var value = values[i];
    buf[i] = {v: value, c: func(value)};
  }
  buf.sort(isArray(buf[0].c) ? sortBy.compareArray : sortBy.compareValue);

  var result = [];
  for (var i = 0; i < n; ++i)
    result[i] = buf[i].v;
  return result;
}
sortBy.compareValue = function(a, b) {
  var ac = a.c, bc = b.c;
  return ac < bc ? -1 : ac > bc ? 1 : 0;
}
sortBy.compareArray = function(a, b) {
  var n = a.length;
  for (var i = 0; i < n; ++i) {
    var d = sortBy.compareValue(a[i], b[i]);
    if (d)
      return d;
  }
  return 0;
}

var toString = Object.prototype.toString;
function isArray(obj) {
  return toString.call(obj) === '[object Array]';
}

function bind(func, obj) {
  return function() {
    return func.apply(obj, arguments);
  }
}

/*
 * Multi-dimensional (2, 3, ...) immutable vector.
 * - components: a vector whose dimension >= 2
 * - x, y      : a vector whose dimension = 2
 * - x, y, z   : a vector whose dimension = 3
 */
function Vector(x, y, z) {
  switch (arguments.length) {
  case 1:
    if (!isArray(x) || x.length < 2)
      throw new Error('An component array with 2 or greater dimension is expected.');
    this.components = x.concat();
    break;
  case 2:
    this.components = [x, y];
    break;
  case 3:
    this.components = [x, y, z];
    break;
  }
}
Vector.fromArray2d = function(array) {
  var n = array.length, vectors = [], i;
  for (i = 0; i < n; ++i)
    vectors.push(new Vector(array[i]));
  return vectors;
};
Vector.sum = function(vectors) {
  var vecCount = vectors.length, comps = [];
  for (var i = 0, dim = vectors[0].dimension(); i < dim; ++i) {
    var ithComps = [];
    for (var j = 0; j < vecCount; ++j)
      ithComps.push(vectors[j].component(i));
    comps.push(sum(ithComps));
  }
  return new Vector(comps);
};
Vector.polynomial = function(t, coefficientVectors) {
  var terms = [], n = coefficientVectors.length - 1,
      i = 0, ti = 1;
  while (true) {
    terms.push(coefficientVectors[i].clone().scalarMult(ti));
//console.log('i=' + i + ', coefficientVectors[i]=' + coefficientVectors[i]);
    if (++i > n)
      break;
    ti *= t;
  }
  return Vector.sum(terms);
};
Vector.zero = function(dimension) {
  var components = [], i;
  for (i = 0; i < dimension; ++i)
    components[i] = 0;
  return new Vector(components);
};
var proto = Vector.prototype;
proto.dimension = function() {
  return this.components.length;
}
proto.component = function(i) {
  return this.components[i];
}
proto.x = function() {
  return this.components[0];
};
proto.y = function() {
  return this.components[1];
};
proto.z = function() {
  return this.components[2];
};
proto.clone = function() {
  return new Vector(this.components.concat());
};
proto.equals = function(vec, epsilon) {
  var dim = this.dimension();
  if (vec.dimension() !== dim)
    return false;
  for (var i = 0; i < dim; ++i) {
    if (!numberEquals(this.component(i), vec.component(i), epsilon))
      return false;
  }
  return true;
};
proto.length = function() {
  return Math.sqrt(this.dotProduct(this));
};
proto.forEach = function(callback) {
  var n = this.dimension(),
      a = this.components,
      ret;
  for (var i = 0; i < n; ++i) {
    if (callback(a[i], i, a))
      return;
  }
};
proto.isSameDimension = function(vectorB) {
  return this.dimension() === vectorB.dimension();
}
proto.add = function(vectorB) {
  if (!this.isSameDimension(vectorB))
    throw new Error('Vector dimension unmatch.');
  var b = vectorB.components, c =[];
  this.forEach(function(ai, i, a) {
    c[i] = ai + b[i];
  });
  return new Vector(c);
};
proto.subtract = function(vectorB) {
  if (!this.isSameDimension(vectorB))
    throw new Error('Vector dimension unmatch.');
  var b = vectorB.components, c = [];
  this.forEach(function(ai, i, a) {
    c[i] = ai - b[i];
  });
  return new Vector(c);
};
proto.scalarMult = function(factor) {
  var b = [];
  this.forEach(function(ai, i, a) {
    b[i] = a[i] * factor;
  });
  return new Vector(b);
};
proto.scalarDiv = function(factor) {
  var b = [];
  this.forEach(function(ai, i, a) {
    b[i] = a[i] / factor;
  });
  return new Vector(b);
};
proto.dotProduct = function(vectorB) {
  var dim = this.dimension();
  if (vectorB.dimension() !== dim)
    throw new Error('Vector dimension unmatch.');
  var terms = [];
  for (var i = 0; i < dim; ++i)
    terms.push(this.component(i) * vectorB.component(i));
  return sum(terms);
};
proto.crossProduct = function(vectorB) {
  // Calculate cross product of this vector and vectorB.
  // Two dimensional vectors is treated as three dimensional vectors with
  // z-component = 0.
  //
  // http://en.wikipedia.org/wiki/Cross_product
  // http://en.wikipedia.org/wiki/Seven-dimensional_cross_product
  var dim = this.dimension(),
      bDim = vectorB.dimension();
  if ((dim === 2 || dim === 3) && (bDim === 2 || bDim === 3)) {
    var ax = this.x(), ay = this.y(), az = this.z() || 0,
        bx = vectorB.x(), by = vectorB.y(), bz = vectorB.z() || 0;
    this.components = [ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx];
  }
  else if (dim === 7 && bDim === 7) {
    var a1 = this.component(0), a2 = this.component(1),
        a3 = this.component(2), a4 = this.component(3),
        a5 = this.component(4), a6 = this.component(5),
        a7 = this.component(6),
        b1 = vectorB.component(0), b2 = vectorB.component(1),
        b3 = vectorB.component(2), b4 = vectorB.component(3),
        b5 = vectorB.component(4), b6 = vectorB.component(5),
        b7 = vectorB.component(6);
    this.components = [
      sum([x2 * y4, -x4 * y2, x3 * y7, -x7 * y3, x5 * y6, -x6 * y5]),
      sum([x3 * y5, -x5 * y3, x4 * y1, -x1 * y4, x6 * y7, -x7 * y6]),
      sum([x4 * y6, -x6 * y4, x5 * y2, -x2 * y5, x7 * y1, -x1 * y7]),
      sum([x5 * y7, -x7 * y5, x6 * y3, -x3 * y6, x1 * y2, -x2 * y1]),
      sum([x6 * y1, -x1 * y6, x7 * y4, -x4 * y7, x2 * y3, -x3 * y2]),
      sum([x7 * y2, -x2 * y7, x1 * y5, -x5 * y1, x3 * y4, -x4 * y3]),
      sum([x1 * y3, -x3 * y8, x2 * y6, -x6 * y2, x4 * y5, -x5 * y4])
    ];
  }
  else
    throw new Error('Cross product is supported for 2, 3 or 7 dimensional vectors only.');
  return this;
};
proto.swapComponents = function(i, j) {
  var c = this.components, b = c.concat();
      tmp = b[i];
  b[i] = b[j];
  b[j] = tmp;
  return new Vector(b);
};
proto.toString = function() {
  return '[' + this.components.join(' ') + ']';
};

// immutable three dimensional vector.
function Vector3d(x, y, z) {
  this.x = x;
  this.y = y;
  this.z = z;
}
Vector3d.sum = function(vectors) {
  var xs = [], ys = [], zs = [], n = vectors.length, i, v;
  for (i = 0; i < n; ++i) {
    v = vectors[i];
    xs.push(v.x);
    ys.push(v.y);
    zs.push(v.z);
  }
  return new Vector3d(sum(xs), sum(ys), sum(zs));
};
proto = Vector3d.prototype;
proto.length = function() {
  return Math.sqrt(this.dotProduct(this));
};
proto.dotProduct = function(vector) {
  return sum([this.x * vector.x, this.y * vector.y, this.z * vector.z]);
};
proto.crossProduct = function(vector) {
  return new Vector3d(
    this.y * vector.z - this.z * vector.y,
    this.z * vector.x - this.x * vector.z,
    this.x * vector.y - this.y * vector.x
  );
};
proto.negate = function() {
  return new Vector3d(-this.x, -this.y, -this.z);
};
proto.scalarMult = function(factor) {
  return new Vector3d(this.x * factor, this.y * factor, this.z * factor);
};
proto.scalarDiv = function(factor) {
  return new Vector3d(this.x / factor, this.y / factor, this.z / factor);
};


function Matrix() {
  var rowSize, columnSize, elements;
  if (isArray(arguments[0])) {
    elements = arguments[0];
    rowSize = elements.length;
    columnSize = elements[0].length;
  }
  else {
    rowSize = arguments[0];
    columnSize = arguments[1] || rowSize;
  }

  var m = this.rowSize = rowSize,
      n = this.columnSize = columnSize,
      a = this.elements = [];
  for (var i = 0; i < m; ++i) {
    a[i] = [];
    for (var j = 0; j < n; ++j)
      a[i][j] = elements ? elements[i][j] : 0;
  }
}
Matrix.identity = function(dimension) {
  var A = new Matrix(dimension, dimension),
      a = A.elements;
  for (var i = 0; i < dimension; ++i)
    a[i][i] = 1;
  return A;
};
Matrix.fromRowVector = function(vector) {
  return new Matrix([vector.components]);
};
Matrix.fromColumnVector = function(vector) {
  var n = vector.dimension(), v = vector.components,
      matrixA = new Matrix(n, 1), a = matrixA.elements, i;
  for (i = 0; i < n; ++i)
    a[i][0] = v[i];
  return matrixA;
};
Matrix.det2by2 = function(a, b, c, d) {
  return a * d - b * c;
};
Matrix.det3by3 = function(a, b, c, d, e, f, g, h, i) {
  return sum(
    [a * e * i, -a * f * h, b * f * g, -b * d * i, c * d * h, -c * e * g]);
};
var proto = Matrix.prototype;
proto.isSameSize = function(matrix) {
  return this.rowSize === matrix.rowSize &&
      this.columnSize === matrix.columnSize;
};
proto.isSquare = function() {
  return this.rowSize === this.columnSize;
};
proto.rowVector = function(i) {
  return new Vector(this.elements[i]);
};
proto.columnVector = function(j) {
  var m = this.rowSize, a = this.elements, components = [], i;
  for (i = 0; i < m; ++i)
    components[i] = a[i][j];
  return new Vector(components);
};
proto.clone = function() {
  return new Matrix(this.elements);
};
proto.forEach = function(callback) {
  var m = this.rowSize,
      n = this.columnSize,
      a = this.elements;
  for (var i = 0; i < m; ++i) {
    for (var j = 0; j < n; ++j) {
      if (callback(a[i][j], i, j, a))
        return;
    }
  }
};
proto.equals = function(matrixB) {
  if (!this.isSameSize(matrixB))
    return false;
  var b = matrixB.elements, ret = true;
  this.forEach(function(aij, i, j, a) {
    if (aij !== b[i][j]) {
      ret = false;
      return true;
    }
  });
  return ret;
};
proto.add = function(matrix) {
  if (!this.isSameSize(matrix))
    throw new Error('Cannot add a matrix of different size.');
  var A = this.clone(),
      b = matrix.elements;
  A.forEach(function(aij, i, j, a) {
    a[i][j] += b[i][j];
  });
  return A;
};
proto.subtract = function(matrix) {
  if (!this.isSameSize(matrix))
    throw new Error('Cannot subtract a matrix of different size.');
  var A = this.clone(),
      b = matrix.elements;
  A.forEach(function(aij, i, j, a) {
    a[i][j] -= b[i][j];
  });
  return A;
};
proto.scaleMult = function(factor) {
  if (!this.isSameSize(matrix))
    throw new Error('Cannot subtract a matrix of different size.');
  var A = this.clone(),
      b = matrix.elements;
  A.forEach(function(aij, i, j, a) {
    a[i][j] *= factor;
  });
  return A;
};
proto.scaleDiv = function(factor) {
  if (!this.isSameSize(matrix))
    throw new Error('Cannot subtract a matrix of different size.');
  var A = this.clone(),
      b = matrix.elements;
  A.forEach(function(aij, i, j, a) {
    a[i][j] /= factor;
  });
  return A;
};
proto.transpose = function() {
  var m = this.rowSize,
      n = this.columnSize,
      B = new Matrix(n, m),
      b = B.elements;
  this.forEach(function(aij, i, j, a) {
    b[j][i] = aij;
  });
  return B;
};
proto.multiply = function() {
  var arg = arguments[0],
      m = this.rowSize,
      n = this.columnSize,
      a = this.elements,
      terms;
  if (arg instanceof Matrix) {
    var matrix = arg;
    if (this.columnSize !== matrix.rowSize)
      throw new Error('Cannot multiply a matrix whose row size does not equal to this matrix\'s column size.');
    var p = matrix.columnSize,
        b = matrix.elements,
        C = new Matrix(m, p),
        c = C.elements;
    for (var i = 0; i < m; ++i) {
      for (var j = 0; j < p; ++j) {
        terms = [];
        for (var r = 0; r < n; ++r)
          terms.push(a[i][r] * b[r][j]);
        c[i][j] = sum(terms);
      }
    }
    return C;
  }
  else if (arg instanceof Vector) {
    return this.multiply(Matrix.fromColumnVector(arg)).columnVector(0);
//    var vector = arg, b = vector.components;
//    if (this.columnSize !== vector.dimension())
//      throw new Error('Cannot multiply a vector whose dimension does not equal to this matrix\'s column size.');
//    var c = [];
//    for (var i = 0; i < m; ++i) {
//      terms = [];
//      for (var j = 0; j < n; ++j)
//        terms.push(a[i][j] * b[j]);
//      c[i] = sum(terms);
//    }
//    return new Vector(c);
  }
};
proto.toString = function() {
  var m = this.rowSize,
      n = this.columnSize,
      a = this.elements,
      rows = [],
      columns;
  for (var i = 0; i < m; ++i) {
    columns = [];
    for (var j = 0; j < n; ++j) {
      columns.push(String(a[i][j]));
    }
    rows.push('[' + columns.join(' ') + ']');
  }
  return rows.join('\n');
};
proto.swapRows = function(i1, i2) {
  var a = this.elements,
      tmp = a[i2];
  a[i2] = a[i1];
  a[i1] = tmp;
};
proto.swapColumns = function(j1, j2) {
  var a = this.elements,
      n = this.rowCount,
      i, tmp;
  for (i = 0; i < n; ++i) {
    tmp = a[i][j1];
    a[i][j1] = a[i][j2];
    a[i][j2] = tmp;
  }
};

function GaussElimination(matrixA, vectorB) {
  this.matrixA = matrixA.clone();
  this.vectorB = vectorB.clone();
  if (vectorB.dimension() !== matrixA.columnSize)
    throw new Error('Column size of the matrix should equal to the vector dimension.');

  var m = this.matrixA.rowSize, n = this.matrixA.columnSize,
      a = this.matrixA.elements, b = this.vectorB.components, i, j, k, l, iMax;

  i = 0;
  j = 0;
  this.singular = false;
  while (i < m && j < n) {
    // Find pivot in column j, starting in row i.
    iMax = i;
    for (k = i + 1; k < m; ++k) {
      if (Math.abs(a[k][j]) > Math.abs(a[iMax][j]))
        iMax = k;
    }

    if (a[iMax][j] !== 0) {
      matrixA.swapRows(i, iMax);
      this.vectorB = this.vectorB.swapComponents(i, iMax);
      b = this.vectorB.components;
      for (l = j + 1; l < n; ++l)
        a[i][l] /= a[i][j];
      b[i] /= a[i][j];
      a[i][j] = 1;
      for (k = i + 1; k < m; ++k) {
        for (l = j + 1; l < n; ++l)
          a[k][l] -= a[k][j] * a[i][l];
        b[k] -= a[k][j] * b[i];
        a[k][j] = 0;
      }
      ++i;
    }
    else
      this.singular = true;
    ++j;
  }
};
var proto = GaussElimination.prototype;
proto.isSingular = function() {
  return this.singular;
};
proto.solution = function() {
  var matrixA = this.matrixA, vectorB = this.vectorB,
      m = matrixA.rowSize, n = matrixA.columnSize, a = matrixA.elements,
      b = vectorB.components, i, j;
  var x = [];
  i = m - 1;
  x[i] = b[i];
  for (--i; i >= 0; --i) {
    var terms = [b[i]];
    for (j = i + 1; j < n; ++j)
      terms.push(-a[i][j] * x[j]);
    x[i] = sum(terms) / a[i][i];
  }
  return new Vector(x);
};

function LUDecomposition(matrixA, singularityThreshold) {
  if (singularityThreshold === undefined)
    singularityThreshold = LUDecomposition.defaultSingularityThreashold;
  if (!matrixA.isSquare())
    throw new Error('Square matrix needed.');

  this.matrixLU = matrixA.clone();
  var n = this.matrixLU.columnSize, lu = this.matrixLU.elements, pivot,
      row, col, i, rowMax, valMax, luRow, terms, absSum, luDiag;
  pivot = [];
  for (row = 0; row < n; ++row)
    pivot[row] = row;
  this.pivotVector = new Vector(pivot);
  this.parity = 1;
  this.singular = false;

  for (col = 0; col < n; ++col) {
    // upper
    for (row = 0; row < col; ++row) {
      luRow = lu[row];
      terms = [luRow[col]];
      for (i = 0; i < row; ++i)
        terms.push(-luRow[i] * lu[i][col]);
      luRow[col] = sum(terms);
    }

    // lower
    rowMax = col;
    valMax = -Infinity;
    for (row = col; row < n; ++row) {
      luRow = lu[row];
      terms = [luRow[col]];
      for (i = 0; i < col; ++i)
        terms.push(-luRow[i] * lu[i][col]);
      luRow[col] = sum(terms);

      absSum = Math.abs(luRow[col]);
      if (absSum > valMax) {
        valMax = absSum;
        rowMax = row;
      }
    }
//log('col=' + col + ', rowMax=' + rowMax);

    // Singularity check
    if (Math.abs(lu[rowMax][col]) <= singularityThreshold) {
      this.singular = true;
      return;
    }
    // Pivot if necessary
    if (rowMax !== col) {
      this.matrixLU.swapRows(rowMax, col);
      this.pivotVector = this.pivotVector.swapComponents(rowMax, col);
      this.parity = -this.parity;
    }

    // Divide the lower elements by the "winning" diagonal element.
    luDiag = lu[col][col];
    for (row = col + 1; row < n; ++row)
      lu[row][col] /= luDiag;
  }
}
LUDecomposition.defaultSingularityThreashold = 1e-12;
var proto = LUDecomposition.prototype;
proto.isSingular = function() {
  return this.singular;
};
proto.matrixL = function() {
  if (!this._matrixL && !this.singular) {
    var n = this.pivotVector.dimension(), lu = this.matrixLU.elements;
    this._matrixL = new Matrix(n, n);
    var l = this._matrixL.elements, i, j, luI;
    for (i = 0; i < n; ++i) {
      luI = lu[i];
      for (j = 0; j < i; ++j)
        l[i][j] = luI[j];
      l[i][i] = 1;
    }
  }
  return this._matrixL;
};
proto.matrixU = function() {
  if (!this._matrixU && !this.singular) {
    var n = this.pivotVector.dimension(), lu = this.matrixLU.elements;
    this._matrixU = new Matrix(n, n);
    var u = this._matrixU.elements, i, j, luI;
    for (i = 0; i < n; ++i) {
      luI = lu[i];
      for (j = i; j < n; ++j)
        u[i][j] = luI[j];
    }
  }
  return this._matrixU;
};
proto.matrixP = function() {
  if (!this._matrixP && !this.singular) {
    var n = this.pivotVector.dimension(), pivot = this.pivotVector.components;
    this._matrixP = new Matrix(n, n);
    var p = this._matrixP.elements, i;
    for (i = 0; i < n; ++i) {
      p[i][pivot[i]] = 1;
    }
  }
  return this._matrixP;
};
proto.determinant = function() {
  var n = this.matrixLU.rowSize, lu = this.matrixLU.elements,
      det = this.parity, i;
  for (i = 0; i < n; ++i)
    det *= lu[i][i];
  return det;
};
proto.solve = function(vectorB) {
  var n = this.pivotVector.dimension(), pivot = this.pivotVector.components,
      lu = this.matrixLU.elements, bp = [], row, col, i, j, bpCol;
  if (!this.pivotVector.isSameDimension(vectorB))
    throw new Error('Vector dimension mismatch.');
  if (this.isSingular())
    throw new Error('Singular matrix.');

  // Apply permutations to b
  for (row = 0; row < n; ++row)
    bp[row] = vectorB.component(pivot[row]);

  // Solve LY = b
  for (col = 0; col < n; ++col) {
    bpCol = bp[col];
    for (i = col + 1; i < n; ++i)
      bp[i] -= bpCol * lu[i][col];
  }

  // Solve UX = Y
  for (col = n - 1; col >= 0; --col) {
    bp[col] /= lu[col][col];
    bpCol = bp[col];
    for (i = 0; i < col; ++i)
      bp[i] -= bpCol * lu[i][col];
  }

  return new Vector(bp);
};

function uniq(array) {
  var hash = {}, n = array.length, ret = [], i;
  for (i = 0; i < n; ++i) {
    var elem = array[i];
    if (!(elem in hash)) {
      hash[elem] = elem;
      ret.push(elem);
    }
  }
  return ret;
}

function uniqAndSort(array) {
  array = uniq(array);
  array.sort(function(a, b) {
    return a < b ? -1 : a > b ? 1 : 0;
  });
  return array;
}

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

function sgn(x) {
  if (isNaN(x))
    return NaN;
  else if (x > 0)
    return 1;
  else if (x < 0)
    return -1;
  else
    return 0;
}

function log() {
  if (console && console.log)
    console.log.apply(console, arguments);
}

function Set(values) {
  this.values = {};
  for (var i = 0, n = values.length; i < n; ++i) {
    var value = values[i];
    if (!(value in this.values))
      this.values[value] = value;
  }
}
mix(Set.prototype, {
  intersect: function intersect(set) {
    var values = this.values;
    if (!(value in values))
      values[value] = value;
    return this;
  }
});

function mix(dest, src) {
  for (var k in src)
    dest[k] = src[k];
  return dest;
}

return {
  MACHINE_EPSILON: MACHINE_EPSILON,
  log: log,
  bind: bind,
  numberEquals: numberEquals,
  sum: sum,
  sgn: sgn,
  uniq: uniq,
  uniqAndSort: uniqAndSort,
  filterValues: filterValues,
  integral: integral,
  sortBy: sortBy,
  OneRootFinder: OneRootFinder,
  QuadraticEquation: QuadraticEquation,
  CubicEquation: CubicEquation,
  Polynomial: Polynomial,
  Vector: Vector,
  Vector3d: Vector3d,
  Matrix: Matrix,
  GaussElimination: GaussElimination,
  LUDecomposition: LUDecomposition
};

}();
