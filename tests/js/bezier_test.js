function testBezier() {

  module('bezier');

  test('bezier polynominal', function() {
    var cp1 = [[300, 200], [1000, 500], [0, 500], [600, 200]],
        bezier1 = new geom2d.Bezier(numeric.Vector.fromArray2d(cp1));
    ok(bezier1.pointAtT(0).equals(new numeric.Vector(cp1[0])), 'start point');
    ok(bezier1.pointAtT(1).equals(new numeric.Vector(cp1[3])), 'end point');
//console.log('bezier1.pointAtT2(0.5)=' + bezier1.pointAtT2(0.5));
    ok(bezier1.pointAtT2(0.5).equals(bezier1.pointAtT(0.5)), 'point at t = 0.5');
//console.log('t=' + 0.3 + ', pointAtT2=' + bezier1.pointAtT2(0.3) + ', pointAtT=' + bezier1.pointAtT(0.3));
    ok(bezier1.pointAtT2(0.3).equals(bezier1.pointAtT(0.3), 1e-12), 'point at t = 0.3');
//console.log('t=' + 0.8 + ', pointAtT2=' + bezier1.pointAtT2(0.8) + ', pointAtT=' + bezier1.pointAtT(0.8));
    ok(bezier1.pointAtT2(0.8).equals(bezier1.pointAtT(0.8), 1e-12), 'point at t = 0.8');
  });

  test('two bezier intersections', function() {
    var cp1 = [[300, 200], [1000, 500], [0, 500], [600, 200]],
        cp2 = [[400, 500], [450, -300], [500, 1100], [550, 100]],
        bezier1 = new geom2d.Bezier(numeric.Vector.fromArray2d(cp1)),
        bezier2 = new geom2d.Bezier(numeric.Vector.fromArray2d(cp2));
    geom2d.Bezier.intersections(bezier1, bezier2);
  });
}
