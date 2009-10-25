function Stopwatch() {}
Stopwatch.prototype.start = function start() {
  this.startTime = +new Date;
};
Stopwatch.prototype.stop = function stop() {
  this.stopTime = +new Date;
};
Stopwatch.prototype.elapsedTimeInSeconds = function elapsedTimeInSeconds() {
  return (this.stopTime - this.startTime) / 1000;
};
Stopwatch.prototype.measure = function measure(count, func) {
  this.start();
  try {
    for (var i = 0; i < count; ++i)
      func();
  }
  finally {
    this.stop();
  }
  return this.elapsedTimeInSeconds();
};
