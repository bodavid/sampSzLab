/**
 * Visualizing sample sizes dependency on even Fourier transforms with noise.
 *
 * Basis function periodicity can be 2*pi or rational. Rational might seem easier but has the disadvantage of being
 * specific on borders on a computer architecture. Basis functions can be unspecified in and around discontinuities.
 */

var SampSzLab = function(defaults) {

//  "use strict";

  var cryptoObj = window.crypto || window.msCrypto;

  this.random = function() {
    var arr = new Uint32Array(2);
    cryptoObj.getRandomValues(arr);
    return arr[0]/4294967296 + arr[1]/18446744073709551616;
  }

  this.basisFunc = Math.cos;
  this.normedBasis = true;
  this.basisMatrix = null;
  this.sampSz = 1007;
  this.period = 2 * Math.PI;
  this.interval = 2 * Math.PI;

  this.rndBMt = function() {
    // Based on  http://www.protonfish.com/jslib/boxmuller.shtml
    // modified for more uniform random()
      var x = 0, y = 0, rds, c;

      // Get two random numbers from -1 to 1.
      // If the radius is zero or greater than 1, throw them out and pick two new ones
      // Rejection sampling throws away about 20% of the pairs.
      do {
        x = this.random()*2-1;
        y = this.random()*2-1;
        rds = x*x + y*y;
      }
      while (rds == 0 || rds > 1)

      // This magic is the Box-Muller Transform
      c = Math.sqrt(-2*Math.log(rds)/rds);

      // It always creates a pair of numbers. I'll return them in an array.
      // This function is quite efficient so don't be afraid to throw one away if you don't need both.
      return [x*c, y*c];
  }

  this.boxMuller = function(retArr, func, stdDev) {
    if (typeof n === 'undefined') n = retArr.length;
    if (typeof func === 'undefined') func = function(x){return 0};
    if (typeof stdDev === 'undefined') stdDev = 1;

    var arr = new Uint32Array(3 * retArr.length + 3);
    cryptoObj.getRandomValues(arr);
    var i = 0, j = 0, rndR;
    var normDist = new Array(n+1);
    while (j < n )  {
//      normDist[j++] = func(j-1) + stdDev * (Math.sqrt(-2 * Math.log((arr[i++]/4294967296 + arr[i++])/18446744073709551616)) * Math.cos(2 * Math.PI * (arr[i++]/4294967296 + arr[i++])/18446744073709551616));
      rndR = -2 * Math.log(arr[i++]/4294967296 + arr[i++]/18446744073709551616);
      normDist[j++] =  func(j) + stdDev * (Math.sqrt(rndR * Math.sin(2 * Math.PI * (arr[i++]/4294967296 + arr[i++]/18446744073709551616))));
      normDist[j++] =  func(j) + stdDev * (Math.sqrt(rndR * Math.cos(2 * Math.PI * (arr[i++]/4294967296 + arr[i++]/18446744073709551616))));
    }
    normDist.length = n;
    return normDist;
  }

  // make a basis vector, based on func
  this.frequencyBasis = function() {
//    var freqVec = new Float64Array(this.sampSz);
    var freqVec = new Array(this.sampSz);
    var dx = this.interval / this.sampSz;
    var l2sum = 0;
    for (var i=0; i < this.sampSz; i++) {
      freqVec[i] = this.basisFunc( dx * i);
      l2sum += freqVec[i] * freqVec[i];
    }
    if (this.normedBasis) {
        var l2lengthInv = 1 / Math.sqrt(l2sum);
        for (var i=0; i < this.sampSz; i++) {
          freqVec[i] = l2lengthInv * freqVec[i];
        }
    }
//    this.freqVec = freqVec;
    return freqVec;
  }

  //make a basis matrix
  this.makeBasisSpectrumMatrix = function(dimensions){
    var basisSpectrumMatrix = new Array(this.sampSz);

    for (var i=0; i < this.sampSz; i++) {
      this.interval = this.period * i;
      basisSpectrumMatrix[i] = this.frequencyBasis();
    }
    this.BasisSpectrumMatrix = BasisSpectrumMatrix;
    return BasisSpectrumMatrix;
  }

  this.makeSpectrum = function() {

  }

}

