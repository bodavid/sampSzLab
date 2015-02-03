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
  this.multiplier = 1;

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
      freqVec[i] = this.basisFunc( this.multiplier * dx * i);
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

/** SteepTrans is a library for analysing rapidly changing processes in noisy environments.
 *  It is designed to be able to optimize due to many zeros and many shared values in common in the basis functions.
 *
 *  Readability of the code has had higher priority than optimization. The code would need optimization for larger sample sizes.
 *
 *  You better know what you are doing when using this lib otherwise you will get inconclusive results and unknown precision.
 *
 *  The SteepTrans library is made as an experiment by David Jonsson 2014.
 */

 this.SteepTrans = function(defaults) {

   /** Default parameters */
   var defaults = defaults || {},
   sampleSize = defaults.sampleSize || 10007,
   dimensions = defaults.dimensions || 23,
   periodicity = defaults.periodicity || 2 * Math.PI;
   func = defaults.func || null;
   basisFunction = defaults.basisFunction || null;

   var stThis = this;

/** The steep(x) function is periodic and looks like
 *
 *  __
 * _  __  _
 *      __
 *
 * over one periodic (2*pi). It is 0, 1 or -1. Using a transcendental number as periodic lowers the risks of aliasing and
 * keeps some features of the sine function.
 */

  prototype.steep = function(t) {
    var tmod = (t / periodicity) % 1;
     if (tmod < 1/8)
       return 0
     else if (tmod < 3/8)
       return 1
     else if (tmod < 5/8)
       return 0
     else if (tmod < 7/8)
       return -1
     else return 0;
   }

  /** The dotProduct between two functions can be used as a way to evaluate the orthogonality of two
  *  functions.
  */
  prototype.dotProduct = function(func, basisFunction, dimensions) {
    if (!(func) {
      func = stThis.steep;
    }
    if (!(basisFunction) {
      basisFunction = stThis.steep;
    }
    if (!(dimensions) {
      dimensions = defaultDimensions;
    }

    var product = Array[dimensions];
    var prod = 0;
    var dt = periodicity / sampleSize;

    for (var dim = 1; dim =< dimensions; dim++) {
      for (var samp = 0.5; samp < sampleSize; samp++) {
       prod += func(samp * dt) * basisFunction(dim * samp * dt);
      }
      product[dim - 1] = prod;
    }

    return product;
  }

  /** The magnitude of a function is the length, size or norm of the function in the sense of functional analysis or
   *  linear algebra.
   */
 prototype.magnitude = function() {
  return Math.sqrt((this.dotProduct(null, null, 1))[0]);
 }

  /** The test object can be used as unit testing and as a way to evaluate the
   *  functions on various platforms, settings and biases.
   */
  prototype.tests = {
   orthogonality: function() {
    return stThis.dotProduct()
  }

 }

}

}

