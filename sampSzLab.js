/**
 * Visualizing sample sizes dependency on some vector spaces (Fourier transforms) with noise.
 *
 * This code is not yet finished and does not run.
 *
 * Basis function periodicity can be 2*pi or rational. Rational might seem easier but has the disadvantage of being
 * specific on borders on a computer architecture. Basis functions can be unspecified in and around discontinuities.
 * The value should be chosen according to the continuous transform and not according to definition.
 */

var SampSzLab = function(defaults) {

//  "use strict";
  var szThis = this;

  var cryptoObj = window.crypto || window.msCrypto;

  /** Default parameters */
  var defaults = defaults || {};
  this.sampleSize = defaults.sampleSize || 10007;
  this.dimensions = defaults.dimensions || 23;
  this.func = defaults.func || null;
  this.basisFunction = defaults.basisFunction || Math.cos;
  this.interval = defaults.interval || 2 * Math.PI;
  this.normedBasis = defaults.normedBasis || true;
  this.basisMatrix = defaults.basisMatrix || null;
  this.period = defaults.period || 2 * Math.PI;
  this.multiplier = defaults.multiplier || 1;
  this.epsilon = defaults.epsilon || 5.56e-17;

  //Somewhat stronger or more uniform than Math.random()
  this.random = function() {
    var arr = new Uint32Array(2);
    cryptoObj.getRandomValues(arr);
    return arr[0]/4294967296 + arr[1]/18446744073709551616;
  }

  this.findMachineEpsilon = function() {
    var eps = 1e-9;
    while ((1 - eps) != 1) {
      eps *= 0.9999
    };
    return eps
  }

    // Based on  http://www.protonfish.com/jslib/boxmuller.shtml
    // modified for more uniform random() with window.crypto.
    this.boxMullerVec = function(n) {
      var i, x, y, rds, c;
      var n = n || 2;
      var bmArr = Array(n+1);

      for (i = 0; i < n; i+=2) {
        // Get two random numbers from -1 to 1.
        // If the radius is zero or greater than 1, throw them out and pick two new ones
        // Rejection sampling throws away about 20% of the pairs.
        do {
          x = this.random() * 2 - 1;
          y = this.random() * 2 - 1;
          rds = x*x + y*y;
        }
        while (rds == 0 || rds > 1)

        // The Box-Muller Transform
        c = Math.sqrt(-2 * Math.log(rds)/rds);
        bmArr[i] = x * c;
        bmArr[i+1] = y * c;
      }
      bmArr.length = n;
      return bmArr;
  }

  // Make a normed basis vector, based on func
  this.frequencyBasis = function() {
    var freqVec = new Array(this.sampleSize);
    var dx = this.interval / this.sampleSize;
    var l2sum = 0;
    for (var i=0; i < this.sampleSize; i++) {
      freqVec[i] = this.basisFunction( this.multiplier * dx * i);
      l2sum += freqVec[i] * freqVec[i];
    }
    if (this.normedBasis) {
        var l2lengthInv = 1 / Math.sqrt(l2sum);
        for (var i=0; i < this.sampleSize; i++) {
          freqVec[i] = l2lengthInv * freqVec[i];
        }
    }
//    this.freqVec = freqVec;
    return freqVec;
  }

  // Make a basis matrix
  this.makeBasisSpectrumMatrix = function(dimensions){
    var basisSpectrumMatrix = new Array(this.sampleSize);

    for (var i=0; i < this.sampleSize; i++) {
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
   *  Where the steep functions changes values are machine dependent.
   *
   *  The SteepTrans library is made as an experiment by David Jonsson 2014, 2015.
   */

  /** The sinSteep(x) function is very similar to Math.round(Math.sin(x)), periodic and looks like
   *
   *  __
   * _  __  _
   *      __
   *
   * over one periodic (2*pi). It is 0, 1 or -1. Using a transcendental number as periodic lowers the risks of aliasing and
   * keeps some features of the sine function.
   */

    this.sinSteep = function(t) {
      var tmod = (t / this.period) % 1;
      if (tmod < 1/8)
       return 0
      else if (tmod < 3/8)
       return Math.sign(tmod)
      else if (tmod < 5/8)
       return 0
      else if (tmod < 7/8)
       return -Math.sign(tmod)
      else return 0;
     }

  /** The cosSteep(x) ( ≈ Math.round(Math.cos(x)) ) function is periodic and looks like
   *
   * _      _
   *  __  __
   *    __
   *
   * over one periodic (2*pi). It is 1, 0 or -1. Using a transcendental number as periodic lowers the risks of aliasing and
   * keeps some features of the cosine function.
   */

    this.cosSteep = function(t) {
      var tmod = Math.abs((t / this.period) % 1);
      if (tmod < 1/8)
        return 1
      else if (tmod >= 3/8)
        return -1
      return 0;
     }

    /** The dotProduct between two functions can be used as a way to evaluate the orthogonality of two
    *  functions. Trapeziod ( mid interval 0.5 ).
    */
    this.dotProduct = function(comparefunc) {

      var product = Array[szThis.dimensions];
      var prod = 0;
      var dt = szThis.period / szThis.sampleSize;

      for (var dim = 1; dim <= szThis.dimensions; dim++) {
        for (var samp = 0.5; samp < szThis.sampleSize; samp++) {
         prod += this.func(samp * dt) * this.func(dim * samp * dt);
        }
        product[dim - 1] = prod;
      }

      return product;
    }

    /** The magnitude of a function is the length, size or norm of the function in the sense of functional analysis or
     *  linear algebra.
     */
    this.magnitude = function() {
      return Math.sqrt((this.dotProduct(null, null, 1))[0]);
    }

    this.plot2D = function(arr3) {
      var c = document.getElementById(szThis.plotDOMid);
      var ctx = c.getContext("2d");
      var imgData = ctx.getImageData(0, 0, c.width, c.height);
      var j = 0; // Image pixel color index
      for (var i = 0; i < arr3[szThis.RGBbuffer].length; i += 3) {
        imgData.data[j++] = (szThis.extremesRGB.rMax === szThis.extremesRGB.rMin) ? 0 : 255 * (arr3[szThis.GXYbuffer][i]     - szThis.extremesRGB.rMin) / (szThis.extremesRGB.rMax - szThis.extremesRGB.rMin);
        imgData.data[j++] = (szThis.extremesRGB.gMax === szThis.extremesRGB.gMin) ? 0 : 255 * (arr3[szThis.GXYbuffer][i + 1] - szThis.extremesRGB.gMin) / (szThis.extremesRGB.gMax - szThis.extremesRGB.gMin);
        imgData.data[j++] = (szThis.extremesRGB.bMax === szThis.extremesRGB.bMin) ? 0 : 255 * (arr3[szThis.GXYbuffer][i + 2] - szThis.extremesRGB.bMin) / (szThis.extremesRGB.bMax - szThis.extremesRGB.bMin);
        imgData.data[j++] = 255;
      }
      ctx.putImageData(imgData, 0, 0);
      
      this.calcStepIterationDOM.innerHTML = this.iteration;
      this.calcStepDurationDOM.innerHTML = this.calcStepDuration;
    };

    this.plotN2D = function(command) {
      var sideLength = 1024, i, j, pixs = 3 * sideLength * sideLength;
      var rndPair = [];
      var plot2Darr = Array(pixs);
      szThis.extremesRGB = {};
      szThis.extremesRGB.rMax = 255;
      szThis.extremesRGB.gMax = 255;
      szThis.extremesRGB.bMax = 255;

      if (defaults && defaults.start) {
        this.draw = function() {
          while  (i < pixs) {
            rndPDair = szThis.boxMullerVec();
            plot2Darr[i++] = rndPair[0];
            plot2Darr[i++] = rndPair[1];
          }
          szThis.plot2D(plot2Darr);
        };
        szThis.nextPlotEvent = requestAnimationFrame(function() {
          //this.
        });
      };
      if (defaults && defaults.stop) {
        cancelAnimationFrame(szThis.nextPlotEvent);
      }



    }

    /** The test object can be used as unit testing and as a way to evaluate the
     *  functions on various platforms, settings and biases.
     */
    this.tests = {
      orthogonality: function() {
        return szThis.dotProduct()
      }

    }

}

