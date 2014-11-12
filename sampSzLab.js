/**
 * Visualizing sample sizes dependency on even Fourier transforms with noise
 */

var SampSzLab = function() {

  "use strict";

  var cryptoObj = window.crypto || window.msCrypto;

  this.boxMuller = function(func, stdDev, n) {
    if (!n) {
      n = 1;
    } 
    var arr = new Uint32Array(4*n);
    cryptoObj.getRandomValues(arr);
    var i = 0, j = 0, rnd;
    var normDist = Float64Array(n);
    while (i < n - 1)  {
      normDist[j++] = mean + stdDev * (Math.sqrt(-2 * Math.log((arr[i++]/4294967296 + arr[i++])/18446744073709551616)) * Math.cos(2 * Math.PI * (arr[i++]/4294967296 + arr[i++])/18446744073709551616));
    }
    return normDist;
  }

  // make a basis vector, based on func
  // func must have period 2pi
  this.frequencyBasis = function(func, normed) {
       if(typeof(normed)==='undefined') normed = true;

  }

  //make a basis matrix
  this.makeBasisSpectrum = function(dimensions){

  }

  this.makeSpectrum = function() {

  }

}

