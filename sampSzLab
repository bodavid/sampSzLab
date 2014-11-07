var sampSzLab = function() {
  var cryptoObj = window.crypto || window.msCrypto;
  this.boxMuller = function(mean, stdDev, n) {
    if (!n) {
      n = 1;
    } 
    var arr = new Uint32Array(4*n);
    cryptoObj.getRandomValues(arr);
    var i = 0, j = 0, rnd;
    var normDist = Float64Array(n);
    while (i < n - 1)  {
      normDist[j++] = mean + stdDev * (Math.sqrt(-2 * Math.log((arr[i++] + 4294967296 * arr[i++])/18446744073709551616)) * Math.cos(2 * Math.PI * (arr[i++] + 4294967296 * arr[i++])/18446744073709551616));
    }
    return normDist;
  }
}

