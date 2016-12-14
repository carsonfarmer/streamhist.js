# StreamHist.js

Javascript implementation of a streaming approximate histogram.

## Description

The streaming histogram data structure is a modified version of the
structure presented in Ben-Haim & Tom-Tov's [Streaming Parallel Decision Tree
Algorithm](http://jmlr.org/papers/volume11/ben-haim10a/ben-haim10a.pdf) (and
based on work by
[Guedalia *et al*](http://www.cs.huji.ac.il/~werman/Papers/guedalia_etal99.pdf)
and [Zhang *et al.*](http://dl.acm.org/citation.cfm?id=1044332)), which reduces
the computational complexity from $O(K_{max}N)$ to $O(log{K_{max}}N)$ without
increasing the space complexity, which is $O(K_{max})$. Additionally, the
histogram can track _multiple_ statistics (in addition to the mean) to help
guide the classification process. This _does_ increase the memory requirements
slightly, but improves the ability to find 'natural breaks' in the data.

The basic algorithm is essentially 'parameter-free', in that it doesn't
require the user to specify any information about the underlying data
distribution explicitly. There are, however, a few parameters that help
control the accuracy of the sketch, including the maximum number of bins,
and whether to use weighting when computing bin means. A larger maximum
number of bins yields a more accurate approximation of the underlying
distribution, at the cost of increased memory utilization and performance.
There is no 'optimal' bin count, but less than 100 bins should generally be
sufficient. Bin weighting is useful for downplaying outliers in the dataset,
as well as providing a more accurate fit to the data in higher density areas
(i.e., more bins in higher density areas).

This Javascript implementation is based on a reading of the paper, with some
inspiration drawn from a
[Javascript implementation](https://github.com/welch/tdigest) of Dunning's
T-Digest for streaming quantile approximation by Will Welch.

## Quick Start

**Warning**: This is total alpha version stuff here... very fragile!

### Node

```bash
npm install streamhist
```

```javascript
var StreamHist = require('streamhist').StreamHist;
var x=[], N = 100000;
for (var i = 0 ; i < N ; i += 1) {
    x.push(Math.random() * 10 - 5);
};
hist = new StreamHist();
hist.push(x);
console.log(td.summary());
console.log("median ~ "+hist.percentile(0.5));
```

### In the browser

The `grunt dist` task has been configured to generate a self-contained
UMD-wrapped version of `streamhist` in `dist/streamhist.js`.

Embed it in HTML like this:

```html
<script src="dist/streamhist.js"></script>
<script>
    var hist = new this.streamhist.StreamHist();
    for (var i=0; i < 1000000; i++) {
        hist.push(Math.random());
    }
    document.write(hist.summary())
</script>
```

## Dependencies

* `bintrees`: https://www.npmjs.com/package/bintrees
* `underscore`: https://www.npmjs.com/package/underscore
* `argminmax`: https://www.npmjs.com/package/argminmax
