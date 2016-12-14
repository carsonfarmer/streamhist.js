/**
 * @fileoverview StreamHist:
 * Javascript implementation of a streaming approximate histogram.
 *
 * The streaming histogram data structure is a modified version of the
 * structure presented in Ben-Haim & Tom-Tov's [Streaming Parallel Decision Tree
 * Algorithm]{@link http://jmlr.org/papers/volume11/ben-haim10a/ben-haim10a.pdf}
 * (and based on work by [Guedalia et al.]{@link http://www.cs.huji.ac.il/
 * ~werman/Papers/guedalia_etal99.pdf} and [Zhang *et al.*]{@link
 * http://dl.acm.org/citation.cfm?id=1044332}), which reduces the computational
 * complexity from $O(K_{max}N)$ to $O(log{K_{max}}N)$ without increasing the
 * space complexity, which is $O(K_{max})$. This is acheived by using an
 * internal Red-Black Search Tree for storing histogram bins.
 *
 * The basic algorithm is essentially 'parameter-free', in that it doesn't
 * require the user to specify any information about the underlying data
 * distribution explicitly. There are, however, a few parameters that help
 * control the accuracy of the sketch, including the maximum number of bins,
 * and whether to use gap weighting when merging bins. A larger maximum
 * number of bins yields a more accurate approximation of the underlying
 * distribution, at the cost of increased memory utilization and performance.
 * There is no 'optimal' bin count, but less than 100 bins should generally be
 * sufficient. Gap weighting is useful for encouraging the histogram to spend
 * more of its bins capturing the densest areas of the distribution. For a
 * Gaussian distribution this means better resolution near the mean and less
 * resolution near the tails.
 *
 * This Javascript implementation is based on a reading of the paper, with some
 * inspiration drawn from a [Javascript implementation]{@link
 * https://github.com/welch/tdigest} of Dunning's T-Digest for streaming
 * quantile approximation by Will Welch.
 *
 * @author carsonfarmer@gmail.com (Carson Farmer)
 * @license Released under the MIT license.
 * Copyright (c) 2016 Carson Farmer.
 * Portions inspired by tdigest (https://github.com/welch/tdigest)
 * Copyright (c) 2015 Will Welch
 */

// Additional information:
// http://druid.io/blog/2013/09/12/the-art-of-approximating-distributions.html

var RBTree = require('bintrees').RBTree,
    math = require('mathjs');
// Adds a helper function (which isn't efficient) to RBTree
// Useful really only for infrequent searches by index...
RBTree.prototype.itemByIndex = function(i) {
    var item;
    var j = 0;
    this.each(function(data) {
        if (j === i) {
            item = data;
            return false;
        }
        j++;
    });
    return item;
};

/**
 * @typedef {Object} Bin
 * @property {number} mean The mean (or center) of the bin.
 * @property {number} count The total count of points within the bin.
 * @property {number} cumn The cumulative count of points within this bin and
 *      all previous bins. This property is not generally required.
 * @property {number} tss The total sum of squared deviations from the mean.
 *      This is essentially a recursive residual metric, and is *not* really
 *      an approximation.
 */

/**
 * Order two bins by their means.
 * @param {Bin} a A histogram bin.
 * @param {Bin} b A histogram bin.
 * @return {number} res Difference between bin means.
 */
function compareBins(a, b) {
    // Part of Step 5 in Algorithm 1 from Ben-Haim & Tom-Tov (2010) p. 851
    return (a.mean > b.mean) ? 1 : (a.mean < b.mean) ? -1 : 0;
};

/**
 * Order two bins by their cumulative counts.
 * @param {Bin} a A histogram bin.
 * @param {Bin} b A histogram bin.
 * @return {number} diff Difference between bin cumulative means.
 */
function compareBinCumns(a, b) {
    // Part of Step 3 in Algorithm 4 from Ben-Haim & Tom-Tov (2010) p. 853
    return (a.cumn - b.cumn);
};

/**
 * Compute the (possibly weighted) difference between two bin means.
 * @param {Bin} a A histogram bin.
 * @param {Bin} b A histogram bin.
 * @param {boolean} weighted Whether to use gap weighting when comaring means.
 * @return {number} diff The (possibly weighted) difference between bin means.
 */
function diffBins(a, b, weighted) {
    // Part of Step 6 in Algorithm 1 from Ben-Haim & Tom-Tov (2010) p. 851
    var diff = Math.pow(b.mean - a.mean, 2)
    if (weighted) {
        diff *= math.log(math.E + math.min(a.count, b.count));
    }
    return diff
};

/**
 * Merges two bins based on a weighted average of their means.
 * @param {Bin} a A histogram bin.
 * @param {Bin} b A histogram bin.
 * @return {Bin} a Bin a merged with Bin b.
 */
function combineBins(a, b) {
    // Part of Step 7 in Algorithm 1 from Ben-Haim & Tom-Tov (2010) p. 851
    var numerator = (b.mean * b.count) + (a.mean * a.count),
        mean = numerator / (a.count + b.count),
        tss = (Math.pow(a.mean - mean, 2) * a.count) +
              (Math.pow(b.mean - mean, 2) * b.count);
    a.mean = mean
    a.count += b.count
    a.tss += b.tss + tss
    return a
};

/**
 * StreamHist class for building a streaming approximate histogram.
 * @param {number} [maxbins=100] The maximum number of bins used to approximate
 *      the data. This should be a positive integer.
 * @param {boolean} [weighted=false] When weighted is true, the histogram is
 *      encouraged to spend more of its bins capturing the densest areas of the
 *      distribution. For a Gaussian distribution this means better resolution
 *      near the mean and less resolution near the tails.
 * @param {number} [freeze=null] After the number of inserts into the histogram
 *      have exceeded the freeze parameter, the histogram bins are locked into
 *      place. As the bin means no longer shift, inserts become computationally
 *      cheap. However the quality of the histogram can suffer if the freeze
 *      parameter is too small. This should be a positive integer, or null to
 *      disable.
 * @constructor
 */
function StreamHist(maxbins, weighted, freeze, warmUp) {
    if (!(this instanceof StreamHist)) // Protect the global namespace!
        return new StreamHist(maxbins, weighted, freeze, warmUp)
    /**
     * Internal Red-Black Search Tree for storing histogram bins.
     * @type {RBTree}
     * @protected
     */
    this._bins = new RBTree(compareBins);
    // Primary parameters, defined in class doc.
    this.maxBins(maxbins || 100);
    this.weighted(weighted || false);
    this.freeze(freeze || 0);
    this.warmUp(warmUp || 0);
    this.reset();
};

/**
 * Adds a point or array of points to the histogram.
 * @param {(number|Array.<number>)} p Point or points to add to the histogram.
 * @param {number} [count=1] The 'weight' to use for the input point(s).
 * @return {StreamHist} this This histogram instance.
 */
StreamHist.prototype.push = function(p, count) {
    // Algorithm 1: Update Procedure
    // Ben-Haim & Tom-Tov (2010) p 851
    // Algorithm input requires a histogram h, and a point p
    // NOTE: Unlike in the paper, we allow for (positive?) integer counts
    count = count || 1;
    p = Array.isArray(p) ? p : [p];
    for (var i = 0 ; i < p.length ; i++) {
        this._insert(p[i], count);
        this._compress()
        if (this.count() === this.warmUp()) {  // Should only happen once...
            this.bins().each(b => {
                b.count = 1;
                b.tss = 0.0;
            })
        }
    }
    // Return a histogram with B bins that represents the set S ∪ {p}, where
    // S is the set represented by h.
    return this
};

/** @protected */
StreamHist.prototype._insert = function(x, count) {
    // Algorithm 1: Update Procedure
    // Ben-Haim & Tom-Tov (2010) p 851
    // NOTE: Unlike in the paper, we track min, max, and count separately
    this._min = math.min(this._min !== null ? this._min : Infinity,
                         x !== null ? x : Infinity);
    this._max = math.max(this._max !== null ? this._max : -Infinity,
                         x !== null ? x : -Infinity);
    this._count += count;
    this._tss = 0.0;  // Reset this because things may have changed...?
    // Steps 1-2: if p == p_i for some i then...
    var nearest = this.findNearest(x);
    if (nearest !== null && (nearest.mean === x ||
        (this.isFrozen() && this.size() === this.maxBins()))) {
        nearest.count += count  // m_i = m_i + 1 (or count in our case)
    } else {
        // Steps 3-5: Add the bin (p, 1) to the histogram, & sort the sequence
        this._newBin(x, count);
    }
};

StreamHist.prototype.tss = function() {
    if (!this._tss) {
        var tss = 0.0;
        this.bins().each(b => tss += b.tss)
        this._tss = tss
    }
    return this._tss
};

/**
 * Find the histogram bin nearest to an input point.
 * @param {number} p Input point for which to find the nearest histogram bin.
 * @return {Bin} nearest The nearest histogram bin to point p.
 */
StreamHist.prototype.findNearest = function(p) {
    // Algorithm 1: Update Procedure
    // Ben-Haim & Tom-Tov (2010) p 851
    // Steps 1-2: if p == p_i for some i then...
    if (this.size() === 0) {
        return null;
    }
    var iter = this.bins().lowerBound({mean:p});  // x <= iter || iter==null
    var nearest = (iter.data() === null) ? iter.prev() : iter.data();
    if (nearest.mean === p) {
        return nearest;
    }
    var prev = iter.prev();
    if (prev && math.abs(prev.mean - p) < math.abs(nearest.mean - p)) {
        return prev;
    } else {
        return nearest;
    }
};

/** @protected */
StreamHist.prototype._compress = function() {
    // Algorithm 1: Update Procedure
    // Ben-Haim & Tom-Tov (2010) p 851
    // Step 6: Find a point q_i that minimizes q_i+1 −q_i.
    var it, b, bi, bi1, a, diff, maxDiff;
    while (this.size() > this.maxBins()) {
        it = this.bins().iterator();
        a = it.next(), maxDiff = Infinity, diff;
        while((b = it.next()) !== null) {
            diff = diffBins(a, b, this.weighted());
            if (diff < maxDiff) {
                maxDiff = diff;
                bi = a, bi1 = b;
            }
            a = b;
        }
        // Step 7: Replace the bins (q_i, k_i), (q_i+1, k_i+1) by the bin:
        // (q_i*k_i + q_{i+1}*k_{i+1} / k_i + K_{i+1)}, k_i + k_{i+1})
        a = combineBins(bi, bi1);
        this.bins().remove(bi1);
    }
};

/** @protected */
StreamHist.prototype._newBin = function(p, count) {
    // Algorithm 1: Update Procedure
    // Ben-Haim & Tom-Tov (2010) p 851
    // Steps 3-5: Add the bin (p, 1) to the histogram, & sort the sequence
    var b = {mean:p, count:count, tss:0.0};
    this.bins().insert(b);
    return b;
};

/**
 * Merge/union two histograms to produce a single output histogram.
 * @param {StreamHist} that The other histogram object to merge into this one.
 * @param {number} maxBins The maximum number of bins to keep after merging.
 * @return {StreamHist} this This histogram instance.
 */
StreamHist.prototype.merge = function(that, maxBins) {
    // Algorithm 2: Merge Procedure
    // Ben-Haim & Tom-Tov (2010) p 852
    // Algorithm input requires two histograms (this and that) (h_1, h_2),
    // and an integer (maxBins) (B)
    // NOTE: Unlike in the paper, we track min, max, and count separately
    this._count += that.count();
    this._min = math.min(this.min() !== null ? this.min() : Infinity,
                         that.min() !== null ? that.min() : Infinity);
    this._max = math.max(this.max() !== null ? this.max() : -Infinity,
                         that.max() !== null ? that.max() : -Infinity);
    // NOTE: Unlike in the paper, B must be <= B_1 + B_2. Default is to use
    // min(B_1, B_2)
    var minBins = math.min(this.maxBins(), that.maxBins());
    var maxBins = Number.isInteger(maxBins) ? maxBins : minBins;
    this.maxBins(math.min(maxBins, this.maxBins() + that.maxBins()));
    // Steps 1-2: Add all bins from h_2 to h_1 and sort the sequence
    var it = that.bins().iterator(), item;
    while((item = it.next()) !== null) {
        this.bins().insert(item);
    }
    // Steps 3-6: Find a point q_i that minimizes q_i+1 − q_i.
    // Replace the bins (q_i, k_i), (q_i+1, k_i+1) by the bin:
    // (q_i*k_i + q_{i+1}*k_{i+1} / k_i + K_{i+1)}, k_i + k_{i+1})
    this._compress();
    // Return a histogram with B bins that represents the set S_1 ∪ S_2, where
    // S_1 and S_2 are the sets represented by h_1 and h_2, respectively.
    return this
};

/**
 * Compute the estimated data value for the given quantile(s).
 * The requested quantile(s) must be between 0 and 1.
 * @param {(number|Array.<number>)} p_or_plist The quantile or array of
 *      quantiles for which to estimate data values.
 * @return {(number|Array.<number>)} qs The data value(s) at the given
 *      quantile(s).
 */
StreamHist.prototype.quantile = function(p_or_plist) {
    // (A variation on) Algorithm 4: Uniform Procedure
    // Ben-Haim & Tom-Tov (2010) p 853
    // Algorithm input requires a Histogram h, and a quantile/percentage p
    // NOTE: This triggers _cumulate() if cumn's are out of date.
    var ps = Array.isArray(p_or_plist) ? p_or_plist : [p_or_plist];
    // Steps 1-6: ∀ j = 1, ..., B (or in our case, ∀ p ∈ ps)
    var qs = ps.map(this._quantile, this);
    // For percentage p {0, ..., 1}, or ∀ p ∈ ps, return
    // a real number u with the property that the number of points <= u
    // is p * ∑_{i=1}^B m_i (total count of points in h)
    return Array.isArray(p_or_plist) ? qs : qs[0];
};

/** @protected */
StreamHist.prototype._quantile = function(q) {
    // (A variation on) Algorithm 4: Uniform Procedure
    // Ben-Haim & Tom-Tov (2010) p 853
    // Some shortcuts to speed up (literal) edge cases
    if (this.size() === 0)
        return null;
    else if (q <= 0.0)
        return this.min();
    else if (q >= 1.0)
        return this.max();
    // Step 2: Set s = j/B ∑_{i=1}^B m_i (or in our case, p * ∑_{i=1}^B m_i)
    var s = this.count() * q;
    // Step 3: Find i such that sum([−∞, p_i ]) < s < sum([−∞, p_i+1 ]).
    // In other words, find bins that bracket s and interpolate s's cumulative
    // sum from their cumulative sums's.
    this._cumulate();
    var bound = this._boundCumn(s);
    var lower = bound[0], upper = bound[1];
    // Step 4: Set d to be the difference between s and sum([−∞, p_i]).
    var d = (lower.cumn - s);
    // Step 5
    // Search for u_j such that
    // d = ((m_i + m_{uj}) / 2) * ((u_j − p_i) / p_{i+1} - p_i)
    // where
    // m_{uj} = m_i + ((m_{i+1} − m_i) / (p_{i+1} - p_i))* (u_j − p_i).
    var a = upper.count - lower.count,  // a = m_{i+1} − m_i
        z = 0.0;
    // NOTE: Unlike in the paper, we handle some edge cases where a == 0,
    // and/or we can do a simpler interpolation...
    if (a == 0) {
        var b = upper.cumn - lower.cumn;
        if (b != 0) {
            z = (s - lower.cumn) / b;
        }  // Otherwise, stick with z = 0 (i.e., no interpolation)
    } else {
        // substituting z = (u_j − p_i) / (p_{i+1} − p_i)
        // we obtain a quadratic equation a * z^2 + b* z + c = 0, with...
        // a as above, and...
        var b = 2 * lower.count;  // b = 2_m
        var c = 2 * d;  // c = -2_d
        // z = -b + √(b^2 - 4ac) / 2a
        z = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a);
    }
    // Hence set u_j = p_i + (p_{i+1} − p_i)* z, where z is above...
    return lower.mean + (upper.mean - lower.mean) * z;
};

/**
 * Estimate values from this histogram's empirical cumulative distribution.
 * @param {(number|Array.<number>)} p_or_plist The value or array of
 *      values at which to estimate the cumulative *count*.
 * @return {(number|Array.<number>)} qs Cumulative count(s) at the given value(s).
 */
StreamHist.prototype.sum = function(p_or_plist, normed) {
    // (A variation on) Algorithm 3: Sum Procedure
    // Ben-Haim & Tom-Tov (2010) p 852
    // Algorithm input requires a Histogram h, and a value p
    var normed = normed == null ? true : normed === true;
    // NOTE: This triggers _cumulate() if cumn's are out of date.
    var ps = Array.isArray(p_or_plist) ? p_or_plist : [p_or_plist];
    // Steps 1-6: ∀ j = 1, ..., B (or in our case, ∀ p ∈ ps)
    var qs = ps.map(this._sum, this);
    // For percentage p {0, ..., 1}, or ∀ p ∈ ps, return
    // a real number u with the property that the number of points <= u
    // is p * ∑_{i=1}^B m_i (total count of points in h)
    // if (normed)  // If false, this is B-H & T-T's Algorithm 3 (sum)
    //     qs = math.divide(qs, this.count());
    return Array.isArray(p_or_plist) ? qs : qs[0];
};

/** @protected */
StreamHist.prototype._sum = function(b) {
    // (A variation on) Algorithm 3: Sum Procedure
    // Ben-Haim & Tom-Tov (2010) p 852
    // Some shortcuts to speed up (literal) edge cases
    if (this.size() === 0)
        return null;
    else if (b < this.min())
        return 0.0;
    else if (b >= this.max())
        return this.count();
    // Step 1: Find i such that p_i < b < p_{i+1}
    // In other words, find bins that bracket b and interpolate b's cumulative
    // sum from their cumulative sums's.
    this._cumulate();
    var bound = this._bound(b),
        lower = bound[0], upper = bound[1];
    // Step 2: Set s = ((m_i + m_b) / 2) * ((b − p_i) / (p_{i+1} − p_i))
    // where m_b = m_i + (m_{i+1} - m_i) / (p_{i+1} - p_i) * (b − p_i).
    var pdiff = (upper.mean - lower.mean),
        bdiff = (b - lower.mean);
    if (bdiff === 0)
        return lower.cumn
    var mdiff = (upper.count - lower.count),
        mb = lower.count + (mdiff / pdiff)  * bdiff,
        s = (lower === this.bins().min()) ? lower.count/2: lower.cumn;
    // Steps 3-6: s = ∑_{j=i-1}^0 m_j + m_i/2
    return s + ((lower.count + mb) / 2) * (bdiff / pdiff)
};

/**
 * Estimate values from this histogram's empirical probability distribution.
 * NOTE: This is based on interpolating a PDF between discrete bins, we
 * *aren't* using Dirac delta functions, which might be a better idea.
 * @see {@link https://www.probabilitycourse.com/chapter4/4_3_2
 *      _delta_function.php}
 * @param {(number|Array.<number>)} p_or_plist The value or array of
 *      values at which to estimate the density.
 * @return {(number|Array.<number>)} qs Mass/density for the given value(s).
 */
StreamHist.prototype.density = function(b) {
    // Despite quite a few tests, this should be treated as experimental
    // Some shortcuts to speed up (literal) edge cases
    if (this.size() === 0)
        return null;
    else if (b < this.min() || b > this.max())
        return 0.0;
    else if (b === this.min() && b === this.max())  // Single point mass
        return Infinity;
    var bound = this._bound(b),
        lower = bound[0], upper = bound[1];
    if (lower.mean === b && upper.mean === b) { // p is exactly a bin mean
        // We use 2* EPSILON because... well I'm not sure, it seems to work.
        var lo = b - Number.EPSILON*10, hi = b + Number.EPSILON*10,
            res = (this.density(lo) + this.density(hi)) / 2;
    } else {
        var pdiff = (upper.mean - lower.mean),
            bdiff = (b - lower.mean),
            mdiff = (upper.count - lower.count),
            ratio = bdiff / pdiff,
            mb = mdiff * ratio,
            res = (lower.count + mb) / pdiff;
    }
    return res >= 0 ? res/this.count() : 0.0
};

/** @protected */
StreamHist.prototype._cumulate = function() {
    // (A variation on) Algorithm 4: Uniform Procedure
    // Step 3: Find i such that sum([−∞, pi ]) < s < sum([−∞, pi+1 ]).
    // Update cumulative counts for each bin for use in above Step 3 search
    if (this._count === this._cumn) {
        return;
    }
    var cumn = 0;
    this.bins().each(c => {
        c.cumn = cumn + c.count / 2;  // half of count at the mean
        cumn += c.count;
    });
    this._count = this._cumn = cumn;
};

/** @protected */
StreamHist.prototype._bound = function(x) {
    // Find bins lower and upper such that lower.mean < x <
    // upper.mean or lower.mean === x === upper.mean.
    // NOTE: Don't call this for x out of bounds.
    var iter = this.bins().lowerBound({mean:x});  // x < iter
    // lower <= x
    var lo = (iter.data() === this.bins().min()) ? iter.data() : iter.prev();
    // upper > x
    var up = (lo.mean === x || lo === this.bins().max()) ? lo : iter.next();
    return [lo, up];
};

/** @protected */
StreamHist.prototype._boundCumn = function(x) {
    // Find bins lower and upper such that lower.cumn < x <
    // upper.cumn or lower.cumn === x === upper.cumn. Don't call
    // this for cumn out of bounds.
    // NOTE: because mean and cumn give rise to the same sort order
    // (up to identical means), use the mean rbtree for our search.
    this.bins()._comparator = compareBinCumns;  // Update comparator
    var iter = this.bins().upperBound({cumn: x});  // cumn < iter
    this.bins()._comparator = compareBins;  // Put it back!
    var lo = iter.prev();  // lower <= cumn
    var up = (lo && lo.mean === x || lo === this.bins().max()) ? lo : iter.next();
    return [lo, up];
};

/** @protected */
StreamHist.prototype._addWeight = function(nearest, x, count) {
    // Add weight at location x to nearest bin. Adding x to
    // nearest will not shift its relative position in the tree and
    // require reinsertion.
    // This isn't required for the Ben-Haim & Tom-Tov version, but _is_ part
    // of the Guedalia et al. paper.
    if (x !== nearest.mean) {
        nearest.mean += count * (x - nearest.mean) / (nearest.count + count);
    }
    // nearest.count += count;
};

/**
 * Reset this histogram so that all counts, bins, etc. are set to zero/null.
 * @return {StreamHist} this This (now reset) histogram instance.
 */
StreamHist.prototype.reset = function() {
    // Reset and prepare to consume new points
    this._bins.clear();
    this._count = 0;
    this._cumn = 0;
    this._min = null;
    this._max = null;
    this._tss = 0.0;
    return this;
};

/**
 * Return this histogram's bins as an array of bin objects.
 * @return {Array.<Bin>} array Array of bin objects ordered by their means.
 */
StreamHist.prototype.toArray = function() {
    var array = [];
    this.bins().each(b => array.push({mean:b.mean, count:b.count, tss:b.tss}));
    return array;
};

/**
 * Return a JSON reprentation of this histogram useful for communicating
 * histogram information across processes or networks.
 * @return {object} json A JSON reprentation of this histogram.
 */
StreamHist.prototype.toJSON = function() {
    var json = {};
    json._bins = this.toArray();
    for (var name in this) {
        if (this.hasOwnProperty(name) &&
            (name !== "_bins" && name !== "_cumn")) {
            json[name] = this[name];
        }
    }
    return json;
};

/**
 * Create a histogram from a JSON data structure.
 * @param {object} json A JSON representation of a histogram.
 * @return {StreamHist} json A JSON reprentation of this histogram.
 * @see StreamHist.toJSON
 */
StreamHist.fromJSON = function(json) {
    var obj = new StreamHist();
    for (var name in json) {
        if (name !== "_bins") {
            obj[name] = json[name];
        }
    }
    // Insert bins from plain array to RBTree
    json._bins.forEach((bin) => obj._bins.insert(bin));
    return obj;
};

/**
 * Return boolean indicating whether or not this histogram is currently frozen.
 * @return {boolean} frozen Is this histogram frozen?
 */
StreamHist.prototype.isFrozen = function() {
    return this.freeze() !== 0 && this.count() > this.freeze();
}

/**
 * Set or get whether this histogram uses gap weighting.
 * @return {boolean} weighted Is gap weighting turned on?
 */
StreamHist.prototype.weighted = function(weighted) {
    if (weighted != null)
        this._weighted = (weighted === true);
    return this._weighted === true;
};

/**
 * Set or get the maximum number of bins this histogram will use.
 * @return {number} maxBins The current maximum number of bins to use.
 */
StreamHist.prototype.maxBins = function(maxBins) {
    if (maxBins != null) {
        this._maxBins = Number.isInteger(maxBins) ? maxBins : 100;
        if (this.size() > this._maxBins)
            this._compress();
    }
    return this._maxBins;
};

/**
 * Set or get the number of inserts into the histogram before the bins are
 * locked into place.
 * @return {number} freeze The current freeze threshold.
 */
StreamHist.prototype.freeze = function(freeze) {
    if (freeze != null)
        this._freeze = Number.isInteger(freeze) ? freeze : 0;
    return this._freeze || 0;
};

/**
 * Set or get the number of warm-up inserts into the histogram.
 * After warmUp inserts have been made, the counts etc. are set to zero.
 * @return {number} warmUp The current warmUp threshold.
 */
StreamHist.prototype.warmUp = function(warmUp) {
    if (warmUp != null)
        this._warmUp = Number.isInteger(warmUp) ? warmUp : 0;
    return this._warmUp || 0;
};

/**
 * Return this histogram's bins.
 * @return {RBTree} bins The internal RBTree containing the histogram's bins.
 */
StreamHist.prototype.bins = function() {
    return this._bins;
};

/**
 * Return this histogram's total count/number of inserts.
 * @return {number} count The total count/number of inserts.
 */
StreamHist.prototype.count = function() {
    // h_c = ∑_{i=1}^B m_i) (using notation from paper)
    return this._count;
};

/**
 * Return the current number of bins in this histogram.
 * @return {number} size The size of the internal RBTree of histogram bins.
 */
StreamHist.prototype.size = function() {
    // This is a wrapper around the histogram's RBTree size
    return this.bins().size;
};

/**
 * Return the minumum and maximum values that this histogram has seen.
 * @return {Array.<number>} limits The extremes of the incoming data stream.
 */
StreamHist.prototype.limits = function () {
    // Min and max are stored outside the histogram bins
    return [this.min(), this.max()];
};

/**
 * Return the maximum value this histogram has seen.
 * @return {number} max The maximum value in the incoming data stream.
 */
StreamHist.prototype.max = function () {
    // Min and max are stored outside the histogram bins
    return this._max;
};

/**
 * Return the minimum value this histogram has seen.
 * @return {number} min The minimum value in the incoming data stream.
 */
StreamHist.prototype.min = function () {
    // Min and max are stored outside the histogram bins
    return this._min;
};

/**
 * Return the estimated mean of the histogram's underlying distribution.
 * @return {number} mean The estimated mean.
 */
StreamHist.prototype.mean = function () {
    if (this.count() === 0)
        return null
    var it = this.bins().iterator(), sum = 0.0, item;
    while((item = it.next()) !== null) {
        sum += (item.mean * item.count);
    }
    return sum / this.count();
};

/**
 * Return the estimated variance of the histogram's underlying distribution.
 * @return {number} variance The estimated variance.
 */
StreamHist.prototype.variance = function () {
    if (this.count() < 2)
        return null
    var it = this.bins().iterator(), sum = 0.0, mean = this.mean(), item;
    while((item = it.next()) !== null) {
        sum += (math.pow(item.mean - mean, 2) * item.count);
    }
    return sum / this.count();
};

/**
 * Return the estimated stand. dev. of the histogram's underlying distribution.
 * @return {number} std The estimated standard deviation.
 */
StreamHist.prototype.std = function() {
    var variance = this.variance()
    return variance ? math.sqrt(variance) : null;
}

/**
 * Return the estimated median of the histogram's underlying distribution.
 * @return {number} median The estimated median.
 */
StreamHist.prototype.median = function(p_or_plist) {
    if (this.count() > this.size()) {
        return this.quantile(0.5);
    } else {  // Return the 'exact' median when possible
        var iter = this.bins().lowerBound({mean:this.count()/2});
        if (this.count() % 2 === 0) {
            var upper = iter.data(), lower = iter.prev();
            return combineBins(lower, upper).mean;
        } else {
            return iter.prev().mean;
        }
    }
};

/**
 * Return a summary of the histogram's underlying distribution.
 * @return {object} summary An object with properties describing various
 * summary statistics. Properties include the count, mean, std, min, 1st, 2nd,
 * and 3rd quartiles, and the max.
 */
StreamHist.prototype.summary = function() {
    var summary = {
        "count": this.count(),
        "mean": this.mean(),
        "std": this.std(),
        "min": this.min(),
        "Q1": this.quantile(0.25),
        "Q2": this.median(),
        "Q3": this.quantile(0.75),
        "max": this.max()
    };
    return summary;
};

StreamHist.prototype.toString = function() {
    var total = this.count();
    var string = "";
    this.bins().each(b => {
        string += b.mean.toFixed(4) + '\t' + '.'.repeat(math.trunc(b.count/total*200)) + '\n';
    });
    return string;
};

/**
 * Compute a streaming histogram with k bins 'on the fly' from an input array.
 * @param {Array.<number>} input An array of input data points.
 * @param {number} k The maximum number of bins to use to approximate the data
 *      distribution. This should be a positive integer.
 * @return {Array.<Bin>} array Array of bin objects ordered by their means.
 */
function fastHist(input, k) {
    return new StreamHist(k).push(input).toArray();
}

module.exports = {'StreamHist': StreamHist, 'fastHist': fastHist,
                  'diffBins': diffBins, 'combineBins': combineBins}
