var chai = require('chai'),
    should = chai.should(),
    expect = chai.expect,
    StreamHist = require('./index.js').StreamHist,
    rand = require('randgen'),
    present = require('present'),
    math = require('mathjs');
chai.use(require("chai-deep-closeto"));

// Array prototype function for aiding in split/apply/combine workflow
Array.prototype.split = function (groupsize) {
    var sets = [];
    var chunks = this.length/groupsize;
    for (var i = 0, j = 0; i < chunks; i++, j += groupsize) {
        sets[i] = this.slice(j, j + groupsize);
    }
    return sets;
};

describe('A basic StreamHist object', function() {

    var hist, data = rand.rvnorm(10000),
        mean = math.mean(data), median = math.median(data),
        std = math.std(data), min = math.min(data),
        max = math.max(data), maxBins = 50, length = data.length,
        unif = rand.rvunif(length);
    // http://stackoverflow.com/a/33352604 (for following range trick)
    var range = Array.from(Array(15).keys());
    // or... Array.from(Array(15)).map((e, i) => i);

    before(function() {
        hist = new StreamHist(maxBins);
    });

    afterEach(function() {
        // Reset to defaults...
        hist.reset();
        hist.maxBins(50);
        hist.weighted(false);
    });

    it('should have reasonable defaults', function() {
        var empty = new StreamHist()
        empty.maxBins().should.equal(100)
        empty.weighted().should.be.false;
        empty.bins().should.have.property('_comparator');  // RBTree property
        empty.bins().size.should.equal(0);
        empty.count().should.equal(0);
        expect(empty.min()).to.be.null;
        expect(empty.min()).to.be.null;
        empty.size().should.equal(0);
    });
    it('should be able to be reset', function() {
        expect(hist.min()).to.be.null;
        expect(hist.max()).to.be.null;
        hist.size().should.equal(0);
        hist.count().should.equal(0);
        hist.push(data);
        hist.min().should.equal(min);
        hist.max().should.equal(max);
        hist.size().should.equal(maxBins);
        hist.count().should.equal(length);
    });
    it('should store the max/min/limits', function() {
        hist.push(range);
        hist.limits().should.deep.equal([0, 14]);
        hist.min().should.equal(0);
        hist.max().should.equal(14);
    });
    it('should store the (current) total count', function() {
        hist.push(data);
        hist.count().should.equal(length);
    });
    it('should store the number of bins', function() {
        hist.push(data);
        hist.size().should.equal(50);
        hist.size().should.equal(hist.maxBins());
        hist.size().should.equal(hist.bins().size)
    });
    it('should compute the mean', function() {
        hist.push(data).mean().should.be.closeTo(mean, 0.05);
        hist.reset();
        hist.push(range).mean().should.equal(7)
    });
    it('should compute the median', function() {
        // Maybe too liberal with the margin or error?
        hist.push(data).median().should.be.closeTo(median, 0.08);
    });
    it('should compute an exact median for n < maxbins', function() {
        // Even number of points
        hist.push(range).median().should.equal(7)
        // Odd number of points
        hist.push(15).median().should.equal(7.5)
    });
    it('should compute variance/std dev unless n < 2', function() {
        expect(hist.push(0).variance()).to.be.null;
        expect(hist.std()).to.be.null;
        hist.reset();
        // Maybe too liberal with the margin or error?
        hist.push(data).variance().should.be.closeTo(Math.pow(std, 2), 0.05);
        hist.std().should.be.closeTo(std, 0.05);
    });
    it('should produce an Array via toArray', function() {
        hist.push(range);
        var array = range.map(d => Object({mean:d, count:1, tss:0}));
        hist.toArray().should.deep.equal(array);
    });
    it('should produce a JSON object via toJSON', function() {
        hist.push(range);
        var obj = {_maxBins: hist.maxBins(),
                   _min: hist.min(), _max: hist.max(),
                   _weighted: false, _freeze: 0, _tss: 0,
                   _count: hist.count(), _warmUp: 0}
        obj._bins = range.map(d => Object({mean:d, count:1, tss:0}));
        var json = hist.toJSON();
        json.should.deep.equal(obj);
        json.should.be.instanceof(Object);
    });
    it('should convert to and from (round trip) JSON', function () {
        var fromJSON = StreamHist.fromJSON(hist.toJSON());
        // First, check bins... since the structure of the RBTree might vary
        // depending on build order, the arrays should be the same, but we
        // can't deep compare the RBTrees directly...
        fromJSON.toArray().should.deep.equal(hist.toArray())
        // Next, check all other properties...
        for (var name in fromJSON) {
            if (fromJSON.hasOwnProperty(name) && name !== "_bins") {
                expect(fromJSON[name]).to.deep.equal(hist[name]);
            }
        }
    });
    it('should push points and arrays into the structure', function() {
        expect(hist.push(1).size()).to.be.equal(1);
        hist.reset();
        hist.push(range).size().should.equal(range.length);
        hist.reset();
        hist.push(data)  // Push array
            .push(0.0)   // Push points...
            .push(0.5)
            .push(1.0).mean().should.be.closeTo(mean, 0.05);
    });
    it('should allow different "weights" or counts to be used', function() {
        hist.push(data)  // Push array
            .push(0.0)   // Push points...
            .push(0.5, 10)
            .push(1.0).count().should.equal(length + 12);
    });
    it('should be able to find the nearest bin to an input point', function() {
        hist.push(data);
        hist.findNearest(min).should.equal(hist.bins().min());
        hist.findNearest(max).should.equal(hist.bins().max());
        // Need a better version of the 'middle' bin...
        // hist.findNearest(mean).mean.should.be.closeTo(mean, 0.05);
    });
    it('should "compress" automatically when adding points', function() {
        this.timeout(2500);
        hist.size().should.equal(0);
        hist.maxBins(80);
        hist.push(data).size().should.equal(80);
        hist.size().should.be.below(length);
        // Now we'll do it manually
        hist.reset();
        // This is comparatively slow...
        for (var i = 0 ; i < length ; i++) {
            hist._insert(data[i], 1);  // Don't try this at home!
        }
        hist.size().should.equal(length);
        hist._compress();
        hist.size().should.equal(hist.maxBins());

    });
    it('should have larger counts at tails with weighted gaps', function() {
        // Histograms using weighted gaps are less eager to merge bins with
        // large counts. This test builds weighted and non-weighted
        // histograms using samples from a normal distribution. The
        // non-weighted histogram should spend more of its bins capturing
        // the tails of the distribution. With that in mind this test makes
        // sure the bins bracketing the weighted histogram have larger
        // counts than the bins bracketing the non-weighted histogram.

        var w = new StreamHist(16, weighted=true).push(data);
        var c = new StreamHist(16, weighted=false).push(data);
        (w.bins().min().count +
         w.bins().max().count).should.be.at.least(c.bins().min().count +
                                                  c.bins().max().count)
    });
    it('should reduce the number of bin merges when using freeze', function() {
        var frozen = new StreamHist(50, false, 500);
        var points = 10000;  // Large number needed to see difference
        frozen.isFrozen().should.be.false;  // First it should be false
        frozen.push(rand.rvnorm(points));
        // We still want to sketch to be mostly accurate
        frozen.median().should.be.closeTo(0.0, 0.05);
        frozen.mean().should.be.closeTo(0.0, 0.05);
        frozen.std().should.be.closeTo(1.0, 0.05);
        // frozen.sum().should.be.closeTo(points/2, points/50)
        frozen.isFrozen().should.be.true;  // Then it should be true!
    });
    it('should be fast to update, especially when frozen', function() {
        var frozen = new StreamHist(50, false, 500);
        var classic = new StreamHist(50);
        var size = 10000,
            points = rand.rvnorm(size);
        var start = present();
        frozen.push(points);
        var one = present() - start;
        start = present();
        classic.push(points);
        var two = present() - start;
        one.should.be.below(two);
    });
    it('should produce a summary of the underlying distribution', function(){
        hist.push(data);
        var summary = hist.summary();
        summary.should.have.property("count", length);
        summary.should.have.property("min", min);
        summary.should.have.property("max", max);
        summary.should.have.property("Q2");
        summary.Q2.should.be.closeTo(median, 0.05);
        summary.should.be.instanceof(Object)
    });
    it('should estimate arbitrary quantiles', function() {
        hist.push(unif).quantile(0.5).should.be.closeTo(0.5, 0.05)
        hist.reset();
        hist.push(data).quantile([0.25, 0.5, 0.75])
            .should.be.deep.closeTo([-0.66, 0.0, 0.66], 0.05)
        hist.quantile(0).should.equal(hist.min());
        hist.quantile(1).should.equal(hist.max());
    });
    it('should allow merging of multiple histograms', function() {
        var initial = new StreamHist(50);  // Empty histogram
        // Also illustrates the split/apply/combine or map/reduce
        // capabilities of the StreamHist objects...
        var merged = data
                        .split(100)
                        .map((arr) => new StreamHist(50).push(arr))
                        .reduce((prev, curr) => prev.merge(curr), initial);
        var single = new StreamHist(50).push(data);
        // This is really slick: http://stackoverflow.com/a/33352604
        var quantiles = Array.from(Array(9)).map((e, i) => (i+1) /10);
        merged.quantile(quantiles)
            .should.be.deep.closeTo(single.quantile(quantiles), 0.05);
    });
    it('should estimate reasonable CDF (sum) values at the limits', function() {
        hist.push(0).push(10);
        hist.sum(5).should.equal(1.0);
        hist.sum(0).should.equal(0.5);
        hist.sum(10).should.equal(2.0);
    });
    it('should estimate reasonable values for the CDF (sum)', function() {
        var scale = data.length/50, middle = data.length/2
        hist.push(data);
        hist.sum(median).should.be.closeTo(middle, scale);
        hist.sum(mean).should.be.closeTo(middle, scale);
        // 'Round-trips'... should be off due to rounding/float errors...
        hist.sum(hist.quantile(0.6)).should.be.closeTo(data.length*0.6, scale);
        hist.quantile(hist.sum(1.5)/hist.count()).should.be.closeTo(1.5, 0.05);
    });

    it('should estimate reasonable PDF values at the limits', function() {
        hist.reset().push([-1, 0, 1]).density(0).should.equal(1/3);
        hist.reset().push(0).density(0).should.equal(Infinity);
    });
    it('should estimate reasonable values for the density', function() {
        hist.reset().push([1, 2, 2, 3]);
        hist.density(0.0).should.be.closeTo(0.0, 1e-10);
        hist.density(0.5).should.be.closeTo(0.0, 1e-10);
        hist.density(1.5).should.be.closeTo(0.375, 1e-10);
        hist.density(2.0).should.be.closeTo(0.5, 1e-10);
        hist.density(2.5).should.be.closeTo(0.375, 1e-10);
        hist.density(3.0).should.be.closeTo(0.25, 1e-10);
        hist.density(3.5).should.be.closeTo(0.0, 1e-10);
        hist.density(4.0).should.be.closeTo(0.0, 1e-10);
    });
});
