(function (factory) {
    if (typeof define === 'function' && define.amd) {
        // AMD. Register as an anonymous module.
        define([], factory);
    } else if (typeof exports === 'object') {
        // Node. Does not work with strict CommonJS, but
        // only CommonJS-like enviroments that support module.exports,
        // like Node.
        module.exports = factory();
    } else {
        // Browser globals
        this.streamhist = factory();
    }
}(function () {
    var global = this, define;
    function _require(id) {
        var module = _require.cache[id];
        if (!module) {
            var exports = {};
            module = _require.cache[id] = {
                id: id,
                exports: exports
            };
            _require.modules[id].call(exports, module, exports);
        }
        return module.exports;
    }
    _require.cache = [];
    _require.modules = [
        function (module, exports) {
            //
            // StreamHist
            //
            // Streaming approximate historgrams for continuous univariate data.
            //
            // This module implements a streaming approximate classifier for continuous
            // univariate data which uses a streaming histogram data structure (sketch)
            // similar to the one proposed in [1] (and based on work from [2] and [3]), to
            // approximate the underlying data distribution. Results are similar to those
            // from Jenk's Natural Breaks algorithm [4], though at _significantly_ reduced
            // computational complexity and memory. The SAC classifier also provides an
            // algorithm for estimating the number of 'natural' breaks in the data that is
            // to be robust to variations in the scale and range of the underlying
            // distribution.
            //
            // The streaming historgram data structure is a modified version of the
            // structure presented in [1], which reduces the computational complexity from
            // O(K_{max}N) to O(log{K_{max}}N) without increasing the space complexity,
            // which is O(K_{max}). Additionally, the histogram can track _multiple_
            // statistics (in addition to the mean) to help guide the classification
            // process. This _does_ increase the memory requirements slightly, but improves
            // the ability to find 'natural breaks' in the data.
            //
            // The basic algorithm is essentially 'parameter-free', in that it doesn't
            // require the user to specify any information about the underlying data
            // distribution explicitly. There are, however, a few parameters that help
            // control the accuracy of the sketch, including the maximum number of bins,
            // and whether to use weighting when computing bin means. A larger maximum
            // number of bins yields a more accurate approximation of the underlying
            // distribution at the cost of increased memory utilization and performance.
            // There is no 'optimal' bin count, but less than 100 bins should generally be
            // sufficient. Bin weighting is useful for downplaying outliers in the dataset,
            // as well as providing a more accurate fit to the data in higher density areas
            // (i.e., more bins in higher density areas).
            // References
            // ----------
            // [1] http://jmlr.org/papers/volume11/ben-haim10a/ben-haim10a.pdf
            // [2] http://www.cs.huji.ac.il/~werman/Papers/guedalia_etal99.pdf
            // [3] http://dl.acm.org/citation.cfm?id=1044332
            // [4] https://en.wikipedia.org/wiki/Jenks_natural_breaks_optimization
            //
            // Copyright (c) 2016 Carson Farmer <carsonfarmer@gmail.com>
            // Portions inspired by tdigest (https://github.com/welch/tdigest)
            // Copyright (c) 2015 Will Welch
            //
            // TODO: Docs need to be updated to reflect the fact that this is the
            //       StreamHist module, _not_ a Classifier, which should be a separate
            //       module.
            // module.paths.push('/usr/local/Cellar/node/4.1.2/libexec/npm/lib/node_modules/');
            var RBTree = _require(2).RBTree, _ = _require(6);
            _.mixin(_require(1));
            function StreamHist(maxbins, weighted) {
                // Allocate a StreamHist structure.
                //
                // maxbins:  The maximum number of bins used to approximate the data
                // weighted: If falsey, a simple minimum distance is used to merge nearby
                //           bins, otherwise, an exponentially weighted distance is used.
                //
                this.maxbins = maxbins || 100;
                this.weighted = weighted || false;
                this.bins = new RBTree(compare_bins);
                this.reset();
            }
            StreamHist.prototype.reset = function () {
                // Reset and prepare to consume new points
                this.bins.clear();
                this.n = 0;
            };
            StreamHist.prototype.size = function () {
                return this.bins.size;
            };
            StreamHist.prototype.toArray = function (everything) {
                // return {mean, n} of centroids as an array ordered by mean.
                var result = [];
                this.bins.each(function (c) {
                    result.push({
                        mean: c.mean,
                        n: c.n
                    });
                });
                return result;
            };
            function compare_bins(a, b) {
                // Order two centroids by mean.
                return a.mean > b.mean ? 1 : a.mean < b.mean ? -1 : 0;
            }
            function diff_bins(a, b, weighted) {
                var diff = b.mean - a.mean;
                if (weighted) {
                    diff *= Math.log(Math.E + Math.min(a.n, b.n));
                }
                return diff;
            }
            function combine_bins(a, b) {
                a.mean = (a.mean * a.n + b.mean * b.n) / (a.n + b.n);
                a.n = a.n + b.n;
                return a;
            }
            StreamHist.prototype.push = function (x, n) {
                // Incorporate value or array of values x, having count n into the
                // histogram. n defaults to 1.
                n = n || 1;
                x = Array.isArray(x) ? x : [x];
                for (var i = 0; i < x.length; i++) {
                    this._insert(x[i], n);
                }
                return this;
            };
            StreamHist.prototype.find_nearest = function (x) {
                // Find the bin closest to x.
                if (this.size() === 0) {
                    return null;
                }
                var iter = this.bins.lowerBound({ mean: x });    // x <= iter || iter==null
                // x <= iter || iter==null
                var c = iter.data() === null ? iter.prev() : iter.data();
                if (c.mean === x) {
                    return c;
                }
                var prev = iter.prev();
                if (prev && Math.abs(prev.mean - x) < Math.abs(c.mean - x)) {
                    return prev;
                } else {
                    return c;
                }
            };
            StreamHist.prototype._new_bin = function (x, n) {
                // Create and insert a new centroid into the digest.
                var c = {
                        mean: x,
                        n: n
                    };
                this.bins.insert(c);
                this.n += n;
                return c;
            };
            StreamHist.prototype._add_weight = function (nearest, x, n) {
                // Add weight at location x to nearest bin. Adding x to
                // nearest will not shift its relative position in the tree and
                // require reinsertion.
                if (x !== nearest.mean) {
                    nearest.mean += n * (x - nearest.mean) / (nearest.n + n);
                }
                nearest.n += n;
                this.n += n;
            };
            StreamHist.prototype._insert = function (x, n) {
                // Incorporate value x, having count n into the histogram.
                var min = this.bins.min();
                var max = this.bins.max();
                var nearest = this.find_nearest(x);
                if (nearest && nearest.mean === x) {
                    // Accumulate exact matches into the bin without limit.
                    this._add_weight(nearest, x, n);
                } else {
                    this._new_bin(x, n);
                }
                this._compress();
            };
            StreamHist.prototype._compress = function () {
                // Compress histogram (back) down to maxbins bins.
                if (this.bins.size > this.maxbins) {
                    var it = this.bins.iterator(), b, bi, bi1, a = it.next(), max_diff = Infinity, diff;
                    while ((b = it.next()) !== null) {
                        diff = diff_bins(a, b, this.weighted);
                        if (diff < max_diff) {
                            max_diff = diff;
                            bi = a, bi1 = b;
                        }
                        a = b;
                    }
                    a = combine_bins(bi, bi1);
                    this.bins.remove(bi1);
                }
            };
            function fasthist(input, k) {
                var hist = new sh.StreamHist(k);
                return hist.push(input).toArray();
            }
            module.exports = {
                'StreamHist': StreamHist,
                'fasthist': fasthist
            };
        },
        function (module, exports) {
            var _ = _require(6);    // Internal function: creates a callback bound to its context if supplied
                                    // Copied from underscore.js: 
                                    // https://github.com/jashkenas/underscore/blob/master/underscore.js
            // Internal function: creates a callback bound to its context if supplied
            // Copied from underscore.js: 
            // https://github.com/jashkenas/underscore/blob/master/underscore.js
            var createCallback = function (func, context, argCount) {
                if (!context)
                    return func;
                switch (argCount == null ? 3 : argCount) {
                case 1:
                    return function (value) {
                        return func.call(context, value);
                    };
                case 2:
                    return function (value, other) {
                        return func.call(context, value, other);
                    };
                case 3:
                    return function (value, index, collection) {
                        return func.call(context, value, index, collection);
                    };
                case 4:
                    return function (accumulator, value, index, collection) {
                        return func.call(context, accumulator, value, index, collection);
                    };
                }
                return function () {
                    return func.apply(this, arguments);
                };
            };    // An internal function to generate lookup iterators. 
                  // Copied from underscore.js: 
                  // https://github.com/jashkenas/underscore/blob/master/underscore.js
            // An internal function to generate lookup iterators. 
            // Copied from underscore.js: 
            // https://github.com/jashkenas/underscore/blob/master/underscore.js
            var lookupIterator = function (value, context, argCount) {
                if (value == null)
                    return _.identity;
                if (_.isFunction(value))
                    return createCallback(value, context, argCount);
                if (_.isObject(value))
                    return _.matches(value);
                return _.property(value);
            };
            module.exports = {
                // return the index of the smallest element in the array "obj".
                // If "iterator" is given, it is called on each element and the result is used for comparison. 
                argmin: function (obj, iterator, context) {
                    var min = null, argmin = null, value;
                    if (iterator)
                        iterator = lookupIterator(iterator, context);
                    for (var i = 0, length = obj.length; i < length; i++) {
                        value = iterator ? iterator(obj[i]) : obj[i];
                        if (min == null || value < min) {
                            min = value;
                            argmin = i;
                        }
                    }
                    ;
                    return argmin;
                },
                // return the index of the largest element in the array "obj".
                // If "iterator" is given, it is called on each element and the result is used for comparison. 
                argmax: function (obj, iterator, context) {
                    var max = null, argmax = null, value;
                    if (iterator)
                        iterator = lookupIterator(iterator, context);
                    for (var i = 0, length = obj.length; i < length; i++) {
                        value = iterator ? iterator(obj[i]) : obj[i];
                        if (max == null || value > max) {
                            max = value;
                            argmax = i;
                        }
                    }
                    ;
                    return argmax;
                }
            };
        },
        function (module, exports) {
            module.exports = {
                RBTree: _require(4),
                BinTree: _require(3)
            };
        },
        function (module, exports) {
            var TreeBase = _require(5);
            function Node(data) {
                this.data = data;
                this.left = null;
                this.right = null;
            }
            Node.prototype.get_child = function (dir) {
                return dir ? this.right : this.left;
            };
            Node.prototype.set_child = function (dir, val) {
                if (dir) {
                    this.right = val;
                } else {
                    this.left = val;
                }
            };
            function BinTree(comparator) {
                this._root = null;
                this._comparator = comparator;
                this.size = 0;
            }
            BinTree.prototype = new TreeBase();    // returns true if inserted, false if duplicate
            // returns true if inserted, false if duplicate
            BinTree.prototype.insert = function (data) {
                if (this._root === null) {
                    // empty tree
                    this._root = new Node(data);
                    this.size++;
                    return true;
                }
                var dir = 0;    // setup
                // setup
                var p = null;    // parent
                // parent
                var node = this._root;    // search down
                // search down
                while (true) {
                    if (node === null) {
                        // insert new node at the bottom
                        node = new Node(data);
                        p.set_child(dir, node);
                        ret = true;
                        this.size++;
                        return true;
                    }    // stop if found
                    // stop if found
                    if (this._comparator(node.data, data) === 0) {
                        return false;
                    }
                    dir = this._comparator(node.data, data) < 0;    // update helpers
                    // update helpers
                    p = node;
                    node = node.get_child(dir);
                }
            };    // returns true if removed, false if not found
            // returns true if removed, false if not found
            BinTree.prototype.remove = function (data) {
                if (this._root === null) {
                    return false;
                }
                var head = new Node(undefined);    // fake tree root
                // fake tree root
                var node = head;
                node.right = this._root;
                var p = null;    // parent
                // parent
                var found = null;    // found item
                // found item
                var dir = 1;
                while (node.get_child(dir) !== null) {
                    p = node;
                    node = node.get_child(dir);
                    var cmp = this._comparator(data, node.data);
                    dir = cmp > 0;
                    if (cmp === 0) {
                        found = node;
                    }
                }
                if (found !== null) {
                    found.data = node.data;
                    p.set_child(p.right === node, node.get_child(node.left === null));
                    this._root = head.right;
                    this.size--;
                    return true;
                } else {
                    return false;
                }
            };
            module.exports = BinTree;
        },
        function (module, exports) {
            var TreeBase = _require(5);
            function Node(data) {
                this.data = data;
                this.left = null;
                this.right = null;
                this.red = true;
            }
            Node.prototype.get_child = function (dir) {
                return dir ? this.right : this.left;
            };
            Node.prototype.set_child = function (dir, val) {
                if (dir) {
                    this.right = val;
                } else {
                    this.left = val;
                }
            };
            function RBTree(comparator) {
                this._root = null;
                this._comparator = comparator;
                this.size = 0;
            }
            RBTree.prototype = new TreeBase();    // returns true if inserted, false if duplicate
            // returns true if inserted, false if duplicate
            RBTree.prototype.insert = function (data) {
                var ret = false;
                if (this._root === null) {
                    // empty tree
                    this._root = new Node(data);
                    ret = true;
                    this.size++;
                } else {
                    var head = new Node(undefined);    // fake tree root
                    // fake tree root
                    var dir = 0;
                    var last = 0;    // setup
                    // setup
                    var gp = null;    // grandparent
                    // grandparent
                    var ggp = head;    // grand-grand-parent
                    // grand-grand-parent
                    var p = null;    // parent
                    // parent
                    var node = this._root;
                    ggp.right = this._root;    // search down
                    // search down
                    while (true) {
                        if (node === null) {
                            // insert new node at the bottom
                            node = new Node(data);
                            p.set_child(dir, node);
                            ret = true;
                            this.size++;
                        } else if (is_red(node.left) && is_red(node.right)) {
                            // color flip
                            node.red = true;
                            node.left.red = false;
                            node.right.red = false;
                        }    // fix red violation
                        // fix red violation
                        if (is_red(node) && is_red(p)) {
                            var dir2 = ggp.right === gp;
                            if (node === p.get_child(last)) {
                                ggp.set_child(dir2, single_rotate(gp, !last));
                            } else {
                                ggp.set_child(dir2, double_rotate(gp, !last));
                            }
                        }
                        var cmp = this._comparator(node.data, data);    // stop if found
                        // stop if found
                        if (cmp === 0) {
                            break;
                        }
                        last = dir;
                        dir = cmp < 0;    // update helpers
                        // update helpers
                        if (gp !== null) {
                            ggp = gp;
                        }
                        gp = p;
                        p = node;
                        node = node.get_child(dir);
                    }    // update root
                    // update root
                    this._root = head.right;
                }    // make root black
                // make root black
                this._root.red = false;
                return ret;
            };    // returns true if removed, false if not found
            // returns true if removed, false if not found
            RBTree.prototype.remove = function (data) {
                if (this._root === null) {
                    return false;
                }
                var head = new Node(undefined);    // fake tree root
                // fake tree root
                var node = head;
                node.right = this._root;
                var p = null;    // parent
                // parent
                var gp = null;    // grand parent
                // grand parent
                var found = null;    // found item
                // found item
                var dir = 1;
                while (node.get_child(dir) !== null) {
                    var last = dir;    // update helpers
                    // update helpers
                    gp = p;
                    p = node;
                    node = node.get_child(dir);
                    var cmp = this._comparator(data, node.data);
                    dir = cmp > 0;    // save found node
                    // save found node
                    if (cmp === 0) {
                        found = node;
                    }    // push the red node down
                    // push the red node down
                    if (!is_red(node) && !is_red(node.get_child(dir))) {
                        if (is_red(node.get_child(!dir))) {
                            var sr = single_rotate(node, dir);
                            p.set_child(last, sr);
                            p = sr;
                        } else if (!is_red(node.get_child(!dir))) {
                            var sibling = p.get_child(!last);
                            if (sibling !== null) {
                                if (!is_red(sibling.get_child(!last)) && !is_red(sibling.get_child(last))) {
                                    // color flip
                                    p.red = false;
                                    sibling.red = true;
                                    node.red = true;
                                } else {
                                    var dir2 = gp.right === p;
                                    if (is_red(sibling.get_child(last))) {
                                        gp.set_child(dir2, double_rotate(p, last));
                                    } else if (is_red(sibling.get_child(!last))) {
                                        gp.set_child(dir2, single_rotate(p, last));
                                    }    // ensure correct coloring
                                    // ensure correct coloring
                                    var gpc = gp.get_child(dir2);
                                    gpc.red = true;
                                    node.red = true;
                                    gpc.left.red = false;
                                    gpc.right.red = false;
                                }
                            }
                        }
                    }
                }    // replace and remove if found
                // replace and remove if found
                if (found !== null) {
                    found.data = node.data;
                    p.set_child(p.right === node, node.get_child(node.left === null));
                    this.size--;
                }    // update root and make it black
                // update root and make it black
                this._root = head.right;
                if (this._root !== null) {
                    this._root.red = false;
                }
                return found !== null;
            };
            function is_red(node) {
                return node !== null && node.red;
            }
            function single_rotate(root, dir) {
                var save = root.get_child(!dir);
                root.set_child(!dir, save.get_child(dir));
                save.set_child(dir, root);
                root.red = true;
                save.red = false;
                return save;
            }
            function double_rotate(root, dir) {
                root.set_child(!dir, single_rotate(root.get_child(!dir), !dir));
                return single_rotate(root, dir);
            }
            module.exports = RBTree;
        },
        function (module, exports) {
            function TreeBase() {
            }    // removes all nodes from the tree
            // removes all nodes from the tree
            TreeBase.prototype.clear = function () {
                this._root = null;
                this.size = 0;
            };    // returns node data if found, null otherwise
            // returns node data if found, null otherwise
            TreeBase.prototype.find = function (data) {
                var res = this._root;
                while (res !== null) {
                    var c = this._comparator(data, res.data);
                    if (c === 0) {
                        return res.data;
                    } else {
                        res = res.get_child(c > 0);
                    }
                }
                return null;
            };    // returns iterator to node if found, null otherwise
            // returns iterator to node if found, null otherwise
            TreeBase.prototype.findIter = function (data) {
                var res = this._root;
                var iter = this.iterator();
                while (res !== null) {
                    var c = this._comparator(data, res.data);
                    if (c === 0) {
                        iter._cursor = res;
                        return iter;
                    } else {
                        iter._ancestors.push(res);
                        res = res.get_child(c > 0);
                    }
                }
                return null;
            };    // Returns an iterator to the tree node at or immediately after the item
            // Returns an iterator to the tree node at or immediately after the item
            TreeBase.prototype.lowerBound = function (item) {
                var cur = this._root;
                var iter = this.iterator();
                var cmp = this._comparator;
                while (cur !== null) {
                    var c = cmp(item, cur.data);
                    if (c === 0) {
                        iter._cursor = cur;
                        return iter;
                    }
                    iter._ancestors.push(cur);
                    cur = cur.get_child(c > 0);
                }
                for (var i = iter._ancestors.length - 1; i >= 0; --i) {
                    cur = iter._ancestors[i];
                    if (cmp(item, cur.data) < 0) {
                        iter._cursor = cur;
                        iter._ancestors.length = i;
                        return iter;
                    }
                }
                iter._ancestors.length = 0;
                return iter;
            };    // Returns an iterator to the tree node immediately after the item
            // Returns an iterator to the tree node immediately after the item
            TreeBase.prototype.upperBound = function (item) {
                var iter = this.lowerBound(item);
                var cmp = this._comparator;
                while (iter.data() !== null && cmp(iter.data(), item) === 0) {
                    iter.next();
                }
                return iter;
            };    // returns null if tree is empty
            // returns null if tree is empty
            TreeBase.prototype.min = function () {
                var res = this._root;
                if (res === null) {
                    return null;
                }
                while (res.left !== null) {
                    res = res.left;
                }
                return res.data;
            };    // returns null if tree is empty
            // returns null if tree is empty
            TreeBase.prototype.max = function () {
                var res = this._root;
                if (res === null) {
                    return null;
                }
                while (res.right !== null) {
                    res = res.right;
                }
                return res.data;
            };    // returns a null iterator
                  // call next() or prev() to point to an element
            // returns a null iterator
            // call next() or prev() to point to an element
            TreeBase.prototype.iterator = function () {
                return new Iterator(this);
            };    // calls cb on each node's data, in order
            // calls cb on each node's data, in order
            TreeBase.prototype.each = function (cb) {
                var it = this.iterator(), data;
                while ((data = it.next()) !== null) {
                    cb(data);
                }
            };    // calls cb on each node's data, in reverse order
            // calls cb on each node's data, in reverse order
            TreeBase.prototype.reach = function (cb) {
                var it = this.iterator(), data;
                while ((data = it.prev()) !== null) {
                    cb(data);
                }
            };
            function Iterator(tree) {
                this._tree = tree;
                this._ancestors = [];
                this._cursor = null;
            }
            Iterator.prototype.data = function () {
                return this._cursor !== null ? this._cursor.data : null;
            };    // if null-iterator, returns first node
                  // otherwise, returns next node
            // if null-iterator, returns first node
            // otherwise, returns next node
            Iterator.prototype.next = function () {
                if (this._cursor === null) {
                    var root = this._tree._root;
                    if (root !== null) {
                        this._minNode(root);
                    }
                } else {
                    if (this._cursor.right === null) {
                        // no greater node in subtree, go up to parent
                        // if coming from a right child, continue up the stack
                        var save;
                        do {
                            save = this._cursor;
                            if (this._ancestors.length) {
                                this._cursor = this._ancestors.pop();
                            } else {
                                this._cursor = null;
                                break;
                            }
                        } while (this._cursor.right === save);
                    } else {
                        // get the next node from the subtree
                        this._ancestors.push(this._cursor);
                        this._minNode(this._cursor.right);
                    }
                }
                return this._cursor !== null ? this._cursor.data : null;
            };    // if null-iterator, returns last node
                  // otherwise, returns previous node
            // if null-iterator, returns last node
            // otherwise, returns previous node
            Iterator.prototype.prev = function () {
                if (this._cursor === null) {
                    var root = this._tree._root;
                    if (root !== null) {
                        this._maxNode(root);
                    }
                } else {
                    if (this._cursor.left === null) {
                        var save;
                        do {
                            save = this._cursor;
                            if (this._ancestors.length) {
                                this._cursor = this._ancestors.pop();
                            } else {
                                this._cursor = null;
                                break;
                            }
                        } while (this._cursor.left === save);
                    } else {
                        this._ancestors.push(this._cursor);
                        this._maxNode(this._cursor.left);
                    }
                }
                return this._cursor !== null ? this._cursor.data : null;
            };
            Iterator.prototype._minNode = function (start) {
                while (start.left !== null) {
                    this._ancestors.push(start);
                    start = start.left;
                }
                this._cursor = start;
            };
            Iterator.prototype._maxNode = function (start) {
                while (start.right !== null) {
                    this._ancestors.push(start);
                    start = start.right;
                }
                this._cursor = start;
            };
            module.exports = TreeBase;
        },
        function (module, exports) {
            (function () {
                //     Underscore.js 1.8.3
                //     http://underscorejs.org
                //     (c) 2009-2015 Jeremy Ashkenas, DocumentCloud and Investigative Reporters & Editors
                //     Underscore may be freely distributed under the MIT license.
                // Baseline setup
                // --------------
                // Establish the root object, `window` in the browser, or `exports` on the server.
                var root = this;    // Save the previous value of the `_` variable.
                // Save the previous value of the `_` variable.
                var previousUnderscore = root._;    // Save bytes in the minified (but not gzipped) version:
                // Save bytes in the minified (but not gzipped) version:
                var ArrayProto = Array.prototype, ObjProto = Object.prototype, FuncProto = Function.prototype;    // Create quick reference variables for speed access to core prototypes.
                // Create quick reference variables for speed access to core prototypes.
                var push = ArrayProto.push, slice = ArrayProto.slice, toString = ObjProto.toString, hasOwnProperty = ObjProto.hasOwnProperty;    // All **ECMAScript 5** native function implementations that we hope to use
                                                                                                                                                 // are declared here.
                // All **ECMAScript 5** native function implementations that we hope to use
                // are declared here.
                var nativeIsArray = Array.isArray, nativeKeys = Object.keys, nativeBind = FuncProto.bind, nativeCreate = Object.create;    // Naked function reference for surrogate-prototype-swapping.
                // Naked function reference for surrogate-prototype-swapping.
                var Ctor = function () {
                };    // Create a safe reference to the Underscore object for use below.
                // Create a safe reference to the Underscore object for use below.
                var _ = function (obj) {
                    if (obj instanceof _)
                        return obj;
                    if (!(this instanceof _))
                        return new _(obj);
                    this._wrapped = obj;
                };    // Export the Underscore object for **Node.js**, with
                      // backwards-compatibility for the old `require()` API. If we're in
                      // the browser, add `_` as a global object.
                // Export the Underscore object for **Node.js**, with
                // backwards-compatibility for the old `require()` API. If we're in
                // the browser, add `_` as a global object.
                if (typeof exports !== 'undefined') {
                    if (typeof module !== 'undefined' && module.exports) {
                        exports = module.exports = _;
                    }
                    exports._ = _;
                } else {
                    root._ = _;
                }    // Current version.
                // Current version.
                _.VERSION = '1.8.3';    // Internal function that returns an efficient (for current engines) version
                                        // of the passed-in callback, to be repeatedly applied in other Underscore
                                        // functions.
                // Internal function that returns an efficient (for current engines) version
                // of the passed-in callback, to be repeatedly applied in other Underscore
                // functions.
                var optimizeCb = function (func, context, argCount) {
                    if (context === void 0)
                        return func;
                    switch (argCount == null ? 3 : argCount) {
                    case 1:
                        return function (value) {
                            return func.call(context, value);
                        };
                    case 2:
                        return function (value, other) {
                            return func.call(context, value, other);
                        };
                    case 3:
                        return function (value, index, collection) {
                            return func.call(context, value, index, collection);
                        };
                    case 4:
                        return function (accumulator, value, index, collection) {
                            return func.call(context, accumulator, value, index, collection);
                        };
                    }
                    return function () {
                        return func.apply(context, arguments);
                    };
                };    // A mostly-internal function to generate callbacks that can be applied
                      // to each element in a collection, returning the desired result — either
                      // identity, an arbitrary callback, a property matcher, or a property accessor.
                // A mostly-internal function to generate callbacks that can be applied
                // to each element in a collection, returning the desired result — either
                // identity, an arbitrary callback, a property matcher, or a property accessor.
                var cb = function (value, context, argCount) {
                    if (value == null)
                        return _.identity;
                    if (_.isFunction(value))
                        return optimizeCb(value, context, argCount);
                    if (_.isObject(value))
                        return _.matcher(value);
                    return _.property(value);
                };
                _.iteratee = function (value, context) {
                    return cb(value, context, Infinity);
                };    // An internal function for creating assigner functions.
                // An internal function for creating assigner functions.
                var createAssigner = function (keysFunc, undefinedOnly) {
                    return function (obj) {
                        var length = arguments.length;
                        if (length < 2 || obj == null)
                            return obj;
                        for (var index = 1; index < length; index++) {
                            var source = arguments[index], keys = keysFunc(source), l = keys.length;
                            for (var i = 0; i < l; i++) {
                                var key = keys[i];
                                if (!undefinedOnly || obj[key] === void 0)
                                    obj[key] = source[key];
                            }
                        }
                        return obj;
                    };
                };    // An internal function for creating a new object that inherits from another.
                // An internal function for creating a new object that inherits from another.
                var baseCreate = function (prototype) {
                    if (!_.isObject(prototype))
                        return {};
                    if (nativeCreate)
                        return nativeCreate(prototype);
                    Ctor.prototype = prototype;
                    var result = new Ctor();
                    Ctor.prototype = null;
                    return result;
                };
                var property = function (key) {
                    return function (obj) {
                        return obj == null ? void 0 : obj[key];
                    };
                };    // Helper for collection methods to determine whether a collection
                      // should be iterated as an array or as an object
                      // Related: http://people.mozilla.org/~jorendorff/es6-draft.html#sec-tolength
                      // Avoids a very nasty iOS 8 JIT bug on ARM-64. #2094
                // Helper for collection methods to determine whether a collection
                // should be iterated as an array or as an object
                // Related: http://people.mozilla.org/~jorendorff/es6-draft.html#sec-tolength
                // Avoids a very nasty iOS 8 JIT bug on ARM-64. #2094
                var MAX_ARRAY_INDEX = Math.pow(2, 53) - 1;
                var getLength = property('length');
                var isArrayLike = function (collection) {
                    var length = getLength(collection);
                    return typeof length == 'number' && length >= 0 && length <= MAX_ARRAY_INDEX;
                };    // Collection Functions
                      // --------------------
                      // The cornerstone, an `each` implementation, aka `forEach`.
                      // Handles raw objects in addition to array-likes. Treats all
                      // sparse array-likes as if they were dense.
                // Collection Functions
                // --------------------
                // The cornerstone, an `each` implementation, aka `forEach`.
                // Handles raw objects in addition to array-likes. Treats all
                // sparse array-likes as if they were dense.
                _.each = _.forEach = function (obj, iteratee, context) {
                    iteratee = optimizeCb(iteratee, context);
                    var i, length;
                    if (isArrayLike(obj)) {
                        for (i = 0, length = obj.length; i < length; i++) {
                            iteratee(obj[i], i, obj);
                        }
                    } else {
                        var keys = _.keys(obj);
                        for (i = 0, length = keys.length; i < length; i++) {
                            iteratee(obj[keys[i]], keys[i], obj);
                        }
                    }
                    return obj;
                };    // Return the results of applying the iteratee to each element.
                // Return the results of applying the iteratee to each element.
                _.map = _.collect = function (obj, iteratee, context) {
                    iteratee = cb(iteratee, context);
                    var keys = !isArrayLike(obj) && _.keys(obj), length = (keys || obj).length, results = Array(length);
                    for (var index = 0; index < length; index++) {
                        var currentKey = keys ? keys[index] : index;
                        results[index] = iteratee(obj[currentKey], currentKey, obj);
                    }
                    return results;
                };    // Create a reducing function iterating left or right.
                // Create a reducing function iterating left or right.
                function createReduce(dir) {
                    // Optimized iterator function as using arguments.length
                    // in the main function will deoptimize the, see #1991.
                    function iterator(obj, iteratee, memo, keys, index, length) {
                        for (; index >= 0 && index < length; index += dir) {
                            var currentKey = keys ? keys[index] : index;
                            memo = iteratee(memo, obj[currentKey], currentKey, obj);
                        }
                        return memo;
                    }
                    return function (obj, iteratee, memo, context) {
                        iteratee = optimizeCb(iteratee, context, 4);
                        var keys = !isArrayLike(obj) && _.keys(obj), length = (keys || obj).length, index = dir > 0 ? 0 : length - 1;    // Determine the initial value if none is provided.
                        // Determine the initial value if none is provided.
                        if (arguments.length < 3) {
                            memo = obj[keys ? keys[index] : index];
                            index += dir;
                        }
                        return iterator(obj, iteratee, memo, keys, index, length);
                    };
                }    // **Reduce** builds up a single result from a list of values, aka `inject`,
                     // or `foldl`.
                // **Reduce** builds up a single result from a list of values, aka `inject`,
                // or `foldl`.
                _.reduce = _.foldl = _.inject = createReduce(1);    // The right-associative version of reduce, also known as `foldr`.
                // The right-associative version of reduce, also known as `foldr`.
                _.reduceRight = _.foldr = createReduce(-1);    // Return the first value which passes a truth test. Aliased as `detect`.
                // Return the first value which passes a truth test. Aliased as `detect`.
                _.find = _.detect = function (obj, predicate, context) {
                    var key;
                    if (isArrayLike(obj)) {
                        key = _.findIndex(obj, predicate, context);
                    } else {
                        key = _.findKey(obj, predicate, context);
                    }
                    if (key !== void 0 && key !== -1)
                        return obj[key];
                };    // Return all the elements that pass a truth test.
                      // Aliased as `select`.
                // Return all the elements that pass a truth test.
                // Aliased as `select`.
                _.filter = _.select = function (obj, predicate, context) {
                    var results = [];
                    predicate = cb(predicate, context);
                    _.each(obj, function (value, index, list) {
                        if (predicate(value, index, list))
                            results.push(value);
                    });
                    return results;
                };    // Return all the elements for which a truth test fails.
                // Return all the elements for which a truth test fails.
                _.reject = function (obj, predicate, context) {
                    return _.filter(obj, _.negate(cb(predicate)), context);
                };    // Determine whether all of the elements match a truth test.
                      // Aliased as `all`.
                // Determine whether all of the elements match a truth test.
                // Aliased as `all`.
                _.every = _.all = function (obj, predicate, context) {
                    predicate = cb(predicate, context);
                    var keys = !isArrayLike(obj) && _.keys(obj), length = (keys || obj).length;
                    for (var index = 0; index < length; index++) {
                        var currentKey = keys ? keys[index] : index;
                        if (!predicate(obj[currentKey], currentKey, obj))
                            return false;
                    }
                    return true;
                };    // Determine if at least one element in the object matches a truth test.
                      // Aliased as `any`.
                // Determine if at least one element in the object matches a truth test.
                // Aliased as `any`.
                _.some = _.any = function (obj, predicate, context) {
                    predicate = cb(predicate, context);
                    var keys = !isArrayLike(obj) && _.keys(obj), length = (keys || obj).length;
                    for (var index = 0; index < length; index++) {
                        var currentKey = keys ? keys[index] : index;
                        if (predicate(obj[currentKey], currentKey, obj))
                            return true;
                    }
                    return false;
                };    // Determine if the array or object contains a given item (using `===`).
                      // Aliased as `includes` and `include`.
                // Determine if the array or object contains a given item (using `===`).
                // Aliased as `includes` and `include`.
                _.contains = _.includes = _.include = function (obj, item, fromIndex, guard) {
                    if (!isArrayLike(obj))
                        obj = _.values(obj);
                    if (typeof fromIndex != 'number' || guard)
                        fromIndex = 0;
                    return _.indexOf(obj, item, fromIndex) >= 0;
                };    // Invoke a method (with arguments) on every item in a collection.
                // Invoke a method (with arguments) on every item in a collection.
                _.invoke = function (obj, method) {
                    var args = slice.call(arguments, 2);
                    var isFunc = _.isFunction(method);
                    return _.map(obj, function (value) {
                        var func = isFunc ? method : value[method];
                        return func == null ? func : func.apply(value, args);
                    });
                };    // Convenience version of a common use case of `map`: fetching a property.
                // Convenience version of a common use case of `map`: fetching a property.
                _.pluck = function (obj, key) {
                    return _.map(obj, _.property(key));
                };    // Convenience version of a common use case of `filter`: selecting only objects
                      // containing specific `key:value` pairs.
                // Convenience version of a common use case of `filter`: selecting only objects
                // containing specific `key:value` pairs.
                _.where = function (obj, attrs) {
                    return _.filter(obj, _.matcher(attrs));
                };    // Convenience version of a common use case of `find`: getting the first object
                      // containing specific `key:value` pairs.
                // Convenience version of a common use case of `find`: getting the first object
                // containing specific `key:value` pairs.
                _.findWhere = function (obj, attrs) {
                    return _.find(obj, _.matcher(attrs));
                };    // Return the maximum element (or element-based computation).
                // Return the maximum element (or element-based computation).
                _.max = function (obj, iteratee, context) {
                    var result = -Infinity, lastComputed = -Infinity, value, computed;
                    if (iteratee == null && obj != null) {
                        obj = isArrayLike(obj) ? obj : _.values(obj);
                        for (var i = 0, length = obj.length; i < length; i++) {
                            value = obj[i];
                            if (value > result) {
                                result = value;
                            }
                        }
                    } else {
                        iteratee = cb(iteratee, context);
                        _.each(obj, function (value, index, list) {
                            computed = iteratee(value, index, list);
                            if (computed > lastComputed || computed === -Infinity && result === -Infinity) {
                                result = value;
                                lastComputed = computed;
                            }
                        });
                    }
                    return result;
                };    // Return the minimum element (or element-based computation).
                // Return the minimum element (or element-based computation).
                _.min = function (obj, iteratee, context) {
                    var result = Infinity, lastComputed = Infinity, value, computed;
                    if (iteratee == null && obj != null) {
                        obj = isArrayLike(obj) ? obj : _.values(obj);
                        for (var i = 0, length = obj.length; i < length; i++) {
                            value = obj[i];
                            if (value < result) {
                                result = value;
                            }
                        }
                    } else {
                        iteratee = cb(iteratee, context);
                        _.each(obj, function (value, index, list) {
                            computed = iteratee(value, index, list);
                            if (computed < lastComputed || computed === Infinity && result === Infinity) {
                                result = value;
                                lastComputed = computed;
                            }
                        });
                    }
                    return result;
                };    // Shuffle a collection, using the modern version of the
                      // [Fisher-Yates shuffle](http://en.wikipedia.org/wiki/Fisher–Yates_shuffle).
                // Shuffle a collection, using the modern version of the
                // [Fisher-Yates shuffle](http://en.wikipedia.org/wiki/Fisher–Yates_shuffle).
                _.shuffle = function (obj) {
                    var set = isArrayLike(obj) ? obj : _.values(obj);
                    var length = set.length;
                    var shuffled = Array(length);
                    for (var index = 0, rand; index < length; index++) {
                        rand = _.random(0, index);
                        if (rand !== index)
                            shuffled[index] = shuffled[rand];
                        shuffled[rand] = set[index];
                    }
                    return shuffled;
                };    // Sample **n** random values from a collection.
                      // If **n** is not specified, returns a single random element.
                      // The internal `guard` argument allows it to work with `map`.
                // Sample **n** random values from a collection.
                // If **n** is not specified, returns a single random element.
                // The internal `guard` argument allows it to work with `map`.
                _.sample = function (obj, n, guard) {
                    if (n == null || guard) {
                        if (!isArrayLike(obj))
                            obj = _.values(obj);
                        return obj[_.random(obj.length - 1)];
                    }
                    return _.shuffle(obj).slice(0, Math.max(0, n));
                };    // Sort the object's values by a criterion produced by an iteratee.
                // Sort the object's values by a criterion produced by an iteratee.
                _.sortBy = function (obj, iteratee, context) {
                    iteratee = cb(iteratee, context);
                    return _.pluck(_.map(obj, function (value, index, list) {
                        return {
                            value: value,
                            index: index,
                            criteria: iteratee(value, index, list)
                        };
                    }).sort(function (left, right) {
                        var a = left.criteria;
                        var b = right.criteria;
                        if (a !== b) {
                            if (a > b || a === void 0)
                                return 1;
                            if (a < b || b === void 0)
                                return -1;
                        }
                        return left.index - right.index;
                    }), 'value');
                };    // An internal function used for aggregate "group by" operations.
                // An internal function used for aggregate "group by" operations.
                var group = function (behavior) {
                    return function (obj, iteratee, context) {
                        var result = {};
                        iteratee = cb(iteratee, context);
                        _.each(obj, function (value, index) {
                            var key = iteratee(value, index, obj);
                            behavior(result, value, key);
                        });
                        return result;
                    };
                };    // Groups the object's values by a criterion. Pass either a string attribute
                      // to group by, or a function that returns the criterion.
                // Groups the object's values by a criterion. Pass either a string attribute
                // to group by, or a function that returns the criterion.
                _.groupBy = group(function (result, value, key) {
                    if (_.has(result, key))
                        result[key].push(value);
                    else
                        result[key] = [value];
                });    // Indexes the object's values by a criterion, similar to `groupBy`, but for
                       // when you know that your index values will be unique.
                // Indexes the object's values by a criterion, similar to `groupBy`, but for
                // when you know that your index values will be unique.
                _.indexBy = group(function (result, value, key) {
                    result[key] = value;
                });    // Counts instances of an object that group by a certain criterion. Pass
                       // either a string attribute to count by, or a function that returns the
                       // criterion.
                // Counts instances of an object that group by a certain criterion. Pass
                // either a string attribute to count by, or a function that returns the
                // criterion.
                _.countBy = group(function (result, value, key) {
                    if (_.has(result, key))
                        result[key]++;
                    else
                        result[key] = 1;
                });    // Safely create a real, live array from anything iterable.
                // Safely create a real, live array from anything iterable.
                _.toArray = function (obj) {
                    if (!obj)
                        return [];
                    if (_.isArray(obj))
                        return slice.call(obj);
                    if (isArrayLike(obj))
                        return _.map(obj, _.identity);
                    return _.values(obj);
                };    // Return the number of elements in an object.
                // Return the number of elements in an object.
                _.size = function (obj) {
                    if (obj == null)
                        return 0;
                    return isArrayLike(obj) ? obj.length : _.keys(obj).length;
                };    // Split a collection into two arrays: one whose elements all satisfy the given
                      // predicate, and one whose elements all do not satisfy the predicate.
                // Split a collection into two arrays: one whose elements all satisfy the given
                // predicate, and one whose elements all do not satisfy the predicate.
                _.partition = function (obj, predicate, context) {
                    predicate = cb(predicate, context);
                    var pass = [], fail = [];
                    _.each(obj, function (value, key, obj) {
                        (predicate(value, key, obj) ? pass : fail).push(value);
                    });
                    return [
                        pass,
                        fail
                    ];
                };    // Array Functions
                      // ---------------
                      // Get the first element of an array. Passing **n** will return the first N
                      // values in the array. Aliased as `head` and `take`. The **guard** check
                      // allows it to work with `_.map`.
                // Array Functions
                // ---------------
                // Get the first element of an array. Passing **n** will return the first N
                // values in the array. Aliased as `head` and `take`. The **guard** check
                // allows it to work with `_.map`.
                _.first = _.head = _.take = function (array, n, guard) {
                    if (array == null)
                        return void 0;
                    if (n == null || guard)
                        return array[0];
                    return _.initial(array, array.length - n);
                };    // Returns everything but the last entry of the array. Especially useful on
                      // the arguments object. Passing **n** will return all the values in
                      // the array, excluding the last N.
                // Returns everything but the last entry of the array. Especially useful on
                // the arguments object. Passing **n** will return all the values in
                // the array, excluding the last N.
                _.initial = function (array, n, guard) {
                    return slice.call(array, 0, Math.max(0, array.length - (n == null || guard ? 1 : n)));
                };    // Get the last element of an array. Passing **n** will return the last N
                      // values in the array.
                // Get the last element of an array. Passing **n** will return the last N
                // values in the array.
                _.last = function (array, n, guard) {
                    if (array == null)
                        return void 0;
                    if (n == null || guard)
                        return array[array.length - 1];
                    return _.rest(array, Math.max(0, array.length - n));
                };    // Returns everything but the first entry of the array. Aliased as `tail` and `drop`.
                      // Especially useful on the arguments object. Passing an **n** will return
                      // the rest N values in the array.
                // Returns everything but the first entry of the array. Aliased as `tail` and `drop`.
                // Especially useful on the arguments object. Passing an **n** will return
                // the rest N values in the array.
                _.rest = _.tail = _.drop = function (array, n, guard) {
                    return slice.call(array, n == null || guard ? 1 : n);
                };    // Trim out all falsy values from an array.
                // Trim out all falsy values from an array.
                _.compact = function (array) {
                    return _.filter(array, _.identity);
                };    // Internal implementation of a recursive `flatten` function.
                // Internal implementation of a recursive `flatten` function.
                var flatten = function (input, shallow, strict, startIndex) {
                    var output = [], idx = 0;
                    for (var i = startIndex || 0, length = getLength(input); i < length; i++) {
                        var value = input[i];
                        if (isArrayLike(value) && (_.isArray(value) || _.isArguments(value))) {
                            //flatten current level of array or arguments object
                            if (!shallow)
                                value = flatten(value, shallow, strict);
                            var j = 0, len = value.length;
                            output.length += len;
                            while (j < len) {
                                output[idx++] = value[j++];
                            }
                        } else if (!strict) {
                            output[idx++] = value;
                        }
                    }
                    return output;
                };    // Flatten out an array, either recursively (by default), or just one level.
                // Flatten out an array, either recursively (by default), or just one level.
                _.flatten = function (array, shallow) {
                    return flatten(array, shallow, false);
                };    // Return a version of the array that does not contain the specified value(s).
                // Return a version of the array that does not contain the specified value(s).
                _.without = function (array) {
                    return _.difference(array, slice.call(arguments, 1));
                };    // Produce a duplicate-free version of the array. If the array has already
                      // been sorted, you have the option of using a faster algorithm.
                      // Aliased as `unique`.
                // Produce a duplicate-free version of the array. If the array has already
                // been sorted, you have the option of using a faster algorithm.
                // Aliased as `unique`.
                _.uniq = _.unique = function (array, isSorted, iteratee, context) {
                    if (!_.isBoolean(isSorted)) {
                        context = iteratee;
                        iteratee = isSorted;
                        isSorted = false;
                    }
                    if (iteratee != null)
                        iteratee = cb(iteratee, context);
                    var result = [];
                    var seen = [];
                    for (var i = 0, length = getLength(array); i < length; i++) {
                        var value = array[i], computed = iteratee ? iteratee(value, i, array) : value;
                        if (isSorted) {
                            if (!i || seen !== computed)
                                result.push(value);
                            seen = computed;
                        } else if (iteratee) {
                            if (!_.contains(seen, computed)) {
                                seen.push(computed);
                                result.push(value);
                            }
                        } else if (!_.contains(result, value)) {
                            result.push(value);
                        }
                    }
                    return result;
                };    // Produce an array that contains the union: each distinct element from all of
                      // the passed-in arrays.
                // Produce an array that contains the union: each distinct element from all of
                // the passed-in arrays.
                _.union = function () {
                    return _.uniq(flatten(arguments, true, true));
                };    // Produce an array that contains every item shared between all the
                      // passed-in arrays.
                // Produce an array that contains every item shared between all the
                // passed-in arrays.
                _.intersection = function (array) {
                    var result = [];
                    var argsLength = arguments.length;
                    for (var i = 0, length = getLength(array); i < length; i++) {
                        var item = array[i];
                        if (_.contains(result, item))
                            continue;
                        for (var j = 1; j < argsLength; j++) {
                            if (!_.contains(arguments[j], item))
                                break;
                        }
                        if (j === argsLength)
                            result.push(item);
                    }
                    return result;
                };    // Take the difference between one array and a number of other arrays.
                      // Only the elements present in just the first array will remain.
                // Take the difference between one array and a number of other arrays.
                // Only the elements present in just the first array will remain.
                _.difference = function (array) {
                    var rest = flatten(arguments, true, true, 1);
                    return _.filter(array, function (value) {
                        return !_.contains(rest, value);
                    });
                };    // Zip together multiple lists into a single array -- elements that share
                      // an index go together.
                // Zip together multiple lists into a single array -- elements that share
                // an index go together.
                _.zip = function () {
                    return _.unzip(arguments);
                };    // Complement of _.zip. Unzip accepts an array of arrays and groups
                      // each array's elements on shared indices
                // Complement of _.zip. Unzip accepts an array of arrays and groups
                // each array's elements on shared indices
                _.unzip = function (array) {
                    var length = array && _.max(array, getLength).length || 0;
                    var result = Array(length);
                    for (var index = 0; index < length; index++) {
                        result[index] = _.pluck(array, index);
                    }
                    return result;
                };    // Converts lists into objects. Pass either a single array of `[key, value]`
                      // pairs, or two parallel arrays of the same length -- one of keys, and one of
                      // the corresponding values.
                // Converts lists into objects. Pass either a single array of `[key, value]`
                // pairs, or two parallel arrays of the same length -- one of keys, and one of
                // the corresponding values.
                _.object = function (list, values) {
                    var result = {};
                    for (var i = 0, length = getLength(list); i < length; i++) {
                        if (values) {
                            result[list[i]] = values[i];
                        } else {
                            result[list[i][0]] = list[i][1];
                        }
                    }
                    return result;
                };    // Generator function to create the findIndex and findLastIndex functions
                // Generator function to create the findIndex and findLastIndex functions
                function createPredicateIndexFinder(dir) {
                    return function (array, predicate, context) {
                        predicate = cb(predicate, context);
                        var length = getLength(array);
                        var index = dir > 0 ? 0 : length - 1;
                        for (; index >= 0 && index < length; index += dir) {
                            if (predicate(array[index], index, array))
                                return index;
                        }
                        return -1;
                    };
                }    // Returns the first index on an array-like that passes a predicate test
                // Returns the first index on an array-like that passes a predicate test
                _.findIndex = createPredicateIndexFinder(1);
                _.findLastIndex = createPredicateIndexFinder(-1);    // Use a comparator function to figure out the smallest index at which
                                                                     // an object should be inserted so as to maintain order. Uses binary search.
                // Use a comparator function to figure out the smallest index at which
                // an object should be inserted so as to maintain order. Uses binary search.
                _.sortedIndex = function (array, obj, iteratee, context) {
                    iteratee = cb(iteratee, context, 1);
                    var value = iteratee(obj);
                    var low = 0, high = getLength(array);
                    while (low < high) {
                        var mid = Math.floor((low + high) / 2);
                        if (iteratee(array[mid]) < value)
                            low = mid + 1;
                        else
                            high = mid;
                    }
                    return low;
                };    // Generator function to create the indexOf and lastIndexOf functions
                // Generator function to create the indexOf and lastIndexOf functions
                function createIndexFinder(dir, predicateFind, sortedIndex) {
                    return function (array, item, idx) {
                        var i = 0, length = getLength(array);
                        if (typeof idx == 'number') {
                            if (dir > 0) {
                                i = idx >= 0 ? idx : Math.max(idx + length, i);
                            } else {
                                length = idx >= 0 ? Math.min(idx + 1, length) : idx + length + 1;
                            }
                        } else if (sortedIndex && idx && length) {
                            idx = sortedIndex(array, item);
                            return array[idx] === item ? idx : -1;
                        }
                        if (item !== item) {
                            idx = predicateFind(slice.call(array, i, length), _.isNaN);
                            return idx >= 0 ? idx + i : -1;
                        }
                        for (idx = dir > 0 ? i : length - 1; idx >= 0 && idx < length; idx += dir) {
                            if (array[idx] === item)
                                return idx;
                        }
                        return -1;
                    };
                }    // Return the position of the first occurrence of an item in an array,
                     // or -1 if the item is not included in the array.
                     // If the array is large and already in sort order, pass `true`
                     // for **isSorted** to use binary search.
                // Return the position of the first occurrence of an item in an array,
                // or -1 if the item is not included in the array.
                // If the array is large and already in sort order, pass `true`
                // for **isSorted** to use binary search.
                _.indexOf = createIndexFinder(1, _.findIndex, _.sortedIndex);
                _.lastIndexOf = createIndexFinder(-1, _.findLastIndex);    // Generate an integer Array containing an arithmetic progression. A port of
                                                                           // the native Python `range()` function. See
                                                                           // [the Python documentation](http://docs.python.org/library/functions.html#range).
                // Generate an integer Array containing an arithmetic progression. A port of
                // the native Python `range()` function. See
                // [the Python documentation](http://docs.python.org/library/functions.html#range).
                _.range = function (start, stop, step) {
                    if (stop == null) {
                        stop = start || 0;
                        start = 0;
                    }
                    step = step || 1;
                    var length = Math.max(Math.ceil((stop - start) / step), 0);
                    var range = Array(length);
                    for (var idx = 0; idx < length; idx++, start += step) {
                        range[idx] = start;
                    }
                    return range;
                };    // Function (ahem) Functions
                      // ------------------
                      // Determines whether to execute a function as a constructor
                      // or a normal function with the provided arguments
                // Function (ahem) Functions
                // ------------------
                // Determines whether to execute a function as a constructor
                // or a normal function with the provided arguments
                var executeBound = function (sourceFunc, boundFunc, context, callingContext, args) {
                    if (!(callingContext instanceof boundFunc))
                        return sourceFunc.apply(context, args);
                    var self = baseCreate(sourceFunc.prototype);
                    var result = sourceFunc.apply(self, args);
                    if (_.isObject(result))
                        return result;
                    return self;
                };    // Create a function bound to a given object (assigning `this`, and arguments,
                      // optionally). Delegates to **ECMAScript 5**'s native `Function.bind` if
                      // available.
                // Create a function bound to a given object (assigning `this`, and arguments,
                // optionally). Delegates to **ECMAScript 5**'s native `Function.bind` if
                // available.
                _.bind = function (func, context) {
                    if (nativeBind && func.bind === nativeBind)
                        return nativeBind.apply(func, slice.call(arguments, 1));
                    if (!_.isFunction(func))
                        throw new TypeError('Bind must be called on a function');
                    var args = slice.call(arguments, 2);
                    var bound = function () {
                        return executeBound(func, bound, context, this, args.concat(slice.call(arguments)));
                    };
                    return bound;
                };    // Partially apply a function by creating a version that has had some of its
                      // arguments pre-filled, without changing its dynamic `this` context. _ acts
                      // as a placeholder, allowing any combination of arguments to be pre-filled.
                // Partially apply a function by creating a version that has had some of its
                // arguments pre-filled, without changing its dynamic `this` context. _ acts
                // as a placeholder, allowing any combination of arguments to be pre-filled.
                _.partial = function (func) {
                    var boundArgs = slice.call(arguments, 1);
                    var bound = function () {
                        var position = 0, length = boundArgs.length;
                        var args = Array(length);
                        for (var i = 0; i < length; i++) {
                            args[i] = boundArgs[i] === _ ? arguments[position++] : boundArgs[i];
                        }
                        while (position < arguments.length)
                            args.push(arguments[position++]);
                        return executeBound(func, bound, this, this, args);
                    };
                    return bound;
                };    // Bind a number of an object's methods to that object. Remaining arguments
                      // are the method names to be bound. Useful for ensuring that all callbacks
                      // defined on an object belong to it.
                // Bind a number of an object's methods to that object. Remaining arguments
                // are the method names to be bound. Useful for ensuring that all callbacks
                // defined on an object belong to it.
                _.bindAll = function (obj) {
                    var i, length = arguments.length, key;
                    if (length <= 1)
                        throw new Error('bindAll must be passed function names');
                    for (i = 1; i < length; i++) {
                        key = arguments[i];
                        obj[key] = _.bind(obj[key], obj);
                    }
                    return obj;
                };    // Memoize an expensive function by storing its results.
                // Memoize an expensive function by storing its results.
                _.memoize = function (func, hasher) {
                    var memoize = function (key) {
                        var cache = memoize.cache;
                        var address = '' + (hasher ? hasher.apply(this, arguments) : key);
                        if (!_.has(cache, address))
                            cache[address] = func.apply(this, arguments);
                        return cache[address];
                    };
                    memoize.cache = {};
                    return memoize;
                };    // Delays a function for the given number of milliseconds, and then calls
                      // it with the arguments supplied.
                // Delays a function for the given number of milliseconds, and then calls
                // it with the arguments supplied.
                _.delay = function (func, wait) {
                    var args = slice.call(arguments, 2);
                    return setTimeout(function () {
                        return func.apply(null, args);
                    }, wait);
                };    // Defers a function, scheduling it to run after the current call stack has
                      // cleared.
                // Defers a function, scheduling it to run after the current call stack has
                // cleared.
                _.defer = _.partial(_.delay, _, 1);    // Returns a function, that, when invoked, will only be triggered at most once
                                                       // during a given window of time. Normally, the throttled function will run
                                                       // as much as it can, without ever going more than once per `wait` duration;
                                                       // but if you'd like to disable the execution on the leading edge, pass
                                                       // `{leading: false}`. To disable execution on the trailing edge, ditto.
                // Returns a function, that, when invoked, will only be triggered at most once
                // during a given window of time. Normally, the throttled function will run
                // as much as it can, without ever going more than once per `wait` duration;
                // but if you'd like to disable the execution on the leading edge, pass
                // `{leading: false}`. To disable execution on the trailing edge, ditto.
                _.throttle = function (func, wait, options) {
                    var context, args, result;
                    var timeout = null;
                    var previous = 0;
                    if (!options)
                        options = {};
                    var later = function () {
                        previous = options.leading === false ? 0 : _.now();
                        timeout = null;
                        result = func.apply(context, args);
                        if (!timeout)
                            context = args = null;
                    };
                    return function () {
                        var now = _.now();
                        if (!previous && options.leading === false)
                            previous = now;
                        var remaining = wait - (now - previous);
                        context = this;
                        args = arguments;
                        if (remaining <= 0 || remaining > wait) {
                            if (timeout) {
                                clearTimeout(timeout);
                                timeout = null;
                            }
                            previous = now;
                            result = func.apply(context, args);
                            if (!timeout)
                                context = args = null;
                        } else if (!timeout && options.trailing !== false) {
                            timeout = setTimeout(later, remaining);
                        }
                        return result;
                    };
                };    // Returns a function, that, as long as it continues to be invoked, will not
                      // be triggered. The function will be called after it stops being called for
                      // N milliseconds. If `immediate` is passed, trigger the function on the
                      // leading edge, instead of the trailing.
                // Returns a function, that, as long as it continues to be invoked, will not
                // be triggered. The function will be called after it stops being called for
                // N milliseconds. If `immediate` is passed, trigger the function on the
                // leading edge, instead of the trailing.
                _.debounce = function (func, wait, immediate) {
                    var timeout, args, context, timestamp, result;
                    var later = function () {
                        var last = _.now() - timestamp;
                        if (last < wait && last >= 0) {
                            timeout = setTimeout(later, wait - last);
                        } else {
                            timeout = null;
                            if (!immediate) {
                                result = func.apply(context, args);
                                if (!timeout)
                                    context = args = null;
                            }
                        }
                    };
                    return function () {
                        context = this;
                        args = arguments;
                        timestamp = _.now();
                        var callNow = immediate && !timeout;
                        if (!timeout)
                            timeout = setTimeout(later, wait);
                        if (callNow) {
                            result = func.apply(context, args);
                            context = args = null;
                        }
                        return result;
                    };
                };    // Returns the first function passed as an argument to the second,
                      // allowing you to adjust arguments, run code before and after, and
                      // conditionally execute the original function.
                // Returns the first function passed as an argument to the second,
                // allowing you to adjust arguments, run code before and after, and
                // conditionally execute the original function.
                _.wrap = function (func, wrapper) {
                    return _.partial(wrapper, func);
                };    // Returns a negated version of the passed-in predicate.
                // Returns a negated version of the passed-in predicate.
                _.negate = function (predicate) {
                    return function () {
                        return !predicate.apply(this, arguments);
                    };
                };    // Returns a function that is the composition of a list of functions, each
                      // consuming the return value of the function that follows.
                // Returns a function that is the composition of a list of functions, each
                // consuming the return value of the function that follows.
                _.compose = function () {
                    var args = arguments;
                    var start = args.length - 1;
                    return function () {
                        var i = start;
                        var result = args[start].apply(this, arguments);
                        while (i--)
                            result = args[i].call(this, result);
                        return result;
                    };
                };    // Returns a function that will only be executed on and after the Nth call.
                // Returns a function that will only be executed on and after the Nth call.
                _.after = function (times, func) {
                    return function () {
                        if (--times < 1) {
                            return func.apply(this, arguments);
                        }
                    };
                };    // Returns a function that will only be executed up to (but not including) the Nth call.
                // Returns a function that will only be executed up to (but not including) the Nth call.
                _.before = function (times, func) {
                    var memo;
                    return function () {
                        if (--times > 0) {
                            memo = func.apply(this, arguments);
                        }
                        if (times <= 1)
                            func = null;
                        return memo;
                    };
                };    // Returns a function that will be executed at most one time, no matter how
                      // often you call it. Useful for lazy initialization.
                // Returns a function that will be executed at most one time, no matter how
                // often you call it. Useful for lazy initialization.
                _.once = _.partial(_.before, 2);    // Object Functions
                                                    // ----------------
                                                    // Keys in IE < 9 that won't be iterated by `for key in ...` and thus missed.
                // Object Functions
                // ----------------
                // Keys in IE < 9 that won't be iterated by `for key in ...` and thus missed.
                var hasEnumBug = !{ toString: null }.propertyIsEnumerable('toString');
                var nonEnumerableProps = [
                        'valueOf',
                        'isPrototypeOf',
                        'toString',
                        'propertyIsEnumerable',
                        'hasOwnProperty',
                        'toLocaleString'
                    ];
                function collectNonEnumProps(obj, keys) {
                    var nonEnumIdx = nonEnumerableProps.length;
                    var constructor = obj.constructor;
                    var proto = _.isFunction(constructor) && constructor.prototype || ObjProto;    // Constructor is a special case.
                    // Constructor is a special case.
                    var prop = 'constructor';
                    if (_.has(obj, prop) && !_.contains(keys, prop))
                        keys.push(prop);
                    while (nonEnumIdx--) {
                        prop = nonEnumerableProps[nonEnumIdx];
                        if (prop in obj && obj[prop] !== proto[prop] && !_.contains(keys, prop)) {
                            keys.push(prop);
                        }
                    }
                }    // Retrieve the names of an object's own properties.
                     // Delegates to **ECMAScript 5**'s native `Object.keys`
                // Retrieve the names of an object's own properties.
                // Delegates to **ECMAScript 5**'s native `Object.keys`
                _.keys = function (obj) {
                    if (!_.isObject(obj))
                        return [];
                    if (nativeKeys)
                        return nativeKeys(obj);
                    var keys = [];
                    for (var key in obj)
                        if (_.has(obj, key))
                            keys.push(key);    // Ahem, IE < 9.
                    // Ahem, IE < 9.
                    if (hasEnumBug)
                        collectNonEnumProps(obj, keys);
                    return keys;
                };    // Retrieve all the property names of an object.
                // Retrieve all the property names of an object.
                _.allKeys = function (obj) {
                    if (!_.isObject(obj))
                        return [];
                    var keys = [];
                    for (var key in obj)
                        keys.push(key);    // Ahem, IE < 9.
                    // Ahem, IE < 9.
                    if (hasEnumBug)
                        collectNonEnumProps(obj, keys);
                    return keys;
                };    // Retrieve the values of an object's properties.
                // Retrieve the values of an object's properties.
                _.values = function (obj) {
                    var keys = _.keys(obj);
                    var length = keys.length;
                    var values = Array(length);
                    for (var i = 0; i < length; i++) {
                        values[i] = obj[keys[i]];
                    }
                    return values;
                };    // Returns the results of applying the iteratee to each element of the object
                      // In contrast to _.map it returns an object
                // Returns the results of applying the iteratee to each element of the object
                // In contrast to _.map it returns an object
                _.mapObject = function (obj, iteratee, context) {
                    iteratee = cb(iteratee, context);
                    var keys = _.keys(obj), length = keys.length, results = {}, currentKey;
                    for (var index = 0; index < length; index++) {
                        currentKey = keys[index];
                        results[currentKey] = iteratee(obj[currentKey], currentKey, obj);
                    }
                    return results;
                };    // Convert an object into a list of `[key, value]` pairs.
                // Convert an object into a list of `[key, value]` pairs.
                _.pairs = function (obj) {
                    var keys = _.keys(obj);
                    var length = keys.length;
                    var pairs = Array(length);
                    for (var i = 0; i < length; i++) {
                        pairs[i] = [
                            keys[i],
                            obj[keys[i]]
                        ];
                    }
                    return pairs;
                };    // Invert the keys and values of an object. The values must be serializable.
                // Invert the keys and values of an object. The values must be serializable.
                _.invert = function (obj) {
                    var result = {};
                    var keys = _.keys(obj);
                    for (var i = 0, length = keys.length; i < length; i++) {
                        result[obj[keys[i]]] = keys[i];
                    }
                    return result;
                };    // Return a sorted list of the function names available on the object.
                      // Aliased as `methods`
                // Return a sorted list of the function names available on the object.
                // Aliased as `methods`
                _.functions = _.methods = function (obj) {
                    var names = [];
                    for (var key in obj) {
                        if (_.isFunction(obj[key]))
                            names.push(key);
                    }
                    return names.sort();
                };    // Extend a given object with all the properties in passed-in object(s).
                // Extend a given object with all the properties in passed-in object(s).
                _.extend = createAssigner(_.allKeys);    // Assigns a given object with all the own properties in the passed-in object(s)
                                                         // (https://developer.mozilla.org/docs/Web/JavaScript/Reference/Global_Objects/Object/assign)
                // Assigns a given object with all the own properties in the passed-in object(s)
                // (https://developer.mozilla.org/docs/Web/JavaScript/Reference/Global_Objects/Object/assign)
                _.extendOwn = _.assign = createAssigner(_.keys);    // Returns the first key on an object that passes a predicate test
                // Returns the first key on an object that passes a predicate test
                _.findKey = function (obj, predicate, context) {
                    predicate = cb(predicate, context);
                    var keys = _.keys(obj), key;
                    for (var i = 0, length = keys.length; i < length; i++) {
                        key = keys[i];
                        if (predicate(obj[key], key, obj))
                            return key;
                    }
                };    // Return a copy of the object only containing the whitelisted properties.
                // Return a copy of the object only containing the whitelisted properties.
                _.pick = function (object, oiteratee, context) {
                    var result = {}, obj = object, iteratee, keys;
                    if (obj == null)
                        return result;
                    if (_.isFunction(oiteratee)) {
                        keys = _.allKeys(obj);
                        iteratee = optimizeCb(oiteratee, context);
                    } else {
                        keys = flatten(arguments, false, false, 1);
                        iteratee = function (value, key, obj) {
                            return key in obj;
                        };
                        obj = Object(obj);
                    }
                    for (var i = 0, length = keys.length; i < length; i++) {
                        var key = keys[i];
                        var value = obj[key];
                        if (iteratee(value, key, obj))
                            result[key] = value;
                    }
                    return result;
                };    // Return a copy of the object without the blacklisted properties.
                // Return a copy of the object without the blacklisted properties.
                _.omit = function (obj, iteratee, context) {
                    if (_.isFunction(iteratee)) {
                        iteratee = _.negate(iteratee);
                    } else {
                        var keys = _.map(flatten(arguments, false, false, 1), String);
                        iteratee = function (value, key) {
                            return !_.contains(keys, key);
                        };
                    }
                    return _.pick(obj, iteratee, context);
                };    // Fill in a given object with default properties.
                // Fill in a given object with default properties.
                _.defaults = createAssigner(_.allKeys, true);    // Creates an object that inherits from the given prototype object.
                                                                 // If additional properties are provided then they will be added to the
                                                                 // created object.
                // Creates an object that inherits from the given prototype object.
                // If additional properties are provided then they will be added to the
                // created object.
                _.create = function (prototype, props) {
                    var result = baseCreate(prototype);
                    if (props)
                        _.extendOwn(result, props);
                    return result;
                };    // Create a (shallow-cloned) duplicate of an object.
                // Create a (shallow-cloned) duplicate of an object.
                _.clone = function (obj) {
                    if (!_.isObject(obj))
                        return obj;
                    return _.isArray(obj) ? obj.slice() : _.extend({}, obj);
                };    // Invokes interceptor with the obj, and then returns obj.
                      // The primary purpose of this method is to "tap into" a method chain, in
                      // order to perform operations on intermediate results within the chain.
                // Invokes interceptor with the obj, and then returns obj.
                // The primary purpose of this method is to "tap into" a method chain, in
                // order to perform operations on intermediate results within the chain.
                _.tap = function (obj, interceptor) {
                    interceptor(obj);
                    return obj;
                };    // Returns whether an object has a given set of `key:value` pairs.
                // Returns whether an object has a given set of `key:value` pairs.
                _.isMatch = function (object, attrs) {
                    var keys = _.keys(attrs), length = keys.length;
                    if (object == null)
                        return !length;
                    var obj = Object(object);
                    for (var i = 0; i < length; i++) {
                        var key = keys[i];
                        if (attrs[key] !== obj[key] || !(key in obj))
                            return false;
                    }
                    return true;
                };    // Internal recursive comparison function for `isEqual`.
                // Internal recursive comparison function for `isEqual`.
                var eq = function (a, b, aStack, bStack) {
                    // Identical objects are equal. `0 === -0`, but they aren't identical.
                    // See the [Harmony `egal` proposal](http://wiki.ecmascript.org/doku.php?id=harmony:egal).
                    if (a === b)
                        return a !== 0 || 1 / a === 1 / b;    // A strict comparison is necessary because `null == undefined`.
                    // A strict comparison is necessary because `null == undefined`.
                    if (a == null || b == null)
                        return a === b;    // Unwrap any wrapped objects.
                    // Unwrap any wrapped objects.
                    if (a instanceof _)
                        a = a._wrapped;
                    if (b instanceof _)
                        b = b._wrapped;    // Compare `[[Class]]` names.
                    // Compare `[[Class]]` names.
                    var className = toString.call(a);
                    if (className !== toString.call(b))
                        return false;
                    switch (className) {
                    // Strings, numbers, regular expressions, dates, and booleans are compared by value.
                    case '[object RegExp]':    // RegExps are coerced to strings for comparison (Note: '' + /a/i === '/a/i')
                    // RegExps are coerced to strings for comparison (Note: '' + /a/i === '/a/i')
                    case '[object String]':
                        // Primitives and their corresponding object wrappers are equivalent; thus, `"5"` is
                        // equivalent to `new String("5")`.
                        return '' + a === '' + b;
                    case '[object Number]':
                        // `NaN`s are equivalent, but non-reflexive.
                        // Object(NaN) is equivalent to NaN
                        if (+a !== +a)
                            return +b !== +b;    // An `egal` comparison is performed for other numeric values.
                        // An `egal` comparison is performed for other numeric values.
                        return +a === 0 ? 1 / +a === 1 / b : +a === +b;
                    case '[object Date]':
                    case '[object Boolean]':
                        // Coerce dates and booleans to numeric primitive values. Dates are compared by their
                        // millisecond representations. Note that invalid dates with millisecond representations
                        // of `NaN` are not equivalent.
                        return +a === +b;
                    }
                    var areArrays = className === '[object Array]';
                    if (!areArrays) {
                        if (typeof a != 'object' || typeof b != 'object')
                            return false;    // Objects with different constructors are not equivalent, but `Object`s or `Array`s
                                             // from different frames are.
                        // Objects with different constructors are not equivalent, but `Object`s or `Array`s
                        // from different frames are.
                        var aCtor = a.constructor, bCtor = b.constructor;
                        if (aCtor !== bCtor && !(_.isFunction(aCtor) && aCtor instanceof aCtor && _.isFunction(bCtor) && bCtor instanceof bCtor) && ('constructor' in a && 'constructor' in b)) {
                            return false;
                        }
                    }    // Assume equality for cyclic structures. The algorithm for detecting cyclic
                         // structures is adapted from ES 5.1 section 15.12.3, abstract operation `JO`.
                         // Initializing stack of traversed objects.
                         // It's done here since we only need them for objects and arrays comparison.
                    // Assume equality for cyclic structures. The algorithm for detecting cyclic
                    // structures is adapted from ES 5.1 section 15.12.3, abstract operation `JO`.
                    // Initializing stack of traversed objects.
                    // It's done here since we only need them for objects and arrays comparison.
                    aStack = aStack || [];
                    bStack = bStack || [];
                    var length = aStack.length;
                    while (length--) {
                        // Linear search. Performance is inversely proportional to the number of
                        // unique nested structures.
                        if (aStack[length] === a)
                            return bStack[length] === b;
                    }    // Add the first object to the stack of traversed objects.
                    // Add the first object to the stack of traversed objects.
                    aStack.push(a);
                    bStack.push(b);    // Recursively compare objects and arrays.
                    // Recursively compare objects and arrays.
                    if (areArrays) {
                        // Compare array lengths to determine if a deep comparison is necessary.
                        length = a.length;
                        if (length !== b.length)
                            return false;    // Deep compare the contents, ignoring non-numeric properties.
                        // Deep compare the contents, ignoring non-numeric properties.
                        while (length--) {
                            if (!eq(a[length], b[length], aStack, bStack))
                                return false;
                        }
                    } else {
                        // Deep compare objects.
                        var keys = _.keys(a), key;
                        length = keys.length;    // Ensure that both objects contain the same number of properties before comparing deep equality.
                        // Ensure that both objects contain the same number of properties before comparing deep equality.
                        if (_.keys(b).length !== length)
                            return false;
                        while (length--) {
                            // Deep compare each member
                            key = keys[length];
                            if (!(_.has(b, key) && eq(a[key], b[key], aStack, bStack)))
                                return false;
                        }
                    }    // Remove the first object from the stack of traversed objects.
                    // Remove the first object from the stack of traversed objects.
                    aStack.pop();
                    bStack.pop();
                    return true;
                };    // Perform a deep comparison to check if two objects are equal.
                // Perform a deep comparison to check if two objects are equal.
                _.isEqual = function (a, b) {
                    return eq(a, b);
                };    // Is a given array, string, or object empty?
                      // An "empty" object has no enumerable own-properties.
                // Is a given array, string, or object empty?
                // An "empty" object has no enumerable own-properties.
                _.isEmpty = function (obj) {
                    if (obj == null)
                        return true;
                    if (isArrayLike(obj) && (_.isArray(obj) || _.isString(obj) || _.isArguments(obj)))
                        return obj.length === 0;
                    return _.keys(obj).length === 0;
                };    // Is a given value a DOM element?
                // Is a given value a DOM element?
                _.isElement = function (obj) {
                    return !!(obj && obj.nodeType === 1);
                };    // Is a given value an array?
                      // Delegates to ECMA5's native Array.isArray
                // Is a given value an array?
                // Delegates to ECMA5's native Array.isArray
                _.isArray = nativeIsArray || function (obj) {
                    return toString.call(obj) === '[object Array]';
                };    // Is a given variable an object?
                // Is a given variable an object?
                _.isObject = function (obj) {
                    var type = typeof obj;
                    return type === 'function' || type === 'object' && !!obj;
                };    // Add some isType methods: isArguments, isFunction, isString, isNumber, isDate, isRegExp, isError.
                // Add some isType methods: isArguments, isFunction, isString, isNumber, isDate, isRegExp, isError.
                _.each([
                    'Arguments',
                    'Function',
                    'String',
                    'Number',
                    'Date',
                    'RegExp',
                    'Error'
                ], function (name) {
                    _['is' + name] = function (obj) {
                        return toString.call(obj) === '[object ' + name + ']';
                    };
                });    // Define a fallback version of the method in browsers (ahem, IE < 9), where
                       // there isn't any inspectable "Arguments" type.
                // Define a fallback version of the method in browsers (ahem, IE < 9), where
                // there isn't any inspectable "Arguments" type.
                if (!_.isArguments(arguments)) {
                    _.isArguments = function (obj) {
                        return _.has(obj, 'callee');
                    };
                }    // Optimize `isFunction` if appropriate. Work around some typeof bugs in old v8,
                     // IE 11 (#1621), and in Safari 8 (#1929).
                // Optimize `isFunction` if appropriate. Work around some typeof bugs in old v8,
                // IE 11 (#1621), and in Safari 8 (#1929).
                if (typeof /./ != 'function' && typeof Int8Array != 'object') {
                    _.isFunction = function (obj) {
                        return typeof obj == 'function' || false;
                    };
                }    // Is a given object a finite number?
                // Is a given object a finite number?
                _.isFinite = function (obj) {
                    return isFinite(obj) && !isNaN(parseFloat(obj));
                };    // Is the given value `NaN`? (NaN is the only number which does not equal itself).
                // Is the given value `NaN`? (NaN is the only number which does not equal itself).
                _.isNaN = function (obj) {
                    return _.isNumber(obj) && obj !== +obj;
                };    // Is a given value a boolean?
                // Is a given value a boolean?
                _.isBoolean = function (obj) {
                    return obj === true || obj === false || toString.call(obj) === '[object Boolean]';
                };    // Is a given value equal to null?
                // Is a given value equal to null?
                _.isNull = function (obj) {
                    return obj === null;
                };    // Is a given variable undefined?
                // Is a given variable undefined?
                _.isUndefined = function (obj) {
                    return obj === void 0;
                };    // Shortcut function for checking if an object has a given property directly
                      // on itself (in other words, not on a prototype).
                // Shortcut function for checking if an object has a given property directly
                // on itself (in other words, not on a prototype).
                _.has = function (obj, key) {
                    return obj != null && hasOwnProperty.call(obj, key);
                };    // Utility Functions
                      // -----------------
                      // Run Underscore.js in *noConflict* mode, returning the `_` variable to its
                      // previous owner. Returns a reference to the Underscore object.
                // Utility Functions
                // -----------------
                // Run Underscore.js in *noConflict* mode, returning the `_` variable to its
                // previous owner. Returns a reference to the Underscore object.
                _.noConflict = function () {
                    root._ = previousUnderscore;
                    return this;
                };    // Keep the identity function around for default iteratees.
                // Keep the identity function around for default iteratees.
                _.identity = function (value) {
                    return value;
                };    // Predicate-generating functions. Often useful outside of Underscore.
                // Predicate-generating functions. Often useful outside of Underscore.
                _.constant = function (value) {
                    return function () {
                        return value;
                    };
                };
                _.noop = function () {
                };
                _.property = property;    // Generates a function for a given object that returns a given property.
                // Generates a function for a given object that returns a given property.
                _.propertyOf = function (obj) {
                    return obj == null ? function () {
                    } : function (key) {
                        return obj[key];
                    };
                };    // Returns a predicate for checking whether an object has a given set of
                      // `key:value` pairs.
                // Returns a predicate for checking whether an object has a given set of
                // `key:value` pairs.
                _.matcher = _.matches = function (attrs) {
                    attrs = _.extendOwn({}, attrs);
                    return function (obj) {
                        return _.isMatch(obj, attrs);
                    };
                };    // Run a function **n** times.
                // Run a function **n** times.
                _.times = function (n, iteratee, context) {
                    var accum = Array(Math.max(0, n));
                    iteratee = optimizeCb(iteratee, context, 1);
                    for (var i = 0; i < n; i++)
                        accum[i] = iteratee(i);
                    return accum;
                };    // Return a random integer between min and max (inclusive).
                // Return a random integer between min and max (inclusive).
                _.random = function (min, max) {
                    if (max == null) {
                        max = min;
                        min = 0;
                    }
                    return min + Math.floor(Math.random() * (max - min + 1));
                };    // A (possibly faster) way to get the current timestamp as an integer.
                // A (possibly faster) way to get the current timestamp as an integer.
                _.now = Date.now || function () {
                    return new Date().getTime();
                };    // List of HTML entities for escaping.
                // List of HTML entities for escaping.
                var escapeMap = {
                        '&': '&amp;',
                        '<': '&lt;',
                        '>': '&gt;',
                        '"': '&quot;',
                        '\'': '&#x27;',
                        '`': '&#x60;'
                    };
                var unescapeMap = _.invert(escapeMap);    // Functions for escaping and unescaping strings to/from HTML interpolation.
                // Functions for escaping and unescaping strings to/from HTML interpolation.
                var createEscaper = function (map) {
                    var escaper = function (match) {
                        return map[match];
                    };    // Regexes for identifying a key that needs to be escaped
                    // Regexes for identifying a key that needs to be escaped
                    var source = '(?:' + _.keys(map).join('|') + ')';
                    var testRegexp = RegExp(source);
                    var replaceRegexp = RegExp(source, 'g');
                    return function (string) {
                        string = string == null ? '' : '' + string;
                        return testRegexp.test(string) ? string.replace(replaceRegexp, escaper) : string;
                    };
                };
                _.escape = createEscaper(escapeMap);
                _.unescape = createEscaper(unescapeMap);    // If the value of the named `property` is a function then invoke it with the
                                                            // `object` as context; otherwise, return it.
                // If the value of the named `property` is a function then invoke it with the
                // `object` as context; otherwise, return it.
                _.result = function (object, property, fallback) {
                    var value = object == null ? void 0 : object[property];
                    if (value === void 0) {
                        value = fallback;
                    }
                    return _.isFunction(value) ? value.call(object) : value;
                };    // Generate a unique integer id (unique within the entire client session).
                      // Useful for temporary DOM ids.
                // Generate a unique integer id (unique within the entire client session).
                // Useful for temporary DOM ids.
                var idCounter = 0;
                _.uniqueId = function (prefix) {
                    var id = ++idCounter + '';
                    return prefix ? prefix + id : id;
                };    // By default, Underscore uses ERB-style template delimiters, change the
                      // following template settings to use alternative delimiters.
                // By default, Underscore uses ERB-style template delimiters, change the
                // following template settings to use alternative delimiters.
                _.templateSettings = {
                    evaluate: /<%([\s\S]+?)%>/g,
                    interpolate: /<%=([\s\S]+?)%>/g,
                    escape: /<%-([\s\S]+?)%>/g
                };    // When customizing `templateSettings`, if you don't want to define an
                      // interpolation, evaluation or escaping regex, we need one that is
                      // guaranteed not to match.
                // When customizing `templateSettings`, if you don't want to define an
                // interpolation, evaluation or escaping regex, we need one that is
                // guaranteed not to match.
                var noMatch = /(.)^/;    // Certain characters need to be escaped so that they can be put into a
                                         // string literal.
                // Certain characters need to be escaped so that they can be put into a
                // string literal.
                var escapes = {
                        '\'': '\'',
                        '\\': '\\',
                        '\r': 'r',
                        '\n': 'n',
                        '\u2028': 'u2028',
                        '\u2029': 'u2029'
                    };
                var escaper = /\\|'|\r|\n|\u2028|\u2029/g;
                var escapeChar = function (match) {
                    return '\\' + escapes[match];
                };    // JavaScript micro-templating, similar to John Resig's implementation.
                      // Underscore templating handles arbitrary delimiters, preserves whitespace,
                      // and correctly escapes quotes within interpolated code.
                      // NB: `oldSettings` only exists for backwards compatibility.
                // JavaScript micro-templating, similar to John Resig's implementation.
                // Underscore templating handles arbitrary delimiters, preserves whitespace,
                // and correctly escapes quotes within interpolated code.
                // NB: `oldSettings` only exists for backwards compatibility.
                _.template = function (text, settings, oldSettings) {
                    if (!settings && oldSettings)
                        settings = oldSettings;
                    settings = _.defaults({}, settings, _.templateSettings);    // Combine delimiters into one regular expression via alternation.
                    // Combine delimiters into one regular expression via alternation.
                    var matcher = RegExp([
                            (settings.escape || noMatch).source,
                            (settings.interpolate || noMatch).source,
                            (settings.evaluate || noMatch).source
                        ].join('|') + '|$', 'g');    // Compile the template source, escaping string literals appropriately.
                    // Compile the template source, escaping string literals appropriately.
                    var index = 0;
                    var source = '__p+=\'';
                    text.replace(matcher, function (match, escape, interpolate, evaluate, offset) {
                        source += text.slice(index, offset).replace(escaper, escapeChar);
                        index = offset + match.length;
                        if (escape) {
                            source += '\'+\n((__t=(' + escape + '))==null?\'\':_.escape(__t))+\n\'';
                        } else if (interpolate) {
                            source += '\'+\n((__t=(' + interpolate + '))==null?\'\':__t)+\n\'';
                        } else if (evaluate) {
                            source += '\';\n' + evaluate + '\n__p+=\'';
                        }    // Adobe VMs need the match returned to produce the correct offest.
                        // Adobe VMs need the match returned to produce the correct offest.
                        return match;
                    });
                    source += '\';\n';    // If a variable is not specified, place data values in local scope.
                    // If a variable is not specified, place data values in local scope.
                    if (!settings.variable)
                        source = 'with(obj||{}){\n' + source + '}\n';
                    source = 'var __t,__p=\'\',__j=Array.prototype.join,' + 'print=function(){__p+=__j.call(arguments,\'\');};\n' + source + 'return __p;\n';
                    try {
                        var render = new Function(settings.variable || 'obj', '_', source);
                    } catch (e) {
                        e.source = source;
                        throw e;
                    }
                    var template = function (data) {
                        return render.call(this, data, _);
                    };    // Provide the compiled source as a convenience for precompilation.
                    // Provide the compiled source as a convenience for precompilation.
                    var argument = settings.variable || 'obj';
                    template.source = 'function(' + argument + '){\n' + source + '}';
                    return template;
                };    // Add a "chain" function. Start chaining a wrapped Underscore object.
                // Add a "chain" function. Start chaining a wrapped Underscore object.
                _.chain = function (obj) {
                    var instance = _(obj);
                    instance._chain = true;
                    return instance;
                };    // OOP
                      // ---------------
                      // If Underscore is called as a function, it returns a wrapped object that
                      // can be used OO-style. This wrapper holds altered versions of all the
                      // underscore functions. Wrapped objects may be chained.
                      // Helper function to continue chaining intermediate results.
                // OOP
                // ---------------
                // If Underscore is called as a function, it returns a wrapped object that
                // can be used OO-style. This wrapper holds altered versions of all the
                // underscore functions. Wrapped objects may be chained.
                // Helper function to continue chaining intermediate results.
                var result = function (instance, obj) {
                    return instance._chain ? _(obj).chain() : obj;
                };    // Add your own custom functions to the Underscore object.
                // Add your own custom functions to the Underscore object.
                _.mixin = function (obj) {
                    _.each(_.functions(obj), function (name) {
                        var func = _[name] = obj[name];
                        _.prototype[name] = function () {
                            var args = [this._wrapped];
                            push.apply(args, arguments);
                            return result(this, func.apply(_, args));
                        };
                    });
                };    // Add all of the Underscore functions to the wrapper object.
                // Add all of the Underscore functions to the wrapper object.
                _.mixin(_);    // Add all mutator Array functions to the wrapper.
                // Add all mutator Array functions to the wrapper.
                _.each([
                    'pop',
                    'push',
                    'reverse',
                    'shift',
                    'sort',
                    'splice',
                    'unshift'
                ], function (name) {
                    var method = ArrayProto[name];
                    _.prototype[name] = function () {
                        var obj = this._wrapped;
                        method.apply(obj, arguments);
                        if ((name === 'shift' || name === 'splice') && obj.length === 0)
                            delete obj[0];
                        return result(this, obj);
                    };
                });    // Add all accessor Array functions to the wrapper.
                // Add all accessor Array functions to the wrapper.
                _.each([
                    'concat',
                    'join',
                    'slice'
                ], function (name) {
                    var method = ArrayProto[name];
                    _.prototype[name] = function () {
                        return result(this, method.apply(this._wrapped, arguments));
                    };
                });    // Extracts the result from a wrapped and chained object.
                // Extracts the result from a wrapped and chained object.
                _.prototype.value = function () {
                    return this._wrapped;
                };    // Provide unwrapping proxy for some methods used in engine operations
                      // such as arithmetic and JSON stringification.
                // Provide unwrapping proxy for some methods used in engine operations
                // such as arithmetic and JSON stringification.
                _.prototype.valueOf = _.prototype.toJSON = _.prototype.value;
                _.prototype.toString = function () {
                    return '' + this._wrapped;
                };    // AMD registration happens at the end for compatibility with AMD loaders
                      // that may not enforce next-turn semantics on modules. Even though general
                      // practice for AMD registration is to be anonymous, underscore registers
                      // as a named module because, like jQuery, it is a base library that is
                      // popular enough to be bundled in a third party lib, but not be part of
                      // an AMD load request. Those cases could generate an error when an
                      // anonymous define() is called outside of a loader request.
                // AMD registration happens at the end for compatibility with AMD loaders
                // that may not enforce next-turn semantics on modules. Even though general
                // practice for AMD registration is to be anonymous, underscore registers
                // as a named module because, like jQuery, it is a base library that is
                // popular enough to be bundled in a third party lib, but not be part of
                // an AMD load request. Those cases could generate an error when an
                // anonymous define() is called outside of a loader request.
                if (typeof define === 'function' && define.amd) {
                    define('underscore', [], function () {
                        return _;
                    });
                }
            }.call(this));
        }
    ];
    return _require(0);
}));