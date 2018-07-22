(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
	typeof define === 'function' && define.amd ? define(['exports'], factory) :
	(factory((global.geometryExtrude = {})));
}(this, (function (exports) { 'use strict';

var earcut_1 = earcut;
var default_1 = earcut;

function earcut(data, holeIndices, dim) {

    dim = dim || 2;

    var hasHoles = holeIndices && holeIndices.length,
        outerLen = hasHoles ? holeIndices[0] * dim : data.length,
        outerNode = linkedList(data, 0, outerLen, dim, true),
        triangles = [];

    if (!outerNode) return triangles;

    var minX, minY, maxX, maxY, x, y, invSize;

    if (hasHoles) outerNode = eliminateHoles(data, holeIndices, outerNode, dim);

    // if the shape is not too simple, we'll use z-order curve hash later; calculate polygon bbox
    if (data.length > 80 * dim) {
        minX = maxX = data[0];
        minY = maxY = data[1];

        for (var i = dim; i < outerLen; i += dim) {
            x = data[i];
            y = data[i + 1];
            if (x < minX) minX = x;
            if (y < minY) minY = y;
            if (x > maxX) maxX = x;
            if (y > maxY) maxY = y;
        }

        // minX, minY and invSize are later used to transform coords into integers for z-order calculation
        invSize = Math.max(maxX - minX, maxY - minY);
        invSize = invSize !== 0 ? 1 / invSize : 0;
    }

    earcutLinked(outerNode, triangles, dim, minX, minY, invSize);

    return triangles;
}

// create a circular doubly linked list from polygon points in the specified winding order
function linkedList(data, start, end, dim, clockwise) {
    var i, last;

    if (clockwise === (signedArea(data, start, end, dim) > 0)) {
        for (i = start; i < end; i += dim) last = insertNode(i, data[i], data[i + 1], last);
    } else {
        for (i = end - dim; i >= start; i -= dim) last = insertNode(i, data[i], data[i + 1], last);
    }

    if (last && equals(last, last.next)) {
        removeNode(last);
        last = last.next;
    }

    return last;
}

// eliminate colinear or duplicate points
function filterPoints(start, end) {
    if (!start) return start;
    if (!end) end = start;

    var p = start,
        again;
    do {
        again = false;

        if (!p.steiner && (equals(p, p.next) || area(p.prev, p, p.next) === 0)) {
            removeNode(p);
            p = end = p.prev;
            if (p === p.next) break;
            again = true;

        } else {
            p = p.next;
        }
    } while (again || p !== end);

    return end;
}

// main ear slicing loop which triangulates a polygon (given as a linked list)
function earcutLinked(ear, triangles, dim, minX, minY, invSize, pass) {
    if (!ear) return;

    // interlink polygon nodes in z-order
    if (!pass && invSize) indexCurve(ear, minX, minY, invSize);

    var stop = ear,
        prev, next;

    // iterate through ears, slicing them one by one
    while (ear.prev !== ear.next) {
        prev = ear.prev;
        next = ear.next;

        if (invSize ? isEarHashed(ear, minX, minY, invSize) : isEar(ear)) {
            // cut off the triangle
            triangles.push(prev.i / dim);
            triangles.push(ear.i / dim);
            triangles.push(next.i / dim);

            removeNode(ear);

            // skipping the next vertice leads to less sliver triangles
            ear = next.next;
            stop = next.next;

            continue;
        }

        ear = next;

        // if we looped through the whole remaining polygon and can't find any more ears
        if (ear === stop) {
            // try filtering points and slicing again
            if (!pass) {
                earcutLinked(filterPoints(ear), triangles, dim, minX, minY, invSize, 1);

            // if this didn't work, try curing all small self-intersections locally
            } else if (pass === 1) {
                ear = cureLocalIntersections(ear, triangles, dim);
                earcutLinked(ear, triangles, dim, minX, minY, invSize, 2);

            // as a last resort, try splitting the remaining polygon into two
            } else if (pass === 2) {
                splitEarcut(ear, triangles, dim, minX, minY, invSize);
            }

            break;
        }
    }
}

// check whether a polygon node forms a valid ear with adjacent nodes
function isEar(ear) {
    var a = ear.prev,
        b = ear,
        c = ear.next;

    if (area(a, b, c) >= 0) return false; // reflex, can't be an ear

    // now make sure we don't have other points inside the potential ear
    var p = ear.next.next;

    while (p !== ear.prev) {
        if (pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
            area(p.prev, p, p.next) >= 0) return false;
        p = p.next;
    }

    return true;
}

function isEarHashed(ear, minX, minY, invSize) {
    var a = ear.prev,
        b = ear,
        c = ear.next;

    if (area(a, b, c) >= 0) return false; // reflex, can't be an ear

    // triangle bbox; min & max are calculated like this for speed
    var minTX = a.x < b.x ? (a.x < c.x ? a.x : c.x) : (b.x < c.x ? b.x : c.x),
        minTY = a.y < b.y ? (a.y < c.y ? a.y : c.y) : (b.y < c.y ? b.y : c.y),
        maxTX = a.x > b.x ? (a.x > c.x ? a.x : c.x) : (b.x > c.x ? b.x : c.x),
        maxTY = a.y > b.y ? (a.y > c.y ? a.y : c.y) : (b.y > c.y ? b.y : c.y);

    // z-order range for the current triangle bbox;
    var minZ = zOrder(minTX, minTY, minX, minY, invSize),
        maxZ = zOrder(maxTX, maxTY, minX, minY, invSize);

    var p = ear.prevZ,
        n = ear.nextZ;

    // look for points inside the triangle in both directions
    while (p && p.z >= minZ && n && n.z <= maxZ) {
        if (p !== ear.prev && p !== ear.next &&
            pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
            area(p.prev, p, p.next) >= 0) return false;
        p = p.prevZ;

        if (n !== ear.prev && n !== ear.next &&
            pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) &&
            area(n.prev, n, n.next) >= 0) return false;
        n = n.nextZ;
    }

    // look for remaining points in decreasing z-order
    while (p && p.z >= minZ) {
        if (p !== ear.prev && p !== ear.next &&
            pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
            area(p.prev, p, p.next) >= 0) return false;
        p = p.prevZ;
    }

    // look for remaining points in increasing z-order
    while (n && n.z <= maxZ) {
        if (n !== ear.prev && n !== ear.next &&
            pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) &&
            area(n.prev, n, n.next) >= 0) return false;
        n = n.nextZ;
    }

    return true;
}

// go through all polygon nodes and cure small local self-intersections
function cureLocalIntersections(start, triangles, dim) {
    var p = start;
    do {
        var a = p.prev,
            b = p.next.next;

        if (!equals(a, b) && intersects(a, p, p.next, b) && locallyInside(a, b) && locallyInside(b, a)) {

            triangles.push(a.i / dim);
            triangles.push(p.i / dim);
            triangles.push(b.i / dim);

            // remove two nodes involved
            removeNode(p);
            removeNode(p.next);

            p = start = b;
        }
        p = p.next;
    } while (p !== start);

    return p;
}

// try splitting polygon into two and triangulate them independently
function splitEarcut(start, triangles, dim, minX, minY, invSize) {
    // look for a valid diagonal that divides the polygon into two
    var a = start;
    do {
        var b = a.next.next;
        while (b !== a.prev) {
            if (a.i !== b.i && isValidDiagonal(a, b)) {
                // split the polygon in two by the diagonal
                var c = splitPolygon(a, b);

                // filter colinear points around the cuts
                a = filterPoints(a, a.next);
                c = filterPoints(c, c.next);

                // run earcut on each half
                earcutLinked(a, triangles, dim, minX, minY, invSize);
                earcutLinked(c, triangles, dim, minX, minY, invSize);
                return;
            }
            b = b.next;
        }
        a = a.next;
    } while (a !== start);
}

// link every hole into the outer loop, producing a single-ring polygon without holes
function eliminateHoles(data, holeIndices, outerNode, dim) {
    var queue = [],
        i, len, start, end, list;

    for (i = 0, len = holeIndices.length; i < len; i++) {
        start = holeIndices[i] * dim;
        end = i < len - 1 ? holeIndices[i + 1] * dim : data.length;
        list = linkedList(data, start, end, dim, false);
        if (list === list.next) list.steiner = true;
        queue.push(getLeftmost(list));
    }

    queue.sort(compareX);

    // process holes from left to right
    for (i = 0; i < queue.length; i++) {
        eliminateHole(queue[i], outerNode);
        outerNode = filterPoints(outerNode, outerNode.next);
    }

    return outerNode;
}

function compareX(a, b) {
    return a.x - b.x;
}

// find a bridge between vertices that connects hole with an outer ring and and link it
function eliminateHole(hole, outerNode) {
    outerNode = findHoleBridge(hole, outerNode);
    if (outerNode) {
        var b = splitPolygon(outerNode, hole);
        filterPoints(b, b.next);
    }
}

// David Eberly's algorithm for finding a bridge between hole and outer polygon
function findHoleBridge(hole, outerNode) {
    var p = outerNode,
        hx = hole.x,
        hy = hole.y,
        qx = -Infinity,
        m;

    // find a segment intersected by a ray from the hole's leftmost point to the left;
    // segment's endpoint with lesser x will be potential connection point
    do {
        if (hy <= p.y && hy >= p.next.y && p.next.y !== p.y) {
            var x = p.x + (hy - p.y) * (p.next.x - p.x) / (p.next.y - p.y);
            if (x <= hx && x > qx) {
                qx = x;
                if (x === hx) {
                    if (hy === p.y) return p;
                    if (hy === p.next.y) return p.next;
                }
                m = p.x < p.next.x ? p : p.next;
            }
        }
        p = p.next;
    } while (p !== outerNode);

    if (!m) return null;

    if (hx === qx) return m.prev; // hole touches outer segment; pick lower endpoint

    // look for points inside the triangle of hole point, segment intersection and endpoint;
    // if there are no points found, we have a valid connection;
    // otherwise choose the point of the minimum angle with the ray as connection point

    var stop = m,
        mx = m.x,
        my = m.y,
        tanMin = Infinity,
        tan;

    p = m.next;

    while (p !== stop) {
        if (hx >= p.x && p.x >= mx && hx !== p.x &&
                pointInTriangle(hy < my ? hx : qx, hy, mx, my, hy < my ? qx : hx, hy, p.x, p.y)) {

            tan = Math.abs(hy - p.y) / (hx - p.x); // tangential

            if ((tan < tanMin || (tan === tanMin && p.x > m.x)) && locallyInside(p, hole)) {
                m = p;
                tanMin = tan;
            }
        }

        p = p.next;
    }

    return m;
}

// interlink polygon nodes in z-order
function indexCurve(start, minX, minY, invSize) {
    var p = start;
    do {
        if (p.z === null) p.z = zOrder(p.x, p.y, minX, minY, invSize);
        p.prevZ = p.prev;
        p.nextZ = p.next;
        p = p.next;
    } while (p !== start);

    p.prevZ.nextZ = null;
    p.prevZ = null;

    sortLinked(p);
}

// Simon Tatham's linked list merge sort algorithm
// http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
function sortLinked(list) {
    var i, p, q, e, tail, numMerges, pSize, qSize,
        inSize = 1;

    do {
        p = list;
        list = null;
        tail = null;
        numMerges = 0;

        while (p) {
            numMerges++;
            q = p;
            pSize = 0;
            for (i = 0; i < inSize; i++) {
                pSize++;
                q = q.nextZ;
                if (!q) break;
            }
            qSize = inSize;

            while (pSize > 0 || (qSize > 0 && q)) {

                if (pSize !== 0 && (qSize === 0 || !q || p.z <= q.z)) {
                    e = p;
                    p = p.nextZ;
                    pSize--;
                } else {
                    e = q;
                    q = q.nextZ;
                    qSize--;
                }

                if (tail) tail.nextZ = e;
                else list = e;

                e.prevZ = tail;
                tail = e;
            }

            p = q;
        }

        tail.nextZ = null;
        inSize *= 2;

    } while (numMerges > 1);

    return list;
}

// z-order of a point given coords and inverse of the longer side of data bbox
function zOrder(x, y, minX, minY, invSize) {
    // coords are transformed into non-negative 15-bit integer range
    x = 32767 * (x - minX) * invSize;
    y = 32767 * (y - minY) * invSize;

    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;

    y = (y | (y << 8)) & 0x00FF00FF;
    y = (y | (y << 4)) & 0x0F0F0F0F;
    y = (y | (y << 2)) & 0x33333333;
    y = (y | (y << 1)) & 0x55555555;

    return x | (y << 1);
}

// find the leftmost node of a polygon ring
function getLeftmost(start) {
    var p = start,
        leftmost = start;
    do {
        if (p.x < leftmost.x) leftmost = p;
        p = p.next;
    } while (p !== start);

    return leftmost;
}

// check if a point lies within a convex triangle
function pointInTriangle(ax, ay, bx, by, cx, cy, px, py) {
    return (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0 &&
           (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0 &&
           (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0;
}

// check if a diagonal between two polygon nodes is valid (lies in polygon interior)
function isValidDiagonal(a, b) {
    return a.next.i !== b.i && a.prev.i !== b.i && !intersectsPolygon(a, b) &&
           locallyInside(a, b) && locallyInside(b, a) && middleInside(a, b);
}

// signed area of a triangle
function area(p, q, r) {
    return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
}

// check if two points are equal
function equals(p1, p2) {
    return p1.x === p2.x && p1.y === p2.y;
}

// check if two segments intersect
function intersects(p1, q1, p2, q2) {
    if ((equals(p1, q1) && equals(p2, q2)) ||
        (equals(p1, q2) && equals(p2, q1))) return true;
    return area(p1, q1, p2) > 0 !== area(p1, q1, q2) > 0 &&
           area(p2, q2, p1) > 0 !== area(p2, q2, q1) > 0;
}

// check if a polygon diagonal intersects any polygon segments
function intersectsPolygon(a, b) {
    var p = a;
    do {
        if (p.i !== a.i && p.next.i !== a.i && p.i !== b.i && p.next.i !== b.i &&
                intersects(p, p.next, a, b)) return true;
        p = p.next;
    } while (p !== a);

    return false;
}

// check if a polygon diagonal is locally inside the polygon
function locallyInside(a, b) {
    return area(a.prev, a, a.next) < 0 ?
        area(a, b, a.next) >= 0 && area(a, a.prev, b) >= 0 :
        area(a, b, a.prev) < 0 || area(a, a.next, b) < 0;
}

// check if the middle point of a polygon diagonal is inside the polygon
function middleInside(a, b) {
    var p = a,
        inside = false,
        px = (a.x + b.x) / 2,
        py = (a.y + b.y) / 2;
    do {
        if (((p.y > py) !== (p.next.y > py)) && p.next.y !== p.y &&
                (px < (p.next.x - p.x) * (py - p.y) / (p.next.y - p.y) + p.x))
            inside = !inside;
        p = p.next;
    } while (p !== a);

    return inside;
}

// link two polygon vertices with a bridge; if the vertices belong to the same ring, it splits polygon into two;
// if one belongs to the outer ring and another to a hole, it merges it into a single ring
function splitPolygon(a, b) {
    var a2 = new Node(a.i, a.x, a.y),
        b2 = new Node(b.i, b.x, b.y),
        an = a.next,
        bp = b.prev;

    a.next = b;
    b.prev = a;

    a2.next = an;
    an.prev = a2;

    b2.next = a2;
    a2.prev = b2;

    bp.next = b2;
    b2.prev = bp;

    return b2;
}

// create a node and optionally link it with previous one (in a circular doubly linked list)
function insertNode(i, x, y, last) {
    var p = new Node(i, x, y);

    if (!last) {
        p.prev = p;
        p.next = p;

    } else {
        p.next = last.next;
        p.prev = last;
        last.next.prev = p;
        last.next = p;
    }
    return p;
}

function removeNode(p) {
    p.next.prev = p.prev;
    p.prev.next = p.next;

    if (p.prevZ) p.prevZ.nextZ = p.nextZ;
    if (p.nextZ) p.nextZ.prevZ = p.prevZ;
}

function Node(i, x, y) {
    // vertice index in coordinates array
    this.i = i;

    // vertex coordinates
    this.x = x;
    this.y = y;

    // previous and next vertice nodes in a polygon ring
    this.prev = null;
    this.next = null;

    // z-order curve value
    this.z = null;

    // previous and next nodes in z-order
    this.prevZ = null;
    this.nextZ = null;

    // indicates whether this is a steiner point
    this.steiner = false;
}

// return a percentage difference between the polygon area and its triangulation area;
// used to verify correctness of triangulation
earcut.deviation = function (data, holeIndices, dim, triangles) {
    var hasHoles = holeIndices && holeIndices.length;
    var outerLen = hasHoles ? holeIndices[0] * dim : data.length;

    var polygonArea = Math.abs(signedArea(data, 0, outerLen, dim));
    if (hasHoles) {
        for (var i = 0, len = holeIndices.length; i < len; i++) {
            var start = holeIndices[i] * dim;
            var end = i < len - 1 ? holeIndices[i + 1] * dim : data.length;
            polygonArea -= Math.abs(signedArea(data, start, end, dim));
        }
    }

    var trianglesArea = 0;
    for (i = 0; i < triangles.length; i += 3) {
        var a = triangles[i] * dim;
        var b = triangles[i + 1] * dim;
        var c = triangles[i + 2] * dim;
        trianglesArea += Math.abs(
            (data[a] - data[c]) * (data[b + 1] - data[a + 1]) -
            (data[a] - data[b]) * (data[c + 1] - data[a + 1]));
    }

    return polygonArea === 0 && trianglesArea === 0 ? 0 :
        Math.abs((trianglesArea - polygonArea) / polygonArea);
};

function signedArea(data, start, end, dim) {
    var sum = 0;
    for (var i = start, j = end - dim; i < end; i += dim) {
        sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
        j = i;
    }
    return sum;
}

// turn a polygon in a multi-dimensional array form (e.g. as in GeoJSON) into a form Earcut accepts
earcut.flatten = function (data) {
    var dim = data[0][0].length,
        result = {vertices: [], holes: [], dimensions: dim},
        holeIndex = 0;

    for (var i = 0; i < data.length; i++) {
        for (var j = 0; j < data[i].length; j++) {
            for (var d = 0; d < dim; d++) result.vertices.push(data[i][j][d]);
        }
        if (i > 0) {
            holeIndex += data[i - 1].length;
            result.holes.push(holeIndex);
        }
    }
    return result;
};

earcut_1.default = default_1;

/*
 (c) 2017, Vladimir Agafonkin
 Simplify.js, a high-performance JS polyline simplification library
 mourner.github.io/simplify-js
*/

// to suit your point format, run search/replace for '.x' and '.y';
// for 3D version, see 3d branch (configurability would draw significant performance overhead)

// square distance between 2 points
function getSqDist(p1, p2) {

    var dx = p1[0] - p2[0],
        dy = p1[1] - p2[1];

    return dx * dx + dy * dy;
}

// square distance from a point to a segment
function getSqSegDist(p, p1, p2) {

    var x = p1[0],
        y = p1[1],
        dx = p2[0] - x,
        dy = p2[1] - y;

    if (dx !== 0 || dy !== 0) {

        var t = ((p[0] - x) * dx + (p[1] - y) * dy) / (dx * dx + dy * dy);

        if (t > 1) {
            x = p2[0];
            y = p2[1];

        } else if (t > 0) {
            x += dx * t;
            y += dy * t;
        }
    }

    dx = p[0] - x;
    dy = p[1] - y;

    return dx * dx + dy * dy;
}
// rest of the code doesn't care about point format

// basic distance-based simplification
function simplifyRadialDist(points, sqTolerance) {

    var prevPoint = points[0],
        newPoints = [prevPoint],
        point;

    for (var i = 1, len = points.length; i < len; i++) {
        point = points[i];

        if (getSqDist(point, prevPoint) > sqTolerance) {
            newPoints.push(point);
            prevPoint = point;
        }
    }

    if (prevPoint !== point) newPoints.push(point);

    return newPoints;
}

function simplifyDPStep(points, first, last, sqTolerance, simplified) {
    var maxSqDist = sqTolerance,
        index;

    for (var i = first + 1; i < last; i++) {
        var sqDist = getSqSegDist(points[i], points[first], points[last]);

        if (sqDist > maxSqDist) {
            index = i;
            maxSqDist = sqDist;
        }
    }

    if (maxSqDist > sqTolerance) {
        if (index - first > 1) simplifyDPStep(points, first, index, sqTolerance, simplified);
        simplified.push(points[index]);
        if (last - index > 1) simplifyDPStep(points, index, last, sqTolerance, simplified);
    }
}

// simplification using Ramer-Douglas-Peucker algorithm
function simplifyDouglasPeucker(points, sqTolerance) {
    var last = points.length - 1;

    var simplified = [points[0]];
    simplifyDPStep(points, 0, last, sqTolerance, simplified);
    simplified.push(points[last]);

    return simplified;
}

// both algorithms combined for awesome performance
function simplify(points, tolerance, highestQuality) {

    if (points.length <= 2) return points;

    var sqTolerance = tolerance !== undefined ? tolerance * tolerance : 1;

    points = highestQuality ? points : simplifyRadialDist(points, sqTolerance);
    points = simplifyDouglasPeucker(points, sqTolerance);

    return points;
}

function dot(v1, v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
function v2Dot(v1, v2) {
    return v1[0] * v2[0] + v1[1] * v2[1];
}

function normalize(out, v) {
    const x = v[0];
    const y = v[1];
    const z = v[2];
    const d = Math.sqrt(x * x + y * y + z * z);
    out[0] = x / d;
    out[1] = y / d;
    out[2] = z / d;
    return out;
}

function v2Normalize(out, v) {
    const x = v[0];
    const y = v[1];
    const d = Math.sqrt(x * x + y * y);
    out[0] = x / d;
    out[1] = y / d;
    return out;
}

function scale(out, v, s) {
    out[0] = v[0] * s;
    out[1] = v[1] * s;
    out[2] = v[2] * s;
    return out;
}



function scaleAndAdd(out, v1, v2, s) {
    out[0] = v1[0] + v2[0] * s;
    out[1] = v1[1] + v2[1] * s;
    out[2] = v1[2] + v2[2] * s;
    return out;
}



function v2Add(out, v1, v2) {
    out[0] = v1[0] + v2[0];
    out[1] = v1[1] + v2[1];
    return out;
}





function v3Sub(out, v1, v2) {
    out[0] = v1[0] - v2[0];
    out[1] = v1[1] - v2[1];
    out[2] = v1[2] - v2[2];
    return out;
}

function v3Normalize(out, v) {
    const x = v[0];
    const y = v[1];
    const z = v[2];
    const d = Math.sqrt(x * x + y * y + z * z);
    out[0] = x / d;
    out[1] = y / d;
    out[2] = z / d;
    return out;
}

function v3Cross(out, v1, v2) {
    var ax = v1[0], ay = v1[1], az = v1[2],
        bx = v2[0], by = v2[1], bz = v2[2];

    out[0] = ay * bz - az * by;
    out[1] = az * bx - ax * bz;
    out[2] = ax * by - ay * bx;
    return out;
}

const rel = [];
// start and end must be normalized
function slerp(out, start, end, t) {
    // https://keithmaggio.wordpress.com/2011/02/15/math-magician-lerp-slerp-and-nlerp/
    const cosT = dot(start, end);
    const theta = Math.acos(cosT) * t;

    scaleAndAdd(rel, end, start, -cosT);
    normalize(rel, rel);// start and rel Orthonormal basis

    scale(out, start, Math.cos(theta));
    scaleAndAdd(out, out, rel, Math.sin(theta));

    return out;
}



function area$1(points, start, end) {
    // Signed polygon area
    const n = end - start;
    if (n < 3) {
        return 0;
    }
    let area = 0;
    for (let i = (end - 1) * 2, j = start * 2; j < end * 2;) {
        const x0 = points[i];
        const y0 = points[i + 1];
        const x1 = points[j];
        const y1 = points[j + 1];
        i = j;
        j += 2;
        area += x0 * y1 - x1 * y0;
    }

    return area;
}

// TODO fitRect x, y are negative?
// TODO Dimensions
// TODO bevel="top"|"bottom"

function triangulate(vertices, holes, dimensions=2) {
    return earcut_1(vertices, holes, dimensions);
}

function flatten(data) {
    return earcut_1.flatten(data);
}

const v1 = [];
const v2 = [];
const v = [];

function innerOffsetPolygon(
    vertices, out, start, end, outStart, offset, miterLimit, close
) {
    const checkMiterLimit = miterLimit != null;
    let outOff = outStart;
    let indicesMap = null;
    if (checkMiterLimit) {
        indicesMap = new Uint32Array(end - start);
    }
    for (let i = start; i < end; i++) {
        const nextIdx = i === end - 1 ? start : i + 1;
        const prevIdx = i === start ? end - 1 : i - 1;
        const x1 = vertices[prevIdx * 2];
        const y1 = vertices[prevIdx * 2 + 1];
        const x2 = vertices[i * 2];
        const y2 = vertices[i * 2 + 1];
        const x3 = vertices[nextIdx * 2];
        const y3 = vertices[nextIdx * 2 + 1];

        v1[0] = x2 - x1;
        v1[1] = y2 - y1;
        v2[0] = x3 - x2;
        v2[1] = y3 - y2;

        v2Normalize(v1, v1);
        v2Normalize(v2, v2);

        checkMiterLimit && (indicesMap[i] = outOff);
        if (!close && i === start) {
            v[0] = v2[1];
            v[1] = -v2[0];
            v2Normalize(v, v);
            out[outOff * 2] = x2 + v[0] * offset;
            out[outOff * 2 + 1] = y2 + v[1] * offset;
            outOff++;
        }
        else if (!close && i === end - 1) {
            v[0] = v1[1];
            v[1] = -v1[0];
            v2Normalize(v, v);
            out[outOff * 2] = x2 + v[0] * offset;
            out[outOff * 2 + 1] = y2 + v[1] * offset;
            outOff++;
        }
        else {
            // PENDING Why using sub will lost the direction info.
            v2Add(v, v2, v1);
            const tmp = v[1];
            v[1] = -v[0];
            v[0] = tmp;

            v2Normalize(v, v);

            const cosA = v2Dot(v, v2);
            const sinA = Math.sqrt(1 - cosA * cosA);
            // PENDING
            const miter = offset * Math.min(10, 1 / sinA);

            const isCovex = offset * cosA < 0;

            if (checkMiterLimit && (1 / sinA) > miterLimit && isCovex) {
                const mx = x2 + v[0] * offset;
                const my = y2 + v[1] * offset;
                const halfA = Math.acos(sinA) / 2;
                const dist = Math.tan(halfA) * Math.abs(offset);
                out[outOff * 2] = mx + v[1] * dist;
                out[outOff * 2 + 1] = my - v[0] * dist;
                outOff++;
                out[outOff * 2] = mx - v[1] * dist;
                out[outOff * 2 + 1] = my + v[0] * dist;
                outOff++;
            }
            else {
                out[outOff * 2] = x2 + v[0] * miter;
                out[outOff * 2 + 1] = y2 + v[1] * miter;
                outOff++;
            }
        }
    }

    return indicesMap;
}

function offsetPolygon(vertices, holes, offset, miterLimit, close) {
    const offsetVertices = miterLimit != null ? [] : new Float32Array(vertices.length);
    const exteriorSize = (holes && holes.length) ? holes[0] : vertices.length / 2;

    innerOffsetPolygon(
        vertices, offsetVertices, 0, exteriorSize, 0, offset, miterLimit, close, false
    );

    if (holes) {
        for (let i = 0; i < holes.length; i++) {
            const start = holes[i];
            const end = holes[i + 1] || vertices.length / 2;
            innerOffsetPolygon(
                vertices, offsetVertices, start, end,
                miterLimit != null ? offsetVertices.length / 2 : start,
                offset, miterLimit, close
            );
        }
    }

    return offsetVertices;
}

function reversePoints(points, stride, start, end) {
    for (let i = 0; i < Math.floor((end - start) / 2); i++) {
        for (let j = 0; j < stride; j++) {
            const a = (i + start) * stride + j;
            const b = (end - i - 1) * stride + j;
            const tmp = points[a];
            points[a] = points[b];
            points[b] = tmp;
        }
    }

    return points;
}

function convertToClockwise(vertices, holes) {
    let polygonVertexCount = vertices.length / 2;
    let start = 0;
    let end = holes && holes.length ? holes[0] : polygonVertexCount;
    if (area$1(vertices, start, end) > 0) {
        reversePoints(vertices, 2, start, end);
    }
    for (let h = 1; h < (holes ? holes.length : 0) + 1; h++) {
        start = holes[h - 1];
        end = holes[h] || polygonVertexCount;
        if (area$1(vertices, start, end) < 0) {
            reversePoints(vertices, 2, start, end);
        }
    }
}

function normalizeOpts(opts) {

    opts.depth = opts.depth || 1;
    opts.bevelSize = opts.bevelSize || 0;
    opts.bevelSegments = opts.bevelSegments == null ? 2 : opts.bevelSegments;
    opts.smoothSide = opts.smoothSide || false;
    opts.smoothBevel = opts.smoothBevel || false;
    opts.simplify = opts.simplify || 0;

    // Normalize bevel options.
    if (typeof opts.depth === 'number') {
        opts.bevelSize = Math.min(!(opts.bevelSegments > 0) ? 0 : opts.bevelSize, opts.depth / 2);
    }
    if (!(opts.bevelSize > 0)) {
        opts.bevelSegments = 0;
    }
    opts.bevelSegments = Math.round(opts.bevelSegments);

    const boundingRect = opts.boundingRect;
    opts.translate = opts.translate || [0, 0];
    opts.scale = opts.scale || [1, 1];
    if (opts.fitRect) {
        let targetX = opts.fitRect.x == null
            ? (boundingRect.x || 0)
            : opts.fitRect.x;
        let targetY = opts.fitRect.y == null
            ? (boundingRect.y || 0)
            : opts.fitRect.y;
        let targetWidth = opts.fitRect.width;
        let targetHeight = opts.fitRect.height;
        if (targetWidth == null) {
            if (targetHeight != null) {
                targetWidth = targetHeight / boundingRect.height * boundingRect.width;
            }
            else {
                targetWidth = boundingRect.width;
                targetHeight = boundingRect.height;
            }
        }
        else if (targetHeight == null) {
            targetHeight = targetWidth / boundingRect.width * boundingRect.height;
        }
        opts.scale = [
            targetWidth / boundingRect.width,
            targetHeight / boundingRect.height
        ];
        opts.translate = [
            (targetX - boundingRect.x) * opts.scale[0],
            (targetY - boundingRect.y) * opts.scale[1]
        ];
    }
}

function generateNormal(indices, position) {

    function v3Set(p, a, b, c) {
        p[0] = a; p[1] = b; p[2] = c;
    }

    const p1 = [];
    const p2 = [];
    const p3 = [];

    const v21 = [];
    const v32 = [];

    const n = [];

    const len = indices.length;
    const normals = new Float32Array(position.length);
    for (let f = 0; f < len;) {
        const i1 = indices[f++] * 3;
        const i2 = indices[f++] * 3;
        const i3 = indices[f++] * 3;

        v3Set(p1, position[i1], position[i1 + 1], position[i1 + 2]);
        v3Set(p2, position[i2], position[i2 + 1], position[i2 + 2]);
        v3Set(p3, position[i3], position[i3 + 1], position[i3 + 2]);

        v3Sub(v21, p1, p2);
        v3Sub(v32, p2, p3);
        v3Cross(n, v21, v32);
        // Already be weighted by the triangle area
        for (let i = 0; i < 3; i++) {
            normals[i1 + i] = normals[i1 + i] + n[i];
            normals[i2 + i] = normals[i2 + i] + n[i];
            normals[i3 + i] = normals[i3 + i] + n[i];
        }
    }

    for (var i = 0; i < normals.length;) {
        v3Set(n, normals[i], normals[i+1], normals[i+2]);
        v3Normalize(n, n);
        normals[i++] = n[0];
        normals[i++] = n[1];
        normals[i++] = n[2];
    }

    return normals;
}
// 0,0----1,0
// 0,1----1,1
const quadToTriangle = [
    [0, 0], [1, 0], [1, 1],
    [0, 0], [1, 1], [0, 1]
];

// Add side vertices and indices. Include bevel.
function addExtrudeSide(
    out, {vertices, topVertices, depth, rect}, start, end,
    cursors, opts
) {
    const ringVertexCount = end - start;
    const splitSide = opts.smoothSide ? 1 : 2;
    const splitRingVertexCount = ringVertexCount * splitSide;

    const splitBevel = opts.smoothBevel ? 1 : 2;
    const bevelSize = Math.min(depth / 2, opts.bevelSize);
    const bevelSegments = opts.bevelSegments;
    const vertexOffset = cursors.vertex;
    const size = Math.max(rect.width, rect.height);

    // Side vertices
    if (bevelSize > 0) {

        const v0 = [0, 0, 1];
        const v1 = [];
        const v2 = [0, 0, -1];
        const v = [];

        let ringCount = 0;
        let vLen = new Float32Array(ringVertexCount);
        for (let k = 0; k < 2; k++) {
            const z = (k === 0 ? (depth - bevelSize) : bevelSize);
            for (let s = 0; s <= bevelSegments * splitBevel; s++) {
                let uLen = 0;
                let prevX;
                let prevY;
                for (let i = 0; i < ringVertexCount; i++) {

                    for (let j = 0; j < splitSide; j++) {
                        // TODO Cache and optimize
                        let idx = ((i + j) % ringVertexCount + start) * 2;
                        v1[0] = vertices[idx] - topVertices[idx];
                        v1[1] = vertices[idx + 1] - topVertices[idx + 1];
                        v1[2] = 0;
                        const l = Math.sqrt(v1[0] * v1[0] + v1[1] * v1[1]);
                        v1[0] /= l;
                        v1[1] /= l;

                        const t = (Math.floor(s / splitBevel) + (s % splitBevel)) / bevelSegments;
                        k === 0 ? slerp(v, v0, v1, t)
                            : slerp(v, v1, v2, t);

                        const t2 = k === 0  ? t : 1 - t;
                        const a = bevelSize * Math.sin(t2 * Math.PI / 2);
                        const b = l * Math.cos(t2 * Math.PI / 2);

                        // ellipse radius
                        const r = bevelSize * l / Math.sqrt(a * a + b * b);

                        const x = v[0] * r + topVertices[idx];
                        const y = v[1] * r + topVertices[idx + 1];
                        const zz = v[2] * r + z;
                        out.position[cursors.vertex * 3] = x;
                        out.position[cursors.vertex * 3 + 1] = y;
                        out.position[cursors.vertex * 3 + 2] = zz;

                        // TODO Cache and optimize
                        if (i > 0 || j > 0) {
                            uLen += Math.sqrt((prevX - x) * (prevX - x) + (prevY - y) * (prevY - y));
                        }
                        if (s > 0 || k > 0) {
                            let tmp = (cursors.vertex - splitRingVertexCount) * 3;
                            let prevX2 = out.position[tmp];
                            let prevY2 = out.position[tmp + 1];
                            let prevZ2 = out.position[tmp + 2];

                            vLen[i] += Math.sqrt(
                                (prevX2 - x) * (prevX2 - x)
                                + (prevY2 - y) * (prevY2 - y)
                                + (prevZ2 - zz) * (prevZ2 - zz)
                            );
                        }
                        out.uv[cursors.vertex * 2] = uLen / size;
                        out.uv[cursors.vertex * 2 + 1] = vLen[i] / size;

                        prevX = x;
                        prevY = y;
                        cursors.vertex++;
                    }

                    if ((splitBevel > 1 && (s % splitBevel)) || (splitBevel === 1 && s >= 1)) {
                        for (let f = 0; f < 6; f++) {
                            const m = (quadToTriangle[f][0] + i * splitSide) % splitRingVertexCount;
                            const n = quadToTriangle[f][1] + ringCount;
                            out.indices[cursors.index++] = (n - 1) * splitRingVertexCount + m + vertexOffset;
                        }
                    }
                }

                ringCount++;
            }
        }
    }
    else {
        for (let k = 0; k < 2; k++) {
            const z = k === 0 ? depth - bevelSize : bevelSize;
            let uLen = 0;
            let prevX;
            let prevY;
            for (let i = 0; i < ringVertexCount; i++) {
                for (let m = 0; m < splitSide; m++) {
                    const idx = ((i + m) % ringVertexCount + start) * 2;
                    const x = vertices[idx];
                    const y = vertices[idx + 1];
                    out.position[cursors.vertex * 3] = x;
                    out.position[cursors.vertex * 3 + 1] = y;
                    out.position[cursors.vertex * 3 + 2] = z;
                    if (i > 0 || m > 0) {
                        uLen += Math.sqrt((prevX - x) * (prevX - x) + (prevY - y) * (prevY - y));
                    }
                    out.uv[cursors.vertex * 2] = uLen / size;
                    out.uv[cursors.vertex * 2 + 1] = z / size;
                    prevX = x;
                    prevY = y;

                    cursors.vertex++;
                }
            }
        }
    }
    // Connect the side
    const sideStartRingN = bevelSize > 0 ? (bevelSegments * splitBevel + 1) : 1;
    for (let i = 0; i < ringVertexCount; i++) {
        for (let f = 0; f < 6; f++) {
            const m = (quadToTriangle[f][0] + i * splitSide) % splitRingVertexCount;
            const n = quadToTriangle[f][1] + sideStartRingN;
            out.indices[cursors.index++] = (n - 1) * splitRingVertexCount + m + vertexOffset;
        }
    }
}

function addTopAndBottom({indices, vertices, topVertices, rect, depth}, out, cursors, opts) {
    if (vertices.length <= 4) {
        return;
    }

    const vertexOffset = cursors.vertex;
    // Top indices
    const indicesLen = indices.length;
    for (let i = 0; i < indicesLen; i++) {
        out.indices[cursors.index++] = vertexOffset + indices[i];
    }
    const size = Math.max(rect.width, rect.height);
    // Top and bottom vertices
    for (let k = 0; k < (opts.excludeBottom ? 1 : 2); k++) {
        for (let i = 0; i < topVertices.length; i += 2) {
            const x = topVertices[i];
            const y = topVertices[i + 1];
            out.position[cursors.vertex * 3] = x;
            out.position[cursors.vertex * 3 + 1] = y;
            out.position[cursors.vertex * 3 + 2] = (1 - k) * depth;

            out.uv[cursors.vertex * 2] = (x - rect.x) / size;
            out.uv[cursors.vertex * 2 + 1] = (y - rect.y) / size;
            cursors.vertex++;
        }
    }
    // Bottom indices
    if (!opts.excludeBottom) {
        const vertexCount = vertices.length / 2;
        for (let i = 0; i < indicesLen; i += 3) {
            for (let k = 0; k < 3; k++) {
                out.indices[cursors.index++] = vertexOffset + vertexCount + indices[i + 2 - k];
            }
        }
    }
}


function innerExtrudeTriangulatedPolygon(preparedData, opts) {
    let indexCount = 0;
    let vertexCount = 0;

    for (let p = 0; p < preparedData.length; p++) {
        const {indices, vertices, holes, depth} = preparedData[p];
        const polygonVertexCount = vertices.length / 2;
        const bevelSize = Math.min(depth / 2, opts.bevelSize);
        const bevelSegments = !(bevelSize > 0) ? 0 : opts.bevelSegments;

        indexCount += indices.length * (opts.excludeBottom ? 1 : 2);
        vertexCount += polygonVertexCount * (opts.excludeBottom ? 1 : 2);
        const ringCount = 2 + bevelSegments * 2;

        let start = 0;
        let end = 0;
        for (let h = 0; h < (holes ? holes.length : 0) + 1; h++) {
            if (h === 0) {
                end = holes && holes.length ? holes[0] : polygonVertexCount;
            }
            else {
                start = holes[h - 1];
                end = holes[h] || polygonVertexCount;
            }

            indexCount += (end - start) * 6 * (ringCount - 1);

            const sideRingVertexCount = (end - start) * (opts.smoothSide ? 1 : 2);
            vertexCount += sideRingVertexCount * ringCount
                // Double the bevel vertex number if not smooth
                + (!opts.smoothBevel ? bevelSegments * sideRingVertexCount * 2 : 0);
        }
    }

    const data = {
        position: new Float32Array(vertexCount * 3),
        indices: new (vertexCount > 0xffff ? Uint32Array : Uint16Array)(indexCount),
        uv: new Float32Array(vertexCount * 2)
    };

    const cursors = {
        vertex: 0, index: 0
    };

    for (let d = 0; d < preparedData.length; d++) {
        addTopAndBottom(preparedData[d], data, cursors, opts);
    }

    for (let d = 0; d < preparedData.length; d++) {
        const {holes, vertices} = preparedData[d];
        const topVertexCount = vertices.length / 2;

        let start = 0;
        let end = (holes && holes.length) ? holes[0] : topVertexCount;
        // Add exterior
        addExtrudeSide(data, preparedData[d], start, end, cursors, opts);
        // Add holes
        if (holes) {
            for (let h = 0; h < holes.length; h++) {
                start = holes[h];
                end = holes[h + 1] || topVertexCount;
                addExtrudeSide(data, preparedData[d], start, end, cursors, opts);
            }
        }
    }

    // Wrap uv
    for (let i = 0; i < data.uv.length; i++) {
        const val = data.uv[i];
        if (val > 0 && Math.round(val) === val) {
            data.uv[i] = 1;
        }
        else {
            data.uv[i] = val % 1;
        }
    }

    data.normal = generateNormal(data.indices, data.position);
    // PENDING
    data.boundingRect = preparedData[0] && preparedData[0].rect;

    return data;
}

function convertPolylineToTriangulatedPolygon(polyline, polylineIdx, opts) {
    const lineWidth = opts.lineWidth;
    const pointCount = polyline.length;
    const points = new Float32Array(pointCount * 2);
    const translate = opts.translate || [0, 0];
    const scale$$1 = opts.scale || [1, 1];
    for (let i = 0, k = 0; i < pointCount; i++) {
        points[k++] = polyline[i][0] * scale$$1[0] + translate[0];
        points[k++] = polyline[i][1] * scale$$1[1] + translate[1];
    }

    if (area$1(points, 0, pointCount) < 0) {
        reversePoints(points, 2, 0, pointCount);
    }

    const insidePoints = [];
    const outsidePoints = [];
    const miterLimit = opts.miterLimit;
    const outsideIndicesMap = innerOffsetPolygon(
        points, outsidePoints, 0, pointCount, 0, -lineWidth / 2, miterLimit, false
    );
    reversePoints(points, 2, 0, pointCount);
    const insideIndicesMap = innerOffsetPolygon(
        points, insidePoints, 0, pointCount, 0, -lineWidth / 2, miterLimit, false
    );

    const polygonVertexCount = (insidePoints.length + outsidePoints.length) / 2;
    const polygonVertices = new Float32Array(polygonVertexCount * 2);

    let offset = 0;
    const outsidePointCount = outsidePoints.length / 2;
    for (let i = 0; i < outsidePoints.length; i++) {
        polygonVertices[offset++] = outsidePoints[i];
    }
    for (let i = 0; i < insidePoints.length; i++) {
        polygonVertices[offset++] = insidePoints[i];
    }

    // Built indices
    const indices = new (polygonVertexCount > 0xffff ? Uint32Array : Uint16Array)(
        ((pointCount - 1) * 2 + (polygonVertexCount - pointCount * 2)) * 3
    );
    let off = 0;
    for (let i = 0; i < pointCount - 1; i++) {
        const i2 = i + 1;
        indices[off++] = outsidePointCount - 1 - outsideIndicesMap[i];
        indices[off++] = outsidePointCount - 1 - outsideIndicesMap[i] - 1;
        indices[off++] = insideIndicesMap[i] + 1 + outsidePointCount;

        indices[off++] = outsidePointCount - 1 - outsideIndicesMap[i];
        indices[off++] = insideIndicesMap[i] + 1 + outsidePointCount;
        indices[off++] = insideIndicesMap[i] + outsidePointCount;

        if (insideIndicesMap[i2] - insideIndicesMap[i] === 2) {
            indices[off++] = insideIndicesMap[i] + 2 + outsidePointCount;
            indices[off++] = insideIndicesMap[i] + 1 + outsidePointCount;
            indices[off++] = outsidePointCount - outsideIndicesMap[i2] - 1;
        }
        else if (outsideIndicesMap[i2] - outsideIndicesMap[i] === 2) {
            indices[off++] = insideIndicesMap[i2] + outsidePointCount;
            indices[off++] = outsidePointCount - 1 - (outsideIndicesMap[i] + 1);
            indices[off++] = outsidePointCount - 1 - (outsideIndicesMap[i] + 2);
        }
    }

    const topVertices = opts.bevelSize > 0
        ? offsetPolygon(polygonVertices, [], opts.bevelSize, null, true) : polygonVertices;
    const boundingRect = opts.boundingRect;
    return {
        vertices: polygonVertices,
        indices,
        topVertices,
        rect: {
            x: boundingRect.x * scale$$1[0] + translate[0],
            y: boundingRect.y * scale$$1[1] + translate[1],
            width: boundingRect.width * scale$$1[0],
            height: boundingRect.height * scale$$1[1],
        },
        depth: typeof opts.depth === 'function' ? opts.depth(polylineIdx) : opts.depth,
        holes: []
    };
}

function removeClosePointsOfPolygon(polygon, epsilon) {
    const newPolygon = [];
    for (let k  = 0; k < polygon.length; k++) {
        const points = polygon[k];
        const newPoints = [];
        const len = points.length;
        let x1 = points[len - 1][0];
        let y1 = points[len - 1][1];
        let dist = 0;
        for (let i = 0; i < len; i++) {
            let x2 = points[i][0];
            let y2 = points[i][1];
            const dx = x2 - x1;
            const dy = y2 - y1;
            dist += Math.sqrt(dx * dx + dy * dy);
            if (dist > epsilon) {
                newPoints.push(points[i]);
                dist = 0;
            }
            x1 = x2;
            y1 = y2;
        }
        if (newPoints.length >= 3) {
            newPolygon.push(newPoints);
        }
    }
    return newPolygon.length > 0 ? newPolygon : null;
}

function simplifyPolygon(polygon, tolerance) {
    const newPolygon = [];
    for (let k  = 0; k < polygon.length; k++) {
        let points = polygon[k];
        points = simplify(points, tolerance, true);
        if (points.length >= 3) {
            newPolygon.push(points);
        }
    }
    return newPolygon.length > 0 ? newPolygon : null;
}
/**
 *
 * @param {Array} polygons Polygons array that match GeoJSON MultiPolygon geometry.
 * @param {Object} [opts]
 * @param {number|Function} [opts.depth]
 * @param {number} [opts.bevelSize = 0]
 * @param {number} [opts.bevelSegments = 2]
 * @param {number} [opts.simplify = 0]
 * @param {boolean} [opts.smoothSide = false]
 * @param {boolean} [opts.smoothBevel = false]
 * @param {boolean} [opts.excludeBottom = false]
 * @param {Object} [opts.fitRect] translate and scale will be ignored if fitRect is set
 * @param {Array} [opts.translate]
 * @param {Array} [opts.scale]
 *
 * @return {Object} {indices, position, uv, normal, boundingRect}
 */
function extrudePolygon(polygons, opts) {

    opts = Object.assign({}, opts);

    const min = [Infinity, Infinity];
    const max = [-Infinity, -Infinity];
    for (let i = 0; i < polygons.length; i++) {
        updateBoundingRect(polygons[i][0], min, max);
    }
    opts.boundingRect = opts.boundingRect || {
        x: min[0], y: min[1], width: max[0] - min[0], height: max[1] - min[1]
    };

    normalizeOpts(opts);

    const preparedData = [];
    const translate = opts.translate || [0, 0];
    const scale$$1 = opts.scale || [1, 1];
    const boundingRect = opts.boundingRect;
    const transformdRect = {
        x: boundingRect.x * scale$$1[0] + translate[0],
        y: boundingRect.y * scale$$1[1] + translate[1],
        width: boundingRect.width * scale$$1[0],
        height: boundingRect.height * scale$$1[1],
    };

    const epsilon = Math.min(
        boundingRect.width, boundingRect.height
    ) / 1e5;
    for (let i = 0; i < polygons.length; i++) {
        let newPolygon = removeClosePointsOfPolygon(polygons[i], epsilon);
        if (!newPolygon) {
            continue;
        }
        const simplifyTolerance = opts.simplify / Math.max(scale$$1[0], scale$$1[1]);
        if (simplifyTolerance > 0) {
            newPolygon = simplifyPolygon(newPolygon, simplifyTolerance);
        }
        if (!newPolygon) {
            continue;
        }

        const {vertices, holes, dimensions} = earcut_1.flatten(newPolygon);

        for (let k = 0; k < vertices.length;) {
            vertices[k] = vertices[k++] * scale$$1[0] + translate[0];
            vertices[k] = vertices[k++] * scale$$1[1] + translate[1];
        }

        convertToClockwise(vertices, holes);

        if (dimensions !== 2) {
            throw new Error('Only 2D polygon points are supported');
        }
        const topVertices = opts.bevelSize > 0
            ? offsetPolygon(vertices, holes, opts.bevelSize, null, true) : vertices;
        const indices = triangulate(topVertices, holes, dimensions);
        preparedData.push({
            indices,
            vertices,
            topVertices,
            holes,
            rect: transformdRect,
            depth: typeof opts.depth === 'function' ? opts.depth(i) : opts.depth
        });
    }
    return innerExtrudeTriangulatedPolygon(preparedData, opts);
}

/**
 *
 * @param {Array} polylines Polylines array that match GeoJSON MultiLineString geometry.
 * @param {Object} [opts]
 * @param {number} [opts.depth]
 * @param {number} [opts.bevelSize = 0]
 * @param {number} [opts.bevelSegments = 2]
 * @param {number} [opts.simplify = 0]
 * @param {boolean} [opts.smoothSide = false]
 * @param {boolean} [opts.smoothBevel = false]
 * @param {boolean} [opts.excludeBottom = false]
 * @param {boolean} [opts.lineWidth = 1]
 * @param {boolean} [opts.miterLimit = 2]
 * @param {Object} [opts.fitRect] translate and scale will be ignored if fitRect is set
 * @param {Array} [opts.translate]
 * @param {Array} [opts.scale]
 * @param {Object} [opts.boundingRect]
 * @return {Object} {indices, position, uv, normal, boundingRect}
 */
function extrudePolyline(polylines, opts) {

    opts = Object.assign({}, opts);

    const min = [Infinity, Infinity];
    const max = [-Infinity, -Infinity];
    for (let i = 0; i < polylines.length; i++) {
        updateBoundingRect(polylines[i], min, max);
    }
    opts.boundingRect = opts.boundingRect || {
        x: min[0], y: min[1], width: max[0] - min[0], height: max[1] - min[1]
    };

    normalizeOpts(opts);
    const scale$$1 = opts.scale || [1, 1];

    if (opts.lineWidth == null) {
        opts.lineWidth = 1;
    }
    if (opts.miterLimit == null) {
        opts.miterLimit = 2;
    }
    const preparedData = [];
    // Extrude polyline to polygon
    for (let i = 0; i < polylines.length; i++) {
        let newPolyline = polylines[i];
        const simplifyTolerance = opts.simplify / Math.max(scale$$1[0], scale$$1[1]);
        if (simplifyTolerance > 0) {
            newPolyline = simplify(newPolyline, simplifyTolerance, true);
        }
        preparedData.push(convertPolylineToTriangulatedPolygon(newPolyline, i, opts));
    }

    return innerExtrudeTriangulatedPolygon(preparedData, opts);
}

function updateBoundingRect(points, min, max) {
    for (let i = 0; i < points.length; i++) {
        min[0] = Math.min(points[i][0], min[0]);
        min[1] = Math.min(points[i][1], min[1]);
        max[0] = Math.max(points[i][0], max[0]);
        max[1] = Math.max(points[i][1], max[1]);
    }
}

/**
 *
 * @param {Object} geojson
 * @param {Object} [opts]
 * @param {number} [opts.depth]
 * @param {number} [opts.bevelSize = 0]
 * @param {number} [opts.bevelSegments = 2]
 * @param {number} [opts.simplify = 0]
 * @param {boolean} [opts.smoothSide = false]
 * @param {boolean} [opts.smoothBevel = false]
 * @param {boolean} [opts.excludeBottom = false]
 * @param {boolean} [opts.lineWidth = 1]
 * @param {boolean} [opts.miterLimit = 2]
 * @param {Object} [opts.fitRect] translate and scale will be ignored if fitRect is set
 * @param {Array} [opts.translate]
 * @param {Array} [opts.scale]
 * @param {Object} [opts.boundingRect]
 * @return {Object} {polyline: {indices, position, uv, normal}, polygon: {indices, position, uv, normal}}
 */

 // TODO Not merge feature
function extrudeGeoJSON(geojson, opts) {

    opts = Object.assign({}, opts);

    const polylines = [];
    const polygons = [];

    const polylineFeatureIndices = [];
    const polygonFeatureIndices = [];

    const min = [Infinity, Infinity];
    const max = [-Infinity, -Infinity];

    for (let i = 0; i < geojson.features.length; i++) {
        const feature = geojson.features[i];
        const geometry = feature.geometry;
        if (geometry && geometry.coordinates) {
            switch (geometry.type) {
                case 'LineString':
                    polylines.push(geometry.coordinates);
                    polylineFeatureIndices.push(i);
                    updateBoundingRect(geometry.coordinates, min, max);
                    break;
                case 'MultiLineString':
                    for (let k = 0; k < geometry.coordinates.length; k++) {
                        polylines.push(geometry.coordinates[k]);
                        polylineFeatureIndices.push(i);
                        updateBoundingRect(geometry.coordinates[k], min, max);
                    }
                    break;
                case 'Polygon':
                    polygons.push(geometry.coordinates);
                    polygonFeatureIndices.push(i);
                    updateBoundingRect(geometry.coordinates[0], min, max);
                    break;
                case 'MultiPolygon':
                    for (let k = 0; k < geometry.coordinates.length; k++) {
                        polygons.push(geometry.coordinates[k]);
                        polygonFeatureIndices.push(i);
                        updateBoundingRect(geometry.coordinates[k][0], min, max);
                    }
                    break;
            }
        }
    }

    opts.boundingRect = opts.boundingRect || {
        x: min[0], y: min[1], width: max[0] - min[0], height: max[1] - min[1]
    };

    const originalDepth = opts.depth;
    return {
        polyline: extrudePolyline(polylines, Object.assign(opts, {
            depth: function (idx) {
                if (typeof originalDepth === 'function') {
                    return originalDepth(
                        geojson.features[polylineFeatureIndices[idx]]
                    );
                }
                return originalDepth;
            }
        })),
        polygon: extrudePolygon(polygons, Object.assign(opts, {
            depth: function (idx) {
                if (typeof originalDepth === 'function') {
                    return originalDepth(
                        geojson.features[polygonFeatureIndices[idx]]
                    );
                }
                return originalDepth;
            }
        }))
    };
}

exports.triangulate = triangulate;
exports.flatten = flatten;
exports.offsetPolygon = offsetPolygon;
exports.extrudePolygon = extrudePolygon;
exports.extrudePolyline = extrudePolyline;
exports.extrudeGeoJSON = extrudeGeoJSON;

Object.defineProperty(exports, '__esModule', { value: true });

})));
