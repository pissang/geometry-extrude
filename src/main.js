// TODO fitRect x, y are negative?
// TODO Extrude dimensions
// TODO bevel="top"|"bottom"
// TODO Not add top and bottom vertices if area is 0
// TODO Calculate max bevel size by calculating intersection of offset vector with all edges.

import earcut from 'earcut';
import doSimplify from './simplify';
import {
    slerp, v2Normalize, v2Dot, v2Add, area,
    v3Normalize, v3Sub, v3Cross, lineIntersection, v2Sub, v2Dist, dot
} from './math';

export function triangulate(vertices, holes, dimensions=2) {
    return earcut(vertices, holes, dimensions);
};

export function flatten(data) {
    return earcut.flatten(data);
}

const v1 = [];
const v2 = [];
const v = [];

function innerAppendArc(out, start, end, center, segments, radius) {
    const PI2 = Math.PI * 2;
    let startAngle = Math.atan2(start[1] - center[1], start[0] - center[0]);
    let endAngle = Math.atan2(end[1] - center[1], end[0] - center[0]);

    if (startAngle > endAngle) {
        endAngle += PI2;
    }
    // Make sure use the smaller arc.
    if (endAngle - startAngle > Math.PI) {
        endAngle -= PI2;
    }

    for (let i = 0; i < segments; i++) {
        let angle = startAngle + i / (segments - 1) * (endAngle - startAngle);
        out.push(center[0] + Math.cos(angle) * radius);
        out.push(center[1] + Math.sin(angle) * radius);
    }
}

function innerOffsetContour(
    vertices, out, start, end, outStart, offset,
    join,   // Default to be miter, can be round.
    miterLimit, close
    // offsetLines
) {
    const needsIndicesMap = miterLimit != null;
    const isComplexJoin = (join && join !== 'miter') || miterLimit;
    const isDynamicArray = isComplexJoin;
    if (isDynamicArray && !out.push) {
        throw new Error('output must be a dynamic array if miterLimit is given');
    }

    let indicesMap = needsIndicesMap ? new Uint32Array(end - start) : null;

    let edges = [];
    // Polygon points have been removed duplication.
    for (let i = start; i < (close ? end : end - 1); i++) {
        const nextIdx = i === end - 1 ? start : i + 1;
        const x1 = vertices[i * 2];
        const y1 = vertices[i * 2 + 1];
        const x2 = vertices[nextIdx * 2];
        const y2 = vertices[nextIdx * 2 + 1];

        let dx1 = y2 - y1;
        let dy1 = x1 - x2;

        const l1 = Math.sqrt(dx1 * dx1 + dy1 * dy1) / offset;

        dx1 /= l1;
        dy1 /= l1;

        edges.push([x1 + dx1, y1 + dy1, x2 + dx1, y2 + dy1])
    }
    // Save raw offseted edges for calculating lineJoin.
    const offsetedEdges = edges.slice();

    const p1 = [];  // Intersection p1
    const p2 = [];  // Intersection p2
    const p0 = [];  // Edge end point
    const p3 = [];  // Edge end point

    let edgeIndices = [];
    for (let i = 0; i < edges.length; i++) {
        const edge = edges[i];
        const dx = edge[2] - edge[0];
        const dy = edge[3] - edge[1];
        const dist = dx * dx + dy * dy;
        // Remove too short or invalid edge.
        if (dist < 1e-8 || isNaN(dist)) {
            edges[i] = null;
            continue;
        }

        edgeIndices.push(i);
    }
    let hasRemoved = false;


    function removeEdgeIntersections() {
        hasRemoved = false;

        const sortedEdges = [];
        const edgesCount = edgeIndices.length;
        for (let i = 0; i < edgesCount; i++) {
            const curr = edges[edgeIndices[i]];
            const dx = curr[2] - curr[0];
            const dy = curr[3] - curr[1];
            const dist = dx * dx + dy * dy;
            sortedEdges.push([i, dist]);
        }
        // Sort edges from short to long.
        // Reemove shorter edges at first.
        sortedEdges.sort((a, b) => a[1] - b[1]);

        for (let i = 0; i < edgesCount; i++) {
            const edgeInfo = sortedEdges[i];
            const idx = edgeInfo[0];
            const curr = edges[edgeIndices[idx]];
            let nextIdx = idx;
            let prevIdx = idx;
            do {
                nextIdx = (nextIdx + 1) % edgesCount;
            } while (edges[edgeIndices[nextIdx]] == null && nextIdx !== idx);
            do {
                prevIdx = prevIdx - 1;
                if (prevIdx < 0) {
                    prevIdx = edgesCount - 1;
                }
            } while (edges[edgeIndices[prevIdx]] == null && prevIdx !== idx);

            if (prevIdx === idx || nextIdx === idx || prevIdx === nextIdx) {
                break;
            }

            const next = edges[edgeIndices[nextIdx]];
            const prev = edges[edgeIndices[prevIdx]];

            lineIntersection(
                prev[0], prev[1], prev[2], prev[3],
                curr[0], curr[1], curr[2], curr[3], p1, 0
            );
            lineIntersection(
                curr[0], curr[1], curr[2], curr[3],
                next[0], next[1], next[2], next[3], p2, 0
            );

            p0[0] = curr[0];
            p0[1] = curr[1];
            p3[0] = curr[2];
            p3[1] = curr[3];

            // Always keep first and last edge if not close
            if (!close && (idx === 0 || idx === edgesCount - 1)) {
                if (i === 0) {
                    // Keep the start point.
                    edges[edgeIndices[idx]] = [p0[0], p0[1], p2[0], p2[1]];
                }
                else {
                    // Keep the end point.
                    edges[edgeIndices[idx]] = [p1[0], p1[1], p3[0], p3[1]];
                }
            }
            else {
                // Check if two intersection points are crossed.
                const t0 = (p1[0] - p0[0]) / (p3[0] - p0[0]);
                const t1 = (p2[0] - p3[0]) / (p0[0] - p3[0]);

                if (
                    // On same side
                    t0 > 1 && t1 < 0
                    || t0 < 0 && t1 > 1
                    // Crossed
                    || t0 + t1 > 1
                ) {
                    hasRemoved = true;
                    edges[edgeIndices[idx]] = null;
                    // Empty edge.
                    continue;
                }

                edges[edgeIndices[idx]] = [p1[0], p1[1], p2[0], p2[1]];
            }
        }

        edgeIndices = [];
        for (let i = 0; i < edges.length; i++) {
            if (edges[i]) {
                edgeIndices.push(i);
            }
        };
    }

    do {
        removeEdgeIntersections();
    } while (edgeIndices.length > 3 && hasRemoved);

    let idx = edges.length - 1;
    let prevEdge;

    if (close) {
        do {
            prevEdge = edges[idx--];
        } while (!prevEdge && idx >= 0)

        if (!prevEdge) {
            // Degenerate
            return;
        }
    }

    let writeIdx = outStart;
    const isRoundJoin = join === 'round';
    const isBevelJoin = join === 'bevel';
    const isMiterJoin = !join || join === 'miter';
    // Temp variables.
    const pt1 = [];
    const pt2 = [];
    const pt3 = [];
    const pt4 = [];
    const v = [];
    const v2 = [];

    for (let i = 0; i < edges.length; i++) {
        const edge = edges[i] || prevEdge;
        const isEnd = edge === prevEdge ? 2 : 0;

        needsIndicesMap && (indicesMap[i] = writeIdx);

        if (isComplexJoin && edge !== prevEdge && i > 0) {
            // Offset vertex
            pt2[0] = prevEdge[2];
            pt2[1] = prevEdge[3];
            // Raw vertex
            pt1[0] = vertices[(start + i) * 2];
            pt1[1] = vertices[(start + i) * 2 + 1];
            v2[0] = edge[2] - edge[0];
            v2[1] = edge[3] - edge[1];
            v2Sub(v, pt2, pt1);
            const isOuterior = v2Dot(v, v2) < 0;

            const d = v2Dist(pt1, pt2);
            if (isOuterior) {
                if (isRoundJoin) {
                    pt3[0] = offsetedEdges[i - 1][2];
                    pt3[1] = offsetedEdges[i - 1][3];

                    pt4[0] = offsetedEdges[i][0];
                    pt4[1] = offsetedEdges[i][1];
                    // TODO not fixed segments
                    innerAppendArc(out, pt3, pt4, pt1, 5, Math.abs(offset));
                    writeIdx += 5;
                }
                else {
                    const limit = isMiterJoin ? miterLimit : Math.abs(offset);
                    // Bevel join or miter limit join
                    if (d > limit) {
                        v2Normalize(v, v);
                        // Calculate the bevel line
                        pt1[0] = pt1[0] + v[0] * limit;
                        pt1[1] = pt1[1] + v[1] * limit;
                        pt2[0] = pt1[0] - v[1];
                        pt2[1] = pt1[1] + v[0];

                        lineIntersection(
                            prevEdge[0], prevEdge[1], prevEdge[2], prevEdge[3],
                            pt1[0], pt1[1], pt2[0], pt2[1], pt3, 0
                        );
                        lineIntersection(
                            pt1[0], pt1[1], pt2[0], pt2[1],
                            edge[0], edge[1], edge[2], edge[3], pt4, 0
                        );

                        out[writeIdx * 2] = pt3[0];
                        out[writeIdx * 2 + 1] = pt3[1];
                        writeIdx++;
                        out[writeIdx * 2] = pt4[0];
                        out[writeIdx * 2 + 1] = pt4[1];
                        writeIdx++;
                    }
                    else {
                        out[writeIdx * 2] = edge[isEnd];
                        out[writeIdx * 2 + 1] = edge[isEnd + 1];
                        writeIdx++;
                    }
                }
            }
            else {
                out[writeIdx * 2] = edge[isEnd];
                out[writeIdx * 2 + 1] = edge[isEnd + 1];
                writeIdx++;
            }
        }
        else {
            out[writeIdx * 2] = edge[isEnd];
            out[writeIdx * 2 + 1] = edge[isEnd + 1];
            writeIdx++;
        }

        prevEdge = edge;
    }

    if (!close) {
        const len = edges.length;
        const lastEdge = edges[len - 1] || prevEdge;
        out[writeIdx * 2] = lastEdge[2];
        out[writeIdx * 2 + 1] = lastEdge[3];

        needsIndicesMap && (indicesMap[len] = writeIdx);
    }

    // if (typeof ctx !== 'undefined') {
    //     offsetedEdges.forEach((edge, idx) => {
    //         ctx.beginPath();
    //         ctx.strokeStyle = '#0f0';
    //         ctx.globalAlpha = idx / 100 + 0.1;
    //         ctx.moveTo(edge[0], edge[1]);
    //         ctx.lineTo(edge[2], edge[3]);
    //         ctx.stroke();
    //     })
    // }

    return indicesMap;
}

export function offsetPolyline(vertices, offset, join, miterLimit) {
    const offsetVertices = [];

    innerOffsetContour(
        vertices, offsetVertices, 0, vertices.length / 2, 0, offset, join, miterLimit, false
    );

    return offsetVertices;
}

export function offsetPolygon(vertices, holes, offset) {
    const offsetVertices = new Float32Array(vertices.length);
    const exteriorSize = (holes && holes.length) ? holes[0] : vertices.length / 2;

    innerOffsetContour(
        vertices, offsetVertices, 0, exteriorSize, 0, offset, null, null, true
    );

    if (holes) {
        for (let i = 0; i < holes.length; i++) {
            const start = holes[i];
            const end = holes[i + 1] || vertices.length / 2;
            innerOffsetContour(
                vertices, offsetVertices, start, end,
                start, offset, null, null, close
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
    if (area(vertices, start, end) > 0) {
        reversePoints(vertices, 2, start, end);
    }
    for (let h = 1; h < (holes ? holes.length : 0) + 1; h++) {
        start = holes[h - 1];
        end = holes[h] || polygonVertexCount;
        if (area(vertices, start, end) < 0) {
            reversePoints(vertices, 2, start, end);
        }
    }
}

function normalizeOpts(opts) {

    opts.depth = opts.depth || 1;
    opts.bevelSize = opts.bevelSize || 0;
    opts.bevelSegments = opts.bevelSegments == null ? 2 : opts.bevelSegments;
    opts.smoothBevel = opts.smoothBevel || false;
    opts.simplify = opts.simplify || 0;

    if (opts.smoothSide == null) {
        opts.smoothSide = 'auto'
    }
    if (opts.smoothSideThreshold == null) {
        opts.smoothSideThreshold = 0.9
    }

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
    out, {vertices, topVertices, splittedMap, depth, rect}, start, end,
    cursors, opts
) {
    const ringVertexCount = end - start;

    const splitBevel = opts.smoothBevel ? 1 : 2;
    const bevelSize = Math.min(depth / 2, opts.bevelSize);
    const bevelSegments = opts.bevelSegments;
    const vertexOffset = cursors.vertex;
    const size = Math.max(rect.width, rect.height, depth);

    function isDuplicateVertex(idx) {
        const nextIdx = (idx + 1) % ringVertexCount;
        const x0 = vertices[idx * 2];
        const y0 = vertices[idx * 2 + 1];
        const x1 = vertices[nextIdx * 2];
        const y1 = vertices[nextIdx * 2 + 1];
        return x0 === x1 && y0 === y1;
    }

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
                    const idx = (i % ringVertexCount + start) * 2;
                    const rawIdx = splittedMap ? splittedMap[idx / 2] * 2 : idx;
                    v1[0] = vertices[idx] - topVertices[rawIdx];
                    v1[1] = vertices[idx + 1] - topVertices[rawIdx + 1];
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

                    const x = v[0] * r + topVertices[rawIdx];
                    const y = v[1] * r + topVertices[rawIdx + 1];
                    const zz = v[2] * r + z;
                    out.position[cursors.vertex * 3] = x;
                    out.position[cursors.vertex * 3 + 1] = y;
                    out.position[cursors.vertex * 3 + 2] = zz;

                    // TODO Cache and optimize
                    if (i > 0) {
                        uLen += Math.sqrt((prevX - x) * (prevX - x) + (prevY - y) * (prevY - y));
                    }
                    if (s > 0 || k > 0) {
                        let tmp = (cursors.vertex - ringVertexCount) * 3;
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

                    // Just ignore this face if vertex are duplicted in `splitVertices`
                    if (isDuplicateVertex(i)) {
                        continue;
                    }
                    if ((splitBevel > 1 && (s % splitBevel)) || (splitBevel === 1 && s >= 1)) {
                        for (let f = 0; f < 6; f++) {
                            const m = (quadToTriangle[f][0] + i) % ringVertexCount;
                            const n = quadToTriangle[f][1] + ringCount;
                            out.indices[cursors.index++] = (n - 1) * ringVertexCount + m + vertexOffset;
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
                const idx = (i % ringVertexCount + start) * 2;
                const x = vertices[idx];
                const y = vertices[idx + 1];
                out.position[cursors.vertex * 3] = x;
                out.position[cursors.vertex * 3 + 1] = y;
                out.position[cursors.vertex * 3 + 2] = z;
                if (i > 0) {
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
    // Connect the side
    const sideStartRingN = bevelSize > 0 ? (bevelSegments * splitBevel + 1) : 1;
    for (let i = 0; i < ringVertexCount; i++) {
        // Just ignore this face if vertex are duplicted in `splitVertices`
        if (isDuplicateVertex(i)) {
            continue;
        }
        for (let f = 0; f < 6; f++) {
            const m = (quadToTriangle[f][0] + i) % ringVertexCount;
            const n = quadToTriangle[f][1] + sideStartRingN;
            out.indices[cursors.index++] = (n - 1) * ringVertexCount + m + vertexOffset;
        }
    }
}

function addTopAndBottom({indices, topVertices, rect, depth}, out, cursors, opts) {
    if (topVertices.length <= 4) {
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
        const vertexCount = topVertices.length / 2;
        for (let i = 0; i < indicesLen; i += 3) {
            for (let k = 0; k < 3; k++) {
                out.indices[cursors.index++] = vertexOffset + vertexCount + indices[i + 2 - k];
            }
        }
    }
}

/**
 * Split vertices for sharp side.
 */
 function splitVertices(vertices, holes, smoothSide, smoothSideThreshold) {
    const isAutoSmooth = smoothSide == null || smoothSide === 'auto';
    if (smoothSide === true) {
        return {vertices, holes};
    }
    const newVertices = [];
    const newHoles = holes && [];
    const count = vertices.length / 2;
    const v1 = [];
    const v2 = [];

    // Map of splitted index to raw index
    const splittedMap = [];

    let start = 0;
    let end = 0;

    const polysCount = (holes ? holes.length : 0) + 1;
    for (let h = 0; h < polysCount; h++) {
        if (h === 0) {
            end = holes && holes.length ? holes[0] : count;
        }
        else {
            start = holes[h - 1];
            end = holes[h] || count;
        }

        for (let i = start; i < end; i++) {
            const x2 = vertices[i * 2];
            const y2 = vertices[i * 2 + 1];
            const nextIdx = i === end - 1 ? start : i + 1;
            const x3 = vertices[nextIdx * 2];
            const y3 = vertices[nextIdx * 2 + 1];

            if (isAutoSmooth) {
                const prevIdx = i === start ? end - 1 : i - 1;
                const x1 = vertices[prevIdx * 2];
                const y1 = vertices[prevIdx * 2 + 1];

                v1[0] = x1 - x2;
                v1[1] = y1 - y2;
                v2[0] = x3 - x2;
                v2[1] = y3 - y2;

                v2Normalize(v1, v1);
                v2Normalize(v2, v2);

                const angleCos = v2Dot(v1, v2) * 0.5 + 0.5;

                if ((1 - angleCos) > smoothSideThreshold) {
                    newVertices.push(x2, y2);
                    splittedMap.push(i);
                }
                else {
                    newVertices.push(x2, y2, x2, y2);
                    splittedMap.push(i, i);
                }
            }
            else {
                newVertices.push(x2, y2, x2, y2);
                splittedMap.push(i, i);
            }
        }

        if (h < polysCount - 1 && newHoles) {
            newHoles.push(newVertices.length / 2);
        }
    }

    return {
        vertices: new Float32Array(newVertices),
        splittedMap,
        holes: newHoles
    };
}

function innerExtrudeTriangulatedPolygon(preparedData, opts) {
    let indexCount = 0;
    let vertexCount = 0;

    for (let p = 0; p < preparedData.length; p++) {
        const {indices, vertices, splittedMap, topVertices, holes, depth} = preparedData[p];
        const bevelSize = Math.min(depth / 2, opts.bevelSize);
        const bevelSegments = !(bevelSize > 0) ? 0 : opts.bevelSegments;

        indexCount += indices.length * (opts.excludeBottom ? 1 : 2);
        vertexCount += topVertices.length / 2 * (opts.excludeBottom ? 1 : 2);
        const ringCount = 2 + bevelSegments * 2;

        let start = 0;
        let end = 0;
        for (let h = 0; h < (holes ? holes.length : 0) + 1; h++) {
            if (h === 0) {
                end = (holes && holes.length) ? holes[0] : vertices.length / 2;
            }
            else {
                start = holes[h - 1];
                end = holes[h] || vertices.length / 2;
            }

            const faceEnd = splittedMap ? splittedMap[end - 1] + 1 : end;
            const faceStart = splittedMap ? splittedMap[start] : start;
            indexCount += (faceEnd - faceStart) * 6 * (ringCount - 1);

            const sideRingVertexCount = end - start;
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
        const vertexCount = vertices.length / 2;

        let start = 0;
        let end = (holes && holes.length) ? holes[0] : vertexCount;
        // Add exterior
        addExtrudeSide(data, preparedData[d], start, end, cursors, opts);
        // Add holes
        if (holes) {
            for (let h = 0; h < holes.length; h++) {
                start = holes[h];
                end = holes[h + 1] || vertexCount;
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
    const scale = opts.scale || [1, 1];
    for (let i = 0, k = 0; i < pointCount; i++) {
        points[k++] = polyline[i][0] * scale[0] + translate[0];
        points[k++] = polyline[i][1] * scale[1] + translate[1];
    }

    if (area(points, 0, pointCount) < 0) {
        reversePoints(points, 2, 0, pointCount);
    }

    const insidePoints = [];
    const outsidePoints = [];
    const miterLimit = opts.miterLimit;
    const lineJoin = opts.lineJoin;
    const outsideIndicesMap = innerOffsetContour(
        points, outsidePoints, 0, pointCount, 0, -lineWidth / 2, lineJoin, miterLimit, false
    );
    reversePoints(points, 2, 0, pointCount);
    const insideIndicesMap = innerOffsetContour(
        points, insidePoints, 0, pointCount, 0, -lineWidth / 2, lineJoin, miterLimit, false
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

        if (insideIndicesMap[i2] - insideIndicesMap[i] > 1) {
            for (let k = 1; k < insideIndicesMap[i2] - insideIndicesMap[i]; k++) {
                indices[off++] = insideIndicesMap[i] + k + 1 + outsidePointCount;
                indices[off++] = insideIndicesMap[i] + k + outsidePointCount;
                indices[off++] = outsidePointCount - outsideIndicesMap[i2] - 1;
            }
        }
        else if (outsideIndicesMap[i2] - outsideIndicesMap[i] > 1) {
            for (let k = 1; k < outsideIndicesMap[i2] - outsideIndicesMap[i]; k++) {
                indices[off++] = insideIndicesMap[i2] + outsidePointCount;
                indices[off++] = outsidePointCount - 1 - (outsideIndicesMap[i] + k);
                indices[off++] = outsidePointCount - 1 - (outsideIndicesMap[i] + k + 1);
            }
        }
    }

    const topVertices = opts.bevelSize > 0
        ? offsetPolygon(polygonVertices, [], opts.bevelSize, null, true) : polygonVertices;
    const boundingRect = opts.boundingRect;

    const res = splitVertices(polygonVertices, null, opts.smoothSide, opts.smoothSideThreshold);
    return {
        vertices: res.vertices,
        rawVertices: polygonVertices,
        splittedMap: res.splittedMap,
        indices,
        topVertices,
        rect: {
            x: boundingRect.x * scale[0] + translate[0],
            y: boundingRect.y * scale[1] + translate[1],
            width: boundingRect.width * scale[0],
            height: boundingRect.height * scale[1],
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
        points = doSimplify(points, tolerance, true);
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
 * @param {boolean} [opts.smoothSide = 'auto']
 * @param {boolean} [opts.smoothSideThreshold = 0.9]    // Will not smooth sharp side.
 * @param {boolean} [opts.smoothBevel = false]
 * @param {boolean} [opts.excludeBottom = false]
 * @param {Object} [opts.fitRect] translate and scale will be ignored if fitRect is set
 * @param {Array} [opts.translate]
 * @param {Array} [opts.scale]
 *
 * @return {Object} {indices, position, uv, normal, boundingRect}
 */
export function extrudePolygon(polygons, opts) {

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
    const scale = opts.scale || [1, 1];
    const boundingRect = opts.boundingRect;
    const transformdRect = {
        x: boundingRect.x * scale[0] + translate[0],
        y: boundingRect.y * scale[1] + translate[1],
        width: boundingRect.width * scale[0],
        height: boundingRect.height * scale[1],
    };

    const epsilon = Math.min(
        boundingRect.width, boundingRect.height
    ) / 1e5;
    for (let i = 0; i < polygons.length; i++) {
        let newPolygon = removeClosePointsOfPolygon(polygons[i], epsilon);
        if (!newPolygon) {
            continue;
        }
        const simplifyTolerance = opts.simplify / Math.max(scale[0], scale[1]);
        if (simplifyTolerance > 0) {
            newPolygon = simplifyPolygon(newPolygon, simplifyTolerance);
        }
        if (!newPolygon) {
            continue;
        }

        const {vertices, holes, dimensions} = earcut.flatten(newPolygon);

        for (let k = 0; k < vertices.length;) {
            vertices[k] = vertices[k++] * scale[0] + translate[0];
            vertices[k] = vertices[k++] * scale[1] + translate[1];
        }

        convertToClockwise(vertices, holes);

        if (dimensions !== 2) {
            throw new Error('Only 2D polygon points are supported');
        }
        const topVertices = opts.bevelSize > 0
            ? offsetPolygon(vertices, holes, opts.bevelSize, null, true) : vertices;
        const indices = triangulate(topVertices, holes, dimensions);
        const res = splitVertices(vertices, holes, opts.smoothSide, opts.smoothSideThreshold)

        preparedData.push({
            indices,
            vertices: res.vertices,
            rawVertices: vertices,
            topVertices,
            holes: res.holes,
            splittedMap: res.splittedMap,
            rect: transformdRect,
            depth: typeof opts.depth === 'function' ? opts.depth(i) : opts.depth
        });
    }
    return innerExtrudeTriangulatedPolygon(preparedData, opts);
};

/**
 *
 * @param {Array} polylines Polylines array that match GeoJSON MultiLineString geometry.
 * @param {Object} [opts]
 * @param {number} [opts.depth]
 * @param {number} [opts.bevelSize = 0]
 * @param {number} [opts.bevelSegments = 2]
 * @param {number} [opts.simplify = 0]
 * @param {boolean} [opts.smoothSide = 'auto']
 * @param {boolean} [opts.smoothSideThreshold = 0.9]    // Will not smooth sharp side.
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
export function extrudePolyline(polylines, opts) {

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
    const scale = opts.scale || [1, 1];

    if (opts.lineWidth == null) {
        opts.lineWidth = 1;
    }
    if (opts.miterLimit == null) {
        opts.miterLimit = 2;
    }
    opts.bevelSize = Math.min(opts.bevelSize, opts.lineWidth / 2 - 1e-5);   // Add a tiny threshold.

    const preparedData = [];
    // Extrude polyline to polygon
    for (let i = 0; i < polylines.length; i++) {
        let newPolyline = polylines[i];
        const simplifyTolerance = opts.simplify / Math.max(scale[0], scale[1]);
        if (simplifyTolerance > 0) {
            newPolyline = doSimplify(newPolyline, simplifyTolerance, true);
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
 * @param {boolean} [opts.smoothSide = 'auto']
 * @param {boolean} [opts.smoothSideThreshold = 0.9]    // Will not smooth sharp side.
 * @param {boolean} [opts.smoothBevel = false]
 * @param {boolean} [opts.excludeBottom = false]
 * @param {boolean} [opts.lineWidth = 1]
 * @param {boolean} [opts.lineWidth = 'miter']
 * @param {boolean} [opts.miterLimit = 2]
 * @param {Object} [opts.fitRect] translate and scale will be ignored if fitRect is set
 * @param {Array} [opts.translate]
 * @param {Array} [opts.scale]
 * @param {Object} [opts.boundingRect]
 * @return {Object} {polyline: {indices, position, uv, normal}, polygon: {indices, position, uv, normal}}
 */

 // TODO Not merge feature
export function extrudeGeoJSON(geojson, opts) {

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