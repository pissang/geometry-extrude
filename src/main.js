import earcut from 'earcut';
import {slerp, scale, normalize, lineIntersection, area, triangleArea} from './math';

export function triangulate(vertices, holes, dimensions=2) {
    return earcut(vertices, holes, dimensions);
};

export function flatten(data) {
    return earcut.flatten(data);
}

function offsetLine(vertices, out, start, offset) {
    const i = start;
    const j = start + 1;
    const x1 = vertices[i * 2];
    const y1 = vertices[i * 2 + 1];
    const x2 = vertices[j * 2];
    const y2 = vertices[j * 2 + 1];
    let dx = y2 - y1;
    let dy = x1 - x2;
    const l = Math.sqrt(dx * dx + dy * dy) / offset;

    dx /= l;
    dy /= l;

    out[i * 2] = x1 + dx;
    out[i * 2 + 1] = y1 + dy;
    out[j * 2] = x2 + dx;
    out[j * 2 + 1] = y2 + dy;
}

export function offsetPolygon(vertices, out, start, end, outStart, offset, miterLimit) {
    const checkMiterLimit = miterLimit != null;
    let outOff = outStart;
    for (let i = start + 1, j = start, k = end - 1; j < end;) {
        const nextIdx = i === end ? start : i;
        const x1 = vertices[k * 2];
        const y1 = vertices[k * 2 + 1];
        const x2 = vertices[j * 2];
        const y2 = vertices[j * 2 + 1];
        const x3 = vertices[nextIdx * 2];
        const y3 = vertices[nextIdx * 2 + 1];

        let dx1 = y2 - y1;
        let dy1 = x1 - x2;

        let dx2 = y3 - y2;
        let dy2 = x2 - x3;

        const l1 = Math.sqrt(dx1 * dx1 + dy1 * dy1) / offset;
        const l2 = Math.sqrt(dx2 * dx2 + dy2 * dy2) / offset;

        dx1 /= l1;
        dy1 /= l1;
        dx2 /= l2;
        dy2 /= l2;

        lineIntersection(
            x1 + dx1, y1 + dy1, x2 + dx1, y2 + dy1,
            x2 + dx2, y2 + dy2, x3 + dx2, y3 + dy2,
            out, outOff * 2
        );

        if (checkMiterLimit) {
            const isConvex = triangleArea(x1, y1, x2, y2, x3, y3) > 0;
            const x = out[outOff * 2];
            const y = out[outOff * 2 + 1];
            const dx = x - x2;
            const dy = y - y2;
            const miter = Math.sqrt(dx * dx + dy * dy);
            const miterInUnit = Math.abs(miter / offset);
            if (miterInUnit > miterLimit && isConvex) {
                const mx = x2 + dx / miterInUnit;
                const my = y2 + dy / miterInUnit;

                lineIntersection(
                    x1 + dx1, y1 + dy1, x2 + dx1, y2 + dy1,
                    mx, my, mx - dy, my + dx,
                    out, outOff * 2
                );
                lineIntersection(
                    x2 + dx2, y2 + dy2, x3 + dx2, y3 + dy2,
                    mx, my, mx - dy, my + dx,
                    out, (outOff + 1) * 2
                );
                outOff++;
            }
        }

        outOff++;

        k = j;
        j = i;
        i++;
    }
}

export function offsetPolygonWithHole(vertices, holes, offset) {
    const offsetVertices = new Float32Array(vertices.length);
    const exteriorSize = (holes && holes.length) ? holes[0] : vertices.length / 2;

    offsetPolygon(vertices, offsetVertices, 0, exteriorSize, 0, offset);

    if (holes) {
        for (let i = 0; i < holes.length; i++) {
            const start = holes[i];
            const end = holes[i + 1] || vertices.length / 2;
            offsetPolygon(vertices, offsetVertices, start, end, start, offset);
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

// 0,0----1,0
// 0,1----1,1
const quadToTriangle = [
    [0, 0], [1, 0], [1, 1],
    [0, 0], [1, 1], [0, 1]
];

// Add side vertices and indices. Include bevel.
function addExtrudeSide(
    out, vertices, topVertices, start, end,
    cursors, opts
) {
    const depth = opts.depth;
    const ringVertexCount = end - start;
    const splitSide = opts.smoothSide ? 1 : 2;
    const splitRingVertexCount = ringVertexCount * splitSide;

    const splitBevel = opts.smoothBevel ? 1 : 2;
    const bevelSize = opts.bevelSize;
    const bevelSegments = opts.bevelSegments;
    const vertexOffset = cursors.vertex;
    // Side vertices
    if (bevelSize > 0) {

        const v0 = [0, 0, 1];
        const v1 = [];
        const v2 = [0, 0, -1];
        const v = [];

        let ringCount = 0;
        for (let k = 0; k < 2; k++) {
            const z = (k === 0 ? (depth - bevelSize) : bevelSize);
            for (let s = 0; s <= bevelSegments * splitBevel; s++) {
                for (let i = 0; i < ringVertexCount; i++) {

                    for (let j = 0; j < splitSide; j++) {
                        // TODO Cache and optimize
                        let idx = ((i + j) % ringVertexCount + start) * 2;
                        v1[0] = vertices[idx] - topVertices[idx];
                        v1[1] = vertices[idx + 1] - topVertices[idx + 1];
                        v1[2] = 0;
                        normalize(v1, v1);

                        const t = (Math.floor(s / splitBevel) + (s % splitBevel)) / bevelSegments;
                        k === 0 ? slerp(v, v0, v1, t)
                            : slerp(v, v1, v2, t);

                        out.position[cursors.vertex * 3] = v[0] * bevelSize + topVertices[idx];
                        out.position[cursors.vertex * 3 + 1] = v[1] * bevelSize + topVertices[idx + 1];
                        out.position[cursors.vertex * 3 + 2] = v[2] * bevelSize + z;
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
            for (let i = 0; i < ringVertexCount; i++) {
                for (let m = 0; m < splitSide; m++) {
                    const idx = ((i + m) % ringVertexCount + start) * 2;
                    out.position[cursors.vertex * 3] = vertices[idx];
                    out.position[cursors.vertex * 3 + 1] = vertices[idx + 1];
                    out.position[cursors.vertex * 3 + 2] = z;
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

function addTopAndBottom({indices, vertices, topVertices}, out, cursors, opts) {
    const depth = opts.depth;
    if (vertices.length <= 2) {
        return;
    }

    const vertexOffset = cursors.vertex;
    // Top indices
    const indicesLen = indices.length;
    for (let i = 0; i < indicesLen; i++) {
        out.indices[cursors.index++] = vertexOffset + indices[i];
    }
    // Top vertices
    for (let i = 0; i < topVertices.length; i += 2) {
        out.position[cursors.vertex * 3] = topVertices[i];
        out.position[cursors.vertex * 3 + 1] = topVertices[i + 1];
        out.position[cursors.vertex * 3 + 2] = depth;
        cursors.vertex++;
    }

    // Bottom indices
    for (let i = 0; i < indicesLen; i += 3) {
        for (let k = 0; k < 3; k++) {
            out.indices[cursors.index++] = cursors.vertex + indices[i + 2 - k];
        }
    }
    // Bottom vertices
    for (let i = 0; i < topVertices.length; i += 2) {
        out.position[cursors.vertex * 3] = topVertices[i];
        out.position[cursors.vertex * 3 + 1] = topVertices[i + 1];
        out.position[cursors.vertex * 3 + 2] = 0;
        cursors.vertex++;
    }
}

function normalizeOpts(opts) {

    opts.depth = opts.depth || 1;
    opts.bevelSize = opts.bevelSize || 0;
    opts.bevelSegments = opts.bevelSegments == null ? 2 : opts.bevelSegments;
    opts.smoothSide = opts.smoothSide || false;
    opts.smoothBevel = opts.smoothBevel || false;

    // Normalize bevel options.
    opts.bevelSize = Math.min(!(opts.bevelSegments > 0) ? 0 : opts.bevelSize, opts.depth / 2);
    if (!(opts.bevelSize > 0)) {
        opts.bevelSegments = 0;
    }
    opts.bevelSegments = Math.round(opts.bevelSegments);
}

function convertToAnticlockwise(vertices, holes) {
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


function extrudeFlattenPolygon(flaternPolygons, opts) {
    const preparedData = [];
    let indexCount = 0;
    let vertexCount = 0;
    for (let p = 0; p < flaternPolygons.length; p++) {
        const {vertices, holes, dimensions} = flaternPolygons[p];

        convertToAnticlockwise(vertices, holes);

        if (dimensions !== 2) {
            throw new Error('Only 2D polygon points are supported');
        }
        let topVertices = vertices;
        if (opts.bevelSize > 0) {
            topVertices = offsetPolygonWithHole(vertices, holes, opts.bevelSize);
        }
        const indices = triangulate(topVertices, holes, dimensions);
        const polygonVertexCount = vertices.length / 2;
        preparedData.push({
            indices,
            vertices,
            topVertices,
            holes
        });
        indexCount += indices.length * 2;
        vertexCount += polygonVertexCount * 2;
        const ringCount = 2 + opts.bevelSegments * 2;

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
                + (!opts.smoothBevel ? opts.bevelSegments * sideRingVertexCount * 2 : 0);
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
        const {holes, vertices, topVertices} = preparedData[d];
        const topVertexCount = vertices.length / 2;

        let start = 0;
        let end = (holes && holes.length) ? holes[0] : topVertexCount;
        // Add exterior
        addExtrudeSide(data, vertices, topVertices, start, end, cursors, opts);
        // Add holes
        if (holes) {
            for (let h = 0; h < holes.length; h++) {
                start = holes[h];
                end = holes[h + 1] || topVertexCount;
                addExtrudeSide(data, vertices, topVertices, start, end, cursors, opts);
            }
        }

    }

    return data;
}
/**
 *
 * @param {Array} polygons Polygons array that match GeoJSON MultiPolygon geometry.
 * @param {Object} [opts]
 * @param {number} [opts.depth]
 * @param {number} [opts.bevelSize = 0]
 * @param {number} [opts.bevelSegments = 2]
 * @param {boolean} [opts.smoothSide = false]
 * @param {boolean} [opts.smoothBevel = false]
 */
// TODO Dimensions
// TODO UV, normal
// TODO If smooth connection between side and bevel.
// TODO anticlockwise
// TODO Ignore bottom, bevel="top"|"bottom"
export function extrudePolygon(polygons, opts) {

    opts = opts || {};
    normalizeOpts(opts);

    const flattenPolygons = [];
    for (let i = 0; i < polygons.length; i++) {
        flattenPolygons.push(earcut.flatten(polygons[i]));
    }
    return extrudeFlattenPolygon(flattenPolygons, opts);
};

function convertPolylineToFlattenPolygon(polyline, lineWidth) {
    const pointCount = polyline.length;
    const points = new Float32Array(pointCount * 2);
    for (let i = 0, k = 0; i < pointCount; i++) {
        points[k++] = polyline[i][0];
        points[k++] = polyline[i][1];
    }

    const polygon = new Float32Array(points.length * 2);
    if (pointCount === 2) {
        offsetLine(points, polygon, 0, lineWidth / 2);
    }
    else {
        offsetPolygon(points, polygon, 0, pointCount, lineWidth / 2);
    }
    for (let i = 0; i < pointCount; i++) {
        polygon[(pointCount * 2 - i - 1) * 2] = polygon[i * 2];
        polygon[(pointCount * 2 - i - 1) * 2 + 1] = polygon[i * 2 + 1];
    }
    if (pointCount === 2) {
        offsetLine(points, polygon, 0, -lineWidth / 2);
    }
    else {
        offsetPolygon(points, polygon, 0, pointCount, -lineWidth / 2);
    }

    console.log(polygon);

    return {
        vertices: polygon,
        holes: [],
        dimensions: 2
    };
}

/**
 *
 * @param {Array} polylines Polylines array that match GeoJSON MultiLineString geometry.
 * @param {Object} [opts]
 * @param {number} [opts.depth]
 * @param {number} [opts.bevelSize = 0]
 * @param {number} [opts.bevelSegments = 2]
 * @param {boolean} [opts.smoothSide = false]
 * @param {boolean} [opts.smoothBevel = false]
 * @param {boolean} [opts.lineWidth = 1]
 */
export function extrudePolyline(polylines, opts) {
    normalizeOpts(opts);
    if (opts.lineWidth == null) {
        opts.lineWidth = 1;
    }
    const polygons = [];
    // Extrude polyline to polygon
    for (let i = 0; i < polylines.length; i++) {
        polygons.push(convertPolylineToFlattenPolygon(polylines[i], opts.lineWidth));
    }

    return extrudeFlattenPolygon(polygons, opts);
}