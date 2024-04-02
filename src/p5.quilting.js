/**
 * (c) JohnC 2024
 */
import { Delaunay } from "d3-delaunay";
import { kdTree } from "kd-tree-javascript";
const KdTree = kdTree;
/**
 * 
 * @param {number} x center x
 * @param {number} y center y
 * @param {number} r radius
 * @returns a function to check if a point is with in the circle
 */
export function circleShape(x, y, r) {
    return (p) => {
        const dx = p[0] - x, dy = p[1] - y;
        return dx * dx + dy * dy < r * r;
    }
}

/**
 * 
 * @param {number} x center x
 * @param {number} y center y
 * @param {number} w width of the rectangle
 * @param {number} h height of the rectangle
 * @returns a function checking if a point is inside a rectangle
 */
export function rectShape(x, y, w, h) {
    return (p) => {
        const dx = p[0] - x, dy = p[1] - y;
        return Math.abs(dx) < w / 2 && Math.abs(dy) < h / 2;
    }
}

/**
 * 
 * @param {number} x1 
 * @param {number} y1 
 * @param {number} x2 
 * @param {number} y2 
 * @param {number} x3 
 * @param {number} y3 
 * @returns a function to check if a point is inside a triangle
 */
export function triShape(x1, y1, x2, y2, x3, y3) {
    return (p) => {
        const sign = (p1x, p1y, p2x, p2y, p3x, p3y) => {
            return (p1x - p3x) * (p2y - p3y) - (p2x - p3x) * (p1y - p3y);
        }

        let d1, d2, d3;
        let b1, b2;
        d1 = sign(p[0], p[1], x1, y1, x2, y2);
        d2 = sign(p[0], p[1], x2, y2, x3, y3);
        d3 = sign(p[0], p[1], x1, y1, x3, y3);
        b1 = (d1 < 0) || (d2 < 0) || (d3 < 0);
        b2 = (d1 > 0) || (d2 > 0) || (d3 > 0);
        return !(has_neg && has_pos);
    }
}

/**
 * 
 * @param {Array} points points of the close polygon in order
 */
export function polygonShape(points){
    return (p) => {
        const ps = points.slice();
        if (ps[0] !== ps[ps.length - 2] || ps[1] !== ps[ps.length - 1]) {
            // 
            ps.push(ps[0], ps[1]);
        }

        let x = p[0], y = p[1];
        let intersection = 0;
        for (let i = 0; i < ps.length - 2; i += 2){
            let x1 = ps[i], y1 = ps[i + 1], x2 = ps[i + 2], y2 = ps[i + 3];
            if (((y < y1) !== (y < y2)) && ((x < x1) !== (x < x2)) && (x2 - x1) * (y - y1) === (x - x1) * (y2 - y1)) return true;
            if (((y < y1) !== (y < y2)) && x < (x2 - x1) * (y - y1) / (y2 - y1) + x1) intersection++;
        }
        if (intersection % 2 === 1) return true;
        return false;
    }
}

/**
 * 
 * @param {Array} center [x,y] center (starting point) of the point distribution
 * @param {Function} shapeFunction a function to check if a point is with in the shape, it should accept [x,y] as input and return a boolean value
 * @param {Array} cellSize a or [w, h], one or two number denoting the cellSize of the grid 
 * @param {number} maxDist stop after the distribution reach this size
 * @param {number} maxPointNo stop after this number of points are placed
 * @returns [x1, y1, x2, y2, ...] points distributed in the shape
 */
export function gridDistribution(center, shapeFunction, cellSize, maxDist, maxPointNo) {
    const allowedPointNo = maxPointNo || Infinity;
    const allowedSize = maxDist || Infinity;
    if (allowedPointNo === Infinity && allowedSize === Infinity) throw new Error("one of the end condition must be provided");
    let res = [];
    const cellW = typeof cellSize === 'number' ? cellSize : cellSize[0];
    const cellH = typeof cellSize === 'number' ? cellSize : cellSize[1];
    const cx = center[0], cy = center[1];
    let n = 0;
    let fail = 0;
    while (res.length / 2 < allowedPointNo) {
        let offset = ulam(n);
        let p = [cx + offset[0] * cellW, cy + offset[1] * cellH];
        if (shapeFunction(p)) {
            res.push(...p);
            fail = 0;
        } else {
            fail++;
        }
        n++;

        // check maxDist
        let d = Math.min(cellW * Math.abs(offset[0]), cellH * Math.abs(offset[1]));
        if (d > allowedSize) break;
        if (fail > 999) break;
    }
    // 
    return res;
}

/**
 * 
 * @param {Array} center [x,y] center (starting point) of the point distribution
 * @param {Function} shapeFunction a function to check if a point is with in the shape, it should accept [x,y] as input and return a boolean value
 * @param {number} a the size of the tessellation unit
 * @param {number} maxDist stop after the distribution reach this size
 * @param {number} maxPointNo stop after this number of points are placed
 * @returns [x1, y1, x2, y2, ...] points distributed in the shape
 */
export function tessellation33434Distribution(center, shapeFunction, a, maxDist, maxPointNo) {
    const u = a / 2;
    const cx = center[0], cy = center[1];
    const allowedPointNo = maxPointNo || Infinity;
    const allowedSize = maxDist || Infinity;
    if (allowedPointNo === Infinity && allowedSize === Infinity) throw new Error("one of the end condition must be provided");
    let res = [];
    let n = 0;
    let fail = 0;
    while (res.length / 2 < allowedPointNo) {
        let offset = ulam(n);
        const xx = cx + offset[0] * u * 6, yy = cy + offset[1] * u * 6;
        let candidates = [[xx, yy + u], [xx + 2 * u, yy + u], [xx + u, yy - u], [xx + u, yy - u * 3],
        [xx - u, yy - 2 * u], [xx - u * 3, yy - u * 2], [xx - u * 2, yy], [xx - u * 2, yy + u * 2]];
        for (const p of candidates) {
            if (shapeFunction(p)) {
                res.push(...p);
                fail = 0;
            } else {
                fail++;
            }
        }
        n++;

        // check maxDist
        let d = a * Math.min(Math.abs(offset[0]), Math.abs(offset[1]));
        d -= 3 * u
        if (d > allowedSize) break;
        if (fail > 999) break;
    }
    // 
    return res;
}

// growing with Ulam spiral
const ulam = (n) => {
    const cycleNumber = (n) => {
        const x = Math.floor(Math.sqrt(n));
        return Math.ceil(x / 2);
    }

    const firstNumberInCycle = (k) => {
        return 4 * k * (k - 1) + 1;
    }

    const sideLengthInCycle = (k) => {
        return 2 * k;
    }

    const cornerVectors = [[1, -1], [-1, -1], [-1, 1], [1, 1]];
    const sideDirections = [[-1, 0], [0, 1], [1, 0], [0, -1]];
    if (n === 0) return [0, 0];
    const k = cycleNumber(n);
    const distFromStart = n - firstNumberInCycle(k);
    const sideLength = sideLengthInCycle(k);
    const side = Math.floor(distFromStart / sideLength);
    const distAlongSide = 1 + distFromStart % sideLength;
    let pos = cornerVectors[side].slice();
    pos[0] *= k, pos[1] *= k;
    pos[0] += distAlongSide * sideDirections[side][0], pos[1] += distAlongSide * sideDirections[side][1];
    return pos;
}

/**
 * 
 * @param {Array} center [x,y] center (starting point) of the point distribution
 * @param {Function} shapeFunction a function to check if a point is with in the shape, it should accept [x,y] as input and return a boolean value
 * @param {number} gap control the gap between points
 * @param {number} maxDist stop after the distribution reach this size
 * @param {number} maxPointNo stop after this number of points are placed
 * @returns [x1, y1, x2, y2, ...] points distributed in the shape
 */
export function phyllotaxisSpiralDistribution(center, shapeFunction, gap, maxDist, maxPointNo) {
    const baseAngle = Math.PI * (1 + Math.sqrt(5));
    const cx = center[0], cy = center[1];
    const allowedPointNo = maxPointNo || Infinity;
    const allowedSize = maxDist || Infinity;
    if (allowedPointNo === Infinity && allowedSize === Infinity) throw new Error("one of the end condition must be provided");
    let res = [];
    let n = 0, fail = 0;
    while (res.length < allowedPointNo) {
        let a = baseAngle * n;
        let r = gap * Math.sqrt(n);
        let x = cx + r * Math.cos(a), y = cy + r * Math.sin(a);
        if (shapeFunction([x, y])) {
            res.push(x, y);
            fail = 0;
        } else {
            fail++;
        }
        n++;
        if (r > maxDist) break;
        if (fail > 999) break;
    }
    return res;
}

/**
 * 
 * @param {Array} center [x,y] center (starting point) of the point distribution
 * @param {Function} shapeFunction a function to check if a point is with in the shape, it should accept [x,y] as input and return a boolean value
 * @param {number} numOfSeed number of point to start with
 * @param {Array} seedingArea [w, h] area to place the random seed
 * @param {number} relaxRound number of times of voronoi relaxation
 * @param {Function} weightMapFunction a function taking a point [x,y] and returning a weight (0-1)
 */
export function weightedVoronoiDistributionAndGraph(center, shapeFunction, numOfSeed, seedingArea, relaxRound, weightMapFunction) {
    const cx = center[0], cy = center[1];

    const weightedVoronoiRelaxation = function (voronoi) {
        let width = voronoi.xmax - voronoi.xmin;
        let seed = voronoi.delaunay.points;
        let weights = new Array(seed.length / 2).fill(0), centroids = new Array(seed.length).fill(0), counts = new Array(seed.length / 2).fill(0);
        let seedIndex = 0;
        for (let x = voronoi.xmin; x < voronoi.xmax; x++) {
            for (let y = voronoi.ymin; y < voronoi.ymax; y++) {
                let weight = weightMapFunction([x, y]);
                //With d3.Delaunay, we can use voronoi.delaunay.find() to find the closest Voronoi cell seed to an location
                seedIndex = voronoi.delaunay.find(x, y, seedIndex);
                centroids[seedIndex * 2] += x * weight;
                centroids[seedIndex * 2 + 1] += y * weight;
                weights[seedIndex] += weight;
                counts[seedIndex]++;
            }
        }
        for (let i = 0; i < centroids.length; i += 2) {
            if (weights[i / 2] > 0) {
                centroids[i] /= weights[i / 2];
                centroids[i + 1] /= weights[i / 2];
            } else {
                centroids[i] = seed[i];
                centroids[i + 1] = seed[i + 1];
            }
        }
        for (let i = 0; i < seed.length; i++) {
            seed[i] += (centroids[i] - seed[i]);
        }
        return new Delaunay(seed).voronoi([voronoi.xmin, voronoi.ymin, voronoi.xmax, voronoi.ymax]);
    }

    let seed = [];
    for (let i = 0; i < numOfSeed; i++) {
        let xx = cx - seedingArea[0] / 2 + Math.random() * seedingArea[0];
        let yy = cy - seedingArea[1] / 2 + Math.random() * seedingArea[1];
        seed.push(xx, yy);
    }
    let voronoi = new Delaunay(seed).voronoi([cx - seedingArea[0] / 2, cy - seedingArea[1] / 2, cx + seedingArea[0] / 2, cy + seedingArea[1] / 2]);
    for (let i = 0; i < relaxRound; i++) {
        voronoi = weightedVoronoiRelaxation(voronoi);
    }
    let ids = [], ps = [], idMap = {};
    let count = 0;
    let graph = new Set();
    for (let i = 0; i < voronoi.delaunay.points.length; i += 2) {
        const p = [voronoi.delaunay.points[i], voronoi.delaunay.points[i + 1]];
        const id = i / 2;
        if (shapeFunction(p)) {
            ids.push(id);
            ps.push(...p);
            idMap[id] = count;
            count++;
        }
    }
    for (let i = 0; i < ids.length; ++i) {
        let nei = voronoi.neighbors(ids[i]);
        for (const id of nei) {
            if (ids.includes(id)) {
                let rId = idMap[id];
                let str = rId > i ? i + " " + rId : rId + " " + i;
                graph.add(str);
            }
        }
    }

    return [ps, graph];
}

/**
 * 
 * @param {Array} points [x1, y1, x2, y2 ...]
 * @param {number} nearestN return links to the nearest N points, set it to Infinite or negative number will disable this limit
 * @param {number} maxDist max distance between points allowed
 * @returns a Set ["p1 p2", ...]. In each string p1 and p2 are the indices of the origin points array
 */
export function distanceGraph(points, nearestN, maxDist) {
    if (nearestN === Infinity || nearestN < 0) nearestN = points.length / 2;
    if (!Array.isArray(points)) throw new Error("points need to be an array");
    let tree = new KdTree([], (a, b) => (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2, [0, 1]);
    for (let i = 0; i < points.length; i += 2) {
        tree.insert([points[i], points[i + 1], i / 2]);
    }
    let lines = new Set();
    for (let i = 0; i < points.length / 2; ++i) {
        const currentP = [points[i * 2], points[i * 2 + 1]];
        let list = tree.nearest(currentP, nearestN + 1, maxDist ** 2);
        list = list.filter(item => item[1] > 0);
        for (let j = 0; j < list.length; ++j) {
            let k = list[j][0][2];
            const id = k > i ? i + " " + k : k + " " + i;
            lines.add(id);
        }
    }
    return lines;
}

/**
 * 
 * @param {Array} points [x1, y1, x2, y2...] 
 * @returns array of line segments [[x1, y1, x2, y2], ...] 
 */
export function primsSpanningTree(points) {
    let p0 = [points[0], points[1]];
    let root = { x: p0[0], y: p0[1], children: [] };
    let unreachedTree = new KdTree([], (a, b) => (a.x - b.x) ** 2 + (a.y - b.y) ** 2, ["x", "y"]);
    for (let i = 2; i < points.length; i += 2) {
        unreachedTree.insert({ x: points[i], y: points[i + 1], children: [] });
    }
    let reached = [root], remaining = points.length / 2 - 1;
    while (remaining > 0) {
        let re = Infinity, start, end;
        for (let i = 0; i < reached.length; ++i) {
            let tt = unreachedTree.nearest(reached[i], 1)[0];
            if (tt[1] < re) {
                re = tt[1], start = reached[i], end = tt[0];
            }
        }
        start.children.push(end)
        unreachedTree.remove(end);
        reached.push(end);
        remaining--;
    }

    return root;
}

/**
 * 
 * @param {Array} points [x1,y1,x2,y2,...]
 * @param {Array} graph 
 * @returns 
 */
export function kruskalsSpanningTree(points, graph){
    const v = points.length / 2;
    
    const union = (subsets, x, y) => {
      let rootX = findRoot(subsets, x); 
      let rootY = findRoot(subsets, y); 
      if (subsets[rootY].rank < subsets[rootX].rank) { 
          subsets[rootY].parent = rootX; 
      } 
      else if (subsets[rootX].rank < subsets[rootY].rank) { 
          subsets[rootX].parent = rootY; 
      } 
      else { 
          subsets[rootY].parent = rootX; 
          subsets[rootX].rank++; 
      } 
    }
  
    const findRoot = (subsets, i) => {
      if (subsets[i].parent === i) return subsets[i].parent; 
      subsets[i].parent = findRoot(subsets, subsets[i].parent); 
      return subsets[i].parent; 
    }
    let nodes = [];
    let edges = [];
    for (let i = 0; i < points.length; i += 2) {
      nodes.push({x: points[i], y: points[i + 1], id: i / 2, children: []});
    }
    for (let pair of graph) {
      const [a,b] = pair.split(" ").map(it => parseInt(it));
      const w = (nodes[a].x - nodes[b].x) ** 2 + (nodes[a].y - nodes[b].y) ** 2;
      edges.push({start:nodes[a], end:nodes[b], weight:w});
    }
    edges.sort((a,b) => a.weight - b.weight);
    console.log(edges)
  
    let j = 0, numberOfEdge = 0;
    let subsets = [];
    let res = [];
    for (let i = 0; i < v; ++i) {
      subsets.push({parent:i, rank:0});
    }
  
    while (numberOfEdge < v-1) {
      let nextEdge = edges[j];
      let x = findRoot(subsets, nextEdge.start.id);
      let y = findRoot(subsets, nextEdge.end.id);
      if (x != y) { 
          res[numberOfEdge] = nextEdge; 
          union(subsets, x, y); 
          numberOfEdge++; 
      } 
      j++;
    }
    return res
  }

/**
 * 
 * @param {Array} points [x1,y1,x2,y2,...]
 * @param {Array} graph corresponding graph from points
 * @param {number} mode 0,1,2,oe 3
 * @returns 
 */
export function turtleSpanningTree(points, graph, mode = 0) {
    let tem = [];
    for (let i = 0; i < points.length; i += 2) {
        tem.push({ x: points[i], y: points[i + 1], children: [], id: i / 2 });
    }
    let root = tem[0];

    const getNei = (i) => {
        let res = [];
        for (const pair of graph) {
            let arr = pair.split(" ").map(it => parseInt(it));
            if (arr[0] === i || arr[1] === i) res.push(arr[0] === i ? arr[1] : arr[0]);
        }
        res = res.filter(nei => !reached.includes(nei));
        return res.map(o => tem[o]);
    }

    let current = root;
    let reached = [0];
    let queue = [root];
    switch (mode) {
        default: break;
        case 0:
            //depth first
            while (queue.length > 0) {
                let options = getNei(current.id);
                if (options.length > 0) {
                    let choice = options[Math.floor(Math.random() * options.length)];
                    current.children.push(choice);
                    reached.push(choice.id);
                    queue.unshift(choice);
                    current = choice;
                } else {
                    current = queue.pop();
                }
            }
            break;
        case 1:
            //breadth first
            while (queue.length > 0) {
                let options = getNei(current.id);
                //sortPointsCCW(current, options);
                options.forEach(choice => {
                    reached.push(choice.id);
                    current.children.push(choice);
                    queue.push(choice);
                })
                current = queue.shift();
            }
            break;
        case 2:
            //mixing depth and breadth first
            while (queue.length > 0) {
                let options = getNei(current.id);
                //sortPointsCCW(current, options);
                options.forEach(choice => {
                    reached.push(choice.id);
                    current.children.push(choice);
                    if (Math.random() > 0.5) queue.push(choice);
                    else queue.unshift(choice);
                })
                current = queue.shift();
            }
            break;
        case 3:
            //3-3
            let bool = false, count = 0;
            //depth first
            while (queue.length > 0) {
                count++;
                if (count > 3) {
                    bool = !bool;
                    count = 0;
                }
                let options = getNei(current.id);
                sortPointsCCW(current, options);
                if (options.length > 0) {
                    let choice = options[bool ? 0 : options.length - 1];
                    current.children.push(choice);
                    reached.push(choice.id);
                    queue.unshift(choice);
                    current = choice;
                } else {
                    current = queue.pop();
                }
            }
            break;
    }

    return root;
}

function sortPointsCCW (center, points) {
    let startAng;
    points.forEach(point => {
        let ang = Math.atan2(point.y - center.y, point.x - center.x);
        if (!startAng) { startAng = ang }
        else {
            if (ang < startAng) {
                ang += Math.PI * 2;
            }
        }
        point.angle = ang;
    });
    points.sort((a, b) => b.angle - a.angle);
    points.forEach(p => delete p.angle)
}

