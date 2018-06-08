import nodeResolve from 'rollup-plugin-node-resolve';
import commonjs from 'rollup-plugin-commonjs';

export default {
    input: __dirname + '/src/main.js',
    plugins: [
        nodeResolve(),
        commonjs()
    ],
    // sourceMap: true,
    output: [
        {
            format: 'umd',
            name: 'geometryExtrude',
            file: 'dist/geometry-extrude.js'
        }
    ]
};