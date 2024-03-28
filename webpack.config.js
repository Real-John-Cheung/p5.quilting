const path = require("path");

module.exports = {
    entry: "./src/p5.quilting.js",
    output:{
        filename: "p5.quilting.min.js",
        globalObject: 'window',
        library: {
            name:"Quilting",
            type: "umd"
        },
        path: path.resolve(__dirname, "dist"),
    }
}