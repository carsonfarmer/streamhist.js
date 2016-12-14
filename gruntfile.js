module.exports = function(grunt) {
    grunt.initConfig({
        pkg: grunt.file.readJSON('package.json'),
        pure_cjs: {
            options: {
                exports: 'streamhist',
                comments: true
            },
            'dist/streamhist.js': './index.js'
        }
    });

    grunt.task.loadNpmTasks('grunt-pure-cjs');
    // grunt.registerTask('dist', ['pure_cjs']);
    grunt.registerTask('default', ['pure_cjs']);
};
