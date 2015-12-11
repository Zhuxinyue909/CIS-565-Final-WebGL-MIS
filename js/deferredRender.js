(function() {
    'use strict';
    // deferredSetup.js must be loaded first

    R.deferredRender = function(state) {
        if (!aborted && (
            !R.progCopy ||
            !R.progRed ||
            !R.progClear ||
            !R.prog_BlinnPhong_PointLight ||
            !R.prog_Debug ||
            !R.prog_pathtrace||
            !R.prog_path_trace_debug)) {
            console.log('waiting for programs to load...');
            return;
        }

        // Move the R.lights
        for (var i = 0; i < R.lights.length; i++) {
            // OPTIONAL TODO: Edit if you want to change how lights move
            var mn = R.light_min[1];
            var mx = R.light_max[1];
           // R.lights[i].pos[1] = (R.lights[i].pos[1] + R.light_dt - mn + mx) % mx + mn;
        }

        // Execute deferred shading pipeline



      
        R.pass_copy.render(state);

        if (cfg && cfg.debugView >= 0) {
          
            R.pass_debug.render(state);
        }
         else if(!cfg.scene2){ 
            R.pass_pathtracing.render(state);
         
        }
        else 
        {
            R.pass_path_trace_debug.render(state);
        }
    };

    /**
     * 'copy' pass: Render into g-buffers
     */
    R.pass_copy.render = function(state) {
        // * Bind the framebuffer R.pass_copy.fbo
        gl.bindFramebuffer(gl.FRAMEBUFFER,R.pass_copy.fbo);
        renderFullScreenQuad(R.progClear);
        gl.clearDepth(1.0);
        gl.clear(gl.DEPTH_BUFFER_BIT);

        // * "Use" the program R.progCopy.prog
        gl.useProgram(R.progCopy.prog); 
        var m = state.cameraMat.elements;
        // * Upload the camera matrix m to the uniform R.progCopy.u_cameraMat
        //   using gl.uniformMatrix4fv
        gl.uniformMatrix4fv(R.progCopy.u_cameraMat,false,m);
        gl.uniform3f(R.progCopy.u_spec, cfg.specular,0.0,0.0);
        drawScene(state);
        drawScene(state);
    };

    var drawScene = function(state) {
        for (var i = 0; i < state.models.length; i++) {
            var m = state.models[i];

            // If you want to render one model many times, note:
            // readyModelForDraw only needs to be called once.
            readyModelForDraw(R.progCopy, m);

            drawReadyModel(m);
        }
    };

    R.pass_debug.render = function(state) {
        // * Unbind any framebuffer, so we can write to the screen
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);

        // * Bind/setup the debug "lighting" pass
        // * Tell shader which debug view to use
        bindTexturesForLightPass(R.prog_Debug);
        gl.uniform1i(R.prog_Debug.u_debug, cfg.debugView);

        // * Render a fullscreen quad to perform shading on
        renderFullScreenQuad(R.prog_Debug);
    };

    /**
     * 'deferred' pass: Add lighting results for each individual light
     */
     R.pass_path_trace_debug.render=function(state) {
        // * Bind R.pass_deferred.fbo to write into for later postprocessing
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);

        // * Clear depth to 1.0 and color to black
        gl.clearColor(0.0, 0.0, 0.0, 0.0);
        gl.clearDepth(1.0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

        gl.enable(gl.BLEND);
        gl.blendEquation( gl.FUNC_ADD );
        gl.blendFunc(gl.ONE,gl.ONE);


        bindTexturesForLightPass(R.prog_BlinnPhong_PointLight);
  for (var i = 0; i < R.lights.length; i++) 
  {
        // TODO: add a loop here, over the values in R.lights, which sets the
       bindTexturesForLightPass(R.prog_path_trace_debug);
       gl.uniform1f(R.prog_path_trace_debug.iGlobalTime, state.iGlobalTime);
    
       renderFullScreenQuad(R.prog_path_trace_debug);
        }
        
        // Disable blending so that it doesn't affect other code
        gl.disable(gl.BLEND);
    };
    R.pass_pathtracing.render = function(state) {
        // * Bind R.pass_deferred.fbo to write into for later postprocessing
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);

        // * Clear depth to 1.0 and color to black
        gl.clearColor(0.0, 0.0, 0.0, 0.0);
        gl.clearDepth(1.0);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

        gl.enable(gl.BLEND);
        gl.blendEquation( gl.FUNC_ADD );
        gl.blendFunc(gl.ONE,gl.ONE);


        bindTexturesForLightPass(R.prog_BlinnPhong_PointLight);
  for (var i = 0; i < R.lights.length; i++) {
        // TODO: add a loop here, over the values in R.lights, which sets the
       bindTexturesForLightPass(R.prog_pathtrace);
       gl.uniform1f(R.prog_pathtrace.iGlobalTime, state.iGlobalTime);
       gl.uniform3f(R.prog_pathtrace.u_sphere_pos,cfg.pos_x,cfg.pos_y,cfg.pos_z);
       gl.uniform1f(R.prog_pathtrace.u_intensity,cfg.intensity);
       if(cfg.GammaCorrection){
                   gl.uniform1f(R.prog_pathtrace.if_gamma,1.0);
       }
       else
        gl.uniform1f(R.prog_pathtrace.if_gamma,-1.0);
       
        renderFullScreenQuad(R.prog_pathtrace);
        }
        
        // Disable blending so that it doesn't affect other code
        gl.disable(gl.BLEND);
    };

    var bindTexturesForLightPass = function(prog) {
        gl.useProgram(prog.prog);

        // * Bind all of the g-buffers and depth buffer as texture uniform
        //   inputs to the shader
        for (var i = 0; i < R.NUM_GBUFFERS; i++) {
            gl.activeTexture(gl['TEXTURE' + i]);
            gl.bindTexture(gl.TEXTURE_2D, R.pass_copy.gbufs[i]);
            gl.uniform1i(prog.u_gbufs[i], i);
        }
        gl.activeTexture(gl['TEXTURE' + R.NUM_GBUFFERS]);
        gl.bindTexture(gl.TEXTURE_2D, R.pass_copy.depthTex);
        gl.uniform1i(prog.u_depth, R.NUM_GBUFFERS);
    };

    /**
     * 'post1' pass: Perform (first) pass of post-processing
     */
    

    var renderFullScreenQuad = (function() {

        var positions = new Float32Array([
            -1.0, -1.0, 0.0,
             1.0, -1.0, 0.0,
            -1.0,  1.0, 0.0,
             1.0,  1.0, 0.0
        ]);

        var vbo = null;

        var init = function() {
            // Create a new buffer with gl.createBuffer, and save it as vbo.
            vbo=gl.createBuffer();
            gl.bindBuffer(gl.ARRAY_BUFFER,vbo);
            gl.bufferData(gl.ARRAY_BUFFER,positions,gl.STATIC_DRAW);
        };

        return function(prog) {
            if (!vbo) {
                // If the vbo hasn't been initialized, initialize it.
                init();
            }

            // Bind the program to use to draw the quad
            gl.useProgram(prog.prog);

            gl.bindBuffer(gl.ARRAY_BUFFER,vbo);
            gl.enableVertexAttribArray(prog.a_position);
            gl.vertexAttribPointer(prog.a_position, 3, gl.FLOAT, false, 0, 0);
            gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

            // Unbind the array buffer.
            gl.bindBuffer(gl.ARRAY_BUFFER, null);
        };
    })();
})();
