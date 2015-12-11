##WebGL Multiple Importance Sampling Path Tracer
###xinyue Zhu
<p>Personal Website: https://seas.upenn.edu/~zhuxiny</p>
![](img/thumb.png)
<p>address: http://zhuxinyue909.github.io/CIS-565-Final-WebGL-MIS/</p>
<p>Desbug Mode:https://www.shadertoy.com/view/Xs33WM</p>
Technique:
1.multiple importance sampling
The Monte Carlo path tracer attempts to solve the light equations for all the visible points in the scene which ask people to recursively trace the ray until it hits the depth. However it may take a long time to converge to a usable image. To optimize the monte carlo path tracer poeple often use parallelization and multiple importance sampling methods. 
<p>BRDF:A function that evaluates the energy emitted along ray given the intersection point of the scene and the direction from which the incoming light emits which is entirely dependent on the attributes of the material sampled at the intersection point. </p>
<p>So In order to reduce the number of the samples that needed to produce the converged scene, we use multiple importace sampling method.</p>
<img src="img/shadertoy.png">
![](img/shadertoy.png)
<img src="img/s1.png"  width="330" height="200">
<img scr="img/shadertoy.png" width="330" height="200">
<p>These following are debug view posted on: https://www.shadertoy.com/view/Xs33WM</p>


<p>Direct Light Sampling: </p>
<p>the following is the debug scene, the material from left to right is blinn-microface(exponent=20),blinn-microface(exponent=10),perfect reflection,blinn-microface(exponent=50),blinn-microface(exponent=100)</p>
<p>when the radius of light is quals 0.3</p>
![](img/light_brdf_r0.3.png)
<p>when the radius of light is quals 0.5</p>
![](img/light_brdf_r0.5.png)
<p>when the radius of light is quals 1.0</p>
![](img/light_brdf_r1.0.png)
<p>Since the light sources are the most important elements in a rendered scene, for some subset of the rays ωi , select directions such that each ωi intersects a given light source at some point</p>
<p>BRDF Sampling: </p>
[![](img/radiance.png)](https://www.youtube.com/watch?v=TKP8JBcbNN8&feature=youtu.be)
<p>Sample ray directions that have a higher contribution to the color reflected along ωo, which is extremely useful in the case where the BRDF has a very narrow set of contributing rays like perfectly reflection case.</p>
<p>the radiance of light when sampling brdf</p>
https://youtu.be/TKP8JBcbNN8
<p>the following image shows the pdf when sampling the brdf of each material.
![](img/brdf_pdf.png)
<p>When a BRDF is more specular, sampling only the light’s PDF makes it less likely that a large light will contribute to its color </p>
<p>When a BRDF is more diffuse, sampling only the BRDF’s PDF makes it less likely that small lights will contribute to its color </p>


<p>2.Gamma correction</p>
In CRT displays, the light intensity varies nonlinearly with the electron-gun voltage. Altering the input signal by gamma compression can cancel this nonlinearity, such that the output picture has the intended luminance. In this circumstance, use gammma correction to modify the color to make it look much prettier., When the GAMMA powers larger than 1 it makes the shadows darker, while Gamma is  smaller than 1 make dark regions lighter.
[![](img/gamma0.png)]
[![](img/gamma1.png)]





















