##WebGL Multiple Importance Sampling Path Tracer
###xinyue Zhu
<p>Personal Website: https://seas.upenn.edu/~zhuxiny
[![](img/xxx.png)]
<p>address: http://zhuxinyue909.github.io/CIS-565-Final-WebGL-MIS/</p>
Technique:
1.multiple importance sampling
The Monte Carlo path tracer attempts to solve the light equations for all the visible points in the scene which ask people to recursively trace the ray until it hits the depth. However it may take a long time to converge to a usable image. To optimize the monte carlo path tracer poeple often use parallelization and multiple importance sampling methods. 
<p>BRDF:A function that evaluates the energy emitted along ray given the intersection point of the scene and the direction from which the incoming light emits which is entirely dependent on the attributes of the material sampled at the intersection point. </p>
<p>So In order to reduce the number of the samples that needed to produce the converged scene, we use multiple importace sampling method.</p> 
<p>Direct Light Sampling: </p>
<p>Since the light sources are the most important elements in a rendered scene, for some subset of the rays ωi , select directions such that each ωi intersects a given light source at some point</p>
<p>BRDF Sampling: </p>
<p>Sample ray directions that have a higher contribution to the color reflected along ωo, which is extremely useful in the case where the BRDF has a very narrow set of contributing rays like perfectly reflection case.</p>

<p>When a BRDF is more specular, sampling only the light’s PDF makes it less likely that a large light will contribute to its color </p>
[![](img/pic1.png)]
<p>When a BRDF is more diffuse, sampling only the BRDF’s PDF makes it less likely that small lights will contribute to its color </p>

[![](img/pic2.png)]
[![](img/pic3.png)]
2.materials






















