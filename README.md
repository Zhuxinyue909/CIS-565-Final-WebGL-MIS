## WebGL path Tracer
###1. Implementation
<p>1.1 WebGL</p>
<p>use webGL directly. </p>
<p>Math library: 1)glMatrix: https://github.com/toji/gl-matrix or</p>
<p>2)Sylvester:https://github.com/jcoglan/sylvester</p>
<p>Base code: WebGL deffered shader:</p>
<p>https://github.com/CIS565-Fall-2015/Project6-WebGL-Deferred-Shading</p>

<p>1.2 method</p>
<p>LD(p,ωo) = 1/N ΣNj (w(ωj)f(p,ωo,ωj)Ld(p,ωj)absdot(ωj, N)/p(ωj))</p>
<p>Using the multiple importance sampling method instead of naïve Monte Carlo method which can accelerate the converging process make the picture becomes more accurate. The formula writes above. Adding the direct lighting and indirect lighting results together. </p>
<p>To reduce the number of samples:</p>
<p>1)	Direct Light sampling: only samples the part that intersect with light </p>
<p>2)	BRDF sampling: sample ray directions that have higher contribution to the reflect color(it is not calculate randomly it’s generate according to the brdf)</p>
<p>Finally multiply of the sampling result with the weight parameter(heuristic power) and add the result together.</p>
<p>1.3 rendering effects/fucntion</p>
<p>Render the flowing effects:</p> 
<p>perfect specular, lambert, blinn phong/blinn phong microfacet method(if have time I will try different material by using different term like: Cook Torrance, Walter)</p>
<p>1.4 performance analysis</p>
<p>Performance contrast of with MIS(multiple importance sampling) and naïve method by comparing the FPS, rendering result, time. </p>
[![](img/thumb.png)]

