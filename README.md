### WebGL path Tracer
##1. Implementation
1.1 WebGL
use webGL directly. 
Math library: 1)glMatrix: https://github.com/toji/gl-matrix or 2)Sylvester:https://github.com/jcoglan/sylvester
Base code: WebGL deffered shader 
1.2 method
LD(p,ωo) = 1/N ΣNj (w(ωj)f(p,ωo,ωj)Ld(p,ωj)absdot(ωj, N)/p(ωj))
Using the multiple importance sampling method instead of naïve Monte Carlo method which can accelerate the converging process make the picture becomes more accurate. The formula writes above. Adding the direct lighting and indirect lighting results together. To reduce the number of samples:
1)	Direct Light sampling: only samples the part that intersect with light 
2)	BRDF sampling: sample ray directions that have higher contribution to the reflect color(it is not calculate randomly it’s generate according to the brdf)
Finally multiply of the sampling result with the weight parameter(heuristic power) and add the result together.
1.3 rendering effects/fucntion
Render the flowing effects: perfect specular, lambert, blinn phong/blinn phong microfacet method(if have time I will try different material by using different term like: Cook Torrance, Walter)
1.4 performance analysis
Performance contrast of with MIS(multiple importance sampling) and naïve method by comparing the FPS, rendering result, time. 
