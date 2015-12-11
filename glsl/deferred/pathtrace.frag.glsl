#version 100
precision highp float;
precision highp int;

#define SAMPLES 		    24
#define EPSILON 			0.00001 
#define GAMMA 				3.2			
#define DEPTH		        5
#define INDIRECT_CLAMP 		20.0	
#define INFINITY            9999.9
#define backgroundColor     vec3( 0.0 )
#define PI 					3.1415926
#define TWO_PI 				6.2831852

#define OBJ_NULL 0
#define OBJ_CUBE 1
#define OBJ_SPHERE 2
#define OBJ_PLANE 3
uniform float if_gamma;
uniform float u_intensity;
uniform float iGlobalTime;
uniform vec3  u_sphere_pos;
varying vec2 v_uv;



float seed;	//seed initialized in main
float rnd() { return fract(sin(seed++)*43758.5453123); }

struct Ray {
    vec3 origin; 
    vec3 dir; 
};
struct Intersection
{
    vec3 point;
    vec3 normal;
    float t;
    int hit_obj;
    //vec3 color;
};
struct Sphere { int materialId; vec3 pos; float r;  float area; };
struct Plane { int materialId; vec4 abcd; };
struct Material { vec3 diff_; vec3 spec_; float roughness_; };
struct RayHit { vec3 pos; vec3 N; vec3 E; vec2 uv; int materialId; };
struct Camera { mat3 rotate; vec3 pos; vec3 target; float fovV; float lensSize; float focusDist; };
struct LightPathNode { vec3 pos_; vec3 N; vec3 L; vec3 Li; int materialId; };
struct Box { vec3 Bl; vec3 Bh; int materialId; };
    
// ************ SCENE ***************
Box boxes[1];
Plane walls[5];
Sphere spherelight;
Sphere sphereGeometry;
Camera camera;
//***********************************



#define WHITE_DIFFUSE			0
#define MTL_WALL			1
#define MTL_WALL_BOTTOM		2
#define MTL_WALL_LEFT		3
#define MTL_WALL_RIGHT		4
#define MTL_GLOSSY	        5
#define MTL_MIRRO	        6
#define MTL_COUNT			7

Material materialLibrary[MTL_COUNT];
#define INIT_MTL(i,diff,spec,roughness) { materialLibrary[i].diff_=diff; materialLibrary[i].spec_=spec; materialLibrary[i].roughness_=roughness; }
void initMaterial() {
   // float intensity = 50.0 ;
   float intensity = u_intensity ;
    INIT_MTL( WHITE_DIFFUSE, vec3( intensity ), vec3( 1.0 ), 1.0 );
    INIT_MTL( MTL_WALL, vec3( 1.0 ), vec3( 0.2 ), 0.7 );
    INIT_MTL( MTL_WALL_BOTTOM, vec3( 0.5, 0.6, 0.8 ), vec3( 0.8 ), 0.3 );
    INIT_MTL( MTL_WALL_LEFT, vec3( 1.0, 0.1, 0.15 ), vec3( 0.2 ), 0.7 );
    INIT_MTL( MTL_WALL_RIGHT, vec3( 0.1, 0.8, 0.1 ), vec3( 0.2 ), 0.7 );
    INIT_MTL( MTL_GLOSSY, vec3( 0.7,0.1,0.4 ), vec3( 1.0 ), 0.3 );
}

Material getMaterialFromLibrary( int index ) {

    for(int i=0; i<MTL_COUNT; i++ ) { 
    	if( index == i ) {
            return materialLibrary[i];
        }
    }
    
    return materialLibrary[0];
}
float iSphere( in vec3 ro, in vec3 rd, in vec4 sph ) {
    vec3 oc = ro - sph.xyz;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - sph.w * sph.w;
    float h = b * b - c;
    if (h < 0.0) return -1.0;

	float s = sqrt(h);
	float t1 = -b - s;
	float t2 = -b + s;
	
	return t1 < 0.0 ? t2 : t1;
}
float iPlane( in vec3 ro, in vec3 rd, in vec4 pla ) {
    return (-pla.w - dot(pla.xyz,ro)) / dot( pla.xyz, rd );
}
void initLightSphere( float time ) {
	spherelight.pos = vec3( 2.0+2.*sin(time*0.9),5.5,-4.0);
}
void initScene() {
    float time = iGlobalTime;
    
    //init lights
    float r = 0.5;
	float r1=0.2;
    spherelight = Sphere( WHITE_DIFFUSE, vec3( 2.0, 5.5, -4.0 ), r, r*r*4.0*PI );

    //init walls
    //bottom
    walls[0].abcd = vec4( 0.0, 1.0, 0.0, 1.0 );
    walls[0].materialId = MTL_WALL_BOTTOM;
    
    //top
    walls[1].abcd = vec4( 0.0, -1.0, 0.0, 5.5 );
    walls[1].materialId = MTL_WALL;
    
    //front
    walls[2].abcd = vec4( 0.0, 0.0, 1.0, 8.0 );
    walls[2].materialId = MTL_WALL;
    
    //left
    walls[3].abcd = vec4( 1.0, 0.0, 0.0, 4.0 );
    walls[3].materialId = MTL_WALL_LEFT;
    
    //right
    walls[4].abcd = vec4( -1.0, 0.0, 0.0, 4.0 );
    walls[4].materialId = MTL_WALL_RIGHT;
    
    //box
    boxes[0].Bl = vec3( 0.2, -1.0, -2.5 );
    boxes[0].Bh = vec3( 2.0, 0.8, -1.6 );
    boxes[0].materialId = MTL_WALL_LEFT;
    
    r = 1.0;
   // sphereGeometry = Sphere( MTL_GLOSSY, vec3( -0.5, 0.0, -3.0 ), r, r*r*4.0*PI  );
   sphereGeometry = Sphere( MTL_MIRRO, u_sphere_pos, r, r*r*4.0*PI  );
}


// intersection function
Intersection raySphereIntersection( Ray r, in Sphere sphr) {
    
    Intersection isx;
    isx.t=-1.0;
   // isx.hit_obj=OBJ_NULL;
    vec3 os=r.origin-sphr.pos;
    //remember to normalize the r.direction
    float A=1.0;
    float B=2.0*(r.dir.x*os.x+r.dir.y*os.y+r.dir.z*os.z);
    float C = pow(os.x,2.0) + pow(os.y,2.0) +pow(os.z,2.0)-sphr.r*sphr.r;
    float Discre=pow(B,2.0)-4.0*A*C;
    
    if(Discre < 0.0) {
        isx.t=-1.0;
        return isx;
    }
    else{
        float t0= (-B - sqrt(Discre))/(2.0*A);
        if(t0>0.0)//intersecton
        {
            isx.point = r.dir*t0+r.origin;//to Local
            isx.normal= normalize(isx.point-sphr.pos);//normalize!!!!!
            isx.t=t0;
           // isx.hit_obj=OBJ_SPHERE;
            //isx.color=sphr.base_color;
            return isx;
        }
        else{
            float t1= (-B + sqrt(Discre))/(2.0*A);
            if(t1>0.0)//t1 intersect at front in sphere
            {
                 isx.point = r.dir*t1+r.origin;
                 isx.normal =normalize(isx.point-sphr.pos);//normalize!!!!!
                 isx.t = t1;
                 //isx.hit_obj=OBJ_SPHERE;
                 // isx.color=sphr.base_color;
                 return isx;
            }
            else {
                isx.t=-1.0;
                return isx;
            }
        }
    }
    
}
Intersection rayCubeIntersection(in Ray r,in Box cube)
{
 float Tnear=-INFINITY;
 float Tfar=INFINITY;
    
 vec3 R=r.dir;
 vec3 Position=r.origin;
 vec3 point;
 
 Intersection isx;
 isx.t=-1.0;
 isx.hit_obj=OBJ_NULL;
 if(abs(R.x)<=EPSILON&&(Position.x < cube.Bl.x||Position.x > cube.Bh.x))return isx;
 else if(abs(R.y)<=EPSILON&&(Position.y<cube.Bl.y||Position.y>cube.Bh.y))return isx;
 else if(abs(R.z)<=EPSILON&&(Position.z<cube.Bl.z||Position.z>cube.Bh.z))return isx;
 else{
          float t1, t2;
          float a;
          int face=-1;
          t1 = (cube.Bl.z -Position.z)/R.z;
          t2 = (cube.Bh.z -Position.z)/R.z;
          if(t1>t2){
              a = t1;t1 = t2;t2 = a;
          }
          if(t1>Tnear){
              Tnear = t1; face=1;
          }
          if(t2<Tfar) Tfar = t2;
    //*********************x
          t1 = (cube.Bl.x - Position.x)/R.x;
          t2 = (cube.Bh.x - Position.x)/R.x;
          if(t1>t2){a = t1;t1 = t2;t2 = a;}
          if(t1>Tnear){Tnear = t1;face=2;}
          if(t2<Tfar) Tfar = t2;
          point.x=Tnear;
   //**********************Y
          t1 = (cube.Bl.y - Position.y)/R.y;
          t2 = (cube.Bh.y - Position.y)/R.y;
          if(t1>t2){a = t1;t1 = t2;t2 = a;}
          if(t1>Tnear) {Tnear = t1;face=3;}
          if(t2<Tfar) Tfar = t2;
          if(Tnear<Tfar){//swap -1
              vec3 normal;
              point = Position+ Tnear*R;//local
              if(face==1){normal=vec3(0.0,0.0,1.0);}
              if(face==2){normal=vec3(1.0,0.0,0.0);}
              if(face==3){normal=vec3(0.0,1.0,0.0);}
              vec3 _sign = sign(point);
              return Intersection(point,vec3(normal.x*_sign.x, normal.y*_sign.y, normal.z*_sign.z),Tnear,OBJ_CUBE);
          }         
          else
              return isx;
       }
}

bool rayPlaneIntersection( Ray ray, Plane plane, out float t ){
    float dotVN = dot( ray.dir, plane.abcd.xyz );
   
    if ( abs( dotVN ) > EPSILON ) {
		t = -(dot( ray.origin, plane.abcd.xyz ) + plane.abcd.w)/dotVN;
    }
    
    return ( t > 0.0 );
}


 
//-------Geometry functions--------
vec2 uniformPointWithinCircle( in float radius, in float Xi1, in float Xi2 ) {
    float r = radius*sqrt(Xi1);
    float theta = Xi2;
	return vec2( r*cos(theta), r*sin(theta) );
}

vec3 uniformDirectionWithinCone( in vec3 d, in float phi, in float sina, in float cosa ) {    
	vec3 w = normalize(d);
    vec3 u = normalize(cross(w.yzx, w));
    vec3 v = cross(w, u);
	return (u*cos(phi) + v*sin(phi)) * sina + w * cosa;
}

vec3 localToWorld( in vec3 localDir, in vec3 normal )
{
    vec3 binormal = normalize( ( abs(normal.x) > abs(normal.z) )?vec3( -normal.y, normal.x, 0.0 ):vec3( 0.0, -normal.z, normal.y ) );
	vec3 tangent = cross( binormal, normal );
    
	return localDir.x*tangent + localDir.y*binormal + localDir.z*normal;
}

vec3 sphericalToCartesian( in float rho, in float phi, in float theta ) {
    float sinTheta = sin(theta);
    return vec3( sinTheta*cos(phi), sinTheta*sin(phi), cos(theta) )*rho;
}

vec3 sampleHemisphereCosWeighted( in vec3 n, in float Xi1, in float Xi2 ) {
    float theta = acos(sqrt(1.0-Xi1));
    float phi = TWO_PI * Xi2;
    return localToWorld( sphericalToCartesian( 1.0, phi, theta ), n );
}

vec3 randomHemisphereDirection( const vec3 n, in float Xi1, in float Xi2 ) {
    vec2 r = vec2(Xi1,Xi2)*TWO_PI;
	vec3 dr=vec3(sin(r.x)*vec2(sin(r.y),cos(r.y)),cos(r.x));
	return dot(dr,n) * dr;
}

vec3 randomDirection( in float Xi1, in float Xi2 ) {
    float theta = acos(1.0 - 2.0*Xi1);
    float phi = TWO_PI * Xi2;
    
    return sphericalToCartesian( 1.0, phi, theta );
}


Ray shootRay( in vec2 pixel, in float rand1, in float rand2 )
{
    Ray ray;
    vec3 eye = vec3(0, 2.0, 12.00);
    vec3 ref = vec3(0, 2.0,  0.00);
    vec3 F = normalize( ref - eye );
    vec3 U = normalize( cross(F,vec3(0.0,1.0,0.0) ) );
    vec3 R = normalize( cross(U,F)); 
	vec3 rd = normalize( (pixel.x+rand1)*U + (pixel.y+rand2)*R +4.* F );
    ray.origin=eye;
    ray.dir=rd;  
	return ray;
}


bool raySceneIntersection( in Ray ray, out RayHit hit, out int objId, out float dist ) {
    float nearest_dist = INFINITY;
    
    //check lights
    Intersection isx0=raySphereIntersection( ray, spherelight);
    if( (isx0.t>0.0) && ( isx0.t < nearest_dist ) ) {
        nearest_dist = isx0.t;
        hit.pos = isx0.point;
        hit.N = isx0.normal;
        hit.materialId = WHITE_DIFFUSE;
        objId = 0;
    }
    
    //sphere
    Intersection isx1=raySphereIntersection( ray, sphereGeometry);
    if( (isx1.t>0.0) && ( isx1.t < nearest_dist ) ) 
    {
        nearest_dist = isx1.t;
        hit.pos =isx1.point;
        hit.N = isx1.normal;
        hit.materialId = MTL_MIRRO;
        objId = 1;
    }
    
    //walls
        for( int i=0; i<5; i++ ){
        float distToPlane;
        if( rayPlaneIntersection( ray, walls[i], distToPlane ) && (distToPlane > 0.0) && (distToPlane < nearest_dist ) ){
            nearest_dist = distToPlane;
            hit.pos = ray.origin + ray.dir*nearest_dist;
    		hit.N = walls[i].abcd.xyz;
    		hit.materialId = walls[i].materialId;
            objId = 2 + i;
        }
    }
    //cube
    float box_t0, box_t1;
    Intersection isx_cube = rayCubeIntersection(ray,boxes[0]);
    if( (isx_cube.t>0.0)  && (isx_cube.t < nearest_dist) ) 
    {
        nearest_dist = isx_cube.t;
        hit.pos = isx_cube.point;
        hit.N = isx_cube.normal;
        hit.materialId = boxes[0].materialId;
        
        objId = 6;
    }
    
    if( nearest_dist < INFINITY) {
        if( hit.pos.z > 0.0 ) {
            return false;
        }
        
    	hit.E = ray.dir*(-1.0);
    	dist = nearest_dist;
        return true;
            
    }
    
    return false;
}

vec3 blin_microfacet(vec3 N,vec3 wo,vec3 wi, vec3 cdiff) {

   //------ Geometry term-------:
    vec3 reflect_color=vec3(1.0,1.0,1.0);
    N=normalize(N);
    vec3 wh=normalize(wo+wi);
    float a=abs(dot(wh,N));
    float M=abs(2.0*dot(N,wh)*dot(N,wo)/dot(wo,wh));
    float S=abs(2.0*dot(N,wh)*dot(N,wi)/dot(wo,wh));
    float _G=min(min(M,S),1.0);


    //------Distribution term------
    float exponent=10.0;
    float _D=( exponent+2.0)/(2.0*PI)*pow(a, exponent);

    float costheta_o=dot(N,wo);//normalize!!!
    float costheta_i=dot(N,wi);
    float brdf=_D*_G/(4.0*costheta_o*costheta_i);
    
    //----color--calculation--
    
    return brdf*cdiff;
}

vec3 brdf_evaluate ( in vec3 N,	in vec3 wo,in vec3 wi,	int matId, vec3 cDiffuse,vec3 cSpecular ) {
    return  blin_microfacet( N, wo, wi, cDiffuse );
}

vec3 brdf_sample (in vec3 N,in vec3 E,out vec3 L,int matId,vec3 cDiffuse,vec3 cSpecular,float Xi1,float Xi2 ) 
{
    if(matId<5)L = sampleHemisphereCosWeighted( N, Xi1, Xi2 );
    else L=reflect(normalize (E), N);
    vec3 brdf = brdf_evaluate( N, E, L, matId, cDiffuse, cSpecular ); 
    float pdf = ( clamp( dot( N, L ), 0.0, 1.0 ) ) / PI;
    
    return brdf*pdf;
}

void sampleDirectLight( vec3 hitpoint, out vec3 dir, out float pdf,Ray r ) {
    
   /*float r=glm::distance2(isx.point,ray.origin);

   float cos_theta=glm::dot(-ray.direction,isx.normal);
   return r/(cos_theta*area);*/
   	vec3 rlen = spherelight.pos - hitpoint;
    float r2 = dot(rlen, rlen);
    float cos_a_max = sqrt( 1.0 - clamp( spherelight.r*spherelight.r / r2, 0.0, 1.0 ) );
    float omega = TWO_PI * (1.0 - cos_a_max);	
    float cosa = mix(cos_a_max, 1.0, rnd());
    float sina = sqrt(1.0 - cosa*cosa);

    dir = uniformDirectionWithinCone( rlen, TWO_PI*rnd(), sina, cosa );
    pdf = 1.0/omega;
    ////////////////
    float distancer=distance(hitpoint,r.origin);
    float dis=pow(distancer,2.0);
    //float cos_theta=dot(-r.dir,)
}

bool ShadowTest( Ray shadowRay ) {
    float distToHit;
    RayHit tmpHit;
    int tmpObjId;
    raySceneIntersection( shadowRay, tmpHit, tmpObjId, distToHit );
    
    return ( tmpObjId == 0 );
}

vec3 calcDirectLightOnSurface( RayHit hit, Material surfMtl,Ray r ) {
    vec3 Li = materialLibrary[WHITE_DIFFUSE].diff_;
    vec3 Lo = vec3( 0.0 );
    
    vec3 wi;
    float pWi;
    sampleDirectLight( hit.pos, wi, pWi,r );
    
    if ( dot( wi, hit.N ) > 0.0 ) {
        if ( ShadowTest( Ray( hit.pos + hit.N*EPSILON, wi ) ) ) {
    		vec3 brdf = brdf_evaluate( hit.N, hit.E, wi, hit.materialId, surfMtl.diff_, surfMtl.spec_ );
    		Lo += (Li*brdf)/pWi;
        }
    }   
    return Lo;
}
vec3 mirro_surface (Ray r) {
	vec3 color = vec3 (0);
	
	float depth = 0.0;
	
	float l = 0.0;
	
	for (int i = 0; i < 3; i++) {
		RayHit h;int objId;float dist;
		//raySceneIntersection (r,h, objId, dist);
		if(!raySceneIntersection (r,h, objId, dist))
		return vec3(0.0);
		else{
		float c = dot (h.N, r.dir);
		depth += (1.0 / 100000.0) + dist;
		
		r = Ray (h.pos + h.N * 0.0001, reflect (normalize (r.dir), h.N));
		
		float d = 1.0 / (depth * depth);
		
		Material surfMtl = getMaterialFromLibrary( h.materialId );
		vec3 hit_color=surfMtl.diff_;
		color = (color + c*c*d) * (1.0 - hit_color * d * 100000.0)/20.0;
		}
	}
	
	return color;
}
vec3 Radiance( in Ray ray, float rand1 ) 
{
    vec3 color = vec3( 0.0 );
    Ray currentRay = ray;
    vec3 weight = vec3( 1.0 );
    
    for( int i=0; i<1 ; i++ ) 
	{
        RayHit hit;
        int objId;
        float dist = INFINITY;
        
        if( raySceneIntersection( currentRay,hit, objId, dist ) ) {
            Material surfMtl = getMaterialFromLibrary( hit.materialId );
           if( hit.materialId == WHITE_DIFFUSE ) {
                color += (i==0)?surfMtl.diff_:vec3(0.0);
                break;
            } 
			else if(hit.materialId ==MTL_MIRRO)
			{
			  color+=mirro_surface (ray);//return color;//mirro_surface (ray);return color;
			 
			  break;
			}
            else {
                vec3 col = weight*calcDirectLightOnSurface( hit, surfMtl,ray);
                if( i==0 ) {
            		color += col;
                } 
                else {
                    float value = max(col.x,max(col.y,col.z));
                    if( value > INDIRECT_CLAMP ) {
                        col *= INDIRECT_CLAMP/value;
                    }
                    color += col;
                }
            }
            
            
            vec3 L;
            vec3 eval = brdf_sample ( hit.N, hit.E, L, hit.materialId, surfMtl.diff_, surfMtl.spec_, rnd(), rnd() );
            
            float dotNL = dot( hit.N, L );
          
            
            weight *= eval*dotNL;

          
            currentRay.origin = hit.pos + hit.N*EPSILON;
            currentRay.dir = L;
        } else {
            break;
        }
   
        }
    return color;

}

void main( )
{

    vec3 iResolution;
	iResolution.xy=gl_FragCoord.xy/v_uv.xy;
    seed = iResolution.y * gl_FragCoord.x / iResolution.x + gl_FragCoord.y / iResolution.y;
    float fov = radians(45.0);
    initMaterial();
    initScene();
	vec3 accumulatedColor = vec3( 0.0 );
    vec2 p = -1.0 + 2.0 * (gl_FragCoord.xy) / iResolution.xy;
    p.x *= iResolution.x/iResolution.y*1.7;
	p.y *= iResolution.y/iResolution.x*1.7;
	p.y-=0.3;
    Ray ray;
    
	for( int i=0; i<SAMPLES; i++){
        initLightSphere( iGlobalTime ); 
        ray = shootRay(p, rnd()/float(SAMPLES), rnd()/float(SAMPLES) );
        accumulatedColor += Radiance( ray, ( float(i) + rnd() )/float(SAMPLES) );
	}
	accumulatedColor = accumulatedColor/float(SAMPLES);
	if(if_gamma>0.0)
    accumulatedColor = pow( accumulatedColor, vec3( 1.0 / GAMMA ) );

	gl_FragColor = vec4( accumulatedColor,1.0 );
	//gl_FragColor = vec4(1.0,0.0,0.0,1.0 );
}