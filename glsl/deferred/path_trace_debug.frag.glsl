#version 100
precision highp float;
precision highp int;
//----------choose the thing u want to sample--
//#define SAMPLE_BRDF
#define COMBINE
//#define SAMPLE_LIGHT
//------------------------------

//#define SCENE_DEFAULT
#define SCENE_SMPL_LIHGT
//#define SCENE_SMPL_BRDF
#define LIGHTCOLOR          vec3(1.0,1.0,1.0)
#define BLINNCOLOR          vec3( 0.1,0.1,0.7 )
#define SAMPLES 		    1
#define N                   5
#define EPSILON 			0.001 
#define GAMMA 				2.2			
#define DEPTH		        5
#define INDIRECT_CLAMP 		20.0	
#define INFINITY            9999.9
#define FOG_STRENGTH 1.0
#define BLACK               vec3( 0.0 )
#define PI 					3.1415926
#define TWO_PI 				6.2831852
#define REFLECTION_STRENGTH 100000.0
float INTENSITY=   50.0;

#define OBJ_NULL 0
#define OBJ_CUBE 1
#define OBJ_SPHERE 2
#define OBJ_PLANE 3

#define COLOR_STRENGTH 10.0
//varying vec2 v_uv;

uniform float u_intensity;
uniform float iGlobalTime;

varying vec2 v_uv;
float seed;	//seed initialized in main
//reference:https://www.shadertoy.com/view/MtfGR4
float rand1() { return fract(sin(seed++)*43758.5453123); }
vec2  rand2() {
    return fract(sin(vec2(seed+=0.1,seed+=0.1))*vec2(43758.5453123,22578.1459123));
}
struct Ray {vec3 origin; vec3 dir;};
struct Intersection
{
    vec3 point;
    vec3 normal;
    float t;
    int hit_obj;
    int hit_mat;
   // vec3 hit_color;
};
//-------------geometry-----------------//
struct Sphere { int materialId; vec3 pos; float r;  float area; };
struct Plane { int materialId; vec4 abcd; };
struct Box { vec3 Bl; vec3 Bh; int materialId; };
//----------------------------------------------    
struct Material { vec3 base; float spec; };
struct Camera { mat3 rotate; vec3 pos; vec3 target; float fovV; float lensSize; float focusDist; };


    
//------------ SCENE--------------
Box boxes[5];
Plane walls[5];
Sphere spherelight;
Sphere spheres[5];
Camera camera;
vec3 light_color=vec3(1.0);
//data structure:https://www.shadertoy.com/view/lts3Dr
//-----------------------------------
//--------------------material-------
//-----------------------------------
//---1)diffuse-----------------------
#define M_DIF_WHITE	        0
#define M_DIF_BLACK 		1
#define M_DIF_GREY	        2
#define M_DIF_RED   		3
#define M_DIF_GREEN 		4
//---2)blinn-------------------------
#define M_BLINN_1	        5
#define M_BLINN_2           6
#define M_BLINN_3           7
#define M_BLINN_4           8
//----3)mirro------------------------
#define M_MIRRO             9
#define M_LIGHT             10
#define M_COUNT			    11

Material materials[M_COUNT];
//interesting... cant use A.base=base;
#define SET_MAT(i,_diff,_spec) { materials[i].base=_diff; materials[i].spec=_spec; }

void initMaterial() {
    float intensity = 50.0 ;
  
    SET_MAT( M_DIF_WHITE, vec3( 1.0 ), 0.0 );
    SET_MAT( M_DIF_BLACK, vec3( 0.0 ),  0.2 );
    SET_MAT( M_DIF_GREY, vec3( 0.5, 0.6, 0.8 ),  0.8 );
    SET_MAT( M_DIF_RED, vec3( 1.0, 0.1, 0.15 ),  0.2 );
    SET_MAT( M_DIF_GREEN, vec3( 0.1, 0.8, 0.1 ), 0.2 );
   
    
    SET_MAT( M_BLINN_1,vec3( 0.1,0.1,0.7 ), 10.0 );
    SET_MAT( M_BLINN_2,vec3( 0.1,0.1,0.7 ), 20.0 );
    SET_MAT( M_BLINN_3,vec3( 0.1,0.1,0.7 ), 50.0 );
    SET_MAT( M_BLINN_4,vec3( 0.1,0.1,0.7 ), 100.0 );
    
    SET_MAT( M_MIRRO, vec3( 0.7,0.1,0.4 ),  1.0 );
    SET_MAT( M_LIGHT, vec3( 1.0 ), INTENSITY );
}

Material GetMaterial(int index ) 
{   
    for(int i=0; i<M_COUNT; i++ ) { 
    	if( index == i ) {
            return materials[i];
        }
    }  
    return materials[0];
}
void initLightSphere( float time ) {
	spherelight.pos = vec3( 3.5*sin(time*0.9),2.0,-2.0);
}

void initScene() 
{

    float time = iGlobalTime;

    //lights
    float rlgiht=0.5;
    float r = 0.5;
    
    spherelight = Sphere( M_LIGHT, vec3( 2.0, 2.0, -2.0 ), rlgiht, rlgiht*rlgiht*4.0*PI );
    SET_MAT( M_LIGHT, vec3( 1.0 ), 50.0 );
    
   // int scene=0;
#ifdef SCENE_DEFAULT
   // if(scene==0){
    //cornelbox scene---
    //walls
    //bottom adc:normal
    walls[0].abcd = vec4( 0.0, 1.0, 0.0, 1.0 );
    walls[0].materialId = M_DIF_GREY;
  
    //top
    walls[1].abcd = vec4( 0.0, -1.0, 0.0, 5.5 );
    walls[1].materialId = M_DIF_BLACK;
    
    //front
    walls[2].abcd = vec4( 0.0, 0.0, 1.0, 6.0 );
    walls[2].materialId = M_DIF_WHITE;
    
    //left
    walls[3].abcd = vec4( 1.0, 0.0, 0.0, 4.0 );
    walls[3].materialId = M_DIF_GREEN;
    
    //right
    walls[4].abcd = vec4( -1.0, 0.0, 0.0, 4.0 );
    walls[4].materialId = M_DIF_RED;
    
    //box
    boxes[0].Bl = vec3( 0.2, -1.0, -2.5 );
    boxes[0].Bh = vec3( 2.0, 0.8, -1.6 );
    boxes[0].materialId = M_BLINN_1;
    
    //sphere
    r = 0.5;
    spheres[1] = Sphere( M_BLINN_1, vec3( -2, 0.0, -3.0 ), r, r*r*4.0*PI  );
    spheres[2] = Sphere( M_BLINN_2, vec3( -0.5, 3.0, -3.0 ), r, r*r*4.0*PI  );
        
    spheres[3] = Sphere( M_BLINN_3, vec3( 3.5, 3.0, -3.0 ), r, r*r*4.0*PI  );
    spheres[4] = Sphere( M_BLINN_4, vec3( -2.5, 3.0, -3.0 ), r, r*r*4.0*PI  );
     
  //  }
#endif    
#ifdef SCENE_SMPL_LIHGT
  //  if(scene==1)
   // {
        //4 plane:
   //  walls[0].abcd = vec4( 0.0, 1.0, 0.0,2.0 );
  //   walls[0].materialId =M_DIF_GREY;
    walls[0].abcd = vec4( 0.0, 1.0, 0.0, 1.0 );
    walls[0].materialId = M_DIF_GREY;
  
    //top
   /* walls[1].abcd = vec4( 0.0, -1.0, 0.0, 5.5 );
    walls[1].materialId = M_DIF_BLACK;
    
    //front
    walls[2].abcd = vec4( 0.0, 0.0, 1.0, 6.0 );
    walls[2].materialId = M_DIF_WHITE;
    
    //left
    walls[3].abcd = vec4( 1.0, 0.0, 0.0, 4.0 );
    walls[3].materialId = M_DIF_GREEN;
    
    //right
    walls[4].abcd = vec4( -1.0, 0.0, 0.0, 4.0 );
    walls[4].materialId = M_DIF_RED;*/
    
    spheres[0] = Sphere( M_MIRRO, vec3( -0.5,1.0, -3.0 ), r, r*r*4.0*PI  );        
    spheres[1] = Sphere( M_BLINN_1, vec3( -2.0, 1.0, -3.0 ), r, r*r*4.0*PI  );
    spheres[2] = Sphere( M_BLINN_2, vec3( -3.5, 1.0, -3.0 ), r, r*r*4.0*PI  );  
        
    spheres[3] = Sphere( M_BLINN_3, vec3( 1.0, 1.0, -3.0 ), r, r*r*4.0*PI  );       
    spheres[4] = Sphere( M_BLINN_4, vec3( 3.0, 1.0, -3.0 ), r, r*r*4.0*PI  );
    
//  }
#endif
    
}

Intersection SetNull()
{
    Intersection isx;
    isx.t=-1.0;
    isx.hit_obj=OBJ_NULL;
    return isx;
}
// intersection function
Intersection raySphereIntersection(in Ray r, in Sphere sphr) 
{    
    vec3 os=r.origin-sphr.pos;
    //remember to normalize the r.direction
    float A=1.0;
    float B=2.0*(r.dir.x*os.x+r.dir.y*os.y+r.dir.z*os.z);
    float C = pow(os.x,2.0) + pow(os.y,2.0) +pow(os.z,2.0)-sphr.r*sphr.r;
    float Discre=pow(B,2.0)-4.0*A*C;
    
    if(Discre < 0.0) {
       
        return SetNull();
    }
    else{
        
        float t0= (-B - sqrt(Discre))/(2.0*A);
         
        if(t0>0.0)//intersecton
        {
            Intersection isx;
            isx.point = r.dir*t0+r.origin;//to Local
            isx.normal= normalize(isx.point-sphr.pos);//normalize!!!!!
            isx.t=t0;
            isx.hit_obj=OBJ_SPHERE;
            
            return isx;
        }
        else{
            float t1= (-B + sqrt(Discre))/(2.0*A);
            
            if(t1>0.0)//t1 intersect at front in sphere
            {
                 Intersection isx;
                 isx.point = r.dir*t1+r.origin;
                 isx.normal =normalize(isx.point-sphr.pos);//normalize!!!!!
                 isx.t = t1;
                 isx.hit_obj=OBJ_SPHERE;
                 
                 return isx;
            }
            else {
                return SetNull();
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
 
   Intersection isx=SetNull();

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
              isx.point=point;
              isx.normal=normal*_sign;
              isx.hit_obj=OBJ_CUBE;
              isx.t=Tnear;
                  return isx;
             // return Intersection(point,vec3(normal.x*_sign.x, normal.y*_sign.y, normal.z*_sign.z),Tnear,OBJ_CUBE);
          }         
          else
              return isx;
       }
}
Intersection rayPlaneIntersection( Ray ray, Plane plane)
{
    Intersection isx=SetNull();

    float dotVN = dot( ray.dir, plane.abcd.xyz );
   //dot(N,(S - R0)):S is random point on plane
    if ( abs( dotVN ) > EPSILON ) {
		isx.t = -(dot(plane.abcd.xyz, ray.origin) + plane.abcd.w)/dotVN;
        isx.hit_obj=OBJ_PLANE;
        isx.point=ray.origin+ray.dir*isx.t;
        isx.normal=plane.abcd.xyz;
    }
    
    return isx;
}


 
//-------Geometry functions--------
float UniformConePdf(float cosThetaMax)
{
    return 1.0 / (2.0 * PI * (1.0 - cosThetaMax));
}
vec3 localToWorld( in vec3 localDir, in vec3 normal )
{
    vec3 binormal = normalize( ( abs(normal.x) > abs(normal.z) )?vec3( -normal.y, normal.x, 0.0 ):vec3( 0.0, -normal.z, normal.y ) );
	vec3 tangent = cross( binormal, normal );
    
	return localDir.x*tangent + localDir.y*binormal + localDir.z*normal;
}
Intersection GetSphereSurfaceSample(in float u1,in float u2)
{
    float z = 1.0 - 2.0 * u1;
    float r = sqrt(max(0.0, 1.0 - z*z));
    float phi = 2.0 * PI * u2;
    float x = r * cos(phi);
    float y = r * sin(phi);
    vec3 normal = normalize(vec3(x,y,z));
    
    vec4 pointL=vec4( x/2.0, y/2.0, z/2.0, 1.0);
    vec3 ppoint=vec3(x/2.0, y/2.0, z/2.0);
    vec4 normalL=vec4(normal,0.0);
   
    vec3 T = normalize(cross(vec3(0,1,0),vec3(normalL)));
    vec3 B = cross(vec3(normalL), T);

    Intersection result;
   
    result.point = spherelight.pos+ppoint;
    result.normal = normalize(result.point-spherelight.pos);

    return result;
}
Ray shootRay( in vec2 pixel,vec2 mo )
{
    Ray ray;
    vec3 m=vec3(mo.x,mo.y,mo.x);
    vec3 eye = vec3(0, 2.0, 12.00);
    vec3 ref = vec3(0, 2.0,  0.00)+m;
    vec3 F = normalize( ref - eye );
    vec3 U = normalize( cross(F,vec3(0.0,1.0,0.0) ) );
    vec3 R = normalize( cross(U,F)); 
	//vec3 rd = normalize( (pixel.x+rand1)*U + (pixel.y+rand2)*R +4.* F );
    vec3 rd = normalize( pixel.x*U + pixel.y*R +4.* F );
    ray.origin=eye;
    ray.dir=rd;  
	return ray;
}

Intersection raySceneIntersection( in Ray ray)
{    
    float nearest_t = INFINITY;
    Intersection isx_f=SetNull();
    
    //check lights
    Intersection isx0=raySphereIntersection( ray, spherelight);
    if(isx0.t>0.0){
        nearest_t=isx0.t;
        isx_f=isx0; 
        isx_f.hit_mat=M_LIGHT;}
    
    //sphere  
    for(int i=0;i<5;i++)
    {
       Intersection isx_sphere=raySphereIntersection( ray, spheres[i]);
       float ts= isx_sphere.t; 
       if(ts<0.0)continue;
       else
       {
           if(ts<nearest_t) 
           {
               nearest_t=ts;
               isx_f=isx_sphere;
               isx_f.hit_mat=spheres[i].materialId;
           }
       }
    }    
     for( int i=0; i<5; i++ )
     {
       
        Intersection isx_p=rayPlaneIntersection(ray, walls[i]);
        float ps=isx_p.t;
        if(ps<0.0)continue;
        else
         {
             if(ps<nearest_t){
                 nearest_t=ps;
                 isx_f=isx_p;
                 isx_f.hit_mat=walls[i].materialId;
             }
         }
     }
    //cube
   
      Intersection isx_cube = rayCubeIntersection(ray,boxes[0]);
      float cs= isx_cube.t;
      if(cs<nearest_t&&cs>0.0)
      {
         nearest_t=cs;isx_f=isx_cube;
         isx_f.hit_mat=boxes[0].materialId;
      }
    
   
      return isx_f;
}

vec3 blin_microfacet(in vec3 normal,in vec3 wo,in vec3 wi, in vec3 diffuse, in float exponent)
{

   //------ Geometry term-------:
    
    normal=normalize(normal);
    vec3 wh=normalize(wo+wi);
    float a=abs(dot(wh,normal));
    float M=abs(2.0*dot(normal,wh)*dot(normal,wo)/dot(wo,wh));
    float S=abs(2.0*dot(normal,wh)*dot(normal,wi)/dot(wo,wh));
    float _G=min(min(M,S),1.0);


    //------Distribution term------
    //float exponent=10.0;
    float _D=( exponent+2.0)/(2.0*PI)*pow(a, exponent);

    float costheta_o=dot(normal,wo);//normalize!!!
    float costheta_i=dot(normal,wi);
    float brdf=_D*_G/(4.0*costheta_o*costheta_i);
    
    //----color--calculation--
    
    return brdf*diffuse;
}

bool ShadowTest( Ray shadowRay ) {
    float distToHit;
   
    Intersection isx=raySceneIntersection( shadowRay);
    if(isx.hit_mat==M_LIGHT)return true;
    else return false;
}
vec3 EvaluateLightBRDF(in vec3 wo, in vec3 wi,in vec3 normal, in int matId)
{
   Material mat=GetMaterial(matId) ;
   //diffuse
   if(matId<=4&&matId>=0){ return mat.base/PI;}
   if(matId>=5&&matId<=8)
    {
       return blin_microfacet(normal,wo, wi, mat.base,mat.spec);
    }
    if(matId==9) return BLACK;
    else return BLACK;
}

float LightRayPDF(in Intersection light_isx,in Ray light_ray)
{
    vec3 Pcenter = spherelight.pos;
    float radius = spherelight.r;
    float a_temp=distance(light_isx.point, Pcenter);
    float dis=pow(a_temp,2.0);
  //  return dis;
    // Return uniform weight if point inside sphere
    if (dis - radius*radius <EPSILON)
    {
      //if(light_isx.hit_obj == OBJ_NULL){return 0.0;}
        
      float cos_theta=dot(-light_ray.dir,light_isx.normal);
      return radius/(cos_theta*spherelight.area);
   }

    // Compute general sphere weight
    float sinThetaMax2 = radius*radius / dis;
    float cosThetaMax = sqrt(max(0.0, 1.0 - sinThetaMax2));
    return UniformConePdf(cosThetaMax);
}
vec3 SampleLambert(in Intersection isx,in vec3 ray_dir,out vec3 new_dir,out float brdf_pdf)
{
      vec2 rnd=rand2();
 
    
       float t = 2.0*PI*rnd.x;
       float u = rnd.y-0.5;//[-0.5,+0.5]*/

       float x=u*cos(t);
       float y=u*sin(t);
       float z=sqrt(0.25-pow(x,2.0)-pow(y,2.0));
       new_dir=normalize(vec3(x,y,z));
       brdf_pdf=dot(new_dir,isx.normal)/PI;
       return EvaluateLightBRDF(ray_dir, new_dir, isx.normal, isx.hit_mat);
     
}
bool samehemi(vec3 x, vec3 y)
{
    return x.z*y.z>0.0;
}
vec3 SampleBlinn(in Intersection isx,in vec3 ray_dir,out vec3 new_dir,out float brdf_pdf)
{
    vec2 rnd2=rand2();
    float rand1=rnd2.x;float rand2=rnd2.y;
    
    Material mat=GetMaterial(isx.hit_mat);
    float exponent=mat.spec;
    float costheta=pow(rand1,1.0/(exponent+1.0));
    float sintheta=sqrt(max(0.0,1.0-costheta*costheta));
    float phi=rand2*2.0*PI;
    
    vec3 wh;
    wh.x=cos(phi)*sintheta;
    wh.y=sin(phi)*sintheta;
    wh.z=costheta;
    wh=localToWorld( wh, isx.normal );
  //  if(!samehemi(ray_dir, wh))wh=-wh;

    new_dir=-ray_dir+2.0*dot(ray_dir,wh)*wh;
    
    vec3 result=EvaluateLightBRDF(ray_dir, new_dir, isx.normal, isx.hit_mat);

    vec3 wm=normalize(ray_dir+new_dir);
    float cst=abs(wm.z);
    brdf_pdf=(exponent+1.0)*pow(cst,exponent)/(2.0*PI*4.0*dot(ray_dir,wm));

    return result;
}
vec3 SampleMirro(in Intersection isx,in vec3 ray_dir,out vec3 new_dir,out float brdf_pdf)
{
   // new_dir=normalize(vec3(-ray_dir.x,-ray_dir.y,ray_dir.z));
  //  new_dir=localToWorld(new_dir, isx.normal );
    new_dir=normalize(reflect(-ray_dir,isx.normal));
    brdf_pdf=1.0;
    return BLINNCOLOR;//to make the compare color he same..so not using base_color
}
vec3 EvaluateSampleBRDF(in Intersection isx,in vec3 ray_dir,out vec3 new_dir,out float brdf_pdf)
{
    vec2 rnd=rand2();
    if(isx.hit_mat<=4)
    { 
        return SampleLambert(isx, ray_dir,new_dir, brdf_pdf);
    }
    if(isx.hit_mat<=8&&isx.hit_mat>=5)
    {
        return SampleBlinn(isx, ray_dir,new_dir, brdf_pdf);
    }
    if(isx.hit_mat==9)
    {
        return SampleMirro(isx, ray_dir,new_dir, brdf_pdf);
    }
    else return vec3(0.0);
}
vec3 Light_Radiance(in vec3 dir, in vec3 normal)
{
   return dot(dir,normal)>0.0? INTENSITY*vec3(1,1,1):vec3(0.0);
}
vec3 Directlight( in Ray ray) {
    
    vec3 color = vec3( 0.0 );
    vec3 weight = vec3( 1.0 );
    //for(int i=0;i<2;i++){
    Intersection isx=raySceneIntersection(ray);
    if(isx.t<0.0) {return vec3(0.0,0.0,0.0);}
    else
    {
        if(isx.hit_mat==M_LIGHT)return LIGHTCOLOR;
        //sample light
      
        vec2 rnd2=rand2();
        vec3 light_color;
        vec3 brdf_color;
        Intersection light_isx= GetSphereSurfaceSample(rnd2.x,rnd2.y);
        
        Ray light_ray;
        light_ray.origin=light_isx.point;
        light_ray.dir=normalize(light_isx.point-isx.point);
        light_ray.origin += light_ray.dir*EPSILON;
        
        //if(ShadowTest(light_ray))
        //{
           //EvaluateScatteredEnergy of light matrial
            vec3 radiance= Light_Radiance(-light_ray.dir,light_isx.normal);
            //vec3 radiance=dot(-light_ray.dir,light_isx.normal)>0.0?50.0*vec3(1,1,1):vec3(0.0);
            float light_pdf=LightRayPDF(light_isx,light_ray);
            vec3 light_brdf=EvaluateLightBRDF(-ray.dir,light_ray.dir,isx.normal,isx.hit_mat);
            vec3 sample_light_color=radiance*light_brdf*abs(dot(isx.normal,light_ray.dir))/light_pdf;
            light_color=light_brdf;
            
        //}
        
        Ray new_ray;
        new_ray.origin=isx.point;
        float brdf_pdf;
        vec3 brdf_brdf=EvaluateSampleBRDF(isx,-ray.dir,new_ray.dir,brdf_pdf);
        //update ray
        new_ray.dir=normalize(new_ray.dir);
        new_ray.origin+=new_ray.dir*EPSILON;
        
        
        //------
        vec3 Ld;
        Intersection isx_new=raySceneIntersection(new_ray);
        if(isx_new.t<0.0){Ld=vec3(0.0);}
        else if(isx_new.hit_mat==M_LIGHT){Ld=vec3(1.0);}//Ld= Light_Radiance(new_ray.dir,isx.normal);}
        else Ld=vec3(0.0);
    
        vec3 test_color=Ld;
       // vec3 test_color=vec3(brdf_pdf,0.0,0.0);
        if(abs(brdf_pdf)>0.0)
        {
            float absd=abs(dot(isx.normal,new_ray.dir));
            vec3 sample_brdf_color=Ld*brdf_brdf*absd/brdf_pdf;
            brdf_color=sample_brdf_color;
        }
        
         float w_light=pow(light_pdf,2.0)/(pow(light_pdf,2.0)+pow(brdf_pdf,2.0));
         float w_brdf=pow(brdf_pdf,2.0)/(pow(light_pdf,2.0)+pow(brdf_pdf,2.0));
        
         vec3 com_color =(light_color*w_light+brdf_color*w_brdf);
        #ifdef SAMPLE_BRDF
         return brdf_color;
        #endif
        #ifdef SAMPLE_LIGHT
        return light_color;
        #endif
        #ifdef COMBINE
        return com_color;
        #endif
        return com_color;
       // return light_color;
        
  //   }
   
       
     }

}       
 

void main()
{
	vec3 iResolution;
	iResolution.xy=gl_FragCoord.xy/v_uv.xy;
    seed = iResolution.y * gl_FragCoord.x / iResolution.x + gl_FragCoord.y / iResolution.y;
    float fov = radians(45.0);
    initMaterial();
    initScene();
    vec2 mo = vec2(0.0);//iMouse.xy/iResolution.xy;
	vec3 accumulatedColor = vec3( 0.0 );
    vec2 p = -1.0 + 2.0 * (gl_FragCoord.xy) / iResolution.xy;
    p.x *= iResolution.x/iResolution.y;//*1.7;
	p.y *= iResolution.y/iResolution.x;//*1.7;
	p.y-=0.3;
    Ray ray;
    
	for( int i=0; i<SAMPLES; i++){
        
        initLightSphere( iGlobalTime );
        ray = shootRay(p,mo);
    //    ray = shootRay(p, rnd()/float(SAMPLES), rnd()/float(SAMPLES) );
        accumulatedColor += Directlight( ray );
       // accumulatedColor += Directlight( ray, ( float(i) + rnd() )/float(SAMPLES) );
	}
	accumulatedColor = accumulatedColor/float(SAMPLES);
	//if(if_gamma>0.0)
    accumulatedColor = pow( accumulatedColor, vec3( 1.0 / GAMMA ) );

	gl_FragColor = vec4( accumulatedColor,1.0 );
	
}