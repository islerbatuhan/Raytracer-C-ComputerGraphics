/////////////////////////////////////////////////////////////////////////////
////////////////// Batuhan Ýþler/1922020
////////////////// recursive raytracer
////////////////// light source(30,10,10)
////////////////// ambient = 10, intensity = 200
////////////////////////////////////////////////////////////////////////////
















#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define ABS(a) ((a)>0?(a):-1*(a))

typedef struct
{
	double x,y,z;
} Vec3;

Vec3 light;

typedef struct
{
	double r,g,b;
} Color;

typedef struct
{
	Vec3 a;
	Vec3 b;
} Ray;

typedef struct
{
	Vec3 center;
	double r;
	Color color;
} Sphere;

typedef struct
{
	Vec3 a,b,c;
	Color color;
} Triangle;

int numSpheres;
int numTriangles;

Sphere *spheres;
Triangle *triangles;

typedef struct 
{
	Vec3 pos;
	Vec3 gaze;
	Vec3 v;
	Vec3 u;
	double l,r,b,t;
	double d;
} Camera;



Camera cam;

int sizeX, sizeY;

double pixelW, pixelH;

double halfPixelW, halfPixelH;

char outFileName[80];

Color **image;

Vec3 cross(Vec3 a, Vec3 b)
{
	Vec3 tmp;
	
	tmp.x = a.y*b.z-b.y*a.z;
	tmp.y = b.x*a.z-a.x*b.z;
	tmp.z = a.x*b.y-b.x*a.y;
	
	return tmp;
}

double dot(Vec3 a, Vec3 b)
{
		return a.x*b.x+a.y*b.y+a.z*b.z;
}

double length2(Vec3 v)
{
	return (v.x*v.x+v.y*v.y+v.z*v.z);
}

double length(Vec3 v)
{
	return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

Vec3 normalize(Vec3 v)
{
	Vec3 tmp;
	double d;
	
	d=length(v);
	tmp.x = v.x/d;
	tmp.y = v.y/d;
	tmp.z = v.z/d;
	
	return tmp;
}

Vec3 add(Vec3 a, Vec3 b)
{
	Vec3 tmp;
	tmp.x = a.x+b.x;
	tmp.y = a.y+b.y;
	tmp.z = a.z+b.z;
	
	return tmp;
}

Vec3 mult(Vec3 a, double c)
{
	Vec3 tmp;
	tmp.x = a.x*c;
	tmp.y = a.y*c;
	tmp.z = a.z*c;
	
	return tmp;
}

double distance(Vec3 a, Vec3 b)
{
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

int equal(Vec3 a, Vec3 b)
{
	double e = 0.000000001;
	
	//printf("%lf %lf %f ----",ABS((a.x-b.x)),ABS((a.y-b.y)),ABS((a.z-b.z)));
	if ((ABS((a.x-b.x))<e) && (ABS((a.y-b.y))<e) && (ABS((a.z-b.z))<e))
		{ return 1;}
	else { return 0;}
}

void readScene(char *fName)
{
    FILE *f;
    int i;
    char line[80]={};
    f = fopen(fName,"r");
    
    fscanf(f,"%d",&numSpheres);
    spheres = (Sphere *)malloc(sizeof(Sphere)*numSpheres);

    for (i=0;i<numSpheres;i++)
    {
	    fscanf(f,"%lf %lf %lf", &(spheres[i].center.x),&(spheres[i].center.y),&(spheres[i].center.z));
	    fscanf(f,"%lf", &(spheres[i].r));
	    fscanf(f,"%lf %lf %lf", &(spheres[i].color.r),&(spheres[i].color.g),&(spheres[i].color.b));
    }

    fscanf(f,"%d",&numTriangles);
    triangles = (Triangle *)malloc(sizeof(Triangle)*numTriangles);

    for (i=0;i<numTriangles;i++)
    {
	    fscanf(f,"%lf %lf %lf", &(triangles[i].a.x),&(triangles[i].a.y),&(triangles[i].a.z));
	    fscanf(f,"%lf %lf %lf", &(triangles[i].b.x),&(triangles[i].b.y),&(triangles[i].b.z));
	    fscanf(f,"%lf %lf %lf", &(triangles[i].c.x),&(triangles[i].c.y),&(triangles[i].c.z));
	    fscanf(f,"%lf %lf %lf", &(triangles[i].color.r),&(triangles[i].color.g),&(triangles[i].color.b));
    }
}

void readCamera(char *fName)
{
    FILE *f;
    int i;
    char line[80]={};
    f = fopen(fName,"r");
    
    fscanf(f,"%lf %lf %lf",&cam.pos.x,&cam.pos.y,&cam.pos.z);
    fscanf(f,"%lf %lf %lf",&cam.gaze.x,&cam.gaze.y,&cam.gaze.z);
    fscanf(f,"%lf %lf %lf",&cam.v.x,&cam.v.y,&cam.v.z);
    
    cam.u = cross(cam.gaze,cam.v);
    
    cam.u = normalize(cam.u);
    
    cam.v = cross(cam.u,cam.gaze);
    cam.v = normalize(cam.v);
    
    cam.gaze = normalize(cam.gaze);
    
    fscanf(f,"%lf %lf %lf %lf",&cam.l,&cam.r,&cam.b,&cam.t);
    fscanf(f,"%lf",&cam.d);
    fscanf(f,"%d %d",&sizeX,&sizeY);
    
    fscanf(f,"%s",outFileName);
    
    // camera reading complete
    
    pixelW = (cam.r-cam.l)/(double)sizeX;
    halfPixelW = pixelW*0.5;
    
    pixelH = (cam.t-cam.b)/(double)sizeY;
    halfPixelH = pixelH*0.5;
    
    image = (Color **)malloc(sizeof(Color *)*sizeX);
    
    if (image==NULL)
    {
        printf("Cannot allocate memory for image.");
        exit(1);
    }
    
    for (i=0;i<sizeX;i++)
    {
        image[i] = (Color *)malloc(sizeof(Color)*sizeY);
        if (image[i]==NULL)
        {
            printf("Cannot allocate memory for image.");
            exit(1);
        }
    }
}

void initImage()
{
	int i,j;

    for (i=0;i<sizeX;i++)
		for (j=0;j<sizeY;j++)
		{
			image[i][j].r = 0;
			image[i][j].g = 0;
			image[i][j].b = 0;
		}
}

Ray generateRay(int i, int j)
{
	Ray tmp;
	
	Vec3 su,sv,s;
	
	tmp.a = cam.pos;
	
	su = mult(cam.u,cam.l+(i*pixelW)+halfPixelW);
	sv = mult(cam.v,cam.b+(j*pixelH)+halfPixelH);
	
	s = add(su,sv);
	
	tmp.b = add(mult(cam.gaze,cam.d),s);
    
	//printVec(add(tmp.a,tmp.b));
	
	return tmp;
}

    
double intersectSphere(Ray r, Sphere s)
{
	double A,B,C; //constants for the quadratic function
	
	double delta;
	
	Vec3 scenter;
    double sradius;
    
    Vec3 p;
		
	double t,t1,t2;
    int i;
	
    scenter = s.center;
    sradius = s.r;
    
    
	C = (r.a.x-scenter.x)*(r.a.x-scenter.x)+(r.a.y-scenter.y)*(r.a.y-scenter.y)+(r.a.z-scenter.z)*(r.a.z-scenter.z)-sradius*sradius;

	B = 2*r.b.x*(r.a.x-scenter.x)+2*r.b.y*(r.a.y-scenter.y)+2*r.b.z*(r.a.z-scenter.z);
	
	A = r.b.x*r.b.x+r.b.y*r.b.y+r.b.z*r.b.z;
	
	delta = B*B-4*A*C;
	
	if (delta<0) return -1;
	else if (delta==0)
	{
		t = -B / (2*A);
	}
	else
	{
		double tmp;
		delta = sqrt(delta);
		A = 2*A;
		t1 = (-B + delta) / A;
		t2 = (-B - delta) / A;
	
		if (t2<t1){
			tmp = t2;
			t2 = t1;
			t1 = tmp;
		}	
	
		t = t1;	
	}
	return t;
}

double intersectTriangle(Ray r, Triangle tr)
{
	double  a,b,c,d,e,f,g,h,i,j,k,l;
	double beta,gamma,t;
	
	double eimhf,gfmdi,dhmeg,akmjb,jcmal,blmkc;

	double M;
	
	double dd;
	Vec3 ma,mb,mc;

	ma = tr.a;
	mb = tr.b;
	mc = tr.c;
	
	a = ma.x-mb.x;
	b = ma.y-mb.y;
	c = ma.z-mb.z;

	d = ma.x-mc.x;
	e = ma.y-mc.y;
	f = ma.z-mc.z;
	
	g = r.b.x;
	h = r.b.y;
	i = r.b.z;
	
	j = ma.x-r.a.x;
	k = ma.y-r.a.y;
	l = ma.z-r.a.z;
	
	eimhf = e*i-h*f;
	gfmdi = g*f-d*i;
	dhmeg = d*h-e*g;
	akmjb = a*k-j*b;
	jcmal = j*c-a*l;
	blmkc = b*l-k*c;

	M = a*eimhf+b*gfmdi+c*dhmeg;
    if (M==0) return -1;
	
	t = -(f*akmjb+e*jcmal+d*blmkc)/M;
	
	gamma = (i*akmjb+h*jcmal+g*blmkc)/M;
	
	if (gamma<0 || gamma>1) return -1;
	
	beta = (j*eimhf+k*gfmdi+l*dhmeg)/M;
	
	if (beta<0 || beta>(1-gamma)) return -1;
	
	return t;
}

int convert (double d)
{
	if (d>255) return 255;
	if (d<0) return 0;
	return round(d);
}

void writeImage(char *fileName)
{
	FILE *outFile;
	int i,j;
	
	outFile = fopen(fileName,"w");
	
	fprintf(outFile,"P3\n");
	fprintf(outFile,"# %s\n",outFileName);
	
	fprintf(outFile,"%d %d\n",sizeX,sizeY);
	
	fprintf(outFile,"255\n");
	
	for (j=sizeY-1;j>=0;j--)
	{
		for (i=0;i<sizeX;i++)
		{
			fprintf(outFile,"%d %d %d ",convert(image[i][j].r),convert(image[i][j].g),convert(image[i][j].b));
		}
		fprintf(outFile,"\n");
	}
	fclose(outFile);
}






int ambient = 10;
int intensity = 200;



Vec3 ReflectRay(Vec3 r,Vec3 n){
	return add( mult( mult( n,dot(r,n)) , 2) , mult(r,-1) ) ; ;
}


Color rayrecursive(Ray r,double tmin,int i,int j,int depth){
	
	Color result;
	result.b=0;result.g=0;result.r=0;
	
	
	if(depth == -1){
		return result;
		
	}
	depth--;
	
	
	int closestObj=-1;
	int isSphere=0;
	
	int k;
	
	
	
	
	
	for (k=0;k<numSpheres;k++)	//foreach object o(sphere)
		    {
			    double t;

			    t = intersectSphere(r,spheres[k]); //if ray intersect o

			    if (t>=1)
			    {
				    if (t<tmin)
				    {
					    tmin = t;		//tmin = t, obj = o
					    closestObj = k;
					    isSphere = 1;
				    }
			    }
		    }
		    for (k=0;k<numTriangles;k++)	//foreach object o(triangle)
		    {
			    double t;

			    t = intersectTriangle(r,triangles[k]);	//if ray intersect o(triangle)

			    if (t>=1)
			    {
				    if (t<tmin)
				    {
					    tmin = t;		//tmin = t, obj = o
					    closestObj = k;
					    isSphere = 0;
				    }
			    }
		    }
		    if (closestObj!=-1)		//if obj not NULL
		    {
			    	Vec3 p;
				Vec3 normal;
				Vec3 toLight;
				Vec3 half;
				Vec3 view;
				Ray sr;
				double diffuse;
				double specular;
				int shadow = 0;

				p = add(r.a,mult(r.b,tmin));

				sr.a = p;


				if (isSphere)
					normal = add(p,mult(spheres[closestObj].center,-1));
				else
					normal = cross(add(triangles[closestObj].b,mult(triangles[closestObj].a,-1)),add(triangles[closestObj].c,mult(triangles[closestObj].a,-1)));

				normal = normalize(normal);

				toLight = add(light,mult(p,-1));

				sr.b = toLight;		//computed shadow ray sr
				
				
		    for (k=0;k<numSpheres;k++)		//for each object p(sphere)
		    {
			    double t;

			    t = intersectSphere(sr,spheres[k]);	//?sr intersects p(sphere)


			    if (t>0.0001 && t<1)		// if sr intersects p(sphere)
			    {
				shadow = 1;
				break;
										//return pixel color(points in shadow)
			    }
		    }
		    if (!shadow) {
		    for (k=0;k<numTriangles;k++)		//for each object p(triangle)
		    {
			    double t;

			    t = intersectTriangle(sr,triangles[k]);		//?sr intersects p(tri)

			    if (t>0.0001 && t<1)		//if sr intersects p(triangle)
			    {
				shadow = 1;					
				break;
										//return pixel color(points in shadow)
			    }
		    }
		    }

			if (!shadow) {					//(part for) add diffuse and specular components

				toLight = normalize(toLight);

				view = mult(r.b,-1);
				
				view = normalize(view);

				half = add(view,toLight);

				half = normalize(half);

				specular = dot(half,normal);

				if (specular<0) specular = 0;
				else
					specular = pow(specular,100);

				diffuse = dot(normal,toLight);
				if (diffuse<0) diffuse=0;
				
				
				
				////////deneme reflectray////
				
				
				///////version1//////////////////////////
				Ray reflectRay;
				//reflectRay.a=r.a;
				//reflectRay.b = add( mult(view,-1) , mult(normal,2*dot(normal,view)) );
				///////////////////////////////////////
				
				
				///////verison2/////////////////
				reflectRay.a=p;
				reflectRay.b = add( mult(view,-1) , mult(normal,2*dot(normal,view)) );
				////////////////////////////////////////////////
				
				
				
				/////version3///////////
				//reflectRay.a=r.a;
				//reflectRay.b= add( mult(r.b,-1) , mult(normal,2*dot(normal,r.b)) );
				/////////////////////////////
				
				
				
				

			if (isSphere)
			{
				result.r= spheres[closestObj].color.r*intensity*diffuse+ambient*spheres[closestObj].color.r+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).r;
				result.g= spheres[closestObj].color.g*intensity*diffuse+ambient*spheres[closestObj].color.g+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).g;
				result.b= spheres[closestObj].color.b*intensity*diffuse+ambient*spheres[closestObj].color.b+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).b;
				
		    //image[i][j].r = spheres[closestObj].color.r*intensity*diffuse+ambient*spheres[closestObj].color.r+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).r;
		    //image[i][j].g = spheres[closestObj].color.g*intensity*diffuse+ambient*spheres[closestObj].color.g+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).g;
		    //image[i][j].b = spheres[closestObj].color.b*intensity*diffuse+ambient*spheres[closestObj].color.b+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).b;
			}
			else
			{
				result.r= triangles[closestObj].color.r*intensity*diffuse+ambient*triangles[closestObj].color.r+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).r;
				result.g= triangles[closestObj].color.g*intensity*diffuse+ambient*triangles[closestObj].color.g+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).g;
				result.b= triangles[closestObj].color.b*intensity*diffuse+ambient*triangles[closestObj].color.b+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).b;
		    //image[i][j].r = triangles[closestObj].color.r*intensity*diffuse+ambient*triangles[closestObj].color.r+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).r;
		    //image[i][j].g = triangles[closestObj].color.g*intensity*diffuse+ambient*triangles[closestObj].color.g+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).g;
		    //image[i][j].b = triangles[closestObj].color.b*intensity*diffuse+ambient*triangles[closestObj].color.b+intensity*specular + 0.3*rayrecursive(reflectRay,40000,i,j,depth).b;
			}
		    }

		   else
			{
				
				
				///////version1//////////////////////////
				Ray reflectRay;
			//	reflectRay.a=r.a;
			//	reflectRay.b = add( mult(view,-1) , mult(normal,2*dot(normal,view)) );
				///////////////////////////////////////
				
				
				///////verison2/////////////////
				reflectRay.a = p;
				reflectRay.b = add( mult(view,-1) , mult(normal,2*dot(normal,view)) );
				////////////////////////////////////////////////
				
				
				///////version3//////////////
			//	reflectRay.a=r.a;
			//	reflectRay.b= add( mult(r.b,-1) , mult(normal,2*dot(normal,r.b)) );
				//////////////////////////////////
				
				
				
				
				
			if (isSphere)
			{
				result.r= ambient*spheres[closestObj].color.r +0.3*rayrecursive(reflectRay,40000,i,j,depth).r;
				result.g= ambient*spheres[closestObj].color.g +0.3*rayrecursive(reflectRay,40000,i,j,depth).g;
				result.b= ambient*spheres[closestObj].color.b +0.3*rayrecursive(reflectRay,40000,i,j,depth).b;
		    //image[i][j].r = ambient*spheres[closestObj].color.r +0.3*rayrecursive(reflectRay,40000,i,j,depth).r;
		    //image[i][j].g = ambient*spheres[closestObj].color.g +0.3*rayrecursive(reflectRay,40000,i,j,depth).g;
		    //image[i][j].b = ambient*spheres[closestObj].color.b +0.3*rayrecursive(reflectRay,40000,i,j,depth).b;
		    
		    
			}
			else
			{
				result.r= ambient*triangles[closestObj].color.r +0.3*rayrecursive(reflectRay,40000,i,j,depth).r;
				result.g= ambient*triangles[closestObj].color.g +0.3*rayrecursive(reflectRay,40000,i,j,depth).g;
				result.b= ambient*triangles[closestObj].color.b +0.3*rayrecursive(reflectRay,40000,i,j,depth).b;
		    //image[i][j].r = ambient*triangles[closestObj].color.r +0.3*rayrecursive(reflectRay,40000,i,j,depth).r;
		    //image[i][j].g = ambient*triangles[closestObj].color.g +0.3*rayrecursive(reflectRay,40000,i,j,depth).g;
		    //image[i][j].b = ambient*triangles[closestObj].color.b +0.3*rayrecursive(reflectRay,40000,i,j,depth).b;
			}
			} 
		}
	
	
	return result;
	
}













int main(int argc, char **argv)
{
    
    char comm[100];
    int i,j,k;
    if (argc<3) {
	printf("Usage: raytracer <camera file> <scene file>\n");
	return 1;
    }
    readCamera(argv[1]);
    readScene(argv[2]);
    
    initImage();

    light.x = 30;
    light.y = 10;
    light.z = 10;

    for (i=0;i<sizeX;i++)	// start foreachpixel
	    for (j=0;j<sizeY;j++) //foreach pixel
	    {
		    Ray r;
		    double tmin = 40000;
		    //int closestObj = -1;
		    //int isSphere = 0;
		    int depth = 1;
		    
		    
		    
		    
		    

		    r = generateRay(i,j);	//compute the viewing ray r(from e to s)
			
		    //rayrecursive(r,tmin,closestObj,isSphere,i,j,k,depth);
		    Color tmp= rayrecursive(r,tmin,i,j,depth);
		    image[i][j].r = tmp.r;
		    image[i][j].g= tmp.g;
		    image[i][j].b= tmp.b;
		    
		    
		       
		    
		    
		    
		    
		    
		}

    writeImage("batuhanIsler.ppm");
    
    return 0;
}
