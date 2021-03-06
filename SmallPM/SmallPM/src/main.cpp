/*********************************************************************************
Copyright (C) 2014 Adrian Jarabo (ajarabo@unizar.es)
Copyright (C) 2014 Diego Gutierrez (diegog@unizar.es)
All rights reserved.

This is an educational Ray Tracer developed for the course 'Informatica Grafica'
(Computer Graphics) tought at Universidad de Zaragoza (Spain). As such, it does not 
intend to be fast or general, but just to provide an educational tool for undergraduate
students. 

This software is provided as is, and any express or implied warranties are disclaimed.
In no event shall copyright holders be liable for any damage.
**********************************************************************************/

#include <iostream>
#include "globals.h"

#include "World.h"
#include "RenderEngine.h"
#include "PhotonMapping.h"
#include "Film.h"
#include "Sphere.h"
#include "Triangle.h"
#include "Plane.h"
#include "Mesh.h"
#include "BSDF.h"
#include "Lambertian.h"
#include "Specular.h"
#include "Transmissive.h"
#include <time.h>

using namespace std;

namespace {
	Film *film;
	RenderEngine *engine;
	PhotonMapping *pm;
	// Ajuste manual del tama�o de la imagen
	int sizex = 100, sizey = 100;
}


// main function - parse arguments and begin rendering.
int main(int argc, char* argv[])
{
	///////////// AJUSTE MANUAL DE VARIABLES //////////////////
	// Nombre
	char *name_file = NULL, *default_name_file = "name_file";
	name_file = default_name_file;

	// Escena
	unsigned int scene = 1;

	// Fotones
	unsigned int photons_global = 10000, 
				 photons_caustic = 10000, 
				 photons_volumetric = 10000,
				 photons_causticvol = 10000,
				 max_shots = 100000, 
				 nb_nearest_photons = 55;

	// Scattering!
	double sigmaScattering = 0.5;
	double sigmaAbsorcion = 0.2;
	double lambda = 0.04;
	bool participative = true;


	//////////////////////////////////////////////////////////
	
	//srand(time(0));
	Real focal_distance = 2.6;

	// ---------------------------------------------------------------------
	// Parse input
	for( int i=1; i<argc; ++i)
	{
		if( !strcmp("-film-name", argv[i]) ) {++i; name_file = argv[i];}
		if( !strcmp("-film-size-x",argv[i]) ) {++i; sizex = atoi(argv[i]); }
		if( !strcmp("-film-size-y",argv[i]) ) {++i; sizey = atoi(argv[i]);}
		
		if (!strcmp("-scene", argv[i])) {++i; scene = atoi(argv[i]);}

		if( !strcmp("-pm-total-photons",argv[i]) ) {++i; photons_global = photons_caustic = atoi(argv[i]); }
		if( !strcmp("-pm-photons-global",argv[i]) ) {++i; photons_global = atoi(argv[i]); }
		if( !strcmp("-pm-photons-caustic",argv[i]) ) {++i; photons_caustic = atoi(argv[i]); }
		if( !strcmp("-pm-max-photons-shot",argv[i]) ) {++i; max_shots = atoi(argv[i]); }
		if( !strcmp("-pm-nb-nearest-neighbor",argv[i]) ) {++i; nb_nearest_photons = atoi(argv[i]); }

	}
  
	max_shots = (max_shots > (photons_global + photons_caustic))? max_shots:(photons_global + photons_caustic);
  
  	// ----------------------------------------------------------------------
	// Set up viewing parameters and create a camera.
	//
	Vector3 camera_position(0,1,3.0);
	Vector3 looking_at(0,1,0);
	Vector3 up(0,1,0);
	Real view_plane_dist = 1;
	Camera camera(camera_position,looking_at,up,view_plane_dist);


	// ----------------------------------------------------------------------
	// Create a world 
	//
	World*  w = new World;
  
	w->set_ambient(Vector3(0.1,0.1,0.1));

	BSDF* glass = new Transmissive(w, 1.5);
	BSDF* mirror = new Specular(w);

	BSDF* white = new Lambertian(w, Vector3(.85,.85,.85));
	BSDF* red = new Lambertian(w, Vector3(.85,.085,.085));
	BSDF* green = new Lambertian(w, Vector3(.085,.85,.085));

	// Suelo
	Triangle* floor1 = new Triangle( Vector3(-1.5,0,1.5),Vector3(1.5,0.,1.5),
									 Vector3(-1.5,0.,-1.5), white);
	w->add_object(floor1); 
	Triangle* floor2 = new Triangle(Vector3(-1.5,0,-1.5), Vector3(1.5,0.,1.5), 
		Vector3(1.5,0.,-1.5), white);
	w->add_object(floor2); 

	// Techo
	Triangle* ceil1 = new Triangle(Vector3(1.5,2.,1.5), Vector3(-1.5,2,1.5), 
		Vector3(-1.5,2.,-1.5), white);
	w->add_object(ceil1); 
	Triangle* ceil2 = new Triangle(Vector3(1.5,2.,1.5), Vector3(-1.5,2,-1.5), 
		Vector3(1.5,2.,-1.5), white);
	w->add_object(ceil2); 

	// Pared fondo
	Triangle* back1 = new Triangle(Vector3(1.5,2.5,-1), Vector3(-1.5,2.5,-1), 
		Vector3(-1.5,-.5,-1), white);
	w->add_object(back1); 
	Triangle* back2 = new Triangle(Vector3(1.5,2.5,-1), Vector3(-1.5,-.5,-1), 
		Vector3(1.5,-.5,-1), white);
	w->add_object(back2); 

	// Pared izquierda
	Triangle* left1 = new Triangle(Vector3(-1,2.5,-1.5),Vector3(-1,2.5,1.5), 
		Vector3(-1,-0.5,-1.5), red);
	w->add_object(left1); 
	Triangle* left2 = new Triangle(Vector3(-1,0,-1.5), Vector3(-1,2.5,1.5),
		Vector3(-1,-.5,1.5), red);
	w->add_object(left2); 

	// Pared derecha
	Object3D* right1 = new Triangle(Vector3(1,2.5,1.5), Vector3(1,2.5,-1.5), 
		Vector3(1,-.5,-1.5), green);
	w->add_object(right1); 
	Object3D* right2 = new Triangle(Vector3(1,2.5,1.5), Vector3(1,-.5,-1.5), 
		Vector3(1,-.5,1.5), green);
	w->add_object(right2);
	
	switch(scene)
	{
	case 1:
	{	
		Object3D* sphere1 = new Sphere(Vector3(0.5,1,0.2), 0.3, glass);
		w->add_object(sphere1);

		Object3D* sphere2 = new Sphere(Vector3(-0.5,0.5,0.3), 0.3, glass);
		w->add_object(sphere2);

		Object3D* sphere3 = new Sphere(Vector3(0,0.5,-0.5), 0.3, mirror);
		w->add_object(sphere3);
	}
	break;
	case 2:
	{	
		Object3D* sphere1 = new Sphere(Vector3(0.5,0.3,.5), 0.3, white);
		w->add_object(sphere1);

		Object3D* sphere2 = new Sphere(Vector3(-0.5,0.5,1.5), 0.3, red);
		w->add_object(sphere2);

		Object3D* sphere3 = new Sphere(Vector3(0.,0.3,.0), 0.3, white);
		w->add_object(sphere3);
	}
	break;
	case 3:
	{	
		Object3D* sphere1 = new Sphere(Vector3(0.5,0.3,.5), 0.3, glass);
		w->add_object(sphere1);


		Mesh* bunny = new Mesh("C:/Users/Jaime/git/PhotonMapping/SmallPM/data/bunny.obj", glass);
		w->add_object(bunny);
	}
	break;
	default:
	{
		Object3D* sphere1 = new Sphere(Vector3(0.5,0.3,.5), 0.3, white);
		w->add_object(sphere1);

		Object3D* sphere2 = new Sphere(Vector3(-0.5,0.5,1.5), 0.3, red);
		w->add_object(sphere2);
	}
	}
	LightSource* ls = new PointLightSource(w, Vector3(0,1.9,0), Vector3(5,5,5));
	w->add_light(ls);

	w->fix();
	// ----------------------------------------------------------------------
	// Create Film and rendering engine
	//
	film = new Film(sizex,sizey);
	pm = new PhotonMapping(w, photons_global, photons_caustic, photons_volumetric, photons_causticvol, max_shots, nb_nearest_photons,
		sigmaScattering, sigmaAbsorcion, lambda, participative );
	engine = new RenderEngine(w, film, &camera, pm);
	
	engine->render(name_file);	

}
