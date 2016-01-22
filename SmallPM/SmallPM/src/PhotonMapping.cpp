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
#include "PhotonMapping.h"
#include "World.h"
#include "Intersection.h"
#include "Ray.h"
#include "BSDF.h"
#include <math.h>

double fRand(double fMin, double fMax)
{
	double f = (double) rand() / RAND_MAX;
	return fMin + f* (fMax - fMin);
}

//*********************************************************************
// Compute the photons by tracing the Ray 'r' from the light source
// through the scene, and by storing the intersections with matter
// in the lists 'xx_photons', storing the diffuse (global) and caustic
// photons respectively. For efficiency, both are computed at the same
// time, since computing them separately would result into a lost of
// several samples marked as caustic or diffuse.
// Same goes with the boolean 'direct', that specifies if direct 
// photons (from light to surface) are being stored or not. 
// The initial traced photon has energy defined by the tristimulus
// 'p', that accounts for the emitted power of the light source.
// The function will return true when there are more photons (caustic
// or diffuse) to be shot, and false otherwise.
//---------------------------------------------------------------------
bool PhotonMapping::trace_ray(const Ray& r, const Vector3 &p, 
			   std::list<Photon> &global_photons, 
			   std::list<Photon> &caustic_photons, 
			   std::list<Photon> &volumetric_photons, 
			   bool participative, bool direct)
{

	//Check if max number of shots done...
	if( ++m_nb_current_shots > m_max_nb_shots )
	{
		return false;
	}
	
	// Compute irradiance photon's energy
	Vector3 energy(p);
	
	Ray photon_ray(r);
	photon_ray.shift();

	bool is_caustic_particle = false;

	//Iterate the path
	while(1)
	{
		// Throw ray and update current_it
		Intersection it;
		world->first_intersection(photon_ray, it);

		if( !it.did_hit() )
			break;

		////////////////// MEDIO PARTICIPATIVO CODE ///////////////////
		if(participative)
		{
			// Coeficiente de extincion
			double sigmaT = 0.4;

			// Coeficiente de scattering
			double sigmaS = fRand(0,sigmaT);

			// Coeficiente de absorcion
			double sigmaA = sigmaT - sigmaS;

			// Booleano para saber si se absorbio
			bool absorbido = false;

			// Lambda
			double lambda = 1.0/sigmaT;

			// Inicio
			Vector3 x(photon_ray.get_origin());

			// Interseccion
			Vector3 xs(it.get_position());

			// Direccion del rayo
			Vector3 dir(photon_ray.get_direction());

			// Nuevo punto a comprobar
			Vector3 xp(x + dir*lambda);

			// Direccion de comprobar
			Vector3 dirComp(xs - xp);

			// Mientras no se haya pasado del punto de interseccion
			while(dirComp.dot(dir)>=0 && !absorbido)
			{
				// Ruleta rusa para saber si el foton continua avanzando
				double ruletitaRusa = fRand(0,1);

				if(ruletitaRusa <= sigmaT)
				{
					// El foton se ve alterado, se ha producido un evento
					if(ruletitaRusa >= sigmaS)
					{
						// Evento de scattering
						// Se guarda
						if( volumetric_photons.size() < m_nb_volumetric_photons )
							volumetric_photons.push_back( Photon(xp, dirComp, energy));
						// Calculamos nueva direccion aleatoria
						double xd,yd,zd;
						xd = fRand(-1,1);
						yd = fRand(-1,1);
						zd = fRand(-1,1);
						Vector3 photonDir(xd,yd,zd);
						Ray* tempRay = new Ray(xp, photonDir);

						Intersection temp;
						world->first_intersection(*tempRay, temp);

						// Si el nuevo rayo se va al infinito lo almacenamos y ya
						if( !temp.did_hit() )
						{
							absorbido = true;
						}
						else
						{
							// Actualizamos valores
							dir = photonDir;
							x = Vector3(xp);
							xs = temp.get_position();
							xp = Vector3(x + dir*lambda);
							dirComp = Vector3(xs - xp);
						}

					}
					else
					{
						// Evento de absorcion
						absorbido = true;
						// Si ha sido absorbido, se guarda
						if( volumetric_photons.size() < m_nb_volumetric_photons )
							volumetric_photons.push_back( Photon(xp, dirComp, energy));
					}
				}
				else
				{
					// El foton sigue avanzando, recalculamos
					x = Vector3(xp);
					xp = Vector3(x + dir*lambda);
					dirComp = Vector3(xs - xp);
				}
			}
		}
		////////////////// FIN MEDIO PARTICIPATIVO ///////////////////

		//Check if has hit a delta material...
		if( it.intersected()->material()->is_delta() )
		{
			// If delta material, then is caustic...
			// Don't store the photon!
			is_caustic_particle = true;
		}
		else if (photon_ray.get_level() > 0 || direct)
		{
			//If non-delta material, store the photon!
			if( is_caustic_particle )	
			{				
				//If caustic particle, store in caustics
				if( caustic_photons.size() < m_nb_caustic_photons )
					caustic_photons.push_back( Photon(it.get_position(), photon_ray.get_direction(), energy ));
			}
			else						
			{
				//If non-caustic particle, store in global
				if( global_photons.size() < m_nb_global_photons )
					global_photons.push_back( Photon(it.get_position(), photon_ray.get_direction(), energy ));
			}
			is_caustic_particle = false;
		}	
		
		// INICIO RULETA RUSA /////////////
		Real pdf;

		Vector3 surf_albedo = it.intersected()->material()->get_albedo(it);
		Real avg_surf_albedo = surf_albedo.avg();

		Real epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
		while (epsilon2 < 0.)
			epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
		
		if (epsilon2 > avg_surf_albedo || photon_ray.get_level() > 20 ) 
			break;
			
		// Random walk's next step
		// Get sampled direction plus pdf, and update attenuation
		it.intersected()->material()->get_outgoing_sample_ray(it, photon_ray, pdf );

		// FIN RULETA RUSA ///////////////

		// Shade...
		energy = energy*surf_albedo;
		if( !it.intersected()->material()->is_delta() )
			energy *= dot_abs(it.get_normal(), photon_ray.get_direction())/3.14159;		

		energy = energy /(pdf*avg_surf_albedo);
	}
	
	if( caustic_photons.size() == m_nb_caustic_photons && 
		global_photons.size() == m_nb_global_photons )
	{
		m_max_nb_shots = m_nb_current_shots-1;
		return false;
	}

	return true;
}

//*********************************************************************
// TODO: Implement the preprocess step of photon mapping,
// where the photons are traced through the scene. To do it,
// you need to follow these steps for each shoot:
//  1 - Sample a world's light source in the scene to create
//		the initial direct photon from the light source.
//	2 - Trace the photon through the scene storing the inter-
//		sections between the photons and matter. You can use
//		the function 'trace_ray' for this purpose.
//	3 - Finally, once all the photons have been shot, you'll
//		need to build the photon maps, that will be used later
//		for rendering. 
//		NOTE: Careful with function
//---------------------------------------------------------------------
void PhotonMapping::preprocess()
{
	int gp = 0;
	int cp = 0;
	int vp = 0;
	bool participativeRoom = true;
	// Muestrea las fuentes de luz de la escena
	for(int i = 0; i < world->nb_lights(); i++){
		
		// Obtiene la fuente de luz i-esima
		Vector3 lightPos = world->light(i).get_position();
		Vector3 lightIntensity = world->light(i).get_intensities();
		LightSource* lt = new PointLightSource(world, lightPos, lightIntensity);
		//Vector3 photonFlux(lightIntensity / lightIntensity);	// energia foton = lightIntensity
		Vector3 photonFlux(lightIntensity / m_max_nb_shots);	// energia foton = lightIntensity / total fotones

		// Muestreo de una esfera, se lanza un rayo en una direccion aleatoria
		// de la esfera. El numero de fotones lanzados es el maximo definido por
		// la variable 'm_max_nb_shots'
		while (m_nb_current_shots < m_max_nb_shots)
		{
			
			// REJECTION SAMPLING
			double x,y,z;
			do {
				x = fRand(-1,1);
				y = fRand(-1,1);
				z = fRand(-1,1);
			} while (pow(x,2) + pow(y,2) + pow(z,2) > 1);
			
			Vector3 photonDir(x,y,z);
			
			//cout << "Omega: " << omega << ", Theta: " << theta << "\n";
			//cout << "(" << x << ", " << y << ", " << z << ")\n";
			//cout << "====================\n";

			// Crea el rayo (foton) a lanzar
			Ray* photonRay = new Ray(lightPos, photonDir);

			// Lanza los fotones muestreados
			std::list<Photon> globalPhotons;
			std::list<Photon> causticPhotons;
			std::list<Photon> volumetricPhotons;

			trace_ray(*photonRay, photonFlux, globalPhotons, causticPhotons, volumetricPhotons, participativeRoom, false);

			// Almacena las colisiones de los fotones difusos
			int k;
			for (k = 0; k < globalPhotons.size(); k++) {
				gp++;

				// Obtiene el foton, lo guarda en el KDTree y lo borra de la lista
				Photon photon = globalPhotons.front();

				std::vector<Real> photonPosition = std::vector<Real>();
				photonPosition.push_back(photon.position.getComponent(0));
				photonPosition.push_back(photon.position.getComponent(1));
				photonPosition.push_back(photon.position.getComponent(2));
				
				m_global_map.store(photonPosition, photon);

				globalPhotons.pop_front(); // elimina el foton almacenado de la lista
			}

			// Almacena las colisiones de los fotones causticos
			for (k = 0; k < causticPhotons.size(); k++) {
				cp++;

				// Obtiene el foton, lo guarda en el KDTree y lo borra de la lista
				Photon photon = causticPhotons.front();

				std::vector<Real> photonPosition = std::vector<Real>();
				photonPosition.push_back(photon.position.getComponent(0));
				photonPosition.push_back(photon.position.getComponent(1));
				photonPosition.push_back(photon.position.getComponent(2));
				
				m_caustics_map.store(photonPosition, photon);

				causticPhotons.pop_front(); // elimina el foton almacenado de la lista
			}

			// Almacena las colisiones de los fotones volumetricos
			for (k = 0; k < volumetricPhotons.size(); k++) {
				vp++;

				// Obtiene el foton, lo guarda en el KDTree y lo borra de la lista
				Photon photon = volumetricPhotons.front();

				std::vector<Real> photonPosition = std::vector<Real>();
				photonPosition.push_back(photon.position.getComponent(0));
				photonPosition.push_back(photon.position.getComponent(1));
				photonPosition.push_back(photon.position.getComponent(2));
				
				m_volumetric_map.store(photonPosition, photon);

				volumetricPhotons.pop_front(); // elimina el foton almacenado de la lista
			}
		}
	}
	// FOTONES ALMACENADOS - PREPROCESO COMPLETADO
	if(gp > 0){
		m_global_map.balance();
	}
	//cout << m_global_map.nb_elements() << "\n";
	if(cp > 0){
		m_caustics_map.balance();
	}

	if(vp > 0){
		m_volumetric_map.balance();
		cout << m_volumetric_map.nb_elements() << endl;
	}
}

//*********************************************************************
// TODO: Implement the function that computes the rendering equation 
// using radiance estimation with photon mapping, using the photon
// maps computed as a proprocess. Note that you will need to handle
// both direct and global illumination, together with recursive the 
// recursive evaluation of delta materials. For an optimal implemen-
// tation you should be able to do it iteratively.
// In principle, the class is prepared to perform radiance estimation
// using k-nearest neighbors ('m_nb_photons') to define the bandwidth
// of the kernel.
//---------------------------------------------------------------------
Vector3 PhotonMapping::shade(Intersection &it0)const
{
	// TODO
	// debug = 1 LUZ DIRECTA
	// debug = 2 LUZ INDIRECTA
	// debug = 3 LUZ DIRECTA + INDIRECTA
	int debug = 4;

	// ESTRUCTURA
	// -----------------------------------------------------------
	// 1.- Calcular iluminacion directa en el punto (a cada PL)
	// 2.- Calcular iluminacion indirecta en el punto a traves
	//		de la estimacion de radiancia
	// -----------------------------------------------------------
	// Debugueo: poner por separado ID e II

	Vector3 L(world->get_background());
	Intersection it(it0);
	Vector3 pN = it.get_normal(); // normal en el punto de interseccion
	Vector3 pI = it.get_position();	// punto de interseccion (x,y,z)

	// REBOTAR MIENTRAS EL OBJETO SEA DELTA (hay que llegar a un solido)
	int MAX_REB = 3;
	int rebotes = 0;
	Ray newRay;

	while (it.intersected()->material()->is_delta() && rebotes < MAX_REB) {

		// Rayo rebotado
		Real pdf;
		it.intersected()->material()->get_outgoing_sample_ray(it, newRay, pdf );

		// Nueva interseccion
		newRay.shift();
		world->first_intersection(newRay, it);
		rebotes++;
	}

	pI = it.get_position();	// punto de interseccion (x,y,z)
	pN = it.get_normal();
	//////////////// FIN DE REBOTES //////////////////
	
	if(debug == 1 || debug == 3 || debug == 4){
		// TERMINO AMBIENTAL
		L += world->get_ambient() * it.intersected()->material()->get_albedo(it);
	
		// LUZ DIRECTA //
		for(int i = 0; i < world->nb_lights(); i++){
		
			// Obtiene la fuente de luz i-esima
			Vector3 lightPos = world->light(i).get_position();
			Vector3 lightIntensity = world->light(i).get_intensities();
			LightSource* lt = new PointLightSource(world, lightPos, lightIntensity);

			Vector3 shadowRay = Vector3() - lt->get_incoming_direction(pI); // (0,0,0) - lightRay = shadowRay

			// Si el objeto es visible se calcula la influencia de la luz
			if (lt->is_visible(pI)) {
				// TERMINO DIFUSO = Kd x Id x (L . N)
				Vector3 Id = lt->get_incoming_light(pI);
				Vector3 Kd = it.intersected()->material()->get_albedo(it);
				float cos = shadowRay.dot(pN);
				L += Kd * Id * cos;

				// TERMINO ESPECULAR: no hay porque son superficies lambertianas
			}
		}
	}
	if(debug == 2 || debug == 3 || debug == 4){
		// LUZ INDIRECTA (Estimacion de radiancia) //
		// pI = punto de interseccion (x,y,z)
		// pN = normal en el punto de interseccion

		std::vector<Real> intersection = std::vector<Real>();
		intersection.push_back(pI.getComponent(0));
		intersection.push_back(pI.getComponent(1));
		intersection.push_back(pI.getComponent(2));

		// FOTONES DIFUSOS
		// Busca los fotones cercanos y los guarda en nearest_photons
		std::vector<const KDTree<Photon, 3>::Node*> nearest_photons;
		Real max_distance;
		Vector3 sumatorio = Vector3(0);
		int i;
		Real area;
	
		if(!m_global_map.is_empty()){
			m_global_map.find(intersection, m_nb_photons, nearest_photons, max_distance);

			for (i = 0; i < nearest_photons.size(); i++) {

				// Obtiene la informacion de un foton
				Photon photon = nearest_photons.at(i)->data();
				Vector3 position = photon.position;
				Vector3 inters = it.get_position();

				// Calcular vector entre posicion y interseccion
				Vector3 vectorcito = position - inters;
				Real distance = vectorcito.length();

				/// ECUACION DE RENDER (suma de flujos de fotones) ///
				sumatorio += photon.flux * it.intersected()->material()->get_albedo(it) * ((max_distance - distance) / max_distance);
			}
			// Calcula el area de un circulo de radio la distancia del foton
			// mas lejano (de los cercanos) respecto al punto de interseccion
			Real area = 3.14 * std::pow(max_distance, 2);
			L += sumatorio / area;
		}
		/////////////////////////////////////////////////////////////////////////
		// FOTONES CAUSTICOS
		if(!m_caustics_map.is_empty()){
			m_caustics_map.find(intersection, m_nb_photons, nearest_photons, max_distance);
			sumatorio = Vector3(0);
			for (i = 0; i < nearest_photons.size(); i++) {

				// Obtiene la informacion de un foton
				Photon photon = nearest_photons.at(i)->data();
				Vector3 position = photon.position;

				/// ECUACION DE RENDER (suma de flujos de fotones) ///
				sumatorio += photon.flux * it.intersected()->material()->get_albedo(it);
			}
			// Calcula el area de un circulo de radio la distancia del foton
			// mas lejano (de los cercanos) respecto al punto de interseccion
			area = 3.14 * std::pow(max_distance, 2);
			L += sumatorio / area;
		}
		////////////////////////////////////////////////////////////////////////

		//cout << "Area: " << area << "\n";
		//cout << "Nearest photons: " << nearest_photons.size() << "\n";

		//cout << "Luz final(rojo): " << L.getComponent(0) << "\n";
		//cout << "Luz final(verde): " << L.getComponent(1) << "\n";
		//cout << "Luz final(azul): " << L.getComponent(2) << "\n";
	}
	if(debug == 4)
	{
		/////////////////////////////////////////////////////////////////////////
		// FOTONES VOLUMETRICOS
		if(!m_volumetric_map.is_empty()){
			// Ray Marching, si, for real, confirmed 2016
			// Copiamos la interseccion por si acaso
			Intersection itV(it0);
			// Calculamos el marching ray
			Ray marchingRay(itV.get_ray());
			// Sacamos el punto de origen
			Vector3 origen(marchingRay.get_origin());
			// Sacamos el punto de destino
			Vector3 destino(itV.get_position());
			// Sacamos el vector de direccion
			Vector3 dir(marchingRay.get_direction());
			// Sacamos el nuevo punto
			double sigmaT = 0.1;
			double landa = 1.0/sigmaT;
			Vector3 nuevo(origen + dir*landa);
			// Sacamos el vector entre destino y nuevo
			Vector3 dirComp(destino - nuevo);
			// Establecemos el radio para beam estimate
			float radius = 3.0;
			// Creamos nearest photons
			std::vector<const KDTree<Photon, 3>::Node*> nearest_photons;
			// Funcion de fase
			double phaseFunction = 1.0 / (4.0*3.14159);

			while(dirComp.dot(dir) >= 0)
			{
				Vector3 sumatorio = Vector3(0);
				std::vector<Real> nuevoVR = std::vector<Real>();
				nuevoVR.push_back(pI.getComponent(0));
				nuevoVR.push_back(pI.getComponent(1));
				nuevoVR.push_back(pI.getComponent(2));
				m_volumetric_map.find(nuevoVR, m_nb_photons, nearest_photons, radius);
				int i;
				for (i = 0; i < nearest_photons.size(); i++) {

					// Obtiene la informacion de un foton
					Photon photon = nearest_photons.at(i)->data();
					Vector3 position = photon.position;

					/// ECUACION DE RENDER (suma de flujos de fotones) ///
					sumatorio += photon.flux * phaseFunction;
				}
				// Calcula el area de un circulo de radio la distancia del foton
				// mas lejano (de los cercanos) respecto al punto de interseccion
				double volumen = (4.0 / 3.0) * 3.14159 * std::pow(radius,3);
				L += sumatorio / volumen;
			}
		}
		////////////////////////////////////////////////////////////////////////
	}


	return L;
}