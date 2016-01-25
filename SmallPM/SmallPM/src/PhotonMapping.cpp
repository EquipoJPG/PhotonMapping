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
			   bool participative, bool direct,
			   double sT,
			   double sS,
			   double sA,
			   double landa)
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

		bool hasHit = it.did_hit();

		////////////////// CODIGO MEDIO PARTICIPATIVO ///////////////////
		// Se considera que toda la sala es medio participativo
		if(participative)
		{

			// Coeficientes (caracteristicas medio participativo)
			// Vienen dados por variables globales que se establecen en preprocess
			// Aunque a trace_ray se le pasan por parametro
			double sigmaT = sT;	// Coeficiente de extincion
			double sigmaS = sS;	// Coeficiente de scattering
			double sigmaA = sA;	// Coeficiente de absorcion
			double e = 2.7189;
			bool absorbido = false;
			double lambda = landa;			// Lambda (mean-free path = 1 / sigmaT)

			// Para cada paso
			Vector3 x(photon_ray.get_origin());			// Inicio
			Vector3 xs(it.get_position());				// Interseccion
			Vector3 w(photon_ray.get_direction());		// Direccion del rayo

			// Variables que indican el nuevo punto de cada paso y la direccion de cada paso
			Vector3 xp(x + w*lambda);				// Nuevo punto a comprobar
			Vector3 dirComp(xs - xp);				// Direccion a comprobar

			// Estas variables se usan para el caso en el que el rayo no interseccione con nada de la escena
			// De tal forma que iterara un maximo de 'max_pasitos' veces (Nota esto es una solucion propia que no
			// figura en ningun paper
			int pasitos = 0;
			int max_pasitos = 30;

			// Mientras no se haya pasado del punto de interseccion y no se haya absorbido y los pasitos sean menores
			// que max_pasitos (esto ultimo solo si no ha habido interseccion)
			while(dirComp.dot(w) >= 0 && !absorbido &&
				pasitos < max_pasitos)
			{
				// Ruleta rusa para saber si el foton continua avanzando
				double ruletitaRusa = fRand(0,1);
				// Calculo de la transmitancia en este paso
				Real Tt = 1.0 - pow(e, (-1) * (Vector3(xp - x).length() * sigmaT));

				if(ruletitaRusa <= Tt)
				{
					// Otra generacion aleatoria para ver si es scattering o absorcion
					ruletitaRusa = fRand(0,1);

					// El foton se ve alterado, se ha producido un evento
					if(ruletitaRusa <= sigmaS/sigmaT)
					{
						// Evento de scattering

						// Se guarda el foton si queda espacio en el mapa de volumen
						if( volumetric_photons.size() < m_nb_volumetric_photons )
							volumetric_photons.push_back( Photon(xp, dirComp, energy));

						// Calcula nueva direccion aleatoria
						// Funcion de fase es constante (1/4PI),
						// asi que la direccion es aleatoria
						double xd,yd,zd;
						xd = fRand(-1,1);
						yd = fRand(-1,1);
						zd = fRand(-1,1);
						Vector3 photonDir(xd,yd,zd);

						// Como nos hemos desviado necesitamos una nueva interseccion
						// con la que comprobar si nos hemos pasado de la escena
						Intersection temp;
						Ray* tempRay = new Ray(xp, photonDir);								
						world->first_intersection(*tempRay, temp);
						
						if( !temp.did_hit() )
						{
							// El nuevo rayo va al infinito
							// El medio participativo ocupa toda la escena,
							// por lo tanto si no hay ningun objeto delante
							// se siguen dando algunos pasos mas
							pasitos++;
						}

						// Si ya hemos dado muchos pasos
						if (pasitos >= max_pasitos) {
							absorbido = true;
						}
						else
						{
							// Actualiza los valores para el nuevo paso
							w = photonDir;
							x = Vector3(xp);
							xs = temp.get_position();
							xp = Vector3(x + w*lambda);
							dirComp = Vector3(xs - xp);
						}
					}
					else
					{
						// Evento de absorcion
						absorbido = true;

						// El foton se guarda si el mapa de volumen no esta lleno
						if( volumetric_photons.size() < m_nb_volumetric_photons )
							volumetric_photons.push_back( Photon(xp, dirComp, energy));
					}
				}
				else
				{
					// El foton sigue avanzando, recalculamos
					x = Vector3(xp);
					xp = Vector3(x + w*lambda);
					dirComp = Vector3(xs - xp);
				}
			}
		}
		////////////////// FIN MEDIO PARTICIPATIVO ///////////////////

		if( !it.did_hit() )
			break;

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
		global_photons.size() == m_nb_global_photons &&
		volumetric_photons.size() == m_nb_volumetric_photons)
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

	//////////////////// AJUSTE DE VARIABLES PARA SCATTERING ////////////////////
	double sigmaS = 0.3;
	double sigmaA = 0.1;
	double sigmaT = sigmaS + sigmaA;
	double lambda = 0.04;
	// Las globales las usa la funcion de shade para no pasarlas por parametro
	globalST = sigmaT;
	globalSS = sigmaS;
	globalSA = sigmaA;
	globalLambda = lambda;
	bool participativeRoom = true;
	////////////////////// FIN AJUSTE VARIABLES SCATTERING //////////////////////

	// Muestrea las fuentes de luz de la escena
	for(int i = 0; i < world->nb_lights(); i++){
		
		// Obtiene la fuente de luz i-esima
		Vector3 lightPos = world->light(i).get_position();
		Vector3 lightIntensity = world->light(i).get_intensities();
		LightSource* lt = new PointLightSource(world, lightPos, lightIntensity);
		Vector3 photonFlux(lightIntensity);	// energia foton = lightIntensity

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

			// Crea el rayo (foton) a lanzar
			Ray* photonRay = new Ray(lightPos, photonDir);

			// Lanza los fotones muestreados
			std::list<Photon> globalPhotons;
			std::list<Photon> causticPhotons;
			std::list<Photon> volumetricPhotons;

			trace_ray(*photonRay, photonFlux, globalPhotons, causticPhotons, volumetricPhotons, participativeRoom, false, sigmaT,
				sigmaS, sigmaA, lambda);

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

	cout << m_nb_current_shots << "/" << m_max_nb_shots << endl;
	cout << "GP: " << gp << endl;
	cout << "CP: " << cp << endl;
	cout << "VP: " << vp << endl;

	// FOTONES ALMACENADOS - PREPROCESO COMPLETADO
	if(gp > 0){
		m_global_map.balance();
	}

	if(cp > 0){
		m_caustics_map.balance();
	}

	if(vp > 0){
		m_volumetric_map.balance();
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
	// debug = 4 LUZ DIRECTA + INDIRECTA + SCATTERING
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
	Vector3 LDirecta = Vector3(0);
	Vector3 LDifusa = Vector3(0);
	Vector3 LCaustica = Vector3(0);
	Vector3 LScattering = Vector3(0);

	////////////////////////////// LUZ DIRECTA /////////////////////////////////////////
	if(debug == 1 || debug == 3 || debug == 4){
		// TERMINO AMBIENTAL
		LDirecta += world->get_ambient() * it.intersected()->material()->get_albedo(it);
	
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
				LDirecta += Kd * Id * cos;

				// TERMINO ESPECULAR: no hay porque son superficies lambertianas
			}
		}
	}
	/////////////////////////////////////////////////////////////////////////////

	///////////////////////// LUZ INDIRECTA /////////////////////////////////////
	if(debug == 2 || debug == 3 || debug == 4){
		/// LUZ INDIRECTA ///
		// pI = punto de interseccion (x,y,z)
		// pN = normal en el punto de interseccion

		std::vector<Real> intersection = std::vector<Real>();
		intersection.push_back(pI.getComponent(0));
		intersection.push_back(pI.getComponent(1));
		intersection.push_back(pI.getComponent(2));

		/// FOTONES DIFUSOS ///
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
				sumatorio += (photon.flux/(double) m_global_map.nb_elements()) * 
									it.intersected()->material()->get_albedo(it) * 
									((max_distance - distance) / max_distance);
			}
			// Calcula el area de un circulo de radio la distancia del foton
			// mas lejano (de los cercanos) respecto al punto de interseccion
			area = 3.14 * std::pow(max_distance, 2);
			LDifusa += sumatorio / area;
		}
		/////////////////////////////////////////////////////////////////////////

		///////////////////////// COMPONENTE CAUSTICA //////////////////////////////////////
		if(!m_caustics_map.is_empty()){
			m_caustics_map.find(intersection, m_nb_photons, nearest_photons, max_distance);
			sumatorio = Vector3(0);

			for (i = 0; i < nearest_photons.size(); i++) {

				// Obtiene la informacion de un foton
				Photon photon = nearest_photons.at(i)->data();
				Vector3 position = photon.position;
				Vector3 inters = it.get_position();

				// Calcular vector entre posicion y interseccion
				Vector3 vectorcito = position - inters;
				Real distance = vectorcito.length();

				/// ECUACION DE RENDER (suma de flujos de fotones) ///
				sumatorio += (photon.flux/(double) m_caustics_map.nb_elements()) * 
									it.intersected()->material()->get_albedo(it) * 
									((max_distance - distance) / max_distance);
			}
			// Calcula el area de un circulo de radio la distancia del foton
			// mas lejano (de los cercanos) respecto al punto de interseccion
			area = 3.14 * std::pow(max_distance, 2);
			LCaustica += sumatorio / area;
		}
		////////////////////////////////////////////////////////////////////////
	}

	/////////////////////// FOTONES VOLUMETRICOS //////////////////////////////
	if(debug == 4 || debug == 5)
	{
		
		if(!m_volumetric_map.is_empty()){

			double e = 2.7189;

			// Coeficientes (caracteristicas medio participativo)
			double sigmaT = globalST;		// Coeficiente de extincion
			double sigmaS = globalSS;		// Coeficiente de scattering
			double sigmaA = globalSA;		// Coeficiente de absorcion
			double landa = 0.5;				// Asignado arbitrariamente. Puede ser 1 / sigmaT? (mean-free path?)

			Intersection itV(it0);

			// Marching ray
			Ray marchingRay(itV.get_ray());
			Vector3 xOriginal(marchingRay.get_origin());
			Vector3 x = Vector3(xOriginal);
			Vector3 w(marchingRay.get_direction());
			
			// Punto destino (xs)
			Vector3 xs(itV.get_position());
			
			// Obtiene el nuevo punto (xt)
			Vector3 xt = Vector3(x + w*landa);
			Vector3 dirComp = Vector3(xs - xt);
			
			// Funcion de fase isotropica (medio participativo homogeneo => 1/(4*PI) )
			double phaseFunction = 1.0 / (4.0*3.14159);

			// ECUACION VOLUMETRICA DE RENDER //
			// e^ -(||xs - x|| * sigmaT)
			Real Ts = pow(e, (-1) * (Vector3(xs - xOriginal).length() * sigmaT));

			// Calcula ray marching
			Vector3 sumLi = Vector3(0);
			Vector3 sumTt = Vector3(0);
			Vector3 sumInScattering = Vector3(0);
			Vector3 LK = Vector3(0);
			Vector3 LKPrev = Vector3(0);
			while(dirComp.dot(w) >= 0 && itV.did_hit())
			{

				// NUEVO PASO EN RAY MARCHING
				std::vector<Real> xp = std::vector<Real>();
				xp.push_back(xt.getComponent(0));
				xp.push_back(xt.getComponent(1));
				xp.push_back(xt.getComponent(2));

				// Transmitancia hasta punto xt
				// e^ -(||xt - x|| * sigmaT)
				Real Tt = pow(e, (-1) * (Vector3(xt - xOriginal).length() * sigmaT));
				
				Vector3 Lid = Vector3(0);
				// Luz directa. Este bucle debe estar mal hecho porque provoca una niebla muy densa
				// que solo se aprecia reduciendo la luz de la imagen en Photoshop
				for(int i = 0; i < world->nb_lights(); i++){
					Vector3 lightPos = world->light(i).get_position();
					double Tr = 1.0 - pow(e, (-1) * (Vector3(lightPos - x).length() * sigmaT));
					Vector3 intensity = world->light(i).get_intensities();
					Lid += Tr * intensity;
				}

				// Obtiene los k fotones cercanos al paso t (xp = xt)
				float radio = 0;
				std::vector<const KDTree<Photon, 3>::Node*> nearest_photons;
				m_volumetric_map.find(xp, m_nb_photons * 2, nearest_photons, radio);
				double volumen = (4.0 / 3.0) * 3.14159 * std::pow(radio,3);

				// Calcula Lii (luz proveniente del evento in-scattering)
				Vector3 Lii = Vector3(0);

				for (int i = 0; i < nearest_photons.size(); i++) {
					Photon photon = nearest_photons.at(i)->data();
					Lii += phaseFunction * ((photon.flux / m_volumetric_map.nb_elements()) / volumen);
				}
				
				Lii = Lii * (1.0/sigmaS);
				Lii = Lii * (sigmaS/sigmaT);

				// Calculo Li
				Vector3 Li = Vector3(0);
				Li = Lid + Lii;

				// Luz paso k esimo
				Vector3 deltaXk = Vector3(xt - x);
				float epsilon = std::pow(e,(-1.0) * sigmaT * deltaXk.length());

				LK = (deltaXk.length() * sigmaS * Li) + LKPrev * epsilon;
				LKPrev = Vector3(LK);
				
				// Calcula la radiancia ganada por in-scattering
				//sumInScattering += Tt * sigmaS * Li;
				
				//sumLi += (1 - Tt) / sigma;
				//sumTt += (1 - Tt) / sigmaT;

				// Obtiene el nuevo punto (xt)
				x = Vector3(xt);
				xt = Vector3(x + w*landa);
				dirComp = Vector3(xs - xt);
			}

			/// ECUACION VOLUMETRICA DE RENDER FINAL ///
			LScattering += LK;
		}
		////////////////////////////////////////////////////////////////////////
	}

	double e = 2.7189;
	Vector3 origen = it.get_ray().get_origin();
	Vector3 destino = it.get_position();
	Real transmit = pow(e, (-1) * (Vector3(destino - origen).length() * globalST));
	L += ((LDirecta) *transmit) + LDifusa + LCaustica + LScattering;

	return L;
}