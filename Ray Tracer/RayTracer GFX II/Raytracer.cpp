#include <math.h>
/***************************************************************************
*                                                                          *
* This is the source file for a ray tracer. It defines most of the		   *
* fundamental functions using in ray tracing.  In particular, "Trace" and  *
* "Shade".  It also defines the MakeImage function, which casts all the    *
* rays from the eye, and several object intersection methods.  Some of the *
* functions are left blank, or with "place holders" that don't do very     *
* much.  You need to fill these in.                                        *
*                                                                          *
*                                                                          *
***************************************************************************/

static const int tree_depth = 2;		// Number of recursions to compute indirect illumination
static const int rays_pixel = 128;		// Quality on deferred rendering, Power of 2, less value = more noise less render time(128); more value = less noise = huge amount of render time (8192)

#include "Raytracer.h"

// Draw image on the screen
void Raytracer::draw( void )
{
	glDrawPixels( resolutionX, resolutionY, GL_RGB, GL_UNSIGNED_BYTE, &(*I)( 0 , 0 ) );
}

// Cast_line casts all the initial rays starting from the eye for a single
//  raster line. Copies pixels to image object.
void Raytracer::cast_line( World world )
{
    Ray ray;
	Ray fRay; // Raig percial pel depth of field
	Color color = Color(0,0,0);	// Color computed when using multiple rays per pixel
	float k;

	ray.origin = world.getCamera().eye; // All initial rays originate from the eye.
	ray.no_emitters = false;

    Vec3 G  = Unit( world.getCamera().lookat - world.getCamera().eye );	// Gaze direction.
    Vec3 U  = Unit( world.getCamera().up / G );							// Up vector.
    Vec3 R  = Unit( G ^ U );											// Right vector.
    Vec3 O  = ( world.getCamera().vpdist * G ) - R + U;					// "Origin" for the raster.
    Vec3 dU = U * ( 2.0 / ( resolutionY - 1 ) );						// Up increments.
	Vec3 dR = R * ( 2.0 / ( resolutionX - 1 ) );						// Right increments.
	float r = world.getCamera().lens_radius;
	float d = world.getCamera().focaldist;

	if( currentLine % 10 == 0 ) cout << "line " << currentLine << endl;
    for( int i = 0; i < resolutionX; i++ )
    {
		if( rays_pixel == 1 )
		{
			// One ray per pixel
			ray.direction = Unit( O + i * dR - currentLine * dU  );
			color = Trace( ray, world.getScene(), tree_depth );
		}
		else if( world.getCamera().lens_radius < Epsilon )
		{
			// Multisampling
			for( int n = 0 ; n < rays_pixel ; n++ )
			{
				ray.direction = Unit( O + ( i + rand( 0.0 , 1.0 ) - 0.5 ) * dR - ( currentLine + rand( 0.0 , 1.0 ) - 0.5 ) * dU  );
				color += Trace( ray, world.getScene(), tree_depth );
			}
		}
		else
		{
			//Depth of Field
			for( int n = 0 ; n < rays_pixel ; n++ )
			{
				ray.direction = Unit( O + ( i + rand( 0.0 , 1.0 ) - 0.5 ) * dR - ( currentLine + rand( 0.0 , 1.0 ) - 0.5 ) * dU  );
				k = d/(ray.direction*G); //Calculamos la distancia en la que se encuentra el punto nítido en dirección en la que lanzamos el rayo
				Vec3 c = world.getCamera().eye; //Punto en el que se encuentra la cámara
				Vec3 F = c + (ray.direction*k); //F es el punto nítido, calculado con la dirección en que se traza el rayo, multiplicado por la distancia k y desplazado al centro de la cámara
				Sample s = SampleDisk(ray.direction); //Sacamos una muestra del espacio donde nos moveremos para crear el desenfoque de movimiento
				ray.origin = (s.P*r)+ world.getCamera().eye; // El origen tendrá en el punto espacial encontrado de la función anterior, desplazado al punto de cámara.
				ray.direction = Unit(F-ray.origin); //La dirección tomada será la del punto F al origen de la cámara.
				color += Trace( ray, world.getScene(), tree_depth ); //Finalmente recogemos el color.
			}
		}

		(*I)( resolutionY-currentLine-1, i ) = ToneMap( color / rays_pixel );

		color.blue = 0;
		color.green = 0;
		color.red = 0;

    }

	if (++currentLine == resolutionY)
	{
		// Image computation done, save it to file
		cout << "done." << endl;
	    I->Write( "Resultat.ppm" );
		isDone = true;
	}
}


// This is a trivial tone mapper; it merely maps values that are
// in [0,1] and maps them to integers between 0 and 255.  If the
// real value is above 1, it merely truncates.  A true tone mapper
// would attempt to handle very large values nicely, without
// truncation; that is, it would try to compensate for the fact that
// displays have a very limited dynamic range.
Pixel Raytracer::ToneMap( const Color &color )
{
	int red   = (int)floor( 256 * color.red   );
    int green = (int)floor( 256 * color.green );
    int blue  = (int)floor( 256 * color.blue  );
    channel r = (channel)( red   >= 255 ? 255 : red   ); 
    channel g = (channel)( green >= 255 ? 255 : green ); 
    channel b = (channel)( blue  >= 255 ? 255 : blue  );
    return Pixel( r, g, b );
}

// Trace is the most fundamental of all the ray tracing functions.  It
// answers the query "What color do I see looking along the given ray
// in the current scene?"  This is an inherently recursive process, as
// trace may again be called as a result of the ray hitting a reflecting
// object.  To prevent the possibility of infinite recursion, a maximum
// depth is placed on the resulting ray tree.
Color Raytracer::Trace( const Ray &ray, const Scene &scene, int max_tree_depth  )
{
    Color   color;                    // The color to return.
    HitInfo hitinfo;                  // Holds info to pass to shader.

	// Intitallizes hit distance to infinity to allow finding intersections in all ray length
	hitinfo.geom.distance = Infinity;

	if (Cast( ray, scene, hitinfo ) && max_tree_depth > 0 )
	{
        // The ray hits an object, so shade the point that the ray hit.
        // Cast has put all necessary information for Shade in "hitinfo".
		
		// If the ray has no_emitters activated and the first hit is an emitter
		//  this ray shouldn't contribute to the color of the current pixel
		if( hitinfo.material.Emitter() && ray.no_emitters == true ) 
			color = Color ();

		// The ray hits an object, so shade the point that the ray hit.
        // Cast has put all necessary information for Shade in "hitinfo".
		else 
			color = Shade( hitinfo, scene, max_tree_depth - 1  );
    }
    else
    {
        // Either the ray has failed to hit anything, or
        // the recursion has bottomed out.
        color = scene.bgcolor;
    }
    
    return color;
}

// Cast finds the first point of intersection (if there is one)
// between a ray and a list of geometric objects.  If no intersection
// exists, the function returns false.  Information about the
// closest object hit is returned in "hitinfo". 
bool Raytracer::Cast( const Ray &ray, const Scene &scene, HitInfo &hitinfo, Object *ignore )
{
    bool hit = false;

    // Each intersector is ONLY allowed to write into the "HitGeom"
    // structure if it has determined that the ray hits the object
    // at a CLOSER distance than currently recorded in HitGeom.distance.
    // When a closer hit is found, the material fields of the "HitInfo"
    // structure are updated to hold the material of the object that 
    // was just hit.

    for( Object *object = scene.first; object != NULL; object = object->next )
    {
        if( object != ignore && object->Intersect( ray, hitinfo.geom ) )
            {
            hitinfo.material = object->material;  // Material of closest surface.
            hit = true;                           // We have hit an object.
            }
    }
    return hit;
}


// Shade assigns a color to a point on a surface, as it is seen
// from another point.  The coordinates of these points, the normal
// of the surface, and the surface material are all recorded in the
// HitInfo structure.  The shader will typically make calls to Trace
// to handle shadows and reflections.
Color Raytracer::Shade( const HitInfo &hit, const Scene &scene, int max_tree_depth )
{
	
	if (hit.material.Emitter())
	{
		return hit.material.m_Emission;
	}

	//Definimos todos los colores que calcularemos (Difusa + Especular = Directa, Reflexiva Difusa y Especular = Indirecta, suma de ambas = Final)
	Color colorDifuso = Color(0,0,0);
	Color colorEspecular = Color (0,0,0);
	Color colorReflectionDifusa = Color(0,0,0);
	Color colorReflectionEspecular = Color(0,0,0);

	Color colorDirecto = Color(0,0,0);
	Color colorIndirecto = Color(0,0,0);
	Color colorFinal = Color(0,0,0);

	//Calculamos la direccion del rayo y se la asignamos a CameraRay
	Vec3 cameraRay = hit.geom.point - hit.geom.origin;
	cameraRay = cameraRay / Length(cameraRay);

	//Calculamos r como el vector reflejado sobre la superficie impactada.
	Vec3 r = Reflection(cameraRay,hit.geom.normal);
	r = r/Length(r);

	//Creamos un rayo auxiliar L para calcular la iluminación directa.
	//Este rayo será trazado desde el punto de colisión hasta la luz, para saber si el punto está en sombra o no.
	Ray	L;
	L.origin = hit.geom.point;
	L.no_emitters = false;

	HitInfo newHit;
	newHit = hit;
	newHit.geom.origin = hit.geom.point;

	//En estas sentencias buscaremos sobre todos los objetos de la escena los objetos que emitan luz.
	Object *o;
	for( Object *object = scene.first; object != NULL; object = object->next )
    {
		if (object->material.Emitter())
		{
			o = object;
		}
	}

	//Llamamos a la función GetSample para recojer una muestra hacia la luz.
	Sample sample = o->GetSample(hit.geom.point, Vec3());

	L.direction = (sample.P - hit.geom.point) / Length(sample.P - hit.geom.point);
	bool intersect = Cast(L , scene, newHit );
	
	//En caso que la muestra generada anteriormente intersecte con un objeto que no sea la luz, estaremos en punto de sombra, en otro caso, deberemos calcular la iluminación
	if(intersect && newHit.geom.distance > Length(sample.P - hit.geom.point) - Epsilon)
		intersect = false;

	if (!intersect)
	{
		//Color Difuso
		Color Difuso = (hit.geom.normal * L.direction) * o->material.m_Emission * hit.material.m_Diffuse * sample.w / Pi;

		//Nos aseguramos que ningún valor sea negativo.
		if(Difuso.blue < 0) Difuso.blue = 0;
		if(Difuso.red < 0) Difuso.red = 0;
		if(Difuso.green < 0) Difuso.green = 0;

		//Sumamos el total de color
		colorDifuso += Difuso;

		//Calculamos RL como el ángulo diferencia entre la reflexión y el vector que va a la luz.
		double RL = r * L.direction;
		
		//Especular color
		//En caso que el objeto tenga coeficiente reflectivo y un angulo de reflexión superior a cero, calculamos la iluminación especular directa.
		if(hit.material.m_Phong_exp != 0 && RL > 0){
			
			//Calculamos el color especular, cerciorándonos que no hayan valores negativos.
			colorEspecular += pow(RL, (double)hit.material.m_Phong_exp) * o->material.m_Emission * hit.material.m_Specular * sample.w * (hit.material.m_Phong_exp + 2)/(2*Pi);
			if(colorEspecular.blue < 0) colorEspecular.blue = 0;
			if(colorEspecular.red < 0) colorEspecular.red = 0;
			if(colorEspecular.green < 0) colorEspecular.green = 0;
		}

		//Finalmente acabamos con el color directo como la suma de la difusa con la especular.
		colorDirecto = colorDifuso + colorEspecular;

	}

	//Reflexion difusa
	//Cogemos una muestra de entre un radio semiesférico sobre la normal en la que el rayo ha impactado.
	Sample S = SampleProjectedHemisphere(hit.geom.normal);
	Ray PS; //Creamos un nuevo rayo del punto de origen al de la muestra.
	PS.origin = hit.geom.point;
	PS.direction = S.P;
	PS.no_emitters = true;
	//Finalmente calculamos la reflexión difusa con la fórmula.
	colorReflectionDifusa = hit.material.m_Diffuse * (S.w/Pi) * Trace(PS,scene,max_tree_depth);

	//Reflexion Especular
	//En caso que el objeto tenga coeficiente especular procedemos a calcular.
	if(hit.material.m_Phong_exp != 0){
		S = SampleSpecularLobe(r,hit.material.m_Phong_exp); //Cogemos una muestra Sobre un lóbulo de la normal en que se ha intersectado
		PS.direction = S.P; //De nuevo y reutilizando el rayo anterior ponemos el nuevo punto de muestra como director.
		//Finalmente aplicamos la formula para sacar el elemento reflectivo especular.
		colorReflectionEspecular = hit.material.m_Specular * (S.w * (hit.material.m_Phong_exp + 2)/(2*Pi)) * Trace(PS,scene,max_tree_depth);
	}

	//Sumamos los coeficientes de difusa y especular reflectivos para sonsacar la iluminación indirecta
	colorIndirecto = colorReflectionDifusa + colorReflectionEspecular;

	//Finalmente sumamos la iluminación directa e indirecta para el color final.
	colorFinal = colorDirecto + colorIndirecto;

	//Devolvemos el color final.
	return colorFinal;
}



// Returns a sample over an oriented disk. The sample is obtained from a uniform distributed random variable.
Sample Raytracer::SampleDisk( const Vec3 &N )
{
	Sample sample; //Creamos una muestra

	//Sacamos s y t como variables aleatorias sobre un cuadrado de 1x1
	float s = rand(0.0,1.0);
	float t = rand(0.0,1.0);

	//Creamos los coeficientes Lx,Ly y Lz para convertir el cuadrado en una muestra circular, y con Lz le damos profundidad, 
	//de esta manera la muestra estará comprendida en una semiesfera.
	float Lx = sqrt(t)*cos(2*Pi*s);
	float Ly = sqrt(t)*sin(2*Pi*s);
	float Lz = sqrt(1-pow(Lx,2)-pow(Ly,2));

	Vec3 vReflect = Unit((N + Vec3(0,0,1))/2);
	//Creamos un vector L que comprende los coeficientes calculador anteriormente.
	Vec3 L = Unit(Vec3(Lx,Ly,Lz));
	//Reorientamos la muestra de la semiesfera hacia la dirección de la normal del objeto intersectado.
	sample.P = Reflection(-1*L,vReflect);

	return sample;
}

// Returns a sample into the projected hemisphere. This is a type of importance sampling.
Sample Raytracer::SampleProjectedHemisphere( const Vec3 &N )
{
	Sample sample;//Creamos una muestra

	//Sacamos s y t como variables aleatorias sobre un cuadrado de 1x1
	float s = rand(0.0,1.0);
	float t = rand(0.0,1.0);

	//Creamos los coeficientes Lx,Ly y Lz para convertir el cuadrado en una muestra circular, y con Lz le damos profundidad, 
	//de esta manera la muestra estará comprendida en una semiesfera.
	float Lx = sqrt(t)*cos(2*Pi*s);
	float Ly = sqrt(t)*sin(2*Pi*s);
	float Lz = sqrt(1-pow(Lx,2)-pow(Ly,2));

	Vec3 vReflect = Unit((N + Vec3(0,0,1))/2);
	Vec3 L = Unit(Vec3(Lx,Ly,Lz));//Creamos un vector L que comprende los coeficientes calculador anteriormente.
	sample.P = Reflection(-1*L,vReflect);//Reorientamos la muestra de la semiesfera hacia la dirección de la normal del objeto intersectado.
	sample.w = Pi;//Damos un peso a la muestra para el cálculo de la luz indirecta difusa.

	return sample;
}

// Returns a sample into the specular lobe. This is a type of importance sampling.
Sample Raytracer::SampleSpecularLobe( const Vec3 &R, float phong_exp  )
{
	Sample sample;//Creamos una muestra

	//Sacamos s y t como variables aleatorias sobre un cuadrado de 1x1
	float s = rand(0.0,1.0);
	float t = rand(0.0,1.0);

	//Creamos los coeficientes Lx,Ly y Lz para convertir el cuadrado en una muestra circular, y con Lz le damos profundidad, 
	//de esta manera la muestra estará comprendida en una semiesfera.
	float Lx = sqrt(1-pow(t,(2/(phong_exp+1))))*cos(2*Pi*s);
	float Ly = sqrt(1-pow(t,(2/(phong_exp+1))))*sin(2*Pi*s);
	float Lz = sqrt(1-pow(Lx,2)-pow(Ly,2));

	Vec3 vReflect = Unit((R+Vec3(0,0,1))/2);
	Vec3 L = Unit(Vec3(Lx,Ly,Lz));//Creamos un vector L que comprende los coeficientes calculador anteriormente.
	sample.P = Reflection(-1*L,vReflect);//Reorientamos la muestra de la semiesfera hacia la dirección de la normal del objeto intersectado.
	sample.w = 2*Pi/(phong_exp + 2);//Damos un peso a la muestra para el cálculo de la luz indirecta especular.

	return sample;
}