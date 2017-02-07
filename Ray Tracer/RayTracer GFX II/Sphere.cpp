#include "Sphere.h"

Sphere::Sphere( const Vec3 &cent, float rad )
{
    center = cent;
    radius = rad;
}


Object *Sphere::ReadString( const char *params ) // Reads params from a string.
{
    float x, y, z, r;
    if( sscanf( params, "sphere (%f,%f,%f) %f", &x, &y, &z, &r ) == 4 )
        return new Sphere( Vec3( x, y, z ), r );
    return NULL;
}

Box3 Sphere::GetBounds() const // Returns a bounding box.
{
    Box3 box;
    box.X.min = center.x - radius;  box.X.max = center.x + radius;
    box.Y.min = center.y - radius;  box.Y.max = center.y + radius;
    box.Z.min = center.z - radius;  box.Z.max = center.z + radius;
    return box;
}

bool Sphere::Intersect( const Ray &ray, HitGeom &hitgeom ) const
{
	//Computar los coeficientes A, B, C
	float a = (ray.direction * ray.direction);
	float b = 2 * (ray.direction * (ray.origin - center));
	float c = ((ray.origin - center)*(ray.origin - center) - (radius * radius));

	//Encontrar el discriminante de la funcion
	float discriminante = (b * b) - (4.0 * a * c);
	
	float t;
	
	//Si el discriminante es negativo, el rayo no atraviesa la geometria de la esfera
	if (discriminante < 0)
		return false;

	float distSqrt = sqrtf(discriminante);
    float q;
    if (b < 0)
        q = (-b - distSqrt)/2.0;
    else
        q = (-b + distSqrt)/2.0;

	float t0 = q / a;
    float t1 = c / q;

	if (discriminante == 0)
		t = -0.5 * b / a;
	else{

		if (t1 < 0.1)
			return false;

		if (t0 < 0.1)
			t = t1;
		else
			t = t0;

	}

		t -= Epsilon;

		if(hitgeom.distance > t){
			hitgeom.distance = t;
			hitgeom.origin = ray.origin;
			hitgeom.point = Vec3(ray.origin.x + (ray.direction.x * t),ray.origin.y + (ray.direction.y * t),ray.origin.z + (ray.direction.z * t));

			Vec3 normal = (hitgeom.point - center);
			normal /= Length(normal);

			hitgeom.normal = normal;

			return true;
		} else {
			return false;
		}
}

Sample Sphere::GetSample( const Vec3 &P, const Vec3 &N ) const
{
	Sample sample; //Creamos una muestra

	float d = Length(center-P);//Calculamos la distancia de la muestra real.
	float h = cos(asin(radius/d) - Epsilon); //Sacamos h para muestrear el peso de la muestra.
	sample.w = 2 * Pi * (1 - h);

	Vec3 PC = (center - P) / d; //Vector director
	Vec3 reflect = (Vec3(0,0,1) + PC) / 2; //Reorientamos la muestra
	reflect = reflect / Length(reflect);

	//calculamos s y t como variables para generar un cuadrado de 4x4
	float s = rand(0.0, 1.0);
	float t = rand(0.0, 1.0);
	float k = sqrt(1 - pow(h + s * (1 - h), 2));//Delimitamos los márgenes de la muestra de la esfera gracias a k.
	//Creamos un vector L como dirección hacia la esfera que emita luz.
	Vec3 L = Vec3(k * cos(2 * Pi * t), k * sin(2 * Pi * t), (1 - h) * s + h);

	Vec3 L2 = Reflection(-1*L, reflect);//Reorientamos la muestra hacia la luz.

	//Generamos un nuevo rayo con dirección L2 y punto de origen en P.
	Ray	ray;
	ray.direction = L2;
	ray.origin = P;
	HitGeom hitgeom;
	hitgeom.distance = Infinity;

	//En caso que intersecte con la geometría ponemos el punto P como punto de intersección
	if (Intersect(ray, hitgeom)){
		sample.P = hitgeom.point;
	}

	return sample;
}