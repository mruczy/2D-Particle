#ifndef PARTICLE_H
#define	PARTICLE_H

#include "ofMain.h"
#include "v2.h"

class Particle
{
private:
	int id;			//particle id
	float r;		//radius
	float m;		//mass
	int cr;			//color r
	int cg;			//color g
	int cb;			//color b
	v2 pos;			//position vector
	v2 vel;			//veloticy vector
	v2 acc;			//acceleration vector

public:
	Particle();

	//position , velocity, acceleration, radius, mass, color 
	Particle(v2 posvec, v2 velvec, v2 accvec, float radius, float mass,int colorR, int colorG, int colorB);
	virtual ~Particle();

	int getID();
	void setID(int id);

	float getRadius();
	void setRadius(float radius);

	float getMass();
	void randomMass(float minMass, float maxMass);
	void setMass(float mass);

	v2 getPosition();
	void setPosition(v2 posvec);

	v2 getVelocity();
	void setVelocity(v2 velvec);

	v2 getAcceleration();
	void setAcceleration(v2 accvec);

	int getColorR();
	int getColorG();
	int getColorB();
	void setColor(int colorR, int colorG, int colorB);

	void update(float dt);
	void draw();
};

#endif