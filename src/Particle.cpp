#include "Particle.h"

Particle::Particle()
{
	this->setPosition(v2(0,0));
	this->setVelocity(v2(0, 0));
	this->setAcceleration(v2(0, 0));
	this->setRadius(0);
	this->setMass(0);
	this->setColor(0, 0, 0);
}

Particle::Particle(v2 posvec, v2 velvec, v2 accvec, float radius, float mass, int colorR, int colorG, int colorB)
{
	this->setPosition(posvec);
	this->setVelocity(velvec);
	this->setAcceleration(accvec);
	this->setRadius(radius);
	this->setMass(mass);
	this->setColor(colorR, colorG, colorB);
}

Particle::~Particle()
{
}

int Particle::getID()
{
	return id;
}

void Particle::setID(int idIn)
{
	id = idIn;
}

float Particle::getRadius()
{
	return r;
}

void Particle::setRadius(float radius)
{
	r = radius;
}

float Particle::getMass()
{
	return m;
}

void Particle::randomMass(float minMass, float maxMass)
{
	if (maxMass > minMass && minMass > 0.0)
	{
		float tempMass = minMass + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (maxMass - minMass)));
		this->setMass(tempMass);
	}
}

void Particle::setMass(float mass)
{
	m = mass;
}

v2 Particle::getPosition()
{
	return pos;
}

void Particle::setPosition(v2 posvec)
{
	pos.x = posvec.x;
	pos.y = posvec.y;
}

v2 Particle::getVelocity()
{
	return vel;
}

void Particle::setVelocity(v2 velvec)
{
	vel.x = velvec.x;
	vel.y = velvec.y;
}

v2 Particle::getAcceleration()
{
	return acc;
}

void Particle::setAcceleration(v2 accvec)
{
	acc.x = accvec.x;
	acc.y = accvec.y;
}


int Particle::getColorR()
{
	return cr;
}

int Particle::getColorG()
{
	return cg;
}

int Particle::getColorB()
{
	return cb;
}

void Particle::setColor(int colorR, int colorG, int colorB)
{
	cr = colorR;
	cg = colorG;
	cb = colorB;
}

void Particle::update(float dt)
{
	pos.x = pos.x + vel.x * dt;
	pos.y = pos.y + vel.y * dt;

	vel.x = vel.x + acc.x * dt;
	vel.y = vel.y + acc.y * dt;
	this->setPosition(v2(pos.x, pos.y));
	this->setVelocity(v2(vel.x, vel.y));
}


void Particle::draw()
{
	if (m > 0)
	{
		ofSetColor(this->getColorR(), this->getColorG(), this->getColorB());
		ofDrawCircle(this->getPosition().x, this->getPosition().y, this->getRadius());
	}
}
