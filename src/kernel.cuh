#ifndef CUDA_H
#define CUDA_H

#include <cuda_runtime.h>
#include <cuda.h>
#include "device_launch_parameters.h"
#include "Particle.h"
#include <stdio.h>


class Cuda
{
public:
	int* id;
	float* r;
	float* m;
	float* posX;
	float* posY;
	float* velX;
	float* velY;
	float* accX;
	float* accY;
	float* fx;
	float* fy;
	float* dist;
	float* rayX;
	float* rayY;
	float* pointX;
	float* pointY;
	float* pointM;

	int* idHost;
	float* rHost;
	float* mHost;
	float* posXHost;
	float* posYHost;
	float* velXHost;
	float* velYHost;
	float* accXHost;
	float* accYHost;
	float* distHost;
	float* rayXHost;
	float* rayYHost;
	float* pointXHost;
	float* pointYHost;
	float* pointMHost;

	int particlesN = 0;
	int pointsN = 0;
	int maxPointsN = 0;

	Cuda(int cudaParticlesN, int cudaPointsN, int cudaMaxPointsN, vector<Particle>* particles, vector<Particle*>* points);
	virtual ~Cuda();

	void run(vector<Particle>* particles, float particleMassMultiplier, float pointMassMultiplier, float viscosity1, float viscosity2, float dt);

	void UpdatePoint(vector<Particle*>* points);
private:

};

#endif