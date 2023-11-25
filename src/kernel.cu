#include "Particle.h"
#include "kernel.cuh"

//calculate force points
__global__ void PointsForce(int particlesN, int maxPointsN, float pointMassMultiplier, float* posX, float* posY, float* fx, float* fy, float* pointX, float* pointY, float* pointM)
{
	int i = threadIdx.x * blockIdx.x;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < particlesN; i += stride)
	{
		float dist;

		for (int j = 0; j < maxPointsN; j++)
		{
			if (pointM[j] != 0.0f)
			{
				dist = sqrtf(pow(posX[i] - pointX[j], 2) + pow(posY[i] - pointY[j], 2));

				fx[i] += pointMassMultiplier * 0.1 * pointM[j] * (posX[i] - pointX[j]) / pow(dist, 3.0);
				fy[i] += pointMassMultiplier * 0.1 * pointM[j] * (posY[i] - pointY[j]) / pow(dist, 3.0);
			}
		}
	}
}

//calculate force between particles
__global__ void ParticlesForce(int particlesN, float particleMassMultiplier, float* posX, float* posY, float* m, float* fx, float* fy)
{
	int i = threadIdx.x * blockIdx.x;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < particlesN; i += stride)
	{
		float dist;

		for (int j = 0; j < particlesN; j++)
		{
			if (i != j)
			{
				dist = sqrtf(pow(posX[i] - posX[j], 2) + pow(posY[i] - posY[j], 2));
				fx[i] += particleMassMultiplier * 0.1 * m[j] * (posX[i] - posX[j]) / pow(dist, 3.0);
				fy[i] += particleMassMultiplier * 0.1 * m[j] * (posY[i] - posY[j]) / pow(dist, 3.0);
			}

		}
	}
}

//Update particles
__global__ void UpdateParticlesPosition(int particlesN, float particleMassMultiplier, float pointMassMultiplier, float viscosity1, float viscosity2, float dt, float* posX, float* posY, float* velX, float* velY, float* accX, float* accY, float* fx, float* fy, float* m, float* r, float* pointX, float* pointY, float* pointM)
{
	int i = threadIdx.x * blockIdx.x;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x;

	for (int i = index; i < particlesN; i += stride)
	{
		float addfx;
		float addfy;
		float dx = 0;
		float dy = 0;
		float q1 = viscosity1;
		float q2 = viscosity2;
		float dist;

		//apply velocity
		velX[i] = velX[i] + accX[i] * dt;
		velY[i] = velY[i] + accY[i] * dt;

		//max vel
		float max = 100;

		if (velX[i] > max) velX[i] = max;
		if (velX[i] < -max) velX[i] = -max;
		if (velY[i] > max) velY[i] = max;
		if (velY[i] < -max) velY[i] = -max;

		//apply possition
		posX[i] = posX[i] + velX[i] * dt;
		posY[i] = posY[i] + velY[i] * dt;

		//apply aerodynamic drag
		float gr = 0.1;
		if (velX[i] > gr || velX[i] < -gr || velY[i] > gr || velY[i] < -gr)
		{
			if (pow(posX[i] - 512.0, 2) + pow(posY[i] - 512.0, 2) >= pow(400.0, 2))
			{
				dx = -6 * 3.14 * velX[i] * q2 * r[i];
				dy = -6 * 3.14 * velY[i] * q2 * r[i];
			}
			else
			{
				dx = -6 * 3.14 * velX[i] * q1 * r[i];
				dy = -6 * 3.14 * velY[i] * q1 * r[i];
			}
		}
		else
		{
			velX[i] = 0.0;
			velY[i] = 0.0;
		}

		fx[i] -= dx;
		fy[i] -= dy;

		//apply acceleration
		accX[i] = -fx[i] / m[i];
		accY[i] = -fy[i] / m[i];

		//apply boundaries
		if (posX[i] - r[i] < 0) velX[i] = -velX[i];
		if (posX[i] + r[i] > 1024) velX[i] = -velX[i];

		if (posY[i] - r[i] < 0) velY[i] = -velY[i];
		if (posY[i] + r[i] > 1024) velY[i] = -velY[i];

		fx[i] = 0.0;
		fy[i] = 0.0;
	}
}

Cuda::Cuda(int cudaParticlesN, int cudaPointsN, int cudaMaxPointsN, vector<Particle>* particles, vector<Particle*>* points)
{
	particlesN = cudaParticlesN;
	pointsN = cudaPointsN;
	maxPointsN = cudaMaxPointsN;

	int* idHost = new int[particlesN];
	float* rHost = new float[particlesN];
	float* mHost = new float[particlesN];
	float* posXHost = new float[particlesN];
	float* posYHost = new float[particlesN];
	float* velXHost = new float[particlesN];
	float* velYHost = new float[particlesN];
	float* accXHost = new float[particlesN];
	float* accYHost = new float[particlesN];
	float* distHost = new float[particlesN];
	float* rayXHost = new float[particlesN];
	float* rayYHost = new float[particlesN];
	float* pointXHost = new float[maxPointsN];
	float* pointYHost = new float[maxPointsN];
	float* pointMHost = new float[maxPointsN];

	cudaMalloc(&id, particlesN * sizeof(int));
	cudaMalloc(&r, particlesN * sizeof(float));
	cudaMalloc(&m, particlesN * sizeof(float));
	cudaMalloc(&posX, particlesN * sizeof(float));
	cudaMalloc(&posY, particlesN * sizeof(float));
	cudaMalloc(&velX, particlesN * sizeof(float));
	cudaMalloc(&velY, particlesN * sizeof(float));
	cudaMalloc(&accX, particlesN * sizeof(float));
	cudaMalloc(&accY, particlesN * sizeof(float));
	cudaMalloc(&fx, particlesN * sizeof(float));
	cudaMalloc(&fy, particlesN * sizeof(float));
	cudaMalloc(&dist, particlesN * sizeof(float));
	cudaMalloc(&rayX, particlesN * sizeof(float));
	cudaMalloc(&rayY, particlesN * sizeof(float));
	cudaMalloc(&pointX, maxPointsN * sizeof(float));
	cudaMalloc(&pointY, maxPointsN * sizeof(float));
	cudaMalloc(&pointM, maxPointsN * sizeof(float));

	for (int i = 0; i < particlesN; i++)
	{
		idHost[i] = particles->at(i).getID();
		rHost[i] = particles->at(i).getRadius();
		mHost[i] = particles->at(i).getMass();
		posXHost[i] = particles->at(i).getPosition().x;
		posYHost[i] = particles->at(i).getPosition().y;
		velXHost[i] = particles->at(i).getVelocity().x;
		velYHost[i] = particles->at(i).getVelocity().y;
		accXHost[i] = particles->at(i).getAcceleration().x;
		accYHost[i] = particles->at(i).getAcceleration().y;
	}

	for (int i = 0; i < maxPointsN; i++)
	{
		pointXHost[i] = 0.0f;
		pointYHost[i] = 0.0f;
		pointMHost[i] = 0.0f;
	}

	for (int i = 0; i < pointsN; i++)
	{
		pointXHost[i] = points->at(i)->getPosition().x;
		pointYHost[i] = points->at(i)->getPosition().y;
		pointMHost[i] = points->at(i)->getMass();
	}

	cudaMemcpy(r, rHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(m, mHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(posX, posXHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(posY, posYHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(velX, velXHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(velY, velYHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(accX, accXHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(accY, accYHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dist, distHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(rayX, rayXHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(rayY, rayYHost, particlesN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(pointX, pointXHost, maxPointsN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(pointY, pointYHost, maxPointsN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(pointM, pointMHost, maxPointsN * sizeof(float), cudaMemcpyHostToDevice);

	free(rHost);
	free(mHost);
	free(posXHost);
	free(posYHost);
	free(velXHost);
	free(velYHost);
	free(accXHost);
	free(accYHost);
	free(distHost);
	free(rayXHost);
	free(rayYHost);
	free(pointXHost);
	free(pointYHost);
	free(pointMHost);
}

Cuda::~Cuda()
{

}

void Cuda::run(vector<Particle>* particles, float particleMassMultiplier, float pointMassMultiplier, float viscosity1, float viscosity2, float dt)
{
	int blockSize = 128;
	int numBlocks = (particlesN + blockSize - 1) / blockSize;
	PointsForce << < numBlocks, blockSize >> > (particlesN, maxPointsN, pointMassMultiplier, posX, posY, fx, fy, pointX, pointY, pointM);
	cudaDeviceSynchronize();
	ParticlesForce << < numBlocks, blockSize >> > (particlesN, particleMassMultiplier, posX, posY, m, fx, fy);
	cudaDeviceSynchronize();
	//Collision << < numBlocks, blockSize >> > (particlesN, pointsN, posX, posY, velX, velY, r, m, pointX, pointY);
	//cudaDeviceSynchronize();
	UpdateParticlesPosition << < numBlocks, blockSize >> > (particlesN, particleMassMultiplier, pointMassMultiplier, viscosity1, viscosity2, dt, posX, posY, velX, velY, accX, accY, fx, fy, m, r, pointX, pointY, pointM);
	cudaDeviceSynchronize();

	float* posXHost = new float[particlesN];
	float* posYHost = new float[particlesN];
	float* velXHost = new float[particlesN];
	float* velYHost = new float[particlesN];
	float* accXHost = new float[particlesN];
	float* accYHost = new float[particlesN];

	cudaMemcpy(posXHost, posX, particlesN * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(posYHost, posY, particlesN * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(velXHost, velX, particlesN * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(velYHost, velY, particlesN * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(accXHost, accX, particlesN * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(accYHost, accY, particlesN * sizeof(float), cudaMemcpyDeviceToHost);


	for (int i = 0; i < particlesN; i++)
	{
		particles->at(i).setPosition(v2(posXHost[i], posYHost[i]));
		particles->at(i).setVelocity(v2(velXHost[i], velYHost[i]));
		particles->at(i).setAcceleration(v2(accXHost[i], accYHost[i]));
	}

	free(posXHost);
	free(posYHost);
	free(velXHost);
	free(velYHost);
	free(accXHost);
	free(accYHost);
}

void Cuda::UpdatePoint(vector<Particle*>* points)
{
	pointsN = points->size();

	float* pointXHost = new float[maxPointsN];
	float* pointYHost = new float[maxPointsN];
	float* pointMHost = new float[maxPointsN];

	for (int i = 0; i < maxPointsN; i++)
	{
		pointXHost[i] = 0.0f;
		pointYHost[i] = 0.0f;
		pointMHost[i] = 0.0f;
	}

	for (int j = 0; j < pointsN; j++)
	{
		pointXHost[j] = points->at(j)->getPosition().x;
		pointYHost[j] = points->at(j)->getPosition().y;
		pointMHost[j] = points->at(j)->getMass();
	}

	int blockSize = 256;
	int numBlocks = (particlesN + blockSize - 1) / blockSize;

	cudaMemcpy(pointX, pointXHost, maxPointsN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(pointY, pointYHost, maxPointsN * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(pointM, pointMHost, maxPointsN * sizeof(float), cudaMemcpyHostToDevice);

	free(pointXHost);
	free(pointYHost);
	free(pointMHost);
}
