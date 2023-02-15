#include "ofApp.h"
#include "kernel.cuh"
#include <random>

float rd(float min, float max)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distr(min, max);
	return distr(gen);
}

//--------------------------------------------------------------
void ofApp::setup() {
	int counter = 0;
	for (int iy = 16; iy < 1024; iy += 32)
	{
		for (int ix = 16; ix < 1024; ix += 32)
		{

			float temp = rd(1.0, 1.0);
			//                    v2(posvec)  v2(velvec)                  v2(accvec)    r     m     color    r g b   
			Particle p = Particle(v2(ix, iy), v2(rd(1.0, 3.0), rd(1.0, 3.0)), v2(0.0, 0.0), temp * 3, temp * 100, rd(1, 1), rd(1, 255), rd(1, 255));
			particles.push_back(p);
			counter++;
		}
	}

	std::cout << counter << std::endl;

	particlesN = particles.size();

	this->points.push_back(new Particle(v2(512.0, 512.0), v2(0.0, 0.0), v2(0.0, 0.0), 10.0, 100000.0, 225, 225, 225));

	this->pointsN = points.size();

	Cuda c = Cuda(particlesN, pointsN, maxPointsN, &particles, &points);
	cuda.push_back(c);

	gui.setup();
	gui.getGroup("mass multiplier").add(particleMassMultiplier.setup("particle", 0, 0, 10000));
	gui.getGroup("mass multiplier").add(pointMassMultiplier.setup("point", 0, 0, 10000));
	gui.getGroup("viscosity").add(viscosity1.setup("1", 0.1, 0.0, 2.0));
	gui.getGroup("viscosity").add(viscosity2.setup("2", 0.2, 0.0, 2.0));
	gui.add(dt.setup("dt", 0.1, 0.01, 1));
}

//--------------------------------------------------------------
void ofApp::update() {
	if (updateP)
	{
		cuda[0].UpdatePoint(&points);
		updateP = false;
	}
	cuda[0].run(&particles, particleMassMultiplier, pointMassMultiplier, viscosity1, viscosity2, dt);
}

//--------------------------------------------------------------
void ofApp::draw() {
	ofSetColor(255);

	for (int i = 0; i < particlesN; i++)
	{
		particles[i].draw();
	}

	for (int j = 0; j < pointsN; j++)
	{
		points[j]->draw();
	}
	gui.draw();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button) {
	//add point by LMB
	if (button == 0)
	{
		if (points.size() < 10)
		{
			this->points.push_back(new Particle(v2(x, y), v2(0.0, 0.0), v2(0.0, 0.0), 10.0, 100000.0, 225, 225, 225));
			pointsN = this->points.size();

			updateP = true;
		}
		else
		{
			addP = true;
			for (int i = 0; i < points.size(); i++)
			{
				if (this->points[i]->getMass() == 0.0)
				{
					if (addP == true)
					{
						points[i] = new Particle(v2(x, y), v2(0.0, 0.0), v2(0.0, 0.0), 10.0, 100000.0, 225, 225, 225);
					}
					addP = false;
				}
			}
		}
	}

	//erase point by RMB
	if (button == 2)
	{
		for (int i = 0; i < points.size(); i++)
		{
			float dist = sqrt(pow(this->points[i]->getPosition().x - x, 2) + pow(this->points[i]->getPosition().y - y, 2));
			if (dist <= this->points[i]->getRadius())
			{
				delete this->points[i];
				this->points[i]->setMass(0.0);
				updateP = true;
			}
		}
	}
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h) {

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg) {

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo) {

}