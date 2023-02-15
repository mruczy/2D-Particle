#pragma once

#include "ofMain.h"
#include "Particle.h"
#include <vector>
#include "kernel.cuh"
#include "ofxGui.h"

class ofApp : public ofBaseApp {

public:

	int particlesN = 0;
	int pointsN = 0;
	int maxPointsN = 10;
	bool updateP = false;
	bool addP = false;

	std::vector<Cuda> cuda;

	std::vector<Particle> particles;
	std::vector<Particle*> points;

	ofxPanel gui;
	ofxFloatSlider particleMassMultiplier;
	ofxFloatSlider pointMassMultiplier;
	ofxFloatSlider viscosity1;
	ofxFloatSlider viscosity2;
	ofxFloatSlider dt;

	void setup();
	void update();
	void draw();

	void keyPressed(int key);
	void keyReleased(int key);
	void mouseMoved(int x, int y);
	void mouseDragged(int x, int y, int button);
	void mousePressed(int x, int y, int button);
	void mouseReleased(int x, int y, int button);
	void mouseEntered(int x, int y);
	void mouseExited(int x, int y);
	void windowResized(int w, int h);
	void dragEvent(ofDragInfo dragInfo);
	void gotMessage(ofMessage msg);
};