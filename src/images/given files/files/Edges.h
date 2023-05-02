#ifndef EDGES_H
#define EDGES_H

#include <iostream>
#include "pixel.h"

class Edges {

	pixel *first;
	pixel *second;

	float weight;
	float geodist;
	float colorDist;
	float signDist;

public:
	Edges() {};
	Edges(pixel *p1, pixel *p2) {
		first = p1;
		second = p2;

		weight = 0.0f;
		colorDist = 0.0f;
		geodist = 0.0f;

		calcWeight();
		calcColorDist();
		calcGeodist();//  geodistance
		//calcsignDist();// calc difference of signs 
	};


	pixel *getTheOtherSide(pixel *p) {
		return p == first ? second : first;
	}


	void calcWeight() // spacial distance: P
	{ // spacial distance 
		weight = (float)first->distTo(second);
	}
	void calcColorDist() { // color distance 
		colorDist = first->colorDistTo(second);
	}
	void calcsignDist() { // sign distance 
		signDist = first->SignDistTo(second);
	}
	void calcGeodist() // sqrt(P2 + C2)
	{
		geodist = sqrt(0.000001f + colorDist); // sqrt (X2 + C2)
	}
	float getWeight() const {
		return (geodist  );
	}
	float colordistance()  {
		return (colorDist);
	}

	float getEclidWeight()  {
		return (weight);
	}

};
struct Comparator {

	bool operator()(const Edges* e1, const Edges* e2)
	{
		return (e1->getWeight() > e2->getWeight());
	}
};

#endif