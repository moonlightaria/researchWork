#pragma once
#ifndef REGION_H
#define REGION_H

#include "crack.h"
#include <iostream>
#include <vector>
#include "pixel.h"
#include "Image.h"


int compress(int below, int tot);
int lookoutbelow(int *hist, int me);
int whereinhist(int *hist, int numbelow);
static bool x_compare(pixel* a, pixel* b) { // return true if a < b
	return a->getPos().x < b->getPos().x;
}
static bool y_compare(pixel* a, pixel* b) {
	return a->getPos().y < b->getPos().y;
}

namespace abstraction {
	struct threeObj {
		cv::Vec3b first;
		cv::Vec3b second;
		float dif;
	};
	static bool smallestDiff(threeObj  pair1, threeObj  pair2)
	{
		return  pair1.dif <  pair2.dif;
	}
	static bool smallestweight(std::pair<float, int>  w1, std::pair<float, int> w2  )
	{
		return  w1.first <  w2.first;
	}
	static bool largestweight(std::pair<float, int>  w1, std::pair<float, int> w2)
	{
		return  w1.first >  w2.first;
	}
	struct RegionEdge {
		int first;
		int second;
		float weight;
		//int getTheOtherSide(Region *p) {
		//	return p->getId() == first ? second : first;
		//}
	};

	class Region {
		friend class pixel;
		friend class Image;

	private:
		int regId = NULL;
	public:
		std::vector<pixel*> regionPix;
		std::vector<pixel*> contourPix;
		std::vector<int> neighbours; // neighbour regions
		std::vector<cv::Vec3b>assigned_colors; // list of assigned colors
		Region* parent;
		std::vector<Region*> children;
		int finalRegId;
		int gray;
		int prior;
		cv::Vec3b median;
		cv::Vec3b mode;
		cv::Vec3b aveColor;
		cv::Vec3b newColor;
		cv::Vec3b BdVecotrs;
		std::vector<cv::Point3d> r__eigen_vecs;
		std::vector<double> r__eigen_val;
		cv::Point3d cntr;
		double theta;
		double phi;
		bool colined;
		bool colored;
		bool visited;
		float dist;
		int minx, maxx, miny, maxy = NULL;
		std::vector<RegionEdge*> edges;

		Region();
		Region(std::vector<pixel*> regPix);
		~Region();
		void setNewColor(cv::Vec3b newCol);
		cv::Vec3b getNewColor();
		int getId();
		void setRegId(int _rid);
		void calcBoundingbox();
		void addElements(pixel* p);
		bool DoBoxesIntersect(Region *B);
		void removeDuplicates(std::vector<int> &v); // remove duplicates
		cv::Vec3b calcClosestToMode();
		cv::Vec3b calcMedianColor();
		cv::Vec3b calcAveColor(/*cv::Mat &img*/);
		int ReturnGray(/*cv::Mat &image*/);
		void SetGray(int gr);
		void CalcPriority();
		void sortNeighbours();
		void RecolorTo(cv::Vec3b _newC);
		void checkVisibility(std::vector<pixel*> allPix);
		
	};

	
} // namespace abstraction
#endif // REGION_H