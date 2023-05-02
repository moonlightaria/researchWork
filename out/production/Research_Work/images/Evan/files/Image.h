#pragma once
#ifndef  IMAGE_H
#define  IMAGE_H

#include "pixel.h"
#include "opencv/cv.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <queue>
#include <map>
#include <set>
#include "slic.h"
#include "simple_svg_1.0.0.hpp"
#include "Region.h"

#include <vector>
#include "ColorSpace.h"
#include "Comparison.h"
#include "Conversion.h"


using namespace svg;
using namespace cv;
// the following functions are based off of the pseudocode
// found on www.easyrgb.com

static cv::Vec3b lab2rgb(cv::Vec3f lab) {
	float y = (lab.val[0] + 16) / 116;
	float x = lab.val[1] / 500 + y;
	float z = y - lab.val[2] / 200;
	float	r; float g; float b;

	x = 0.95047 * ((x * x * x > 0.008856) ? x * x * x : (x - 16 / 116) / 7.787);
	y = 1.00000 * ((y * y * y > 0.008856) ? y * y * y : (y - 16 / 116) / 7.787);
	z = 1.08883 * ((z * z * z > 0.008856) ? z * z * z : (z - 16 / 116) / 7.787);

	r = x * 3.2406 + y * -1.5372 + z * -0.4986;
	g = x * -0.9689 + y * 1.8758 + z * 0.0415;
	b = x * 0.0557 + y * -0.2040 + z * 1.0570;

	r = (r > 0.0031308) ? (1.055 * pow(r, 1 / 2.4) - 0.055) : 12.92 * r;
	g = (g > 0.0031308) ? (1.055 * pow(g, 1 / 2.4) - 0.055) : 12.92 * g;
	b = (b > 0.0031308) ? (1.055 * pow(b, 1 / 2.4) - 0.055) : 12.92 * b;

	return cv::Vec3b (max(0, min(1,(int) r)) * 255, max(0, min(1, (int)g)) * 255, max(0, min(1, (int)b)) * 255);
}

static cv::Vec3f rgb2lab(cv::Vec3b rgb) {
	float r = rgb.val[0] / 255;
	float g = rgb.val[1] / 255;
	float b = rgb.val[2] / 255;
	float x; float y; float z;

	r = (r > 0.04045) ? pow((r + 0.055) / 1.055, 2.4) : r / 12.92;
	g = (g > 0.04045) ? pow((g + 0.055) / 1.055, 2.4) : g / 12.92;
	b = (b > 0.04045) ? pow((b + 0.055) / 1.055, 2.4) : b / 12.92;

	x = (r * 0.4124 + g * 0.3576 + b * 0.1805) / 0.95047;
	y = (r * 0.2126 + g * 0.7152 + b * 0.0722) / 1.00000;
	z = (r * 0.0193 + g * 0.1192 + b * 0.9505) / 1.08883;

	x = (x > 0.008856) ? pow(x, 1 / 3) : (7.787 * x) + 16 / 116;
	y = (y > 0.008856) ? pow(y, 1 / 3) : (7.787 * y) + 16 / 116;
	z = (z > 0.008856) ? pow(z, 1 / 3) : (7.787 * z) + 16 / 116;

	return cv::Vec3f((116 * y) - 16, 500 * (x - y), 200 * (y - z));
}

// calculate the perceptual distance between colors in CIELAB
// https://github.com/THEjoezack/ColorMine/blob/master/ColorMine/ColorSpaces/Comparisons/CmcComparison.cs
static float deltaE(cv::Vec3f labA, cv::Vec3f labB) {
	float deltaL = labA.val[0] - labB.val[0];
	float deltaA = labA.val[1] - labB.val[1];
	float deltaB = labA.val[2] - labB.val[2];
	float c1 = sqrt(labA.val[1] * labA.val[1] + labA.val[2] * labA.val[2]);
	float c2 = sqrt(labB.val[1] * labB.val[1] + labB.val[2] * labB.val[2]);
	float deltaC = c1 - c2;
	float deltaH = deltaA * deltaA + deltaB * deltaB - deltaC * deltaC;
	deltaH = deltaH < 0 ? 0 : sqrt(deltaH);
	float sc = 1.0 + 0.045 * c1;
	float sh = 1.0 + 0.015 * c1;
	float deltaLKlsl = deltaL / (1.0);
	float deltaCkcsc = deltaC / (sc);
	float deltaHkhsh = deltaH / (sh);
	float i = deltaLKlsl * deltaLKlsl + deltaCkcsc * deltaCkcsc + deltaHkhsh * deltaHkhsh;
	return i < 0 ? 0 : sqrt(i);
}

static int compress(int below, int tot)
{
	// compress "below" count toward 50%
	float compressratio = 0.5f;
	float rat = below / (tot + 0.5f);

	if (rat > 0.5f) rat = 0.5f + (rat - 0.5f)*compressratio;
	else
		rat = 0.5f - (0.5f - rat)*compressratio;

	int nrat = (int)(tot*rat);
	return nrat;
}

static int lookoutbelow(int *hist, int me)
{
	// find out how many are lte slot "me".
	// note, counts number "equal"
	// further note, returns absolute count, not percentile

	int i;
	int numbelow = 0;
	for (i = 0; i <= me; i++)
	{
		numbelow += hist[i];
	}

	return numbelow;
}

static int whereinhist(int *hist, int numbelow)
{
	// takes in histogram (256-element array of counts) and target, returns
	// index where numbelow elements are below
	//
	// note, target is an absolute count, not a percentile
	//

	int cbelow = 0;
	int i = 0;
	while ((cbelow <= numbelow) && (i < 255))
	{
		cbelow += hist[i];
		i++;
	}
	return i;
}



namespace abstraction {

	class Image {

		friend class Region;

	public:

		Image();
		Image(cv::Mat &image);
		cv::Mat img;
		cv::Mat baseImage;
		cv::Mat gray;
		cv::Mat LabImage;

		IplImage* image2;// keep slic

		std::vector<Region*> AllReions,AllRegions;
		std::vector <std::vector <int> > regionsIDs;
		std::vector<Region*> BackGround;

		pixel *allpixels[ROWS*COLS];
		std::vector<pixel*> allPix;
		int rows;
		int cols;
		int numnodes;
		int mask;

		int numSlic;
		int maxLev;
		int maxRes = 0;

		cv::Mat filtered,filtered10;

		std::string path, filteredPath,filtered10Path;
		cv::Mat  nres;
		cv::Mat  pres;
		cv::Mat colorRes;

		cv::Mat signedRes960, signedRes480,signedRes160, signedRes40, signedRes;
		cv::Mat  gr_res;

		std::vector<std::vector<int>> slicIDs; // slic id per pixel

		std::vector < std::vector<pixel*>> blobsP,blobTheRest, nutBlob960, blobsN;

		std::vector < std::vector<pixel*>> slicBlobs;

		std::vector<pixel*> seeds,seeds2,seeds3; //seeds = bodycolor seeds/ seeds2 = maxvalue /seeds3 = medianloc

		std::vector<std::vector<pixel*>> regions, collectedRegions;
		std::vector<pixel*>  extremas;


		void initialize(cv::Mat &image, bool colour);
		void FindBlobs1(const cv::Mat &binary, std::vector < std::vector<pixel*> > &blobs, cv::Mat &blobdraw);
		void calcSLIC(int mask, int numsp,cv::Mat &draw);
		void FindCC(const cv::Mat & visible,/* std::vector<std::vector<pixel*>>& blobs,*/ int & maxLabel);
		void floodFill(cv::Mat & image, int sr, int sc, int newColor, cv::Vec3b oldColor);
		void Flood_fill(cv::Mat & image, int sr, int sc, int newColor);
		void separateConnectedRegions(cv::Mat &visible, int &maxLabel);
		void DFSUtil(int v,/* bool visited[],*/ int label);
		void separateConnectedRegions2(int &label);
		void createRegions(int mask,  std::string str);
		void makeGraphonRegions(cv::Mat& img, std::string str);
		void colorRegions();
		void mergingSmallRegionsOfSize(int sz);
		void mergingSmallRegions();
		void assignSLICid(std::vector<std::vector<int>> SLICids);		
		void calcKmean_reColor(std::vector<pixel*>& regionpix,/* std::vector<cv::Point > contours,*/ cv::Mat &draw, int K= 2);
		void calcMahalanobisAB(std::vector<cv::Vec3f> initDistribution, std::vector<cv::Vec3f> distribution, std::vector<float> &covarElement, double & m_dist);
		void calcMahalanobis(cv::Vec3f point,std::vector<cv::Vec3f> distribution,double &m_dist);
		void colorMode(std::vector<pixel*> &region, cv::Mat &screen, cv::Vec3b &mode);
		void colorAverage(std::vector<pixel*> &region, cv::Mat &screen, cv::Vec3b &ave);
	
		// assign slic ids from slic operator for each pixel
		void calcResiduals(cv::Mat &f,cv::Mat &colorRes, cv::Mat &pres, cv::Mat &nres, cv::Mat &signedres);
		void getConnectedComponents(int val, void*);
		void SeedGeneration(int mask, Image * image, const char * filename);	
		void generateSeeds(std::vector<pixel*> colorBin[100], bool flag); // falg = pos/neg res
		void calculate_cdf(std::vector<int> hist, std::vector<float>& normalized_cdf);
		void HistogramMatching(std::vector<float> grays ,std::vector<int>& lookUp_table);
		
		void FindBlobs(const cv::Mat &binary, std::vector < std::vector<pixel*> > &blobs, cv::Mat &blobdraw,int Mask, bool pos);
		void generateSeedsWithMaxvalue(std::vector<pixel*> &pixels);
		void generateMedianLocSeeds(std::vector<pixel*>& pixels);

		void makegraph();
			
		void setResiduals();	

		cv::Vec3b returnAverageColor(std::vector<pixel*>& region);

		void calcCovarMat(std::vector<cv::Vec3f> distribution, cv::Mat & covarMat, cv::Mat &mean);
		void calcKmean(int K, cv::Mat &gy);
		void lineSimpilification(std::vector<std::vector<pixel*>> &regions, cv::Mat &draw);
		void assignRegionToPixels(int mask, int &maxreg);
		cv::Vec3b PerturbInHSL(cv::Vec3b col);
		void smoothBoundary(std::vector< std::vector<cv::Point>> &contours, std::vector<pixel*> &region, cv::Mat &draw);

	};
}//namespace abstraction
#endif