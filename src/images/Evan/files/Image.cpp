#include "Image.h"
#include "Edges.h"
#include <iostream>
#include <random>
#include <string>
#include <iterator>
#include <algorithm>
#include <list>
#include "heap.h"
#include <opencv2/opencv.hpp>
#include "Region_Graph.h"

namespace abstraction {

	int xx = 0; int yy = 0;

	static float Min(float a, float b) {
		return a <= b ? a : b;
	}
	static float Max(float a, float b) {
		return a >= b ? a : b;
	}
	static float HueToRGB(float v1, float v2, float vH) {
		if (vH < 0)
			vH += 1;

		if (vH > 1)
			vH -= 1;

		if ((6 * vH) < 1)
			return (v1 + (v2 - v1) * 6 * vH);

		if ((2 * vH) < 1)
			return v2;

		if ((3 * vH) < 2)
			return (v1 + (v2 - v1) * ((2.0f / 3) - vH) * 6);

		return v1;
	}
	
	bool wayToSort(std::vector<cv::Point3d> i, std::vector<cv::Point3d> j)
		{	return i.size() > j.size(); }

	bool ToSort(std::vector<std::vector<pixel*>> a, std::vector<std::vector<pixel*>> b)
	{
		for (int i = 0; i < a.size();i++) {
			if(a[i].size() > b[i].size())
				return a[i] > b[i];
			else
				return a[i] < b[i];
		}
	}

	bool findId(pixel *pix, std::vector<pixel*> pixels) {
		bool flag = false;
		for (int i = 0; i < pixels.size(); i++) {
			if (pix->getId() == pixels[i]->getId())
				flag = true;
		}
		return flag;
	}

	bool checkALLIDS(std::vector<int> newIds) {
		bool flag = false;
		int count = 0;
		for (int i = 0; i < newIds.size(); i++) {
			for (int j = 0; j < newIds.size(); j++) {
				if (newIds[j] == newIds[i])
				{
					count++;
				}
			}
		}
		if (count == newIds.size() * newIds.size()) {
			flag = true;
		}
		return flag;
	}
	
	int findMode(std::vector<int> data) {

		int number = data[0];
		int mode = number;
		int count = 1;
		int countMode = 1;

		for (int i = 1; i < data.size(); i++)
		{
			if (data[i] == number)
			{ // count occurrences of the current number
				++count;
			}
			else
			{ // now this is a different number
				if (count > countMode)
				{
					countMode = count; // mode is the biggest ocurrences
					mode = number;
				}
				count = 1; // reset count for the new number
				number = data[i];
			}
		}
		return mode;
	}

	static double CalcMHWScore(std::vector<int> scores){
		size_t size = scores.size();
		if (size == 0)
		{
			return 0;  // Undefined, really.
		}
		else
		{
			sort(scores.begin(), scores.end());
			if (size % 2 == 0)
			{
				return (scores[size / 2 - 1] + scores[size / 2]) / 2;
			}
			else
			{
				return scores[size / 2];
			}
		}
	}

	int Calcmode(std::vector<int> scores, int K) {
		std::vector<int> histogram(K);
		for (int i = 0; i < K; ++i)
			for (int j = 0; j < scores.size(); j++) {
				if (scores[j] == i)
				{
					++histogram[i];
				}
			}
		// return bin id to find the size of the region
		return int(std::max_element(histogram.begin(), histogram.end()) -histogram.begin());
	}
	
	cv::Mat mkKernel(int ks, double sig, double th, double lm, double ps)
	{
		int hks = (ks - 1) / 2;
		double theta = th*CV_PI / 180;
		double psi = ps*CV_PI / 180;
		double del = 2.0 / (ks - 1);
		double lmbd = lm;
		double sigma = sig / ks;
		double x_theta;
		double y_theta;
		cv::Mat kernel(ks, ks, CV_32F);
		for (int y = -hks; y <= hks; y++)
		{
			for (int x = -hks; x <= hks; x++)
			{
				x_theta = x*del*cos(theta) + y*del*sin(theta);
				y_theta = -x*del*sin(theta) + y*del*cos(theta);
				kernel.at<float>(hks + y, hks + x) = (float)exp(-0.5*(pow(x_theta, 2) + pow(y_theta, 2)) / pow(sigma, 2))* cos(2 * CV_PI*x_theta / lmbd + psi);
			}
		}
		return kernel;
	}
	

	float distance(int x, int y, int i, int j) {
		return float(sqrt(pow(x - i, 2) + pow(y - j, 2)));
	}
	float gaussian(float x, double sigma) {
		return exp(-(pow(x, 2)) / (2 * pow(sigma, 2))) / (2 * CV_PI * pow(sigma, 2));
	}
	
	struct CompareDist
	{
		bool operator()(std::pair<float, pixel*> pairs1, std::pair<float, pixel*> pairs2) 
		{
			return pairs1.first > pairs2.first; // > for smallest on top. < for biggest on top
		}
	};
	struct CompareRegions
	{
		bool operator()(std::vector<pixel*> &region1, std::vector<pixel*> &region2)
		{
			return region1.size() < region2.size(); // > for smallest on top. < for biggest on top
		}
	};
	struct Compare
	{
		bool operator()(pixel *n1,  pixel *n2) const
		{
			return n1->dist > n2->dist;
		}
	};

	Image::Image() {

	}

	Image::Image(cv::Mat &image) {
		/****initialization***/
		bool colour;
		if (image.channels() == 1) {
			img = cv::Mat::zeros(image.rows, image.cols, CV_8UC1);
			colour = false;
		}
		rows = image.rows;
		cols = image.cols;

		if (image.channels() == 3) {
			img = cv::Mat::zeros(image.rows, image.cols, CV_8UC3);
			colour = true;
		}
		gray = cv::Mat::zeros(image.rows, image.cols, CV_8UC1);
		LabImage = cv::Mat::zeros(image.rows, image.cols, CV_8UC3);

		cvtColor(image, gray, cv::COLOR_BGR2GRAY);
		cvtColor(image, LabImage, CV_BGR2Lab);
		
		img = image.clone();
	
		initialize(image,colour);
	}

	void Image::makegraph() {
		/********make a 8-connectivity graph on all image pixels******/
		Edges *temp;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (j + 1 < cols)
				{		
					temp =  new Edges(allPix[j + i* img.cols], allPix[(j+1) + i* img.cols]);
					allPix[j + i* img.cols]->edges[RighT] = temp; //right =6
					allPix[(j+1) + i* img.cols]->edges[LEFT] = temp; //left = 2
				}
				if (j - 1 >= 0 )
				{
					if (i - 1 >= 0) {
						temp = new Edges(allPix[j + i* img.cols], allPix[(j - 1) + (i - 1)* img.cols]);
						allPix[j + i* img.cols]->edges[TOPL] = temp;
						allPix[(j - 1) + (i - 1)* img.cols]->edges[DOWNR] = temp;
					}
					if (i + 1 < rows) {
						temp = new Edges(allPix[j + i* img.cols], allPix[(j - 1) + (i + 1)* img.cols]);
						allPix[j + i* img.cols]->edges[DOWNL] = temp;
						allPix[(j - 1) + (i + 1)* img.cols]->edges[TOPR] = temp;
					}
				}
				if (i + 1 < rows)
				{
					temp = new Edges(allPix[j + i* img.cols], allPix[j + (i+1)* img.cols]);
					allPix[j + i* img.cols]->edges[DOWN] = temp; //down = 0
					allPix[j + (i + 1)* img.cols]->edges[TOP] = temp; //top =4
				}

			}
		}

	 
		for (int ii = 0; ii < allPix.size();ii++) {
			for (unsigned int i = 0; i < 8; i++) {
				if (allPix[ii]->edges[i] != 0) {
					int nid = allPix[ii]->edges[i]->getTheOtherSide(allPix[ii])->id;
					allPix[ii]->neighbours.push_back(allPix[nid]);
				}
			}
		}
	}

	void Image::setResiduals() 
	{
		/*****both residuals and intensities are set herefor both 1 channel nd 3 channels images******/

		cv::Mat sres = cv::Mat::zeros(img.size(), CV_32F);
		sres = signedRes.clone();
		int max_res = -500;
		if (!img.data)
		{
			std::cout << "Failed to load file " << std::endl;
			assert(false);
		}
		if (img.channels() == 1) {
			for (int r = 0; r < img.rows; r++)
				for (int c = 0; c < img.cols; c++)
				{
					uchar color;
					color = img.at<uchar>(r, c);

					int intensity = color;
					allPix[c + r* img.cols]->setIntensity(intensity);
					allPix[c + r* img.cols]->setResidual(sres.at<float>(r, c));
				}
		}
		else {

			for (int r = 0; r < img.rows; r++)
				for (int c = 0; c < img.cols; c++)
				{
					cv::Vec3b color;
					color = img.at<cv::Vec3b>(r, c);
					allPix[c + r* img.cols]->setColor(color);
					allPix[c + r* img.cols]->setIntensity(gray.at<uchar>(r, c));
					allPix[c + r* img.cols]->setResidual(sres.at<float>(r, c));
					if (abs(sres.at<float>(r, c)) > max_res)
						max_res = sres.at<float>(r, c);
				}
		}
		maxRes = max_res;
		YD = img.rows;
		XD = img.cols;
	}

	void Image::initialize(cv::Mat &image, bool colour) {
		/****image position, color and indexes are set here *****/

		filtered = cv::Mat::zeros(image.rows, image.cols, CV_8UC3);
		colorRes = cv::Mat::zeros(image.rows, image.cols, CV_8UC3);
	
		nres = cv::Mat::zeros(image.rows, image.cols, CV_8U);
		pres = cv::Mat::zeros(image.rows, image.cols, CV_8UC1);

		gr_res = cv::Mat::zeros(image.rows, image.cols, CV_8UC1);
		signedRes = cv::Mat::zeros(image.size(), CV_32F);

		numnodes = 0;
		for (int r = 0; r < image.rows; r++) {
			for (int c = 0; c < image.cols; c++) {
				pixel *pix = new pixel();
				if (image.channels() == 1) {

					pix->setPos(cv::Point2f(c,r));
					pix->setIntensity(gray.at<uchar>(r, c));
					pix->setId(c + r * gray.cols);
				}
				if (image.channels() == 3) {

					pix->setPos(cv::Point2f(c, r));
					pix->setColor(image.at<cv::Vec3b>(r, c));
					pix->setLabColor(LabImage.at<cv::Vec3b>(r, c));
					pix->setId(c + r * image.cols);
					pix->setIntensity(gray.at<uchar>(r, c));
					
				}


				allPix.push_back(pix);
				allpixels[numnodes] = pix;
				numnodes++;
			}
		}
	}
	
	double PerpendicularDistance(const cv::Point2d &pt, const cv::Point2d &lineStart, const cv::Point2d &lineEnd)
	{
		double dx = lineEnd.x - lineStart.x;
		double dy = lineEnd.y - lineStart.y;

		//Normalise
		double mag = pow(pow(dx, 2.0) + pow(dy, 2.0), 0.5);
		if (mag > 0.0)
		{
			dx /= mag; dy /= mag;
		}

		double pvx = pt.x - lineStart.x;
		double pvy = pt.y - lineStart.y;

		//Get dot product (project pv onto normalized direction)
		double pvdot = dx * pvx + dy * pvy;

		//Scale line direction vector
		double dsx = pvdot * dx;
		double dsy = pvdot * dy;

		//Subtract this from pv
		double ax = pvx - dsx;
		double ay = pvy - dsy;

		return pow(pow(ax, 2.0) + pow(ay, 2.0), 0.5);
	}

	void RamerDouglasPeucker(const std::vector<cv::Point> &pointList, double epsilon, std::vector<cv::Point> &out)
	{
		if (pointList.size() < 2) {
			//throw invalid_argument("Not enough points to simplify");			
		}
		else {
			// Find the point with the maximum distance from line between start and end
			double dmax = 0.0;
			size_t index = 0;
			size_t end = pointList.size() - 1;
			for (size_t i = 1; i < end; i++)
			{
				double d = PerpendicularDistance(pointList[i], pointList[0], pointList[end]);
				if (d > dmax)
				{
					index = i;
					dmax = d;
				}
			}

			// If max distance is greater than epsilon, recursively simplify
			if (dmax > epsilon)
			{
				// Recursive call
				std::vector<cv::Point> recResults1;
				std::vector<cv::Point> recResults2;
				std::vector<cv::Point> firstLine(pointList.begin(), pointList.begin() + index + 1);
				std::vector<cv::Point> lastLine(pointList.begin() + index, pointList.end());
				RamerDouglasPeucker(firstLine, epsilon, recResults1);
				RamerDouglasPeucker(lastLine, epsilon, recResults2);

				// Build the result list
				out.assign(recResults1.begin(), recResults1.end() - 1);
				out.insert(out.end(), recResults2.begin(), recResults2.end());
				if (out.size() < 2)
					throw runtime_error("Problem assembling output");
			}
			else
			{
				//Just return start and end points
				out.clear();
				out.push_back(pointList[0]);
				out.push_back(pointList[end]);
			}
		}
	}

	//Ramer-Douglas-Peucker line simplification
	//https://rosettacode.org/wiki/Ramer-Douglas-Peucker_line_simplification#C.2B.2B

	void Image::lineSimpilification(std::vector<std::vector<pixel*>> &regions, cv::Mat &draw)
	{
	//	cv::Vec3b mode;
		cv::Mat draw0 = cv::Mat::zeros(img.size(), CV_8UC3); //draw0 = img.clone();
		cv::Mat draw1 = cv::Mat::zeros(img.size(), CV_8UC3); draw.copyTo(draw1);
	//	cv::Mat draw1(img.size(), CV_8UC3, Scalar(255, 255, 255));//= cv::Mat::zeros(img.size(), img.type());;//
		std::vector<std::vector<cv::Point>> Out;
		for (int i = 0; i < regions.size(); i++) {

			if (regions[i].size()> 0) 
			{
			    cv::Vec3b mode;
				colorMode(regions[i], draw0, mode);
				//colorAverage(regions[i], draw0, mode);
				cv::Mat contourImage = cv::Mat::zeros(img.size(), img.type());
				cv::Mat gbp = cv::Mat::zeros(img.size(), CV_8UC1);
				std::vector<cv::Vec4i> hierarchy;
				std::vector< std::vector<cv::Point>> contours, contours3;
				std::vector<cv::Point > contours2;
				for (int j = 0; j < regions[i].size(); j++) {
					int r = regions[i][j]->getPos().y;
					int c = regions[i][j]->getPos().x;
					contourImage.at<cv::Vec3b>(r, c) = regions[i][j]->getColor();
				}
				cvtColor(contourImage, gbp, cv::COLOR_BGR2GRAY);
				findContours(gbp, contours, hierarchy, CV_RETR_EXTERNAL /*CV_RETR_CCOMP*/, CV_CHAIN_APPROX_NONE, cv::Point(0, 0));

				if (contours[0].size() > 5) {
					contours2.resize(contours[0].size());
					approxPolyDP(Mat(contours[0]), contours2, 4, true); //Douglas-Peucker algorithm
					if (contours2.size() >5) {
						contours3.push_back(contours2);
						polylines(draw1, contours2, true, Scalar(mode.val[0], mode.val[1], mode.val[2]), 1);
						if(contours3[0].size()> 5)
							fillPoly(draw1, contours3, Scalar(mode.val[0], mode.val[1], mode.val[2]));
					}
				}
				contours.clear(); contours2.clear(); contours3.clear();
				imwrite(path +"SimplifiledContours.png",draw1);
			}
		}
		draw = draw1.clone();			
	}

	void Image::assignRegionToPixels(int m, int &maxreg)
	{
		//Note: most pixels already have a reggionid
		// finf the CC of pixels with no id
		maxreg++;
		cv::Mat mask = cv::Mat::zeros(img.size(),CV_8UC1);
		cv::Mat draw = cv::Mat::zeros(img.size(), CV_8UC3);
		cv::Mat draw2 = cv::Mat::zeros(img.size(), CV_8UC3);

		std::vector<std::vector<pixel*>> holes;
		for (int i = 0; i < allPix.size();i++) {
			if (allPix[i]->getRegId() < 0) {
				int r = allPix[i]->pos.y; int c = allPix[i]->pos.x;
				mask.at<uchar>(r,c) = 1;
				draw2.at<cv::Vec3b>(r, c) = cv::Vec3b(0,255,00);
			}
		}
		FindBlobs(mask,holes,draw,m,NULL);
		for (int i = 0; i < holes.size();i++) {
			for (int j = 0; j < holes[i].size(); j++) {
				holes[i][j]->setRegionId( maxreg);
			}
			maxreg++;
		}
		maxreg--;
		std::cout << "maxreg on the background layer :" << maxreg << std::endl;
		//imwrite(path + "holes.png",draw2);
	}

	cv::Vec3b Image::PerturbInHSL(cv::Vec3b col) {
		float perturb = (1)*(0.1);
		cv::Vec3f hsl = cv::Vec3f(0.0, 0.0, 0.0); //H,S,L

		float r = (col.val[0] / 255.0f);
		float g = (col.val[1] / 255.0f);
		float b = (col.val[2] / 255.0f);

		float min = Min(Min(r, g), b);
		float max = Max(Max(r, g), b);
		float delta = max - min;

		hsl.val[2] = ((max + min) / 2);

		if (delta == 0)
		{
			hsl.val[0] = 0;
			hsl.val[1] = 0.0f;
		}
		else
		{
			hsl.val[1] = (hsl.val[2] <= 0.5) ? (delta / (max + min))+perturb : (delta / (2 - max - min));

			float hue;

			if (r == max)
			{
				hue = ((g - b) / 6) / delta;
			}
			else if (g == max)
			{
				hue = (1.0f / 3) + ((b - r) / 6) / delta;
			}
			else
			{
				hue = (2.0f / 3) + ((r - g) / 6) / delta;
			}

			if (hue < 0)
				hue += 1;
			if (hue > 1)
				hue -= 1;

			hsl.val[0] = (int)(hue * 360);
		}
		/*****change HSL  to RGB****/
		cv::Vec3b newcol;

		if (hsl.val[1] == 0)
		{
			newcol.val[0] = newcol.val[1] = newcol.val[2] = (unsigned char)(hsl.val[2] * 255);
		}
		else
		{
			float v1, v2;
			float hue = (float)hsl.val[0] / 360;

			v2 = (hsl.val[2] < 0.5) ? (hsl.val[2] * (1 + hsl.val[1])) : ((hsl.val[2] + hsl.val[1]) - (hsl.val[2] * hsl.val[1]));
			v1 = 2 * hsl.val[2] - v2;

			newcol.val[0] = (unsigned char)(255 * HueToRGB(v1, v2, hue + (1.0f / 3)));
			newcol.val[1] = (unsigned char)(255 * HueToRGB(v1, v2, hue));
			newcol.val[2] = (unsigned char)(255 * HueToRGB(v1, v2, hue - (1.0f / 3)));
		}
	
		return newcol;
	}

	void Image::calcResiduals(cv::Mat &f, cv::Mat &colorRes, cv::Mat &pres, cv::Mat &nres, cv::Mat &signedres) {
		int thresh = 3;

		cv::Mat fg = cv::Mat::zeros(f.rows, f.cols, CV_8UC1);

		cv::Mat pn(f.rows, f.cols, CV_8UC1);
		cv::Mat gr(f.rows, f.cols, CV_8UC1);

		cvtColor(f, fg, cv::COLOR_BGR2GRAY);

		for (int r = 0; r < fg.rows; r++)
			for (int c = 0; c < fg.cols; c++) {
				pn.at<uchar>(r, c) = 128;
				gr.at<uchar>(r, c) = 128;
			}

		for (int r = 0; r < f.rows; r++) {
			for (int c = 0; c < f.cols; c++) {
				colorRes.at<cv::Vec3b>(r, c) = img.at<cv::Vec3b>(r, c) - f.at<cv::Vec3b>(r, c);
			}
		}


		for (int r = 0; r < fg.rows; r++) {
			for (int c = 0; c < fg.cols; c++) {
				signedres.at<float>(r, c) = gray.at<uchar>(r, c) - fg.at<uchar>(r, c); // I changed it for a test to abs values
				if (gray.at<uchar>(r, c) - fg.at<uchar>(r, c) == 0)
				{
					allPix[c + r*img.cols]->white = false;
					allPix[c + r*img.cols]->black = false;
					allPix[c + r*img.cols]->nut = true;
				}

				if (gray.at<uchar>(r, c) - fg.at<uchar>(r, c) > 0)
				{
					allPix[c + r*img.cols]->white = true; allPix[c + r*img.cols]->black = false; allPix[c + r*img.cols]->nut = false;

					gr.at<uchar>(r, c) = 128 + abs(gray.at<uchar>(r, c) - fg.at<uchar>(r, c));

					if (128 + abs(gray.at<uchar>(r, c) - fg.at<uchar>(r, c)) > 255) gr.at<uchar>(r, c) = 255;

					pres.at<uchar>(r, c) = abs(gray.at<uchar>(r, c) - fg.at<uchar>(r, c));
					if (abs(gray.at<uchar>(r, c) - fg.at<uchar>(r, c)) > thresh)
					{
						pn.at<uchar>(r, c) = 255;
					}
					else {
						pn.at<uchar>(r, c) = 128;
					}
				}
				else
				{
					pres.at<uchar>(r, c) = 0;
				}
			}
		}
		for (int r = 0; r < fg.rows; r++) {
			for (int c = 0; c < fg.cols; c++) {
				if (gray.at<uchar>(r, c) - fg.at<uchar>(r, c) < 0)
				{
					allPix[c + r*img.cols]->white = false; allPix[c + r*img.cols]->black = true; allPix[c + r*img.cols]->nut = false;

					gr.at<uchar>(r, c) = 128 - abs(gray.at<uchar>(r, c) - fg.at<uchar>(r, c));
					if (128 - abs(gray.at<uchar>(r, c) - fg.at<uchar>(r, c)) < 0) gr.at<uchar>(r, c) = 0;
					nres.at<uchar>(r, c) = abs(gray.at<uchar>(r, c) - fg.at<uchar>(r, c));
					if (abs(gray.at<uchar>(r, c) - fg.at<uchar>(r, c)) > thresh) {
						pn.at<uchar>(r, c) = 0;					
					}
					else {
						pn.at<uchar>(r, c) = 128;				
					}
				}
				else {
					nres.at<uchar>(r, c) = 0;
				}
			}
		}

		cv::imwrite(path + "pos_neg_gray" + std::to_string(mask) + ".png", gr);
		cv::imwrite(path + "colorRes" + std::to_string(mask) + ".png", colorRes);

		gr_res = gr.clone();

		cv::Mat close = nres.clone();
		cv::Mat element = getStructuringElement(2 /*MORPH_RECT 0*/, cv::Size(3, 3), cv::Point(0, 0));   //cv::MORPH_RECT = 0,   cv::MORPH_ELLIPSE = 2  
		morphologyEx(nres, close, 3 /*MORPH_CLOSE*/, element); //closing
		//cv::imwrite(path + "neg_close_" + mask + ".png", close);

	}

	void Image::getConnectedComponents(int val, void*userdata)
	{
		Image *image = (Image*)userdata;
		cv::Mat bn = cv::Mat::zeros(image->img.size(), CV_8UC1); 
		cv::Mat binary = cv::Mat::zeros(image->img.size(), CV_8UC1);	
		cv::Mat blobdraw = cv::Mat::zeros(image->img.size(), CV_8UC3);
		cv::Mat blobNeg = cv::Mat::zeros(image->img.size(), CV_8UC3);
		cv::Mat theRest = cv::Mat::zeros(image->img.size(), CV_8UC1);
		cv::Mat Rest = cv::Mat::zeros(image->img.size(), CV_8UC3);

		for (int r = 0; r < image->img.rows; r++) {
			for (int c = 0; c < image->img.cols; c++) {

				if (image->pres.at<uchar>(r, c) > val) {
					binary.at<uchar>(r, c) = 1;
				}
				if (image->nres.at<uchar>(r, c) > val) {
					bn.at<uchar>(r, c) = 1;
				}
				if (image->pres.at<uchar>(r, c) <= val && image->pres.at<uchar>(r, c) != 0)
				{
					theRest.at<uchar>(r, c) = 1;
				}
				if (image->nres.at<uchar>(r, c) <= val && image->nres.at<uchar>(r, c) != 0)
				{
					theRest.at<uchar>(r, c) = 1;
				}
				if (image->nres.at<uchar>(r, c) == 0 && image->pres.at<uchar>(r, c) == 0)
				{
					theRest.at<uchar>(r, c) = 1;
				}
			}
		}

		/**************find connected components in positive and negative residuals******************/
		image->FindBlobs(binary, image->blobsP, blobdraw, 160, true);
		image->FindBlobs(bn, image->blobsN, blobNeg, 160, false);
		image->FindBlobs(theRest, image->blobTheRest, Rest, 160, NULL);

		cv::imwrite(path + "blobs"+std::to_string(image->mask)+".png", blobdraw);
		cv::imwrite(path + "blobsN"+std::to_string(image->mask)+".png", blobNeg);
		cv::imwrite(path + "blobs_Rest" + std::to_string(image->mask) + ".png", Rest);

	}

	void Image::SeedGeneration(int mask, Image *image, const char *filename)
	{
		int ct = 0;
		for (int i = 0; i < image->blobsP.size(); i++)
		{

			generateSeedsWithMaxvalue(image->blobsP[i]);
			generateMedianLocSeeds(image->blobsP[i]);
			ct++;
		}
		std::cout << "count of PCC: " << ct << std::endl;
		for (int i = 0; i < image->blobsN.size(); i++)
		{
			generateSeedsWithMaxvalue(image->blobsN[i]); // negative blobs
			generateMedianLocSeeds(image->blobsN[i]);
			ct++;
		}
		std::cout << "count of PCC + NCC: " << ct << std::endl;
	}

	void Image::generateSeedsWithMaxvalue(std::vector<pixel*>& pixels)
	{
		float max_gray = -500; int pid;
		for (int i = 0; i < pixels.size();i++) {
			int res = pixels[i]->res;
			if (abs(res) > max_gray) {
				max_gray = res;
				pid = pixels[i]->getId();
			}
		}
		seeds2.push_back(allPix[pid]);
	}

	void Image::generateMedianLocSeeds(std::vector<pixel*>& pixels)
	{
		std::vector<int> xx; std::vector<int> yy;
		for (int i = 0; i < pixels.size(); i++) {
			// calc median of pixels in this cc;
			xx.push_back(pixels[i]->getPos().x);
			yy.push_back(pixels[i]->getPos().y);
			
		}
		int id = CalcMHWScore(xx) + CalcMHWScore(yy) * img.cols;
		seeds3.push_back(allPix[id]);
		std::sort(seeds3.begin(), seeds3.end(), [](const pixel* a, const pixel* b) {
			return a->res < b->res;// sort small to large
		});

	}

	void Image::generateSeeds(std::vector<pixel*> colorBin[100], bool flag) {

		std::vector<int> neighbours[100]; // keep neighbours for each colorBin id (maximjm  8neighbours)

		for (int i = 0; i < 100; i++) {
			if (i - 1 % 10 >= 0 && i % 10 != 0)
				neighbours[i].push_back(i - 1);
			if (i + 1 % 10 != 0)
				neighbours[i].push_back(i + 1);
			if (i - 10 % 10 >= 0 && i > 9)
				neighbours[i].push_back(i - 10);
			if (i < 90)
				neighbours[i].push_back(i + 10);
			if (i> 10 && i % 10 != 0)
				neighbours[i].push_back(i - 11);
			if (i % 10 != 9)
				neighbours[i].push_back(i + 11);
			if (i - 9 > 0 && i % 10 != 9)
				neighbours[i].push_back(i - 9);
			if (i % 10 != 0 && i < 90)
				neighbours[i].push_back(i + 9);
		}
		std::vector<std::vector<pixel*>> maxima;
		for (int i = 0; i < 100; i++) {
			if (colorBin[i].size() > 0) {
				int count = 0;
		//		std::cout << " i: " << i << std::endl;
				for (int k = 0; k < neighbours[i].size(); k++) {
					if (colorBin[i].size() > colorBin[neighbours[i][k]].size())
						count++;
				}
				if (count == neighbours[i].size())
					maxima.push_back(colorBin[i]);
			}
		}
		//std::cout << " num of local maxima: " << maxima.size() << std::endl;
		for (int i = 0; i < maxima.size(); i++) {
			// calc median of pixels in this bin;
			cv::Point2f mean = cv::Point2f(0, 0);
			std::vector<int> xx; std::vector<int> yy;
			for (int j = 0; j < maxima[i].size(); j++) {
				xx.push_back(maxima[i][j]->getPos().x);
				yy.push_back(maxima[i][j]->getPos().y);
			}
			int id = CalcMHWScore(xx) + CalcMHWScore(yy) * img.cols;
		//	std::cout << "median: " << cv::Point2f(CalcMHWScore(xx), CalcMHWScore(yy)) << " , pix id:" << id << std::endl;
			//if (seeds.size() > 1)
			//{
			//	int ccc = 0;
			//	for (int i = 0; i < seeds.size();i++) {
			//		if (sqrt(seeds[i]->distTo(allPix[id])) < 4)
			//			ccc++;
			//	}
			//	if(ccc == 0)
			//		seeds.push_back(allPix[id]);
			//}else
			seeds.push_back(allPix[id]);

		}

	}

	void Image::calculate_cdf(std::vector<int> hist, std::vector<float> &normalized_cdf) {
		normalized_cdf.resize(hist.size());
		int sum = 0;
		std::vector<float> cdf(hist.size());
		for (int i = 0; i < hist.size(); ++i)
		{
			sum += hist[i];
			cdf[i] = sum;
			//lut[i] = sum * MAX_INTENSITY / pixelCount;
		}
		float max_cdf = *max_element(cdf.begin(), cdf.end());
		float min_cdf = *min_element(cdf.begin(), cdf.end());

		for (int i = 0; i < hist.size(); ++i)
		{
			normalized_cdf[i] = round((cdf[i] - min_cdf) / (img.rows * img.cols - 1) * (255 - 1));
		}
	}

	void Image::HistogramMatching(std::vector<float> grays,  std::vector<int> &lookUp_table) {
		lookUp_table.resize(grays.size());
		//grays hist		
		sort(grays.begin(), grays.end());
		std::vector<int> rArr(grays.size());
		int c = 0;
		for (int i = 0; i < grays.size(); i++) {
			rArr[c] = (int)grays[i];
			//std::cout << rArr[c] << std::endl;
			c++;
		}
		int m = *max_element(rArr.begin(), rArr.end()); //find max value of data points
		int frsize;
		if (m > c - 1)
			frsize = m + 1;
		else frsize = c - 1;
		std::vector<int> rfreq(frsize); //declare frequency array with an appropriate size
		for (int i = 0; i < frsize; i++) //initialize frequency array
			rfreq[i] = 0;
		//compute frequencies
		for (int i = 0; i < c - 1; i++)
			rfreq[rArr[i]]++;

		//calculate_cdf(reg_hist)
		std::vector<float> normalized_rcdf;
		calculate_cdf(rfreq, normalized_rcdf);

		// lookup table
		
		for (int i= 0; i < grays.size();i++) {
			lookUp_table[i] = normalized_rcdf[i];
		}
		//for (int i = 0; i < edgeweights.size();i++) {
		//	std::cout << edgeweights[i] << "  -->  "<< lookUp_table[edgeweights[i]] << std::endl;
		//}

	}

	void Image::FindBlobs(const cv::Mat &binary, std::vector < std::vector<pixel*> > &blobs, cv::Mat &blobdraw, int Mask, bool pos)
	{
		blobs.clear();
		cv::Mat draw = cv::Mat::zeros(img.rows, img.cols, CV_8UC3);

		// Fill the label_image with the blobs
		// 0  - background
		// 1  - unlabelled foreground
		// 2+ - labelled foreground

		cv::Mat label_image;
		binary.convertTo(label_image, CV_32SC1);

		int label_count = 2; // sstump2s at 2 because 0,1 are used already

		for (int y = 0; y < label_image.rows; y++) {
			int *row = (int*)label_image.ptr(y);
			for (int x = 0; x < label_image.cols; x++) {
				if (row[x] != 1) {
					continue;
				}

				cv::Rect rect;
				cv::floodFill(label_image, cv::Point(x, y), label_count, &rect, 0, 0, 4);

				std::vector <pixel*> blob;

				for (int i = rect.y; i < (rect.y + rect.height); i++) {
					int *row2 = (int*)label_image.ptr(i);
					for (int j = rect.x; j < (rect.x + rect.width); j++) {
						if (row2[j] != label_count) {
							continue;
						}
						blob.push_back(allPix[j + i * img.cols]);
					}
				}
				if (blob.size() > 10) {
					blobs.push_back(blob);
					label_count++;
				}
			}
		}
		//======================draw blobs ============
		if (pos != NULL){
			for (int i = 0; i < blobs.size(); i++) {
				cv::Vec3b c = cv::Vec3b((rand() & 255), (rand() & 200), (rand() & 255));// rand() % 100 + 100;
				for (int j = 0; j < blobs[i].size(); j++) {
					int rr = blobs[i][j]->pos.y; int cc = blobs[i][j]->pos.x;
					allPix[cc + rr * img.cols]->setPosNegBlobId(i + 1, pos);
					draw.at<cv::Vec3b>(rr, cc) = c;
				}
			}
		}
		if(pos == NULL)
		{
			for (int i = 0; i < blobs.size(); i++) {
				cv::Vec3b c = cv::Vec3b((rand() & 255), (rand() & 200), (rand() & 255));
				for (int j = 0; j < blobs[i].size(); j++) {
					int rr = blobs[i][j]->pos.y; int cc = blobs[i][j]->pos.x;
					draw.at<cv::Vec3b>(rr, cc) = c;
				}
			}
		}
		blobdraw = draw.clone();
	}

	void Image::FindBlobs1(const cv::Mat &binary, std::vector < std::vector<pixel*> > &blobs, cv::Mat &blobdraw)
	{
		blobs.clear();
		cv::Mat draw = cv::Mat::zeros(img.rows, img.cols, CV_8UC3);

		// Fill the label_image with the blobs
		// 0  - background
		// 1  - unlabelled foreground
		// 2+ - labelled foreground

		cv::Mat label_image;
		binary.convertTo(label_image, CV_32SC1);
		
		int label_count = 2; // sstump2s at 2 because 0,1 are used already

		for (int y = 0; y < label_image.rows; y++) {
			int *row = (int*)label_image.ptr(y);
			for (int x = 0; x < label_image.cols; x++) {
				if (row[x] != 1) {
					continue;
				}
				cv::Rect rect;
				cv::floodFill(label_image, cv::Point(x, y), label_count, &rect, 0, 0, 4);

				std::vector <pixel*> blob;

				for (int i = rect.y; i < (rect.y + rect.height); i++) {
					int *row2 = (int*)label_image.ptr(i);
					for (int j = rect.x; j < (rect.x + rect.width); j++) {
						if (row2[j] != label_count) {
							continue;
						}
						blob.push_back(allPix[j + i * img.cols]);
					}
				}

				blobs.push_back(blob);
				label_count++;		

			}
		}
		//======================draw blobs ============

		for (int i = 0; i < blobs.size(); i++) {
			cv::Vec3b c = cv::Vec3b((rand() & 255), (rand() & 200), (rand() & 255));
			for (int j = 0; j < blobs[i].size(); j++) {
				int rr = blobs[i][j]->pos.y; int cc = blobs[i][j]->pos.x;
				draw.at<cv::Vec3b>(rr, cc) = c;
			}
		}
		blobdraw = draw.clone();
	}

	void Image::calcSLIC(int mask, int numsp, cv::Mat &draw) {

		/**********delete slicBlobs**********/
		for (int i = 0; i < slicBlobs.size(); i++) {
			slicBlobs[i].clear();
		}
		slicBlobs.clear();
		/***********************************/
		IplImage* input = cvCloneImage(&(IplImage)img); //image2; // image2 filtered image
		IplImage* lab_image = cvCloneImage(input);

		// apply pre-processing functions
		cvCvtColor(input, lab_image, CV_BGR2Lab);

		/* Yield the number of superpixels and weight-factors from the user. */
		int w = input->width, h = input->height;

		int nr_superpixels = numsp /*atoi(argv[2])*/;
		int nc = 30;// 30 /*atoi(argv[3])*/; smaller nc more details
		if (numsp == 10) nc = 100;
		double step = sqrt((w * h) / (double)nr_superpixels);

		/* Perform the SLIC superpixel algorithm. */
		Slic slic;
		slic.generate_superpixels(lab_image, step, nc);
		slic.create_connectivity(lab_image);

		/* Display the contours and show the result. */
	//	slic.display_contours(input, CV_RGB(255, 0, 0));

		numSlic = 0;
		slicIDs = slic.PixelSLICid();
		assignSLICid(slicIDs);
		numSlic = slic.centersSize();

	//	cv::Mat draw = cv::Mat::zeros(img.size(), CV_8UC3);
		//std::vector<std::vector<double>> cent = slic.getCenters();

		for (int k = 0; k < numSlic;k++) {
			std::vector<pixel*> blob;
			for (int i = 0; i < allPix.size(); i++) {
				if (allPix[i]->slicId == k) {
					blob.push_back(allPix[i]);

					//int idcol = (int)cent[k][3] + (int)cent[k][4] * img.cols;
					//int r = allPix[i]->getPos().y; int c = allPix[i]->getPos().x;
					//draw.at<cv::Vec3b>(r, c) = allPix[idcol]->getColor();
				}
			}
			slicBlobs.push_back(blob);
		}
		//cv::imwrite(path + "SLIC_blobs" + std::to_string(mask) + "_" + std::to_string(numsp) + "sp.png", draw);
		/********draw ave of slic blobs*************************/
		//if (numsp == 10) {
			cv::Mat draw1 = cv::Mat::zeros(img.size(), CV_8UC3);

			std::vector<std::vector<double>> cent = slic.getCenters();
			for (int i = 0; i < slicBlobs.size(); i++) {
				cv::Vec3b col(0,0,0);
				int sum1 = 0,sum2=0,sum3=0;
				for (int j = 0; j < slicBlobs[i].size(); j++) {
					sum1= sum1 + slicBlobs[i][j]->getColor().val[0];
					sum2 = sum2 + slicBlobs[i][j]->getColor().val[1];
					sum3 = sum3 + slicBlobs[i][j]->getColor().val[2];
					//std::cout << "sum1: "<< sum1 << "sum2:" <<sum2 << "sum3: " << sum3 << std::endl;
				}
				col.val[0] = sum1 / slicBlobs[i].size();
				col.val[1] = sum2 / slicBlobs[i].size();
				col.val[2] = sum3 / slicBlobs[i].size();
				//std::cout << "ave col" << col << std::endl;

				for (int j = 0; j < slicBlobs[i].size(); j++) {
				
					int r = slicBlobs[i][j]->getPos().y; int c = slicBlobs[i][j]->getPos().x;
					draw.at<cv::Vec3b>(r, c) = col;
				}
			}
		//	cv::imwrite(path + "SLIC_blobs" + std::to_string(mask) + "_" + std::to_string(numsp) + ".png", draw1);
		//}
		/********************************/
	
		cent.clear();
		for (int i = 0; i < slicIDs.size();i++) {
			slicIDs[i].clear();
		}
		cv::Mat m = cv::cvarrToMat(input);
		/*******black and white *******/
		cv::Mat draw2 = cv::Mat::zeros(img.size(), CV_8UC3);
		for (int i = 0; i < slicBlobs.size();++i) {
			cv::Vec3b col(0,0,0);
			cv::Vec3i sum(0, 0, 0);
			for (int j = 0; j < slicBlobs[i].size(); ++j) {
				sum += slicBlobs[i][j]->getColor();
			}
			col = cv::Vec3b(sum.val[0]/slicBlobs[i].size(), sum.val[1] / slicBlobs[i].size(), sum.val[2] / slicBlobs[i].size());
		
			if ( (col.val[0]+ col.val[1]+ col.val[2])/3.0> 120) {
				for (int j = 0; j < slicBlobs[i].size(); ++j) {
					int r = slicBlobs[i][j]->getPos().y; int c = slicBlobs[i][j]->getPos().x;
					draw2.at<cv::Vec3b>(r, c) = cv::Vec3b(255,255,255);
				}
			}
			else {
				for (int j = 0; j < slicBlobs[i].size(); ++j) {
					int r = slicBlobs[i][j]->getPos().y; int c = slicBlobs[i][j]->getPos().x;
					draw2.at<cv::Vec3b>(r, c) = cv::Vec3b(0, 0, 0);
				}
			}
		}
		cv::imwrite(path + "SLIC_BW_" + std::to_string(numsp) + ".png", draw2);
		//cv::imwrite(path + "SLIC" + std::to_string(mask) + "_"+std::to_string(numsp)+ "sp.png", m);
	}

	void Image::FindCC(const cv::Mat &visible, int &maxLabel)
	{
		maxLabel = 0;
		cv::Mat draw = cv::Mat::zeros(img.rows, img.cols, CV_8UC3);

		for (int v = 0; v < img.rows * img.cols; v++) {
			allPix[v]->visited = false;
			allPix[v]->CCid = -1;
		}

		cv::Mat  grayVisible = cv::Mat::zeros(visible.size(), CV_8UC1);
		cvtColor(visible, grayVisible, cv::COLOR_BGR2GRAY);

		cv::Mat label_image;
		grayVisible.convertTo(label_image, CV_32SC1);
		int label_count = -1; // sstump2s at 2 because 0,1 are used already

		for (int y = 0; y < label_image.rows; y++) {
			int *row = (int*)label_image.ptr(y);
			for (int x = 0; x < label_image.cols; x++) {
				//if (row[x] != 1) {
				//	continue;
				//}
				if (allPix[x + y * img.cols]->visited == true || allPix[x + y * img.cols]->CCid != -1) {
					continue;
				}
				cv::Rect rect;
				cv::floodFill(label_image, cv::Point(x, y), label_count++, &rect, 8);

				for (int i = rect.y; i < (rect.y + rect.height); i++) {
					int *row2 = (int*)label_image.ptr(i);
					for (int j = rect.x; j < (rect.x + rect.width); j++) {
						if (j >= 0 && j < img.cols && i >= 0 && i < img.rows) {
							if (row2[j] != label_image.at<int>(y, x)) {
								continue;
							}
							allPix[j + i * img.cols]->visited = true;
							allPix[j + i * img.cols]->setCCId(label_count);
						}
					}
				}
				//label_count++;
			}
		}	
		maxLabel = label_count;

		regions.clear(); regions.resize(maxLabel + 1);
		for (int i = 0; i < allPix.size(); i++) {
			regions[allPix[i]->CCid].push_back(allPix[i]);
		}
		//======================draw blobs ============

		//for (int i = 0; i < maxLabel; i++) {
		//	cv::Vec3b c = cv::Vec3b((rand() & 255), (rand() & 200), (rand() & 255));
		//	for (int j = 0; j < allPix.size(); j++) {
		//		if (allPix[j]->CCid == i) {
		//			int rr = allPix[j]->pos.y; int cc = allPix[j]->pos.x;
		//			draw.at<cv::Vec3b>(rr, cc) = c;
		//		}
		//	}
		//}
		//imwrite(path + "labeled_mask1.png", draw);
	}

	void Image::floodFill(cv::Mat &image, int sr, int sc, int newcolor, cv::Vec3b oldColor)
	{
		if (sr < 0 || sr >= image.rows || sc < 0 || sc >= image.cols || image.at<cv::Vec3b>(sr,sc) != oldColor )
			return;
		allPix[sc + sr * img.cols]->setCCId(newcolor);
		floodFill(image, sr - 1, sc, newcolor, oldColor);
		floodFill(image, sr + 1, sc, newcolor, oldColor);
		floodFill(image, sr, sc - 1, newcolor, oldColor);
		floodFill(image, sr, sc + 1, newcolor, oldColor);
	}

	void Image::Flood_fill(cv::Mat &image, int sr, int sc, int newcolor) // newcolor = new id
	{	
		if (/*image.at<int>(sr, sc) != newcolor*/ allPix[sc + sr * img.cols]->CCid != newcolor) // 只有 oldColor 是应当被替换的颜色，其他颜色不应被替换
			floodFill(image, sr, sc, newcolor, image.at<cv::Vec3b>(sr, sc));
	}

	void Image::separateConnectedRegions(cv::Mat &visible, int &maxLabel)
	{
		cv::Mat label,stats,centroids;
		cv::Mat  grayVisible = cv::Mat::zeros(visible.size(),CV_8UC1);

		cvtColor(visible, grayVisible, cv::COLOR_BGR2GRAY);

		cv::Mat draw = cv::Mat::zeros(img.size(), CV_8UC3);

		int k = 0;
		int count = 0;
		while(k<256) {
			cv::Mat  mask = cv::Mat::zeros(visible.size(), CV_8UC1);
			for (int r = 0; r < img.rows; r++) {
				for (int c = 0; c < img.cols; c++) {
					if(grayVisible.at<uchar>(r,c) == k )
						mask.at<uchar>(r, c) = 1;
				}
			}
			std::vector<std::vector<pixel*>> holes;
			FindBlobs1(mask, holes, draw);
			for (int i = 0; i < holes.size(); i++) {
				for (int j = 0; j < holes[i].size(); j++) {
					int c = holes[i][j]->getPos().x; int r = holes[i][j]->getPos().y;
					allPix[c + r * img.cols]->visited = true;
					allPix[c + r * img.cols]->setCCId(count);
					allPix[c + r * img.cols]->setRegionId(count);
				}
				count++;	
			}
			k++;
		}
		maxLabel= count-1;
	}

	void Image::DFSUtil(int v, /*bool visited[],*/int label)
	{
		// Mark the current node as visited and print it 
		//visited[v] = true;
		allPix[v]->visited = true;
		//cout << v << " ";
		allPix[v]->setCCId(label);
		// Recur for all the vertices 
		// adjacent to this vertex 
		//std::vector<pixel>::iterator i;
		//for (i = allPix[v]->neighbours.begin(); i != allPix[v]->neighbours.end(); ++i)
		//{
		//	if (allPix[v]->getRegId() == (i)->getRegId())
		//	{
		//		int id = (i)->getId();
		//		if (/*!visited[id]*/(i)->visited == false)
		//			DFSUtil(id,/* visited,*/ label);
		//	}
		//}
		for (int j = 0; j < allPix[v]->neighbours.size(); j++) {
			//int regid = allPix[v]->edges[j]->getTheOtherSide(allPix[v])->getRegId();
		//	if (allPix[v]->edges[j] != 0 && allPix[v]->getRegId() == allPix[v]->edges[j]->getTheOtherSide(allPix[v])->getRegId())
			if ( allPix[v]->getRegId() == allPix[v]->neighbours[j]->getRegId())
			{
				int id = allPix[v]->neighbours[j]->getId();
				if (allPix[id]->visited == false/* visited[id] == false*/)
					DFSUtil(id, /*visited,*/ label);
			}
		}
	}
	
	void Image::separateConnectedRegions2(int &label)
	{
		// Mark all the vertices as not visited 
		bool *visited = new bool[img.rows * img.cols];
		for (int v = 0; v < img.rows * img.cols; v++) {
			//visited[v] = false;
			allPix[v]->visited = false;
			allPix[v]->CCid = -1;
		}
		label = -1;
		for (int v = 0; v < img.rows * img.cols; v++)
		{
			if (allPix[v]->CCid != -1)// since ccid is -1 it has not been visited but it's also a new cluster -> increase the label
				continue;
			else 
				label++;
			if (allPix[v]->visited == false  /*visited[v] == false*/)
			{
				// print all reachable vertices from v
				DFSUtil(v, /*visited ,*/ label);
			//	cout << "\n";
			}
	
		}

	}

	void Image::assignSLICid(std::vector<std::vector<int>> SLICids) {
		for (int c = 0; c < img.cols;c++) {
			for (int r = 0; r < img.rows;r++) {
				int id = c + r * img.cols;
				allPix[id]->slicId = SLICids[c][r];
			}
		}
	}

	void Image::calcKmean_reColor(std::vector<pixel*> &regionpix/*, std::vector<cv::Point> contours*/ ,cv::Mat &draw3, int K) {
	
	//	if (regionpix.size() > contours.size() + 10 ) {
			// exclude contours from region

			std::vector<pixel*> cents;
			std::vector<cv::Point3f> dataCol;
			std::vector<pixel*> pixx(regionpix.size());
			pixx.clear();

			for (int i = 0; i < regionpix.size(); i++) {
				bool flag = false;
				cv::Point3f col;
				col.x = regionpix[i]->getColor().val[0];
				col.y = regionpix[i]->getColor().val[1];
				col.z = regionpix[i]->getColor().val[2];
				//for (int j = 0; j < contours.size(); j++) {
				//	if (contours[j].x == regionpix[i]->getPos().x && contours[j].y == regionpix[i]->getPos().y)
				//		flag = true; // it's a contour
				//}
			//	if (flag == false) {
					dataCol.push_back(col);
					pixx.push_back(regionpix[i]);
			//	}
			}
			cv::Mat labels, centers;
			std::vector<std::vector<pixel*>> clustered_colors(K);
			std::vector<std::vector<cv::Point3f>> clustered_centers(K);

			if (dataCol.size() > 10)
			{
				kmeans(dataCol, K, labels, cv::TermCriteria(), 1, cv::KMEANS_PP_CENTERS, centers);
				for (int i = 0; i < labels.rows; i++)
				{
					int idx = labels.at<int>(i);
					cv::Point3f original_color = dataCol[i];
					cv::Point3f clustered_center;
					clustered_center.x = centers.at<float>(idx, 0);
					clustered_center.y = centers.at<float>(idx, 1);
					clustered_center.z = centers.at<float>(idx, 2);
					//	cerr << i << " " << idx << " " << original_color << " " << clustered_center << endl;
					clustered_colors[idx].push_back(pixx[i]);
					clustered_centers[idx].push_back(clustered_center);
				}

				if (K > 2) sort(clustered_colors.begin(), clustered_colors.end());

				//	std::cout << "size of the clusters: " << clustered_centers[0].size() << ", " << clustered_centers[1].size() << std::endl;

					// return just center pixel of first two largest bins
				for (int i = 0; i < 2; i++) {
					cv::Point2f sum(0, 0);
					std::vector<pixel*> closeCol;
					int count = 0;
					cv::Vec3b color1 = cv::Vec3b(
						(uchar)floor(clustered_centers[i][0].x),
						(uchar)floor(clustered_centers[i][0].y),
						(uchar)floor(clustered_centers[i][0].z));
					//std::cout << "color:  " << color1 << std::endl;
					for (int j = 0; j < clustered_colors[i].size(); j++) {
						//sum = sum + clustered_colors[i][j]->getPos();
						if (clustered_colors[i][j]->sqrtColorDistTo(color1) < 3.0f) {
							//closeCol.push_back(clustered_colors[i][j]);
							sum = sum + clustered_colors[i][j]->getPos();
							count++;
						}
					}
					if (count == 0)
					{
						sum = cv::Point2f(0, 0);
						for (int j = 0; j < clustered_colors[i].size(); j++) {
							sum = sum + clustered_colors[i][j]->getPos();
							count++;
						}
					}
					//		std::cout << "count: " << count << std::endl;
					if (count > 0) {
						cv::Point2f c = cv::Point2f((sum.x / count), (sum.y / count));
						int id = floor(c.x) + floor(c.y) * img.cols;
						cents.push_back(allPix[id]);
					}
				}
				if (cents.size() == 2) {
					float var0 = 0, var1 = 0;
					for (int i = 0; i < clustered_colors[0].size(); i++) {
						var0 += clustered_colors[0][i]->distTo(cents[0]);
					}
					for (int i = 0; i < clustered_colors[1].size(); i++) {
						var1 += clustered_colors[1][i]->distTo(cents[1]);
					}
					var0 /= (clustered_colors[0].size() - 1);
					var1 /= (clustered_colors[1].size() - 1);
					float var = (var0 + var1) / 2;

					for (int i = 0; i < regionpix.size(); i++) {
						float W0, W1; cv::Point3i C0, C1;

						//float d1 = sqrt(regionpix[i]->distTo(cents[0]));
						//float d2 = sqrt(regionpix[i]->distTo(cents[1]));

						W0 = exp(-(regionpix[i]->distTo(cents[0])) / (2 * var0));
						W1 = exp(-(regionpix[i]->distTo(cents[1])) / (2 * var1));
						float dis = sqrt(cents[1]->distTo(cents[0]));

						C0 = clustered_centers[0][0];
						C1 = clustered_centers[1][0];

						cv::Vec3b p = cv::Vec3b(int((W0 * C0.x + W1 * C1.x) / (W0 + W1)),
							int((W0 * C0.y + W1 * C1.y) / (W0 + W1)),
							int((W0 * C0.z + W1 * C1.z) / (W0 + W1)));
						//if ( d1 == dis/(2.0f) || d1 == dis / (3.0f) ||d1 == dis / (4.0f) || d1 == dis / (5.0f))
						//	p = cv::Vec3b(0, 0, 0);
						//if (d2 == dis / (2.0f) || d2 == dis / (3.0f) || d2 == dis / (4.0f) || d2 == dis / (5.0f))
						//	p = cv::Vec3b(0, 100, 0);

						regionpix[i]->setreColor(p);
						int c = regionpix[i]->getPos().x;	int r = regionpix[i]->getPos().y;
						draw3.at<cv::Vec3b>(r, c) = p;// regionpix[i]->reColor;
					}
				}
				else {
					cv::Point3f sumc = cv::Point3f(0.0, 0.0, 0.0);
					for (int i = 0; i < regionpix.size(); i++) {
						int c = regionpix[i]->getPos().x;	int r = regionpix[i]->getPos().y;
						sumc += cv::Point3f(regionpix[i]->getColor().val[0], regionpix[i]->getColor().val[1], regionpix[i]->getColor().val[2]);
					}
					for (int i = 0; i < regionpix.size(); i++) {
						int c = regionpix[i]->getPos().x;	int r = regionpix[i]->getPos().y;
						draw3.at<cv::Vec3b>(r, c) = cv::Vec3b((int)(sumc.x/ regionpix.size()),(int)(sumc.y/regionpix.size()),
							(int)(sumc.z/ regionpix.size()));// regionpix[i]->reColor;
					}
				}
			}
			else{
					cv::Point3f sumc = cv::Point3f(0.0,0.0,0.0);
					for (int i = 0; i < regionpix.size(); i++) {
						int c = regionpix[i]->getPos().x;	int r = regionpix[i]->getPos().y;
						sumc += cv::Point3f(regionpix[i]->getColor().val[0], regionpix[i]->getColor().val[1], regionpix[i]->getColor().val[2]);
					}
					for (int i = 0; i < regionpix.size(); i++) {
						int c = regionpix[i]->getPos().x;	int r = regionpix[i]->getPos().y;
						draw3.at<cv::Vec3b>(r, c) = cv::Vec3b((int)(sumc.x / regionpix.size()), (int)(sumc.y / regionpix.size()),
							(int)(sumc.z / regionpix.size()));// regionpix[i]->reColor;
					}
				}
		//}

		cv::imwrite(path + "recoloredReg.png", draw3);

		//cv::circle(draw3,cents[0]->getPos(), 1, cvScalar(0,0,0), 1, 8);
		//cv::circle(draw3, cents[1]->getPos(), 1, cvScalar(0, 0, 0), 1, 8);
	//	cv::line(draw3, cents[0]->getPos(), cents[1]->getPos(), cvScalar(0, 0, 0), 1, 8);
	}

	//*******Bhattacharyya measurement********/
	void Image::calcMahalanobisAB(std::vector<cv::Vec3f> distribution, std::vector<cv::Vec3f> initDistribution,std::vector<float> &covarElement, double &min_dist)
	{
		double m_distAB, m_distBA;
		cv::Mat covarMat(3, 3, CV_32FC1); cv::Mat covarMatB(3, 3, CV_32FC1);;
		cv::Mat mean(1, 3, CV_32FC1); cv::Mat meanB(1, 3, CV_32FC1);
		/*********************upcoming distribtion**************************/
		cv::Mat samplesB(distribution.size(), 3, CV_32FC1);
		for (int i = 0; i < distribution.size(); i++) {
			samplesB.at<float>(i, 0) = distribution[i].val[0];
			samplesB.at<float>(i, 1) = distribution[i].val[1];
			samplesB.at<float>(i, 2) = distribution[i].val[2];
		}
		cv::calcCovarMatrix(samplesB, covarMatB, meanB, CV_COVAR_NORMAL | CV_COVAR_ROWS, 5);
		 /*********************initial distribtion**************************/
		cv::Mat samples(initDistribution.size(), 3, CV_32FC1);
		for (int i = 0; i < initDistribution.size(); i++) {
			samples.at<float>(i, 0) = initDistribution[i].val[0];
			samples.at<float>(i, 1) = initDistribution[i].val[1];
			samples.at<float>(i, 2) = initDistribution[i].val[2];
		}
		cv::calcCovarMatrix(samples, covarMat, mean, CV_COVAR_NORMAL | CV_COVAR_ROWS, 5);

		m_distBA = cv::Mahalanobis(meanB, mean, covarMat);
		m_distAB = cv::Mahalanobis(mean, meanB, covarMatB);
		covarElement.push_back(covarMat.at<float>(0, 0));	
		covarElement.push_back(covarMat.at<float>(1, 1));
		covarElement.push_back(covarMat.at<float>(2, 2));
		covarElement.push_back(covarMatB.at<float>(0, 0));
		covarElement.push_back(covarMatB.at<float>(1, 1));
		covarElement.push_back(covarMatB.at<float>(2, 2));
		//min_dist = (mean-meanB)/sqrt(covarMat.at<float>(0, 0) + covarMatB.at<float>(0, 0));// min(m_distBA, m_distAB);
		//Bhattacharya distance
		cv::Mat trans;// = cv::Mat(3, 1, mean.type());
		cv::transpose(mean - meanB, trans);
		cv::Mat tomultiply = cv::Mat(3,1,mean.type());
		cv::Mat inv_mat = cv::Mat(3, 3, mean.type()); inv_mat = ((covarMat + covarMatB) / 2.0).inv(1);
		tomultiply = (inv_mat) * (trans); // 3x3 3x1
		//std::cout << "mean size: "<< mean.rows << mean.cols << ", trans rows: "<< trans.rows << ", tomultiply size: " << tomultiply.size() << std::endl;
		//std::cout << tomultiply << std::endl;
		float BD; //Bhattacharya distance
		BD = 0.125f * cv::determinant((mean - meanB) * tomultiply) + 0.5 *log((cv::determinant((covarMat + covarMatB) / 2))/(sqrt(cv::determinant(covarMat)*cv::determinant(covarMatB))));
		//Bhattacharya Coeff
		min_dist = 1.0f / exp(BD); //Bhattacharya coeff
		double stddev0 = covarMat.at<float>(0, 0);
		double stddev1 = covarMat.at<float>(1, 1);
		double stddev2 = covarMat.at<float>(2, 2);
		//std::cout << "std0 : " << stddev0 << ", std1 : " << stddev1 << " ,std2 : " << stddev2 << std::endl;
	}

	void Image::calcMahalanobis(cv::Vec3f newPoint, std::vector<cv::Vec3f> distribution, double &m_dist)
	{
		cv::Mat newP(1,3,CV_32FC1); 
		newP.at<float>(0, 0) = newPoint.val[0];
		newP.at<float>(0, 1) = newPoint.val[1];
		newP.at<float>(0, 2) = newPoint.val[2];

		cv::Mat covarMat;
		cv::Mat mean(1,3, CV_32FC1);

		cv::Mat samples(distribution.size(), 3, CV_32FC1);
		for (int i = 0; i < distribution.size();i++) {
			samples.at<float>(i, 0) = distribution[i].val[0];
			samples.at<float>(i, 1) = distribution[i].val[1];
			samples.at<float>(i, 2) = distribution[i].val[2];
		}
		cv::calcCovarMatrix(samples, covarMat, mean,CV_COVAR_NORMAL | CV_COVAR_ROWS,5);

		m_dist = cv::Mahalanobis(newP, mean, covarMat);
		double stddev0 = covarMat.at<float>(0, 0);
		double stddev1 = covarMat.at<float>(1, 1);
		double stddev2 = covarMat.at<float>(2, 2);
		std::cout << "std0 : " << stddev0 << ", std1 : " << stddev1 << " ,std2 : " << stddev2 << std::endl;
	}

	void Image::colorMode(std::vector<pixel*>& region, cv::Mat &screen, cv::Vec3b &mode)
	{
		std::vector<int>count(64, 0);
		int maxIndex = -1, maxCount = 0;
		int index;
		for (int i = 0; i < region.size();i++) {
			int r = region[i]->getPos().y;
			int c = region[i]->getPos().x;

			index = floor(gray.at<uchar>(r, c)/4.0f);
			count[index]++;
			if (count[index] > maxCount)
			{			
				maxCount = count[index];
				maxIndex = index;
			}
		}
		cv::Vec3d ave;
		std::vector<cv::Vec3b> modes;

		for(int i = 0; i < region.size(); i++) {
			int r = region[i]->getPos().y;
			int c = region[i]->getPos().x;
			if (maxIndex == floor(gray.at<uchar>(r, c) / 4.0f)) {
				modes.push_back( img.at<cv::Vec3b>(r, c));
			}
		}
		if (modes.size() > 0 && maxIndex !=-1) {
			for (int i = 0; i < modes.size(); i++) {
				ave.val[0] += modes[i].val[0]; ave.val[1] += modes[i].val[1]; ave.val[2] += modes[i].val[2];
			}
			mode = cv::Vec3b(ave.val[0] / modes.size(), ave.val[1] / modes.size(), ave.val[2] / modes.size());

			//	std::cout << "mode color: "<< mode << std::endl;
			for (int i = 0; i < region.size(); i++) {
				int r = region[i]->getPos().y;
				int c = region[i]->getPos().x;
				screen.at<cv::Vec3b>(r, c) = mode;
			}
		}
		else {
			for (int i = 0; i < region.size(); i++) {
				int r = region[i]->getPos().y;
				int c = region[i]->getPos().x;
				screen.at<cv::Vec3b>(r, c) = img.at<cv::Vec3b>(r, c);
			}
		}

	}
	
	void Image::colorAverage(std::vector<pixel*> &region, cv::Mat &screen, cv::Vec3b &ave)
	{
		cv::Vec3d sum;
		if (region.size() > 0) {
			for (int p = 0; p < region.size(); p++) {

				int x = region[p]->getPos().x; int y = region[p]->getPos().y;
				sum.val[0] += region[p]->getColor().val[0];
				sum.val[1] += region[p]->getColor().val[1];
				sum.val[2] += region[p]->getColor().val[2];
			}
			ave = cv::Vec3b((float)sum.val[0] / region.size(), (float)sum.val[1] / region.size(), (float)sum.val[2] / region.size());
			//-----------average color-----------------
			for (int p = 0; p < region.size(); p++) {
				int x = region[p]->getPos().x; 
				int y = region[p]->getPos().y;
				screen.at<cv::Vec3b>(y, x) = ave;
			}
		}
	}


	inline cv::Vec3b Image::returnAverageColor(std::vector<pixel*>& region)
	{
		cv::Vec3b ave(0, 0, 0);
		cv::Vec3d sum;
		//--- for color averaging----/*
		for (int p = 0; p < region.size(); p++) {

			int x = region[p]->getPos().x; int y = region[p]->getPos().y;
			sum.val[0] += region[p]->getColor().val[0];
			sum.val[1] += region[p]->getColor().val[1];
			sum.val[2] += region[p]->getColor().val[2];
		}
		ave = cv::Vec3b((float)sum.val[0] / region.size(), (float)sum.val[1] / region.size(), (float)sum.val[2] / region.size());

		return ave;
	}
		
	void Image::calcCovarMat(std::vector<cv::Vec3f> distribution,cv::Mat &covarMat, cv::Mat &mean) {
		//mean = cv::Mat::zeros(1, 3, CV_32FC1);
		//covarMat = cv::Mat::zeros(3, 3, CV_32FC1);
		cv::Mat samples(distribution.size(), 3, CV_32FC1);
		for (int i = 0; i < distribution.size(); i++) {
			samples.at<float>(i, 0) = distribution[i].val[0];
			samples.at<float>(i, 1) = distribution[i].val[1];
			samples.at<float>(i, 2) = distribution[i].val[2];
		}
		cv::calcCovarMatrix(samples, covarMat, mean, CV_COVAR_NORMAL | CV_COVAR_ROWS, 5);
	}

	void Image::calcKmean(int K, cv::Mat &gy) {

		int n = img.rows * img.cols;
		cv::Mat data = img.reshape(1, n);
		data.convertTo(data, CV_32F);
		std::vector<int> labels;
		cv::Mat colors = cv::Mat(img.size(), CV_32F);
		kmeans(data, K, labels, cv::TermCriteria(), 1, cv::KMEANS_PP_CENTERS, colors);
		for (int i = 0; i < n; ++i)
		{
			data.at<float>(i, 0) = colors.at<float>(labels[i], 0);
			data.at<float>(i, 1) = colors.at<float>(labels[i], 1);
			data.at<float>(i, 2) = colors.at<float>(labels[i], 2);
		}
		cv::Mat reduced = data.reshape(3, img.rows);
		reduced.convertTo(reduced, CV_8U);
		cvtColor(reduced, gy, cv::COLOR_BGR2GRAY);
		cv::imwrite(path + "Reduced" + std::to_string(K) + ".png", gy);
	}
	
	/*********************generate regions *********************************/	
	void Image::createRegions(int mask,  std::string str) { // find extremas and put them in pq
		struct CompareRes
		{
			bool operator()(pixel *n1, pixel *n2) const
			{
				return abs(n1->res) < abs(n2->res); // < larger to smaller values - > small to large priority
			}
		};
		struct CompareRegion
		{
			bool operator()(pixel *n1, pixel *n2) const
			{
				return abs(n1->regionsize) < abs(n2->regionsize); // < larger to smaller values - > small to large priority
			}
		};
		int rcount[256];
		int gcount[256];
		int bcount[256];
		int icount[256];
		cv::Mat screen = cv::Mat::zeros(img.size(), CV_8UC3);
		cv::Mat screen3 = cv::Mat::zeros(img.size(), CV_8UC3);
		cv::Mat gy = cv::Mat::zeros(img.size(), CV_8UC1); 

		/************************************/
		for (int r = 0; r < screen.rows; r++) {
			for (int c = 0; c < screen.cols; c++) {
				screen.at<cv::Vec3b>(r, c) = cv::Vec3b(128, 128, 128);
			}
		}

		/*********k-mean quantization****************/
		const int K = 64;
		calcKmean(K, gy);
		/***************************************************************************/
		std::priority_queue< pixel*, std::vector< pixel*>, CompareRes > pixel_queue, pixel_queue2;
		std::priority_queue< pixel*, std::vector< pixel*>, CompareRegion > pixels;
		std::priority_queue<std::vector< pixel*>, std::vector< std::vector< pixel*>>, CompareRegions > Regions_Q;
		std::pair<float, pixel*> pairs;
		pixel_queue.empty();
		/**********************seeds = body color-- seeds2 = maxvalue -- seeds3 = median(sorted) ***********************/
		
		for (int i = 0; i <seeds3.size();i++) {
			//if(i > 0.5 * seeds3.size())
			//pixel_queue.push(seeds3[(seeds3.size()-1)%rand()]);
			pixel_queue.push(seeds3[i]);
		}

		/****************** Seed in grid structure**************/
		//std::vector<int> num_s ; // seed id
		//for (int r= 10; r < img.rows;r+=30) {
		//	for (int c = 10; c < img.cols;c+= 30) {
		//		int id = c + r * img.cols; 
		//		pixel_queue.push(allPix[id]);
		//		num_s.push_back(id);				
		//	}
		//}
		/********************************************************/
		std::cout << "size of the seeds: " << seeds3.size()/* seeds.size()*/  << std::endl;
		std::cout << "size of the pixel_queue: " << pixel_queue.size() << std::endl;
	
		int regsize;
		int mycost;
		int whereami;

		double totred = 0;
		double totgreen = 0;
		double totblue = 0;
		int s, t, curcen;
		int  p, q;

		rgbvector neighcol, origcol, herecol;

		int TARGREGSIZE = 10000;//accum; stop when the pixel number is accum

		// initialization pass here:
		for (int i = 0; i <allPix.size(); i++)
		{
			allPix[i]->partcost = allPix.size();// XD*YD;
			allPix[i]->dist = FLT_MAX;
			allPix[i]->isVisited = false;
			allPix[i]->parent = NULL;
			allPix[i]->ancesId = NULL;
			allPix[i]->ancesParent = NULL;
			allPix[i]->setNewColor(128);
			allPix[i]->sn = NULL; // blob id
			allPix[i]->parentsize = 0;
			allPix[i]->sumColdist = 0.0f;
			allPix[i]->lock = false;
		}
		cv::Mat colorimage = img.clone();
		std::cout << "Starting the semi" << std::endl;
		int regId =  0;
		int step = 0; int st = 0;
		/*************************************************/

		while (!pixel_queue.empty()) // large to small priority
		{
			pixel* up = pixel_queue.top();
			pixel_queue.pop();
			pixels.push(up);
		}

		/*******base layer ******/
		cv::Mat base = cv::Mat::zeros(screen.size(), screen.type());
		cv::Mat mask_background = cv::Mat::zeros(screen.size(), CV_8UC1);
		cv::Mat mask_draw = cv::Mat::zeros(screen.size(), CV_8UC3);

		//-----for backgroud using slic image------------
		cv::Mat slicimage = cv::Mat::zeros(img.size(), img.type());
		calcSLIC(160, 580,slicimage);
		for (int i = 0; i< allPix.size(); i++) {
			allPix[i]->slicId = -1;
		}
		slicBlobs.clear();

		cv::imwrite(path + "slicimage.png", slicimage);
		screen = base.clone();
		
		clock_t t1, t2, t3, t4,t5,t6;
		t1 = clock();
		/***************************/
		float thresh = FLT_MAX; // initial error : used for thresholding

		/*******************************************************/
		while (!pixels.empty()) // large to small priority
		{
			cv::Mat growth = cv::Mat::zeros(base.size(), CV_8UC3);
		
			for (int i = 0; i < allPix.size(); i++)
			{
				allPix[i]->lock = false;
			}
			st++;

			pixel *up = pixels.top();
			pixels.pop();

			up->parentsize = 0;

			TARGREGSIZE = /*min(up->regionsize, 1000);*/  up->regionsize;
			//TARGREGSIZE = 5000;
			// a new building process centred on current pixel
			curcen = up->id;

			allPix[curcen]->dist = 0;
			allPix[curcen]->setRegionId( regId);

			origcol.red = allPix[curcen]->getColor().val[0];
			origcol.green = allPix[curcen]->getColor().val[1];
			origcol.blue = allPix[curcen]->getColor().val[2];

			//clearhist();
			for (int i = 0; i < 255; i++)
			{
				rcount[i] = 0;
				gcount[i] = 0;
				bcount[i] = 0;
				icount[i] = 0;
			}

			std::priority_queue<std::pair<float, pixel*>, std::vector<std::pair<float, pixel*>>, CompareDist > region_queue;

			region_queue.empty();
			up->dist = 0;
			region_queue.push(std::make_pair(up->dist, up));

			regsize = 0;

			totred = 0;
			totblue = 0;
			totgreen = 0;
			float dist = 0;
			/*********************************************************/
			int totwhite = 0;
			int totblack = 0;
			bool grow = true;
			up->grow = true;

			std::vector<pixel*> region;
			up->parent = up; // seed pparent  = seed

			bool flag = false;
			thresh = 1000000000000.0; // initial error : used for stoppingg growth

			while (!region_queue.empty()   /* regsize < TARGREGSIZE*/)
			{
				std::pair<float, pixel*> pair = region_queue.top();
				region_queue.pop();

				int idu = pair.second->id;

				dist = (pair.second)->dist;// pair.first;

				allPix[idu]->isVisited = true;
				if (allPix[idu]->lock == true) continue;

			//	if (dist == allPix[idu]->dist){
				
					allPix[idu]->setRegionId(regId);
					// this node is good, pix[s] can become part of current region

					p = allPix[idu]->getPos().y;
					q = allPix[idu]->getPos().x;

					herecol.red = allPix[idu]->getColor().val[0];
					herecol.green = allPix[idu]->getColor().val[1];
					herecol.blue = allPix[idu]->getColor().val[2];

					totred += allPix[idu]->getColor().val[0];
					totgreen += allPix[idu]->getColor().val[1];
					totblue += allPix[idu]->getColor().val[2];

					rcount[herecol.red]++;
					gcount[herecol.green]++;
					bcount[herecol.blue]++;

					// increment region size
					regsize++;

					// run through all neighbours and add to heap:
					std::vector<int> neighbours; // index of neighbours
					neighbours.empty();
					for (unsigned int i = 0; i < 8; i++) {
						if (allPix[idu]->edges[i] != 0)
							neighbours.push_back(i);
					}
				
					for (int j = 0; j < neighbours.size(); j++)
					{
						//pixel *v = allPix[idu]->edges[neighbours[j]]->getTheOtherSide(allPix[idu]);

						int nbId = allPix[idu]->edges[neighbours[j]]->getTheOtherSide(allPix[idu])->id;

						// cgrf from center to pixel(v) [changed up->distTo(v) to u->distTo(v)]			
						float nWeight = 1.0+  sqrt(allPix[idu]->distTo(allPix[nbId])) + up->colorDistTo(allPix[nbId]) /*up->angular_colorDistTo(allPix[nbId])*/;
						float alt = dist + nWeight;
						if (alt < allPix[nbId]->dist /* && (allPix[nbId]->getRegId() < allPix[idu]->getRegId())*/ )
						{
							allPix[nbId]->dist = alt;
							float pred_dist = 0.0;
							/**************set the flag when 10 parents has been passed **************/
							if (allPix[idu]->parentsize == 10 && flag == false  ) {
								thresh = allPix[nbId]->dist;// allPix[idu]->returnParent(10)->dist; // initial error
								if (thresh == 0) std::cout << "thresh is 0 ..."  << std::endl;
								//std::cout << "parent size: " << allPix[idu]->parentsize << std::endl;
								flag = true;		
							}
							if (allPix[idu]->parentsize > 10 ) {
								pred_dist = allPix[idu]->returnParent(10)->dist;
							}
							//if the difference of current pixel and the predessesor (10th parent) is over a threshold stops to grow
							if (abs(allPix[nbId]->dist - pred_dist) > abs( 100 * thresh )  ) { 
								continue;
							}
							else  {
								allPix[nbId]->parent = allPix[idu];
								allPix[nbId]->parentsize = allPix[idu]->parentsize + 1;

								allPix[nbId]->regionIds.push_back(allPix[nbId]->getRegId());
								region_queue.push(std::make_pair(alt, allPix[nbId]));
							}
									
						}
						else {
							// do nothing
						}
					}
			//	}
			

			} // end "big enough region"
		//	std::cout<< "reigon grow stopped -----" << std::endl;
	
			// now, region completed, compute average and put into output image
			// mean

			for (int k = 0; k < allPix.size(); k++)
			{
				if (allPix[k]->getRegId() == regId) {
					int pp = allPix[k]->getPos().y;
					int qq = allPix[k]->getPos().x;
					allPix[k]->finalized = true; // for drawing
					region.push_back(allPix[k]);
				}
			}
			Regions_Q.push(region);
			regions.push_back(region);

			region.clear();
			/************************************************************************************/
			//// progress meter to console
			if (regId % 5 == 0 || regId % 2 == 0 || regId % 3 == 0)
			{
				//std::cout << "region id:   " << regId << std::endl;
				//cv::imwrite(path + "inProgress_" + std::to_string(regId) + ".png", screen);
			//	cv::imwrite(path + "growth2/growth" + std::to_string(regId) + ".png", growth);
			}
			region_queue.empty();
			regId++;
		}

		t2 = clock();
		float diff((float)t2 - (float)t1);
		std::cout << diff / CLOCKS_PER_SEC << "  seconds to generate original regions " << std::endl;
		/************sort vector of regions based on size******************************************/
		std::sort(regions.begin(), regions.end(), [](const std::vector<pixel*>& a, const std::vector<pixel*>& b) {
			return a.size() > b.size();
		});

		cv::Mat simply = cv::Mat::zeros(img.size(), img.type()); //averagecol.copyTo(simply);
		simply = img.clone();
		slicimage.copyTo(simply);


		std::vector<int> sortedRegId;
		/*****draw regions after sort by size**********************************/
		cv::Vec3b mode;
		slicimage.copyTo(screen3);
		for (int ii = 0; ii < regions.size(); ii++) {
			if (regions[ii].size() > 0) {
				sortedRegId.push_back(regions[ii][0]->getRegId()); // sorted based on region size . larger to smaller																//---------mode color----------------
				colorMode(regions[ii], screen3, mode);
				if (ii % 150 == 0) // progress meter to console
				{
					std::cout << " number of regions:   " << regions[ii].size() << std::endl;
				}
			}
		}

		cv::imwrite(path + "inOrder_v1.png", screen3);
		std::cout << " number of regions before find CC:   " << regions.size() << std::endl;

		/***new region id assigned to visible pixels*****/
		int maxleb = 0;
		t3 = clock();
		FindCC(screen3, maxleb);// (works well )

		t4 = clock();
		float diff2((float)t4 - (float)t3);
		std::cout << diff2 / CLOCKS_PER_SEC << "  seconds to find all connected components (flattening) " << std::endl;

		std::cout << "number of over segmented regions before merging:" << regions.size() << std::endl;
		/**********merging small regions- first round********/
		mergingSmallRegions();
		
		maxleb = 0; cv::Vec3b ave;
		screen3 = cv::Mat::zeros(img.size(), img.type());
		for (int i = 0; i < regions.size(); i++) {
			if (regions.size() != 0)
				colorAverage(regions[i], screen3, ave);
		}
		
		int count = 0;
		for (int i = 0; i < regions.size(); i++)
			if (regions[i].size() > 0)
				count++;
		std::cout << "number of over segmented regions after merging single pixels:" << count << std::endl;
		/**********merging small regions- second round********/
		int sz = 0;
		for (sz = 1; sz < 10; sz++) {
			mergingSmallRegionsOfSize(sz); // based on neighbors color (can have more than one neighbor)	
		}
		screen3 = cv::Mat::zeros(img.size(), img.type());

		for (int i = 0; i < regions.size(); i++) {
			if (regions.size() != 0)
				colorAverage(regions[i], screen3, ave);
		}
		/******find Connected Components after second round of merging**********/
		FindCC(screen3, maxleb);// (works well )
		/*******************/
		count = 0;
		for (int i = 0; i < regions.size(); i++)
			if (regions[i].size() != 0)
				count++;
		std::cout << "number of over segmented regions after merging regions of size smaller than "<<sz<<":" << count << std::endl;


		colorRegions();

		//**** clear ******
		AllRegions.clear();
		//regions.clear();
		//allPix.clear();

		fprintf(stderr, "All done!\n");

	}
	
	/*************make graph on regions and recolor them *************/

	void Image::makeGraphonRegions(cv::Mat& img, std::string str) {
		clock_t t1, t2;
		t1 = clock();
		Region_Graph* Graph = new  Region_Graph(path, img, allPix, regions);
		t2 = clock();
		float difft((float)t2 - (float)t1);
		std::cout << difft / CLOCKS_PER_SEC << "  seconds to create graph and neighbours " << std::endl;
		std::cout << " number of edges of the graph: " << Graph->edgeweights.size() / 2 << std::endl;
		cv::Mat draw = cv::Mat::zeros(img.size(), CV_8UC3);
		
		/***********************************/
		// assign random colors
		cv::Mat randCol = cv::Mat::zeros(img.size(), CV_8UC3);
		Graph->RecoloringRandomColors(randCol);
		cv::imwrite(path + "RandColorRegs.png", randCol);
		/***********************************/
		std::vector<cv::Vec3b> palette;
		std::vector<std::vector<cv::Vec3b>> palettes;
		//  wpap_palette
		//palette.push_back(cv::Vec3b(255, 255, 255));   //1 white
		palette.push_back(cv::Vec3b(244, 229, 192));  //  2 light blue
		palette.push_back(cv::Vec3b(224, 181, 12));   //  3 blue
		palette.push_back(cv::Vec3b(142, 76, 5));   //  4 dark blue
		palette.push_back(cv::Vec3b(29, 16, 14));   //  5  darkest blue. almost black
		palette.push_back(cv::Vec3b(26, 236, 242));   //  6 yellow
		palette.push_back(cv::Vec3b(60, 205, 161));  //  7 light green
		palette.push_back(cv::Vec3b(69, 145, 27));   //  8  green
		palette.push_back(cv::Vec3b(74, 28, 235));   //  9  red-pink
		palette.push_back(cv::Vec3b(67, 18, 122));   //  10  purple
		palette.push_back(cv::Vec3b(64, 14, 34));   //  11  dark purple
		palettes.push_back(palette); palette.clear();
		//calmwater
		//http://colrd.com/palette/23719/
		palette.push_back(cv::Vec3b(230, 183, 191));  //rgb(191,183,230)   
		palette.push_back(cv::Vec3b(245, 163, 145));  //rgb(145,163,245)  
		palette.push_back(cv::Vec3b(193, 134, 125));  //rgb(125,134,193)
		palette.push_back(cv::Vec3b(116, 56, 64));  //rgb(64,56,116)   
		palette.push_back(cv::Vec3b(78, 28, 38));  //rgb(38,28,78)  
		palette.push_back(cv::Vec3b(55, 9, 31));  //rgb(31,9,55)
		palette.push_back(cv::Vec3b(49, 67, 87));  //rgb(87,67,49)   
		palette.push_back(cv::Vec3b(33, 145, 157));  //rgb(157,145,33)   
		palette.push_back(cv::Vec3b(89, 153, 164));  //rgb(164,153,89)
		palette.push_back(cv::Vec3b(126, 179, 182));  //rgb(182,179,126)
		palettes.push_back(palette); palette.clear();
		//  colorful palette
		palette.push_back(cv::Vec3b(207, 116, 158));   //  1 light purple
		palette.push_back(cv::Vec3b(132, 52, 87));  //  2  darker purple
		palette.push_back(cv::Vec3b(84, 32, 41));   //  3  darkest purple
		palette.push_back(cv::Vec3b(112, 54, 162));   //  4 dark pink 
		palette.push_back(cv::Vec3b(117, 68, 237));   //  5  light pink
		palette.push_back(cv::Vec3b(142, 210, 254));   //  6  yellow
		palette.push_back(cv::Vec3b(116, 170, 255));  //  7  yello-orange-light
		palette.push_back(cv::Vec3b(89, 127, 250));   //  8  orange
		palette.push_back(cv::Vec3b(75, 86, 227));   //  9  darker orange
		palette.push_back(cv::Vec3b(75, 38, 169));   //  10  almost red
		palettes.push_back(palette); palette.clear();

		std::vector<int> indexes;
		int mid_reg;
		std::vector<int> compression_params;
		compression_params.push_back(IMWRITE_PNG_COMPRESSION);
		compression_params.push_back(9);
		std::vector < std::pair<int, int>> values_p;
		std::vector<cv::Vec3b> newpalet = palettes[0];
		std::string palettename[3] = { "wpap","calmwater","colorful"};

		for (int p = 0; p < palettes.size(); p++) {
			/*********start from a region in the middle of the image*******/
			int pid = (img.cols / 3) + (img.rows / 3) * img.cols;
			mid_reg = Graph->RegionPixels[pid]->getRegId();
			/************************/
			draw = cv::Mat::zeros(img.size(), CV_8UC3);
			Graph->RecoloringByBottleneck(draw, /*newpalet*/ palettes[p], palettename[p], indexes, values_p, mid_reg); // increasekey
			imwrite(path + "RecoloredBN_" + palettename[p] + "_Palette_Euclid_lut_" + str + "_r-mid-tableW-midw_1+b.png", draw, compression_params);

			/**********Recoloring without path planning***************/
			//draw = cv::Mat::zeros(img.size(), CV_8UC3); values_p.empty();
			//Graph->RecoloringPathFree(draw, palettes[p], values_p);
			//imwrite(path + "RecoloredPathFree_" + palettename[p] + "_Palette_Euclid_" + str + ".png", draw);

			/**********Random recoloring***************/
			//draw = cv::Mat::zeros(img.size(), CV_8UC3); values_p.empty();
			//Graph->Recoloring_random(draw, palettes[p], values_p);
			//imwrite(path + "RecoloredRandom_" + palettename[p] + "_Palette_Euclid_" + str + ".png", draw);

			/***********Recoloring without lookup table**************/
			//draw = cv::Mat::zeros(img.size(), CV_8UC3); indexes.clear(); values_p.clear();
			//Graph->RecoloringByBottleneck_noLUT(draw, palettes[p], palettename[p], indexes, values_p, mid_reg); //no Histogram Matching
			//imwrite(path + "RecoloredBN_" + palettename[p] + "_Palette_Euclid_Nolut_" + str + "_r-mid-OrgW-midw.png", draw);
		}

		delete Graph;
	}

	void Image::colorRegions(){
		cv::Vec3b ave,mode;
		cv::Mat draw1 = cv::Mat::zeros(img.size(), img.type());
		cv::Mat draw2 = cv::Mat::zeros(img.size(), img.type());
		for (int i = 0; i < regions.size(); i++)
		{
			if (regions[i].size() > 0) {
				colorAverage(regions[i], draw1, ave);
				colorMode(regions[i], draw2, mode);
			}
		}
		cv::imwrite(path + "AverageColors.png", draw1);
		cv::imwrite(path + "ModeColors.png", draw2);
	}

	void Image::mergingSmallRegionsOfSize(int sz) {

		for (int i = 0; i < regions.size(); i++) {
			bool checkequality = false;
			if (regions[i].size() <= sz && regions[i].size() > 0) {
				std::vector<int> newids;
				for (int j = 0; j < regions[i].size(); j++) {
					int id = regions[i][j]->getId();
					//std::cout << "id: " << id << std::endl;
					for (int n = 0; n < allPix[id]->neighbours.size(); n++) {
						if (allPix[id]->neighbours[n] != NULL) {
							if (allPix[id]->neighbours[n]->CCid != allPix[id]->CCid) {
								newids.push_back(allPix[id]->neighbours[n]->CCid);
							}
						}
					}
				}
				//std::cout << "newids size: " << newids.size() << std::endl;
				checkequality = checkALLIDS(newids); // JUST ONE NEIGHBOR
				if (checkequality) {
					for (int j = 0; j < regions[i].size(); j++) {
						//regions[i][j]->setCCId(newids[0]);
						allPix[regions[i][j]->getId()]->setCCId(newids[0]);
						regions[newids[0]].push_back(allPix[regions[i][j]->getId()]);
					}
					newids.clear();
				}
				else {
					/****remove duplicate neighbour ids****/
					std::vector<int>::iterator ip;
					// Sorting the array 
					std::sort(newids.begin(), newids.end());
					// Using std::unique 
					ip = std::unique(newids.begin(), newids.begin() + newids.size());
					// Resizing the vector so as to remove the undefined terms 
					newids.resize(std::distance(newids.begin(), ip));

					/*******assign same CCid to the neighbor region thathas closest color******/
					float mincol = 10000; int newid;
					for (int n = 0; n < newids.size(); n++) {
						float distcolor = colorDifferance(returnAverageColor(regions[newids[n]]), returnAverageColor(regions[i]));
						if (distcolor < mincol) {
							mincol = distcolor;
							newid = newids[n];
						}
					}
					for (int j = 0; j < regions[i].size(); j++) {
						//regions[i][j]->setCCId(newid);
						allPix[regions[i][j]->getId()]->setCCId(newid);
						regions[newid].push_back(allPix[regions[i][j]->getId()]);
					}
					newids.clear();
				}
				regions[i].clear();
			}

		}
		int maxleb = regions.size();
		regions.clear();
		regions.resize(maxleb);
		for (int i = 0; i < allPix.size(); ++i) {
			regions[allPix[i]->CCid].push_back(allPix[i]);
		}
	}


	void Image::mergingSmallRegions() {
		std::vector<std::vector<pixel*>> smallregions;
		for (int i = 0; i < regions.size(); i++) {
			if (regions[i].size() == 1)
			{
				smallregions.push_back(regions[i]);
			}
		}
		std::random_shuffle(smallregions.begin(), smallregions.end());
		for (int i = 0; i < smallregions.size(); i++) {
			//find neighbours regionid, if more than 80% of them belong to one region assign same id to current small region
			int newid; int regid = smallregions[i][0]->CCid;
			for (int j = 0; j < smallregions[i].size(); j++) {
				int id = smallregions[i][j]->getId();
				for (int k = 0; k < allPix[id]->neighbours.size(); k++) {
					if (allPix[id]->neighbours[k] != NULL) {
						if (allPix[id]->neighbours[k]->CCid != allPix[id]->CCid) {
							newid = allPix[id]->neighbours[k]->CCid;// allPix[id]->getRegId();
						}
					}
				}
			}
			for (int j = 0; j < smallregions[i].size(); j++) {
				smallregions[i][j]->setCCId(newid);
			}
		}
		smallregions.clear();
		int maxleb = regions.size();
		regions.clear();
		regions.resize(maxleb);
		for (int i = 0; i < allPix.size(); ++i) {
			regions[allPix[i]->CCid].push_back(allPix[i]);
		}
	}

	void Image::smoothBoundary(std::vector<std::vector<cv::Point>>& contours, std::vector<pixel*> &region,cv::Mat &draw)
	{
		cv::Mat draw2 = draw.clone();
		cv::Vec3b mode;
		colorMode(region, draw2, mode);
	//	draw2.convertTo(draw, -1, 1.0, 20);

		for (int i = 0; i < region.size(); i++) {
			int y = region[i]->getPos().y;	int x = region[i]->getPos().x;
			for (int c = 0; c < img.channels(); c++) {
				draw.at<cv::Vec3b>(y, x)[c] = saturate_cast<uchar>(1.0*draw2.at<Vec3b>(y, x)[c] + 20);
			}
		}

		int maxindex = 0; int max = -10;
		for (int i = 0; i < contours.size();i++) {
			if (contours[i].size() > max) {
				maxindex = i;
				max = contours[i].size();
			}
		}

		int N = contours[maxindex].size();
		if (N >= 5) {

			// store complex number in a two channel image
			// y axis is the imaginary axis (channel 1)
			cv::Mat in(1, N, CV_32FC2);
			//for each point i in the contour
			for (int i = 0; i < N; i++) {
				in.at<cv::Vec2f>(0, i)[0] = float(contours[maxindex][i].x);
				in.at<cv::Vec2f>(0, i)[1] = float(contours[maxindex][i].y);
			}

			// forward transform
			// trans will store the real and imaginary components of the 
			// transformed contour in channels 0 and 1

			cv::Mat trans, ctrs;
			dft(in, trans, cv::DFT_COMPLEX_OUTPUT);

			int ncoeff = 5; //7 
			int start = (ncoeff / 2) + 1;
			int end = N - ((ncoeff / 2)) - 1;
			for (int i = start; i < end; i++) {
				trans.at<cv::Vec2f>(0, i)[0] = 0;
				trans.at<cv::Vec2f>(0, i)[1] = 0;
			}

			idft(trans, ctrs, cv::DFT_SCALE);

			std::vector<cv::Point> contour; std::vector<std::vector<cv::Point>> allcontour;
		//	if (N > 6) {
				
				for (int i = 0; i < N; i++) {
					if (ctrs.at<cv::Vec2f>(0, i).val[1] >= 0 && ctrs.at<cv::Vec2f>(0, i).val[1] < img.rows
						&&  ctrs.at<cv::Vec2f>(0, i).val[0] >= 0 && ctrs.at<cv::Vec2f>(0, i).val[0] < img.cols) {
						//draw.at<cv::Vec3b>(ctrs.at<cv::Vec2f>(0, i).val[1], ctrs.at<cv::Vec2f>(0, i).val[0]) = cv::Vec3b(255, 255, 255);
						contour.push_back(cv::Point(ctrs.at<cv::Vec2f>(0, i).val[0], ctrs.at<cv::Vec2f>(0, i).val[1]));
					}
					else if (ctrs.at<cv::Vec2f>(0, i).val[1] < 0) {
						//draw.at<cv::Vec3b>(0, ctrs.at<cv::Vec2f>(0, i).val[0]) = cv::Vec3b(255, 255, 255);
						contour.push_back(cv::Point(ctrs.at<cv::Vec2f>(0, i).val[0], 0));
					}
					else if (ctrs.at<cv::Vec2f>(0, i).val[1] >= img.rows) {
						//draw.at<cv::Vec3b>(img.rows - 1, ctrs.at<cv::Vec2f>(0, i).val[0]) = cv::Vec3b(255, 255, 255);
						contour.push_back(cv::Point(ctrs.at<cv::Vec2f>(0, i).val[0], img.rows - 1));
					}
					else if (ctrs.at<cv::Vec2f>(0, i).val[0] < 0) {
						//draw.at<cv::Vec3b>(ctrs.at<cv::Vec2f>(0, i).val[1], 0) = cv::Vec3b(255, 255, 255);
						contour.push_back(cv::Point(0, ctrs.at<cv::Vec2f>(0, i).val[1]));
					}
					else if (ctrs.at<cv::Vec2f>(0, i).val[0] >= img.cols) {
						//draw.at<cv::Vec3b>(ctrs.at<cv::Vec2f>(0, i).val[1], img.cols - 1) = cv::Vec3b(255, 255, 255);
						contour.push_back(cv::Point(img.cols - 1, ctrs.at<cv::Vec2f>(0, i).val[1]));
					}
				}
				if (contour.size() > 0) {
					allcontour.push_back(contour);					
					cv::Scalar color = cv::Scalar(cv::saturate_cast<uchar>(mode.val[0]-5), cv::saturate_cast<uchar>(mode.val[1]-5),
						cv::saturate_cast<uchar>(mode.val[2]-5));
					drawContours(draw, allcontour,-1, color, CV_FILLED);
					drawContours(draw, allcontour, -1, cvScalar(0));// , 0.5, 8);
				}
				//drawContours(draw, contour, cv::Scalar(region[0]->getColor()), 0.5, 8);// , hierarchy, 0);
		//	}

		}
	}
	

} // namespace abstraction
