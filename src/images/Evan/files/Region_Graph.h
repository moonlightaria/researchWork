#pragma once
#ifndef REGION_GRAPH_H
#define REGION_GRAPH_H

#include "Region.h"

template <class myType>
static myType CalcMHWScore(std::vector<myType> scores) {
	int size = scores.size();
	//if (size == 0)
	//{
	//	return 0;  // Undefined, really.
	//}
	//else
	//{
		sort(scores.begin(), scores.end());
		if (size % 2 == 0)
		{
			return (myType)(scores[size / 2 - 1] + scores[size / 2]) / 2;
		}
		else
		{
			return (myType)scores[size / 2];
		}
	//}
}
template <class myType>
static myType CalcAverage(std::vector<myType> scores) {
	myType avesize=0.0;
	for (int i = 0; i < scores.size(); i++) {
		avesize += scores[i];
	}
	return (avesize / scores.size());
}

static int counter = 0;

namespace abstraction {

	static int findmax(int vals[], int n) {
		int max = vals[0];
		for (int i = 0; i<n; i++)
			if (max<vals[i]) max = vals[i];
		return max;
	}
	static void calculate_cdf(float MaxValue,std::vector<int> hist, std::vector<float> &normalized_cdf ) {

		int maxHist = *max_element(hist.begin(), hist.end());
		normalized_cdf.resize(hist.size());
		float sum = 0.0f;
		float sH = 0.0;// some of all elements
		for (int i = 0; i < hist.size(); ++i)
		{
			sH += hist[i];
		}
		std::vector<float> cdf(hist.size());
		for (int i = 0; i < hist.size(); i++)
		{
			sum += (float)(hist[i]/sH);
			//std::cout << "cdf[i]: " << sum << std::endl;
			cdf[i] = sum ;
		}
		float max_cdf = *max_element(cdf.begin(),cdf.end());
		float min_cdf = *min_element(cdf.begin(), cdf.end());

		for (int i = 0; i < hist.size(); i++)
		{
			normalized_cdf[i] = cdf[i];// (max_cdf * (cdf[i] - min_cdf) / (max_cdf - min_cdf));// / max_cdf;	
			//std::cout << "normalized_cdf[i]:  " << normalized_cdf[i] << std::endl;
		}
		std::cout << "normalized_cdf[final]:  " << normalized_cdf[hist.size()-1] << std::endl;
		cdf.clear();
	}
	static std::vector<int> lookupTable(std::vector<float> &normalized_rcdf, std::vector<float> &normalized_tcdf, std::vector<float> &tArr) {
		int size = normalized_rcdf.size();
		std::vector<int> lookup_table(size);// = np.zeros(256);
		for (int i = 0; i < size;i++) {
			lookup_table[i] = 0;
		}

		int lookup_val = 0;
		int j = 0;
		for (int i = 0; i <size;i++) {
			while (normalized_tcdf[j] < normalized_rcdf[i] 
				/*(abs(normalized_tcdf[j + 1] - normalized_rcdf[i]) <= abs(normalized_tcdf[j] - normalized_rcdf[i]))*/ && j < normalized_tcdf.size()) {
				j++;
			}
			lookup_table[i] = j;// (int)floor((normalized_rcdf.size() - 1) * normalized_tcdf[j]);// j;
		}
		return lookup_table;
	}

	static void SortedlookupTable(std::vector<float> rArr, std::vector<float> tArr, std::vector<int> &lookup_table) {
		lookup_table.resize(rArr.size());

		float m = *max_element(rArr.begin(), rArr.end());
		float mm = *max_element(tArr.begin(), tArr.end());
		float min = *min_element(tArr.begin(), tArr.end());
		std::cout << "rArr size : " << rArr.size() << "tArr size : " << tArr.size() << std::endl;
		for (int i = 0; i < rArr.size(); i++) {
			int id = (float(/*rArr[i] / m*/ (i*1.0f)/rArr.size())) * (tArr.size());// (int) floor((float(i / rArr.size())) * tArr.size());
			//std::cout << "id: " << id << std::endl;
			if (id >= tArr.size()) id = tArr.size() - 1;
			if (id < 0) id = 0;
			lookup_table[i] = tArr[id];
		}

	}

	static float colorDifferance(cv::Vec3b col1, cv::Vec3b col2) {
	
		ColorSpace::Rgb Rgb1(col1.val[2],col1.val[1], col1.val[0]);
		ColorSpace::Rgb Rgb2(col2.val[2], col2.val[1], col2.val[0]);
		ColorSpace::Lab Lab1, Lab2;
		double diff = ColorSpace::EuclideanComparison::Compare(&Rgb1, &Rgb2);

		return (float) diff;
	}

	class Region_Graph {

		friend struct RegionEdge;
	public:

		Region_Graph();
		~Region_Graph();
		int numReg;
		std::vector<Region*> graph_regions;
		std::vector<pixel*> RegionPixels;
		std::vector<std::vector<int>> adj;
		std::vector<RegionEdge*> edges[8];
		cv::Mat image;
		std::string path;
		std::vector<float> edgeweights;
		std::pair<int, cv::Vec3b> in_colors[256];

		std::vector<float> treeEdges; // edge weights in the tree
		std::vector<float> normalized_rcdf;
		std::vector<float> normalized_tcdf;

		void make();

		Region_Graph(std::string path,cv::Mat &img, std::vector<pixel*> &allPix, std::vector<std::vector<pixel*>> &_regions );

		Region * returnRegion(int _regid);

		cv::Vec3b colorAverage(std::vector<int> ids);

		float colorDiff(Region * Reg1, Region * Reg2);
		void DrawBarCOLOR(std::vector<int> bin, std::vector<cv::Vec3b> palette, cv::Mat &screen , std::vector<int> &indexes, std::vector < std::pair<int, int>> &values_p);
		cv::Point retunMidPoint(Region *reg);
		void DrawChildren(Region * root, cv::Mat & draw);
		void AssignColors(Region * source, std::vector<cv::Vec3b> palette, cv::Mat &table, std::vector<int> lookUp_table, int binsize, Region *node);
		void recolorWithModeofAssignedColors(cv::Mat &draw, std::vector<cv::Vec3b> palette);
		void DrawTree(Region *root, cv::Mat & abstraction,int count);
		float getMidWeight(std::vector<float> weights, std::vector<int> &lookUp_table);
		void HistogramMatching(std::vector<float> differences, cv::Mat &table, std::vector<int> &lookUp_table, int &binsize);
		bool checkOverlaps(Region* reg1, Region* reg2);
		void addEdge(std::vector<std::vector<int>> &adj, int u, int v);
		void addEdge(Region * Reg, int v);
		void findNeighbours(std::vector<pixel*> &allPix);
		void findNeighbours2(std::vector<pixel*>& allPix);
		void quaitizeImageColorTo(std::vector<cv::Vec3b> palette, cv::Mat &coloredImage);
		void colorRegion(Region *reg, cv::Mat &contourImage);
		void colorRegion(Region * reg, cv::Mat &contourImage, cv::Vec3b palettecolor);
		void drawEdges(cv::Mat &screen);
		std::pair<int, int> ReturnClosestPalette(RegionEdge * edge, cv::Mat & table, std::vector<cv::Vec3b> palette, std::vector<int> lookUp_table);
		int ReturnClosestPaletteIndex(RegionEdge * edge, std::vector<cv::Vec3b> palette, std::vector<int> lookUp_table);
		int ReturnClosestPaletteIndexForLargeRegions(Region * source, std::vector<cv::Vec3b> palette, cv::Mat & table, std::vector<int> lookUp_table, Region * curReg);
		int ReturnClosestPaletteIndex1(Region *source, std::vector<cv::Vec3b> palette, cv::Mat & table, std::vector<int> lookUp_table, int binsize, Region * reg);
		int ReturnClosestPaletteIndex2(RegionEdge * edge, std::vector<cv::Vec3b> palette, cv::Mat &table, std::vector<int> lookUp_table, cv::Vec3b C);
		void ClosestPaletteIndexes(RegionEdge * edge, cv::Mat & table, std::vector<cv::Vec3b> palette, std::vector<int> lookUp_table, std::vector<threeObj> &pairs);
		void RecoloringByColorPalette(cv::Mat & draw);
		void reducePalette(std::vector<int>& indexes, std::vector<cv::Vec3b>& palette, std::vector<cv::Vec3b>& newpalet);
		void BFS(Region * curr, cv::Mat &draw);
		void RecoloringByBottleneck(cv::Mat & screen, std::vector<cv::Vec3b>& palette, std::string palettename, 
			std::vector<int> &indexes, std::vector < std::pair<int, int>> &values_p, int &mid_reg);
		void RecoloringByBottleneck_noLUT(cv::Mat & screen, std::vector<cv::Vec3b>& palette0, std::string palettename, std::vector<int>& indexes, std::vector<std::pair<int, int>>& values_p, int & mid_reg);
		void RecoloringRandomColors(cv::Mat& screen);
		void RecoloringPathFree(cv::Mat & screen, std::vector<cv::Vec3b>& palette, std::vector < std::pair<int, int>> &values_p);
		void Recoloring_random(cv::Mat & screen, std::vector<cv::Vec3b>& palette, std::vector<std::pair<int, int>>& values_p);
		void Histogram(std::vector<float> differences,std::string);
		void printGraph();// V= num regions //A utility function to print the adjacency list representation of graph 
	};
} //namespace abstraction
#endif //