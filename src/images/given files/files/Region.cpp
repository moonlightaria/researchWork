#include "Region.h"

namespace abstraction {

	Region::Region() {
		regId = NULL;
	}
	Region::Region(std::vector<pixel*> regPix) {

		for (int i = 0; i < regPix.size(); i++) {
			regionPix.push_back(regPix[i]);
		}
		setRegId(regPix[0]->CCid);
	}
	Region::~Region() {
		for (int i = 0; i < regionPix.size();i++) {
			delete regionPix[i];
		}
	}

	void Region::setNewColor(cv::Vec3b newCol) {
		newColor = newCol;
	}
	cv::Vec3b Region::getNewColor()
	{
		return newColor;
	}
	int Region::getId() {
		return regId;
	}

	void Region::setRegId(int _rid) {
		regId = _rid;
	}

	void Region::calcBoundingbox() {
		maxx = (*std::max_element(regionPix.begin(), regionPix.end(), x_compare))->getPos().x;
		maxy = (*std::max_element(regionPix.begin(), regionPix.end(), y_compare))->getPos().y;
		minx = (*std::min_element(regionPix.begin(), regionPix.end(), x_compare))->getPos().x;
		miny = (*std::min_element(regionPix.begin(), regionPix.end(), y_compare))->getPos().y;

		std::cout << "regId: " << regId << ", " << minx << "," << maxx << "," << miny << "," << maxy << std::endl;
	}

	void Region::addElements(pixel* p) {
		regionPix.push_back(p);
	}

	bool Region::DoBoxesIntersect(Region *B) {

		return (this->minx < B->minx + abs(B->minx - B->maxx) &&
			this->minx + abs(this->minx - this->maxx) > B->minx &&
			this->miny < B->miny + abs(B->miny - B->maxy) &&
			this->miny + abs(this->miny - this->maxy) > B->miny);
		// collision detected!

		//return (abs(this->minx - B->minx) * 2 < (abs(this->minx -this->maxx) + abs(B->minx - B->maxx))) &&
		//	(abs(this->miny - B->miny) * 2 < (abs(this->miny - this->maxy) + abs(B->miny - B->maxy)));
	}
	void Region::removeDuplicates(std::vector<int> &v) // remove duplicates
	{
		auto end = v.end();
		for (auto it = v.begin(); it != end; ++it) {
			end = std::remove(it + 1, end, *it);
		}
		v.erase(end, v.end());
		for (int i = 0; i < v.size(); i++) {
			neighbours.push_back(v[i]);
		}
	}

	cv::Vec3b Region::calcClosestToMode() {
		// find the mode of the colors in region pixels, 
		//average the colors of fthe mode, then find the closest pixel to the mode color
		int regsize = regionPix.size();
		if (regsize == 0) std::cout << "region size is zero!!!" << std::endl;
		//std::vector<int> histogram(256);
		int histogram[256];
		std::vector<std::vector<pixel*>> pixelss(256);
		std::vector<pixel* > binPixels;
		for (int i = 0; i < 256; i++)
		{
			histogram[i] = 0;
			pixelss[i].clear();
		}

		for (int j = 0; j < regsize; j++) {
			int bin = regionPix[j]->getIntensity();
			histogram[bin]++;
			//pixelss[bin].push_back(regionPix[j]);
		}
		//int maxHist = std::max_element(histogram.begin(), histogram.end()) -histogram.begin();
		int maxHist = -10;
		int maxHistValue = -1;
		for (int i = 0; i < 256; i++)
		{
			if (histogram[i] > maxHistValue) {
				maxHistValue = histogram[i];
				maxHist = i;
			}
		}

		for (int j = 0; j < regsize; j++) {
			int bin = regionPix[j]->getIntensity();
			if (bin == maxHist)
				binPixels.push_back(regionPix[j]);
		}

		//	std::cout << " maxHist:  " << maxHist << std::endl;
		// calc mode from maxhHist
		cv::Vec3i sum(0, 0, 0); int ave = 0; int min = 10000;
		int totred = 0; int totgreen = 0; int totblue = 0;
		float diff = 0;
		if (binPixels.size() <= 0 || maxHist <= 20 || maxHist >= 235)
		{
			mode = calcMedianColor();
		}
		else
		{
			for (int k = 0; k < binPixels.size(); k++) {
				totred += binPixels[k]->getColor().val[0];
				totgreen += binPixels[k]->getColor().val[1];
				totblue += binPixels[k]->getColor().val[2];
			}

			// ind the closest pxel to ave value in mode
			int red = (totred / binPixels.size());
			if (red > 255) red = 255; if (red < 0) red = 0;
			int green = (totgreen / binPixels.size());
			if (green > 255) green = 255; if (green < 0) green = 0;
			int blue = (totblue / binPixels.size());
			if (blue > 255) blue = 255; if (blue < 0) blue = 0;
			mode.val[0] = red; mode.val[1] = green; mode.val[2] = blue;
		}
		return mode;
	}
	
	cv::Vec3b Region::calcMedianColor() {
		// median
		int whereami = 0;

		int regsize = regionPix.size();
		int rcount[256];
		int gcount[256];
		int bcount[256];
		int icount[256];
		for (int i = 0; i < 256; i++)
		{
			rcount[i] = 0;
			gcount[i] = 0;
			bcount[i] = 0;
			icount[i] = 0;
		}
		cv::Vec3b herecol(0, 0, 0);

		pixel* randpix = regionPix[rand() % regionPix.size()];
		int originalRed = randpix->getColor().val[0];
		int originalGreen = randpix->getColor().val[1];
		int originalBlue = randpix->getColor().val[2];

		int totred = 0; int totgreen = 0; int totblue = 0;
		for (int i = 0; i < regionPix.size(); i++) {
			herecol.val[0] = regionPix[i]->getColor().val[0];
			herecol.val[1] = regionPix[i]->getColor().val[1];
			herecol.val[2] = regionPix[i]->getColor().val[2];
			totred += regionPix[i]->getColor().val[0];
			totgreen += regionPix[i]->getColor().val[1];
			totblue += regionPix[i]->getColor().val[2];

			rcount[herecol.val[0]]++;
			gcount[herecol.val[1]]++;
			bcount[herecol.val[2]]++;
		}
		whereami = compress(lookoutbelow(rcount, originalRed), regsize);
		whereami = regsize / 2;
		median.val[0] = whereinhist(rcount, whereami);
		whereami = compress(lookoutbelow(gcount, originalGreen), regsize);
		whereami = regsize / 2;
		median.val[1] = whereinhist(gcount, whereami);
		whereami = compress(lookoutbelow(bcount, originalBlue), regsize);
		whereami = regsize / 2;
		median.val[2] = whereinhist(bcount, whereami);
		return median;
	}
	
	cv::Vec3b Region::calcAveColor(/*cv::Mat &img*/) {
		int totred = 0; int totgreen = 0; int totblue = 0;
		for (int i = 0; i < regionPix.size(); i++) {
			int x = regionPix[i]->getPos().x; int y = regionPix[i]->getPos().y;
			totred += regionPix[i]->getColor().val[0];
			totgreen += regionPix[i]->getColor().val[1];
			totblue += regionPix[i]->getColor().val[2];
		}
		return cv::Vec3b( uchar(totred/ regionPix.size()),uchar( totgreen/ regionPix.size()) , uchar(totblue/ regionPix.size()));
	}

	int Region::ReturnGray(/*cv::Mat &img*/) {
		cv::Vec3b ave = calcAveColor(/*img*/);	
		return int(uchar((0.3 *ave.val[2]) + (0.59 *ave.val[1]) + (0.11 * ave.val[0])));
	}
	
	void Region::SetGray(int gr) {
		this->gray = gr;
	}
	void Region::CalcPriority() {
		//prior = (abs(255 - gray) > gray) ? 255 - gray : gray;

		if (abs(255 - gray) > gray)
			prior = 255 - gray;
		else
			prior = gray;		
	}
	void Region::sortNeighbours()
	{
		float mindif = 10000;
		std::vector<std::pair<float, int>> neighbors;
		for (unsigned int e = 0; e < neighbours.size(); e++) {
			int nid = neighbours[e]; // neighbour region
			float nWeight = edges[e]->weight;
			if (nWeight < mindif) {
				mindif = nWeight;
				neighbors.push_back(std::make_pair(nWeight, nid));
			}
			sort(neighbors.begin(), neighbors.end(), largestweight);
		}
		for (unsigned int e = 0; e < neighbours.size(); e++) {
			neighbours[e] = neighbors[e].second;
		}
	}
	void Region::RecolorTo(cv::Vec3b _newC) {

	}
	void Region::checkVisibility(std::vector<pixel*> allPix)
	{
		for (int i = 0; i < regionPix.size(); i++) {
			int id = regionPix[i]->getId();
			if (regId >= 0 && allPix[id]->getRegId() != regId)
				regionPix[i]->setVisibility(false);
			if (regId >= 0 && allPix[id]->getRegId() == regId)
				regionPix[i]->setVisibility(true);
		}
	}

}
