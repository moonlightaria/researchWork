#include "Region_Graph.h"
//#include "Matching.cpp"


bool findId(pixel *pix, std::vector<pixel*> pixels) {
	bool flag = false;
	for (int i = 0; i < pixels.size(); i++) {
		if (pix->getId() == pixels[i]->getId())
			flag = true;
	}
	return flag;
}
cv::Vec3b findMode(std::vector<cv::Vec3b> data) {
	cv::Vec3b number = data[0];
	cv::Vec3b mode = number;
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

std::vector<cv::Vec3b> findModesCol(std::vector<cv::Vec3b> data, std::vector<cv::Vec3b> palette) {
	std::vector<cv::Vec3b> modes;
	std::vector<int> votes;
	int ind; int max = -10; int ind1; int ind2;
	std::vector<int> bin;
	for (int p = 0; p < palette.size(); p++)
		bin.push_back(0);
	for (int i = 0; i <data.size(); i++) {
		cv::Vec3b pcol = data[i];
		for (int p = 0; p < palette.size(); p++) {
			if (palette[p] == pcol)
				ind = p;
		}
		bin[ind]++;
	}
	for (int i = 0; i < bin.size(); i++)
	{
		if (bin[i] > max)
		{
			max = bin[i];
			ind1 = i;
		}
	}
	max = -10;
	votes.push_back(bin[ind1]);
	for (int i = 0; i < bin.size(); i++)
	{
		if (i != ind1) {
			if (bin[i] > max)
			{
				max = bin[i];
				ind2 = i;
			}
		}
	}
	votes.push_back(bin[ind2]);
	//if (bin[ind1] >= 15) {
	//	modes.push_back(cv::Vec3b(0, 0, 0));
	//	modes.push_back(cv::Vec3b(0, 0, 0));
	//}
	//else {
	modes.push_back(palette[ind1]);
	modes.push_back(palette[ind2]);
	//}

	return modes;
}
std::vector<int> findModes(std::vector<cv::Vec3b> data, std::vector<cv::Vec3b> palette) {
	std::vector<cv::Vec3b> modes;
	std::vector<int> votes;
	int ind; int max = -10; int ind1; int ind2;
	std::vector<int> bin;
	for (int p = 0; p < palette.size(); p++)
		bin.push_back(0);
	for (int i = 0; i <data.size(); i++) {
		cv::Vec3b pcol = data[i];
		for (int p = 0; p < palette.size(); p++) {
			if (palette[p] == pcol)
				ind = p;
		}
		bin[ind]++;
	}
	for (int i = 0; i < bin.size(); i++)
	{
		if (bin[i] > max)
		{
			max = bin[i];
			ind1 = i;
		}
	}
	max = -10;
	votes.push_back(bin[ind1]);
	for (int i = 0; i < bin.size(); i++)
	{
		if (i != ind1) {
			if (bin[i] > max)
			{
				max = bin[i];
				ind2 = i;
			}
		}
	}
	votes.push_back(bin[ind2]);
	//if (bin[ind1] >= 15) {
	//	modes.push_back(cv::Vec3b(0, 0, 0));
	//	modes.push_back(cv::Vec3b(0, 0, 0));
	//}
	//else {
		modes.push_back(palette[ind1]);
		modes.push_back(palette[ind2]);
	//}
		
	return votes/*modes*/;
}
typedef struct {
	long int distance;
	long int temp;
	int id;
} strippednode;

namespace abstraction {

	//struct RegionEdge {
	//	int first;
	//	int second;
	//	float weight;
	//	int getTheOtherSide(Region *p) {
	//		return p->getId() == first ? second : first;
	//	}
	//};
	//struct RegionEdge* newEdge(Region *p1, Region *p2) {
	//	// declare and allocate new edge  
	//	struct RegionEdge* edge = new struct RegionEdge();
	//	edge->weight = 0.0f;   // Assign weight to this edge
	//						   // Initialize left and right side of the edge  
	//	edge->first = p1->getId();
	//	edge->second = p2->getId();
	//	return(edge);
	//};

	Region_Graph::Region_Graph() {
		numReg = NULL;
		edgeweights.clear();
		edges->clear();
		adj.clear();
		RegionPixels.clear();
		graph_regions.clear();
	}
	Region_Graph::~Region_Graph() {

		for (int i = 0; i < graph_regions.size();i++) {

			graph_regions[i]->neighbours.clear();

			delete graph_regions[i];
		}
		edges->clear();
	}

	Region_Graph::Region_Graph(std::string _path ,cv::Mat &img, std::vector<pixel*> &allPix,
		std::vector<std::vector<pixel*>> &_regions ) {
		numReg = NULL;
		edgeweights.clear();
		edges->clear();
		adj.clear();
		RegionPixels.clear();
		graph_regions.clear();

		RegionPixels = allPix;
		path = _path;

		std::vector<cv::Vec3b> palette;
		cv::Mat coloredImage = cv::Mat::zeros(img.size(), CV_8UC3);;

		image = img.clone();
		/*************for each vector of pixels create a Region object ***************************/
		for (int i = 0; i < _regions.size(); ++i) {
			if (_regions[i].size() != 0) {
				Region *R = new Region(_regions[i]);
				R->colored = false;
				R->edges.clear();
				R->SetGray(R->ReturnGray(/*img*/));		
				graph_regions.push_back(R);
			}
		}
		findNeighbours2(allPix);
	}

	Region *Region_Graph::returnRegion(int _regid) {
		for (int i = 0; i < graph_regions.size();i++) {
			if (graph_regions[i]->getId() == _regid)
				return graph_regions[i];
		}
		return NULL;
	}
	
	cv::Vec3b Region_Graph::colorAverage(/*std::vector<pixel*> &region*/std::vector<int> ids)
	{
		cv::Vec3d sum(0.0, 0.0, 0.0);
		for (int p = 0; p < ids.size(); p++) {

			int x = RegionPixels[ids[p]]->getPos().x; int y = RegionPixels[ids[p]]->getPos().y;
			sum.val[0] += RegionPixels[ids[p]]->getColor().val[0];
			sum.val[1] += RegionPixels[ids[p]]->getColor().val[1];
			sum.val[2] += RegionPixels[ids[p]]->getColor().val[2];
		}
		return cv::Vec3b((float)sum.val[0] / ids.size(), (float)sum.val[1] / ids.size(), (float)sum.val[2] / ids.size());
	}

	float Region_Graph::colorDiff(Region *Reg1, Region *Reg2) {
		cv::Vec3b col1 = Reg1->calcAveColor();// Reg1->calccalcAveColor()(/*image*/);
		cv::Vec3b col2 = Reg2->calcAveColor(); //Reg2->calcAveColor(/*image*/);
		//float diff = sqrt((col1.val[0] - col2.val[0]) * (col1.val[0] - col2.val[0]) +
		//	(col1.val[1] - col2.val[1]) * (col1.val[1] - col2.val[1]) + (col1.val[2] - col2.val[2]) * (col1.val[2] - col2.val[2]));
		

		ColorSpace::Rgb Rgb1(col1.val[2], col1.val[1], col1.val[0]);
		ColorSpace::Rgb Rgb2(col2.val[2], col2.val[1], col2.val[0]);
		ColorSpace::Lab Lab1, Lab2;
		//Rgb1.To<ColorSpace::Lab>(&Lab1); Rgb2.To<ColorSpace::Lab>(&Lab2);
		double diffLab = ColorSpace::EuclideanComparison::Compare(&Rgb1, &Rgb2);

		return (float)diffLab;
	}

	void Region_Graph::DrawBarCOLOR(std::vector<int> bin, std::vector<cv::Vec3b> palette , cv::Mat &screen ,
		std::vector<int> &indexes, std::vector < std::pair<int, int>> &values_p) { // portions of the largest
		cv::Mat draw = cv::Mat(50, screen.cols, CV_8UC3, Scalar(255, 255, 255));
		int xpos1 = 0;
		//std::vector < std::pair<int, int>> values_p; // porportion, palette_id
		for (int p = 0; p < palette.size();p++) {
			int xpos2 =  (bin[p] * (screen.cols) / RegionPixels.size()) + xpos1;
			values_p.push_back(std::make_pair(abs(xpos2 - xpos1), p));
			cv::Scalar col = cv::Scalar(palette[p].val[0], palette[p].val[1], palette[p].val[2]);
			cv::rectangle(draw, cv::Point(xpos1, 10), cv::Point(xpos2, draw.rows), col, -1);
			xpos1 = xpos2+1;			
		}
		screen.push_back(draw);
		//cv::imwrite(path + "BarColor_chipmunk_cmc.png", draw);

		// Sorting the vector elements on the basis of first element of pairs in decending order.
		sort(values_p.begin(), values_p.end());

		//if (values_p.size() == 4)
		//{
		//	indexes.push_back(values_p[0].second); // the first two smallest 
		//	float dist1 = colorDifferance(palette[values_p[3].second],palette[values_p[1].second]);
		//	float dist2 = colorDifferance(palette[values_p[3].second],palette[values_p[2].second]);
		//	if(dist1 < dist2)
		//		indexes.push_back(values_p[1].second);
		//	else 
		//		indexes.push_back(values_p[2].second);
		//	
		//}
		//else {
			indexes.push_back(values_p[0].second); // the first two smallest 
			indexes.push_back(values_p[1].second);
		//}
		//std::cout << "indexes size :" << indexes.size() << std::endl;
		//for (int i = 0; i < values_p.size();i++) {
		//	std::cout <<  "val : "<< values_p[i].first << std::endl;
		//}
	}
	
	cv::Point Region_Graph::retunMidPoint(Region *reg) {
		float x = 0; float y = 0;
		int regsize = reg->regionPix.size();
		if (regsize > 1) {
			for (int j = 0; j < regsize; j++) {
				x += reg->regionPix[j]->getPos().x;
				y += reg->regionPix[j]->getPos().y;
			}
			x /= regsize; y /= regsize;
		}
		else {
			//if (regsize == 0)
			//	std::cout << "region size is 0 ...." << std::endl;
			x = reg->regionPix[0]->getPos().x;
			y = reg->regionPix[0]->getPos().y;
		}
		if (x > image.cols - 1) x = image.cols - 1;
		if (y > image.rows - 1) y = image.rows - 1;
		return cv::Point((int)x, (int)y);
	}
	
	void Region_Graph::DrawChildren(Region *root, cv::Mat &draw ) {
		struct CompareRegion
		{
			bool operator()(Region *n1, Region *n2) const
			{
				return (n1->dist) < (n2->dist); // < larger to smaller values - > small to large priority
			}
		};
		std::priority_queue< Region*, std::vector< Region*>, CompareRegion > Regions;
		Regions.empty();
		for (int i = 0; i < graph_regions.size(); i++) {
			graph_regions[i]->visited = false;
		}
		Regions.push(root);
		while (!Regions.empty()) {
			Region *up = Regions.top();
			Regions.pop();
			up->visited = true;
			for (int j = 0; j < up->children.size(); j++) {
				int nid = up->children[j]->getId();
				if (up->children[j]->visited == false &&  up == up->children[j]->parent) {
					cv::Point p1 = retunMidPoint(up); cv::Point p2 = retunMidPoint(up->children[j]);
					cv::line(draw, cv::Point(4 * p1.x, 4 * p1.y), cv::Point(4 * p2.x, 4 * p2.y), cv::Scalar(0, 0, 0), 2, 8);
					Regions.push(up->children[j]);
				}
				
			}
		}
	}

	void Region_Graph::AssignColors(Region *root, std::vector<cv::Vec3b> palette, cv::Mat &table, std::vector<int> lookUp_table,int binsize ,Region *node) {
		
		//int indx = ReturnClosestPaletteIndex1(root, palette, table, lookUp_table, node);
		int indx=0;
		//if (node->regionPix.size() > 200 || node->parent->regionPix.size()> 200 )
		//	indx = ReturnClosestPaletteIndexForLargeRegions(root, palette, table, lookUp_table, node);
		//else
			indx = ReturnClosestPaletteIndex1(root, palette, table, lookUp_table, binsize, node);

		node->setNewColor(palette[indx]);
		node->assigned_colors.push_back(palette[indx]);
		node->colored = true;

		if (node != root)
			treeEdges.push_back(colorDiff(node, node->parent));

		for (int j = 0; j < node->children.size(); j++) {
			Region *neighb = node->children[j];
			if (neighb->colored == false) {
				AssignColors(root, palette, table, lookUp_table, binsize, neighb);
			}
		}


		/*
		struct CompareRegion
		{
			bool operator()(Region *n1, Region *n2) const
			{
				return (n1->dist) < (n2->dist); // < larger to smaller values - > small to large priority
			}
		};
		std::priority_queue< Region*, std::vector< Region*>, CompareRegion > Regions;
		Regions.empty();
		for (int i = 0; i < graph_regions.size(); i++) {
			graph_regions[i]->visited = false;
			graph_regions[i]->colored = false;
			if (graph_regions[i]->getId() == root->getId()) {
				int indx = ReturnClosestPaletteIndex1(root, palette, table, lookUp_table, graph_regions[i]);
				graph_regions[i]->setNewColor(palette[indx]);
				graph_regions[i]->colored = true;
			}
			Regions.push(graph_regions[i]);
		}

		//Regions.push(root);
		while (!Regions.empty()) {
			Region *up = Regions.top();
			Regions.pop();
			if (up->colored) {
				//int indx = ReturnClosestPaletteIndex1(root, palette, table, lookUp_table, up);
				//up->setNewColor(palette[indx]);
				//up->colored = true;
				for (int j = 0; j < up->edges.size(); j++) {	
					int nid = up->edges[j]->second;
					Region *neighb = returnRegion(nid);
					if (neighb->colored == false && up->getId() == neighb->parent->getId() ) {
						int indx = ReturnClosestPaletteIndex1(root, palette, table, lookUp_table, neighb);
						neighb->setNewColor(palette[indx]);
						neighb->colored = true;
						//Regions.push(up->children[j]);
					}
				}
			}
			else {
				Regions.push(up);
			}
		}
		*/
	}

	void Region_Graph::recolorWithModeofAssignedColors(cv::Mat &draw, std::vector<cv::Vec3b> palette) {
		std::vector<cv::Vec3b> mapcolors; 
		mapcolors.push_back(cv::Vec3b(0,0,230)); mapcolors.push_back(cv::Vec3b(77, 77, 255));
		mapcolors.push_back(cv::Vec3b(153, 153, 255)); mapcolors.push_back(cv::Vec3b(77, 195, 255));
		mapcolors.push_back(cv::Vec3b(0, 170, 255)); mapcolors.push_back(cv::Vec3b(153, 255, 153));
		mapcolors.push_back(cv::Vec3b(26, 255, 26)); mapcolors.push_back(cv::Vec3b(60, 179, 0));
		ofstream color_votes1(path + "/color_votes1.csv"); ofstream color_votes2(path + "/color_votes2.csv");
		for (int i = 0; i < graph_regions.size(); i++) {
			std::vector<cv::Vec3b>  Col_modes = findModesCol(graph_regions[i]->assigned_colors, palette);
			std::vector<int>  vote_modes = findModes(graph_regions[i]->assigned_colors, palette);
			//cv::Vec3b mode = findMode(graph_regions[i]->assigned_colors);
			cv::Vec3b col;
			if ((vote_modes[0] - vote_modes[1]) <= 2) col = mapcolors[0];
			if ((vote_modes[0] - vote_modes[1]) > 2 && (vote_modes[0] - vote_modes[1]) <= 4)  col = mapcolors[1];
			if ((vote_modes[0] - vote_modes[1]) > 4 && (vote_modes[0] - vote_modes[1]) <= 6)  col = mapcolors[2];
			if ((vote_modes[0] - vote_modes[1]) > 6 && (vote_modes[0] - vote_modes[1]) <= 8)  col = mapcolors[3];
			if ((vote_modes[0] - vote_modes[1]) > 8 && (vote_modes[0] - vote_modes[1]) <= 10)  col = mapcolors[4];
			if ((vote_modes[0] - vote_modes[1]) > 10 && (vote_modes[0] - vote_modes[1]) <= 12)  col = mapcolors[5];
			if ((vote_modes[0] - vote_modes[1]) > 12 && (vote_modes[0] - vote_modes[1]) <= 14)  col = mapcolors[6];
			if ((vote_modes[0] - vote_modes[1]) > 14 && (vote_modes[0] - vote_modes[1]) <= 16)  col = mapcolors[7];
			for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
			{
				int x = graph_regions[i]->regionPix[j]->getPos().x;
				int y = graph_regions[i]->regionPix[j]->getPos().y;				
				draw.at<cv::Vec3b>(y, x) = col;

				if ( (vote_modes[0] - vote_modes[1]) >= 15) {
					graph_regions[i]->setNewColor(Col_modes[0]);
					graph_regions[i]->colored = true;
				}
			}			
			if (color_votes1.is_open())
			{
				color_votes1 << vote_modes[0] << std::endl;
			}
			if (color_votes2.is_open())
			{
				color_votes2 << vote_modes[1] << std::endl;
			}
			graph_regions[i]->assigned_colors.clear();
		}
		color_votes1.close();
		color_votes2.close();
	}

	void Region_Graph::DrawTree(Region *root,cv::Mat &draw, int count) {
		float px = 0.0; float py = 0;
		float x = 0; float y = 0;

		for (int i = 0; i < root->neighbours.size(); i++) {
			uchar colorb = rand() % 255; uchar colorg = rand() % 255; uchar colorr = rand() % 255;
			int nid = root->neighbours[i];
			if (root == returnRegion(nid)->parent) {
				int regsize = root->regionPix.size();
				
				if (regsize > 0) {
					for (int j = 0; j < regsize; j++) {
						x += root->regionPix[j]->getPos().x;
						y += root->regionPix[j]->getPos().y;
					}
					x /= regsize; y /= regsize;
				}

				int childsize = returnRegion(nid)->regionPix.size();
				if (childsize > 0) {
					for (int p = 0; p < childsize; p++) {
						px += returnRegion(nid)->regionPix[p]->getPos().x;
						py += returnRegion(nid)->regionPix[p]->getPos().y;
					}
					px /= childsize; py /= childsize;
				}
				cv::line(draw, cv::Point(4 * px,4*  py), cv::Point(4 *x,4* y), cv::Scalar(0, 0, 0), 2, 8);

				//DrawTree(returnRegion(nid), draw,count++);
			}
		}
		for (int i = 0; i < root->neighbours.size(); i++) {
			int nid = root->neighbours[i];
			if (root == returnRegion(nid)->parent) {
				DrawTree(returnRegion(nid), draw, count);
			}
		}

		

	/*	for (int i = 0; i < graph_regions.size(); i++) {
			uchar colorb = rand() % 255; uchar colorg = rand() % 255; uchar colorr = rand() % 255;
			for (int j = 0; j < graph_regions[i]->neighbours.size(); j++) {
				int regsize = graph_regions[i]->regionPix.size();
				int nid = graph_regions[i]->neighbours[j];
				if (graph_regions[i] == returnRegion(nid)->parent) {	
					if (regsize > 0) {
						for (int p = 0; p < regsize; p++) {
							x += graph_regions[i]->regionPix[p]->getPos().x;
							y += graph_regions[i]->regionPix[p]->getPos().y;
						}
						x /= regsize; y /= regsize;
					}
					int childsize = returnRegion(nid)->regionPix.size();
					if (childsize > 0) {
						for (int p = 0; p <  childsize; p++) {
							px += returnRegion(nid)->regionPix[p]->getPos().x;
							py += returnRegion(nid)->regionPix[p]->getPos().y;
						}
						px /= childsize; py /= childsize;
					}
					if (regsize > 0 && childsize > 0)
						cv::line(draw, cv::Point(px, py), cv::Point(x, y), cv::Scalar(colorb,colorg,colorr), 1, 8);
				}
			}
		}*/

	}
	
	float Region_Graph::getMidWeight(std::vector<float> weights, std::vector<int> &lookUp_table)
	{
		std::vector<int> nweights;
		if (lookUp_table.empty()) {
			for (int i = 0; i < weights.size(); i++) {
				nweights.push_back(weights[i]);
			}
		}
		else {
			for (int i = 0; i < weights.size(); i++) {
				nweights.push_back(weights[i]/*lookUp_table[weights[i]]*/);
			}
		}
		float midweight = 0.0f;
		if(nweights.size() <= 1){
			midweight = 0.0f;
		}
		else {
			int maxw = *max_element(nweights.begin(), nweights.end());
			int minw = *min_element(nweights.begin(), nweights.end());
			midweight = (maxw + minw) / 2.0f;
		}
		return midweight;
	}
	
	void Region_Graph::make(/*std::vector<Region*> &_regions*/) {
		numReg = graph_regions.size();
		for (int i = 0; i < numReg; i++) {
			std::vector<int> vec;
			adj.push_back(vec);
		}
		for (int u = 0; u < numReg; ++u)
		{
			for (int v = u; v < numReg; ++v)
			{
				if (u != v) {
					/******check if they are neighbor regions********/
					bool neighb = false;
					//overl = checkOverlaps(_regions[u], _regions[v]);

					if (neighb == true) {
						addEdge(adj, graph_regions[u]->getId(), graph_regions[v]->getId());
					}
				}
			}
		}
	}

	bool Region_Graph::checkOverlaps(Region* reg1, Region* reg2) {

		if (reg1->DoBoxesIntersect(reg2) == false)
			return false;

		//bool overlap = false;
		//std::vector<pixel*> overlaped;
		for (int i = 0; i < reg1->regionPix.size(); ++i) {
			for (int j = 0; j < reg2->regionPix.size(); ++j) {
				if (reg1->regionPix[i]->id == reg2->regionPix[j]->id) {
					//overlaped.push_back(reg2->regionPix[i]);//
					//overlap = true;//
					return true;
				}
			}
		}
		//if (overlaped.size() > 0) {
		//	overlap = true;
		//}
		return false;
	}
	
	void Region_Graph::addEdge(std::vector<std::vector<int>> &adj, int u, int v)
	{
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	
	void Region_Graph::addEdge(Region *Reg, int v)
	{
		if (returnRegion(v)->regionPix.size() > 0) {

			Reg->neighbours.push_back(v);
			RegionEdge *temp = new RegionEdge();
			temp->first = Reg->getId();
			temp->second = v;

			float w = colorDiff(Reg, returnRegion(v));
			temp->weight = w;
			Reg->edges.push_back(temp);
		}
	}

	void Region_Graph::findNeighbours(std::vector<pixel*> &allPix) {
		for (int i = 0; i < graph_regions.size(); i++) {
			if (graph_regions[i]->regionPix.size() > 0) {
				std::vector<int> vec;
				adj.push_back(vec);
			}
		}
		// look for the different regions ids around this region // first find the contours
		//std::vector<float> edgeweights;
		for (int i = 0; i < graph_regions.size(); i++) {

			cv::Mat gbp = cv::Mat::zeros(image.size(), CV_8UC1);
			cv::Mat contourImage = cv::Mat::zeros(image.size(), CV_8UC3);
			std::vector<cv::Vec4i> hierarchy;
			std::vector<std::vector<cv::Point>> contours;
			colorRegion(graph_regions[i], contourImage);
			cvtColor(contourImage, gbp, cv::COLOR_BGR2GRAY);
			findContours(gbp, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cv::Point(0, 0));
			if (contours.size() == 0) {
				for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
					graph_regions[i]->contourPix.push_back(graph_regions[i]->regionPix[j]);
			}
			else
			{
				if (contours[0].size() > 4) {
					for (int j = 0; j < contours[0].size(); j++) {
						int id = contours[0][j].x + contours[0][j].y * image.cols;
						graph_regions[i]->contourPix.push_back(allPix[id]);
					}
				}
				else {
					for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
						graph_regions[i]->contourPix.push_back(graph_regions[i]->regionPix[j]);
				}
			}

			int regid = graph_regions[i]->getId();
			std::vector<int> neighbId; // neighbour regions

			for (int p = 0; p < graph_regions[i]->contourPix.size(); p++) {
				int pid = graph_regions[i]->contourPix[p]->id;
				for (int n = 0; n < allPix[pid]->neighbours.size(); n++) {
					if (allPix[pid]->neighbours[n]->CCid != allPix[pid]->CCid)
					{
						neighbId.push_back(allPix[pid]->neighbours[n]->CCid);
					}
				}
			}

			//remove duplicate neighbour ids
			std::vector<int>::iterator ip;
			// Sorting the array 
			std::sort(neighbId.begin(), neighbId.end());
			// Using std::unique 
			ip = std::unique(neighbId.begin(), neighbId.begin() + neighbId.size());
			// Resizing the vector so as to remove the undefined terms 
			neighbId.resize(std::distance(neighbId.begin(), ip));

			for (int n = 0; n < neighbId.size(); n++) {
				//addEdge(adj, regid, neighbId[n]);
				addEdge(graph_regions[i], neighbId[n]);
				edgeweights.push_back(graph_regions[i]->edges[n]->weight);
			}
			neighbId.clear();
			
		}
		/******histogram of edge weights = color differences***********/
		Histogram(edgeweights, path);
	
		//remove duplicates from adj list
		for (int i = 0; i < adj.size(); i++) {
			std::vector<int>::iterator ip;
			// Sorting the array 
			std::sort(adj[i].begin(), adj[i].end());
			// Using std::unique 
			ip = std::unique(adj[i].begin(), adj[i].begin() + adj[i].size());
			// Resizing the vector so as to remove the undefined terms 
			adj[i].resize(std::distance(adj[i].begin(), ip));
		}
	}

	void Region_Graph::findNeighbours2(std::vector<pixel*> &allPix) {

		/*****************************************************************/
		for (int i = 0; i < graph_regions.size(); i++) {
			if (graph_regions[i]->regionPix.size() != 0) {
				std::vector<int> vec;
				adj.push_back(vec);
			}
		}
		// look for the different regions ids around this region // first find the contours
		/******************find neighbours by contours*******************/
		//for (int i = 0; i < graph_regions.size(); i++) {
		//	cv::Mat gbp = cv::Mat::zeros(image.size(), CV_8UC1);
		//	cv::Mat contourImage = cv::Mat::zeros(image.size(), CV_8UC3);
		//	std::vector<cv::Vec4i> hierarchy;
		//	std::vector<std::vector<cv::Point>> contours;
		//	colorRegion(graph_regions[i], contourImage);
		//	cvtColor(contourImage, gbp, cv::COLOR_BGR2GRAY);
		//	findContours(gbp, contours, hierarchy,CV_RETR_TREE,/*CV_RETR_EXTERNAL,*/ CV_CHAIN_APPROX_NONE, cv::Point(0, 0));
		//	
		//	graph_regions[i]->ListOfcontourPix.resize(contours.size());
		//	for (int a = 0; a < contours.size(); a++) {
		//		if (contours[a].size() > 4) {
		//			for (int j = 0; j < contours[a].size(); j++) {
		//				int id = contours[a][j].x + contours[a][j].y * image.cols;
		//				graph_regions[i]->ListOfcontourPix[a].push_back(allPix[id]);
		//			}
		//		}
		//		else {
		//			for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
		//				graph_regions[i]->ListOfcontourPix[a].push_back(graph_regions[i]->regionPix[j]);
		//		}
		//	}
		//	int regid = graph_regions[i]->getId();
		//	std::vector<int> neighbId; // neighbour regions
		//	/*****************CCId********************/
		//	for (int j = 0; j < graph_regions[i]->ListOfcontourPix.size(); j++) {
		//		for (int p = 0; p < graph_regions[i]->ListOfcontourPix[j].size(); p++) {
		//			int pid = graph_regions[i]->ListOfcontourPix[j][p]->id;
		//			for (int n = 0; n < allPix[pid]->neighbours.size(); n++) {
		//				if (allPix[pid]->neighbours[n]->CCid != allPix[pid]->CCid)
		//				{
		//					neighbId.push_back(allPix[pid]->neighbours[n]->CCid);
		//				}
		//			}
		//		}
		//	}
		//	/////*************slicid**************/
		//	//for (int j = 0; j < graph_regions[i]->ListOfcontourPix.size(); j++) {
		//	//	for (int p = 0; p < graph_regions[i]->ListOfcontourPix[j].size(); p++) {
		//	//		int pid = graph_regions[i]->ListOfcontourPix[j][p]->id;
		//	//		for (int n = 0; n < allPix[pid]->neighbours.size(); n++) {
		//	//			if (allPix[pid]->neighbours[n]->slicId != allPix[pid]->slicId)
		//	//			{
		//	//				neighbId.push_back(allPix[pid]->neighbours[n]->slicId);
		//	//			}
		//	//		}
		//	//	}
		//	//}
		//	//for (int a = 0; a < contours.size();a++) {
		//	//	if (contours[a].size() > 4) {
		//	//		for (int j = 0; j < contours[a].size(); j++) {
		//	//			int id = contours[a][j].x + contours[a][j].y * image.cols;
		//	//			graph_regions[i]->contourPix.push_back(allPix[id]);
		//	//		}
		//	//	}
		//	//	else {
		//	//		for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
		//	//			graph_regions[i]->contourPix.push_back(graph_regions[i]->regionPix[j]);
		//	//	}
		//	//}
		//	//int regid = graph_regions[i]->getId();
		//	//std::vector<int> neighbId; // neighbour regions
		//	//for (int p = 0; p < graph_regions[i]->contourPix.size(); p++) {
		//	//	int pid = graph_regions[i]->contourPix[p]->id;
		//	//	for (int n = 0; n < allPix[pid]->neighbours.size(); n++) {
		//	//		if (allPix[pid]->neighbours[n]->CCid != allPix[pid]->CCid)
		//	//		{
		//	//			neighbId.push_back(allPix[pid]->neighbours[n]->CCid);
		//	//		}
		//	//	}
		//	//}
		//	//remove duplicate neighbour ids
		//	std::vector<int>::iterator ip;
		//	// Sorting the array 
		//	std::sort(neighbId.begin(), neighbId.end());
		//	// Using std::unique 
		//	ip = std::unique(neighbId.begin(), neighbId.begin() + neighbId.size());
		//	// Resizing the vector so as to remove the undefined terms 
		//	neighbId.resize(std::distance(neighbId.begin(), ip));
		//	for (int n = 0; n < neighbId.size(); n++) {
		//		addEdge(graph_regions[i], neighbId[n]);
		//		edgeweights.push_back(graph_regions[i]->edges[n]->weight);
		//	}
		//	neighbId.clear();
		//}
		/*********************find neigbor regions by pixels neighbors*******************/
		for (int i = 0; i < graph_regions.size(); i++) {
			if (graph_regions[i]->regionPix.size() != 0) {
				std::vector<int> neighbId; // neighbour regions
				for (int j = 0; j < graph_regions[i]->regionPix.size(); j++) {
					int pid = graph_regions[i]->regionPix[j]->getId();
					for (int n = 0; n < allPix[pid]->neighbours.size(); n++) {
						if (allPix[pid]->neighbours[n] != NULL && allPix[pid]->neighbours[n]->CCid != allPix[pid]->CCid)
						{
							neighbId.push_back(allPix[pid]->neighbours[n]->CCid);
						}
					}
				}
				//remove duplicate neighbour ids
				std::vector<int>::iterator ip;
				// Sorting the array 
				std::sort(neighbId.begin(), neighbId.end());
				// Using std::unique 
				ip = std::unique(neighbId.begin(), neighbId.begin() + neighbId.size());
				// Resizing the vector so as to remove the undefined terms 
				neighbId.resize(std::distance(neighbId.begin(), ip));

				for (int n = 0; n < neighbId.size(); n++) {
					addEdge(graph_regions[i], neighbId[n]);
					edgeweights.push_back(graph_regions[i]->edges[n]->weight);
				}
				neighbId.clear();
			}
		}

		/******histogram of edge weights = color differences***********/
		Histogram(edgeweights, path);

		//remove duplicates from adj list
		for (int i = 0; i < adj.size(); i++) {
			std::vector<int>::iterator ip;
			// Sorting the array 
			std::sort(adj[i].begin(), adj[i].end());
			// Using std::unique 
			ip = std::unique(adj[i].begin(), adj[i].begin() + adj[i].size());
			// Resizing the vector so as to remove the undefined terms 
			adj[i].resize(std::distance(adj[i].begin(), ip));
		}
	}


	void Region_Graph::quaitizeImageColorTo(std::vector<cv::Vec3b> palette, cv::Mat &coloredImage) {
	
		//barn color
		//palette.push_back(cv::Vec3b(72, 189, 237));/* yellow*/ palette.push_back(cv::Vec3b(17, 46, 121));/*brown*/
		//palette.push_back(cv::Vec3b(39, 102, 148)); /*sienna*/
		//palette.push_back(cv::Vec3b(230, 197, 165)); /*light blue*/ palette.push_back(cv::Vec3b(33, 34, 41)); /*black*/
	
		//art deco color palette
		palette.push_back(cv::Vec3b(174, 210, 240)); 
		palette.push_back(cv::Vec3b(30, 114, 166));
		palette.push_back(cv::Vec3b(90, 133, 90)); 
		palette.push_back(cv::Vec3b(49, 83, 43));
		palette.push_back(cv::Vec3b(66, 66, 214));
		palette.push_back(cv::Vec3b(51, 59, 32));

		for (int i = 0; i < graph_regions.size(); i++) {
			float mindif = 1000.0f;
			int index;
			cv::Vec3b col = graph_regions[i]->calcAveColor();// graph_regions[i]->calcAveColor(/*image*/);
			for (int j = 0; j < palette.size(); j++) {
				float dif = colorDifferance(col, palette[j]);
				if (dif < mindif) {
					mindif = dif;
					index = j;
				}
			}
			colorRegion(graph_regions[i], coloredImage, palette[index]);
		}

	}

	void Region_Graph::colorRegion(Region *reg, cv::Mat  &contourImage) {
		contourImage = cv::Mat::zeros(image.size(), CV_8UC3);
		for (int i = 0; i<reg->regionPix.size(); i++)
		{
			int x = reg->regionPix[i]->getPos().x; int y = reg->regionPix[i]->getPos().y;
			contourImage.at<cv::Vec3b>(y, x) = reg->regionPix[i]->getColor();
		}
	}
	
	void Region_Graph::colorRegion(Region *reg, cv::Mat  &contourImage, cv::Vec3b palettecolor) {
		
		for (int i = 0; i<reg->regionPix.size(); i++)
		{
			int x = reg->regionPix[i]->getPos().x; int y = reg->regionPix[i]->getPos().y;
			reg->regionPix[i]->setColor(palettecolor);
			contourImage.at<cv::Vec3b>(y, x) = palettecolor;
		}
	}

	void Region_Graph::drawEdges(cv::Mat &screen)
	{
		for (int i = 0; i < graph_regions.size();i++) {
			float x = 0; float y = 0;
			for (int j = 0; j < graph_regions[i]->regionPix.size(); j++) {
				x += graph_regions[i]->regionPix[j]->getPos().x;
				y += graph_regions[i]->regionPix[j]->getPos().y;
			}
			x /= graph_regions[i]->regionPix.size();
			y /= graph_regions[i]->regionPix.size();
			
			for (int e = 0; e < graph_regions[i]->neighbours.size(); e++) {
				float xx = 0; float yy = 0;
				int nid = graph_regions[i]->edges[e]->second;
				std::vector<pixel*> regpix = returnRegion(nid)->regionPix;
				for (int n = 0; n <regpix.size(); n++) {
					xx += regpix[n]->getPos().x;
					yy += regpix[n]->getPos().y;
				}
				xx /= regpix.size();
				yy /= regpix.size();

				cv::line(screen, cv::Point(2 * x, 2 * y), cv::Point(2 * xx, 2 * yy), cvScalar(255, 0, 0), 1);
			}
			cv::circle(screen, cv::Point(2 * x, 2 * y), 1, cvScalar(0, 0, 255), 1);
		}
	}
	
	std::pair<int, int> Region_Graph::ReturnClosestPalette(RegionEdge * edge,cv::Mat &table,
		std::vector<cv::Vec3b> palette, std::vector<int> lookUp_table) {
		float mindif = 1000.0f;
		std::pair<int, int> index;
		bool signr = false;
		bool signp = false;
		float colorDif =  lookUp_table[edge->weight];
		/*************************************/
		cv::Vec3b col1 = returnRegion(edge->first)->calcAveColor();// calcAveColor();
		cv::Vec3b col2 = returnRegion(edge->second)->calcAveColor();// calcAveColor();
		int g1 = int(uchar((0.3 *col1.val[2]) + (0.59 *col1.val[1]) + (0.11 * col1.val[0])));	// return gray of average color	
		int g2 = int(uchar((0.3 *col2.val[2]) + (0.59 *col2.val[1]) + (0.11 * col2.val[0])));	// return gray of average color	
	
		if (g1 > g2) signr = true;//positive
		else signr = false;//negative

		for (int i = 0; i < table.rows; i++) {
			for (int j = 0; j < table.cols; j++) {
				 float dif = abs(colorDif - table.at<float>(i, j));
				 cv::Vec3b pal1 = palette[i];  cv::Vec3b pal2 = palette[j];
				 //std::cout <<"diff:"<< dif << " , edge weight: " << colorDif << " , table: "<< table.at<float>(i, j) << std::endl;
				 int p1 = int(uchar((0.3 *pal1.val[2]) + (0.59 *pal1.val[1]) + (0.11 * pal1.val[0])));	// return gray of average color	
				 int p2 = int(uchar((0.3 *pal2.val[2]) + (0.59 *pal2.val[1]) + (0.11 * pal2.val[0])));	// return gray of average color
				 if (p1 > p2) signp = true;
				 else signp = false;
				 if (dif < mindif /*&& signr == signp*/) {
					 mindif = dif;
					 index.first = i;
					 index.second = j;
				 }	
				 //if (dif < mindif && signr != signp) {
					// mindif = dif;
					// index.first = j;
					// index.second = i;
				 //}
			}
		}
		return index;
	}
	
	int Region_Graph::ReturnClosestPaletteIndex(RegionEdge *edge, std::vector<cv::Vec3b> palette, std::vector<int> lookUp_table) {
		float mindif = 1000.0f;
		//********************
		int index; cv::Vec3b col, col2;
		float colorDif = lookUp_table[edge->weight];
		if (returnRegion(edge->first)->colored == false) {
			col = returnRegion(edge->first)->calcAveColor();// calccalcAveColor()(/*image*/);
			col2 = returnRegion(edge->second)->calcAveColor();// calcAveColor(/*image*/);
		}
		if (returnRegion(edge->second)->colored == false) {
			col = returnRegion(edge->second)->calcAveColor();// calcAveColor(/*image*/);
			col2 = returnRegion(edge->first)->calcAveColor();//calcAveColor(/*image*/);
		}
		int g = int(uchar((0.3 *col.val[2]) + (0.59 *col.val[1]) + (0.11 * col.val[0])));	// return gray of average color	
		std::vector<int> pals;
		for (int i = 0; i < palette.size(); i++) { 
			int p = int(uchar((0.3 *palette[i].val[2]) + (0.59 *palette[i].val[1]) + (0.11 * palette[i].val[0])));
			float diff = colorDifferance(col, palette[i]);// abs(g - p);// lookUp_table[colorDifferance(col, palette[i])];//
			//float dif = abs(colorDif - diff);
			if (diff < mindif) {
				mindif = diff;
				index = i;
			}
		}
		return index;
	}

	int Region_Graph::ReturnClosestPaletteIndexForLargeRegions(Region * source, std::vector<cv::Vec3b> palette, cv::Mat &table, std::vector<int> lookUp_table, Region * curReg) {
		float mindif = 1000.0f;
		int index = -1;
		bool signr = false;
		bool signp = false;
		float g1, g2; float p1, p2;
		int pidex = -1; //parent color index 
		float colordiff; float colorDif;
		cv::Vec3b col1, col2;
		int regid = curReg->getId();
		col2 = returnRegion(regid)->calcAveColor();// calcAveColor();
		/**could be the first region or the ones didn't get a parent***/
		if (pidex == -1 || curReg->parent == NULL || curReg->getId() == source->getId()) {
			for (int i = 0; i < curReg->neighbours.size(); i++) {
				int id_n = curReg->neighbours[i];
				if (returnRegion(id_n)->colored) {
					for (int j = 0; j < palette.size(); j++) {
						if (returnRegion(id_n)->getNewColor() == palette[j])
							pidex = j;
					}
				}
			}
			if (pidex == -1) { // or find the closest pair
				for (int i = 0; i < palette.size(); i++) {
					float diff = colorDifferance(col2, palette[i]);
					if (diff < mindif) {
						mindif = diff;
						index = i;
					}
				}
			}
		}
		if (curReg->parent != NULL) {
			cv::Mat gbp = cv::Mat::zeros(image.size(), CV_8UC1);
			cv::Mat contourImage = cv::Mat::zeros(image.size(), CV_8UC3);
			std::vector<cv::Vec4i> hierarchy;
			std::vector<std::vector<cv::Point>> contours;
			colorRegion(curReg, contourImage);
			cvtColor(contourImage, gbp, cv::COLOR_BGR2GRAY);
			findContours(gbp, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
			sort(contours.begin(), contours.end(), [](const std::vector<cv::Point> & a, const std::vector<cv::Point> & b) { return a.size() > b.size(); });
			std::vector<pixel*> boundary;
			for (int a = 0; a < contours[0].size(); a++) {
				if (contours[0].size() > 4) {
					for (int j = 0; j < contours[0].size(); j++) {
						int idd = contours[0][j].x + contours[0][j].y * image.cols;
						boundary.push_back(RegionPixels[idd]);
					}
				}
				else {
					for (int j = 0; j < curReg->regionPix.size(); j++)
						boundary.push_back(curReg->regionPix[j]);
				}
			}
			
			int parent_id = curReg->parent->getId();
			std::vector<int> cur_pix; // current region pixels
			std::vector<int> neighbId; // neighbour region pixels
			/*****************neighbor pixels to shared contour********************/
			for (int j = 0; j < boundary.size(); j++) {
				int x = boundary[j]->getPos().x; int y = boundary[j]->getPos().y;
				for (int xx = x -10; xx < x + 11; xx++) {
					for (int yy = y - 10; yy < y +11; yy++) {
						if (xx >= 0 && xx < image.cols && yy >=0 && yy < image.rows) {
							int nid = xx + yy * image.cols;
							if (RegionPixels[nid]->CCid == regid /*&& findId(RegionPixels[nid], cur_pix) == false*/) {
								cur_pix.push_back(nid);
							}
							if (RegionPixels[nid]->CCid == parent_id /*&& findId(RegionPixels[nid], neighbId) == false*/) {
								neighbId.push_back(nid);
							}
						}
					}
				}			
			}
			//remove duplicate neighbour ids
			std::vector<int>::iterator ip, ipp;
			// Sorting the array 
			std::sort(neighbId.begin(), neighbId.end());
			// Using std::unique 
			ip = std::unique(neighbId.begin(), neighbId.begin() + neighbId.size());
			// Resizing the vector so as to remove the undefined terms 
			neighbId.resize(std::distance(neighbId.begin(), ip));
			//remove duplicate neighbour ids
			// Sorting the array 
			std::sort(cur_pix.begin(), cur_pix.end());
			// Using std::unique 
			ipp = std::unique(cur_pix.begin(), cur_pix.begin() + cur_pix.size());
			// Resizing the vector so as to remove the undefined terms 
			cur_pix.resize(std::distance(cur_pix.begin(), ipp));
			/******************/
			int pid = curReg->parent->getId(); // parent id
			col1 = colorAverage(neighbId);// boundary in parent side
			col2 = colorAverage(cur_pix); //boundary in current region side
			for (int i = 0; i < palette.size(); i++) {
				if (curReg->parent->getNewColor() == palette[i])
					pidex = i;
			}
			p1 = (((0.2126 *palette[pidex].val[2]) + (0.7152 *palette[pidex].val[1]) + (0.0722 * palette[pidex].val[0])));
			g1 = (((0.2126 *col1.val[2]) + (0.7152 *col1.val[1]) + (0.0722 * col1.val[0])));	// return luminance 
			g2 = (((0.2126 *col2.val[2]) + (0.7152 *col2.val[1]) + (0.0722 * col2.val[0])));	// return luminance 
			float regiondif = colorDifferance(col1, col2);
			float colorDif = 0.0;
			if (!lookUp_table.empty())
				colorDif = lookUp_table[regiondif];
			else
				colorDif = regiondif;

			if (g1 >= g2) signr = true;//positive
			else signr = false;//negative
			mindif = 1000.0f;
			for (int j = 0; j < palette.size(); j++) {
				float difp = colorDifferance(palette[pidex], palette[j]);
				p2 = (((0.2126 *palette[j].val[2]) + (0.7152 *palette[j].val[1]) + (0.0722 * palette[j].val[0])));
				if (p1 >= p2) signp = true;
				else signp = false;
				float dif = abs(colorDif - difp);
				if (signr == signp) {
					if (dif < mindif) {
						mindif = dif;
						index = j;
					}
				}
			}
			if (index == -1)
			{
				counter++;
				float mindif = 1000.0f;
				for (int j = 0; j < palette.size(); j++) {
					float difp = colorDifferance(palette[pidex], palette[j]);
					float dif = abs(colorDif - difp);
					if (dif < mindif) {
						mindif = dif;
						index = j;
					}
				}
			}
		}
		return index;
	}

	int Region_Graph::ReturnClosestPaletteIndex1(Region * source, std::vector<cv::Vec3b> palette,
		cv::Mat &table, std::vector<int> lookUp_table, int binsize, Region * curReg) {
		float mindif = 1000.0f;
		float difp;
		int index = -1;
		bool signr = false;
		bool signp = false;
		float g1, g2; float p1, p2;
		int pidex = -1; //parent color index 
		float colordiff; float colorDif;
		cv::Vec3b col1, col2;
		double miPaleDif, maxPaleDif;
		cv::minMaxLoc(table, &miPaleDif, &maxPaleDif);
		float maxEdgeDif = *max_element(edgeweights.begin(), edgeweights.end());
		//********************
		std::vector<float> tArr(palette.size()*palette.size());
		int count = 0;
		for (int i = 0; i < table.rows; i++) {
			for (int j = 0; j < table.cols; j++) {
				tArr[count] = table.at<float>(i, j);
				//std::cout << tArr[count] << std::endl;
				count++;
			}
		}

		/************/
		//int id = curReg->getId();
		col2 = curReg->calcAveColor();// calcAveColor();
		/**could be the first region or the ones didn't get a parent***/
		if (pidex == -1 || curReg->parent == NULL || curReg->getId() == source->getId()) {
			for (int i = 0; i < curReg->neighbours.size(); i++) {
				int id_n = curReg->neighbours[i];
				Region *neighb = returnRegion(id_n);
				if (neighb->colored) {
					for (int j = 0; j < palette.size(); j++) {
						if (neighb->getNewColor() == palette[j])
							pidex = j;
					}
				}
			}
			if (pidex == -1) { // or find the closest pair
				for (int i = 0; i < palette.size(); i++) {
					float diff = colorDifferance(col2, palette[i]);
					if (diff < mindif) {
						mindif = diff;
						index = i;
					}
				}
			}
		}
		if (curReg->parent != NULL) {
			//int pid = curReg->parent->getId(); // parent id
			Region *parent = curReg->parent;
			col1 = parent->calcAveColor();// calcAveColor(); // parent
			for (int i = 0; i < palette.size(); i++) {
				if (parent->getNewColor() == palette[i]) {
					pidex = i;
					break;
				}
			}

			p1 = (((0.2126 *palette[pidex].val[2]) + (0.7152 *palette[pidex].val[1]) + (0.0722 * palette[pidex].val[0])));
			g1 = (((0.2126 *col1.val[2]) + (0.7152 *col1.val[1]) + (0.0722 * col1.val[0])));	// return luminance 
			g2 = (((0.2126 *col2.val[2]) + (0.7152 *col2.val[1]) + (0.0722 * col2.val[0])));	// return luminance 
			float regiondif = colorDifferance(col1, col2);
			float colorDif = 0.0;
			if (!lookUp_table.empty())
			{
				int id = (int) round(regiondif);
				//if (id < 0) id = 0; if (id >= lookUp_table.size()) id = lookUp_table.size() - 1;
				if (id >= lookUp_table.size() || id < 0)
					std::cout << "id: " << id << std::endl;
				colorDif = lookUp_table[id];

				//colorDif = lookUp_table[(int) ((regiondif /*/ maxEdgeDif * maxPaleDif*/))];
			}
			else
				colorDif = regiondif;
			
			/**********Design1***********/
			//float mindif = 1000.0f;
			//for (int j = 0; j < palette.size(); j++) {
			//	float difp = colorDifferance(palette[pidex], palette[j]);
			//	float dif = abs(colorDif - difp);
			//	if (dif < mindif) {
			//		mindif = dif;
			//		index = j;
			//	}
			//}
			/*********Design 2************/
			if (g1 >= g2) signr = true;//positive
			else signr = false;//negative
			mindif = 1000.0f;
			for (int j = 0; j < palette.size(); j++) {
				auto palette_j = palette[j];
				difp = colorDifferance(palette[pidex], palette_j);
				p2 = (((0.2126 *palette_j.val[2]) + (0.7152 *palette_j.val[1]) + (0.0722 * palette_j.val[0])));
				if (p1 >= p2) signp = true;
				else signp = false;
				float dif = abs(colorDif - difp);
				if (signr == signp) {
					if (dif < mindif) {
						mindif = dif;
						index = j;
					}
				}
			}
			if (index == -1)
			{
				counter++;
				float mindif = 1000.0f;
				for (int j = 0; j < palette.size(); j++) {
					difp = colorDifferance(palette[pidex], palette[j]);
					float dif = abs(colorDif - difp);
					if (dif < mindif) {
						mindif = dif;
						index = j;
					}
				}
			}
			//std::cout <<"colorDif: "<< colorDif << " , difp: "<< difp << std::endl;
		}
		return index;
	}

	int Region_Graph::ReturnClosestPaletteIndex2(RegionEdge *edge, std::vector<cv::Vec3b> palette, cv::Mat &table, std::vector<int> lookUp_table, cv::Vec3b C) {
		float mindif = 1000.0f;
		int index= -1;
		bool signr = false;
		bool signp = false;
		float g1, g2;
		//********************

		cv::Vec3b col1 = returnRegion(edge->first)->calcAveColor();// calcAveColor();
		cv::Vec3b col2 = returnRegion(edge->second)->calcAveColor();// calcAveColor();
		float colorDif = lookUp_table[edge->weight];
		int pidex =-1;
		for (int i = 0; i < palette.size(); i++) {
			if (C == palette[i])
				pidex = i;
		}
		if (pidex == -1) {
			std::cout << " C is not equal to any palette:" << C << std::endl;
			std::cout << "col1: " << col1 << ", col2: " << col2 << ", C: " << C << std::endl;
		}
		//if (returnRegion(edge->first)->getNewColor() == C)
		//{
		//	g1 = (((0.2126 *col1.val[2]) + (0.7152 *col1.val[1]) + (0.0722 * col1.val[0])));	// return luminance 
		//	g2 = (((0.2126 *col2.val[2]) + (0.7152 *col2.val[1]) + (0.0722 * col2.val[0])));	// return luminance 
		//}
		//if (returnRegion(edge->second)->getNewColor() == C) {
		//	g1 = (((0.2126 *col2.val[2]) + (0.7152 *col2.val[1]) + (0.0722 * col2.val[0])));	// return luminance 
		//	g2 = (((0.2126 *col1.val[2]) + (0.7152 *col1.val[1]) + (0.0722 * col1.val[0])));	// return luminance 
		//}
		cv::Vec3b pal1 = palette[pidex];	
		float p1 = (((0.2126 *pal1.val[2]) + (0.7152 *pal1.val[1]) + (0.0722 * pal1.val[0])));	// return luminance 							
		g1 = (((0.2126 *col1.val[2]) + (0.7152 *col1.val[1]) + (0.0722 * col1.val[0])));	// return luminance 
		g2 = (((0.2126 *col2.val[2]) + (0.7152 *col2.val[1]) + (0.0722 * col2.val[0])));	// return luminance 
		/**********Design1***********/
		//if (g1 >= g2) signr = true;//positive
		//else signr = false;//negative
		//for (int j = 0; j < palette.size(); j++) {
		//	float difp = colorDifferance(palette[pidex], palette[j]);
		//	//int p2 = int(uchar((0.3 *palette[j].val[2]) + (0.59 *palette[j].val[1]) + (0.11 * palette[j].val[0])));
		//	float p2 = (((0.2126 *palette[j].val[2]) + (0.7152 *palette[j].val[1]) + (0.0722 * palette[j].val[0])));	// return gray of average color	
		//	if (p1 >= p2) signp = true;
		//	else signp = false;
		//	float dif = abs(colorDif - difp);
		//	if (dif < mindif && signr == signp) {
		//		mindif = dif;
		//		index = j;
		//	}
		//}
		//if (index == -1) {
		//	mindif = 1000.0f;
		//	for (int j = 0; j < palette.size(); j++) {
		//		float difp = colorDifferance(palette[pidex], palette[j]);
		//		float dif = abs(colorDif - difp);
		//		if (dif < mindif) {
		//			mindif = dif;
		//			index = j;
		//		}
		//	}	
		//}
		/**********Design2***********/
		//if (g1 >= g2) signr = true;//positive
		//else signr = false;//negative
		//for (int j = 0; j < palette.size(); j++) {
		//	float difp = colorDifferance(palette[pidex], palette[j]);
		//	//int p2 = int(uchar((0.3 *palette[j].val[2]) + (0.59 *palette[j].val[1]) + (0.11 * palette[j].val[0])));
		//	float p2 = (((0.2126 *palette[j].val[2]) + (0.7152 *palette[j].val[1]) + (0.0722 * palette[j].val[0])));	// return gray of average color	
		//	if (p1 >= p2) signp = true;
		//	else signp = false;
		//	float dif = abs(colorDif - difp);
		//	if (dif < mindif && signr == signp) {
		//		mindif = dif;
		//		index = j;
		//	}
		//}
		//if (index == -1) {
		//	index = pidex;
		//	mindif = 1000.0f;
		//	if (g1 >= g2)
		//	{
		//		for (int j = 0; j < palette.size(); j++) {
		//			float p2 = (((0.2126 *palette[j].val[2]) + (0.7152 *palette[j].val[1]) + (0.0722 * palette[j].val[0])));
		//			float difp = colorDifferance(palette[pidex], palette[j]);
		//			float dif = abs(colorDif - difp);
		//			if (p1 >= p2 && dif < mindif) {
		//				mindif = dif;
		//				index = j;
		//			}
		//		}
		//	}
		//	if (g1 < g2)
		//	{
		//		for (int j = 0; j < palette.size(); j++) {
		//			float p2 = (((0.2126 *palette[j].val[2]) + (0.7152 *palette[j].val[1]) + (0.0722 * palette[j].val[0])));
		//			float difp = colorDifferance(palette[pidex], palette[j]);
		//			float dif = abs(colorDif - difp);
		//			if (p1 <= p2 && dif < mindif) {
		//				mindif = dif;
		//				index = j;
		//			}
		//		}
		//	}
		//}
		/*********Design 3************/
		//if (g1 >= g2) signr = true;//positive
		//else signr = false;//negative
		//for (int j = 0; j < palette.size(); j++) {
		//	float difp = colorDifferance(palette[pidex], palette[j]);
		//	//int p2 = int(uchar((0.3 *palette[j].val[2]) + (0.59 *palette[j].val[1]) + (0.11 * palette[j].val[0])));
		//	float p2 = (((0.2126 *palette[j].val[2]) + (0.7152 *palette[j].val[1]) + (0.0722 * palette[j].val[0])));	// return gray of average color	
		//	if (p1 >= p2) signp = true;
		//	else signp = false;
		//	float dif = abs(colorDif - difp);
		//	if (dif < mindif && signr == signp) {
		//		mindif = dif;
		//		index = j;
		//	}
		//}
		//if (index == -1) { // find a new pair
		//	index = pidex;
		//	mindif = 1000.0f;
		//	std::pair<int, int> indexes;
		//	float colorDif = lookUp_table[edge->weight];
		//	for (int i = 0; i < table.rows; i++) {
		//		for (int j = 0; j < table.cols; j++) {
		//			float dif = abs(colorDif - table.at<float>(i, j));
		//			if (dif < mindif) {
		//				mindif = dif;
		//				indexes.first = i;
		//				indexes.second = j;
		//			}
		//		}
		//	}
		//	float p1 = (((0.2126 *palette[indexes.first].val[2]) + (0.7152 *palette[indexes.first].val[1]) +(0.0722 * palette[indexes.first].val[0])));	// return luminance
		//	float p2 = (((0.2126 *palette[indexes.second].val[2]) + (0.7152 *palette[indexes.second].val[1]) +(0.0722 * palette[indexes.second].val[0])));	// return luminance	
		//	if (p1 >= p2) signp = true;
		//	else signp = false;
		//	if (signr == signp) {
		//		returnRegion(edge->first)->setNewColor(palette[indexes.first]);
		//		returnRegion(edge->second)->setNewColor(palette[indexes.second]);
		//	}
		//	else {
		//		returnRegion(edge->first)->setNewColor(palette[indexes.second]);
		//		returnRegion(edge->second)->setNewColor(palette[indexes.first]);
		//	}
		//	returnRegion(edge->first)->colored = true;
		//	returnRegion(edge->second)->colored = true;		
		//}

		/*********Design 4************/
		if (g1 >= g2) signr = true;//positive
		else signr = false;//negative
		mindif = 1000.0f;
		for (int j= 0; j < palette.size(); j++) {
			float difp = colorDifferance(palette[pidex], palette[j]);
			float p2 = (((0.2126 *palette[j].val[2]) + (0.7152 *palette[j].val[1]) + (0.0722 * palette[j].val[0])));	
			if (p1 >= p2) signp = true;
			else signp = false;
			float dif = abs(colorDif - difp);
			if (signr == signp ) {
				if (dif <= mindif) {
					mindif = dif;
					index = j;
				}
		    }
		}	
		if (index == -1)
		{
			float mindif = 1000.0f;
			for (int j = 0; j < palette.size(); j++) {
				float difp = colorDifferance(palette[pidex], palette[j]);
				float dif = abs(colorDif - difp);
				if (dif < mindif) {
					mindif = dif;
					index = j;
				}
			}
		}
		return index;
	}

	void Region_Graph::ClosestPaletteIndexes(RegionEdge * edge, cv::Mat &table,	std::vector<cv::Vec3b> palette, 
		std::vector<int> lookUp_table, std::vector<threeObj> &pairs) {
		float mindif = 1000.0f;
		std::pair<int, int> index;
		float colorDif = lookUp_table[edge->weight];
		for (int i = 0; i < table.rows; i++) {
			for (int j = 0; j < table.cols; j++) {
				float dif = abs(colorDif - table.at<float>(i, j));
				//std::cout << "diff:" << dif << " , edge weight: " << colorDif << " , table: " << table.at<float>(i, j) << std::endl;
				if (dif < mindif) {
					mindif = dif;
					index.first = i;
					index.second = j;
					threeObj obj ;
					obj.first = palette[i];
					obj.second = palette[j];
					obj.dif = dif;
					pairs.push_back(obj);
				}
			}
		}
		sort(pairs.begin(), pairs.end(),smallestDiff);
		//std::cout << "1 pair: " << pairs[0].dif << std::endl;
		//std::cout << "2 pair: " << pairs[1].dif << std::endl;
		//std::cout << "3 pair: " << pairs[2].dif << std::endl;
		//std::cout << "4 pair: " << pairs[3].dif << std::endl;
	}

	void Region_Graph::RecoloringByColorPalette(cv::Mat &screen) {
		std::vector<cv::Vec3b> palette;
		std::pair<int, int> index;
		struct CompareEdgeReg
		{
			bool operator()(RegionEdge *n1, RegionEdge *n2) const
			{
				return (n1->weight) < (n2->weight); // < larger to smaller values - > small to large priority
			}
		};
		std::priority_queue< RegionEdge*, std::vector< RegionEdge*>, CompareEdgeReg > Reg_edges;
		for (int i = 0; i < graph_regions.size(); i++) {
			graph_regions[i]->colored = false;
			graph_regions[i]->setNewColor(cv::Vec3b(0,0,0));
			for (int e = 0; e < graph_regions[i]->neighbours.size(); e++) {
				//int nid = graph_regions[i]->edges[e]->second;
				Reg_edges.push(graph_regions[i]->edges[e]);
			}
		}
		//  colorfull palette
		palette.push_back(cv::Vec3b(207, 116, 158));   //  1
		palette.push_back(cv::Vec3b(132, 52, 87));  //  2
		palette.push_back(cv::Vec3b(84, 32, 41));   //  3
		palette.push_back(cv::Vec3b(112, 54, 162));   //  4
		palette.push_back(cv::Vec3b(117, 68, 237));   //  5
		palette.push_back(cv::Vec3b(142, 210, 254));   //  6
		palette.push_back(cv::Vec3b(116, 170, 255));  //  7
		palette.push_back(cv::Vec3b(89, 127, 250));   //  8
		palette.push_back(cv::Vec3b(75, 86, 227));   //  9
		palette.push_back(cv::Vec3b(75, 38, 169));   //  10
		
		cv::Mat table = cv::Mat::zeros(palette.size(), palette.size(), CV_32F);
		//*** table keeps the pairwise differences between color palettes
		for (int i = 0; i < palette.size();i++) {
			for (int j = 0; j < palette.size(); j++) {
				ColorSpace::Rgb Rgb1(palette[i].val[2], palette[i].val[1],palette[i].val[0]);
				ColorSpace::Rgb Rgb2(palette[j].val[2], palette[j].val[1], palette[j].val[0]);
				double diff = ColorSpace::EuclideanComparison::Compare(&Rgb1,&Rgb2);
				
				table.at<float>(i, j) = float(diff) ;// deltaE(palette1, palette2);
			}
		}
		/***************palette histogram**************/
		////std::string path = "D:/Code/BW_details/results/" + std::to_string(160) + "/12003/";
		//ofstream table_p(path + "diff_colors_palette.csv");
		////cv::Mat table2 = cv::Mat::zeros(256, 256, CV_32F);
		////*** table keeps the pairwise differences between region colors
		//for (int i = 0; i < palette.size(); i++) {
		//	for (int j = 0; j <  palette.size(); j++) {
		//		//hist[table2.at<int>(i, j)]++;
		//		std::cout << table.at<float>(i, j) << std::endl;
		//		if (table_p.is_open())
		//		{
		//			table_p << (int) table.at<float>(i, j) << std::endl;
		//		}
		//	}
		//}
		//table_p.close();
		/**********Histogram matching**********************/
		std::vector<int> lookUp_table;
		int binsize = 20;
		HistogramMatching(edgeweights, table, lookUp_table, binsize); // between color differences of graph_regions and table 
		//**** assign the closest colors to regions 
		while (!Reg_edges.empty()) {
			RegionEdge *edge = Reg_edges.top();
			Reg_edges.pop();
			//Region *one = new Region();  //one = returnRegion(edge->first);
			//Region *two = new Region();  two = returnRegion(edge->second);
			/***********case 2: both regions are colored ******************/
			if (returnRegion(edge->second)->colored == true && returnRegion(edge->first)->colored == true) {
				continue;
			}
			/***********Case 1 : both regions are not colored***************/
			if (returnRegion(edge->second)->colored == false && returnRegion(edge->first)->colored == false) {
				//***** return closest to palette**** 
				index = ReturnClosestPalette(edge, table, palette, lookUp_table);
				cv::Vec3b pcol1 = palette[index.first];
				cv::Vec3b pcol2 = palette[index.second];
				//****regions colors****
				cv::Vec3b col1 = returnRegion(edge->first)->calcAveColor();// calcAveColor(/*image*/);
				cv::Vec3b col2 = returnRegion(edge->second)->aveColor;// calccalcAveColor()(/*image*/);
				/***********************compute differences in Lab******************/
				ColorSpace::Rgb Rgb1(col1.val[2], col1.val[1], col1.val[0]);
				ColorSpace::Rgb Rgb2(col2.val[2], col2.val[1], col2.val[0]);
				ColorSpace::Rgb Rgbp1(pcol1.val[2], pcol1.val[1], pcol1.val[0]);
				ColorSpace::Rgb Rgbp2(pcol2.val[2], pcol2.val[1], pcol2.val[0]);
				double diffLab1 = ColorSpace::EuclideanComparison::Compare(&Rgbp1, &Rgb1);
				double diffLab2 = ColorSpace::EuclideanComparison::Compare(&Rgbp1, &Rgb2);
				int p1 = int(uchar((0.3 *pcol1.val[2]) + (0.59 *pcol1.val[1]) + (0.11 * pcol1.val[0])));	// return gray of average color	
				int p2 = int(uchar((0.3 *pcol2.val[2]) + (0.59 *pcol2.val[1]) + (0.11 * pcol2.val[0])));	// return gray of average color	
				int r1 = int(uchar((0.3 *col1.val[2]) + (0.59 *col1.val[1]) + (0.11 * col1.val[0])));	// return gray of average color	
				int r2 = int(uchar((0.3 *col2.val[2]) + (0.59 *col2.val[1]) + (0.11 * col2.val[0])));	// return gray of average color	
				/*****************************************************************/
				//bool signr = false; bool signp = false;
				//if (r1 > r2) signr = true;//positive
				//else signr = false;//negative
				//if (p1 > p2) signp = true;//positive
				//else signp = false;//negative
				//if (signp == signr) {
				//	for (int i = 0; i < returnRegion(edge->first)->regionPix.size(); i++)
				//	{
				//		int x = returnRegion(edge->first)->regionPix[i]->getPos().x;
				//		int y = returnRegion(edge->first)->regionPix[i]->getPos().y;
				//		screen.at<cv::Vec3b>(y, x) = pcol1;
				//	}
				//	returnRegion(edge->first)->colored = true;
				//	returnRegion(edge->first)->setNewColor( pcol1);
				//	for (int i = 0; i < returnRegion(edge->second)->regionPix.size(); i++)
				//	{
				//		int x = returnRegion(edge->second)->regionPix[i]->getPos().x;
				//		int y = returnRegion(edge->second)->regionPix[i]->getPos().y;
				//		screen.at<cv::Vec3b>(y, x) = pcol2;
				//	}
				//	returnRegion(edge->second)->colored = true;
				//	returnRegion(edge->second)->setNewColor(pcol2);
				//}
				//else {
				//	for (int i = 0; i < returnRegion(edge->first)->regionPix.size(); i++)
				//	{
				//		int x = returnRegion(edge->first)->regionPix[i]->getPos().x;
				//		int y = returnRegion(edge->first)->regionPix[i]->getPos().y;
				//		screen.at<cv::Vec3b>(y, x) = pcol2;
				//	}
				//	returnRegion(edge->first)->colored = true;
				//	returnRegion(edge->first)->setNewColor( pcol2);
				//	for (int i = 0; i < returnRegion(edge->second)->regionPix.size(); i++)
				//	{
				//		int x = returnRegion(edge->second)->regionPix[i]->getPos().x;
				//		int y = returnRegion(edge->second)->regionPix[i]->getPos().y;
				//		screen.at<cv::Vec3b>(y, x) = pcol1;
				//	}
				//	returnRegion(edge->second)->colored = true;
				//	returnRegion(edge->second)->setNewColor( pcol1);
				//}

				if (diffLab1 > diffLab2 ) {
					
					for (int i = 0; i < returnRegion(edge->second)->regionPix.size(); i++)
					{
						int x = returnRegion(edge->second)->regionPix[i]->getPos().x;
						int y = returnRegion(edge->second)->regionPix[i]->getPos().y;
						screen.at<cv::Vec3b>(y, x) = pcol1;
					}
					returnRegion(edge->second)->colored = true;
					returnRegion(edge->second)->setNewColor( pcol1);
					for (int i = 0; i < returnRegion(edge->first)->regionPix.size(); i++)
					{
						int x = returnRegion(edge->first)->regionPix[i]->getPos().x;
						int y = returnRegion(edge->first)->regionPix[i]->getPos().y;
						screen.at<cv::Vec3b>(y, x) = pcol2;
					}
					returnRegion(edge->first)->colored = true;
					returnRegion(edge->first)->setNewColor(pcol2);
				}
				else {
					for (int i = 0; i < returnRegion(edge->second)->regionPix.size(); i++)
					{
						int x = returnRegion(edge->second)->regionPix[i]->getPos().x;
						int y = returnRegion(edge->second)->regionPix[i]->getPos().y;
						screen.at<cv::Vec3b>(y, x) = pcol2;
					}
					returnRegion(edge->second)->colored = true;
					returnRegion(edge->second)->setNewColor(pcol2);
					for (int i = 0; i < returnRegion(edge->first)->regionPix.size(); i++)
					{
						int x = returnRegion(edge->first)->regionPix[i]->getPos().x;
						int y = returnRegion(edge->first)->regionPix[i]->getPos().y;
						screen.at<cv::Vec3b>(y, x) = pcol1;
					}
					returnRegion(edge->first)->colored = true;
					returnRegion(edge->first)->setNewColor( pcol1);
				}
			}
			
			if (returnRegion(edge->second)->colored == false || returnRegion(edge->first)->colored == false) {
				////***** return closest to palette**** 
				//int indx = ReturnClosestPaletteIndex(edge, palette, lookUp_table);
				//cv::Vec3b pcol ;//= palette[indx];
				
				if (returnRegion(edge->second)->colored == false && returnRegion(edge->first)->colored == true) {
					////***** return closest to palette**** 
					int indx = ReturnClosestPaletteIndex2(edge, palette,table, lookUp_table, returnRegion(edge->first)->getNewColor());
					cv::Vec3b pcol = palette[indx];

					for (int i = 0; i < returnRegion(edge->second)->regionPix.size(); i++)
					{
						int x = returnRegion(edge->second)->regionPix[i]->getPos().x;
						int y = returnRegion(edge->second)->regionPix[i]->getPos().y;
						screen.at<cv::Vec3b>(y, x) = pcol;
					}
					returnRegion(edge->second)->colored = true;
					returnRegion(edge->second)->setNewColor(pcol);
				}
				
				if (returnRegion(edge->first)->colored == false && returnRegion(edge->second)->colored == true) {
					////***** return closest to palette**** 
					int indx = ReturnClosestPaletteIndex2(edge, palette,table, lookUp_table, returnRegion(edge->second)->getNewColor());
					cv::Vec3b pcol = palette[indx];

					for (int i = 0; i < returnRegion(edge->first)->regionPix.size(); i++)
					{
						int x = returnRegion(edge->first)->regionPix[i]->getPos().x;
						int y = returnRegion(edge->first)->regionPix[i]->getPos().y;
						screen.at<cv::Vec3b>(y, x) = pcol;
					}
					returnRegion(edge->first)->colored = true;
					returnRegion(edge->first)->setNewColor( pcol);
				}
			}
			
		}
	}
	
	void Region_Graph::reducePalette(std::vector<int> &indexes,std::vector<cv::Vec3b> &palette, std::vector<cv::Vec3b> &newpalet) {
		std::list<cv::Vec3b> newpal;
		std::list<cv::Vec3b>::iterator it1;
		for (int i = 0; i < newpalet.size(); i++) {
			if(i != indexes[0] && i != indexes[1])
				newpal.push_back(newpalet[i]);
		}
		//list<cv::Vec3b>::iterator iter = newpal.begin();

		//list<cv::Vec3b>::iterator t = newpal.begin();
		//list<cv::Vec3b>::iterator t1 = newpal.begin();
		//std::advance(t, indexes[0]);
		//newpal.erase(t);
		//t = newpal.begin();
		//std::advance(t, indexes[1]);
		//newpal.erase(t);

		std::cout << "palette size: "<< newpal.size() << std::endl;

		std::cout << "palette size after reduction: " << newpal.size() << std::endl;
		int size = newpal.size();

		newpalet.clear(); newpalet.reserve(size);
		for (it1 = newpal.begin(); it1 != newpal.end(); ++it1) {
			newpalet.push_back(*it1);
		}
		std::cout << "palette size after reduction: " << newpalet.size() << std::endl;

		indexes.resize(0);
	}
	
	void Region_Graph::BFS(Region *curr,cv::Mat &etretat)
	{
		// Mark the current node as visited and print it 

		curr->visited = true;

		//allPix[v]->setCCId(label);

		for (int j = 0; j < curr->children.size(); j++) {
			int nid = curr->children[j]->getId();
			//int regid = allPix[v]->edges[j]->getTheOtherSide(allPix[v])->getRegId();
			if (curr->getId() == graph_regions[nid]->parent->getId())
			{
				cv::Point p1 = retunMidPoint(graph_regions[curr->getId()]); cv::Point p2 = retunMidPoint(graph_regions[nid]);
				cv::line(etretat, cv::Point(4 * p1.x, 4 * p1.y), cv::Point(4 * p2.x, 4 * p2.y), cv::Scalar(255, 255, 20), 2, 8);
				
				if (graph_regions[nid]->visited == false)
					BFS(graph_regions[nid],etretat);
			}
		}
	}
	
	/*******Recoloring with lookup table******/
	void Region_Graph::RecoloringByBottleneck(cv::Mat &screen, std::vector<cv::Vec3b> &palette, std::string palettename, 
		std::vector<int> &indexes, std::vector < std::pair<int, int>> &values_p, int &mid_reg) {

		/**************SVG file ****************/
		//Dimensions dimensions(screen.cols, screen.rows);
		//Document doc(path + "svg_regions_"+ palettename +".svg", Layout(dimensions, Layout::TopLeft));

		cv::Mat screen2 = screen.clone();
		cv::Mat contourImage = cv::Mat::zeros(image.size(), CV_8UC3);
		std::vector<int> regsizes;
		for (int i = 0; i < graph_regions.size(); i++) {
			if (graph_regions[i]->regionPix.size() > 0) {
				regsizes.push_back(graph_regions[i]->regionPix.size());
				cv::Vec3b color = graph_regions[i]->calcAveColor(); //graph_regions[i]->aveColor;// 
				
				for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
				{
					int x = graph_regions[i]->regionPix[j]->getPos().x;
					int y = graph_regions[i]->regionPix[j]->getPos().y;
					contourImage.at<cv::Vec3b>(y, x) = color;
				}
			}
		}
		cv::imwrite(path + "AveColorRegs.png", contourImage);

		/**********Median region size********/
		double medsize = CalcMHWScore(regsizes);
		double avesize = CalcAverage(regsizes);
		std::cout << "median region size: "<< medsize << ", number of regions above 0: "<< regsizes.size() << std::endl;
		clock_t tt1, tt2, tt3,tt4, tt5,tt6;
		strippednode  sn, newnode;
		int s;
		
		struct CompareRegion2
		{
			bool operator()(strippednode n1, strippednode n2) const
			{
				return (n1.distance) < (n2.distance); // < larger to smaller values - > small to large priority
			}
		};
		struct CompareRegion
		{
			bool operator()(Region *n1, Region *n2) const
			{

				return (n1->dist) < (n2->dist); // < larger to smaller values - > small to large priority
			}
		};
		
		std::priority_queue< Region*, std::vector< Region*>, CompareRegion > Regions;
		Regions.empty();
		std::priority_queue< strippednode, std::vector< strippednode>, CompareRegion2 > Regions2;
		Regions2.empty();
		//std::vector<cv::Vec3b> palette;
		std::pair<int, int> index;

		//int pid = (image.cols / 2) + (image.rows / 2) * image.cols;// RegionPixels.size();  
		   //int pid = (image.rows / 4) + (image.cols / 2) * image.cols; // top
		   //int pid = (image.rows * 3/4) + (image.cols / 2) * image.cols; //down
		//int mid_reg = RegionPixels[pid]->getRegId();
		std::vector< RegionEdge*> Reg_edges;
		for (int i = 0; i < graph_regions.size(); i++) {
			graph_regions[i]->children.clear();
			graph_regions[i]->setNewColor(cv::Vec3b(-1, -1, -1));
			graph_regions[i]->dist = FLT_MIN;
			graph_regions[i]->parent = NULL;
			graph_regions[i]->colored = false;
			graph_regions[i]->visited = false;
			if (i == mid_reg) { // push to q
				graph_regions[i]->dist = FLT_MAX;
				graph_regions[i]->visited = true;
				graph_regions[i]->parent = graph_regions[i];
				sn.id = graph_regions[i]->getId();
				sn.distance = FLT_MAX;
				Regions2.push(sn);
				Regions.push(graph_regions[i]);
			}
		}

		cv::Mat table = cv::Mat::zeros(palette.size(), palette.size(), CV_32F);
		//*** table keeps the pairwise differences between color palettes
		for (int i = 0; i < palette.size(); i++) {
			for (int j = 0; j < palette.size(); j++) {

				ColorSpace::Rgb Rgb1(palette[i].val[2], palette[i].val[1], palette[i].val[0]);
				ColorSpace::Rgb Rgb2(palette[j].val[2], palette[j].val[1], palette[j].val[0]);
				double diff = ColorSpace::EuclideanComparison::Compare(&Rgb1, &Rgb2);

				table.at<float>(i, j) = float(diff);
			}
		}
		/***************palette histogram**************/
		ofstream table_p(path + "diff_colors_palette.csv");
		//*** table keeps the pairwise differences between region colors
		for (int i = 0; i < palette.size(); i++) {
			for (int j = 0; j <  palette.size(); j++) {
				if (table_p.is_open())
				{
					table_p << (int) table.at<float>(i, j) << std::endl;
				}
			}
		}
		table_p.close();
		/*******************/
		tt5 = clock();
		std::vector<int> lookUp_table; 
		int binsize = palette.size()*palette.size();

		HistogramMatching(edgeweights, table, lookUp_table,binsize);// between color differences of graph_regions and table 
		float midweight = getMidWeight(edgeweights, lookUp_table);
		tt6 = clock();
		/******************/
		ofstream table_in(path + "diff_colors_Before_After_Mapping.csv");

		/******starting region in the middle of the image*****/
		Region* source = NULL;// = graph_regions[0];// NULL
		//	graph_regions[mid_reg]->dist = 0;
		//	graph_regions[mid_reg]->parent = graph_regions[mid_reg];
		//	Regions.push(graph_regions[mid_reg]);
		int count = 0;
		int b = 0; int n = 0;
		cv::Mat draw = cv::Mat::zeros(screen.size(), screen.type());
		std::vector<float> ratio_cp; //ratio of child to parent size
		tt1 = clock();

		while (!Regions2.empty()) {
			sn = Regions2.top();
			Regions2.pop();

			s = sn.id;
			Region* cur = returnRegion(s);
			if (count == 0)
			{
				source = cur;
			}
			if (cur->colored == false)
			{
				/*****recolor*********/
				if (cur->parent != NULL)
					(cur->parent)->children.push_back(cur);
				//int indx = ReturnClosestPaletteIndex1(source, palette, table, lookUp_table, cur);
				cur->colored = true;
				//cur->setNewColor(palette[indx]);

				float colordist = cur->dist;

				int nsize = cur->edges.size();
				for (int e = 0; e < nsize; e++) {
					int nid = cur->edges[e]->second; // neighbour region
					float nWeight0 = cur->edges[e]->weight;// lookUp_table[cur->edges[e]->weight  /*+ colorDiff(source, up)*/];// /*+ colorDiff(source, up)*/;
					Region *neighbReg = returnRegion(nid);
					//float sizefactor = returnRegion(nid)->regionPix.size() / (returnRegion(nid)->regionPix.size() + up->regionPix.size());
					float nWeight = abs(nWeight0 - midweight) + (float)(neighbReg->regionPix.size() / avesize) *  abs(nWeight0 - midweight); // (1+b) Delta
					//float nWeight = abs(nWeight0 - midweight) * (min(1.0, (neighbReg->regionPix.size() / avesize))); // min(1,b)
					//float nWeight = abs(nWeight0 - midweight);
					float alt = min(colordist, nWeight); //max(neighbReg->dist, min(colordist, up->edges[e]->weight));
					/******************/
					if (alt > neighbReg->dist)
					{
						if (neighbReg->colored == false) {
							neighbReg->dist = alt;
							neighbReg->parent = cur;
							
							newnode.id = nid;
							newnode.distance = alt;
							Regions2.push(newnode);
						}
					}
					else {
					}
				}
			}
			count++;
		}

		table_in.close();

		//std::cout << "num of times recoloring couldn't find color based on luminance: " << counter << std::endl;
		counter = 0;
		tt2 = clock();
		/*********fill the regions that are not colored***************/
		int c = 0; int d = 0; int v = 0;
		for (int i = 0; i < graph_regions.size(); i++) {
			graph_regions[i]->colored = false;
		}
		std::cout << "num regions not colored initially: " << c << " , no parent set: " << d << " , not visited: " << v << std::endl;
		/**************** render the colors******************************/
		tt3 = clock();
		/***********set the colors *******/	
		AssignColors(source, palette, table, lookUp_table,binsize, source);
		
		tt4 = clock();
		/*********Mapped colors cdf**********/
		//double miPaleDif, maxPaleDif;
		//cv::minMaxLoc(table, &miPaleDif, &maxPaleDif);
		//float mmm = *max_element(treeEdges.begin(), treeEdges.end());
		//std::vector<int> Mappedfreq(500);
		//for (int i = 0; i < Mappedfreq.size(); i++) //initialize frequency array
		//	Mappedfreq[i] = 0;
		////compute frequencies
		//std::vector<float> mapArr(treeEdges.size());
		//for (int i = 0; i < treeEdges.size();i++) {
		//	mapArr.push_back(lookUp_table[(int)floor((treeEdges[i] / mmm )* (maxPaleDif))]);
		//}
		//for (int i = 0; i < mapArr.size(); i++)
		//	Mappedfreq[(int)mapArr[i]]++;
		////calculate cdf after mapping
		//float nn = *max_element(mapArr.begin(), mapArr.end());
		//std::vector<float> normalized_Mappedtcdf;
		//calculate_cdf(nn, Mappedfreq, normalized_Mappedtcdf);
		//ofstream table_Mappedcdf(path + "normalized_Mappedtcdf_colorful_Euclid.csv");
		//if (table_Mappedcdf.is_open() )
		//{
		//	for (int i = 0; i <normalized_Mappedtcdf.size(); i++) {
		//		table_Mappedcdf << normalized_rcdf[i] << "," << normalized_Mappedtcdf[i] << std::endl;
		//	}
		//}
		//table_Mappedcdf.close();
		/********************/
		int ind;
		std::vector<int> bin;
		for (int p = 0; p < palette.size(); p++)
			bin.push_back(0);
		for (int i = 0; i < graph_regions.size(); i++) {
			cv::Vec3b pcol = graph_regions[i]->getNewColor();
			for (int p = 0; p < palette.size(); p++) {
				if (palette[p] == pcol)
					ind = p;
			}
			for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
			{
				int x = graph_regions[i]->regionPix[j]->getPos().x;
				int y = graph_regions[i]->regionPix[j]->getPos().y;
				screen.at<cv::Vec3b>(y, x) = pcol;
				bin[ind]++;
			}

			/************ SVG ****************/
			//for (int j = 0; j < graph_regions[i]->ListOfcontourPix.size(); j++) {
			//	if (graph_regions[i]->ListOfcontourPix[j].size() > 0) {
			//		Polygon region(Color(pcol.val[2], pcol.val[1], pcol.val[0]), Stroke(.5, Color(0, 0, 0)));
			//		for (int p = 0; p < graph_regions[i]->ListOfcontourPix[j].size(); p++) {
			//			svg::Point svgp((double)graph_regions[i]->ListOfcontourPix[j][p]->getPos().x, (double)graph_regions[i]->ListOfcontourPix[j][p]->getPos().y);
			//			region.operator<< (svgp);
			//		}
			//		doc << region;
			//	}
			//}
			//if (graph_regions[i]->contourPix.size() > 0) {
			//	Polygon region(Color(pcol.val[2], pcol.val[1], pcol.val[0]), Stroke(.5, Color(0, 0, 0)));
			//	for (int p = 0; p < graph_regions[i]->contourPix.size(); p++) {
			//		svg::Point svgp((double)graph_regions[i]->contourPix[p]->getPos().x, (double)graph_regions[i]->contourPix[p]->getPos().y);
			//		region.operator<< (svgp);
			//	}
			//	doc << region;
			//}
			
		}
		//doc.save(); // save svg file

		float difft1((float)tt2 - (float)tt1);
		float difft2((float)tt4 - (float)tt3);
		float difft3((float)tt6 - (float)tt5);
		std::cout << difft1 / CLOCKS_PER_SEC << "  seconds to create the tree  " << std::endl;
		std::cout << difft2 / CLOCKS_PER_SEC << "  seconds to assign colors  " << std::endl;
		std::cout << difft3 / CLOCKS_PER_SEC << "  seconds to apply HM & calc midweight  " << std::endl;
		/********num of pixels per color in the palette******/
		ofstream lumHist(path + "palette_distribution.csv");
		if (lumHist.is_open())
		{
			for (int i = 0; i <bin.size(); i++) {
				lumHist << (int)bin[i] << std::endl;
			}
		}
		lumHist.close();		
		/********Draw tree *********/
		//source->visited = false;
		//int countt = 0;
		//draw = screen.clone();
		//DrawTree(source, draw, countt);
		//imwrite(path + "RecoloredBN_tree_tableW-midw.png", draw);
		//for (int i = 0; i < graph_regions.size(); i++) {
		//	if (graph_regions[i]->visited == false) {
		//		BFS(source, etretat);
		//	}
		//}

		/****************************************/
		cv::Mat lines;
		cv::resize(screen, lines, cv::Size(), 4.0, 4.0);
		DrawChildren(source,lines);
		//cv::circle(etretat,8* retunMidPoint(source), 2, cv::Scalar(0, 200, 0), 2, 8);
		//imwrite(path + "RecoloredBN_tree_" + palettename +"_EuclideanComparison.png", lines);

		/***********add bar color on the bottom of the result images **************/
		DrawBarCOLOR(bin, palette, screen,indexes,values_p);
		lookUp_table.clear();
	}
	
	/*******Recoloring without lookup table******/
	void Region_Graph::RecoloringByBottleneck_noLUT(cv::Mat &screen, std::vector<cv::Vec3b> &palette, std::string palettename, 
		std::vector<int> &indexes , std::vector < std::pair<int, int>> &values_p, int & mid_reg) {
		/**************SVG file ****************/
		//Dimensions dimensions(screen.cols, screen.rows);
		//Document doc(path + "svg_regions_"+ palettename +".svg", Layout(dimensions, Layout::TopLeft));
		//cv::Mat etretat = cv::imread(path + "/RecoloredBN_colorful_Palette_Euclid_lut_etretat1024_r-mid-tableW-midw_4x.png", 1);
		cv::Mat screen2 = screen.clone();
		cv::Mat contourImage = cv::Mat::zeros(image.size(), CV_8UC3);
		std::vector<int> regsizes;
		for (int i = 0; i < graph_regions.size(); i++) {
			cv::Mat draw5 = cv::Mat::zeros(image.size(), CV_8UC3);
			if (graph_regions[i]->regionPix.size() > 10)
				regsizes.push_back(graph_regions[i]->regionPix.size());
			cv::Vec3b color = graph_regions[i]->calcAveColor();// calcAveColor();
			for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
			{
				int x = graph_regions[i]->regionPix[j]->getPos().x;
				int y = graph_regions[i]->regionPix[j]->getPos().y;
				contourImage.at<cv::Vec3b>(y, x) = color;
				draw5.at<cv::Vec3b>(y, x) = color;
			}
			cv::Mat gbbp = cv::Mat::zeros(screen.size(), CV_8UC1);
			cvtColor(draw5, gbbp, cv::COLOR_BGR2GRAY);
			std::vector<cv::Vec4i> hierarchy;
			std::vector< std::vector<cv::Point>> contours;
			findContours(gbbp, contours, hierarchy, CV_RETR_EXTERNAL /*CV_RETR_CCOMP*/, CV_CHAIN_APPROX_NONE, cv::Point(0, 0));
			for (int j = 0; j < contours.size(); j++) {
				if(contours[j].size() >10)
					drawContours(contourImage, contours, j, cvScalar(0), 1, 8);// , hierarchy, 0);
			}
		}
		//cv::imwrite(path + "AveColorRegs.png", contourImage);
		/**************************************/
		cv::imwrite(path + "contours.png", contourImage);
		/**********MEdian region size********/
		double medsize = CalcMHWScore(regsizes);
		double avesize = CalcAverage(regsizes);
		std::cout << "median region size: " << medsize << ", number of regions above 10: " << regsizes.size() << std::endl;
		clock_t tt1, tt2, tt3, tt4, tt5,tt6;
		strippednode  sn, newnode;
		int s;
		struct CompareRegion2
		{
			bool operator()(strippednode n1, strippednode n2) const
			{
				return (n1.distance) < (n2.distance); // < larger to smaller values - > small to large priority
			}
		};
		struct CompareRegion
		{
			bool operator()(Region *n1, Region *n2) const
			{

				return (n1->dist) < (n2->dist); // < larger to smaller values - > small to large priority
			}
		};
		std::priority_queue< Region*, std::vector< Region*>, CompareRegion > Regions;
		Regions.empty();
		std::priority_queue< strippednode, std::vector< strippednode>, CompareRegion2 > Regions2;
		Regions2.empty();
		//std::vector<cv::Vec3b> palette;
		//palette.push_back(palette0[values_p[values_p.size() - 1].second]); // two largets palette from previouse palette
		//palette.push_back(palette0[values_p[values_p.size() - 2].second]);

		std::pair<int, int> index;

		//int pid = RegionPixels.size();// (image.rows / 2) + (image.cols / 2) * image.cols; ;// 
		//int mid_reg = RegionPixels[pid/2]->getRegId();
		std::vector< RegionEdge*> Reg_edges;
		for (int i = 0; i < graph_regions.size(); i++) {
			graph_regions[i]->children.clear();
			graph_regions[i]->setNewColor(cv::Vec3b(-1, -1, -1));
			graph_regions[i]->dist = FLT_MIN;
			graph_regions[i]->parent = NULL;
			graph_regions[i]->colored = false;
			graph_regions[i]->visited = false;
			if (i == mid_reg) { // push to q
				graph_regions[i]->dist = FLT_MAX;
				graph_regions[i]->visited = true;
				graph_regions[i]->parent = graph_regions[i];
				sn.id = graph_regions[i]->getId();
				sn.distance = FLT_MAX;
				Regions2.push(sn);
				Regions.push(graph_regions[i]);
			}
		}

		cv::Mat table = cv::Mat::zeros(palette.size(), palette.size(), CV_32F);
		//*** table keeps the pairwise differences between color palettes
		for (int i = 0; i < palette.size(); i++) {
			for (int j = 0; j < palette.size(); j++) {

				ColorSpace::Rgb Rgb1(palette[i].val[2], palette[i].val[1], palette[i].val[0]);
				ColorSpace::Rgb Rgb2(palette[j].val[2], palette[j].val[1], palette[j].val[0]);
				double diff = ColorSpace::EuclideanComparison::Compare(&Rgb1, &Rgb2);

				table.at<float>(i, j) = float(diff);
			}
		}
		/*******midweight************/
		tt5 = clock();
		std::vector<int> lookUp_table;
		lookUp_table.empty();
		float midweight = getMidWeight(edgeweights, lookUp_table);
		lookUp_table.clear();
		int binsize = 20;
	//	HistogramMatching(edgeweights, table, lookUp_table, binsize);
		//// between color differences of graph_regions and table 
		//float midweight = getMidWeight(edgeweights, lookUp_table);
		tt6 = clock();
		/******************/
		/******starting region in the middle of the image*****/
		Region* source = NULL;// = graph_regions[0];// NULL

		int count = 0;
		int b = 0; int n = 0;
		cv::Mat draw = cv::Mat::zeros(screen.size(), screen.type());
		tt1 = clock();


		while (!Regions2.empty()) {
			sn = Regions2.top();
			Regions2.pop();

			s = sn.id;
			Region* cur = returnRegion(s);
			if (count == 0)
			{
				source = cur;
			}
			if (cur->colored == false)
			{
				/*****recolor*********/
				if (cur->parent != NULL)
					(cur->parent)->children.push_back(cur);
				int indx = ReturnClosestPaletteIndex1(source, palette, table, lookUp_table,binsize, cur);				
				cur->setNewColor(palette[indx]);
				cur->colored = true;

				float colordist = cur->dist;

				int nsize = cur->edges.size();
				for (int e = 0; e < nsize; e++) {
					int nid = cur->edges[e]->second; // neighbour region
					float nWeight0 = cur->edges[e]->weight;
					Region *neighbReg = returnRegion(nid);
					//float sizefactor = returnRegion(nid)->regionPix.size() / (returnRegion(nid)->regionPix.size() + up->regionPix.size());
					float nWeight = abs(nWeight0 - midweight) + (float)(neighbReg->regionPix.size() / avesize) *  abs(nWeight0 - midweight); // (1+b) Delta
					//float nWeight = abs(nWeight0 - midweight) * (min(1.0, (neighbReg->regionPix.size() / avesize))); // min(1,b)
					//float nWeight = abs(nWeight0 - midweight);
					float alt = min(colordist, nWeight); //max(neighbReg->dist, min(colordist, up->edges[e]->weight));
					//std::cout << "neighbReg->dist: " << neighbReg->dist<< ", alt: " << alt << ", nWeight0: " << nWeight0 << std::endl;
					/******************/
					if (alt > neighbReg->dist)
					{
						if (neighbReg->colored == false) {
							neighbReg->dist =  alt;
							neighbReg->parent = cur;

							newnode.id = nid;
							newnode.distance = alt;
							Regions2.push(newnode);
						}
					}
					else {
					}
				}
			}
			count++;
		}

		tt2 = clock();

		/**************** render the colors******************************/
		tt3 = clock();
		//************set the colors *****	
		//HistogramMatching(edgeweights, table, lookUp_table);
		//AssignColors(source, palette, table, lookUp_table, source);

		tt4 = clock();

		int ind;
		std::vector<int> bin;
		for (int p = 0; p < palette.size(); p++)
			bin.push_back(0);
		for (int i = 0; i < graph_regions.size(); i++) {
			cv::Vec3b pcol = graph_regions[i]->getNewColor();
			for (int p = 0; p < palette.size(); p++) {
				if (palette[p] == pcol)
					ind = p;
			}
			for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
			{
				int x = graph_regions[i]->regionPix[j]->getPos().x;
				int y = graph_regions[i]->regionPix[j]->getPos().y;
				screen.at<cv::Vec3b>(y, x) = pcol;
				bin[ind]++;
			}

			/************ SVG ****************/
			//for (int j = 0; j < graph_regions[i]->ListOfcontourPix.size(); j++) {
			//	if (graph_regions[i]->ListOfcontourPix[j].size() > 0) {
			//		Polygon region(Color(pcol.val[2], pcol.val[1], pcol.val[0]), Stroke(.5, Color(0, 0, 0)));
			//		for (int p = 0; p < graph_regions[i]->ListOfcontourPix[j].size(); p++) {
			//			svg::Point svgp((double)graph_regions[i]->ListOfcontourPix[j][p]->getPos().x, (double)graph_regions[i]->ListOfcontourPix[j][p]->getPos().y);
			//			region.operator<< (svgp);
			//		}
			//		doc << region;
			//	}
			//}
			//if (graph_regions[i]->contourPix.size() > 0) {
			//	Polygon region(Color(pcol.val[2], pcol.val[1], pcol.val[0]), Stroke(.5, Color(0, 0, 0)));
			//	for (int p = 0; p < graph_regions[i]->contourPix.size(); p++) {
			//		svg::Point svgp((double)graph_regions[i]->contourPix[p]->getPos().x, (double)graph_regions[i]->contourPix[p]->getPos().y);
			//		region.operator<< (svgp);
			//	}
			//	doc << region;
			//}

		}
		//doc.save(); // save svg file
		
		float difft1((float)tt2 - (float)tt1);
		float difft2((float)tt4 - (float)tt3);
		float difft3((float)tt6 - (float)tt5);
		std::cout << difft1 / CLOCKS_PER_SEC << "  seconds to create the tree  " << std::endl;
		std::cout << difft2 / CLOCKS_PER_SEC << "  seconds to assign colors  " << std::endl;
		std::cout << difft3 / CLOCKS_PER_SEC << "  seconds to apply HM & calc midweight  " << std::endl;

		/********num of pixels per color in the palette******/
		ofstream lumHist(path + "palette_distribution.csv");
		if (lumHist.is_open())
		{
			for (int i = 0; i <bin.size(); i++) {
				lumHist << (int)bin[i] << std::endl;
			}
		}
		lumHist.close();

		/****************/
		cv::Mat gray = cv::Mat::zeros(image.rows, image.cols, CV_8UC3);
		cv::Mat lab = cv::Mat::zeros(image.rows, image.cols, CV_8UC3);
		cv::Mat fin_img = cv::Mat::zeros(image.rows, image.cols, CV_8UC3);
		cv::cvtColor(image, gray, cv::COLOR_BGR2Lab/*cv::COLOR_BGR2GRAY*/);
		cv::cvtColor(screen, lab, cv::COLOR_BGR2Lab/*cv::COLOR_BGR2GRAY*/);
		std::vector<cv::Mat> lab_channels;
		std::vector<cv::Mat> channels, mixchannels;
		cv::split(gray, lab_channels);// original image
		cv::split(lab, channels); // recolored
		cv::imwrite(path + "L_image.png", lab_channels[0]);
		cv::imwrite(path + "A_image.png", lab_channels[1]);
		cv::imwrite(path + "B_image.png", lab_channels[2]);
		/*******************Draw children*****************/
		cv::Mat lines;
		cv::resize(screen, lines, cv::Size(), 4.0, 4.0);
		DrawChildren(source, lines);
		//cv::circle(etretat,8* retunMidPoint(source), 2, cv::Scalar(0, 200, 0), 2, 8);
		//imwrite(path + "RecoloredBN_noLUT_tree_" + palettename + "_EuclideanComparison.png", lines);
		/*****merge L-channel and recolored*********/

		mixchannels.push_back(lab_channels[0]);
		mixchannels.push_back(channels[1]);
		mixchannels.push_back(channels[2]);
		cv::merge(mixchannels, fin_img);
		cv::cvtColor(fin_img, fin_img, cv::COLOR_Lab2BGR);
		cv::imwrite(path + "mereged.png", fin_img);
		/***********add bar color**************/
		DrawBarCOLOR(bin, palette, screen, indexes, values_p);
	}
	/*********RecoloringRandomColors***********/
	void Region_Graph::RecoloringRandomColors(cv::Mat& screen)
	{
		for (int i = 0; i < graph_regions.size(); i++) {
			if (graph_regions[i]->regionPix.size() > 0) {
				cv::Vec3b color = cv::Vec3b(rand() * 255, rand() * 255, rand() * 255);// graph_regions[i]->aveColor;// calcAveColor();
				for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
				{
					int x = graph_regions[i]->regionPix[j]->getPos().x;
					int y = graph_regions[i]->regionPix[j]->getPos().y;
					screen.at<cv::Vec3b>(y, x) = color;
				}
			}
		}
	}
	/*******Recoloring without path planning******/
	void Region_Graph::RecoloringPathFree(cv::Mat & screen, std::vector<cv::Vec3b> &palette, std::vector < std::pair<int, int>> &values_p)
	{
		std::vector<int> indexes;
		int ind;
		std::vector<int> bin;
		for (int p = 0; p < palette.size(); p++)
			bin.push_back(0);

		for (int i = 0; i < graph_regions.size(); i++) {
			float mindif = 1000.0; int index = -1;
			cv::Vec3b color = graph_regions[i]->calcAveColor();// calcAveColor();
			//int indx = ReturnClosestPaletteIndex(edge, palette, lookUp_table);
			for (int p = 0; p < palette.size();p++) {
				float diff = colorDifferance(color, palette[p]);															  
				if (diff < mindif) {
					mindif = diff;
					index = p;
				}
			}
			ind = index;
			cv::Vec3b pcol = palette[index];
			for (int j = 0; j < graph_regions[i]->regionPix.size();j++)
			{
				int x = graph_regions[i]->regionPix[j]->getPos().x;
				int y = graph_regions[i]->regionPix[j]->getPos().y;
				screen.at<cv::Vec3b>(y, x) = pcol;
				bin[ind]++;
			}
		}
		/***********add bar color**************/
		DrawBarCOLOR(bin, palette, screen, indexes, values_p);
	}
	
	/*******Radom color assignments ******/
	void Region_Graph::Recoloring_random(cv::Mat & screen, std::vector<cv::Vec3b> &palette, std::vector < std::pair<int, int>> &values_p)
	{
		std::vector<int> indexes;
		int ind;
		std::vector<int> bin;
		for (int p = 0; p < palette.size(); p++)
			bin.push_back(0);

		for (int i = 0; i < graph_regions.size(); i++) {
			int index = -1;
			cv::Vec3b color = graph_regions[i]->calcAveColor();// calcAveColor();
														 //int indx = ReturnClosestPaletteIndex(edge, palette, lookUp_table);
			for (int p = 0; p < palette.size(); p++) {
				index = rand() % palette.size();			
			}
			ind = index;
			cv::Vec3b pcol = palette[index];
			for (int j = 0; j < graph_regions[i]->regionPix.size(); j++)
			{
				int x = graph_regions[i]->regionPix[j]->getPos().x;
				int y = graph_regions[i]->regionPix[j]->getPos().y;
				screen.at<cv::Vec3b>(y, x) = pcol;
				bin[ind]++;
			}
		}
		/***********add bar color**************/
		DrawBarCOLOR(bin, palette, screen, indexes, values_p);
	}

	void Region_Graph::Histogram(std::vector<float> differences, std::string path) {
		//int hist[443];
		ofstream table_in(path + "diff_colors_cmc.csv");
		//cv::Mat table2 = cv::Mat::zeros(256,256, CV_32F);
		
		//*** table keeps the pairwise differences between region colors
		//for (int i = 0; i <256; i++) {
		//	for (int j = i; j < 256; j++) {			
		//		table2.at<int>(i, j) = sqrt(
		//			(in_colors[i].second.val[0] - in_colors[j].second.val[0])
		//			*(in_colors[i].second.val[0] - in_colors[j].second.val[0])+
		//			(in_colors[i].second.val[1] - in_colors[j].second.val[1])
		//			*(in_colors[i].second.val[1] - in_colors[j].second.val[1])+
		//			(in_colors[i].second.val[2] - in_colors[j].second.val[2])
		//			*(in_colors[i].second.val[2] - in_colors[j].second.val[2]));
		//		//hist[table2.at<int>(i, j)]++;
		//		std::cout << table2.at<int>(i, j) << std::endl;
		//		if (table_in.is_open())
		//		{
		//			table_in << table2.at<int>(i, j) << std::endl;
		//		}		
		//	}
		//}
		if (table_in.is_open())
		{
			for (int i = 0; i <differences.size(); i++) {
				//std::cout << differences[i] << std::endl;
				table_in << (int) differences[i] << std::endl;
			}	
		}
		table_in.close();
	}

	void Region_Graph::HistogramMatching(std::vector<float> differences, cv::Mat &table, std::vector <int> &lookUp_table, int &binsize) {
	
		//regions hist		
		sort(differences.begin(), differences.end());
		std::vector<float> rArr(differences.size());
		int c = 0;
		for (int i = 0; i < differences.size(); i++) {
			rArr[c] =  differences[i];
			//std::cout << rArr[c] << std::endl;
			c++;
		}
		float m = *max_element(rArr.begin(), rArr.end()); //find max value of data points
		int frsize;
		if (m > c - 1)
			frsize = m + 1;
		else frsize = c - 1;

		////removed uplicates
		//auto end = rArr.end();
		//for (auto it = rArr.begin(); it != end; ++it) {
		//	end = std::remove(it + 1, end, *it);
		//}
		//rArr.erase(end, rArr.end());
		//
		////
		//if (rArr.size() < frsize) frsize = rArr.size()+1;

		//table hist
		int numbins =  table.rows * table.cols;
		std::vector<float> tArr(numbins);
		int count = 0;
		for (int i = 0; i < table.rows; i++) {
			for (int j = 0; j < table.cols; j++) {
				tArr[count] =  (float)table.at<float>(i, j); //count ;
				//std::cout << tArr[count] << std::endl;
				count++;
			}
		}
		float mm = *max_element(tArr.begin(), tArr.end());// findmax(tArr, count - 1); //find max value of data points
		int fsize;
		if (mm > count - 1)
			fsize = mm + 1;
		else fsize = count - 1;

		//**************
		//sort(tArr.begin(), tArr.end());
		////removed uplicates
		//auto end = tArr.end();
		//for (auto it = tArr.begin(); it != end; ++it) {
		//	end = std::remove(it + 1, end, *it);
		//}
		//tArr.erase(end, tArr.end());
		//if (tArr.size() < fsize) fsize = tArr.size();
		/***************/
		int bsize = 500;// (int)max(m, mm) + 1;  //m+1;//
		float ratio = (float)(max(m, mm)+1) / bsize;
		std::cout << "bsize: " << bsize << ", ratio:" << ratio << std::endl;
		std::vector<int> rfreq(bsize); //declare frequency array with an appropriate size
		for (int i = 0; i < rfreq.size(); i++) //initialize frequency array
			rfreq[i] = 0;
		//compute frequencies
		for (int i = 0; i < rArr.size(); i++)
			rfreq[(int)(rArr[i]/*/((m+1)/bsize) */)]++;
		std::cout << "rfreq bin size: " << rfreq.size() << std::endl;
		std::vector<int> tfreq(bsize); //declare frequency array with an appropriate size
		for (int i = 0; i < tfreq.size(); i++) //initialize frequency array
			tfreq[i] = 0;
		//compute frequencies
		for (int i = 0; i < tArr.size(); i++)
			tfreq[(int)(tArr[i]/*/((mm+1)/bsize)*/ )]++;
		
		std::cout << "tfreq bin size: " << tfreq.size() << std::endl;

		//calculate_cdf(reg_hist)

		calculate_cdf(m,rfreq, normalized_rcdf);
		//calculate_cdf(table_hist)

		calculate_cdf(mm,tfreq, normalized_tcdf);

		// lookup table
		lookUp_table = lookupTable(normalized_rcdf, normalized_tcdf, tArr);

		//**************mapped*******************

		std::vector<int> Mappedfreq(bsize);
		for (int i = 0; i < Mappedfreq.size(); i++) //initialize frequency array
			Mappedfreq[i] = 0;
		//compute frequencies
		std::vector<int> mapArr(differences.size());
		for (int i = 0; i < differences.size();i++) {
			int id = round(differences[i] /*/ ((m) / bsize)*/);
			//std::cout <<"differences[i]: " << round(differences[i]) <<  ", id: " << id << std::endl;
			mapArr.push_back(lookUp_table[id]);
		}
		float mmm = *max_element(mapArr.begin(), mapArr.end());

		for (int i = 0; i < differences.size(); i++)
		{
			int id =round(differences[i]);
			Mappedfreq[/*mapArr[i]*/  lookUp_table[id]]++;
		}
		
		//calculate cdf after mapping		
		std::vector<float> normalized_Mappedtcdf;
		calculate_cdf(mmm, Mappedfreq, normalized_Mappedtcdf);

		//for (int i = 0; i < edgeweights.size();i++) {
		//	if (lookUp_table[int(edgeweights[i])] == -1)
		//	{
		//		std::cout <<  "null at " << i << std::endl;
		//	}
		//	//std::cout << edgeweights[i] << "  -->  "<< lookUp_table[(int)floor(edgeweights[i] / m * (mm))] << std::endl;
		//}
		std::cout << "lookup table size:" << lookUp_table.size() << " , normalized_rcdf size: " << normalized_rcdf.size() 
			<< " , normalized_tcdf size: " << normalized_tcdf.size() << std::endl;
		std::cout << "m: " << m << ", mm: " << mm << ", mmm: " << mmm << std::endl;
		ofstream table_in(path + "diff_colors_Before_After_HM.csv");
		//ofstream table_HM(path + "diff_colors_cmc_Before_HM.csv");
		if (table_in.is_open() /*|| table_HM.is_open()*/)
		{
			for (int i = 0; i <differences.size(); i++) {
				table_in << (int)differences[i] << "," << lookUp_table[round(differences[i] /*/((m) /bsize)*/ )]<< std::endl;
			}
		}
		table_in.close();
		//table_HM.close();
		///************************/
		ofstream table_rcdf(path + "normalized_rcdf_tcdf_colorful_Euclid.csv");
		ofstream table_Mappedcdf(path + "normalized_Mappedtcdf_colorful_Euclid.csv");
		if (table_rcdf.is_open() )
		{
			for (int i = 0; i <normalized_rcdf.size(); i++) {
				table_rcdf << normalized_rcdf[i] << "," << normalized_tcdf[i] << std::endl;
			}
		}
		if (table_Mappedcdf.is_open() )
		{
			for (int i = 0; i <normalized_rcdf.size(); i++) {
				table_Mappedcdf << normalized_rcdf[i] << "," << normalized_Mappedtcdf[i] << std::endl;
			}
		}
		table_rcdf.close();
		table_Mappedcdf.close();

		/****Histograms******/
		ofstream table_hist1(path + "Histogram_image.csv");
		ofstream table_hist2(path + "Histogram_palette.csv"); ofstream table_hist3(path + "Histogram_mapped.csv");
		for (int i = 0; i <rfreq.size(); i++) {
			table_hist1 << rfreq[i] << std::endl;
		}
		for (int i = 0; i <tfreq.size(); i++) {
			table_hist2 << tfreq[i] << std::endl;
		}
		for (int i = 0; i <Mappedfreq.size(); i++) {
			table_hist3 << Mappedfreq[i] << std::endl;
		}
		table_hist1.close(); table_hist2.close(); table_hist3.close();
		std::cout << "files are closed." << std::endl;
	}

	void Region_Graph::printGraph() // V= num regions
	{
		for (int v = 0; v < adj.size(); ++v)
		{
			cout << "\n Adjacency list of region "
				<< v << "\n head ";
			for (auto x : adj[v])
				std::cout << "-> " << x;
			std::cout << std::endl;
		}
	}

} //namespace abstraction

