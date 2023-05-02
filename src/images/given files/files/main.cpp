/***********************************************************************************************/
/*     You can change the colorspaces for diffrent color diffrence methods;                    */
/*     Euclidean to CMC method: ColorSpace::EuclideanComparison to ColorSpace::CmcComparison   */ 
/*     Change paths for results, images and filtered images before running the code            */
/*     Make a new folder in results folder for any new image you want to test                  */
/***********************************************************************************************/

#include <iostream>
#include "opencv\cv.hpp"
#include "opencv2\core\core.hpp"
#include "opencv2\highgui\highgui.hpp"
#include "opencv2\imgproc\imgproc.hpp"
#include "Image.h"
#include <stdlib.h>
#include "ColorSpace.h"
#include <time.h>

using namespace cv;
using namespace abstraction;


int main() {

	char filenames[80][20] = { "metal2", "dog1","bird","fall3","indian1024", "bike1","12003", "65019","175083",
		"PARS_Export","deer2", "yemeni1024",
		//"snow1024","city1024","cabbage1024","barn1024","cat1024","mac1024","angel1024",  "berries1024","arch1024",
		// "mountains1024","daisy1024","athletes1024","headlight1024",  "tomato1024","darkwoods1024","oparara1024", "rim1024","toque1024","desert1024",
		
	};

	char curfile[80];
	clock_t t1, t2;

	/**********process all images************************/
 	for (int k = 0; k <5; k++)
	{
		t1 = clock();
		int mask = 160;// 160;// CRGF mask size 
		sprintf_s(curfile, "E:/myApp/geofilter/images/%s.jpg", filenames[k]);
		std::string str(filenames[k]);
		cv::Mat image = cv::imread(curfile, 1);
		//cv::Mat image_down;
		//resize(image, image_down, cv::Size(6500 / 4, 4000 / 4), CV_INTER_LINEAR);

		//image = image_down.clone();

		abstraction::Image *in_img = new abstraction::Image(image);
		
		std::string path = "D:/Git codes/ReColorRegions/results/"+std::to_string(mask)+"/"+str+"/";
		in_img->path = path;
		in_img->filteredPath = "E:/myApp/geofilter/results/mask"+std::to_string(mask)+"/filtered_"+str+".png";
		in_img->filtered = imread(in_img->filteredPath, 1); // CRGF: edge-aware smoothing filter
		in_img->mask = mask;

		/*******make 8-connectivity graph on pixels**************/
		in_img->makegraph();
		/*******get the pos-neg residuals**************/
		in_img->calcResiduals(in_img->filtered, in_img->colorRes, in_img->pres, in_img->nres, in_img->signedRes);
		/*****set intensity or colors and set signed residuals***********/
		in_img->setResiduals();		
		/******get onnected components from residuals : we use these CC to select the seeds **********/
		in_img->getConnectedComponents( 3, (void*)in_img); //threshval = 3 ..  (Dijkstra is disabled)
	
		/*****seed generation (2 methods :max value / median position) *************/
		in_img->SeedGeneration(mask, in_img, curfile);

		/******** region generation ******  larger regions of smaller residuals **************/
		std::cout << "   Start create regions for  " << str << std::endl;
		
    	in_img->createRegions(mask, str); //based on cumulative error

		/****************Make a graph on the regions and recolor them ************************/
		in_img->makeGraphonRegions(image, str);

		delete in_img;
		t2 = clock();
		float diff((float)t2 - (float)t1);
		std::cout << diff / CLOCKS_PER_SEC << "  seconds to process "<< filenames[k] << std::endl;

	}
	std::system("PAUSE");
	cv::waitKey(0);
	return 0;
}