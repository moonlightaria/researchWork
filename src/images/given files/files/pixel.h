#pragma once
#include "opencv\cv.hpp"
#include "crack.h"
#include <list>

#define DOWN 0
#define DOWNL 1 //
#define LEFT 2 // 1  
#define TOPL 3 //
#define TOP   4 //2
#define TOPR  5
#define RighT 6 //3
#define DOWNR 7


//class Image;

class pixel {
	friend class Edges;
private:
	int regionid;
public:
	pixel() {
		max = NULL;
		min = NULL;
		edges[8] = { 0 };
		finalized = false;
		contour = false;
		newCol = cv::Vec3b(0, 0, 0);
		chosen = false;
		pblobId = 0;
		nblobId = 0.0f;
		visibility = false;
		regionid = -50;
		tempRegion = -50;
		finalAssignedRegion = -50;
		regColor = cv::Vec3b(0, 0, 0);
		parentsize = 0;
		difWeight = 0;
		ilegal = false;
		state = 0;
	};

	pixel(cv::Point2f _pos, uchar _color) {
		pos = _pos;
		//xx = _pos.x; yy = pos.y;
		gray = _color;
		max = NULL;
		min = NULL;
		edges[8] = { 0 };
		finalized = false;
		contour = false;
		newCol = cv::Vec3b(0,0,0);
		chosen = false;
		pblobId = 0;
		nblobId = 0.0f;
		visibility = false;
		regionid = -50;
		tempRegion = -50;
		finalAssignedRegion = -50;
		regColor = cv::Vec3b(0, 0, 0);
		parentsize = 0;
		difWeight = 0;
		ilegal = false;
		state = 0;
	}

	pixel(cv::Point2f _pos, cv::Vec3b _color) {
		pos = _pos;
		color = _color;
		purturbed = color;
		gray = 0.299* _color.val[3] + 0.587*_color.val[1] + 0.114*_color.val[0];
		max = NULL;
		min = NULL;
		edges[8] = { 0 };
		finalized = false;
		contour = false;
		newCol = cv::Vec3b(0, 0, 0);
		reColor = _color;
		pblobId = 0;
		nblobId = 0.0f;
		visibility = false;
		regionid = -50;
		tempRegion = -50;
		finalAssignedRegion = -50;
		regColor = cv::Vec3b(0, 0, 0);
		parentsize = 0;
		difWeight = 0;
		ilegal = false;
		state = 0;
		regionIds.resize(2000);
	}
	~pixel() {
	}

	float dist;
	int parentsize;
	float sumColdist;
	int nb; // id in the path
	int sn; // source num
	float score; // score is a normalized value based on the distance from ancesstor 

	float difWeight; // difference of curent pix weight from 10 away parent
	float seedrange; // for switch between different masks 

	int numedges;
	bool finalized;
	bool contour;
	bool grow;
	bool ilegal;
	bool finalReg, tempReg;

	bool visibility;

	int tempRegion, finalAssignedRegion;

	std::vector<int> regionIds;
	int partcost;
	int regionsize;

	bool black, white, nut; 
	bool edge;
	bool chosen;// chosen for the range of residual values
	bool lock; // lock from growing
	bool belongToreg;
	bool dead;
	int state; // = 0 (unknown/not visited); = 5 (on the heap) ; = 10 (temp) ; = 15 (final) ; = 20 (rejected)

	cv::Point2f pos;

	Edges *edges[8];
	pixel *parent;
	std::vector<int> childrenId;
	std::vector<pixel*> childrens;
	std::vector<pixel*> parents;
	pixel *ancesParent;
	int CCid;
	int ancesId;
	int pblobId, nutblobId;
	float nblobId;

	float W12,W21; // probablit that pixel blong to the inner bound

	int slicId;
	int level;
	std::vector<std::pair<int, cv::Vec3b >> leveles_values; // region Id , color

	cv::Vec3b color,oldcolor,newCol, reColor,purturbed;
	cv::Vec3f LabColor;
	cv::Vec3f LabColors[4];
	uchar newColor;


	int gray;
	float res;
	float coherency;
	double orientation;
	bool max; bool min;
	bool visited; // for seperating CC
	cv::Vec3b regColor;// region color

	int id;
	bool idc480, idc960;
	bool idc40, idc160;
	bool idcn480, idcn960;
	bool idcn40, idcn160;

	bool isVisited;
	std::vector<pixel*> neighbours;
	float pdf,pdf_s;

	void addChild(pixel *_child) {
		//childrenId.push_back(_childid);
		childrens.push_back(_child);
	}
	void addRegionId(int _regid) {
		regionIds.push_back(_regid);
	}
	void setRegionId(int _regid) {
		regionid = _regid;
	}
	void setCCId(int _regid) {
		CCid = _regid;
	}
	int getRegId() {
		return regionid;
	}
	void setVisibility(bool flag) {
		visibility = flag;
	}

	void setScore(int _sn,int _nb) { // sn = source number, nb = id in the path
		sn = _sn; // contour id
		nb = _nb;
		ancesId = parent->ancesId;
		score = sqrt(colorDistTo(ancesParent));// sqrt(distTo(ancesParent));
	}

	void setNewColor(uchar newCol) {
		newColor = newCol;
	}
	void setreColor(cv::Vec3b reCol) {
		reColor = reCol;
	}
	float distTo(pixel *p)
	{
		return ((pos.x - p->pos.x) * (pos.x - p->pos.x) + (pos.y - p->pos.y) * (pos.y - p->pos.y)); // XxX
	}
	float colorDistTo(pixel *p)
	{
		//setLabColor(seedrange);
		//p->setLabColor(seedrange);
		//Lab color difference
		return pow( sqrt((2*(float)(LabColor.val[0] - p->LabColor.val[0]) * (float)(LabColor.val[0] - p->LabColor.val[0]) +
			(float)(LabColor.val[1] - p->LabColor.val[1]) * (float)(LabColor.val[1] - p->LabColor.val[1]) +
			(float)(LabColor.val[2] - p->LabColor.val[2]) * (float)(LabColor.val[2] - p->LabColor.val[2]))),4.5) ; // CxC
		//// rgb color differences
		//return sqrt(((float)(color.val[0] - p->color.val[0]) * (float)(color.val[0] - p->color.val[0]) + 
		//	(float)(color.val[1] - p->color.val[1]) * (float)(color.val[1] - p->color.val[1]) +
		//	(float)(color.val[2] - p->color.val[2]) * (float)(color.val[2] - p->color.val[2]) )); // CxC
	}
	float angular_colorDistTo(pixel *p)
	{

		cv::Vec3f a = cv::Vec3f(color.val[0], color.val[1], color.val[2]); 
		cv::Vec3f b = cv::Vec3f(p->color.val[0], p->color.val[1], p->color.val[2]);
		float La = sqrt((color.val[0] * color.val[0] + color.val[1] * color.val[1] + color.val[2] * color.val[2]));
		float Lb = sqrt((p->color.val[0] * p->color.val[0] + p->color.val[1] * p->color.val[1] + p->color.val[2] * p->color.val[2]));
		float b1 =  a.dot(b) / La;
		float h = sqrt(Lb*Lb - b1*b1);
		float area = (La * h) /2.0;

		float anglDist = pow(( 1.0f - a.dot(b)/ (La * Lb)),0.05);// (1- cost)^ 1/20
		//float anglDist = 1.0f - (1.0/(area+0.001))* a.dot(b) / (La * Lb);// (1- (1/area)cost)

		return  1000*anglDist; 
	}
	float sqrtColorDistTo(pixel *p)
	{
		return  sqrt(((float)(color.val[0] - p->color.val[0]) * (float)(color.val[0] - p->color.val[0]) +
			(float)(color.val[1] - p->color.val[1]) * (float)(color.val[1] - p->color.val[1]) +
			(float)(color.val[2] - p->color.val[2]) * (float)(color.val[2] - p->color.val[2]))); // CxC
	}
	float sqrtColorDistTo(cv::Vec3b p)
	{
		return  sqrt(((float)(color.val[0] - p.val[0]) * (float)(color.val[0] - p.val[0]) +
			(float)(color.val[1] - p.val[1]) * (float)(color.val[1] - p.val[1]) +
			(float)(color.val[2] - p.val[2]) * (float)(color.val[2] - p.val[2])) ); 
	}
	float SignDistTo(pixel *p)
	{
		return (expf(-1 *abs(res - p->res)));
	}

	cv::Vec3b getColor() {
		return color;
	}
	cv::Point2f getPos() {
		return pos;
	}
	void setId(int _id) {
		id = _id;
	}
	int getId() {
		return id;
	}
	int getPosBlobid() {
		return pblobId;
	}
	float getNegBlobid() {
		return nblobId;
	}
	int getNutBlobid() {
		return nutblobId;
	}
	void setNegBlobId(float nid) {
		nblobId =  nid;
	}
	void setNutBlobId(int nutid) {
		nutblobId = nutid;
	}
	void setPosNegBlobId(int pid, bool flag) {
		if(flag == true)
			pblobId = pid;
		if (flag == false)
			nblobId -=  pid;
		if (flag == NULL)
			nutblobId = pid;
	}

	void setPos(cv::Point2f p) {
		pos = p;
	}
	void setColor(cv::Vec3b c) {
		color = c;
	}
	void setLabColor(cv::Vec3f _c) {
		LabColor = _c;
	}
	void setLabColor(float range) {

		if (range < 0.25f)
		{
			this->LabColor = LabColors[3];
		}
		if (range < 0.5f && range >= 0.25f)
		{
			this->LabColor = LabColors[3];
		}
		if (range < 0.75f && range >= 0.5f)
		{
			this->LabColor = LabColors[3];
		}
		//if (range < 1.0f && range >= 0.75f)
		//{
		//	this->LabColor =  LabColors[3];
		//}
	}
	void setIntensity(int g) {
		gray = g;
	}
	int getIntensity() {
		return int(gray);
	}
	void setResidual(float _res) {
		res = _res;
	}
	void setLevel(int _lev) {
		level = _lev;
	}
	void initialLevelValues(int maxLev) {
		for (int i = 0; i < maxLev;i++) {
			leveles_values[i].first = -1;
			leveles_values[i].second = cv::Vec3b(-1,-1,-1);
		}
	}
	void setValueforLevel(int _lev, cv::Vec3b _col, int reg_Id) {
		leveles_values[_lev].first = reg_Id;
		leveles_values[_lev].second = _col;
	}
	void setParent(pixel *pp) {
		parent = pp;
		parentsize++;
	}
	pixel* returnParent(int _pid) {// return parent number _pid
		int count = 0;
		pixel* p = this->parent;
		while (count < _pid) {
			p = p->parent;
			count++;
		}
		return p;
	}

	std::vector<pixel*> returnChilderns() {// return all children that has the same regionid
		int regid = this->getRegId();
		std::vector<pixel*> allchildren;
		allchildren.push_back(this);

		std::list<pixel*> childs;
		std::list<pixel*>::iterator it1;
		for (int i = 0; i < this->childrens.size();i++) {
			childs.push_back(this->childrens[i]);
		} 
		while (childs.size() > 0) {
			for(it1 = childs.begin(); it1 != childs.end(); ++it1){
				if ((*it1)->getRegId() == regid) {
					allchildren.push_back((*it1));
					for (int j = 0;j < (*it1)->childrens.size();j++) {
						childs.push_back((*it1)->childrens[j]);
					}
				}
				//else
				childs.erase(it1);
			}
		}
		return allchildren;
	}
};

