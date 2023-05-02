#ifndef CRACK_H
#define CRACK_H

#include "defs.h"

#define bXD 2000
#define bYD 2000
#define MAXPATHLEN 512

#define COLS 2000
#define ROWS 2000



static int XD = 512;
static int YD = 512;

typedef struct {
	int nearend;
	int farend;
	int cost;
} edgetype;


typedef struct {
	intvector loc;
	int id;
	int partcost;
	int dist;
	int regionid;
	edgetype edges[9];
	int whichpass;
	int numedges;

} nodetype;


//nodetype pix[ROWS*COLS];




//rgbvector screen[ROWS][COLS];


//rgbvector m[YD][XD];

//int BASE[5][5]; //[bXD][bYD];
//rgbvector postfilter_c[ROWS][COLS];

#define NOTEXIST -99

//rgbvector BACO[bXD][bYD];



//void makegraph(void);

//void processmany(void);

//int TARGREGSIZE = 40;


/*------------------------------------------------------------------------*/

#endif