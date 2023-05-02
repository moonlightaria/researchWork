#ifndef DEFS_H
#define DEFS_H

typedef struct {
	float x;
	float y;
	float z;
} vector;

typedef struct {
	double red;
	double green;
	double blue;
} rgbtotal;

typedef struct {
	unsigned char red;
	unsigned char green;
	unsigned char blue;
} rgbvector;

typedef struct {
	int x;
	int y;
	int z;
} intvector;


#define PI 3.1415926536
#define TRUE 1
#define FALSE 0

#endif
