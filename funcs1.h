#ifndef FUNCS1_H
#define FUNCS1_H

// header file for funcs1


/*
 * Params[in]:
 * pointer to mesh m
 */
extern "C" {
extern int initABC( int downward_counts[4][4], int numberEntities[4]);
extern void resetVertIt();
extern void resetEdgeIt();
extern void resetFaceIt();
extern void incrementVertIt();
extern void incrementEdgeIt();
extern void incrementFaceIt();
extern void incrementElIt();
extern void resetElIt();
extern void setGlobalVertNumber(int val);
extern int getGlobalVertNumber();
extern void checkVars();
extern void checkNums();
extern void getVertCoords(double coords[][3], int sx, int sy);
extern int getEdgeCoords(double coords[2][3], int sx, int sy);
extern int getFaceCoords(double coords[][3], int sx, int sy);
extern int getElCoords(double coords[][3], int sx, int sy);


}

#endif
