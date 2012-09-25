// -*- c++ -*-
// Copyright 2009 Isis Innovation Limited
//
// EyeGame.h
// Declares the EyeGame class
// EyeGame is a trivial AR app which draws some 3D graphics
// Draws a bunch of 3d eyeballs remniscient of the 
// AVL logo
//
#ifndef __EYEGAME_H
#define __EYEGAME_H

#include <TooN/TooN.h>
#include "OpenGL.h"
#include "Game.h"
#include "MapPoint.h"
#include "Map.h"
#include "Box.h"

namespace PTAMM {

using namespace TooN;

class EyeGame : public Game
{
  public:
    EyeGame( );
    virtual void Draw3D(GLWindow2 &glWindow, Map &map, SE3<> se3CfromW );
    virtual void Reset();
    virtual void Init();
    virtual void HandleKeyPress( std::string sKey );
    virtual void HandleClick(Vector<2> v2VidCoords, Vector<2> v2UFB, Vector<3> v3RayDirnW, Vector<2> v2Plane, int nButton);
    
private:
    
    Vector<3> _SelectClosestMapPoint(Vector<3> v3RayDirnW);
    void _DrawOrientPoints();
    void _HandleOrient(Vector<3> v3RayDirnW, Vector<2> v2Plane, int nButton);
    void _DrawBlocks();
    void _UpdateDistanceField();
    void _DrawWireframeWalls();
    void _DrawWall();
    void _DrawGrid();
    void _UpdateBoxes();
    void _DrawBoxes();
    
    void SetupModelView(SE3<> se3WorldFromCurrent = SE3<>());
    void SetupFrustum();
    void DrawCamera(SE3<> se3, bool bSmall=false);

    vector< vector<double> > distanceField;
    vector< vector< Vector<3> > > fieldCoordinates;
    vector< Box > mvbpBoxes;
    vector< Box > mvbpBoxesWithPts;
    
    bool drawFrame;
    
    // create a window to visualize what's going on
    //GLWindow2 mGLWindow2;
    
                               
    int mvDivs;
    int mwDivs;
    int mnFrameCounter;
    
    Map * mpMap;                            // The associated map
    SE3<> mse3CfW;                          // The current camera position
    
    // for third person view
    SE3<> mse3ViewerFromWorld;
    Vector<3> mv3MassCenter;
   
     
    SE3<> mse3MfromW;                       // Model's pose and location in the world.
    double mdScaleMult;                     // the scale multiplier for the model
    double mdScale;                         // the model's scale;
    double mdFieldMultiplier;               // a constant multiplier on the field
    double mdFieldWidth;
    double mdFieldHeight;
    double mdMinHeight;
    
    Vector<3> mvTranslation;

    bool planeIntersect1, planeIntersect2, planeIntersect3;

    bool blockVisible;
    Vector<3> mp1, mp2, mp3;
    Vector<3> wp1, wp2, wp3;
    double mdNoisePower;
    double mdMasterHeight;
    
    // Perlin noise
        double noise( double x, double y, double z );
        double octaveNoise( const Vector<3>& pt, int octaves );
        double fade( double t );
        double lerp( double t, double a, double b );
        double grad( int hash, double x, double y, double z );
        static int p[512];

};
    
}

#endif
