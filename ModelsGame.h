// Copyright 2009 Isis Innovation Limited
//
// C++ Interface: ModelsGame
//
// Description: An AR game for placing 3ds models in a map
//
//
// Author: Robert Castle <bob@robots.ox.ac.uk>, (C) 2009
//
//
#ifndef MODELSGAME_H
#define MODELSGAME_H

#include "Game.h"
#include "Model3ds.h"
#include "ModelBrowser.h"
#include "ModelControls.h"
#include "ModelsGameData.h"

#include <cvd/image.h>

namespace PTAMM {

/**
 * AR game for displaying 3DS models in a map
 *
 * @author Robert Castle <bob@robots.ox.ac.uk>
 */
class ModelsGame : public Game
{
  public:
    ModelsGame();
    ~ModelsGame();

    void Draw3D( GLWindow2 &glWindow, Map &map, SE3<> se3CfromW);
    void Draw2D( GLWindow2 &glWindow, Map &map);
    
    void Reset();
    void Init();
    
    void HandleClick(Vector<2> v2VidCoords, Vector<2> v2UFB, Vector<3> v3RayDirnW, Vector<2> v2Plane, int nButton);

    std::string Save(std::string sMapPath);
    void Load(std::string sDataPath);

  private:
    Vector<3> _SelectClosestMapPoint(Vector<3> v3RayDirnW);
    bool _HandleModelBrowserActions( Vector<2> v2VidCoords );
    void _DrawSelectedTarget();
    void _DrawAxes();
    bool _GetVertexIntersect(Vector<3>& position, Vector<2> v2Plane, Vector<3> v3RayDirnW);
    void _HandleOrient(Vector<3> v3RayDirnW, Vector<2> v2Plane, int nButton);
    void _DrawOrientPoints();
    
  private:
    Map * mpMap;                    // The associated map
    SE3<> mse3CfW;                  // The current camera position

    ModelsGameData mData;           // The game data
    ModelBrowser mModelBrowser;     // The model browser
    ModelControls mModelControls;   // The model controls

    GLuint mnTargetTex;             // The current model target texture
    
    bool modelIntersect1, modelIntersect2, modelIntersect3;
    bool planeIntersect1, planeIntersect2, planeIntersect3;
    
    Vector<3> mp1, mp2, mp3;
    Vector<3> wp1, wp2, wp3;
       
};

}

#endif
