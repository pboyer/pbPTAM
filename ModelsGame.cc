// Copyright 2009 Isis Innovation Limited
//
// C++ Implementation: ModelsGame
//
// Description: An AR game for placing 3ds models in a map
//
// Author: Robert Castle <bob@robots.ox.ac.uk>, (C) 2009
//
#include "ModelsGame.h"
#include "OpenGL.h"
#include "MapPoint.h"
#include "Map.h"

#include <cvd/gl_helpers.h>
#include <cvd/image_io.h>
#include <gvars3/instances.h>


namespace PTAMM {

using namespace CVD;
using namespace TooN;
using namespace GVars3;
             
class Map;
             
/**
 * Constructor
 */
ModelsGame::ModelsGame()
 : Game( "Models" ),
   mData( Name() ),
   mModelControls( mData )
{
  Reset();
}

/**
 * Destructor
 */
ModelsGame::~ModelsGame()
{
}



/**
 * Reset the game.
 */
void ModelsGame::Reset()
{
  mData.Reset();
}


/**
 * Initialize the game.
 */
void ModelsGame::Init()
{
  if( mbInitialised ) {
    return;
  }

    modelIntersect1 = false;
    modelIntersect2 = false;
    modelIntersect3 = false;
    planeIntersect1 = false;
    planeIntersect2 = false;
    planeIntersect3 = false;
    
  bool bOk = false;

  // initialize the model browser
  bOk = mModelBrowser.Init( ImageRef(BROWSER_X, BROWSER_Y), "ARData/Overlays/model_browser.png");

  if( !bOk ) {
    return;
  }

  // initialize the model controls
  mModelControls.Init();
  

  // load the current model target
  try {
    Image<Rgba<CVD::byte> > imTarget;
    CVD::img_load( imTarget, "ARData/Overlays/selected.png" );
    glGenTextures( 1, &mnTargetTex);
    glBindTexture( GL_TEXTURE_2D, mnTargetTex );
    glTexImage2D(imTarget);
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );  
  }
  catch(CVD::Exceptions::All err) {
    cerr << "Failed to load target image " << "ARData/Overlays/selected.png" << ": " << err.what << endl;
    return;
  }

 
  mbInitialised = true;
}



/**
 * Draw the 3D components (the 3DS models)
 * @param glWindow The GL Window
 * @param map the current map
 * @param se3CfromW The current camera position
 */
void ModelsGame::Draw3D( GLWindow2 &glWindow, Map &map, SE3<> se3CfromW)
{
  if( !mbInitialised ) {
    Init();
  }

  if( mData.mbHideAR )  {
    return;
  }

  mpMap = &map;
  mse3CfW = se3CfromW;
  
  glEnable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  GLfloat tmp[]={0.51, 0.5, 0.5, 0.5};
  glLightfv(GL_LIGHT0, GL_AMBIENT, tmp);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, tmp);
  glLightfv(GL_LIGHT0, GL_SPECULAR, tmp);

  GLfloat tmp2[]={1.0, 0.0, 1.0, 0.5};
  glLightfv(GL_LIGHT0, GL_POSITION, tmp2);
  GLfloat tmp3[]={1.0, 1.0, 1.0, 1.0};
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, tmp3);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 70.0);

  glMatrixMode(GL_MODELVIEW);
  
  glLoadIdentity();
  glMultMatrix(SE3<>());

  // Display the map points?
  if( !mData.mbHidePoints )
  {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPointSize(5);
    glColor3f(1,1,1);
    glBegin(GL_POINTS);
    for( size_t ii = 0; ii < map.vpPoints.size(); ++ii )
    {
      MapPoint &p = *(map.vpPoints[ii]);
      if( p.nSourceLevel != 0 ) {
        continue;
      }

      glVertex( p.v3WorldPos );
    }
    glEnd();
  }  

    _DrawOrientPoints();

  // Draw the origin axes ?
  if( GV3::get<int>("ModelsGame.DrawOrigin", "1", SILENT) ) {
    _DrawAxes();
  }


  // Draw the models
  std::vector<Model3DS*>::const_iterator itr;
  for ( itr = mData.Models().begin(); itr != mData.Models().end(); ++itr )
  {
    (*itr)->Draw();
  }


  // Find the current model and draw the target
  Model3DS * m = mData.CurrentModel();
  if( m && !mData.mbHideControls )
  {
    glPushMatrix();
    glMultMatrix( m->Pose() );
    glScaled( m->Diameter(), m->Diameter(), m->Diameter() );
    _DrawSelectedTarget();
    glPopMatrix();
  }
  
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  // A rather nasty hack to get around the window event limitations.
  // Check to see if the user is holding the mouse button down and if
  // it is on a repeating button, call it again.
  if( mData.mDisplayState == ModelsGameData::ControlState &&
      glWindow.IsMouseButtonPressed( CVD::GLWindow::BUTTON_LEFT ) ) {
    mModelControls.HandlePressAndHold();
  }
}

void ModelsGame::_DrawOrientPoints()
{    
    
    if (modelIntersect1) {
        
        // draw a point
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glPointSize(20);
        glColor3f(1,0,0);
        glBegin(GL_POINTS);
        
        glVertex(mp1);
        
        glEnd();
        
    }
    
    if (modelIntersect2) {
        
        // draw a point
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glPointSize(20);
        glColor3f(0,0,1);
        glBegin(GL_POINTS);
        
        glVertex(mp2);
        
        glEnd();
        
        // draw a line
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glLineWidth(5);
        glColor3f(1,1,1);
        glBegin(GL_LINES);
        
        glVertex(mp1);
        glVertex(mp2);
        
        glEnd();
        
    }
    
    if (modelIntersect3) {
        
        // draw a point
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glPointSize(20);
        glColor3f(1,0,0);
        glBegin(GL_POINTS);
        
        glVertex(mp3);
        //glVertex3d( mi1[0],mi1[1],mi[2] );
        
        glEnd();
        
        // draw all of the marker points
        std::vector<MapPoint*> vpPoints = mpMap->vpPoints;
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glPointSize(8);
        glColor3f(0.8, 0.8, 0.8);
        glBegin(GL_POINTS);
        
        for (int i = 0; i < vpPoints.size(); i++) {
            glVertex(vpPoints[i]->v3WorldPos);
        }
        
        glEnd();
        
    }
    
    if (planeIntersect1) {
        
            glPointSize(20);
        
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            glColor3f(1,0,1);
            glBegin(GL_POINTS);
            
            glVertex(wp1);
            
            glEnd();

  
        
        if (planeIntersect2) {
            
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            glColor3f(1,0,1);
            glBegin(GL_POINTS);
            
            glVertex(wp2);
            
            glEnd();
            
        }
        
        if (planeIntersect3) {
            
            // draw a point
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            glColor3f(1,0,1);
            glBegin(GL_POINTS); 
            
            glVertex(wp3);
            
            glEnd();
            
        }
    }
}

/**
 * Draw the 2D components (menus etc)
 * @param glWindow The GL Window
 * @param map The current map
 */
void ModelsGame::Draw2D( GLWindow2 &glWindow, Map &map)
{
  glEnable(GL_BLEND);
  
  switch(mData.mDisplayState)
  {
    case ModelsGameData::BrowserState:
      mModelBrowser.Draw(glWindow);       // draw the model browser
      break;
    case ModelsGameData::ControlState:
      mModelControls.Draw(glWindow);      // draw the controls
      break;
    case ModelsGameData::HiddenState:     // nothing to draw
    default:
      //do nothing
      break;
  }
  
  glDisable(GL_BLEND);
}
    
void ModelsGame::_HandleOrient(Vector<3> v3RayDirnW, Vector<2> v2Plane, int nButton) {
    
    //select appropriate reponse for each overlay
    cout << "button: " << nButton << endl;
    
    if (nButton == 0) {
        cout << nButton << endl;
    }
    
    if (!modelIntersect1) {
        Vector<3> position = makeVector(0,0,0);
        bool foundIntersection = this->_GetVertexIntersect(position, v2Plane, v3RayDirnW);
        if (foundIntersection) {
            cout << "Found model point 1" << endl;
            cout << "Supply model point 2 to define scale and orient coordinates" << endl;
            modelIntersect1 = true;
            mp1 = position;
        } else {
            cout << "Missed initial point, try again!" << endl;
        }
        return;
    }
    
    if (modelIntersect1 && !modelIntersect2) {
        Vector<3> position = makeVector(0,0,0);
        bool foundIntersection = this->_GetVertexIntersect(position, v2Plane, v3RayDirnW);
        
        if (foundIntersection) {
            cout << "Found model point 2" << endl;
            cout << "Supply model point 3 to completely define the plane of points you want to orient with" << endl;
            modelIntersect2 = true;
            mp2 = position;
        } else {
            cout << "missed model point 2, restarting orient" << endl;
            modelIntersect1 = false;
        }
        return;
    }
    
    if (modelIntersect1 && modelIntersect2 && !modelIntersect3) {
        Vector<3> position = makeVector(0,0,0);
        bool foundIntersection = this->_GetVertexIntersect(position, v2Plane, v3RayDirnW);
        
        if (foundIntersection) {
            cout << "Found model point 3" << endl;
            cout << "Now define the target points in world space.  Here's the model positions." << endl;
            modelIntersect3 = true;
            mp3 = position;
        } else {
            cout << "Missed model point 3, restarting orient" << endl;
            modelIntersect1 = false;
            modelIntersect2 = false;
        }
        return;
    }
    

    if (modelIntersect1 && modelIntersect2 && modelIntersect3 && !planeIntersect1) {
        cout << "Obtained target point 1"  << endl;
        
        wp1 = this->_SelectClosestMapPoint( v3RayDirnW );
        
        //wp1 = v2Plane;   
        planeIntersect1 = true;
        return;
    }
    
    if (planeIntersect1 && !planeIntersect2) {
        cout << "Obtained target point 2"  << endl;
        cout << "This defines the scale and mostly defines the transformation."  << endl;
        wp2 = this->_SelectClosestMapPoint( v3RayDirnW );
        
        //wp2 = v2Plane;
        planeIntersect2 = true;
        return;
    }
    
    if (planeIntersect1 && planeIntersect2) {
        
        wp3 = this->_SelectClosestMapPoint( v3RayDirnW );
        
        planeIntersect3 = true;
    
        
        // PERFORM REORIENTATION OF MODEL HERE
        
        std::vector< Model3DS* > models = mData.Models();
        
        if (models.size() > 0) {
            
            // The scale comes from the ratio between the lengths of the first two origin points
            // and the first two target points
            double scaleRatio = norm(wp1 - wp2) / norm(mp1 - mp2); 
            
            Vector<3> a = mp1;
            Vector<3> ap = mp2;
            Vector<3> app = mp3;
            
            Vector<3> b = wp1; 
            Vector<3> bp = wp2; 
            Vector<3> bpp = wp3; 

            Vector<3> a1 = ap - a;
            Vector<3> a2 = app - a;
            
            Vector<3> b1 = bp - b;
            Vector<3> b2 = bpp - b; 
            
            a1 = unit(a1);
            a2 = unit(a2);
            Vector<3> a3 = Zeros;
            
            a3 = a1 ^ a2;
            a2 = a3 ^ a1;
            a3 = unit(a3);
            a2 = unit(a2);
            
            // a2apnorm, a2appnorm, a2apppnorm makes up orthonormal coordinate system from model points
            // use rotation to express world points (bs) in model point coordinate system (as)
            Matrix<3> aRot = Zeros;
            aRot[0] = a1;
            aRot[1] = a2;
            aRot[2] = a3;
            
            // now we express the b points in terms of this coordinate system
            // then the three columns b1, b2, b3 define our 3d rotation matrix
            b1 = unit(b1);
            b2 = unit(b2);
            Vector<3> b3 = Zeros;
            
            b3 = b1 ^ b2;
            b2 = b3 ^ b1;
            b3 = unit(b3);
            b2 = unit(b2);
            
            Vector<3> b1a = aRot * b1;
            Vector<3> b2a = aRot * b2;
            Vector<3> b3a = aRot * b3;
            
            Matrix<3> bRotIna = Zeros;
            bRotIna[0] = b1a;
            bRotIna[1] = b2a;
            bRotIna[2] = b3a;
            
            // just the orthonormal columns obtained from the world positions
            // expressed in coordinate system defined by 3 points on model
        
            Matrix<3> bRot = bRotIna.T();
            SO3<> rotation = bRot; // put in nicer format
            
            // order of transformations is
            // 1) rotate at origin
            // 2) scale at origin
            // 3) translate to final position
            
            // so, our translation needs to ta
            Vector<3> ca = models[0]->mse3MfromW.get_translation();
            
            Vector<3> amca = mp1 - ca;
            Vector<3> bmcb = rotation*(amca);
            bmcb = scaleRatio * bmcb;
            Vector<3> bmca = wp1 - ca;
            
            Vector<3> newPos = bmca - bmcb + ca;
            
            
            models[0]->mdScale = models[0]->mdScale * scaleRatio;
            // add the new translation
            models[0]->MoveTo(newPos);  
            // this rotation is the last one we will perform
            models[0]->mse3MfromW.get_rotation() = rotation * models[0]->mse3MfromW.get_rotation();  
            
            
            
        }
        
        modelIntersect1 = false;
        modelIntersect2 = false;
        modelIntersect3 = false;
        planeIntersect1 = false;
        planeIntersect2 = false;
        planeIntersect3 = false;
        return;
    }
    
    
}


/*
void ModelsGame::_HandleOrient(Vector<3> v3RayDirnW, Vector<2> v2Plane, int nButton) {
        
    //select appropriate reponse for each overlay
    cout << "button: " << nButton << endl;
    
    if (nButton == 0) {
        cout << nButton << endl;
    }
    
    if (!modelIntersect1) {
        Vector<3> position = makeVector(0,0,0);
        bool foundIntersection = this->_GetVertexIntersect(position, v2Plane, v3RayDirnW);
        if (foundIntersection) {
            cout << "found model point 1" << endl;
            modelIntersect1 = true;
            mp1 = position;
        } else {
            cout << "missed initial point" << endl;
        }
        return;
    }
    
    if (modelIntersect1 && !modelIntersect2) {
        Vector<3> position = makeVector(0,0,0);
        bool foundIntersection = this->_GetVertexIntersect(position, v2Plane, v3RayDirnW);
        
        if (foundIntersection) {
            cout << "found model point 2" << endl;
            modelIntersect2 = true;
            mp2 = position;
        } else {
            cout << "missed model point 2, restarting orient" << endl;
            modelIntersect1 = false;
        }
        return;
    }
    
    // get third point for plane
    
    
    
    
    
    if (modelIntersect1 && modelIntersect2 && !planeIntersect1) {
        cout << "obtained out point1"  << endl;
        
        wp1 = this->_SelectClosestMapPoint( v3RayDirnW );
        
        //wp1 = v2Plane;   
        planeIntersect1 = true;
        return;
    }
    
    if (modelIntersect1 && modelIntersect2 && planeIntersect1) {
        
        wp2 = this->_SelectClosestMapPoint( v3RayDirnW );
        
        //wp2 = v2Plane;
        planeIntersect2 = true;
        
        // PERFORM REORIENTATION OF MODEL HERE
        
        std::vector< Model3DS* > models = mData.Models();
        
        // need a 3 point method
        // 
        
        if (models.size() > 0) {
     
            double scaleRatio = norm(wp1 - wp2) / norm(mp1 - mp2);
            cout << "scaleRatio" << scaleRatio << endl;
            
            Vector<3> a = mp1;
            Vector<3> ap = mp2;
            
            Vector<3> b = wp1; //makeVector(wp1[0],wp1[1],0);
            Vector<3> bp = wp2; //makeVector(wp2[0],wp2[1],0);
            
            
            // get angle between the vectors
            Vector<3> a2ap = ap - a;
            Vector<3> b2bp = bp - b;
            
            Vector<3> a2apnorm = unit(a2ap);
            Vector<3> b2bpnorm = unit(b2bp);
            
            // calculate rotation
            double yawToAdd = acos( a2apnorm*b2bpnorm ); // always makes positive angle, need to check cross product
            
            Vector<3> planeNorm = a2apnorm ^ b2bpnorm;
            
            Vector<3> z = makeVector(0,0,1);
            
            if (planeNorm * z < 0) {
                yawToAdd = -yawToAdd;
            }
            
            SE3<> rotate;
            Vector<3> rot = makeVector(0.0, 0.0, yawToAdd);
            rotate.get_rotation() *= SO3<>::exp(rot);
            
            Vector<3> ca = models[0]->mse3MfromW.get_translation();

            Vector<3> amca = a - ca;
            Vector<3> bmcb = rotate*(amca);
            bmcb = scaleRatio * bmcb;
            
            Vector<3> bmca = b - ca;
            
            //Vector<3> newTranslation = wp13 - mc2wp1;
            Vector<3> newPos = bmca - bmcb + ca;
            //newPos[2] = ca[2]; // maintain plane
            
            cout << "a2ap" << a2ap<< endl;
            cout << "b2bp" << b2bp << endl;
            cout << "rotationAngle" << yawToAdd << endl;
            cout << "new Position" << newPos << endl;
            cout << "current translation" << models[0]->mse3MfromW.get_translation() << endl;
            
            
            // SCALE THE MODEL
            models[0]->mdScale = models[0]->mdScale * scaleRatio;
            models[0]->MoveTo(newPos);
            models[0]->mse3MfromW.get_rotation() *= SO3<>::exp(rot); // add to model
            
            
            
            

            
        }
        
//        cen = models[i]->mv3Offset + cen;
//        cen = models[i]->mdScale * cen;
//        cen = models[i]->mse3MfromW * models[i]->mse3ModelOffset * cen;
        
        modelIntersect1 = false;
        modelIntersect2 = false;
        planeIntersect1 = false;
        planeIntersect2 = false;
        return;
    }
    
        
}
*/

/**
 * Handle a click from the user. The user will either be clicking
 * on a button overlay or on the scene.
 * @param v2VidCoords location in the image coords.
 * @param v2UFB Location in the undistorted framebuffer coords
 * @param v3RayDirnW Direction of the click ray
 * @param v2Plane The in-plane location of the click
 * @param nButton The button clicked
 */
void ModelsGame::HandleClick(Vector<2> v2VidCoords, Vector<2> v2UFB, Vector<3> v3RayDirnW, Vector<2> v2Plane, int nButton)
{
    
  switch(mData.mDisplayState)
  { 
    case ModelsGameData::BrowserState:
      _HandleModelBrowserActions( v2VidCoords );
      break;
    default:
        if( !mModelControls.HandleClick( v2VidCoords, nButton )
           && !mData.mbHideAR) 
          {
              _HandleOrient(v3RayDirnW, v2Plane, nButton);   
              
          }
          
        
      if( !mModelControls.HandleClick( v2VidCoords, nButton ) && mData.mbSnapTo
           && !mData.mbHideControls  && !mData.mbHideAR)
      {
        Model3DS *p = mData.CurrentModel();
        if( p )
        {
          Vector<3> v3Loc = _SelectClosestMapPoint( v3RayDirnW );
          p->MoveTo( v3Loc );
        }
      }
      break;
  }
}
    
/**
* Handle and act on the results from the user clicking on the model browser.
* @param v2VidCoords the video coordinates
* @return true if used the click
*/
bool ModelsGame::_GetVertexIntersect(Vector<3>& position, Vector<2> v2Plane, Vector<3> v3RayDirnW) {
    
    std::vector< Model3DS* > models = mData.Models(); // get all the models in the scene
    
    Vector<3> closestPt = makeVector(0,0,0);
    double closestDist = 100000000;
    
    Vector<3> o = mse3CfW.inverse().get_translation();     // camera origin?
    Vector<3> d = v3RayDirnW;
    double r = 3;  // default radius for sphere
    
    bool foundIntersection = false;
    
    for (int i = 0; i < models.size(); i++) {
        
        // need to take into account the actual transformation of the points
        vector< Vector<3> > vertices = models[i]->vertices;
        
        for (int j = 0; j < vertices.size(); j++) {
            
            // get vertex
            Vector<3> cen = vertices[j];
            
            cen = models[i]->mv3Offset + cen;
            cen = models[i]->mdScale * cen;
            cen = models[i]->mse3MfromW * models[i]->mse3ModelOffset * cen;
            
            //Compute A, B and C coefficients
            float a = d*d;
            float b = 2 * d * (o - cen);
            float c = (o - cen) * (o-cen) - (r * r);
            
            //Find discriminant
            float disc = b * b - 4 * a * c;
            
            // if discriminant is negative there are no real roots, so return 
            // false as ray misses sphere
            if (disc < 0)
                continue;
            
            double q;
            
            if (b < 0) {
                q = (-b + sqrt(disc))/2;
            } else {
                q = (-b - sqrt(disc))/2;
            }
            
           
            double t0 = q/a;
            double t1 = c/q;
            
            if (t0 < closestDist && t0 > 0) {
                foundIntersection = true;
                closestDist = t0;
                closestPt = cen;
            }
            
            if (t1 < closestDist && t1 > 0) {
                foundIntersection = true;
                closestDist = t1;
                closestPt = cen;
            }
            
        }
    }
    
    position = closestPt;
    
    return foundIntersection;
    
}



/**
 * Handle and act on the results from the user clicking on the model browser.
 * @param v2VidCoords the video coordinates
 * @return true if used the click
 */
bool ModelsGame::_HandleModelBrowserActions( Vector<2> v2VidCoords )
{
  ModelBrowser::CLICK_STATUS cs = mModelBrowser.HandleClick( v2VidCoords );

  if(cs == ModelBrowser::MODEL)
  {
    const ModelData * pModelData = mModelBrowser.GetSelectedModelData();

    if( pModelData ) {
      mData.AddModel( pModelData->sLocation, pModelData->sModelFile, pModelData->sName, pModelData->v3Rotation );
    }
  }
  else if( cs == ModelBrowser::CANCEL )
  {
    if( mData.NumModels() != 0 ) {
      mData.mDisplayState = ModelsGameData::ControlState;
    }
  }
  else if( cs == ModelBrowser::OK ) {
  }
  else {
    return false;
  }

  return true;

}


/**
 * Select the closest map point to the ray from click
 */
Vector<3> ModelsGame::_SelectClosestMapPoint( Vector<3> v3RayDirnW )
{
  if(!mpMap) {
    return makeVector( 0, 0, 0 );
  }

  Vector<3> v3CamPos = mse3CfW.inverse().get_translation();
  
  Vector<3> v3Best;
  MapPoint* pBest = NULL;
  double dBest = -999999.8;
  for( size_t i = 0; i < mpMap->vpPoints.size(); i++ )
  {
    MapPoint &p = *(mpMap->vpPoints[i]);
    if( p.nSourceLevel != 0 ) {
      continue;
    }
    Vector<3> v3RayToPoint = p.v3WorldPos - v3CamPos;
    normalize(v3RayToPoint);
    if( (v3RayToPoint * v3RayDirnW) > dBest )
    {
      dBest = v3RayToPoint * v3RayDirnW;
      pBest = &p;
    }
  }

  return pBest->v3WorldPos;
}


/**
 * Draw a coordinate frame
 */
void ModelsGame::_DrawAxes()
{
  float w = 0;
  glGetFloatv(GL_LINE_WIDTH, &w);

  glLineWidth(5);
  const float len = 0.5;

  //draw axis
  glBegin(GL_LINES);
  glColor4f(1.0f,0.0f,0.0f,1.0f);
  glVertex3f(0.0f,0.0f,0.0f);
  glVertex3f(len,0.0f,0.0f);

  glColor4f(0.0f,1.0f,0.0f,1.0f);
  glVertex3f(0.0f,0.0f,0.0f);
  glVertex3f(0.0f,len,0.0f);

  glColor4f(0.0f,0.0f,1.0f,1.0f);
  glVertex3f(0.0f,0.0f,0.0f);
  glVertex3f(0.0f,0.0f,len);

  glEnd();

  glLineWidth(w);
}



/**
 * Draw the selection target
 */
void ModelsGame::_DrawSelectedTarget()
{
  glPushMatrix();

  glEnable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_NORMALIZE);
  glEnable(GL_TEXTURE_2D);
  glShadeModel(GL_SMOOTH);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glDisable(GL_LIGHTING);
  
  //draw the selected icon
  glColor4f(1,1,1,1);
  glBindTexture ( GL_TEXTURE_2D , mnTargetTex );
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glBegin(GL_QUADS);
  glTexCoord2f(0.0,0.0);  glVertex3f( 0.5,  0.5, 0.0);
  glTexCoord2f(1.0,0.0);  glVertex3f( 0.5, -0.5, 0.0);
  glTexCoord2f(1.0,1.0);  glVertex3f(-0.5, -0.5, 0.0);
  glTexCoord2f(0.0,1.0);  glVertex3f(-0.5,  0.5, 0.0);
  glEnd();
  
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_DEPTH_TEST);
  
  glPopMatrix();
}


/**
 * Save the game to disk
 * @param sMapPath Path to the map
 * @return The file name
 */
std::string ModelsGame::Save(std::string sMapPath)
{
  std::string sFileName = "ModelsGame.xml";
  std::string sFilePath = sMapPath + "/" + sFileName;
  
  mData.Save( sFilePath );
  
  return sFileName;
}

/**
 * Load a game from disk
 * @param sDataPath The path to the game file.
 */
void ModelsGame::Load(std::string sDataPath)
{
  mData.Load( sDataPath );
}



}


