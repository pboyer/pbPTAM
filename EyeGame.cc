// Copyright 2009 Isis Innovation Limited
#include "EyeGame.h"


namespace PTAMM {

using namespace CVD;

/**
 * Constructor
 */
EyeGame::EyeGame()
  : Game( "Blocks" ), mvDivs(20), mwDivs(20), drawFrame(false)
{
    mvbpBoxes = vector<Box>();
    mvbpBoxesWithPts = vector<Box>();

    mdFieldWidth = 5;
    mdFieldHeight = 5;
    mbInitialised = false;
    mp1 = makeVector(0,0,0);
    mp2 = makeVector(0,0,0);
    mp3 = makeVector(0,0,0);
    
    wp1 = makeVector(0,0,0);
    wp2 = makeVector(0,0,0);
    wp3 = makeVector(0,0,0);
    distanceField = vector< vector<double> >(mwDivs);
    fieldCoordinates = vector< vector< Vector<3> > >(mvDivs);
    mdFieldMultiplier = 1.0;
    mdScale = 1.0;
    mse3MfromW = SE3<>();
    mnFrameCounter = 0;
    mvTranslation = makeVector(0,0,0);
    mdNoisePower = 0.9;
    mdMinHeight = 0.3;
    mdMasterHeight = 0.4;
    mse3ViewerFromWorld = 
    SE3<>::exp(makeVector(0,0,2,0,0,0)) * SE3<>::exp(makeVector(0,0,0,0.8 * M_PI,0,0));
    mv3MassCenter = Zeros;
};
    
void EyeGame::HandleClick(Vector<2> v2VidCoords, Vector<2> v2UFB, Vector<3> v3RayDirnW, Vector<2> v2Plane, int nButton)
{
    if (!drawFrame)
        _HandleOrient(v3RayDirnW, v2Plane, nButton);  
}
    
void EyeGame::HandleKeyPress( std::string sKey ) {
    cout << "Unhandled key press: " << sKey << endl;
    if (strcmp ( sKey.c_str(), "j") ) {
        cout << "K captured" << endl;
        drawFrame = !drawFrame;
    }
}  

vector< Vector<3> > getPtsInBox(Box b, vector< Vector<3> > ptsToContain) {
    
    vector< Vector<3> > containedPts;
    
    for (int i = 0; i < ptsToContain.size(); i++) {
        Vector<3> p = ptsToContain[i];
        if ( b.Contains(p) ) 
            containedPts.push_back(p);
    }
    
    return containedPts;
    
}

void Octree(Box rootBox, vector< Vector<3> > pointsInside, vector<Box>& finalBoxes, vector<Box>& boxesContainingPts, int maxDepth, int currentDepth) {
    
    if (pointsInside.size() <= 1 || currentDepth == maxDepth) {
        //cout << "adding root box" << endl;
        boxesContainingPts.push_back(rootBox);
        if (currentDepth == 0) 
            finalBoxes.push_back(rootBox);
        return;
    }
    
    if (currentDepth == 0) 
        finalBoxes.push_back(rootBox);
    
    vector<Box> babyBoxes = rootBox.GetBabies();
    
    for (int i = 0; i < babyBoxes.size(); i++) {
        Box babyBox = babyBoxes[i];
        vector< Vector<3> > ptsInside = getPtsInBox(babyBox, pointsInside);
        Octree(babyBox, ptsInside, finalBoxes, boxesContainingPts, maxDepth, currentDepth + 1);
        finalBoxes.push_back(babyBox);
    }
    
}

    
void EyeGame::_UpdateBoxes() {
    
    // get the bounding box on all of the points
    // go through all points, z value bigger than current
    // need zmax, zmin, etc

    double xmin = 1e12;
    double ymin = 1e12;
    double zmin = 1e12;
    double xmax = -1e12;
    double ymax = -1e12;
    double zmax = -1e12;
    
    std::vector<MapPoint*> vpPoints = mpMap->vpPoints;
    std::vector< Vector<3> > pointsInside = std::vector< Vector<3> >();
    
    Vector<3> camPos = mse3CfW.inverse().get_translation();
    pointsInside.push_back(camPos);
    
    double x = camPos[0];
    double y = camPos[1];
    double z = camPos[2];
    
    //cout << "Num pts: " << vpPoints.size() << endl;
    
    for (int i = 0; i < vpPoints.size(); i++) {

        if (abs(z) < 1e3 && abs(y) < 1e3 && abs(x) < 1e3) {
        
            Vector<3> p = makeVector(x,y,z);
            pointsInside.push_back(p);
            
            if ( x < xmin) xmin = x;
            if ( y < ymin) ymin = y;
            if ( z < zmin) zmin = z;
            if ( x > xmax) xmax = x;
            if ( y > ymax) ymax = y;
            if ( z > zmax) zmax = z;
       
        }
        
        x = vpPoints[i]->v3WorldPos[0];
        y = vpPoints[i]->v3WorldPos[1];
        z = vpPoints[i]->v3WorldPos[2];

    }
    
    Vector<3> tr = makeVector(xmax, ymax, zmax);
    Vector<3> bl = makeVector(xmin, ymin, zmin);

//    cout << "bl: " << bl << endl;
//    cout << "tr: " << tr << endl;
    
    Box rootBox(bl, tr);
    
    mvbpBoxes = vector<Box>();
    mvbpBoxesWithPts = vector<Box>();
    
    Octree(rootBox, pointsInside, mvbpBoxes, mvbpBoxesWithPts, 8, 0);

    
}

void EyeGame::_DrawBoxes() {
    
    glLineWidth(1);
    
    for (unsigned i = 0; i < mvbpBoxesWithPts.size(); i++) {
        
        vector< Vector<3> > corners = mvbpBoxesWithPts[i].GetCorners();
//        cout << "num corners: " << corners.size() << endl;
//        cout << "Corner 0: " << corners[0] << endl;
//        cout << "Corner 1: " << corners[1] << endl;
//        cout << "Corner 2: " << corners[2] << endl;
//        cout << "Corner 3: " << corners[3] << endl;
//        cout << "Corner 4: " << corners[4] << endl;
//        cout << "Corner 5: " << corners[5] << endl;
//        cout << "Corner 6: " << corners[6] << endl;
//        cout << "Corner 7: " << corners[7] << endl;
        
//        cout << "X: " << mvbpBoxes[i].diff << endl;
        
        if (corners.size() == 8) {
            glDisable(GL_LIGHTING);
            glColor4f(0.6,0.5,0.5,0.4);
            
            // bottom face
            glBegin(GL_LINES);
            // top
            glVertex( corners[0] );
            glVertex( corners[1] );
            
            glVertex( corners[1] );
            glVertex( corners[2] );
            
            glVertex( corners[2] );
            glVertex( corners[3] );
            
            glVertex( corners[3] );
            glVertex( corners[0] );
            
            glVertex( corners[4] );
            glVertex( corners[0] );
            
            glVertex( corners[1] );
            glVertex( corners[5] );
            
            glVertex( corners[6] );
            glVertex( corners[2] );

            glVertex( corners[3] );
            glVertex( corners[7] );
            
            glVertex( corners[4] );
            glVertex( corners[5] );
            
            glVertex( corners[5] );
            glVertex( corners[6] );
            
            glVertex( corners[6] );
            glVertex( corners[7] );
            
            glVertex( corners[7] );
            glVertex( corners[4] );
            
            glEnd();
            
            glEnable(GL_LIGHTING);
            glEnable(GL_BLEND);
            
            GLfloat af[4]; 
            
            double a = ( 0.5 + 0.5 * sin( 2 * M_PI * ( (mnFrameCounter%800)/((double) 800) )  ) );
            double b = ( 0.5 + 0.5 * cos( 2 * M_PI * ( (mnFrameCounter%1600)/((double) 1600) ) ) );
            double c = ( 0.5 + 0.5 * sin( 2 * M_PI * ( corners[3][0] * (mnFrameCounter%3200)/((double) 3200) )  ) );
            
            af[0]= 0.5 + 0.5 * noise( a, b, c);
            af[1]= 0.5 + 0.5 * noise( b, c, a);
            af[2]= 0.5 + 0.5 * noise( c, b, a);
            
            af[3]= 0.2f;
            
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, af);
            
            glBegin(GL_QUADS);
            
            glVertex( corners[3] );
            glVertex( corners[2] );
            glVertex( corners[1] );
            glVertex( corners[0] );
            
            glVertex( corners[0] );
            glVertex( corners[1] );
            glVertex( corners[2] );
            glVertex( corners[3] );
            
            glVertex( corners[4] );
            glVertex( corners[5] );
            glVertex( corners[6] );
            glVertex( corners[7] );
            
            glVertex( corners[0] );
            glVertex( corners[1] );
            glVertex( corners[5] );
            glVertex( corners[4] );
            
            glVertex( corners[1] );
            glVertex( corners[2] );
            glVertex( corners[6] );
            glVertex( corners[5] );
            
            glVertex( corners[2] );
            glVertex( corners[3] );
            glVertex( corners[7] );
            glVertex( corners[6] );
            
            glVertex( corners[3] );
            glVertex( corners[0] );
            glVertex( corners[4] );
            glVertex( corners[7] );
            
            glEnd();
            
            
        }
    }
    
  

}
    
/**
 * Draw the game
 * @param glWindow the gl window
 * @param map the current map
 * @param se3CfromW the location of the camera
 */
void EyeGame::Draw3D( GLWindow2 &glWindow, Map &map, SE3<> se3CfromW)
{
    if (!mbInitialised)
        Init();
    
    mpMap = &map;
    mse3CfW = se3CfromW; 
    
    //this->_UpdateBoxes();
    
    this->_UpdateDistanceField();

    mnFrameCounter++;
    
    if (!drawFrame) {

        this->_DrawOrientPoints();
        
        glEnable(GL_BLEND);
        glShadeModel(GL_SMOOTH);                        // Enable Smooth Shading
        //glClearDepth(1.0f);                             // Depth Buffer Setup
        glEnable(GL_DEPTH_TEST);                        // Enables Depth Testing
        glDepthFunc(GL_LEQUAL);                         // The Type Of Depth Testing To Do
        glDisable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);

        GLfloat af[4]; 
        af[0]=0.5; af[1]=0.5; af[2]=0.5; af[3]=1.0;
        glLightfv(GL_LIGHT0, GL_DIFFUSE, af);
        glLightfv(GL_LIGHT0, GL_AMBIENT, af);
        glLightfv(GL_LIGHT0, GL_SPECULAR, af);
        af[0]=1.0; af[1]=0.0; af[2]=1.0; af[3]=1.0;
        glLightfv(GL_LIGHT0, GL_POSITION, af);
               
        glMatrixMode(GL_MODELVIEW);
        
        glPushMatrix();
        
        glMultMatrix( mse3MfromW );
        glScaled( mdScale, mdScale, mdScale );

        //this->_DrawBoxes();
        this->_DrawBlocks();
        
        glPopMatrix();

        glDisable(GL_LIGHTING);

        glLoadIdentity();

        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        
    } else {
        
        {
            std::pair<TooN::Vector<6>, TooN::Vector<6> >  pv6 = glWindow.GetMousePoseUpdate();
            SE3<> se3CamFromMC;
            se3CamFromMC.get_translation() = mse3ViewerFromWorld * mv3MassCenter;
            mse3ViewerFromWorld = SE3<>::exp(pv6.first) * 
            se3CamFromMC * SE3<>::exp(pv6.second) * se3CamFromMC.inverse() * mse3ViewerFromWorld;
        }
        
        //glWindow.SetupViewport();
        
        glClearColor(0,0,0,0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_BLEND);
        glShadeModel(GL_SMOOTH);                        // Enable Smooth Shading
        glClearDepth(1.0f);                             // Depth Buffer Setup
        glEnable(GL_DEPTH_TEST);                        // Enables Depth Testing
        glDepthFunc(GL_LEQUAL);                         // The Type Of Depth Testing To Do
        glDisable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        
        //this->_DrawOrientPoints();

        GLfloat af[4]; 
        af[0]=0.5; af[1]=0.5; af[2]=0.5; af[3]=1.0;
        glLightfv(GL_LIGHT0, GL_DIFFUSE, af);
        glLightfv(GL_LIGHT0, GL_AMBIENT, af);
        glLightfv(GL_LIGHT0, GL_SPECULAR, af);
        af[0]=1.0; af[1]=0.0; af[2]=1.0; af[3]=1.0;
        glLightfv(GL_LIGHT0, GL_POSITION, af);
        
        this->SetupFrustum();
        this->SetupModelView();
        
        this->_DrawGrid();
        this->_DrawBlocks();
        
        glDisable(GL_LIGHTING);
        
        this->DrawCamera(mse3CfW);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        
        
    }
    
    
};
    
/**
 * Draw the camera / keyframe
 * @param se3CfromW Camera's location in thr world
 * @param bSmall Draw the small camera (keyframe)
 */
void EyeGame::DrawCamera(SE3<> se3CfromW, bool bSmall)
{
    
    SetupModelView(se3CfromW.inverse());
    SetupFrustum();
    
    if(bSmall)
        glLineWidth(1);
    else
        glLineWidth(7);
    
    glBegin(GL_LINES);
    glColor4f(1,0,0,1);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.1f, 0.0f, 0.0f);
    glColor4f(0,1,0,1);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.1f, 0.0f);
    glColor4f(1,1,1,1);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 20.0f);
    glEnd();
    
    if(!bSmall)
    {
        glLineWidth(1);
        glColor3f(0.5,0.5,0.5);
        SetupModelView();
        Vector<2> v2CamPosXY = se3CfromW.inverse().get_translation().slice<0,2>();
        glBegin(GL_LINES);
        glColor3f(1,1,1);
        glVertex2d(v2CamPosXY[0] - 0.04, v2CamPosXY[1] + 0.04);
        glVertex2d(v2CamPosXY[0] + 0.04, v2CamPosXY[1] - 0.04);
        glVertex2d(v2CamPosXY[0] - 0.04, v2CamPosXY[1] - 0.04);
        glVertex2d(v2CamPosXY[0] + 0.04, v2CamPosXY[1] + 0.04);
        glEnd();
    }
    
}

//    List<Box> boxesInRecursion = new List<Box>();
//    List<Box> boxesContainingPts = new List<Box>();
//    Octree(rootBox, points, boxesInRecursion, boxesContainingPts, recursionLimits, 0);
//    
//    boxesFinal = boxesInRecursion;
//    boxesWithPts = boxesContainingPts;



/**
 * set up the viewer frustrum
 */
void EyeGame::SetupFrustum()
{
    glMatrixMode(GL_PROJECTION);  
    glLoadIdentity();
    double zNear = 0.03;
    glFrustum(-zNear, zNear, 0.75*zNear,-0.75*zNear,zNear,50);
    glScalef(1,1,-1);
    return;
}

/**
 * Set up the model view.
 * @param se3WorldFromCurrent se3 that converts from
 * world frame to camera frame
 */
void EyeGame::SetupModelView(SE3<> se3WorldFromCurrent)
{
    glMatrixMode(GL_MODELVIEW);  
    glLoadIdentity();
    glMultMatrix(mse3ViewerFromWorld * se3WorldFromCurrent);
    return;
}
    
/**
 * Draw the Grid
 */
void EyeGame::_DrawGrid()
{
    glLineWidth(1);
    
    glBegin(GL_LINES);
    
    // Draw a larger grid around the outside..
    double dGridInterval = 0.1;
    
    double dMin = -100.0 * dGridInterval;
    double dMax =  100.0 * dGridInterval;
    
    for(int x=-10;x<=10;x+=1)
    {
        if(x==0)
            glColor3f(1,1,1);
        else
            glColor3f(0.3,0.3,0.3);
        glVertex3d((double)x * 10 * dGridInterval, dMin, 0.0);
        glVertex3d((double)x * 10 * dGridInterval, dMax, 0.0);
    }
    for(int y=-10;y<=10;y+=1)
    {
        if(y==0)
            glColor3f(1,1,1);
        else
            glColor3f(0.3,0.3,0.3);
        glVertex3d(dMin, (double)y * 10 *  dGridInterval, 0.0);
        glVertex3d(dMax, (double)y * 10 * dGridInterval, 0.0);
    }
    
    glEnd();
    
    glBegin(GL_LINES);
    dMin = -10.0 * dGridInterval;
    dMax =  10.0 * dGridInterval;
    
    for(int x=-10;x<=10;x++)
    {
        if(x==0)
            glColor3f(1,1,1);
        else
            glColor3f(0.5,0.5,0.5);
        
        glVertex3d((double)x * dGridInterval, dMin, 0.0);
        glVertex3d((double)x * dGridInterval, dMax, 0.0);
    }
    for(int y=-10;y<=10;y++)
    {
        if(y==0)
            glColor3f(1,1,1);
        else
            glColor3f(0.5,0.5,0.5);
        glVertex3d(dMin, (double)y * dGridInterval, 0.0);
        glVertex3d(dMax, (double)y * dGridInterval, 0.0);
    }
    
    glColor3f(1,0,0);
    glVertex3d(0,0,0);
    glVertex3d(1,0,0);
    glColor3f(0,1,0);
    glVertex3d(0,0,0);
    glVertex3d(0,1,0);
    glColor3f(1,1,1);
    glVertex3d(0,0,0);
    glVertex3d(0,0,1);
    glEnd();
    
    //   glColor3f(0.8,0.8,0.8);
    //   glRasterPos3f(1.1,0,0);
    //   mGLWindow.PrintString("x");
    //   glRasterPos3f(0,1.1,0);
    //   mGLWindow.PrintString("y");
    //   glRasterPos3f(0,0,1.1);
    //   mGLWindow.PrintString("z");
}


/**
 * Reset the game
 */
void EyeGame::Reset()
{

    // get rid of all your orient points
    mdFieldMultiplier = 5.0;
    mdScale = 0.5;
    mse3MfromW = SE3<>();
    mnFrameCounter = 0;

    
};


/**
 * Initialize the game
 */
void EyeGame::Init()
{
    if(mbInitialised) return;
    mbInitialised = true;

    planeIntersect1 = false;
    planeIntersect2 = false;
    planeIntersect3 = false;

    double fieldHeight = mdFieldHeight;
    double fieldWidth = mdFieldWidth;
    
    for (int i = 0; i < mvDivs; i++) {
        
        double y = ((double) -fieldHeight)/2 + fieldWidth * ((double) i) /((double)mvDivs - 1);
        
        distanceField[i] =  vector<double>( mwDivs ) ;
        fieldCoordinates[i] = vector< Vector<3> >( mwDivs );
        
        for (int j = 0; j < mwDivs; j++) {
            
            double x = ((double) -fieldHeight)/2 + fieldHeight * ((double) j) /((double)mwDivs - 1);
            Vector<3> position = makeVector(x, y, 0.0);
            
            fieldCoordinates[i][j] = position;
            distanceField[i][j] = 0.0;

        }
        
    }

    Reset();
    

};
    
 
void EyeGame::_DrawOrientPoints()
{    
    
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

};
    
double EyeGame::noise( double x, double y, double z )
{
    int X = (int)floor(x) & 255;                  // FIND UNIT CUBE THAT
    int Y = (int)floor(y) & 255;                  // CONTAINS POINT.
    int Z = (int)floor(z) & 255;
    x -= floor(x);                                // FIND RELATIVE X,Y,Z
    y -= floor(y);                                // OF POINT IN CUBE.
    z -= floor(z);
    double u = fade(x);                                // COMPUTE FADE CURVES
    double v = fade(y);                                // FOR EACH OF X,Y,Z.
    double w = fade(z);
    int A = p[X  ]+Y; int AA = p[A]+Z; int AB = p[A+1]+Z;      // HASH COORDINATES OF
    int B = p[X+1]+Y; int BA = p[B]+Z; int BB = p[B+1]+Z;      // THE 8 CUBE CORNERS,
    return lerp(w, lerp(v, lerp(u, grad(p[AA  ], x  , y  , z   ),  // AND ADD
                                grad(p[BA  ], x-1, y  , z   )), // BLENDED
                        lerp(u, grad(p[AB  ], x  , y-1, z   ),  // RESULTS
                             grad(p[BB  ], x-1, y-1, z   ))),// FROM  8
                lerp(v, lerp(u, grad(p[AA+1], x  , y  , z-1 ),  // CORNERS
                             grad(p[BA+1], x-1, y  , z-1 )), // OF CUBE
                     lerp(u, grad(p[AB+1], x  , y-1, z-1 ),
                          grad(p[BB+1], x-1, y-1, z-1 ))));
};


double EyeGame::octaveNoise( const Vector<3>& pt, int octaves )
{
    double answer = 0;
    for (int i = 0; i < octaves; i++)
    {
        float tmp = pow(2.0f,i);
        answer += noise(tmp*pt[0],tmp*pt[1],tmp*pt[2]) / float(tmp);
    }
    return answer;
};

double EyeGame::fade( double t )
{ 
    return t * t * t * (t * (t * 6 - 15) + 10);
}

double EyeGame::lerp( double t, double a, double b )
{ 
    return a + t * (b - a);
}

double EyeGame::grad( int hash, double x, double y, double z )
{
    int h = hash & 15;                      // CONVERT LO 4 BITS OF HASH CODE
    double u = h<8 ? x : y;                 // INTO 12 GRADIENT DIRECTIONS.
    double v = h<4 ? y : h==12||h==14 ? x : z;
    return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v);
    
}
    
/**
 * Updates set of shifting vertices based on position of person
 */
void EyeGame::_UpdateDistanceField()
{
    // draw a rectangular array of blocks that shift in scale with time
    
    Vector<3> camPos = mse3CfW.inverse().get_translation();
    Vector<3> camDir = mse3CfW.inverse().get_rotation() * makeVector(0,0,1);
    normalize(camDir);
    //cout << camDir << endl;

//    cout << camPos << endl;
////    cout << mvDivs << " " << mwDivs << endl;
//    
//    Vector<3> pos = fieldCoordinates[0][0];
//    cout << pos << endl;
//    cout << mse3MfromW * (mdScale * pos) << endl;
    
    for (int i = 0; i < mvDivs; i++) {
        for (int j = 0; j < mwDivs; j++) {
            
            Vector<3> pos = fieldCoordinates[i][j]; 
            pos = mse3MfromW * pos;
            pos = mdScale * pos;
            
            // project point onto line and get norm from that
            Vector<3> toPoint = pos - camPos;
            Vector<3> projectedPoint = (toPoint * camDir) * camDir;
                
            double distToViewLine = norm(projectedPoint - pos);
            double val = norm(toPoint);
            
            distanceField[i][j] = mdMasterHeight * ( mdMinHeight * val * distToViewLine + mdNoisePower * noise( fieldCoordinates[i][j][0] * val + 0.03 * sin(2 * M_PI * ( (mnFrameCounter%800)/((double) 800) ) ), fieldCoordinates[i][j][1] * val + 0.03 * sin(((double) mnFrameCounter)/800), val  ) );
            
        }
    }
    
};
    
void _DrawBlock(double height) {
    
    glBegin(GL_QUADS);	
    // Front Face
    glNormal3f( 0.0f, 0.0f, 1.0f);                  // Normal Pointing Towards Viewer
    glTexCoord2f(0.0f, 0.0f); glVertex3f(-1.0f, -1.0f,  height);  // Point 1 (Front)
    glTexCoord2f(1.0f, 0.0f); glVertex3f( 1.0f, -1.0f,  height);  // Point 2 (Front)
    glTexCoord2f(1.0f, 1.0f); glVertex3f( 1.0f,  1.0f,  height);  // Point 3 (Front)
    glTexCoord2f(0.0f, 1.0f); glVertex3f(-1.0f,  1.0f,  height);  // Point 4 (Front)
    // Back Face
    glNormal3f( 0.0f, 0.0f,-1.0f);                  // Normal Pointing Away From Viewer
    glTexCoord2f(1.0f, 0.0f); glVertex3f(-1.0f, -1.0f, -1.0f);  // Point 1 (Back)
    glTexCoord2f(1.0f, 1.0f); glVertex3f(-1.0f,  1.0f, -1.0f);  // Point 2 (Back)
    glTexCoord2f(0.0f, 1.0f); glVertex3f( 1.0f,  1.0f, -1.0f);  // Point 3 (Back)
    glTexCoord2f(0.0f, 0.0f); glVertex3f( 1.0f, -1.0f, -1.0f);  // Point 4 (Back)
    // Top Face
    glNormal3f( 0.0f, 1.0f, 0.0f);                  // Normal Pointing Up
    glTexCoord2f(0.0f, 1.0f); glVertex3f(-1.0f,  1.0f, -1.0f);  // Point 1 (Top)
    glTexCoord2f(0.0f, 0.0f); glVertex3f(-1.0f,  1.0f,  height);  // Point 2 (Top)
    glTexCoord2f(1.0f, 0.0f); glVertex3f( 1.0f,  1.0f,  height);  // Point 3 (Top)
    glTexCoord2f(1.0f, 1.0f); glVertex3f( 1.0f,  1.0f, -1.0f);  // Point 4 (Top)
    // Bottom Face
    glNormal3f( 0.0f,-1.0f, 0.0f);                  // Normal Pointing Down
    glTexCoord2f(1.0f, 1.0f); glVertex3f(-1.0f, -1.0f, -1.0f);  // Point 1 (Bottom)
    glTexCoord2f(0.0f, 1.0f); glVertex3f( 1.0f, -1.0f, -1.0f);  // Point 2 (Bottom)
    glTexCoord2f(0.0f, 0.0f); glVertex3f( 1.0f, -1.0f,  height);  // Point 3 (Bottom)
    glTexCoord2f(1.0f, 0.0f); glVertex3f(-1.0f, -1.0f,  height);  // Point 4 (Bottom)
    // Right face
    glNormal3f( 1.0f, 0.0f, 0.0f);                  // Normal Pointing Right
    glTexCoord2f(1.0f, 0.0f); glVertex3f( 1.0f, -1.0f, -1.0f);  // Point 1 (Right)
    glTexCoord2f(1.0f, 1.0f); glVertex3f( 1.0f,  1.0f, -1.0f);  // Point 2 (Right)
    glTexCoord2f(0.0f, 1.0f); glVertex3f( 1.0f,  1.0f,  height);  // Point 3 (Right)
    glTexCoord2f(0.0f, 0.0f); glVertex3f( 1.0f, -1.0f,  height);  // Point 4 (Right)
    // Left Face
    glNormal3f(-1.0f, 0.0f, 0.0f);                  // Normal Pointing Left
    glTexCoord2f(0.0f, 0.0f); glVertex3f(-1.0f, -1.0f, -1.0f);  // Point 1 (Left)
    glTexCoord2f(1.0f, 0.0f); glVertex3f(-1.0f, -1.0f,  height);  // Point 2 (Left)
    glTexCoord2f(1.0f, 1.0f); glVertex3f(-1.0f,  1.0f,  height);  // Point 3 (Left)
    glTexCoord2f(0.0f, 1.0f); glVertex3f(-1.0f,  1.0f, -1.0f);  // Point 4 (Left)
    glEnd();										// Done Drawing The Quad
    
    
    glDisable(GL_LIGHTING);
    glLineWidth(1);
    
    glColor4f(1.0f, 1.0f, 1.0f, 0.8f);
    
    glBegin(GL_LINES);									// Draw A Quad
    glVertex3f( 1.0f, 1.0f,-1.0f);					// Top Right Of The Quad (Top)
    glVertex3f(-1.0f, 1.0f,-1.0f);					// Top Left Of The Quad (Top)
    glVertex3f(-1.0f, 1.0f, height);					// Bottom Left Of The Quad (Top)
    glVertex3f( 1.0f, 1.0f, height);	
    glEnd();
    
    glBegin(GL_LINES);	
    glVertex3f( 1.0f,-1.0f, height);					// Top Right Of The Quad (Bottom)
    glVertex3f(-1.0f,-1.0f, height);					// Top Left Of The Quad (Bottom)
    glVertex3f(-1.0f,-1.0f,-1.0f);					// Bottom Left Of The Quad (Bottom)
    glVertex3f( 1.0f,-1.0f,-1.0f);					// Bottom Right Of The Quad (Bottom)
    glEnd();
    
    glBegin(GL_LINES);	
    //glColor3f(1.0f,0.0f,0.0f);						// Set The Color To Red
    glVertex3f( 1.0f, 1.0f, height);					// Top Right Of The Quad (Front)
    glVertex3f(-1.0f, 1.0f, height);					// Top Left Of The Quad (Front)
    glVertex3f(-1.0f,-1.0f, height);					// Bottom Left Of The Quad (Front)
    glVertex3f( 1.0f,-1.0f, height);					// Bottom Right Of The Quad (Front)
    glEnd();
    
    glBegin(GL_LINES);	
    glVertex3f( 1.0f,-1.0f,-1.0f);					// Top Right Of The Quad (Back)
    glVertex3f(-1.0f,-1.0f,-1.0f);					// Top Left Of The Quad (Back)
    glVertex3f(-1.0f, 1.0f,-1.0f);					// Bottom Left Of The Quad (Back)
    glVertex3f( 1.0f, 1.0f,-1.0f);					// Bottom Right Of The Quad (Back)
    glEnd();

    glVertex3f(-1.0f, 1.0f, height);					// Top Right Of The Quad (Left)
    glVertex3f(-1.0f, 1.0f,-1.0f);					// Top Left Of The Quad (Left)
    glVertex3f(-1.0f,-1.0f,-1.0f);					// Bottom Left Of The Quad (Left)
    glVertex3f(-1.0f,-1.0f, height);					// Bottom Right Of The Quad (Left)
    
    glEnd();
    glBegin(GL_LINES);	
    glVertex3f( 1.0f, 1.0f,-1.0f);					// Top Right Of The Quad (Right)
    glVertex3f( 1.0f, 1.0f, height);					// Top Left Of The Quad (Right)
    glVertex3f( 1.0f,-1.0f, height);					// Bottom Left Of The Quad (Right)
    glVertex3f( 1.0f,-1.0f,-1.0f);					// Bottom Right Of The Quad (Right)
    glEnd();
    
    glEnable(GL_LIGHTING);
    
};

    
/**
 * Draw shifty blocks
 */
void EyeGame::_DrawBlocks()
{
    // draw a rectangular array of blocks that shift in scale with time
    double blockScaleW = mdFieldWidth / (2 * (mwDivs - 1));
    double blockScaleV = mdFieldHeight / (2 * (mvDivs - 1));
    
    for (int i = 0; i < mvDivs; i++) {
        for (int j = 0; j < mwDivs; j++) {
            
            GLfloat af[4]; 
            
            double a = ( 0.5 + sin( 2 * M_PI * ( (mnFrameCounter%800)/((double) 800) )  ) );
            double b = ( 0.5 + cos( 2 * M_PI * ( (mnFrameCounter%1600)/((double) 1600) ) ) );
            double c = ( 0.5 + sin( distanceField[i][j] * 2 * M_PI * ( (mnFrameCounter%3200)/((double) 3200) )  ) );
            
            af[0]= 0.5 + 0.5 * noise( a, b, c);
            af[1]= 0.5 + 0.5 * noise( b, c, a);
            af[2]= 0.5 + 0.5 * noise( c, b, a);
            
            af[3]= 0.9f;
            
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, af);
            
//            af[0]=1.0; af[1]=1.0; af[2]=1.0f; af[3]=1.0;
//            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, af);
            //glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0);
            
            Vector<3> translation = fieldCoordinates[i][j];
            double h = distanceField[i][j];

            glTranslatef( translation[0], translation[1], translation[2] );	  
            // scale locally
            glScaled(blockScaleW, blockScaleV, 1.0);
            _DrawBlock(h);
            glScaled(1/blockScaleW, 1/blockScaleV, 1.0);
            
            glTranslatef( -translation[0], -translation[1], -translation[2] );	
            
        }
        
    }
    
    // blocks transparency will be augmented
    // blocks color will be augmented
    
};
    
/**
 * Draw shifty blocks
 */
void EyeGame::_DrawWireframeWalls()
{

    
};
    
void EyeGame::_DrawWall()
{

    
};
    
    

    
void EyeGame::_HandleOrient(Vector<3> v3RayDirnW, Vector<2> v2Plane, int nButton) {
    
    
    if (!planeIntersect1) {
        cout << "Obtained target point 1"  << endl;
        
        wp1 = this->_SelectClosestMapPoint( v3RayDirnW );
           
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
                
         if (blockVisible) {
            
            // The scale comes from the ratio between the lengths of the first two origin points
            // and the first two target points
             
             mp1 = mse3MfromW * makeVector(0, 0, 0);
             mp2 = mse3MfromW * makeVector(mdFieldWidth * mdScale, 0, 0);
             mp3 = mse3MfromW * makeVector(0, mdFieldHeight * mdScale, 0);
             
             cout << "M to W mat" << endl;
             cout << mse3MfromW << endl;
             
             cout << "model points" << endl;
             cout << mp1 << endl;
             cout << mp2 << endl;
             cout << mp3 << endl;
             
//             mp1 = makeVector(0, 0, 0);
//             mp2 = makeVector(1, 0, 0);
//             mp3 = makeVector(0, 1, 0);
             
             cout << "world points" << endl;
             cout << wp1 << endl;
             cout << wp2 << endl;
             cout << wp3 << endl;
             
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
        
// this should be the models current translation
        
            Vector<3> ca = this->mse3MfromW.get_translation();
            
            Vector<3> amca = mp1 - ca;
            Vector<3> bmcb = rotation*(amca);
            bmcb = scaleRatio * bmcb;
            Vector<3> bmca = wp1 - ca;
            
            Vector<3> newPos = bmca - bmcb + ca;
            
// update here

            this->mdScale = this->mdScale * scaleRatio;
//            // add the new translation
            mse3MfromW.get_translation() = newPos;
//            // this rotation is the last one we will perform
            this->mse3MfromW.get_rotation() = rotation * this->mse3MfromW.get_rotation();  
            
        }
        
        planeIntersect1 = false;
        planeIntersect2 = false;
        planeIntersect3 = false;
        return;
    }
};
    
/**
 * Select the closest map point to the ray from click
 */
Vector<3> EyeGame::_SelectClosestMapPoint( Vector<3> v3RayDirnW )
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
};
    
    // permutation
    int EyeGame::p[512] = 
    { 151,160,137, 91, 90, 15,131, 13,201, 95, 96, 53,194,233,  7,225,
    140, 36,103, 30, 69,142,  8, 99, 37,240, 21, 10, 23,190,  6,148,
    247,120,234, 75,  0, 26,197, 62 ,94,252,219,203,117, 35, 11, 32,
    57,177, 33, 88,237,149, 56, 87,174, 20,125,136,171,168, 68,175,
    74,165, 71,134,139, 48, 27,166, 77,146,158,231, 83,111,229,122,
    60,211,133,230,220,105, 92, 41, 55, 46,245, 40,244,102,143, 54,
    65, 25, 63,161,  1,216, 80, 73,209, 76,132,187,208, 89, 18,169,
    200,196,135,130,116,188,159, 86,164,100,109,198,173,186,  3, 64,
    52,217,226,250,124,123,  5,202, 38,147,118,126,255, 82, 85,212,
    207,206, 59,227, 47, 16, 58, 17,182,189, 28, 42,223,183,170,213,
    119,248,152,  2, 44,154,163, 70,221,153,101,155,167, 43,172,  9,
    129, 22, 39,253, 19, 98,108,110, 79,113,224,232,178,185,112,104,
    218,246, 97,228,251, 34,242,193,238,210,144, 12,191,179,162,241,
    81, 51,145,235,249, 14,239,107, 49,192,214, 31,181,199,106,157,
    184, 84,204,176,115,121, 50, 45,127,  4,150,254,138,236,205, 93,
    222,114, 67, 29, 24, 72,243,141,128,195, 78, 66,215, 61,156,180,
    
    // repeat
    151,160,137, 91, 90, 15,131, 13,201, 95, 96, 53,194,233,  7,225,
    140, 36,103, 30, 69,142,  8, 99, 37,240, 21, 10, 23,190,  6,148,
    247,120,234, 75,  0, 26,197, 62 ,94,252,219,203,117, 35, 11, 32,
    57,177, 33, 88,237,149, 56, 87,174, 20,125,136,171,168, 68,175,
    74,165, 71,134,139, 48, 27,166, 77,146,158,231, 83,111,229,122,
    60,211,133,230,220,105, 92, 41, 55, 46,245, 40,244,102,143, 54,
    65, 25, 63,161,  1,216, 80, 73,209, 76,132,187,208, 89, 18,169,
    200,196,135,130,116,188,159, 86,164,100,109,198,173,186,  3, 64,
    52,217,226,250,124,123,  5,202, 38,147,118,126,255, 82, 85,212,
    207,206, 59,227, 47, 16, 58, 17,182,189, 28, 42,223,183,170,213,
    119,248,152,  2, 44,154,163, 70,221,153,101,155,167, 43,172,  9,
    129, 22, 39,253, 19, 98,108,110, 79,113,224,232,178,185,112,104,
    218,246, 97,228,251, 34,242,193,238,210,144, 12,191,179,162,241,
    81, 51,145,235,249, 14,239,107, 49,192,214, 31,181,199,106,157,
    184, 84,204,176,115,121, 50, 45,127,  4,150,254,138,236,205, 93,
    222,114, 67, 29, 24, 72,243,141,128,195, 78, 66,215, 61,156,180};
    
    

 
} // END PTAMM Namespace

