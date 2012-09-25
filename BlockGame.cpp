////
////  BlockGame.cpp
////  pbPTAM
////
////  Created by Peter Boyer on 4/5/12.
////  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
////
//
//#include "BlockGame.h"
//
//namespace PTAMM {
//    
//    /**
//     * Constructor
//     */
//    BlockGame::BlockGame()
//    : Game( "BlockGame" )
//    {
//        
//    }
//    
//    /**
//     * Draw the game
//     * @param glWindow the gl window
//     * @param map the current map
//     * @param se3CfromW the location of the camera
//     */
//    void BlockGame::Draw3D( const GLWindow2 &glWindow, Map &map, SE3<> se3CfromW)
//    {
////        Vector<3> v3CameraPos = se3CfromW.inverse().get_translation();
////        
////        if(!mbInitialised) {
////            Init();
////        }
////        
////        mnFrameCounter ++;
////        
////        glDisable(GL_BLEND);
////        glEnable(GL_CULL_FACE);
////        glFrontFace(GL_CW);
////        glEnable(GL_DEPTH_TEST);
////        glDepthFunc(GL_LEQUAL);
////        glEnable(GL_LIGHTING);
////        glEnable(GL_LIGHT0);
////        glEnable(GL_NORMALIZE);
////        glEnable(GL_COLOR_MATERIAL);
////        
////        GLfloat af[4]; 
////        af[0]=0.5; af[1]=0.5; af[2]=0.5; af[3]=1.0;
////        glLightfv(GL_LIGHT0, GL_AMBIENT, af);
////        glLightfv(GL_LIGHT0, GL_DIFFUSE, af);
////        af[0]=1.0; af[1]=0.0; af[2]=1.0; af[3]=0.0;
////        glLightfv(GL_LIGHT0, GL_POSITION, af);
////        af[0]=1.0; af[1]=1.0; af[2]=1.0; af[3]=1.0;
////        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, af);
////        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0);
////        
////        glMatrixMode(GL_MODELVIEW);
////        
////        glDisable(GL_LIGHTING);
////        
////        glLoadIdentity();
//
//        
//    };
//    
//    
//    /**
//     * Reset the game
//     */
//    void BlockGame::Reset()
//    {
//
//    };
//    
//
//    
//    /**
//     * Initialize the game
//     */
//    void BlockGame::Init()
//    {
//        if(mbInitialised) return;
//        mbInitialised = true;
//        
//        // Set up the display list for the eyeball.
//        //mnEyeDisplayList = glGenLists(1);
//        
////        glNewList(mnEyeDisplayList,GL_COMPILE);
////        DrawEye();
////        glEndList();
//
//        Reset();
//    };
//    
//}
//    
//
//
