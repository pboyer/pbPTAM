//
//  Box.h
//  pbPTAM
//
//  Created by Peter Boyer on 4/25/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//
#include "Box.h"

using namespace TooN;
using namespace std;
using namespace PTAMM;

Box::Box(Vector<3> blin, Vector<3> trin) {
    
    bl = blin;
    tr = trin;
    diff = trin - blin;
    
//    cout << "diff:" << diff << endl;
    
    double X = diff[0];
    double Y = diff[1];
    double Z = diff[2];
    
    //        cout << "X:" << X << endl;
    //        cout << "Y:" << Y << endl;
    //        cout << "Z:" << Z << endl;
    
    //        if (X > 1e8 || Y > 1e8 || Z > 1e8)
    //            cout << "BIG VALUE ENCOUNTERED" << endl;
    
    // counter clockwise winding of bottom
    // counter clockwise winding of top
    corners = vector< Vector<3> >(8);
    
    Vector<3> x_z = makeVector(X,0,Z);
    Vector<3> xy_ = makeVector(X,Y,0);
    Vector<3> x__ = makeVector(X,0,0);
    Vector<3> _y_ = makeVector(0,Y,0);
    Vector<3> __z = makeVector(0,0,Z);
    
    corners[0] = bl;
    corners[1] = bl + x__;
    corners[2] = bl + xy_;
    corners[3] = bl + _y_;
    corners[4] = bl + __z;
    corners[5] = bl + x_z;
    corners[6] = tr;
    corners[7] = tr - x__;
    
}
        
//    Box() {
//        bl = makeVector(0,0,0);
//        tr = makeVector(0,0,0); // bl tr points
//        X = 10;
//        Y = 10;
//        Z = 10;  // lengths
//        corners = vector< Vector<3> >(8);
//        
//        Vector<3> x_z = makeVector(X,0,Z);
//        Vector<3> xy_ = makeVector(X,Y,0);
//        Vector<3> x__ = makeVector(X,0,0);
//        Vector<3> _y_ = makeVector(0,Y,0);
//        Vector<3> __z = makeVector(0,0,Z);
//        
//        corners[0] = bl;
//        corners[1] = bl + x__;
//        corners[2] = bl + xy_;
//        corners[3] = bl + _y_;
//        corners[4] = bl + __z;
//        corners[5] = bl + x_z;
//        corners[6] = tr;
//        corners[7] = tr - x__;
//        
//    }
        
bool Box::Contains(Vector<3> p) {
    
    if ( p[0] > bl[0] && p[0] < tr[0] ) {
        if ( p[1] > bl[1] && p[1] < tr[1] ) {
            if ( p[2] > bl[2] && p[2] < tr[2] ) {
                
                return true;
                
            }
        }
    }
    
    return false;
    
}
        
vector< Vector<3> > Box::GetCorners() {

    return corners;

}
    
vector<Box> Box::GetBabies() {
    
    vector<Box> babies;
    
    double X = diff[0];
    double Y = diff[1];
    double Z = diff[2];
    
    Vector<3> xyz = makeVector(X/2, Y/2, Z/2);
    Vector<3> x_z = makeVector(X/2,0,Z/2);
    Vector<3> xy_ = makeVector(X/2,Y/2,0);
    Vector<3> x__ = makeVector(X/2,0,0);
    Vector<3> _y_ = makeVector(0,Y/2,0);
    Vector<3> __z = makeVector(0,0,Z/2);
    Vector<3> _yz = makeVector(0,Y/2,Z/2);
    
    //        babies[0] = Box( bl                        , bl + makeVector(X/2,Y/2,Z/2) );
    //        babies[1] = Box( bl + makeVector(X/2, 0, 0), bl + makeVector(X,Y/2,Z/2) );
    //        babies[2] = Box( bl + makeVector(X/2, Y/2, 0), tr - makeVector(0,0,Z/2) );
    //        babies[3] = Box( bl + makeVector(0, Y/2, 0), bl + makeVector(X/2,Y,Z/2) );
    //        babies[4] = Box( bl + makeVector(0, 0, Z/2), bl + makeVector(X/2,Y/2,Z) );
    //        babies[5] = Box( bl + makeVector(X/2, 0, Z/2), bl + makeVector(X,Y/2,Z) );
    //        babies[6] = Box( bl + makeVector(X/2, Y/2, Z/2), tr );
    //        babies[7] = Box( bl + makeVector(0, Y/2, Z/2), bl + makeVector(X/2,Y,Z) );
    
    babies.push_back( Box( bl, bl + xyz ) );
    babies.push_back( Box( bl + x__, bl + xyz + x__ ) );
    babies.push_back( Box( bl + xy_, tr - __z ) );
    babies.push_back( Box( bl + _y_, bl + xyz + _y_ ) );
    babies.push_back( Box( bl + __z, bl + xyz + __z ) );
    babies.push_back( Box( bl + x_z, bl + xyz + x__ + __z ) );
    babies.push_back( Box( bl + xyz, tr ) );
    babies.push_back( Box( bl + _yz, bl + xyz + _y_ + __z ) );
    
    return babies;
    
}


