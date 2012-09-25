//
//  Box.h
//  pbPTAM
//
//  Created by Peter Boyer on 4/25/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef pbPTAM_Box_h
#define pbPTAM_Box_h
#include <TooN/TooN.h>
#include <vector>

    
// OCTREE STUFF


namespace PTAMM {
    using namespace TooN;
    using namespace std;
        
    class Box {
        
    public:

        Box(Vector<3> blin, Vector<3> trin);
        bool Contains(Vector<3> p);
        vector< Vector<3> > GetCorners();
        vector<Box> GetBabies();
        
        Vector<3> bl, tr, diff;
        vector< Vector<3> > corners;
        
    };
    
}

#endif
