#ifndef _MYWORLD_
#define _MYWORLD_

#include <vector>
#include <algorithm>
#include "dart/dart.h"


class MyWorld {
 public:
    MyWorld();
    virtual ~MyWorld();
    dart::dynamics::SkeletonPtr getSkel() {
        return mSkel;
    }

    void solve();
    void createConstraint(int _index);
    void modifyConstraint(int _index, Eigen::Vector3d _deltaP);
    void removeConstraint(int _index);
    dart::dynamics::Marker* getMarker(int _index);

 protected:
    Eigen::VectorXd updateGradients();
    void createMarkers();

    dart::dynamics::SkeletonPtr mSkel;
    std::vector<dart::dynamics::Marker*> mMarkers;
    Eigen::VectorXd mC;
    Eigen::MatrixXd mJ;
    Eigen::Vector3d mTarget[100]; // The target location of the constriant
    std::vector <int> mConstrainedMarker; // The indices of the constrained markers
};

#endif
