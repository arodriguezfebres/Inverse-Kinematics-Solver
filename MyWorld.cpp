#include "MyWorld.h"
#include <iostream>

using namespace Eigen;
using namespace dart::dynamics;

MyWorld::MyWorld() {
  
	//std::cout << "CHECKPOINT 0" << std::endl;

	// Load a skeleton from file
  mSkel = dart::utils::SkelParser::readSkeleton(DART_DATA_PATH"skel/human.skel");

  // Create markers
  createMarkers();
  
  // Initialize Jacobian assuming that there is only one constraint
  mJ = MatrixXd::Zero(300, mSkel->getNumDofs());
  //mTarget = Vector3d::Zero(300); This is getting frustrating.
  for (int index = 0; index < 100; index++) {
	  mTarget[index] = Vector3d(0.0, 0.0, 0.0);
  }
  mC = VectorXd::Zero(300); //Fix arbitrary size decision.

}

MyWorld::~MyWorld() {
}

void MyWorld::solve() {
  if (mConstrainedMarker.empty())
    return; 
  int numIter = 50;
  double alpha = 0.01;
  int nDof = mSkel->getNumDofs();
  VectorXd gradients(nDof);
  VectorXd newPose(nDof);
  //std::cout << "CHECKPOINT 1" << std::endl;
  for (int i = 0; i < numIter; i++) {
    gradients = updateGradients();
    newPose = mSkel->getPositions() - alpha * gradients;
    mSkel->setPositions(newPose); 
    mSkel->computeForwardKinematics(true, false, false); // DART updates all the transformations based on newPose
  }
  
}

// Current code only works for the left leg with only one constraint
VectorXd MyWorld::updateGradients() {

	Vector4d offset;
	BodyNode *node;
	Matrix4d worldToParent;// = node->getParentBodyNode()->getTransform().matrix();
	Matrix4d parentToJoint;// = joint->getTransformFromParentBodyNode().matrix();
	Matrix4d dR;// = joint->getTransformDerivative(0); // Doesn't need .matrix() because it returns a Matrix4d instead of Isometry3d
	Matrix4d R;// = joint->getTransform(1).matrix();
	Matrix4d jointToChild;// = joint->getTransformFromChildBodyNode().inverse().matrix();
	Vector4d jCol;// = worldToParent * parentToJoint * dR * R * jointToChild * offset;
	int colIndex;// = joint->getIndexInSkeleton(0);
	Matrix4d R2;




	for (int index = 0; index < mConstrainedMarker.size(); index++) {
		
		//std::cout << "CHECKPOINT 2" << std::endl;

		int markerNum = mConstrainedMarker[index];
		//Vector3d posChange = getMarker(markerNum)->getWorldPosition() - Vector3d(mTarget(markerNum, 0), mTarget(markerNum, 1), mTarget(markerNum, 2) );
		Vector3d posChange = getMarker(markerNum)->getWorldPosition() - mTarget[markerNum];


		mC[(index * 3) + 0] = posChange(0);
		mC[(index * 3) + 1] = posChange(1);
		mC[(index * 3) + 2] = posChange(2);

		offset << getMarker(markerNum)->getLocalPosition(), 1;
		node = getMarker(markerNum)->getBodyNode();

		//std::cout << "CHECKPOINT 3" << std::endl;

		while (node != NULL) {
			//Do we have to check if the joint is valid? Abdomen?
			Joint *joint = node->getParentJoint();
			if (joint != NULL) {
				int myDOF = joint->getNumDofs();
				bool isRoot = (node->getParentBodyNode() == NULL);
				switch (myDOF) {
					case 1:
						//logic 1

						//1st and only DOF. This comment is a little unnecessary.
						if (!isRoot) {
							worldToParent = node->getParentBodyNode()->getTransform().matrix();
							parentToJoint = joint->getTransformFromParentBodyNode().matrix();
							jointToChild = joint->getTransformFromChildBodyNode().inverse().matrix();
							dR = joint->getTransformDerivative(0);

							jCol = worldToParent * parentToJoint * dR * jointToChild * offset;
						}
						else {
							jointToChild = joint->getTransformFromChildBodyNode().inverse().matrix();
							dR = joint->getTransformDerivative(0);
							
							jCol = dR * jointToChild * offset;
						}
						
						colIndex = joint->getIndexInSkeleton(0);
						mJ.col(colIndex)[(index * 3) + 0] = jCol.x(); mJ.col(colIndex)[(index * 3) + 1] = jCol.y(); mJ.col(colIndex)[(index * 3) + 2] = jCol.z();

						offset = ((!isRoot) ? parentToJoint * joint->getTransform(0).matrix() * jointToChild * offset : joint->getTransform(0).matrix() * jointToChild * offset); //I could set p2j to identity mat?
						break;
					case 2:
						//logic 2 - alternate dR and R for both Dofs

						//1st DOF
						if (!isRoot) {
							worldToParent = node->getParentBodyNode()->getTransform().matrix();
							parentToJoint = joint->getTransformFromParentBodyNode().matrix();
							jointToChild = joint->getTransformFromChildBodyNode().inverse().matrix();
							dR = joint->getTransformDerivative(0);
							R = joint->getTransform(1).matrix();

							jCol = worldToParent * parentToJoint * dR * R * jointToChild * offset;
						}
						else {
							jointToChild = joint->getTransformFromChildBodyNode().inverse().matrix();
							dR = joint->getTransformDerivative(0);
							R = joint->getTransform(1).matrix();

							jCol = dR * R * jointToChild * offset;
						}

						colIndex = joint->getIndexInSkeleton(0);
						mJ.col(colIndex)[(index * 3) + 0] = jCol.x(); mJ.col(colIndex)[(index * 3) + 1] = jCol.y(); mJ.col(colIndex)[(index * 3) + 2] = jCol.z();

						//2nd DOF
						//std::cout << "CHECKPOINT 4." << std::endl;
						if (!isRoot) {
							//std::cout << "CHECKPOINT 5." << std::endl;
							dR = joint->getTransformDerivative(1);
							R = joint->getTransform(0).matrix();

							jCol = worldToParent * parentToJoint * dR * R * jointToChild * offset;
						}
						else {
							dR = joint->getTransformDerivative(1);
							R = joint->getTransform(0).matrix();

							jCol = dR * R * jointToChild * offset;
						}

						colIndex = joint->getIndexInSkeleton(1);
						mJ.col(colIndex)[(index * 3) + 0] = jCol.x(); mJ.col(colIndex)[(index * 3) + 1] = jCol.y(); mJ.col(colIndex)[(index * 3) + 2] = jCol.z();

						offset = ((!isRoot) ? parentToJoint * joint->getTransform(0).matrix() * joint->getTransform(1).matrix() * jointToChild * offset 
							: joint->getTransform(0).matrix() * joint->getTransform(1).matrix() * jointToChild * offset);
						break;
					case 3:
						//logic 3

						//1st DOF
						if (!isRoot) {
							worldToParent = node->getParentBodyNode()->getTransform().matrix();
							parentToJoint = joint->getTransformFromParentBodyNode().matrix();
							jointToChild = joint->getTransformFromChildBodyNode().inverse().matrix();
							dR = joint->getTransformDerivative(0);
							R = joint->getTransform(1).matrix();
							R2 = joint->getTransform(2).matrix();

							jCol = worldToParent * parentToJoint * dR * R * R2 * jointToChild * offset;
						}
						else {
							jointToChild = joint->getTransformFromChildBodyNode().inverse().matrix();
							dR = joint->getTransformDerivative(0);
							R = joint->getTransform(1).matrix();
							R2 = joint->getTransform(2).matrix();

							jCol = dR * R * R2 * jointToChild * offset;
						}

						colIndex = joint->getIndexInSkeleton(0);
						mJ.col(colIndex)[(index * 3) + 0] = jCol.x(); mJ.col(colIndex)[(index * 3) + 1] = jCol.y(); mJ.col(colIndex)[(index * 3) + 2] = jCol.z();

						//2nd DOF
						if (!isRoot) {
							dR = joint->getTransformDerivative(1);
							R = joint->getTransform(0).matrix();
							R2 = joint->getTransform(2).matrix();

							jCol = worldToParent * parentToJoint * dR * R * R2 * jointToChild * offset;
						}
						else {
							dR = joint->getTransformDerivative(1);
							R = joint->getTransform(0).matrix();
							R2 = joint->getTransform(2).matrix();

							jCol = dR * R * R2 * jointToChild * offset;
						}

						colIndex = joint->getIndexInSkeleton(1);
						mJ.col(colIndex)[(index * 3) + 0] = jCol.x(); mJ.col(colIndex)[(index * 3) + 1] = jCol.y(); mJ.col(colIndex)[(index * 3) + 2] = jCol.z();

						//3rd DOF
						if (!isRoot) {
							dR = joint->getTransformDerivative(2);
							R = joint->getTransform(0).matrix();
							R2 = joint->getTransform(1).matrix();

							jCol = worldToParent * parentToJoint * dR * R * R2 * jointToChild * offset;
						}
						else {
							dR = joint->getTransformDerivative(2);
							R = joint->getTransform(0).matrix();
							R2 = joint->getTransform(1).matrix();

							jCol = dR * R * R2 * jointToChild * offset;
						}

						colIndex = joint->getIndexInSkeleton(2);
						mJ.col(colIndex)[(index * 3) + 0] = jCol.x(); mJ.col(colIndex)[(index * 3) + 1] = jCol.y(); mJ.col(colIndex)[(index * 3) + 2] = jCol.z();

						offset = ((!isRoot) ? parentToJoint * joint->getTransform(0).matrix() * joint->getTransform(1).matrix() * joint->getTransform(2).matrix() * jointToChild * offset 
							: joint->getTransform(0).matrix() * joint->getTransform(1).matrix() * joint->getTransform(2).matrix() * jointToChild * offset);

						break;
					default:
						//Pointer failure.
						break;
				}
			}
			node = node->getParentBodyNode();
		}
	}



  // compute gradients
  VectorXd gradients = 2 * mJ.transpose() * mC;
  return gradients;
}

// Current code only handles one constraint on the left foot.
void MyWorld::createConstraint(int _index) {
  /*if (_index == 0) {
    mTarget = getMarker(_index)->getWorldPosition();
    mConstrainedMarker = _index;
  } else {
    mConstrainedMarker = -1;
  }*/

	Vector3d myPos = getMarker(_index)->getWorldPosition();

	//mTarget(_index, 0) = myPos.x();
	//mTarget(_index, 1) = myPos.y();
	//mTarget(_index, 2) = myPos.z();
	
	mTarget[_index] = myPos;
	
	mConstrainedMarker.push_back(_index);
	
	for (int index = 0; index < mConstrainedMarker.size(); index++) {
		std::cout << mConstrainedMarker[index] << std::endl;
	}
}

void MyWorld::modifyConstraint(int _index, Vector3d _deltaP) {
  //if (mConstrainedMarker == 0)
  //  mTarget += _deltaP;
	mTarget[_index] += _deltaP;
}

void MyWorld::removeConstraint(int _index) {
	mConstrainedMarker.erase(std::remove(mConstrainedMarker.begin(), mConstrainedMarker.end(), _index), mConstrainedMarker.end());
}

Marker* MyWorld::getMarker(int _index) {
  return mMarkers[_index];
}

void MyWorld::createMarkers() {
  Vector3d offset(0.2, 0.0, 0.0);
  BodyNode* bNode = mSkel->getBodyNode("h_heel_right");
  Marker* m = new Marker("right_foot", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.2, 0.0, 0.0);
  bNode = mSkel->getBodyNode("h_heel_left");
  m = new Marker("left_foot", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.065, -0.3, 0.0);
  bNode = mSkel->getBodyNode("h_thigh_right");
  m = new Marker("right_thigh", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.065, -0.3, 0.0);
  bNode = mSkel->getBodyNode("h_thigh_left");
  m = new Marker("left_thigh", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.0, 0.0, 0.13);
  bNode = mSkel->getBodyNode("h_pelvis");
  m = new Marker("pelvis_right", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.0, 0.0, -0.13);
  bNode = mSkel->getBodyNode("h_pelvis");
  m = new Marker("pelvis_left", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.075, 0.1, 0.0);
  bNode = mSkel->getBodyNode("h_abdomen");
  m = new Marker("abdomen", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.0, 0.18, 0.075);
  bNode = mSkel->getBodyNode("h_head");
  m = new Marker("head_right", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.0, 0.18, -0.075);
  bNode = mSkel->getBodyNode("h_head");
  m = new Marker("head_left", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.0, 0.22, 0.0);
  bNode = mSkel->getBodyNode("h_scapula_right");
  m = new Marker("right_scapula", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.0, 0.22, 0.0);
  bNode = mSkel->getBodyNode("h_scapula_left");
  m = new Marker("left_scapula", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.0, -0.2, 0.05);
  bNode = mSkel->getBodyNode("h_bicep_right");
  m = new Marker("right_bicep", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.0, -0.2, -0.05);
  bNode = mSkel->getBodyNode("h_bicep_left");
  m = new Marker("left_bicep", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.0, -0.1, 0.025);
  bNode = mSkel->getBodyNode("h_hand_right");
  m = new Marker("right_hand", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);

  offset = Vector3d(0.0, -0.1, -0.025);
  bNode = mSkel->getBodyNode("h_hand_left");
  m = new Marker("left_hand", offset, bNode);
  mMarkers.push_back(m);
  bNode->addMarker(m);
}
