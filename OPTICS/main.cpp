/***************************
Main cpp for testing OPTICS
***************************/
#define NOMINMAX

#include <conio.h>
#include <iostream>
#include <random>
#include <opencv2/opencv.hpp>

#include "OPTICS/optics.hpp"

#include <barn_common.hpp>
#include <barn_open_cv_common.hpp>
#include <Persistence1D/src/persistence1d/persistence1d.hpp>

#include "OPTICS_test.hpp"

using namespace cv;
using namespace std;


const string image_file = "nested.png";


void main() {
    
    Mat3b testset = imread( image_file);
    //noise( testset, 0.05f, Vec3b(255,255,255));
    //testset = testset.t();

    while(1) {

        float eps= -1;
        unsigned int min_pts;
        float persistence = -1;
        unsigned int n_clusters;
        bool use_n_clusters;
        float outlier_threshold;
        cout << "epsilon : "; cin >> eps;
        cout << "min_pts : "; cin >> min_pts;
        cout << "oose n_clusters instead of persistence? : "; cin >> use_n_clusters;
        if( use_n_clusters) {
            cout << "n_clusters : "; cin >> n_clusters;
        } else {
            cout << "persistence : "; cin >> persistence;
        }
        cout << "outlier threshold : "; cin >> outlier_threshold;

        cout << endl;
        

        test_optics( testset, 
                  false,
                  eps, 
                  min_pts, 
                  persistence, 
                  n_clusters, 
                  use_n_clusters, 
                  outlier_threshold);
        cout << "===============================================================================\n";

    } // END while(1)
}