/***************************************
Functions and classes for testing OPTICS
***************************************/
#define NOMINMAX

#include <conio.h>
#include <iomanip>
#include <iostream>
#include <random>
#include <opencv2/opencv.hpp>

#include "OPTICS/optics.hpp"

#include <barn_common.hpp>
#include <barn_open_cv_common.hpp>
#include <Persistence1D/src/persistence1d/persistence1d.hpp>

using namespace cv;
using namespace std;

const unsigned int max_hist_height = 8000;

const string hist_file_name    = "hist.txt";

const string winname_hist      = "hist";
const string winname_testset   = "testset";
const string winname_resultset = "resultset";



const Vec3b color_background            ( 0  ,   0,   0 );
const Vec3b color_normal_point          ( 224, 224, 224 );
const Vec3b color_hist_bar              ( 224, 224, 224 );
const Vec3b color_hist_cluster_border   ( 255,   0, 255 );
const Vec3b color_hist_unreachable      ( 0,   255,   0 );
const Vec3b color_per_se_reachable      ( 0,     0,  96 );
const Vec3b color_reached               ( 0,     0, 255 );
const Vec3b color_unreachable           ( 255,   0,   0 );
const Vec3b color_marked                ( 0,   255, 255 );
const Vec3b color_marked_reachable      ( 0,   128, 128 );
const Vec3b color_marked_unreachable    ( 255, 128, 128 );
const Vec3b color_invalid               ( 255, 128, 255 ); // u know something is wrong if bright pink is the color you see


// callback class handling all the callback stuff for the hist_mouse_callback function
struct Callback {
    const OPTICS::DataVector& result;
    Mat3b* orig_hist; // backup hist
    Mat3b show_hist;
    const Mat3b* testset;
    float max_reach_dist;
    unsigned int mark_begin;
    unsigned int mark_end;
    int hist_thresh_row;

    /*
    */
    Callback( const float max_r_dist, Mat3b* hist, const Mat3b* test_set, const OPTICS::DataVector& rslt) 
        : max_reach_dist( max_r_dist), testset( test_set), result(rslt), hist_thresh_row(0) {
            
            setHist( hist);
            
    }


    /*
    */
    void setHist( Mat3b* hist) {
        orig_hist = hist;
        hist->copyTo( show_hist);
    }

    /*
    */
    void set_reachability_line( int y) {

        hist_thresh_row = y;

        // draw reachability line in histogram
        orig_hist->copyTo(show_hist);
        Vec3b* show_hist_r_ptr = show_hist.ptr<Vec3b>(y);
        for( int c=0; c < show_hist.cols; ++c) {
            show_hist_r_ptr[c] = Vec3b(255, 0, 0);
        }

        // draw resultset
        Mat3b resultset( testset->rows, testset->cols, Vec3b(0,0,0));
        for( unsigned int c=0; c<result.size(); ++c) {
            
            const OPTICS::DataPoint* dp = result[c];
            const float rdist = dp->reachability_distance();
    
            Vec3b color = color_per_se_reachable;

            if( 1 - (float)y/show_hist.rows >= rdist/max_reach_dist)
                // reached
                color = color_reached;
            if( rdist == OPTICS::UNDEFINED)
                // unreachable
                color = color_unreachable;
            
            resultset( (int)dp->data()[0], (int)dp->data()[1]) = color;
        }

        namedWindow( winname_resultset,  WINDOW_NORMAL);
        imshow( winname_hist,      show_hist);
        imshow( winname_resultset, resultset);
    }
    
    /*
    */
    void mark() {

        // draw histogram
        orig_hist->copyTo(show_hist);
        const unsigned int start = (std::min)(mark_begin, mark_end);
        const unsigned int end = (std::max)(mark_begin, mark_end);
        // draw histogram marking bar
        for( int r=0; r < show_hist.rows; ++r) {
        for( unsigned int c=start; c < end; ++c) {
            Vec3b* p = show_hist.ptr<Vec3b>(r);
            p[c] += Vec3b(0, 0, 192);
        }
        }
        // draw histogram threshold line
        Vec3b* p = show_hist.ptr<Vec3b>( this->hist_thresh_row);
        for( int c=0; c < show_hist.cols; ++c) {
            p[c] = Vec3b(255, 0, 0);
        }


        // draw resultset
        Mat3b resultset( testset->rows, testset->cols, Vec3b(0,0,0));
        unsigned int resultsize = result.size();
        for( unsigned int i=0; i<resultsize; ++i) {
        
            const OPTICS::DataPoint* dp = result[i];
            const float rd = dp->reachability_distance();

            Vec3b color = color_per_se_reachable;
            if( rd == OPTICS::UNDEFINED)
                // unreachable
                color = color_unreachable;

            if( start <= i && end >= i) {
                // marked
                color = color_marked;

                if( rd > max_reach_dist - (float)this->hist_thresh_row/show_hist.rows * max_reach_dist)
                    // marked,not reachable with given threshold
                    color = color_marked_reachable;

                if( rd == OPTICS::UNDEFINED)
                    // unreachable
                    color = color_marked_unreachable;
            }
            
            resultset( (int)dp->data()[0], (int)dp->data()[1]) = color;
        }

        namedWindow( winname_resultset,  WINDOW_NORMAL);
        imshow( winname_hist,      show_hist);
        imshow( winname_resultset, resultset);
    }
};

// histogram window callback function
void hist_mouse_callback( int evt, int x, int y, int flags, void* userdata) {
    
    if( x < 0 || y < 0)
        return;

    Callback* cb = (Callback*)userdata;

    if( evt == cv::EVENT_RBUTTONDOWN) {
        cb->set_reachability_line(y);
    

    } else if( evt == cv::EVENT_LBUTTONDOWN) {
        cb->mark_begin = x;
    
    } else if( evt == cv::EVENT_LBUTTONUP) {
        cb->mark_end = x;
        cb->mark();
    }
}


void test_optics( const Mat3b& testset,
                  float eps, 
                  const unsigned int min_pts, 
                  const float persistence, 
                  const unsigned int n_clusters, 
                  const bool use_n_clusters,
                  const float outlier_threshold);
OPTICS::DataVector scan_testset( const Mat3b& testset);
Mat3b build_histogram( const float rows, const vector<float>& reachabilities);
std::vector<unsigned int> find_k_histogram_peaks( const vector<OPTICS::real>& reachabilities,
                                                  const uint n_clusters);
std::vector<unsigned int> find_histogram_peaks( const vector<OPTICS::real>& reachabilities, 
                                                const OPTICS::real persistence);
vector<Mat3b> create_cluster_images( const vector<OPTICS::DataVector>& clusters, unsigned int rows, unsigned int cols);


/*
*/
void test_optics( const Mat3b& testset,
                  const bool shuffle,
                  float eps, 
                  const unsigned int min_pts, 
                  const float persistence, 
                  const unsigned int n_clusters, 
                  const bool use_n_clusters,
                  const float outlier_threshold) {

    // adjust epsilon
    if( eps < 0) 
        eps = std::numeric_limits<float>::max();

    // print parameters
    cout << ">>> epsilon    : " << eps << endl;
    cout << ">>> min_pts    : " << min_pts << endl;
    if( !use_n_clusters)
        cout << ">>> persistence: " << persistence << endl;
    else
        cout << ">>> n_clusters: " << n_clusters << endl;
    cout << ">>> outlier threshold : " << outlier_threshold << endl;

    // scan test set
    OPTICS::DataVector db =  scan_testset(testset);

    // shuffle data
    if( shuffle) {
        cout << "Shuffling...\n";
        std::random_shuffle( db.begin(), db.end() );
    }

    cout << fixed;

    // run optics
    unsigned int n_processed = 0;
    cout << "\nRunning OPTICS with " << db.size() << " samples...\n";
    OPTICS::DataVector result = OPTICS::optics( db, 
                                                eps, 
                                                min_pts, 
                                                [&n_processed, &db](const OPTICS::DataPoint* p){ 
                                                     
                                                    if(n_processed % 100 == 0)
                                                        cout << setprecision(2) << 100.0f * n_processed / db.size() << "% processed"<< "\n";
                                                    ++n_processed;
                                                });

    cout << "done. Found " << result.size() << " results.\n";

    // extract reachability distances
    vector<float> reachabilities;
    std::for_each(result.begin(), result.end(), [&reachabilities]( const OPTICS::DataPoint* p) { reachabilities.push_back( p->reachability_distance()); });

    // write reachabilities to text file
    to_file( hist_file_name, reachabilities);
    

    // count # unreachables
    const unsigned int n_unreachables = std::count( reachabilities.begin(), reachabilities.end(), OPTICS::UNDEFINED);
    cout << "# unreachables: " << n_unreachables << std::endl;

    // find max maximum reachability distance     (filter out OPTICS::UNDEFINED)
    const float max_r_dist = *std::max_element(reachabilities.begin(), reachabilities.end(), []( float a, float b){ return a == OPTICS::UNDEFINED ? true : a<b; });

    // build histogram
    Mat3b hist = build_histogram( max_r_dist, reachabilities);

    // find histogram maximum peaks
    std::vector<unsigned int> cluster_borders;
    if( use_n_clusters) {
        cluster_borders = find_k_histogram_peaks( reachabilities, n_clusters);
    } else {
        cluster_borders = find_histogram_peaks( reachabilities, persistence);
    }
    std::sort( cluster_borders.begin(), cluster_borders.end());

    // draw cluster borders into histogram
    for( auto it=cluster_borders.begin(); it!=cluster_borders.end(); ++it)
        for( int r=0; r<hist.rows; ++r)
            hist(r,*it) = color_hist_cluster_border;

    // create separate image for each cluster
    vector<OPTICS::DataVector> clusters = OPTICS::extract_clusters( result, cluster_borders, outlier_threshold);
    vector<Mat3b> cluster_images = create_cluster_images( clusters, testset.rows, testset.cols);


    // show images
    namedWindow( winname_testset, WINDOW_NORMAL);
    namedWindow( winname_hist,    WINDOW_NORMAL);
    
    
    // setup the callback function
    Callback callback( max_r_dist, &hist, &testset, result);

    setMouseCallback(winname_hist,  hist_mouse_callback, &callback);
    
    imshow( winname_testset, testset);
    imshow( winname_hist, hist);
    

    // show cluster images
    for( unsigned int i=0; i<cluster_images.size(); ++i) {
        namedWindow( itos(i), WINDOW_NORMAL);
        imshow( itos(i), cluster_images[i]);
        stringstream ss;
        ss << "cluster_ " << i << ".png";
        imwrite( ss.str(), cluster_images[i]);
    }

    

    // write histogram to image file
    imwrite("hist.png", hist);
    
    waitKey();
    destroyAllWindows();

    // cleanup
    for( auto it = db.begin(); it!=db.end(); ++it) {
        delete *it;
    }
}

/*
*/
OPTICS::DataVector scan_testset( const Mat3b& testset) {

    cout << "Scanning " << testset.rows << " x " << testset.cols << " test set... ";
    OPTICS::DataVector db;
    for( int r=0; r<testset.rows; ++r)
    for( int c=0; c<testset.cols; ++c) {
        if( r%50 == 0 && c==0)
            cout << r << "   ";

        if( testset( r, c)[0] > 128) {
            OPTICS::DataPoint* p = new OPTICS::DataPoint();
            p->data().push_back( (float)r);
            p->data().push_back( (float)c);
            db.push_back( p);
        }
    }
    cout << endl;
    return db;
}

/*
*/
Mat3b build_histogram( const float max_r_dist, const vector<float>& reachabilities) {
    
    const float frows = std::min( static_cast<float>(max_hist_height), max_r_dist);
    const unsigned int rows = static_cast<unsigned int>(frows);

    Mat3b ret( rows, reachabilities.size(), color_background);
    for( int c=0; c<ret.cols; ++c) {
        for( int r=0; r<ret.rows; ++r) {
            
            if( (float)reachabilities[c]/rows > (float)(rows-r)/rows)
                ret(r,c) = color_hist_bar;
            if( reachabilities[c] == OPTICS::UNDEFINED)
                ret(r,c) = color_hist_unreachable;
        }
    }
    return ret;
}

/** Given the OPTICS ordered output, finds the k most persistent maxima peaks 
 * of the reachability distances, which are presumably cluster-borders.
 * @param reachabilities The OPTICS ordered reachability distances of the DataPoints 
 *        that where the input of the OPTICS function.
 * @param n_clusters the number of clusters that shall is we want to extract.
 *        The n_clusters-1 most persistent histogram peaks are assumed to be their borders.
 *        If less than n_clusters-1 histogram peaks exist, all peak incides will be returned.
 * @return All the n_clusters-1 biggest histogram peak indices .
 *         The indices are ordered in descending order to the persistence of the peaks at these positions.
 * @see optics()
 */
std::vector<unsigned int> find_k_histogram_peaks( const vector<OPTICS::real>& reachabilities,
                                                  const uint n_clusters) {
    std::vector<unsigned int> ret;

    p1d::Persistence1D p;
    p.RunPersistence(reachabilities);
    vector<p1d::TPairedExtrema> extrema;
    p.GetPairedExtrema(extrema, 0);

    unsigned int n=1;
    for( int i=extrema.size()-1; i>=0 && n<n_clusters; --i, ++n)
        ret.push_back( extrema[i].MaxIndex);

    return ret;
}


/** Given the OPTICS ordered output, finds all maxima peaks with a persistence 
 * greater than a given threshold. These are presumably cluster-borders.
 * @param reachabilities The OPTICS ordered reachability distances of the DataPoints 
 *        that where the input of the OPTICS function.
 * @param persistence The persistence of the histogram peaks to retain.
 * @return All histogram peak indices that are above the given persistence threshold.
 *         The indices are ordered in ascending order to the persistence of the peaks at these positions.
 * @see optics()
 */
std::vector<unsigned int> find_histogram_peaks( const vector<OPTICS::real>& reachabilities, 
                                                const OPTICS::real persistence) {
    std::vector<unsigned int> ret;

    p1d::Persistence1D p;
    p.RunPersistence(reachabilities);
    vector<p1d::TPairedExtrema> extrema;
    p.GetPairedExtrema(extrema, persistence);

    std::for_each(extrema.begin(), extrema.end(), [&ret]( const p1d::TPairedExtrema& e) { ret.push_back( e.MaxIndex); });

    return ret;
}

/*
*/
vector<Mat3b> create_cluster_images( const vector<OPTICS::DataVector>& clusters, unsigned int rows, unsigned int cols) {
    vector<Mat3b> ret;
    
    for( auto it=clusters.begin(); it!=clusters.end(); ++it) {
        const OPTICS::DataVector& cluster = *it;
        Mat3b mat( rows, cols, color_background);

        for( unsigned int i=0; i<it->size(); ++i) {
            const OPTICS::DataPoint& p = *(cluster[i]);
            mat( (int)p[0], (int)p[1]) = color_normal_point;
        }
        ret.push_back( mat);
    }
    return ret;
}