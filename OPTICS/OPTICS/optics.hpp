/******************************************************************************
/* @file Contains the OPTICS algorithm implementation based on the paper
/*       "OPTICS: Ordering Points To Identify the Clustering Structure"
/*       by Ankerst, Breunig, Kriegel & Sander.
/*       (http://fogo.dbs.ifi.lmu.de/Publikationen/Papers/OPTICS.pdf)
/*
/* 
/* The design & implementation is based on
/*    - readability
/*    - ease of use
/*    - small weight
/*    - zero dependencies (except for the STL)
/*
/*
/* @author langenhagen
/* @version 150605
/******************************************************************************/
#pragma once

///////////////////////////////////////////////////////////////////////////////
// INCLUDES project headers

#include "DataPoint.hpp"

///////////////////////////////////////////////////////////////////////////////
//INCLUDES C/C++ standard library (and other external libraries)

#include <algorithm> // nth_element
#include <functional>

///////////////////////////////////////////////////////////////////////////////
// NAMESPACE, CONSTANTS and TYPE DECLARATIONS/IMPLEMENTATIONS


/// Namespace of the OPTICS module.
namespace OPTICS {



    // FUNCTION DECLARATIONS ######################################################################

    // non-callback version
    DataVector optics( DataVector& db, const real eps, const unsigned int min_pts);
    void expand_cluster_order( DataVector& db, DataPoint* p, const real eps, const unsigned int min_pts, DataVector& o_ordered_vector);
    
    // callback version
    DataVector optics( DataVector& db, 
                       const real eps, 
                       const unsigned int min_pts,
                       std::function<void(const DataPoint* p)> point_processed_callback);
    void expand_cluster_order( DataVector& db, 
                               DataPoint* p, 
                               const real eps, 
                               const unsigned int min_pts, 
                               DataVector& o_ordered_vector, 
                               std::function<void(const DataPoint* p)> point_processed_callback);
    
    // utility functions
    std::vector<DataVector> extract_clusters( const DataVector& result, const std::vector<unsigned int>& cluster_borders, real outlier_threshold);

    // helpers
    void update_seeds( const DataVector& N_eps, const DataPoint* center_object, const real c_dist, DataSet& o_seeds);
    DataVector get_neighbors( const DataPoint* p, const real eps, DataVector& db);
    real squared_core_distance( const DataPoint* p, const unsigned int min_pts, DataVector& N_eps);
    real squared_distance( const DataPoint* a, const DataPoint* b);
    


    // NON-CALLBACK VERSION #######################################################################


    /** Performs the classic OPTICS algorithm.
     * @param db All data points that are to be considered by the algorithm. Changes their values.
     * @param eps The epsilon representing the radius of the epsilon-neighborhood.
     * @param min_pts The minimum number of points to be found within an epsilon-neigborhood.
     * @return Return the OPTICS ordered list of Data points with reachability-distances set.
     */
    DataVector optics( DataVector& db, const real eps, const unsigned int min_pts) {
        assert( eps >= 0 && "eps must not be negative");
        assert( min_pts > 0 && "min_pts must be greater than 0");
        DataVector ret;

        for( auto p_it = db.begin(); p_it != db.end(); ++p_it) {
            DataPoint* p = *p_it;
            
            if( p->is_processed())
                continue;
            
            expand_cluster_order( db, p, eps, min_pts, ret);
        }
        return ret;
    }


    /** Expands the cluster order while adding new neighbor points to the order.
     * @param db All data points that are to be considered by the algorithm. Changes their values.
     * @param p The point to be examined.
     * @param eps The epsilon representing the radius of the epsilon-neighborhood.
     * @param min_pts The minimum number of points to be found within an epsilon-neigborhood.
     * @param o_ordered_vector The ordered vector of data points. Elements will be added to this vector.
     */
    void expand_cluster_order( DataVector& db, DataPoint* p, const real eps, const unsigned int min_pts, DataVector& o_ordered_vector) {
        assert( eps >= 0 && "eps must not be negative");
        assert( min_pts > 0 && "min_pts must be greater than 0");
        
        DataVector N_eps = get_neighbors( p, eps, db);
        p->reachability_distance( OPTICS::UNDEFINED);
        const real core_dist_p = squared_core_distance( p, min_pts, N_eps);
        p->processed( true);
        o_ordered_vector.push_back( p);
    
        if( core_dist_p == OPTICS::UNDEFINED)
            return;

        DataSet seeds;
        update_seeds( N_eps, p, core_dist_p, seeds);
        
        while( !seeds.empty()) {
            DataPoint* q = *seeds.begin();
            seeds.erase( seeds.begin()); // remove first element from seeds

            DataVector N_q = get_neighbors( q, eps, db);
            const real core_dist_q = squared_core_distance( q, min_pts, N_q);
            q->processed( true);
            o_ordered_vector.push_back( q);
            if( core_dist_q != OPTICS::UNDEFINED) {
                // *** q is a core-object ***
                update_seeds( N_q, q, core_dist_q, seeds);
            }
        }
    }



    // CALLBACK VERSION ###########################################################################


    /** Performs the classic OPTICS algorithm.
     * Because OPTICS can take a while on big data sets or when working with high dimensions,
     * a callback function informs you when a new point is inserted into the OPTICS ordering.
     * @param db All data points that are to be considered by the algorithm. Changes their values.
     * @param eps The epsilon representing the radius of the epsilon-neighborhood.
     * @param min_pts The minimum number of points to be found within an epsilon-neigborhood.
     * @param point_processed_callback Callback function that is called when one point is 
     *        added to the ordered output list. It takes the pointer to the data point as an argument.
     * @return Return the OPTICS ordered list of Data points with reachability-distances set.
     */
    DataVector optics( DataVector& db, 
                       const real eps, 
                       const unsigned int min_pts, 
                       std::function<void(const DataPoint* p)> point_processed_callback) {
        assert( eps >= 0 && "eps must not be negative");
        assert( min_pts > 0 && "min_pts must be greater than 0");
        DataVector ret;

        for( auto p_it = db.begin(); p_it != db.end(); ++p_it) {
            DataPoint* p = *p_it;
            
            if( p->is_processed())
                continue;
            
            expand_cluster_order( db, p, eps, min_pts, ret, point_processed_callback);
        }
        return ret;
    }


    /** Expands the cluster order while adding new neighbor points to the order.
     * Because OPTICS can take a while on big data sets or when working with high dimensions,
     * a callback function informs you when a new point is inserted into the OPTICS ordering.
     * @param db All data points that are to be considered by the algorithm. Changes their values.
     * @param p The point to be examined.
     * @param eps The epsilon representing the radius of the epsilon-neighborhood.
     * @param min_pts The minimum number of points to be found within an epsilon-neigborhood.
     * @param o_ordered_vector The ordered vector of data points. Elements will be added to this vector.
     * @param point_processed_callback Callback function that is called when one point is 
     *        added to the ordered output list. It takes the pointer to the data point as an argument.
     */
    void expand_cluster_order( DataVector& db,
                               DataPoint* p, 
                               const real eps,
                               const unsigned int min_pts,
                               DataVector& o_ordered_vector,
                               std::function<void(const DataPoint* p)> point_processed_callback) {
        assert( eps >= 0 && "eps must not be negative");
        assert( min_pts > 0 && "min_pts must be greater than 0");
        
        DataVector N_eps = get_neighbors( p, eps, db);
        p->reachability_distance( OPTICS::UNDEFINED);
        const real core_dist_p = squared_core_distance( p, min_pts, N_eps);
        p->processed( true);
        o_ordered_vector.push_back( p);
        point_processed_callback( p);
    
        if( core_dist_p == OPTICS::UNDEFINED)
            return;

        DataSet seeds;
        update_seeds( N_eps, p, core_dist_p, seeds);
        
        while( !seeds.empty()) {
            DataPoint* q = *seeds.begin();
            seeds.erase( seeds.begin()); // remove first element from seeds

            DataVector N_q = get_neighbors( q, eps, db);
            const real core_dist_q = squared_core_distance( q, min_pts, N_q);
            q->processed( true);
            o_ordered_vector.push_back( q);
            point_processed_callback( p);
            if( core_dist_q != OPTICS::UNDEFINED) {
                // *** q is a core-object ***
                update_seeds( N_q, q, core_dist_q, seeds);
            }
        }
    }


    
    // HELPERS ####################################################################################


    /** Updates the seeds priority queue with new neighbors or neighbors that now have a better 
     * reachability distance than before.
     * @param N_eps All points in the the epsilon-neighborhood of the center_object, including p itself.
     * @param center_object The point on which to start the update process.
     * @param c_dist The core distance of the given center_object.
     * @param o_seeds The seeds priority queue (aka set with special comparator function) that will be modified.
     */
    void update_seeds( const DataVector& N_eps, const DataPoint* center_object, const real c_dist, DataSet& o_seeds) {
        assert( c_dist != OPTICS::UNDEFINED && "the core distance must be set <> UNDEFINED when entering update_seeds");
        
        for( DataVector::const_iterator it=N_eps.begin(); it!=N_eps.end(); ++it) {
            DataPoint* o = *it;

            if( o->is_processed())
                continue;

            const real new_r_dist = std::max( c_dist, squared_distance( center_object, o));
            // *** new_r_dist != UNDEFINED ***
        
            if( o->reachability_distance() == OPTICS::UNDEFINED) {
                // *** o not in seeds ***
                o->reachability_distance( new_r_dist);
                o_seeds.insert( o);

            } else if( new_r_dist < o->reachability_distance()) {
                // *** o already in seeds & can be improved ***
                o_seeds.erase( o);
                o->reachability_distance( new_r_dist);
                o_seeds.insert( o);
            }
        }
    }


    /** Retrieves all points in the epsilon-surrounding of the given data point, including the point itself.
     * @param p The datapoint which represents the center of the epsilon surrounding.
     * @param eps The epsilon value that represents the radius for the neigborhood search.
     * @param db The database consisting of all datapoints that are checked for neighborhood.
     * @param A vector of pointers to datapoints that lie within the epsilon-neighborhood 
     *        of the given point p, including p itself.
     */
    DataVector get_neighbors( const DataPoint* p, const real eps, DataVector& db) {
        assert( eps >= 0 && "eps must not be negative");
        DataVector ret;

        const real eps_sq = eps*eps;

        for( auto q_it=db.begin(); q_it!=db.end(); ++q_it) {
            DataPoint* q = *q_it;
            if( squared_distance( p, q) <= eps_sq) {
                ret.push_back( q);
            }
        }
        return ret;
    }


    /** Finds the squared core distance of one given point.
     * @param p The point to be examined.
     * @param min_pts The minimum number of points to be found within an epsilon-neigborhood.
     * @param N_eps All points in the the epsilon-neighborhood of p, including p itself.
     * @return The squared core distance of p.
     */
    real squared_core_distance( const DataPoint* p, const unsigned int min_pts, DataVector& N_eps) {    
        assert( min_pts > 0 && "min_pts must be greater than 0");
        real ret( OPTICS::UNDEFINED);
    
        if( N_eps.size() > min_pts) {
            std::nth_element( N_eps.begin(), 
                              N_eps.begin()+min_pts, 
                              N_eps.end(), 
                              [p]( const DataPoint* a, const DataPoint* b){ return squared_distance( p, a) < squared_distance( p, b); } );

            ret = squared_distance( p, N_eps[min_pts]);
        }
        return ret;
    }


    /** Retrieves the squared euclidean distance of two DataPoints.
     * @param a The first DataPoint.
     * @param b The second DataPoint. Both data points must have the same dimensionality.
     */
    real squared_distance( const DataPoint* a, const DataPoint* b) {
        const std::vector<real>& a_data = a->data();
        const std::vector<real>& b_data = b->data();
        const unsigned int vec_size = static_cast<unsigned int>(a_data.size());
        assert( vec_size == b_data.size() && "Data-vectors of both DataPoints must have same dimensionality");
        real ret(0);
    
        for( unsigned int i=0; i<vec_size; ++i) {
            ret += std::pow( a_data[i]-b_data[i], 2);
        }
        //return std::sqrt( ret);
        return ret;
    }


    /// A Comp_DataPoint_f comparison functor ()-operator implementation.
    bool Comp_DataPoint_Ptr_f::operator() (const DataPoint* lhs, const DataPoint* rhs) const {
        assert( lhs != nullptr && "nullptr objects are not allowed");
        assert( rhs != nullptr && "nullptr objects are not allowed");
        assert( lhs->data().size() == rhs->data().size() && "Comparing DataPoints requires them to have same dimensionality");

        //return lhs->reachability_distance() < rhs->reachability_distance();    
        if( lhs->reachability_distance() < rhs->reachability_distance())
            return true;
        else if( lhs->reachability_distance() == rhs->reachability_distance() && lhs < rhs)
            return true;
        else /*lhs->reachability_distance() == rhs->reachability_distance() && lhs >= rhs || lhs->reachability_distance() > rhs->reachability_distance()*/
            return false;
    }


    
    // UTILITY FUNCTIONS ##########################################################################


    /** Partitions the specified OPTICS ordered data points along the given cluster borders.
     * Points that lie above a specified threshold are put into a separate outlier cluster.
     * @param result The OPTICS ordered result vector of the optics function.
     * @param cluster_borders A vector of indices specifiying the cluster borders.
     *        IMPORTANT: The vector must be sorted in ascending order.
     * @param outlier_threshold All values above that outlier_threshold are considered outliers
     *        and will be put in a special outlier cluster. Is the threshold value set 
     *        to 0 or negative no point will be considered as an outlier.
     * @return A vector of different disjoint data point containers, each making up one cluster. 
     *         The first container stores the data points that are considered outliers.
     * @see optics()
     */
    std::vector<DataVector> extract_clusters( const DataVector& result, const std::vector<unsigned int>& cluster_borders, real outlier_threshold) {
        std::vector<DataVector> ret;
        ret.push_back( DataVector()); // outlier container
    
        if( outlier_threshold <= 0)
            outlier_threshold = std::numeric_limits<real>::max();
    
    
        for( unsigned int i=0; i<=cluster_borders.size(); ++i) {
            
            const unsigned int lower_idx = i == 0                        ? 0                                        : cluster_borders[i-1];
            const unsigned int upper_idx = i == cluster_borders.size()   ? static_cast<unsigned int>(result.size()) : cluster_borders[i];
        
            DataVector cluster_i;

            for( unsigned int j=lower_idx; j<upper_idx; ++j) {
               DataPoint* p = result[j];
           
               if( p->reachability_distance() > outlier_threshold) {
                   ret[0].push_back( p);
               } else {
                   cluster_i.push_back( p);
               }
            }
            ret.push_back( cluster_i);
        }
        return ret;
    }

} // END namespace OPTICS
