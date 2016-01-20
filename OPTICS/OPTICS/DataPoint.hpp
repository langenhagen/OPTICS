/******************************************************************************
/* @file Contains the DataPoint class that represents multi-dimensional points
/*       in the OPTICS framework.
/*
/*
/* @author langenhagen
/* @version 150520
/******************************************************************************/
#pragma once

///////////////////////////////////////////////////////////////////////////////
// INCLUDES project headers

#include "common.hpp"

///////////////////////////////////////////////////////////////////////////////
//INCLUDES C/C++ standard library (and other external libraries)

#include <assert.h>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
// NAMESPACE, CONSTANTS and TYPE DECLARATIONS/IMPLEMENTATIONS

namespace OPTICS {

    /// Implements multi-dimensional numeric points.
    class DataPoint {

    private: // vars

        std::vector<real> _data;        ///< The data elements.
        real _reachability_distance;    ///< The reachability distance.
        bool _is_processed;             ///< A flag indicating if the object is already processed.
    
    public: // ctor & dtor

        /** Main constructor.
         * Sets the reachability distance to OPTICS::UNDEFINED and sets the processed-flag to false.
         */
        DataPoint() : _data( std::vector<real>()), _reachability_distance( UNDEFINED), _is_processed( false) 
        {}

        //
        // Copy construction done by the compiler generated copy constructor.
        //

        /// Destructor.
        virtual ~DataPoint() 
        {}

    public: // methods

        /** Sets the reachability distance.
         * @param d The new reachability distance. The value must not be negative.
         */
        inline void reachability_distance( real d) {
            assert( d>=0 && "Reachability distance must not be negative.");
            _reachability_distance = d;
        }

        /** Retrieves the current reachability distance.
         * @return The reachability distance. Can be OPTICS::UNDEFINED.
         */
        inline real reachability_distance() const { return _reachability_distance; }
    
        /** Sets the processed flag.
         * @param b The new processed flag.
         */
        inline void processed( bool b) { _is_processed = b; }
    
        /** Retrieves the processed flag.
         * @return Returns either TRUE or FALSE.
         */
        inline bool is_processed() const { return _is_processed; }
    
        /** Retrieves a reference to the data vector.
         * @return A reference to the data vector that stores the data elements.
         */
        inline std::vector<real>& data() { return _data; }

        /** Retrieves a const reference to a data vector.
         * Constant method.
         * @return A const reference to the data vector that stores the data elements.
         */
        inline const std::vector<real>& data() const { return _data; }

    public: // operators

        /** An index operator that retrieves the idx-th element from the data vector.
         * @param idx The index pointing to the position of the data element to retrieve. Must be within the range of the data-vectors dimensionality.
         * @return Returns the element at the idx-th position of the DataPoint.
         * @see data()
         */
        inline real operator[]( const std::size_t idx) const { 
            assert( _data.size()>idx && "Index must be within OPTICS::DataPoint::_data's range.");
            return _data[idx];
        }
    };






    /// Implements multi-dimensional numeric points that can carry a label.
    template<typename T=int>
    class LabelledDataPoint : public DataPoint {

    private: // vars

        T _label; ///< The object's individual label. Can carry anything you want, e.g. a pointer.

    public: // ctor & dtor

        /** Main constructor.
         * Sets the reachability distance to OPTICS::UNDEFINED and sets the processed-flag to false.
         * @param label The label that will be stored within the object.
         */
        LabelledDataPoint( T label) : DataPoint(), _label(label)
        {}

        //
        // Copy construction done by the compiler generated copy constructor.
        //

        /// Destructor.
        ~LabelledDataPoint() 
        {}

    public: // methods

        /** Sets the label.
         * @param l The new label.
         */
        inline void label( const T& l) { _label = l; }

        /** Retrieves the current label.
         * @return The label.
         */
        inline const T& label() const { return _label; }
    };

} // END namespace OPTICS