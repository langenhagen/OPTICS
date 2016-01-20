/******************************************************************************
/* @file Contains common elements, constants and typedefs of the OPTICS module.
/*
/*
/* @author langenhagen
/* @version 150515
/******************************************************************************/
#pragma once

///////////////////////////////////////////////////////////////////////////////
//INCLUDES C/C++ standard library (and other external libraries)

#include <limits>
#include <set>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
// NAMESPACE, CONSTANTS and TYPE DECLARATIONS/IMPLEMENTATIONS


/// Namespace of the OPTICS module.
namespace OPTICS {

    /// typedef for abstracting single/double precision. Change at will.
    typedef float real;
    
    /// "Undefined" value for distance measures (which are always >= 0 by nature).
    const real UNDEFINED = std::numeric_limits<real>::max();

    /// The DataPoint class.
    class DataPoint;
    
    /** A comparison functor for comparing DataPoints according to their reachability distance.
     * Reachability distance values must not be UNDEFINED for both left hand side and right hand side operands.
     */
    struct Comp_DataPoint_Ptr_f { 
        bool operator() (const DataPoint* lhs, const DataPoint* rhs) const; 
    };
    
    /// A set of data points equipped with a Comp_DataPoint_Ptr_f comparison functor.
    typedef std::set<DataPoint*, Comp_DataPoint_Ptr_f> DataSet;

    /// A vector of Pointers to DataPoints.
    typedef std::vector<DataPoint*> DataVector;

} // END namespace OPTICS