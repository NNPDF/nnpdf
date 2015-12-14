// $Id: thpredictions.h 3044 2015-06-18 19:06:34Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "common.h"
#include "experiments.h"
#include "pdfset.h"
#include "fastkernel.h"
#include "fkset.h"

namespace NNPDF
{
    /**
     * \class ThPredictions
     * \brief Computes the predictions and statistics for a PDFSet and DataSet
     */
    class ThPredictions {

      private:
        ThPredictions();                                //!< Disable default constructor
        ThPredictions& operator=(const ThPredictions&); //!< Disable copy-assignment

        real *fObs;          //!< Theory predictions for Observables
        real fTconv;         //!< time required for the convolution in msecs

        const int fNpdf;     //!< Number of PDFs used in the computation
        const int fNData;    //!< Number of data points

        std::string fPDFName; //!< Name of PDF used to do the calculation
        std::string fSetName; //!< Name of dataset in the calculation

        const int  fEtype;   //!< Uncertainty type

        static void GetNZPDF(const PDFSet*, const FKTable*, real* pdf); //!< Form up PDF array
        static void GetNZPDF(const PDFSet*, const PDFSet*, const FKTable*, real* pdf); //!< Form up PDF array (differing beams)

      public:
        ThPredictions(const PDFSet*, const Experiment*); //!< Constructor
        ThPredictions(const PDFSet*, const FKTable*);    //!< Constructor
        ThPredictions(const PDFSet*, const FKSet*);      //!< Constructor

        ThPredictions(const PDFSet*, const PDFSet*, const FKTable*);    //!< Different-beam constructor

        ThPredictions(const ThPredictions&);  //!< Copy-constructor
        ~ThPredictions();  //!< Destructor

        // Operators
        ThPredictions& operator+=(const ThPredictions&);
        ThPredictions& operator-=(const ThPredictions&);
        ThPredictions& operator*=(const ThPredictions&);
        ThPredictions& operator/=(const ThPredictions&);

        const ThPredictions operator+(const ThPredictions &) const;
        const ThPredictions operator-(const ThPredictions &) const;
        const ThPredictions operator/(const ThPredictions &) const;
        const ThPredictions operator*(const ThPredictions &) const;

        // Fast static convolution functions
        static void Convolute(const PDFSet*,const FKTable*,real*);
        static void Convolute(const PDFSet*,const FKSet*,real*);
        static void Convolute(const PDFSet*,const PDFSet*,const FKTable*,real*);
            
        // Get Methods
        real* GetObs() const {return fObs;}; //!< Return Obs array
      
        real GetObsCV     (int const& idat) const;   //!< Get Central Value
        real GetObsError  (int const& idat) const;   //!< Get Error value (fEtype)
        real GetTConv()   const { return fTconv; }   //!< Get convolution time

        int  GetNPdf()    const { return fNpdf;    } //!< Return the number of pdfs
        int  GetNData()   const { return fNData;   } //!< Return the number of datapoints

	    std::string GetSetName() const { return fSetName; } //!< Return the set name
	    std::string GetPDFName() const { return fPDFName; } //!< Return the set name

        // Set Methods
        void SetPDFName( std::string const& newname ) {fPDFName = newname;};
        void SetDataName( std::string const& newname ) {fSetName = newname;};

        // Print thpredictions to file
        void Print(std::ostream&, bool latex = false) const;

        // Verbosity
        static bool Verbose;
    };
}
