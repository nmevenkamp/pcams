#ifndef __SEGMENTATION_H
#define __SEGMENTATION_H

#include <quoc.h>
#include <configurators.h>
#include <convolution.h>
#include <firstOrderTVAlgos.h>
#include <ctrlCCatcher.h>
#include <progressBar.h>
#include <gnuplotter.h>
#include <connectedComponents.h>
#include <clustering.h>
#ifdef USE_LIB_HUNGARIAN
#include <Hungarian.h>
#endif


template <typename RealType>
static inline void cropSegment ( qc::ScalarArray<RealType, qc::QC_2D> &ArgDest, const qc::ScalarArray<int, qc::QC_2D> &Segmentation, const int SegmentIdx ) {
  const int nx = ArgDest.getNumX ( ), ny = ArgDest.getNumY ( );
  
  if ( Segmentation.getNumX ( ) != nx || Segmentation.getNumY ( ) != ny ) throw aol::Exception ( "Input and segmentation dimensions do not match!", __FILE__, __LINE__ );
  if ( SegmentIdx < 0 || SegmentIdx > Segmentation.getMaxValue ( ) ) throw aol::Exception ( "Segment index invalid!", __FILE__, __LINE__ );
  
  // Determine crop start & end
  aol::Vec<2,int> cropStart; cropStart[0] = -1; cropStart[1] = -1;
  aol::Vec<2,int> cropEnd; cropEnd[0] = nx; cropEnd[1] = ny;
  qc::ScalarArray<int, qc::QC_1D> lineX ( nx ), lineY ( ny );
  do { ++cropStart[0]; Segmentation.getLine ( qc::QC_X, cropStart[0], lineY ); } while ( cropStart[0] < nx-1 && lineY.numOccurence ( SegmentIdx ) == 0 );
  do { ++cropStart[1]; Segmentation.getLine ( qc::QC_Y, cropStart[1], lineX ); } while ( cropStart[1] < ny-1 && lineX.numOccurence ( SegmentIdx ) == 0 );
  do { --cropEnd[0]; Segmentation.getLine ( qc::QC_X, cropEnd[0], lineY ); } while ( cropEnd[0] > 0 && lineY.numOccurence ( SegmentIdx ) == 0 );
  do { --cropEnd[1]; Segmentation.getLine ( qc::QC_Y, cropEnd[1], lineX ); } while ( cropEnd[1] > 0 && lineX.numOccurence ( SegmentIdx ) == 0 );

  aol::Vec<2,int> cropSize; cropSize[0] = cropEnd[0] - cropStart[0] + 1; cropSize[1] = cropEnd[1] - cropStart[1] + 1;
  ArgDest.crop ( cropStart, cropSize );
}

static inline void transformGrayValuesToLabels ( qc::ScalarArray<int, qc::QC_2D> &Segmentation ) {
  int label = -1;
  while ( Segmentation.getMaxValue ( ) >= 0 ) {
    int val = Segmentation.getMaxValue ( );
    for ( int i=0; i<Segmentation.size ( ) ; ++i ) {
      if ( Segmentation[i] == val )
        Segmentation[i] = label;
    }
    --label;
  }
  Segmentation *= -1;
  Segmentation.addToAll ( -1 );
}

template <typename RealType, typename PictureType>
static inline void plotSegmentationBoundariesOntoImage ( const char* path, const PictureType &Input, const qc::ScalarArray<int, qc::QC_2D> &Segmentation ) {
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > ConfiguratorType;
  
  qc::ScalarArray<int, qc::QC_2D> seg ( Segmentation );
  transformGrayValuesToLabels ( seg );
  aol::Plotter<RealType> plotter;
  aol::PlotDataFileHandler<RealType> plotDataFileHandler;
  plotDataFileHandler.generateBackgroundPNGData ( Input );
  plotter.addPlotCommandsFromHandler ( plotDataFileHandler );
  
  aol::MultiVector<RealType> levelsetFunctions ( seg.getMaxValue ( ) + 1, seg.size ( ) );
  for ( int i=0; i<seg.size ( ) ; ++i ) levelsetFunctions[seg[i]][i] = 1.0;
  typename ConfiguratorType::InitType grid ( qc::GridSize<qc::QC_2D> ( Input.getNumX ( ), Input.getNumY ( ) ) );
  ConfiguratorType confDummy ( grid );
  plotDataFileHandler.generateIsolineData ( levelsetFunctions, grid, confDummy, 0.5, false, true );
  const RealType rgb[3] = {1, 0, 0};
  for ( unsigned int i=1; i<plotDataFileHandler.getDataFileNames().size ( ) ; ++i )
    plotter.addLinePlot ( plotDataFileHandler.getDataFileNames ( )[i].c_str ( ), NULL, 3, rgb );
  
  plotter.set_outfile_base_name ( path );
  plotter.setBackgroundPNGSettings ( Input.getNumX ( ), Input.getNumY ( ) );
  plotter.genPlot ( aol::GNUPLOT_PDF );
}


/**
 *  \brief Hungarian Method (Kuhn-Munkres algorithm) for finding a perfect (edge weight minimizing) matching in a bipartite graph
 *
 *  Input: adjacency matrix of a bipartite graph (edge weights)
 *         use the value NEG_INF to designate non-existing edges
 *  Output: a vector that assigns each node in the left part of the graph (index of the vector entry) at most one node in the right part of the graph (value of the vector entry)
 *          vector entries equal to -1 indicate that the corresponding node in the left part of the graph is not connected to any node in the right part of the graph
 *
 *  This is just a wrapper for the Hungarian library found at: http://www.frc.ri.cmu.edu/~lantao/codes/hungarian.php
 *
 *  Lantao Liu, Dylan Shell. "Assessing Optimal Assignment under Uncertainty: An Interval-based Algorithm".
 *  International Journal of Robotics Research (IJRR). vol. 30, no. 7, pp 936-953. Jun 2011.
 *
 *  \author mevenkamp
 */
template <typename RealType>
static inline void HungarianMethod ( const aol::FullMatrix<RealType> &AdjacencyMatrix, aol::Vector<int> &Memberships ) {
#ifdef USE_LIB_HUNGARIAN
  const int rows = AdjacencyMatrix.getNumRows ( ), cols = AdjacencyMatrix.getNumCols ( );
  const int k = aol::Max<int> ( rows, cols );
  HungarianMatrix m ( k, std::vector<Edge> ( k ) );
  for ( int i=0; i<k ; ++i )
    for ( int j=0; j<k ; ++j )
      m[i][i].SetWeight ( NEG_INF );
  for ( int i=0; i<rows; ++i )
    for ( int j=0; j<cols ; ++j )
      m[i][j].SetWeight ( -AdjacencyMatrix.get ( i, j ) );

  BipartiteGraph bg(m);
  Hungarian h(bg);
  h.HungarianAlgo();

  Memberships.setAll ( -1 );
  BipartiteGraph* bgFinal = h.GetBG();
  for(unsigned int i=0; i<bgFinal->GetNumAgents(); i++){
    for(unsigned int j=0; j<bgFinal->GetNumTasks(); j++){
      if(bgFinal->GetMatrix(i,j)->GetMatchedFlag() && m[i][j].GetWeight ( ) != NEG_INF ) Memberships[i] = j;
    }
  }
#else
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( AdjacencyMatrix );
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( Memberships );
  throw aol::Exception ( "Hungarian method requires library hungarian. Compile with -BUILD_AND_USE_HUNGARIAN", __FILE__, __LINE__ );
#endif
}

template <typename RealType>
class SegmentationQualityQuantifier {
  typedef qc::ScalarArray<int, qc::QC_2D> PictureType;
protected:
  PictureType _gt, _seg;
  const int _numPixels;
  int _N, _M, _K;
  aol::FullMatrix<int> _errorMatrix;
  aol::Vector<int> _iHat, _i;
  const RealType _overlapAcceptanceThreshold;
public:
  SegmentationQualityQuantifier ( const PictureType &GroundTruth, const PictureType &Segmentation )
    : _gt ( GroundTruth, aol::DEEP_COPY ), _seg ( Segmentation, aol::DEEP_COPY ), _numPixels ( GroundTruth.size ( ) ),
      _overlapAcceptanceThreshold ( 0.75 ) {

    transformGrayValuesToLabels ( _gt );
    transformGrayValuesToLabels ( _seg );
      
    _N = _gt.getMaxValue ( ) + 1;
    _M = _seg.getMaxValue ( ) + 1;
    _K = aol::Max<int> ( _N, _M );
        
    // Compute error matrix (n_{i,j}) (i.e. number of pixels interpreted as i-th class, but belonging to j-th class
    _errorMatrix.reallocate ( _K, _K );
    for ( int i=0; i<_M ; ++i ) {
      for ( int j=0; j<_N ; ++j ) {
        for ( int k=0; k<_numPixels ; ++k ) {
          if ( _seg[k] == i && _gt[k] == j )
            _errorMatrix.add ( i, j, 1 );
        }
      }
    }
        
    // Compute memberships (i.e. mapping of ground truth segments to classified labels)
    aol::FullMatrix<int> adjacencyMatrix ( _N, _M );
    for ( int gtLabel=0; gtLabel<_N ; ++gtLabel ) {
      for ( int segLabel=0; segLabel<_M ; ++segLabel ) {
        for ( int k=0; k<_numPixels ; ++k ) {
          if ( ( _gt[k] == gtLabel && _seg[k] != segLabel ) || ( _seg[k] == segLabel && _gt[k] != gtLabel ) )
            adjacencyMatrix.add ( gtLabel, segLabel, 1 );
        }
      }
    }
    aol::Vector<int> memberships ( _K );
    HungarianMethod ( adjacencyMatrix, memberships );

    _iHat.reallocate ( _K ), _i.reallocate ( _K );
    for ( int i=0; i<_K ; ++i ) {
      if ( memberships[i] >= 0 ) {
        _iHat[i] = memberships[i];
        _i[_iHat[i]] = i;
      }
    }
  }
  
  void saveSegmentationsWithMatchingLabels ( const char* outputdir ) {
    PictureType seg ( _seg );
    for ( int k=0; k<_numPixels ; ++k ) {
      if ( _i[_seg[k]] < _N )
        seg[k] = _i[_seg[k]];
    }
    
    _gt.setOverflowHandlingToCurrentValueRange ( );
    _gt.savePNG ( aol::strprintf ( "%s/gt.png", outputdir ).c_str ( ) );
    
    seg.setOverflowHandlingToCurrentValueRange ( );
    seg.savePNG ( aol::strprintf ( "%s/segmented.png", outputdir ).c_str ( ) );
  }
  
  void printReport ( const char* path ) {
    std::ofstream txtFile ( path );
    txtFile << "Report on segmentation quality" << std::endl << std::endl;
    
    txtFile << "Pixel-wise criteria" << std::endl;
    txtFile << "CS  = " << getCorrectDetection ( ) * 100.0 << "%" << std::endl;
    txtFile << "O   = " << getOmissionError ( ) * 100.0 << "%" << std::endl;
    txtFile << "C   = " << getCommissionError ( ) * 100.0 << "%" << std::endl;
    txtFile << "CA  = " << getClassAccuracy ( ) * 100 << "%" << std::endl;
    txtFile << "CO  = " << getCorrectAssignment ( ) * 100 << "%" << std::endl;
    txtFile << "CC  = " << getObjectAccuracy ( ) * 100 << "%" << std::endl;
    txtFile << std::endl;
  }
  
  /* Begin: region-wise criteria */
  
  RealType getCorrectDetection ( ) const {
    RealType cs = 0;
    for ( int i=0; i<_N ; ++i ) {
      const int numOverlap = getNumOverlap ( i );
      if ( numOverlap >= _overlapAcceptanceThreshold * getNIDot ( _iHat[i] ) && numOverlap >= _overlapAcceptanceThreshold * getNDotI ( i ) )
        ++cs;
    }
    cs /= static_cast<RealType> ( _N );
    return cs;
  }
  
  /* End: region-wise criteria */
  
  /* Begin: pixel-wise criteria */
  
  RealType getOmissionError ( ) const {
    aol::Vector<RealType> omissionErrors ( _N );
    for ( int i=0; i<_N ; ++i ) {
      RealType ndoti = getNDotI ( i );
      omissionErrors[i] = ( ndoti - _errorMatrix.get ( _iHat[i], i ) ) / static_cast<RealType> ( ndoti );
    }
    return omissionErrors.getMedianValue ( );
  }
  
  RealType getCommissionError ( ) const {
    aol::Vector<RealType> commissionErrors ( _M );
    for ( int iHat=0; iHat<_M ; ++iHat ) {
      const RealType nidot = getNIDot ( iHat );
      commissionErrors[iHat] = ( nidot - _errorMatrix.get ( iHat, _i[iHat] ) ) / static_cast<RealType> ( nidot );
    }
    return commissionErrors.getMedianValue ( );
  }
  
  RealType getClassAccuracy ( ) const {
    RealType ca = 0;
    for ( int i=0; i<_K ; ++i ) {
      const RealType ndoti = getNDotI ( i );
      const RealType nidot = getNIDot ( _iHat[i] );
      ca += _errorMatrix.get ( _iHat[i], i ) * ndoti / ( ndoti + nidot - _errorMatrix.get ( _iHat[i], i ) );
    }
    return ca / static_cast<RealType> ( _numPixels );
  }
  
  RealType getCorrectAssignment ( ) const {
    RealType co = 0;
    for ( int i=0; i<_K ; ++i ) co += _errorMatrix.get ( _iHat[i], i );
    return co / static_cast<RealType> ( _numPixels );
  }
  
  RealType getObjectAccuracy ( ) const {
    RealType cc = 0;
    for ( int i=0; i<_K ; ++i ) {
      const RealType nidot = getNIDot ( _iHat[i] );
      cc += ( nidot > 0 ) ? _errorMatrix.get ( _iHat[i], i ) * getNDotI ( i ) / nidot : 0.0;
    }
    return cc / static_cast<RealType> ( _numPixels );
  }
  
  /* End: pixel-wise criteria */
protected:
  RealType getNIDot ( const int I ) const {
    RealType nidot = 0;
    for ( int j=0; j<_N ; ++j ) nidot += _errorMatrix.get ( I, j );
    return nidot;
  }
  
  RealType getNDotI ( const int I ) const {
    RealType ndoti = 0;
    for ( int j=0; j<_M ; ++j ) ndoti += _errorMatrix.get ( j, I );
    return ndoti;
  }
  
  int getNumOverlap ( const int I ) const {
    int numOverlap = 0;
    for ( int k=0; k<_numPixels ; ++k ) {
      if ( _gt[k] == I && _seg[k] == _iHat[I] )
        ++numOverlap;
    }
    return numOverlap;
  }
};

/*
 * BEGIN: Operators suitable for segmentation
 */

template <typename _RealType, typename _PictureType>
class DOSTOp : public aol::Op<_PictureType, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const int _size, _sizeSqr;
  const int _blockOffset;
public:
  DOSTOp ( const int Size )
    : _size ( Size ), _sizeSqr ( Size * Size ),
      _blockOffset ( ( Size - 1 ) / 2 ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    Dest.resize ( 2 * _sizeSqr, Arg.size ( ) );
    qc::FastILexMapper<qc::QC_2D> mapper ( Arg.getNumX ( ), Arg.getNumY ( ) );
    qc::MultiArray<RealType, 2, 2> function ( _size, _size ), transform ( _size, _size );
    for ( int y=0; y<Arg.getNumY ( ) ; ++y ) {
      for ( int x=0; x<Arg.getNumX ( ) ; ++x ) {
        for ( int dy=-_blockOffset ; dy<=_blockOffset ; ++dy )
          for ( int dx=-_blockOffset ; dx<=_blockOffset ; ++dx )
            function[0].set ( dx + _blockOffset, dy + _blockOffset, Arg.getReflection ( x + dx, y + dy ) );
        qc::StockwellTransform<RealType> ( function, transform );
        transform /= static_cast<RealType> ( _sizeSqr );
        for ( int i=0; i<_sizeSqr ; ++i )
          for ( int k=0; k<2 ; ++k )
            Dest[k * _sizeSqr + i][mapper.getGlobalIndex ( x, y )] = transform[k][i];
      }
    }
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

template <typename _RealType, typename _PictureType>
class DOSTModulusOp : public aol::Op<_PictureType, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const int _size, _sizeSqr;
  const int _blockOffset;
public:
  DOSTModulusOp ( const int Size )
    : _size ( Size ), _sizeSqr ( Size * Size ),
      _blockOffset ( ( Size - 1 ) / 2 ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    Dest.resize ( _sizeSqr, Arg.size ( ) );
    qc::FastILexMapper<qc::QC_2D> mapper ( Arg.getNumX ( ), Arg.getNumY ( ) );
    qc::MultiArray<RealType, 2, 2> function ( _size, _size ), transform ( _size, _size );
    for ( int y=0; y<Arg.getNumY ( ) ; ++y ) {
      for ( int x=0; x<Arg.getNumX ( ) ; ++x ) {
        for ( int dy=-_blockOffset ; dy<=_blockOffset ; ++dy )
          for ( int dx=-_blockOffset ; dx<=_blockOffset ; ++dx )
            function[0].set ( dx + _blockOffset, dy + _blockOffset, Arg.getReflection ( x + dx, y + dy ) );
        qc::StockwellTransform<RealType> ( function, transform );
        transform /= static_cast<RealType> ( _sizeSqr );
        for ( int i=0; i<_sizeSqr ; ++i )
          Dest[i][mapper.getGlobalIndex ( x, y )] = sqrt ( aol::Sqr<RealType> ( transform[0][i] ) + aol::Sqr<RealType> ( transform[1][i] ) );
      }
    }
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

template <typename _RealType, typename _PictureType>
class FFTModulusOp : public aol::Op<_PictureType, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const int _size, _sizeSqr;
  const int _blockOffset;
public:
  FFTModulusOp ( const int Size )
    : _size ( Size ), _sizeSqr ( Size * Size ),
      _blockOffset ( ( Size - 1 ) / 2 ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    Dest.resize ( _sizeSqr, Arg.size ( ) );
    qc::FastILexMapper<qc::QC_2D> mapper ( Arg.getNumX ( ), Arg.getNumY ( ) );
    for ( int y=0; y<Arg.getNumY ( ) ; ++y ) {
      for ( int x=0; x<Arg.getNumX ( ) ; ++x ) {
        PictureType block ( _size ), modulus ( _size );
        for ( int dy=-_blockOffset ; dy<=_blockOffset ; ++dy )
          for ( int dx=-_blockOffset ; dx<=_blockOffset ; ++dx )
            block.set ( dx + _blockOffset, dy + _blockOffset, Arg.getReflection ( x + dx, y + dy ) );
        qc::computeLogFFTModulus<RealType> ( block, modulus, 0, false );
        modulus /= static_cast<RealType> ( _sizeSqr );
        for ( int k=0; k<modulus.size ( ) ; ++k )
          Dest[k][mapper.getGlobalIndex ( x, y )] = modulus[k];
      }
    }
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

template <typename _PictureType>
class FilterConvolutionOp : public aol::Op<_PictureType> {
  typedef _PictureType PictureType;
  typedef typename PictureType::RealType RealType;
protected:
  const int _size, _offset;
  aol::FullMatrix<RealType> _M;
public:
  FilterConvolutionOp ( const int Size ) : _size ( Size ), _offset ( ( Size - 1 ) / 2 ), _M ( Size, Size ) { }
  
  virtual void apply ( const PictureType &Arg, PictureType &Dest ) const {
    Dest.setZero ( );
    for ( int y=0; y<Arg.getNumY ( ) ; ++y )
      for ( int x=0; x<Arg.getNumX ( ) ; ++x )
        for ( int dy=-_offset; dy<=_offset ; ++dy )
          for ( int dx=-_offset; dx<=_offset ; ++dx )
            Dest.add ( x, y, _M.get ( dy + _offset , dx + _offset ) * Arg.getReflection ( x + dx, y + dy ) );
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, PictureType &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  const aol::FullMatrix<RealType>& getFilterKernel ( ) const {
    return _M;
  }
};

template <typename _RealType, typename _PictureType> class FilterResponseLocalHistogramOp;

template <typename _RealType>
class FilterResponseLocalHistogramOp<_RealType, qc::ScalarArray<_RealType, qc::QC_2D> > : public aol::Op<qc::ScalarArray<_RealType, qc::QC_2D>, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> PictureType;
protected:
  const FilterConvolutionOp<PictureType> *_filter;
  const int _integrationScale, _integrationScaleSqr, _offset;
  const int _numBins;
  RealType _minVal, _maxVal, _numBinsInvMinMaxDiff;
public:
  FilterResponseLocalHistogramOp ( const FilterConvolutionOp<PictureType> *Filter, const int IntegrationScale, const int NumBins )
    : _filter ( Filter ),
      _integrationScale ( IntegrationScale ),
      _integrationScaleSqr ( aol::Sqr<RealType> ( IntegrationScale ) ),
      _offset ( ( IntegrationScale - 1 ) / 2 ), _numBins ( NumBins ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    Dest.resize ( _numBins, Arg.size ( ) );
    qc::FastILexMapper<qc::QC_2D> mapper ( Arg.getNumX ( ), Arg.getNumY ( ) );
    PictureType response ( Arg, aol::STRUCT_COPY );
    _filter->apply ( Arg, response );
    RealType minVal = response.getMinValue ( ), maxVal = response.getMaxValue ( );
    RealType numBinsInvMinMaxDiff = static_cast<RealType> ( _numBins ) / ( maxVal - minVal );
    aol::Vector<RealType> histogram ( _numBins );
    for ( int y=0; y<response.getNumY ( ) ; ++y ) {
      for ( int x=0; x<response.getNumX ( ) ; ++x ) {
        histogram.setZero ( );
        int binIdx;
        for ( int dy=-_offset ; dy<=_offset ; ++dy ) {
          for ( int dx=-_offset ; dx<=_offset ; ++dx ) {
            binIdx = floor ( ( response.getClip ( x + dx, y + dy ) - minVal ) * numBinsInvMinMaxDiff );
            if ( binIdx < 0 ) binIdx = 0;
            if ( binIdx >= _numBins ) binIdx = _numBins - 1;
            histogram[binIdx] += 1.0;
          }
        }
        histogram /= static_cast<RealType> ( _integrationScaleSqr );
        for ( int k=0; k<_numBins ; ++k )
          Dest[k][mapper.getGlobalIndex ( x, y )] = histogram[k];
      }
    }
  }
  
  int getNumBins ( ) const {
    return _numBins;
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

template <typename _RealType>
class FilterResponseLocalHistogramOp<_RealType, qc::MultiArray<_RealType, qc::QC_2D, 3> > : public aol::Op<qc::MultiArray<_RealType, qc::QC_2D, 3>, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> PictureType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> ScalarPictureType;
protected:
  const FilterResponseLocalHistogramOp<RealType, ScalarPictureType> _scalarFilterResponceLocalHistogramOp;
  const bool _convertToGray;
public:
  FilterResponseLocalHistogramOp ( const FilterConvolutionOp<ScalarPictureType> *Filter, const int IntegrationScale, const int NumBins,
                                   const bool ConvertToGray = true )
    : _scalarFilterResponceLocalHistogramOp ( Filter, IntegrationScale, NumBins ),
      _convertToGray ( ConvertToGray ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    const int numBins = _scalarFilterResponceLocalHistogramOp.getNumBins ( );
    if ( _convertToGray ) {
      ScalarPictureType gray ( Arg[0] ); // Assume CIELab color space, where first coordinate is luminance/gray-scale
      Dest.resize ( numBins, gray.size ( ) );
      _scalarFilterResponceLocalHistogramOp.apply ( gray, Dest );
    } else {
      Dest.resize ( 3 * numBins, Arg[0].size ( ) );
      aol::MultiVector<RealType> tmp ( numBins, Arg[0].size ( ) );
      for ( int c=0; c<3 ; ++ c ) {
        _scalarFilterResponceLocalHistogramOp.apply ( Arg[c], tmp );
        for ( int i=0; i<numBins ; ++i ) Dest[c * numBins + i] = tmp[i];
      }
    }
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

template <typename _RealType, typename _PictureType>
class BaseLocalSpectralHistogramOp : public aol::Op<_PictureType, aol::MultiVector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  std::vector<FilterResponseLocalHistogramOp<_RealType, _PictureType>*> _histOps;
public:
  BaseLocalSpectralHistogramOp ( ) { }
  
  void apply ( const PictureType &Arg, aol::MultiVector<RealType> &Dest ) const {
    for ( int i=0; i<static_cast<int>( _histOps.size ( ) ) ; ++i ) {
      aol::MultiVector<RealType> *dest = new aol::MultiVector<RealType> ( );
      _histOps[i]->apply ( Arg, *dest );
      Dest.appendReference ( *dest, true );
    }
  }
  
  virtual void applyAdd ( const PictureType &/*Arg*/, aol::MultiVector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void pushBack ( const FilterConvolutionOp<PictureType> *Filter, const int IntegrationScale, const int NumBins ) {
    _histOps.push_back ( new FilterResponseLocalHistogramOp<RealType, PictureType> ( Filter, IntegrationScale, NumBins ) );
  }
};

template <typename _RealType, typename _PictureType>
class LocalSpectralHistogramOp : public BaseLocalSpectralHistogramOp<_RealType, _PictureType> { };

template <typename _RealType>
class LocalSpectralHistogramOp<_RealType, qc::MultiArray<_RealType, qc::QC_2D, 3> > : public BaseLocalSpectralHistogramOp<_RealType, qc::MultiArray<_RealType, qc::QC_2D, 3> > {
public:
  void pushBack ( const FilterConvolutionOp<qc::ScalarArray<_RealType, qc::QC_2D> > *Filter, const int IntegrationScale, const int NumBins,
                  const bool ConvertToGray = true ) {
    BaseLocalSpectralHistogramOp<_RealType, qc::MultiArray<_RealType, qc::QC_2D, 3> >::_histOps.push_back ( new FilterResponseLocalHistogramOp<_RealType, qc::MultiArray<_RealType, qc::QC_2D, 3> > ( Filter, IntegrationScale, NumBins, ConvertToGray ) );
  }
};

template <typename RealType>
static void histogramEqualization ( qc::ScalarArray<RealType, qc::QC_2D> &ArgDest ) {
  const int argMax = ArgDest.getMaxValue ( );
  aol::Vector<int> frequencies ( argMax + 1 ), lookUp ( argMax + 1 );
  for ( int k=0; k<ArgDest.size ( ) ; ++k ) ++frequencies[static_cast<int> ( ArgDest[k] )];
  for ( int i=0; i<=argMax ; ++i ) {
    for ( int j=0; j<=i ; ++j ) lookUp[i] += frequencies[j];
    lookUp[i] *= argMax / static_cast<RealType> ( ArgDest.size ( ) );
  }
  for ( int k=0; k<ArgDest.size ( ) ; ++k ) ArgDest[k] = lookUp[static_cast<int> ( ArgDest[k] )];
}

template <typename _PictureType>
class IntensityFilterOp : public FilterConvolutionOp<_PictureType> {
  typedef _PictureType PictureType;
  typedef typename PictureType::RealType RealType;
public:
  IntensityFilterOp ( ) : FilterConvolutionOp<PictureType> ( 1 ) {
    this->_M.set ( 0, 0, 1.0 );
  }
};

template <typename _PictureType>
class GaborFilterOp : public FilterConvolutionOp<_PictureType> {
  typedef _PictureType PictureType;
  typedef typename PictureType::RealType RealType;
protected:
  const RealType _s;
  const RealType _theta;
  const RealType _phase;
  const RealType _gamma;
  const RealType _f;
public:
  GaborFilterOp ( const int Size, const RealType Theta )
    : FilterConvolutionOp<PictureType> ( Size ),
      _s ( 0.5 * this->_offset ), _theta ( Theta ), _phase ( 0.0 ), _gamma ( 1.0 ),
      _f ( 1.0 / ( 2.0 * _s ) ) {
    RealType sSqr = aol::Sqr<RealType> ( _s );
    RealType xPrime, yPrime;
    for ( int y=-this->_offset; y<=this->_offset ; ++y ) {
      for ( int x=-this->_offset; x<=this->_offset ; ++x ) {
        xPrime = x * cos(_theta) + y * sin(_theta);
        yPrime = y * cos(_theta) - x * sin(_theta);
        this->_M.set ( x+this->_offset, y+this->_offset, 1.0 / ( 2.0 * aol::NumberTrait<RealType>::pi * sSqr ) * exp ( -0.5 * ( aol::Sqr<RealType> ( xPrime )
                                                         + aol::Sqr<RealType> ( yPrime * _gamma ) ) / sSqr ) * cos ( 2.0 * aol::NumberTrait<RealType>::pi * _f * xPrime + _phase ) );
      }
    }
  }
};

template <typename _PictureType>
class GaussianFilterOp : public FilterConvolutionOp<_PictureType> {
  typedef _PictureType PictureType;
  typedef typename PictureType::RealType RealType;
protected:
  const RealType _sigma, _minusInvTwoSigmaSqr;
public:
  GaussianFilterOp ( const int Size, const RealType Sigma )
    : FilterConvolutionOp<PictureType> ( Size ),
      _sigma ( Sigma ), _minusInvTwoSigmaSqr ( -1.0 / ( 2.0 * aol::Sqr<RealType> ( Sigma ) ) ) {
    RealType sum = 0;
    for ( int y=-this->_offset; y<=this->_offset ; ++y ) {
      for ( int x=-this->_offset; x<=this->_offset ; ++x ) {
        this->_M.set ( x+this->_offset, y+this->_offset, exp ( ( aol::Sqr<RealType> ( x ) + aol::Sqr<RealType> ( y ) ) * _minusInvTwoSigmaSqr ) );
        sum += this->_M.get ( x+this->_offset, y+this->_offset );
      }
    }
    this->_M *= 1.0 / sum;
  }
  
  RealType getMinusInvTwoSigmaSqr ( ) const {
    return _minusInvTwoSigmaSqr;
  }
};

template <typename _PictureType>
class LaplacianOfGaussianFilterOp : public FilterConvolutionOp<_PictureType> {
  typedef _PictureType PictureType;
  typedef typename PictureType::RealType RealType;
protected:
  const RealType _sigma, _sigmaSqr;
public:
  LaplacianOfGaussianFilterOp ( const int Size, const RealType Sigma )
    : FilterConvolutionOp<PictureType> ( Size ),
      _sigma ( Sigma ), _sigmaSqr ( aol::Sqr<RealType> ( Sigma ) ) {
    GaussianFilterOp<PictureType> gaussianFilterOp ( this->_size, _sigma );
    for ( int y=-this->_offset; y<=this->_offset ; ++y )
      for ( int x=-this->_offset; x<=this->_offset ; ++x )
        this->_M.set ( x+this->_offset, y+this->_offset, gaussianFilterOp.getFilterKernel ( ).get ( x+this->_offset, y+this->_offset )
                                                         * ( 1 + ( aol::Sqr<RealType> ( x ) + aol::Sqr<RealType> ( y ) ) * gaussianFilterOp.getMinusInvTwoSigmaSqr ( ) ) );
    this->_M *= -1.0 / ( aol::NumberTrait<RealType>::pi * aol::Sqr<RealType> ( _sigmaSqr ) );
  }
};

/*
 * END: Local spectral histogram ops
 */


template <typename RealType>
static void saveOptimalRegionIndicators ( const char* baseName, const qc::ScalarArray<int, qc::QC_2D> &GroundTruth, const aol::MultiVector<RealType> &Indicator ) {
  qc::ScalarArray<int, qc::QC_2D> groundTruth ( GroundTruth, aol::DEEP_COPY );
  transformGrayValuesToLabels ( groundTruth );
  for ( int l=0; l<=groundTruth.getMaxValue ( ) ; ++l ) {
    int n = 0;
    aol::Vector<RealType> mean ( Indicator.numComponents ( ) );
    for ( int i=0; i<groundTruth.size ( ) ; ++i ) {
      if ( groundTruth[i] == l ) {
        for ( int c=0; c<Indicator.numComponents ( ) ; ++c )
          mean[c] += Indicator[c][i];
        ++n;
      }
    }
    mean /= static_cast<RealType> ( n );

    qc::ScalarArray<RealType, qc::QC_2D> errorImg ( groundTruth.getNumX ( ), groundTruth.getNumY ( ) );
    aol::Vector<RealType> error ( Indicator.numComponents ( ) );
    for ( int i=0; i<groundTruth.size ( ) ; ++i ) {
      error = mean;
      for ( int c=0; c<Indicator.numComponents ( ) ; ++c )
        error[c] -= Indicator[c][i];
      errorImg[i] = error.norm ( );
    }
    errorImg.setOverflowHandlingToCurrentValueRange ( );
    errorImg.savePNG ( aol::strprintf ( "%s_%d.png", baseName, l ).c_str ( ) );
  }
}


/**
 * A class for PCA-based segmentation of high-dimensional piece-wise constant vector fields
 * e.g. arising in texture segmentation based on local spectral histograms or unitary transforms (e.g. Stockwell)
 *
 * As BaseClass, use either PiecewiseConstantMultiPhaseMSSegmentor<ConfiguratorType, FirstOrderPrimalDualMultiPhaseMSSegmentor<ConfiguratorType> > (multi phase)
 *                       or PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, FirstOrderPrimalDualTwoPhaseMSSegmentor<ConfiguratorType> > (two phase)
 *
 * \author Mevenkamp
 */
template <typename ConfiguratorType, typename BaseClass>
class PCAMSSegmentor : public BaseClass {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  const aol::MultiVector<RealType> _highDimImageMVec;
  aol::MultiVector<RealType> _eigenVecs;
  aol::Vector<RealType> _eigenVals;
  const int _numEvals;
  const RealType _omega;
  const bool _orthonormalizeMeanValues;
  std::string _outputDir;
  qc::ScalarArray<int, qc::QC_2D> _groundTruth;
public:
  PCAMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                   const RealType Gamma,
                   aol::MultiVector<RealType> &ImageMVec,
                   const bool InitializeGrayValues = true,
                   const int NumSegments = 0,
                   const int NumEvals = 0,
                   const RealType Omega = 0.05,
                   const bool OrthonormalizeMeanValues = false,
                   const int NumGhostCells = 0,
                   const std::string &OutputDir = "" )
    : BaseClass ( Initializer, Gamma, ImageMVec, false ),
      _highDimImageMVec ( ImageMVec ), _numEvals ( ( NumEvals == 0 ) ? NumSegments : NumEvals ), _omega ( Omega ),
      _orthonormalizeMeanValues ( OrthonormalizeMeanValues ), _outputDir ( OutputDir ) {
    this->setNumGhostCells ( NumGhostCells );
    setEigenVecs ( );
    this->setNumSegments ( ( NumSegments > 0 ) ? NumSegments : _eigenVecs.numComponents ( ) );
    if ( NumSegments == 0 && !this->_quietMode ) std::cerr << "PCA detected " << this->_numSegments << " segments." << std::endl;
    setProjectedCoefficients ( );
    this->_meanValues.resize ( this->_numSegments, this->_imageMVec.numComponents ( ) );
    if ( InitializeGrayValues ) initializeGrayValues ( );
        
    if ( _outputDir != "" && _groundTruth.size ( ) > 0 )
      saveOptimalRegionIndicators ( aol::strprintf ( "%s/pcaRegionIndicator", _outputDir.c_str ( ) ).c_str ( ), _groundTruth, this->_imageMVec );
  }
  
  static void getEigenVecs ( aol::MultiVector<RealType> &EigenVecs, aol::Vector<RealType> &EigenVals, const aol::MultiVector<RealType> &HighDimInput, const int NumEvals, const RealType Omega = 0 ) {
    // Setup system matrix
    const int K = HighDimInput.numComponents ( );
    const int N = HighDimInput[0].size ( ) * log ( static_cast<RealType> ( K ) ) / static_cast<RealType> ( K );
    
    // Compute sample mean from data
    aol::Vector<RealType> sampleMean ( K );
    HighDimInput.getMeanComponents ( sampleMean );
    
    // Assemble covariance matrix from N random samples
    aol::MultiVector<RealType> samples ( K, N );
    aol::RandomGenerator randomGenerator;
//    randomGenerator.randomize ( );    // if "true" randomness is required, uncomment this line
    aol::Vector<int> indices ( N );
    randomGenerator.rIntVecPairwiseDifferent ( indices, 0, HighDimInput[0].size ( ) );
    for ( int row=0; row<N ; ++row )
      for ( int col=0; col<K ; ++col )
        samples[col][row] = HighDimInput[col][indices[row]] - sampleMean[col];
    
    // Assemble AA^T
    aol::ProgressBar<> progressBar ( "Assembling ATA" );
    progressBar.start ( ceil ( K * ( K + 1 ) / 2 ) );
    aol::FullMatrix<RealType> AAT ( K, K );
    for ( int i=0; i<K ; ++i ) {
      for ( int j=0; j<=i ; ++j ) {
        AAT.set ( i, j, samples[i].dotProduct ( samples[j] ) );
        progressBar++;
      }
    }
    progressBar.finish ( );
    for ( int i=0; i<K ; ++i )
      for ( int j=i+1 ; j<K ; ++j )
        AAT.set ( i, j, AAT.get ( j, i ) );
    
    // Compute eigenvalues
    aol::MultiVector<RealType> eig;
    aol::DeflationEigenvectorOp<aol::FullMatrix<RealType> > evalOp ( NumEvals );
    RealType sumSingularValsSqr = 0;
    for ( int row = 0; row < N; ++row )
      for ( int col = 0; col < K; ++col )
        sumSingularValsSqr += aol::Sqr<RealType> ( samples[col][row] );
    if ( NumEvals == 0 )
      evalOp.setThreshold ( Omega * N, sumSingularValsSqr );
    evalOp.apply ( AAT, eig );
    EigenVecs.resize ( eig.numComponents ( ) - 1, K );
    for ( int ev = 0; ev<eig.numComponents ( ) - 1; ++ev )
      for ( int i = 0; i<K; ++i )
        EigenVecs[ev][i] = eig[ev + 1][i];
    EigenVals.resize ( eig.numComponents ( ) ); // DeflationEigenvectorOp does not return eigenvalues
    EigenVals[0] = sumSingularValsSqr;
    for ( int ev = 1; ev < eig.numComponents ( ); ++ev )
      EigenVals[ev] = eig[0][ev];

    //aol::EigenLibraryInterfaceOp<RealType, aol::FullMatrix<RealType> >::getEigenVectorsRealSymmetric ( AAT, EigenVecs, EigenVals );
    //
    //// Determine number of eigen values / vectors to be used (if necessary)
    //if ( NumEvals > 0 ) EigenVecs.resize ( NumEvals, K );
    //else {
    //  int numEvals = EigenVals.size ( );
    //  RealType lse = 0;
    //  while ( numEvals > 1 && lse < Omega * N ) {
    //    --numEvals;
    //    lse += EigenVals[numEvals];
    //  }
    //  EigenVecs.resize ( numEvals+1, K );

    //  std::cerr << lse << std::endl;
    //  std::cerr << N << std::endl;

    //  std::cerr << "PCAMSSegmentor: detected " << numEvals + 1 << " segments." << std::endl;
    //}
  }
  
  static void getProjectedCoefficients ( aol::MultiVector<RealType> &LowDimOutput, const aol::MultiVector<RealType> &HighDimInput, const aol::MultiVector<RealType> &EigenVecs ) {
    const int K = HighDimInput.numComponents ( );
    
    // Compute sample mean from data
    aol::Vector<RealType> sampleMean ( K );
    HighDimInput.getMeanComponents ( sampleMean );
    
    // Project mean-centralized high-dimensional input onto specified eigen-vector basis
    aol::ProgressBar<> progressBar ( "Calculating projected coefficients" );
    progressBar.start ( HighDimInput[0].size ( ) );
    LowDimOutput.reallocate ( EigenVecs.numComponents ( ), HighDimInput[0].size ( ) );
    aol::Vector<RealType> highDimDatum ( HighDimInput.numComponents ( ) );
    for ( int k=0; k<HighDimInput[0].size ( ) ; ++k ) {
      HighDimInput.getTo ( k, highDimDatum );
      highDimDatum -= sampleMean;
      for ( int ev=0; ev<EigenVecs.numComponents ( ) ; ++ev )
        LowDimOutput[ev][k] = highDimDatum.dotProduct ( EigenVecs[ev] );
      progressBar++;
    }
    progressBar.finish ( );
  }
protected:
  void setEigenVecs ( ) {
    getEigenVecs ( _eigenVecs, _eigenVals, _highDimImageMVec, _numEvals, _omega );
  }
  
  void setProjectedCoefficients ( ) const {
    getProjectedCoefficients ( this->_imageMVec, _highDimImageMVec, _eigenVecs );
  }
  
  void orthoNormalizeGrayValues ( ) {
    if ( !this->_quietMode ) std::cerr << "Orthonormalizing basis representation of mean values.." << std::endl;
    // Assemble matrix A containing the coefficient representation of mean values in eigen vector basis
    aol::FullMatrix<RealType> A ( _eigenVecs.numComponents ( ), this->_numSegments );
    for ( int l=0; l<this->_numSegments ; ++l )
      for ( int ev=0; ev<_eigenVecs.numComponents ( ) ; ++ev )
        A.set ( ev, l, this->_meanValues[l][ev] );
    
    // Perform orthonormalization of mean values
    aol::FullMatrix<RealType> Q ( _eigenVecs.numComponents ( ), this->_numSegments ), R ( this->_numSegments, this->_numSegments );
    aol::QRDecomposeModifiedGramSchmidt<RealType> qrGramSchmidt;
    qrGramSchmidt.transform ( A, R, Q );
    for ( int l=0; l<this->_numSegments ; ++l )
      for ( int ev=0; ev<_eigenVecs.numComponents ( ) ; ++ev )
        this->_meanValues[l][ev] = Q.get ( ev, l );
    
    // Apply change of basis to eigen vectors
    aol::FullMatrix<RealType> QT ( Q.getNumCols ( ), Q.getNumRows ( ) );
    Q.transposeTo ( QT );
    aol::FullMatrix<RealType> AQT ( A.getNumRows ( ), A.getNumRows ( ) );
    AQT.makeProduct ( A, QT );
    aol::QRInverse<RealType> qrInv ( AQT );
    aol::FullMatrix<RealType> B = qrInv.getFM ( );
    aol::Vector<RealType> vec ( _eigenVecs.numComponents ( ) ), transformedVec ( _eigenVecs.numComponents ( ) );
    for ( int k=0; k<_eigenVecs[0].size ( ) ; ++k ) {
      _eigenVecs.getTo ( k, vec );
      B.apply ( vec, transformedVec );
      _eigenVecs.set ( k, transformedVec );
    }
    
    getProjectedCoefficients ( this->_imageMVec, _highDimImageMVec, _eigenVecs );
  }

  virtual void initializeGrayValues ( ) {
    aol::KMeansClusterer<RealType> kMeansClusterer;
    aol::MultiVector<RealType> clusters;
    kMeansClusterer.applyMultipleRNG ( this->_imageMVec, clusters, this->_numSegments );
    const int imageDim = this->_imageMVec.numComponents ( );
    for ( int l=0; l<this->_numSegments ; ++l )
      for ( int j=0; j<imageDim ; ++j )
        this->_meanValues[l][j] = clusters[j][l];
    if ( !this->_quietMode ) cerr << this->_meanValues << endl;
    
    if ( _orthonormalizeMeanValues ) orthoNormalizeGrayValues ( );
  }
  virtual void updateGrayValues ( const typename BaseClass::ArrayType &CurrentSegmentation ) {
    BaseClass::updateGrayValues ( CurrentSegmentation );
    if ( _orthonormalizeMeanValues ) orthoNormalizeGrayValues ( );
  }
public:
  void setMeanValues ( const aol::Vector<int> &Indices ) {
    BaseClass::setMeanValues ( Indices );
    if ( _orthonormalizeMeanValues ) orthoNormalizeGrayValues ( );
  }
  
  void setGroundTruth ( const qc::ScalarArray<int, qc::QC_2D> &GroundTruth ) {
    _groundTruth.reallocate ( GroundTruth.getSize ( ) );
    _groundTruth = GroundTruth;
  }
};


class SEGMENTATION_BOUNDARY_CONDITIONS {
  
};


template<typename _RealType, typename _PictureType>
class LocalFeatureMSSegmentationOptions {
public:
  // I/O
  const _PictureType &input;
  std::string outputDir;
  bool quietMode;
  aol::ProgressBar<> *progressBar;
  qc::ScalarArray<int, qc::QC_2D> groundTruth;
  
  // General parameters
  short resampleFactor;
  int numSegments;
  int numGhostCells;
  bool includeBoundary;
  _RealType lowEdgenessThreshold;
  aol::MultiVector<_RealType> seedPositions;
  
  // Parameters for PCA
  int numEvals;
  _RealType omega;
  bool orthonormalizeMeanValues;

  // Parameters for Clustering
  int numDominantOperators;
  int transformScale;
  bool forceMultiPhaseSegmentation;
  
  bool forceGroundTruthClusters; // for testing purposes
  
  // Parameters for Mumford-Shah
  _RealType gamma;
  int maxIt;
  _RealType epsilon;
  int numOuterIterations;
  
  // Operators
  std::vector<aol::Op<_PictureType, aol::MultiVector<_RealType> >*> operators;
  aol::Vector<_RealType> operatorWeights;
  
  LocalFeatureMSSegmentationOptions ( const _PictureType &Input, const std::string OutputDir = "", const bool Quiet = true )
    : input ( Input ), outputDir ( OutputDir ), quietMode ( Quiet ), progressBar ( NULL ),
      resampleFactor ( 1 ), numSegments ( 0 ), numGhostCells ( 0 ),
      includeBoundary ( true ), lowEdgenessThreshold ( 0.5 ), seedPositions ( ),
      numEvals ( 0 ), omega ( 0.05 ), orthonormalizeMeanValues ( false ),
      numDominantOperators ( 1 ), transformScale ( 21 ), forceMultiPhaseSegmentation ( false ),
      forceGroundTruthClusters ( false ),
      gamma ( gammaDefault ), maxIt ( maxItDefault ), epsilon ( epsilonDefault ), numOuterIterations ( 3 ) { }
  
  LocalFeatureMSSegmentationOptions ( const LocalFeatureMSSegmentationOptions<_RealType, _PictureType> &Options )
    : input ( Options.input ), outputDir ( Options.outputDir ), quietMode ( Options.quietMode ), progressBar ( Options.progressBar ),
      resampleFactor ( Options.resampleFactor ), numSegments ( Options.numSegments ), numGhostCells ( Options.numGhostCells ),
      includeBoundary ( Options.includeBoundary ), lowEdgenessThreshold ( Options.lowEdgenessThreshold ), seedPositions ( Options.seedPositions ),
      numEvals ( Options.numEvals ), omega ( Options.omega ), orthonormalizeMeanValues ( Options.orthonormalizeMeanValues ),
      numDominantOperators ( Options.numDominantOperators ), transformScale ( Options.transformScale ), forceMultiPhaseSegmentation ( Options.forceMultiPhaseSegmentation ),
      forceGroundTruthClusters ( Options.forceGroundTruthClusters ),
      gamma ( Options.gamma ), maxIt ( Options.maxIt ), epsilon ( Options.epsilon ), numOuterIterations ( Options.numOuterIterations ) { }
  
  static _RealType gammaDefault, epsilonDefault;
  static int maxItDefault;
};
template<typename _RealType, typename _PictureType> _RealType LocalFeatureMSSegmentationOptions<_RealType, _PictureType>::gammaDefault = 1e-4;
template<typename _RealType, typename _PictureType> _RealType LocalFeatureMSSegmentationOptions<_RealType, _PictureType>::epsilonDefault = 1e-3;
template<typename _RealType, typename _PictureType> int LocalFeatureMSSegmentationOptions<_RealType, _PictureType>::maxItDefault = 10000;

template <typename _RealType, typename _PictureType>
class LocalFeatureMSSegmentor {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> SegmentationType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfiguratorType;
  typedef typename ConfiguratorType::InitType InitType;
  typedef typename qc::ComponentsCollection<RealType>::NonEmptyComponentsIterator NonEmptyComponentsIterator;
  typedef LocalFeatureMSSegmentationOptions<RealType, PictureType> OptionsType;
  typedef qc::PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, 0, qc::FirstOrderPrimalTwoPhaseMSSegmentor<ConfiguratorType> > TwoPhaseSegmentorType;
  typedef qc::PiecewiseConstantMultiPhaseMSSegmentor<ConfiguratorType, qc::FirstOrderPrimalDualMultiPhaseMSSegmentor<ConfiguratorType> > MultiPhaseSegmentorType;
protected:
  std::string _outputDir;
  aol::ProgressBar<> *_progressBar;
  bool _quietMode;
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
  int _numDominantComponents;
public:
  LocalFeatureMSSegmentor ( ) : _outputDir ( "" ), _quietMode ( false ), _catchCtrlC ( false ) { }
  
  virtual ~LocalFeatureMSSegmentor ( ) { }
  
  void apply ( OptionsType &Options, qc::ScalarArray<int, qc::QC_2D> &Segmentation ) {
    Segmentation.reallocate ( Options.input.getNumX ( ), Options.input.getNumY ( ) );
    initialize ( Options );
    getSegmentationFromLocalFeatures ( Options, Segmentation );
  }
  
  void getLocalFeatures ( const OptionsType &Options, aol::MultiVector<RealType> &LocalFeatures ) {
    InitType *grid = NULL;
    getLocalFeatures ( Options, LocalFeatures, grid );
  }
  
  void setQuietMode ( bool Quiet ) {
    _quietMode = Quiet;
  }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
protected:
  void initialize ( const OptionsType &Options ) {
    _outputDir = Options.outputDir;
    _quietMode = Options.quietMode;
    _progressBar = Options.progressBar;
    
    // Copy input to destination
    if ( !_quietMode && _outputDir != "" ) {
      std::stringstream ss;
      ss << _outputDir << "/input.png";
      PictureType u ( Options.input );
      u.scaleValuesTo01 ( );
      u.setOverflowHandling ( aol::SCALE, 0, 1 );
      u.savePNG ( ss.str ( ).c_str ( ) );
    }
  }
  
  void getSegmentationFromLocalFeatures ( OptionsType &Options, qc::ScalarArray<int, qc::QC_2D> &Dest ) {
    aol::MultiVector<RealType> localFeatures;
    InitType *grid = NULL;
    getLocalFeatures ( Options, localFeatures, grid );
    getHardSegmentation ( Options, Dest, *grid, localFeatures );
    upsampleHardSegmentation ( Dest, aol::Vec2<short> ( Options.input.getNumX ( ), Options.input.getNumY ( ) ) );
    
    if ( !_quietMode && _outputDir != "" )
      plotSegmentationBoundariesOntoImage<RealType, PictureType> ( aol::strprintf ( "%s/boundaries", _outputDir.c_str ( ) ).c_str ( ), Options.input, Dest );
  }
  
  void getLocalFeatures ( const OptionsType &Options, aol::MultiVector<RealType> &LocalFeatures, InitType *&Grid ) {
    // Resample image
    PictureType u ( Options.input.getNumX ( ) / Options.resampleFactor, Options.input.getNumY ( ) / Options.resampleFactor );
    u.resampleFrom ( Options.input );
    
    // Post-process operator weights
    aol::Vector<RealType> operatorWeights ( static_cast<int>( Options.operators.size ( ) ) );
    RealType sum = Options.operatorWeights.sum ( ) + Options.operators.size ( ) - Options.operatorWeights.size ( );
    operatorWeights.setAll ( 1.0 / sum );
    for ( int k=0; k<Options.operatorWeights.size ( ) ; ++k ) operatorWeights[k] = Options.operatorWeights[k] / sum;
    
    // Apply specified operators and assemble MultiVector of local features
    qc::FastILexMapper<qc::QC_2D> mapper ( u.getNumX ( ), u.getNumY ( ) );
    int offset = ( Options.includeBoundary ) ? 0 : ( Options.transformScale - 1 ) / 2;
    Grid = new InitType ( qc::GridSize<qc::QC_2D> ( u.getNumX ( ) - 2 * offset, u.getNumY ( ) - 2 * offset ) );
    aol::MultiVector<RealType> localFeatures;
    if ( this->_progressBar != NULL ) {
      this->_progressBar->setText ( "Computing local descriptors" );
      this->_progressBar->start ( Options.operators.size ( ) );
    }
    _numDominantComponents = 0;
    for ( int k=0; k<static_cast<int>(Options.operators.size ( )) ; ++k ) {
      Options.operators[k]->apply ( u, localFeatures );
      localFeatures *= operatorWeights[k];
      if ( k < Options.numDominantOperators ) _numDominantComponents += localFeatures.numComponents ( );
      int c0 = LocalFeatures.numComponents ( );
      LocalFeatures.resize ( c0 + localFeatures.numComponents ( ), Grid->getNumberOfNodes ( ) );
      int j = 0;
      for ( int y=offset; y<u.getNumY ( )-offset ; ++y ) {
        for ( int x=offset; x<u.getNumX ( )-offset ; ++x ) {
          for ( int c=0; c<localFeatures.numComponents ( ) ; ++c )
            LocalFeatures[c0+c][j] = localFeatures[c][mapper.getGlobalIndex ( x, y )];
          ++j;
        }
      }
      if ( this->_progressBar != NULL ) (*this->_progressBar)++;
    }
    if ( this->_progressBar != NULL ) this->_progressBar->finish ( );
    
    if ( !_quietMode && _outputDir != "" ) {
      // Save resampled input image
      std::stringstream ss;
      ss << _outputDir << "/input_resampled.png";
      u.scaleValuesTo01 ( );
      u.setOverflowHandling ( aol::SCALE, 0, 1 );
      u.savePNG ( ss.str ( ).c_str ( ) );
    }
    
    if ( !_quietMode && _outputDir != "" && Options.groundTruth.size ( ) > 0 )
      saveOptimalRegionIndicators ( aol::strprintf ( "%s/highDimFeatureIndicator", _outputDir.c_str ( ) ).c_str ( ), Options.groundTruth, LocalFeatures );
  }
  
  void getHardSegmentation ( OptionsType &Options, qc::ScalarArray<int, qc::QC_2D> &HardSegmentation, const InitType &Grid, aol::MultiVector<RealType> &Array ) const {
    const bool initialMeanPositions = Options.numSegments > 0 && Options.seedPositions.numComponents ( ) == Options.numSegments && Options.seedPositions[0].size ( ) == 2;
    aol::Vector<int> indices ( Options.numSegments );
    int offset = ( Options.includeBoundary ) ? 0 : ( Options.transformScale - 1 ) / 2;
    if ( initialMeanPositions ) {
      qc::FastILexMapper<qc::QC_2D> mapper ( Grid );
      int x, y;
      for ( int l=0; l<Options.numSegments ; ++l ) {
        x = Options.seedPositions[l][0] / Options.resampleFactor - offset;
        y = Options.seedPositions[l][1] / Options.resampleFactor - offset;
        if ( x < 0 || y < 0 ) throw aol::Exception ( "Seed position outside of support! If boundary is not included, choose seed positions away from boundary.", __FILE__, __LINE__ );
        indices[l] = mapper.getGlobalIndex ( x, y );
      }
    } else if ( Options.seedPositions.numComponents ( ) > 0 ) std::cerr << "WARNING: seed positions were specified but have wrong dimensions! Will use automatic seeding instead.." << std::endl;
    
    HardSegmentation.reallocate ( Grid );
    if ( Options.numSegments == 2 && !Options.forceMultiPhaseSegmentation ) {
      // Compute soft segmentation
      SegmentationType segmentation ( Grid );
      PCAMSSegmentor<ConfiguratorType, TwoPhaseSegmentorType> segmentor ( Grid, Options.gamma, Array, !initialMeanPositions, Options.numSegments,
                                                                          Options.numEvals, Options.omega, Options.orthonormalizeMeanValues, Options.numGhostCells );
      segmentor.setQuietMode ( this->_quietMode );
      segmentor.setCatchCtrlC ( _catchCtrlC );
      segmentor.setMaxIterations ( Options.maxIt );
      segmentor.setStopEpsilon ( Options.epsilon );
      segmentor.setOuterIterations ( Options.numOuterIterations );
      if ( initialMeanPositions ) segmentor.setMeanValues ( indices );
      segmentor.segmentAndAdjustGrayValues ( segmentation );
      
      if ( !_quietMode && _outputDir != "" ) {
        // Save soft segmentation
        std::stringstream ss;
        ss << _outputDir << "/softSegmentation" << qc::getDefaultArraySuffix ( qc::QC_2D );
        segmentation.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      }
      
      // Threshold to receive hard segmentation
      segmentation.threshold ( 0.5, 0, 1 );
      for ( int k=0; k<segmentation.size ( ) ; ++k ) HardSegmentation[k] = static_cast<int> ( segmentation[k] );
    } else {
      aol::VectorContainer<SegmentationType> initialSegmentation ( aol::Max<int> ( 1, Options.numSegments ), SegmentationType ( Grid ) );
      if ( !initialMeanPositions ) getInitialSegmentationFromLowEdgenessDominantFeatureClustering ( Options, initialSegmentation, Grid, Array );
      
      // Compute soft segmentation
      aol::VectorContainer<SegmentationType> segmentation ( Options.numSegments, SegmentationType ( Grid ) );
      PCAMSSegmentor<ConfiguratorType, MultiPhaseSegmentorType> segmentor ( Grid, Options.gamma, Array, false, Options.numSegments,
                                                                            Options.numEvals, Options.omega, Options.orthonormalizeMeanValues, Options.numGhostCells, _outputDir );
      segmentor.setQuietMode ( this->_quietMode );
      segmentor.setOutputDirectory ( _outputDir );
      segmentor.setCatchCtrlC ( _catchCtrlC );
      segmentor.setMaxIterations ( Options.maxIt );
      segmentor.setStopEpsilon ( Options.epsilon );
      segmentor.setOuterIterations ( Options.numOuterIterations );
      if ( initialMeanPositions ) segmentor.setMeanValues ( indices );
      else segmentor.setMeanValuesFromInitialSegmentation ( initialSegmentation );
      segmentor.segmentAndAdjustGrayValues ( segmentation );
      
      if ( !_quietMode && _outputDir != "" ) {
        // Save soft segmentations
        for ( int l=0; l<Options.numSegments ; ++l ) {
          std::stringstream ss;
          ss << _outputDir << "/softSegmentation_" << l << qc::getDefaultArraySuffix ( qc::QC_2D );
          segmentation[l].save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
        }
      }
      
      segmentor.getHardSegmentation ( HardSegmentation, segmentation );
    }
    
    // Extend segmentation to boundary (if necessary)
    if ( !Options.includeBoundary ) {
      qc::ScalarArray<int, qc::QC_2D> hardSegmentation ( HardSegmentation );
      HardSegmentation.reallocate ( Options.input.getNumX ( ) / Options.resampleFactor, Options.input.getNumY ( ) / Options.resampleFactor );
      for ( int y=0; y<HardSegmentation.getNumY ( ) ; ++y )
        for ( int x=0; x<HardSegmentation.getNumX ( ) ; ++x )
          HardSegmentation.set ( x, y, hardSegmentation.getClip ( x - offset, y - offset ) );
    }
  }
  
  void upsampleHardSegmentation ( qc::ScalarArray<int, qc::QC_2D> &HardSegmentation, const aol::Vec2<short> &Size ) const {
    qc::ScalarArray<int, qc::QC_2D> hardSegmentation ( HardSegmentation );
    HardSegmentation.reallocate ( Size[0], Size[1] );
    for ( int y=0; y<HardSegmentation.getNumY ( ) ; ++y )
      for ( int x=0; x<HardSegmentation.getNumX ( ) ; ++x )
        HardSegmentation.set ( x, y, hardSegmentation.get ( round ( x * hardSegmentation.getNumX ( ) / Size[0] ), round ( y * hardSegmentation.getNumY ( ) / Size[1] ) ) );
  }
  
  void getInitialSegmentationFromLowEdgenessDominantFeatureClustering ( OptionsType &Options, aol::VectorContainer<SegmentationType> &Segmentation, const InitType &Grid, aol::MultiVector<RealType> &Array ) const {
    // Reduce array to the components created by the dominant operators
    // Then perform PCA on it
    aol::MultiVector<RealType> reducedArray ( _numDominantComponents, Array[0].size ( ) );
    for ( int i=0; i<Array[0].size ( ) ; ++i )
      for ( int c=0; c<reducedArray.numComponents ( ); ++c )
        reducedArray[c][i] = Array[c][i];
    aol::MultiVector<RealType> eigenVecs, reducedPCAArray;
    aol::Vector<RealType> eigenVals;
    PCAMSSegmentor<ConfiguratorType, MultiPhaseSegmentorType>::getEigenVecs ( eigenVecs, eigenVals, reducedArray, Options.numEvals, Options.omega );
    PCAMSSegmentor<ConfiguratorType, MultiPhaseSegmentorType>::getProjectedCoefficients ( reducedPCAArray, reducedArray, eigenVecs );
    if ( Options.numSegments == 0 ) {
      Options.numSegments = eigenVecs.numComponents ( );
      if ( !this->_quietMode ) std::cerr << "PCA detected " << Options.numSegments << " segments." << std::endl;
    }
    
    if ( !_quietMode && _outputDir != "" ) {
      eigenVals.saveASCII ( aol::strprintf ( "%s/eigenValues_full.txt", this->_outputDir.c_str ( ) ).c_str ( ) );
      std::ofstream txtFile ( aol::strprintf ( "%s/numSamples_full.txt", this->_outputDir.c_str ( ) ).c_str ( ) );
      txtFile << reducedArray[0].size ( ) << std::endl;
      txtFile.close ( );
    }
    
    // Determine edge-ness indicators
    std::cerr << Grid.getNumberOfNodes ( ) << std::endl;
    std::cerr << reducedPCAArray[0].size ( ) << std::endl;
    if ( Grid.getNumberOfNodes ( ) != reducedPCAArray[0].size ( ) ) throw aol::Exception ( "Grid dimensions and PCA array size do not match!", __FILE__, __LINE__ );
    qc::FastILexMapper<qc::QC_2D> mapper ( Grid );
    int x, y, h = ( Options.transformScale - 1 ) / 2;
    int i1, i2;
    aol::Vector<RealType> f1 ( reducedPCAArray.numComponents ( ) ), f2 ( f1, aol::STRUCT_COPY );
    aol::Vector<RealType> indicators ( reducedPCAArray[0].size ( ) );
    RealType meanIndicator = 0;
    int numValidIndicators = 0;
    for ( int i=0; i<reducedPCAArray[0].size ( ) ; ++i ) {
      mapper.splitGlobalIndex ( i, x, y );
      if ( x - h >= 0 && x + h < Grid.getNumX ( ) && y - h >= 0 && y + h < Grid.getNumY ( ) ) {
        i1 = mapper.getGlobalIndex ( x - h, y );
        i2 = mapper.getGlobalIndex ( x + h, y );
        reducedPCAArray.getTo ( i1, f1 );
        reducedPCAArray.getTo ( i2, f2 );
        f1 -= f2;
        indicators[i] += f1.normSqr ( );
        i1 = mapper.getGlobalIndex ( x, y - h );
        i2 = mapper.getGlobalIndex ( x, y + h );
        reducedPCAArray.getTo ( i1, f1 );
        reducedPCAArray.getTo ( i2, f2 );
        f1 -= f2;
        indicators[i] += f1.normSqr ( );
        
        meanIndicator += indicators[i];
        ++numValidIndicators;
      } else indicators[i] = -1;
    }
    meanIndicator /= static_cast<RealType> ( numValidIndicators );
    
    // Construct array with low edge-ness
    aol::MultiVector<RealType> reducedPCAArrayLowEdgeness ( reducedPCAArray.numComponents ( ), reducedPCAArray[0].size ( ) );
    aol::Vector<int> indexCorrespondences ( Array[0].size ( ) );
    int j = 0;
    for ( int i = 0; i<reducedPCAArray[0].size ( ); ++i ) {
      if ( indicators[i] > 0 && indicators[i] < Options.lowEdgenessThreshold * meanIndicator ) {
        for ( int c = 0; c<reducedPCAArray.numComponents ( ); ++c )
          reducedPCAArrayLowEdgeness[c][j] = reducedPCAArray[c][i];
        indexCorrespondences[j] = i;
        ++j;
      }
    }
    reducedPCAArrayLowEdgeness.resize ( reducedPCAArray.numComponents ( ), j );
    
    // Perform clustering on reduced array with low edge-ness
    aol::KMeansClusterer<RealType> kMeansClusterer;
    aol::MultiVector<RealType> clusters;
    aol::Vector<int> clusterLabels;
    
    RealType residual = 0;
    if ( Options.forceGroundTruthClusters && Options.groundTruth.size ( ) > 0 ) {
      // If clustering result is bad, use this to check if it could be improved by a better initialization of k-means
      qc::ScalarArray<int, qc::QC_2D> gt ( Options.groundTruth, aol::DEEP_COPY );
      transformGrayValuesToLabels ( gt );
      aol::MultiVector<RealType> trueCenters ( reducedPCAArrayLowEdgeness.numComponents ( ), Options.numSegments );
      for ( int l=0; l<=gt.getMaxValue ( ) ; ++l ) {
        qc::BitArray<qc::QC_2D> bitMask ( gt.getNumX ( ), gt.getNumY ( ) );
        for ( int y=1; y<gt.getNumY ( )-1 ; ++y )
          for ( int x=1; x<gt.getNumX ( )-1 ; ++x )
            if ( gt.get ( x, y ) == l ) bitMask.set ( x, y, true );
        bitMask.invert ( );
        for ( int j=0; j<12 ; ++j ) bitMask.dilateByOne ( );
        bitMask.invert ( );
        
        int n = 0;
        for ( int i=0; i<bitMask.size ( ) ; ++i ) {
          if ( indicators[i] > 0 && indicators[i] < Options.lowEdgenessThreshold * meanIndicator && bitMask[i] ) {
            for ( int c=0; c<reducedPCAArray.numComponents ( ) ; ++c )
              trueCenters[c][l] += reducedPCAArray[c][i];
            ++n;
          }
        }
        for ( int c=0; c<reducedPCAArray.numComponents ( ) ; ++c )
          trueCenters[c][l] /= static_cast<RealType> ( n );
      }
      kMeansClusterer.setInitialCenters ( trueCenters );
      residual = kMeansClusterer.apply ( reducedPCAArrayLowEdgeness, clusters, Options.numSegments, clusterLabels, aol::KMEANS_INIT_METHOD::MAN );
    } else
      residual = kMeansClusterer.applyMultipleRNG ( reducedPCAArrayLowEdgeness, clusters, Options.numSegments, clusterLabels, 1000 );
    if ( !_quietMode ) {
      std::cerr << "Residual after clustering: " << residual << std::endl;
      std::cerr << clusters << std::endl;
    }
    
    // Create segmentation labels based on clustering
    for ( int l=Segmentation.size ( ) ; l<Options.numSegments ; ++l ) Segmentation.pushBack ( Segmentation[0] ); // add required number of segments if automatically detected
    for ( int j=0; j<clusterLabels.size ( ) ; ++j ) Segmentation[clusterLabels[j]][indexCorrespondences[j]] = 1.0;
    
    if ( !_quietMode && _outputDir != "" ) {
      // Save initial segmentations
      for ( int l=0; l<Options.numSegments ; ++l ) {
        std::stringstream ss;
        ss << _outputDir << "/initialSegmentation_" << l << qc::getDefaultArraySuffix ( qc::QC_2D );
        Segmentation[l].save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      }
    }
  }
  
  void setCtrlCHandler () const {
    if (_catchCtrlC)
      _previousCtrlCHandler = signal ( InterruptSignal, aol::ctrlCHandler );
  }
  
  void unsetCtrlCHandler () const {
    if (_catchCtrlC)
      signal ( InterruptSignal, _previousCtrlCHandler );
  }
  
  bool wantsInterrupt() const {
    if (!_catchCtrlC || !aol::getCtrlCState())
      return false;
    else
      return true;
  }
};

template<typename _RealType, typename _PictureType>
class PiecewisePeriodicSegmentationOptions : public LocalFeatureMSSegmentationOptions<_RealType, _PictureType> {
public:
  PiecewisePeriodicSegmentationOptions ( const _PictureType &Input, const int BlockSize, const std::string OutputDir = "", const bool Quiet = true )
    : LocalFeatureMSSegmentationOptions<_RealType, _PictureType> ( Input, OutputDir, Quiet ) {
    this->operators.push_back ( new FFTModulusOp<_RealType, _PictureType> ( BlockSize ) );
    this->includeBoundary = false;
    this->operatorWeights.pushBack ( 1.0 );
    this->numDominantOperators = 1;
    this->transformScale = BlockSize;
  }
};

#endif
