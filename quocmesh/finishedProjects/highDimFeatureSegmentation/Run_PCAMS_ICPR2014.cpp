#include <aol.h>
#include <quoc.h>
#include <segmentation.h>


typedef double RealType;
typedef qc::ScalarArray<RealType, qc::QC_2D> ScalarPictureType;
typedef qc::MultiArray<RealType, qc::QC_2D, 3> PictureType;


int main ( int argc, char** argv ) {
  try {
    if ( argc < 3 ) std::cerr << "Usage: Run_PCAMS_PragueICPR2014 <sourceDirectory> <outputDirectory> [verbose]" << std::endl;
    else {
      const std::string srcDir = argv[1];
      const std::string outputDir = argv[2];
      aol::makeDirectory ( outputDir.c_str ( ) );
      
      const int N = 80;
      for ( int n=0; n<N ; ++n ) {
        std::cerr << "n = " << n+1 << std::endl;
        
        // Load texture mosaic and convert to CIELab color space
        const int n1 = std::floor ( n / 4 ) + 1, n2 = std::floor ( ( n % 4 ) / 2 ) + 1, n3 = ( n % 4 ) % 2 + 1;
        const std::string srcPath = aol::strprintf ( "%s/tm%d_%d_%d.png", srcDir.c_str ( ), n1, n2, n3 );
        PictureType u ( srcPath.c_str ( ) );
        RealType rgb[3], CIELab[3];
        for ( int j=0; j<u[0].size ( ) ; ++j ) {
          for ( int c=0; c<3 ; ++c ) rgb[c] = u[c][j];
          aol::RGBColorMap<RealType>::rgb2CIELab ( rgb, CIELab );
          for ( int c=0; c<3 ; ++c ) u[c][j] = CIELab[c];
        }
        u.scaleValuesTo01 ( );
        
        // Include possibility to produce verbose output (initial segmentation after clustering & segmentations after each outer iteration)
        std::string verboseOutputDir = ( argc == 4 && ( std::strcmp ( argv[3], "verbose" ) == 0 ) ) ? aol::strprintf ( "%s/tm%d_%d_%d", outputDir.c_str ( ), n1, n2, n3 ) : "";
        if ( verboseOutputDir != "" ) aol::makeDirectory ( verboseOutputDir.c_str ( ) );
      
        // Setup general options
        LocalFeatureMSSegmentationOptions<RealType, PictureType> options ( u, verboseOutputDir, false );
        options.numGhostCells = 15;
        options.omega = 0.05;
        options.gamma = 0.01;
        options.maxIt = 10000;
        options.epsilon = 0.001;
        options.numOuterIterations = 3;
        options.lowEdgenessThreshold = 0.5;
        options.forceMultiPhaseSegmentation = true;
        
        // Setup local spectral histograms
        const int numBins = 11;
        aol::Vector<int> integrationScales;
        integrationScales.pushBack ( 61 );
        integrationScales.pushBack ( 31 );
        options.operatorWeights.pushBack ( 0.2 );
        options.operatorWeights.pushBack ( 0.8 );
        options.numDominantOperators = 2;
        options.transformScale = integrationScales.getMaxValue ( );
        
        for ( int i=0; i<integrationScales.size ( ) ; ++i ) {
          LocalSpectralHistogramOp<RealType, PictureType> *spectralHistOp = new LocalSpectralHistogramOp<RealType, PictureType> ( );
          spectralHistOp->pushBack ( new IntensityFilterOp<ScalarPictureType> ( ), integrationScales[i], numBins, false );
          spectralHistOp->pushBack ( new GaborFilterOp<ScalarPictureType> ( 5, 0.5 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<ScalarPictureType> ( 5, 0.0 ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<ScalarPictureType> ( 5, 0.25 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<ScalarPictureType> ( 5, -0.25 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<ScalarPictureType> ( 7, 0.5 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<ScalarPictureType> ( 7, 0.0 ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<ScalarPictureType> ( 7, 0.25 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<ScalarPictureType> ( 7, -0.25 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          options.operators.push_back ( spectralHistOp );
        }
        
        // Apply PCA-based Mumford-Shah segmentation
        qc::ScalarArray<int, qc::QC_2D> segmentation;
        LocalFeatureMSSegmentor<RealType, PictureType> segmentor;
        segmentor.apply ( options, segmentation );
        segmentation.addToAll ( 1 );
        segmentation.savePNG ( aol::strprintf ( "%s/seg%d_%d_%d.png", outputDir.c_str ( ), n1, n2, n3 ).c_str ( ) );
      }
    }
  } catch ( aol::Exception& ex ) {
    ex.dump();
  }
}
