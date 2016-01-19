#include <aol.h>
#include <quoc.h>
#include <segmentation.h>


typedef double RealType;
typedef qc::ScalarArray<RealType, qc::QC_2D> PictureType;


int main ( int argc, char** argv ) {
  try {
    if ( argc < 3 ) std::cerr << "Usage: Run_PCAMS_Outex <sourceDirectory> <outputDirectory> [verbose]"<< std::endl;
    else {
      const std::string srcDir = argv[1];
      const std::string outputDir = argv[2];
      aol::makeDirectory ( outputDir.c_str ( ) );

      const int N = 100;
      for ( int n=0; n<N ; ++n ) {
        std::cerr << "n = " << n+1 << std::endl;
        
        // Load texture mosaic
        const std::string srcPath = aol::strprintf ( "%s/%.3d.png", srcDir.c_str ( ), n );
        PictureType u ( srcPath.c_str ( ) );
        u.scaleValuesTo01 ( );

        // Include possibility to produce verbose output (initial segmentation after clustering & segmentations after each outer iteration)
        std::string verboseOutputDir = ( argc == 4 && ( std::strcmp ( argv[3], "verbose" ) == 0 ) ) ? aol::strprintf ( "%s/%.3d", outputDir.c_str ( ), n ) : "";
        if ( verboseOutputDir != "" ) aol::makeDirectory ( verboseOutputDir.c_str ( ) );
        
        // Setup general options
        LocalFeatureMSSegmentationOptions<RealType, PictureType> options ( u, verboseOutputDir, false );
        options.numSegments = 5;
        options.numEvals = 5;
        options.numGhostCells = 30;
        options.gamma = 0.005;
        options.maxIt = 10000;
        options.epsilon = 0.001;
        options.numOuterIterations = 3;
        options.lowEdgenessThreshold = 0.25;
        options.forceMultiPhaseSegmentation = true;
      
        // Setup local spectral histograms
        const int numBins = 11;
        aol::Vector<int> integrationScales;
        integrationScales.pushBack ( 61 );
        integrationScales.pushBack ( 31 );
        options.operatorWeights.pushBack ( 0.2 );
        options.operatorWeights.pushBack ( 0.8 );
        options.numDominantOperators = 1;
        options.transformScale = integrationScales.getMaxValue ( );
        
        for ( int i=0; i<integrationScales.size ( ) ; ++i ) {
          LocalSpectralHistogramOp<RealType, PictureType> *spectralHistOp = new LocalSpectralHistogramOp<RealType, PictureType> ( );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 5, 0.5 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 5, 0.0 ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 5, 0.25 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 5, -0.25 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 7, 0.5 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 7, 0.0 ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 7, 0.25 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 7, -0.25 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 9, 0.5 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 9, 0.0 ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 9, 0.25 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 9, -0.25 * aol::NumberTrait<RealType>::pi ), integrationScales[i], numBins );
          options.operators.push_back ( spectralHistOp );
        }
      
        // Apply PCA-based Mumford-Shah segmentation
        qc::ScalarArray<int, qc::QC_2D> segmentation;
        LocalFeatureMSSegmentor<RealType, PictureType> segmentor;
        segmentor.apply ( options, segmentation );
        segmentation.addToAll ( 1 );
        segmentation.savePNG ( aol::strprintf ( "%s/seg%.3d.png", outputDir.c_str ( ), n ).c_str ( ) );
      }
    }
  } catch ( aol::Exception& ex ) {
    ex.dump();
  }
}
