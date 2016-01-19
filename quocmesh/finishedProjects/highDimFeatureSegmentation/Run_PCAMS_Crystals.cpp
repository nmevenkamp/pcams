#include <aol.h>
#include <quoc.h>
#include <segmentation.h>
#include <randomGenerator.h>


typedef double RealType;
typedef qc::ScalarArray<RealType, qc::QC_2D> PictureType;
typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType,qc::QC_2D,3> > ConfiguratorType;


int main ( int argc, char** argv ) {
  try {
    if ( argc < 3 ) std::cerr << "Usage: Run_PCAMS_PragueICPR2014 <sourceDirectory> <outputDirectory> [verbose]" << std::endl;
    else {
      const std::string srcDir = argv[1];
      const std::string outputDir = argv[2];
      aol::makeDirectory ( outputDir.c_str ( ) );
      
      const int N = 4;
      const std::string fileNames[N] = { "grain", "grains3", "grains5", "CMS-GaAs" };
      const int numSegments[N] = { 2, 3, 5, 2 };
      const int blockSizes[N] = { 31, 31, 31, 41 };
      const RealType noiseStdDevs[N] = { 255.0, 168.3, 255.0, 255.0 };
      
      for ( int n=0; n<N ; ++n ) {
        for ( int addNoise=0; addNoise<=1 ; ++addNoise ) {
          std::cerr << "n = " << n+1;
          
          // Load crystal image
          const std::string srcPath = aol::strprintf ( "%s/%s.png", srcDir.c_str ( ), fileNames[n].c_str ( ) );
          PictureType u ( srcPath.c_str ( ) );

          if ( addNoise ) {
            std::cerr << " (added AGWN with sigma = " << noiseStdDevs[n] / u.getMaxValue ( ) * 100.0 << "% of maximum intensity)";
            aol::RandomGenerator rng;
            for ( int k=0; k<u.size ( ) ; ++k )
              u[k] += rng.normalrReal<RealType> ( 0.0, noiseStdDevs[n] );
          }
          
          std::cerr << std::endl;

          // Include possibility to produce verbose output (initial segmentation after clustering & segmentations after each outer iteration)
          std::string verboseOutputDir = ( argc == 4 && ( std::strcmp ( argv[3], "verbose" ) == 0 ) ) ? aol::strprintf ( "%s/%s%s", outputDir.c_str ( ), fileNames[n].c_str ( ), ( addNoise ) ? "_noisy" : "" ) : "";
          if ( verboseOutputDir != "" ) aol::makeDirectory ( verboseOutputDir.c_str ( ) );

          // Setup options
          PiecewisePeriodicSegmentationOptions<RealType, PictureType> options ( u, blockSizes[n], verboseOutputDir, false );
          options.numSegments = numSegments[n];
          options.numEvals = numSegments[n];
          options.gamma = 25;
          options.maxIt = 10000;
          options.epsilon = 0.001;
          options.numOuterIterations = 3;
          options.lowEdgenessThreshold = 1.0;
          options.forceMultiPhaseSegmentation = true;

          // Apply PCA-based Mumford-Shah segmentation
          qc::ScalarArray<int, qc::QC_2D> segmentation;
          LocalFeatureMSSegmentor<RealType, PictureType> segmentor;
          segmentor.apply ( options, segmentation );
          segmentation.addToAll ( 1 );
          segmentation.savePNG ( aol::strprintf ( "%s/seg_%s%s.png", outputDir.c_str ( ), fileNames[n].c_str ( ), ( addNoise ) ? "_noisy" : "" ).c_str ( ) );
          plotSegmentationBoundariesOntoImage<RealType, PictureType> ( aol::strprintf ( "%s/bdry_%s%s", outputDir.c_str ( ), fileNames[n].c_str ( ), ( addNoise ) ? "_noisy" : "" ).c_str ( ), u, segmentation );
        }
      }
    }
  } catch ( aol::Exception& ex ) {
    ex.dump();
  }
}
