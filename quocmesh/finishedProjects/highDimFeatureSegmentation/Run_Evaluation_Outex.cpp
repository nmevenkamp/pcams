#include <aol.h>
#include <segmentation.h>

typedef double RType;
typedef qc::ScalarArray<int, qc::QC_2D> PicType;

typedef RType RealType;
typedef PicType PictureType;


int main( int argc, char **argv ) {
  try {
#ifndef USE_LIB_HUNGARIAN
    throw aol::Exception ( "This executable requires the Hungarian external library. Please refer to the README.txt located inside the same source code folder as this executable for further instructions.", __FILE__, __LINE__ );
#endif
    if ( argc < 3 ) std::cerr << "Usage: Run_Evaluation_Outex <segmentationDirectory> <pathToGroundTruthSegments>"<< std::endl;
    else {
      const std::string segDir = argv[1];
      const std::string gtPath = argv[2];
      
      const int N = 100;
      aol::Vector<RealType> o, c, ca, co, cc, cs;
      aol::ProgressBar<> progressBar ( "Analyzing segmentation error" );
      progressBar.start ( N );
      for ( int i=0; i<N ; ++i, progressBar++ ) {
        PictureType seg ( aol::strprintf ( "%s/seg%03d.png", segDir.c_str ( ), i ).c_str ( ) );
        PictureType gt ( gtPath.c_str ( ) );
        
        SegmentationQualityQuantifier<RealType> segQualityQuantifier ( gt, seg );
        o.pushBack ( segQualityQuantifier.getOmissionError ( ) );
        c.pushBack ( segQualityQuantifier.getCommissionError ( ) );
        ca.pushBack ( segQualityQuantifier.getClassAccuracy ( ) );
        co.pushBack ( segQualityQuantifier.getCorrectAssignment ( ) );
        cc.pushBack ( segQualityQuantifier.getObjectAccuracy ( ) );
        cs.pushBack ( segQualityQuantifier.getCorrectDetection ( ) );
      }
      progressBar.finish ( );
      
      std::ofstream txtFile ( aol::strprintf ( "%s/quantitativeEvaluation.txt", segDir.c_str ( ) ).c_str ( ) );
      txtFile << "O  = " << o.getMeanValue ( ) << " +- " << o.getStdDev ( ) << std::endl;
      txtFile << "C  = " << c.getMeanValue ( ) << " +- " << c.getStdDev ( ) << std::endl;
      txtFile << "CA = " << ca.getMeanValue ( ) << " +- " << ca.getStdDev ( ) << std::endl;
      txtFile << "CO = " << co.getMeanValue ( ) << " +- " << co.getStdDev ( ) << std::endl;
      txtFile << "CC = " << cc.getMeanValue ( ) << " +- " << cc.getStdDev ( ) << std::endl;
      txtFile << "CS = " << cs.getMeanValue ( ) << " +- " << cs.getStdDev ( ) << std::endl;
    }
  } catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
