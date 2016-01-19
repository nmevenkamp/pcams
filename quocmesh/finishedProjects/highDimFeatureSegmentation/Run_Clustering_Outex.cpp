#include <aol.h>
#include <quoc.h>
#include <segmentation.h>


typedef double RealType;
typedef qc::ScalarArray<RealType, qc::QC_2D> PictureType;


int main ( int argc, char** argv ) {
  try {
    if ( argc != 3 ) std::cerr << "Usage: Run_Clustering_Outex <sourceDirectory> <outputDirectory>"<< std::endl;
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
        
        // Setup local spectral histograms
        LocalFeatureMSSegmentationOptions<RealType, PictureType> options ( u, "", false );
        const int numBins = 11;
        options.transformScale = 61;
        options.operatorWeights.pushBack ( 1.0 );
        LocalSpectralHistogramOp<RealType, PictureType> *spectralHistOp = new LocalSpectralHistogramOp<RealType, PictureType> ( );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 5, 0.5 * aol::NumberTrait<RealType>::pi ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 5, 0.0 ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 5, 0.25 * aol::NumberTrait<RealType>::pi ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 5, -0.25 * aol::NumberTrait<RealType>::pi ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 7, 0.5 * aol::NumberTrait<RealType>::pi ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 7, 0.0 ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 7, 0.25 * aol::NumberTrait<RealType>::pi ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 7, -0.25 * aol::NumberTrait<RealType>::pi ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 9, 0.5 * aol::NumberTrait<RealType>::pi ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 9, 0.0 ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 9, 0.25 * aol::NumberTrait<RealType>::pi ), options.transformScale, numBins );
        spectralHistOp->pushBack ( new GaborFilterOp<PictureType> ( 9, -0.25 * aol::NumberTrait<RealType>::pi ), options.transformScale, numBins );
        options.operators.push_back ( spectralHistOp );
        
        // Extract features using local spectral histograms
        aol::MultiVector<RealType> features;
        LocalFeatureMSSegmentor<RealType, PictureType> segmentor;
        segmentor.getLocalFeatures ( options, features );
        
        // Apply K-means clustering
        const int numClusters = 5;
        aol::KMeansClusterer<RealType> kMeansClusterer;
        aol::MultiVector<RealType> clusters;
        aol::Vector<int> clusterLabels;
        kMeansClusterer.apply ( features, clusters, numClusters, clusterLabels );
        
        // Save cluster labels as segmentation
        qc::ScalarArray<int, qc::QC_2D> segmentation ( u.getNumX ( ), u.getNumY ( ) );
        for ( int j=0; j<clusterLabels.size ( ) ; ++j ) segmentation[j] = clusterLabels[j]+1;
        segmentation.savePNG ( aol::strprintf ( "%s/seg%.3d.png", outputDir.c_str ( ), n ).c_str ( ) );
      }
    }
  } catch ( aol::Exception& ex ) {
    ex.dump();
  }
}
