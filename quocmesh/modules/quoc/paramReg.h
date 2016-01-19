#ifndef __PARAMREG_H
#define __PARAMREG_H

#include <registration.h>
#include <convolution.h>
#include <mutualInformation.h>
#include <tensor.h>

namespace qc {

/** 
 * \author Berkels
 */
template <bool AllowScaling>
class ParametricRigidBodyMotion2DHelper {};

template <>
class ParametricRigidBodyMotion2DHelper<false> {
public:
  static const int NumberOfScalingParameters = 0;
  static aol::Vec<2, int> getDeformParametersSize ( ) {
    return aol::Vec2<int> ( 2, 1 );
  }
};

template <>
class ParametricRigidBodyMotion2DHelper<true> {
public:
  static const int NumberOfScalingParameters = 1;
  static aol::Vec<3, int> getDeformParametersSize ( ) {
    return aol::Vec3<int> ( 2, 1, 1 );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class ParametricDeformationBase {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
  const typename ConfiguratorType::InitType &_grid;
public:
  ParametricDeformationBase ( const typename ConfiguratorType::InitType &Initializer )
    : _grid ( Initializer ) { }

  //! \note Inefficient implementation, don't use for anything where performance is critical.
  template <typename ParametricDeformationType>
  void evaluateDeformationOn01 ( const ParametricDeformationType &ParDef, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    typename ConfiguratorType::ElementType el;
    typename ConfiguratorType::DomVecType localCoord;
    qc::getLocalCoordsRegularRectangularGrid<ConfiguratorType> ( Position, _grid, el, localCoord );
    typename ConfiguratorType::ElementType transformedEl;
    typename ConfiguratorType::VecType transformedLocalCoord;
    ParDef.template evaluateDeformation<false> ( el, localCoord, transformedEl, transformedLocalCoord );

    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      DeformedPosition[i] = ( transformedEl[i] + transformedLocalCoord[i] ) * _grid.H();
  }

  const typename ConfiguratorType::InitType &getGridRef () const {
    return _grid;
  }
};

/**
 * \note With AllowScaling == true the name is misleading: The class is a rigid body motion plus scaling in this case.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, bool AllowScaling = false>
class ParametricRigidBodyMotion2D {
public:
  static const int NumberOfDeformParameters = 1 + ConfiguratorType::Dim + ParametricRigidBodyMotion2DHelper<AllowScaling>::NumberOfScalingParameters;
  static const int NumOfDeformParametersComponents = 2 + ParametricRigidBodyMotion2DHelper<AllowScaling>::NumberOfScalingParameters;
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const RealType _rotationAngle;
  aol::Matrix22<RealType> _rotationMatrix;
  typename ConfiguratorType::VecType _translation;
  RealType _scaling;
public:
  ParametricRigidBodyMotion2D ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
    : _grid ( Initializer ),
      _rotationAngle ( DeformParameters[1][0] ),
      _scaling ( AllowScaling ? DeformParameters[2][0] : 1 ) {
    _rotationMatrix.makeRotation ( _rotationAngle );
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      _translation[i] = DeformParameters[0][i] + 0.5;
  }

  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
    if ( AllowScaling )
      DeformParameters[ConfiguratorType::Dim][0] = 1;
  }

  // DestDeform := ArgDeform \circ DestDeform
  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &ArgDeformPars, aol::MultiVector<RealType> &DestDeformPars ) {
    aol::Matrix22<RealType > _rotationMatrix;
    aol::Vec2<RealType> tmp ( ArgDeformPars[0][0], ArgDeformPars[0][1] );
    _rotationMatrix.makeRotation ( ArgDeformPars[1][0] );
    _rotationMatrix.multAdd ( DestDeformPars[0], tmp );
    if ( AllowScaling ) {
      tmp *= ArgDeformPars[2][0];
      DestDeformPars[2][0] *= ArgDeformPars[2][0];
    }

    for ( int i = 0; i < 2; ++i )
      DestDeformPars[0][i] = tmp[i];

    DestDeformPars[1][0] += ArgDeformPars[1][0];
  }

  static aol::Vec<NumOfDeformParametersComponents, int> getDeformParametersSize ( ) {
    return ParametricRigidBodyMotion2DHelper<AllowScaling>::getDeformParametersSize();
  }

  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {

    typename ConfiguratorType::VecType coord, transformedCoord;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      coord[i] = El[i] + RefCoord[i] - 0.5 / _grid.H();
    }
    _rotationMatrix.mult ( coord, transformedCoord );
    if ( AllowScaling )
      transformedCoord *= _scaling;
    return qc::transformCoord<ConfiguratorType, ClipCoord> ( _grid, transformedCoord, _translation, TransformedEl, TransformedLocalCoord );
  }

  void evaluateDeformationOn01 ( const ParametricRigidBodyMotion2D<ConfiguratorType, AllowScaling> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    typename ConfiguratorType::VecType tmp ( Position );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      tmp[i] -= 0.5;
    _rotationMatrix.mult ( tmp, DeformedPosition );
    if ( AllowScaling )
      DeformedPosition *= _scaling;
    DeformedPosition += _translation;
  }

  void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &Position,
                                           aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    const RealType cosAlpha = _rotationMatrix[0][0];
    const RealType sinAlpha = _rotationMatrix[1][0];
    Jacobian[ConfiguratorType::Dim][0] = ( AllowScaling ? _scaling : 1 ) * ( -sinAlpha * ( Position[0] - 0.5 ) - cosAlpha * ( Position[1] - 0.5 ) );
    Jacobian[ConfiguratorType::Dim][1] = ( AllowScaling ? _scaling : 1 ) * (  cosAlpha * ( Position[0] - 0.5 ) - sinAlpha * ( Position[1] - 0.5 ) );

    Jacobian[0][0] = Jacobian[1][1] = 1;
    Jacobian[1][0] = Jacobian[0][1] = 0;

    if ( AllowScaling ) {
      Jacobian[ConfiguratorType::Dim+1][0] = ( cosAlpha * ( Position[0] - 0.5 ) - sinAlpha * ( Position[1] - 0.5 ) );
      Jacobian[ConfiguratorType::Dim+1][1] = ( sinAlpha * ( Position[0] - 0.5 ) + cosAlpha * ( Position[1] - 0.5 ) );
    }
  }

  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &El,
                                       const typename ConfiguratorType::VecType &RefCoord,
                                       aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    typename ConfiguratorType::VecType x;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      x[i] = ( El[i] + RefCoord[i] ) * _grid.H();
    }
    evaluateDerivativeDeformationOn01 ( x, Jacobian );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class ParametricRigidBodyMotion3D : public ParametricDeformationBase<ConfiguratorType> {
public:
  static const int NumberOfDeformParameters = 6;
  static const int NumOfDeformParametersComponents = 2;
  typedef typename ConfiguratorType::RealType RealType;
private:
  const RealType _yaw;
  const RealType _pitch;
  const RealType _roll;
  aol::Matrix33<RealType> _rotationYaw;
  aol::Matrix33<RealType> _rotationPitch;
  aol::Matrix33<RealType> _rotationRoll;
  typename ConfiguratorType::VecType _translation;
public:
  ParametricRigidBodyMotion3D ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ),
      _yaw ( DeformParameters[1][0] ),
      _pitch ( DeformParameters[1][1] ),
      _roll ( DeformParameters[1][2] ) {
    _rotationYaw.setRotationAboutZ ( _yaw );
    _rotationPitch.setRotationAboutY ( _pitch );
    _rotationRoll.setRotationAboutX ( _roll );
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      _translation[i] = DeformParameters[0][i] + 0.5;
  }

  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
  }

  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &/*ArgDeformPars*/, aol::MultiVector<RealType> &/*DestDeformPars*/ ) {
    throw aol::UnimplementedCodeException ( "not implemented", __FILE__, __LINE__ );
  }

  static aol::Vec<NumOfDeformParametersComponents, int> getDeformParametersSize ( ) {
    return aol::Vec2<int> ( 3, 3 );
  }

  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {

    typename ConfiguratorType::VecType coord, transformedCoord;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      coord[i] = El[i] + RefCoord[i] - 0.5 / this->_grid.H();
    }
    _rotationYaw.mult ( coord, transformedCoord );
    _rotationPitch.mult ( transformedCoord, coord );
    _rotationRoll.mult ( coord, transformedCoord );
    return qc::transformCoord<ConfiguratorType, ClipCoord> ( this->_grid, transformedCoord, _translation, TransformedEl, TransformedLocalCoord );
  }
  void evaluateDeformationOn01 ( const ParametricRigidBodyMotion3D<ConfiguratorType> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    typename ConfiguratorType::VecType tmp ( Position );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      tmp[i] -= 0.5;
    _rotationYaw.mult ( tmp, DeformedPosition );
    _rotationPitch.mult ( DeformedPosition, tmp );
    _rotationRoll.mult ( tmp, DeformedPosition );
    DeformedPosition += _translation;
  }
  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &El,
                                       const typename ConfiguratorType::VecType &RefCoord,
                                       aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    typename ConfiguratorType::VecType temp;
    typename ConfiguratorType::VecType x;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      x[i] = ( El[i] + RefCoord[i] ) * this->_grid.H() - 0.5;
    }

    // Yaw derivative
    {
      const RealType cosYaw = _rotationYaw[0][0];
      const RealType sinYaw = _rotationYaw[1][0];
      Jacobian[ConfiguratorType::Dim][0] = -sinYaw * x[0] - cosYaw * x[1];
      Jacobian[ConfiguratorType::Dim][1] =  cosYaw * x[0] - sinYaw * x[1];
      Jacobian[ConfiguratorType::Dim][2] =  0;
      _rotationPitch.mult ( Jacobian[ConfiguratorType::Dim], temp );
      _rotationRoll.mult ( temp, Jacobian[ConfiguratorType::Dim] );
    }

    // Pitch derivative
    {
      const RealType cosPitch = _rotationPitch[0][0];
      const RealType sinPitch = _rotationPitch[2][0];
      _rotationYaw.mult ( x, Jacobian[ConfiguratorType::Dim+1] );
      temp[0] = -sinPitch * x[0] - cosPitch * x[2];
      temp[1] = 0;
      temp[2] = cosPitch * x[0] - sinPitch * x[2];
      _rotationRoll.mult ( temp, Jacobian[ConfiguratorType::Dim+1] );
    }

    // Roll derivative
    {
      const RealType cosRoll = _rotationRoll[1][1];
      const RealType sinRoll = _rotationRoll[2][1];
      _rotationYaw.mult ( x, Jacobian[ConfiguratorType::Dim+2] );
      _rotationPitch.mult ( Jacobian[ConfiguratorType::Dim+2], temp );
      Jacobian[ConfiguratorType::Dim+2][0] = 0;
      Jacobian[ConfiguratorType::Dim+2][1] = -sinRoll * x[1] - cosRoll * x[2];
      Jacobian[ConfiguratorType::Dim+2][2] =  cosRoll * x[1] - sinRoll * x[2];
    }

    Jacobian[0][0] = Jacobian[1][1] = Jacobian[2][2] = 1;
    Jacobian[1][0] = Jacobian[0][1] = Jacobian[2][0] = Jacobian[0][2] = Jacobian[2][1] = Jacobian[1][2] = 0;
  }
};

/** 
 * \author Berkels
 */
template <typename ConfiguratorType>
class ParametricTranslation : public ParametricDeformationBase<ConfiguratorType> {
public:
  static const int NumberOfDeformParameters = ConfiguratorType::Dim;
  static const int NumOfDeformParametersComponents = 1;
  typedef typename ConfiguratorType::RealType RealType;
private:
  typename ConfiguratorType::VecType _translation;
public:
  ParametricTranslation ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ) {
    DeformParameters.copyTo ( _translation );
  }

  ParametricTranslation ( const typename ConfiguratorType::InitType &Initializer, const aol::Vec<ConfiguratorType::Dim, RealType> &Translation )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ) {
    _translation = Translation;
  }

  ParametricTranslation ( const typename ConfiguratorType::InitType &Initializer )
    : ParametricDeformationBase<ConfiguratorType> ( Initializer ) {
    _translation.setZero();
  }

  void setTranslation ( const aol::Vec<ConfiguratorType::Dim, RealType> &Translation ) {
    _translation = Translation;
  }

  void setTranslation ( const aol::Vector<RealType> &Translation ) {
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      _translation[i] = Translation[i];
  }

  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
  }

  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &ArgDeformPars, aol::MultiVector<RealType> &DestDeformPars ) {
    DestDeformPars += ArgDeformPars;
  }

  static aol::Vec<1, int> getDeformParametersSize ( ) {
    return aol::Vec<1, int> ( ConfiguratorType::Dim );
  }

  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                             const typename ConfiguratorType::VecType &RefCoord,
                             qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {
    return qc::transformCoord<ConfiguratorType, ClipCoord> ( this->_grid, El, RefCoord, _translation, TransformedEl, TransformedLocalCoord );
  }
  void evaluateDeformationOn01 ( const aol::Vec<ConfiguratorType::Dim, RealType> &Position, typename aol::Vec<ConfiguratorType::Dim, RealType> &DeformedPosition ) const {
    DeformedPosition = Position;
    DeformedPosition += _translation;
  }
  void evaluateDeformationOn01 ( const ParametricTranslation<ConfiguratorType> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    evaluateDeformationOn01 ( Position, DeformedPosition );
  }
  void evaluateDeformationOn0N ( typename aol::Vec<ConfiguratorType::Dim, RealType> &DeformedPosition ) const {
    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      DeformedPosition[i] += _translation[i] / this->_grid.H();
  }
  void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &/*Position*/,
                                           aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian.setIdentity();
  }
  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &/*El*/,
                                       const typename ConfiguratorType::VecType &/*RefCoord*/,
                                       aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian.setIdentity();
  }
  void evaluateDerivativeDeformationOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                           aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian.setZero();
    Jacobian.addToDiagonal ( aol::ZOTrait<RealType>::one / this->_grid.H() );
  }

  void evaluateSecondDerivativeDeformationOn0N ( const typename ConfiguratorType::VecType &/*Position*/,
                                                 aol::Tensor3rdOrder< NumberOfDeformParameters, NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Hessian ) const {
    Hessian.setZero();
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricDeformedPosition {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::VecType _pos01;
  typename ConfiguratorType::VecType _deformedPos;
  bool _transformPositionInDomain;
public:
  ParametricDeformedPosition ( const int X, const int Y,
                              const ParametricDeformationType &ParDef,
                              const aol::Vec3<int> &GridSize,
                              const RealType H )
  : _pos01 ( X*H, Y*H ),
  _transformPositionInDomain ( true ) {

    ParDef.evaluateDeformationOn01 ( ParDef, _pos01, _deformedPos );
    _deformedPos /= H;
    for ( int c = 0; c < ConfiguratorType::Dim; c++ ) {
      const RealType temp = aol::Clamp ( _deformedPos[c], aol::ZOTrait<RealType>::zero, static_cast<RealType> ( GridSize[c] - 1 ) );
      if ( temp != _deformedPos[c] )
        _transformPositionInDomain = false;

      _deformedPos[c] = temp;
    }
  }

  const typename ConfiguratorType::VecType &getUndefPos01 ( ) const {
    return _pos01;
  }

  const typename ConfiguratorType::VecType &getDefPos ( ) const {
    return _deformedPos;
  }

  bool isDefPosInDomain ( ) const {
    return _transformPositionInDomain;
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType, qc::Dimension Dim>
struct doFastParametricDeformImage {
static void apply ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                    const typename ConfiguratorType::InitType &Grid,
                    aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                    const ParametricDeformationType &ParDef );
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
struct doFastParametricDeformImage<ConfiguratorType, ParametricDeformationType, qc::QC_2D> {
  static void apply ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                      const typename ConfiguratorType::InitType &Grid,
                      aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                      const ParametricDeformationType &ParDef ) {
    typedef typename ConfiguratorType::RealType RealType;
    const typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );

    const aol::Vec3<int> gridSize = Grid.getSize();
    const RealType h = static_cast<RealType> ( Grid.H() );

    for ( int j = 0; j < gridSize[1]; ++j ) {
      for ( int i = 0; i < gridSize[0]; ++i ) {
        const int index = qc::ILexCombine2( i, j, Grid.getNumX() );
        qc::ParametricDeformedPosition<ConfiguratorType, ParametricDeformationType> defPos ( i, j, ParDef, gridSize, h );
        DeformedImage[index] = imageArray.interpolate ( defPos.getDefPos() );
      }
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
struct doFastParametricDeformImage<ConfiguratorType, ParametricDeformationType, qc::QC_3D> {
  static void apply ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                      const typename ConfiguratorType::InitType &Grid,
                      aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                      const ParametricDeformationType &ParDef ) {
    typedef typename ConfiguratorType::RealType RealType;
    const typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );

    const aol::Vec3<int> gridSize = Grid.getSize();
    const RealType h = static_cast<RealType> ( Grid.H() );

    for ( int k = 0; k < gridSize[2]; ++k ) {
      for ( int j = 0; j < gridSize[1]; ++j ) {
        for ( int i = 0; i < gridSize[0]; ++i ) {
          const int index = qc::ILexCombine3( i, j, k, Grid.getNumX(), Grid.getNumY() );
          const typename ConfiguratorType::VecType pos ( i*h, j*h, k*h );
          typename ConfiguratorType::VecType ds;
          ParDef.evaluateDeformationOn01 ( ParDef, pos, ds );
          for ( int c = 0; c < ConfiguratorType::Dim; c++ )
            ds[c] = aol::Clamp ( ds[c] / h, aol::ZOTrait<RealType>::zero, static_cast<RealType> ( gridSize[c] - 1 ) );
          DeformedImage[index] = imageArray.interpolate ( ds );
        }
      }
    }
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType, qc::Dimension Dim>
void FastParametricDeformImage ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                                 const typename ConfiguratorType::InitType &Grid,
                                 aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                                 const ParametricDeformationType &ParDef ) {
  doFastParametricDeformImage<ConfiguratorType, ParametricDeformationType, Dim>::apply ( Image, Grid, DeformedImage, ParDef );
}

/**
 * \f$ \frac{1}{2}\int ((T\circ\Phi)(x)-R(x))^2 dx \f$ where \f$\Phi\f$ is a parametric deformation
 * whose parameters are given in the constructor and T is given as argument in apply(Add).
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricSSDEnergyFromImage : public aol::FENonlinIntegrationScalarInterface<ConfiguratorType, ParametricSSDEnergyFromImage<ConfiguratorType, ParametricDeformationType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const ParametricDeformationType &_parDef;
public:

  ParametricSSDEnergyFromImage ( const typename ConfiguratorType::InitType &Grid, 
                                 const aol::Vector<RealType> &ImR,
                                 const ParametricDeformationType &ParDef )
  : aol::FENonlinIntegrationScalarInterface<ConfiguratorType, ParametricSSDEnergyFromImage<ConfiguratorType, ParametricDeformationType> > ( Grid ),
    _r( Grid, ImR ),
    _parDef ( ParDef ) {}

  RealType evaluateIntegrand ( const aol::DiscreteFunctionDefault<ConfiguratorType> &DiscFuncs,
                               const typename ConfiguratorType::ElementType &El,
                               int QuadPoint, const typename ConfiguratorType::DomVecType &RefCoord ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    if ( !_parDef.template evaluateDeformation<true> ( El, RefCoord, transformedEl, transformedLocalCoord) ) {
      return 0;
    }

    return 0.5 * aol::Sqr ( DiscFuncs.evaluate(transformedEl, transformedLocalCoord) - _r.evaluateAtQuadPoint(El, QuadPoint) );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class VariationOfParametricSSDEnergyWRTParamsFromImage
  : public aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                           VariationOfParametricSSDEnergyWRTParamsFromImage<ConfiguratorType, ParametricDeformationType>,
                                                           1, ParametricDeformationType::NumberOfDeformParameters > {
public:
  typedef typename ConfiguratorType::RealType RealType;
protected:
private:
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const ParametricDeformationType &_parDef;

public:
  VariationOfParametricSSDEnergyWRTParamsFromImage ( const typename ConfiguratorType::InitType &Initializer,
                                                     const aol::Vector<RealType> &ImR,
                                                     const ParametricDeformationType &ParDef )
    : aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                      VariationOfParametricSSDEnergyWRTParamsFromImage<ConfiguratorType, ParametricDeformationType>,
                                                      1, ParametricDeformationType::NumberOfDeformParameters > ( Initializer ),
      _r( Initializer, ImR ),
      _parDef ( ParDef ) {}

  void evaluateIntegrand ( const aol::auto_container<1, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs, // DiscrFuncs[0] is supposed to contain T.
                           const typename ConfiguratorType::ElementType &El,
                           int QuadPoint, const typename ConfiguratorType::VecType &RefCoord, aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> &Integrand ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    if ( !_parDef.template evaluateDeformation<true> ( El, RefCoord, transformedEl, transformedLocalCoord) ) {
      Integrand.setZero();
      return;
    }

    aol::Mat< ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType > jacobian;
    _parDef.evaluateDerivativeDeformation ( El, RefCoord, jacobian );

    typename ConfiguratorType::VecType gradT;
    DiscrFuncs[0].evaluateGradient ( transformedEl, transformedLocalCoord, gradT );

    jacobian.mult ( gradT, Integrand );
    Integrand *= ( DiscrFuncs[0].evaluate ( transformedEl, transformedLocalCoord ) - _r.evaluateAtQuadPoint ( El, QuadPoint ) );
  }
};

/**
 * \brief  \f$ \frac{1}{2}\int ((T\circ\Phi)(x)-R(x))^2 dx \f$ where \f$\Phi\f$ is a parametric deformation
 * whose parameters are given as argument in apply(Add).
 *
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricSSDEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const aol::Vector<RealType> _r;
  const aol::Vector<RealType> _t;
public:
  ParametricSSDEnergy ( const typename ConfiguratorType::InitType &Grid, 
                        const aol::Vector<RealType> &ImR,
                        const aol::Vector<RealType> &ImT )
  : _grid ( Grid ),
    _r( ImR ), 
    _t( ImT ) {}

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    const ParametricDeformationType parDef ( _grid, MArg );
    ParametricSSDEnergyFromImage<ConfiguratorType, ParametricDeformationType> E ( _grid, _r, parDef );
    E.applyAdd ( _t, Dest );
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    const ParametricDeformationType parDef ( _grid, MArg );

    aol::Vector<RealType> dest ( ParametricDeformationType::NumberOfDeformParameters );
    VariationOfParametricSSDEnergyWRTParamsFromImage<ConfiguratorType, ParametricDeformationType> E ( _grid, _r, parDef );
    E.applyAdd ( _t, dest );
    MDest.copySplitFrom ( dest );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricCrossCorrelationForce
  : public aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                           ParametricCrossCorrelationForce<ConfiguratorType, ParametricDeformationType>,
                                                           1, ParametricDeformationType::NumberOfDeformParameters > {
public:
  typedef typename ConfiguratorType::RealType RealType;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const ParametricDeformationType &_parDef;

  ParametricCrossCorrelationForce ( const typename ConfiguratorType::InitType &Initializer,
                                    const aol::Vector<RealType> &ImR,
                                    const ParametricDeformationType &ParDef )
    : aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                      ParametricCrossCorrelationForce<ConfiguratorType, ParametricDeformationType>,
                                                      1, ParametricDeformationType::NumberOfDeformParameters > ( Initializer ),
      _r( Initializer, ImR ),
      _parDef ( ParDef ) {}

  void evaluateIntegrand ( const aol::auto_container<1, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs, // DiscrFuncs[0] is supposed to contain T.
                           const typename ConfiguratorType::ElementType &El,
                           int QuadPoint, const typename ConfiguratorType::VecType &RefCoord, aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> &Integrand ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    if ( !_parDef.template evaluateDeformation<true> ( El, RefCoord, transformedEl, transformedLocalCoord) ) {
      Integrand.setZero();
      return;
    }

    aol::Mat< ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType > jacobian;
    _parDef.evaluateDerivativeDeformation ( El, RefCoord, jacobian );

    typename ConfiguratorType::VecType gradT;
    DiscrFuncs[0].evaluateGradient ( transformedEl, transformedLocalCoord, gradT );

    jacobian.mult ( gradT, Integrand );
    Integrand *= - _r.evaluateAtQuadPoint ( El, QuadPoint );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricNCCEnergy : public qc::NormalizedCrossCorrelationEnergy<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;

  ParametricNCCEnergy ( const typename ConfiguratorType::InitType &Grid, 
                        const aol::Vector<RealType> &ImR,
                        const aol::Vector<RealType> &ImT )
  : qc::NormalizedCrossCorrelationEnergy<ConfiguratorType> ( Grid, ImR, ImT ) {}

  // By overloading deformImage, we can use apply(Add) from the base class to calculate the energy.
  void deformImage ( const typename ConfiguratorType::ArrayType &Image, typename ConfiguratorType::ArrayType &DeformedImage , const aol::MultiVector<RealType> &Phi) const {
    ParametricDeformationType parDef ( this->_grid, Phi );
    qc::FastParametricDeformImage<ConfiguratorType, ParametricDeformationType, ConfiguratorType::Dim> ( Image, this->_grid, DeformedImage, parDef );
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    aol::Scalar<RealType> energy;
    this->apply ( MArg, energy );
    aol::Vector<RealType> temp ( this->_normalizedR );
    temp.addMultiple ( this->_lastNormalizedDeformedT, energy[0] );
    temp /= this->_varOfLastDeformedT;

    aol::Vector<RealType> destVec ( ParametricDeformationType::NumberOfDeformParameters );
    ParametricDeformationType parDef ( this->_grid, MArg );
    ParametricCrossCorrelationForce<ConfiguratorType, ParametricDeformationType> force ( this->_grid, temp, parDef );
    force.apply ( this->_t, destVec );
    MDest.copySplitFrom ( destVec );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class DeformationMatrix : protected aol::SparseMatrix<typename ConfiguratorType::RealType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const ParametricDeformationType _parDef;
public:
  DeformationMatrix ( const typename ConfiguratorType::InitType &Grid,
                     const aol::MultiVector<RealType> &DeformParams )
  : aol::SparseMatrix<RealType> ( Grid ),
  _grid ( Grid ),
  _parDef ( Grid, DeformParams ) {

    const aol::Vec3<int> gridSize = _grid.getSize();
    const RealType h = static_cast<RealType> ( _grid.H() );

    for ( int j = 0; j < gridSize[1]; ++j ) {
      for ( int i = 0; i < gridSize[0]; ++i ) {
        const int index = qc::ILexCombine2( i, j, _grid.getNumX() );
        qc::ParametricDeformedPosition<ConfiguratorType, ParametricDeformationType> defPos ( i, j, _parDef, gridSize, h );

        RealType X = defPos.getDefPos()[0];
        RealType Y = defPos.getDefPos()[1];
        int xL = static_cast<int> ( X ), yL = static_cast<int> ( Y );

        if ( X >= static_cast<RealType> ( gridSize[0] - 1 ) ) xL = gridSize[0] - 2;
        if ( Y >= static_cast<RealType> ( gridSize[1] - 1 ) ) yL = gridSize[1] - 2;
        if ( X < 0. ) xL = 0;
        if ( Y < 0. ) yL = 0;

        X -= xL;
        Y -= yL;

        const int xU = xL + 1, yU = yL + 1;

        this->set ( index, qc::ILexCombine2( xL, yL, _grid.getNumX() ), ( 1 - X ) * ( 1 - Y ) );
        this->set ( index, qc::ILexCombine2( xU, yL, _grid.getNumX() ), X * ( 1 - Y ) );
        this->set ( index, qc::ILexCombine2( xL, yU, _grid.getNumX() ), ( 1 - X ) * Y );
        this->set ( index, qc::ILexCombine2( xU, yU, _grid.getNumX() ), X * Y );
      }
    }
  }

  const aol::SparseMatrix<typename ConfiguratorType::RealType> &getSparseMatrixRef ( ) const {
    return *this;
  }

  using aol::SparseMatrix<RealType>::apply;
  using aol::SparseMatrix<RealType>::applyAddTransposed;
  using aol::SparseMatrix<RealType>::transposeTo;
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
class SmoothedImage {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  typename ConfiguratorType::ArrayType _image;
  typename ConfiguratorType::ArrayType _imageDx;
  typename ConfiguratorType::ArrayType _imageDy;
  const RealType _sigma;
  const int _kernelSize;
public:
  SmoothedImage ( const typename ConfiguratorType::InitType &Grid,
                  const aol::Vector<RealType> &Image,
                  const bool ComputeSmoothDerivative,
                  const RealType Sigma = 1.5,
                  const int KernelSize = 5 )
    : _grid ( Grid ),
      _image ( Grid ),
      _imageDx ( ),
      _imageDy ( ),
      _sigma ( Sigma ),
      _kernelSize ( KernelSize ) {
    const typename ConfiguratorType::ArrayType image ( Image, Grid, aol::FLAT_COPY );
    qc::GaussKernel2d<RealType> kernel ( _kernelSize, _sigma );
    image.applyLinearFilterTo ( kernel, _image );
    if ( ComputeSmoothDerivative ) {
      _imageDx.reallocate ( Grid );
      _imageDy.reallocate ( Grid );
      qc::GaussDiffKernel2d<RealType> kernelX ( _kernelSize, _sigma, qc::DIFF_X );
      qc::GaussDiffKernel2d<RealType> kernelY ( _kernelSize, _sigma, qc::DIFF_Y );
      image.applyLinearFilterTo ( kernelX, _imageDx );
      image.applyLinearFilterTo ( kernelY, _imageDy );
    }
  }

  const typename ConfiguratorType::ArrayType &getSmoothImRef ( ) const {
    return _image;
  }

  const typename ConfiguratorType::ArrayType &getSmoothDXImRef ( ) const {
    return _imageDx;
  }

  const typename ConfiguratorType::ArrayType &getSmoothDYImRef ( ) const {
    return _imageDy;
  }

  RealType interpolate ( const typename ConfiguratorType::VecType &Pos ) const {
    return _image.interpolate ( Pos );
  }

  RealType interpolateDX ( const typename ConfiguratorType::VecType &Pos ) const {
    return _imageDx.interpolate ( Pos );
  }

  RealType interpolateDY ( const typename ConfiguratorType::VecType &Pos ) const {
    return _imageDy.interpolate ( Pos );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class SimpleParametricSSDEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
public:
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const qc::SmoothedImage<ConfiguratorType> _r;
  const qc::SmoothedImage<ConfiguratorType> _t;
public:
  SimpleParametricSSDEnergy ( const typename ConfiguratorType::InitType &Grid,
                              const aol::Vector<RealType> &ImR,
                              const aol::Vector<RealType> &ImT )
    : _grid ( Grid ),
      _r ( Grid, ImR, false ),
      _t ( Grid, ImT, true ) { }

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    const ParametricDeformationType parDef ( _grid, MArg );
    typename ConfiguratorType::ArrayType deformedT ( _grid );
    qc::FastParametricDeformImage<ConfiguratorType, ParametricDeformationType, ConfiguratorType::Dim>
      ( _t.getSmoothImRef(), _grid, deformedT, parDef );

    deformedT -= _r.getSmoothImRef();
    Dest[0] += 0.5*deformedT.normSqr();
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    const ParametricDeformationType parDef ( _grid, MArg );

    aol::Mat< ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType > jacobian;

    const aol::Vec3<int> gridSize = _grid.getSize();
    const RealType h = static_cast<RealType> ( _grid.H() );

    aol::Vector<RealType> dest ( ParametricDeformationType::NumberOfDeformParameters );

    for ( int j = 0; j < gridSize[1]; ++j ) {
      for ( int i = 0; i < gridSize[0]; ++i ) {
        const int index = qc::ILexCombine2( i, j, _grid.getNumX() );
        ParametricDeformedPosition<ConfiguratorType, ParametricDeformationType> defPos ( i, j, parDef, gridSize, h );

        if ( defPos.isDefPosInDomain() == false )
          continue;

        parDef.evaluateDerivativeDeformationOn01 ( defPos.getUndefPos01(), jacobian );
        const typename ConfiguratorType::VecType gradT ( _t.interpolateDX ( defPos.getDefPos() ), _t.interpolateDY ( defPos.getDefPos() ) );
        aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> tmp;

        jacobian.mult ( gradT, tmp );
        tmp *= ( _t.interpolate ( defPos.getDefPos() ) - _r.getSmoothImRef()[index] );
        for ( int k = 0; k < ParametricDeformationType::NumberOfDeformParameters; ++k )
          dest[k] += tmp[k];
      }
    }
    MDest.copySplitFrom ( dest );
  }
};

/**
 * \author Berkels
 */
template <typename RealType>
class SquaredL2DistanceEnergy : public aol::Op<aol::MultiVector<RealType>, aol::Scalar<RealType> > {
private:
  const aol::MultiVector<RealType> &_origin;
public:
  SquaredL2DistanceEnergy ( const aol::MultiVector<RealType> &Origin )
    : _origin ( Origin ) {}

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    Dest[0] *= 2;
    for ( int i = 0; i < _origin.numComponents(); ++i )
      for ( int j = 0; j < _origin[i].size(); ++j )
        Dest[0] += aol::Sqr ( MArg[i][j] - _origin[i][j] );
    Dest[0] /= 2;
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    MDest.setSum ( MArg, _origin, -1 );
  }
};
  

/**
 * \author Berkels
 */
template <typename ConfiguratorType>
typename ConfiguratorType::RealType updateDeformParameters ( const typename ConfiguratorType::InitType &Grid,
                              const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > &E,
                              const aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::MultiVector<typename ConfiguratorType::RealType> > &DE,
                              aol::MultiVector<typename ConfiguratorType::RealType> &DeformParameters,
                              const int MaxGradientDescentSteps = 1000,
                              const bool UseComponentWiseTimestep = true ) {
  typedef typename ConfiguratorType::RealType RealType;

  /*
  aol::FirstDerivativeValidator<aol::MultiVector<RealType> > tester ( E, DE, Grid.H(), aol::FirstDerivativeValidator<aol::MultiVector<RealType> >::LINEAR, 0.01 );
  tester.testDirection ( DeformParameters, "test/test" );
  tester.testAllDirections( DeformParameters, "test/test" );
  */

  aol::MultiVector<RealType> temp ( DeformParameters );
  aol::DeleteFlagPointer< aol::GradientDescentBase<RealType, aol::MultiVector<RealType> > > pGradientDescentSolver;
  typedef aol::GridlessGradientDescent<RealType, aol::MultiVector<RealType> > GradientDescentType;
  if ( UseComponentWiseTimestep ) {
    // Since GridlessGradientDescent doesn't need a grid, the Grid argument here is completely ignored
    // It's only necessary for the syntax.
    // Furthermore, we use a negative stopping epsilon here since GradientDescentComponentWiseTimestepControlled
    // only calculates tau for the different components approximately, possibly leading to a premature stopping.
    pGradientDescentSolver.reset ( new aol::GradientDescentComponentWiseTimestepControlled<ConfiguratorType, GradientDescentType> ( Grid, E, DE, MaxGradientDescentSteps, 1, -1000 ), true );
  }
  else {
    pGradientDescentSolver.reset ( new GradientDescentType ( Grid, E, DE, MaxGradientDescentSteps ), true );
  }
  pGradientDescentSolver->apply ( DeformParameters, temp );
  DeformParameters = temp;
  return pGradientDescentSolver->getEnergyAtLastPosition();
}

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
void ParametricDeformImage ( const aol::Vector<typename ConfiguratorType::RealType> &Image,
                             const typename ConfiguratorType::InitType &Grid,
                             aol::Vector<typename ConfiguratorType::RealType> &DeformedImage,
                             const aol::MultiVector<typename ConfiguratorType::RealType> &DeformParameters,
                             const bool ExtendWithConstant = true,
                             const typename ConfiguratorType::RealType ExtensionConstant = 0,
                             const bool NearestNeighborInterpolation = false ) {
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::ArrayType imageArray ( Image, Grid, aol::FLAT_COPY );

  qc::DataGenerator<ConfiguratorType> generator ( Grid );
  qc::MultiArray<RealType, ConfiguratorType::Dim> deformation ( Grid );
  ParametricDeformationType parDef ( Grid, DeformParameters );
  generator.template generateDeformationFromParametricDeformation<ParametricDeformationType, false> ( parDef, deformation );
  qc::DeformImage<ConfiguratorType> ( Image, Grid, DeformedImage, deformation, ExtendWithConstant, ExtensionConstant, NearestNeighborInterpolation );
}
  
/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricMIForce
: public aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                         ParametricMIForce<ConfiguratorType, ParametricDeformationType>,
                                                         1, ParametricDeformationType::NumberOfDeformParameters > {
  typedef typename ConfiguratorType::RealType RealType;
  const int _numberOfIntensityValues;
  const qc::ScalarArray<RealType, qc::QC_2D> &_mutualInformation;
  const typename ConfiguratorType::InitType &_grid;
  const aol::DiscreteFunctionDefault<ConfiguratorType> _r;
  const ParametricDeformationType &_parDef;
public:
  ParametricMIForce ( const typename ConfiguratorType::InitType &Grid,
                      JointHistogram<RealType> &Hist,
                      const ParametricDeformationType &ParDef )
    : aol::VectorFENonlinIntegrationVectorInterface < ConfiguratorType,
                                                      ParametricMIForce<ConfiguratorType, ParametricDeformationType>,
                                                      1, ParametricDeformationType::NumberOfDeformParameters >  ( Grid ),
      _numberOfIntensityValues ( Hist.getNumberOfIntensityValues() ),
      _mutualInformation ( Hist.getMutualInformation() ),
      _grid ( Grid ),
      _r ( Grid, Hist.getReference() ),
      _parDef ( ParDef ) {}

  void evaluateIntegrand ( const aol::auto_container<1, aol::DiscreteFunctionDefault<ConfiguratorType> > &DiscrFuncs, // DiscrFuncs[0] is supposed to contain T.
                          const typename ConfiguratorType::ElementType &El,
                          int QuadPoint, const typename ConfiguratorType::VecType &RefCoord,
                          aol::Vec<ParametricDeformationType::NumberOfDeformParameters, RealType> &Integrand ) const {

    typename ConfiguratorType::VecType transformedLocalCoord;
    qc::Element transformedEl;
    if ( !_parDef.template evaluateDeformation<false> ( El, RefCoord, transformedEl, transformedLocalCoord) ) {
      Integrand.setZero();
      return;
    }

    aol::Mat< ParametricDeformationType::NumberOfDeformParameters, ConfiguratorType::Dim, RealType > jacobian;
    _parDef.evaluateDerivativeDeformation ( El, RefCoord, jacobian );

    typename ConfiguratorType::VecType gradT;
    DiscrFuncs[0].evaluateGradient ( transformedEl, transformedLocalCoord, gradT );
    jacobian.mult ( gradT, Integrand );

    const RealType Td = DiscrFuncs[0].evaluate ( transformedEl, transformedLocalCoord );
    const RealType R = _r.evaluateAtQuadPoint ( El, QuadPoint );
    const int ir = aol::Min ( static_cast<int> ( ( _numberOfIntensityValues - 1 ) * R ), _numberOfIntensityValues - 1 );
    const int it = aol::Min ( static_cast<int> ( ( _numberOfIntensityValues - 1 ) * Td ), _numberOfIntensityValues - 1 );
    Integrand *= _mutualInformation.get ( ir, it );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricDeformationType>
class ParametricMIEnergy : public aol::Op<aol::MultiVector<typename ConfiguratorType::RealType>, aol::Scalar<typename ConfiguratorType::RealType> > {
  const typename ConfiguratorType::ArrayType _r;
  const typename ConfiguratorType::ArrayType _t;
  const qc::MIRegistrationConfigurator<ConfiguratorType> _regisConfig;
  qc::MIRegistrationEnergyWithRegardToPhi<ConfiguratorType> _MIEnergy;
public:
  typedef typename ConfiguratorType::RealType RealType;

  ParametricMIEnergy ( const typename ConfiguratorType::InitType &Grid,
                       const aol::Vector<RealType> &ImR,
                       const aol::Vector<RealType> &ImT )
    : _r ( ImR, Grid, aol::FLAT_COPY ),
      _t ( ImT, Grid, aol::FLAT_COPY ),
      _regisConfig ( 6, 1 ),
      _MIEnergy ( Grid, _r, _t, _regisConfig ) {}

  void applyAdd ( const aol::MultiVector<RealType> &MArg, aol::Scalar<RealType> &Dest ) const {
    // This implementation is not very efficient but is hopefully sufficient to test how well
    // MI works in the parametric case.
    qc::DataGenerator<ConfiguratorType> generator ( _MIEnergy.getGridReference() );
    qc::MultiArray<RealType, ConfiguratorType::Dim> deformation ( _MIEnergy.getGridReference() );
    ParametricDeformationType parDef ( _MIEnergy.getGridReference(), MArg );
    generator.template generateDeformationFromParametricDeformation<ParametricDeformationType, true> ( parDef, deformation );
    _MIEnergy.applyAdd ( deformation, Dest );
  }

  void applyDerivative ( const aol::MultiVector<RealType> &MArg, aol::MultiVector<RealType> &MDest ) const {
    // This implementation is not very efficient but is hopefully sufficient to test how well
    // MI works in the parametric case.
    typename ConfiguratorType::ArrayType deformedTemplate ( _MIEnergy.getGridReference() );
    qc::ParametricDeformImage<ConfiguratorType, ParametricDeformationType> ( _t, _MIEnergy.getGridReference(), deformedTemplate, MArg );
    JointHistogram<RealType> hist ( _r, deformedTemplate, _regisConfig.getNumberOfIntensityLevels(), _regisConfig.getBeta() );

    aol::Vector<RealType> destVec ( ParametricDeformationType::NumberOfDeformParameters );
    ParametricDeformationType parDef (  _MIEnergy.getGridReference(), MArg );
    ParametricMIForce<ConfiguratorType, ParametricDeformationType> force ( _MIEnergy.getGridReference(), hist, parDef );
    force.apply ( _t, destVec );
    MDest.copySplitFrom ( destVec );
  }
};

/**
 * \author Berkels
 */
template <typename ConfiguratorType, typename ParametricEnergyType>
typename ConfiguratorType::RealType updateDeformParameters ( const typename ConfiguratorType::InitType &Grid,
                              const aol::Vector<typename ConfiguratorType::RealType> &ImR,
                              const aol::Vector<typename ConfiguratorType::RealType> &ImT,
                              aol::MultiVector<typename ConfiguratorType::RealType> &DeformParameters,
                              const int MaxGradientDescentSteps = 1000,
                              const bool UseComponentWiseTimestep = true,
                              const typename ConfiguratorType::RealType ParameterPenaltyWeight = 0 ) {
  typedef typename ConfiguratorType::RealType RealType;
  ParametricEnergyType dataE ( Grid, ImR, ImT );
  aol::DerivativeWrapper<RealType, ParametricEnergyType, aol::MultiVector<RealType> > dataDE( dataE );
  if ( ParameterPenaltyWeight <= 0 )
    return updateDeformParameters<ConfiguratorType> ( Grid, dataE, dataDE, DeformParameters, MaxGradientDescentSteps, UseComponentWiseTimestep );
  else {
    SquaredL2DistanceEnergy<RealType> penaltyE ( DeformParameters );
    aol::DerivativeWrapper<RealType, SquaredL2DistanceEnergy<RealType>, aol::MultiVector<RealType> > penaltyDE( penaltyE );
    aol::LinCombOp<aol::MultiVector<RealType>, aol::Scalar<RealType> > E;
    E.appendReference ( dataE );
    E.appendReference ( penaltyE, ParameterPenaltyWeight );
    aol::LinCombOp<aol::MultiVector<RealType> > DE;
    DE.appendReference ( dataDE );
    DE.appendReference ( penaltyDE, ParameterPenaltyWeight );
    return updateDeformParameters<ConfiguratorType> ( Grid, E, DE, DeformParameters, MaxGradientDescentSteps, UseComponentWiseTimestep );
  }
}

/**
 * \author Berkels
 */
template <typename _ConfiguratorType, typename ParametricDeformationType, typename ParametricEnergyType>
class ParametricRegistrationMultilevelDescent : public qc::RegistrationMultilevelDescentInterfaceBase<_ConfiguratorType> {
public:
  typedef _ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef aol::MultiVector<RealType> TransformationDOFType;
  static const bool IsParametric = true;

protected:
  aol::MultiVector<RealType> _deformParameters;
  int _checkboxWidth, _maxGradientDescentSteps;
  bool _saveSteps;
  bool _useComponentWiseTimestep;
  RealType _parameterPenaltyWeight;
  RealType _energyOfLastSolution;

public:
  ParametricRegistrationMultilevelDescent ( const aol::ParameterParser &Parser, const bool TryToInitializeTranslation = false )
    : qc::RegistrationMultilevelDescentInterfaceBase<ConfiguratorType> ( Parser ),
      _deformParameters ( ParametricDeformationType::getDeformParametersSize() ),
      _checkboxWidth( Parser.getInt("checkboxWidth") ),
      _maxGradientDescentSteps( Parser.getInt("MaxGradientDescentSteps") ),
      _saveSteps( true ),
      _useComponentWiseTimestep ( Parser.getInt ( "UseComponentWiseTimestep" ) != 0 ),
      _parameterPenaltyWeight ( Parser.getRealOrDefault<RealType> ( "parameterPenaltyWeight", 0 ) ),
      _energyOfLastSolution ( aol::NumberTrait<RealType>::NaN ) {
    ParametricDeformationType::setIdentityDeformationParameters ( _deformParameters );

    if ( TryToInitializeTranslation ) {
      const aol::Vec<ConfiguratorType::Dim, RealType> centerRef = qc::getCenterOfMassOfArray<RealType, ConfiguratorType::Dim> ( this->getRefImageReference() ) / this->_grid.getNumX();
      const aol::Vec<ConfiguratorType::Dim, RealType> centerTem = qc::getCenterOfMassOfArray<RealType, ConfiguratorType::Dim> ( this->getTemplImageReference() ) / this->_grid.getNumX();
      const aol::Vec<ConfiguratorType::Dim, RealType> offset = ( centerTem - centerRef );

      for ( int i = 0; i < ConfiguratorType::Dim; ++i )
        _deformParameters[0][i] = offset[i]; 
    }
  }

  ParametricRegistrationMultilevelDescent ( const int MaxDepth, bool saveSteps = true )
    : qc::RegistrationMultilevelDescentInterfaceBase<ConfiguratorType> ( MaxDepth ),
      _deformParameters ( ParametricDeformationType::getDeformParametersSize() ),
      _checkboxWidth( 32 ),
      _maxGradientDescentSteps( 1000 ),
      _saveSteps( saveSteps ),
      _useComponentWiseTimestep ( true ),
      _parameterPenaltyWeight ( 0 ) {
    ParametricDeformationType::setIdentityDeformationParameters ( _deformParameters );
  }

  virtual ~ParametricRegistrationMultilevelDescent( ) {}

  void prolongate( ) {
    if ( this->_curLevel < this->getMaxGridDepth() )
      this->setLevel ( this->_curLevel + 1 );
  }

  aol::Vec<ParametricDeformationType::NumOfDeformParametersComponents, int> getTransformationDOFInitializer ( ) const {
    return ParametricDeformationType::getDeformParametersSize();
  }

  void setTransformation ( const aol::MultiVector<RealType> &DeformParameters ) {
    _deformParameters = DeformParameters;
  }

  void getTransformation ( aol::MultiVector<RealType> &DeformParameters ) const {
    DeformParameters = _deformParameters;
  }

  RealType getTransformationNorm ( ) const {
    return _deformParameters.norm();
  }

  void addTransformationTo ( aol::MultiVector<RealType> &DeformParameters ) const {
    ParametricDeformationType::concatenateDeformationParameters ( _deformParameters, DeformParameters );
  }

  void setTransformationToZero ( ) {
    ParametricDeformationType::setIdentityDeformationParameters ( _deformParameters );
  }

  void composeDeformations ( const aol::MultiVector<RealType> &DeformParameters1, const aol::MultiVector<RealType> &DeformParameters2, aol::MultiVector<RealType> &CompositionDeformParameters ) {
    CompositionDeformParameters = DeformParameters2;
    ParametricDeformationType::concatenateDeformationParameters ( DeformParameters1, CompositionDeformParameters );
  }

  void setTransformationToComposition ( const aol::MultiVector<RealType> &/*DeformParameters1*/, const aol::MultiVector<RealType> &/*DeformParameters2*/ ) {
    throw aol::UnimplementedCodeException ( "not implemented", __FILE__, __LINE__ );
  }

  void applyTransformation ( const aol::MultiVector<RealType> &DeformParameters, const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage, const RealType ExtensionConstant = aol::NumberTrait<RealType>::Inf, const bool NearestNeighborInterpolation = false ) const {
    qc::ParametricDeformImage<ConfiguratorType, ParametricDeformationType> ( InputImage, this->_grid, DeformedImage, DeformParameters, true, ExtensionConstant, NearestNeighborInterpolation );
  }

  void applySavedTransformation ( const char *DefBaseName, const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage, const RealType ExtensionConstant = aol::NumberTrait<RealType>::Inf, const bool NearestNeighborInterpolation = false ) const {
    aol::MultiVector<RealType> deformParameters;
    loadTransformationTo ( DefBaseName, deformParameters );
    applyTransformation ( deformParameters, InputImage, DeformedImage, ExtensionConstant, NearestNeighborInterpolation );
  }

  void loadTransformationTo ( const char *DefBaseName, aol::MultiVector<RealType> &DeformParameters ) const {
    DeformParameters.load ( aol::strprintf ( "%s%s", DefBaseName, getDeformationFileNameSuffix().c_str() ).c_str() );
  }

  void saveTransformationTo ( const aol::MultiVector<RealType> &DeformParameters, const char *DefBaseName ) const {
    DeformParameters.save ( aol::strprintf ( "%s%s", DefBaseName, getDeformationFileNameSuffix().c_str() ).c_str() );
  }

  int getMaxGradientSteps ( ) const {
    return _maxGradientDescentSteps;
  }

  void setMaxGradientSteps ( int maxGradientDescentSteps ) {
    _maxGradientDescentSteps = maxGradientDescentSteps;
  }

  void descentOnCurrentGrid( ) {
    cerr << "Registration on level " << this->_curLevel << " started\n";

    _energyOfLastSolution = updateDeformParameters<ConfiguratorType, ParametricEnergyType> ( this->getCurrentGrid(), this->_org_reference[ this->_curLevel ], this->_org_template[ this->_curLevel ], _deformParameters, _maxGradientDescentSteps, _useComponentWiseTimestep, _parameterPenaltyWeight );

    cerr << "Detected parameters are \n" << _deformParameters << endl;

    if( _saveSteps )
      saveCurrentDeformation ( );
  }


  void setCheckboxWidth ( int checkboxWidth ) {
    _checkboxWidth = checkboxWidth;
  }

  string getDeformationFileNameSuffix ( ) const {
    return ".dat";
  }

  void saveCurrentDeformation ( ) const {
    qc::DataGenerator<ConfiguratorType> generator ( this->_grid );
    qc::MultiArray<RealType, ConfiguratorType::Dim> deformation ( this->_grid );
    ParametricDeformationType parDef ( this->_grid, _deformParameters );
    generator.template generateDeformationFromParametricDeformation<ParametricDeformationType, false> ( parDef, deformation );

    qc::RegistrationStepSaver<ConfiguratorType, typename ConfiguratorType::ArrayType, false, false>
      stepSaver( this->_grid, this->getRefImageReference(), this->getTemplImageReference(), _checkboxWidth  );
    stepSaver.setSaveName ( aol::strprintf ( "_%02d", this->_curLevel ).c_str() );
    stepSaver.setSaveDirectory ( this->getSaveDirectory() );
    if ( this->getParserReference().checkAndGetBool ( "onlySaveDisplacement" ) == false )
      stepSaver.saveStep ( deformation, -1 );
    _deformParameters.save ( stepSaver.createSaveName ( "deformation", ".dat", -1 ).c_str() );
  }

  RealType getEnergyOfLastSolution ( ) const {
    return _energyOfLastSolution;
  }

  void applyCurrentTransformation ( const aol::Vector<RealType> &InputImage, aol::Vector<RealType> &DeformedImage ) const {
    applyTransformation ( _deformParameters, InputImage, DeformedImage );
  }

  void saveTransformation ( const char *FileName ) const {
    _deformParameters.save ( FileName );
  }

  void loadTransformation ( const char *FileName ) {
    _deformParameters.load ( FileName );
  }

};

/**
 * \brief Registers two images using phase correlation to find the optimal integer shift between the two images.
 *
 * \author Berkels
 */
template<typename ConfiguratorType>
class PhaseCorrelationRegistration {
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
  static const qc::Dimension Dim = ConfiguratorType::Dim;
public:
  static qc::CoordType registerImages ( const ArrayType &Reference, const ArrayType &Template ) {
    const typename ConfiguratorType::InitType grid ( Reference.getSize() );
    ArrayType phaseCorrelation ( grid );

    qc::Convolution<Dim, RealType> conv ( qc::GridSize<Dim> ( grid ).getSizeAsVecDim() );
    conv.phaseCorrelation  ( Reference, Template, phaseCorrelation );

    qc::OTFILexMapper<Dim> mapper ( grid );
    return mapper.splitGlobalIndex ( phaseCorrelation.getMaxIndexAndValue().first );
  }
};
  
/**
 *\brief Class for affine transformation. It allows translation, rotation, scaling and shearing
 *\author Tatano
 *
 */
template <typename ConfiguratorType>
class ParametricAffineTransformation2D {
public:
  static const int NumberOfDeformParameters = 2 + 2 * ConfiguratorType::Dim ;
  static const int NumOfDeformParametersComponents = 4;
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  const RealType _rotationAngle;
  aol::Matrix22<RealType> _rotationMatrix;
  aol::Matrix22<RealType> _shearingMatrix;
  typename ConfiguratorType::VecType _translation;
  typename ConfiguratorType::VecType _shearing;
  RealType _scaling;
public:
  ParametricAffineTransformation2D ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
  : _grid ( Initializer ),
  _rotationAngle ( DeformParameters[1][0] ),
  _scaling ( DeformParameters[2][0] ){
    _rotationMatrix.makeRotation ( _rotationAngle );
    _shearingMatrix.setIdentity();
    for ( int i = 0; i < ConfiguratorType::Dim; ++i ){
      _translation[i] = DeformParameters[0][i] + 0.5;
      _shearing[i] = DeformParameters[3][i];
    }
    _shearingMatrix[1][0] = _shearing[1];
    _shearingMatrix[0][1] = _shearing[0];
  }
    
  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
    DeformParameters[ConfiguratorType::Dim][0] = 1;
  }
    
  // DestDeform := ArgDeform \circ DestDeform
  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &ArgDeformPars, aol::MultiVector<RealType> &DestDeformPars ) {
    throw aol::Exception("Not implemented yet!", __FILE__, __LINE__);
  }
  
  static aol::Vec<NumOfDeformParametersComponents, int> getDeformParametersSize ( ) {
    return aol::Vec<NumOfDeformParametersComponents, int> (2, 1, 1, 2);
  }
    
  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                              const typename ConfiguratorType::VecType &RefCoord,
                              qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {
      
    typename ConfiguratorType::VecType coord, transformedCoord, tmp;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      coord[i] = El[i] + RefCoord[i] - 0.5 / _grid.H();
    }
    _shearingMatrix.mult ( coord, tmp );
    _rotationMatrix.mult ( tmp, transformedCoord );
    transformedCoord *= _scaling;
      
    return qc::transformCoord<ConfiguratorType, ClipCoord> ( _grid, transformedCoord, _translation, TransformedEl, TransformedLocalCoord );
  }
    
  void evaluateDeformationOn01 ( const ParametricAffineTransformation2D<ConfiguratorType> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    typename ConfiguratorType::VecType tmp ( Position ), tmp2(tmp, aol::STRUCT_COPY);
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      tmp[i] -= 0.5;
    _shearingMatrix.apply ( tmp, tmp2 );
    _rotationMatrix.apply ( tmp2, DeformedPosition );
    DeformedPosition *= _scaling;
    DeformedPosition += _translation;
  }
    
  void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &Position, aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    const RealType cosAlpha = _rotationMatrix[0][0];
    const RealType sinAlpha = _rotationMatrix[1][0];
      
    const RealType tempX = (Position[0] - 0.5) + _shearing[0] * (Position[1] - 0.5);
    const RealType tempY = (Position[1] - 0.5) + _shearing[1] * (Position[0] - 0.5);
      
    Jacobian[ConfiguratorType::Dim][0] = _scaling * ( -sinAlpha * tempX - cosAlpha * tempY );
    Jacobian[ConfiguratorType::Dim][1] = _scaling * (  cosAlpha * tempX - sinAlpha * tempY );
      
    Jacobian[0][0] = Jacobian[1][1] = 1;
    Jacobian[1][0] = Jacobian[0][1] = 0;
      
    Jacobian[ConfiguratorType::Dim+1][0] = ( cosAlpha * tempX - sinAlpha * tempY );
    Jacobian[ConfiguratorType::Dim+1][1] = ( sinAlpha * tempX + cosAlpha * tempY );
      
    Jacobian[ConfiguratorType::Dim+2][0] =   _scaling * cosAlpha * (Position[1] - 0.5);
    Jacobian[ConfiguratorType::Dim+2][1] = - _scaling * sinAlpha * (Position[0] - 0.5);
    Jacobian[ConfiguratorType::Dim+3][1] =   _scaling * sinAlpha * (Position[1] - 0.5);
    Jacobian[ConfiguratorType::Dim+3][0] =   _scaling * cosAlpha * (Position[0] - 0.5);
  }
    
  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &El,
                                        const typename ConfiguratorType::VecType &RefCoord,
                                        aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    typename ConfiguratorType::VecType x;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      x[i] = ( El[i] + RefCoord[i] ) * _grid.H();
    }
    evaluateDerivativeDeformationOn01 ( x, Jacobian );
  }
};
  
/**
 *\brief Class for fully affine transformation
 *\author Tatano
 *
 */
template <typename ConfiguratorType>
class FullyAffineTransformation2D {
public:
  static const int NumberOfDeformParameters = 6 ;
  static const int NumOfDeformParametersComponents = 2;
  typedef typename ConfiguratorType::RealType RealType;
private:
  const typename ConfiguratorType::InitType &_grid;
  aol::Matrix22<RealType> _affineMatrix;
  typename ConfiguratorType::VecType _translation;
public:
  FullyAffineTransformation2D ( const typename ConfiguratorType::InitType &Initializer, const aol::MultiVector<RealType> &DeformParameters )
  : _grid ( Initializer ){
    _affineMatrix[0][0] = DeformParameters[1][0];
    _affineMatrix[0][1] = DeformParameters[1][1];
    _affineMatrix[1][0] = DeformParameters[1][2];
    _affineMatrix[1][1] = DeformParameters[1][3];
    _translation[0] = DeformParameters[0][0] + 0.5;
    _translation[1] = DeformParameters[0][1] + 0.5;
  }
    
  static void setIdentityDeformationParameters ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
    DeformParameters[1][0] = 1;
    DeformParameters[1][3] = 1;
  }
  
  static void setDeformationParametersToShrinked ( aol::MultiVector<RealType> &DeformParameters ) {
    DeformParameters.setZero();
    DeformParameters[1][0] = 1./2;
    DeformParameters[1][3] = 1./2;
  }
    
  // DestDeform := ArgDeform \circ DestDeform
  static void concatenateDeformationParameters ( const aol::MultiVector<RealType> &ArgDeformPars, aol::MultiVector<RealType> &DestDeformPars ) {

    aol::Matrix22<RealType> argMat, destMat;
    argMat.fill ( ArgDeformPars[1] );
    destMat.fill ( DestDeformPars[1] );

    aol::Vec2<RealType> b, destB;
    for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
      b[i] = ArgDeformPars[0][i];
      destB[i] = DestDeformPars[0][i];
    }

    argMat.multAdd ( destB, b );
    argMat *= destMat;

    DestDeformPars[1][0] = argMat[0][0];
    DestDeformPars[1][1] = argMat[0][1];
    DestDeformPars[1][2] = argMat[1][0];
    DestDeformPars[1][3] = argMat[1][1];

    for ( int i = 0; i < ConfiguratorType::Dim; ++i )
      DestDeformPars[0][i] = b[i];
  }
  
  static aol::Vec<NumOfDeformParametersComponents, int> getDeformParametersSize ( ) {
    return aol::Vec2<int> (2, 4);
  }
    
  template <bool ClipCoord>
  bool evaluateDeformation ( const typename ConfiguratorType::ElementType &El,
                              const typename ConfiguratorType::VecType &RefCoord,
                              qc::Element &TransformedEl, typename ConfiguratorType::VecType &TransformedLocalCoord ) const {
      
    typename ConfiguratorType::VecType coord, transformedCoord;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      coord[i] = El[i] + RefCoord[i] - 0.5 / _grid.H();
    }
    _affineMatrix.apply(coord, transformedCoord);
    return qc::transformCoord<ConfiguratorType, ClipCoord> ( _grid, transformedCoord, _translation, TransformedEl, TransformedLocalCoord );
  }
    
  void evaluateDeformationOn01 ( const FullyAffineTransformation2D<ConfiguratorType> &/*ParDef*/, const typename ConfiguratorType::VecType &Position, typename ConfiguratorType::VecType &DeformedPosition ) const {
    typename ConfiguratorType::VecType tmp ( Position );
    for ( int i = 0; i < ConfiguratorType::Dim; i++ )
      tmp[i] -= 0.5;
    _affineMatrix.apply ( tmp, DeformedPosition );
    DeformedPosition += _translation;
  }
    
  void evaluateDerivativeDeformationOn01 ( const typename ConfiguratorType::VecType &Position, aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    Jacobian[0][0] = 1;                       Jacobian[0][1] = 0;
    Jacobian[1][0] = 0;                       Jacobian[1][1] = 1;
    Jacobian[2][0] = Position[0] - 0.5;       Jacobian[2][1] = 0;
    Jacobian[3][0] = Position[1] - 0.5;       Jacobian[3][1] = 0;
    Jacobian[4][0] = 0;                       Jacobian[4][1] = Position[0] - 0.5;
    Jacobian[5][0] = 0;                       Jacobian[5][1] = Position[1] - 0.5;
  }
    
  void evaluateDerivativeDeformation ( const typename ConfiguratorType::ElementType &El,
                                        const typename ConfiguratorType::VecType &RefCoord,
                                        aol::Mat< NumberOfDeformParameters, ConfiguratorType::Dim, RealType > &Jacobian ) const {
    typename ConfiguratorType::VecType x;
    for ( int i = 0; i < ConfiguratorType::Dim; i++ ) {
      x[i] = ( El[i] + RefCoord[i] ) * _grid.H();
    }
    evaluateDerivativeDeformationOn01 ( x, Jacobian );
  }
};

} // end namespace qc

#endif // __PARAMREG_H
