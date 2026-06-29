foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        Core/WarpXInsert.cpp
        Background/BackgroundCoupledDensity.cpp
        Background/InsertBackgroundDensity.cpp
        Boundary/InsertBoundaryParticles.cpp
        Boundary/InsertBoundaryPhi.cpp
        Boundary/ZMinWallCharge.cpp
        Collisions/CoupledBackgroundMCCCollision.cpp
        Collisions/IonizationSourceTable.cpp
        Diagnostics/InsertRuntimeDiagnostics.cpp
        Fields/ECDIChargeFilter.cpp
        Fields/SpectralBoundarySchur.cpp
        Injection/DistributionSampler1D.cpp
        Injection/HallCoordinateDistribution.cpp
        Injection/HallCoordinateTransform.cpp
        Injection/HallDistribution1D.cpp
        Injection/HallInjectionSource.cpp
        Injection/HallInjector.cpp
        Injection/HallRateModel.cpp
        Injection/IonizationSourceSampler.cpp
        Injection/InsertInjection.cpp
        Utils/InsertUtils.cpp
    )
endforeach()
