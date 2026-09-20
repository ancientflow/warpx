foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    if(D STREQUAL "3")
        target_compile_definitions(lib_${SD} PUBLIC WARPX_USE_HALL_ELECTROSTATIC_MATERIALS)
    endif()
    target_sources(lib_${SD}
      PRIVATE
        Core/WarpXInsert.cpp
        Background/BackgroundCoupledDensity.cpp
        Background/InsertBackgroundDensity.cpp
        Boundary/AnalyticBoundaryGeometry.cpp
        Boundary/AnalyticBoundaryInteraction.cpp
        Boundary/InsertBoundaryParticles.cpp
        Boundary/InsertBoundaryPhi.cpp
        Boundary/ZMinWallCharge.cpp
        Collisions/CoupledBackgroundMCCCollision.cpp
        Collisions/IonizationSourceTable.cpp
        Diagnostics/InsertRuntimeDiagnostics.cpp
        Diagnostics/ZmaxRadialExitStats.cpp
        Fields/ECDIChargeFilter.cpp
        Fields/HallElectrostaticMaterial.cpp
        Fields/SpectralBoundarySchur.cpp
        Geometry/TruncatedConeSurfaceSampler.cpp
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
