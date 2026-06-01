foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        Core/WarpXInsert.cpp
        Background/BackgroundCoupledDensity.cpp
        Background/InsertBackgroundDensity.cpp
        Boundary/InsertBoundaryParticles.cpp
        Boundary/InsertBoundaryPhi.cpp
        Collisions/CoupledBackgroundMCCCollision.cpp
        Diagnostics/InsertRuntimeDiagnostics.cpp
        Injection/InsertInjection.cpp
    )
endforeach()
