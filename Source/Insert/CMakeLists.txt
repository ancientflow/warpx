foreach(D IN LISTS WarpX_DIMS)
    warpx_set_suffix_dims(SD ${D})
    target_sources(lib_${SD}
      PRIVATE
        WarpXInsert.cpp
        BackgroundCoupledDensity.cpp
        InsertBackgroundDensity.cpp
        InsertBoundaryParticles.cpp
        InsertBoundaryPhi.cpp
        InsertInjection.cpp
        InsertRuntimeDiagnostics.cpp
    )
endforeach()
