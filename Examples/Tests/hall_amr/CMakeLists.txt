# Small coupled checks for the direct-density Insert configuration.
foreach(case IN ITEMS single multi mpi)
    set(levels 1)
    set(ranks 1)
    if(case STREQUAL "single")
        set(levels 0)
    elseif(case STREQUAL "mpi")
        set(ranks 2)
    endif()
    set(mesh_args "amr.max_level=${levels}")
    if(levels GREATER 0)
        string(APPEND mesh_args
            " amr.ref_ratio_vect='2 2 2'"
            " warpx.fine_tag_lo='-0.002 -0.002 -0.002'"
            " warpx.fine_tag_hi='0.002 0.002 0.002'"
        )
    endif()
    set(test_name test_3d_hall_amr_${case})
    add_warpx_test(
        ${test_name}
        3
        ${ranks}
        "inputs_test_3d_hall_amr ${mesh_args}"
        "analysis.py ${levels}"
        OFF
        OFF
    )
    if(TEST ${test_name}.run)
        configure_file(
            elastic.dat
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_name}/elastic.dat
            COPYONLY
        )
    endif()
endforeach()
