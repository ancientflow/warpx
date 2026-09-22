# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_hall_rotating_axis  # name
    3  # dims
    1  # nprocs
    inputs_test_3d_hall_rotating_axis  # inputs
    "analysis_rotating_axis.py diags/diag000000"  # analysis
    OFF  # checksum
    OFF  # dependency
)
