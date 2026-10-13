/* Copyright 2016-2020 Andrew Myers, Ann Almgren, Axel Huebl
 *                     David Grote, Jean-Luc Vay, Remi Lehe
 *                     Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "Initialization/WarpXInit.H"

#include <ablastr/profiler/ProfilerWrapper.H>
#include <ablastr/utils/timer/Timer.H>

#include <AMReX_Print.H>

#include <cstring>
#include <memory>
#include <string>

int
main (int argc, char* argv[]) {
    std::unique_ptr<char[]> relocated_input;
    if (argc >= 2) {
        std::string const relocated = std::string("../../Script/") + argv[1];
        relocated_input = std::make_unique<char[]>(relocated.size() + 1);
        std::memcpy(relocated_input.get(), relocated.c_str(), relocated.size() + 1);
        argv[1] = relocated_input.get();
    }

    warpx::initialization::initialize_external_libraries(argc, argv);
    {
        ABLASTR_PROFILE_VAR("main()", pmain);

        auto timer = ablastr::utils::timer::Timer{};
        timer.record_start_time();

        auto& warpx = WarpX::GetInstance();
        warpx.InitData();
        warpx.Evolve();
        const auto is_warpx_verbose = warpx.Verbose();
        WarpX::Finalize();

        timer.record_stop_time();
        if (is_warpx_verbose) {
            amrex::Print() << "Total Time                     : "
                           << timer.get_global_duration() << '\n';
        }

        ABLASTR_PROFILE_VAR_STOP(pmain);
    }
    warpx::initialization::finalize_external_libraries();
}
