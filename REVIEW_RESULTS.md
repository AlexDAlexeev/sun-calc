# Review results

## Findings

### High — One CI workflow invokes an unsupported CMake option

`/home/runner/work/sun-calc/sun-calc/.github/workflows/c-cpp.yml:30` runs `cmake --test`. CMake has no `--test` mode, so this step exits with an argument error instead of running tests. The other workflow correctly invokes `ctest` at `/home/runner/work/sun-calc/sun-calc/.github/workflows/cmake-multi-platform.yml:75`.

### Medium — Polar fallback ignores the configured timezone offset

When an event does not occur on the requested date, `/home/runner/work/sun-calc/sun-calc/sun-calc/sun_calc.cpp:57-59` converts the fallback Julian date with `CalcDateTimeFromJD`. That helper calls `GetTimeLocal` in `/home/runner/work/sun-calc/sun-calc/sun-calc/noaa_calc.cpp:395-397`, which applies the machine's local timezone rather than `tz_offset_minutes`. Results therefore depend on the host timezone and can disagree with the fixed offset supplied to `SunCalc`. The existing fallback test uses offset zero and derives its expected result from the same helper (`/home/runner/work/sun-calc/sun-calc/test/main.cpp:125-154`).

### Medium — UTC event results do not account for UTC date rollover

`GetSunriseSetUTC` in `/home/runner/work/sun-calc/sun-calc/sun-calc/sun_calc.cpp:37-40` always uses the input Julian date for the returned calendar date and stores the calculated minutes directly as `tm_hour` and `tm_min`. For valid longitudes where the event falls before 00:00 or at/after 24:00 UTC, the returned `tm` can have out-of-range time fields and the wrong date, despite the method promising a UTC time.

### Medium — Local event rollover adjusts the day but discards it

`NOAACalc::CalcSunriseSet` adjusts its returned Julian day when local event minutes cross midnight (`/home/runner/work/sun-calc/sun-calc/sun-calc/noaa_calc.cpp:573-580`). The caller then constructs the result's date from the original input Julian day instead of the adjusted day (`/home/runner/work/sun-calc/sun-calc/sun-calc/sun_calc.cpp:49-55`). The returned local time can consequently be paired with the wrong calendar date.

### Medium — Invalid coordinates can make local-time calls loop indefinitely

`CalcHourAngleSunrise` applies `acos` without validating its domain (`/home/runner/work/sun-calc/sun-calc/sun-calc/noaa_calc.cpp:228-234`). For invalid coordinates, such as a latitude outside `[-90, 90]`, this can yield a non-finite event time. The no-event fallback loop in `calcJDofNextPrevRiseSet` keeps advancing days until it gets a finite result, with no bound or invalid-input exit (`/home/runner/work/sun-calc/sun-calc/sun-calc/noaa_calc.cpp:536-547`). The public constructor does not validate coordinates (`/home/runner/work/sun-calc/sun-calc/sun-calc/sun_calc.cpp:6-9`).

## Coverage and maintainability gaps

- The solar API tests calculate expected values by calling the same NOAA calculation functions used by production code (`/home/runner/work/sun-calc/sun-calc/test/main.cpp:95-111`, `136-141`, `167-172`). They verify wrapper behavior but do not independently validate numerical accuracy. Add reference fixtures and cases for UTC/local midnight rollover and nonzero-offset polar fallback.
- GoogleTest is fetched from the moving `main` branch (`/home/runner/work/sun-calc/sun-calc/CMakeLists.txt:7-11`), so identical project revisions can configure against different dependency revisions.
- The library target does not publish its header directory as a usage requirement (`/home/runner/work/sun-calc/sun-calc/CMakeLists.txt:19-24`); only the test target adds it explicitly (`/home/runner/work/sun-calc/sun-calc/CMakeLists.txt:32-35`). Consumers linking the target do not automatically inherit the public include path.

## Scope and limitations

Reviewed the C++ implementation and headers, tests, CMake configuration, workflows, README, and repository configuration. This was a static review; no independent numerical comparison against NOAA reference outputs or cross-timezone execution was performed.
