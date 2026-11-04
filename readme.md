# Sun Calc

[![CMake on multiple platforms](https://github.com/AlexDAlexeev/sun-calc/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/AlexDAlexeev/sun-calc/actions/workflows/cmake-multi-platform.yml)

A C++20 library for calculating sunrise and sunset times from a date and geographic coordinates. Its solar calculations are based on NOAA's [Solar Calculation details](https://gml.noaa.gov/grad/solcalc/calcdetails.html).

## Features

- Sunrise and sunset as local wall-clock times for a fixed UTC offset.
- Sunrise and sunset as UTC times, with an empty result when the event does not occur on the requested date.
- Local-time methods return the nearest occurrence when the event does not occur on the requested date.

## Requirements

- CMake 3.25 or newer
- A C++20 compiler
- Internet access when configuring the project for the first time, to fetch GoogleTest

## Build and test

From the repository root, configure and build the library and tests, then run CTest:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

The static library target is named `sun-calc`; its public headers are in `sun-calc/`.

## Usage

Pass latitude and longitude in degrees, followed by a fixed timezone offset in minutes east of UTC. For example, `60` means UTC+01:00 and `-300` means UTC−05:00. The library does not determine a location's timezone or apply daylight-saving rules.

The methods take a Unix timestamp (`std::time_t`) identifying the date to calculate. For a specific UTC calendar date, use a timestamp at the start of that date in UTC.

```cpp
#include "sun_calc.h"

#include <ctime>
#include <iomanip>
#include <iostream>

int main()
{
    constexpr double latitude = 37.7749;
    constexpr double longitude = -122.4194;
    constexpr int timezoneOffsetMinutes = 0;

    const SolarCalc::SunCalc sunCalc(latitude, longitude, timezoneOffsetMinutes);
    const std::time_t date = std::time(nullptr);

    const auto sunriseUtc = sunCalc.GetSunriseUTC(date);
    const auto sunsetUtc = sunCalc.GetSunsetUTC(date);

    if (sunriseUtc)
    {
        std::cout << "Sunrise (UTC): " << std::put_time(&*sunriseUtc, "%Y-%m-%d %H:%M") << '\n';
    }
    else
    {
        std::cout << "No sunrise on this date\n";
    }

    if (sunsetUtc)
    {
        std::cout << "Sunset (UTC): " << std::put_time(&*sunsetUtc, "%Y-%m-%d %H:%M") << '\n';
    }
    else
    {
        std::cout << "No sunset on this date\n";
    }

    const auto sunriseLocal = sunCalc.GetSunrise(date);
    const auto sunsetLocal = sunCalc.GetSunset(date);
    std::cout << "Sunrise (fixed-offset local time): "
              << std::put_time(&sunriseLocal, "%Y-%m-%d %H:%M") << '\n';
    std::cout << "Sunset (fixed-offset local time): "
              << std::put_time(&sunsetLocal, "%Y-%m-%d %H:%M") << '\n';
}
```

### `SolarCalc::SunCalc` API

| Method | Return type | Behavior |
| --- | --- | --- |
| `GetSunriseUTC(timestamp)` | `std::optional<tm>` | UTC sunrise; empty if there is no sunrise on the date. |
| `GetSunsetUTC(timestamp)` | `std::optional<tm>` | UTC sunset; empty if there is no sunset on the date. |
| `GetSunrise(timestamp)` | `tm` | Local sunrise; returns the nearest sunrise when none occurs on the date. |
| `GetSunset(timestamp)` | `tm` | Local sunset; returns the nearest sunset when none occurs on the date. |

The returned `tm` values contain the calendar date and time; `tm_year` is years since 1900 and `tm_mon` is zero-based, following the standard C/C++ `tm` convention.

## License

No license file is currently included in this repository. Check with the maintainers before redistributing the library.
