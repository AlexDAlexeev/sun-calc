[![CMake on multiple platforms](https://github.com/AlexDAlexeev/sun-calc/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/AlexDAlexeev/sun-calc/actions/workflows/cmake-multi-platform.yml)

# Sunset/Sunrise calculator
Translated to C++ from the https://gml.noaa.gov/grad/solcalc/index.html.

# Usage

```cpp
const auto longitude = 93; // longitude in degrees
const auto latitude = 74; // latitude in degrees
const int timezone = +1*60; // timezone in minutes

SolarCalc::SunCalc sun_calc(latitude, longitude, timezone);

auto sunrise = sun_calc.GetSunrise(); // returns a std::optional<tm> filled with sunrise time and date
auto sunset = sun_calc.GetSunset(); // returns a std::optional<tm> filled with sunset time and date

auto sunriseUTC = sun_calc.GetSunriseUTC(); // returns a time_t filled with sunrise time in UTC, if  there is no sunrise, returns empty
auto sunsetUTC = sun_calc.GetSunsetUTC(); // returns a time_t filled with sunset time in UTC, if  there is no sunset, returns empty
```

