# Sunset/Sunrise calculator
Translated to C++ from the https://gml.noaa.gov/grad/solcalc/index.html.

# Usage

```cpp
const auto longitude = 93; // longitude in degrees
const auto lattitude = 74; // latitude in degrees
const int timezone = +1*60; // timezone in minutes

SolarCalc::SunCalc sun_calc(lattitude, longitude, timezone);

auto sunrise = sun_calc.getSunrise(); // returns a tm struct filled with sunrise time and date
auto sunset = sun_calc.getSunset(); // returns a tm struct filled with sunset time and date

auto sunriseUTC = sun_calc.getSunriseUTC(); // returns a time_t filled with sunrise time in UTC, if  there is no sunrise, returns empty
auto sunsetUTC = sun_calc.getSunsetUTC(); // returns a time_t filled with sunset time in UTC, if  there is no sunset, returns empty
```

