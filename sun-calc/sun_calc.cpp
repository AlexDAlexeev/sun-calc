#include "sun_calc.h"
#include "noaa_calc.h"

namespace SolarCalc
{
SunCalc::SunCalc(const double lat, const double lon, const int timezone_offset_minutes) :
    latitude(lat), longitude(lon), tz_offset_minutes(timezone_offset_minutes)
{
}

std::optional<tm> SunCalc::GetSunriseUTC(time_t timestamp) const
{
    return GetSunriseSetUTC(true, timestamp);
}

std::optional<tm> SunCalc::GetSunsetUTC(time_t timestamp) const
{
    return GetSunriseSetUTC(false, timestamp);
}

std::optional<tm> SunCalc::GetSunrise(time_t timestamp) const
{
    return GetSunriseSet(true, timestamp);
}

std::optional<tm> SunCalc::GetSunset(time_t timestamp) const
{
    return GetSunriseSet(false, timestamp);
}

std::optional<tm> SunCalc::GetSunriseSetUTC(bool calc_sunrise, time_t timestamp) const
{
    const auto jd = NOAACalc::UnixTimeToJulianTime(timestamp);
    const auto sunrise = NOAACalc::CalcSunriseSetUTC(calc_sunrise, jd, latitude, longitude);
    if (NOAACalc::IsNumber(sunrise))
    {
        auto date = NOAACalc::CalcDateFromJD(jd);
        date.tm_hour = static_cast<int>(sunrise) / 60;
        date.tm_min = static_cast<int>(sunrise) % 60;
        return date;
    }
    return std::nullopt;
}

tm SunCalc::GetSunriseSet(bool calc_sunrise, time_t timestamp) const
{
    const auto jd = NOAACalc::UnixTimeToJulianTime(timestamp);
    const auto sunset = NOAACalc::CalcSunriseSet(calc_sunrise, jd, latitude, longitude, tz_offset_minutes);
    const auto time_jd = std::get<1>(sunset);
    if (std::get<2>(sunset) != -1.0)
    {
        auto date = NOAACalc::CalcDateFromJD(time_jd);
        date.tm_hour = static_cast<int>(time_jd) / 60;
        date.tm_min = static_cast<int>(time_jd) % 60;
        return date;
    }
    else
    {
        return NOAACalc::CalcDateTimeFromJD(std::get<0>(sunset));
    }
}
}
