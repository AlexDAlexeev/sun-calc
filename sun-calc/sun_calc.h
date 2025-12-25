#pragma once
#include <optional>

namespace SolarCalc
{
/**
 * Sunrise/sunset calculator helper class
 */
class SunCalc
{
public:
    /**
     *
     * @param lat desired latitude
     * @param lon desired longitude
     * @param timezone_offset_minutes Timezone offset in minutes
     */
    SunCalc(double lat, double lon, int timezone_offset_minutes);

    /**
     * Determines the UTC time of sunrise. In certain situations, a sunrise may not occur for the given date and coordinates.
     * @param timestamp Date timestamp
     * @return Empty if there is no sunrise, tm struct filled with hours and minutes
     */
    [[nodiscard]] std::optional<tm> GetSunriseUTC(time_t timestamp) const;
    /**
     * Determines the UTC time of sunset. In certain situations, a sunset may not occur for the given date and coordinates.
     * @param timestamp Date timestamp
     * @return Empty if there is no sunset, tm struct filled with hours and minutes
     */
    [[nodiscard]] std::optional<tm> GetSunsetUTC(time_t timestamp) const;

    /**
     * Determines the local time of sunrise. In case if there is no sunrise on a given date returns the nearest sunrise date and time
     * @param timestamp Date timestamp
     * @return tm filled with the sunrise date and time
     */
    [[nodiscard]] tm GetSunrise(time_t timestamp) const;
    /**
     * Determines the local time of sunset. In case if there is no sunset on a given date returns the nearest sunset date and time
     * @param timestamp Date timestamp
     * @return tm filled with the sunset date and time
     */
    [[nodiscard]] tm GetSunset(time_t timestamp) const;

private:
    double latitude = 0.0;
    double longitude = 0.0;
    int tz_offset_minutes = 0;

    [[nodiscard]] std::optional<tm> GetSunriseSetUTC(bool calc_sunrise, time_t timestamp) const;
    [[nodiscard]] tm GetSunriseSet(bool calc_sunrise, time_t timestamp) const;
};
}
