#pragma once

#include <cmath>
#include <ctime>
#include <optional>
#include <tuple>

namespace NOAACalc
{
    tm GetTimeLocal(time_t now = std::time(nullptr));
    tm GetTimeUTC(time_t now = std::time(nullptr));

    /**
     * Finds numerical day-of-year from month, day and 4-digit year info
     * @param year 4-digit year
     * @param month month number, January is 1
     * @param day  day in month number, starting from 1
     * @return The numerical day of year
     */
    int CalcDayOfYear(int year, int month, int day);
    /**
     * Returns number of days in a month for a given year
     * @param year 4-digit year
     * @param month month number, January is 1
     * @return Number of days in the month
     */
    int DaysInMonth(int year, int month);
    /**
     * Returns day of week for the given date
     * @param year 4-digit year
     * @param month month number, January is 1
     * @param day day in month number, starting from 1
     * @return Day of the week, Sunday is 0
     */
    int CalcDayOfWeek(int year, int month, int day);
    /**
     * Returns day of week for the given date in Julian format
     * @param julianTime julian date
     * @return Day of the week, Sunday is 0
     */
    int CalcDayOfWeek(double julianTime);

    /**
     * Julian day from Unix timestamp. Number is returned for the start of the day.
     * @param unixTime Unix time
     * @return Julian time
     */
    double UnixTimeToJulianTime(time_t unixTime);

    /**
    * Julian day from calendar day. Number is returned for the start of the day.  Fractional days should be	added later.
     * @param year 4-digit year
     * @param month month number, January is 1
     * @param day day in month number, starting from 1
     * @return Julian date
     */
    double DateToJulianDate(int year, int month, int day);

    /**
     * Convert Julian date to Unix timestamp
     * @param julianTime julian time
     * @return Unix timestamp
     */
    time_t JulianTimeToUnixTime(double julianTime);

    /**
     * Calendar date from Julian date
     * @param jd Julian date
     * @return tm filled with year, month, day
     */
    tm CalcDateFromJD(double jd);

    /**
     * Calendar local date and time from Julian date
     * @param jd Julian date
     * @return tm filled time
     */
    tm CalcDateTimeFromJD(double jd);

    /**
     * Calendar UTC date and time from Julian date
     * @param jd Julian date
     * @return tm filled time
     */
    tm CalcDateTimeUTCFromJD(double jd);

    /**
     * Convert Julian Day to centuries since J2000.0.
     * @param jd Julian date
     * @return century
     */
    double CalcTimeJulianCent(double jd);

    /**
     * Convert centuries since J2000.0 to Julian Day
     * @param jt number of Julian centuries since J2000.0
     * @return Julian date
     */
    double CalcJDFromJulianCent(double jt);

    inline bool IsNumber(const double value)
    {
        return !std::isnan(value) && std::isfinite(value);
    }

    /**
     * Calculate the Universal Coordinated Time (UTC) of sunset for the given day at the given location on earth
     * If there is no sunrise or sunset for the given date returns NaN
     * @param rise true for sunrise, false for sunset
     * @param jd julian day
     * @param latitude latitude of observer in degrees
     * @param longitude longitude of observer in degrees
     * @return time in minutes from zero Z
     */
    double CalcSunriseSetUTC(bool rise, double jd, double latitude, double longitude);

    /**
     * Calculate the local time of sunset or sunrise for the given day at the given location on earth
     * @param rise true for sunrise, false for sunset
     * @param jd julian day
     * @param latitude latitude of observer in degrees
     * @param longitude longitude of observer in degrees
     * @param timezone observer timezone offset
     * @return
     */
    std::tuple<double, double, double> CalcSunriseSet(bool rise, double jd, double latitude, double longitude,
                                                      double timezone);
}
