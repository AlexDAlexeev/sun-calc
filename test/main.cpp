#include <gtest/gtest.h>
#include <iomanip>
#include <tuple>

#include "noaa_calc.h"
#include "sun_calc.h"


static tm GetTimeLocal(const time_t now = std::time(nullptr))
{
    struct tm tm_buf = {};
#ifdef WIN32
    localtime_s(&tm_buf, &now);
#else
    localtime_r(&now, &tm_buf);
#endif
    return tm_buf;
}


TEST(Dates, DayOfYear)
{
    // Non-Leap year
    const int nl_year = 2019;
    int num_day = 1;
    for (int m = 1; m <= 12; ++m)
    {
        for (int d = 1; d <= NOAACalc::DaysInMonth(nl_year, m); ++d)
        {
            ASSERT_EQ(NOAACalc::CalcDayOfYear(nl_year, m, d), num_day);
            num_day++;
        }
    }
    const int l_year = 2020;
    num_day = 1;
    for (int m = 1; m <= 12; ++m)
    {
        for (int d = 1; d <= NOAACalc::DaysInMonth(l_year, m); ++d)
        {
            ASSERT_EQ(NOAACalc::CalcDayOfYear(l_year, m, d), num_day);
            num_day++;
        }
    }
}

TEST(Dates, DayOfWeek)
{
    ASSERT_EQ(NOAACalc::CalcDayOfWeek(1970, 1, 1), 4);
    ASSERT_EQ(NOAACalc::CalcDayOfWeek(2440587.5), 4);
    ASSERT_EQ(NOAACalc::CalcDayOfWeek(2025, 12, 23), 2);
}

TEST(Dates, GregorianToJulian)
{
    EXPECT_EQ(2440587.5, NOAACalc::UnixTimeToJulianTime(0));
    EXPECT_EQ(2440587.5, NOAACalc::DateToJulianDate(1970, 1, 1));
    EXPECT_EQ(0, NOAACalc::JulianTimeToUnixTime(2440587.5));
}

TEST(Dates, JulianToDate)
{
    const auto date = NOAACalc::CalcDateFromJD(2440587.5);
    EXPECT_EQ(70, date.tm_year);
    EXPECT_EQ(0, date.tm_mon);
    EXPECT_EQ(1, date.tm_mday);

    const auto now = time(nullptr);
    const auto expected_now = GetTimeLocal(now);
    const auto nowJT = NOAACalc::UnixTimeToJulianTime(now);
    const auto nowDate = NOAACalc::CalcDateFromJD(nowJT);
    EXPECT_EQ(nowDate.tm_year, expected_now.tm_year);
    EXPECT_EQ(nowDate.tm_mon, expected_now.tm_mon);
    EXPECT_EQ(nowDate.tm_mday, expected_now.tm_mday);
}

static time_t MakeUtcDateTimestamp(const int year, const int month, const int day)
{
    const auto jd = NOAACalc::DateToJulianDate(year, month, day);
    return NOAACalc::JulianTimeToUnixTime(jd);
}

TEST(SunCalc, GetSunriseUTCAndSunsetUTC)
{
    constexpr double lat = 37.7749;
    constexpr double lon = -122.4194;
    const auto ts = MakeUtcDateTimestamp(2024, 6, 21);

    const SolarCalc::SunCalc calc(lat, lon, 0);
    const auto sunriseUtc = calc.GetSunriseUTC(ts);
    const auto sunsetUtc = calc.GetSunsetUTC(ts);

    ASSERT_TRUE(sunriseUtc.has_value());
    ASSERT_TRUE(sunsetUtc.has_value());

    const auto jd = NOAACalc::UnixTimeToJulianTime(ts);
    const auto expectedDate = NOAACalc::CalcDateFromJD(jd);

    const auto expectedSunriseMinutes = static_cast<int>(NOAACalc::CalcSunriseSetUTC(true, jd, lat, lon));
    const auto expectedSunsetMinutes = static_cast<int>(NOAACalc::CalcSunriseSetUTC(false, jd, lat, lon));

    EXPECT_EQ(expectedDate.tm_year, sunriseUtc->tm_year);
    EXPECT_EQ(expectedDate.tm_mon, sunriseUtc->tm_mon);
    EXPECT_EQ(expectedDate.tm_mday, sunriseUtc->tm_mday);
    EXPECT_EQ(expectedSunriseMinutes / 60, sunriseUtc->tm_hour);
    EXPECT_EQ(expectedSunriseMinutes % 60, sunriseUtc->tm_min);

    EXPECT_EQ(expectedDate.tm_year, sunsetUtc->tm_year);
    EXPECT_EQ(expectedDate.tm_mon, sunsetUtc->tm_mon);
    EXPECT_EQ(expectedDate.tm_mday, sunsetUtc->tm_mday);
    EXPECT_EQ(expectedSunsetMinutes / 60, sunsetUtc->tm_hour);
    EXPECT_EQ(expectedSunsetMinutes % 60, sunsetUtc->tm_min);
}

TEST(SunCalc, GetSunriseUTCAndSunsetUTCWhenNoEvent)
{
    constexpr double lat = 89.0;
    constexpr double lon = 0.0;
    const auto ts = MakeUtcDateTimestamp(2024, 12, 21);

    const SolarCalc::SunCalc calc(lat, lon, 0);
    EXPECT_FALSE(calc.GetSunriseUTC(ts).has_value());
    EXPECT_FALSE(calc.GetSunsetUTC(ts).has_value());
}

TEST(SunCalc, GetSunriseAndSunsetLocalWithFallback)
{
    constexpr double lat = 89.0;
    constexpr double lon = 0.0;
    constexpr int tzOffsetMinutes = 0;
    const auto ts = MakeUtcDateTimestamp(2024, 12, 21);

    const SolarCalc::SunCalc calc(lat, lon, tzOffsetMinutes);
    const auto sunriseLocal = calc.GetSunrise(ts);
    const auto sunsetLocal = calc.GetSunset(ts);

    const auto jd = NOAACalc::UnixTimeToJulianTime(ts);
    const auto expectedSunriseTuple = NOAACalc::CalcSunriseSet(true, jd, lat, lon, tzOffsetMinutes);
    const auto expectedSunsetTuple = NOAACalc::CalcSunriseSet(false, jd, lat, lon, tzOffsetMinutes);

    const auto expectedSunrise = NOAACalc::CalcDateTimeFromJD(std::get<0>(expectedSunriseTuple));
    const auto expectedSunset = NOAACalc::CalcDateTimeFromJD(std::get<0>(expectedSunsetTuple));

    EXPECT_EQ(expectedSunrise.tm_year, sunriseLocal.tm_year);
    EXPECT_EQ(expectedSunrise.tm_mon, sunriseLocal.tm_mon);
    EXPECT_EQ(expectedSunrise.tm_mday, sunriseLocal.tm_mday);
    EXPECT_EQ(expectedSunrise.tm_hour, sunriseLocal.tm_hour);
    EXPECT_EQ(expectedSunrise.tm_min, sunriseLocal.tm_min);

    EXPECT_EQ(expectedSunset.tm_year, sunsetLocal.tm_year);
    EXPECT_EQ(expectedSunset.tm_mon, sunsetLocal.tm_mon);
    EXPECT_EQ(expectedSunset.tm_mday, sunsetLocal.tm_mday);
    EXPECT_EQ(expectedSunset.tm_hour, sunsetLocal.tm_hour);
    EXPECT_EQ(expectedSunset.tm_min, sunsetLocal.tm_min);
}

TEST(SunCalc, GetSunriseAndSunsetLocal)
{
    constexpr double lat = 51.5074;
    constexpr double lon = -0.1278;
    constexpr int tzOffsetMinutes = 60;
    const auto ts = MakeUtcDateTimestamp(2024, 3, 21);

    const SolarCalc::SunCalc calc(lat, lon, tzOffsetMinutes);
    const auto sunriseLocal = calc.GetSunrise(ts);
    const auto sunsetLocal = calc.GetSunset(ts);

    const auto jd = NOAACalc::UnixTimeToJulianTime(ts);
    const auto expectedSunriseTuple = NOAACalc::CalcSunriseSet(true, jd, lat, lon, tzOffsetMinutes);
    const auto expectedSunsetTuple = NOAACalc::CalcSunriseSet(false, jd, lat, lon, tzOffsetMinutes);

    const auto expectedSunriseDate = NOAACalc::CalcDateFromJD(std::get<1>(expectedSunriseTuple));
    const auto expectedSunsetDate = NOAACalc::CalcDateFromJD(std::get<1>(expectedSunsetTuple));

    EXPECT_EQ(expectedSunriseDate.tm_year, sunriseLocal.tm_year);
    EXPECT_EQ(expectedSunriseDate.tm_mon, sunriseLocal.tm_mon);
    EXPECT_EQ(expectedSunriseDate.tm_mday, sunriseLocal.tm_mday);
    EXPECT_EQ(static_cast<int>(std::get<1>(expectedSunriseTuple)) / 60, sunriseLocal.tm_hour);
    EXPECT_EQ(static_cast<int>(std::get<1>(expectedSunriseTuple)) % 60, sunriseLocal.tm_min);

    EXPECT_EQ(expectedSunsetDate.tm_year, sunsetLocal.tm_year);
    EXPECT_EQ(expectedSunsetDate.tm_mon, sunsetLocal.tm_mon);
    EXPECT_EQ(expectedSunsetDate.tm_mday, sunsetLocal.tm_mday);
    EXPECT_EQ(static_cast<int>(std::get<1>(expectedSunsetTuple)) / 60, sunsetLocal.tm_hour);
    EXPECT_EQ(static_cast<int>(std::get<1>(expectedSunsetTuple)) % 60, sunsetLocal.tm_min);
}
