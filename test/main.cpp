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

TEST(Dates, GrigorianToJulian)
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
