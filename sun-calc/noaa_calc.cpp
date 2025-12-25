#include "noaa_calc.h"

#include <cmath>
#include <tuple>
#include <utility>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace
{
double radToDeg(const double rad)
{
    return rad * 180.0 / M_PI;
}

double degToRad(const double deg)
{
    return M_PI * deg / 180.0;
}

/**
 * Check the leap year
 * @param year 4-digit year
 * @return true if a year is a leap year
 */
bool IsLeapYear(const int year)
{
    return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
}

/**
 * Calculate the Geometric Mean Longitude of the Sun
 * @param t number of Julian centuries since J2000.0
 * @return The Geometric Mean Longitude of the Sun in degrees
 */
double CalcGeomMeanLongSun(const double t)
{
    auto L0 = 280.46646 + t * (36000.76983 + 0.0003032 * t);
    while (L0 > 360.0)
    {
        L0 -= 360.0;
    }
    while (L0 < 0.0)
    {
        L0 += 360.0;
    }
    return L0; // in degrees
}

/**
 * Calculate the Geometric Mean Anomaly of the Sun
 * @param t number of Julian centuries since J2000.0
 * @return The Geometric Mean Anomaly of the Sun in degrees
 */
double CalcGeomMeanAnomalySun(const double t)
{
    return 357.52911 + t * (35999.05029 - 0.0001537 * t); // in degrees
}

/**
 * Calculate the eccentricity of earth's orbit
 * @param t Number of Julian centuries since J2000.0
 * @return The unitless eccentricity
 */
double CalcEccentricityEarthOrbit(const double t)
{
    return 0.016708634 - t * (0.000042037 + 0.0000001267 * t); // unitless
}

/**
 * Calculate the equation of a center for the sun
 * @param t number of Julian centuries since J2000.0
 * @return degrees
 */
double CalcSunEqOfCenter(const double t)
{
    const auto m = CalcGeomMeanAnomalySun(t);

    const auto mrad = degToRad(m);
    const auto sinm = sin(mrad);
    const auto sin2m = sin(mrad + mrad);
    const auto sin3m = sin(mrad + mrad + mrad);

    return sinm * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin2m * (0.019993 - 0.000101 * t) + sin3m * 0.000289;
    // in degrees
}

/**
 * Calculate the true longitude of the sun
 * @param t number of Julian centuries since J2000.0
 * @return sun's true longitude in degrees
 */
double CalcSunTrueLong(const double t)
{
    const auto l0 = CalcGeomMeanLongSun(t);
    const auto c = CalcSunEqOfCenter(t);
    return l0 + c; // in degrees
}

/**
 * Calculate the true anomaly of the sun
 * @param t number of Julian centuries since J2000.0
 * @return  sun's true anomaly in degrees
 */
double CalcSunTrueAnomaly(const double t)
{
    const auto m = CalcGeomMeanAnomalySun(t);
    const auto c = CalcSunEqOfCenter(t);
    return m + c; // in degrees
}

/**
 * Calculate the distance to the sun in AU
 * @param t number of Julian centuries since J2000.0
 * @return sun radius vector in AUs
 */
double CalcSunRadVector(const double t)
{
    const auto v = CalcSunTrueAnomaly(t);
    const auto e = CalcEccentricityEarthOrbit(t);

    return 1.000001018 * (1 - e * e) / (1 + e * cos(degToRad(v))); // in AUs
}

/**
 * Calculate the apparent longitude of the sun
 * @param t number of Julian centuries since J2000.0
 * @return sun's apparent longitude in degrees
 */
double CalcSunApparentLong(const double t)
{
    const auto o = CalcSunTrueLong(t);

    const auto omega = 125.04 - 1934.136 * t;
    const auto lambda = o - 0.00569 - 0.00478 * sin(degToRad(omega));
    return lambda; // in degrees
}

/**
 * Calculate the mean obliquity of the ecliptic
 * @param t number of Julian centuries since J2000.0
 * @return mean obliquity in degrees
 */
double CalcMeanObliquityOfEcliptic(const double t)
{
    const auto seconds = 21.448 - t * (46.8150 + t * (0.00059 - t * 0.001813));
    const auto e0 = 23.0 + (26.0 + seconds / 60.0) / 60.0;
    return e0; // in degrees
}

/**
 * Calculate the corrected obliquity of the ecliptic
 * @param t number of Julian centuries since J2000.0
 * @return corrected obliquity in degrees
 */
double CalcObliquityCorrection(const double t)
{
    const auto e0 = CalcMeanObliquityOfEcliptic(t);

    const auto omega = 125.04 - 1934.136 * t;
    const auto e = e0 + 0.00256 * cos(degToRad(omega));
    return e; // in degrees
}

/**
 * Calculate the right ascension of the sun
 * @param t number of Julian centuries since J2000.0
 * @return sun's right ascension in degrees
 */
double CalcSunRtAscension(const double t)
{
    const auto e = CalcObliquityCorrection(t);
    const auto lambda = CalcSunApparentLong(t);
    const auto tananum = cos(degToRad(e)) * sin(degToRad(lambda));
    const auto tanadenom = cos(degToRad(lambda));
    const auto alpha = radToDeg(atan2(tananum, tanadenom));
    return alpha; // in degrees
}

/**
 * Calculate the declination of the sun
 * @param t number of Julian centuries since J2000.0
 * @return sun's declination in degrees
 */
double CalcSunDeclination(const double t)
{
    auto e = CalcObliquityCorrection(t);
    auto lambda = CalcSunApparentLong(t);
    auto sint = sin(degToRad(e)) * sin(degToRad(lambda));
    auto theta = radToDeg(asin(sint));
    return theta; // in degrees
}

/**
 * Calculate the difference between true solar time and mean solar time
 * @param t number of Julian centuries since J2000.0
 * @return equation of time in minutes of time
 */
double CalcEquationOfTime(const double t)
{
    const auto epsilon = CalcObliquityCorrection(t);
    const auto l0 = CalcGeomMeanLongSun(t);
    const auto e = CalcEccentricityEarthOrbit(t);
    const auto m = CalcGeomMeanAnomalySun(t);

    auto y = tan(degToRad(epsilon) / 2.0);
    y *= y;

    const auto sin2l0 = sin(2.0 * degToRad(l0));
    const auto sinm = sin(degToRad(m));
    const auto cos2l0 = cos(2.0 * degToRad(l0));
    const auto sin4l0 = sin(4.0 * degToRad(l0));
    const auto sin2m = sin(2.0 * degToRad(m));

    const auto Etime = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0 - 0.5 * y * y * sin4l0 - 1.25 * e *
        e * sin2m;
    return radToDeg(Etime) * 4.0; // in minutes of time
}

/**
 * Calculate the hour angle of the sun at sunrise for the latitude
 * @param lat latitude of observer in degrees
 * @param solarDec declination angle of sun in degrees
 * @return hour angle of sunrise in radians
 */
double CalcHourAngleSunrise(const double lat, const double solarDec)
{
    const auto latRad = degToRad(lat);
    const auto sdRad = degToRad(solarDec);
    const auto HAarg = cos(degToRad(90.833)) / (cos(latRad) * cos(sdRad)) - tan(latRad) * tan(sdRad);
    const auto HA = acos(HAarg);
    return HA; // in radians (for sunset, use -HA)
}

/**
 * Calculate the Universal Coordinated Time (UTC) of solar noon for the given day at the given location on earth
 * @param jd number of Julian centuries since J2000.0
 * @param longitude longitude of observer in degrees
 * @param timezone timezone offset, minutes
 * @return time in minutes from zero Z
 */
double CalcSolNoon(const double jd, const double longitude, const double timezone)
{
    const auto tnoon = NOAACalc::CalcTimeJulianCent(jd - longitude / 360.0);
    auto eqTime = CalcEquationOfTime(tnoon);
    const auto solNoonOffset = 720.0 - longitude * 4 - eqTime; // in minutes
    const auto newt = NOAACalc::CalcTimeJulianCent(jd - 0.5 + solNoonOffset / 1440.0);
    eqTime = CalcEquationOfTime(newt);
    auto solNoonLocal = 720 - longitude * 4 - eqTime + timezone; // in minutes
    while (solNoonLocal < 0.0)
    {
        solNoonLocal += 1440.0;
    }
    while (solNoonLocal >= 1440.0)
    {
        solNoonLocal -= 1440.0;
    }
    return solNoonLocal;
}
}

namespace NOAACalc
{
double UnixTimeToJulianTime(const time_t unixTime)
{
    return static_cast<double>(unixTime) / 86400.0 + 2440587.5;
}

double DateToJulianDate(int year, int month, const int day)
{
    if (month <= 2)
    {
        year -= 1;
        month += 12;
    }
    const auto A = floor(year / 100.0);
    const auto B = 2 - A + floor(A / 4.0);

    const auto JD = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + day + B - 1524.5;
    return JD;
}

time_t JulianTimeToUnixTime(const double julianTime)
{
    return static_cast<time_t>((julianTime - 2440587.5) * 86400);
}


tm GetTimeLocal(const time_t now)
{
    tm tm_buf = {};
#ifdef WIN32
    localtime_s(&tm_buf, &now);
#else
    localtime_r(&now, &tm_buf);
#endif
    return tm_buf;
}

tm GetTimeUTC(time_t now)
{
    tm tm_buf = {};
#ifdef WIN32
    gmtime_s(&tm_buf, &now);
#else
    gmtime_r(&now, &tm_buf);
#endif
    return tm_buf;
}

int CalcDayOfYear(const int year, const int month, const int day)
{
    const auto k = IsLeapYear(year) ? 1 : 2;
    const auto doy = floor(275.0 * month / 9) - k * floor((month + 9) / 12) + day - 30;
    return static_cast<int>(doy);
}

int CalcDayOfYear(const double jd)
{
    auto date = CalcDateFromJD(jd);
    const auto k = IsLeapYear(date.tm_year + 1900) ? 1 : 2;
    const auto doy = floor(275.0 * date.tm_mon / 9) - k * floor((date.tm_mon + 9) / 12) + date.tm_mday - 30;
    return static_cast<int>(doy);
}

int DaysInMonth(const int year, const int month)
{
    switch (month)
    {
        case 1:
        case 3:
        case 5:
        case 7:
        case 8:
        case 10:
        case 12:
            return 31;
        case 4:
        case 6:
        case 9:
        case 11:
            return 30;
        case 2:
            return IsLeapYear(year) ? 29 : 28;
        default:
            return 0;
    }
}

int CalcDayOfWeek(const int year, const int month, const int day)
{
    return CalcDayOfWeek(DateToJulianDate(year, month, day));
}

int CalcDayOfWeek(const double julianTime)
{
    const auto A = static_cast<int>(floor(julianTime + 1.5)) % 7;
    return A;
}

tm CalcDateFromJD(const double jd)
{
    const auto z = floor(jd + 0.5);
    const auto f = jd + 0.5 - z;

    double A;
    if (z < 2299161)
    {
        A = z;
    }
    else
    {
        const auto alpha = floor((z - 1867216.25) / 36524.25);
        A = z + 1 + alpha - floor(alpha / 4);
    }

    const auto B = A + 1524;
    const auto C = floor((B - 122.1) / 365.25);
    const auto D = floor(365.25 * C);
    const auto E = floor((B - D) / 30.6001);

    const int day = static_cast<int>(B - D - floor(30.6001 * E) + f);
    const int month = static_cast<int>(E < 14 ? E - 1 : E - 13);
    const int year = static_cast<int>(month > 2 ? C - 4716 : C - 4715);
    tm date = {};
    date.tm_year = year - 1900;
    date.tm_mon = month - 1;
    date.tm_mday = day;
    return date;
}

tm CalcDateTimeFromJD(double jd)
{
    return GetTimeLocal(JulianTimeToUnixTime(jd));
}

tm CalcDateTimeUTCFromJD(double jd)
{
    return GetTimeUTC(JulianTimeToUnixTime(jd));
}

double CalcTimeJulianCent(const double jd)
{
    return (jd - 2451545.0) / 36525.0;
}

double CalcJDFromJulianCent(const double jt)
{
    return jt * 36525.0 + 2451545.0;
}


double CalcRefraction(const double elev)
{
    double correction = 0.0;
    if (elev > 85.0)
    {
        correction = 0.0;
    }
    else
    {
        auto te = tan(degToRad(elev));
        if (elev > 5.0)
        {
            correction = 58.1 / te - 0.07 / (te * te * te) + 0.000086 / (te * te * te * te * te);
        }
        else if (elev > -0.575)
        {
            correction = 1735.0 + elev * (-518.2 + elev * (103.4 + elev * (-12.79 + elev * 0.711)));
        }
        else
        {
            correction = -20.774 / te;
        }
        correction = correction / 3600.0;
    }
    return correction;
}

std::pair<double, double> calcAzEl(const double T, const double localtime, const double latitude,
                                   const double longitude, const double timezone)
{
    const auto eqTime = CalcEquationOfTime(T);
    const auto theta = CalcSunDeclination(T);

    const auto solarTimeFix = eqTime + 4.0 * longitude - timezone;
    auto earthRadVec = CalcSunRadVector(T);
    auto trueSolarTime = localtime + solarTimeFix;
    while (trueSolarTime > 1440)
    {
        trueSolarTime -= 1440;
    }
    auto hourAngle = trueSolarTime / 4.0 - 180.0;
    if (hourAngle < -180)
    {
        hourAngle += 360.0;
    }
    const auto haRad = degToRad(hourAngle);
    auto csz = sin(degToRad(latitude)) * sin(degToRad(theta)) + cos(degToRad(latitude)) * cos(degToRad(theta)) *
        cos(haRad);
    if (csz > 1.0)
    {
        csz = 1.0;
    }
    else if (csz < -1.0)
    {
        csz = -1.0;
    }
    const auto zenith = radToDeg(acos(csz));
    const auto azDenom = cos(degToRad(latitude)) * sin(degToRad(zenith));
    auto azimuth = 0.0;
    if (abs(azDenom) > 0.001)
    {
        auto azRad = (sin(degToRad(latitude)) * cos(degToRad(zenith)) - sin(degToRad(theta))) / azDenom;
        if (abs(azRad) > 1.0)
        {
            if (azRad < 0)
            {
                azRad = -1.0;
            }
            else
            {
                azRad = 1.0;
            }
        }
        azimuth = 180.0 - radToDeg(acos(azRad));
        if (hourAngle > 0.0)
        {
            azimuth = -azimuth;
        }
    }
    else
    {
        if (latitude > 0.0)
        {
            azimuth = 180.0;
        }
        else
        {
            azimuth = 0.0;
        }
    }
    if (azimuth < 0.0)
    {
        azimuth += 360.0;
    }
    const auto exoatmElevation = 90.0 - zenith;

    // Atmospheric Refraction correction
    const auto refractionCorrection = CalcRefraction(exoatmElevation);

    const auto solarZen = zenith - refractionCorrection;
    auto elevation = 90.0 - solarZen;

    return {azimuth, elevation};
}


double CalcSunriseSetUTC(const bool rise, const double jd, const double latitude, const double longitude)
{
    const auto t = CalcTimeJulianCent(jd);
    const auto eqTime = CalcEquationOfTime(t);
    const auto solarDec = CalcSunDeclination(t);
    auto hourAngle = CalcHourAngleSunrise(latitude, solarDec);
    if (!rise)
    {
        hourAngle = -hourAngle;
    }
    const auto delta = longitude + radToDeg(hourAngle);
    const auto timeUTC = 720 - 4.0 * delta - eqTime; // in minutes
    return timeUTC;
}

double calcJDofNextPrevRiseSet(const bool next, const bool rise, const double JD, const double latitude,
                               const double longitude, const double tz)
{
    auto julianday = JD;
    auto increment = next ? 1.0 : -1.0;
    auto time = CalcSunriseSetUTC(rise, julianday, latitude, longitude);

    while (!IsNumber(time))
    {
        julianday += increment;
        time = CalcSunriseSetUTC(rise, julianday, latitude, longitude);
    }
    auto timeLocal = time + tz * 60.0;
    while (timeLocal < 0.0 || timeLocal >= 1440.0)
    {
        auto incr = timeLocal < 0 ? 1 : -1;
        timeLocal += incr * 1440.0;
        julianday -= incr;
    }

    return julianday;
}

std::tuple<double, double, double> CalcSunriseSet(const bool rise, const double jd, const double latitude,
                                                  const double longitude, const double timezone)
{
    auto timeUTC = CalcSunriseSetUTC(rise, jd, latitude, longitude);
    auto newTimeUTC = CalcSunriseSetUTC(rise, jd + timeUTC / 1440.0, latitude, longitude);
    auto jday = jd;
    auto timeLocal = 0.0;
    auto azimuth = -1.0;
    if (IsNumber(newTimeUTC))
    {
        timeLocal = newTimeUTC + timezone;
        auto riseT = CalcTimeJulianCent(jd + newTimeUTC / 1440.0);
        auto riseAzEl = calcAzEl(riseT, timeLocal, latitude, longitude, timezone);
        azimuth = riseAzEl.second;
        if (timeLocal < 0.0 || timeLocal >= 1440.0)
        {
            auto increment = timeLocal < 0 ? 1 : -1;
            while (timeLocal < 0.0 || timeLocal >= 1440.0)
            {
                timeLocal += increment * 1440.0;
                jday -= increment;
            }
        }
    }
    else
    {
        // no sunrise/set found

        azimuth = -1.0;
        auto doy = CalcDayOfYear(jd);
        if ((latitude > 66.4 && doy > 79 && doy < 267) ||
            (latitude < -66.4 && (doy < 83 || doy > 263)))
        {
            //previous sunrise/next sunset
            jday = calcJDofNextPrevRiseSet(!rise, rise, jd, latitude, longitude, timezone);
        }
        else
        {
            //previous sunset/next sunrise
            jday = calcJDofNextPrevRiseSet(rise, rise, jd, latitude, longitude, timezone);
        }
    }

    return std::make_tuple(jday, timeLocal, azimuth);
}
}
