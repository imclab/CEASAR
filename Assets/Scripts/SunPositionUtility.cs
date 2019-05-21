using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public class SunPositionUtility : MonoBehaviour
{


    double radToDeg(double angleRad)
    {
        return Mathf.Rad2Deg * angleRad; // (180.0f * angleRad / Math.PI);
    }

    double degToRad(double angleDeg)
    {
        return Mathf.Deg2Rad * angleDeg; // (Math.PI * angleDeg / 180.0f);
    }

    double calcGeomMeanLongSun(double t)
    {
        var L0 = 280.46646f + t * (36000.76983f + t * (0.0003032f));
        while (L0 > 360.0f)
        {
            L0 -= 360.0f;
        }
        while (L0 < 0.0f)
        {
            L0 += 360.0f;
        }
        return L0;     // in degrees
    }

    double calcGeomMeanAnomalySun(double t)
    {
        var M = 357.52911f + t * (35999.05029f - 0.0001537f * t);
        return M;     // in degrees
    }

    double calcEccentricityEarthOrbit(double t)
    {
        var e = 0.016708634f - t * (0.000042037f + 0.0000001267f * t);
        return e;     // unitless
    }

    double calcSunEqOfCenter(double t)
    {
        var m = calcGeomMeanAnomalySun(t);
        var mrad = degToRad(m);
        var sinm = Math.Sin(mrad);
        var sin2m = Math.Sin(mrad + mrad);
        var sin3m = Math.Sin(mrad + mrad + mrad);
        var C = sinm * (1.914602f - t * (0.004817f + 0.000014f * t)) + sin2m * (0.019993f - 0.000101f * t) + sin3m * 0.000289f;
        return C;     // in degrees
    }

    double calcSunTrueLong(double t)
    {
        var l0 = calcGeomMeanLongSun(t);
        var c = calcSunEqOfCenter(t);
        var O = l0 + c;
        return O;     // in degrees
    }

    double calcSunTrueAnomaly(double t)
    {
        var m = calcGeomMeanAnomalySun(t);
        var c = calcSunEqOfCenter(t);
        var v = m + c;
        return v;     // in degrees
    }

    double calcSunRadVector(double t)
    {
        var v = calcSunTrueAnomaly(t);
        var e = calcEccentricityEarthOrbit(t);
        var R = (1.000001018f * (1 - e * e)) / (1 + e * Math.Cos(degToRad(v)));
        return R;     // in AUs
    }

    double calcSunApparentLong(double t)
    {
        var o = calcSunTrueLong(t);
        var omega = 125.04f - 1934.136f * t;
        var lambda = o - 0.00569f - 0.00478f * Math.Sin(degToRad(omega));
        return lambda;        // in degrees
    }

    double calcMeanObliquityOfEcliptic(double t)
    {
        var seconds = 21.448f - t * (46.8150f + t * (0.00059f - t * (0.001813f)));
        var e0 = 23.0f + (26.0f + (seconds / 60.0f)) / 60.0f;
        return e0;        // in degrees
    }

    double calcObliquityCorrection(double t)
    {
        var e0 = calcMeanObliquityOfEcliptic(t);
        var omega = 125.04f - 1934.136f * t;
        var e = e0 + 0.00256f * Math.Cos(degToRad(omega));
        return e;     // in degrees
    }

    double calcSunRtAscension(double t)
    {
        var e = calcObliquityCorrection(t);
        var lambda = calcSunApparentLong(t);
        var tananum = (Math.Cos(degToRad(e)) * Math.Sin(degToRad(lambda)));
        var tanadenom = (Math.Cos(degToRad(lambda)));
        var alpha = radToDeg(Math.Atan2(tananum, tanadenom));
        return alpha;     // in degrees
    }

    double calcSunDeclination(double t)
    {
        var e = calcObliquityCorrection(t);
        var lambda = calcSunApparentLong(t);

        var sint = Math.Sin(degToRad(e)) * Math.Sin(degToRad(lambda));
        var theta = radToDeg(Math.Asin(sint));
        return theta;     // in degrees
    }

    double calcEquationOfTime(double t)
    {
        var epsilon = calcObliquityCorrection(t);
        var l0 = calcGeomMeanLongSun(t);
        var e = calcEccentricityEarthOrbit(t);
        var m = calcGeomMeanAnomalySun(t);

        var y = Math.Tan(degToRad(epsilon) / 2.0f);
        y *= y;

        var sin2l0 = Math.Sin(2.0f * degToRad(l0));
        var sinm = Math.Sin(degToRad(m));
        var cos2l0 = Math.Cos(2.0f * degToRad(l0));
        var sin4l0 = Math.Sin(4.0f * degToRad(l0));
        var sin2m = Math.Sin(2.0f * degToRad(m));

        var Etime = y * sin2l0 - 2.0f * e * sinm + 4.0f * e * y * sinm * cos2l0 - 0.5f * y * y * sin4l0 - 1.25f * e * e * sin2m;
        return radToDeg(Etime) * 4.0f;  // in minutes of time
    }

    double calcHourAngleSunrise(double lat, double solarDec)
    {
        var latRad = degToRad(lat);
        var sdRad = degToRad(solarDec);
        var HAarg = (Math.Cos(degToRad(90.833f)) / (Math.Cos(latRad) * Math.Cos(sdRad)) - Math.Tan(latRad) * Math.Tan(sdRad));
        var HA = Math.Acos(HAarg);
        return HA;        // in radians (for sunset, use -HA)
    }

    double calcAzEl(bool output, double T, double localtime, double latitude, double longitude, double zone)
    {
        var eqTime = calcEquationOfTime(T);
        var theta = calcSunDeclination(T);
        if (output)
        {
            Debug.Log("eqtbox " + Math.Floor(eqTime * 100 + 0.5f) / 100.0f);
            Debug.Log("sdbox " + Math.Floor(theta * 100 + 0.5f) / 100.0f);
        }
        var solarTimeFix = eqTime + 4.0f * longitude - 60.0f * zone;
        var earthRadVec = calcSunRadVector(T);
        var trueSolarTime = localtime + solarTimeFix;
        while (trueSolarTime > 1440)
        {
            trueSolarTime -= 1440;
        }
        var hourAngle = trueSolarTime / 4.0f - 180.0f;
        if (hourAngle < -180)
        {
            hourAngle += 360.0f;
        }
        var haRad = degToRad(hourAngle);
        var csz = Math.Sin(degToRad(latitude)) * Math.Sin(degToRad(theta)) + Math.Cos(degToRad(latitude)) * Math.Cos(degToRad(theta)) * Math.Cos(haRad);
        if (csz > 1.0f)
        {
            csz = 1.0f;
        }
        else if (csz < -1.0f)
        {
            csz = -1.0f;
        }
        var zenith = radToDeg(Math.Acos(csz));
        var azDenom = (Math.Cos(degToRad(latitude)) * Math.Sin(degToRad(zenith)));
        double azRad = 0;
        double azimuth;
        if (Math.Abs(azDenom) > 0.001f)
        {
            azRad = ((Math.Sin(degToRad(latitude)) * Math.Cos(degToRad(zenith))) - Math.Sin(degToRad(theta))) / azDenom;
            if (Math.Abs(azRad) > 1.0f)
            {
                azRad = azRad < 0 ? -1.0f : 1.0f;

            }
            azimuth = 180.0f - radToDeg(Math.Acos(azRad));
            if (hourAngle > 0.0f)
            {
                azimuth = -azimuth;
            }
        }
        else
        {
            if (latitude > 0.0f)
            {
                azimuth = 180.0f;
            }
            else
            {
                azimuth = 0.0f;
            }
        }
        if (azimuth < 0.0f)
        {
            azimuth += 360.0f;
        }
        var exoatmElevation = 90.0f - zenith;

        // Atmospheric Refraction correction
        double refractionCorrection = 0.0f;
        if (exoatmElevation > 85.0f)
        {
            refractionCorrection = 0.0f;
        }
        else
        {
            var te = Math.Tan(degToRad(exoatmElevation));
            if (exoatmElevation > 5.0f)
            {
                refractionCorrection = 58.1f / te - 0.07f / (te * te * te) + 0.000086f / (te * te * te * te * te);
            }
            else if (exoatmElevation > -0.575f)
            {
                refractionCorrection = 1735.0f + exoatmElevation * (-518.2f + exoatmElevation * (103.4f + exoatmElevation * (-12.79f + exoatmElevation * 0.711f)));
            }
            else
            {
                refractionCorrection = -20.774f / te;
            }
            refractionCorrection = refractionCorrection / 3600.0f;
        }

        var solarZen = zenith - refractionCorrection;

        if ((output) && (solarZen > 108.0f))
        {
            Debug.Log("azbox " + "dark");
            Debug.Log("elbox " + "dark");
        }
        else if (output)
        {
            Debug.Log("azbox " + Math.Floor(azimuth * 100 + 0.5f) / 100.0f);
            Debug.Log("elbox " + Math.Floor((90.0f - solarZen) * 100 + 0.5f) / 100.0f);
        }
        return (azimuth);
    }

    double calcSolNoon(double jd, double longitude, double timezone, double dst)
    {
        var tnoon = TimeConverter.JulianDateToJulianCentuary(jd - longitude / 360.0f);
        var eqTime = calcEquationOfTime(tnoon);
        var solNoonOffset = 720.0f - (longitude * 4) - eqTime; // in minutes
        var newt = TimeConverter.JulianDateToJulianCentuary(jd + solNoonOffset / 1440.0f);
        eqTime = calcEquationOfTime(newt);
        var solNoonLocal = 720 - (longitude * 4) - eqTime + (timezone * 60.0f);// in minutes
        if (dst != 0) solNoonLocal += 60.0f;
        while (solNoonLocal < 0.0f)
        {
            solNoonLocal += 1440.0f;
        }
        while (solNoonLocal >= 1440.0f)
        {
            solNoonLocal -= 1440.0f;
        }
        return solNoonLocal;
    }

    double calcSunriseSetUTC(bool rise, double JD, double latitude, double longitude)
    {
        var t = TimeConverter.JulianDateToJulianCentuary(JD);
        var eqTime = calcEquationOfTime(t);
        var solarDec = calcSunDeclination(t);
        var hourAngle = calcHourAngleSunrise(latitude, solarDec);
        //alert("HA = " + radToDeg(hourAngle));
        if (!rise) hourAngle = -hourAngle;
        var delta = longitude + radToDeg(hourAngle);
        var timeUTC = 720 - (4.0f * delta) - eqTime;    // in minutes
        return timeUTC;
    }

    void calcSunriseSet(bool rise, double JD, double latitude, double longitude, double timezone, bool dst)
    {
        var timeUTC = calcSunriseSetUTC(rise, JD, latitude, longitude);
        var newTimeUTC = calcSunriseSetUTC(rise, JD + timeUTC / 1440.0f, latitude, longitude);
        var timeLocal = newTimeUTC + (timezone * 60.0f);
        if (rise)
        {
            var riseT = TimeConverter.JulianDateToJulianCentuary(JD + newTimeUTC / 1440.0f);
            var riseAz = calcAzEl(false, riseT, timeLocal, latitude, longitude, timezone);
            //if (rise)
            //    {
            // showLineGeodesic2("sunrise", "#00aa00", riseAz, latitude, longitude);
            //    }
            //    else
            //    {
            // showLineGeodesic2("sunset", "#ff0000", riseAz, latitude, longitude);
            //    }
            //}
            timeLocal += ((dst) ? 60.0f : 0.0f);
            if ((timeLocal >= 0.0f) && (timeLocal < 1440.0f))
            {
                Debug.Log(timeLocal);
            }
            else
            {
                var jday = JD;
                var increment = ((timeLocal < 0) ? 1 : -1);
                while ((timeLocal < 0.0f) || (timeLocal >= 1440.0f))
                {
                    timeLocal += increment * 1440.0f;
                    jday -= increment;
                }
                Debug.Log(TimeConverter.JulianToCalendarDate(jday).AddMinutes(timeLocal));
            }
        }
        else
        { // no sunrise/set found
            var doy = TimeConverter.CalcDoyFromJD(JD);
            double jdy = 0;
            if (((latitude > 66.4f) && (doy > 79) && (doy < 267)) ||
          ((latitude < -66.4f) && ((doy < 83) || (doy > 263))))
            {   //previous sunrise/next sunset
                if (rise)
                { // find previous sunrise
                    jdy = calcJDofNextPrevRiseSet(false, rise, JD, latitude, longitude, timezone, dst);
                }
                else
                { // find next sunset
                    jdy = calcJDofNextPrevRiseSet(true, rise, JD, latitude, longitude, timezone, dst);
                }
                Debug.Log(TimeConverter.JulianToCalendarDate(jdy));
            }
            else
            {   //previous sunset/next sunrise
                if (rise)
                { // find previous sunrise
                    jdy = calcJDofNextPrevRiseSet(true, rise, JD, latitude, longitude, timezone, dst);
                }
                else
                { // find next sunset
                    jdy = calcJDofNextPrevRiseSet(false, rise, JD, latitude, longitude, timezone, dst);
                }
                Debug.Log(TimeConverter.JulianToCalendarDate(jdy));
            }
        }

    }
    double calcJDofNextPrevRiseSet(bool next, bool rise, double JD, double latitude, double longitude, double tz, bool dst)
    {
        var julianday = JD;
        var increment = ((next) ? 1.0f : -1.0f);

        var time = calcSunriseSetUTC(rise, julianday, latitude, longitude);
        while (Double.IsNaN(time))
        {
            julianday += increment;
            time = calcSunriseSetUTC(rise, julianday, latitude, longitude);
        }
        var timeLocal = time + tz * 60.0f + ((dst) ? 60.0f : 0.0f);
        while ((timeLocal < 0.0f) || (timeLocal >= 1440.0f))
        {
            var incr = ((timeLocal < 0) ? 1 : -1);
            timeLocal += (incr * 1440.0f);
            julianday -= incr;
        }
        return julianday;
    }

    // void calculate()
    // {
    //     var jday = getJD(DateTime.Now);
    //     var tl = getTimeLocal();
    //     var tz = -5;
    //     var dst = true;
    //     var total = jday + tl / 1440.0f - tz / 24.0f;
    //     var T = calcTimeJulianCent(total);
    //     var lat = 51;
    //     var lng = 47;
    //     calcAzEl(1, T, tl, lat, lng, tz);
    //     calcSolNoon(jday, lng, tz, dst);
    //     var rise = calcSunriseSet(1, jday, lat, lng, tz, dst);
    //     var set = calcSunriseSet(0, jday, lat, lng, tz, dst);
    //     //alert("JD " + jday + "  " + rise + "  " + set + "  ")
    // }
}
