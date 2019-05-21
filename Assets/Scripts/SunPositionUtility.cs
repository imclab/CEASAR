using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public class SunPositionUtility
{
    public SunPositionUtility(){}
    public double calcGeomMeanLongSun(double t)
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

    public double calcGeomMeanAnomalySun(double t)
    {
        var M = 357.52911f + t * (35999.05029f - 0.0001537f * t);
        return M;     // in degrees
    }

    public double calcEccentricityEarthOrbit(double t)
    {
        var e = 0.016708634f - t * (0.000042037f + 0.0000001267f * t);
        return e;     // unitless
    }

    public double calcSunEqOfCenter(double t)
    {
        var m = calcGeomMeanAnomalySun(t);
        var mrad = Mathf.Deg2Rad *(m);
        var sinm = Math.Sin(mrad);
        var sin2m = Math.Sin(mrad + mrad);
        var sin3m = Math.Sin(mrad + mrad + mrad);
        var C = sinm * (1.914602f - t * (0.004817f + 0.000014f * t)) + sin2m * (0.019993f - 0.000101f * t) + sin3m * 0.000289f;
        return C;     // in degrees
    }

    public double calcSunTrueLong(double t)
    {
        var l0 = calcGeomMeanLongSun(t);
        var c = calcSunEqOfCenter(t);
        var O = l0 + c;
        return O;     // in degrees
    }

    public double calcSunTrueAnomaly(double t)
    {
        var m = calcGeomMeanAnomalySun(t);
        var c = calcSunEqOfCenter(t);
        var v = m + c;
        return v;     // in degrees
    }

    public double calcSunRadVector(double t)
    {
        var v = calcSunTrueAnomaly(t);
        var e = calcEccentricityEarthOrbit(t);
        var R = (1.000001018f * (1 - e * e)) / (1 + e * Math.Cos(Mathf.Deg2Rad *(v)));
        return R;     // in AUs
    }

    public double calcSunApparentLong(double t)
    {
        var o = calcSunTrueLong(t);
        var omega = 125.04f - 1934.136f * t;
        var lambda = o - 0.00569f - 0.00478f * Math.Sin(Mathf.Deg2Rad *(omega));
        return lambda;        // in degrees
    }

    public double calcMeanObliquityOfEcliptic(double t)
    {
        var seconds = 21.448f - t * (46.8150f + t * (0.00059f - t * (0.001813f)));
        var e0 = 23.0f + (26.0f + (seconds / 60.0f)) / 60.0f;
        return e0;        // in degrees
    }

    public double calcObliquityCorrection(double t)
    {
        var e0 = calcMeanObliquityOfEcliptic(t);
        var omega = 125.04f - 1934.136f * t;
        var e = e0 + 0.00256f * Math.Cos(Mathf.Deg2Rad *(omega));
        return e;     // in degrees
    }

    public double calcSunRtAscension(double t)
    {
        var e = calcObliquityCorrection(t);
        var lambda = calcSunApparentLong(t);
        var tananum = (Math.Cos(Mathf.Deg2Rad *(e)) * Math.Sin(Mathf.Deg2Rad *(lambda)));
        var tanadenom = (Math.Cos(Mathf.Deg2Rad *(lambda)));
        var alpha = Mathf.Rad2Deg * Math.Atan2(tananum, tanadenom);
        return alpha;     // in degrees
    }

    public double calcSunDeclination(double t)
    {
        var e = calcObliquityCorrection(t);
        var lambda = calcSunApparentLong(t);

        var sint = Math.Sin(Mathf.Deg2Rad *(e)) * Math.Sin(Mathf.Deg2Rad *(lambda));
        var theta = Mathf.Rad2Deg * (Math.Asin(sint));
        return theta;     // in degrees
    }

    public double calcEquationOfTime(double t)
    {
        var epsilon = calcObliquityCorrection(t);
        var l0 = calcGeomMeanLongSun(t);
        var e = calcEccentricityEarthOrbit(t);
        var m = calcGeomMeanAnomalySun(t);

        var y = Math.Tan(Mathf.Deg2Rad *(epsilon) / 2.0f);
        y *= y;

        var sin2l0 = Math.Sin(2.0f * Mathf.Deg2Rad *(l0));
        var sinm = Math.Sin(Mathf.Deg2Rad *(m));
        var cos2l0 = Math.Cos(2.0f * Mathf.Deg2Rad *(l0));
        var sin4l0 = Math.Sin(4.0f * Mathf.Deg2Rad *(l0));
        var sin2m = Math.Sin(2.0f * Mathf.Deg2Rad *(m));

        var Etime = y * sin2l0 - 2.0f * e * sinm + 4.0f * e * y * sinm * cos2l0 - 0.5f * y * y * sin4l0 - 1.25f * e * e * sin2m;
        return  Mathf.Rad2Deg * (Etime) * 4.0f;  // in minutes of time
    }

    public double calcHourAngleSunrise(double lat, double solarDec)
    {
        var latRad = Mathf.Deg2Rad *(lat);
        var sdRad = Mathf.Deg2Rad *(solarDec);
        var HAarg = (Math.Cos(Mathf.Deg2Rad *(90.833f)) / (Math.Cos(latRad) * Math.Cos(sdRad)) - Math.Tan(latRad) * Math.Tan(sdRad));
        var HA = Math.Acos(HAarg);
        return HA;        // in radians (for sunset, use -HA)
    }

    public double calcAzEl(bool output, double T, double localtime, double latitude, double longitude, double zone)
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
        var haRad = Mathf.Deg2Rad *(hourAngle);
        var csz = Math.Sin(Mathf.Deg2Rad *(latitude)) * Math.Sin(Mathf.Deg2Rad *(theta)) + Math.Cos(Mathf.Deg2Rad *(latitude)) * Math.Cos(Mathf.Deg2Rad *(theta)) * Math.Cos(haRad);
        if (csz > 1.0f)
        {
            csz = 1.0f;
        }
        else if (csz < -1.0f)
        {
            csz = -1.0f;
        }
        var zenith =  Mathf.Rad2Deg *(Math.Acos(csz));
        var azDenom = (Math.Cos(Mathf.Deg2Rad *(latitude)) * Math.Sin(Mathf.Deg2Rad *(zenith)));
        double azRad = 0;
        double azimuth;
        if (Math.Abs(azDenom) > 0.001f)
        {
            azRad = ((Math.Sin(Mathf.Deg2Rad *(latitude)) * Math.Cos(Mathf.Deg2Rad *(zenith))) - Math.Sin(Mathf.Deg2Rad *(theta))) / azDenom;
            if (Math.Abs(azRad) > 1.0f)
            {
                azRad = azRad < 0 ? -1.0f : 1.0f;

            }
            azimuth = 180.0f -  Mathf.Rad2Deg *(Math.Acos(azRad));
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
            var te = Math.Tan(Mathf.Deg2Rad *(exoatmElevation));
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

    public double calcSolNoon(double jd, double longitude, double timezone, double dst)
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

    public double calcSunriseSetUTC(bool rise, double JD, double latitude, double longitude)
    {
        var t = TimeConverter.JulianDateToJulianCentuary(JD);
        var eqTime = calcEquationOfTime(t);
        var solarDec = calcSunDeclination(t);
        var hourAngle = calcHourAngleSunrise(latitude, solarDec);
        //alert("HA = " +  Mathf.Rad2Deg *(hourAngle));
        if (!rise) hourAngle = -hourAngle;
        var delta = longitude +  Mathf.Rad2Deg *(hourAngle);
        var timeUTC = 720 - (4.0f * delta) - eqTime;    // in minutes
        return timeUTC;
    }

    public void calcSunriseSet(bool rise, double JD, double latitude, double longitude, double timezone, bool dst)
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
    public double calcJDofNextPrevRiseSet(bool next, bool rise, double JD, double latitude, double longitude, double tz, bool dst)
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
