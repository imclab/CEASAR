using System;
using System.Text;

public struct SkyPosition
{
    public double Azimuth;
    public double Elevation;
    public SkyPosition(double azimuth, double elevation)
    {
        Azimuth = azimuth;
        Elevation = elevation;
    }
}
public class SunPositionUtility
{
    public SunPositionUtility()
    {
        outputMessages = new StringBuilder();
    }
    private StringBuilder outputMessages;
    public StringBuilder OutputMessages
    {
        get { return outputMessages; }
    }


    public double Deg2Rad(double deg)
    {
        return Math.PI * deg / 180;
    }
    public double Rad2Deg(double rad)
    {
        return 180.0 * rad / Math.PI;
    }
    public double calcGeomMeanLongSun(double t)
    {
        double L0 = 280.46646 + t * (36000.76983 + t * (0.0003032));
        while (L0 > 360.0)
        {
            L0 -= 360.0;
        }
        while (L0 < 0.0)
        {
            L0 += 360.0;
        }
        return L0;     // in degrees
    }

    public double calcGeomMeanAnomalySun(double t)
    {
        double M = 357.52911 + t * (35999.05029 - 0.0001537 * t);
        return M;     // in degrees
    }

    public double calcEccentricityEarthOrbit(double t)
    {
        double e = 0.016708634 - t * (0.000042037 + 0.0000001267 * t);
        return e;     // unitless
    }

    public double calcSunEqOfCenter(double t)
    {
        double m = calcGeomMeanAnomalySun(t);
        double mrad = Deg2Rad(m);
        double sinm = Math.Sin(mrad);
        double sin2m = Math.Sin(mrad + mrad);
        double sin3m = Math.Sin(mrad + mrad + mrad);
        double C = sinm * (1.914602 - t * (0.004817 + 0.000014 * t)) + sin2m * (0.019993 - 0.000101 * t) + sin3m * 0.000289;
        return C;     // in degrees
    }

    public double calcSunTrueLong(double t)
    {
        double l0 = calcGeomMeanLongSun(t);
        double c = calcSunEqOfCenter(t);
        double O = l0 + c;
        return O;     // in degrees
    }

    public double calcSunTrueAnomaly(double t)
    {
        double m = calcGeomMeanAnomalySun(t);
        double c = calcSunEqOfCenter(t);
        double v = m + c;
        return v;     // in degrees
    }

    public double calcSunRadVector(double t)
    {
        double v = calcSunTrueAnomaly(t);
        double e = calcEccentricityEarthOrbit(t);
        double R = (1.000001018 * (1 - e * e)) / (1 + e * Math.Cos(Deg2Rad(v)));
        return R;     // in AUs
    }

    public double calcSunApparentLong(double t)
    {
        double o = calcSunTrueLong(t);
        double omega = 125.04 - 1934.136 * t;
        double lambda = o - 0.00569 - 0.00478 * Math.Sin(Deg2Rad(omega));
        return lambda;        // in degrees
    }

    public double calcMeanObliquityOfEcliptic(double t)
    {
        double seconds = 21.448 - t * (46.8150 + t * (0.00059 - t * (0.001813)));
        double e0 = 23.0 + (26.0 + (seconds / 60.0)) / 60.0;
        return e0;        // in degrees
    }

    public double calcObliquityCorrection(double t)
    {
        double e0 = calcMeanObliquityOfEcliptic(t);
        double omega = 125.04 - 1934.136 * t;
        double e = e0 + 0.00256 * Math.Cos(Deg2Rad(omega));
        return e;     // in degrees
    }

    public double calcSunRtAscension(double t)
    {
        double e = calcObliquityCorrection(t);
        double lambda = calcSunApparentLong(t);
        double tananum = (Math.Cos(Deg2Rad(e)) * Math.Sin(Deg2Rad(lambda)));
        double tanadenom = (Math.Cos(Deg2Rad(lambda)));
        double alpha = Rad2Deg(Math.Atan2(tananum, tanadenom));
        return alpha;     // in degrees
    }

    public double calcSunDeclination(double t)
    {
        double e = calcObliquityCorrection(t);
        double lambda = calcSunApparentLong(t);

        double sint = Math.Sin(Deg2Rad(e)) * Math.Sin(Deg2Rad(lambda));
        double theta = Rad2Deg(Math.Asin(sint));
        return theta;     // in degrees
    }

    public double calcEquationOfTime(double t)
    {
        double epsilon = calcObliquityCorrection(t);
        double l0 = calcGeomMeanLongSun(t);
        double e = calcEccentricityEarthOrbit(t);
        double m = calcGeomMeanAnomalySun(t);

        double y = Math.Tan(Deg2Rad(epsilon) / 2.0);
        y *= y;

        double sin2l0 = Math.Sin(2.0 * Deg2Rad(l0));
        double sinm = Math.Sin(Deg2Rad(m));
        double cos2l0 = Math.Cos(2.0 * Deg2Rad(l0));
        double sin4l0 = Math.Sin(4.0 * Deg2Rad(l0));
        double sin2m = Math.Sin(2.0 * Deg2Rad(m));

        double Etime = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0 - 0.5 * y * y * sin4l0 - 1.25 * e * e * sin2m;
        return Rad2Deg(Etime) * 4.0;  // in minutes of time
    }

    public double calcHourAngleSunrise(double lat, double solarDec)
    {
        double latRad = Deg2Rad(lat);
        double sdRad = Deg2Rad(solarDec);
        double HAarg = (Math.Cos(Deg2Rad(90.833)) / (Math.Cos(latRad) * Math.Cos(sdRad)) - Math.Tan(latRad) * Math.Tan(sdRad));
        double HA = Math.Acos(HAarg);
        return HA;        // in radians (for sunset, use -HA)
    }

    public SkyPosition calcAzEl(bool output, double T, double localtime, double latitude, double longitude, double zone)
    {
        double eqTime = calcEquationOfTime(T);
        double theta = calcSunDeclination(T);
        if (output)
        {
            outputMessages.Append("Equation of Time: " + Math.Floor(eqTime * 100 + 0.5) / 100.0 + "\n");
            outputMessages.Append("Sun Declination " + Math.Floor(theta * 100 + 0.5) / 100.0 + "\n");
        }
        double solarTimeFix = eqTime + 4.0 * longitude - 60.0 * zone;
        double earthRadVec = calcSunRadVector(T);
        double trueSolarTime = localtime + solarTimeFix;
        outputMessages.Append("true solar time " + trueSolarTime + "\n");
        while (trueSolarTime > 1440)
        {
            trueSolarTime -= 1440;
        }
        double hourAngle = trueSolarTime / 4.0 - 180.0;
        if (hourAngle < -180)
        {
            hourAngle += 360.0;
        }
        double haRad = Deg2Rad(hourAngle);
        double csz = Math.Sin(Deg2Rad(latitude)) * Math.Sin(Deg2Rad(theta)) + Math.Cos(Deg2Rad(latitude)) * Math.Cos(Deg2Rad(theta)) * Math.Cos(haRad);
        if (csz > 1.0)
        {
            csz = 1.0;
        }
        else if (csz < -1.0f)
        {
            csz = -1.0;
        }
        double zenith = Rad2Deg(Math.Acos(csz));
        double azDenom = (Math.Cos(Deg2Rad(latitude)) * Math.Sin(Deg2Rad(zenith)));
        double azRad = 0;
        double azimuth;
        if (Math.Abs(azDenom) > 0.001)
        {
            azRad = ((Math.Sin(Deg2Rad(latitude)) * Math.Cos(Deg2Rad(zenith))) - Math.Sin(Deg2Rad(theta))) / azDenom;
            if (Math.Abs(azRad) > 1.0)
            {
                azRad = azRad < 0 ? -1.0 : 1.0;

            }
            azimuth = 180.0 - Rad2Deg(Math.Acos(azRad));
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
        double exoatmElevation = 90.0 - zenith;

        // Atmospheric Refraction correction
        double refractionCorrection = 0.0;
        if (exoatmElevation > 85.0)
        {
            refractionCorrection = 0.0;
        }
        else
        {
            double te = Math.Tan(Deg2Rad(exoatmElevation));
            if (exoatmElevation > 5.0)
            {
                refractionCorrection = 58.1 / te - 0.07 / (te * te * te) + 0.000086 / (te * te * te * te * te);
            }
            else if (exoatmElevation > -0.575)
            {
                refractionCorrection = 1735.0 + exoatmElevation * (-518.2 + exoatmElevation * (103.4 + exoatmElevation * (-12.79 + exoatmElevation * 0.711)));
            }
            else
            {
                refractionCorrection = -20.774 / te;
            }
            refractionCorrection = refractionCorrection / 3600.0;
        }

        double solarZen = zenith - refractionCorrection;

        double az = Math.Floor(azimuth * 100 + 0.5) / 100.0;
        double el = Math.Floor((90.0 - solarZen) * 100 + 0.5) / 100.0;

        if (output)
        {
            string dark = (solarZen > 108.0) ? " (dark)" : "";
            outputMessages.Append("Azimuth " + az + dark + "\n");
            outputMessages.Append("Elevation " + el + dark + "\n");
        }
        return new SkyPosition(az, el);
    }

    public double calcSolNoon(double jd, double longitude, double timezone, double dst)
    {
        double tnoon = TimeConverter.calcTimeJulianCent(jd - longitude / 360.0);
        double eqTime = calcEquationOfTime(tnoon);
        double solNoonOffset = 720.0 - (longitude * 4) - eqTime; // in minutes
        double newt = TimeConverter.calcTimeJulianCent(jd + solNoonOffset / 1440.0);
        eqTime = calcEquationOfTime(newt);
        double solNoonLocal = 720 - (longitude * 4) - eqTime + (timezone * 60.0);// in minutes
        if (dst > 0 || dst < 0) solNoonLocal += 60.0;
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

    public double calcSunriseSetUTC(bool rise, double JD, double latitude, double longitude)
    {
        double t = TimeConverter.calcTimeJulianCent(JD);
        double eqTime = calcEquationOfTime(t);
        double solarDec = calcSunDeclination(t);
        double hourAngle = calcHourAngleSunrise(latitude, solarDec);
        //alert("HA = " +  Rad2Deg(hourAngle));
        if (!rise) hourAngle = -hourAngle;
        double delta = longitude + Rad2Deg(hourAngle);
        double timeUTC = 720 - (4.0 * delta) - eqTime;    // in minutes
        return timeUTC;
    }

    public void calcSunriseSet(bool rise, double JD, double latitude, double longitude, double timezone, bool dst)
    {
        double timeUTC = calcSunriseSetUTC(rise, JD, latitude, longitude);
        double newTimeUTC = calcSunriseSetUTC(rise, JD + timeUTC / 1440.0, latitude, longitude);
        double timeLocal = newTimeUTC + (timezone * 60.0);
        if (rise)
        {
            double riseT = TimeConverter.calcTimeJulianCent(JD + newTimeUTC / 1440.0);
            SkyPosition pos = calcAzEl(false, riseT, timeLocal, latitude, longitude, timezone);
            double riseAz = pos.Azimuth;
            //if (rise)
            //    {
            // showLineGeodesic2("sunrise", "#00aa00", riseAz, latitude, longitude);
            //    }
            //    else
            //    {
            // showLineGeodesic2("sunset", "#ff0000", riseAz, latitude, longitude);
            //    }
            //}
            timeLocal += ((dst) ? 60.0 : 0.0);
            if ((timeLocal >= 0.0) && (timeLocal < 1440.0))
            {
                outputMessages.Append(timeLocal + "\n");
            }
            else
            {
                double jday = JD;
                double increment = ((timeLocal < 0) ? 1 : -1);
                while ((timeLocal < 0.0) || (timeLocal >= 1440.0))
                {
                    timeLocal += increment * 1440.0;
                    jday -= increment;
                }
                outputMessages.Append(TimeConverter.JulianToCalendarDate(jday).AddMinutes(timeLocal) + "\n");
            }
        }
        else
        { // no sunrise/set found
            double doy = TimeConverter.CalcDoyFromJD(JD);
            double jdy = 0;
            if (((latitude > 66.4) && (doy > 79) && (doy < 267)) || ((latitude < -66.4) && ((doy < 83) || (doy > 263))))
            {   //previous sunrise/next sunset
                if (rise)
                { // find previous sunrise
                    jdy = calcJDofNextPrevRiseSet(false, rise, JD, latitude, longitude, timezone, dst);
                }
                else
                { // find next sunset
                    jdy = calcJDofNextPrevRiseSet(true, rise, JD, latitude, longitude, timezone, dst);
                }
                outputMessages.Append(TimeConverter.JulianToCalendarDate(jdy) + "\n");
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
                outputMessages.Append(TimeConverter.JulianToCalendarDate(jdy) + "\n");
            }
        }

    }
    public double calcJDofNextPrevRiseSet(bool next, bool rise, double JD, double latitude, double longitude, double tz, bool dst)
    {
        double julianday = JD;
        double increment = ((next) ? 1.0 : -1.0);

        double time = calcSunriseSetUTC(rise, julianday, latitude, longitude);
        while (Double.IsNaN(time))
        {
            julianday += increment;
            time = calcSunriseSetUTC(rise, julianday, latitude, longitude);
        }
        double timeLocal = time + tz * 60.0 + ((dst) ? 60.0 : 0.0);
        while ((timeLocal < 0.0) || (timeLocal >= 1440.0))
        {
            double incr = ((timeLocal < 0) ? 1 : -1);
            timeLocal += (incr * 1440.0);
            julianday -= incr;
        }
        return julianday;
    }

    // void calculate()
    // {
    //     double jday = getJD(DateTime.Now);
    //     double tl = getTimeLocal();
    //     double tz = -5;
    //     double dst = true;
    //     double total = jday + tl / 1440.0f - tz / 24.0f;
    //     double T = calcTimeJulianCent(total);
    //     double lat = 51;
    //     double lng = 47;
    //     calcAzEl(1, T, tl, lat, lng, tz);
    //     calcSolNoon(jday, lng, tz, dst);
    //     double rise = calcSunriseSet(1, jday, lat, lng, tz, dst);
    //     double set = calcSunriseSet(0, jday, lat, lng, tz, dst);
    //     //alert("JD " + jday + "  " + rise + "  " + set + "  ")
    // }
}
