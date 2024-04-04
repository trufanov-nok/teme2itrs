#include "teme2itrs.h"

#include <erfa.h>

/*
 * Получаем время UT1 и TT через преобразование UTC в TT в TAI (атомарные часы) в TT.
 * На входе требуется время и dUT1 - delta_utc_ut1, параметр получаемый из EOP (earth orientation parameters)
 *
 * На выходе UT1 и TT представлены парой чисел, для сохранения точности
 * Источником может быть:
 * - https://celestrak.org/SpaceData/
 * - https://datacenter.iers.org/eop.php - Международная служба вращения Земли (IERS)
 * - http://pvz.vniiftri.ru/bull_A_Q.php - Центр обработки и анализа данных о параметрах вращения Земли Главного метрологического центра ГСВЧ (ЦОАД ПВЗ ГМЦ ГСВЧ)
 */

void time_to_UT1_and_TT(int y, int m, int d,
                        int hh, int mm, double sec,
                        double dUT1,
                        double& jdut1, double& jdut2,
                        double& jdtt1, double& jdtt2)
{
    double jdutc1, jdutc2;
    eraDtf2d("UTC", y, m, d, hh, mm, sec, &jdutc1, &jdutc2);
    eraUtcut1(jdutc1, jdutc2, dUT1, &jdut1, &jdut2);

    double tai1, tai2;
    eraUtctai(jdutc1, jdutc2, &tai1, &tai2);
    eraTaitt(tai1, tai2, &jdtt1, &jdtt2);
}

/*
 * Матрица поворота TEME в конкретный момент времени
 */

Eigen::Matrix3d TEME_rotation_at(double jdut1, double jdut2,
                                 double jdtt1, double jdtt2,
                                 double dUT1)
{
    /* Classical nutation x precession x bias matrix, IAU 2000A. */
    double rnpb[3][3];
    eraPnm06a(jdut1, jdut2, rnpb);


    /*
     * Выбрана модель прецесси и нутации Земли IAU 2006,
     * т.к. именно ее коэффициенты я нашел в python пакете skyfield,
     * с которым конверсию и сравнивал.
     */

    /* Greenwich apparent sidereal time. */
    const double gast = eraGst06(jdut1, jdut2, jdtt1, jdtt2, rnpb);
    const Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> NPB(*rnpb);


    const double theta = eraGmst82(jdut1, jdut2);

    const double angle = theta - gast;
    const double cos_angle = cos(angle);
    const double sin_angle = sin(angle);
    return Eigen::Matrix3d({
               {cos_angle, -sin_angle, 0.},
               {sin_angle,  cos_angle, 0.},
               {       0.,         0., 1.}
           }) * NPB;
}

void teme2gcrs(const Eigen::Vector3d& teme_r,
               const Eigen::Vector3d& teme_v,
               double jdut1, double jdut2,
               double jdtt1, double jdtt2,
               double dUT1,
               Eigen::Vector3d& gcrs_r,
               Eigen::Vector3d& gcrs_v)
{
    const Eigen::Matrix3d teme2gcrs =
        TEME_rotation_at(jdut1, jdut2, jdtt1, jdtt2, dUT1).transpose();
    gcrs_r = teme2gcrs * teme_r;
    gcrs_v = teme2gcrs * teme_v;
}

void gcrs2teme(const Eigen::Vector3d& gcrs_r,
               const Eigen::Vector3d& gcrs_v,
               double jdut1, double jdut2,
               double jdtt1, double jdtt2,
               double dUT1,
               Eigen::Vector3d& teme_r,
               Eigen::Vector3d& teme_v)
{
    const Eigen::Matrix3d gcrs2teme =
        TEME_rotation_at(jdut1, jdut2, jdtt1, jdtt2, dUT1);
    teme_r = gcrs2teme * gcrs_r;
    teme_v = gcrs2teme * gcrs_v;
}


/*
 *  Переводим GCRS в ITRF
 */

Eigen::Matrix3d
GCRStoGreenwich(double jdut1, double jdut2, double jdtt1, double jdtt2)
{
    /* Classical nutation x precession x bias matrix, IAU 2000A. */
    double rnpb[3][3];
    eraPnm06a(jdut1, jdut2, rnpb);

    /* Greenwich apparent sidereal time. */
    const double gast = eraGst06(jdut1, jdut2, jdtt1, jdtt2, rnpb);

    const Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> NPB(*rnpb);

    const double cos_gast = cos(gast);
    const double sin_gast = sin(gast);
    const Eigen::Matrix3d ST_transp = Eigen::Matrix3d({
        { cos_gast, sin_gast, 0.},
        {-sin_gast, cos_gast, 0.},
        {0.,              0., 1.}
    });

    return  ST_transp*NPB;
}

// ω_з - угловая скорость вращения Земли (документ ПЗ-90.11) [рад/с]
const double OMEGA_EARTH = 7.2921150e-5;
const Eigen::Matrix3d _ITRS_angvel_matrix ({
                                           { 0.0,  OMEGA_EARTH, 0.0},
                                           { -OMEGA_EARTH, 0.0, 0.0},
                                           { 0.0,          0.0, 0.0},
                                           });

void gcrs2itrf(const Eigen::Vector3d& gcrs_r,
               const Eigen::Vector3d& gcrs_v,
               double jdut1, double jdut2, double jdtt1, double jdtt2,
               Eigen::Vector3d& itrf_r,
               Eigen::Vector3d& itrf_v)
{
    const Eigen::Matrix3d gcrs2itrf = GCRStoGreenwich(jdut1, jdut2, jdtt1, jdtt2);
    itrf_r = gcrs2itrf * gcrs_r;
    itrf_v = gcrs2itrf * gcrs_v + _ITRS_angvel_matrix * itrf_r;
}


void teme2itrs(const Eigen::Vector3d& teme_r,
               const Eigen::Vector3d& teme_v,
               double unixtime, double dUT1,
               Eigen::Vector3d& itrf_r,
               Eigen::Vector3d& itrf_v)
{
    time_t unix_time_ = unixtime;
    struct tm time_;
#ifndef _WIN32
    gmtime_r((time_t *)&unix_time_, &time_);
#else
    gmtime_s(&time_, (time_t *)&unix_time_);
#endif

    const int y = time_.tm_year + 1900;
    const int m = time_.tm_mon + 1;
    const int d = time_.tm_mday;

    const double unix_time_at_midight = ((int)unixtime / 86400/*SECONDS_IN_DAY_INT*/)*86400;
    const double time = unixtime - unix_time_at_midight; // seconds_since_midnight

    const int hh     = (int) time / 3600 /*sec in hour*/;
    const double remainingSec = fmod(time, 3600.0);
    const int mm     = (int) remainingSec / 60;
    const double sec = fmod(remainingSec, 60.);



    double jdut1, jdut2, jdtt1, jdtt2;
    time_to_UT1_and_TT(y, m, d, hh, mm, sec, dUT1,
                       jdut1, jdut2, jdtt1, jdtt2);

    Eigen::Vector3d r;
    Eigen::Vector3d v;
    teme2gcrs(teme_r, teme_v, jdut1, jdut2, jdtt1, jdtt2, dUT1, r, v);
    gcrs2itrf(r, v, jdut1, jdut2, jdtt1, jdtt2, itrf_r, itrf_v);
}
