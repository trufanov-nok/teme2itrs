#include <iomanip>
#include <iostream>
#include "teme2itrs.h"

using namespace std;

int print_help()
{
    cout <<
        "teme2itrs Rx Ry Rz Vx Vy Vz unixtime dUT1\n"
        "Rx - x component of radius vector in TEME (km)\n"
        "Ry - y component of radius vector in TEME (km)\n"
        "Rz - z component of radius vector in TEME (km)\n"
        "Vx - x component of velocity vector in TEME (km/sec)\n"
        "Vy - y component of velocity vector in TEME (km/sec)\n"
        "Vz - z component of velocity vector in TEME (km/sec)\n"
        "unixtime - time since 1970/01/01 (sec)\n"
        "dUT1 - optional time correction (sec) (see https://en.wikipedia.org/wiki/DUT1)\n";
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 8) {
        return print_help();
    }

    const Eigen::Vector3d teme_r({stod(argv[1]), stod(argv[2]), stod(argv[3])});
    const Eigen::Vector3d teme_v({stod(argv[4]), stod(argv[5]), stod(argv[6])});
    const double unixtime = stod(argv[7]);
    const double dUT1 = (argc > 8) ? stod(argv[8]) : 0.0;

    cout << std::setprecision(16);
    cout << "TEME R: " << teme_r.x() << ", " << teme_r.y() << ", " << teme_r.z() << " km\n"
         << "TEME V: " << teme_v.x() << ", " << teme_v.y() << ", " << teme_v.z()  << " km\\sec\n"
         << "unixtime: " << unixtime << " sec\n"
         << "dUT1: " << dUT1 << " sec\n\n";

    Eigen::Vector3d itrs_r;
    Eigen::Vector3d itrs_v;
    teme2itrs(teme_r, teme_v, unixtime, dUT1, itrs_r, itrs_v);
    cout << "ITRS R: " << itrs_r.x() << ", " << itrs_r.y() << ", " << itrs_r.z()  << " km\n"
         << "ITRS V: " << itrs_v.x() << ", " << itrs_v.y() << ", " << itrs_v.z() << " km\\sec\n";

    return 0;
}
