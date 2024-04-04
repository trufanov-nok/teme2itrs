This is an example of TEME to ITRS conversion (via GCRS) written with C++ using IAU's [ERFA](https://github.com/liberfa/erfa) library.

**TEME**: the [ECI frame](https://en.wikipedia.org/wiki/Earth-centered_inertial) used for the NORAD [two-line elements](https://en.wikipedia.org/wiki/Two-line_element_set) is sometimes called _true equator, mean equinox_ (TEME) although it does not use the conventional mean equinox.   
**GCRS**: ([Geocentric Celestial Reference System](https://en.wikipedia.org/wiki/Barycentric_and_geocentric_celestial_reference_systems)) is an astronomical coordinate system with the earth as its origin, which is relative to the earth. It's also [ECI](https://en.wikipedia.org/wiki/Earth-centered_inertial) as TEME.   
**ITRS**: ([International Terrestrial Reference System](https://en.wikipedia.org/wiki/International_Terrestrial_Reference_System_and_Frame) an [ECEF](https://en.wikipedia.org/wiki/Earth-centered,_Earth-fixed_coordinate_system) coordinate system.   
**ERFA**: ([Essential Routines for Fundamental Astronomy](https://github.com/liberfa/erfa)) is a C library containing key algorithms for astronomy, and is based on the [SOFA library](http://www.iausofa.org/) published by the International Astronomical Union (IAU).

# How to build

```bash
sudo apt install libeigen3-dev liberfa-dev

git clone https://github.com/trufanov-nok/teme2itrs.git
cd teme2itrs
mkdir build
cd build
cmake ..
make
```

# CLI
teme2itrs Rx Ry Rz Vx Vy Vz unixtime dUT1   

Rx - x component of radius vector in TEME (km)   
Ry - y component of radius vector in TEME (km)   
Rz - z component of radius vector in TEME (km)   
Vx - x component of velocity vector in TEME (km/sec)   
Vy - y component of velocity vector in TEME (km/sec)   
Vz - z component of velocity vector in TEME (km/sec)   
unixtime - time since 1970/01/01 (sec)   
[dUT1](https://en.wikipedia.org/wiki/DUT1) - optional time correction (sec) (see https://en.wikipedia.org/wiki/DUT1)   

# Test
Use skyfield or astropy python packages [to test](https://github.com/trufanov-nok/teme2itrs/wiki).

Example calls:    
```bash
./teme2itrs 3469.9479844480247 -2690.388430365502 5175.8319246510355 5.810229142098143 4.802261184575433 -1.3882803330121878  1575909509.3634412 -0.1724485450587226

TEME R: 3469.947984448025, -2690.388430365502, 5175.831924651035 km
TEME V: 5.810229142098143, 4.802261184575433, -1.388280333012188 km\sec
unixtime: 1575909509.363441 sec
dUT1: -0.1724485450587226 sec

ITRS R: 4370.210055865812, -424.2558056512487, 5175.831924651035 km
ITRS V: 2.321254144507617, 6.842860358413781, -1.388280333012187 km\sec
```

Astropy for TLE   
1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991   
2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482   

```bash
TEME R: (3469.94808421, -2690.38834791, 5175.83190082) km 
TEME V: (5.81022907, 4.80226124, -1.38828045) km / s 
unixtime: 1575909509.3634412 
dUT1: -0.1724485450587226

ITRS R: (4370.2127611, -424.26248152, 5175.82909345) km 
ITRS V: (2.32125334, 6.84286215, -1.38827265) km / s
```

Skyfield for same TLE:
```bash
unixtime: 1575909509.3634338 
dUT1: -0.1724535139325667 
ITRS R: [4370.21007939 -424.25573126 5175.83191099] 
ITRS V: [ 2.32125409  6.84286036 -1.3882804 ]
```
