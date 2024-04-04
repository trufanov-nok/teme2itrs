#ifndef TEME2ITEF_H
#define TEME2ITEF_H

#include <Eigen/Dense>

/*
 * Конвертация векторов RV из СК TEME в СК ITRS
 * через СК GCRS (тоже инарционная как и TEME)
 * при помощи библиотеки SOFA (её opensource вариант - ERFA) IAU (Международный астрономический соююз).
 */


void teme2itrs(const Eigen::Vector3d& teme_r, // радиус-вектор в TEME
               const Eigen::Vector3d& teme_v, // скорость в TEME
               double unixtime,               // время в unixtime
               /*
                * https://en.wikipedia.org/wiki/DUT1
                * Источником может быть:
                * - https://celestrak.org/SpaceData/
                * - https://datacenter.iers.org/eop.php - Международная служба вращения Земли (IERS)
                * - http://pvz.vniiftri.ru/bull_A_Q.php - Центр обработки и анализа данных о параметрах вращения Земли Главного метрологического центра ГСВЧ (ЦОАД ПВЗ ГМЦ ГСВЧ)
                */
               double dUT1,
               Eigen::Vector3d& itrs_r,       // [OUT] радиус-вектор в ITRS
               Eigen::Vector3d& itrs_v);      // [OUT] радиус-вектор в ITRS

#endif // TEME2ITEF_H
