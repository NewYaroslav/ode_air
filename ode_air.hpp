#ifndef ODE_AIR_HPP_INCLUDED
#define ODE_AIR_HPP_INCLUDED

#include <ode/ode.h>

/** \brief Применить к телу сопротивление воздуха
 * Определение характерной площади зависит от формы тела:
 * в простейшем случае (шар) — площадь поперечного сечения
 * для крыльев и оперения — площадь крыла/оперения в плане
 * для пропеллеров и несущих винтов вертолётов — либо площадь лопастей, либо ометаемая площадь винта
 * для подводных объектов обтекаемой формы — площадь смачиваемой поверхности
 * для продолговатых тел вращения, ориентированных вдоль потока (фюзеляж, оболочка дирижабля) — приведённая волюметрическая площадь, равная V2/3, где V — объём тела
 * Примеры Cxo КСФ для некоторых форм:
 * Сфера 0.47
 * Конус 2:1 (острием к потоку) 0.5
 * Куб (поверхностью к потоку) 1.05
 * Цилиндр (длина равна двум диаметрам, торцом к потоку) 0.82
 * Вытянутое каплевидное тело 0,04
 * \param body тело
 * \param Cxo безразмерный аэродинамический коэффициент сопротивления
 * \param S характерной площадь
 * \param p плотность воздуха (1.22 при высоте 0м над уровнем моря)
 */
void dBodyCalcAerodynamicDrag (dBodyID body, dReal Cxo, dReal S, dReal p = 1.22d);

/** \brief Применить к телу сопротивление воздуха (для сферического тела)
 * \param body тело
 * \param Cxo безразмерный аэродинамический коэффициент сопротивления
 * \param r радиус сферы
 * \param p плотность воздуха (1.22 при высоте 0м над уровнем моря)
 */
void dBodyCalcAerodynamicDragSphere (dBodyID body, dReal r, dReal p = 1.22d);

/** \brief Применить к телу силу ветра
 * \param body тело
 * \param Cxo безразмерный аэродинамический коэффициент сопротивления
 * \param S характерной площадь
 * \param p плотность воздуха
 */
void dBodyCalcWindStrength (dBodyID body, const dReal velWind[3], const dReal Cxo, const dReal S, const dReal p = 1.22d);

void dBodyCalcWindStrength3D (dBodyID body, const dReal velWind[3], const dReal Cxo[3], const dReal S[3], const dReal p = 1.22d);

/** \brief Применить к телу силу ветра (для сферического тела)
 * \param body тело
 * \param Cxo безразмерный аэродинамический коэффициент сопротивления
 * \param S характерной площадь
 * \param p плотность воздуха
 */
void dBodyCalcWindStrengthSphere (dBodyID body, const dReal velWind[3], const dReal r, const dReal p = 1.22d);

dReal getAirMolarMass(const dReal P, const dReal Pv);

/** \brief Получить ускорение свободного падения в зависимости от широты
 * \param psy широта (0 - 90)
 * \param h высота над уровнем моря
 * \return ускорение свободного падения
 */
dReal getAccelerationGravityLatitude (const dReal psy, const dReal h = 0.0d);

/** \brief Получить плотность воздуха
 * \param T абсолютная температура (в Кельвинах)
 * \param P абсолютное давление (в Паскалях)
 * \param RH относительная влажность (0 - 100)
 * \return плотность воздуха (кг/м3)
 */
dReal getAirDensity (const dReal T, const dReal P, const dReal RH);

/** \brief Давление от высоты над уровнем моря
 * Зависимость давления газа от высоты (барометрическая формула)
 * \param h высота
 * \param T абсолютная температура (в Кельвинах)
 * \param P0 давление на высоте уровня моря (в Паскалях)
 * \param g ускорение свободного падения (м/сс)
 * \return давление
 */
dReal getAirPressureFromAltitude (const dReal h, const dReal T, const dReal P0, const dReal RH, const dReal g = 9.81d);

class OdeAir {
    private:
    dReal T, RH, P0, latitude;
    dReal g, density0;
    public:

    OdeAir() {};

    OdeAir(dReal T, dReal RH, dReal P0, dReal latitude);

    /** \brief Получить ускорение свободного падения
     * \return ускорение свободного падения
     */
    inline dReal getAccelerationGravity() {return g;};

    /** \brief Получить давление
     * \param h высота над уровнем моря
     * \return давление
     */
    inline dReal getPressure(const dReal h) {return getAirPressureFromAltitude(h, T, P0, RH, g);};

    /** \brief Получить плотность воздуха
     * \param h высота над уровнем моря
     * \return плотность воздуха
     */
    inline dReal getDensity(const dReal h) {return getAirDensity(T, getAirPressureFromAltitude(h, T, P0, RH, g), RH);};

};

#endif // ODE_AIR_RESISTANCE_HPP_INCLUDED
