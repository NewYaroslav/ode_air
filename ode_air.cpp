#include "ode_air.hpp"

void dBodyCalcAerodynamicDrag (dBodyID body, dReal Cxo, dReal S, dReal p) {
    const dReal minVel = 0.0d;
    const dReal* vel = dBodyGetLinearVel (body);
    dReal fAir[3];
    for(int i = 0; i < 3; ++i) {
        if(abs(vel[i]) > minVel) {
            fAir[i] = vel[i] * vel[i] * Cxo * p * S / 2.0d;
            if(vel[i] > 0) fAir[i] = -fAir[i];
        } else fAir[i] = 0.0d;
    } // for
    dBodyAddForce (body, fAir[0], fAir[1], fAir[2]);
}

void dBodyCalcWindStrength (dBodyID body, const dReal velWind[3], const dReal Cxo, const dReal S, const dReal p) {
    const dReal minVel = 0.0d;
    const dReal* vel = dBodyGetLinearVel (body);
    dReal fAir[3];
    for(int i = 0; i < 3; ++i) {
        dReal dvel = vel[i] - velWind[i];
        if(abs(dvel) > minVel) {
            fAir[i] = dvel * dvel * Cxo * p * S / 2.0d;
            if(dvel > 0) fAir[i] = -fAir[i];
        } else fAir[i] = 0.0d;
    } // for
    dBodyAddForce (body, fAir[0], fAir[1], fAir[2]);
}

void dBodyCalcWindStrength3D (dBodyID body, const dReal velWind[3], const dReal Cxo[3], const dReal S[3], const dReal p) {
    const dReal minVel = 0.0d;
    const dReal* vel = dBodyGetLinearVel (body);
    dReal fAir[3];
    for(int i = 0; i < 3; ++i) {
        dReal dvel = vel[i] - velWind[i];
        if(abs(dvel) > minVel) {
            fAir[i] = dvel * dvel * Cxo[i] * p * S[i] / 2.0d;
            if(dvel > 0) fAir[i] = -fAir[i];
        } else fAir[i] = 0.0d;
    } // for
    dBodyAddForce (body, fAir[0], fAir[1], fAir[2]);
}

void dBodyCalcAerodynamicDragSphere (dBodyID body, dReal r, dReal p) {
    const dReal minVel = 0.0;
    const dReal Cxo = 0.47d;
    const dReal PI = 3.1415926535897932384626433832795d;
    dReal S = r * r * PI;
    const dReal* vel = dBodyGetLinearVel (body);
    dReal fAir[3];
    for(int i = 0; i < 3; ++i) {
        if(abs(vel[i]) > minVel) {
            fAir[i] = vel[i] * vel[i] * Cxo * p * S / 2.0d;
            if(vel[i] > 0) fAir[i] = -fAir[i];
        } else fAir[i] = 0.0d;
    } // for
    dBodyAddForce (body, fAir[0], fAir[1], fAir[2]);
}

void dBodyCalcWindStrengthSphere (dBodyID body, const dReal velWind[3], const dReal r, const dReal p) {
    const dReal minVel = 0.0d;
    const dReal Cxo = 0.47d;
    const dReal PI = 3.1415926535897932384626433832795d;
    dReal S = r * r * PI;
    const dReal* vel = dBodyGetLinearVel (body);
    dReal fAir[3];
    for(int i = 0; i < 3; ++i) {
        dReal dvel = vel[i] - velWind[i];
        if(abs(dvel) > minVel) {
            fAir[i] = dvel * dvel * Cxo * p * S / 2.0d;
            if(dvel > 0) fAir[i] = -fAir[i];
        } else fAir[i] = 0.0d;
    } // for
    dBodyAddForce (body, fAir[0], fAir[1], fAir[2]);
}

dReal getAirMolarMass(const dReal P, const dReal Pv) {
    const dReal uv = 0.02896; // молярная масса сухого воздуха
    return uv - 0.010944 * (Pv / P);
}

dReal getAccelerationGravityLatitude (const dReal psy, const dReal h) {
    const dReal PI = 3.1415926535897932384626433832795d;
    dReal rad = psy * PI / 180.0d;
    dReal coeff1 = sin(rad);
    coeff1 = coeff1 * coeff1;
    dReal coeff2 = sin(2 * rad);
    coeff2 = coeff2 * coeff2;
    return 9.780318d * (1.0d + 0.005302d * coeff1 - 0.000006d * coeff2) - 0.000003086d * h;
}

dReal getAirDensity (dReal T, dReal P, dReal RH) {
    dReal Psat = 6.1078 * pow(10, (7.5 * T - 2048.625)/(T - 35.85)); // парциальное давление насыщенного пара (в миллибарах)
    Psat *= 100; // Приводим к ПА
    dReal Pv = RH * Psat; // давление водяного пара
    dReal Pd = P - Pv; // давление сухого воздуха
    const dReal Rv = 461.495; // постоянная для пара (461,495 Дж/кг·К)
    const dReal Rd = 287.058; // газовая постоянная для сухого воздуха (287,058 Дж/кг·К)
    return (Pd / (Rd * T)) + (Pv / (Rv * T));
}

dReal getAirPressureFromAltitude (const dReal h, const dReal T, const dReal P0, const dReal RH, const dReal g) {
    dReal Psat = 6.1078 * pow(10, (7.5 * T - 2048.625)/(T - 35.85)); // парциальное давление насыщенного пара (в миллибарах)
    Psat *= 100; // Приводим к ПА
    dReal Pv = RH * Psat; // давление водяного пара
    const dReal u = getAirMolarMass(P0, Pv); // молярная масса воздуха
    const dReal R = 8.31;  // универсальная газовая постоянная
    return P0 * exp((-u * g * h) / (R * T));
}

OdeAir::OdeAir(dReal T, dReal RH, dReal P0, dReal latitude) {
    OdeAir::T = T;
    OdeAir::RH = RH;
    OdeAir::P0 = P0;
    OdeAir::latitude =latitude;
    g = getAccelerationGravityLatitude(latitude, 0);
    density0 = getAirDensity(T, P0, RH);
}
