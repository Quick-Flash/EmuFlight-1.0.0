/*
 * This file is part of Cleanflight and Betaflight and EmuFlight.
 *
 * Cleanflight and Betaflight and EmuFlight are free software. You can redistribute
 * this software and/or modify this software under the terms of the
 * GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Cleanflight and Betaflight and EmuFlight are distributed in the hope that they
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.
 *
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "sensors/gyro.h"
#include "filter.h"

#define MAX_KALMAN_WINDOW_SIZE 300

#define VARIANCE_SCALE 0.67f


typedef struct kalman
{
    float Q00, Q11, Q22, Q33;     //process noise covariance
    float r;     //measurement noise covariance
    float P[3][3];     //estimation error covariance matrix
    float A01, A02, A03;
    float K0, K1, K2, K3;     //kalman gain
    float e;
    float ak, vk, xk, jk;
    float dT, dT2, dT3;
    float axisVar;
    uint16_t windex;
    float axisWindow[MAX_KALMAN_WINDOW_SIZE];
    float axisSumMean;
    float axisMean;
    float axisSumVar;
    float inverseN;
    uint16_t w;
} kalman_t;

extern void kalman_init(void);
extern float kalman_update(float input, int axis);
