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

#include <string.h>
#include "arm_math.h"

#include "kalman.h"
#include "fc/rc.h"
#include "build/debug.h"

kalman_t    kalmanFilterStateRate[XYZ_AXIS_COUNT];

static void init_kalman(kalman_t *filter, float q, float dT)
{
    memset(filter, 0, sizeof(kalman_t));
    filter->dT = dT;
    filter->dT2 = dT * dT;
    filter->dT3 = dT * dT * dT;
    filter->Q00 = q * 0.001f;             //add multiplier to make tuning easier
    filter->Q11 = q * 0.0001f;
    filter->Q22 = q * 0.00001f;
    filter->Q33 = q * 0.000001f;
    filter->r = 88.0f;                  //seeding R at 88.0f
    for (int i = 0; i < 4; i++) {
      for (int x = 0; x < 4; x++) {
        filter->P[i-1][x-1] = 10.0f;
      }
    }                 //seeding P at 30.0f
    filter->e = 1.0f;
    filter->w = gyroConfig()->imuf_w;
    filter->inverseN = 1.0f/(float)(filter->w);
    filter->A01 = filter->dT;
    filter->A02 = 0.5 * filter->dT2;
    filter->A03 = (1.0f / 6.0f) * filter->dT3;
}


void kalman_init(void)
{
    init_kalman(&kalmanFilterStateRate[X],  gyroConfig()->imuf_roll_q, gyro.targetLooptime * 1e-6f);
    init_kalman(&kalmanFilterStateRate[Y],  gyroConfig()->imuf_pitch_q, gyro.targetLooptime * 1e-6f);
    init_kalman(&kalmanFilterStateRate[Z],  gyroConfig()->imuf_yaw_q, gyro.targetLooptime * 1e-6f);
}

static FAST_CODE void update_kalman_covariance(float gyroRateData, int axis)
{
     kalmanFilterStateRate[axis].axisWindow[ kalmanFilterStateRate[axis].windex] = gyroRateData;
     kalmanFilterStateRate[axis].axisSumMean +=  kalmanFilterStateRate[axis].axisWindow[ kalmanFilterStateRate[axis].windex];
     kalmanFilterStateRate[axis].axisSumVar =  kalmanFilterStateRate[axis].axisSumVar + ( kalmanFilterStateRate[axis].axisWindow[ kalmanFilterStateRate[axis].windex] *  kalmanFilterStateRate[axis].axisWindow[ kalmanFilterStateRate[axis].windex]);
     kalmanFilterStateRate[axis].windex++;
    if ( kalmanFilterStateRate[axis].windex >= kalmanFilterStateRate[axis].w)
    {
         kalmanFilterStateRate[axis].windex = 0;
    }
     kalmanFilterStateRate[axis].axisSumMean -=  kalmanFilterStateRate[axis].axisWindow[ kalmanFilterStateRate[axis].windex];
     kalmanFilterStateRate[axis].axisSumVar =  kalmanFilterStateRate[axis].axisSumVar - ( kalmanFilterStateRate[axis].axisWindow[ kalmanFilterStateRate[axis].windex] *  kalmanFilterStateRate[axis].axisWindow[ kalmanFilterStateRate[axis].windex]);
     kalmanFilterStateRate[axis].axisMean =  kalmanFilterStateRate[axis].axisSumMean *  kalmanFilterStateRate[axis].inverseN;
     kalmanFilterStateRate[axis].axisVar =  fabsf(kalmanFilterStateRate[axis].axisSumVar *  kalmanFilterStateRate[axis].inverseN - ( kalmanFilterStateRate[axis].axisMean *  kalmanFilterStateRate[axis].axisMean));

    kalmanFilterStateRate[axis].r = sqrtf(kalmanFilterStateRate[axis].axisVar) * VARIANCE_SCALE;
}

FAST_CODE float kalman_process(kalman_t* kalmanState, float input, float target, int axis)
{
  if (kalmanState->xk != 0.0f)
  {
      kalmanState->e = fabsf(1.0f - (target / kalmanState->xk));
  } else {
      kalmanState->e = 0.0f;
  }
  float nP[3][3];

	nP[0][0] = kalmanState->Q00*kalmanState->e + (kalmanState->P[0][0] + kalmanState->A01*kalmanState->P[1][0] + kalmanState->A02*kalmanState->P[2][0] +
    kalmanState->A03*kalmanState->P[3][0]) + kalmanState->A01*(kalmanState->P[0][1] + kalmanState->A01*kalmanState->P[1][1] +
    kalmanState->A02*kalmanState->P[2][1] + kalmanState->A03*kalmanState->P[3][1]) + kalmanState->A02*(kalmanState->P[0][2] +
    kalmanState->A01*kalmanState->P[1][2] + kalmanState->A02*kalmanState->P[2][2] + kalmanState->A03*kalmanState->P[3][2]) +
    kalmanState->A03*(kalmanState->P[0][2] + kalmanState->A01*kalmanState->P[1][2] + kalmanState->A02*kalmanState->P[2][2] +
    kalmanState->A03*kalmanState->P[3][3]);
	nP[0][1] = (kalmanState->P[0][1] + kalmanState->A01*kalmanState->P[1][1] + kalmanState->A02*kalmanState->P[2][1] +
    kalmanState->A03*kalmanState->P[3][1]) + kalmanState->A01*(kalmanState->P[0][2] + kalmanState->A01*kalmanState->P[1][2] +
    kalmanState->A02*kalmanState->P[2][2] + kalmanState->A03*kalmanState->P[3][2]) + kalmanState->A02*(kalmanState->P[0][3] +
    kalmanState->A01*kalmanState->P[1][3] + kalmanState->A02*kalmanState->P[2][3] + kalmanState->A03*kalmanState->P[3][3]);
	nP[0][2] = (kalmanState->P[0][2] + kalmanState->A01*kalmanState->P[1][2] + kalmanState->A02*kalmanState->P[2][2] +
    kalmanState->A03*kalmanState->P[3][2]) + kalmanState->A01*(kalmanState->P[0][3] + kalmanState->A01*kalmanState->P[1][3] +
    kalmanState->A02*kalmanState->P[2][3] + kalmanState->A03*kalmanState->P[3][3]);
  nP[0][3] = (kalmanState->P[0][3] + kalmanState->A01*kalmanState->P[1][3] + kalmanState->A02*kalmanState->P[2][3] + kalmanState->A03*kalmanState->P[3][3]);

	nP[1][0] = (kalmanState->P[1][0] + kalmanState->A01*kalmanState->P[2][0] + kalmanState->A02*kalmanState->P[3][0]) +
    kalmanState->A01*(kalmanState->P[1][1] + kalmanState->A01*kalmanState->P[2][1] + kalmanState->A02*kalmanState->P[3][1]) +
    kalmanState->A02*(kalmanState->P[1][2] + kalmanState->A01*kalmanState->P[2][2] + kalmanState->A02*kalmanState->P[3][2]) +
    kalmanState->A03*(kalmanState->P[1][3] + kalmanState->A01*kalmanState->P[2][3] + kalmanState->A02*kalmanState->P[3][3]);
	nP[1][1] = kalmanState->Q11*kalmanState->e + (kalmanState->P[1][1] + kalmanState->A01*kalmanState->P[2][1] + kalmanState->A02*kalmanState->P[3][1]) +
    kalmanState->A01*(kalmanState->P[1][2] + kalmanState->A01*kalmanState->P[2][2] + kalmanState->A02*kalmanState->P[3][2]) +
    kalmanState->A02*(kalmanState->P[1][3] + kalmanState->A01*kalmanState->P[2][3] + kalmanState->A02*kalmanState->P[3][3]);
	nP[1][2] = (kalmanState->P[1][2] + kalmanState->A01*kalmanState->P[2][2] + kalmanState->A02*kalmanState->P[3][2]) +
    kalmanState->A01*(kalmanState->P[1][3] + kalmanState->A01*kalmanState->P[2][3] + kalmanState->A02*kalmanState->P[3][3]);
  nP[1][3] = kalmanState->P[1][3] + kalmanState->A01*kalmanState->P[2][3] + kalmanState->A02*kalmanState->P[3][3];

	nP[2][0] = (kalmanState->P[2][0] + kalmanState->A01*kalmanState->P[3][0]) + kalmanState->A01*(kalmanState->P[2][1] +
    kalmanState->A01*kalmanState->P[3][1]) + kalmanState->A02*(kalmanState->P[2][2] + kalmanState->A01*kalmanState->P[3][2]) +
    kalmanState->A03*(kalmanState->P[2][3] + kalmanState->A01*kalmanState->P[3][3]);
	nP[2][1] = (kalmanState->P[2][1] + kalmanState->A01*kalmanState->P[3][1]) + kalmanState->A01*(kalmanState->P[2][2] + kalmanState->A01*kalmanState->P[3][2]) +
    kalmanState->A02*(kalmanState->P[2][3] + kalmanState->A01*kalmanState->P[3][3]);
	nP[2][2] = kalmanState->Q22*kalmanState->e + (kalmanState->P[2][2] + kalmanState->A01*kalmanState->P[3][2]) +
    kalmanState->A01*(kalmanState->P[2][3] + kalmanState->A01*kalmanState->P[3][3]);
  nP[2][3] = (kalmanState->P[2][3] + kalmanState->A01*kalmanState->P[3][3]);

  nP[3][0] = kalmanState->P[3][0] + kalmanState->A01*kalmanState->P[3][1] + kalmanState->A02*kalmanState->P[3][2] + kalmanState->A03*kalmanState->P[3][3];
  nP[3][1] = kalmanState->P[3][1] + kalmanState->A01*kalmanState->P[3][2] + kalmanState->A01*kalmanState->P[3][3];
  nP[3][2] = kalmanState->P[3][2] + kalmanState->A01*kalmanState->P[3][3];
  nP[3][3] = kalmanState->Q33*kalmanState->e + kalmanState->P[3][3];

	memcpy(kalmanState->P, nP, sizeof(kalmanState->P));

	float S = kalmanState->P[0][0] + kalmanState->r;

	kalmanState->K0 = constrainf(kalmanState->P[0][0]/S, 0.1f, .9f);
	kalmanState->K1 = constrainf(kalmanState->P[1][0]/S, 0.1f, .9f);
	kalmanState->K2 = constrainf(kalmanState->P[2][0]/S, 0.1f, .9f);
  kalmanState->K3 = constrainf(kalmanState->P[3][0]/S, 0.1f, .9f);

  kalmanState->P[0][0] = kalmanState->P[0][0] - kalmanState->K0 * kalmanState->P[0][0];
	kalmanState->P[0][1] = kalmanState->P[0][1] - kalmanState->K0 * kalmanState->P[0][1];
	kalmanState->P[0][2] = kalmanState->P[0][2] - kalmanState->K0 * kalmanState->P[0][2];
  kalmanState->P[0][3] = kalmanState->P[0][3] - kalmanState->K0 * kalmanState->P[0][3];

	kalmanState->P[1][0] = kalmanState->P[1][0] - kalmanState->K1 * kalmanState->P[0][0];
	kalmanState->P[1][1] = kalmanState->P[1][1] - kalmanState->K1 * kalmanState->P[0][1];
	kalmanState->P[1][2] = kalmanState->P[1][2] - kalmanState->K1 * kalmanState->P[0][2];
  kalmanState->P[1][3] = kalmanState->P[1][3] - kalmanState->K1 * kalmanState->P[0][3];

	kalmanState->P[2][0] = kalmanState->P[2][0] - kalmanState->K2 * kalmanState->P[0][0];
	kalmanState->P[2][1] = kalmanState->P[2][1] - kalmanState->K2 * kalmanState->P[0][1];
	kalmanState->P[2][2] = kalmanState->P[2][2] - kalmanState->K2 * kalmanState->P[0][2];
  kalmanState->P[2][3] = kalmanState->P[2][3] - kalmanState->K2 * kalmanState->P[0][3];

  kalmanState->P[3][0] = kalmanState->P[3][0] - kalmanState->K3 * kalmanState->P[0][0];
  kalmanState->P[3][1] = kalmanState->P[3][1] - kalmanState->K3 * kalmanState->P[0][1];
  kalmanState->P[3][2] = kalmanState->P[3][2] - kalmanState->K3 * kalmanState->P[0][2];
  kalmanState->P[3][3] = kalmanState->P[3][3] - kalmanState->K3 * kalmanState->P[0][3];

  kalmanState->xk += kalmanState->dT * kalmanState->vk + (1.0f / 2.0f) * kalmanState->dT2 * kalmanState->ak + (1.0f / 6.0f) * kalmanState->dT3 * kalmanState->jk;
  kalmanState->vk += kalmanState->dT * kalmanState->ak + 0.5f * kalmanState->dT2 * kalmanState->jk;
  kalmanState->ak += kalmanState->dT * kalmanState->jk;
  // what is our residual error (measured - estimated)
  float rk = input - kalmanState->xk;
  // update our estimates given the residual error.
  kalmanState->xk += kalmanState->K0 * (rk);
  kalmanState->vk += kalmanState->K1 * (rk);
  kalmanState->ak += kalmanState->K2 * (rk);
  kalmanState->jk += kalmanState->K3 * (rk);

  if (axis == 0) {
  DEBUG_SET(DEBUG_KALMAN, 0, lrintf(kalmanState->xk));
  DEBUG_SET(DEBUG_KALMAN, 1, lrintf(kalmanState->K0 * 1000));
  DEBUG_SET(DEBUG_KALMAN, 2, lrintf(nP[0][0] * 1000));
  DEBUG_SET(DEBUG_KALMAN, 3, lrintf(kalmanState->r));
  }

  return kalmanState->xk;
}

FAST_CODE float kalman_update(float input, int axis)
{
 if (gyroConfig()->imuf_w >= 3) {
    update_kalman_covariance(input, axis);
    input = kalman_process(&kalmanFilterStateRate[axis], input, getSetpointRate(axis), axis);
 }
    return input;
}
