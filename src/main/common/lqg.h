/*
 ******************************************************************************
 * @addtogroup Libraries Libraries
 * @{
 * @addtogroup FlightMath math support libraries
 * @{
 *
 * @file       lqg.h
 * @author     dRonin, http://dronin.org, Copyright (C) 2018
 * @brief      LQG Control algorithm
 *
 * @see        The GNU Public License (GPL) Version 3
 *
 *****************************************************************************/
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see <http://www.gnu.org/licenses/>
 */

#pragma once

#define LQG_SOLVER_FAILED -1
#define LQG_SOLVER_RUNNING 0
#define LQG_SOLVER_DONE 1

typedef struct lqr_state_s {
	int16_t solver_iterations;

	float A[2][2];
	float B[2];
	float K[2];
	float P[2][2];
	float R;
	float Q[2][2];
	float u;

	float beta;
	float tau;
} lqr_state_t;

typedef struct rtkf_state_s {
	int solver_iterations;

	float A[3][3];
	float B[3];
	float K[3];
	float P[3][3];
	float Q[3][3];
	float R;
	float X[3];

	float biaslim;
} rtkf_state_t;

extern void rtkf_create(rtkf_state_t *rtkf, float beta, float tau, float Ts, float R, float Q1, float Q2, float Q3, float biaslim);
extern void rtkf_stabilize_covariance(rtkf_state_t *rtkf);
// extern int rtkf_solver_status(rtkf_state_t *rtkf);

extern void lqr_create(lqr_state_t *lqr, float beta, float tau, float Ts, float q1, float q2, float r);
extern void lqr_stabilize_covariance(lqr_state_t *lqr);
// extern int lqr_solver_status(lqr_state_t *lqr);

extern void lqr_update(lqr_state_t *lqr, float q1, float q2, float r);
extern void lqr_get_gains(lqr_state_t *lqr, float K[2]);

// extern int lqg_solver_status(lqg_t lqg);
// extern lqr_t lqg_get_lqr(lqg_t lqg);
// extern rtkf_t lqg_get_rtkf(lqg_t lqg);

extern void lqg_run_covariance(rtkf_state_t *rtkf, lqr_state_t *lqg);

extern float lqg_controller(rtkf_state_t *rtkf, lqr_state_t *lqg, float signal, float setpoint);

// extern void lqg_set_x0(lqg_t lqq, float x0);
// void lqg_get_rtkf_state(lqg_t lqg, float *rate, float *torque, float *bias);

/**
 * @}
 * @}
 */
