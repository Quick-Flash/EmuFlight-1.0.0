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

#include "common/axis.h"
#include "common/maths.h"

#include "pg/pg.h"

typedef struct boardAlignment_s {
    int32_t rollDegrees;
    int32_t pitchDegrees;
    int32_t yawDegrees;
} boardAlignment_t;

PG_DECLARE(boardAlignment_t, boardAlignment);

void alignSensorViaMatrix(float *dest, fp_rotationMatrix_t* rotationMatrix);
void alignSensorViaRotation(float *dest, uint8_t rotation);

void initBoardAlignment(const boardAlignment_t *boardAlignment);
bool isBoardAlignmentStandard(const boardAlignment_t *boardAlignment);