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

#include <stdint.h>

#include "pg/pg.h"

typedef enum {
    RATES_TYPE_BETAFLIGHT = 0,
    RATES_TYPE_RACEFLIGHT,
    RATES_TYPE_KISS,
    RATES_TYPE_ACTUAL,
} ratesType_e;

typedef enum {
    THROTTLE_LIMIT_TYPE_OFF = 0,
    THROTTLE_LIMIT_TYPE_SCALE,
    THROTTLE_LIMIT_TYPE_CLIP,
    THROTTLE_LIMIT_TYPE_COUNT   // must be the last entry
} throttleLimitType_e;

typedef struct rateDynamics_s {
    uint8_t rateSensCenter;
    uint8_t rateSensEnd;
    uint8_t rateCorrectionCenter;
    uint8_t rateCorrectionEnd;
    uint8_t rateWeightCenter;
    uint8_t rateWeightEnd;
} rateDynamics_t;

#define MAX_RATE_PROFILE_NAME_LENGTH 8u

typedef struct controlRateConfig_s {
    uint8_t thrMid8;
    uint8_t thrExpo8;
    uint8_t rates_type;
    uint8_t rcRates[3];
    uint8_t rcExpo[3];
    uint8_t rates[3];
    rateDynamics_t rateDynamics;
    uint8_t throttle_limit_type;            // Sets the throttle limiting type - off, scale or clip
    uint8_t throttle_limit_percent;         // Sets the maximum pilot commanded throttle limit
    uint8_t addRollToYawRc;
    uint8_t addYawToRollRc;
    uint8_t rollPitchMagExpo;       // expo applied when pitch and roll are both high
    char profileName[MAX_RATE_PROFILE_NAME_LENGTH + 1]; // Descriptive name for rate profile
} controlRateConfig_t;

PG_DECLARE_ARRAY(controlRateConfig_t, CONTROL_RATE_PROFILE_COUNT, controlRateProfiles);

extern controlRateConfig_t *currentControlRateProfile;

void loadControlRateProfile(void);
void changeControlRateProfile(uint8_t controlRateProfileIndex);

void copyControlRateProfile(const uint8_t dstControlRateProfileIndex, const uint8_t srcControlRateProfileIndex);
