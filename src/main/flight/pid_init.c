/*
 * This file is part of Cleanflight and Betaflight.
 *
 * Cleanflight and Betaflight are free software. You can redistribute
 * this software and/or modify this software under the terms of the
 * GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Cleanflight and Betaflight are distributed in the hope that they
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.
 *
 * If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "platform.h"

#include "build/build_config.h"
#include "build/debug.h"

#include "common/axis.h"
#include "common/filter.h"
#include "common/maths.h"

#include "drivers/dshot_command.h"

#include "fc/rc_controls.h"
#include "fc/runtime_config.h"

#include "flight/interpolated_setpoint.h"
#include "flight/pid.h"
#include "flight/rpm_filter.h"

#include "sensors/gyro.h"
#include "sensors/sensors.h"

#include "pid_init.h"

#if defined(USE_D_MIN)
#define D_MIN_RANGE_HZ 80    // Biquad lowpass input cutoff to peak D around propwash frequencies
#define D_MIN_LOWPASS_HZ 10  // PT1 lowpass cutoff to smooth the boost effect
#define D_MIN_GAIN_FACTOR 0.00005f
#define D_MIN_SETPOINT_GAIN_FACTOR 0.00005f
#endif

#define ANTI_GRAVITY_THROTTLE_FILTER_CUTOFF 15  // The anti gravity throttle highpass filter cutoff
#define ANTI_GRAVITY_SMOOTH_FILTER_CUTOFF 3  // The anti gravity P smoothing filter cutoff

static void pidSetTargetLooptime(uint32_t pidLooptime)
{
    targetPidLooptime = pidLooptime;
    pidRuntime.dT = targetPidLooptime * 1e-6f;
    pidRuntime.pidFrequency = 1.0f / pidRuntime.dT;
#ifdef USE_DSHOT
    dshotSetPidLoopTime(targetPidLooptime);
#endif
}

void pidInitFilters(const pidProfile_t *pidProfile)
{
    STATIC_ASSERT(FD_YAW == 2, FD_YAW_incorrect); // ensure yaw axis is 2

    if (targetPidLooptime == 0) {
        // no looptime set, so set all the filters to null
        pidRuntime.dtermNotchApplyFn = nullFilterApply;
        pidRuntime.dtermLowpassApplyFn = nullFilterApply;
        pidRuntime.dtermLowpass2ApplyFn = nullFilterApply;
        pidRuntime.ptermYawLowpassApplyFn = nullFilterApply;
        pidRuntime.dtermABGApplyFn = nullFilterApply;
        return;
    }

    const uint32_t pidFrequencyNyquist = pidRuntime.pidFrequency / 2; // No rounding needed

    uint16_t dTermNotchHz;
    if (pidProfile->dterm_notch_hz <= pidFrequencyNyquist) {
        dTermNotchHz = pidProfile->dterm_notch_hz;
    } else {
        if (pidProfile->dterm_notch_cutoff < pidFrequencyNyquist) {
            dTermNotchHz = pidFrequencyNyquist;
        } else {
            dTermNotchHz = 0;
        }
    }

    if (dTermNotchHz != 0 && pidProfile->dterm_notch_cutoff != 0) {
        pidRuntime.dtermNotchApplyFn = (filterApplyFnPtr)biquadFilterApply;
        const float notchQ = filterGetNotchQ(dTermNotchHz, pidProfile->dterm_notch_cutoff);
        for (int axis = FD_ROLL; axis <= FD_YAW; axis++) {
            biquadFilterInit(&pidRuntime.dtermNotch[axis], dTermNotchHz, targetPidLooptime, notchQ, FILTER_NOTCH);
        }
    } else {
        pidRuntime.dtermNotchApplyFn = nullFilterApply;
    }

    //1st Dterm Lowpass Filter
    uint16_t dterm_lowpass_hz = pidProfile->dterm_lowpass_hz;

#ifdef USE_DYN_LPF
    if (pidProfile->dyn_lpf_dterm_min_hz) {
        dterm_lowpass_hz = pidProfile->dyn_lpf_dterm_min_hz;
    }
#endif

    if (dterm_lowpass_hz > 0 && dterm_lowpass_hz < pidFrequencyNyquist) {
        switch (pidProfile->dterm_filter_type) {
        case FILTER_PT1:
            pidRuntime.dtermLowpassApplyFn = (filterApplyFnPtr)pt1FilterApply;
            for (int axis = FD_ROLL; axis <= FD_YAW; axis++) {
                pt1FilterInit(&pidRuntime.dtermLowpass[axis].pt1Filter, pt1FilterGain(dterm_lowpass_hz, pidRuntime.dT));
            }
            break;
        case FILTER_BIQUAD:
#ifdef USE_DYN_LPF
            pidRuntime.dtermLowpassApplyFn = (filterApplyFnPtr)biquadFilterApplyDF1;
#else
            pidRuntime.dtermLowpassApplyFn = (filterApplyFnPtr)biquadFilterApply;
#endif
            for (int axis = FD_ROLL; axis <= FD_YAW; axis++) {
                biquadFilterInitLPF(&pidRuntime.dtermLowpass[axis].biquadFilter, dterm_lowpass_hz, targetPidLooptime);
            }
            break;
        default:
            pidRuntime.dtermLowpassApplyFn = nullFilterApply;
            break;
        }
    } else {
        pidRuntime.dtermLowpassApplyFn = nullFilterApply;
    }

    //2nd Dterm Lowpass Filter
    if (pidProfile->dterm_lowpass2_hz == 0 || pidProfile->dterm_lowpass2_hz > pidFrequencyNyquist) {
        pidRuntime.dtermLowpass2ApplyFn = nullFilterApply;
    } else {
        switch (pidProfile->dterm_filter2_type) {
        case FILTER_PT1:
            pidRuntime.dtermLowpass2ApplyFn = (filterApplyFnPtr)pt1FilterApply;
            for (int axis = FD_ROLL; axis <= FD_YAW; axis++) {
                pt1FilterInit(&pidRuntime.dtermLowpass2[axis].pt1Filter, pt1FilterGain(pidProfile->dterm_lowpass2_hz, pidRuntime.dT));
            }
            break;
        case FILTER_BIQUAD:
            pidRuntime.dtermLowpass2ApplyFn = (filterApplyFnPtr)biquadFilterApply;
            for (int axis = FD_ROLL; axis <= FD_YAW; axis++) {
                biquadFilterInitLPF(&pidRuntime.dtermLowpass2[axis].biquadFilter, pidProfile->dterm_lowpass2_hz, targetPidLooptime);
            }
            break;
        default:
            pidRuntime.dtermLowpass2ApplyFn = nullFilterApply;
            break;
        }
    }

    if (pidProfile->yaw_lowpass_hz == 0 || pidProfile->yaw_lowpass_hz > pidFrequencyNyquist) {
        pidRuntime.ptermYawLowpassApplyFn = nullFilterApply;
    } else {
        pidRuntime.ptermYawLowpassApplyFn = (filterApplyFnPtr)pt1FilterApply;
        pt1FilterInit(&pidRuntime.ptermYawLowpass, pt1FilterGain(pidProfile->yaw_lowpass_hz, pidRuntime.dT));
    }

    if (pidProfile->dtermAlpha == 0) {
        pidRuntime.dtermABGApplyFn = nullFilterApply;
    } else {
        pidRuntime.dtermABGApplyFn = (filterApplyFnPtr)alphaBetaGammaApply;
        for (int axis = FD_ROLL; axis <= FD_YAW; axis++) {
            ABGInit(&pidRuntime.dtermABG[axis], pidProfile->dtermAlpha, pidRuntime.dT);
        }
    }

#if defined(USE_THROTTLE_BOOST)
    pt1FilterInit(&throttleLpf, pt1FilterGain(pidProfile->throttle_boost_cutoff, pidRuntime.dT));
#endif
#if defined(USE_ITERM_RELAX)
    if (pidRuntime.itermRelaxCutoff || pidRuntime.itermRelaxCutoffYaw) {
        for (int i = 0; i < XYZ_AXIS_COUNT; i++) {
            if (i != FD_YAW) {
                pt1FilterInit(&pidRuntime.windupLpf[i], pt1FilterGain(pidRuntime.itermRelaxCutoff, pidRuntime.dT));
            } else {
                pt1FilterInit(&pidRuntime.windupLpf[i], pt1FilterGain(pidRuntime.itermRelaxCutoffYaw, pidRuntime.dT));
            }
        }
    }
#endif
#if defined(USE_D_MIN)

    // Initialize the filters for all axis even if the d_min[axis] value is 0
    // Otherwise if the pidProfile->d_min_xxx parameters are ever added to
    // in-flight adjustments and transition from 0 to > 0 in flight the feature
    // won't work because the filter wasn't initialized.
    for (int axis = FD_ROLL; axis <= FD_YAW; axis++) {
        biquadFilterInitLPF(&pidRuntime.dMinRange[axis], D_MIN_RANGE_HZ, targetPidLooptime);
        pt1FilterInit(&pidRuntime.dMinLowpass[axis], pt1FilterGain(D_MIN_LOWPASS_HZ, pidRuntime.dT));
     }
#endif

    pt1FilterInit(&pidRuntime.antiGravityThrottleLpf, pt1FilterGain(ANTI_GRAVITY_THROTTLE_FILTER_CUTOFF, pidRuntime.dT));
    pt1FilterInit(&pidRuntime.antiGravitySmoothLpf, pt1FilterGain(ANTI_GRAVITY_SMOOTH_FILTER_CUTOFF, pidRuntime.dT));

    pidRuntime.ffBoostFactor = (float)pidProfile->ff_boost / 10.0f;
}

void pidInit(const pidProfile_t *pidProfile)
{
    pidSetTargetLooptime(gyro.targetLooptime); // Initialize pid looptime
    pidInitFilters(pidProfile);
    pidInitConfig(pidProfile);
#ifdef USE_RPM_FILTER
    rpmFilterInit(rpmFilterConfig());
#endif
}

#ifdef USE_RC_SMOOTHING_FILTER
void pidInitSetpointDerivativeLpf(uint16_t filterCutoff, uint8_t debugAxis, uint8_t filterType)
{
    pidRuntime.rcSmoothingDebugAxis = debugAxis;
    pidRuntime.rcSmoothingFilterType = filterType;
    if ((filterCutoff > 0) && (pidRuntime.rcSmoothingFilterType != RC_SMOOTHING_DERIVATIVE_OFF)) {
        pidRuntime.setpointDerivativeLpfInitialized = true;
        for (int axis = FD_ROLL; axis <= FD_YAW; axis++) {
            switch (pidRuntime.rcSmoothingFilterType) {
                case RC_SMOOTHING_DERIVATIVE_PT1:
                    pt1FilterInit(&pidRuntime.setpointDerivativePt1[axis], pt1FilterGain(filterCutoff, pidRuntime.dT));
                    break;
                case RC_SMOOTHING_DERIVATIVE_BIQUAD:
                    biquadFilterInitLPF(&pidRuntime.setpointDerivativeBiquad[axis], filterCutoff, targetPidLooptime);
                    break;
            }
        }
    }
}

void pidUpdateSetpointDerivativeLpf(uint16_t filterCutoff)
{
    if ((filterCutoff > 0) && (pidRuntime.rcSmoothingFilterType != RC_SMOOTHING_DERIVATIVE_OFF)) {
        for (int axis = FD_ROLL; axis <= FD_YAW; axis++) {
            switch (pidRuntime.rcSmoothingFilterType) {
                case RC_SMOOTHING_DERIVATIVE_PT1:
                    pt1FilterUpdateCutoff(&pidRuntime.setpointDerivativePt1[axis], pt1FilterGain(filterCutoff, pidRuntime.dT));
                    break;
                case RC_SMOOTHING_DERIVATIVE_BIQUAD:
                    biquadFilterUpdateLPF(&pidRuntime.setpointDerivativeBiquad[axis], filterCutoff, targetPidLooptime);
                    break;
            }
        }
    }
}
#endif // USE_RC_SMOOTHING_FILTER

void pidInitConfig(const pidProfile_t *pidProfile)
{
    if (pidProfile->feedForwardTransition == 0) {
        pidRuntime.feedForwardTransition = 0;
    } else {
        pidRuntime.feedForwardTransition = 100.0f / pidProfile->feedForwardTransition;
    }
    for (int axis = FD_ROLL; axis <= FD_YAW; axis++) {
        pidRuntime.pidCoefficient[axis].Kp = PTERM_SCALE * pidProfile->pid[axis].P;
        pidRuntime.pidCoefficient[axis].Ki = ITERM_SCALE * pidProfile->pid[axis].I;
        pidRuntime.pidCoefficient[axis].Kd = DTERM_SCALE * pidProfile->pid[axis].D;
        pidRuntime.pidCoefficient[axis].Kf = FEEDFORWARD_SCALE * (pidProfile->pid[axis].F / 100.0f);
        pidRuntime.dynThr[axis] = pidProfile->dynThr[axis];
        for (int pid = 0; pid <= 2; pid++) {
            pidRuntime.stickPositionTransition[pid][axis] = (pidProfile->stickTransition[pid][axis] / 100.0f) - 1.0f;
        }
    }
    pidRuntime.trueYawFF = YAW_TRUE_FF_SCALE * pidProfile->yaw_true_ff;

    pidRuntime.dtermMeasurementSlider = pidProfile->dtermMeasurementSlider / 100;
    pidRuntime.dtermMeasurementSliderInverse = 1 - (pidProfile->dtermMeasurementSlider / 100);

    pidRuntime.emuBoostPR = (pidProfile->emuBoostPR * pidProfile->emuBoostPR / 1000000) * 0.003;
    pidRuntime.emuBoostY = (pidProfile->emuBoostY * pidProfile->emuBoostY / 1000000) * 0.003;
    pidRuntime.emuBoostLimitPR = powf(pidProfile->emuBoostPR, 0.75f) * 1.4;
    pidRuntime.emuBoostLimitY = powf(pidProfile->emuBoostY, 0.75f) * 1.4;
    pidRuntime.dtermBoost = (pidProfile->dtermBoost * pidProfile->dtermBoost / 1000000) * 0.003;
    pidRuntime.dtermBoostLimit = powf(pidProfile->dtermBoost, 0.75f) * 1.4;
    pidRuntime.iDecay = pidProfile->i_decay;
    pidRuntime.iDecayCutoff = pidProfile->i_decay_cutoff;

    pidRuntime.P_angle_low = pidProfile->pid[PID_LEVEL_LOW].P * 0.1f;
    pidRuntime.D_angle_low = pidProfile->pid[PID_LEVEL_LOW].D * 0.00017f;
    pidRuntime.P_angle_high = pidProfile->pid[PID_LEVEL_HIGH].P * 0.1f;
    pidRuntime.D_angle_high = pidProfile->pid[PID_LEVEL_HIGH].D * 0.00017f;
    pidRuntime.F_angle = pidProfile->pid[PID_LEVEL_LOW].F * 0.00000125f;
    pidRuntime.horizonTransition = (float)pidProfile->horizonTransition;
    pidRuntime.horizonCutoffDegrees = pidProfile->racemode_tilt_effect;
    pidRuntime.horizonGain = pidProfile->horizonGain / 10.0f;
    pidRuntime.horizonTiltExpertMode = pidProfile->racemode_horizon;
    pidRuntime.horizonFactorRatio = (100 - pidProfile->racemode_tilt_effect) * 0.01f;
    pidRuntime.itermWindupPointInv = 0.0f;
    if (pidProfile->itermWindupPointPercent != 0) {
        const float itermWindupPoint = pidProfile->itermWindupPointPercent / 100.0f;
        pidRuntime.itermWindupPointInv = 1.0f / itermWindupPoint;
    }
    pidRuntime.itermAcceleratorGain = pidProfile->itermAcceleratorGain;
    pidRuntime.crashGyroThreshold = pidProfile->crash_gthreshold;
    pidRuntime.crashDtermThreshold = pidProfile->crash_dthreshold;
    pidRuntime.crashSetpointThreshold = pidProfile->crash_setpoint_threshold;
    pidRuntime.itermLimit = pidProfile->itermLimit;
#if defined(USE_THROTTLE_BOOST)
    throttleBoost = pidProfile->throttle_boost * 0.1f;
#endif
    pidRuntime.itermRotation = pidProfile->iterm_rotation;
    pidRuntime.antiGravityMode = pidProfile->antiGravityMode;

    // Calculate the anti-gravity value that will trigger the OSD display.
    // For classic AG it's either 1.0 for off and > 1.0 for on.
    // For the new AG it's a continuous floating value so we want to trigger the OSD
    // display when it exceeds 25% of its possible range. This gives a useful indication
    // of AG activity without excessive display.
    pidRuntime.antiGravityOsdCutoff = 0.0f;
    if (pidRuntime.antiGravityMode == ANTI_GRAVITY_SMOOTH) {
        pidRuntime.antiGravityOsdCutoff += (pidRuntime.itermAcceleratorGain / 1000.0f) * 0.25f;
    }
    pidRuntime.tpaBreakpoint = pidProfile->tpa_breakpoint;

#if defined(USE_ITERM_RELAX)
    pidRuntime.itermRelaxCutoff = pidProfile->iterm_relax_cutoff;
    pidRuntime.itermRelaxCutoffYaw = pidProfile->iterm_relax_cutoff_yaw;
    pidRuntime.itermRelaxThreshold = pidProfile->iterm_relax_threshold;
    pidRuntime.itermRelaxThresholdYaw = pidProfile->iterm_relax_threshold_yaw;
#endif

#ifdef USE_DYN_LPF
    if (pidProfile->dyn_lpf_dterm_min_hz > 0) {
        switch (pidProfile->dterm_filter_type) {
        case FILTER_PT1:
            pidRuntime.dynLpfFilter = DYN_LPF_PT1;
            break;
        case FILTER_BIQUAD:
            pidRuntime.dynLpfFilter = DYN_LPF_BIQUAD;
            break;
        default:
            pidRuntime.dynLpfFilter = DYN_LPF_NONE;
            break;
        }
    } else {
        pidRuntime.dynLpfFilter = DYN_LPF_NONE;
    }
    pidRuntime.dynLpfMin = pidProfile->dyn_lpf_dterm_min_hz;
    pidRuntime.dynLpfMax = pidProfile->dyn_lpf_dterm_max_hz;
    pidRuntime.dynLpfCurveExpo = pidProfile->dyn_lpf_curve_expo;
    pidRuntime.dynLpf2Gain = pidProfile->dterm_dynlpf2_gain;
    pidRuntime.dynLpf2Max = pidProfile->dterm_dynlpf2_fmax;
#endif

#ifdef USE_LAUNCH_CONTROL
    pidRuntime.launchControlMode = pidProfile->launchControlMode;
    if (sensors(SENSOR_ACC)) {
        pidRuntime.launchControlAngleLimit = pidProfile->launchControlAngleLimit;
    } else {
        pidRuntime.launchControlAngleLimit = 0;
    }
    pidRuntime.launchControlKi = ITERM_SCALE * pidProfile->launchControlGain;
#endif

#ifdef USE_THRUST_LINEARIZATION
    pidRuntime.thrustLinearization = pidProfile->thrustLinearization / 100.0f;
    pidRuntime.throttleCompensateAmount = pidRuntime.thrustLinearization - 0.5f * powerf(pidRuntime.thrustLinearization, 2);
#endif
#if defined(USE_D_MIN)
    for (int axis = FD_ROLL; axis <= FD_YAW; ++axis) {
        const uint8_t dMin = pidProfile->d_min[axis];
        if ((dMin > 0) && (dMin < pidProfile->pid[axis].D)) {
            pidRuntime.dMinPercent[axis] = dMin / (float)(pidProfile->pid[axis].D);
        } else {
            pidRuntime.dMinPercent[axis] = 0;
        }
    }
    pidRuntime.dMinGyroGain = pidProfile->d_min_gain * D_MIN_GAIN_FACTOR / D_MIN_LOWPASS_HZ;
    pidRuntime.dMinSetpointGain = pidProfile->d_min_gain * D_MIN_SETPOINT_GAIN_FACTOR * pidProfile->d_min_advance * pidRuntime.pidFrequency / (100 * D_MIN_LOWPASS_HZ);
    // lowpass included inversely in gain since stronger lowpass decreases peak effect
#endif
#ifdef USE_INTERPOLATED_SP
    pidRuntime.ffFromInterpolatedSetpoint = pidProfile->ff_interpolate_sp;
    pidRuntime.ffSmoothFactor = 1.0f - ((float)pidProfile->ff_smooth_factor) / 100.0f;
    interpolatedSpInit(pidProfile);
#endif

    pidRuntime.nfeRaceMode = pidProfile->nfe_racemode;
}

void pidCopyProfile(uint8_t dstPidProfileIndex, uint8_t srcPidProfileIndex)
{
    if (dstPidProfileIndex < PID_PROFILE_COUNT && srcPidProfileIndex < PID_PROFILE_COUNT
        && dstPidProfileIndex != srcPidProfileIndex) {
        memcpy(pidProfilesMutable(dstPidProfileIndex), pidProfilesMutable(srcPidProfileIndex), sizeof(pidProfile_t));
    }
}
