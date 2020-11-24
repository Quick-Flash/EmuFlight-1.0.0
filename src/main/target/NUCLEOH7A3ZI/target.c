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

#include <stdint.h>

#include "platform.h"
#include "drivers/io.h"

#include "drivers/dma.h"
#include "drivers/timer.h"
#include "drivers/timer_def.h"

const timerHardware_t timerHardware[USABLE_TIMER_CHANNEL_COUNT] = {
    // DEF_TIM(TIM2,  CH2, PB3,  TIM_USE_LED,                 0,  7,  0 ), // SPI1_SCK
    DEF_TIM(TIM12, CH2, PB15, TIM_USE_LED,                 0,  0,  0 ),

    DEF_TIM(TIM1,  CH2, PE11, TIM_USE_PWM|TIM_USE_PPM,     0,  0,  0 ),
    DEF_TIM(TIM1,  CH1, PE9,  TIM_USE_PWM,                 0,  0,  0 ),
    DEF_TIM(TIM2,  CH1, PA15, TIM_USE_PWM,                 0,  0,  0 ),

    DEF_TIM(TIM8,  CH1, PC6,  TIM_USE_MOTOR,                0, 0,  4 ),
    DEF_TIM(TIM8,  CH2, PC7,  TIM_USE_MOTOR,                0, 1,  4),
    DEF_TIM(TIM8,  CH3, PC8,  TIM_USE_MOTOR,                0, 2,  4),
    DEF_TIM(TIM8,  CH4, PC9,  TIM_USE_MOTOR,                0, 3,  4),

    DEF_TIM(TIM3,  CH1, PB4,  TIM_USE_MOTOR,               0,  4,  0 ),
    DEF_TIM(TIM3,  CH2, PB5,  TIM_USE_MOTOR,               0,  5,  0 ),

    DEF_TIM(TIM5,  CH3, PA2,  TIM_USE_MOTOR,               0,  0,  0 ),
    DEF_TIM(TIM5,  CH2, PA1,  TIM_USE_MOTOR,               0,  4,  0 ),
    DEF_TIM(TIM5,  CH4, PA3,  TIM_USE_MOTOR,               0,  1,  0 ),
    DEF_TIM(TIM5,  CH1, PA0,  TIM_USE_MOTOR,               0,  2,  0 ),

    // DEF_TIM(TIM4,  CH3, PB8,  TIM_USE_MOTOR,               0,  8,  0 ), // I2C1
    // DEF_TIM(TIM4,  CH4, PB9,  TIM_USE_MOTOR,               0,  0,  0 ), // I2C1
    // DEF_TIM(TIM9,  CH2, PE6,  TIM_USE_MOTOR,               0,  0,  0 ),
};
