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

#include "rx/rx_spi.h"

#define INTERVAL_RX_LOSS_MS 1000
#define INTERVAL_RX_BIND_MS 250

void rxSpiCommonIOInit(const rxSpiConfig_t *rxSpiConfig);

void rxSpiLedOn(void);
void rxSpiLedOff(void);
void rxSpiLedToggle(void);
void rxSpiLedBlink(timeMs_t blinkMs);
void rxSpiLedBlinkRxLoss(rx_spi_received_e result);
void rxSpiLedBlinkBind(void);

void rxSpiBind(void);
bool rxSpiCheckBindRequested(bool reset);
