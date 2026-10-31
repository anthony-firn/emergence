// Add to AM32's Inc/targets.h, then: make ARM_SDK_PREFIX=arm-none-eabi- EMERGENCE12S_F051
// Written against am32-firmware/AM32 at 86cf1e5. Not yet compiled.
#ifdef EMERGENCE12S_F051
// The ESC board of the Emergence inspection quadcopter: group B's pins, an 18:1 bus divider
// for 12S, a 1 mR shunt into an INA240A2 (50 mV/A), and 500 ns of dead time for a UCC27211.
#define FILE_NAME "EMERGENCE12S_F051"
#define FIRMWARE_NAME "Emergence12S"
#define DEAD_TIME 24
#define HARDWARE_GROUP_F0_B
#define MILLIVOLT_PER_AMP 50
#define CURRENT_OFFSET 0
#define TARGET_VOLTAGE_DIVIDER 180
#define VOLTAGE_ADC_CHANNEL LL_ADC_CHANNEL_6
#define VOLTAGE_ADC_PIN LL_GPIO_PIN_6
#define CURRENT_ADC_CHANNEL LL_ADC_CHANNEL_3
#define CURRENT_ADC_PIN LL_GPIO_PIN_3
#define USE_SERIAL_TELEMETRY
#define TARGET_STALL_PROTECTION_INTERVAL 9000
#endif
