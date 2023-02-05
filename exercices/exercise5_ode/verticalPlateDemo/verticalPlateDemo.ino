// A simple demo case on using Arduino and infrared temperature sensor to live plot data
// Modified version of this tutorial
//  http://wiki.seeedstudio.com/Grove-Digital_Infrared_Temperature_Sensor/
// code based on this
//  https://github.com/Seeed-Studio/Digital_Infrared_Temperature_Sensor_MLX90615/blob/master/examples/singleDevice/singleDevice.ino
// MIT license, see original for details.
//
// Antti Mikkonen, a.mikkonen@iki.fi, 2018
// Tampere University, 

// Include libraries
#include "MLX90615.h"
#include <I2cMaster.h>

// Choose pins for data and clock signals
#define SDA_PIN 4    //define the SDA pin, digital
#define SCL_PIN 5    //define the SCL pin, digital

// Initialize infrared 
SoftI2cMaster i2c(SDA_PIN, SCL_PIN);
MLX90615 mlx90615(DEVICE_ADDR, &i2c);

// This runs once at start up
void setup()
{
  // Initializes serial communication, data is red on Python side  
  Serial.begin(9600);

  // Emissivity could be defined here 
  //mlx90615.writeEEPROM(Default_Emissivity); //write data into EEPROM when you need to adjust emissivity.
  //mlx90615.readEEPROM(); //read EEPROM data to check whether it's a default one.
}

// This keeps looping until shutdown
void loop()
{
  // Read infrared temperature and write to serial  
  Serial.print(mlx90615.getTemperature(MLX90615_OBJECT_TEMPERATURE));
  // Write comma to serial
  Serial.print(", ");
  // Read sensor temperature and write to serial
  Serial.println(mlx90615.getTemperature(MLX90615_AMBIENT_TEMPERATURE));
  // Wait 1000ms (1s) 
  delay(1000);
}
