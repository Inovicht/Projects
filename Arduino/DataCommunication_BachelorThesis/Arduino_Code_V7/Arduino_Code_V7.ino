#include <SPI.h>
#include <SdFat.h>
#include <Wire.h>
#include "ICM_20948.h"
#include <SparkFun_I2C_Mux_Arduino_Library.h>
SdFat SDCard;  // Initializing the SD-Card object
SdFile file;
ICM_20948_I2C IMU; // Initializing the IMU object
QWIICMUX MUX;
#define AD0_VAL 1
#define Ser Serial
#define NUMBER_SENSOR 3 // Modulart - Angiv antal Sensorer
#define InterruptPin 2
#define RED_PIN 3
#define BLUE_PIN 7
#define PIN 13
#define SPI_SPEED SD_SCK_MHZ(4)
#define ChipSelect 10
float Period = 0; // 60 Hz?
int Values[NUMBER_SENSOR*6];
bool initialized = false;
String Name;
unsigned long Time_Now = 0;
unsigned long Time_reset = 0;
bool start = 0;

void setup() 
{
  Ser.begin(2000000);                       // Initializing Serial (printing)
  pinMode(RED_PIN,OUTPUT);                 // High when not running
  pinMode(BLUE_PIN,OUTPUT);                // High when running test
  pinMode(InterruptPin, INPUT_PULLUP);     // High when not presssed Pin 2 also suitable for interrupt
  attachInterrupt(digitalPinToInterrupt(InterruptPin),&Interrupt,FALLING);
  Wire.begin(); // Initializing I2C
  MUX.begin();                            // Initializing Multiplexer
  Start_Sensor();                         // Checking Sensor if ready
  Start_SDCard();                         // Checking SDCard if ready
  Wire.setClock(400000); // Setting clock
  Ser.println("Initialization is done");
  digitalWrite(RED_PIN,HIGH);
}
void loop() 
{
  if(start)                            // When interrupt, start test
  {
    File_Maker();                         // Making file on SDCard - Data folder
    Time_reset = millis();                // Taking time at test start - Usable for multiple test without compile
    file.open(Name.c_str(), O_WRONLY | O_APPEND);
    digitalWrite(BLUE_PIN,HIGH);          // Indicate start of test
    digitalWrite(RED_PIN,LOW);
    while(1)
    {
      if(millis() >= Time_Now+Period)     // having frequency for each test
      {
        Time_Now += Period;
        Run();
        if(!start || !file){file.close();digitalWrite(RED_PIN,HIGH);digitalWrite(BLUE_PIN,LOW);break;}
        Data_Store();
      }
    }
  }
  delay(200);
}

void Run()
{
  for(byte x=0; x<NUMBER_SENSOR;x++)
  {
    MUX.setPort(x);                       // Setting port on Multiplexer
    IMU.getAGMT();                        // Loading Data
    if(!IMU.dataReady()){start = !start; Ser.println("Error - Sensor");break;}
    else
    {
      if(file)
      {Data_Collector(&IMU, x);}         // Data is ready - Collecting in array
    }
  }
}
