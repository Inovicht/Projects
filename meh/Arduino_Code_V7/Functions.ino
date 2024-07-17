// Several functions

void Interrupt()                           // Interrupt - Starting test with Switch
{
  static unsigned long last_interrupt_time = 0;
  unsigned long interrupt_time = millis();
  if (interrupt_time - last_interrupt_time > 300) 
  {
    start = !start;
  }
  last_interrupt_time = interrupt_time;
}

void Data_Collector(ICM_20948_I2C *sensor, byte x)
{
  Values[x]                  = (sensor->accX()); // 0-2
  Values[NUMBER_SENSOR+x]    = (sensor->accY()); // 3-5
  Values[2*NUMBER_SENSOR+x]  = (sensor->accZ()); // 6-8
  Values[3*NUMBER_SENSOR+x]  = (sensor->gyrX()*100); // 9-11
  Values[4*NUMBER_SENSOR+x]  = (sensor->gyrY())*100; // 12-14
  Values[5*NUMBER_SENSOR+x]  = (sensor->gyrZ())*100; // 15-17
}

void Data_Store()
{
  file.print(millis()-Time_reset);
  for(byte i=0;i<sizeof(Values)/sizeof(int);i++)
  {
    file.write(",");
    file.print(Values[i]);
  }
  file.println();
}
void Start_Sensor()
{
   // Starting the sensors
  for (byte x=0; x<NUMBER_SENSOR; x++)
  {
    initialized = 0;
    MUX.setPort(x);
    while(!initialized)
    {
      if(IMU.begin() != 0)
      {
        Ser.print("Sensor: ");
        Ser.print(x);
        Ser.println(" Not working");
        digitalWrite(PIN,HIGH);delay(1000);}
      else {initialized = true; 
        Ser.print("Sensor: ");
        Ser.print(x);
        Ser.println(" working");
      }
    }
  }
}
void Start_SDCard()
{
  if(!SDCard.begin(ChipSelect,SPI_SPEED))
  {Ser.println("SDCard Init Failed! - While(1)");while(1);}
  else{Ser.println("SDCard Init Worked");}
}
void File_Maker()
{
  String Dir = "DATA";
  String File_Name = "Test";
  if(!SDCard.mkdir(Dir))
    {SDCard.mkdir(Dir);Ser.println("Directory created");}
    byte m = 0;
    while(1)
    {
      Name = "";
      Name += Dir;
      Name += "/";
      Name += File_Name;
      Name += String(m);
      Name += ".txt";
      if(!SDCard.exists(Name))
      {
        char navn[Name.length()+1];
        Name.toCharArray(navn,sizeof(navn));
        Ser.println("Creating File");
        if(!file.open(Name.c_str(), O_WRONLY | O_CREAT | O_EXCL)) {}
        file.close();
        break;
      }
      else{m++; Ser.println("File_name exists");}
    }
}
