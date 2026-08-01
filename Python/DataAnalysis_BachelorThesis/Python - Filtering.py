import pandas as pd
import numpy as np
import os
from scipy import signal
import matplotlib.pyplot as plt
import math

def Import_Data(Test_Number,Excel_name,drive):
    Dir_name = os.path.dirname(__file__)
    path = "{}:\\DATA\\TEXT{}.txt".format(drive, Test_Number)
    Read_Data= pd.read_csv(path,sep=',')
    N_sensor = int((len(Read_Data.columns)-1)/6)
    Header = np.empty(1+N_sensor*6,dtype=object)
    Header[0] = "Time"
    Dir = np.array(['accX','accY','accZ','gyrX','gyrY','gyrZ'])
    for y in range(6):
        for x in range(N_sensor):
            Header[1+y*N_sensor+x] = Dir[y]+str(x+1)
    Excel_path = 'Data/{}.xlsx'.format(Excel_name)
    Excel_path = os.path.join(Dir_name,Excel_path)
    Read_Data.to_excel(Excel_path)
    return Excel_path


class Complementary_filter():
    def __init__(self, path,rangedata=0):
        self.path = path
        self.Data = pd.read_excel(self.path)
        if rangedata == 0: 
            self.Data = self.Convert_data(self.Data)
        else:
            self.Data = self.Data[rangedata[0]:rangedata[1]]
            self.Data = self.Convert_data(self.Data)
            
    def Convert_data(self,Data):
        Data['Time']  = Data['Time'].div(1000)
        Data['gyrX1'] = Data['gyrX1'].div(100)
        Data['gyrX2'] = Data['gyrX2'].div(100)
        Data['gyrX3'] = Data['gyrX3'].div(100)
        Data['gyrY1'] = Data['gyrY1'].div(100)
        Data['gyrY2'] = Data['gyrY2'].div(100)
        Data['gyrY3'] = Data['gyrY3'].div(100)
        Data['gyrZ1'] = Data['gyrZ1'].div(100)
        Data['gyrZ2'] = Data['gyrZ2'].div(100)
        Data['gyrZ3'] = Data['gyrZ3'].div(100)
        # Indexing each row (data set)
        Data.rename( columns={'Unnamed: 0':'N'}, inplace=True )
        return Data
    
    def get_Fundamental_Data(self):
        Time = self.Data['Time'][self.Data.index[-1]]
        Samples = self.Data.index[-1]
        SamplesPerSecond = Samples/Time
        Frequency = 1/SamplesPerSecond
        print("Time: %s, Samples: %s, SamplesPerSecond: %s, Frequency: %s" % (Time, Samples, SamplesPerSecond, Frequency))
    
    def lfilter(self,n,data):
        b = [1.0 / n] * n
        a = 1
        data_filtered = signal.lfilter(b, a, data)
        return data_filtered
    
    def GetAcc_angle(self, Trans_1,Trans_2,n_filter):
        value1 = "acc{}".format(Trans_1)
        value2 = "acc{}".format(Trans_2)
        acc_angle_Raw = np.arctan2(self.Data[value1],self.Data[value2])*(180/math.pi)
        self.acc_angle_Raw = np.absolute(acc_angle_Raw) # Remove inverted sign (Sudden drops)
        self.acc_angle_Fil = self.lfilter(n_filter, self.acc_angle_Raw)
        
        
    def Gyr_filter(self,Trans_1,n_filter, Inverse = 0):
        value1 = "gyr{}".format(Trans_1)
        self.gyr_Fil = self.lfilter(n_filter,self.Data[value1])
        if Inverse == 0:    
            self.gyr_Fil *= -1
        
    def get_angle(self,state,avg_n,rangedata = 0):
        avg = 0
        avg = np.mean(self.gyr_Fil[:avg_n])*-1
        self.gyr_Fil += avg
        # If no range
        if rangedata == 0:
            rangedata = [0,len(self.Data)]
        #Deciding the start angle with accelerometer
        angle = np.zeros(len(self.Data))
        angle[0] = self.acc_angle_Fil[0]

        #Finding the angle for each iteration (data set)
        for i in range(rangedata[1]-rangedata[0]-1):
            dt = self.Data.iloc[i+1]['Time']-self.Data.iloc[i]['Time']  
            if state == 0:
                angle[i+1] = 0.98*(angle[i]+self.gyr_Fil[i+1]*dt)+0.02*self.acc_angle_Fil[i+1]
            if state == 1:
                angle[i+1] = 1*(angle[i]+self.gyr_Fil[i+1]*dt)
        return angle
    def makeithappen(self,state,Trans_1,Trans_2,Trans_1_gyr, avg_n = 100,n_filter = 1,Inverse = 0):
        #self.get_Fundamental_Data()
        self.Gyr_filter(Trans_1_gyr, n_filter, Inverse)
        self.GetAcc_angle(Trans_1,Trans_2, n_filter)
        angle = self.get_angle(state, avg_n, rangedata)
        return angle
        

    #Importing data
drive = "E"
Excel_name = "Meh"
Test_Number = "6"
file_n = "ExcelDynamic"
Excel_Path = 'Data/{}.xlsx'.format(file_n)
avg_n = 200
#rangedata = [5000,5600]  
rangedata = 0  
n_filter = 1 # the larger n is, the smoothlener curve will be (lfilter
CF = False # 0 = Complementary filter or Gyroscope?
Inverse = 1 # Inverse the result of gyr
File_Name_Export = "ValidationData_Dynamic"

if __name__ == "__main__":
    #Importing data from Txt file
    Excel_Path = Import_Data(Test_Number, Excel_name, drive)
    # Using complementary filter
    uf = Complementary_filter(Excel_Path, rangedata)
    angle1 = uf.makeithappen(CF,"X1","Z1", "Y1",avg_n,n_filter)
    angle2 = uf.makeithappen(CF, "X2", "Z2", "Y2",avg_n, n_filter)
    angle3 = uf.makeithappen(CF, "X3", "Z3", "Y3",avg_n, n_filter,Inverse)
    
    # Deleting useless data from start-index
    start = 0
    end = 0
    range1 = np.arange(start,end)
    Time = uf.Data.drop(uf.Data.index[[range1]])
    angle1 = np.delete(angle1,np.arange(start,end))
    angle2 = np.delete(angle2,np.arange(start,end))
    angle3 = np.delete(angle3,np.arange(start,end))
    
    # Plotting
    plt.plot(Time['Time'],angle1)
    plt.plot(Time['Time'],angle2)
    plt.plot(Time['Time'],angle3)
    
    # # Exporting to Excel
    x = np.array([Time['Time'],angle1, angle2, angle3])
    m = np.asmatrix(x)
    m = m.T
    df = pd.DataFrame(m, columns=['Time','Angle1','Angle2','Angle3'])
    Export = 'Data/{}.xlsx'.format(File_Name_Export)
    df.to_excel(Export)