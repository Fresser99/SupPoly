
# 导入需要的库
import os
import time
import numpy as np
# 创建AP的本地服务器
import win32com.client as win32
# AP8.8版替换下面的36.0为34.0; 9.0替换为35.0; 10.0替换为36.0； 11.0替换为37.0
Application = win32.Dispatch('Apwn.Document.40.0')

# 获取当前文件夹的地址，三效蒸发器文件和本程序文件需放置在同一个文件夹
address = os.getcwd()

# AP的bkp文件的文件名，即三效蒸发器bkp文件的文件名
SimulationName = 'ddd'

# 打开三效蒸发器文件
Application.InitFromArchive2('D:\\poly\\SupPoly\\sd\\ddd$backup2222.bkp')

# 设置AP用户界面的可见性，1为可见，0为不可见
Application.Visible = 1

# 压制对话框的弹出，1为压制；0为不压制
Application.SuppressDialogs = 1
time.sleep(10)
# 试运行三效蒸发器模拟
Application.Engine.Run2(1)

# 因为AP运行比较慢，所以用一个循环语句，每两秒钟检查一次是否运行完毕
while Application.Engine.IsRunning == 1:
    time.sleep(2)

print('run complete')
# 定义输入、输出变量、收敛标签
# FeedFlow = [1200, 1400]
# HeatDuty = []
# SimulationConvergency = []
#
# # 程序核心部分
#
# # 用for循环依次改变流量为1200、1400，并运行模拟，然后存储数据，检查历史文件
# for ii in range(len(FeedFlow)):
#     # 设置输入流量
#     Application.Tree.FindNode("\Data\Streams\BRINE\Input\TOTFLOW\MIXED").Value = FeedFlow[ii]
#
#     # 初始化，用于清除之前的数据
#     Application.Reinit
#
#     # 运行程序
#     Application.Engine.Run2(1)
#
#     # 定义本次的输出数据
#     heatduty = []
#
#     # 每两秒钟检查一次是否运行完毕
#     while Application.Engine.IsRunning == 1:
#         time.sleep(2)
#
#     # 获取三个蒸发器的热负荷
#     heatduty.append(Application.Tree.FindNode("\Data\Blocks\STAGE-1\Output\QCALC").Value)
#     heatduty.append(Application.Tree.FindNode("\Data\Blocks\STAGE-2\Output\QCALC").Value)
#     heatduty.append(Application.Tree.FindNode("\Data\Blocks\STAGE-3\Output\QCALC").Value)
#
#     # 存储到HeatDuty中
#     HeatDuty.append(heatduty)
#
#     # 获取运行ID，并以文本格式打开相应的历史文件，搜索文件中是否有关键词
#     name = Application.Tree.FindNode("\Data\Results Summary\Run-Status\Output\RUNID").Value
#     Filename = address + '\\' + name + '.his'
#
#     # 如果关键词出现，则标记这一次结果为没有收敛
#     with open(Filename, 'r') as f:
#         isError = np.any(np.array([line.find('SEVERE ERROR') for line in f.readlines()]) >= 0)
#         SimulationConvergency.append(not isError)

# 关闭AP
Application.Close
Application.Quit