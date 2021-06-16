# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 08:05:52 2020

@author: BK
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

def par(alpha,bata):
    """
    並聯運算
    """
    par_result = (alpha*bata)/(alpha+bata)       
    return par_result

def circuit_list2mat_diagonal(circuit_list_main,anti_diagonal=0):
    """
    建立矩陣，輸入一list例如[A,B,C,D]
    list內有4個element，為了建立符合這個list的矩陣
    因此建立一常數J獲得list內的element的數量，python的計數從0開始，因此數量要減1
    才能建立符合list的矩陣。
    anti_digonal:如矩陣建立方向由z(4,4)建立到z(1,1),正常不會用到，我寫爽的
    """
    j=np.int(len(circuit_list_main))-1
    circuit_mat=np.zeros([len(circuit_list_main),len(circuit_list_main)],dtype=complex)
    if anti_diagonal == 0:
        for i in range(len(circuit_list_main)):
            circuit_mat[i,i]=circuit_list_main[i]
    elif anti_diagonal == 1:
        for i in range(len(circuit_list_main)):
            circuit_mat[j,i]=circuit_list_main[i]
            j=j-1
    else:
        print("input index/diagonal error,please check!")
    return circuit_mat

def circuit_mat_subelement(circuit_mat,circuit_list_sub):
    """
    將circuit_list2mat_diagonal建立起的對角線矩陣帶入，並加入非對角線的矩陣element
    同時建立矩陣對稱性關聯z(2,1)=z(1,2)
    """
    circuit_mat_whole=np.zeros([circuit_mat.shape[0],circuit_mat.shape[0]],dtype=complex)
    for i in range(len(circuit_mat)-1):
        for j in range(len(circuit_mat)-1):
            circuit_mat_whole[i,j+1]=circuit_list_sub[i,j]
    for i in range(len(circuit_mat)):
        for j in range(len(circuit_mat)):
            circuit_mat_whole[j,i]=circuit_mat_whole[i,j]
    circuit_mat_whole=circuit_mat_whole+circuit_mat
    return circuit_mat_whole

def circuit_impedence_cal (imp_matrix,vol_mat_number,force_mat_number,curve_type):
    """
    計算速度/體積速度:
    ex:[1 0 0 0]*inv(Z)*[1;0;0;0]
    判讀等效電路的長寬後，建立同等元素數的列、行矩陣
    依據目標的迴路，給與行、列矩陣需要的目標元素位置(vol_mat_number,force_mat_number)
    再以上方舉例的運算式算出目標的速度/體積速度(依使用者轉換得domain決定)
    
    下方暫定註解(有誤)
    impedance_type:
    0 for frequencr response
    1 for impedance curve    
    """
    impedance=np.zeros(1,dtype=complex)
    vol_mat=np.zeros([1,len(imp_matrix)])
    vol_mat[0,vol_mat_number-1]=1
    force_mat=np.zeros([len(imp_matrix),1])
    force_mat[force_mat_number-1,0]=1
    if curve_type == 0:
        imp_matrix_inv=np.linalg.inv(imp_matrix)
        impedance=(np.dot(np.dot(vol_mat, imp_matrix_inv),force_mat))
    elif curve_type == 1:
        imp_matrix_inv=np.linalg.inv(imp_matrix)
        impedance=(np.dot(np.dot(vol_mat, imp_matrix_inv),force_mat))
    return impedance, imp_matrix_inv

def complex_value_calculate(complex_number):
    return cmath.sqrt(np.real(complex_number)**2+np.imag(complex_number)**2)
    

