#If using python2 remove comment from ine under
#from __future__ import unicode_literals

import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import json

#Opens the comfigurantion file
with open('config.json') as config_file:
    data = json.load(config_file)


#### Please go to config.json to change the configurations ####
#Configurations that determines the plots and algorithms
alg=data['alg']
res=int(data['res'])
wr=data['wr']
Max=int(data['Max'])
plot_solution_var=data['plot_solution']
plot_number_of_transfos_var=data['plot_number_of_transformations']
run_unit_tests=data['run_unit_tests']

#Function plotting the number of transformations vs the size of the matrix.
def plot_number_of_transfos(res):
  #Data to compare
  count = np.loadtxt("../cpp/output/number_of_transfos.txt",skiprows=1, unpack=True)
  
  #Plotting configurations
  plt.figure()
  plt.xscale('log')
  plt.yscale('log')
  plt.grid(1)
  plt.title("Number of transformations")
  plt.ylabel('Number of transformations')
  plt.xlabel('Size of matrix (n)')

  #Plots
  x=np.arange(3,len(count)+3)
  plt.plot(x,count,'b.', label='Number of transformations')
  
  #Finds the second order polynomial which fits the function and plots it
  p=np.polyfit(x,count,2)
  plt.plot(x,p[2]+p[1]*x+p[0]*x*x, label='Polynomial fit by $%sx^2+%sx%s$'%(str(round(p[0],1)),str(round(p[1],1)),str(round(p[2],1))))
  
  plt.legend()
  plt.savefig("img/number_of_transfos_res_%s.pdf"%(res))


def plot_solution(alg , res ,wr ,Max):
  #Plotting the two electron problem
  if(alg.lower()=='qho_int'): 
    
    #Gets the eigenvals and eigenvecs
    v_0=np.zeros(res+2)
    j_vals, a_vals, v_0[1:-1],v_1,v_2=np.loadtxt("output/%s_res_%d_wr=%s.txt"%(alg,res,str(wr)), skiprows=2, unpack=True)
    x=np.linspace(0,Max,res+2)

    #Plots and plotting configs
    plt.figure(figsize=(8, 6))
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"$\varphi(\rho)$")
    plt.title("Interacting quantum harmonic oscillator ground state. $\omega_r=%s$"%(str(wr)))
    plt.plot(x,v_0, label = 'Ground state. $\lambda = %s$'%(str(round(j_vals[0],1))))
    plt.legend()
    plt.savefig("img/%s_res_%d_wr=%s.pdf"%(alg,res,str(wr)))

  #One electron problem
  elif(alg.lower()=='qho_no_int'):
    #Gets the eigenvals and eigenvecs
    j_vals, a_vals, v_0,v_1,v_2=np.loadtxt("output/%s_res_%d.txt"%(alg,res), skiprows=2, unpack=True)
    x=np.arange(0,Max,Max/len(v_0))
    
    #Plots and plotting configs
    plt.figure(figsize=(8, 6))
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"$\varphi(\rho)$")
    plt.title("Noninteracting quantum harmonic oscillator ground state.")
    plt.plot(x,v_0, label = 'Ground state. $\lambda = %s$'%(str(round(j_vals[0],1))))
    plt.legend()
    plt.savefig("img/%s_res_%d.pdf"%(alg,res))

    #Beam problem
  else:
    #Gets eigenvecs and eigenvals
    v_0 =np.zeros(res+2)
    v_1 =np.zeros(res+2)
    v_2 =np.zeros(res+2)
    j_vals, a_vals, v_0[1:-1],v_1[1:-1],v_2[1:-1]=np.loadtxt("output/%s_res_%d.txt"%(alg,res), skiprows=2, unpack=True)
    x=np.arange(0,1,1/len(v_0))

    #Plots and plotting configs
    plt.figure(figsize=(8, 6))
    plt.xlabel(r"$\rho = x/L$")
    plt.ylabel(r"u($\rho$)") 
    plt.title("First three solutions to the buckeling beam problem.")
    plt.plot(x,v_0, label = 'First eigenvalue. $\lambda = %s$'%(str(round(j_vals[0],1))))
    plt.plot(x,v_1, label = 'Second eigenvalue. $\lambda = %s$'%(str(round(j_vals[1],1))))
    plt.plot(x,v_2, label = 'Third eigenvalue. $\lambda = %s$'%(str(round(j_vals[2],1))))
    plt.legend(prop={'size': 8})
    plt.savefig("img/%s_res_%d.pdf"%(alg,res))


#Determines what runs decided by the configurations
if(plot_solution_var.lower()=='yes'):
  os.system("g++ -O3 -o %s %s.cpp jacobi.cpp -larmadillo && ./%s %s %s %s"%(alg,alg,alg,str(res),str(Max),str(wr)))
  plot_solution(alg,res,wr,Max)

if(plot_number_of_transfos_var.lower()=='yes'):
  os.system("g++ -O3 -o count count.cpp jacobi.cpp -larmadillo && ./count %s"%(str(res)))
  plot_number_of_transfos(str(res))

if(run_unit_tests.lower()=='yes'):
  os.system("g++ -O3 -o test-functions test-functions.cpp tests-main.cpp catch.hpp jacobi.cpp -larmadillo && ./test-functions")