import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
import json

with open('config.json') as config_file:
    data = json.load(config_file)


#### Please go to config.json to change the configurations ####
alg=data['alg']
res=int(data['res'])
wr=data['wr']
Max=int(data['Max'])
plot_solution_var=data['plot_solution']
plot_number_of_transfos_var=data['plot_number_of_transformations']
run_unit_tests=data['run_unit_tests']

def plot_number_of_transfos(res):

  count = np.loadtxt("../cpp/output/number_of_transfos.txt",skiprows=1, unpack=True)

  plt.figure()
  plt.xscale('log')
  plt.yscale('log')
  x=np.arange(3,len(count)+3)
  plt.plot(x,count,'b.', label='Number of transformations')
  
  #Finds the second order polynomial which fits the function
  p=np.polyfit(x,count,2)
  plt.plot(x,p[2]+p[1]*x+p[0]*x*x, label='Polynomial fit by $%sx^2+%sx%s$'%(str(round(p[0],1)),str(round(p[1],1)),str(round(p[2],1))))
  plt.grid(1)
  plt.title("Number of transformations")
  plt.ylabel('Number of transformations')
  plt.xlabel('Size of matrix')
  plt.legend()
  plt.savefig("img/number_of_transfos_res_%s.png"%(res))


def plot_solution(alg , res ,wr ,Max):
  #Must take into account wr if we have interactions
  if(alg.lower()=='qho_int'):
    #Gets the eigenvals and eigenvecs
    j_vals, a_vals, v_0,v_1,v_2=np.loadtxt("output/%s_res_%d_wr=%s.txt"%(alg,res,str(wr)), skiprows=2, unpack=True)
    x=np.arange(0,Max,Max/len(v_0))
    
    plt.figure(figsize=(8, 6))
    plt.title("Interacting quantum harmonic oscillator ground state. $\omega_r=%s$"%(str(wr)))
    plt.plot(x,v_0, label = 'Ground state. $\lambda = %s$'%(str(round(j_vals[0],1))))
    plt.legend()
    plt.savefig("img/%s_res_%d_wr=%s.png"%(alg,res,str(wr)))

  elif(alg.lower()=='qho_no_int'):
    #Gets the eigenvals and eigenvecs
    j_vals, a_vals, v_0,v_1,v_2=np.loadtxt("output/%s_res_%d.txt"%(alg,res), skiprows=2, unpack=True)
    x=np.arange(0,Max,Max/len(v_0))
    
    plt.figure(figsize=(8, 6))
    plt.title("Noninteracting quantum harmonic oscillator ground state.")
    plt.plot(x,v_0, label = 'Ground state. $\lambda = %s$'%(str(round(j_vals[0],1))))
    plt.legend()
    plt.savefig("img/%s_res_%d.png"%(alg,res))
  else:
    j_vals, a_vals, v_0,v_1,v_2=np.loadtxt("output/%s_res_%d.txt"%(alg,res), skiprows=2, unpack=True)
    x=np.arange(0,Max,Max/len(v_0))

    plt.figure(figsize=(8, 6))
    plt.title("First three solutions to the buckeling beam problem.")
    plt.plot(x,v_0, label = 'First eigenvalue. $\lambda = %s$'%(str(round(j_vals[0],1))))
    plt.plot(x,v_1, label = 'Second eigenvalue. $\lambda = %s$'%(str(round(j_vals[1],1))))
    plt.plot(x,v_2, label = 'Third eigenvalue. $\lambda = %s$'%(str(round(j_vals[2],1))))
    plt.legend(prop={'size': 8})
    plt.savefig("img/%s_res_%d.png"%(alg,res))

  
if(plot_solution_var.lower()=='yes'):
  os.system("g++ -o %s %s.cpp jacobi.cpp -larmadillo && ./%s %s %s %s"%(alg,alg,alg,str(res),str(Max),str(wr)))
  plot_solution(alg,res,wr,Max)

if(plot_number_of_transfos_var.lower()=='yes'):
  os.system("g++ -o count count.cpp jacobi.cpp -larmadillo && ./count %s"%(str(res)))
  plot_number_of_transfos(str(res))

if(run_unit_tests.lower()=='yes'):
  os.system("g++ -o test-functions test-functions.cpp tests-main.cpp catch.hpp jacobi.cpp -larmadillo && ./test-functions")