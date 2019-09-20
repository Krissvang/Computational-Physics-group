import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

alg =2
res=100
wr="0.5"
Max=10

def plot_transfos():

  count = np.loadtxt("../cpp/output/number_of_transfos.txt",skiprows=1, unpack=True)

  plt.figure()
  plt.xscale('log')
  plt.yscale('log')
  x=np.arange(3,len(count)+3)
  plt.plot(x,count,'b.')
  
  p=np.polyfit(x,count,2)
  
  plt.plot(x,p[2]+p[1]*x+p[0]*x*x)

  plt.grid(1)
  plt.title("Number of transformations")
  plt.ylabel('Number of transformations')
  plt.xlabel('Size of matrix')
  plt.savefig("number_of_transfos.png")


def plot_eigenpairs(alg , res ,wr ,Max):
  if(alg==0):
    j_vals, a_vals, v_0,v_1,v_2=np.loadtxt("output/beam_res_%d.txt"%(res), skiprows=2, unpack=True)
    x=np.arange(0,1,1/len(v_0))
    
    plt.figure()
    plt.plot(x,v_0)
    plt.plot(x,v_1)
    plt.plot(x,v_2)
    plt.savefig("img/Beam_res_%d.png"%(res))

  elif(alg==1):
    j_vals, a_vals, v_0,v_1,v_2=np.loadtxt("output/qho_no_int_res_%d.txt"%(res), skiprows=2, unpack=True)
    x=np.arange(0,Max,Max/len(v_0))
    
    plt.figure()
    plt.plot(x,v_0)
    plt.plot(x,v_1)
    plt.plot(x,v_2)
    plt.savefig("img/QHO_no_int_res_%d.png"%(res))

  elif(alg==2):
    j_vals, a_vals, v_0,v_1,v_2=np.loadtxt("output/qho_int_res_%d_wr=%s.txt"%(res,wr), skiprows=2, unpack=True)
    x=np.arange(0,Max,Max/len(v_0))
    
    plt.figure()
    plt.plot(x,v_0)
    plt.plot(x,v_1)
    plt.plot(x,v_2)
    plt.savefig("img/QHO_int_res_%d_wr=%s.png"%(res,wr))



print("Algorithm ", alg)
plot_eigenpairs(alg,res,wr,Max)
