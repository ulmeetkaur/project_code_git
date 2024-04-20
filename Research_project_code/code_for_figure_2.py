import numpy as np
import random
import matplotlib.pyplot as plt

start_time = 0.0                            # start time of the simulation
time = [0.0]                                # list to save time for each iteration


initial_conditions = np.array([0,                       # initial no. of molecules ptch_mRNA
                               3010,                    # initial no. of molecules ptch_cyt
                               301,                     # initial no. of molecules ptch_mem
                               6020,                    # initial no. of molecules Shh (20nm)
                               0,                       # initial no. of molecules Shh-ptch
                               12040,                   # initial no. of molecules cholesterol(inactive)
                               12040,                   # initial no. of molecules cholesterol(active)
                               0,                       # initial no. of molecules Smo
                               0,                       # initial no. of molecules Smo-Chol
                               0,                       # initial no. of molecules GliFL
                               0,                       # initial no. of molecules GliA
                               2500])                   # initial no. of molecules GliR
k = np.copy(initial_conditions)
gliA = [0]                                               # to save the result of gli_A after each increment

stoichiometry=np.array([[ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],                 # Reaction 1
                        [-1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        [-1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        [ 0, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        [ 0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        [ 0,  0, -1, -1,  1,  0,  0,  0,  0,  0,  0,  0],
                        [ 0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0],
                        [ 0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0],
                        [ 0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0],
                        [ 0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0],
                        [ 0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0],
                        [ 0,  0,  0,  0,  0,  0, -1, -1,  1,  0,  0,  0],
                        [ 0,  0,  0,  0,  0,  0,  1,  1, -1,  0,  0,  0],
                        [ 0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0],
                        [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0],
                        [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  1,  0],
                        [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0],
                        [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0],
                        [ 0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  1],
                        [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1],                      # Reaction 20
                        [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],])


individual_rates = np.array([40.6845,                              # r1- rate of reaction 1
                  0.00626 * k[0],                                  # r2- rate of reaction 2
                  0.0058 * k[0],                                   # r3- rate of reaction 3
                  0.0036 * k[1],                                   # r4- rate of reaction 4
                  0.0068 * k[2],                                   # r5- rate of reaction 5
                  0.0000143 * k[2] * k[3],                         # r6- rate of reaction 6
                  0.0068 * k[4],                                   # r7- rate of reaction 7
                  0.009 * k[5],                                    # r8- rate of reaction 8
                  0.25 * k[2] * k[6],                              # r9- rate of reaction 9
                  40.6485,                                         # r10- rate of reaction 10
                  0.0058 * k[7],                                   # r11- rate of reaction 11
                  0.0001072 * k[6] * k[7],                         # r12- rate of reaction 12
                  0.004 * k[8],                                    # r13- rate of reaction 13
                  0.0019 * k[8],                                   # r14- rate of reaction 14
                  45.165,                                          # r15- rate of reaction 15
                  0.00008701 * k[8] * k[9],                        # r16- rate of reaction 16
                  0.0029 * k[9],                                   # r17- rate of reaction 17
                  0.0173 * k[10],                                  # r18- rate of reaction 18
                  0.0015 * k[9],                                   # r19- rate of reaction 19
                  0.005 * k[11],                                   # r20- rate of reaction 20
                  0.05 *(k[10]/k[11])])                            # mrna upreg

end_time = 2880                                                    # max time


while time[-1] < end_time:
    Total_rate = np.sum(individual_rates)                          # Total reaction rate
    if Total_rate == 0:
        break
    tau = np.random.exponential(scale=1 / Total_rate)              # time interval generated randomly
    time.append(time[-1] + tau)

    rand = random.uniform(0, Total_rate)
    cumulative_sum = np.cumsum(individual_rates)
    reaction_index = np.argmax(cumulative_sum >= rand)              # determine which reaction will occur
    k += stoichiometry[reaction_index, : ]                          # update number of molecules after each iteration
    result_gliA= k[10]
    gliA.append(result_gliA)


    #to update reaction rate after each iteration
    individual_rates = np.array([40.6845,                           # r1- rate of reaction 1
                                 0.00626 * k[0],                    # r2- rate of reaction 2
                                 0.0058 * k[0],                     # r3- rate of reaction 3
                                 0.0036 * k[1],                     # r4- rate of reaction 4
                                 0.0068 * k[2],                     # r5- rate of reaction 5
                                 0.0000143 * k[2] * k[3],           # r6- rate of reaction 6
                                 0.0068 * k[4],                     # r7- rate of reaction 7
                                 0.009 * k[5],                      # r8- rate of reaction 8
                                 0.25 * k[2] * k[6],                # r9- rate of reaction 9
                                 40.6485,                           # r10- rate of reaction 10
                                 0.0058 * k[7],                     # r11- rate of reaction 11
                                 0.0001072 * k[6] * k[7],           # r12- rate of reaction 12
                                 0.004 * k[8],                      # r13- rate of reaction 13
                                 0.0019 * k[8],                     # r14- rate of reaction 14
                                 45.165,                            # r15- rate of reaction 15
                                 0.00008701 * k[8] * k[9],          # r16- rate of reaction 16
                                 0.0029 * k[9],                     # r17- rate of reaction 17
                                 0.0173 * k[10],                    # r18- rate of reaction 18
                                 0.0015 * k[9],                     # r19- rate of reaction 19
                                 0.005 * k[11],                     # r20- rate of reaction 20
                                 0.05 * (k[10] / k[11])])           # mrna upreg


plt.plot(time, gliA)

plt.xlabel("time (minutes)")
plt.ylabel("number of molecules")
plt.title("GliA profile in response to Shh signalling")
plt.grid()
plt.legend(loc=2)
plt.show()

