import numpy as np
import matplotlib.pyplot as plt
import csv


csv_file_path = '............shh_data_for_figure3.csv'            # file path for the data


ln_co = []
t = []

with open(csv_file_path, 'r') as csv_file:
    csv_reader = csv.DictReader(csv_file)

    for row in csv_reader:
        L_0 = 88.1
        tau = 41.6
        t_value = tau * np.log((float(row['Shh_size'])) / L_0)
        t.append(t_value)

        ln_co_value = np.log(float(row['Shh_CO']))
        ln_co.append(ln_co_value)


plt.scatter(t, ln_co,  marker='o', color='blue', alpha=0.7)

coefficients = np.polyfit(t, ln_co, 1)

slope = coefficients[0]
intercept = coefficients[1]


print("Intercept:", intercept)
print("Slope:", slope)


line_of_best_fit = np.polyval(coefficients, t)
plt.plot(t, line_of_best_fit, color='red')

plt.xlabel("time (hours)")
plt.ylabel("ln C_0 (arbitrary units)")
plt.title("ln C_0 vs time")
plt.grid()
plt.legend(loc=2)
plt.show()


